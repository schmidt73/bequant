from Bio import SeqIO
import argparse
import csv
from tqdm import tqdm
import random
import json
import pandas as pd
import numpy as np
from abc import abstractmethod
from loguru import logger
from collections import defaultdict
import os
import sys

def filter_fastq(input_fastq, whitelist_file, output_fastq, sensor_structure, filter_length, column_name='sensor_amplicon', debug=False):
    """
    Filter sequences from the input FASTQ file based on the list of sequences in the whitelist.
    Only sequences that contain any of the filter sequences as substrings are retained.
    """

    # Scaffold to use for finding matching sgRNA to sensor
    scaffold = sensor_structure['sgrna-right-scaffold']

    filter_sequences = set()

    # Read sequences from the whitelist CSV file
    with open(whitelist_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if column_name in row:
                full_sequence = row[column_name].strip()
                if not full_sequence:
                    # Skip empty sequences
                    continue

                idx = full_sequence.find(scaffold)
                if idx == -1:
                    # Skip sequences that do not contain the sgRNA scaffold
                    continue

                filter_sequence = full_sequence[idx - filter_length:idx+len(scaffold) + filter_length]
                filter_sequences.add(filter_sequence)
            else:
                print(f"Column '{column_name}' not found in row: {row}")

    # Count total records for progress bar
    total_records = sum(1 for _ in SeqIO.parse(input_fastq, 'fastq'))

    # Filter the FASTQ file
    with open(output_fastq, 'w') as out_fh, tqdm(total=total_records, desc="Filtering FASTQ", unit="read") as pbar:
        for record in SeqIO.parse(input_fastq, 'fastq'):
            idx = record.seq.find(scaffold)
            if idx == -1:
                continue

            filter_match_sequence = record.seq[idx - filter_length:idx+len(scaffold) + filter_length]
            if filter_match_sequence in filter_sequences:
                SeqIO.write(record, out_fh, 'fastq')

            pbar.update(1)

class Matcher:
    @abstractmethod
    def match(self, sgRNA):
        pass

class ExactMatcher(Matcher):
    """
    Matcher that performs exact matching of sgRNA sequences against a whitelist.
    """
    def __init__(self, whitelist):
        self.whitelist = whitelist
        self.whitelist_dictionary = {guide['sgRNA']: guide for guide in self.whitelist.values()}

    def match(self, sgRNA):
        return self.whitelist_dictionary.get(sgRNA, {}).get('guide_ID')

class ProbabilisticMatcher(Matcher):
    """
    Matcher that performs probabilistic matching of sgRNA sequences using random hash functions.
    """
    def __init__(self, whitelist, seed, n, m, k):
        self.whitelist = whitelist
        self.n = n
        self.m = m
        self.k = k
        indices = list(range(n))
        random.seed(seed)
        self.random_hash_functions = [random.sample(indices, m) for _ in range(k)]
        self.tables = [defaultdict(list) for _ in range(len(self.random_hash_functions))]
        for idx, positions in enumerate(self.random_hash_functions):
            for guide in self.whitelist.values():
                key = self.random_hash(idx, guide['sgRNA'])
                self.tables[idx][key].append(guide)

    def random_hash(self, idx, s):
        return hash(''.join(s[p] for p in self.random_hash_functions[idx]))

    def match(self, sgRNA):
        if len(sgRNA) != self.n:
            return None

        def hd(s1, s2):
            return sum(s1[i] != s2[i] for i in range(len(s1)))

        candidates = []
        for i in range(len(self.tables)):
            key = self.random_hash(i, sgRNA)
            candidates.extend(self.tables[i].get(key, []))

        if not candidates:
            return None

        closest_candidate, dist = min(((c, hd(sgRNA, c['sgRNA'])) for c in candidates), key=lambda x: x[1])
        return closest_candidate['guide_ID']

def parse_read(sensor_structure, read):
    """
    Parse a read to extract the sgRNA and genomic sensor sequences based on the sensor structure.
    """
    sgrna_start = read.find(sensor_structure['sgrna-left-scaffold'])
    if sgrna_start == -1:
        return None
    sgrna_start += len(sensor_structure['sgrna-left-scaffold'])
    sgrna_end = read.find(sensor_structure['sgrna-right-scaffold'])
    if sgrna_end == -1:
        return None

    genomic_sensor_start = read.find(sensor_structure['sensor-left-scaffold'])
    if genomic_sensor_start == -1:
        return None
    genomic_sensor_end = read.find(sensor_structure['sensor-right-scaffold'])
    if genomic_sensor_end == -1:
        return None

    genomic_sensor_start += len(sensor_structure['sensor-left-scaffold']) + sensor_structure['sensor-left-flank-length']
    genomic_sensor_end -= sensor_structure['sensor-right-flank-length']

    return {
        "sgrna": read[sgrna_start:sgrna_end],
        "genomic_sensor": read[genomic_sensor_start:genomic_sensor_end]
    }

def extract_sensor_sequence(sensor_structure, sensor_amplicon):
    """
    Extract the sensor sequence from the sensor amplicon based on the sensor structure.
    """
    genomic_sensor_start = sensor_amplicon.find(sensor_structure['sensor-left-scaffold'])
    if genomic_sensor_start == -1:
        return None
    genomic_sensor_end = sensor_amplicon.find(sensor_structure['sensor-right-scaffold'])
    if genomic_sensor_end == -1:
        return None

    genomic_sensor_start += len(sensor_structure['sensor-left-scaffold']) + sensor_structure['sensor-left-flank-length']
    genomic_sensor_end -= sensor_structure['sensor-right-flank-length']

    return sensor_amplicon[genomic_sensor_start:genomic_sensor_end]

def match_fastq(sensor_structure, matcher, fastq_file, progress=10000):
    """
    Match reads from the FASTQ file to guides in the whitelist and extract relevant sequences.
    """
    matched_reads = []
    with open(fastq_file, 'r') as fastq:
        counter = 0
        while True:
            if counter % progress == 0:
                logger.info(f"Processed {counter} reads from {fastq_file}.")
            try:
                header, raw_read, sense, qc = next(fastq), next(fastq), next(fastq), next(fastq)
                read = parse_read(sensor_structure, raw_read)
                if read is None:
                    continue
                match = matcher.match(read['sgrna'])
                matched_reads.append({
                    "read": read,
                    "match": match
                })
            except StopIteration:
                break
            counter += 1
    return matched_reads

def quantify_guide_editing(guide, matched_reads, sensor_structure, base, format_type):
    """
    Quantify the editing outcomes for a guide based on the matched reads.
    """
    def hd(s1, s2):
        return sum(s1[i] != s2[i] for i in range(len(s1)))

    sensor_amplicon = extract_sensor_sequence(sensor_structure, guide['sensor_amplicon'])
    if sensor_amplicon is None:
        raise ValueError("Unable to extract sensor sequence from sensor_amplicon.")

    target_position = guide['target_position'] if format_type == 'target' else None

    outcomes = defaultdict(lambda: defaultdict(int))
    non_edit_count = 0
    indel_count = 0

    for read in matched_reads:
        genomic_sensor = read['genomic_sensor']
        if len(genomic_sensor) != len(sensor_amplicon):
            indel_count += 1
            continue

        hamming_distance = hd(sensor_amplicon, genomic_sensor)
        if hamming_distance == 0:
            non_edit_count += 1
            continue

        single_edit = hamming_distance == 1

        for i in range(len(sensor_amplicon)):
            if format_type == 'target' and i != target_position:
                continue

            if sensor_amplicon[i] == base:
                if genomic_sensor[i] == 'A':
                    outcomes[i]['A'] += 1
                    if single_edit:
                        outcomes[i]['A_single'] += 1
                elif genomic_sensor[i] == 'C':
                    outcomes[i]['C'] += 1
                    if single_edit:
                        outcomes[i]['C_single'] += 1
                elif genomic_sensor[i] == 'G':
                    outcomes[i]['G'] += 1
                    if single_edit:
                        outcomes[i]['G_single'] += 1
                elif genomic_sensor[i] == 'T':
                    outcomes[i]['T'] += 1
                    if single_edit:
                        outcomes[i]['T_single'] += 1

    total_reads = len(matched_reads)
    total_minus_indels = total_reads - indel_count
    outcomes['total_reads'] = total_reads
    outcomes['total_minus_indels'] = total_minus_indels
    outcomes['non_edits'] = non_edit_count
    outcomes['indels'] = indel_count

    return outcomes

def postprocess_outcomes(whitelist, guides_to_outcomes, sensor_structure, base, format_type):
    """
    Post-process the quantified editing outcomes and format them into a DataFrame.
    """
    columns = [
        "guide_ID", "sequence", "PAM", "edit_position", 
        "total_reads", "total_minus_indels", "non_edits", "indels"
    ]

    rows = []
    for guide_id, outcomes in guides_to_outcomes.items():
        guide = whitelist[guide_id]
        sequence = guide['sgRNA']
        pam = guide['PAM']
        total_reads = outcomes['total_reads']
        total_minus_indels = outcomes['total_minus_indels']
        non_edits = outcomes['non_edits']
        indels = outcomes['indels']

        sensor_amplicon = extract_sensor_sequence(sensor_structure, guide['sensor_amplicon'])
        if sensor_amplicon is None:
            raise ValueError("Unable to extract sensor sequence from sensor_amplicon.")

        for position in range(len(sensor_amplicon)):
            if format_type == 'target' and position != guide['target_position']:
                continue

            if sensor_amplicon[position] == base:
                row = {
                    "guide_ID": guide_id,
                    "sequence": sequence,
                    "PAM": pam,
                    "edit_position": position,
                    "total_reads": total_reads,
                    "total_minus_indels": total_minus_indels,
                    "non_edits": non_edits,
                    "indels": indels,
                    f"t{base}AN": outcomes[position]['A'],
                    f"t{base}CN": outcomes[position]['C'],
                    f"t{base}GN": outcomes[position]['G'],
                    f"t{base}TN": outcomes[position]['T'],
                    f"t{base}A_single": outcomes[position]['A_single'],
                    f"t{base}C_single": outcomes[position]['C_single'],
                    f"t{base}G_single": outcomes[position]['G_single'],
                    f"t{base}T_single": outcomes[position]['T_single']
                }
                rows.append(row)

    df = pd.DataFrame(rows)
    columns_to_drop = [col for col in df.columns if any(col.startswith(f"t{base}{base}N") or col.startswith(f"t{base}{base}_single") for base in "ACGT")]
    filtered_columns = [col for col in df.columns if col not in columns_to_drop]
    df = df[filtered_columns]

    return df

def process_whitelist(whitelist_file):
    """
    Process the whitelist CSV file and return it as a dictionary.
    """
    whitelist = pd.read_csv(whitelist_file)
    whitelist = whitelist[['guide_ID', 'sgRNA', 'target_position', 'PAM', 'sensor_amplicon']]
    if not all(~whitelist['guide_ID'].duplicated()):
        raise Exception("Not all guide_IDs in whitelist are unique.")
    whitelist = whitelist.to_dict(orient='records')
    whitelist_dict = {guide['guide_ID']: guide for guide in whitelist}
    return whitelist_dict

def process_sensor_structure(sensor_structure_file):
    """
    Process the sensor structure JSON file and return it as a dictionary.
    """
    with open(sensor_structure_file, 'r') as f:
        sensor_structure = json.load(f)
    return sensor_structure

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Quantifies the editing outcomes from a base editing sensor screen.'
    )

    parser.add_argument(
        '-i', '--input_fastq', type=str, required=True,
        help='Path to the input FASTQ file containing the sequencing reads.'
    )
    parser.add_argument(
        '-w', '--whitelist', type=str, required=True,
        help='CSV file containing the list of sequences to keep and guides to quantify.'
    )
    parser.add_argument(
        '-o', '--output_fastq', type=str, required=True,
        help='Path to the output filtered FASTQ file.'
    )
    parser.add_argument(
        '-s', '--sensor_structure', type=str, required=True,
        help='Path to JSON file describing sensor structure.'
    )
    parser.add_argument(
        '-q', '--quant_output', type=str, default='analysis.csv',
        help='Path to output file for quantified editing outcomes (default: analysis.csv).'
    )
    parser.add_argument(
        '-c', '--column_name', type=str, default='sensor_amplicon',
        help='Column name in the CSV file containing sequences (default: sensor_amplicon).'
    )
    parser.add_argument(
        '--filter_length', type=int, default=5,
        help='Length of the filter around the scaffold (default: 5).'
    )
    parser.add_argument(
        '-f', '--format', type=str, choices=['target', 'all'], default='target',
        help='Format for the output: "target" for specific target positions, "all" for all positions.'
    )
    parser.add_argument(
        '-b', '--base', type=str, default='A', choices=['A', 'C', 'G', 'T'],
        help='Base to quantify.'
    )
    parser.add_argument(
        '-m', '--matcher', type=str, choices=['exact', 'probabilistic'], default='probabilistic',
        help='Type of matching algorithm.'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true', help='Enable debugging output'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Check if input FASTQ file exists
    if not os.path.isfile(args.input_fastq):
        print(f"Error: Input FASTQ file '{args.input_fastq}' not found.")
        sys.exit(1)

    # Step 1: Read sensor structure and whitelist
    logger.info(f"Reading sensor structure from {args.sensor_structure}")
    sensor_structure = process_sensor_structure(args.sensor_structure)
    sensor_structure['sensor-left-flank-length'] = 10
    sensor_structure['sensor-right-flank-length'] = 13

    # Step 2: Filter the FASTQ file
    logger.info(f"Filtering FASTQ file {args.input_fastq} using sequences from {args.whitelist}")
    filter_fastq(args.input_fastq, args.whitelist, args.output_fastq, sensor_structure, args.filter_length, args.column_name, args.debug)

    logger.info(f"Reading whitelist from {args.whitelist}")
    whitelist = process_whitelist(args.whitelist)

    logger.info(f"Constructing sgRNA matcher from whitelist.")
    if args.matcher == 'exact':
        matcher = ExactMatcher(whitelist)
    else:
        matcher = ProbabilisticMatcher(whitelist, 73, 20, 7, 20)

    # Step 3: Match reads to guides in the filtered FASTQ file
    logger.info(f"Matching reads in {args.output_fastq} to guides in whitelist...")
    matched_reads = match_fastq(sensor_structure, matcher, args.output_fastq, progress=50000)

    guides_to_reads = defaultdict(list)
    for matched_read in matched_reads:
        if matched_read['match'] is None:
            continue
        guides_to_reads[matched_read['match']].append(matched_read['read'])

    # Step 4: Quantify editing outcomes
    logger.info(f"Quantifying editing outcomes...")
    guides_to_outcomes = {}
    for guide, matched_reads in guides_to_reads.items():
        editing_outcomes = quantify_guide_editing(whitelist[guide], matched_reads, sensor_structure, args.base, args.format)
        guides_to_outcomes[guide] = editing_outcomes

    # Step 5: Postprocess and save the outcomes
    logger.info(f"Formatting editing outcomes...")
    result_df = postprocess_outcomes(whitelist, guides_to_outcomes, sensor_structure, args.base, args.format)

    print(result_df)
    result_df.to_csv(args.quant_output, index=False)

    logger.info(f"Processing complete.")

