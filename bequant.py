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
import gzip


def parse_args():
    parser = argparse.ArgumentParser(description='Quantifies the editing outcomes from a base editing sensor screen.')

    parser.add_argument('-i', '--input_fastq', type=str, required=True)
    parser.add_argument('-w', '--whitelist', type=str, required=True)
    parser.add_argument('-o', '--output_fastq', type=str, required=True)
    parser.add_argument('-s', '--sensor_structure', type=str, required=True)
    parser.add_argument('-q', '--quant_output', type=str, default='analysis.csv')
    parser.add_argument('-c', '--column_name', type=str, default='sensor_amplicon')
    parser.add_argument('--filter_length', type=int, default=5)
    parser.add_argument('-f', '--format', type=str, choices=['target', 'all'], default='target')
    parser.add_argument('-b', '--base', type=str, default='A', choices=['A', 'C', 'G', 'T'])
    parser.add_argument('--matcher', type=str, choices=['exact', 'probabilistic'], default='probabilistic')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--self_editing_analysis', action='store_true', help='Enable sgRNA editing quantification')
    parser.add_argument('--run_hd_summary', action='store_true', help='Enable hamming distance summary output')
    parser.add_argument('--run_mismatch_analysis', action='store_true', help='Enable mismatch scaffold analysis')
    parser.add_argument('--hamming_threshold', type=int, default=9, help='Max Hamming distance to include reads (default: 9)')
    parser.add_argument('--m', type=int, default=7, help='Number of positions per hash (probabilistic)')
    parser.add_argument('--k', type=int, default=20, help='Number of hash tables (probabilistic)')
    

    return parser.parse_args()
    
def filter_fastq(input_fastq, whitelist_file, output_fastq, sensor_structure, filter_length, column_name='sensor_amplicon', debug=False, mismatched_output_fastq=None, qc_stats=None):
    """
    Filter sequences from the input FASTQ file based on the list of sequences in the whitelist.
    Also track and optionally save reads with mismatched scaffolds.
    """
    scaffold = sensor_structure['sgrna-right-scaffold']
    filter_sequences = set()

    with open(whitelist_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if column_name in row:
                full_sequence = row[column_name].strip()
                if not full_sequence:
                    continue

                idx = full_sequence.find(scaffold)
                if idx == -1:
                    continue

                filter_sequence = full_sequence[idx - filter_length:idx + len(scaffold) + filter_length]
                filter_sequences.add(filter_sequence)
            else:
                print(f"Column '{column_name}' not found in row: {row}")

    if input_fastq.endswith('.gz'):
        handle = gzip.open(input_fastq, 'rt')
    else:
        handle = open(input_fastq, 'r')

    total_records = sum(1 for _ in SeqIO.parse(handle, 'fastq'))
    handle.close()

    if input_fastq.endswith('.gz'):
        handle = gzip.open(input_fastq, 'rt')
    else:
        handle = open(input_fastq, 'r')

    with handle, \
         open(output_fastq, 'w') as out_fh, \
         open(mismatched_output_fastq, 'w') if mismatched_output_fastq else open(os.devnull, 'w') as mismatch_fh, \
         tqdm(total=total_records, desc="Filtering FASTQ", unit="read") as pbar:

        for record in SeqIO.parse(handle, 'fastq'):
            idx = record.seq.find(scaffold)
            if idx == -1:
                if qc_stats is not None:
                    qc_stats["mutated_scaffold"] += 1
                SeqIO.write(record, mismatch_fh, 'fastq')
                pbar.update(1)
                continue

            found_seq = record.seq[idx:idx + len(scaffold)]
            if found_seq != scaffold:
                if qc_stats is not None:
                    qc_stats["mutated_scaffold"] += 1
                SeqIO.write(record, mismatch_fh, 'fastq')
                pbar.update(1)
                continue
            else:
                if qc_stats is not None:
                    qc_stats["perfect_scaffold"] += 1

            filter_match_sequence = record.seq[idx - filter_length:idx + len(scaffold) + filter_length]
            if filter_match_sequence in filter_sequences:
                SeqIO.write(record, out_fh, 'fastq')
                if qc_stats is not None:
                    qc_stats["analyzed_reads"] += 1

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

def parse_read(sensor_structure, read, qc_stats=None, allow_scaffold_mismatch=False):
    scaffold_seq = sensor_structure['sgrna-right-scaffold']

    sgrna_end = read.find(scaffold_seq)
    if sgrna_end == -1:
        if qc_stats is not None:
            qc_stats["mutated_scaffold"] += 1
        if not allow_scaffold_mismatch:
            return None

        sgrna_start = read.find(sensor_structure['sgrna-left-scaffold'])
        if sgrna_start == -1:
            return None
        sgrna_start += len(sensor_structure['sgrna-left-scaffold'])
        sgrna_end = sgrna_start + 20  # assume sgRNA is 20bp long when scaffold is missing
    else:
        sgrna_start = read.find(sensor_structure['sgrna-left-scaffold'])
        if sgrna_start == -1:
            return None
        sgrna_start += len(sensor_structure['sgrna-left-scaffold'])

        found_seq = read[sgrna_end : sgrna_end + len(scaffold_seq)]
        if found_seq != scaffold_seq:
            if not allow_scaffold_mismatch:
                return None

    sgrna = read[sgrna_start:sgrna_end]

    genomic_sensor_start = read.find(sensor_structure['sensor-left-scaffold'])
    genomic_sensor_end = read.find(sensor_structure['sensor-right-scaffold'])
    if genomic_sensor_start == -1 or genomic_sensor_end == -1:
        return None

    genomic_sensor_start += len(sensor_structure['sensor-left-scaffold']) + sensor_structure['sensor-left-flank-length']
    genomic_sensor_end -= sensor_structure['sensor-right-flank-length']

    return {
        "sgrna": sgrna,
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

def match_fastq(sensor_structure, matcher, fastq_file, read_qc_stats=None, progress=10000, allow_scaffold_mismatch=False):
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
                read = parse_read(sensor_structure, raw_read, qc_stats=read_qc_stats, allow_scaffold_mismatch=allow_scaffold_mismatch)
                if read is None:
                    continue  # ✅ skip reads that failed parsing
                    
                match = matcher.match(read['sgrna'])
                
                matched_reads.append({
                    "read": read,
                    "match": match
                })
            except StopIteration:
                break
            counter += 1
    return matched_reads

def process_mismatched_reads(sensor_structure, matcher, whitelist, args):
    if not args.run_mismatch_analysis:
        return

    logger.info("Processing mismatched scaffold reads for sensor-only editing...")

    mismatch_fastq = args.output_fastq.replace(".fastq", "_mismatched.fastq")
    if not os.path.exists(mismatch_fastq):
        logger.warning(f"No mismatched FASTQ found at {mismatch_fastq}. Skipping sensor-only analysis.")
        return

    logger.info(f"Found mismatched FASTQ: {mismatch_fastq}")

    mismatched_matched_reads = match_fastq(
        sensor_structure,
        matcher,
        mismatch_fastq,
        read_qc_stats=None,
        progress=50000,
        allow_scaffold_mismatch=True
    )

    mismatch_guides_to_reads = defaultdict(list)
    for matched_read in mismatched_matched_reads:
        if matched_read['match'] is None:
            continue
        mismatch_guides_to_reads[matched_read['match']].append(matched_read['read'])

    mismatch_guides_to_outcomes = {}
    mismatch_hd_summaries = {}
    for guide_id, matched_reads in mismatch_guides_to_reads.items():
        outcomes, hd_summary = quantify_guide_editing(
            whitelist[guide_id], matched_reads, sensor_structure, args.base, args.format, args.hamming_threshold
        )
        mismatch_guides_to_outcomes[guide_id] = outcomes
        mismatch_hd_summaries[guide_id] = hd_summary

    mismatch_result_df = postprocess_outcomes(
        whitelist, mismatch_guides_to_outcomes, sensor_structure, args.base, args.format
    )

    mismatch_output_file = args.quant_output.replace(".csv", "_sensor_only_mismatched.csv")
    mismatch_result_df.to_csv(mismatch_output_file, index=False)
    logger.info(f"Saved mismatch sensor-only analysis to {mismatch_output_file}")

    mismatch_hd_rows = []
    for guide_id, summary in mismatch_hd_summaries.items():
        for distance, (count, excluded) in sorted(summary.items()):
            mismatch_hd_rows.append({
                "guide_ID": guide_id,
                "hamming_distance": distance,
                "read_count": count,
                "excluded": excluded
            })

    mismatch_hd_df = pd.DataFrame(mismatch_hd_rows)
    mismatch_hd_summary_file = args.quant_output.replace(".csv", "_sensor_only_mismatched_hd_summary.csv")
    mismatch_hd_df.to_csv(mismatch_hd_summary_file, index=False)
    logger.info(f"Saved mismatch HD summary to {mismatch_hd_summary_file}")
    
def quantify_guide_editing(guide, matched_reads, sensor_structure, base, format_type, hamming_threshold=9):
    def hd(s1, s2):
        return sum(s1[i] != s2[i] for i in range(len(s1)))

    sensor_amplicon = extract_sensor_sequence(sensor_structure, guide['sensor_amplicon'])
    if sensor_amplicon is None:
        raise ValueError("Unable to extract sensor sequence from sensor_amplicon.")

    target_position = guide['target_position'] if format_type == 'target' else None

    outcomes = defaultdict(lambda: defaultdict(int))
    hd_summary = defaultdict(lambda: [0, False])
    non_edit_count = 0
    indel_count = 0
    above_hd_threshold = 0

    for read in matched_reads:
        genomic_sensor = read['genomic_sensor']
        if len(genomic_sensor) != len(sensor_amplicon):
            indel_count += 1
            continue

        hamming_distance = hd(sensor_amplicon, genomic_sensor)
        excluded = hamming_distance > hamming_threshold
        hd_summary[hamming_distance][0] += 1
        hd_summary[hamming_distance][1] = excluded

        if hamming_distance == 0:
            non_edit_count += 1
            continue

        if excluded:
            above_hd_threshold += 1
            continue

        single_edit = hamming_distance == 1

        for i in range(len(sensor_amplicon)):
            if format_type == 'target' and i != target_position:
                continue
            if sensor_amplicon[i] == base:
                if genomic_sensor[i] in 'ACGT':
                    outcomes[i][genomic_sensor[i]] += 1
                    if single_edit:
                        outcomes[i][f"{genomic_sensor[i]}_single"] += 1

    total_reads = len(matched_reads)
    total_minus_indels = total_reads - indel_count
    outcomes['total_reads'] = total_reads
    outcomes['total_minus_indels'] = total_minus_indels
    outcomes['non_edits'] = non_edit_count
    outcomes['indels'] = indel_count
    outcomes['above_hd_threshold'] = above_hd_threshold

    return outcomes, hd_summary

def quantify_sgrna_editing(guide, matched_reads):
    """
    Quantify the editing outcomes for the sgRNA region.
    """
    def hd(s1, s2):
        return sum(s1[i] != s2[i] for i in range(len(s1)))

    sgrna_reference = guide['sgRNA']  # Expected sgRNA sequence
    sgrna_length = len(sgrna_reference)

    outcomes = defaultdict(lambda: defaultdict(int))
    non_edit_count = 0
    indel_count = 0

    for read in matched_reads:
        sgrna_seq = read['sgrna']  # Extracted sgRNA from read
        
        if len(sgrna_seq) != sgrna_length:
            indel_count += 1
            continue

        hamming_distance = hd(sgrna_reference, sgrna_seq)
        if hamming_distance == 0:
            non_edit_count += 1
            continue

        single_edit = hamming_distance == 1
        
        for i in range(sgrna_length):
            if sgrna_reference[i] != sgrna_seq[i]:
                outcomes[i][sgrna_seq[i]] += 1  # Count observed base at mismatch
                if single_edit:
                    outcomes[i][f"{sgrna_seq[i]}_single"] += 1

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

def postprocess_sgrna_outcomes(whitelist, guides_to_sgrna_outcomes, target_base):
    """
    Post-process the quantified sgRNA editing outcomes into a DataFrame.
    Uses the simple calculation: total reads with an edit = (total reads - indels - non-edits).
    """
    columns = ["guide_ID", "sequence", "total_reads", "total_minus_indels", "non_edits", "indels", "self_edited_sgRNAs"]
    
    rows = []
    for guide_id, outcomes in guides_to_sgrna_outcomes.items():
        guide = whitelist[guide_id]
        sequence = guide["sgRNA"]
        total_reads = outcomes["total_reads"]
        total_minus_indels = outcomes["total_minus_indels"]
        non_edits = outcomes["non_edits"]
        indels = outcomes["indels"]

        # ✅ Simple subtraction approach: Count reads with at least one edit at the target base
        reads_with_target_base_edit = max(0, total_minus_indels - non_edits)  # Ensure no negatives

        row = {
            "guide_ID": guide_id,
            "sequence": sequence,
            "total_reads": total_reads,
            "total_minus_indels": total_minus_indels,
            "non_edits": non_edits,
            "indels": indels,
            "self_edited_sgRNAs": reads_with_target_base_edit
        }
        
        rows.append(row)
    
    return pd.DataFrame(rows)

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

def write_hd_summary(hd_summaries, output_path):
    rows = []
    for guide_id, summary in hd_summaries.items():
        for distance, (count, excluded) in sorted(summary.items()):
            rows.append({
                "guide_ID": guide_id,
                "hamming_distance": distance,
                "read_count": count,
                "excluded": excluded
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    logger.info(f"Saved Hamming distance summary to {output_path}")

def combine_results(sensor_df, sgrna_df):
    sgrna_subset = sgrna_df[["guide_ID", "self_edited_sgRNAs"]]
    return pd.merge(sensor_df, sgrna_subset, on="guide_ID", how="left")

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

    read_qc_stats = {
        "mutated_scaffold": 0,
        "perfect_scaffold": 0,
        "analyzed_reads": 0
    }
    
    # Step 2: Filter the FASTQ file
    logger.info(f"Filtering FASTQ file {args.input_fastq} using sequences from {args.whitelist}")
    filter_fastq(
        args.input_fastq,
        args.whitelist,
        args.output_fastq,
        sensor_structure,
        args.filter_length,
        args.column_name,
        args.debug,
        mismatched_output_fastq=args.output_fastq.replace('.fastq', '_mismatched.fastq'),
        qc_stats=read_qc_stats
    )
    
    logger.info(f"Reading whitelist from {args.whitelist}")
    whitelist = process_whitelist(args.whitelist)

    logger.info(f"Constructing sgRNA matcher from whitelist.")
    if args.matcher == 'exact':
        matcher = ExactMatcher(whitelist)
    else:
        matcher = ProbabilisticMatcher(whitelist, seed=73, n=20, m=args.m, k=args.k)
    
    # Step 3: Match reads to guides in the filtered FASTQ file
    logger.info(f"Matching reads in {args.output_fastq} to guides in whitelist...")
    matched_reads = match_fastq(sensor_structure, matcher, args.output_fastq, read_qc_stats=read_qc_stats, progress=50000)

    guides_to_reads = defaultdict(list)
    for matched_read in matched_reads:
        if matched_read['match'] is None:
            continue
        guides_to_reads[matched_read['match']].append(matched_read['read'])

    # Step 4: Quantify editing outcomes
    logger.info(f"Quantifying editing outcomes...")
    guides_to_outcomes = {}
    guide_hd_summaries = {}
    for guide, matched_reads in guides_to_reads.items():
        editing_outcomes, hd_summary = quantify_guide_editing(
            whitelist[guide], matched_reads, sensor_structure, args.base, args.format, args.hamming_threshold)
        guides_to_outcomes[guide] = editing_outcomes
        guide_hd_summaries[guide] = hd_summary
        
    # Quantifying sgRNA editing outcomes
    if args.self_editing_analysis:
        logger.info(f"Quantifying sgRNA editing outcomes...")
        guides_to_sgrna_outcomes = {}
        for guide, matched_reads in guides_to_reads.items():
            sgrna_editing_outcomes = quantify_sgrna_editing(whitelist[guide], matched_reads)
            guides_to_sgrna_outcomes[guide] = sgrna_editing_outcomes

        logger.info(f"Total guides processed for sgRNA editing: {len(guides_to_sgrna_outcomes)}")

    # Step 5: Postprocess and save the outcomes
    logger.info(f"Formatting editing outcomes...")
    result_df = postprocess_outcomes(whitelist, guides_to_outcomes, sensor_structure, args.base, args.format)

    # Optional: Save Hamming distance summary for each guide
    logger.info(f"Saving Hamming distance mismatch summary...")

    if args.run_hd_summary:
        hd_summary_file = args.quant_output.replace(".csv", "_hd_summary.csv")
        write_hd_summary(guide_hd_summaries, hd_summary_file)
    
    print(result_df)
    result_df.to_csv(args.quant_output, index=False)
    
    # Postprocess sgRNA editing outcomes
    if args.self_editing_analysis:
        logger.info(f"Processing sgRNA editing outcomes...")
        sgrna_result_df = postprocess_sgrna_outcomes(whitelist, guides_to_sgrna_outcomes, args.base)
        sgrna_output_file = args.quant_output.replace(".csv", "_self_editing.csv")
        sgrna_result_df.to_csv(sgrna_output_file, index=False)
        logger.info(f"Saved sgRNA editing analysis to {sgrna_output_file}")
    
        if sgrna_result_df is None:
            logger.error("Error: postprocess_sgrna_outcomes() returned None!")
        elif sgrna_result_df.empty:
            logger.warning("Warning: sgRNA result dataframe is empty! Check quantification logic.")
        else:
            logger.info(f"sgRNA result dataframe shape: {sgrna_result_df.shape}")

    # Merge sensor and sgRNA results
        logger.info("Merging sensor and sgRNA results...")
        final_result_df = combine_results(result_df, sgrna_result_df)

    # Check if merging was successful
        if final_result_df is None:
            logger.error("Error: combine_results() returned None!")
        elif final_result_df.empty:
            logger.warning("Warning: Combined analysis file is empty! Check guide_IDs.")
        else:
            logger.info(f"Final merged dataframe shape: {final_result_df.shape}")

    # Save the combined results
        combined_output_file = args.quant_output.replace(".csv", "_combined.csv")
        final_result_df.to_csv(combined_output_file, index=False)
        logger.info(f"Saved combined analysis to {combined_output_file}")
    
    process_mismatched_reads(sensor_structure, matcher, whitelist, args)

    print("\n🔍 QC Summary:")
    print(f"  Mutated scaffolds (missing or mismatched): {read_qc_stats['mutated_scaffold']:,}")
    perfect_scaffold_reads = sum(1 for _ in open(args.output_fastq)) // 4
    print(f"  Perfect scaffold reads: {perfect_scaffold_reads:,}")
    total_analyzed_reads = sum(out['total_reads'] for out in guides_to_outcomes.values())
    print(f"  Analyzed reads (retained for editing analysis): {total_analyzed_reads:,}")
    sensor_indel_count = sum(out['indels'] for out in guides_to_outcomes.values())
    above_hd_count = sum(out.get('above_hd_threshold', 0) for out in guides_to_outcomes.values())
    print(f"  Reads with sensor indels: {sensor_indel_count:,}")
    print(f"  Reads above Hamming threshold: {above_hd_count:,}")

    logger.info(f"Processing complete.")
