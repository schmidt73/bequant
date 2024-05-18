import argparse
import random
import sys
import os
import json

import pandas as pd
import numpy as np

from abc import abstractmethod
from loguru import logger
from dataclasses import dataclass
from collections import defaultdict

class Matcher:
    """
    A matcher object provides functionality to take a  (potentially 
    edited) sgRNA and match it against  a whitelist of sgRNAs. The 
    matcher object provides the `match` method which returns the ID 
    of the best  matching sgRNA from the whitelist, if any match 
    exists.
    """

    @abstractmethod
    def match(self, sgRNA):
        pass

class ExactMatcher(Matcher):
    """
    An exact matcher, which requires the sgRNA exactly matches a
    sgRNA in the whitelist.
    """

    def __init__(self, whitelist):
        self.whitelist = whitelist
        self.whitelist_dictionary = {}
        for guide in self.whitelist.values():
            self.whitelist_dictionary[guide['sgRNA']] = guide

    def match(self, sgRNA):
        if sgRNA in self.whitelist_dictionary:
            return self.whitelist_dictionary[sgRNA]['guide_ID']

        return None

class ProbabilisticMatcher(Matcher):
    """
    A probabilistic matching algorithm, parameterized by three
    variables:
    - n: the length of the sequence to be matched against
    - m <= n: the number of positions to select
    - k: the number of random hash functions to construct
    The complexity of the matcher grows *linearly* in k, the
    number of hash functions.
    """

    def __init__(self, whitelist, seed, n, m, k):
        self.whitelist = whitelist
        self.n = n
        self.m = m
        self.k = k

        indices = list(range(n))

        random.seed(seed)
        self.random_hash_functions = []
        for _ in range(k):
            positions = random.sample(indices, m) # randomly select m of n positions to use for hashing
            self.random_hash_functions.append(positions)

        self.tables = []
        for idx in range(len(self.random_hash_functions)):
            table = defaultdict(list)
            for guide in self.whitelist.values():
                key = self.random_hash(idx, guide['sgRNA'])
                table[key].append(guide)
            self.tables.append(table)

    def random_hash(self, idx, s):
        positions = self.random_hash_functions[idx]
        return hash(''.join(s[p] for p in positions))

    def match(self, sgRNA):
        if len(sgRNA) != self.n:
            return None

        def hd(s1, s2):
            return sum(s1[i] != s2[i] for i in range(len(s1)))

        candidates = []
        for i in range(len(self.tables)):
            key = self.random_hash(i, sgRNA)
            if key in self.tables[i]:
                candidates += self.tables[i][key]

        if not candidates:
            return None

        closest_candidate, dist = None, -1
        for candidate in candidates:
            if closest_candidate is None:
                dist = hd(sgRNA, candidate['sgRNA']) 
                closest_candidate = candidate
                continue

            candidate_dist = hd(sgRNA, candidate['sgRNA']) 
            if candidate_dist < dist:
                dist = candidate_dist
                closest_candidate = candidate

        return closest_candidate['guide_ID']

def parse_read(sensor_structure, read):
    """
    Parses a read, extracting the sgRNA and the genomic sensor 
    within it, returning None if the read is malformed.
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

    genomic_sensor_start += len(sensor_structure['sensor-left-scaffold'])
    genomic_sensor_start += sensor_structure['sensor-left-flank-length']
    genomic_sensor_end -= sensor_structure['sensor-right-flank-length']

    return {
        "sgrna": read[sgrna_start:sgrna_end],
        "genomic_sensor": read[genomic_sensor_start:genomic_sensor_end]
    }

def match_fastq(sensor_structure, matcher, fastq_file, progress=10000):
    """
    Matches all reads from a FASTq file with the matcher,
    using the specified sensor structure configuration.

    Returns a list mapping each read to its corresponding 
    guide, if such a guide exists.
    """

    matched_reads = []
    with open(fastq_file, 'r') as fastq:
        counter = 0
        while True:
            if counter % progress == 0:
                logger.info(f"Processed {counter} reads from {fastq_file}.")

            try: 
                header, raw_read, sense, qc = next(fastq), next(fastq), next(fastq), next(fastq) # process 4 lines at a time
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

def quantify_guide_editing(guide, matched_reads):
    """
    Quantifies the editing outcomes for all reads
    matching the guide. Returns a dictionary counting
    the different type of edits.
    """

    def hd(s1, s2):
        return sum(s1[i] != s2[i] for i in range(len(s1)))

    sgrna = guide['sgRNA']

    outcomes = {
        'non_edits': 0,
        'indels': 0,
        'edits': defaultdict(int)
    }

    for read in matched_reads:
        genomic_sensor = read['genomic_sensor']
        if len(genomic_sensor) != len(sgrna):
            outcomes['indels'] += 1
            continue

        hamming_distance = hd(sgrna, genomic_sensor)
        if hamming_distance == 0:
            outcomes['non_edits'] += 1
            continue

        for i in range(len(sgrna)):
            if sgrna[i] == genomic_sensor[i]: continue
            edit = (sgrna[i], genomic_sensor[i], hamming_distance, i) # (from, to, hamming-distance, position)
            outcomes['edits'][edit] += 1

    outcomes['total_edits'] = len(matched_reads) - outcomes['non_edits']
    return outcomes

def postprocess_outcomes(whitelist, guides_to_outcomes, base, targeted=True):
    """
    Post-processes the outcomes, formatting them into an easy to 
    understand Pandas dataframe.
    """

    columns = [
        "guide_ID", "sequence", "PAM", "exPAM", "edit_position", 
        "surrounding nucleotide context (NNCNN)", "total",
        f"t{base}AN", f"t{base}CN",  f"t{base}GN", f"t{base}TN", 
        f"t{base}A",  f"t{base}C",  f"t{base}G",  f"t{base}T",  
        "tINDEL"
    ]

    def count_edit(edits, src, dst, target, isolated):
        count = 0
        for edit in edits:
            if src != edit[0] or dst != edit[1] or (target is not None and target != edit[3]):
                continue
            if not isolated:
                count += edits[edit]
                continue
            if edit[2] == 1:
                count += edits[edit]
                continue
        return count

    rows = []
    for guide_id, outcomes in guides_to_outcomes.items():
        guide = whitelist[guide_id]
        target   = guide['target_position'] if targeted else None
        sequence = guide['sgRNA']
        surrounding_context = guide['sgRNA'][target - 2:target + 3] if target else None

        edits = outcomes['edits']
        tAN = count_edit(edits, base, 'A', target, False) if base != 'A' else None
        tCN = count_edit(edits, base, 'C', target, False) if base != 'C' else None
        tGN = count_edit(edits, base, 'G', target, False) if base != 'G' else None
        tTN = count_edit(edits, base, 'T', target, False) if base != 'T' else None
        tA  = count_edit(edits, base, 'A', target, True)  if base != 'A' else None
        tC  = count_edit(edits, base, 'C', target, True)  if base != 'C' else None
        tG  = count_edit(edits, base, 'G', target, True)  if base != 'G' else None
        tT  = count_edit(edits, base, 'T', target, True)  if base != 'T' else None
        tINDEL = outcomes['indels']
        total = outcomes['total_edits']

        row =  [
            guide_id, sequence, guide['PAM'], None, target,
            surrounding_context, total, tAN, tCN, tGN, tTN, 
            tA, tC, tG, tT, tINDEL
        ]

        rows.append(dict(zip(columns, row)))

    return pd.DataFrame(rows)

def process_whitelist(whitelist_file):
    """
    Reads the whitelist and returns a dictionary mapping guide IDs
    to the guide, where each guide is represented as a dictionary.
    """

    whitelist = pd.read_csv(whitelist_file)
    whitelist =  whitelist[['guide_ID', 'sgRNA', 'target_position', 'PAM',]]
    if not all(~whitelist['guide_ID'].duplicated()):
        raise Exception("Not all guide_IDs in whitelist are unique.")

    whitelist = whitelist.to_dict(orient='records')
    whitelist_dict = {}
    for guide in whitelist:
        whitelist_dict[guide['guide_ID']] = guide

    return whitelist_dict

def process_sensor_structure(sensor_structure_file):
    """
    Reads the sensor structure file and returns a dictionary 
    containing the sensor structure.
    """

    with open(sensor_structure_file, 'r') as f:
        sensor_structure = json.load(f)

    return sensor_structure

def parse_args():
    parser = argparse.ArgumentParser(
        description='Quantifies the editing outcomes from a base editing sensor screen.'
    )

    parser.add_argument(
        'whitelist', type=str,
        help='Path to the whitelist (CSV) containing guides to quantify.'
    )

    parser.add_argument(
        'sensor_structure', type=str,
        help='Path to JSON file describing sensor structure.'
    )

    parser.add_argument(
        'fastq', type=str, nargs='+',
        help='Path to FASTq files containing the sequencing reads.'
    )

    parser.add_argument(
        '-o', '--output', type=str, default='analysis.csv',
        help='Path to output file.'
    )

    parser.add_argument(
        '-f', '--format', type=str, choices=['target', 'all'], default='target',
        help='Path to fastq file containing the sequencing reads.'
    )

    parser.add_argument(
        '-b', '--base', type=str, default='A', choices=['A', 'C', 'G', 'T'],
        help='Base to quantify.'
    )

    parser.add_argument(
        '-m', '--matcher', type=str, choices=['exact', 'probabilistic'], default='exact',
        help='Type of matching algorithm.'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    """ Pre-processing whitelist and configuration. """
    logger.info(f"Reading sensor structure from {args.sensor_structure}")
    sensor_structure = process_sensor_structure(args.sensor_structure)
    sensor_structure['sensor-left-flank-length'] = 10
    sensor_structure['sensor-right-flank-length'] = 13

    logger.info(f"Reading whitelist from {args.whitelist}")
    whitelist = process_whitelist(args.whitelist)

    """ Constructing sgRNA matcher object. """
    logger.info(f"Constructing sgRNA matcher from whitelist.")
    if args.matcher == 'exact':
        matcher = ExactMatcher(whitelist)
    else:
        matcher = ProbabilisticMatcher(whitelist, 73, 20, 10, 10)

    """ Matching reads to guides in whitelist. """
    matched_reads = []
    for fastq_file in args.fastq:
        logger.info(f"Matching reads in {fastq_file} to guides in whitelist...")
        matches = match_fastq(sensor_structure, matcher, fastq_file, progress=50000)
        matched_reads += matches

    guides_to_reads = defaultdict(list)
    for matched_read in matched_reads:
        if matched_read['match'] is None: 
            continue
        guides_to_reads[matched_read['match']].append(matched_read['read'])

    """ Quantifying editing outcomes. """
    logger.info(f"Quantifying editing outcomes...")
    guides_to_outcomes = {}
    for guide, matched_reads in guides_to_reads.items():
        editing_outcomes = quantify_guide_editing(whitelist[guide], matched_reads)
        guides_to_outcomes[guide] = editing_outcomes

    """ Formatting editing outcomes. """
    logger.info(f"Formatting editing outcomes...")
    if args.format == 'target':
        result_df = postprocess_outcomes(whitelist, guides_to_outcomes, args.base, True)
        print(result_df)
        result_df.to_csv(args.output)
    else:
        result_df = postprocess_outcomes(whitelist, guides_to_outcomes, args.base, False)
        print(result_df)
        result_df.to_csv(args.output)

    logger.info(f"Processing complete.")
