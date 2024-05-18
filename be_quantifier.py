import argparse
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

def match_fastq(sensor_structure, matcher, fastq_file):
    """
    Matches all reads from a FASTq file with the matcher,
    using the specified sensor structure configuration.

    Returns a list mapping each read to its corresponding 
    guide, if such a guide exists.
    """

    matched_reads = []
    with open(fastq_file, 'r') as fastq:
        while True:
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
        '-f', '--format', type=str, choices=['target', 'all', 'json'], default='target',
        help='Path to fastq file containing the sequencing reads.'
    )

    parser.add_argument(
        '-b', '--base', type=str, default='A', choices=['A', 'C', 'G', 'T'],
        help='Base to quantify.'
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
    matcher = ExactMatcher(whitelist)

    """ Matching reads to guides in whitelist. """
    matched_reads = []
    for fastq_file in args.fastq:
        logger.info(f"Matching reads in {fastq_file} to guides in whitelist...")
        matches = match_fastq(sensor_structure, matcher, fastq_file)
        matched_reads += matches

    guides_to_reads = defaultdict(list)
    for matched_read in matched_reads:
        if matched_read['match'] is None: 
            continue
        guides_to_reads[matched_read['match']].append(matched_read['read'])

    """ Quantifying editing outcomes. """
    guides_to_outcomes = {}
    for guide, matched_reads in guides_to_reads.items():
        editing_outcomes = quantify_guide_editing(whitelist[guide], matched_reads)
        guides_to_outcomes[guide] = editing_outcomes

    print(guides_to_outcomes)


