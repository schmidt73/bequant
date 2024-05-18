import argparse
import sys
import os
import json

import pandas as pd
import numpy as np

from abc import abstractmethod
from loguru import logger
from dataclasses import dataclass

class Matcher:
    """
    A matcher object provides functionality to take 
    a (potentially edited) sgRNA and match it against 
    a whitelist of sgRNAs. The matcher object provides
    the `match` method which returns the best matching
    sgRNA from the whitelist, if any.
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
        for guide in self.whitelist:
            self.whitelist_dictionary[guide['sgRNA']] = guide

    def match(self, sgRNA):
        if sgRNA in self.whitelist_dictionary:
            return self.whitelist_dictionary[sgRNA]

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


def process_whitelist(whitelist_file):
    """
    Reads the whitelist and returns a list of guides to quantify,
    each guide is represented as a dictionary.
    """

    whitelist = pd.read_csv(whitelist_file)
    whitelist =  whitelist[['guide_ID', 'sgRNA', 'target_position', 'PAM',]]
    whitelist = whitelist.to_dict(orient='records')
    return whitelist

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

    logger.info(f"Reading sensor structure from {args.sensor_structure}")
    sensor_structure = process_sensor_structure(args.sensor_structure)
    sensor_structure['sensor-left-flank-length'] = 10
    sensor_structure['sensor-right-flank-length'] = 13

    logger.info(f"Reading whitelist from {args.whitelist}")
    whitelist = process_whitelist(args.whitelist)

    logger.info(f"Constructing sgRNA matcher from whitelist.")
    matcher = ExactMatcher(whitelist)

    matched_reads = []
    for fastq_file in args.fastq:
        logger.info(f"Matching reads in {fastq_file} to guides in whitelist...")
        matches = match_fastq(sensor_structure, matcher, fastq_file)
        matched_reads.append(matches)

    print(matched_reads)
    
    

