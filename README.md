# bequant

`bequant` processes raw data from a sensor based base editing
screen into a mapping from sgRNAs to editing outcomes. 
`bequant` uses raw, paired-end read data to classify the editing
events that resulted in the observed target sites. Complete details 
describing the method can be found at:
```
```

## usage

To use `bequant` simply execute:
```
$ python bequant.py --help
Usage: java -jar analyze-fastqs rep1 ... repN [options]

Options:
  -o, --output FILE                   Output file (REQUIRED).
  -w, --whitelist WHITELIST           Whitelist file for screen (REQUIRED).
  -c, --config CONFIG                 Configuration file describing sensor and whitelist structure (REQUIRED).
  -f, --format OUTPUT_FORMAT  TARGET  Output either only target cytosine or all cytosines.
  -h, --help

Expects FASTq files consisting of the paired end reads
for each replicate.
```

## example

As an example, we will use `bequant` to process a subset of data
in an mKRAS tiling screen of cells treated with Adagrasib. To apply
`bequant` to this example, run the following command:

```
$ python bequant.py -w example/mKras_whitelist.csv -i example/R1_D10_MRTX849_full.fastq \
                    -s example/mKras_scaffold.json -o example/R1_D10_MRTX849_full.processed.fast \
                    -q example/R1_D10_MRTX849_full_quantification.csv
```

Output files...
