# Processing Sensor Outcomes

This script processes raw data from a **sensor based** base editing
screen into a mapping from sgRNAs to editing outcomes. In particular,
the script uses raw, paired-end read data to classify the editing
events that resulted in the observed target sites. It is assumes that 
the sensor screen has the following structure: 
```
left-scaffold --- sgrna --- middle-scaffold --- target-site --- right-scaffold 
``` 
and that the sgRNA region remains unedited during the screen. 
Reads with edited sgRNA regions are thrown out as errors.

## Script execution

To run the script, first compile the program with `lein` using
```
lein with-profile prod uberjar
```

Then, run the following command to analyze the raw read data.
```
$ java -jar analysis.jar 
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
