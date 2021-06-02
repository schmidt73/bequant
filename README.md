# Processing Sensor Outcomes

To analyze sensor outcomes from raw FASTq reads, we first performed
standard paired-end alignment using the PANDAseq tool. This generated
a pair of merged FASTq files, one for each replicate. This pair of
merged FASTq files was then used directly for sensor outcome analysis.

We first ensured that the merged FASTq reads were not defective as an
artifact of sequencing or cloning. Namely, we first threw out reads
with broken linker constructs and non-matching 5' or 3'
scaffolds. Then, we used the 5' scaffold and linker construct to
associate an sgRNA with the read. sgRNAs that did not match those in
our whitelist were thrown out, as they correspond to
cloning/oligo-construction? errors. With this set of reads, we then
aligned the sensor construct in the read to the target sensor
construct found in the whitelist and classified the outcome. This
yielded a set of outcomes for each sgRNA that we screened.

Sensor constructs of length different than that of the target sensor
were classified as indels; we did not distinguish these further. If
the lengths matched, we classified the outcome as a set of pairs
$(position, mutation_type)$ where $mutation\_type$ could be one of the
six single nucleotide mutations. Typically, most of this information
was thrown out in downstream analysis as we were only concerned with
editing at substrate nucleotides.



