Required Arguments
------------------

    -g/--ref-genome

A file containing the reference genome in FASTA format. This file will be used throughout the stages in the pipeline.

    -1/--mate1-file

The name of the file containing the /1 mates of the paired-end reads in FASTQ format.

    -2/--mate2-file

The name of the file containing the /2 mates of the paired-end reads in FASTQ format.

Options
-------

    -o/--output-dir

Sets the name of the directory to which results files will be written. Default is "ip_out" in the current working directory.

    -f/--fragment-length

Sets the fragment length value for use in the filtering and grouping steps. If no distance is provided, it will be calculated
from the mapping pairs in the initial Bowtie alignment. If the distance is known ahead of time, specifying this option will
result in a shorter overall run time.

    -a/--min-mates

The minimum number of mates required to form a group in the grouping step in the pipeline.

    -l/--min-contig-length

The minimum contig length needed for assembline the mates in a group. This value is passed to the Taipan assembler. Default is 16.

    -m/--max-intron-length

The intron length value to use when looking for groups of introns near one another. This value contributes to a distance value
that determines whether mates are close enough together to form a group. Default is 250.

    -t/--tolerance-simpair

The tolerance with which to consider simulated paired-end reads to have mapped to the reference genome. The tolerance represents the number of base positions behind or forward of the expected alignment position allowed to consider the alignment as acceptable. This value is added and subtracted from the expected alignment position, so for example a value of 5 is considered as plus or minus 5 and creates a window of size 10. Default is 5.

    -r/--tolerance-blast

The tolerance with which to consider Blast alignments of unmapped mates as aligned to the reference genome relative to the expected position of the alignment based on the mapped mate in the pair and the insert length. Default is 500.

    --validate-reads

Runs the validate_fastq.pl validation script on the FASTQ reads files.

    --version

Prints the version number and exits.
