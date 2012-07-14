Required Arguments
------------------

    -g/--ref-genome

A file containing the reference genome in FASTA format. This file will be used throughout the pipeline.

    -r/--reads-basename

The basename of the files containing paired-end reads in FASTQ format. Files should end with _1.fq and _2.fq respectively.

Options
-------

    -o/--output-dir

Sets the name of the directory to which results files will be written. Default is "./ip_out".

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
that determines whether mates are close enough together to form a group.

    -I/--minins

The minimum fragment length for valid paired-end alignments. This value is passed to Bowtie for the initial mapping.

    -X/--maxins

The maximum fragment length for valid paired-end alignments. This value is passed to Bowtie for the initial mapping.

    -p/--num-threads

Number of parallel search threads to use for the mapping and filtering steps. Default is to use a number of threads equal to
the number of available cores on the system. 

    -s/--skip-to

Pipeline step to skip ahead to. This value can be used to restart a partially completed run, or to start the pipeline with
existing data. Value must be one of C (skip to collection), F (skip to filtering), A (skip to assembly), or L (skip to
alignment). 

If C is specified, then -w must be used to specify the existing alignment file to use. If F or A are specified,
-w and -e must be used to specify existing alignment and halfmapping reads files must be used. If L is specified, then
-c must be used to specify existing contigs to use.

    -b/--existing-bwt-index-dir

Path to an existing index to use for mapping. If no index is found, a new one will be created.

    -w/--existing-align-file

A SAM format file containing the result of a Bowtie 2 alignment.

    -e/--existing-halfmap-file

A SAM format file containing half-mapping reads to use for the filtering or assembly steps in the pipeline.

    -c/--existing-contigs-file

A multi-FASTA file containing contigs to use for the alignment step in the pipeline.

    --validate-reads

Runs the validate_fastq.pl validation script on the FASTQ reads files.

    --version

Prints the version number and exits.
