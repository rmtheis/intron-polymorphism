# intron-polymorphism
* * *

Discovers intron insertions and deletions from samples of closely related organisms.

## EXTERNAL DEPENDENCIES

The following packages must be installed and added to your PATH:

1. [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)

2. [Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

3. [Taipan](http://sourceforge.net/projects/taipan/)

## INSTALLATION
**_Linux (64-bit)_**

Install required Perl modules:

    cpan IPC::System::Simple
    cpan Sys::CPU

Install project code:

    git clone git://github.com/rmtheis/intron-polymorphism.git

## RUNNING

    cd intron-polymorphism
    perl ip_handler.pl --ref-genome genome.fa -1 myfile_1.fq -2 myfile_2.fq

## PARAMETERS

See [arguments.md](https://github.com/rmtheis/intron-polymorphism/blob/master/arguments.md) 
for a description of the command line parameters.

```text
Usage:
    ip_handler.pl [required arguments] [options]

Required Arguments:
    -g/--ref-genome              <string>
    -1/--mate1-file              <string>
    -2/--mate2-file              <string>

Options:
    -o/--output-dir               <string>     [ default: ./ip_out  ]
    -f/--fragment-length          <int>        [ default: detect    ]
    -a/--min-mates                <int>        [ default: 3         ]
    -l/--min-contig-length        <int>        [ default: 70        ]
    -m/--max-intron-length        <int>        [ default: 250       ]
    -I/--minins                   <int>        [ default: 0         ]
    -X/--maxins                   <int>        [ default: 3000      ]
    -t/--tolerance                <int>        [ default: 10        ]
    -p/--num-threads              <int>        [ default: all cores ]
    -s/--skip-to                  <C,F,A,L>
    -b/--existing-bwt-index-dir   <string>
    -w/--existing-align-file      <string>
    -e/--existing-halfmap-file    <string>
    -c/--existing-contigs-file    <string>
    --bowtie1
    --validate-reads
    --version
```

## LICENSE

[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html)
