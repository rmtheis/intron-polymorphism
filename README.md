# intron-polymorphism
* * *

Discovers intron insertions and deletions from samples of closely related organisms.

## EXTERNAL DEPENDENCIES

The following packages must be installed and added to your PATH:

1. [BWA](http://bio-bwa.sourceforge.net/)

2. [Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

3. [ClustalW](http://www.clustal.org/clustal2/)

4. [Taipan](http://sourceforge.net/projects/taipan/)

## INSTALLATION
**_Linux (64-bit)_**

Install required Perl module:

    sudo cpan IPC::System::Simple

Install project code:

    git clone git://github.com/rmtheis/intron-polymorphism.git

## RUNNING

    cd intron-polymorphism
    perl ip_handler.pl --ref-genome genome.fa -1 myfile_1.fq -2 myfile_2.fq --trim-reads

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
    -o/--output-dir               <string>     [ default: ip_out    ]
    -f/--fragment-length          <int>        [ default: detect    ]
    -a/--min-mates                <int>        [ default: 3         ]
    -l/--min-contig-length        <int>        [ default: 70        ]
    -m/--max-intron-length        <int>        [ default: 250       ]
    -t/--tolerance-simpair        <int>        [ default: 10        ]
    -r/--tolerance-blast          <int>        [ default: 500       ]
    --trim-reads
    --validate-reads
    --version
```

## LICENSE

[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html)
