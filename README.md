# intron-polymorphism
* * *

Discovers intron insertions and deletions from samples of closely related organisms.

This project is under construction--please do not try to use it until it is finished!

## EXTERNAL DEPENDENCIES

The following packages must be installed and added to your PATH:

1. [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)

2. [ClustalW](http://www.ebi.ac.uk/Tools/msa/clustalw2/)

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

## LICENSE

[Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html)
