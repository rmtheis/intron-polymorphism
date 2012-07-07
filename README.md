# intron-polymorphism
* * *

This project is under construction--please do not try to use it until it is finished!

## INSTALLATION
**_Linux (64-bit)_**

Install required Perl modules:

    cpan IPC::System::Simple
    cpan Sys::CPU

Install project code:

    git clone git://github.com/rmtheis/intron-polymorphism.git
    cd intron-polymorphism
    chmod +x intron-polymorphism/*.pl

## EXTERNAL DEPENDENCIES

The following packages must be installed and accessible from the command line:

1. [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)

2. [ClustalW](http://www.ebi.ac.uk/Tools/msa/clustalw2/)

3. [Taipan](http://sourceforge.net/projects/taipan/)

## RUNNING

    cd intron-polymorphism
    ./ip_handler.pl -b reads_basename -g ref_genome.fa 1> run.output 2>&1 &
    tail -f run.output

## LICENSE

[Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0.html)
