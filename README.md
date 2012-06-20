# intron-polymorphism
* * *

This project is under construction--please do not try to use it until it is finished!

## INSTALLATION
**_Linux (64-bit)_**

Install required Perl modules:

    cpan IPC::System::Simple
    cpan Sys::CPU

Install project code:

    git clone git://github.com/rmtheis/intron-polymorphism.git intron-poly

## EXTERNAL DEPENDENCIES

#### Installing [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)

    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta6/bowtie2-2.0.0-beta6-linux-x86_64.zip
    unzip bowtie2-2.0.0-beta6-linux-x86_64.zip

or Download [here](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta6/bowtie2-2.0.0-beta6-linux-x86_64.zip)

#### Installing [Velvet](http://http://www.ebi.ac.uk/~zerbino/velvet/)

    git clone git://github.com/dzerbino/velvet.git velvet
    git checkout aeb11f8058e4ea794a6ec425c168ffcbbfd1bbbc
    cd velvet
    make

## RUNNING

    cd intron-poly
    ./ip_handler.pl

## LICENSE

[Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0.html)
