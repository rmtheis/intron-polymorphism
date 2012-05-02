#intron-polymorphism
* * *

##INTRODUCTION##

This project is under construction--please do not try to use it until it is finished!

##PREREQUISITES##

Bowtie 2
Velvet

##INSTALLATION##
**_Linux (64-bit)_**

#### 1. Install Bowtie 2

Download [here](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta5/bowtie2-2.0.0-beta5-linux-x86_64.zip)

#### 2. Install Velvet

    git clone git://github.com/dzerbino/velvet.git velvet
    git checkout aeb11f8058e4ea794a6ec425c168ffcbbfd1bbbc
    cd velvet
    make

#### 3. Install this project

    git clone git://github.com/rmtheis/intron-polymorphism intron-poly
    cpan File::Basename
    cpan IO::File
    cpan IPC::System::Simple
    cpan Sys::CPU


##RUNNING##

    ip_handler.pl


##MANUAL##


##LICENSE##

Apache 2.0
