#!/usr/bin/perl -w

#  Copyright 2012 Robert Theis
# 
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

use strict;
use IntronPoly;
use Getopt::Long;
#use IPC::System::Simple qw(capture $EXITVAL);

##################################
# CONFIGURE USER RUNTIME OPTIONS #
##################################

# Default general parameters
my $scripts_path = "/home/theis/intron-polymorphism/"; # Location of this and other scripts
my $input_path = "/home/theis/intron-polymorphism/bowtie-0.12.7/reads/"; # Location of reads files
my $input_basename = "e_coli_1000"; # Base of name of reads, without extension ".1.fq" or ".2.fq"
my $input_ref = "/home/theis/intron-polymorphism/bowtie-0.12.7/genomes/NC_008253.fna"; # Bowtie reference genome
my $skipto = ""; # Skip to a step in the pipeline

# Default tools paths
my $bowtie_path = "/home/theis/intron-polymorphism/bowtie-0.12.7"; # Location of bowtie aligner

# Parse command-line options
GetOptions(
  'i=s' => \$input_path,
  'k=s' => \$skipto,
  's=s' => \$scripts_path,
);

# Initialize the project
my $project = IntronPoly->new();
my $workdir = $project->set_workdir( $scripts_path );
$project->set_tooldirs( $bowtie_path );

# Set paths to input data
$project->build_db();

##############################
# JUMP TO THE REQUESTED STEP #
##############################

if ( $skipto eq "C" ) { print "Skipping to collecting\n"; goto COLLECT; }
if ( $skipto eq "F" ) { print "Skipping to filtering\n";  goto FILTER; }
if ( $skipto eq "A" ) { print "Skipping to assembly\n";   goto ASSEMBLE; }
if ( $skipto eq "L" ) { print "Skipping to alignment\n";  goto ALIGN; }
if ( $skipto eq "N" ) { print "Skippint to analysis\n";   goto ANALYZE; }

####################
# RUN THE PIPELINE #
####################

MAP:
$project->build_mapping_db($input_ref);
$project->run_mapping($input_basename);

#COLLECT:


#FILTER:


#ASSEMBLE:


#ALIGN:


#ANALYZE:



