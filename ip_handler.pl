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

##########################################
# CONFIGURE DEFAULT USER RUNTIME OPTIONS #
##########################################

# Skip ahead to a step in the pipeline
my $skip_to = "";

# Previously used output directory name, for restarting a partially completed run 
my $resume_work_dir = "";

# Directory name where this script and associated scripts are located
my $scripts_dir = "/home/theis/intron-polymorphism";

# Reference genome filename (FastA format)
my $ref_genome_filename = "/home/theis/intron-polymorphism/genomes/NC_008253.fna";

# Directory where the unmapped reads are located
my $reads_dir = "/home/theis/intron-polymorphism/reads";

# Base of the unmapped reads filenames, without ".1.fq" or ".2.fq" extension (FastQ format)
my $reads_basename = "testData";

# Directory name where the bowtie executables are located
my $bowtie_dir = "/home/theis/intron-polymorphism/bowtie-0.12.7";

# Directory name where existing bowtie index files may be located
my $bowtie_index_dir = "/home/theis/intron-polymorphism/bowtie-index";

# Number of threads to use when running bowtie
my $bowtie_num_threads = 8;

###########################
# INITIALIZE THE PIPELINE #
###########################

# Parse command-line options, overriding any default options set above
GetOptions(
  'g:s' => \$ref_genome_filename,
  'k:s' => \$skip_to,
  'rd:s' => \$reads_dir,
  'rb:s' => \$reads_basename,
  't:s' => \$bowtie_num_threads,
);

# Initialize the project
my $project = IntronPoly->new();

# Create output directory
my $work_dir = $project->set_work_dir( $scripts_dir, $resume_work_dir );

# Set paths to data used throughout pipeline
$project->build_db( $ref_genome_filename );

##############################
# JUMP TO THE REQUESTED STEP #
##############################

if ( $skip_to eq "C" ) { print "Skipping to collecting\n"; goto COLLECT; }
if ( $skip_to eq "F" ) { print "Skipping to filtering\n";  goto FILTER; }
if ( $skip_to eq "A" ) { print "Skipping to assembly\n";   goto ASSEMBLE; }
if ( $skip_to eq "L" ) { print "Skipping to alignment\n";  goto ALIGN; }
if ( $skip_to eq "N" ) { print "Skippint to analysis\n";   goto ANALYZE; }

####################
# RUN THE PIPELINE #
####################

MAP:
$project->mapping_setup(
  $bowtie_dir,
  $bowtie_index_dir,
  $reads_dir,
  $reads_basename,
);
$project->build_mapping_index();
$project->run_mapping( $bowtie_num_threads );
$project->collect_half_mapping_pairs( $bowtie_num_threads );

COLLECT: # find half mapping mates that are close to one another

FILTER:


ASSEMBLE:


ALIGN:


ANALYZE:


exit;
