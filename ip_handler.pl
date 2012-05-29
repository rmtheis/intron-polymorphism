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
use Sys::CPU;

##########################################
# CONFIGURE DEFAULT USER RUNTIME OPTIONS #
##########################################

# Skip ahead to a step in the pipeline
my $skip_to = "";

# Directory name where this script and associated scripts are located
my $scripts_dir = "/home/theis/intron-polymorphism";

# Reference genome filename (FastA format)
my $ref_genome_filename = "/home/theis/intron-polymorphism/testdata/NC_008253.fna";

# Previously used output directory name, for restarting a partially completed run 
# Leave this option empty ("") to start a new run of the pipeline
my $resume_work_dir = "";

# Directory where the unmapped reads are located
my $reads_dir = "/home/theis/intron-polymorphism/testdata";

# Base of the unmapped reads filenames, without "_1.fq" or "_2.fq" extension (FastQ format)
# This base filename will be used for all data subsequently generated from these reads.
my $data_basename = "E100000";

# Directory name where the bowtie executables are located
my $bowtie1_dir = "/home/theis/intron-polymorphism/bowtie-0.12.7";
my $bowtie2_dir = "/home/theis/intron-polymorphism/bowtie2-2.0.0-beta5";

# Directory name where existing bowtie index files may be located
my $bowtie_index_dir = "/home/theis/intron-polymorphism/bowtie-index";

# Directory name where the Velvet executable is located
#my $velvet_dir = "/home/theis/intron-polymorphism/velvet";

# Default -I/--minins <int> value for Bowtie
my $minins = 100;

# Default -X/--maxins <int> value for Bowtie
my $maxins = 3000;

###########################
# INITIALIZE THE PIPELINE #
###########################

# Check for invalid input
if ($data_basename eq "") {
  print "Error: data_basename not initialized. Set a value for data_basename.\n";
}

# Set default number of threads to number of available cores
my $bowtie_num_threads = Sys::CPU::cpu_count();

# Parse command-line options, overriding any default options set above
GetOptions(
  'b:s' => \$data_basename,
  'g:s' => \$ref_genome_filename,
  'k:s' => \$skip_to,
  'rd:s' => \$reads_dir,
  't:s' => \$bowtie_num_threads,
);

# Initialize the project
my $project = IntronPoly->new();

# Create output directory
my $work_dir = $project->set_work_dir( $scripts_dir, $resume_work_dir );

# Set paths to data used throughout pipeline
$project->build_db( $ref_genome_filename, $data_basename );

# Track whether mapping-related parameters have been set
my $mapping_setup_completed = 0;

##############################
# JUMP TO THE REQUESTED STEP #
##############################

if ( defined $skip_to && $skip_to ne "" && !(-e $resume_work_dir )) {
  die "$0: resume_work_directory not found";
}

if ( $skip_to eq "C" ) { print "Skipping to collecting\n"; goto COLLECTION; }
if ( $skip_to eq "F" ) { print "Skipping to filtering\n";  goto FILTERING; }
if ( $skip_to eq "A" ) { print "Skipping to assembly\n";   goto ASSEMBLY; }
if ( $skip_to eq "L" ) { print "Skipping to alignment\n";  goto ALIGNMENT; }
if ( $skip_to eq "N" ) { print "Skippint to analysis\n";   goto ANALYSIS; }

####################
# RUN THE PIPELINE #
####################

MAPPING:
print "Running MAPPING...\n";

# Set mapping-related parameters
$project->mapping_setup(
  $bowtie2_dir,
  $bowtie_index_dir,
  $reads_dir,
  $data_basename,
);
$mapping_setup_completed = 1;

# Build the index required for Bowtie to run
$project->build_bowtie2_index();

# Map the reads to the reference genome to identify unaligning pairs
$project->run_bowtie2_mapping( $bowtie_num_threads, $minins, $maxins );

COLLECTION:
print "Running COLLECTION...\n";

# Ensure mapping-related parameters are set, in case we skipped to this step
if ($mapping_setup_completed != 1) {
  $project->mapping_setup(
    $bowtie2_dir,
    $bowtie_index_dir,
    $reads_dir,
    $data_basename,
  );
}

# Identify all the half-mapping read pairs
$project->bowtie2_identify();

FILTERING:
print "Running FILTERING...\n";
$project->filter( $bowtie_num_threads );

ASSEMBLY:
print "Running ASSEMBLY...\n";
$project->group();

#ALIGNMENT:


#ANALYSIS:

# Copy latest run to "run-latest" directory for quickly locating the last-run results
system( "rm -rf run-latest/*" );
system( "cp -R $work_dir/* run-latest" );

exit;
