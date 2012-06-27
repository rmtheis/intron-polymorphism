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
use FindBin;
use IntronPoly;
use Getopt::Long;
use Sys::CPU;

##################################
# CONFIGURE USER RUNTIME OPTIONS #
##################################

my $ref_genome_file = undef;            # Reference genome (FastA format)
my $reads_dir = undef;                  # Directory containing read pairs
my $reads_basename = "m5";             # Read pairs base name, without "_1.fq" or "_2.fq" extension (FastQ format)
my $blast_dir = undef;                  # Directory to Blast executable, if not in $PATH
my $bowtie_dir = undef;                 # Directory to Bowtie executable, if not in $PATH
my $velvet_dir = undef;                 # Directory to Velvet executable, if not in $PATH
my $minins = undef;                     # Bowtie -I/--minins <int> value
my $maxins = undef;                     # Bowtie -X/--maxins <int> value
my $intron_length = undef;              # Expected intron length for assembly
my $num_aln = undef;                    # Minimum number of nearby half-mapping mates for local assembly
my $cov_cutoff = undef;                 # Velvet coverage cutoff
my $hash_length = undef;                # Velvet hash length
my $index_dir = undef;                  # Directory containing index files for Bowtie/Blast
my $existing_alignment_file = undef;    # Path to existing SAM data. Leave as undef for new run
my $existing_halfmapping_file = undef;  # Path to existing SAM data. Leave as undef for new run
my $skip_to = undef;                    # Pipeline step to skip ahead to
my $bowtie_num_threads = undef;         # Number of Bowtie parallel search threads

###########################
# INITIALIZE THE PIPELINE #
###########################

# Parse command-line options, overriding any options set above
GetOptions(
  'a:s' => \$existing_alignment_file,
  'b:s' => \$reads_basename,
  'bw:s' => \$bowtie_dir,
  'c:i' => \$cov_cutoff,
  'g:s' => \$ref_genome_file,
  'h:s' => \$existing_halfmapping_file,
  'hl:s' => \$hash_length,
  'idx:s' => \$index_dir,
  'k:s' => \$skip_to,
  'l:i' => \$intron_length,
  'max:i' => \$maxins,
  'min:i' => \$minins,
  'n:i' => \$num_aln,
  'rd:s' => \$reads_dir,
  's:s' => \$blast_dir,
  't:s' => \$bowtie_num_threads,
  'v:s' => \$velvet_dir,
) || die "$0: Bad option";

# Use defaults for undefined values
$ref_genome_file = $ref_genome_file || $ENV{HOME} . "/intron-polymorphism/testdata/m_chr1.fa";
$reads_dir = $reads_dir || $ENV{HOME} . "/intron-polymorphism/testdata";
$reads_basename = $reads_basename || "m_chr1_6";
$blast_dir = $blast_dir || "";
$bowtie_dir = $bowtie_dir || "";
$velvet_dir = $velvet_dir || "";
$minins = $minins || 0;
$maxins = $maxins || 700;
$intron_length = $intron_length || 250;
$num_aln = $num_aln || 3;
$cov_cutoff = $cov_cutoff || 2;
$hash_length = $hash_length || 15;
$index_dir = $index_dir || $ENV{HOME} . "/intron-polymorphism/index";
$existing_alignment_file = $existing_alignment_file || "";
$existing_halfmapping_file = $existing_halfmapping_file || "";
$skip_to = $skip_to || "";
$bowtie_num_threads = $bowtie_num_threads || Sys::CPU::cpu_count();

# Initialize the project
my $project = IntronPoly->new();

# Create output directory
my $scripts_dir = $FindBin::Bin;
my $work_dir = $project->set_work_dir( $scripts_dir );

# Set paths to data used throughout pipeline
$project->build_db(
  $ref_genome_file,
  $reads_basename,
);

##############################
# JUMP TO THE REQUESTED STEP #
##############################

if ( $skip_to eq "C" ) { print "Skipping to collecting\n"; goto COLLECTION; }
if ( $skip_to eq "F" ) { print "Skipping to filtering\n";  goto FILTERING; }
if ( $skip_to eq "A" ) { print "Skipping to assembly\n";   goto ASSEMBLY; }
if ( $skip_to eq "L" ) { print "Skipping to alignment\n";  goto ALIGNMENT; }
if ( $skip_to eq "N" ) { print "Skipping to analysis\n";   goto ANALYSIS; }

####################
# RUN THE PIPELINE #
####################

MAPPING:
print "Running MAPPING...\n";
$project->mapping_setup(
  $bowtie_dir,
  $index_dir,
  $reads_dir,
  $reads_basename,
  $velvet_dir,
  $blast_dir,
);
$project->build_bowtie_index( $index_dir );
$project->run_bowtie_mapping( $bowtie_num_threads, $minins, $maxins );

COLLECTION:
print "Running COLLECTION...\n";
$project->bowtie_identify( $existing_alignment_file );

FILTERING:
print "Running FILTERING...\n";
$project->create_fake_pairs( $existing_alignment_file, $existing_halfmapping_file );
$project->run_fake_pairs_alignment( $bowtie_num_threads );
$project->filter1();

ASSEMBLY:
print "Running ASSEMBLY...\n";
$project->assemble_groups_velvet( $intron_length, $num_aln, $hash_length, $cov_cutoff );

ALIGNMENT:
print "Running ALIGNMENT...\n";
$project->build_blast_index( $index_dir );
$project->align_groups();

#ANALYSIS:

# Copy latest run to "run-latest" directory for quickly locating the last-run results
system( "rm -rf run-latest/*" );
system( "cp -R $work_dir/* run-latest" );

exit;
