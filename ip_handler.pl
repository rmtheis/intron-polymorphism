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
my $reads_basename = "m5";             # Read pairs pathname, without "_1.fq" or "_2.fq" extension (FastQ format)
my $minins = undef;                     # Bowtie -I/--minins <int> value
my $maxins = undef;                     # Bowtie -X/--maxins <int> value
my $intron_length = undef;              # Expected intron length for assembly
my $num_aln = undef;                    # Minimum number of nearby half-mapping mates for local assembly
my $min_contig_length;                  # Minimum contig length for local assembly
my $index_dir = undef;                  # Directory containing index files for Bowtie/Blast
my $existing_alignment_file = undef;    # Path to existing SAM alignment data. Leave as undef for new run
my $existing_halfmapping_file = undef;  # Path to existing SAM halfmapping mate data. Leave as undef for new run
my $skip_to = undef;                    # Pipeline step to skip ahead to
my $num_threads = undef;                # Number of parallel search threads to use for alignment

###########################
# INITIALIZE THE PIPELINE #
###########################

# Parse command-line options, overriding any options set above
GetOptions(
  'a:s' => \$existing_alignment_file,
  'b:s' => \$reads_basename,
  'c:s' => \$min_contig_length,
  'g:s' => \$ref_genome_file,
  'h:s' => \$existing_halfmapping_file,
  'idx:s' => \$index_dir,
  'k:s' => \$skip_to,
  'l:i' => \$intron_length,
  'max:i' => \$maxins,
  'min:i' => \$minins,
  'n:i' => \$num_aln,
  't:s' => \$num_threads,
) || die "$0: Bad option";

# Use defaults for undefined values
$ref_genome_file = $ref_genome_file || "testdata/m_chr1.fa";
$reads_basename = $reads_basename || "m_chr1_6";
$minins = $minins || 0;
$maxins = $maxins || 700;
$intron_length = $intron_length || 250;
$num_aln = $num_aln || 3;
$min_contig_length = $min_contig_length || 16;
$index_dir = $index_dir || $ENV{HOME} . "/intron-polymorphism/index";
$existing_alignment_file = $existing_alignment_file || "";
$existing_halfmapping_file = $existing_halfmapping_file || "";
$skip_to = $skip_to || "";
$num_threads = $num_threads || Sys::CPU::cpu_count();

# Resolve relative pathnames
$reads_basename =~ s/^~/$ENV{HOME}/;
$ref_genome_file =~ s/^~/$ENV{HOME}/;
$existing_alignment_file =~ s/^~/$ENV{HOME}/;
$existing_halfmapping_file =~ s/^~/$ENV{HOME}/;

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
$project->mapping_setup( $index_dir );
$project->build_bowtie_index( $index_dir );
$project->run_bowtie_mapping( $num_threads, $minins, $maxins );

COLLECTION:
$project->bowtie_identify( $existing_alignment_file );

FILTERING:
$project->create_simulated_pairs( $existing_alignment_file, $existing_halfmapping_file );
$project->align_simulated_pairs( $num_threads );
$project->filter1();

ASSEMBLY:
$project->assemble_groups( $intron_length, $num_aln, $min_contig_length );

ALIGNMENT:
$project->build_blast_index( $index_dir );
$project->align_groups_blast();
$project->align_groups_clustal();

# Copy latest run to "run-latest" directory for quickly locating the last-run results
system( "rm -rf run-latest/*" );
system( "cp -R $work_dir/* run-latest" );

exit;
