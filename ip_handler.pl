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
use lib $FindBin::Bin;
use IntronPoly;
use Getopt::Long;
use Sys::CPU;

use constant VERSION => ("v0.9.0");
my $usage_msg = "Usage:\n"
  . "    ip_handler.pl [required arguments] [options]\n"
  . "\n"
  . "Required Arguments:\n"
  . "    -g/--ref-genome              <string>\n"
  . "    -1/--mate1-file              <string>\n"
  . "    -2/--mate2-file              <string>\n"
  . "\n"
  . "Options:\n"
  . "    -o/--output-dir               <string>     [ default: ip_out    ]\n"
  . "    -f/--fragment-length          <int>        [ default: detect    ]\n"
  . "    -a/--min-mates                <int>        [ default: 3         ]\n"
  . "    -l/--min-contig-length        <int>        [ default: 70        ]\n"
  . "    -m/--max-intron-length        <int>        [ default: 250       ]\n"
  . "    -t/--tolerance-simpair        <int>        [ default: 5         ]\n"
  . "    -r/--tolerance-blast          <int>        [ default: 500       ]\n"
  . "    --trim-reads\n"
  . "    --validate-reads\n"
  . "    --version\n";

##################################
# CONFIGURE USER RUNTIME OPTIONS #
##################################

my $ref_genome_file = undef;            # Reference genome (FastA format)
my $mate1s = undef;                     # File containing mate 1s in FastQ format
my $mate2s = undef;                     # File containing mate 2s in FastQ format
my $output_dir = undef;                 # Directory to put output files in
my $output_file = undef;                # Filename to save contig alignment to; Used by Galaxy
my $fragment_length = undef;            # Outer distance between mates
my $min_mates = undef;                  # Minimum number of nearby half-mapping mates for local assembly
my $min_contig_length;                  # Minimum contig length for local assembly
my $intron_length = undef;              # Expected intron length for assembly
my $tolerance_simpair = undef;          # Alignment position +/- tolerance for simulated pairs
my $tolerance_blast = undef;            # Alignment position +/- tolerance for Blast-aligned pairs
my $index_dir = undef;                  # Directory containing index files for BWA/Blast
my $trim_reads = undef;                 # Flag indicating whether to trim reads
my $validate_reads;                     # Flag indicating whether to validate Fastq reads file
my $version;                            # Flag to print version number and exit
my $help;                               # Flag to print usage message and exit

###########################
# INITIALIZE THE PIPELINE #
###########################

# Parse command-line options, overriding any options set above
GetOptions(
  "g|ref-genome=s" => \$ref_genome_file,
  "1|mate1-file=s" => \$mate1s,
  "2|mate2-file=s" => \$mate2s,
  "o|output-dir:s" => \$output_dir,
  "u|output-file:s" => \$output_file,
  "f|fragment-length:i" => \$fragment_length,
  "a|min-mates:i" => \$min_mates,
  "l|min-contig-length:s" => \$min_contig_length,
  "m|max-intron-length:i" => \$intron_length,
  "t|tolerance-simpair:s" => \$tolerance_simpair,
  "r|tolerance-blast:s" => \$tolerance_blast,
  "trim-reads" => \$trim_reads,
  "validate-reads" => \$validate_reads,
  "v|version" => \$version,
  "h|help" => \$help,
) || die "$0: Bad option";
die $usage_msg unless ( (defined $ref_genome_file && defined $mate1s && defined $mate2s)
    || defined $version || defined $help );

# Use defaults for undefined values
$output_dir = $output_dir || "./ip_out/";
$output_file = $output_file || "";
$fragment_length = $fragment_length || -1;
$min_mates = $min_mates || 3;
$min_contig_length = $min_contig_length || 70;
$intron_length = $intron_length || 250;
$tolerance_simpair = $tolerance_simpair || 5;
$tolerance_blast = $tolerance_blast || 500;
$index_dir = $index_dir || "index";
$trim_reads = $trim_reads || 0;
$validate_reads = $validate_reads || 0;

# Print version number if requested
if ($version) {
  print VERSION . "\n";
  exit;
}

# Print usage message if requested
if ($help) {
  print $usage_msg;
  exit;
}

# Resolve relative pathnames
$mate1s =~ s/^~/$ENV{HOME}/;
$mate2s =~ s/^~/$ENV{HOME}/;
$ref_genome_file =~ s/^~/$ENV{HOME}/;
$output_dir =~ s/^~/$ENV{HOME}/;

# Initialize the project
my $project = IntronPoly->new();

# Create output directory
my $scripts_dir = $FindBin::Bin;
my $work_dir = $project->set_work_dir( $output_dir, $scripts_dir );

# Set paths to data used throughout pipeline
$project->build_db(
  $ref_genome_file,
  $mate1s,
  $mate2s
);

# Ensure required executables are available
die "blastall not available on PATH" if system("which blastall >/dev/null 2>/dev/null") != 0;
die "bwa not available on PATH" if system("which bwa >/dev/null 2>/dev/null") != 0;
die "taipan not available on PATH" if system("which taipan >/dev/null 2>/dev/null") != 0;
die "clustalw not available on PATH" if system("which clustalw >/dev/null 2>/dev/null") != 0;

####################
# RUN THE PIPELINE #
####################

$project->mapping_setup( $index_dir, $validate_reads, $trim_reads);
$project->build_bwa_index( $index_dir );
$project->run_bwa_mapping();

#COLLECTION:
$project->bwa_identify();

#FILTERING:
$project->set_fragment_length( $fragment_length );
$project->create_simulated_pairs();
$project->align_simulated_pairs_bwa();
$project->filter1( $tolerance_simpair );
$project->filter2( $tolerance_blast );

#ASSEMBLY:
$project->assemble_groups( $intron_length, $min_mates, $min_contig_length );

#ALIGNMENT:
$project->build_blast_index( $index_dir );
$project->align_contigs_clustal( $output_file );

# Copy latest run to "run-latest" directory for quickly locating the last-run results
#system( "rm -rf run-latest/*" );
#system( "cp -R $work_dir/* run-latest" );

exit;
