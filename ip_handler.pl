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
  . "    -o/--output-dir               <string>     [ default: ./ip_out  ]\n"
  . "    -f/--fragment-length          <int>        [ default: detect    ]\n"
  . "    -a/--min-mates                <int>        [ default: 3         ]\n"
  . "    -l/--min-contig-length        <int>        [ default: 70        ]\n"
  . "    -m/--max-intron-length        <int>        [ default: 250       ]\n"
  . "    -I/--minins                   <int>        [ default: 0         ]\n"
  . "    -X/--maxins                   <int>        [ default: 3000      ]\n"
  . "    -t/--tolerance                <int>        [ default: 5         ]\n"
  . "    -p/--num-threads              <int>        [ default: all cores ]\n"
  . "    --bowtie1\n"
  . "    --validate-reads\n"
  . "    --version\n";

##################################
# CONFIGURE USER RUNTIME OPTIONS #
##################################

my $ref_genome_file = undef;            # Reference genome (FastA format)
my $mate1s = undef;                     # File containing mate 1s in FastQ format
my $mate2s = undef;                     # File containing mate 2s in FastQ format
my $output_dir = undef;                 # Directory to put output files in
my $output_file = undef;                # Filename to save alignment to; for Galaxy
my $fragment_length = undef;            # Outer distance between mates
my $min_mates = undef;                  # Minimum number of nearby half-mapping mates for local assembly
my $min_contig_length;                  # Minimum contig length for local assembly
my $intron_length = undef;              # Expected intron length for assembly
my $minins = undef;                     # Bowtie -I/--minins <int> value
my $maxins = undef;                     # Bowtie -X/--maxins <int> value
my $tolerance = undef;                  # Alignment position +/- tolerance for simulated pairs
my $num_threads = undef;                # Number of parallel search threads to use for alignment
my $skip_to = undef;                    # Pipeline step to skip ahead to
my $index_dir = undef;                  # Directory containing index files for Bowtie/Blast
my $bowtie1 = undef;                    # Flag to use Bowtie 1 instead of Bowtie 2
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
  "X|maxins:i" => \$maxins,
  "I|minins:i" => \$minins,
  "t|tolerance:s" => \$tolerance,
  "p|num-threads:s" => \$num_threads,
  "s|skip-to:s" => \$skip_to,
  "bowtie1" => \$bowtie1,
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
$minins = $minins || 0;
$maxins = $maxins || 3000;
$tolerance = $tolerance || 5;
$num_threads = $num_threads || Sys::CPU::cpu_count();
$skip_to = $skip_to || "";
$bowtie1 = $bowtie1 || 0;
$index_dir = $index_dir || "./index";
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

# Set the Bowtie version to use
my $bowtie_version;
if ($bowtie1) {
  $bowtie_version = 1;
} else {
  $bowtie_version = 2; # Default is Bowtie 2.
}
$project->set_bowtie_version( $bowtie_version );

# Ensure required executables are available
die "bowtie2 not available on PATH" if system("bowtie2 --version >/dev/null 2>/dev/null") != 0;
die "taipan not available on PATH" if system("which taipan >/dev/null 2>/dev/null") != 0;
die "clustalw not available on PATH" if system("which clustalw >/dev/null 2>/dev/null") != 0;

####################
# RUN THE PIPELINE #
####################

$project->mapping_setup( $index_dir, $validate_reads );
#$project->build_bowtie_index( $index_dir );
$project->build_bwa_index( $index_dir );
#$project->run_bowtie_mapping( $num_threads, $minins, $maxins );
$project->run_bwa_mapping();

#COLLECTION:
#$project->bowtie_identify( $existing_alignment_file );
$project->bwa_identify();

#FILTERING:
$project->set_fragment_length( $fragment_length );
$project->create_simulated_pairs();
#$project->align_simulated_pairs_bowtie( $num_threads, $minins, $maxins );
$project->align_simulated_pairs_bwa();
$project->filter1( $tolerance );
$project->filter2( $tolerance );

#ASSEMBLY:
$project->assemble_groups( $intron_length, $min_mates, $min_contig_length );

#ALIGNMENT:
$project->build_blast_index( $index_dir );
$project->align_groups_blast();
#$project->align_groups_clustal( $output_file );

# Copy latest run to "run-latest" directory for quickly locating the last-run results
#system( "rm -rf run-latest/*" );
#system( "cp -R $work_dir/* run-latest" );

exit;
