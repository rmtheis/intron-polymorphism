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

use constant VERSION => ("v0.9.0");
my $usage_msg = "Usage:\n"
  . "    ip_handler.pl [required arguments] [options]\n"
  . "\n"
  . "Required Arguments:\n"
  . "    -g/--ref-genome              <string>\n"
  . "    -r/--reads-basename          <string>\n"
  . "\n"
  . "Options:\n"
  . "    -o/--output-dir               <string>     [ default: ./ip_out  ]\n"
  . "    -f/--fragment-length          <int>        [ default: infer     ]\n"
  . "    -a/--min-mates                <int>        [ default: 3         ]\n"
  . "    -l/--min-contig-length        <int>        [ default: 16        ]\n"
  . "    -m/--max-intron-length        <int>        [ default: 250       ]\n"
  . "    -I/--minins                   <int>        [ default: 0         ]\n"
  . "    -X/--maxins                   <int>        [ default: 3000      ]\n"
  . "    -p/--num-threads              <int>        [ default: all cores ]\n"
  . "    -s/--skip-to                  <C,F,A,L>\n"
  . "    -b/--existing-bwt-index-dir   <string>\n"
  . "    -w/--existing-align-file      <string>\n"
  . "    -e/--existing-halfmap-file    <string>\n"
  . "    -c/--existing-contigs-file    <string>\n"
  . "    --validate-reads\n"
  . "    --version\n";

##################################
# CONFIGURE USER RUNTIME OPTIONS #
##################################

my $ref_genome_file = undef;            # Reference genome (FastA format)
my $reads_basename = undef;             # Read pairs basename, without "_1.fq" or "_2.fq" extension (FastQ format)
my $output_dir = undef;                 # Directory to put output files in
my $fragment_length = undef;            # Outer distance between mates
my $min_mates = undef;                  # Minimum number of nearby half-mapping mates for local assembly
my $min_contig_length;                  # Minimum contig length for local assembly
my $intron_length = undef;              # Expected intron length for assembly
my $minins = undef;                     # Bowtie -I/--minins <int> value
my $maxins = undef;                     # Bowtie -X/--maxins <int> value
my $num_threads = undef;                # Number of parallel search threads to use for alignment
my $skip_to = undef;                    # Pipeline step to skip ahead to
my $index_dir = undef;                  # Directory containing index files for Bowtie/Blast
my $existing_alignment_file = undef;    # Path to existing SAM alignment data. Leave as undef for new run
my $existing_halfmapping_file = undef;  # Path to existing SAM halfmapping mate data. Leave as undef for new run
my $existing_contigs_file = undef;      # Path to existing multi-Fasta contigs data. Leave as undef for new run
my $validate_reads;                     # Flag indicating whether to validate Fastq reads file
my $version;                            # Flag to print version number and exit
my $help;                               # Flag to print usage message and exit

###########################
# INITIALIZE THE PIPELINE #
###########################

# Parse command-line options, overriding any options set above
GetOptions(
  "g|ref-genome=s" => \$ref_genome_file,
  "r|reads-basename=s" => \$reads_basename,
  "o|output-dir:s" => \$output_dir,
  "f|fragment-length:i" => \$fragment_length,
  "a|min-mates:i" => \$min_mates,
  "l|min-contig-length:s" => \$min_contig_length,
  "m|max-intron-length:i" => \$intron_length,
  "X|maxins:i" => \$maxins,
  "I|minins:i" => \$minins,
  "p|num-threads:s" => \$num_threads,
  "s|skip-to:s" => \$skip_to,
  "b|existing-bwt-index-dir:s" => \$index_dir,
  "w|existing-align-file:s" => \$existing_alignment_file,
  "e|existing-halfmap-file:s" => \$existing_halfmapping_file,
  "c|existing-contigs-file:s" => \$existing_contigs_file,
  "validate-reads" => \$validate_reads,
  "v|version" => \$version,
  "h|help" => \$help,
) || die "$0: Bad option";
die $usage_msg unless ( (defined $ref_genome_file && defined $reads_basename)
                       || defined $version || defined $help );

# Use defaults for undefined values
$output_dir = $output_dir || "./ip_out/";
$fragment_length = $fragment_length || -1;
$min_mates = $min_mates || 3;
$min_contig_length = $min_contig_length || 16;
$intron_length = $intron_length || 250;
$minins = $minins || 0;
$maxins = $maxins || 3000;
$num_threads = $num_threads || Sys::CPU::cpu_count();
$skip_to = $skip_to || "";
$index_dir = $index_dir || $ENV{HOME} . "/intron-polymorphism/index";
$existing_alignment_file = $existing_alignment_file || "";
$existing_halfmapping_file = $existing_halfmapping_file || "";
$existing_contigs_file = $existing_contigs_file || "";
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

# If jumping to a later pipeline step, ensure required files are available
if ( $skip_to eq "C" && !defined $existing_alignment_file ) {
  die "Must specify --existing-align-file with '--skip_to C'\n";
}
if ( ($skip_to eq "F" && !(defined $existing_alignment_file && defined $existing_halfmapping_file)) ) {
  die "Must specify --existing-align-file and --existing_halfmap-file with '--skip_to F'\n";
}
if ( ($skip_to eq "A" && !(defined $existing_alignment_file && defined $existing_halfmapping_file)) ) {
  die "Must specify --existing-align-file and --existing_halfmap-file with '--skip_to A'\n";
}
if ( $skip_to eq "L" && !defined $existing_contigs_file ) {
  die "Must specify --existing-contigs-file with '--skip_to L'\n";
}

# Resolve relative pathnames
$reads_basename =~ s/^~/$ENV{HOME}/;
$ref_genome_file =~ s/^~/$ENV{HOME}/;
$output_dir =~ s/^~/$ENV{HOME}/;
$index_dir =~ s/^~/$ENV{HOME}/;
$existing_alignment_file =~ s/^~/$ENV{HOME}/;
$existing_halfmapping_file =~ s/^~/$ENV{HOME}/;
$existing_contigs_file =~ s/^~/$ENV{HOME}/;

# Initialize the project
my $project = IntronPoly->new();

# Create output directory
my $scripts_dir = $FindBin::Bin;
my $work_dir = $project->set_work_dir( $output_dir, $scripts_dir );

# Set paths to data used throughout pipeline
$project->build_db(
  $ref_genome_file,
  $reads_basename,
);

# Ensure required executables are available
die "bowtie2 not available on PATH" if system("bowtie2 --version >/dev/null 2>/dev/null") != 0;
die "taipan not available on PATH" if system("which taipan >/dev/null 2>/dev/null") != 0;
die "clustalw not available on PATH" if system("which clustalw >/dev/null 2>/dev/null") != 0;

##############################
# JUMP TO THE REQUESTED STEP #
##############################

if ( $skip_to eq "C" ) { print "Skipping to collecting\n"; goto COLLECTION; }
if ( $skip_to eq "F" ) { print "Skipping to filtering\n";  goto FILTERING; }
if ( $skip_to eq "A" ) { print "Skipping to assembly\n";   goto ASSEMBLY; }
if ( $skip_to eq "L" ) { print "Skipping to alignment\n";  goto ALIGNMENT; }

####################
# RUN THE PIPELINE #
####################

$project->mapping_setup( $index_dir, $validate_reads );
$project->build_bowtie_index( $index_dir );
$project->run_bowtie_mapping( $num_threads, $minins, $maxins );

COLLECTION:
$project->build_bowtie_index( $index_dir );
$project->bowtie_identify( $existing_alignment_file );

FILTERING:
$project->build_bowtie_index( $index_dir );
$project->set_fragment_length( $fragment_length, $existing_alignment_file );
$project->create_simulated_pairs( $existing_alignment_file, $existing_halfmapping_file );
$project->align_simulated_pairs( $num_threads );
$project->filter1( $existing_halfmapping_file );

ASSEMBLY:
$project->set_fragment_length( $fragment_length, $existing_alignment_file );
$project->assemble_groups( $intron_length, $min_mates, $min_contig_length, $existing_halfmapping_file );

ALIGNMENT:
#$project->build_blast_index( $index_dir );
#$project->align_groups_blast( $existing_contigs_file );
$project->align_groups_clustal( $existing_contigs_file );

# Copy latest run to "run-latest" directory for quickly locating the last-run results
system( "rm -rf run-latest/*" );
system( "cp -R $work_dir/* run-latest" );

exit;
