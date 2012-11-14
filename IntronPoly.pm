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

package IntronPoly;

use strict;
use File::Basename;
use IO::File;
use IPC::System::Simple qw(capture $EXITVAL);

use constant DEBUG => (1);    # Flag to print info for debugging

=head2 new

 Title   : new
 Usage   : $project = IntronPoly->new()
 Function: Initializes a new intron polymorphism analysis project
 Example : $project = IntronPoly->new()
           $project->set_work_dir( "/home/theis/work" )
 Returns : IntronPoly project object
 Args    : None

=cut

sub new {
  my $self = {};
  $self->{"scripts_dir"}    = undef;   # Scalar path to scripts directory
  $self->{"work_dir"}       = undef;   # Scalar path to work directory
  $self->{"reads"}          = undef;   # Hash of reads related values
  $self->{"ref_genome"}     = undef;   # Hash of reference genome related values
  $self->{"bowtie_db"}      = undef;   # Hash of data directories for bowtie
  $self->{"alignment_file"} = undef;   # Scalar path to initial read pairs alignments
  $self->{"read_length"}    = undef;   # Median read length determined from alignment
  $self->{"frag_length"}    = undef;   # Fragment length
  bless($self);
  return $self;
}

=head2 set_work_dir

 Title   : set_work_dir
 Usage:  : $project->set_work_dir( "path_to_scripts" )
 Function: Creates directory for pipeline working data and output files
 Example : my $work_dir = "/home/theis/work";
           $project->set_work_dir( $work_dir )
 Returns : The path to the directory for pipeline working data and output files
 Args    : Scalar of full path to the parent directory of the work directory

=cut

sub set_work_dir {
  my $self        = shift;
  my $output_dir  = shift;
  my $scripts_dir = shift || "";
  $scripts_dir =~ s!/*$!/! if ($scripts_dir ne ""); # Add trailing slash if not already present
  $self->{"scripts_dir"} = $scripts_dir;
  $output_dir =~ s!/*$!/! if ($output_dir ne ""); # Add trailing slash if not already present
  my $work_dir = shift || "${output_dir}run-" . &_datestamp;
  unless ( -e $output_dir ) {
    &_make_dir($output_dir);
  }
  unless ( -e $work_dir ) {
    $work_dir = $1 if ( $work_dir =~ /(.*)\/$/ );
    &_make_dir($work_dir);
    print "Created output directory $work_dir\n";
  }
  else {
    $work_dir = &_check_dir($work_dir);
    print "Using existing output directory $work_dir\n";
  }
  $self->{"work_dir"} = $work_dir;

  my $assembly_dir = "$work_dir/assembly";
  unless ( -d $assembly_dir ) {
    &_make_dir($assembly_dir);
    print "Created assembly data directory $assembly_dir\n";
  }
  $self->{"assembly_dir"} = $assembly_dir;
  return $work_dir;
}

=head2 build_db

 Title   : build_db
 Usage   : $project->build_db( "path_to_reference_genome", "read_pairs_base_name" )
 Function: Sets references to data files used in multiple steps in the pipeline
 Example : my $ref_genome = "/home/theis/genome/mygenome.fna";
           $project->build_db( $ref_genome, "reads_1.fq", "reads_2.fq" );
 Returns : No return value
 Args    : Scalar of full path to the reference genome Fasta file, scalar of the FastQ
           file containing mate 1s, and scalar of the FastQ file containing mate 2s

=cut

sub build_db {
  my $self                 = shift;
  my $ref_genome_filename  = shift;
  my $mate1s               = shift;
  my $mate2s               = shift;

  die "Cannot find $ref_genome_filename" if !(-e $ref_genome_filename);
  die "Cannot find $mate1s" if !(-e $mate1s);
  die "Cannot find $mate2s" if !(-e $mate2s);
  
  # Use the base of the read pairs filenames for saving data throughout pipeline
  {
    my ( $file, $dir, $ext ) = fileparse( $mate1s, qr/\.[^.]*/ );
    $file =~ s/_1$//;
    $self->{"reads"}->{"mate1s"}   = $mate1s;
    $self->{"reads"}->{"mate2s"}   = $mate2s;
    $self->{"reads"}->{"basename"} = $file;
    $self->{"reads"}->{"dir"}      = $dir;
  }
  
  # Split reference genome pathname into components
  {
    my ( $file, $dir, $ext ) = fileparse( $ref_genome_filename, qr/\.[^.]*/ );
    $self->{"ref_genome"}->{"full_pathname"} = $ref_genome_filename;
    $self->{"ref_genome"}->{"basename"}      = $file;
    $self->{"ref_genome"}->{"dir"}           = $dir;
    $self->{"ref_genome"}->{"extension"}     = $ext;
  }
}

=head2 set_bowtie_version

 Title   : set_bowtie_version
 Usage   : $project->set_bowtie_version( $ver )
 Function: Sets major version number of Bowtie to use: Bowtie 2 or Bowtie 1.
 Example : $project->set_bowtie_version( 2 );
 Returns : No return value
 Args    : Use 1 for Bowtie 1. Use 2 for Bowtie 2.

=cut

sub set_bowtie_version {
  my $self     = shift;
  my $version  = shift;
  $self->{"bowtie_version"} = $version;
}

=head2 mapping_setup

 Title   : mapping_setup
 Usage   : $project->mapping_setup( "index_dir", "reads_dir", "read_pairs_base_name")
 Function: Sets references to data files used for mapping reads to the reference genome, and validates
           input files
 Example : $project->mapping_setup( "/tmp" );
 Returns : No return value
 Args    : Scalar full path to the index directory for the alignment program (BWA or Bowtie 1/2),
           flag indicating whether to validate the reads file (optional), flag indicating whether
           to trim the reads using the trim2.pl adaptive trimming script (optional)

=cut

sub mapping_setup {
  my $self                 = shift;
  my $index_dir            = shift;
  my $validate_reads       = shift || 0;
  my $trim_reads           = shift || 0;
  my $scripts_dir          = $self->{"scripts_dir"};
  my $work_dir             = $self->{"work_dir"};
  my $ref_genome           = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename  = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_dir       = $self->{"ref_genome"}->{"dir"};
  my $reads_dir            = $self->{"reads"}->{"dir"};
  my $reads_basename       = $self->{"reads"}->{"basename"};
  my $mate1s               = $self->{"reads"}->{"mate1s"};
  my $mate2s               = $self->{"reads"}->{"mate2s"};
  my $trim_output_file     = "${work_dir}/${reads_basename}_trim.out";
  die "$0: reference genome directory $ref_genome_dir does not exist" unless ( -d $ref_genome_dir );
  
  print "Validating Fasta reference genome data...\n";
  
  # Perform basic validation on reference genome
  capture("${scripts_dir}validate_fasta.pl -i $ref_genome");
  die "$0: validate_fasta.pl exited unsuccessful" if ( $EXITVAL != 0 );

  # Ensure that read pairs files exist
  $reads_dir =~ s!/*$!/! if ($reads_dir ne ""); # Add trailing slash if not already present
  die "$0: reads file $mate1s does not exist" unless ( -e $mate1s );
  die "$0: reads file $mate2s does not exist" unless ( -e $mate2s );
  
  if ($validate_reads) {
    # Perform basic validation on read pairs files
    print "Validating Fastq read pair data...\n";   
    capture("${scripts_dir}validate_fastq.pl -1 ${reads_dir}$mate1s -2 ${reads_dir}$mate2s");
    die "$0: validate_fastq.pl exited unsuccessful" if ( $EXITVAL != 0 );
  }

  if ($trim_reads) {
    # Trim reads by quality scores using the trim2.pl adaptive trim script, t=3
    print "Trimming reads...\n";
    capture("${scripts_dir}trim2.pl $mate1s $mate2s 3 40 > $trim_output_file");
    die "$0: trim2.pl exited unsuccessful" if ( $EXITVAL != 0 );
    $mate1s .= ".trimmed";
    $mate2s .= ".trimmed";
  }
  
  # Create directory for aligner index files if necessary
  $index_dir = $1 if ( $index_dir =~ /(.*)\/$/ );
  unless ( -d $index_dir ) {
    &_make_dir($index_dir);
  }

  # Remember paths needed for running
  $self->{"bowtie_db"} = {
    reads_file_one => $mate1s,
    reads_file_two => $mate2s
  };
  $self->{"index_dir"} = $index_dir;
  print "Using reference genome $ref_genome\n";
  print "Using Fastq file: $mate1s\n";
  print "Using Fastq file: $mate2s\n";
}

=head2 build_bowtie_index

 Title   : build_bowtie_index
 Usage   : $project->build_bowtie_index();
 Function: Runs the Bowtie indexer to create a Bowtie index from the reference genome
 Example : $project->mapping_setup( $index_dir, $reads_dir, "reads" );
           $project->build_bowtie_index();
 Returns : No return value
 Args    : Directory containing Bowtie index files (optional)

=cut

sub build_bowtie_index {
  my $self                = shift;
  my $index_dir           = shift || $self->{"index_dir"};
  my $ref_genome          = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie_version      = $self->{"bowtie_version"};
  $self->{"index_dir"}    = $index_dir;
  
  if ($bowtie_version == 2) {
    # Bowtie 2
    unless ( -e "$index_dir/$ref_genome_basename.1.bt2"
      && "$index_dir/$ref_genome_basename.2.bt2"
      && "$index_dir/$ref_genome_basename.3.bt2"
      && "$index_dir/$ref_genome_basename.4.bt2"
      && "$index_dir/$ref_genome_basename.rev.1.bt2"
      && "$index_dir/$ref_genome_basename.rev.2.bt2" )
    {
      # Create the Bowtie index
      print "Creating bowtie index in $index_dir...\n";
      capture( "bowtie2-build $ref_genome $index_dir/$ref_genome_basename" );
      die "$0: bowtie2-build exited unsuccessful" if ( $EXITVAL != 0 );
    }
  } else {
    # Bowtie 1
    unless ( -e "$index_dir/$ref_genome_basename.1.ebwt"
      && "$index_dir/$ref_genome_basename.2.ebwt"
      && "$index_dir/$ref_genome_basename.3.ebwt"
      && "$index_dir/$ref_genome_basename.4.ebwt"
      && "$index_dir/$ref_genome_basename.rev.1.ebwt"
      && "$index_dir/$ref_genome_basename.rev.2.ebwt" )
    {
      # Create the Bowtie index
      print "Creating bowtie index in $index_dir...\n";
      capture( "bowtie-build $ref_genome $index_dir/$ref_genome_basename" );
      die "$0: bowtie-build exited unsuccessful" if ( $EXITVAL != 0 );
    }
  }
}

=head2 build_bwa_index

 Title   : build_bwa_index
 Usage   : $project->build_bwa_index();
 Function: Runs BWA indexing to create a BWA index from the reference genome
 Example : $project->mapping_setup( $index_dir );
           $project->build_bwa_index();
 Returns : No return value
 Args    : Directory containing BWA index files (optional)

=cut

sub build_bwa_index {
  my $self                 = shift;
  my $index_dir            = shift || $self->{"index_dir"};
  my $ref_genome           = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename  = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_extension = $self->{"ref_genome"}->{"extension"};
  $self->{"index_dir"}     = $index_dir;
  my $mate1s               = $self->{"reads"}->{"mate1s"};
  my $mate2s               = $self->{"reads"}->{"mate2s"};
  my $mate_basename        = $self->{"reads"}->{"basename"};
  my $sai1                 = $mate_basename . ".sai1";
  my $sai2                 = $mate_basename . ".sai2";
  $self->{"sai1"}          = $sai1;
  $self->{"sai2"}          = $sai2;

  # Copy reference genome to index directory so BWA can create index there
  my $copied_genome = "$index_dir/$ref_genome_basename$ref_genome_extension";
  unless ( -e $copied_genome) {
    capture( "cp $ref_genome $copied_genome" );
    die "$0: cp exited unsuccessful" if ( $EXITVAL != 0 );
  }
  
  # Create index for reference genome
  unless ( -e "$copied_genome.amb"
    && "$copied_genome.ann"
    && "$copied_genome.bwt"
    && "$copied_genome.pac"
    && "$copied_genome.sa" )
  {
    # Create the BWA index.
    print "Creating bwa index in $index_dir...\n";
    capture( "bwa index $copied_genome > $index_dir/$ref_genome_basename.out 2> $index_dir/$ref_genome_basename.err" );
    die "$0: bwa index exited unsuccessful" if ( $EXITVAL != 0 );
  }

  # Create SAI files
  unless ( -e "$index_dir/$sai1" && "$index_dir/$sai2") {
    capture( "bwa aln $copied_genome $mate1s > $index_dir/$sai1 2> $mate1s.err" );
    die "$0: bwa aln exited unsuccessful" if ( $EXITVAL != 0 );
    capture( "bwa aln $copied_genome $mate2s > $index_dir/$sai2 2> $mate2s.err" );
    die "$0: bwa aln exited unsuccessful" if ( $EXITVAL != 0 );
  }
}

=head2 run_bowtie_mapping

 Title   : run_bowtie_mapping
 Usage   : $project->run_bowtie_mapping( num_threads, minins, maxins )
 Function: Aligns reads to the reference genome and saves output file
 Example : $project->mapping_setup( $index_dir, $reads_dir, "reads" );
           $project->build_bowtie_index();
           $project->run_bowtie_mapping( 8, 100, 500 );
 Returns : No return value
 Args    : Number of parallel search threads to use for Bowtie, numeric length to use for
           Bowtie -I/--minins parameter, numeric length to use for Bowtie -X/--maxins parameter

=cut

sub run_bowtie_mapping {
  my $self                  = shift;
  my $num_threads           = shift;
  my $minins                = shift;
  my $maxins                = shift;
  my $bowtie_version        = $self->{"bowtie_version"};
  my $mate1s                = $self->{"reads"}->{"mate1s"};
  my $mate2s                = $self->{"reads"}->{"mate2s"};
  my $ref_genome_basename   = $self->{"ref_genome"}->{"basename"};
  my $index_dir             = $self->{"index_dir"};
  my $reads_file_one        = $self->{"bowtie_db"}->{"reads_file_one"};
  my $reads_file_two        = $self->{"bowtie_db"}->{"reads_file_two"};
  my $work_dir              = $self->{"work_dir"};
  my $reads_basename        = $self->{"reads"}->{"basename"};
  my $prefix                = "$work_dir/${reads_basename}";
  my $output_file           = "${prefix}_initial_alignments.sam";
  $self->{"alignment_file"} = $output_file;
  
  if ($bowtie_version == 2) {
    
    # Call Bowtie 2 to run the mapping
    print "Mapping using Bowtie 2, using $num_threads threads...\n";
    capture(
          "bowtie2 -x $index_dir/$ref_genome_basename "
        . "--threads $num_threads --reorder --no-hd "
        . "--maxins $maxins --minins $minins "
        . "--no-discordant "
        . "--no-contain --no-overlap "
        . "-k 3 -1 $reads_file_one -2 $reads_file_two "
        . "-S $output_file 2> $work_dir/${reads_basename}_bowtie2_initial_mapping_stats.txt"
    );
    die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
    print "Saved alignment data to $output_file\n";
    
  } else {
    
    # Call Bowtie 1 to run the mapping
    print "Mapping using Bowtie 1, using $num_threads threads...\n";
    capture(
          "bowtie $index_dir/$ref_genome_basename " .
          "--threads $num_threads " .
          "--maxins $maxins --minins $minins " .
          "--un ${prefix}_initial_unaligned.fq " . 
          "-m 1 -1 $reads_file_one -2 $reads_file_two " .
          "-S $output_file 2> $work_dir/${reads_basename}_bowtie1_initial_mapping_stats.txt"
    );
    die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
  
    # Identify the half-mapping read pairs by running unaligned mates individually.
    # Read pairs are already those that align uniquely to the reference sequence because Bowtie 1 was
    # run with option '-m 1'.
    capture(
          "bowtie $index_dir/$ref_genome_basename " .
          "--threads $num_threads " .
          "-m 1 ${prefix}_initial_unaligned_1.fq " . 
          "--sam ${prefix}_mate1.sam " .
          "2> $work_dir/${reads_basename}_bowtie1_secondary_mapping_stats1.txt"
    );
    die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
    unlink "${prefix}_initial_unaligned_1.fq" or die "Can't delete ${prefix}_initial_unaligned_1.fq: $!";
    
    capture(
          "bowtie $index_dir/$ref_genome_basename " .
          "--threads $num_threads " .
          "-m 1 ${prefix}_initial_unaligned_2.fq " . 
          "--sam ${prefix}_mate2.sam " .
          "2> $work_dir/${reads_basename}_bowtie1_secondary_mapping_stats2.txt"
    );
    die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
    unlink "${prefix}_initial_unaligned_2.fq" or die "Can't delete ${prefix}_initial_unaligned_2.fq: $!";
    
    # Concatenate the SAM files for the individual mates
    capture( "cat ${prefix}_mate1.sam ${prefix}_mate2.sam > ${prefix}_combined.sam");
    die "$0: cat exited unsuccessful" if ( $EXITVAL != 0 );
    unlink "${prefix}_mate1.sam" or die "Can't delete ${prefix}_mate1.sam: $!";
    unlink "${prefix}_mate2.sam" or die "Can't delete ${prefix}_mate2.sam: $!";
    
    # Sort by ID
    capture("sort -V -k1,1 -o ${prefix}_combined_sorted.sam ${prefix}_combined.sam");
    die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );
    unlink "${prefix}_combined.sam" or die "Can't delete ${prefix}_combined.sam: $!";
    
  }
}

=head2 run_bwa_mapping

 Title   : run_bwa_mapping
 Usage   : $project->run_bwa_mapping( num_threads, minins, maxins )
 Function: Aligns reads to the reference genome and saves output file
 Example : $project->mapping_setup( $index_dir );
           $project->build_bwa_index( $index_dir );
           $project->run_bwa_mapping();
 Returns : No return value
 Args    : None

=cut

sub run_bwa_mapping {
  my $self                  = shift;
  my $mate1s                = $self->{"reads"}->{"mate1s"};
  my $mate2s                = $self->{"reads"}->{"mate2s"};
  my $ref_genome_basename   = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_extension  = $self->{"ref_genome"}->{"extension"};
  my $index_dir             = $self->{"index_dir"};
  my $reads_file_one        = $self->{"bowtie_db"}->{"reads_file_one"};
  my $reads_file_two        = $self->{"bowtie_db"}->{"reads_file_two"};
  my $work_dir              = $self->{"work_dir"};
  my $scripts_dir             = $self->{"scripts_dir"};
  my $reads_basename        = $self->{"reads"}->{"basename"};
  my $prefix                = "$work_dir/${reads_basename}";
  my $initial_align_file    = "${prefix}_initial_alignments.sam";
  $self->{"alignment_file"} = $initial_align_file;
  my $sai1                  = $self->{"sai1"};
  my $sai2                  = $self->{"sai2"};
  my $second_align_file     = "${prefix}_unaligned.sam";

  # Call BWA to run the mapping
  print "Mapping using BWA...\n";
  my $ref = "$index_dir/$ref_genome_basename$ref_genome_extension";
  capture("bwa sampe -s $ref $index_dir/$sai1 $index_dir/$sai2 $mate1s $mate2s > $initial_align_file " .
          "2> $work_dir/${reads_basename}_bwa_initial_mapping_stats.txt"
  );
  die "$0: BWA exited unsuccessful" if ( $EXITVAL != 0 );

  # Write the unaligned mates to a file
  my $ifh = IO::File->new( $initial_align_file, 'r' ) or die "Can't open $initial_align_file: $!";
  my $ofh = IO::File->new( $second_align_file, 'w' ) or die "Can't open $second_align_file: $!";
  my $prev_id         = "";
  my $discard_this_id = 0;
  my $mapping         = 0;
  my @print_lines = ();
  while ( my $line = $ifh->getline ) {
    next if $line =~ m/^@/;
    my @fields = split( /\t/, $line );
    my ($mate_id, $flags) = ($fields[0], $fields[1]);
    if (&_isInNonMappingPair($flags)) {
      print $ofh $line;
    }
  }
  close $ifh;
  close $ofh;
}

=head2 bowtie_identify

 Title   : bowtie_identify
 Usage   : $project->bowtie_identify()
 Function: Identifies half-mapping read pairs from the read alignment output file
 Example : $project->mapping_setup( $index_dir );
           $project->build_bowtie_index();
           $project->bowtie_identify();
 Returns : No return value
 Args    : Scalar path to alignment file to use (optional)

=cut

sub bowtie_identify {
  my $self                    = shift;
  my $alignment_file          = shift || $self->{"alignment_file"};
  my $bowtie_version          = $self->{"bowtie_version"};
  my $reads_basename          = $self->{"reads"}->{"basename"};
  my $scripts_dir             = $self->{"scripts_dir"};
  my $work_dir                = $self->{"work_dir"};
  my $bowtie1_infile          = "$work_dir/${reads_basename}_combined_sorted.sam";
  my $outfile                 = "$work_dir/${reads_basename}_halfm_pairs.sam";
  my $outfile2                = "$work_dir/${reads_basename}_halfm_unal_mates.sam";
  my $outfile3                = $outfile;
  $outfile3                   =~ s/(\.[^.]+)$/_sorted.sam/;
  my $outfile4                = "$work_dir/${reads_basename}_halfm_unal_mates_unsorted.sam";
  $self->{"halfmapping_file"} = $outfile2;
  
  if (!-e $alignment_file) {
    print STDERR "Alignment file not found.\n";
    return;
  }
  $scripts_dir =~ s!/*$!/! if ($scripts_dir ne ""); # Add trailing slash if not already present
  print "Identifying half-mapping read pairs...\n";
  my $count = 0;
  if ($bowtie_version == 2) {
    
    # Bowtie 2
    
    my $ifh  = IO::File->new( $alignment_file, 'r' ) or die "Can't open $alignment_file: $!";
    my $ofh  = IO::File->new( $outfile, 'w' ) or die "Can't create $outfile: $!";
    my $ofh2 = IO::File->new( $outfile2, 'w' ) or die "Can't create $outfile2: $!";
    my $prev_id         = "";
    my $discard_this_id = 0;
    my $mapping         = 0;
    my @print_lines = ();
    while ( my $line = $ifh->getline ) {
      next if $line =~ m/^@/;
      my @fields = split( /\t/, $line );
      my ($mate_id, $flags) = ($fields[0], $fields[1]);
  
      # Check if this mate ID is different from the last line
      if ( $mate_id ne $prev_id ) {
        # This line is the first line with a new read ID
        if ( scalar @print_lines == 2 && $mapping == 1 && $prev_id ne "" ) {
          my $line1  = shift(@print_lines);
          my @fd1    = split( /\t/, $line1 );
          my $flags1 = $fd1[1];
          my $line2  = shift(@print_lines);
          my @fd2    = split( /\t/, $line2 );
          my $flags2 = $fd2[1];
  
          # Print the half-mapping mate pairs with the last line's ID
          print $ofh "$line1$line2";
          if (&_isUnalignedMate($flags1)) {
            print $ofh2 $line1;
          } else {
            print $ofh2 $line2;
          }
          $count++;
        }
  
        # Discard mates with multiple alignments here
        $discard_this_id = 0;
        $prev_id         = $mate_id;
      } else {
        # This line is the second or greater line with this read ID
        next if $discard_this_id;
        if ( scalar @print_lines == 2 && $mapping == 0 ) {
          # Found a secondary alignment preceded by two half-mapping mates. Discard all lines with this ID.
          @print_lines     = ();
          $discard_this_id = 1;
          next;
        } elsif ( scalar @print_lines == 3 ) {
          # We have four alignments for this read pair. Discard all lines with this ID.
          @print_lines     = ();
          $discard_this_id = 1;
          next;
        }
      }
  
      # At this point we have an alignment we may want to keep, so look at the sum-of-flags value to classify.
      if ( &_isInMappingPair($flags) ) {
        # Save this full-mapping line for printing
        $mapping = 2;
        push( @print_lines, $line );
      } elsif ( &_isInNonMappingPair($flags) ) {
        # We have a non-mapping ID. Discard all lines with this ID.
        $mapping         = 0;
        $discard_this_id = 1;
      } elsif ( &_isHalfMappingMate($flags) ) {
        # Save this half-mapping line for printing
        $mapping = 1;
        push( @print_lines, $line );
      } else {
        # We have a secondary alignment.
        $mapping = 2;
        push( @print_lines, $line );
      }
    }
    $ifh->close;
    $ofh->close;
    
    capture( "sort -k3,3 -k8n,8 -o $outfile3 $outfile" );
    die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );
    
  } else {
    
    # Bowtie 1
    
    # Get the alignment flags, and remove the cases where both mates fail to align
    # We have unpaired alignments, so possible flag values are:
    # 0 or 16 for alignment to the forward or reverse strand respectively
    # 4 for no reported alignments
    my $ifh = IO::File->new( $bowtie1_infile, 'r' ) or die "Can't open $bowtie1_infile: $!";
    my $ofh = IO::File->new( $outfile4, 'w' ) or die "Can't open $outfile4: $!";
    my $ofh2 = IO::File->new( $outfile, 'w' ) or die "Can't open $outfile: $!"; 
    while ( my $line1 = $ifh->getline ) {
      next if $line1 =~ m/^@/;
  
      # Read in both mates in the pair
      my @a = split( /\t/, $line1 );
      my ($mate1_id, $flags1) = ($a[0], $a[1]);
      my $line2 = $ifh->getline;
      my @b = split( /\t/, $line2 );
      my ($mate2_id, $flags2) = ($b[0], $b[1]);
      
      # Look for half-mapping alignments where one of the two mates has no
      # reported alignments (flags value = 4 for one mate but not for both mates)
      if ( ($flags1 == 4) ^ ($flags2 == 4) ) {
        # Record the unaligned mate along with the position where its paired mate aligned
        if ($flags1 == 4) {
          # Mate 1 is unaligned, so record mate 1 with the alignment position of mate 2, as in Bowtie 2 SAM output:
          # For field 2, use 69 if other mate aligns to forward strand, use 101 otherwise
          # For field 7, record the name of reference seq where mate's alignment occurs
          # For field 8, record the offset into the fwd ref strand where mate alignment occurs
          my $newflags;
          if ($flags2 == 16) {
            $a[1] = 101;
          } else {
            $a[1] = 69;
          }
          print $ofh "$a[0]\t$a[1]\t$b[2]\t$b[3]\t$a[4]\t$a[5]\t$b[2]\t$b[3]\t$a[8]\t$a[9]\t$a[10]\t$a[11]";
        } else {
          if ($flags1 == 16) {
            $b[1] = 101;
          } else {
            $b[1] = 69;
          }
          # Mate 2 is unaligned, so record mate 2 with the alignment position of mate 1
          print $ofh "$b[0]\t$b[1]\t$a[2]\t$a[3]\t$b[4]\t$b[5]\t$a[2]\t$a[3]\t$b[8]\t$b[9]\t$b[10]\t$b[11]";
        }
        $count++;
        print $ofh2 "$line1$line2";
      } 
    }
    $ifh->close;
    $ofh->close;
    $ofh2->close;
    unlink $bowtie1_infile or die "Can't delete $bowtie1_infile: $!";
    capture( "sort -V -k1,1 -o $outfile2 $outfile4" );
    die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );
    unlink $outfile4 or die "Can't delete $outfile4: $!";
  }
  if ($count == 0) {
    print STDERR "No half-mapping read pairs found.\n";
    exit;
  }
  print "Identified $count half-mapping read pair alignments.\n";
}

=head2 bwa_identify

 Title   : bwa_identify
 Usage   : $project->bwa_identify()
 Function: Identifies half-mapping read pairs from the read alignment output file
 Example : $project->mapping_setup( $index_dir );
           $project->build_bwa_index( $index_dir );
           $project->bwa_identify();
 Returns : No return value
 Args    : Scalar path to alignment file to use (optional)

=cut

sub bwa_identify {
  my $self                    = shift;
  my $alignment_file          = shift || $self->{"alignment_file"};
  my $reads_basename          = $self->{"reads"}->{"basename"};
  my $scripts_dir             = $self->{"scripts_dir"};
  my $work_dir                = $self->{"work_dir"};
  my $outfile                 = "$work_dir/${reads_basename}_halfm_pairs.sam";
  my $outfile2                = "$work_dir/${reads_basename}_halfm_unal_mates_unsorted.sam";
  my $outfile4                = "$work_dir/${reads_basename}_halfm_unal_mates.sam";
  $self->{"halfmapping_file"} = $outfile4;
  
  if (!-e $alignment_file) {
    print STDERR "Alignment file not found.\n";
    return;
  }
  $scripts_dir =~ s!/*$!/! if ($scripts_dir ne ""); # Add trailing slash if not already present
  print "Identifying half-mapping read pairs...\n";
  my $count = 0;
  my $ifh  = IO::File->new( $alignment_file, 'r' ) or die "Can't open $alignment_file: $!";
  my $ofh  = IO::File->new( $outfile, 'w' ) or die "Can't create $outfile: $!";
  my $ofh2 = IO::File->new( $outfile2, 'w' ) or die "Can't create $outfile2: $!";
  my $prev_id         = "";
  my $discard_this_id = 0;
  my $mapping         = 0;
  my @print_lines = ();
  while ( my $line1 = $ifh->getline ) {
    next if $line1 =~ m/^@/;
    my @fields1 = split( /\t/, $line1 );
    my ($mate_id1, $flags1) = ($fields1[0], $fields1[1]);

    my $line2 = $ifh->getline;
    my @fields2 = split( /\t/, $line2 );
    my ($mate_id2, $flags2) = ($fields2[0], $fields2[1]);

    if (&_isHalfMappingPair($flags1, $flags2)) {
      print $ofh "$line1$line2";
      if (&_isUnalignedMate($flags1)) {
        print $ofh2 $line1;
      } else {
        print $ofh2 $line2;
      }
    }
  }
  $ifh->close;
  $ofh->close;
  $ofh2->close;
  
  capture( "sort -V -k1,1 -o $outfile4 $outfile2" );
  die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );
  unlink $outfile2 or die "Can't delete $outfile2: $!";
}  

=head2 set_fragment_length

 Title   : set_fragment_length
 Usage   : $project->set_fragment_length()
 Function: Sets fragment length value for use in filtering and grouping pipeline steps. Distance will
           be measured from the given alignment file if one is provided and a length of -1 is used.
 Example : $project->set_fragment_length( -1, "existing_alignment_file.sam" );
 Returns : No return value
 Args    : Fragment length: outer distance between mates, or -1 as a place holder, and path to an
           existing alignment file to use for inferring the median fragment length

=cut

sub set_fragment_length {
  my $self = shift;
  my $fragment_length = shift;
  my $alignment_file = shift || $self->{"alignment_file"};
  my $scripts_dir = $self->{"scripts_dir"};
  
  return if (defined $self->{"frag_length"});
  
  if (!defined $fragment_length || $fragment_length == -1) {
    print "Determining fragment length...\n";
    $fragment_length = &_compute_frag_length( $scripts_dir, $alignment_file );
    $fragment_length =~ s/\s+$//;
    print "Median fragment length: $fragment_length nt\n";
  }
  $self->{"frag_length"} = $fragment_length;
}
  
=head2 create_simulated_pairs

 Title   : create_simulated_pairs
 Usage   : $project->create_simulated_pairs()
 Function: Create simulated paired-end reads from the outer portions of non-aligning mate sequences,
           using the half-mapping read pairs as a starting point
 Example : $project->run_bwa_mapping();
           $project->bwa_identify();
           $project->create_simulated_pairs();
 Returns : No return value
 Args    : Scalar path to SAM format alignment file (optional), scalar path to SAM half-mapping
           alignment file (optional). Note: Specifying arg 2 requires arg 1

=cut

sub create_simulated_pairs {
  my $self           = shift;
  my $alignment_file = shift || $self->{"alignment_file"};
  my $infile         = shift || $self->{"halfmapping_file"};
  my $work_dir       = $self->{"work_dir"};
  my $reads_basename = $self->{"reads"}->{"basename"};
  my $frag_length    = $self->{"frag_length"};
  my $outfile1       = "$work_dir/${reads_basename}_sim_pairs_1.fq";
  my $outfile2       = "$work_dir/${reads_basename}_sim_pairs_2.fq";
  $self->{"sim_pairs_file_1"} = $outfile1;
  $self->{"sim_pairs_file_2"} = $outfile2;
  
  return if (-z $infile || !(-e $infile));
  print "Creating simulated paired-end reads...\n";

  # Define length of each mate for the simulated read pairs
  my $mate_length = 16;

  my $ifh  = IO::File->new( $infile, 'r' ) or die "$0: Can't open $infile: $!";
  my $ofh1 = IO::File->new( $outfile1, 'w' ) or die "$0: Can't open $outfile1: $!";
  my $ofh2 = IO::File->new( $outfile2, 'w' ) or die "$0: Can't open $outfile2: $!";

  # Boolean to track whether we have any half-mapping pairs
  my $have_candidates = 0;

  # Get the mates we want and write them as simulated paired-end reads to FastQ
  while ( my $line = $ifh->getline ) {
    next if $line =~ m/^@/;
    my @fields = split( /\t/, $line );
    my ($id, $flags, $sequence, $quality_scores) = ($fields[0], $fields[1], $fields[9], $fields[10]);
    $quality_scores =~ s/^\s+|\s+$//g; # Remove leading/trailing whitespace
    
    # Simulate a read pair from the unaligned mate in the pair
    if ( &_isUnalignedMate($flags) ) {

      # Set a flag indicating that we got at least one half-mapping pair
      $have_candidates = 1;

      # Reverse complement reads that aligned to the reverse strand, to get original read
      if ( ( $flags & 16 ) == 1 ) {
        $sequence = &_reverseComplement($sequence);
      }

      # Create the simulated mates. Mate 2 gets reverse complemented.
      my $mate_1 = substr( $sequence, 0, $mate_length );
      my $mate_2 =
        &_reverseComplement( substr( $sequence, -1 * $mate_length ) );
      my $quality_scores_1 = substr( $quality_scores, 0, $mate_length );
      my $quality_scores_2 =
        scalar reverse substr( $quality_scores, -1 * $mate_length );

      # Write the simulated mates
      print $ofh1 "\@$id\/1\n$mate_1\n+\n$quality_scores_1\n";
      print $ofh2 "\@$id\/2\n$mate_2\n+\n$quality_scores_2\n";
    }
  }
  $ifh->close;
  $ofh1->close;
  $ofh2->close;

  print STDERR "No half-mapping read pairs found.\n" if ( !$have_candidates );
}

=head2 align_simulated_pairs_bowtie

 Title   : align_simulated_pairs_bowtie
 Usage   : $project->align_simulated_pairs_bowtie( $num_threads )
 Function: Aligns simulated read pairs to the reference genome using Bowtie.
 Example : $project->bowtie_identify();
           $project->align_simulated_pairs_bowtie( 8 );
 Returns : No return value
 Args    : Number of parallel search threads to use for Bowtie

=cut

sub align_simulated_pairs_bowtie {
  my $self                = shift;
  my $num_threads         = shift;
  my $minins              = shift;
  my $maxins              = shift;
  my $pairs_file_1        = shift || $self->{"sim_pairs_file_1"};
  my $pairs_file_2        = shift || $self->{"sim_pairs_file_2"};
  my $work_dir            = $self->{"work_dir"};
  my $bowtie_version      = $self->{"bowtie_version"};
  my $reads_basename      = $self->{"reads"}->{"basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $index_dir           = $self->{"index_dir"};
  my $outfile             = "$work_dir/${reads_basename}_sim_pairs_alignment.sam";
  my $unsorted_outfile    = "$work_dir/${reads_basename}_sim_pairs_alignment_unsorted.sam";
  $self->{"realignment_file"} = $outfile;
  
  return if (-z $pairs_file_1 || !(-e $pairs_file_1));
 
  # Run Bowtie with the simulated read pairs as input.
  if ($bowtie_version == 2) {
    # Bowtie 2
    print "Aligning simulated paired-end reads using Bowtie 2...\n";
    capture(
          "bowtie2 -x $index_dir/$ref_genome_basename "
        . "--threads $num_threads --reorder --no-hd "
        . "--no-mixed --no-discordant --no-contain --no-overlap "
        . "-1 $pairs_file_1 -2 $pairs_file_2 "
        . "-S $outfile 2> $work_dir/${reads_basename}_bowtie2_sim_pairs_stats.txt"
    );
    die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
  } else {
    # Bowtie 1
    print "Aligning simulated paired-end reads using Bowtie 1...\n";
    capture(
          "bowtie $index_dir/$ref_genome_basename "
        . "--threads $num_threads "
	. "--maxins $maxins --minins $minins "
        . "-m 1 -1 $pairs_file_1 -2 $pairs_file_2 "
        . "-S $unsorted_outfile 2> $work_dir/${reads_basename}_bowtie1_sim_pairs_stats.txt"
    );
    die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
    
    # Sort alignment file by ID
    capture( "sort -V -k1,1 -o $outfile $unsorted_outfile" );
    die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );
    unlink $unsorted_outfile or die "$0: Can't delete $unsorted_outfile: $!";
  }
}

=head2 align_simulated_pairs_bwa

 Title   : align_simulated_pairs_bwa
 Usage   : $project->align_simulated_pairs_bwa( $num_threads )
 Function: Aligns simulated read pairs to the reference genome using BWA.
 Example : $project->bwa_identify();
           $project->align_simulated_pairs_bwa();
 Returns : No return value
 Args    : None

=cut

sub align_simulated_pairs_bwa {
  my $self                 = shift;
  my $num_threads          = shift;
  my $minins               = shift;
  my $maxins               = shift;
  my $pairs_file_1         = shift || $self->{"sim_pairs_file_1"};
  my $pairs_file_2         = shift || $self->{"sim_pairs_file_2"};
  my $work_dir             = $self->{"work_dir"};
  my $reads_basename       = $self->{"reads"}->{"basename"};
  my $ref_genome_basename  = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_extension = $self->{"ref_genome"}->{"extension"};
  my $index_dir            = $self->{"index_dir"};
  my $outfile              = "$work_dir/${reads_basename}_sim_pairs_alignment.sam";
  my $unsorted_outfile     = "$work_dir/${reads_basename}_sim_pairs_alignment_unsorted.sam";
  $self->{"realignment_file"} = $outfile;  
  
  return if (-z $pairs_file_1 || !(-e $pairs_file_1));
  print "Aligning simulated paired-end reads using BWA...\n";
  
  # Create the index for the simulated pairs
  my $copied_genome = "$index_dir/$ref_genome_basename$ref_genome_extension";
  my $sai1 = $pairs_file_1;
  $sai1 =~ s/(\.[^.]+)$/.sai/;
  my $sai2 = $pairs_file_2;
  $sai2 =~ s/(\.[^.]+)$/.sai/;
  capture( "bwa aln $copied_genome $pairs_file_1 > $sai1 2> $sai1.err" );
  die "$0: bwa aln exited unsuccessful" if ( $EXITVAL != 0 );
  capture( "bwa aln $copied_genome $pairs_file_2 > $sai2 2> $sai2.err" );
  die "$0: bwa aln exited unsuccessful" if ( $EXITVAL != 0 );
  
  # Align the simulated pairs to the reference genome
  capture("bwa sampe -s $copied_genome $sai1 $sai2 $pairs_file_1 $pairs_file_2 > $unsorted_outfile " .
          "2> $work_dir/${reads_basename}_bwa_sim_pairs_mapping_stats.txt"
  );
  die "$0: BWA exited unsuccessful" if ( $EXITVAL != 0 );
  
  # Sort alignment file by ID
  capture( "sort -V -k1,1 -o $outfile $unsorted_outfile" );
  die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );
  unlink $unsorted_outfile or die "$0: Can't delete $unsorted_outfile: $!";
}

=head2 filter1

 Title   : filter1
 Usage   : $project->filter1();
 Function: Creates a filtered list of alignments, removing half-mapping read pairs whose unaligned
           mates were successfully aligned to the reference genome as a simulated paired-end read
 Example : $project->bwa_identify();
           $project->align_simulated_pairs_bwa();
           $project->filter1( 5 );
 Returns : No return value
 Args    : Alignment tolerance, path to file containing original half-mapping reads (optional),
           path to file containing alignments for simulated paired-end reads

=cut

sub filter1 {
  my $self                = shift;
  my $tolerance           = shift;
  my $infile              = shift || $self->{"halfmapping_file"};
  my $realignment_file    = shift || $self->{"realignment_file"};
  my $work_dir            = $self->{"work_dir"};
  my $reads_basename      = $self->{"reads"}->{"basename"};
  my $frag_length         = $self->{"frag_length"};
  my $debug_file          = "$work_dir/${reads_basename}_sim_pairs_alignment.debug" if DEBUG;
  my $outfile             = "$work_dir/${reads_basename}_filtered.sam";
  
  # Save the filtered results in the place of the half-mapping alignment file
  $self->{"halfmapping_file"} = $outfile;

  return if (-z $realignment_file || !(-e $realignment_file));
  
  print "Filtering using simulated pair alignment...\n";

  my $ifh1 = IO::File->new( $realignment_file, 'r' ) or die "$0: Can't open $realignment_file: $!";
  my $ifh2 = IO::File->new( $infile, 'r' ) or die "$0: Can't open $infile: $!";
  my $ofh1 = IO::File->new( $outfile, 'w' ) or die "$0: Can't open $outfile: $!";
  my $ofh2 = IO::File->new( $debug_file, 'w' ) or die "$0: Can't open $debug_file: $!" if DEBUG;

  my $read_count = 0;
  my $discard_count = 0;

  # Get an alignment from the originally half-mapping pairs
  while ( my $line = $ifh2->getline ) {
    next if $line =~ m/^@/;
    my @f = split( /\t/, $line );
    my ($orig_id, $orig_flags, $other_mate_chr, $other_mate_pos, $seq) = ($f[0], $f[1], $f[2], $f[7], $f[9]);
    
    # Consider only the unaligned mate from the half-mapping read pair
    next if !&_isUnalignedMate($orig_flags);
    $read_count++;
    
    # Get the corresponding alignment from the simulated pairs
    my $sim_line;
    while ( $sim_line = $ifh1->getline ) {
      next if $sim_line =~ m/^@/;
      my @f = split( /\t/, $sim_line );
      my $id = $f[0];
      while ( $id ne $orig_id ) {
        $sim_line  = $ifh1->getline;
        @f = split( /\t/, $sim_line );
        $id = $f[0];
      }
      my ($chr, $pos) = ($f[2], $f[3]);

      print $ofh2 "Original half-mapping pair unaligned mate: $line" if DEBUG;
      print $ofh2 "Simulated pair alignment line 1: $sim_line"       if DEBUG;

      # Calculate the window for an acceptable alignment, based on which strand we expect to align to
      my ($max_offset, $min_offset);
      if ( $pos != 0 && !&_isOtherMateAlignedToReverseStrand($orig_flags)) {
        $min_offset = $other_mate_pos + $frag_length - $tolerance;
        $max_offset = $other_mate_pos + $frag_length + $tolerance;
      } elsif ( $pos != 0 ) {
        $min_offset = $other_mate_pos - $frag_length - $tolerance;
        $max_offset = $other_mate_pos - $frag_length + $tolerance;        
      }

      # Discard simulated pairs that align close to the originally aligning mate in the half-mapping read pair      
      if ( $chr eq $other_mate_chr && $pos != 0 && $pos < $max_offset && $pos > $min_offset ) {
        print $ofh2 "DISCARDED because $pos between $min_offset and $max_offset, and $chr=$other_mate_chr\n\n" if DEBUG;
        $discard_count++;
        last;
      }

      # Print debug info about why the simulated pair was retained
      if ($pos == 0) {
        print $ofh2 "RETAINED because unaligned\n\n" if DEBUG;
      } else {
        print $ofh2 "RETAINED: other_mate_pos = $other_mate_pos; pos is $pos, which is not between $min_offset and $max_offset, or $chr!=$other_mate_chr\n\n" if DEBUG;
      }

      # Keep this mate
      print $ofh1 $line;
      last;
    }
  }
  $ifh1->close;
  $ifh2->close;
  $ofh1->close;
  $ofh2->close if DEBUG;
  print "Discarded $discard_count mates out of $read_count\n";
}

=head2 build_blast_index

 Title   : build_blast_index
 Usage   : $project->build_blast_index()
 Function: Creates the Blast index needed for running a Blast alignment.
 Example : $project->bwa_identify();
           $project->align_simulated_pairs_bwa();
           $project->filter1( 5 );
           $project->build_blast_index();
           $project->filter2( 500 );
 Returns : No return value
 Args    : Directory to use for the index (optional)

=cut

sub build_blast_index() {
  my $self                = shift;
  my $index_dir           = shift || $self->{"index_dir"};
  my $ref_genome          = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $reads_basename      = $self->{"reads"}->{"basename"};
  unless ( -e "$index_dir/$ref_genome_basename.nhr"
    && "$index_dir/$ref_genome_basename.nin"
    && "$index_dir/$ref_genome_basename.nsq" )
  {
    # Call formatdb to create the Blast index
    print "Creating Blast index in $index_dir...\n";
    capture( "formatdb -p F -i $ref_genome -n $index_dir/$ref_genome_basename" );
    die "$0: formatdb exited unsuccessful" if ( $EXITVAL != 0 );
  }
}

=head2 filter2

 Title   : filter2
 Usage   : $project->filter2();
 Function: Creates a filtered list of alignments, removing half-mapping read pairs whose unaligned
           mates were successfully aligned to the reference genome meeting specific criteria
           using Blast
 Example : $project->bwa_identify();
           $project->align_simulated_pairs_bwa();
           $project->filter1( 5 );
           $project->filter2( 500 );
 Returns : No return value
 Args    : Alignment tolerance

=cut

sub filter2 {
  my $self                = shift;
  my $tolerance           = shift;
  my $work_dir            = $self->{"work_dir"};
  my $index_dir           = $self->{"index_dir"};
  my $scripts_dir         = $self->{"scripts_dir"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $reads_basename      = $self->{"reads"}->{"basename"};
  my $frag_length         = $self->{"frag_length"};  
  my $index               = "$index_dir/${ref_genome_basename}";
  my $halfmap_file        = $self->{"halfmapping_file"};
  my $infasta             = "$work_dir/${reads_basename}_blast_input.fa";
  my $blastfile           = "$work_dir/${reads_basename}_halfm_unal_mates.blast";
  my $errfile             = "$work_dir/${reads_basename}_halfm_unal_mates.blast.err";
  my $hitfile             = "$work_dir/${reads_basename}_halfm_unal_mates_to_remove.txt";
  my $hitfile_sorted      = "$work_dir/${reads_basename}_halfm_unal_mates_to_remove_sorted.txt";
  my $debug_file          = "$work_dir/${reads_basename}_halfm_unal_mates_to_remove.debug" if DEBUG;
  my $outfile             = "$work_dir/${reads_basename}_filtered2.sam";
  
  return if (-z $halfmap_file || !(-e $halfmap_file));
  print "Filtering using Blast...\n";
  
  # Convert SAM to Fasta for Blast
  capture( "${scripts_dir}sam2fasta.pl --rev -i $halfmap_file > $infasta" );
  die "$0: sam2fasta.pl exited unsuccessful" if ( $EXITVAL != 0 );

  unless ( -e "$index_dir/$ref_genome_basename.nhr"
    && "$index_dir/$ref_genome_basename.nin"
    && "$index_dir/$ref_genome_basename.nsq" )
  {
    &build_blast_index($self);
  }
  
  # Call Blast to run the mapping
  capture( "blastall -m 8 -i $infasta -d $index -p blastn -r 2 -G 2 -E 2 -e0.00001 -b1 " .
           "> $blastfile 2> $errfile" );
  die "$0: Blast exited unsuccessful" if ( $EXITVAL != 0 );
  
  # Parse Blast output
  my $last_id = "";
  my $mapped_count = 0;
  my $ifh = IO::File->new( $blastfile, 'r' ) or die "$0: Can't open $blastfile: $!";
  my $ofh = IO::File->new( $hitfile, 'w' ) or die "$0: Can't open $hitfile: $!";
  my $ofh2 = IO::File->new( $debug_file, 'w' ) or die "$0: Can't open $debug_file: $!" if DEBUG;
  while ( my $line = $ifh->getline ) {
    # Retrieve the first match per ID
    my @f = split( /\t/, $line );
    my ($id, $chr, $qstart, $qend, $sstart, $send) = ($f[0], $f[1], $f[6], $f[7], $f[8], $f[9]);
    next if ($id eq $last_id);
    $last_id = $id;
    my ($qid, $qflags, $qchr, $qpos, $qlen) = split( ',', $id);
    
    # Look for matches within 5 of the ends, on the same chromosome
    next if ($chr ne $qchr || $qstart > 5 || $qend < ($qlen - 5));
    
    # Check orientation and position
    my ($max_offset, $min_offset);
    my $strand = 0;
    $strand = 1 if (($send - $sstart) > 0);
    if ( !&_isOtherMateAlignedToReverseStrand($qflags)) {
      # Other mate aligns to forward strand; this mate should align to reverse strand
      if (($send - $sstart) > 0) {
        #print $ofh2 "strand mismatch\n\n" if DEBUG;
        next;
      }
      $min_offset = $qpos + $frag_length - $tolerance;
      $max_offset = $qpos + $frag_length + $tolerance;     
      if ($send <= $max_offset && $send >= $min_offset) {
        print $ofh2 "SUCCESSFULLY MAPPED: " if DEBUG;
        print $ofh "$qid\n";
        $mapped_count++;
      }
      print $ofh2 "$qid $qflags $qchr $qlen $qstart $qend QPOS=$qpos " if DEBUG;
      print $ofh2 "aligned SSTART=$sstart to SEND=$send (on reverse strand), " if DEBUG;
      print $ofh2 "window is MIN_OFFSET=$min_offset to MAX_OFFSET=$max_offset\n\n" if DEBUG;
    } else {
      # Other mate aligns to reverse strand; this mate should align to forward strand
      if (($send - $sstart) < 0) {
        #print $ofh2 "strand mismatch\n\n" if DEBUG;
        next;
      }   
      $min_offset = $qpos - $frag_length - $tolerance;
      $max_offset = $qpos - $frag_length + $tolerance;
      if ($sstart <= $max_offset && $sstart >= $min_offset) {
        print $ofh2 "SUCCESSFULLY MAPPED: " if DEBUG;
        print $ofh "$qid\n";
        $mapped_count++;
      }
      print $ofh2 "$qid $qflags $qchr $qlen $qstart $qend QPOS=$qpos " if DEBUG;
      print $ofh2 "aligned SSTART=$sstart to SEND=$send (on forward strand), " if DEBUG;
      print $ofh2 "window is MIN_OFFSET=$min_offset, MAX_OFFSET=$max_offset\n\n" if DEBUG;
    }
  }
  $ifh->close;
  $ofh->close;
  $ofh2->close if DEBUG;

  # If we're not removing any mates, stop here
  if (-z $hitfile) {
    print "No mates discarded\n";
    return;
  }
  
  # Sort the list of successfully mapped IDs
  capture( "sort -V -k1,1 -o $hitfile_sorted $hitfile" );
  die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 ); 
  
  # Remove mapped IDs from half-mapping list
  $ifh = IO::File->new( $hitfile_sorted, 'r' ) or die "$0: Can't open $hitfile_sorted: $!";
  my $ifh2 = IO::File->new( $halfmap_file, 'r' ) or die "$0: Can't open $halfmap_file: $!";
  $ofh = IO::File->new( $outfile, 'w' ) or die "$0: Can't open $outfile: $!";
  while ( my $line = $ifh->getline ) {
    $line =~ s/^\s+|\s+$//g; # Remove leading/trailing whitespace
    next if $line =~ m/^@/;
    my $id = "";
    while ($line ne $id) {
      my $orig_line = $ifh2->getline;
      my @f = split( /\t/, $orig_line );
      $id = $f[0];
      if ($line eq $id) {
        next;
      } else {
        print $ofh $orig_line;
      }
    }
    
  }
  $ifh->close;
  $ifh2->close;
  $ofh->close;
  
  # Print statistics about how many half-mapping read pairs were removed
  print "Discarded $mapped_count mates\n";
  
  # Save the filtered results in the place of the half-mapping alignment file
  $self->{"halfmapping_file"} = $outfile;
}

=head2 assemble_groups

 Title   : assemble_groups
 Usage   : $project->assemble_groups()
 Function: Identifies groups of half-mapping read pairs that align to the reference genome at
           locations near one another, and performs a local assembly on the unaligned mates in each
           group using Taipan
 Example : $project->bwa_identify();
           $project->filter( 8 );
           $project->assemble_groups( 250, 3, 16 );
 Returns : No return value
 Args    : Expected intron length for group identification, Minimum number of unaligned mates
           from half-mapping read pairs needed near one another to consider them to be a group,
           path to file containing half-mapping reads (optional)

=cut

sub assemble_groups {
  my $self                = shift;
  my $intron_length       = shift;
  my $min_mates           = shift;
  my $min_contig_length   = shift;
  my $halfmap_file        = shift || $self->{"halfmapping_file"};
  my $sorted_file         = $halfmap_file;
  $sorted_file            =~ s/(\.[^.]+)$/_sorted.sam/;
  my $work_dir            = $self->{"work_dir"};
  my $reads_basename      = $self->{"reads"}->{"basename"};
  my $frag_length         = $self->{"frag_length"};
  my $assembly_dir        = $self->{"assembly_dir"};
  my $results_file        = "$work_dir/${reads_basename}_contigs.fa";
  $self->{"contigs_file"} = $results_file;
  my $contigs_file_suffix = "_k16_t8_o1.contigs";
  
  return if (-z $halfmap_file || !(-e $halfmap_file));
  print "Performing local assemblies...\n";
  
  # Sort by alignment position
  capture( "sort -k3,3 -k8n,8 -o $sorted_file $halfmap_file" );
  die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );

  # Get the groups one at a time
  my $ifh = IO::File->new( $sorted_file, 'r' ) or die "$0: $sorted_file: $!";
  my $ofh3 = IO::File->new( $results_file, 'w' ) or die "$0: $results_file: $!";
  my @groups = ();
  my $last_pos = 0;
  my $last_chr = "";
  my $last_len = 0;
  my $overlap_dist;
  my $count = 0;
  my $new_group = 1;
  my $left_pos = 0;
  my $unique = 0;
  while ( my $line = $ifh->getline ) {
    next if $line =~ m/^@/;
    my @fields = split( /\t/, $line );
    my ($id, $flags, $chr, $pos, $seq) = ($fields[0], int($fields[1]), $fields[2], $fields[7], $fields[9]);

    # Calculate overlap distance for this alignment as distance = insert_length + intron_length
    $overlap_dist = $frag_length + $intron_length;

    # Record the first position in the group
    if ($new_group == 1) {
      $left_pos = $pos;
      $new_group = 0;
    }

    # Add or discard based on position, treating upstream mates differently from downstream mates
    my $reversed;
    if ( (($flags & 8) == 0) && (($flags & 0x20 ) == 0) ) { # 69 or 133
      if ( (( $pos - ($frag_length - length($seq)) - $last_len ) - $last_pos) < $overlap_dist
          && (($chr eq $last_chr) || $last_chr eq "" )) {
        push( @groups, $line );
        $last_pos = $pos - ($frag_length - length($seq));
        $last_chr = $chr;
        $last_len = length($seq);
        next;
      }
      $reversed = 1;
    } else { # 101 or 165
      if ( (( $pos - $last_pos ) - $last_len) <= $overlap_dist
          && (($chr eq $last_chr) || $last_chr eq "" )) {
        push( @groups, $line );
        $last_pos = $pos;
        $last_chr = $chr;
        $last_len = length($seq);
        next;
      }
      $reversed = 0;
    }

    # Assemble groups containing at least $min_mates reads
    if ( scalar(@groups) >= $min_mates ) {
      ++$count;
      my $samfile = "$work_dir/assembly/${reads_basename}_group$count.sam";
      my $rawfile = "$work_dir/assembly/${reads_basename}_group$count.raw";
      my $ofh = IO::File->new( $samfile, 'w' ) or die "$0: Can't open $samfile: $!";
      my $ofh2 = IO::File->new( $rawfile, 'w' ) or die "$0: Can't open $rawfile: $!";
      while ( scalar(@groups) > 0 ) {
        my $ln = shift(@groups);
        my @fields = split( /\t/, $ln );
        $chr = $fields[2];
        print $ofh $ln; # Record the SAM alignment data to a SAM file for this group
        print $ofh2 "$fields[9]\n"; # Record the raw read to a .raw file for this group
      }
      $ofh->close;
      $ofh2->close;

      capture( "taipan -f $rawfile -c $min_contig_length -k 16 -o 1 -t 8"
           . " >> $work_dir/${reads_basename}_taipan_output.txt 2>&1" );
      die "$0: taipan exited unsuccessful" if ( $EXITVAL != 0 );

      # Append this group's results to the multi-FastA output file containing all contigs
      my $contigs_file = "${rawfile}${contigs_file_suffix}";
      if (!(-z $contigs_file)) {
        my $ifh2 = IO::File->new( $contigs_file, 'r' ) or die "$0: Can't open $contigs_file: $!";
        while ( my $contig_line = $ifh2->getline ) {
          if ($contig_line =~ m/^>/) {
            $unique++;
            print $ofh3 ">Group${count}-${unique}|$chr|${left_pos}|${last_pos}|" . substr($contig_line, 1);
          } else {
            print $ofh3 $contig_line;
          }
        }
        $ifh2->close;
      }
    } else {
      @groups = ();
    }
    $new_group = 1;
    $last_chr = $chr;
    $last_len = length($seq);
    if ($reversed == 1) {
      $last_pos = $pos - ($frag_length - length($seq));
    } else {
      $last_pos = $pos;
    }
  }
  $ifh->close;
  
  # Handle the last group
  if ( scalar(@groups) >= $min_mates ) {
    ++$count;
    my $samfile = "$work_dir/assembly/${reads_basename}_group$count.sam";
    my $rawfile = "$work_dir/assembly/${reads_basename}_group$count.raw";
    my $ofh = IO::File->new( $samfile, 'w' ) or die "$0: Can't open $samfile: $!";
    my $ofh2 = IO::File->new( $rawfile, 'w' ) or die "$0: Can't open $rawfile: $!";
    my $chr;
    while ( scalar(@groups) > 0 ) {
      my $ln = shift(@groups);
      my @fields = split( /\t/, $ln );
      $chr = $fields[2];
      print $ofh $ln; # Record the SAM alignment data to a SAM file for this group
      print $ofh2 "$fields[9]\n"; # Record the raw read to a .raw file for this group
    }
    $ofh->close;
    $ofh2->close;

    capture( "taipan -f $rawfile -c $min_contig_length -k 16 -o 1 -t 8"
           . " >> $work_dir/${reads_basename}_taipan_output.txt 2>&1" );
    die "$0: taipan exited unsuccessful" if ( $EXITVAL != 0 );

    # Append this group's results to the multi-FastA output file containing all contigs
    my $contigs_file = "${rawfile}${contigs_file_suffix}";
    if (!(-z $contigs_file)) {
      my $ifh2 = IO::File->new( $contigs_file, 'r' ) or die "$0: Can't open $contigs_file: $!";
      while ( my $contig_line = $ifh2->getline ) {
        if ($contig_line =~ m/^>/) {
          $unique++;
          print $ofh3 ">Group${count}|${unique}|$chr|${left_pos}|${last_pos}|" . substr($contig_line, 1);
        } else {
          print $ofh3 $contig_line;
        }
      }
      $ifh2->close;
    }
  }
  if ($count == 0) {
    print "No groups identified for assembly.\n";
    return;
  }
  if ($unique > 0) {
    print "Successfully assembled $unique contigs out of $count groups identified\n";
  } else {
    print "No contigs assembled from $count groups identified.\n";
  }
}

=head2 align_contigs_clustal

 Title   : align_contigs_clustal
 Usage   : $project->align_contigs_clustal()
 Function: Aligns contigs assembled from half-mapping reads
 Example : $project->bwa_identify();
           $project->filter( 8 );
           $project->assemble_groups( 250, 3 );
           $project->align_contigs_clustal();
 Returns : No return value
 Args    : Path to multi-FastA format file containing contigs for alignment (optional)

=cut

sub align_contigs_clustal() {
  my $self           = shift;
  my $contigs_file   = shift || $self->{"contigs_file"};
  my $work_dir       = $self->{"work_dir"};
  my $assembly_dir   = $self->{"assembly_dir"};
  my $reads_basename = $self->{"reads"}->{"basename"};
  my $ref_genome     = $self->{"ref_genome"}->{"full_pathname"};
  my $scripts_dir    = $self->{"scripts_dir"};
  my $output_file    = "$work_dir/${reads_basename}_candidates.fa";
  
  return if (-z $contigs_file || !(-e $contigs_file));
  print "Aligning assembled contigs to reference genome using Clustal...\n";

  # Generate the files containing sequences for each group to be sent to Clustal
  my $padding = 200; # Length of additional sequence past bounds to retrieve from reference
  my $chr = "";
  my $group = "";
  my $last_group = "";
  my $header_line = "";
  my $count = 0;
  my @group_lines = ();
  my ($start, $stop); # Left/right bounds for the group
  my $ifh = IO::File->new( $contigs_file, 'r' ) or die "$0: $contigs_file: $!";
  my $seq_file;
  my @seq_files; # All multi-Fasta input files for Clustal
  my $ofh = undef;
  while ( my $line = $ifh->getline ) { 
    # Get the contigs that correspond to each sequence/chromosome of the reference genome
    if ($line =~ m/^>/) {
      $header_line = $line;
      my @fields = split( /\|/, $line );
      $group = $fields[0];
      $chr = $fields[2];
      $start = $fields[3];
      $stop = $fields[4];
      $group =~ s/^\>//g;
      next;
    }
    if ($group ne $last_group) { # Got the next group
      # Define the file for writing the reference genome segment and all contigs for this group
      $seq_file = "$assembly_dir/${reads_basename}_${group}_sequences.fa";
      push(@seq_files, $seq_file);
      
      # Get the portion of the reference genome we want to align to
      $start = $start - $padding;
      $start = 0 if ($start < 0);
      $stop = $stop + $padding;
      capture( "${scripts_dir}trim_fasta.pl -i $ref_genome -s $chr --start " .
              ($start + $padding) . " --stop " . ($stop + $padding) . " > $seq_file");
      die "$0: trim_fasta.pl exited unsuccessful" if ( $EXITVAL != 0 );

      # Write the first sequence for this group
      $ofh->close if (defined $ofh);
      $ofh = IO::File->new( $seq_file, 'a' ) or die "$0: $seq_file: $!";
      print $ofh $header_line . $line;
    } else { # Got another line for this group
      print $ofh $header_line . $line
    }
    $last_group = $group; 
  }
  $ifh->close;
  $ofh->close;

  # Run Clustal for every group
  foreach my $file (@seq_files) {
    my $err_file = $file;
    $err_file =~ s/(\.[^.]+)$/.err/;
    
    # Call Clustal to run the mapping
    capture( "clustalw -infile=$file -gapopen=50 -gapext=0.01 -output=gde > $err_file" );
    die "$0: clustalw exited unsuccessful" if ( $EXITVAL != 0 );
    #capture( "clustalw -infile=$file -gapopen=50 -gapext=0.01" ) if DEBUG;
    #die "$0: clustalw exited unsuccessful" if ( $EXITVAL != 0 && DEBUG);
  }

  # Parse the Clustal results for every group
  my $reg_count = 0;
  foreach my $file (@seq_files) {
    my $gde_file = $file;
    $gde_file =~ s/(\.[^.]+)$/.gde/;
    my $err_file = $file;
    $err_file =~ s/(\.[^.]+)$/.err/;
    my $counts_file = $file if DEBUG;
    $counts_file =~ s/(\.[^.]+)$/.counts/ if DEBUG;
    
    # Use a hash to count the number of aligned sequences at each position
    my %counts = ();
    my $pos = 0;
    my ($chr, $start);
    $ifh = IO::File->new( $gde_file, 'r' ) or die "$0: $gde_file: $!";
    $ofh = IO::File->new( $counts_file, 'w' ) or die "$0: $counts_file: $!" if DEBUG;
    while ( my $line = $ifh->getline ) {
      if ($line =~ m/^#/) {
        if ($line =~ m/^##(\S+)\|(\d+)\|/) { # Trimmed sequence from reference genome
          $chr = $1;
          $start = $2;
        }
        $pos = 0;
        next;
      }
      $line =~ s/\s+$//g;
      foreach my $c (split //, $line) {
        if ($c ne "-") {
          $counts{$pos}++;
        }
        $pos++;
      }
    }
    $ifh->close;
    
    # Convert the number of aligned sequences at each position into a string
    my $count_string = "";
    for (my $i = 0; $i < $pos; $i++) {
      if ($counts{$i} <= 9) {
        $count_string .= "$counts{$i}";
      } else {
        $count_string .= "9"; # Use a single digit for values >9 
      }
    }
    for ( my $i = 0; $i < length($count_string); $i += 20) {
      print $ofh substr($count_string, $i, 20), "\n" if DEBUG;
    }
    $ofh->close if DEBUG;
    
    # Find regions between alignments
    my $saved_sequence = "";
    my $i = $start;
    my $last_c = "1";
    my $aln_start = -1;
    my $aln_stop = -1;
    foreach my $c (split //, $count_string) {
      $i++;
      # Look for the end of an aligned region
      if ($last_c > 1 && $c eq "1") {
        $aln_start = $i;
      }
      
      # Look for the beginning of an aligned region, and a minimum length of 50
      if ($last_c eq "1" && $c > 1 && $aln_start != -1) {
        $aln_stop = $i;
        if ($aln_stop - $aln_start >= 50) {
          $reg_count++;
          # Get the sequence between two alignments directly from the reference genome
          capture("${scripts_dir}trim_fasta.pl -i $ref_genome --start $aln_start " .
                  " -s $chr --stop $aln_stop >> $output_file");
          die "$0: trim_fasta.pl exited unsuccessful" if ( $EXITVAL != 0 );
          $aln_start = -1;
        }
      }
      $last_c = $c;
    }
  }
  print "Identified $reg_count candidate insertion/deletion regions.\n";
}

############ Subroutines for internal use by this module ############

# Runs infer_fraglen.pl script to determine median fragment length from aligned read pairs
sub _compute_frag_length {
  my $scripts_dir    = shift;
  my $alignment_file = shift;
  # Determine length of read pair for use in filtering and grouping steps
  my $frag_length = capture("${scripts_dir}infer_fraglen.pl -i $alignment_file -m");
  die "$0: infer_fraglen.pl exited unsuccessful" if ( $EXITVAL !=0 );
  return $frag_length;
}

# Returns the reverse complement of the given string
sub _reverseComplement {
  my $value = shift;
  $value = scalar reverse $value;
  for ( 0 .. length($value) - 1 ) {
    substr( $value, $_, 1 ) = &_complement( substr( $value, $_, 1 ) );
  }
  return $value;
}

# Returns the complement of the given string
sub _complement {
  my $input = shift;
  my %complementMap = (
    "A" => "T",
    "T" => "A",
    "a" => "t",
    "t" => "a",
    "C" => "G",
    "G" => "C",
    "c" => "g",
    "g" => "c",
    "N" => "N",
    "." => ".",
    "n" => "n"
  );
  my $complementedString = $complementMap{ $input }
    or die "$0: Can't get reverse complement for '$input'";
  return $complementedString;
}

# Returns 1 if the given sum-of-flags value identifies a secondary alignment, otherwise returns 0
sub _isSecondaryAlignment {
  my $flags = int(shift);
  return ( (($flags & 0x100) != 0) );
}

# Returns 1 if the given sum-of-flags value identifies a mate in a half-mapping read pair, otherwise returns 0
# The mate may be either the aligning mate or the non-aligning mate in the half-mapping read pair.
# Note: Also returns 1 for "far-mapping" alignments, in which both mates align independently but do not meet
# fragment length constraints defined by Bowtie's "--minins" and "--maxins" options
sub _isHalfMappingMate {
  my $flags = int(shift);
  return ( (($flags & 4) != 0) ^ (($flags & 8) != 0) );
}

# Returns 1 if the given sum-of-flags value identifies a mate in an aligning read pair, otherwise returns 0
# Note: Also returns 1 for discordant alignments, which must be suppressed using Bowtie's "--no-discordant" option
sub _isInMappingPair {
  my $flags = int(shift);
  return ( (($flags & 4) == 0) && (($flags & 8) == 0) );
}

# Returns 1 if the given sum-of-flags value identifies a mate in a no-alignments read pair, otherwise returns 0
sub _isInNonMappingPair {
  my $flags = int(shift);
  return ( (($flags & 4) != 0) && (($flags & 8) != 0) );
}

# Returns 1 if the given sum-of-flags value identifies a mate in a pair that aligns discordantly, otherwise returns 0
# Discordant alignments are those in which the mates overlap, or both mates map to the same strand
sub _isInDiscordantPair {
  my $flags = int(shift);
  return ( ($flags == 65 || $flags == 81 || $flags == 97 || $flags == 113 || 
            $flags == 129 || $flags == 145 || $flags == 161 || $flags == 177) );
}

# Returns 1 if the given sum-of-flags values represent a pair that aligns discordantly, otherwise returns 0
# Discordant alignments are those in which the mates overlap, or both mates map to the same strand
sub _isDiscordantPair {
  my $flags1 = int(shift);
  my $flags2 = int(shift);
  return ( _isInDiscordantPair($flags1) || _isInDiscordantPair($flags2));
}

# Returns 1 if the given sum-of-flags values represent a half-mapping read pair, otherwise returns 0
# Bowtie 2 and BWA differ on whether they set the 16 (0x10) flag in two half-mapping cases (looks like a bug in BWA)
sub _isHalfMappingPair {
  my $flags1 = int(shift);
  my $flags2 = int(shift);
  return ( ($flags1 == 69 && $flags2 == 137) || ($flags1 == 137 && $flags2 == 69)   # Bowtie 2/BWA
        || ($flags1 == 73 && $flags2 == 133) || ($flags1 == 133 && $flags2 == 73)   # Bowtie 2/BWA
        || ($flags1 == 89 && $flags2 == 165) || ($flags1 == 165 && $flags2 == 89)   # Bowtie 2
        || ($flags1 == 101 && $flags2 == 153) || ($flags1 == 153 && $flags2 == 101) # Bowtie 2
        || ($flags1 == 89 && $flags2 == 181) || ($flags1 == 181 && $flags2 == 89)   # BWA
        || ($flags1 = 117 && $flags2 == 153) || ($flags1 == 153 && $flags2 == 117)  # BWA
         ); 
}

# Returns 1 if the given sum-of-flags values represent a read pair containing mates that align independently but
# do not meet fragment length constraints, otherwise returns 0
sub _isFarMappingPair {
  my $flags1 = int(shift);
  my $flags2 = int(shift);
  return ( ((($flags1 & 4) == 0) && (($flags1 & 8) != 0))
        && ((($flags2 & 4) == 0) && (($flags2 & 8) != 0)) );
}

# Returns 1 if the given sum-of-flags value indicates the other mate aligned to the reverse strand
sub _isOtherMateAlignedToReverseStrand {
  my $flags = int(shift);
  return ( (($flags & 32) != 0) );
}

# Returns 1 if the given sum-of-flags value indicates the mate is mate 1 in the pair, otherwise returns 0
sub _isMateOneInPair {
  my $flags = int(shift);
  return ( (($flags & 64) != 0))
}

# Returns 1 if the given sum-of-flags value identifies a non-aligning mate, otherwise returns 0
sub _isUnalignedMate {
  my $flags = int(shift);
  return ( (($flags & 4) != 0) );
}

# Ensures that the given folder exists, and removes trailing slash from string
sub _check_dir {
  my $dir = shift;
  if ( -d $dir ) {
    $dir = $1 if ( $dir =~ /(.*)\/$/ );
    return $dir;
  }
  die "$0: directory $dir does not exist";
}

# Creates the given directory, and gives an error message on failure
sub _make_dir {
  my $dirname = shift;
  mkdir $dirname, 0755 or die "$0: could not create $dirname";
}

# Returns a date/time string with no whitespace
sub _datestamp {
  my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
    localtime(time);
  my $datestring = sprintf(
    "%4d%02d%02d_%02d_%02d_%02d",
    $year + 1900,
    $mon + 1, $mday, $hour, $min, $sec
  );
}

# _find_mates: Identify corresponding (paired) mates for the given file of mates
# Arguments: 1. Unmatched mates file 2. File with candidate mates to search 3. Output file
# Note: Matches in output file are unsorted, so may be listed in a different order.
sub _find_mates {
  my $unpaired_mates = shift;
  my $candidate_mates = shift;
  my $output_file = shift;

  # Open the unpaired mates file
  my $ifh = new IO::File($unpaired_mates, 'r') or die "$0: Can't open $unpaired_mates: $!";

  # Loop through the fastq unpaired mates file and get all the search IDs
  my $line_count = 0;
  my $line = 0;
  my $combined_search = "";
  while ( my $unpaired_line = $ifh->getline ) {
    my $id;
    $line++;
    if ( $line_count == 0 ) {
      $id = substr( $unpaired_line, 1 );
      chomp( $id );
      $id =~ s/\/[12]$//;

      # Append this ID to the combined search string using alternation metacharacter
      # Combined search string will look like r4|r12|r21|r36|r142|r141|r153
      if ( $combined_search ne "" ) {
        $combined_search = $combined_search . "|$id";
      } else {
        $combined_search = $id;
      }
    }
    $line_count = ($line_count + 1) % 4;
  }
  $ifh->close;

  # Open the file containing all potential mates
  my $ifhb = new IO::File($candidate_mates, 'r') or die "$0: Can't open $candidate_mates: $!";

  # Open the output file for writing
  my $ofh = IO::File->new($output_file, "w") or die "$0: Can't create $output_file: $!";

  # Find and save the entries for all IDs in the combined search string
  $line_count = 0;
  my $save_flag = 0;
  while ( my $candidate_line = $ifhb->getline ) {
    if ( $line_count == 0 ) {
      # Find the pairing mate for that ID in the candidate mates file
      $save_flag = 0;
      if ( $candidate_line =~ m/$combined_search/ && $combined_search ne "") {

        # Set the flag to print all lines for this mate to the output file
        $save_flag = 1;

        # Get the ID that matched
        my $id;
        $id = substr( $candidate_line, 1 );
        chomp( $id );
        $id =~ s/\/[12]$//;

        # Remove the ID that matched from the combined search string
        $combined_search =~ s/\|\Q$id\E\|/\|/g; # Match in middle of line
        $combined_search =~ s/^\Q$id\E\|//g; # Match at beginning of line
        $combined_search =~ s/\|\Q$id\E$//g; # Match at end of line
        $combined_search =~ s/^\Q$id\E$//; # Match whole line
      }
    }
    $line_count = ($line_count + 1) % 4;

    # Write this line to the output file if it's a match
    if ( $save_flag == 1 ) {
      print $ofh $candidate_line;
    }
  }
  $ifhb->close;
  $ofh->close;

  # Make sure we found all the mates we were looking for
  if ( $combined_search ne "") {
    die "$0: mates not found in $candidate_mates: $combined_search\n";
  }
}

# Sorts the entries in a FastQ file by their identifiers using Unix sort
sub _sort_fastq_by_id {
  my $input_file = shift;
  my $output_file = shift;
  my $linearized_file = $input_file . "_linearized.tmp";
  my $linearized_sorted_file = $linearized_file . "_sorted.tmp";
  my $ifh = new IO::File($input_file, 'r') or die "$0: Can't open $input_file: $!";
  my $lfh = IO::File->new($linearized_file, "w") or die "$0: Can't create $linearized_file: $!";

  # Linearize the file by merging multiple lines into one
  my $line_count = 0;
  while ( my $line = $ifh->getline ) {
    chomp($line);
    print $lfh $line;
    $line_count = ($line_count + 1) % 4;
    if ($line_count == 0) {
      print $lfh "\n";
    } else {
      print $lfh "\t";
    }
  }
  $ifh->close;
  $lfh->close;

  # Sort the file
  system( "sort -V -t '\t' -k1,1 $linearized_file -o $linearized_sorted_file" );
  unlink $linearized_file or die "$0: Can't delete $linearized_file: $!";

  # Un-linearize the file by splitting lines into multiple lines, and save
  $lfh = IO::File->new($linearized_sorted_file, "r") or die "$0: Can't open $linearized_sorted_file: $!";
  my $ofh = IO::File->new($output_file, "w") or die "$0: Can't create $output_file: $!";
  $line_count = 0;
  while ( my $line = $lfh->getline ) {
    chomp( $line );
    my @fields = split( /\t/, $line );
    print $ofh $fields[0] . "\n";
    print $ofh $fields[1] . "\n";
    print $ofh $fields[2] . "\n";
    print $ofh $fields[3] . "\n";
  }
  $ofh->close;
  $lfh->close;
  unlink $linearized_sorted_file or die "$0: Can't delete $linearized_sorted_file: $!";
}

1;
