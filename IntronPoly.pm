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

use constant DEBUG => (1);    # Flag to print info for debugging (temporary)

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
  $self->{"reads_basename"} = undef;   # Scalar base filename for read pairs/data
  $self->{"ref_genome"}     = undef;   # Hash of reference genome related values
  $self->{"bowtie_db"}      = undef;   # Hash of data directories for bowtie
  $self->{"alignment_file"} = undef;   # Scalar path to initial read pairs alignments
  $self->{"read_length"}    = undef;   # Median read length determined from alignment
  bless($self);
  return $self;
}

=head2 set_work_dir

 Title   : set_work_dir
 Usage:  : $project->set_work_dir( "path_to_scripts" )
 Function: Creates directory for pipeline working data and output files
 Example : my $work_dir = "/home/theis/work";
           $project->set_work_dir( $scripts_dir )
 Returns : The path to the directory for pipeline working data and output files
 Args    : Scalar of full path to the parent directory of the work directory

=cut

sub set_work_dir {
  my $self        = shift;
  my $scripts_dir = shift;
  $scripts_dir =~ s!/*$!/! if ($scripts_dir ne ""); # Add trailing slash if not already present
  $self->{"scripts_dir"} = $scripts_dir;
  my $work_dir = shift || "${scripts_dir}run-" . &_datestamp;
  unless ( -e $work_dir ) {
    $work_dir = $1 if ( $work_dir =~ /(.*)\/$/ );
    &_make_dir($work_dir);
    print "Created work directory $work_dir\n";
  }
  else {
    $work_dir = &_check_dir($work_dir);
    print "Using existing work directory $work_dir\n";
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
           $project->build_db( $ref_genome, "reads" );
 Returns : No return value
 Args    : Scalar of full path to the reference genome Fasta file, and a scalar of the base
           name of the read pairs files (filename without "_1.fq" or "_2.fq" ending)

=cut

sub build_db {
  my $self                 = shift;
  my $ref_genome_filename  = shift;
  my $reads_basename       = shift;

  # Use the base of the read pairs filenames for saving data throughout pipeline
  $self->{"reads_basename"} = $reads_basename;

  # Split reference genome pathname into components
  my ( $file, $dir, $ext ) = fileparse( $ref_genome_filename, qr/\.[^.]*/ );
  $self->{"ref_genome"}->{"full_pathname"} = $ref_genome_filename;
  $self->{"ref_genome"}->{"basename"}      = $file;
  $self->{"ref_genome"}->{"dir"}           = $dir;
}

=head2 mapping_setup

 Title   : mapping_setup
 Usage   : $project->mapping_setup( "index_dir", "reads_dir", "read_pairs_base_name")
 Function: Sets references to data files used for mapping reads to the reference genome, and validates
           input files
 Example : my $index_dir = "/home/theis/bt2/";
           my $reads_dir = "/home/theis/reads/";
           $project->mapping_setup( $index_dir, $reads_dir, "reads" );
 Returns : No return value
 Args    : Scalar full path to the Bowtie executable directory, scalar full path to the Bowtie
           index directory, scalar full path to the directory containing the Fastq read pairs files,
           and a scalar of the base name of the read pairs files (filename without "_1.fq" or
           "_2.fq" ending)

=cut

sub mapping_setup {
  my $self                 = shift;
  my $index_dir            = shift;
  my $reads_dir            = shift;
  my $reads_basename       = shift;
  my $scripts_dir          = $self->{"scripts_dir"};
  my $work_dir             = $self->{"work_dir"};
  my $ref_genome           = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename  = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_dir       = $self->{"ref_genome"}->{"dir"};
  die "$0: reference genome directory $ref_genome_dir does not exist" unless ( -d $ref_genome_dir );
  
  print "Validating Fasta reference genome data...\n";
  
  # Perform basic validation on reference genome
  capture("${scripts_dir}validate_fasta.pl -i $ref_genome");
  die "$0: validate_fasta.pl exited unsuccessful" if ( $EXITVAL != 0 );

  # Ensure that read pairs files exist
  $reads_dir = $1 if ( $reads_dir =~ /(.*)\/$/ );
  my $reads_file_one = "$reads_dir/${reads_basename}_1.fq";
  my $reads_file_two = "$reads_dir/${reads_basename}_2.fq";
  die "$0: reads file $reads_file_one does not exist" unless ( -e $reads_file_one );
  die "$0: reads file $reads_file_two does not exist" unless ( -e $reads_file_two );

  print "Validating Fastq read pair data...\n";
  
  # Perform basic validation on read pairs files
  capture("${scripts_dir}validate_fastq.pl -i $reads_dir/$reads_basename");
  die "$0: validate_fastq.pl exited unsuccessful" if ( $EXITVAL != 0 );

  # Create directory for Bowtie index files if necessary
  $index_dir = $1 if ( $index_dir =~ /(.*)\/$/ );
  unless ( -d $index_dir ) {
    &_make_dir($index_dir);
  }

  # Remember paths needed for running
  $self->{"bowtie_db"} = {
    reads_file_one => $reads_file_one,
    reads_file_two => $reads_file_two
  };
  $self->{"index_dir"} = $index_dir;
  print "Using reference genome $ref_genome\n";
  print "Using reads file 1: $reads_file_one\n";
  print "Using reads file 2: $reads_file_two\n";
}

=head2 build_bowtie_index

 Title   : build_bowtie_index
 Usage   : $project->build_bowtie_index();
 Function: Runs the bowtie2-build indexer to create a Bowtie index from the reference genome
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
  $self->{"index_dir"}    = $index_dir;
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
  else {
    print "Bowtie index already exists, not re-creating.\n";
    print "Using Bowtie index at $index_dir/$ref_genome_basename.*.bt2\n";
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
  my $reads_basename        = $self->{"reads_basename"};
  my $ref_genome_basename   = $self->{"ref_genome"}->{"basename"};
  my $index_dir             = $self->{"index_dir"};
  my $reads_file_one        = $self->{"bowtie_db"}->{"reads_file_one"};
  my $reads_file_two        = $self->{"bowtie_db"}->{"reads_file_two"};
  my $work_dir              = $self->{"work_dir"};
  my $output_file           = "$work_dir/${reads_basename}_initial_alignments.sam";
  $self->{"alignment_file"} = $output_file;
  print "Running mapping using Bowtie, using $num_threads threads...\n";
  
  # Call Bowtie to run the mapping
  capture(
        "bowtie2 -x $index_dir/$ref_genome_basename "
      . "--threads $num_threads --reorder --no-hd "
      . "--maxins $maxins --minins $minins "
      . "--no-discordant "
      . "--no-contain --no-overlap "
      . "-k 3 -1 $reads_file_one -2 $reads_file_two "
      . "-S $output_file"
  );
  die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
  print "Saved alignment data to $output_file\n";
}

=head2 bowtie_identify

 Title   : bowtie_identify
 Usage   : $project->bowtie_identify()
 Function: Identifies half-mapping read pairs from the read alignment output file
 Example : $project->build_bowtie_index();
           $project->mapping_setup( $index_dir, $reads_dir, "reads" );
           $project->bowtie_identify();
 Returns : No return value
 Args    : Scalar path to alignment file to use (optional)

=cut

sub bowtie_identify {
  my $self                    = shift;
  my $alignment_file          = shift || $self->{"alignment_file"};
  my $reads_basename          = $self->{"reads_basename"};
  my $scripts_dir             = $self->{"scripts_dir"};
  my $work_dir                = $self->{"work_dir"};
  my $outfile                 = "$work_dir/${reads_basename}_halfm_pairs.sam";
  my $outfile2                = "$work_dir/${reads_basename}_halfm_unal_mates.sam";
  my $outfile3                = $outfile;
  $outfile3                   =~ s/(\.[^.]+)$/_sorted.sam/;
  $self->{"halfmapping_file"} = $outfile2;
  
  if (!-e $alignment_file) {
    print STDERR "Alignment file not found.\n";
    return;
  }

  $scripts_dir =~ s!/*$!/! if ($scripts_dir ne ""); # Add trailing slash if not already present

  # Determine length of read pair for use in filtering and grouping steps
  print "Determining fragment length...\n";
  my $fragment_length = &_compute_frag_length( $scripts_dir, $alignment_file );
  $fragment_length =~ s/\s+$//;
  $self->{"frag_length"} = $fragment_length;
  print "Median fragment length: $fragment_length nt\n";

  print "Identifying half-mapping read pairs...\n";
  
  my $ifh  = IO::File->new( $alignment_file, 'r' ) or die "Can't open $alignment_file: $!";
  my $ofh  = IO::File->new( $outfile, 'w' ) or die "Can't create $outfile: $!";
  my $ofh2 = IO::File->new( $outfile2, 'w' ) or die "Can't create $outfile2: $!";
  my $prev_id         = "";
  my $discard_this_id = 0;
  my $mapping         = 0;
  my @print_lines = ();
  my $count = 0;
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
    }
    else {

      # This line is the second or greater line with this read ID
      next if $discard_this_id;

      if ( scalar @print_lines == 2 && $mapping == 0 ) {

        # Found a secondary alignment preceded by two half-mapping mates. Discard all lines with this ID.
        @print_lines     = ();
        $discard_this_id = 1;
        next;
      }
      elsif ( scalar @print_lines == 3 ) {

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
    }
    elsif ( &_isInNonMappingPair($flags) ) {

      # We have a non-mapping ID. Discard all lines with this ID.
      $mapping         = 0;
      $discard_this_id = 1;
    }
    elsif ( &_isHalfMappingMate($flags) ) {

      # Save this half-mapping line for printing
      $mapping = 1;
      push( @print_lines, $line );
    }
    else {

      # We have a secondary alignment.
      $mapping = 2;
      push( @print_lines, $line );
    }
  }
  $ifh->close;
  $ofh->close;
  if ($count == 0) {
    print STDERR "No half-mapping read pairs found.\n";
    return;
  }
  print "Saved $count half-mapping read pair alignments to $outfile\n";
  print "Sorting...\n";
  capture( "sort -k3,3 -k8n,8 -o $outfile3 $outfile" );
  die "$0: sort exited unsuccessful" if ( $EXITVAL != 0 );
}

=head2 create_fake_pairs

 Title   : create_fake_pairs
 Usage   : $project->create_fake_pairs()
 Function: Create simulated paired-end reads from the outer portions of non-aligning mate sequences,
           using the half-mapping read pairs as a starting point
 Example : $project->run_bowtie_mapping( 8, 100, 500 );
           $project->bowtie_identify();
           $project->create_fake_pairs();
 Returns : No return value
 Args    : Scalar path to SAM format alignment file (optional), scalar path to SAM half-mapping
           alignment file (optional). Note: Specifying arg 2 requires arg 1

=cut

sub create_fake_pairs {
  my $self           = shift;
  my $alignment_file = shift || $self->{"alignment_file"};
  my $infile         = shift || $self->{"halfmapping_file"};
  my $work_dir       = $self->{"work_dir"};
  my $reads_basename = $self->{"reads_basename"};
  my $scripts_dir    = $self->{"scripts_dir"};
  my $frag_length    = $self->{"frag_length"};
  my $outfile1       = "$work_dir/${reads_basename}_fake_pairs_1.fq";
  my $outfile2       = "$work_dir/${reads_basename}_fake_pairs_2.fq";
  $self->{"fake_pairs_file_1"} = $outfile1;
  $self->{"fake_pairs_file_2"} = $outfile2;
  
  return if (-z $infile || !(-e $infile));
  print "Creating simulated paired-end reads...\n";

  # Define length of each mate for the fake read pairs
  my $mate_length = 16;

  # Make sure we have computed the fragment length
  $self->{"frag_length"} = &_compute_frag_length( $scripts_dir, $alignment_file ) if ($frag_length eq "");

  my $ifh  = IO::File->new( $infile, 'r' ) or die "$0: Can't open $infile: $!";
  my $ofh1 = IO::File->new( $outfile1, 'w' ) or die "$0: Can't open $outfile1: $!";
  my $ofh2 = IO::File->new( $outfile2, 'w' ) or die "$0: Can't open $outfile2: $!";

  # Boolean to track whether we have any half-mapping pairs
  my $have_candidates = 0;

  # Get the mates we want and write them as fake paired-end reads to FastQ
  while ( my $line = $ifh->getline ) {
    next if $line =~ m/^@/;
    my @fields = split( /\t/, $line );
    my ($id, $flags, $sequence, $quality_scores) = ($fields[0], $fields[1], $fields[9], $fields[10]);

    # Create a fake pair from the unaligned mate in the pair
    if ( &_isUnalignedMate($flags) ) {

      # Set a flag indicating that we got at least one half-mapping pair
      $have_candidates = 1;

      # Reverse complement reads that aligned to the reverse strand, to get original read
      if ( ( $flags & 16 ) == 1 ) {
        $sequence = &_reverseComplement($sequence);
      }

      # Create the fake mates. Mate 2 gets reverse complemented.
      my $mate_1 = substr( $sequence, 0, $mate_length );
      my $mate_2 =
        &_reverseComplement( substr( $sequence, -1 * $mate_length ) );
      my $quality_scores_1 = substr( $quality_scores, 0, $mate_length );
      my $quality_scores_2 =
        scalar reverse substr( $quality_scores, -1 * $mate_length );

      # Write the fake mates
      print $ofh1 "\@$id\/1\n$mate_1\n+\n$quality_scores_1\n";
      print $ofh2 "\@$id\/2\n$mate_2\n+\n$quality_scores_2\n";
    }
  }
  $ifh->close;
  $ofh1->close;
  $ofh2->close;

  print STDERR "No half-mapping read pairs found.\n" if ( !$have_candidates );
}

=head2 run_fake_pairs_alignment

 Title   : run_fake_pairs_alignment
 Usage   : $project->run_fake_pairs_alignment( $num_threads )
 Function: Aligns simulated read pairs to the reference genome.
 Example : $project->bowtie_identify();
           $project->run_fake_pairs_alignment( 8 );
 Returns : No return value
 Args    : Number of parallel search threads to use for Bowtie

=cut

sub run_fake_pairs_alignment {
  my $self                = shift;
  my $num_threads         = shift;
  my $pairs_file_1        = shift || $self->{"fake_pairs_file_1"};
  my $pairs_file_2        = shift || $self->{"fake_pairs_file_2"};
  my $work_dir            = $self->{"work_dir"};
  my $reads_basename      = $self->{"reads_basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $index_dir           = $self->{"index_dir"};
  my $outfile             = "$work_dir/${reads_basename}_fake_pairs_alignment.sam";
  $self->{"realignment_file"} = $outfile;
  
  return if (-z $pairs_file_1 || !(-e $pairs_file_1));
  print "Aligning simulated paired-end reads...\n";
  
  # Run Bowtie with the fake read pairs as input
  capture(
        "bowtie2 -x $index_dir/$ref_genome_basename "
      . "--threads $num_threads --reorder --no-hd "
      . "--no-mixed --no-discordant --no-contain --no-overlap "
      . "-1 $pairs_file_1 -2 $pairs_file_2 "
      . "-S $outfile"
  );
  die "$0: Bowtie exited unsuccessful" if ( $EXITVAL != 0 );
}

=head2 filter1

 Title   : filter1
 Usage   : $project->filter1();
 Function: Creates a filtered list of alignments, removing half-mapping read pairs whose unaligned
           mates were successfully aligned to the reference genome using relaxed mapping criteria
 Example : $project->bowtie_identify();
           $project->run_fake_pairs_alignment( 8 );
           $project->filter1();
 Returns : No return value
 Args    : Path to file containing alignments from a secondary run of Bowtie using relaxed criteria
           (optional), path to file containing original half-mapping reads (optional)
=cut

sub filter1 {
  my $self                = shift;
  my $realignment_file    = shift || $self->{"realignment_file"};
  my $infile              = shift || $self->{"halfmapping_file"};
  my $work_dir            = $self->{"work_dir"};
  my $reads_basename      = $self->{"reads_basename"};
  my $frag_length         = $self->{"frag_length"};
  my $debug_file          = "$work_dir/${reads_basename}_fake_pairs_alignment.debug" if DEBUG;
  my $outfile             = "$work_dir/${reads_basename}_filtered.sam";

  # Save the filtered results in the place of the half-mapping alignment file
  $self->{"halfmapping_file"} = $outfile;

  return if (-z $realignment_file || !(-e $realignment_file));
  
  print "Filtering alignments...\n";

  my $ifh1 = IO::File->new( $realignment_file, 'r' ) or die "$0: Can't open $realignment_file: $!";
  my $ifh2 = IO::File->new( $infile, 'r' ) or die "$0: Can't open $infile: $!";
  my $ofh1 = IO::File->new( $outfile, 'w' ) or die "$0: Can't open $outfile: $!";
  my $ofh2 = IO::File->new( $debug_file, 'w' ) or die "$0: Can't open $debug_file: $!" if DEBUG;

  my $read_count = 0;
  my $discard_count = 0;

  # Get an alignment from the originally half-mapping pairs
  while ( my $line = $ifh2->getline ) {
    next if $line =~ m/^@/;
    my @fields = split( /\t/, $line );
    my ($orig_id, $orig_flags, $other_mate_pos) = ($fields[0], $fields[1], $fields[7]);

    # Consider only the unaligned mate from the half-mapping read pair
    next if !&_isUnalignedMate($orig_flags);
    $read_count++;

    # Get the corresponding alignment from the fake pairs
    my $fake_line;
    while ( $fake_line = $ifh1->getline ) {
      next if $fake_line =~ m/^@/;
      my @fields = split( /\t/, $fake_line );
      my ($id, $pos) = ($fields[0], $fields[3]);
      while ( $id ne $orig_id ) {
        $fake_line = $ifh1->getline;
        my @fields = split( /\t/, $fake_line );
        ($id, $pos) = ($fields[0], $fields[3]);
      }

      print $ofh2 "\nFAKE PAIR ALIGNMENT LINE 1: $fake_line"         if DEBUG;
      print $ofh2 "ORIGINAL HALF-MAPPING PAIR UNALIGNED MATE: $line" if DEBUG;
      print $ofh2 "ORIG_POS: $other_mate_pos POS: $pos\n"            if DEBUG;

      # Discard fake pairs that align close to the originally aligning mate in the half-mapping read pair
      if (
        $pos != 0
        && ( $pos < ( $other_mate_pos + $frag_length + 10 )
          && $pos > ( $other_mate_pos + $frag_length - 10 ) )
        )
      {
        print $ofh2 "DISCARDED because $pos too close to $other_mate_pos plus $frag_length +/- 10\n"
          if DEBUG;
        $discard_count++;
        last;
      }
      print $ofh2 "RETAINED: pos = $pos, other_mate_pos = $other_mate_pos\n"
        if DEBUG;

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

=head2 assemble_groups

 Title   : assemble_groups
 Usage   : $project->assemble_groups()
 Function: Identifies groups of half-mapping read pairs that align to the reference genome at
           locations near one another, and performs a local assembly on the unaligned mates in each
           group using Taipan
 Example : $project->bowtie_identify();
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
  my $num_aln             = shift;
  my $min_contig_length   = shift;
  my $halfmap_file        = shift || $self->{"halfmapping_file"};
  my $sorted_file         = $halfmap_file;
  $sorted_file            =~ s/(\.[^.]+)$/_sorted.sam/;
  my $work_dir            = $self->{"work_dir"};
  my $reads_basename      = $self->{"reads_basename"};
  my $frag_length         = $self->{"frag_length"};
  my $assembly_dir        = $self->{"assembly_dir"};
  my $results_file        = "$work_dir/${reads_basename}_contigs.fa";
  $self->{"contigs_file"} = $results_file;
  my $contigs_file_suffix = "_k${min_contig_length}_t1_o1.contigs";
  
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
  my $last_chr;
  my $last_len = 0;
  my $overlap_dist;
  my $count = 0;
  my $num_met_cutoff = 0;
  my $new_group = 1;
  my $left_pos = 0;
  while ( my $line = $ifh->getline ) {
    next if $line =~ m/^@/;
    my @fields = split( /\t/, $line );
    my ($id, $flags, $chr, $pos, $seq) = ($fields[0], int($fields[1]), $fields[2], $fields[7], $fields[9]);

    # Calculate overlap distance for this alignment as distance = 2 x insert_length + intron_length
    $overlap_dist = 2 * ($frag_length - 2 * length($seq)) + $intron_length;

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

    # Assemble groups containing at least $num_aln reads
    if ( scalar(@groups) >= $num_aln ) {
      ++$count;
      my $samfile = "$work_dir/assembly/${reads_basename}_group$count.sam";
      my $rawfile = "$work_dir/assembly/${reads_basename}_group$count.raw";
      my $ofh = IO::File->new( $samfile, 'w' ) or die "$0: Can't open $samfile: $!";
      my $ofh2 = IO::File->new( $rawfile, 'w' ) or die "$0: Can't open $rawfile: $!";
      while ( scalar(@groups) > 0 ) {
        my $ln = shift(@groups);
        my @fields = split( /\t/, $ln );
        print $ofh $ln; # Record the SAM alignment data to a SAM file for this group
        print $ofh2 "$fields[9]\n"; # Record the raw read to a .raw file for this group
      }
      $ofh->close();
      $ofh2->close();

      capture( "taipan -f $rawfile -c $min_contig_length -k 16 -o 1 -t 1" );
      die "$0: taipan exited unsuccessful" if ( $EXITVAL != 0 );

      # Append this group's results to the multi-FastA output file containing all contigs
      my $contigs_file = "${rawfile}${contigs_file_suffix}";
      if (!(-z $contigs_file)) {
        $num_met_cutoff++;
        my $ifh2 = IO::File->new( $contigs_file, 'r' ) or die "$0: Can't open $contigs_file: $!";
        while ( my $contig_line = $ifh2->getline ) {
          if ($contig_line =~ m/^>/) {
            print $ofh3 ">Group${count}(${left_pos}-${last_pos})|" . substr($contig_line, 1);
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
  if ( scalar(@groups) >= $num_aln ) {
    ++$count;
    my $samfile = "$work_dir/assembly/${reads_basename}_group$count.sam";
    my $rawfile = "$work_dir/assembly/${reads_basename}_group$count.raw";
    my $ofh = IO::File->new( $samfile, 'w' ) or die "$0: Can't open $samfile: $!";
    my $ofh2 = IO::File->new( $rawfile, 'w' ) or die "$0: Can't open $rawfile: $!";
    while ( scalar(@groups) > 0 ) {
      my $ln = shift(@groups);
      my @fields = split( /\t/, $ln );
      print $ofh $ln; # Record the SAM alignment data to a SAM file for this group
      print $ofh2 "$fields[9]\n"; # Record the raw read to a .raw file for this group
    }
    $ofh->close();
    $ofh2->close();

    capture( "taipan -f $rawfile -c $min_contig_length -k 16 -o 1 -t 1" );
    die "$0: taipan exited unsuccessful" if ( $EXITVAL != 0 );

    # Append this group's results to the multi-FastA output file containing all contigs
    my $contigs_file = "${rawfile}${contigs_file_suffix}";
    if (!(-z $contigs_file)) {
      $num_met_cutoff++;
      my $ifh2 = IO::File->new( $contigs_file, 'r' ) or die "$0: Can't open $contigs_file: $!";
      while ( my $contig_line = $ifh2->getline ) {
        if ($contig_line =~ m/^>/) {
          print $ofh3 ">Group${count}(${left_pos}-${last_pos})|" . substr($contig_line, 1);
        } else {
          print $ofh3 $contig_line;
        }
      }
      $ifh2->close;
    }
  }
  print "Identified $count groups for assembly.\n";
  print "Successfully assembled $num_met_cutoff contigs.\n";
}

=head2 build_blast_index

 Title   : build_blast_index
 Usage   : $project->build_blast_index()
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut

sub build_blast_index() {
  my $self                = shift;
  my $ref_genome          = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $index_dir           = $self->{"index_dir"};
  my $reads_basename      = $self->{"reads_basename"};
  unless ( -e "$index_dir/$ref_genome_basename.nhr"
    && "$index_dir/$ref_genome_basename.nin"
    && "$index_dir/$ref_genome_basename.nsq" )
  {
    # Call formatdb to create the Blast index
    print "Creating Blast index in $index_dir...\n";
    capture( "formatdb -p F -i $ref_genome -n $index_dir/$ref_genome_basename" );
    die "$0: formatdb exited unsuccessful" if ( $EXITVAL != 0 );
  }
  else {
    print "Blast index already exists, not re-creating.\n";
    print "Using Blast index at $index_dir/$ref_genome_basename.n*\n";
  }
}

=head2 align_groups_blast

 Title   : align_groups_blast
 Usage   : $project->align_groups_blast()
 Function: Aligns contigs assembled from half-mapping reads
 Example : $project->bowtie_identify();
           $project->filter( 8 );
           $project->assemble_groups( 250, 3 );
           $project->align_groups_blast();
 Returns : No return value
 Args    : Path to multi-FastA format file containing contigs for alignment (optional)

=cut

sub align_groups_blast() {
  my $self                          = shift;
  my $contigs_file                  = shift || $self->{"contigs_file"};
  my $index_dir                     = $self->{"index_dir"};
  my $work_dir                      = $self->{"work_dir"};
  my $reads_basename                = $self->{"reads_basename"};
  my $ref_genome_basename           = $self->{"ref_genome"}->{"basename"};
  my $output_file                   = "$work_dir/${reads_basename}_contigs_aligned_blast";
  my $index                         = "$index_dir/${ref_genome_basename}";
  $self->{"contigs_alignment_file"} = $output_file;

  return if (-z $contigs_file || !(-e $contigs_file));
  print "Aligning assembled contigs to reference genome using Blast...\n";
  
  # Call Blast to run the mapping
  capture( "blastall -p blastn -d $index -i $contigs_file -o $output_file" );
  die "$0: Blast exited unsuccessful" if ( $EXITVAL != 0 );
  
  print "Alignment results saved to $output_file\n";
}

=head2 align_groups_clustal

 Title   : align_groups_clustal
 Usage   : $project->align_groups_clustal()
 Function: Aligns contigs assembled from half-mapping reads
 Example : $project->bowtie_identify();
           $project->filter( 8 );
           $project->assemble_groups( 250, 3 );
           $project->align_groups_clustal();
 Returns : No return value
 Args    : Path to multi-FastA format file containing contigs for alignment (optional)

=cut

sub align_groups_clustal() {
  my $self                          = shift;
  my $contigs_file                  = shift || $self->{"contigs_file"};
  my $reads_basename                = $self->{"reads_basename"};
  my $ref_genome                    = $self->{"ref_genome"}->{"full_pathname"};
  my $work_dir                      = $self->{"work_dir"};
  my $output_file                   = "$work_dir/${reads_basename}_contigs_aligned_clustal";
  $self->{"contigs_alignment_file"} = $output_file;
  my $clustal_input_file            = "$work_dir/${reads_basename}_pre-alignment";

  return if (-z $contigs_file || !(-e $contigs_file));
  print "Aligning assembled contigs to reference genome using Clustal...\n";
  
  # Append reference genome file to assembled contigs for multiple sequence alignment
  capture( "cat $contigs_file $ref_genome > $clustal_input_file");
  die "$0: cat exited unsuccessful" if ( $EXITVAL != 0 );
  
  # Call Clustal to run the mapping
  capture( "clustalw -infile=$clustal_input_file -gapopen=50 -gapext=0.01 -outfile=$output_file" );
  die "$0: clustalw exited unsuccessful" if ( $EXITVAL != 0 );
  
  print "Alignment results saved to $output_file\n";
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

# Returns 1 if the given sum-of-flags values represent a read pair containing mates that align independently but
# do not meet fragment length constraints, otherwise returns 0
sub _isFarMappingPair {
  my $flags1 = int(shift);
  my $flags2 = int(shift);
  return ( ((($flags1 & 4) == 0) && (($flags1 & 8) != 0))
        && ((($flags2 & 4) == 0) && (($flags2 & 8) != 0)) );
}

# Returns 1 if the given sum-of-flags value identifies a non-aligning mate, otherwise returns 0
sub _isUnalignedMate {
  my $flags = int(shift);
  return ( (($flags & 4) != 0) )
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

1;
