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
  $self->{"scripts_dir"}   = undef;   # Scalar path to scripts directory
  $self->{"work_dir"}      = undef;   # Scalar path to work directory
  $self->{"data_basename"} = undef;   # Scalar base filename for read pairs/data
  $self->{"ref_genome"}    = undef;   # Hash of reference genome related values
  $self->{"bowtie_db"}     = undef;   # Hash of data directories for bowtie
  bless($self);
  return $self;
}

=head2 set_work_dir

 Title   : set_work_dir
 Usage:  : $project->set_work_dir( "path_to_scripts" , "" )
 Function: Creates directory for pipeline working data and output files
 Example : my $work_dir = "/home/theis/work";
           my $resume_dir = "";
           $project->set_work_dir( $scripts_dir, $resume_dir )
 Returns : The path to the directory for pipeline working data and output files
 Args    : Scalar of full path to the parent directory of the work directory, and a scalar of 
           full path to work folder from a previous run (or containing existing unmapped read pairs)

=cut

sub set_work_dir {
  my $self        = shift;
  my $scripts_dir = shift;
  $self->{"scripts_dir"} = $scripts_dir;
  my $work_dir = shift || "$scripts_dir/run-" . &_datestamp;
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

  my $velvet_dir = "$work_dir/velvet-data";
  unless ( -d $velvet_dir ) {
    &_make_dir($velvet_dir);
    print "Created velvet directory $velvet_dir\n";
  }
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
  my $self                = shift;
  my $ref_genome_filename = shift;
  my $data_basename       = shift;

  # Use the base of the read pairs filenames for saving data throughout pipeline
  $self->{"data_basename"} = $data_basename;

  # Split reference genome pathname into components
  my ( $file, $dir, $ext ) = fileparse( $ref_genome_filename, qr/\.[^.]*/ );
  $self->{"ref_genome"}->{"full_pathname"} = $ref_genome_filename;
  $self->{"ref_genome"}->{"basename"}      = $file;
  $self->{"ref_genome"}->{"dir"}           = $dir;
}

=head2 mapping_setup

 Title   : mapping_setup
 Usage   : $project->mapping_setup( "bowtie2dir", "bowtie_index_dir", "reads_dir", "read_pairs_base_name")
 Function: Sets references to data files used for mapping reads to the reference genome
 Example : my $bowtie_dir = "/home/theis/bowtie2-2.0.0-beta6/";
           my $index_dir = "/home/theis/bt2/";
           my $reads_dir = "/home/theis/reads/";
           $project->mapping_setup( $bowtie_dir, $index_dir, $reads_dir, "reads" );
 Returns : No return value
 Args    : Scalar full path to the Bowtie 2 executable directory, scalar full path to the Bowtie 2
           index directory, scalar full path to the directory containing the Fastq read pairs files,
           and a scalar of the base name of the read pairs files (filename without "_1.fq" or
           "_2.fq" ending)

=cut

sub mapping_setup {
  my $self                = shift;
  my $bowtie1_dir         = shift;
  my $bowtie2_dir         = shift;
  my $bowtie_index_dir    = shift;
  my $reads_dir           = shift;
  my $data_basename       = shift;
  my $scripts_dir         = $self->{"scripts_dir"};
  my $work_dir            = $self->{"work_dir"};
  my $ref_genome          = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_dir      = $self->{"ref_genome"}->{"dir"};
  unless ( -d $ref_genome_dir ) {
    die "$0: reference genome directory $ref_genome_dir does not exist";
  }

  # Perform basic validation on reference genome
  my $results = capture("$scripts_dir/validate_fasta.pl -i $ref_genome");
  if ( $EXITVAL != 0 ) {
    die "$0: validate_fasta.pl exited unsuccessful";
  }

  # Ensure that read pairs files exist
  $reads_dir = $1 if ( $reads_dir =~ /(.*)\/$/ );
  my $reads_file_one = "$reads_dir/${data_basename}_1.fq";
  my $reads_file_two = "$reads_dir/${data_basename}_2.fq";
  unless ( -e $reads_file_one ) {
    die "$0: reads file $reads_file_one does not exist";
  }
  unless ( -e $reads_file_two ) {
    die "$0: reads file $reads_file_two does not exist";
  }

  # Perform basic validation on read pairs files
  $results =
    capture("$scripts_dir/validate_fastq.pl -i $reads_dir/$data_basename");
  if ( $EXITVAL != 0 ) {
    die "$0: validate_fastq.pl exited unsuccessful";
  }

  # Create directory for bowtie index files if necessary
  $bowtie_index_dir = $1 if ( $bowtie_index_dir =~ /(.*)\/$/ );
  unless ( -d $bowtie_index_dir ) {
    &_make_dir($bowtie_index_dir);
  }

  # Remember paths needed for running bowtie
  $self->{"bowtie_db"} = {
    bowtie1_dir      => $bowtie1_dir,
    bowtie2_dir      => $bowtie2_dir,
    bowtie_index_dir => $bowtie_index_dir,
    reads_file_one   => $reads_file_one,
    reads_file_two   => $reads_file_two
  };
  print "Using reference genome $ref_genome\n";
  print "Using reads file 1: $reads_file_one\n";
  print "Using reads file 2: $reads_file_two\n";
}

=head2 build_bowtie2_index

 Title   : build_bowtie2_index
 Usage   : $project->build_bowtie2_index();
 Function: Runs the bowtie2-build indexer to create a Bowtie index from the reference genome
 Example : $project->mapping_setup( $bowtie_dir, $index_dir, $reads_dir, "reads" );
           $project->build_bowtie2_index();
 Returns : No return value
 Args    : No arguments

=cut

sub build_bowtie2_index {
  my $self                = shift;
  my $ref_genome          = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie2_dir         = $self->{"bowtie_db"}->{"bowtie2_dir"};
  my $bowtie_index_dir    = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  unless ( -e "$bowtie_index_dir/$ref_genome_basename.1.bt2"
    && "$bowtie_index_dir/$ref_genome_basename.2.bt2"
    && "$bowtie_index_dir/$ref_genome_basename.3.bt2"
    && "$bowtie_index_dir/$ref_genome_basename.4.bt2"
    && "$bowtie_index_dir/$ref_genome_basename.rev.1.bt2"
    && "$bowtie_index_dir/$ref_genome_basename.rev.2.bt2" )
  {

    # Call bowtie-build to create the bowtie index
    print "Creating bowtie index in $bowtie_index_dir...\n";
    my $results = capture( "$bowtie2_dir/bowtie2-build "
        . "$ref_genome $bowtie_index_dir/$ref_genome_basename" );
    if ( $EXITVAL != 0 ) {
      die "$0: bowtie2-build exited unsuccessful";
    }
  }
  else {
    print "Bowtie2 index already exists, not re-creating.\n";
    print "Using Bowtie2 index at $bowtie_index_dir/$ref_genome_basename.*\n";
  }
}

=head2 run_bowtie2_mapping

 Title   : run_bowtie2_mapping
 Usage   : $project->run_bowtie2_mapping( num_threads, minins, maxins )
 Function: Aligns reads to the reference genome and saves output file
 Example : $project->mapping_setup( $bowtie_dir, $index_dir, $reads_dir, "reads" );
           $project->build_bowtie2_index();
           $project->run_bowtie2_mapping( 8, 100, 500 );
 Returns : No return value
 Args    : Number of parallel search threads to use for Bowtie 2, numeric length to use for
           Bowtie 2 -I/--minins parameter, numeric length to use for Bowtie 2 -X/--maxins parameter

=cut

sub run_bowtie2_mapping {
  my $self                = shift;
  my $num_threads         = shift;
  my $minins              = shift;
  my $maxins              = shift;
  my $data_basename       = $self->{"data_basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie2_dir         = $self->{"bowtie_db"}->{"bowtie2_dir"};
  my $bowtie_index_dir    = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  my $reads_file_one      = $self->{"bowtie_db"}->{"reads_file_one"};
  my $reads_file_two      = $self->{"bowtie_db"}->{"reads_file_two"};
  my $work_dir            = $self->{"work_dir"};
  print "Running mapping using bowtie2, using $num_threads threads...\n";

  # Call bowtie to run the mapping
  my $results = capture(
        "$bowtie2_dir/bowtie2 -x $bowtie_index_dir/$ref_genome_basename "
      . "--threads $num_threads --reorder --sam-no-hd "
      . "--maxins $maxins --minins $minins "
      . "--no-discordant "
      . "--no-contain --no-overlap "
      . "-k 3 -1 $reads_file_one -2 $reads_file_two "
      .

      #"--al-conc $work_dir/${data_basename}_al-conc.%.fq " .
      #"--un-conc $work_dir/${data_basename}_un-conc.%.fq " .
      "-S $work_dir/${data_basename}_alignment.sam"
  );

  if ( $EXITVAL != 0 ) {
    die "$0: bowtie2 exited unsuccessful";
  }

}

=head2 bowtie2_identify

 Title   : bowtie2_identify
 Usage   : $project->bowtie2_identify( $num_threads )
 Function: Identifies half-mapping read pairs from the read alignment output file
 Example : $project->build_bowtie2_index();
           $project->mapping_setup( $bowtie_dir, $index_dir, $reads_dir, "reads" );
           $project->bowtie2_identify();
 Returns : No return value
 Args    : No argumnets

=cut

sub bowtie2_identify {
  my $self          = shift;
  my $data_basename = $self->{"data_basename"};
  my $work_dir      = $self->{"work_dir"};

  print "Identifying half-mapping read pairs...\n";

# # Remove all rows we're not interested in now
# my $results = capture( "cut -f 1,2,4,10,11 $work_dir/${data_basename}_alignment.sam " .
#        "> $work_dir/${data_basename}_alignment_pruned" );
# if( $EXITVAL != 0 ) {
#   die "$0: cut exited unsuccessful";
# }

  my $ifh = new IO::File( "$work_dir/${data_basename}_alignment.sam", 'r' )
    or die "Can't open $work_dir/${data_basename}_alignment.sam: $!";
  my $ofh = IO::File->new( "$work_dir/${data_basename}_halfmapping1.sam", 'w' )
    or die "Can't create $work_dir/${data_basename}_halfmapping1.sam: $!";
  my $ofh2 = IO::File->new( "$work_dir/${data_basename}_halfmapping2.sam", 'w' )
    or die "Can't create $work_dir/${data_basename}_halfmapping2.sam: $!";
  my $ofh3 = IO::File->new( "$work_dir/${data_basename}_multmapping.sam", 'w' )
    or die "Can't create $work_dir/${data_basename}_multmapping.sam: $!";
  my $prev_id         = "";
  my $discard_this_id = 0;
  my $mapping         = 0;
  my @print_lines;

  while ( my $line = $ifh->getline ) {

    # Ensure we're reading an alignment line and not a header line
    if ( $line =~ m/^@/ ) {
      next;
    }
    my @fields   = split( /\t/, $line );
    my $mate_id  = $fields[0];
    my $flag_sum = $fields[1];

    # Check if this mate ID is different from the last line
    if ( $mate_id ne $prev_id ) {

      # This line is the first line with a new read ID
      if ( scalar @print_lines == 2 && $mapping == 1 && $prev_id ne "" ) {

        # Check if the pair is discordant
        my $line1 = shift(@print_lines);
        my @fd1   = split( /\t/, $line1 );
        my $flag1 = $fd1[1];
        my $line2 = shift(@print_lines);
        my @fd2   = split( /\t/, $line2 );
        my $flag2 = $fd2[1];
        if ( &_isDiscordant( $flag1, $flag2 ) ) {
          print $ofh2 "$line1$line2";
        }
        else {

          # Print the half-mapping mate pairs with the last line's ID
          print $ofh "$line1$line2";
        }
      }
      elsif ( scalar @print_lines == 3 && $prev_id ne "" ) {

        # Print the trio where two mates align concordantly, and one has a secondary alignment
        foreach my $ln (@print_lines) {
          print $ofh3 $ln;
        }
        @print_lines = ();
      } # Note: can print cases of multiple (>3/pair) alignments by extending this if-else block
      $discard_this_id = 0;
      $prev_id         = $mate_id;
    }
    else {

      # This line is the second or greater line with this read ID
      if ($discard_this_id) {
        next;
      }
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
    if ( &_isMapping($flag_sum) ) {

      # Save this full-mapping line for printing
      $mapping = 2;
      push( @print_lines, $line );
    }
    elsif ( &_isNonMapping($flag_sum) ) {

      # We have a non-mapping ID. Discard all lines with this ID.
      $mapping         = 0;
      $discard_this_id = 1;
    }
    elsif ( &_isHalfMapping($flag_sum) ) {

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
  $ofh2->close;
  $ofh3->close;
}

=head2 filter

 Title   : filter
 Usage   : $project->filter( $num_threads )
 Function: Reduces the number of half-mapping read pairs by re-running alignments using looser
           matching criteria.
 Example : $project->bowtie2_identify();
           $project->filter( 8 );
 Returns : No return value
 Args    : Number of parallel search threads to use for Bowtie 2

=cut

sub filter {
  my $self                = shift;
  my $num_threads         = shift;
  my $work_dir            = $self->{"work_dir"};
  my $data_basename       = $self->{"data_basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie2_dir         = $self->{"bowtie_db"}->{"bowtie2_dir"};
  my $bowtie_index_dir    = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  my $ifh = new IO::File( "$work_dir/${data_basename}_halfmapping1.sam", 'r' )
    or die "$0: Can't open $work_dir/${data_basename}_halfmapping1.sam: $!";
  my $ofh =
    new IO::File( "$work_dir/${data_basename}_fake_paired_end_1.fq", 'w' )
    or die "$0: Can't open $work_dir/${data_basename}_fake_paired_end_1.fq: $!";
  my $ofh2 =
    new IO::File( "$work_dir/${data_basename}_fake_paired_end_2.fq", 'w' )
    or die "$0 Can't open $work_dir/${data_basename}_fake_paired_end_2.fq: $!";

  # Define length of each mate for the fake read pairs
  my $mate_length = 16;

  # Define length of read pair for fake paired-end alignment
  my $read_length = 500;    # TODO get this from infer_fraglen.pl

  # Boolean to track whether we have any half-mapping pairs
  my $have_candidates = 0;

  # Get the mates we want and write them as fake paired-end reads to FastQ
  while ( my $line = $ifh->getline ) {
    if ( $line =~ m/^@/ ) {
      next;
    }
    my @fields         = split( /\t/, $line );
    my $id             = $fields[0];
    my $flag_sum       = $fields[1];
    my $sequence       = $fields[9];
    my $quality_scores = $fields[10];

    # Create a fake pair from the unaligned mate in the pair
    if ( &_isUnalignedMate($flag_sum) ) {

      # Set a flag indicating that we got at least one half-mapping pair
      $have_candidates = 1;

      # Reverse complement reads that aligned to the reverse strand, to get original read
      if ( ( $flag_sum & 16 ) == 1 ) {
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
      print $ofh "\@$id\/1\n$mate_1\n+\n$quality_scores_1\n";
      print $ofh2 "\@$id\/2\n$mate_2\n+\n$quality_scores_2\n";
    }
  }
  $ifh->close;
  $ofh->close;
  $ofh2->close;

  # Define the maximum/minimum insert length as read length +/- 10
  if ( !$have_candidates ) {
    print STDERR "$0: filter(): No half-mapping pairs available.\n";
    return 1;
  }

  # Run Bowtie 2 with the fake read pairs as input
  my $results = capture(
        "$bowtie2_dir/bowtie2 -x $bowtie_index_dir/$ref_genome_basename "
      . "--threads $num_threads --reorder --sam-no-hd "
      . "--no-mixed --no-discordant --no-contain --no-overlap "
      . "-1 $work_dir/${data_basename}_fake_paired_end_1.fq "
      . "-2 $work_dir/${data_basename}_fake_paired_end_2.fq "
      .

      #"--al-conc $work_dir/${data_basename}_fake_pairs_al-conc.%.fq " .
      "-S $work_dir/${data_basename}_fake_pairs_alignment.sam"
  );
  if ( $EXITVAL != 0 ) {
    die "$0: bowtie2 exited unsuccessful";
  }

  # Find reads representing insertions/deletions by considering all the fake pairs that aligned
  $ifh =
    new IO::File( "$work_dir/${data_basename}_fake_pairs_alignment.sam", 'r' )
    or die
    "$0: Can't open $work_dir/${data_basename}_fake_pairs_alignment.sam: $!";
  my $ifh2 = new IO::File( "$work_dir/${data_basename}_halfmapping1.sam", 'r' )
    or die "$0: Can't open $work_dir/${data_basename}_halfmapping1.sam: $!";
  $ofh = new IO::File( "$work_dir/${data_basename}_filtered.sam", 'w' )
    or die "$0: Can't open $work_dir/${data_basename}_filtered.sam: $!";
  $ofh2 =
    new IO::File( "$work_dir/${data_basename}_fake_pairs_alignment.debug", 'w' )
    or die
    "$0: Can't open $work_dir/${data_basename}_fake_pairs_alignment.debug: $!"
    if DEBUG;

  # Get an alignment from the originally half-mapping pairs
  while ( my $line = $ifh2->getline ) {
    if ( $line =~ m/^@/ ) {
      next;
    }
    my @fields         = split( /\t/, $line );
    my $orig_id        = $fields[0];
    my $orig_flag_sum  = $fields[1];
    my $other_mate_pos = $fields[7];

    # Consider only the unaligned mate from the half-mapping read pair
    if ( !&_isUnalignedMate($orig_flag_sum) ) {
      next;
    }

    # Get the corresponding alignment from the fake pairs
    my $fake_line;
    while ( $fake_line = $ifh->getline ) {
      if ( $fake_line =~ m/^@/ ) {
        next;
      }
      my @fields = split( /\t/, $fake_line );
      my $id     = $fields[0];
      my $pos    = $fields[3];
      while ( $id ne $orig_id ) {
        $fake_line = $ifh->getline;
        my @fields = split( /\t/, $fake_line );
        $id  = $fields[0];
        $pos = $fields[3];
      }

      print $ofh2 "\nFAKE PAIR ALIGNMENT LINE 1: $fake_line"         if DEBUG;
      print $ofh2 "ORIGINAL HALF-MAPPING PAIR UNALIGNED MATE: $line" if DEBUG;
      print $ofh2 "ORIG_POS: $other_mate_pos POS: $pos\n"            if DEBUG;

      # Discard fake pairs that align close to the originally aligning mate in the half-mapping read pair
      if (
        $pos != 0
        && ( $pos < ( $other_mate_pos + $read_length + 10 )
          && $pos > ( $other_mate_pos + $read_length - 10 ) )
        )
      {
        print $ofh2 "DISCARDED because $pos too close to $other_mate_pos plus $read_length plus or minus ten\n"
          if DEBUG;
        last;
      }
      print $ofh2 "RETAINED: pos = $pos, other_mate_pos = $other_mate_pos\n"
        if DEBUG;

      # Keep this mate
      print $ofh $line;
      last;
    }
  }
  $ifh->close;
  $ifh2->close;
  $ofh->close;
  $ofh2->close if DEBUG;
}

=head2 group

 Title   : group
 Usage   : $project->group()
 Function: Gathers groups of half-mapping read pairs that align to the reference genome at locations
           near one another, and performs a local assembly on the unaligned mates in each group
 Example : $project->bowtie2_identify();
           $project->filter( 8 );
           $project->group();
 Returns : No return value
 Args    : No arguments

=cut

sub group {
  my $self          = shift;
  my $work_dir      = $self->{"work_dir"};
  my $data_basename = $self->{"data_basename"};

  # Expected intron length
  my $intron_length = 250;

  # Insert length
  my $ins = 100;

  # Sort by alignment positionn
  my $results = capture(
    "cat $work_dir/${data_basename}_filtered.sam " .

      #"$work_dir/${data_basename}_halfmapping2.sam " .
      "> $work_dir/${data_basename}_halfmapping_all.sam"
  );
  if ( $EXITVAL != 0 ) {
    die "$0: cat exited unsuccessful";
  }
  $results = capture(
    "sort -k 8,8 -n -o $work_dir/${data_basename}_halfmapping_all_sorted.sam "
      . "$work_dir/${data_basename}_halfmapping_all.sam" );
  if ( $EXITVAL != 0 ) {
    die "$0: sort exited unsuccessful";
  }

  # Get the clusters one at a time
  my $ifh =
    new IO::File( "$work_dir/${data_basename}_halfmapping_all_sorted.sam", 'r' )
    or die
    "$0: Can't open $work_dir/${data_basename}_halfmapping_all_sorted.sam: $!";
  my @cluster;
  my $have_cluster = 0;
  my $last_pos     = 0;
  my $overlap_dist = 2 * $ins + $intron_length;

  my $id;
  my $flag;
  my $pos;
  my $seq;
  my $count = 0;
  while ( my $line = $ifh->getline ) {
    if ( $line =~ m/^@/ ) {
      next;
    }
    my @fields = split( /\t/, $line );
    $id   = $fields[0];
    $flag = $fields[1];
    $pos  = $fields[7];
    $seq  = $fields[9];

    # Add or discard based on position, treating upstream mates differently from downstream mates
    if ( ( $flag & 20 ) == 0 ) {
      if ( ( $pos - $last_pos ) < $overlap_dist ) {
        push( @cluster, $line );
        print "before: last_pos=$last_pos, pos=$pos\n";
        $last_pos = $pos;
        print "after: last_pos=$last_pos, pos=$pos\n";
        next;
      }
      $last_pos = $pos;
    }
    else {
      if ( ( $pos - $last_pos ) < ( $ins - length($seq) + $intron_length ) ) {
        push( @cluster, $line );
        $last_pos = $pos;
        next;
      }
      $last_pos = $pos;
    }

    # Run velveth on groups of 2 or more reads using stdout
    if ( scalar(@cluster) >= 2 ) {

      # Open output file
      ++$count;
      my $ofh = new IO::File(
        "$work_dir/velvet-data/${data_basename}_cluster_$count.txt", 'w' )
        or die
        "$0: Can't open $work_dir/velvet-data/${data_basename}_cluster_$count.txt: $!";

      while ( scalar(@cluster) > 0 ) {
        print $ofh shift(@cluster);
      }
      $ofh->close();
      $results =
        capture( "velveth $work_dir/velvet-data/cluster_$count 13 -sam -short "
          . "$work_dir/velvet-data/${data_basename}_cluster_$count.txt" );
      if ( $EXITVAL != 0 ) {
        die "$0: velveth exited unsuccessful";
      }

# Delete the SAM file of read data unless we're debugging
#unlink( "$work_dir/velvet-data/${data_basename}_cluster_$count.txt" ) if !DEBUG
#        or die "$0: cannot delete $work_dir/velvet-data/${data_basename}_cluster_$count.txt";

      # Run velvetg on groups
      $results = capture("velvetg $work_dir/velvet-data/cluster_${count}/");
      if ( $EXITVAL != 0 ) {
        die "$0: velvetg exited unsuccessful";
      }
    }
    else {
      @cluster = ();
    }

  }

}

############ Subroutines for internal use by this module ############

# Returns the reverse complement of the given string
sub _reverseComplement {
  my ($value) = shift;
  $value = scalar reverse $value;
  for ( 0 .. length($value) - 1 ) {
    substr( $value, $_, 1 ) = &_complement( substr( $value, $_, 1 ) );
  }
  return $value;
}

# Returns the complement of the given string
sub _complement {
  my %complementMap = (
    "A" => "T",
    "T" => "A",
    "a" => "t",
    "t" => "a",
    "C" => "G",
    "G" => "C",
    "c" => "g",
    "g" => "c",
    "R" => "Y",
    "Y" => "R",
    "r" => "y",
    "y" => "r",
    "M" => "K",
    "K" => "M",
    "m" => "k",
    "k" => "m",
    "S" => "S",
    "W" => "W",
    "s" => "s",
    "w" => "w",
    "B" => "V",
    "V" => "B",
    "b" => "v",
    "v" => "b",
    "H" => "D",
    "D" => "H",
    "h" => "d",
    "d" => "h",
    "N" => "N",
    "." => ".",
    "n" => "n"
  );
  my $complementedString = $complementMap{ $_[0] }
    or die "$0: Can't get reverse complement for '$_[0]'";
  return $complementedString;
}

# Returns 1 if the given sum-of-flags value identifies a mate in a half-mapping read pair, otherwise returns 0
sub _isHalfMapping {
  my $flag = shift;
  if ( $flag == 69
    || $flag == 73
    || $flag == 89
    || $flag == 101
    || $flag == 133
    || $flag == 137
    || $flag == 153
    || $flag == 165 )
  {
    return 1;
  }
  else {
    return 0;
  }
}

# Returns 1 if the given sum-of-flags value identifies a mate in a concordant alignment, otherwise returns 0
sub _isMapping {
  my $flag = shift;
  if ( $flag == 83
    || $flag == 99
    || $flag == 147
    || $flag == 163 )
  {
    return 1;
  }
  else {
    return 0;
  }
}

# Returns 1 if the given sum-of-flags value identifies a mate in a no-alignments pair, otherwise returns 0
sub _isNonMapping {
  my $flag = shift;
  if ( $flag == 77 || $flag == 141 ) {
    return 1;
  }
  else {
    return 0;
  }
}

# Returns 1 if the given pair of sum-of-flags values represents mates that align discordantly, otherwise returns 0
sub _isDiscordant {
  my $flag1 = shift;
  my $flag2 = shift;

  if ( ( $flag1 == 73 && $flag2 == 153 )
    || ( $flag1 == 153 && $flag2 == 73 )
    || ( $flag1 == 89  && $flag2 == 137 )
    || ( $flag1 == 137 && $flag2 == 89 ) )
  {
    return 1;
  }
  else {
    return 0;
  }
}

# Returns 1 if the given sum-of-flags value identifies an unaligned mate in a half-mapping pair, otherwise returns 0
sub _isUnalignedMate {
  my $flag = shift;
  if ( $flag == 69
    || $flag == 101
    || $flag == 133
    || $flag == 165 )
  {
    return 1;
  }
  else {
    return 0;
  }
}

# Checks if the given folder exists, and removes trailing slash from string
sub _check_dir {
  my $dir = shift;
  if ( -d $dir ) {
    $dir = $1 if ( $dir =~ /(.*)\/$/ );
    return $dir;
  }
  else {
    die("$0: directory $dir does not exist");
  }
}

# Creates the given directory, and gives an error message on failure
sub _make_dir {
  my ($dirname) = @_;
  mkdir $dirname, 0755 || die "$0: could not create $dirname";
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
