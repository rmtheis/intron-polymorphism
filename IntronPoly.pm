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
  $self->{"scripts_dir"} = undef; # Scalar path to scripts directory
  $self->{"work_dir"} = undef; # Scalar path to work directory
  $self->{"data_basename"} = undef; # Scalar base filename for read pairs/data
  $self->{"ref_genome"} = undef; # Hash of reference genome related values
  $self->{"bowtie_db"} = undef; # Hash of data directories for bowtie
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
  my $self = shift;
  my $scripts_dir = shift;
  $self->{"scripts_dir"} = $scripts_dir;
  my $work_dir = shift || "$scripts_dir/run-" . &_datestamp;
  unless (-e $work_dir) {  
    $work_dir = $1 if ($work_dir =~ /(.*)\/$/);
    &_make_dir( $work_dir );
    print "Created work directory $work_dir\n"
  } else {
    $work_dir = &_check_dir( $work_dir );
    print "Using existing work directory $work_dir\n";
  }
  $self->{"work_dir"} = $work_dir;
  return $work_dir;
}

=head2 build_db

 Title   : build_db
 Usage   :
 Function:
 Example : 
 Returns : 
 Args    :

=cut

sub build_db {
  my $self = shift;
  my $ref_genome_filename = shift;
  my $data_basename = shift;

  # Use the base of the read pairs filenames for saving data throughout pipeline
  $self->{"data_basename"} = $data_basename;

  # Split reference genome pathname into components
  my ($file, $dir, $ext) = fileparse($ref_genome_filename, qr/\.[^.]*/);
  $self->{"ref_genome"}->{"full_pathname"} = $ref_genome_filename;
  $self->{"ref_genome"}->{"basename"} = $file;
  $self->{"ref_genome"}->{"dir"} = $dir;
}

=head2 mapping_setup

 Title   : mapping_setup
 Usage   :
 Function:
 Example : 
 Returns : 
 Args    :

=cut

sub mapping_setup {
  my $self = shift;
  my $bowtie1_dir = shift;
  my $bowtie2_dir = shift;
  my $bowtie_index_dir = shift;
  my $reads_dir = shift;
  my $data_basename = shift;
  my $scripts_dir = $self->{"scripts_dir"};
  my $work_dir = $self->{"work_dir"};
  my $ref_genome = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_dir = $self->{"ref_genome"}->{"dir"};
  unless (-d $ref_genome_dir) {
    die "$0: reference genome directory $ref_genome_dir does not exist" ;
  }

  # Perform basic validation on reference genome
  my $results = capture( "$scripts_dir/validate_fasta.pl -i $ref_genome" );
  if ( $EXITVAL != 0 ) {
    die "$0: validate_fasta.pl exited unsuccessful";
  }

  # Ensure that read pairs files exist
  $reads_dir = $1 if ($reads_dir =~ /(.*)\/$/);
  my $reads_file_one = "$reads_dir/$data_basename.1.fq";
  my $reads_file_two = "$reads_dir/$data_basename.2.fq";
  unless (-e $reads_file_one) {
    die "$0: reads file $reads_file_one does not exist";
  }
  unless (-e $reads_file_two) {
    die "$0: reads file $reads_file_two does not exist";
  }

  # Perform basic validation on read pairs files
  $results = capture ( "$scripts_dir/validate_fastq.pl -i $reads_dir/$data_basename" );
  if ( $EXITVAL != 0 ) {
    die "$0: validate_fastq.pl exited unsuccessful";
  }

  # Create directory for bowtie index files if necessary
  $bowtie_index_dir = $1 if ($bowtie_index_dir =~ /(.*)\/$/);
  unless (-d $bowtie_index_dir) {
    &_make_dir( $bowtie_index_dir );
  }

  # Remember paths needed for running bowtie
  $self->{"bowtie_db"} = {
                          bowtie1_dir => $bowtie1_dir,
                          bowtie2_dir => $bowtie2_dir,
                          bowtie_index_dir => $bowtie_index_dir,
                          reads_file_one => $reads_file_one,
                          reads_file_two => $reads_file_two
                         };
  print "Using reference genome $ref_genome\n";
  print "Using reads file 1: $reads_file_one\n";
  print "Using reads file 2: $reads_file_two\n";
}

=head2 build_bowtie1_index

 Title   : build_bowtie1_index
 Usage   :
 Function:
 Example : 
 Returns : 
 Args    :

=cut

sub build_bowtie1_index {
  my $self = shift;
  my $ref_genome = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie1_dir = $self->{"bowtie_db"}->{"bowtie1_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  unless (-e "$bowtie_index_dir/$ref_genome_basename.1.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.2.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.3.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.4.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.rev.1.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.rev.2.ebwt") {
    # Call bowtie-build to create the bowtie index
    print "Creating bowtie index in $bowtie_index_dir...";
    my $results = capture( "$bowtie1_dir/bowtie-build " .
                           "$ref_genome $bowtie_index_dir/$ref_genome_basename" );
    if ( $EXITVAL != 0 ) {
      die "$0: bowtie-build exited unsuccessful";
    }
  } else {
    print "Bowtie1 index already exists, not re-creating.\n";
    print "Using Bowtie1 index at $bowtie_index_dir/$ref_genome_basename.*\n";
  }
}

=head2 build_bowtie2_index

 Title   : build_bowtie2_index
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub build_bowtie2_index {
  my $self = shift;
  my $ref_genome = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie2_dir = $self->{"bowtie_db"}->{"bowtie2_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  unless (-e "$bowtie_index_dir/$ref_genome_basename.1.bt2" &&
             "$bowtie_index_dir/$ref_genome_basename.2.bt2" &&
             "$bowtie_index_dir/$ref_genome_basename.3.bt2" &&
             "$bowtie_index_dir/$ref_genome_basename.4.bt2" &&
             "$bowtie_index_dir/$ref_genome_basename.rev.1.bt2" &&
             "$bowtie_index_dir/$ref_genome_basename.rev.2.bt2") {
    # Call bowtie-build to create the bowtie index
    print "Creating bowtie index in $bowtie_index_dir...";
    my $results = capture( "$bowtie2_dir/bowtie2-build " .
                           "$ref_genome $bowtie_index_dir/$ref_genome_basename" );
    if ( $EXITVAL != 0 ) {
      die "$0: bowtie2-build exited unsuccessful";
    }
    print "finished.\n";
  } else {
    print "Bowtie2 index already exists, not re-creating.\n";
    print "Using Bowtie2 index at $bowtie_index_dir/$ref_genome_basename.*\n";
  }
}

=head2 run_bowtie1_mapping

 Title   : run_bowtie1_mapping
 Usage   : 
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut
 
sub run_bowtie1_mapping {
  my $self = shift;
  my $num_threads = shift;
  my $data_basename = $self->{"data_basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie1_dir = $self->{"bowtie_db"}->{"bowtie1_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  my $reads_file_one = $self->{"bowtie_db"}->{"reads_file_one"};
  my $reads_file_two = $self->{"bowtie_db"}->{"reads_file_two"};
  my $work_dir = $self->{"work_dir"};
  print "Running mapping using bowtie1, using $num_threads threads...\n";

  # Call bowtie to run the mapping
  my $results = capture( "$bowtie1_dir/bowtie $bowtie_index_dir/$ref_genome_basename " .
                         "--threads $num_threads --suppress 3,6 --time " .
                         "--un $work_dir/${data_basename}_unaligned " . 
                         "-m 1 -1 $reads_file_one -2 $reads_file_two " .
                         "--al $work_dir/${data_basename}_aligned " .
                         "$work_dir/${data_basename}_hits.map"
                       );

  if( $EXITVAL != 0 ) {
    die "$0: bowtie1 exited unsuccessful";
  }
  print "Bowtie1 finished.\n";

  # Sort aligning reads by match position
  #print "Sorting aligned reads by hit position...";
  #$results = capture( "sort -k 3,3 --sort=n $work_dir/hits.map -o $work_dir/hits_sorted.map" );
  #if( $EXITVAL != 0 ) {
  #  die "$0: sort exited unsuccessful";
  #}
  #print "finished.\n";
}

=head2 run_bowtie2_mapping

 Title   : run_bowtie2_mapping
 Usage   : 
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut

sub run_bowtie2_mapping {
  my $self = shift;
  my $num_threads = shift;
  my $data_basename = $self->{"data_basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie2_dir = $self->{"bowtie_db"}->{"bowtie2_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  my $reads_file_one = $self->{"bowtie_db"}->{"reads_file_one"};
  my $reads_file_two = $self->{"bowtie_db"}->{"reads_file_two"};
  my $work_dir = $self->{"work_dir"};
  print "Running mapping using bowtie2, using $num_threads threads...\n";

  # Call bowtie to run the mapping
  my $results = capture( "$bowtie2_dir/bowtie2 -x $bowtie_index_dir/$ref_genome_basename " .
                         "--threads $num_threads --reorder --sam-no-hd " .
                         "--no-discordant --no-contain --no-overlap " .
                         "-k 4 -1 $reads_file_one -2 $reads_file_two " .
                         "--al-conc $work_dir/${data_basename}_al-conc.%.fq " .
                         "--un-conc $work_dir/${data_basename}_un-conc.%.fq " .
                         "-S $work_dir/${data_basename}_alignment.sam"
                       );

  if( $EXITVAL != 0 ) {
    die "$0: bowtie2 exited unsuccessful";
  }
  print "Bowtie2 finished.\n";

}

=head2 bowtie1_identify

 Title   : bowtie1_identify
 Usage   : 
 Function: Identifies the half-mapping read pairs from Bowtie 1 run data
 Example : 
 Returns : 
 Args    : 

=cut

sub bowtie1_identify {
  my $self = shift;
  my $num_threads = shift;
  my $data_basename = $self->{"data_basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie1_dir = $self->{"bowtie_db"}->{"bowtie1_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  my $work_dir = $self->{"work_dir"};

  # Identify the half-mapping read pairs. Read pairs are already those that align uniquely to the
  # reference sequence because Bowtie 1 was run with option '-m 1'.
  print "Running bowtie to identify half-mapping read pairs...\n";
  my $results = capture( "$bowtie1_dir/bowtie $bowtie_index_dir/$ref_genome_basename " .
                         "--threads $num_threads --suppress 3,6 --time " .
                         "-m 1 $work_dir/${data_basename}_unaligned_1 " . 
                         "--al $work_dir/${data_basename}_1_halfmapping.1.fq " . 
                         #"--un $work_dir/mate1_unaligned .
                         "$work_dir/${data_basename}_mate1_hits.map"
                       );

  $results = capture( "$bowtie1_dir/bowtie $bowtie_index_dir/$ref_genome_basename " .
                         "--threads $num_threads --suppress 3,6 --time " .
                         "-m 1 $work_dir/${data_basename}_unaligned_2 " . 
                         "--al $work_dir/${data_basename}_2_halfmapping.2.fq " . 
                         #"--un $work_dir/mate2_unaligned " . 
                         "$work_dir/${data_basename}_mate2_hits.map"
                       );

  # Gather the corresponding mates that are paired with our identified half-mapping mates
  print "Gathering mates...\n";
  _find_mates("$work_dir/${data_basename}_1_halfmapping.1.fq",
              "$work_dir/${data_basename}_unaligned_2",
              "$work_dir/${data_basename}_1_halfmapping.2.fq");
  _find_mates("$work_dir/${data_basename}_2_halfmapping.2.fq",
              "$work_dir/${data_basename}_unaligned_1",
              "$work_dir/${data_basename}_2_halfmapping.1.fq");

  # Concatenate the results
  system( "cat $work_dir/${data_basename}_1_halfmapping.1.fq " .
          "$work_dir/${data_basename}_2_halfmapping.1.fq > " .
          "$work_dir/${data_basename}_halfmapping.1.fq_unsorted.tmp" );
  unlink "$work_dir/${data_basename}_1_halfmapping.1.fq" 
          or die "Can't delete $work_dir/${data_basename}_1_halfmapping.1.fq: $!";
  system( "cat $work_dir/${data_basename}_1_halfmapping.2.fq " .
          "$work_dir/${data_basename}_2_halfmapping.2.fq > " .
          "$work_dir/${data_basename}_halfmapping.2.fq_unsorted.tmp" );
  unlink "$work_dir/${data_basename}_2_halfmapping.2.fq"
          or die "Can't delete $work_dir/${data_basename}_2_halfmapping.2.fq: $!";

  # Sort our mates by ID and save the sorted output files
  _sort_fastq_by_id( "$work_dir/${data_basename}_halfmapping.1.fq_unsorted.tmp", 
                     "$work_dir/${data_basename}_halfmapping.1.fq" );
  unlink "$work_dir/${data_basename}_halfmapping.1.fq_unsorted.tmp" 
          or die "Can't delete $work_dir/${data_basename}_halfmapping.1.fq_unsorted.tmp: $!";
  _sort_fastq_by_id( "$work_dir/${data_basename}_halfmapping.2.fq_unsorted.tmp", 
                     "$work_dir/${data_basename}_halfmapping.2.fq" );
  unlink "$work_dir/${data_basename}_halfmapping.2.fq_unsorted.tmp" 
          or die "Can't delete $work_dir/${data_basename}_halfmapping.2.fq_unsorted.tmp: $!";
  print "Half-mapping read pairs are saved:\n";
  print "Created $work_dir/${data_basename}_halfmapping.1.fq\n";
  print "Created $work_dir/${data_basename}_halfmapping.2.fq\n";
}

=head2 bowtie2_identify

 Title   : bowtie2_identify
 Usage   : 
 Function: Using alignment results, identifies half-mapping read pairs with no secondary alignments
 Example : 
 Returns : 
 Args    : 

=cut

sub bowtie2_identify {
  my $self = shift;
  my $data_basename = $self->{"data_basename"};
  my $work_dir = $self->{"work_dir"};

  print "Identifying half-mapping read pairs...\n";

  # # Remove all rows we're not interested in now
  # my $results = capture( "cut -f 1,2,4,10,11 $work_dir/${data_basename}_alignment.sam " .
  #        "> $work_dir/${data_basename}_alignment_pruned" );
  # if( $EXITVAL != 0 ) {
  #   die "$0: cut exited unsuccessful";
  # }

  my $ifh = new IO::File("$work_dir/${data_basename}_alignment.sam", 'r') 
          or die "Can't open $work_dir/${data_basename}_alignment.sam: $!";
  my $ofh = IO::File->new("$work_dir/${data_basename}_halfmapping1.sam", 'w') 
          or die "Can't create $work_dir/${data_basename}_halfmapping1.sam: $!";
  my $ofh2 = IO::File->new("$work_dir/${data_basename}_halfmapping2.sam", 'w')
          or die "Can't create $work_dir/${data_basename}_halfmapping2.sam: $!";
  my $prev_id = "";
  my $discard_this_id = 0;
  my $mapping = 0;
  my @print_lines;
  while( my $line = $ifh->getline ) {
    # Ensure we're reading an alignment line and not a header line
    my $mate_id;
    my $flag_sum;
    if( $line =~ m/^(\S+)\s(\d+).*/ ) {
      $mate_id = $1;
      $flag_sum = $2;
    } else {
      next;
    }
  
    # Check if this mate ID is different from the last line
    if ($mate_id ne $prev_id) {
      # This line is the first line with a new read ID
      if (scalar @print_lines == 2 && $mapping == 1 && $prev_id ne "") {
        # Print the half-mapping mate pairs with the last line's ID
        print $ofh shift(@print_lines); # First mate in pair
        print $ofh shift(@print_lines); # Second mate in pair
      } elsif (scalar @print_lines == 3 && $mapping == 2 && $prev_id ne "") {
        # Print the trio where two mates align concordantly, and one has a secondary alignment
        foreach my $ln (@print_lines) {
          print $ofh2 $ln;
        }
        @print_lines = ();
      }
      $discard_this_id = 0;
      $prev_id = $mate_id;
    } else {
      # This line is the second or greater line with this read ID
      if ($discard_this_id) {
        next;
      }
      if (scalar @print_lines == 2 && $mapping == 0) {
        # Found a secondary alignment preceded by two half-mapping mates. Discard all lines with this ID.
        @print_lines = ();
        $discard_this_id = 1;
        next;
      } elsif (scalar @print_lines == 3) {
        # We have four alignments for this read pair. Discard all lines with this ID.
        @print_lines = ();
        $discard_this_id = 1;
        next;
      }
    }

    # At this point we have an alignment we may want to keep, so look at the sum-of-flags value to classify.
    if ( &_isMapping($flag_sum) ) {
      # Save this full-mapping line for printing
      $mapping = 2;
      push( @print_lines, $line );
    } elsif ( &_isNonMapping($flag_sum) ) {
      # We have a non-mapping ID. Discard all lines with this ID.
      $mapping = 0;
      $discard_this_id = 1;
    } elsif( &_isHalfMapping($flag_sum) ) {
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
  $ofh2->close;
}

=head2 filter

 Title   : filter
 Usage   : 
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut

sub filter {
  my $self = shift;
  my $num_threads = shift;
  my $work_dir = $self->{"work_dir"};
  my $data_basename = $self->{"data_basename"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie2_dir = $self->{"bowtie_db"}->{"bowtie2_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  my $ifh = new IO::File("$work_dir/${data_basename}_halfmapping1.sam", 'r')
          or die "$0: Can't open $work_dir/${data_basename}_halfmapping1.sam: $!";
  my $ofh = new IO::File("$work_dir/${data_basename}_fake_paired_end.1.fq", 'w')
          or die "$0: Can't open $work_dir/${data_basename}_fake_paired_end.1.fq: $!";
  my $ofh2 = new IO::File("$work_dir/${data_basename}_fake_paired_end.2.fq", 'w')
          or die "$0 Can't open $work_dir/${data_basename}_fake_paired_end.2.fq: $!";
  # Define length of each mate for the fake read pairs
  my $mate_length = 16;

  # Get the mates we want and write them as fake paired-end reads to FastQ
  my $read_length;
  while( my $line = $ifh->getline ) {
    my ( $id, $flag_sum, $sequence, $quality_scores );
    if( $line =~ m/^(\S+)\s(\d+)\s\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s(\S+)\s(\S+)\s.*/ ) {
      $id = $1;
      $flag_sum = $2;
      $sequence = $3;
      $quality_scores = $4;
    } else {
      next;
    }

    # Only create fake pairs for the unaligned mate in the pair
    if( &_isUnalignedMate($flag_sum) ) {
      # Reverse complement reads that aligned to the reverse strand, to get original read
      if(($flag_sum & 16) == 1) {
        $sequence = &_reverseComplement($sequence);
      }

      # Create the fake mates. Mate 2 gets reverse complemented.
      my $mate_1 = substr($sequence, 0, $mate_length);
      my $mate_2 = &_reverseComplement(substr($sequence, -1 * $mate_length));
      my $quality_scores_1 = substr($quality_scores, 0, $mate_length);
      my $quality_scores_2 = scalar reverse substr($quality_scores, -1 * $mate_length);

      # Write the fake mates
      print $ofh "\@$id\/1\n$mate_1\n+\n$quality_scores_1\n";
      print $ofh2 "\@$id\/2\n$mate_2\n+\n$quality_scores_2\n";

      # Record read length for fake pair alignment parameters
      $read_length = length($sequence);
    }
  }
  $ifh->close;
  $ofh->close;
  $ofh2->close;

  # Define the maximum/minimum insert length as read length +/- 10
  if (!defined $read_length) {
    print STDERR "$0: filter(): No half-mapping pairs available.\n";
    return 1;
  }
  my $minins = $read_length - 10;
  if( $minins < 0 ) {
    $minins = 0;
  }
  my $maxins = $read_length + 10;
  print "using minins $minins, maxins $maxins\n";

  # Run Bowtie 2 with the fake read pairs as input
  my $results = capture( "$bowtie2_dir/bowtie2 -x $bowtie_index_dir/$ref_genome_basename " .
                         "--threads $num_threads --time --reorder --sam-no-hd " .
                         "--minins $minins --maxins $maxins " .
                         "--no-mixed --no-discordant --no-contain --no-overlap " .
                         "-1 $work_dir/${data_basename}_fake_paired_end.1.fq " . 
                         "-2 $work_dir/${data_basename}_fake_paired_end.2.fq " .
                         "--al-conc $work_dir/${data_basename}_fake_pairs_al-conc.%.fq " .
                         "-S $work_dir/${data_basename}_fake_pairs_aligned.sam"
                       );
  if ( $EXITVAL != 0 ) {
    die "$0: bowtie2 exited unsuccessful";
  }

  # Find reads representing insertions/deletions by considering all the fake pairs that aligned
  $ifh = new IO::File("$work_dir/${data_basename}_halfmapping1.sam", 'r')
          or die "$0: Can't open $work_dir/${data_basename}_halfmapping1.sam: $!"; 
  my $ifh2 = new IO::File("$work_dir/${data_basename}_fake_pairs_aligned.sam", 'r')
          or die "$0: Can't open $work_dir/${data_basename}_fake_pairs_aligned.sam: $!";
  $ofh = new IO::File("$work_dir/${data_basename}_filtered.sam", 'w')
          or die "$0: Can't open $work_dir/${data_basename}_filtered.sam: $!";
  $ofh2 = new IO::File("$work_dir/${data_basename}_fake_pairs_aligned.debug", 'w')
          or die "$0: Can't open $work_dir/${data_basename}_fake_pairs_aligned.debug: $!";
  while( my $line = $ifh2->getline ) {
    my $mate_id;
    my $flag_sum;
    my $pos;
    if( $line =~  m/^(\S+)\s(\d+)\s\S+\s(\S+).*/ ) {
      $mate_id = $1;
      $flag_sum = $2;
      $pos = $3;
    } else {
      next;
    }

    # Discard the fake pairs that aligned discordantly
    if( !&_isMapping($flag_sum) ) {
      next;
    }

    # Get the corresponding alignment from the original half-mapping reads
    while (my $orig_line = $ifh->getline) {
      my $orig_id;
      my $orig_flag_sum;
      my $orig_pos;
      if( $orig_line =~ m/^(\S+)\s(\d+)\s\S+\s(\S+).*/ ) {
        $orig_id = $1;
        $orig_flag_sum = $2;
        $orig_pos = $3;
      } else {
        next;
      }
      while( $orig_id ne $mate_id ) {
        $orig_line = $ifh->getline;
        if ( $orig_line =~ m/(\S+)\s(\d+)\s\S+\s(\S+).*/ ) {
          $orig_id = $1;
          $orig_flag_sum = $2;
          $orig_pos = $3;
        }
      }
      # Get the line with the second mate in the fake pair
      my $line2 = $ifh2->getline;
      # Get second mate from the half-mapping read pair
      my $orig_line2 = $ifh->getline;

      print $ofh2 "\nFAKE PAIR ALIGNMENT LINE 1: $line";
      print $ofh2 "FAKE PAIR ALIGNMENT LINE 2: $line2";
      print $ofh2 "ORIGINAL HALF-MAPPING PAIR LINE 1: $orig_line";
      print $ofh2 "ORIGINAL HALF-MAPPING PAIR LINE 2: $orig_line2";

      # Discard fake pairs that align close to the originally aligning mate in the half-mapping read pair
      if( $pos - $orig_pos < $maxins || $orig_pos - $pos < $maxins ) {
        print $ofh2 "DISCARDED -- $pos - $orig_pos < +/- $maxins\n";
        last;
      }

      # Keep only the unaligned mate from the half-mapping read pair
      if( &_isUnalignedMate($orig_flag_sum) ) {
        print $ofh $orig_line;
      } else {
        print $ofh $orig_line2;
      } 
      last;
    }
  }
  $ifh->close;
  $ifh2->close;
}

############ Subroutines for internal use by this module ############

# _find_mates: Identify corresponding (paired) mates for the given file of mates
# Arguments: 1. Unmatched mates file 2. File with candidate mates to search 3. Output file
# Note: Matches in output file are unsorted, so may be listed in a different order.
sub _find_mates {
  my $unpaired_mates = shift;
  my $candidate_mates = shift;
  my $output_file = shift;

  # Open the unpaired mates file
  my $ifh = new IO::File($unpaired_mates, 'r') or die "$0: Can't open $unpaired_mates: $!";

  # Loop through the unpaired mates file and get all the search IDs
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

      # Append this ID to the combined search string
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

# Returns the reverse complement of the given string
sub _reverseComplement {
  my ($value) = shift;
  $value = scalar reverse $value;
  for(0..length($value)-1) { substr($value, $_, 1) = &_complement(substr($value, $_, 1)); }
  return $value;
}

# Returns the complement of the given string
sub _complement {
  my %complementMap = (
    "A" => "T", "T" => "A", "a" => "t", "t" => "a",
    "C" => "G", "G" => "C", "c" => "g", "g" => "c",
    "R" => "Y", "Y" => "R", "r" => "y", "y" => "r",
    "M" => "K", "K" => "M", "m" => "k", "k" => "m",
    "S" => "S", "W" => "W", "s" => "s", "w" => "w",
    "B" => "V", "V" => "B", "b" => "v", "v" => "b",
    "H" => "D", "D" => "H", "h" => "d", "d" => "h",
    "N" => "N", "." => ".", "n" => "n" );
  my $complementedString = $complementMap{$_[0]} or die "$0: Can't get reverse complement for '$_[0]'";
  return $complementedString;
}

# Returns 1 if the given sum-of-flags value identifies a mate in a half-mapping read pair, otherwise returns 0
sub _isHalfMapping {
  my $flag = shift;
  if( $flag == 69 ||
      $flag == 73 ||
      $flag == 89 ||
      $flag == 101 ||
      $flag == 133 ||
      $flag == 137 ||
      $flag == 153 ||
      $flag == 165 ) {
    return 1;
  } else {
    return 0;
  }
}

# Returns 1 if the given sum-of-flags value identifies a mate in a concordant alignment, otherwise returns 0
sub _isMapping {
  my $flag = shift;
  if( $flag == 83 ||
      $flag == 99 ||
      $flag == 147 ||
      $flag == 163 ) {
    return 1;
  } else {
    return 0;
  }
}

# Returns 1 if the given sum-of-flags value identifies a mate in a no-alignments pair, otherwise returns 0
sub _isNonMapping {
  my $flag = shift;
  if( $flag == 77 || $flag == 141 ) {
    return 1;
  } else {
    return 0;
  }
}

# Returns 1 if the given sum-of-flags value identifies an unaligned mate in a half-mapping pair, otherwise returns 0
sub _isUnalignedMate {
  my $flag = shift;
  if( $flag == 69 ||
      $flag == 101 ||
      $flag == 133 ||
      $flag == 165 ) {
    return 1;
  } else {
    return 0;
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

# Checks if the given folder exists, and removes trailing slash from string
sub _check_dir {
  my $dir = shift;
  if (-d $dir) {
    $dir = $1 if ($dir =~ /(.*)\/$/);
    return $dir;
  } else {
    die("$0: directory $dir does not exist");
  }
}

# Creates the given directory, and gives an error message on failure
sub _make_dir {
  my( $dirname ) = @_;
  mkdir $dirname, 0755 || die "$0: could not create $dirname";
}

# Returns a date/time string with no whitespace
sub _datestamp {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $datestring = sprintf("%4d%02d%02d_%02d_%02d_%02d",$year+1900,$mon+1,
                           $mday,$hour,$min,$sec);
}

1;
