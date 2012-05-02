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

use constant DEBUG => (1); # Flag to print info for debugging (temporary)

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
  my $reads_file_one = "$reads_dir/${data_basename}_1.fq";
  my $reads_file_two = "$reads_dir/${data_basename}_2.fq";
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
    print "Creating bowtie index in $bowtie_index_dir...\n";
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
    print "Creating bowtie index in $bowtie_index_dir...\n";
    my $results = capture( "$bowtie2_dir/bowtie2-build " .
                           "$ref_genome $bowtie_index_dir/$ref_genome_basename" );
    if ( $EXITVAL != 0 ) {
      die "$0: bowtie2-build exited unsuccessful";
    }
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

  # Sort aligning reads by match position
  #print "Sorting aligned reads by hit position...\n";
  #$results = capture( "sort -k 3,3 --sort=n $work_dir/hits.map -o $work_dir/hits_sorted.map" );
  #if( $EXITVAL != 0 ) {
  #  die "$0: sort exited unsuccessful";
  #}
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
                         "--maxins 500 --minins 100 " .
                         "--no-discordant " . 
                         "--no-contain --no-overlap " .
                         "-k 3 -1 $reads_file_one -2 $reads_file_two " .
                         "--al-conc $work_dir/${data_basename}_al-conc.%.fq " .
                         "--un-conc $work_dir/${data_basename}_un-conc.%.fq " .
                         "-S $work_dir/${data_basename}_alignment.sam"
                       );

  if( $EXITVAL != 0 ) {
    die "$0: bowtie2 exited unsuccessful";
  }

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
                         "--threads $num_threads --suppress 3,6 " .
                         "-m 1 $work_dir/${data_basename}_unaligned_1 " . 
                         "--al $work_dir/${data_basename}_1_halfmapping.1.fq " . 
                         #"--un $work_dir/mate1_unaligned .
                         "$work_dir/${data_basename}_mate1_hits.map"
                       );

  $results = capture( "$bowtie1_dir/bowtie $bowtie_index_dir/$ref_genome_basename " .
                         "--threads $num_threads --suppress 3,6 " .
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

# =head2 find_distant_alignments
# 
#  Title   : find_distant_alignments
#  Usage   : 
#  Function: Run Bowtie to find read pairs where both mates align concordantly, but far from one another
#  Example : 
#  Returns : 
#  Args    : 
# 
# =cut
# 
# sub find_distant_alignments {
#   my $self = shift;
#   my $num_threads = shift;
#   my $data_basename = $self->{"data_basename"};
#   my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
#   my $bowtie2_dir = $self->{"bowtie_db"}->{"bowtie2_dir"};
#   my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
#   my $reads_file_one = $self->{"bowtie_db"}->{"reads_file_one"};
#   my $reads_file_two = $self->{"bowtie_db"}->{"reads_file_two"};
#   my $work_dir = $self->{"work_dir"};
#   print "Running mapping using bowtie2 to find distant alignments, using $num_threads threads...\n";
# 
#   # Combine the read pairs files into one file with unique IDs for unpaired mapping
#   my $unpaired_reads_file = "$work_dir/${data_basename}_unpaired.fq";
#   &_combine_fastq( $reads_file_one, $reads_file_two, $unpaired_reads_file );
# 
# 
#   # Call bowtie to run the mapping as unpaired reads
#   my $results = capture( "$bowtie2_dir/bowtie2 -x $bowtie_index_dir/$ref_genome_basename " .
#                          "--threads $num_threads --reorder --sam-no-hd " .
#                          "--no-discordant --no-contain --no-overlap " .
#                          "-k 3 -U $reads_file_one,$reads_file_two " .
#                          #"--al-conc $work_dir/${data_basename}_al-conc.%.fq " .
#                          #"--un-conc $work_dir/${data_basename}_un-conc.%.fq " .
#                          "-S $work_dir/${data_basename}_halfmapping3.sam.tmp"
#                        );
#
#  if( $EXITVAL != 0 ) {
#    die "$0: bowtie2 exited unsuccessful";
#  }
# 
#   # Parse the alignments file to compare distances between reads
# 
#   print "debug: find_distant_alignments() done.\n";
# }

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
  my $ofh = new IO::File("$work_dir/${data_basename}_fake_paired_end_1.fq", 'w')
          or die "$0: Can't open $work_dir/${data_basename}_fake_paired_end_1.fq: $!";
  my $ofh2 = new IO::File("$work_dir/${data_basename}_fake_paired_end_2.fq", 'w')
          or die "$0 Can't open $work_dir/${data_basename}_fake_paired_end_2.fq: $!";
  # Define length of each mate for the fake read pairs
  my $mate_length = 16;

  # Define length of read pair for fake paired-end alignment
  my $read_length = 500;

  # Boolean to track whether we have any half-mapping pairs
  my $have_candidates = 0;

  # Get the mates we want and write them as fake paired-end reads to FastQ
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

    # Create a fake pair from the unaligned mate in the pair
    if( &_isUnalignedMate($flag_sum) ) {
      # Set a flag indicating that we got at least one half-mapping pair
      $have_candidates = 1;

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
    }
  }
  $ifh->close;
  $ofh->close;
  $ofh2->close;

  # Define the maximum/minimum insert length as read length +/- 10
  if (!$have_candidates) {
    print STDERR "$0: filter(): No half-mapping pairs available.\n";
    return 1;
  }

  # Run Bowtie 2 with the fake read pairs as input
  my $results = capture( "$bowtie2_dir/bowtie2 -x $bowtie_index_dir/$ref_genome_basename " .
                         "--threads $num_threads --reorder --sam-no-hd " .
                         "--no-mixed --no-discordant --no-contain --no-overlap " .
                         "-1 $work_dir/${data_basename}_fake_paired_end_1.fq " . 
                         "-2 $work_dir/${data_basename}_fake_paired_end_2.fq " .
                         #"--al-conc $work_dir/${data_basename}_fake_pairs_al-conc.%.fq " .
                         "-S $work_dir/${data_basename}_fake_pairs_alignment.sam"
                       );
  if ( $EXITVAL != 0 ) {
    die "$0: bowtie2 exited unsuccessful";
  }

  # Find reads representing insertions/deletions by considering all the fake pairs that aligned
  $ifh = new IO::File("$work_dir/${data_basename}_fake_pairs_alignment.sam", 'r')
          or die "$0: Can't open $work_dir/${data_basename}_fake_pairs_alignment.sam: $!"; 
  my $ifh2 = new IO::File("$work_dir/${data_basename}_halfmapping1.sam", 'r')
          or die "$0: Can't open $work_dir/${data_basename}_halfmapping1.sam: $!";
  $ofh = new IO::File("$work_dir/${data_basename}_filtered.sam", 'w')
          or die "$0: Can't open $work_dir/${data_basename}_filtered.sam: $!";
  $ofh2 = new IO::File("$work_dir/${data_basename}_fake_pairs_alignment.debug", 'w')
          or die "$0: Can't open $work_dir/${data_basename}_fake_pairs_alignment.debug: $!" if DEBUG;

  # Get an alignment from the originally half-mapping pairs
  while( my $line = $ifh2->getline ) {
    my $orig_id;
    my $orig_flag_sum;
    my $other_mate_pos;
    if( $line =~  m/^(\S+)\s(\d+)\s\S+\s\S+\s\S+\s\S+\s\S+\s(\S+)/ ) {
      $orig_id = $1;
      $orig_flag_sum = $2;
      $other_mate_pos = $3;
    } else {
      next;
    }

    # Consider only the unaligned mate from the half-mapping read pair
    if( !&_isUnalignedMate($orig_flag_sum) ) {
      next;
    }

    # Get the corresponding alignment from the fake pairs
    while (my $fake_line = $ifh->getline) {
      my $id;
      my $pos;
      if( $fake_line =~ m/^(\S+)\s\d+\s\S+\s(\S+).*/ ) {
        $id = $1;
        $pos = $2;
      } else {
        next;
      }
      while( $id ne $orig_id ) {
        $fake_line = $ifh->getline;
        if( $fake_line =~ m/(\S+)\s\d+\s\S+\s(\S+).*/ ) {
          $id = $1;
          $pos = $2;
        }
      }

      print $ofh2 "\nFAKE PAIR ALIGNMENT LINE 1: $fake_line" if DEBUG;
      print $ofh2 "ORIGINAL HALF-MAPPING PAIR UNALIGNED MATE: $line" if DEBUG;
      print $ofh2 "ORIG_POS: $other_mate_pos POS: $pos\n" if DEBUG;

      # Discard fake pairs that align close to the originally aligning mate in the half-mapping read pair
      if( $pos != 0 && ($pos < ($other_mate_pos + $read_length + 10) && $pos > ($other_mate_pos + $read_length - 10))) {
        print $ofh2 "DISCARDED because $pos too close to $other_mate_pos plus $read_length plus or minus ten\n" if DEBUG;
        last;
      }
      print $ofh2 "RETAINED: pos = $pos, other_mate_pos = $other_mate_pos\n" if DEBUG;

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
 Usage   : 
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut

sub group {
  my $self = shift;
  my $work_dir = $self->{"work_dir"};
  my $data_basename = $self->{"data_basename"};

  # Expected intron length
  my $intron_length = 250;

  # Read length
  my $read_length = 500;

  # Mate length
  my $mate_length = 125; # TODO get this from SAM

  # Sort by alignment positionn
  my $results = capture( "cat $work_dir/${data_basename}_filtered.sam " .
                         "$work_dir/${data_basename}_halfmapping2.sam " .
                         "> $work_dir/${data_basename}_halfmapping_all.sam" );
  if ( $EXITVAL != 0 ) {
    die "$0: cat exited unsuccessful";
  }
  $results = capture( "sort -k 4,4 -n -o $work_dir/${data_basename}_halfmapping_all_sorted.sam " .
                         "$work_dir/${data_basename}_halfmapping_all.sam" );
  if ( $EXITVAL != 0 ) {
    die "$0: sort exited unsuccessful";
  }

  # Get the clusters one at a time
  my $ifh = new IO::File("$work_dir/${data_basename}_filtered.sam", 'r')
          or die "$0: Can't open $work_dir/${data_basename}_filtered.sam: $!";
  my @cluster;
  my $have_cluster = 0;
  my $last_pos = 0;
  my $overlap_dist = 2 * (($read_length - 2 * $mate_length) + $intron_length);
 
  my $id;
  my $pos; 
  while( my $line = $ifh->getline ) {
    if( $line =~  m/^(\S+)\s\d+\s\S+\s\S+\s\S+\s\S+\s\S+\s(\S+)/ ) {
      $id = $1;
      $pos = $2;
    } else {
      next;
    }
    
    # Add or discard based on position

    # Run velvet on the reads one at a time using stdout

    # velveth $work_dir/velvet-results-separate-dir 13 -sam -shortPaired 

  }
  $ifh->close;

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

# Combines two FastQ files, appending "-1" and "-2" to the read IDs for uniqueness
sub _combine_fastq {
  my $input_file1 = shift;
  my $input_file2 = shift;
  my $output_file = shift;
  my $ifh1 = new IO::File($input_file1, 'r') or die "$0: Can't open $input_file1: $!";
  my $ifh2 = new IO::File($input_file2, 'r') or die "$0: Can't open $input_file2: $!";
  my $ofh = new IO::File($output_file, 'w') or die "$0: Can't open $output_file: $!";

  while ( my $line = $ifh1->getline ) {
    print $ofh $line; #TODO change the IDs. Put "-1" before any "/1"
  }
  while ( my $line = $ifh2->getline ) {
    print $ofh $line; #TODO change the IDs.
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
