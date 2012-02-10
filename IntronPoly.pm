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
use Bio::Tools::GuessSeqFormat;
use File::Basename;
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
  $self->{"work_dir"} = undef; # Scalar path to work directory
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
           full path to work folder from a previous run (or containing existing unmapped reads)

=cut

sub set_work_dir {
  my $self = shift;
  my $scripts_dir = shift;
  my $work_dir = shift || "$scripts_dir/run-" . &datestamp;
  unless (-e $work_dir) {  
    $work_dir = $1 if ($work_dir =~ /(.*)\/$/);
    &make_dir( $work_dir );
    print "Created work directory $work_dir\n"
  } else {
    $work_dir = &check_dir( $work_dir );
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

  # Split reference genome pathname into components
  my ($file, $dir, $ext) = fileparse($ref_genome_filename, qr/\.[^.]*/);
  $self->{"ref_genome"}->{"full_pathname"} = $ref_genome_filename;
  $self->{"ref_genome"}->{"basename"} = $file;
  $self->{"ref_genome"}->{"dir"} = $dir;
};

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
  my $bowtie_dir = shift;
  my $bowtie_index_dir = shift;
  my $reads_dir = shift;
  my $reads_basename = shift;
  my $work_dir = $self->{"work_dir"};

  # Ensure reference genome appears to be FastA format
  my $ref_genome = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $ref_genome_dir = $self->{"ref_genome"}->{"dir"};
  unless (-d $ref_genome_dir) {
    die "$0: reference genome directory $ref_genome_dir does not exist" ;
  }
  my $guesser = Bio::Tools::GuessSeqFormat->new();
  $guesser->file( $ref_genome );
  my $format = $guesser->guess;
  if ($format ne "fasta") {
    die "$0: reference genome $ref_genome does not appear to be a valid FastA file";
  }

  # Ensure that reads files exist (but don't verify format--GuessSeqFormat can't get it right)
  $reads_dir = $1 if ($reads_dir =~ /(.*)\/$/);
  my $reads_file_one = "$reads_dir/$reads_basename.1.fq";
  my $reads_file_two = "$reads_dir/$reads_basename.2.fq";
  unless (-e $reads_file_one) {
    die "$0: reads file $reads_file_one does not exist";
  }
  unless (-e $reads_file_two) {
    die "$0: reads file $reads_file_two does not exist";
  }

  # TODO Ensure number of reads in pair matches between file 1 and 2

  # Create directory for bowtie index files if necessary
  $bowtie_index_dir = $1 if ($bowtie_index_dir =~ /(.*)\/$/);
  unless (-d $bowtie_index_dir) {
    &make_dir( $bowtie_index_dir );
  }

  # Remember paths needed for running bowtie
  $self->{"bowtie_db"} = {
                          bowtie_dir => $bowtie_dir,
                          bowtie_index_dir => $bowtie_index_dir,
                          reads_file_one => $reads_file_one,
                          reads_file_two => $reads_file_two
                         };
  print "Using reference genome $ref_genome\n";
  print "Using reads file 1: $reads_file_one\n";
  print "Using reads file 2: $reads_file_two\n";
}

=head2 build_mapping_index

 Title   : build_mapping_index
 Usage   :
 Function:
 Example : 
 Returns : 
 Args    :

=cut

sub build_mapping_index {
  my $self = shift;
  my $ref_genome = $self->{"ref_genome"}->{"full_pathname"};
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie_dir = $self->{"bowtie_db"}->{"bowtie_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  unless (-e "$bowtie_index_dir/$ref_genome_basename.1.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.2.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.3.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.4.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.rev.1.ebwt" &&
             "$bowtie_index_dir/$ref_genome_basename.rev.2.ebwt") {
    # Call bowtie-build to create the bowtie index
    print "Creating bowtie index in $bowtie_index_dir...";
    my $results = capture( "$bowtie_dir/bowtie-build " .
                           "$ref_genome $bowtie_index_dir/$ref_genome_basename" );
    if ( $EXITVAL != 0 ) {
      die "$0: bowtie-build exited unsuccessful";
    }
    print "finished.\n";
  } else {
    print "Bowtie index already exists, not re-creating.\n";
  }
}

=head2 run_mapping

 Title   : run_mapping
 Usage   : 
 Function: 
 Example : 
 Returns : 
 Args    : 

=cut
 
sub run_mapping {
  my $self = shift;
#  my $threads = shift; # Note can't use mult threads with --refout
  my $ref_genome_basename = $self->{"ref_genome"}->{"basename"};
  my $bowtie_dir = $self->{"bowtie_db"}->{"bowtie_dir"};
  my $bowtie_index_dir = $self->{"bowtie_db"}->{"bowtie_index_dir"};
  my $reads_file_one = $self->{"bowtie_db"}->{"reads_file_one"};
  my $reads_file_two = $self->{"bowtie_db"}->{"reads_file_two"};
  my $work_dir = $self->{"work_dir"};
  print "Running mapping using bowtie...\n";

  # Call bowtie to run the mapping
  my $results = capture( "$bowtie_dir/bowtie $bowtie_index_dir/$ref_genome_basename " .
                         "-1 $reads_file_one -2 $reads_file_two " .
                         "--al $work_dir/aligned --un $work_dir/unaligned --refout"
                       );

  if( $EXITVAL != 0 ) {
    die "$0: bowtie exited unsuccessful";
  }
  print "Bowtie finished.\n";
}


###################################################################################

# Checks if the given folder exists, and removes trailing slash from string
sub check_dir {
  my $dir = shift;
  if (-d $dir) {
    $dir = $1 if ($dir =~ /(.*)\/$/);
    return $dir;
  } else {
    die("$0: directory $dir does not exist");
  }
}

# Creates the given directory, and gives an error message on failure
sub make_dir {
  my( $dirname ) = @_;
  mkdir $dirname, 0755 || die "$0: could not create $dirname";
}

# Returns a date/time string with no whitespace
sub datestamp {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $datestring = sprintf("%4d%02d%02d_%02d_%02d_%02d",$year+1900,$mon+1,
                           $mday,$hour,$min,$sec);
}

1;
