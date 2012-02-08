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
#use Bio::Tools;
use File::Basename;
#use IPC::System::Simple qw(capture $EXITVAL);

=head2 new

 Title   : new
 Usage   : $project = IntronPoly->new()
 Function: Initializes a new intron polymorphism analysis project
 Example : $project = IntronPoly->new();
           $project->set_workdir("/home/theis/work");
           
 Returns : IntronPoly project object
 Args    : None

=cut

sub new {
  my $self = {};
  $self->{"workdir"} = undef; # working directory
  $self->{"outputdir"} = undef; # output directory
  $self->{"toolpaths"} = undef; # hash of paths to tool directories
  $self->{"db"} = undef; # hash of text file data directories
  bless($self);
  return $self;
}

sub set_workdir {
  my $self = shift;
  my $path = shift;

  # Ensure path has a trailing slash
  unless ( $path =~ m/\/$/ ) { $path = $path . '/'; }

  # Set the working directory
  my $outputdir = $path . &datestamp . "_output";
  $self->{"workdir"} = $path;
  $self->{"outputdir"} = $outputdir;
  &makedir( $outputdir );
  print "Setting working directory to $path\n";
  print "Created output directory $outputdir\n";
  return $self->{"workdir"};
}

sub set_tooldirs {
  my $self = shift;
  my $path = shift;

  # Set the path to bowtie
  $self->{"toolpaths"}->{"bowtie"} = $path;
}

sub build_db {
  my $self = shift;
  my $workdir = $self->{"workdir"};
  $self->{"db"} = {
                   # Path to Bowtie database created by bowtie-build
                   bowtiedb => $workdir . "/bowtie-db/",
                   # Base name of reference genome for bowtie
                   bowtie_ref_base => undef,
                   # Path to reads for mapping step
                   reads => $workdir . "/reads/",
                  };
}

sub build_mapping_db {
  my $self = shift;
  my $input_ref = shift;
  my $workdir = $self->{"workdir"};
  my $bowtie_path = $self->{"toolpaths"}->{"bowtie"};
  my $bowtie_db = $self->{"db"}->{"bowtiedb"};

  # Get the base name, without extension, of the reference genome
  my ($ref_file, $ref_dir, $ref_ext) = fileparse($input_ref, qr/\.[^.]*/);
  $self->{"db"}->{"bowtie_ref_base"} = $ref_file;

  # TODO change this to check for actual ref file name in db dir, not just any files
  unless (-e "$bowtie_db/$ref_file.1.ebwt") {
    unless (-d $bowtie_db) { &makedir($bowtie_db); }

    print "bowtie database will be created in " . $bowtie_db . "\n";
    system("$bowtie_path/bowtie-build $input_ref $bowtie_db/$ref_file");

    # TODO check exit val

  } else {
    print "bowtie database already exists in $bowtie_db/$ref_file.*, not re-creating.\n";
  }

}

sub run_mapping {
  my $self = shift;
  my $input_base = shift;
  my $input = $self->{"db"}->{"reads"};
  my $bowtie_db = $self->{"db"}->{"bowtiedb"};
  my $bowtie_path = $self->{"toolpaths"}->{"bowtie"};
  my $bowtie_ref_base = $self->{"db"}->{"bowtie_ref_base"};

  print "bowtie path is $bowtie_path\n";
  print "input base is $input_base\n";

  system("$bowtie_path/bowtie $bowtie_db/$bowtie_ref_base -1 $input/$input_base.1.fq -2 $input/$input_base.2.fq" );

  #if( $EXITVAL != 0 ) {
  #  warn("Error running mapping\n");
  #  exit(0);
  #}
}

###################################################################################

# Returns 1 if the given folder is empty, and creates it if it does not exist
sub is_folder_empty {
    my $dirname = shift;
    unless (-d $dirname) {
      # Create the folder
      &makedir($dirname);
      return 1;
    }
    opendir(my $dh, $dirname) or die "$dirname: $!";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

# Creates the given directory, and gives an error message on failure
sub makedir {
  my( $dirname ) = @_;
  mkdir $dirname, 0755 || die "$0: could not create $dirname\n";
}

# Returns a date/time string with no whitespace
sub datestamp {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $datestring = sprintf("%4d%02d%02d_%02d_%02d_%02d",$year+1900,$mon+1,
                           $mday,$hour,$min,$sec);
}

1;
