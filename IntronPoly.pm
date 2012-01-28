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
#use IPC::System::Simple qw(capture $EXITVAL);

=head2 new

 Title   : new
 Usage   : $project = Intron->new()
 Function: Initializes a new intron polymorphism analysis project
 Example : $project = Intron->new();
           $project->set_workdir("/home/theis/work");
           $project->set_sample("Reads_Name");
 Returns : Intron project object
 Args    : None

=cut

sub new {
  my $self = {};
  $self->{"workdir"} = undef; # working directory
  $self->{"outputdir"} = undef; # output directory
  bless($self);
  return $self;
}

sub set_workdir {
  my $self = shift;
  my $path = shift;

  # Ensure path has a trailing slash
  unless ( $path =~ m/\/$/ ) { $path = $path . '/'; }

  my $outputdir = $path . &datestamp . "_output";
  $self->{"workdir"} = $path;
  $self->{"outputdir"} = $outputdir;
  &makedir( $outputdir );
  print "Setting working directory to $path\n";
  print "Created output directory $outputdir\n";
  return $self->{"workdir"};
}

sub run_mapping {
  my ( $self ) = @_;

  #if( $EXITVAL != 0 ) {
  #  warn("Error running mapping\n");
  #  exit(0);
  #}
}

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
