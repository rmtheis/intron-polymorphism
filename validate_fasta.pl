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
use Getopt::Long;
use IO::File;

#
# Performs basic validation on files in FastA format.
#

my $usage_msg = "Performs basic validation on a Fasta file.\n"
              . "Usage: validate_fasta.pl -i fasta_file\n";
die $usage_msg unless ( @ARGV );
my $input_file;
GetOptions( "i=s" =>\$input_file ) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );

$input_file =~ s/^~/$ENV{HOME}/;
my $ifh = new IO::File( $input_file, 'r' ) or die "Can't open $input_file: $!";
my $first_line = $ifh->getline;
if ( $first_line !~ m/^>\s*\w/ ) {
  print STDERR "Validation error: First FastA line is not a comment: $first_line";
  print STDERR "File: $input_file\n";
  die;
}

while( my $line = $ifh->getline ) {

    if( !($line =~ /^[A-IK-NP-Z]+$/i || $line =~ /^>\s*\w/ || $line =~ /^\s*$/) ) {
      print STDERR "Validation error: FastA line is not valid sequence or comment: $line";
      print STDERR "File: $input_file\n";
      die;
    }

}
$ifh->close;
exit;
