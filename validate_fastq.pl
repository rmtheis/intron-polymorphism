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
# Performs basic validation on paired read files in FastQ format. Ensures that:
#
#   1. Read IDs and their order match between mate 1 and mate 2
#   2. There are no tabs anywhere in the files
#   3. There are no pipes '|' on read ID lines
#   4. Header lines do not start with '@'
#
# If any problems are found, this script alerts the user and immediately exits. No changes are made
# to the input files.
#

my $usage_msg = "Performs basic validation on a Fastq file.\n"
              . "Usage: validate_fastq.pl -1 fastq_file_1 -2 fastq_file_2\n";
die $usage_msg unless ( @ARGV );
my ( $mate1s, $mate2s );
GetOptions(
  "1=s" =>\$mate1s,
  "2=s" =>\$mate2s,
  ) || die "$0: Bad option";
die $usage_msg unless ( defined $mate1s && defined $mate2s );

if ( `wc -l $mate1s | cut -d ' ' -f 1` != `wc -l $mate2s | cut -d ' ' -f 1` ) {
  print STDERR "Validation error: number of lines does not match between files\n";
  print STDERR "Files: $mate1s, $mate2s\n";
  die;
}

$mate1s =~ s/^~/$ENV{HOME}/;
$mate2s =~ s/^~/$ENV{HOME}/;
my $ifh1 = new IO::File( $mate1s, 'r' ) or die "Can't open $mate1s: $!";
my $ifh2 = new IO::File( $mate2s, 'r' ) or die "Can't open $mate2s: $!";

my $line_count = 0;
while ( my $line_1 = $ifh1->getline ) {

  my $line_2 = $ifh2->getline;

  if ( $line_count == 0 ) {
    if ( substr( $line_1, 0, 1 ) ne "@" ) {
      print STDERR "Validation error: FastQ ID does not begin with \'@\': $line_1";
      print STDERR "File: $mate1s\n";
      die;
    }

    if ( $line_1 =~ m/\|/ ) {
      print STDERR "Validation error: Pipe character (\"\|\") in ID: $line_1";
      print STDERR "File: $mate1s\n";
      die;
    }
    if ( $line_2 =~ m/\|/ ) {
      print STDERR "Validation error: Pipe character (\"\|\") in ID: $line_2";
      print STDERR "File: $mate2s\n";
      die;
    }
  }

  if ( $line_1 =~ m/\t/ ) {
    print STDERR "Validation error: FastQ read file contains a tab: $mate1s\n";
    die;
  }
  if ( $line_2 =~ m/\t/ ) {
    print STDERR "Validation error: FastQ read file contains a tab: $mate2s\n";
    die;
  }
  $line_count = ($line_count + 1) % 4;
}
$ifh1->close;
$ifh2->close;
exit;
