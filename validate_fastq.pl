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
              . "Usage: validate_fastq.pl -i fastq_file_basename\n"
              . "(Basename is filename without '_1.fq' or '_2.fq' extension)\n";
unless ( @ARGV ) {
  print $usage_msg;
  exit;
}
my $base_name;
GetOptions( "i=s" =>\$base_name );
unless ( defined $base_name ) {
  print $usage_msg;
  exit;  
}

my $input_file_1 = $base_name . "_1.fq";
my $input_file_2 = $base_name . "_2.fq";
my $ifh1 = new IO::File( $input_file_1, 'r' ) or die "Can't open $input_file_1: $!";
my $ifh2 = new IO::File( $input_file_2, 'r' ) or die "Can't open $input_file_2: $!";

if ( `wc -l $input_file_1 | cut -d ' ' -f 1` != `wc -l $input_file_2 | cut -d ' ' -f 1` ) {
  print STDERR "Validation error: number of lines does not match between files\n";
  print STDERR "Files: $input_file_1, $input_file_2\n";
  die;
}

my $line_count = 0;
while ( my $line_1 = $ifh1->getline ) {

  my $line_2 = $ifh2->getline;

  if ( $line_count == 0 ) {
    if ( substr( $line_1, 0, 1 ) ne "@" ) {
      print STDERR "Validation error: FastQ ID does not begin with \'@\': $line_1";
      print STDERR "File: $input_file_1\n";
      die;
    }

    if ( $line_1 =~ m/\|/ ) {
      print STDERR "Validation error: Pipe character (\"\|\") in ID: $line_1";
      print STDERR "File: $input_file_1\n";
      die;
    }
    if ( $line_2 =~ m/\|/ ) {
      print STDERR "Validation error: Pipe character (\"\|\") in ID: $line_2";
      print STDERR "File: $input_file_2\n";
      die;
    }
  }

  if ( $line_1 =~ m/\t/ ) {
    print STDERR "Validation error: FastQ read file contains a tab: $input_file_1\n";
    die;
  }
  if ( $line_2 =~ m/\t/ ) {
    print STDERR "Validation error: FastQ read file contains a tab: $input_file_2\n";
    die;
  }
  $line_count = ($line_count + 1) % 4;
}
$ifh1->close;
$ifh2->close;
exit;
