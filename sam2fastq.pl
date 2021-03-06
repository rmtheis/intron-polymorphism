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

use IntronPoly;
use lib '..';
use strict;
use Getopt::Long;
use IO::File;

#
# Converts a SAM format file to Fastq format.
#
# Output is printed to standard out.
#

my $usage_msg =
  "Converts a SAM format file to Fastq format.\n"
  . "Usage: sam2fastq.pl -i sam_file > output.fq\n";
die $usage_msg unless ( @ARGV );
my $input_file;
GetOptions( "i=s" =>\$input_file ) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );

$input_file =~ s/^~/$ENV{HOME}/;
my $ifh = new IO::File( $input_file, 'r' ) or die "Can't open $input_file: $!";
while ( my $line = $ifh->getline ) {
  next if $line =~ m/^@/;
  my @f = split(/\t/, $line);
  my ($id, $flags, $chr, $offset, $seq, $qual) = ($f[0], $f[1], $f[2], $f[7], $f[9], $f[10]);
  $qual =~ s/\s+$//;

  # Write Fastq output
  print "@" . "${id}_${flags}\n";
  print "$seq\n";
  print "+\n";
  print "$qual\n";
}
$ifh->close;
exit;


