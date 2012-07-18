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
# Trims bases from Fastq paired-end reads.
#

my $usage_msg = "Trims bases from Fastq paired-end reads.\n"
              . "Usage: trim_fastq.pl -1 fastq1 -2 fastq2 -o1 outfile1 -o2 outfile2 -5 <int> -3 <int>\n";
die $usage_msg unless ( @ARGV );
my ( $mate1s, $mate2s, $out1, $out2, $trim3, $trim5 );
GetOptions(
  "1=s" => \$mate1s,
  "2=s" => \$mate2s,
  "o1=s" => \$out1,
  "o2=s" => \$out2,
  "3|trim3:i" => \$trim3,
  "5|trim5:i" => \$trim5,
  ) || die "$0: Bad option";
die $usage_msg unless ( defined $mate1s && defined $mate2s && defined $out1 && defined $out2);

my $ifh1 = new IO::File( $mate1s, 'r' ) or die "Can't open $mate1s: $!";
my $ifh2 = new IO::File( $mate2s, 'r' ) or die "Can't open $mate2s: $!";
my $ofh1 = new IO::File( $out1, 'w' ) or die "Can't open $out1: $!";
my $ofh2 = new IO::File( $out2, 'w' ) or die "Can't open $out2: $!";
my $line_count = 0;
while ( my $line_1 = $ifh1->getline ) {
  my $line_2 = $ifh2->getline;
  if ( $line_count == 1 ) {
    my $trimmed1 = $line_1;
    my $trimmed2 = $line_2;
    $trimmed1 =~ s/\s+$//;
    $trimmed2 =~ s/\s+$//;
    if (defined $trim5) {
      $trimmed1 = substr($trimmed1, $trim5);
      $trimmed2 = substr($trimmed2, 0, length($trimmed2) - $trim5);
    }
    if (defined $trim3) {
      $trimmed1 = substr($trimmed1, 0, length($trimmed2) - $trim3);
      $trimmed2 = substr($trimmed2, $trim3);
    }
    print $ofh1 "$trimmed1\n";
    print $ofh2 "$trimmed2\n";
  } else {
    print $ofh1 $line_1;
    print $ofh2 $line_2;
  }
  $line_count = ($line_count + 1) % 4;
}
$ifh1->close;
$ifh2->close;
$ofh1->close;
$ofh2->close;
exit;
