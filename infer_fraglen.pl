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
# Infers fragment length by building a distribution from uniquely aligning mates.
#

unless ( @ARGV ) {
  print "Usage: $0 -i input_file {-m}\n";
  die;
}
my $infile;
my $median;
GetOptions(
  'i:s' => \$infile, # SAM format input file to use for inferring fragment length
  'm' => \$median # Optional flag for reporting only median fragment length. Default: binned output
) || die "$0: Bad option";

my $ifh  = IO::File->new( $infile, 'r' ) or die "Can't open $infile: $!";
my $bin_size = 10;
my %fragments = ();
while ( my $line1 = $ifh->getline ) {
  next if $line1 =~ m/^@/;
  my @fields = split(/\t/, $line1);
  my ($name1, $chr1, $frag1) = ($fields[0], $fields[2], $fields[8]);
  next if ( ($chr1 eq "*") || ($frag1 == 0) );

  # Get the second alignment in the pair
  my $line2 = $ifh->getline;
  @fields = split( /\t/, $line2 );
  my $chr2 = $fields[2];
  next if $chr1 ne $chr2;
  
  # Round to nearest bin
  my $frag = int( (abs($frag1) + ($bin_size / 2)) / $bin_size );
  $fragments{$frag}++;
}
$ifh->close;

if ($median) {
  # Print bin with largest number of fragments
  my $b = 0;
  my $max_fragments = 0;
  for my $key (keys %fragments) {
    if ($fragments{$key} > $max_fragments) {
      $max_fragments = $fragments{$key};
      $b = $key;
    }
  }
  my $bin = $b * $bin_size;
  print "$bin\n";
} else {
  # Print sorted bins and frequencies
  for my $key (sort {$a <=> $b} keys %fragments) {
    my $bin = $key * $bin_size;
    print "$bin, $fragments{$key}\n";
  }  
}
exit;
