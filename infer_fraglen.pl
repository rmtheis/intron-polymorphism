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
# Infers fragment length by building a distribution from alignments in a SAM file.
#

my $usage_msg = "Infers fragment length from SAM format alignments.\n"
              . "Usage: infer_fraglen.pl -i sam_file [-b bin_size (Default: 10)] [-m (Median)]\n";
die $usage_msg unless ( @ARGV );
my $bin_size = 10;
my ( $input_file, $median );
GetOptions(
  "b:i" => \$bin_size,
  "i=s" => \$input_file, # SAM input file
  "m" => \$median # Optional flag for reporting only median fragment length. Default: binned output
) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );

$input_file =~ s/^~/$ENV{HOME}/;
my $ifh = IO::File->new( $input_file, 'r' ) or die "Can't open $input_file: $!";
my %fragments = ();
my $count = 0;
while ( my $line1 = $ifh->getline ) {
  next if $line1 =~ m/^@/;
  my @fields = split(/\t/, $line1);
  my ($chr1, $frag1) = ($fields[2], $fields[8]);
  next if ( ($chr1 eq "*") || ($frag1 == 0) );

  # Get the second alignment in the pair
  my $line2 = $ifh->getline;
  @fields = split( /\t/, $line2 );
  my $chr2 = $fields[2];
  next if $chr1 ne $chr2;
  
  # Round to nearest bin
  my $frag = int( (abs($frag1) + ($bin_size / 2)) / $bin_size );
  $fragments{$frag}++;
  $count++;
  last if ($count == 100000);
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
