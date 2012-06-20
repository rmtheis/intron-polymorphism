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
# Bins the half-mapping reads in a SAM file.
#

unless ( @ARGV ) {
  print "Usage: $0 -i input_file\n";
  exit;
}
my $aligned;
my $infile;
GetOptions(
  #'a' => \$aligned, # Optional flag for also reporting binned read depth for aligning reads.
  'i:s' => \$infile # SAM input file
) || die "$0: Bad option";

my $ifh = IO::File->new( $infile, 'r' ) or die "Can't open $infile: $!";
my $bin_size = 100;
my %fragments = ();
while ( my $line1 = $ifh->getline ) {
  next if $line1 =~ m/^@/;
  my @fields = split(/\t/, $line1);
  my ($flags1, $chr1) = ($fields[1], $fields[2]);
  next if (IntronPoly::_isInMappingPair($flags1) != 0);
  next if (IntronPoly::_isSecondaryAlignment($flags1) != 0);

  # Get the second alignment in the pair
  my $line2 = $ifh->getline;
  @fields = split( /\t/, $line2 );
  my ($flags2, $chr2, $pos) = (int($fields[1]), $fields[2], $fields[7]);
  next if $chr1 ne $chr2;
  next if (IntronPoly::_isSecondaryAlignment($flags2) != 0);
  next if (IntronPoly::_isFarMappingPair($flags1, $flags2) != 0);
  if (IntronPoly::_isHalfMappingMate($flags2) != 0) {
    # Round to nearest bin
    my $frag = int( ($pos + ($bin_size / 2)) / $bin_size );
    $fragments{$frag}++;
  }
}
$ifh->close;

# Print sorted bins and frequencies
for my $key (sort {$a <=> $b} keys %fragments) {
  my $bin = $key * $bin_size;
  print "$bin, $fragments{$key}\n";
}  
exit;
