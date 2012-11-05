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
# Bins the half-mapping reads in a SAM file by position in the reference genome.
#

my $usage_msg = "Bins distribution of half-mapping alignments by position specified in a SAM file.\n"
              . "Usage: find_distrib.pl -i sam_file [-b bin_size (default: 100)] [--bowtie1]\n";
die $usage_msg unless ( @ARGV );
my ( $aligned, $input_file, $bowtie1 );
my $bin_size = 100;
GetOptions(
  "b:i" => \$bin_size,
  "i=s" => \$input_file, # SAM input file
  "bowtie1" => \$bowtie1 # Expect Bowtie 1 flags for unpaired reads
) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );

$input_file =~ s/^~/$ENV{HOME}/;
my $ifh = IO::File->new( $input_file, 'r' ) or die "Can't open $input_file: $!";
my %fragments = ();

if ($bowtie1) {
  # Bowtie 1
  # Use unpaired read flags from Bowtie 1 output
  # For a half-mapping pair, one of the two mates will have 4 as its flag,
  # and the other mate will have 0 or 16 as its flag
  while ( my $line1 = $ifh->getline ) {
    next if $line1 =~ m/^@/;
    my @fields = split(/\t/, $line1);
    my ($flags1, $chr1, $pos1) = (int($fields[1]), $fields[2], $fields[3]);
    
    # Get the second alignment in the pair
    my $line2 = $ifh->getline;
    @fields = split( /\t/, $line2 );
    my ($flags2, $chr2, $pos2) = (int($fields[1]), $fields[2], $fields[3]);
  
    if (($flags1 == 4) ^ ($flags2 == 4)) {
      my $pos;
      if ($pos1 != 0) {
        $pos = $pos1;
      } else {
        $pos = $pos2;
      }
      # Round to nearest bin
      my $frag = int( ($pos + ($bin_size / 2)) / $bin_size );
      $fragments{$frag}++;
    }
  }
} else {
  # Bowtie 2
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
}
$ifh->close;

# Print sorted bins and frequencies
for my $key (sort {$a <=> $b} keys %fragments) {
  my $bin = $key * $bin_size;
  print "$bin, $fragments{$key}\n";
}  
exit;
