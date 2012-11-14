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
# Given sequence name and start/stop positions, return the resulting subsequence from a Fasta file.
#
# Positions are indexed starting with the first character in the given sequence.
# Output goes to standard out.
#

my $usage_msg = "From a Fasta file, returns the subsequence delimited by <start> and <stop> index.\n"
              . "Usage: trim_fasta.pl -s scaffold_1 -i myfile.fa --start <int> --stop <int>\n";
die $usage_msg unless ( @ARGV );
my ( $in, $chr, $start, $stop );
GetOptions(
  "i=s" => \$in,
  "s=s" => \$chr,
  "l|start=s" => \$start,
  "r|stop=s" => \$stop,
  ) || die "$0: Bad option";
die $usage_msg unless ( defined $in && defined $chr && defined $start && defined $stop );
die "Bad start/stop position." if ( $stop <= $start || $start < 0 || $stop < 0 );
$in =~ s/^~/$ENV{HOME}/;
my $ifh = new IO::File( $in, 'r' ) or die "Can't open $in: $!";
my $index = 0;
my $seq = "";
# Find the header line for the requested section
while ( my $line = $ifh->getline ) {
  if ($line =~ m/^>/) {
    $line =~ s/^\>|\s+$//g;
    last if ($line eq $chr);
  }
}

# Find the subsequence from this section
while ( my $line = $ifh->getline ) {
  last if ($line =~ m/^>/); # Stop if we've reached the end of this sequence
  $line =~ s/\s+$//g; # Remove trailing whitespace
  $index += length($line);
  next if ($index < $start);
  
  if ($index >= $stop && abs($stop - $start) <= length($line)) {
    # Get the entire sequence we want, it's all on one line
    $seq = substr($line, length($line) - ($index - $start), $stop - $start);
    last;
  }
  
  if ($index >= $stop) {
    # Get the last portion of the sequence we want
    $seq = $seq . substr($line, 0, length($line) - ($index - $stop));      
    last;
  }

  if (abs($index - $start) < length($line)) {
    # Get the beginning of the sequence we want
    $seq = $seq . substr($line, -($index - $start))
  } else {
    # Get intermediate lines in the sequence we want
    $seq .= $line;
  }
  
}
$ifh->close;
die if (length($seq) == 0);
print ">#$chr|$start|$stop\n";
my $length = 70; # Number of characters per line for sequences
for ( my $pos = 0; $pos < length($seq); $pos += $length) {
  print substr($seq, $pos, $length), "\n";
}
exit;
