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
# Prints a Clustal multiple sequence alignment output file, ignoring unaligned regions.
#
# Output is printed to standard out.
#

my $usage_msg =
  "Prints a Clustal multiple sequence alignment output file, ignoring unaligned regions.\n"
  . "Usage: trim_clustal.pl -i clustal_file.aln {--margin [int, default=100]} {--no-num}\n";
die $usage_msg unless ( @ARGV );
my $input_file;
my $margin_size = 100;
my $suppress_line_numbers = 0;
GetOptions(
   "i|input=s" => \$input_file,
   "m|margin:i" => \$margin_size, # Length of sequence to show before/after aligned regions
   "n|no-num" => \$suppress_line_numbers
) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );

$input_file =~ s/^~/$ENV{HOME}/;
my $ifh = new IO::File( $input_file, 'r' ) or die "Can't open $input_file: $!";
my $first_line = $ifh->getline;
print STDERR "Warning: first line does not begin with 'CLUSTAL'\n" if $first_line !~ m/^CLUSTAL/;
$ifh->getline;
$ifh->getline;
my @buffer = ();
my $buffer_size = 0;
my $cnt = 0;
my $pos = 0;
my $trailing_positions = 0;
while ( my $line = $ifh->getline ) { 
  # Handle a group when we see a blank line or '*' characters indicating aligned positions
  if ($line =~ /^\s*$/ || $line =~ m/\*/) {
    $pos += 50;
    $buffer_size += 50;
    push(@buffer, $line);
    push(@buffer, $ifh->getline);  
    
    # Check for an alignment (two or more sequences are present)
    if ($cnt >= 2) {
      print shift(@buffer) while ( scalar(@buffer) > 0 );
      $buffer_size = 0;
      $trailing_positions = $margin_size;
    } elsif ($trailing_positions > 0) { # Check if this group follows an alignment
      print shift(@buffer) while ( scalar(@buffer) > 0 );
      $trailing_positions -= 50;
      $buffer_size -= 50;
      $cnt = 0;
      push(@buffer, "$pos\n") unless $suppress_line_numbers;
      next;
    }    
    $cnt = 0;
    
    # Check if we need to flush a group from the buffer
    if ($buffer_size > $margin_size) {
      while ( scalar(@buffer) > 0 ) {
        if (shift(@buffer) =~ /^\s*$/) { # A blank line ends the group
          shift(@buffer);
          last;
        }
      }
    }
    push(@buffer, "$pos\n") unless $suppress_line_numbers; 
  } else {
    push(@buffer, $line); # Add the current non-blank line to the buffer
    my @fields = split(/\s+/, $line);
    $cnt++ if ($fields[1] !~ m/^-+$/);
  }
}
$ifh->close;
exit;
