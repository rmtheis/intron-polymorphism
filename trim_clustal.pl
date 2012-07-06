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
# Trims unaligned portions of a clustal multiple sequence alignment output file.
#
# Output is printed to standard out.
#

my $usage_msg =
  "Views a clustal multiple sequence alignment output file, ignoring unaligned regions.\n"
  . "Usage: trim_clustal.pl -i clustal_file.aln {-m margin (default:50)}\n";
die $usage_msg unless ( @ARGV );
my $input_file;
my $margin_size = 50;
GetOptions(
   "i:s" =>\$input_file,
   "m:i" =>\$margin_size
) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );

my $ifh = new IO::File( $input_file, 'r' ) or die "Can't open $input_file: $!";
my $first_line = $ifh->getline;
print STDERR "Warning: first line does not begin with 'CLUSTAL'\n" if $first_line !~ m/^CLUSTAL/;
$ifh->getline;
$ifh->getline;

my $pos = 0;
my $cnt = 0;
my @buffer = ();
my $buffer_size = 0;
my $print_trailing_positions = 0;
while ( my $line = $ifh->getline ) { 
  # Handle a group when we see a blank line or '*' characters indicating aligned positions
  if ($line =~ /^\s*$/ || $line =~ m/\*/) {
    push(@buffer, $line);
    push(@buffer, $ifh->getline);
    $pos += 50;
    $buffer_size += 50;
    
    # Check for an alignment (two or more sequences are present)
    if ($cnt >= 2) {
      print "$pos\n";
      while ( scalar(@buffer) > 0 ) {
        print shift(@buffer);
      }
      $buffer_size = 0;
      $print_trailing_positions = $margin_size;
    }

    # Check if this group follows an alignment
    elsif ($print_trailing_positions > 0) {
      print shift(@buffer) while ( scalar(@buffer) > 0 );
      $print_trailing_positions -= 50;
      $buffer_size -= 50;
      $cnt = 0;
      next;
    }    
    $cnt = 0;
    
    # Check if we need to flush a group from the buffer
    if ($buffer_size > $margin_size) {
      my $buffer_line = "";
      while ( scalar(@buffer) > 0 ) {
        $buffer_line = shift(@buffer);
        if ($buffer_line =~ /^\s*$/) {
          shift(@buffer);
          last;
        }
      }
    }
  } else {
    # Add the current non-blank line to the buffer
    push(@buffer, $line);
    my @fields = split(/\s+/, $line);
    $cnt++ if ($fields[1] !~ m/^-+$/);
  }
}
$ifh->close;
exit;
