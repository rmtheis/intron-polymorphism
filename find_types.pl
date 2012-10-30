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
# Prints the distribution of alignment types in a SAM file.
#

my $usage_msg =
  "Prints the distribution of alignment types in a SAM file.\n"
  . "Usage: find_types.pl -i sam_file\n";
die $usage_msg unless ( @ARGV );
my $input_file;
GetOptions( "i=s" =>\$input_file ) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );
$input_file =~ s/^~/$ENV{HOME}/;

my $non = 0;
my $full = 0;
my $half = 0;
my $disc = 0;
my $total = 0;
my ( %flags_disc, %flags_full, %flags_non, %flags_half );
my $ifh = new IO::File( $input_file, 'r' ) or die "Can't open $input_file: $!";
while ( my $line1 = $ifh->getline ) {
  next if $line1 =~ m/^@/;
  my @f = split(/\t/, $line1);
  my ($id1, $flags1) = ($f[0], $f[1]);
  
  my $line2 = $ifh->getline;
  my @g = split(/\t/, $line2);
  my ($id2, $flags2) = ($g[0], $g[1]); 
  $total++;
  
  # Check for discordant pairs first, because they may also match other categories
  if (IntronPoly::_isDiscordantPair($flags1, $flags2)) { $disc++; $flags_disc{"$flags1,$flags2"}++; } 
  elsif (IntronPoly::_isHalfMappingPair($flags1, $flags2)) { $half++; $flags_half{"$flags1,$flags2"}++; }
  elsif (IntronPoly::_isInMappingPair($flags1)) { $full++; $flags_full{"$flags1,$flags2"}++; }
  elsif (IntronPoly::_isInNonMappingPair($flags1)) { $non++; $flags_non{"$flags1,$flags2"}++; }
  else { print "not categorized: $flags1, $flags2\n"; }
}
$ifh->close;

my $percent_non = sprintf('%.3f', $non / $total * 100);
my $percent_full = sprintf('%.3f', $full / $total * 100);
my $percent_half = sprintf('%.3f', $half / $total * 100);
my $percent_disc = sprintf('%.3f', $disc / $total * 100);
print "Total reads: $total\n";
print "Full-mapping: $full ($percent_full%)\n";
foreach my $key (sort { $flags_full{$b} <=> $flags_full{$a} } keys %flags_full) {
  print "  $key: " . $flags_full{$key} . "\n";
}
print "Non-mapping: $non ($percent_non%)\n";
foreach my $key (sort { $flags_non{$b} <=> $flags_non{$a} } keys %flags_non) {
  print "  $key: " . $flags_non{$key} . "\n";
}
print "Half-mapping: $half ($percent_half%)\n";
foreach my $key (sort { $flags_half{$b} <=> $flags_half{$a} } keys %flags_half) {
  print "  $key: " . $flags_half{$key} . "\n";
}
print "Discordant: $disc ($percent_disc%)\n";
foreach my $key (sort { $flags_disc{$b} <=> $flags_disc{$a} } keys %flags_disc) {
  print "  $key: " . $flags_disc{$key} . "\n";
}

exit;


