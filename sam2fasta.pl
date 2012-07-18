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
# Prints a SAM file to standard out in Fasta format.
#

my $usage_msg = "Prints a SAM file to a Fasta file.\n"
              . "Usage: sam2fasta.pl -i sam_file > output.fa\n";
die $usage_msg unless ( @ARGV );
my $input_file;
GetOptions( "i=s" =>\$input_file ) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file );
$input_file =~ s/^~/$ENV{HOME}/;

my $ifh = new IO::File( $input_file, 'r' ) or die "Can't open $input_file: $!";
while ( my $line = $ifh->getline ) {
  next if $line =~ m/^@/;
  my @fields = split(/\t/, $line);
  my ( $id, $seq ) = ( $fields[0], $fields[9] );
  print ">$id\n";
  print "$seq\n";
}
$ifh->close;
exit;

