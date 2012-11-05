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

my $usage_msg =
  "Converts a SAM format file to paired-end Fastq. SAM flags determine whether mates are mate 1 or mate 2.\n"
  . "Non-mapping pairs and secondary alignments are discarded.\n"
  . "Usage: sam2pefastq.pl -i sam_file -1 outfile1.fq -2 outfile2.fq\n";
die $usage_msg unless ( @ARGV );
my ($input_file, $output_file1, $output_file2);
GetOptions( "i=s" =>\$input_file,
            "1=s" =>\$output_file1,
	    "2=s" =>\$output_file2,
	  ) || die "$0: Bad option";
die $usage_msg unless ( defined $input_file && defined $output_file1 && defined $output_file2 );

$input_file =~ s/^~/$ENV{HOME}/;
my $ifh = IO::File->new( $input_file, 'r' ) or die "Can't open $input_file: $!";
my $ofh1 = IO::File->new( $output_file1, 'w' ) or die "Can't create $output_file1: $!";
my $ofh2 = IO::File->new( $output_file2, 'w' ) or die "Can't create $output_file2: $!";
while ( my $line = $ifh->getline ) {
  next if $line =~ m/^@/;
  my @f = split(/\t/, $line);
  my ($id, $flags, $chr, $offset, $seq, $qual) = ($f[0], $f[1], $f[2], $f[7], $f[9], $f[10]);
  next if (IntronPoly::_isInNonMappingPair($flags) != 0);
  next if (IntronPoly::_isSecondaryAlignment($flags) != 0);

  # Write Fastq output
  if (IntronPoly::_isMateOneInPair($flags) != 0) {
    print $ofh1 "@" . "${id}\n";
    print $ofh1 "$seq\n";
    print $ofh1 "+\n";
    print $ofh1 "$qual\n";
  } else {
    print $ofh2 "@" . "${id}\n";
    print $ofh2 "$seq\n";
    print $ofh2 "+\n";
    print $ofh2 "$qual\n";
  }
}
$ifh->close;
exit;


