#!/usr/bin/perl -w

use strict;
use IntronPoly;

my $project = IntronPoly->new();
$project->set_work_dir( "/home/theis/ccmp-rerun", "/home/theis/intron-polymorphism" );
$project->build_db(
  "/home/theis/testdata/MicromonasCCMP1545.fasta",
  "/home/theis/testdata/CCMP490_7.1.trimmed",
  "/home/theis/testdata/CCMP490_7.2.trimmed",
);

$project->set_fragment_length( 280 );
#$project->assemble_groups( 250, 3, 70, "/home/theis/intron-polymorphism/ip_out/run-20121102_16_09_18/CCMP490_7.1_filtered2_sorted.sam" );
$project->align_contigs_clustal( "/home/theis/ccmp-reassemble/run-20121107_22_51_55/CCMP490_7.1_contigs.fa" );

