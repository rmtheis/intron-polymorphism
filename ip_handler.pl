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
use IntronPoly;
use Getopt::Long;
#use IPC::System::Simple qw(capture $EXITVAL);

##################################
# CONFIGURE USER RUNTIME OPTIONS #
##################################

# Default general parameters
my $scripts_path = "/home/theis/intron-polymorphism/"; # Location of intron-polymorphism scripts
my $input_path = "/home/theis/intron-polymorphism/reads/"; # Location of reads files
my $skipto = ""; # Skip to a step in the pipeline

# Parse command-line options
GetOptions(
  'i=s' => \$input_path,
  'k=s' => \$skipto,
  's=s' => \$scripts_path,
);

##############################
# JUMP TO THE REQUESTED STEP #
##############################

if ( $skipto eq "C" ) { print "Skipping to collecting\n"; goto COLLECT; }
if ( $skipto eq "F" ) { print "Skipping to filtering\n";  goto FILTER; }
if ( $skipto eq "A" ) { print "Skipping to assembly\n";   goto ASSEMBLE; }
if ( $skipto eq "L" ) { print "Skipping to alignment\n";  goto ALIGN; }
if ( $skipto eq "N" ) { print "Skippint to analysis\n";   goto ANALYZE; }

####################
# RUN THE PIPELINE #
####################

# Initialize the project
my $project = IntronPoly->new();
$project->set_workdir( $scripts_path );

MAP:
#$project->set_inputdir( $input_path );
#$project->run_mapping();

#COLLECT:


#FILTER:


#ASSEMBLE:


#ALIGN:


#ANALYZE:



