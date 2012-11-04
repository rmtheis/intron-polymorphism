#!/usr/bin/perl

#  TAKES A PAIR OF MATE-PAIR FILES.
#  FOR EACH READ, TRIMS THEM BY THE STRATEGY LAID OUT IN KELLEY ET AL
#  (PMID: 21114842), EQUATION 1.

($#ARGV == 3) || die "usage: trim2.pl mate1.fq mate2.fq t-score-parameter minimum_output_read_length\n\n";
$t = $ARGV[2];
$min_read_length = $ARGV[3];


#  Input, output files
open (IN1, $ARGV[0]);
open (IN2, $ARGV[1]);
open (OUT1, ">$ARGV[0].trimmed");
open (OUT2, ">$ARGV[1].trimmed");


#  LOOPS THROUGH, TRIMS EACH READ PAIR, ASKS
#  IF EITHER READ IS NOW SHORTER THAN $min_read_length.
#  IF NOT, PRINTS THEM BOTH OUT.

while (!eof(IN1)) {
    $discard = 0;
    (@{$data[0]}) = get_fastq(*IN1);
    (@{$data[1]}) = get_fastq(*IN2);
    #  This is the meet of the algorithm
    for $f (0,1) {
        $score = $max = $trim[$f] = 0;
        #  Quality scores in reverse order
        @b = (reverse $data[$f][3] =~ /./g);
        #  Go through the quality scores one by one
        for $i (0..$#b) {
            #  Increment the score by $t - qual_value
            $score += $t - ((ord $b[$i])-66);

            if ($score > $max) {
                $max = $score;
                $trim[$f] = $i+1;
            }
            #  If the score is way less than the already-established maximum,
            #  don't waste your time.
            elsif ($score < $max - 500) {last};
        }
        #  This variable keeps track of whether either read will be too short
        $discard += ($trim[$f] > ((length $data[0][1])-$min_read_length));
    }
    print "$discard, $trim[0], $trim[1]\n";
    
    #  Trim reads, print 'em
    unless ($discard) {
    $data[0][1] =~ s/.{$trim[$f]}(?=\n)$//;
    $data[0][3] =~ s/.{$trim[$f]}(?=\n)$//;
    $data[1][1] =~ s/.{$trim[$f]}(?=\n)$//;
    $data[1][3] =~ s/.{$trim[$f]}(?=\n)$//;
    print OUT1 join ("", @{$data[0]});
    print OUT2 join ("", @{$data[1]});
    }
}


sub get_fastq () {
    my $handle = $_[0];
    my @l;
    $l[0] = <$handle>;
    $l[1] = <$handle>;
    $l[2] = <$handle>;
    $l[3] = <$handle>;
    @l;
}

