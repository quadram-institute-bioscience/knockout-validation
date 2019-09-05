#!/usr/bin/env perl
use 5.018;
use warnings;
use Getopt::Long;
use FASTX::Reader;

my $opt_fasta;
my $opt_cut;
my $opt_width = 350;
my $opt = GetOptions(
 'f|fasta=s'          => \$opt_fasta,
 'w|flanking-width=i' => \$opt_width,
 'c|cut'              => \$opt_cut,
);

my $program_usage=<<END;
   FASTA GREP

   USAGE:
   grepfa.pl [options] PATTERN1 PATTERN2...
 ---------------------------------------------------------------------------
   -f, --fasta  FILE              FASTA file to scan for patterns 
   -w, --flanking-width INT       Number of flanking bases to report [$opt_width]
   -c, --cut                      If two patterns are found in the same
                                  sequence (contig), report also the inner
                                  region
 ---------------------------------------------------------------------------
END

unless (-e "$opt_fasta") {
  say STDERR " * Missing fasta file (-f FILENAME)";
}
unless ($ARGV[0]) {
  say STDERR " * Specify at least one DNA pattern";
}

my $FASTA = FASTX::Reader->new({ filename => "$opt_fasta" });

say STDERR "# Scanning: $opt_fasta";

while (my $seq = $FASTA->getRead() ) {
 foreach my $pattern (@ARGV) {
   $pattern = uc($pattern);
   my $reverse = rc($pattern);
   my $regex = "($pattern|$reverse)";
   my @starts = ();
   while ($seq->{seq} =~/$regex/g) {
     my $strand = '+';
     $strand    = '-' if ($1 eq $reverse);
     my $left  = substr($seq->{seq}, $-[0] - $opt_width, $opt_width);
     my $right = substr($seq->{seq}, $+[0], $opt_width);
     push(@starts, $-[0]);
     say "# $pattern\t$strand\t$seq->{name}\t$-[0]-$+[0]";
     say ">Upstream_$1 ctg=$seq->{name}\n$left\n>Downstream_$1 ctg=$seq->{name}\n$right";
     if ($opt_cut and $#starts == 1) {
    	say ">", $seq->{name}, " upstream\n",
    		substr($seq->{seq}, $starts[0] - $opt_width, $opt_width);

	    say ">", $seq->{name}, " slice:@starts\n",
		    substr($seq->{seq}, $starts[0], $starts[1] - $starts[0]);

    	say ">", $seq->{name}, " downstream\n",
    		substr($seq->{seq}, $starts[1] , $opt_width);

     }
   }  
 }
}

sub rc  {
 my $s = reverse $_[0];
 $s =~tr/ACGT/TGCA/;
 return $s;
}
