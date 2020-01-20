#!/usr/bin/env perl
use 5.012;
use warnings;
use FindBin qw($RealBin);
use FASTX::Reader;
use Data::Dumper;
use List::Util qw/shuffle/;
my $ref_file = "$RealBin/ref/example.fa";

my $R = FASTX::Reader->new({ filename => "$ref_file" });
my $first = $R->getRead();
my $ref = $first->{seq};

my $total = 20;
my $size = length($ref);

my $seqs_done = 0;
for (my $i = 1; $i <= $total; $i++) {
	my $pos = int(rand($size - 50));
	for (my $oligo_len = 16; $oligo_len <= 20; $oligo_len += 2) {
		for my $strand ('for', 'rev') {
			for my $mismatches (0, 1, 3, 5) {
 				$seqs_done++;
 				my $s = substr($ref, $pos, $oligo_len);
 				my $z = $s;
 				$z = rc($z) if ($strand eq 'rev');
 				
 				$z = mismatches($z, $mismatches);
 				die "MM Seq #$seqs_done oligo_${strand}_${oligo_len}bp_pos${pos}_err$mismatches:\n$s!$z\n"
 					unless $z;
 			
 				say ">oligo${seqs_done}_${strand}_${oligo_len}bp_pos${pos}_err$mismatches";
 				say "$z";
 			}
		}
	}
}


sub mismatches {
	my ($seq, $err) = @_;
	return $seq unless $err;
	my @pos = shuffle 0..length($seq)-1;
	die "Unsupported errors: $err > 6" if ($err > 6);
	for (my $e = 0; $e < $err; $e++) {
		substr($seq, $pos[$e]) = newletter(substr($seq, $pos[$e], 1)) . substr($seq, $pos[$e] + 1);
	}

	return $seq;
}


sub newletter {
	my $letter = shift @_;
	my %seq = (
		'A' => ['C', 'G', 'T'],
		'C' => ['A', 'G', 'T'],
		'G' => ['C', 'A', 'T'],
		'T' => ['C', 'G', 'A'],
	);
	return $seq{$letter}[rand(3)];
}

sub rc {
	my $seq = shift @_;
	$seq = reverse($seq);
	$seq =~tr/ACGTacgt/TGCAtgca/;
	return $seq;
}