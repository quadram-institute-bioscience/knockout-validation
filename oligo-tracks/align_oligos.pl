#!/usr/bin/perl env
use 5.012;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

my ($opt_ref, $opt_oligos, $opt_output);

my $_opt = GetOptions(
	'r|ref=s'      => \$opt_ref,
	'i|input=s'    => \$opt_oligos,
	'o|output=s'   => \$opt_output,
);

if (!defined $opt_ref or !defined $opt_oligos) {
	usage();

}

die "Unable to find reference <$opt_ref>\n" unless (-e "$opt_ref");
die "Unable to find reference <$opt_oligos>\n" unless (-e "$opt_oligos");

say STDERR<<END;
 * Reference    $opt_ref
 * Oligos       $opt_oligos
 * Output       $opt_output
END

if ( ! -e "${opt_ref}.bwt" ) {
	say STDERR " * Indexing reference: $opt_ref";
	run( qq(bwa index "$opt_ref" ) );
} else {
	say STDERR " * Index found: $opt_ref";
}

run('Align', qq(bwa mem -k 8 -r 1 -T 10 "$opt_ref" "$opt_oligos" 2>/dev/null | samtools view -bS | samtools sort -o "${opt_output}.bam" -));
run('Index', qq(samtools index "${opt_output}.bam"));
run('Coverage', qq(covtobed -m 1 "${opt_output}.bam" > "${opt_output}.bed"));

sub run {
	my ($desc,$command) = @_;
	die "Missing command to run for $desc" unless defined $command;
	say STDERR " * ", GREEN, BOLD, $desc, RESET, ": $command";
	system($command);
	die "ERROR running <$command>\n" if ($?);
}

sub usage {
	say<<END;

 align_oligos.pl -r REFERENCE -i OLIGO_FASTA -o OUTPUT_BASENAME
END
	exit;
}