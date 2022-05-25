#
# Compute sims with diamond
#

use strict;
use Getopt::Long::Descriptive;
use File::Path qw(make_path);
use IPC::Run 'run';

my($opt, $usage) = describe_options("%c %o sims-dir",
				    ["threads=i", "Number of threads to use", { default => 4 }],
				    ["keep=i", "Keep this many matches", { default => 200 }],
				    ["diamond-sensitivity=s", "Diamond sensitivity flag to use", { default => "sensitive" }],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 1;

my $sims_dir = shift;

my($tag) = $sims_dir =~ /(\d+)$/;
make_path("$sims_dir/sims.split");

my @cmd = ("diamond", "blastp", "-v",
	      "-b5",
	      "-c1",
	      "-e", "1e-4",
	      "--masking", 0,
	      "--compress", 1,
	      "--threads", $opt->threads,
	      "-k", $opt->keep,
	      "--outfmt", qw(6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen),
	      "-d", "$sims_dir/nr",
	      "-q", "$sims_dir/seqs.added",
	      "-o", "$sims_dir/diamond.out");
my $rc = system(@cmd);

$rc == 0 or die "Error running @cmd\n";

my $ok = run(["zcat", "$sims_dir/diamond.out.gz"],
	  "|", 
	  ["reduce_sims", "$sims_dir/nr-len.btree", $opt->keep],
	  "|",
	  ["split_sims", "$sims_dir/sims.split", "sims.$tag"]);
$ok or die "Reduce/split pipeline failed\n";
