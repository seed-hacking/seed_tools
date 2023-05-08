
#
# Prepare for sims computation.
#
# Compute the differences in the new and old NR, and pull the fasta data.
#
# Create the diamond index for the NR.
#

use strict;
use FIG;
use FIG_Config;
use File::Basename;
use IPC::Run 'run';

use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o job-dir",
	["max-chunk=i", "Maximum chunk size for sims compute", { default => 500_000 }],
	["help|h" => "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 1;

my $jobdir = shift;

-d $jobdir or die "$0: job dir $jobdir does not exist\n";

my @cmd = ("sims-determine-seqs-to-update",
	   "--max-chunk", $opt->max_chunk,
	   "$jobdir/prev_nr",
	   "$jobdir/nr",
	   "$jobdir/seqs");
my $ok = run(\@cmd);
$ok or die "Error running @cmd\n";

$ok = run(["diamond", "makedb", "--in", "$jobdir/nr", "--db", "$jobdir/nr"]);
$ok or die "Error running diamond makedb\n";
    
