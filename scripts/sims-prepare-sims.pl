
#
# Prepare for sims computation.
#
# Compute the differences in the new and old NR, and pull the fasta data.
#

use strict;
use FIG;
use FIG_Config;
use File::Basename;
use IPC::Run 'run';

@ARGV == 1 or die "Usage: $0 job-dir\n";

my $jobdir = shift;

-d $jobdir or die "$0: job dir $jobdir does not exist\n";


#
# First compute the changes in ids.
#

my @cmd = ("$FIG_Config::bin/compute_changed_ids_for_nrs",
	   "$jobdir/prev_nr",
	   "$jobdir/prev_syn",
	   "$jobdir/nr",
	   "$jobdir/ids.added",
	   "$jobdir/ids.changed",
	   "$jobdir/ids.deleted");
my $ok = run(\@cmd);
$ok or die "Error running @cmd\n";

my $cmd = "$FIG_Config::bin/pull_fasta_entries $jobdir/nr < $jobdir/ids.added > $jobdir/seqs.added";

$ok = run(["pull_fasta_entries", "$jobdir/nr"],
	  "<", "$jobdir/ids.added",
	  ">", "$jobdir/seqs.added");
$ok or die "pull_fasta failed\n";

$ok = run(["diamond", "makedb", "--in", "$jobdir/nr", "--db", "$jobdir/nr"]);
$ok or die "Error running diamond makedb\n";
    
