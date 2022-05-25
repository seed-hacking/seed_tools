
#
# Build an NR.
#

use strict;
use FIG;
use FIG_Config;
use File::Basename;
use IPC::Run 'run';

@ARGV == 1 or die "Usage: $0 job-dir\n";

my $jobdir = shift;

-d $jobdir or die "$0: job dir $jobdir does not exist\n";

my @cmd = ("build_nr_md5",
	   "$jobdir/nr.sources",
	   "$jobdir/nr",
	   "$jobdir/peg.synonyms",
	   "$jobdir/nr-len.btree",
	   "$jobdir/nr-fig-ids");


my $ok = run(\@cmd);
$ok or die "Error building NR: build_nr @cmd\n";

#
# And the pegsyn, for reduce_sims use.
#
$ok = run(["make_reduce_sims_index",
	   "$jobdir/peg.synonyms",
	   "$jobdir/peg.synonyms.reduce_sims_index.btree"]);
$ok or die "Error running make_reduce_sims_index\n";

