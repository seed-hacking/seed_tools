#
# Create a SEED-import job.
#
# The job initially has a listing of SEED fasta sources, NR fasta sources, and
# RAST fasta sources. The new NR will not have been built - it will be the first stage
# in the pipeline.
#
# For now we need to pass in the path to the NR and peg.synonyms files we are building from.
#
# We create the following files:
#
#    nr.dirs
#	Tab-delimited data of db-name, source path, size of fasta file
#
#    nr.sources
#	Listing of all fasta source files from which the nr is to be built
#
# We hardcode in the script, for now, the source locations of things. This is an ANL internal
# application at this point.
#
# 
#

use strict;
use Data::Dumper;
use File::Basename;
use DirHandle;
use NRTools;
use File::Copy;
use FIG_Config;
use Getopt::Long::Descriptive;

# my $usage = "create_import_job [-no-rast-jobs] [-lustre] [-extra-jobs jobfile] [-no-old-NR] [-new-nr-data dir] [-import-biodb] [-use-last-NR] [-from-job jobnum] [prev-nr prev-syn prev-sims]";

my($opt, $usage) = describe_options("%c %o sims-dir [prev-nr prev-syn prev-sims]",
				    ["from-job=s", "Initialize from prior job number"],
				    ["help|h", "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if !(@ARGV == 1 || @ARGV == 4);

my $sims_dir = shift;
-d $sims_dir or die "Sims directory $sims_dir does not exist\n";

my $prev_nr_src;
my $prev_syn_src;
my $prev_sim_dir;

if ($opt->from_job)
{
    @ARGV == 0 or die $usage->text;

    my $dir = $sims_dir . "/" . $opt->from_job;
    -d $dir or die "from-job directory $dir does not exist\n";
    $prev_nr_src = "$dir/nr";
    $prev_syn_src = "$dir/peg.synonyms";
    $prev_sim_dir = sprintf("$dir/Sims.%03d", $opt->from_job);
}
else
{
    @ARGV == 3 or die $usage->text;

    $prev_nr_src = shift;
    $prev_syn_src = shift;
    $prev_sim_dir = shift;
}

#
# Validate
#
if (open(F, "<$prev_nr_src"))
{
    $_ = <F>;
    if ($_ && ! /^>/)
    {
	die "$prev_nr_src does not look like a fasta file\n";
    }
    close(F);
}
else
{
    die "Cannot open previous NR file $prev_nr_src: $!\n";
}

if (open(F, "<$prev_syn_src"))
{
    $_ = <F>;
    if ($_ && !/^(gnl\|md5\|[a-f0-9]+|xxx\d+),\d+\t/)
    {
	die "$prev_syn_src does not look like a peg.synonyms file\n";
    }
    close(F);
}
else
{
    die "Cannot open previous synonyms file $prev_syn_src: $!\n";
}

my @sfiles = <$prev_sim_dir/sims*>;
if (not(-d $prev_sim_dir and @sfiles > 0))
{
    warn "previous sim dir $prev_sim_dir does not appear to contain sims\n";
}


print "Creating import job\n";
print "\tprev_nr=$prev_nr_src\n";
print "\tprev_syn=$prev_syn_src\n";
print "\tprev_sim=$prev_sim_dir\n";

#
# Create our jobdir.
#

my ($jobnum, $err) = create_new_job($sims_dir);

if (!$jobnum)
{
    die "Create failed with error: $err\n";
}

my $jobdir = "$sims_dir/$jobnum";

#
# Symlink to prev_nr and prev_syn in the job directory.
#

my $prev_nr = "$jobdir/prev_nr";
my $prev_syn = "$jobdir/prev_syn";
my $prev_sims = "$jobdir/prev_sims";

unlink($prev_nr, $prev_syn, $prev_sims);

symlink($prev_nr_src, $prev_nr) or die "symlimk $prev_nr_src $prev_nr failed: $!";
symlink($prev_syn_src, $prev_syn) or die "symlimk $prev_syn_src $prev_syn failed: $!";
symlink($prev_sim_dir, $prev_sims) or die "symlimk $prev_sim_dir $prev_sims failed: $!";

#
# Build list of NR sources. We start with the directories in the reference
# SEED's NR dir, and override with anything in the biodb NR dir.
#


#
# Scan for SEED organisms.
#

my %NR_dirs;
scan_seed_dir(\%NR_dirs, $FIG_Config::organisms);


open(F, ">$jobdir/all.nr.dirs");
open(F2, ">$jobdir/nr.sources");
for my $d (sort bydb keys %NR_dirs)
{
    print F join("\t", $d, @{$NR_dirs{$d}}{'path', 'size'}), "\n";
    print F2 $NR_dirs{$d}->{fasta_path} . "\n";
}
close(F);
close(F2);

sub bydb
{
    if ($a =~ /^(\d+)\.(\d+)$/)
    {
	my($ga, $ia) = ($1, $2);
	if ($b =~ /^(\d+)\.(\d+)$/)
	{
	    my($gb, $ib) = ($1, $2);
	    return $ga <=> $gb or $ia <=> $ib;
	}
	else
	{
	    return 1;
	}
    }
    elsif ($b =~ /^\d+\.\d+$/)
    {
	return -1;
    }
    else
    {
	return $a cmp $b;
    }
}


sub validate_dirs
{
    my(@dirs) = @_;

    my $err;
    for my $dir (@dirs)
    {
	if (! -d $dir)
	{
	    warn "Required directory $dir is not present\n";
	    $err++;
	}
    }
    exit(1) if $err;
}


sub create_new_job
{
    my ($jobs_dir) = @_;

    # init job counter if necessary
    umask 0000;
    unless (-f "$jobs_dir/JOBCOUNTER") {
	open(FH, ">$jobs_dir/JOBCOUNTER") or die "could not create jobcounter file: $!\n";
	print FH "000\n";
	close FH;
    }
    
    # get new job id from job counter
    open(FH, "$jobs_dir/JOBCOUNTER") or die "could not open jobcounter file: $!\n";
    my $jobnumber = <FH>;
    chomp $jobnumber;
    $jobnumber = sprintf("%03d", $jobnumber + 1);
    close FH;
    while (-d $jobs_dir.'/'.$jobnumber)
    {
	$jobnumber = sprintf("%03d", $jobnumber + 1);
    }

    # create job directory
    my $job_dir = $jobs_dir.'/'.$jobnumber;
    mkdir $job_dir;
    open(FH, ">$jobs_dir/JOBCOUNTER") or die "could not write to jobcounter file: $!\n";
    print FH $jobnumber;
    close FH;
    
    unless (-d $job_dir) {
	return (undef, 'The job directory could not be created.');
    }

    return ($jobnumber,'');
}
