use strict;
use FIG;
use File::Basename;
use Proc::ParallelLoop;
use PerlIO::via::Blockwise;
use DB_File;

# /home/olson/FIGdisk/FIG/bin/update_sims2 nr peg.synonyms 300 prev_sims Sims.005 sims.job.sge/sims.flipped ids.deleted


my $usage = "usage: update_sims nr peg.synonyms MaxSz InDir OutDir SortedNewSims IdsToDelete unused";

my($nr, $added,$deleted,$new_sims,@files,$file,$full_file,$hdr,%in,$sim,$syns,$maxsz);
my(@sims,$curr,$currQ,$x,%new,$remapped,@flds,%remapped,$new);
my($in_dir,$out_dir, $sorted_sims, $to_delete);

@ARGV == 6 or @ARGV == 7 or @ARGV == 8 or die $usage;

$nr      = shift @ARGV;
$syns    = shift @ARGV;
$maxsz   = shift @ARGV;
$in_dir  = shift @ARGV;
$out_dir = shift @ARGV;
$sorted_sims = shift @ARGV;
$to_delete = shift @ARGV;
my $unused = shift @ARGV;

#
# Read the NR so we can eliminate sims for sequences that
# do not occur in the NR.
#

my %nr_seqs;
if ($nr =~ /btree$/)
{
    tie %nr_seqs, 'DB_File', $nr, 0, O_RDONLY, $DB_BTREE or die "Cannot tie $nr: $!";
}
else
{
    print "Reading NR $nr\n";
    open(NR,"<:via(Blockwise)", $nr) || die "could not open $nr";
    
    while (<NR>)
    {
	if ($_ =~ /^>(\S+)/)
	{
	    $nr_seqs{$1}++;
	}
    }
    close(NR);
}

-f $syns or die "Synonyms file $syns not found\n";

my $unused_fh;
if ($unused)
{
    open($unused_fh, ">$unused") or die "cannot open unused file $unused: $!\n";
}

my %deleted;
if ($to_delete && -f $to_delete)
{
    open(DEL, "<$to_delete") or die "Cannot open delete file to_delete: $!";
    while (<DEL>)
    {
	if (/(\S+)/)
	{
	    $deleted{$1}++;
	}
    }
    close(DEL);
}

if (-d $in_dir)
{
    opendir(SIMS,$in_dir)
	|| die "could not open $in_dir";
    @files = map { "$in_dir/$_" }  grep { $_ !~ /^\./ } readdir(SIMS);
    closedir(SIMS);
}
elsif (open(SFILE, "<$in_dir"))
{
    @files = <SFILE>;
    chomp @files;
    close(SFILE);
}
else
{
    die "Bad sims input $in_dir\n";
}

#-d $out_dir and die "Output directory $out_dir already exists\n";


-d $out_dir || mkdir($out_dir,0777) || die "could not make $out_dir";

#
# Read the sorted sims, maintaining an index of peg -> seek & length.
#

open(S, "<$sorted_sims") or die "Cannot open $sorted_sims: $!\n";
$| = 1;

my %idx;
my $offset = tell(S);
my $line = <S>;

my $curr;
my $next_offset;
while (defined($line) and $line =~ /^(\S+)/)
{
    print "." if $. % 100000 == 0;
    $curr = $1;
    while (defined($line) and
	   $line =~ /^(\S+)/ and
	   $1 eq $curr)
    {
	$next_offset = tell(S);
	$line = <S>;
    }
    my $len = $next_offset - $offset;
    $idx{$curr} = [$offset, $len];

    $offset = $next_offset;
}
print "\n";
#
# Outer: Loop over the existing sims files.
#
@files = sort @files;
#foreach $file (@files)
#{
pareach \@files, sub {
    my $file = shift;
    print STDERR "processing $file\n";

    open(OLD,"<:via(Blockwise)", $file)
        || die "could not open $file";

    #
    # Shortcut if the old file is empty; don't want to pay the price
    # of starting reduce_sims.
    #
    my $base = basename($file);

    my @stat = stat(OLD);
    if (@stat and $stat[7] == 0)
    {
	open(NEW, ">$out_dir/$base");
	close(NEW);
	close(OLD);
	return;
    }

    open(NEW,"| reduce_sims $syns $maxsz > $out_dir/$base")
        || die "could not open $out_dir/$base";
#    open(NEW,"| reduce_sims $syns $maxsz | reformat_sims $nr > $out_dir/$base")
#        || die "could not open $out_dir/$base";

    #
    # Rescan the old sims file, adding the new sims.
    #
    my $n_read = 0;
    my $n_written = 0;
    my $n_bad = 0;
    $sim = <OLD>;
    while (defined($sim))
    {
        if ($sim =~ /^(\S+)\t(\S+)/)
        {
	    if (!(sim_ok($1) && sim_ok($2)))
	    {
		$sim = <OLD>;
		$n_bad++;
		next;
	    }
	    $n_read++;
            $curr = $1;
            $currQ = quotemeta $curr;
            @sims = ();
            while (defined($sim) && ($sim =~ /^$currQ\s/))
            {
                if ($sim =~ /^$currQ\t(\S+)/)
                {
		    my $ok = sim_ok($1);
		    if (sim_ok($1))
		    {
			push(@sims,$sim);
			$n_read++;
		    }
		    else
		    {
			$n_bad++;
		    }
		}
                $sim = <OLD>;
            }

	    if (my $where = $idx{$curr})
	    {
		my($seek, $len) = @$where;
		seek(S, $seek, 0);
		my $buf;
		read(S, $buf, $len);
		my @new = split(/\n/, $buf);

#		print "Got sims=@sims\n new=$buf\n";

		#
		# we have new sims available for this id1.
		#
		# Create a list of tuples [id1, id2, rest of simdata] where tup[10] is the e-score
		# Sort by e-score
		# Select the highest score id2 if there are multiple id2s.
		# and map back to tab-separated list to print.
		#
                my %seen;
                @sims = map  { join("\t",@$_) . "\n" }
		grep { if (! $seen{$_->[1]}) { $seen{$_->[1]} = 1; 1 } else { 0 } }
		sort {
                                ($a->[10] <=> $b->[10])
				}
		map {  chomp; [split(/\t/,$_)] }
		(@sims,@new);

		#
		# Remove the index entry for $curr. That way we can print them out
		# at the end of the process, if we have any left over.
		#
		delete $idx{$curr};
            }
	    $n_written += @sims;
            print NEW join("",@sims) or die "$0 cannot write to $out_dir/$base: $!";
        }
        else
        {
            $sim = <OLD>;
        }
    }
    print "$out_dir/$base\t$n_read\t$n_written\t$n_bad\n";
    close(NEW);
    close(OLD);
}, { Max_Workers => 3 };

sub sim_ok
{
    my($id1) = @_;
	
    return !$deleted{$id1} && exists($nr_seqs{$id1});
}
sub rev_sim {
    my($sim) = @_;

    chomp $sim;
    my($id1,$id2,$iden,$ali_ln,$mismatches,$gaps,$b1,$e1,$b2,$e2,$psc,$bsc,$ln1,$ln2) = split(/\t/,$sim);
    return join("\t",($id2,$id1,$iden,$ali_ln,$mismatches,$gaps,$b2,$e2,$b1,$e1,$psc,$bsc,$ln2,$ln1)) . "\n";
}
