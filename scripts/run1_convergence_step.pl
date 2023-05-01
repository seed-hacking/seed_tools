use strict;
use Data::Dumper;
use File::Path 'make_path';
use Proc::ParallelLoop;
use Getopt::Long::Descriptive;

# 1) You need to keep $dir/Seqs up-to-date with coreSEED genomes.
#    Seqs/GENOME must contain the fasta file of peg translations.
# 2) Annotations/0/GENOME must contain the current coreSEED assignments
#

my($opt, $usage) = describe_options("%c %o Step-Number Dir",
				    ["parallel|p=i" => "Number of processes to use in running annotations", { default => 4 }],
				    ["help|h" => "Show this help message"]);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 2;

my $n = shift @ARGV;
my $dir = shift @ARGV;

if (! -d $dir)
{
    die "Data directory $dir not present\n";
}

my $last = $n-1;
print STDERR "going from $last to $n\n";
use FIG;
my $fig = new FIG;

opendir(SEQS,"$dir/Seqs") || die "Cannot opendir $dir/Seqs: $!";
my @genomes = grep { $_ =~ /^\d+\.\d+$/ } readdir(SEQS);
closedir(SEQS);
#goto x;
foreach my $g (@genomes)
{
    make_path("$dir/Organisms/$g/Features/peg");
    if (! -s "$dir/Organisms/$g/Features/peg/fasta")
    {
	symlink("$dir/Seqs/$g", "$dir/Organisms/$g/Features/peg/fasta");
    }
    &SeedUtils::run("cp $dir/Annotations/$last/$g $dir/Organisms/$g/assigned_functions");
}

make_path("$dir/Annotations/$n");
make_path("$dir/Calls/$n");
make_path("$dir/New/$n");

&make_Data($fig,$dir,$last,\@genomes);

x:
pareach \@genomes, sub { 
    my $g = shift;

    &SeedUtils::run("kmer_search -d $dir/Data.$last -a -m 10 -g 50 < $dir/Organisms/$g/Features/peg/fasta > $dir/Calls/$n/$g");

    #
    # Create merged annotations in Annotations/<new>/genome
    # W
    #
    my %anno = map { ($_ =~ /^(\S+)\t(\S[^\t]*\S)/) ? ($1 => $2) : () } `cat $dir/Annotations/$last/$g`;
    foreach $_ (`cat $dir/Calls/$n/$g`)
    {
	if ($_ =~ /^(\S+)\t(\S[^\t]*\S)/)
	{
	    $anno{$1} = $2;
	}
    }
    open(MERGED,">$dir/Annotations/$n/$g") || die "could not open $dir/Annotations/$n/$g";
    foreach $_ (sort { &SeedUtils::by_fig_id($a,$b) } keys(%anno))
    {
	print MERGED join("\t",($_,$anno{$_})),"\n";
    }
    close(MERGED);

    open(NEW,">$dir/New/$n/$g") || die "could not open $dir/New/$n/$g";
    my %anno1 = map { ($_ =~ /^(\S+)\t(\S[^\t]*\S)/) ? ($1 => $2) : () } `cat $dir/Annotations/$last/$g`;
    foreach $_ (`cat $dir/Calls/$n/$g`)
    {
	if ($_ =~ /^(\S+)\t(\S[^\t]*\S)/)
	{
	    my $x = $anno1{$1};
	    if ($x && $2 && ($x ne $2))
	    {
		print NEW join("\t",($1,$x,$2)),"\n";
	    }
	}
    }
    close(NEW);
}, { Max_Workers => $opt->parallel };

if (-f "$dir/New/$n/83333.1")
{
    my @changed_ec = `cat $dir/New/$n/83333.1`;
    if (@changed_ec > 500)
    {
	die "too many changes - check $dir/New directory";
    }
}

sub nxt_iter {
    my($dir) = @_;

    opendir(ANNO,"$dir/Annotations") || die "could not open $dir/Annotations";
    my(@done) = sort { $b <=> $a } grep { ($_ =~ /^\d+$/) &&
					  (-s "$dir/Annotations/$_/83333.1") } readdir(ANNO);
    my $last = (@done > 0) ? $done[0] : 0;
    return $last+1;
}

sub make_Data {
    my($fig,$dir,$last,$genomes) = @_;

    make_path("$dir/Data.$last");
    system("cp", "$dir/subsystem.roles", "$dir/Data.$last/subsystem.roles") if -f "$dir/subsystem.roles";
    system("cp", "$dir/additional.funcs", "$dir/Data.$last/additional.funcs") if -f "$dir/additional.funcs";
    open(GS,">$dir/Data.$last/genomes") || die "could not open $dir/Data.$last/genomes";
    my $agenome;
    my $unk = 1;
    foreach my $g (@$genomes)
    {
	my $gs = $fig->genus_species($g);
	if (!$gs)
	{
	    if (open(GN, "$dir/gnames/$g"))
	    {
		$gs = <GN>;
		chomp $gs;
		close(GN);
	    }
	    else
	    {
		warn "Cannot open $dir/gnames/$g: $!";
		$gs = "unknown genome $unk";
		$unk++;
	    }
	}
	print GS join("\t",($gs,$g)),"\n";
	$agenome = $g;
    }
    close(GS);
#   unlink("$dir/Data.$last/final.kmers");
#   unlink("$dir/Data.$last/kmer.table.mem_map");
    my $cmd = "kmer_search --allow-rebuild -r $dir/Organisms -a -d $dir/Data.$last < $dir/Organisms/$agenome/Features/peg/fasta > $dir/Data.$last/test.out";
    print STDERR "$cmd\n";
    &SeedUtils::run($cmd);
}
