########################################################################
use strict;
use Data::Dumper;
use gjoseqlib;
use File::Path 'make_path';
use File::Basename;
use Getopt::Long::Descriptive;
use FIG_Config;
use FIG;

my %supported_seeds = (core => "/vol/core-seed",
		       open => "/vol/open-seed",
		       bob => '/home/olson',
		       pubseed => "/vol/public-pseed",
		       pseed => "/vol/pseed");

my($opt, $usage) = describe_options("%c %o data-dir",
				    ['skip-pegs=s', 'Skip the pegs in this file'],
				    ['skip-weighted-score-computation', "Skip the pass to computed weighted scores"],
				    ["keep-intermediates", "Keep intermediate files instead of deleting when finished"],
				    ['additional-fasta=s@', "Use this fasta file to mix in additional sequences to build kmers"],
				    ['sam-data=s@', "Use the data formatted from Sam Seaver to mix in additional sequences to build kmers"],
				    ["id-map=s", "File to which additional fasta file entries' id mapping is written"],
				    ["no-strip", "Don't strip comments from functions"],
				    ["protect-subsystem-roles", "Keep kmers for functions that contain a subsystem role"],
				    ["min-reps-required=i", "Number of proteins required to form a family",
				     { default => 5 }],
				    ["seed=s", "Choose which seed to draw genomes from", { default => 'core' }],
				    ["this-seed", "Force use of this SEED"],
				    ["sort-memory", "Use this much sort memory", { default => "4G" }],
				    ["virus-dir=s", "Incorporate virus genomes from here"],
				    ["function-overrides=s", "File of function overrides"],
				    ["parallel|p=i" => "Number of processes to use in running annotations", { default => 4 }],
				    ["otu-reps=s" => "Use the representatives from the OTUs in this file as the genome set"],
				    ['genome=s@' => "Use this genome to create kmers (may be repeated)"],
				    ["help|h" => "Show this help message"]);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

$ENV{KM_MIN_REPS_REQUIRED} = $opt->min_reps_required;
$ENV{KM_SKIP_WEIGHTED_SCORE_COMPUTATION} = $opt->skip_weighted_score_computation;
$ENV{KM_KEEP_INTERMEDIATES} = $opt->keep_intermediates;

my $seed_path = $supported_seeds{$opt->seed};
if (! -d $seed_path && $opt->this_seed)
{
    $seed_path = dirname($FIG_Config::fig_disk);
    print "set to $seed_path\n";
}
if (!$seed_path)
{
    die $opt->seed . " is not a supported SEED. Valid SEEds are " . join(" ", sort keys %supported_seeds) . "\n";
}

if ($FIG_Config::fig_disk ne "$seed_path/FIGdisk")
{
    die "This script must be run with the chosen SEED's configuration sourced\nConfig says $FIG_Config::fig_disk seedpath=$seed_path\n";
}

my $dir = shift @ARGV;
if (! -d $dir)
{
    die "The data directory $dir must already exist\n";
}

my @genomes;
if ($opt->otu_reps)
{
    open(O, "<", $opt->otu_reps) or die "The OTU-reps file " . $opt->otu_reps . " cannot be opened: $!\n";

    my %seen;
    while (<O>)
    {
	chomp;
	my($set, $id, $name) = split(/\t/);
	if (!$seen{$set})
	{
	    push(@genomes, $id);
	    $seen{$set} = 1;
	}
    }
    close(O);
}

if ($opt->genome)
{
    push(@genomes, @{$opt->genome});
}

my %skip_pegs;
if ($opt->skip_pegs)
{
    open(F, "<", $opt->skip_pegs) or die "Cannot open skip-pegs file " . $opt->skip_pegs . ": $!";
    while (<F>)
    {
	chomp;
	my($peg, $fun) = split(/\t/);
	$skip_pegs{$peg} = 1;
    }
    close(F);
}

make_path("$dir/Seqs", "$dir/Annotations/0", "$dir/gnames");

my $cs_orgs = $FIG_Config::organisms;
my $fig = FIG->new();

if (!@genomes)
{
    opendir(CS,$cs_orgs) || die "could not open $cs_orgs: $!";
    while (my $g = readdir(CS))
    {
	if ($g =~ /^\d+\.\d+$/)
	{
	    push(@genomes, $g);
	}
    }
    
    closedir(CS);
}

if ($opt->virus_dir)
{
    my %viral_roles;
    opendir(D, $opt->virus_dir . "/fasta") or die "Cannot opendir " . $opt->virus_dir . "/fasta: $!";
    my @virus_genomes = sort grep { /^\d+\.\d+$/ } readdir(D);
    closedir(D);

    for my $org (@virus_genomes)
    {
	my $fasta_file = $opt->virus_dir . "/fasta/$org";
	my $anno_file = $opt->virus_dir . "/anno/$org";

	if (! -f $anno_file)
	{
	    die "Missing anno file $anno_file\n";
	}
	
	symlink($fasta_file, "$dir/Seqs/$org") or die "cannot symlink $fasta_file $dir/Seqs/$org: $!";
	open(ANN,">$dir/Annotations/0/$org") || die "could not open $dir/Annotations/0/$org";

	open(FROM, "<", $anno_file) or die "Cannot open $anno_file: $!";
	while (<FROM>)
	{
	    chomp;
	    my($id, $func) = split(/\t/);

	    $viral_roles{$func} = 1;
	    print ANN "$id\t$func\n";
	}
	close(ANN);
	close(FROM);
    }
    # if ($opt->keep_viral_roles)
    # {
    # 	open(VRF, ">", "$kmer_dir/viral.roles");
    # 	print VRF "$_\n" foreach sort keys %viral_roles;
    # 	close(VRF);
    # 	push(@good_roles, "--good-roles", "$kmer_dir/viral.roles");
    # }
}

for my $g (@genomes)
{
    if (! -s "$dir/Seqs/$g")
    {
	next if (! -d "$cs_orgs/$g");
	print "Adding $g\n";
	my %fasta = map { ($_->[0] => $_->[2]) } &gjoseqlib::read_fasta("$cs_orgs/$g/Features/peg/fasta");
	open(SEQS,">$dir/Seqs/$g") || die "could not open $dir/Seqs/$g";
	foreach my $peg (keys(%fasta))
	{
	    if ($skip_pegs{$peg})
	    {
		print "Skipping $peg in seqs file\n";
		next;
	    }
	    print SEQS ">$peg\n$fasta{$peg}\n";
	}
	close(SEQS);
    }
    if (1) # (! -s "$dir/Annotations/0/$g")
    {
	open(ANN,">$dir/Annotations/0/$g") || die "could not open $dir/Annotations/0/$g";

	my %suffix;
	my $feats = $fig->all_features_detailed_fast($g);
	my @feats = sort { FIG::by_locus($a->[1], $b->[1]) } @$feats;
	
	my $i = 0;
	my $j;
	my @runs;
	while ($i < @feats)
	{
	    $j = $i + 1;
	    my $elt = $feats->[$i];
	    next if $elt->[3] ne 'peg';

	    # print "$i\t$elt->[0]\t$elt->[6]\n";
	    
	    my ($fun, $comment) = SeedUtils::strip_func_comment($elt->[6]);
	    $comment =~ s/^\s*\#\s*//;

	    next if $comment !~ /(fragment|frameshift)/;

	    while ($j < @feats)
	    {
		my $elt2 = $feats->[$j];
		next if $elt2->[3] ne 'peg';

		my ($fun2, $comment2) = SeedUtils::strip_func_comment($elt2->[6]);
		$comment2 =~ s/^\s*\#\s*//;

		# print "comp\t$i\t$fun\t$j\t$fun2\n";

		last if $fun ne $fun2 || $comment2 !~ /(fragment|frameshift)/;
	    }
	    continue {
		$j++;
	    }
	    if ($j - $i > 1)
	    {
		my($undef, $beg, $end, $strand) = $fig->boundaries_of($elt->[1]);
		print STDERR "$i $elt->[0] $elt->[1] $strand\n";
		my @locs = $i..$j-1;
		if ($strand eq '-')
		{
		    @locs = reverse(@locs);
		}
		my $val = 1;
		for my $x (@locs)
		{
		    $suffix{$x} = $val++;
		}
	    }
	}
	continue
	{
	    $i = $j;
	}
	       
	for my $i (0..$#feats)
	{
	    my $elt = $feats->[$i];
	    next if $elt->[3] ne 'peg';

	    my $func = $elt->[6];
	    
	    if ($suffix{$i})
	    {
		$func .= " $suffix{$i}";
	    }
	    else
	    {
		$func = SeedUtils::strip_func_comment($func);
	    }
	    if ($skip_pegs{$elt->[0]})
	    {
		print "Skipping $elt->[0] in anno file\n";
		next;
	    }
	    print ANN "$elt->[0]\t$func\n";
	}
	close(ANN);
	
# 	my %deleted;
# 	open(ANN,">$dir/Annotations/0/$g") || die "could not open $dir/Annotations/0/$g";
# 	if (-s "$cs_orgs/$g/Features/peg/deleted.features")
# 	{
# 	    %deleted = map { chomp; ($1 => 1) } 
# 	               `cat $cs_orgs/$g/Features/peg/deleted.features`;
# 	}
# 	my %tmp_funcs = map { ($_ =~ /^(\S+)\t(\S[^#\t]*\S)/) ? ($1 => $2) : () }
# 	                `cat $cs_orgs/$g/assigned_functions`;
# 	foreach my $peg (keys(%tmp_funcs))
# 	{
# 	    if (! $deleted{$peg})
# 	    {
# 		print ANN join("\t",($peg,$tmp_funcs{$peg})),"\n";
# 	    }
# 	}
# 	close(ANN);
    }
}

my $gbase = 7777777;
my $gidx = 1;

if ($opt->additional_fasta)
{
    my $mfh;
    if ($opt->id_map)
    {
	open($mfh, ">", $opt->id_map) or die "Cannot open " . $opt->id_map . " for writing: $!";
    }

    my %additional_funcs;

    for my $fn (@{$opt->additional_fasta})
    {
	open(my $fh, "<", $fn) or die "Cannot open fasta file $fn: $!";

	my $g = "$gbase.$gidx";
	$gidx++;

	open(GN, ">", "$dir/gnames/$g") or die "Cannot open $dir/gnames/$g: $!";
	print GN "Extra " . basename($fn) . "\n";
	close(GN);
    
	open(SEQS, ">", "$dir/Seqs/$g") || die "could not open $dir/Seqs/$g";
	open(ANN,">$dir/Annotations/0/$g") || die "could not open $dir/Annotations/0/$g";

	my $next_id = 1;
	
	while (my($id, $def, $seq) = read_next_fasta_seq($fh))
	{
	    $def =~ s/\s+\[[^\]]+\]\s*$//;
	    my $nid = "fig|$g.peg.$next_id";
	    $next_id++;
	    print $mfh "$id\t$nid\n" if $mfh;
	    
	    print SEQS ">$nid\n$seq\n";
	    my ($fun, $comment) = SeedUtils::strip_func_comment($def);
	    $comment =~ s/^\s*\#\s*//;

    	    $additional_funcs{$fun}++;
	    
	    print ANN "$nid\t$fun\n";
	}
	close(ANN);
	close(SEQS);
    }

    open(D, ">", "$dir/additional.funcs") or die "Cannot write $dir/additional.funcs: $!";
    print D "$_\n" foreach sort { $a cmp $b } keys %additional_funcs;
    close(D);
}

if ($opt->sam_data)
{
    my $g = "$gbase.$gidx";
    $gidx++;
    
    my $mfh;
    if ($opt->id_map)
    {
	open($mfh, ">>", $opt->id_map) or die "Cannot open " . $opt->id_map . " for writing: $!";
    }

    open(SEQS, ">", "$dir/Seqs/$g") || die "could not open $dir/Seqs/$g";
    open(ANN,">$dir/Annotations/0/$g") || die "could not open $dir/Annotations/0/$g";

    my %additional_funcs;

    for my $fn (@{$opt->sam_data})
    {
	open(my $fh, "<", $fn) or die "Cannot open fasta file $fn: $!";
	
	my $next_id = 1;

	while (<$fh>)
	{
	    chomp;
	    my($id, $seq, $func) = split(/\t/);
	    
	    my $nid = "fig|$g.peg.$next_id";
	    $next_id++;
	    print $mfh "$id\t$nid\n" if $mfh;
	    
	    print SEQS ">$nid\n$seq\n";
	    my ($fun, $comment) = SeedUtils::strip_func_comment($func);

    	    $additional_funcs{$fun}++;
	    
	    print ANN "$nid\t$fun\n";
	}
    }
    close(ANN);
    close(SEQS);

    open(D, ">", "$dir/additional.funcs") or die "Cannot write $dir/additional.funcs: $!";
    print D "$_\n" foreach sort { $a cmp $b } keys %additional_funcs;
    close(D);
}

#
# If we have a function overrides file, append these to the individual
# per-genome annotation files in $dir/Annotations/0
#
if ($opt->function_overrides)
{
    my %override;
    open(F, "<", $opt->function_overrides) or die "Cannot read " . $opt->function_overrides . ": $!";
    while (<F>)
    {
	chomp;
	my($id, $fn) = split(/\t/);
	push(@{$override{genome_of($id)}}, [$id, $fn]);
    }
    close(F);

    while (my($genome, $list) = each %override)
    {
	my $afile = "$dir/Annotations/0/$genome";
	-f $afile or warn "Annotation file $afile is missing\n";
	open(A, ">>", $afile) or die "Cannot append to $afile: $!";
	print A join("\t", @$_), "\n" foreach @$list;
	close(A);
    }
}

if ($opt->protect_subsystem_roles && ! -s "$dir/subsystem.roles")
{
    my $list = $fig->subsystem_roles();
    open(D, ">", "$dir/subsystem.roles") or die "Cannot write $dir/subsystem.roles: $!";
    print D "$_\n" foreach sort { $a cmp $b } keys %$list;
    close(D);
}

if ($opt->no_strip)
{
    $ENV{KM_BUILD_NO_STRIP} = 1;
}

unlink("$dir/Data.0/final.kmers");
unlink("$dir/Data.1/final.kmers");
unlink("$dir/Data.2/final.kmers");
my $par = "--parallel " . $opt->parallel;
my $mem = "--sort-memory " . $opt->sort_memory;
&SeedUtils::run("run1_convergence_step $par $mem 1 $dir");
#exit;
&SeedUtils::run("run1_convergence_step $par $mem 2 $dir");
&SeedUtils::run("run1_convergence_step $par $mem 3 $dir");
print STDERR "now computing stats to support z-scores on density\n";
&SeedUtils::run("compute_stats_for_kmer_hits -d $dir/Data.2");
print STDERR "$dir/Data.2 should be ready to go\n";
