#
# Use BV-BRC RASTtk API to reannotate proteins in the given genome
#

use strict;
use Data::Dumper;
use GenomeTypeObject;
use FIG;
use IPC::Run qw(run);
use POSIX;
use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o user genome-id",
				    ["save-intermediates" => "Show and save the temp files for genome annotation"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit (0) if $opt->help;
die($usage->text), if @ARGV != 2;

my $user = shift;
my $genome = shift;

my $now = strftime("%Y-%m-%d-%H-%M-%S", localtime);

my $fig = new FIG;

my $name = $fig->genus_species($genome);
$name or die "Could not look up $genome\n";

$name =~ s/^(\S+\s+\S+).*$/$1/;

my $dir = $fig->organism_directory($genome);
print "Saving to $dir\n";

my $gobj = GenomeTypeObject->new;

$gobj->set_metadata({ scientific_name => $name });

my $feats = $fig->all_features_detailed_fast($genome);
my @sorted =  sort { $a->[1] cmp $b->[1] or $a->[4] <=> $b->[4] } map { $_->[1] =~ s/_\d+_\d+$//; $_ } @$feats;

my @pegs = map { [@$_[0,6]] } grep { $_->[3] eq 'peg' } @sorted;

my $trans = $fig->get_translation_bulk([ map { $_->[0] } @pegs]);

my %orig = map { $_->[0] => $_->[1] } @pegs;

for my $ent (@pegs)
{
    my($id, $func) = @$ent;
    my $aa = $trans->{$id};
    $gobj->add_feature({
	-id => $id,
	    -location => [0,0],
	    -type => 'peg',
	    # -function => $func,
	    -protein_translation => $aa,
		       });
}

my $tmp_in = File::Temp->new(UNLINK => $opt->save_intermediates);
my $tmp_out = File::Temp->new(UNLINK => $opt->save_intermediates);

if ($opt->save_intermediates)
{
    print "Genome files: $tmp_in $tmp_out\n";
}

$tmp_in->close();
$tmp_out->close();

$gobj->destroy_to_file("$tmp_in");


if (1)
{
    my $ok = run(["rast-annotate-proteins-kmer-v2", "-i", "$tmp_in"],
		 "|",
		 ["rast-annotate-families-patric", "-o", "$tmp_out"]);
    $ok or die "Cannot annotate\n";
}
else
{
    $tmp_in = "z.json";
}

my $gobj = GenomeTypeObject->new({ file => "$tmp_out" });

open(PF, ">", "$dir/pattyfams.txt") or die "Cannot open $dir/pattyfams.txt: $!";
open(NF, ">", "$dir/reannoated.$now.txt") or die "Cannot open $dir/reannotated.$now.txt: $!";
open(CHANGES, ">", "$dir/changes.$now.txt") or die "Cannot open $dir/changes.$now.txt: $!";

for my $feat ($gobj->features)
{
    my $id = $feat->{id};
    my $fun = $feat->{function};

    if ($fun && $fun ne $orig{$id})
    {
	print NF "$id\t$fun\n";
	print CHANGES "$id\t$orig{$id}\t$fun\n";
    }
    for my $fam (@{$feat->{family_assignments}})
    {
	my($which, $fam, $fun, $ver, $score) = @$fam;
	$score //= 0;
	print PF join("\t", $id, $fam, $score, $fun), "\n";
    }
}
close(PF);
print "Load annotations\n";
my $ok = run(["fig", "assign_functionF", $user, "$dir/reannoated.$now.txt"]);
$ok or die "assignment failed\n";
print "Load families\n";
my $ok = run(["load_pattyfams", $genome]);
$ok or die "assignment failed\n";
