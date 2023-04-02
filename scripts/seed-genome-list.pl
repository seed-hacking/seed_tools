=head1 NAME 

seed-genome-list

=head1 SYNOPSIS

seed-genome-list [--from-file ... ] create name [genome-id...]
seed-genome-list [--from-file ... ] create name [genome-id...]
seed-genome-list get name
seed-genome-list list

=cut

use strict;
use GenomeLists;
use Data::Dumper;
use Getopt::Long::Descriptive;

my %commands = (create => [\&do_create, sub { @ARGV > 1 }],
		get => [\&do_get, sub { @ARGV == 1 }],
		list => [\&do_list, sub { @ARGV == 0}],
    );

my($opt, $usage) = describe_options("%c %o cmd [opts]",
				    ["from-file=s", "Read genome IDs from given file instead of command line"],
				    ["description=s", "Description for a new list"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 1;

my $cmd = shift @ARGV;

my($func, $check) = @{$commands{$cmd}};

$func or die $usage->text;

&$check() or die $usage->text;

&$func();

sub do_list
{
    my @list = GenomeLists->getListsForUser();
    print "$_\n" foreach @list;
}

sub do_get
{
    my $name = shift @ARGV;

    my $list = GenomeLists::load($name);
    if (!$list)
    {
	die "List $name not found\n";
    }
    print "$_\n" foreach @{$list->{genomes}};
}

sub do_create
{
    my $name = shift @ARGV;
    my @genomes;
    if (@ARGV)
    {
	@genomes = @ARGV;
    }
    if ($opt->from_file)
    {
	if (open(my $fh, "<", $opt->from_file))
	{
	    while (<$fh>)
	    {
		if (/^\s*(\d+\.\d+)\s*$/)
		{
		    push(@genomes, $1);
		}
		else
		{
		    die "Invalid genome id at line $. of " . $opt->from_file;
		}
	    }
	    close($fh);
	}
    }
    GenomeLists::create($name, $opt->description, \@genomes);
}

