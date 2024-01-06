#
# Given a list of IDs for which we have no sims, look in the sims files
# for entries where the ID is a sim target and generate flipped sims.
#
use strict;

use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o missing-ids sim-dir",
	["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 2;

my $missing_file = shift;
my $dir = shift;


my %missing;
open(M, "<", $missing_file) or die "Cannot open $missing_file: $!";
while (<M>)
{
    chomp;
    $missing{$_} = 1;
}
close(M);

for my $sfile (<$dir/Sims*/*>)
{
    open(S, "<", $sfile) or die "Cannot open $sfile:$ !";
    print STDERR "$sfile\n";
    while (<S>)
    {
	chomp;
	my @sim = split(/\t/);
	if ($missing{$sim[1]})
	{
	    print join("\t", $sim[1], $sim[0], @sim[2..$#sim]), "\n";
	}
    }
    close(S);
}
