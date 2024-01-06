#
# Find sequences in peg.synonyms that don't appear as id1 in the sims collection.
#
# Initial implementation; scan all sims and retain a hash of the id1. See if we can do this on the 
# more limited memory systems.
#
use strict;

use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o sim-dir",
	["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 1;

my $dir = shift;

my %found;

for my $sfile (<$dir/Sims*/*>)
{
    open(S, "<", $sfile) or die "Cannot open $sfile:$ !";
    print STDERR "$sfile\n";
    while (<S>)
    {
	my($id) = /^([^\t]+)/;
	$found{$id}++;
    }
    close(S);
}

print STDERR "Scan peg.syn\n";

# gnl|md5|000036fa7e8cefb24ccbb49c395e55b6,287    fig|457393.3.peg.343,287;fig|820.10001.peg.2231,287;fig|820.10003.peg.2950,287

open(PS, "<", "$dir/peg.synonyms") or die "Cannot open $dir/peg.synonyms: $!";

while (<PS>)
{
    if (/^([^\t,]+),/)
    {
	if (!$found{$1})
	{
	    print "$1\n";
	}
    }
    else
    {
	warn "Bad peg.synonyms parse $. $_";
    }
}
