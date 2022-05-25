#
# Create a btree containing xxx => FIG id mapping that is used 
# in reduce_sims.
#

my $usage = "make_reduce_sims_index peg.syn peg.syn.btree";

use strict;
use FIG;
use DB_File;

@ARGV == 2 or die $usage;

my $pegsyn_file = shift;
my $btree_file = shift;

open(SYN, "<$pegsyn_file") or die "Cannot open $pegsyn_file: $!";
#
# Unlink existing file to ensure it is clean.
#
if (-f $btree_file)
{
    unlink($btree_file);
}

my %btree;
my $btree_tie = tie %btree, 'DB_File', $btree_file, O_RDWR | O_CREAT, 0666, $DB_BTREE;

$btree_tie or die "Cannot create btree $btree_file: $!\n";

while (defined($_ = <SYN>))
{
    if ($_ =~ /^([^,]+),\d+\t(\S+)/)
    {
	my $maj = $1;
	my $rest = $2;
	if (($maj !~ /^fig/) && ($rest =~ /(fig\|\d+\.\d+\.peg\.\d+)/))
	{
	    $btree{$maj} = $1;
	}
    }
}
close(SYN);
untie %btree;
