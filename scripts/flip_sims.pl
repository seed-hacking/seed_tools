#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

use strict;

my $usage = "usage: flip_sims sims.dir sorted\n";
my($in_dir, $sorted);

(
 ($in_dir  = shift @ARGV) &&
 ($sorted = shift @ARGV) &&
 (@ARGV == 0)
)
    || die $usage;


#
# if in_dir is actually a file, just flip that file.
#

my @files;
if (-f $in_dir)
{
    @files = ($in_dir);
}
else
{
    opendir(SIMS,$in_dir)
	|| die "could not open $in_dir";
    @files = grep { $_ !~ /^\./ } readdir(SIMS);
    closedir(SIMS);
    @files = map { "$in_dir/$_" } @files;
}

open(SORT, "| sort $ENV{SORT_ARGS} > $sorted") or die "cannot open $sorted: $!\n";
foreach my $new_sims (@files)
{
    open(NEWSIMS,"<$new_sims")
	|| die "could not open $new_sims";
    while (defined(my $sim = <NEWSIMS>))
    {
	if ($sim =~ /^\S+\t(\S+)/)
	{
	    my $rev = &rev_sim($sim);
	    print SORT $rev;
        }
    }
    close(NEWSIMS);
}

close(SORT);

sub rev_sim {
    my($sim) = @_;

    chomp $sim;
    my($id1,$id2,$iden,$ali_ln,$mismatches,$gaps,$b1,$e1,$b2,$e2,$psc,$bsc,$ln1,$ln2) = split(/\t/,$sim);
    return join("\t",($id2,$id1,$iden,$ali_ln,$mismatches,$gaps,$b2,$e2,$b1,$e1,$psc,$bsc,$ln2,$ln1)) . "\n";
}
