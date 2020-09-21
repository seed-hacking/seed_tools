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

use FIG;
my $fig = new FIG;

# usage: pegs_in_conflict [FileOfSubsystems] > PEG 

if (@ARGV > 0)
{
    %subsys = map { $_ =~ /^\s*(\S[^\t]+\S)/; $sub = $1; $sub =~ s/ /_/g; $sub => 1 }
              `cat $ARGV[0]`;
}

open(PEGS,"pegs_in_subsystems all | function_of |")
    || die "could not get pegs_in_subsystems | function_of";

while (defined($_ = <PEGS>))
{
    chomp;
    ($sub,$role,$peg,undef,$func) = split(/\t/,$_);
    $sub =~ s/ /_/g;
    next if ((@ARGV > 0) && (! $subsys{$sub}));

    @roles = $fig->roles_of_function($func);
    for ($i=0; ($i < @roles) && ($roles[$i] ne $role); $i++) {}
    if ($i == @roles)
    {
	$bad{$peg} = 1;
    }
}
close(PEGS);

foreach $peg (sort { &FIG::by_fig_id($a,$b) } keys(%bad))
{
    print "$peg\n";
}
