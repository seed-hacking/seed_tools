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

#
use FIG;
$fig = new FIG;

# usage: peg_to_subsystems.pl < PEG > PEG-FUNC-SUBSYSTEM-ROLE

while (defined($_ = <STDIN>))
{
    if ($_ =~ /(fig\|\d+\.\d+\.peg\.\d+)$/)
    {
        $peg = $1;
        $func = $fig->function_of($peg);
        foreach $x ($fig->subsystems_for_peg($peg))
        {
            ($subsystem,$role) = @$x;
	    $curator = $fig->subsystem_curator($subsystem);
            $curator =~ s/^master://;
            print "$peg\t$func\t$subsystem\t$role\t$curator\n";
        }
    }
}
