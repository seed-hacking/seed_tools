# -*- perl -*-
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
$fig = new FIG;

my $usage = "usage: pegs_in_subsystems [ all ] > SUBSYS-ROLE-PEG\nThe 'all' parameter is required to include genomes without an operational variant\n\n";

my $all = ((@ARGV > 0) && ($ARGV[0] =~ /^all/i));

foreach $subsys (grep { $fig->usable_subsystem($_) } $fig->all_subsystems)
{
    @roles    = $fig->subsystem_to_roles($subsys);
    if (! $all)
    {
	@roles = grep { ! $fig->is_aux_role_in_subsystem($subsys,$_) } @roles;
    }

    $genomes  = $fig->subsystem_genomes($subsys,$all);
    foreach $role (@roles)
    {
	foreach $pair (@$genomes)
	{
	    ($genome,$gs) = @$pair;
	    @pegs = $fig->pegs_in_subsystem_cell($subsys,$genome,$role);
	    foreach $peg (@pegs)
	    {
		print "$subsys\t$role\t$peg\n";
	    }
	}
    }
}
