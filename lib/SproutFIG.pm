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

package SproutFIG;

#
# This package returns a FIG object or a Sprout object translated to use the
# FIG API, depending on whether or not the constructor has any operands.
#
# 	my $sprout_or_fig = SproutFig->new()
#
# creates a FIG object.
#
#   my $sprout_or_fig = SproutFig->new('', '')
#
# creates a translated Sprout object.

use strict;
use Carp;

#
# Safely import Sprout.
#
my $sproutAvail;
eval {
    require Sprout;
    import Sprout;
    require SFXlate;
    import SFXlate;
    $sproutAvail = 1;
};

if ($@)
{
    #
    # Croak on errors other than can't find sprout.
    #
    if ($@ !~ /Can\'t locate Sprout.pm/)
    {
	die $@;
    }
}

use FIG_Config;
use SeedDas;

sub new
{
    my($class, $sproutDB, $sproutData) = @_;
    my $fig = FIG->new();

    if (!defined $sproutDB)
    {
	return $fig;
    }
    elsif ($sproutAvail)
    {
    	my $sprout = SFXlate->new($fig, $sproutDB, $sproutData);
	return $sprout;
    }
    else
    {
	confess "SproutFIG: Attempt to use Sprout database when Sprout code is not available\n";
    }
}


1;
