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

use strict;
use DB_File;

# usage: reduce_sims peg.synonyms TargetN < sims > reduced.sims

my($syn, $targetN);
(($syn      = shift @ARGV) &&
 ($targetN  = shift @ARGV)
)
    || die "usage: reduce_sims peg.synonyms TargetN < sims > reduced.sims";

#
# If $syn ends in .btree, it is a precompute btree index for the data we're looking for.
#
# If it ends in .fig, it is a 2-column table of id and SEED organism.

my %to_peg;
my %to_org;
my $have_orgtable;

if ($syn =~ /\.fig$/)
{
    open(F, "<$syn") or die "cannot open $syn: $!";
    while (<F>)
    {
	chomp;
	my($id, $org) = split(/\t/);
	$to_org{$id} = $org;
    }
    $have_orgtable = 1;
}
elsif ($syn =~ /\.btree$/)
{
    my $btree_tie = tie %to_peg, 'DB_File', $syn, O_RDONLY, 0666, $DB_BTREE;

    $btree_tie or die "Cannot open btree $syn: $!\n";
}
else
{
    open(SYN,"<$syn") || die "could not open synonyms file '$syn': $!";
    
    #
    # Parse the peg.synonyms file.
    #
    # For each entry, determine if the major id is a FIG id.
    #
    # If it is not, find the first FIG id in the synonym list and
    # remember it in %to_peg. This hash maps from non-fig major ID
    # to the longest FIG id that it maps to.
    #
    
    while (defined($_ = <SYN>))
    {
	if ($_ =~ /^([^,]+),\d+\t(\S+)/)
	{
	    my $maj = $1;
	    my $rest = $2;
	    if (($maj !~ /^fig/) && ($rest =~ /(fig\|\d+\.\d+\.peg\.\d+)/))
	    {
		$to_peg{$maj} = $1;
	    }
	}
    }
    close(SYN);
}

#
# Read the sims entires from stdin.
#

my $sim = &get_input;
while (defined($sim) && ($sim =~ /^(\S+)/))
{
    #
    # $curr is the id of the sim entry we're currently examining.
    #
    
    my $curr = $1;
    my $currQ = quotemeta $curr;
    my @sims = ();
    my %id2s;

    #
    # Read the rest of the sims for this id, collecting them in
    # @sims.
    #
    # We discard entries where id1 == id2, and where we
    # have already seen an id2 entry for this peg.
    #
    # In other words, we keep the first of any id1/id2 pair.
    # (The above logic is also present in reformat_sims).
    #

    while ($sim && ($sim =~ /^$currQ\t(\S+)/))
    {
	#
	# $1 is id2 for this sim record.
	#
        if (($curr ne $1) && (! $id2s{$1}))
        {
            push(@sims,[$sim,$1]);
	    #
	    # Commenting out this line preserves the behavior incurred
	    # when there was a typo in that line, introduced in vers 1.6.
            #$id2s{$1} = 1;
        }
        $sim = &get_input;
    }

    #
    # If we have $targetN (from the command line) or fewer sims for
    # this peg, print them all.
    #

    my %pegorg;
    if (@sims <= $targetN)
    {
	foreach $_ (@sims)
	{
	    print $_->[0] or die "reduce_sims cannot write to stdout: $!";
	}
    }
    else
    {
	#
	# Otherwise, do this.
	#
	my $tot = 0;
	my %orgs;
	foreach $_ (@sims)
	{
	    my($sim1,$id2) = @$_;

	    my $org;
	    if ($have_orgtable)
	    {
		$org = $to_org{$id2} or "sp";
	    }
	    else
	    {
		
		#
		# Map $id2 to be a FIG id if we can.
		#

		if (my $id2T = $to_peg{$id2})
		{
		    $id2 = $id2T;
		}
		
		#
		# If $id2 is a FIG id, increment the count of
		# hits for that organism.
		#
		# If it is not, increment the hitcount for "sp".
		#
		# In any event, count the number of sims.
		#
		if ($id2 =~ /^fig\|(\d+\.\d+)/)
		{
		    $org = $1;
		}
	    }
	    # print STDERR "$id2 $org\n";
	    
	    if ($org)
	    {
		$tot++;
		$orgs{$org}++;
		$pegorg{$id2} = $org;
	    }
	    else
	    {
		$tot++;
		$orgs{"sp"}++;
	    }
	}

	my %cnt;

	#
	# Determine how many of each organisms pegs we retain.
	#
	# We target a total of $targetN pegs, and assign them
	# in proportion to the number of organisms that showed up
	# in the full set of sims. We ensure sim entry for each
	# organism is kept.
	#
	foreach my $org (keys(%orgs))
	{
	    $cnt{$org} = int(0.99999 + (($orgs{$org}/$tot) * $targetN));
	    $_ = @sims;
	    # print STDERR "$org $orgs{$org} $cnt{$org} $_\n";
	}

	foreach $_ (@sims)
	{
	    my ($sim1,$id2) = @$_;
	    if (my $id2T = $to_peg{$id2})
	    {
		$id2 = $id2T;
	    }
	    my $org = $pegorg{$id2} or "sp";

	    if ($cnt{$org}-- > 0)
	    {
		print $sim1 or die "reduce_sims cannot write to stdout: $!";
	    }
	}
    }
}

#
# Return the next line from standard in that looks like
# a tab-separated set of data.
#
sub get_input {
    my($line);

    $line = <STDIN>;
    while (defined($line) && ($line !~ /^\S+\t\S+/))
    {
	$line = <STDIN>;
    }
    return $line;
}
