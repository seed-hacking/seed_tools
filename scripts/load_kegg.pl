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


# -*- perl -*-

###########################################
use strict;

# KEY RELATIONAL TABLES:
#    
#   1.  compound(Cid,Priority,Name)    % we have a prioritized set of names for a compound
#   2.  comp_to_CAS(Cid,CASid)         % connection to chemical abstract society [optional]
#   3.  reaction(Rid,Reversible)       % nonreversible go from substrates to products
#   4.  reaction_to_compound(Rid,0/1[substrate or product],Cid,Stoich,main_compound [part of major transformation])
#   5.  reaction_to_role(Rid,FunctionalRole)
#

use FIG;
my $fig = new FIG;

my $usage = "usage: load_kegg";
my $base_path = $FIG_Config::kegg || "$FIG_Config::data/KEGG";

use Tracer;
TSetup('1 *', 'WARN');

if (! -f "$base_path/ligand/enzyme/enzyme")
{
    warn "Skipping KEGG load\n";
    exit;
}


&load_ec_and_map_data;
&load_compounds;
&load_reactions_and_catalyzes;

undef $fig;

sub load_ec_and_map_data {

    Open(\*TMPIN, "<$base_path/ligand/enzyme/enzyme");
    Open(\*ECMAP,">$FIG_Config::temp/ec_map.table");

	Trace("Reading KEGG enzymes.") if T(2);
    my($ec,%name,$map);
    $/ = "\n///\n";
    while (defined($_ = <TMPIN>))
    {
		if ($_ =~ /ENTRY\s+EC\s+(\d+\.\d+\.\d+\.\d+)/s)
		{
			$ec = $1;
			while ($_ =~ /PATH:\s+ec(\d+)\s+(\S[^\n]+\S)/sg)
			{
			    # mdj: prefix of 'map' added for backwards compatibility
			    my $map_name = "map".$1;
			    print ECMAP "$ec\t$map_name\n";
			    $name{$map_name} = $2;
			}
		}
    }
    $/ = "\n";
    close(TMPIN);
    close(ECMAP);
	
	Trace("Writing map table.") if T(2);
    Open(\*MAP, ">$FIG_Config::temp/map_name.table");

    foreach $map (keys(%name))
    {
		print MAP "$map\t$name{$map}\n";
    }
    close(MAP);

    $fig->reload_table('all', "ec_map",
					   "ec varchar(100), map varchar(100)",
					   { index_ec_map_ec => "ec", index_ec_map_map => "map" },
					   "$FIG_Config::temp/ec_map.table");
#    unlink("$FIG_Config::temp/ec_map.table");
	$fig->reload_table('all', "map_name",
					   "map varchar(100) UNIQUE NOT NULL, mapname varchar(200), primary key ( map )",
					   { }, "$FIG_Config::temp/map_name.table");
#    unlink("$FIG_Config::temp/map_name.table");
}

sub load_compounds {
    
	Trace("Loading compounds.") if T(2);
	
    Open(\*TMPIN, "<$base_path/ligand/compound/compound");
    Open(\*COMP, ">$FIG_Config::temp/comp_name.table");
    Open(\*CAS, ">$FIG_Config::temp/comp_cas.table");

    my($cid,$name,$cas,$names,$tmp,$n,$entry);
    $/ = "\n///\n";
    while (defined($entry = <TMPIN>))
    {
		if ($entry =~ /ENTRY\s+(C\d+).*\nNAME\s+(\S[^\n]*)\n((\s+(\S[^\n]*\S)\n)*)/s)
		{
			$cid = $1;
			$names = $2;

			if ($3)
			{
				$tmp = $3;
				chop $tmp;
				$tmp =~ s/^\s+/ /;
				$names = $names . $tmp;
				$names =~ s/\n\s+/ /g;
				$names =~ s/- /-/g;
			}
	
			$n = 1;
			foreach $name (map { $_ =~ s/^\s+//; $_ =~ s/\s+$//; $_ } split(/;/,$names))
			{
				print COMP "$cid\t$n\t$name\n";
				$n++;
				if (length $name > 200) { print "$cid, $name\n" }
			}
		}

		if ($entry =~ /DBLINKS\s+CAS:\s+(\S+)/s)
		{
			print CAS "$cid\t$1\n";
		}
    }
    $/ = "\n";
    close(TMPIN);
    close(COMP);
    close(CAS);

	$fig->reload_table('all', "comp_name",
					   "cid varchar(7), pos integer, name varchar(200)",
					   { index_comp_name_cid => "cid",
						 index_comp_name_name => "name" },
					"$FIG_Config::temp/comp_name.table");
#    unlink("$FIG_Config::temp/comp_name.table");

	$fig->reload_table('all', "comp_cas",
					   "cid varchar(7), cas varchar(100)",
					   { index_comp_cas_cid => "cid",
						 index_comp_cas_cas => "cas" },
					   "$FIG_Config::temp/comp_cas.table");
#    unlink("$FIG_Config::temp/comp_cas.table");
}

sub load_reactions_and_catalyzes {

    my($react,$path,$sub,$prod,@sub,@prod,$subs,$prods,$dir);
    my($cid,$n,$main,$paths,%reaction,$x);

	Trace("Loading reactions.") if T(2);
    Open(\*REAC, "<$base_path/ligand/reaction/reaction");
    Open(\*REACTION, "<$base_path/ligand/reaction/reaction.lst");
    Open(\*RMAPFORMULA, "<$base_path/ligand/reaction/reaction_mapformula.lst");
    Open(\*R2C, ">$FIG_Config::temp/reaction_to_compound.table");
    Open(\*REV, ">$FIG_Config::temp/rev.table");
    Open(\*RDIR, ">$FIG_Config::temp/reaction_direction.table");
    Open(\*REAC2ENZ, ">$FIG_Config::temp/reaction_to_enzyme.table");

	Trace("Reading reaction list file.") if T(2);

    while (defined($_ = <REACTION>))
    {
		if ($_ =~ /(R\d+):\s+(\S.*\S)\s+<=>\s+(\S.*\S)/)
		{
			$react = $1;
			$sub   = $2;
			$prod  = $3;
			@sub = split(/\s+\+\s+/,$sub);
			@prod = split(/\s+\+\s+/,$prod);
			@sub   = map { $_ =~ /^(([\dmn\(]\S*)\s+)?([CG]\d+)/; $2 ? [$3,$2,[]] : [$3,1,[]] } @sub;
			@prod  = map { $_ =~ /^(([\dmn\(]\S*)\s+)?([CG]\d+)/; $2 ? [$3,$2,[]] : [$3,1,[]] } @prod;
			$reaction{$react} = [[@sub],[@prod]];
		}
		else
		{
			Trace("Invalid reaction format: $_") if T(1);
		}
    }
    close(REACTION);

	Trace("Reading main reaction file.") if T(2);

    my %reversibility;

    while (defined($_ = <RMAPFORMULA>))
    {
	# ignore overview map 01100
	next if ($_ =~ /\s01100:/);

	if ($_ =~ /^(R\d+):\s+(\d+):\s+(\S.*\S)\s(\<?=\>?)\s(\S.*\S)/)
	{
	    $react = $1;
	    $path  = $2;
	    $sub   = $3;
	    $dir   = $4;
	    $prod  = $5;

	    if (exists($reaction{$react}))
	    {
		$subs  = $reaction{$react}->[0];
		$prods = $reaction{$react}->[1];
		my $rc = &mark_main($sub,$subs,$path);

		if ($rc == 0)
		{
		    $rc = &mark_main($sub,$prods,$path) && &mark_main($prod,$subs,$path);

		    if ($rc != 0 and $dir ne "<=>")
		    {
			if ($dir eq "=>")
			{
			    $dir = "<=";
			}
			else
			{
			    $dir = "=>";
			}
		    }
		}
		else
		{
		    $rc = &mark_main($prod,$prods, $path);
		}

		if ($rc == 0)
		{
		    print "Can't handle $_\n";
		}

		# entry for the reaction in the context of the pathway
		$reversibility{$react.":".$path} = $dir;

		# since there can be multiple entries per reaction, with different
		# reversibility info (in the context of different pathways),
		# reversible trumps non-reversible for the general reaction entry
		if (($dir eq "<=") || ($dir eq "=>"))
		{
		    if (! defined($reversibility{$react}))
		    {
			$reversibility{$react} = $dir;
		    }
		    elsif ($reversibility{$react} ne "<=>" && $reversibility{$react} ne $dir)
		    {
			$reversibility{$react} = "<=>";
		    }
		}
		else
		{
		    $reversibility{$react} = "<=>";
		}
	    }
	}
    }
    close(RMAPFORMULA);

    my($entry);

    Trace("Reading KEGG reaction file.") if T(2);
    my($rid,$ec,@ecs,$ecs);
    $/ = "\n///\n";
    while (defined($entry = <REAC>))
    {
	if ($entry =~ /ENTRY\s+(R\d+).*\nENZYME\s+(\S[^a-zA-Z\/]+)/s)
	{
	    $rid = $1;
	    $ecs = $2;

	    foreach $ec (split(/\s+/,$ecs))
	    {
		print REAC2ENZ "$rid\t$ec\n";
	    }
	}
    }

    $/ = "\n";
    close(REAC);
    close(REAC2ENZ);

    foreach $react (keys %reversibility)
    {
	my ($rid, $mapid) = split ":", $react;
	
	if (defined $mapid)
	{
	    if ($reversibility{$react} eq "<=>")
	    {
		print RDIR "$rid\t$mapid\tB\n";
	    }
	    elsif ($reversibility{$react} eq "<=")
	    {
		print RDIR "$rid\t$mapid\tL\n";
	    }
	    else
	    {
		print RDIR "$rid\t$mapid\tR\n";
	    }
	}
	else
	{
	    if ($reversibility{$react} eq "<=>")
	    {
		print REV "$rid\t1\n";
	    }
	    else
	    {
		print REV "$rid\t0\n";
	    }
	}
    }

    close(REV);
    close(RDIR);

    Trace("Connecting reactions to compounds.") if T(2);

    foreach $react (sort keys(%reaction))
    {
	($subs,$prods) = @{$reaction{$react}};
	foreach $x (@$subs)
	{
	    ($cid,$n,$paths) = @$x;

	    if (scalar @{$paths} > 0)
	    {
		foreach $path (@{$paths})
		{
		    print R2C "$react\t0\t$cid\t$n\t1\t$path\n";
		}
	    }
	    else
	    {
		print R2C "$react\t0\t$cid\t$n\t0\t\n";
	    }
	}

	foreach $x (@$prods)
	{
	    ($cid,$n,$paths) = @$x;

	    if (scalar @{$paths} > 0)
	    {
		foreach $path (@{$paths})
		{
		    print R2C "$react\t1\t$cid\t$n\t1\t$path\n";
		}
	    }
	    else
	    {
		print R2C "$react\t1\t$cid\t$n\t0\t\n";
	    }
	}
    }
    close(R2C);

    $fig->reload_table('all', "reaction_to_compound", 
		       "rid varchar(8), setn char(1), cid varchar(8), stoich char(6), main char(1), path char(5)",
		       { index_reaction_to_compound_rid => "rid",
			 index_reaction_to_compound_cid => "cid" },
		       "$FIG_Config::temp/reaction_to_compound.table"); 
#    unlink("$FIG_Config::temp/reaction_to_compound.table");

    $fig->reload_table('all', "reversible",
		       "rid varchar(8) UNIQUE NOT NULL, reversible char(1), primary key(rid)",
		       { }, "$FIG_Config::temp/rev.table");
#    unlink("$FIG_Config::temp/rev.table");

    $fig->reload_table('all', "reaction_direction",
		       "rid varchar(8) NOT NULL, mapid varchar(8), direction char(1)",
		       { index_reaction_direction_rid => "rid",
			 index_reaction_direction_mapid => "mapid" }, 
		       "$FIG_Config::temp/reaction_direction.table");
#    unlink("$FIG_Config::temp/reaction_direction.table");

    Trace("Reactions processed.") if T(2);

    $fig->reload_table('all', "reaction_to_enzyme",
		       "rid varchar(8), role varchar(100)",
		       { index_reaction_to_enzyme_rid => "rid",
			 index_reaction_to_enzyme_role => "role" },
		       "$FIG_Config::temp/reaction_to_enzyme.table");
#    unlink("$FIG_Config::temp/reaction_to_enzyme.table");
    Trace("Enzyme reactions loaded.") if T(2);
}

sub mark_main {
    my($main,$set,$path) = @_;
    my($cid,$i);

    foreach $cid (split(/\s+\+\s+/,$main))
    {
		for ($i=0; ($i < @$set) && ($set->[$i]->[0] ne $cid); $i++) {}
		if ($i == @$set)
		{
		    # compound id was not found in the reaction
		    return 0;
		}
		else
		{
		    push @{$set->[$i]->[2]}, $path;
		}
    }

    return 1;
}
