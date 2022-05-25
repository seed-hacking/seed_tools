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

package CompareMR;

use Carp;
use strict;
use FIG_Config;
use Tracer;
use Data::Dumper;

use FIG;
use FIGV;

# "arg1" is evidently either a genome id or a pointer to a RAST directory.
# $genome1 is irrelevant, if $fig_services has not been instantiated.
# If it has, then it is the real 'genome1'.  This is pretty clutzy.

sub compare_genomes_MR {
    my($arg1,$genome2, $fig_services, $genome1) = @_;

    unless ($fig_services) {
      if (-d "$FIG_Config::organisms/$arg1")
	{
	  $fig_services = FIG->new();
	  $genome1      = $arg1;
	}
      else 
	{
	  if (-d $arg1)
	    {
	      my $figV = new FIGV($arg1,FIG->new());
	      $fig_services = $figV;
	      if ($arg1 =~ /(\d+\.\d+)$/)
		{
		  $genome1 = $1;
		}
	      else
		{
		  die "invalid first argument: $arg1";
		}
	    }
	  else
	    {
	      die "Invalid first argument;  $arg1";
	    }
	}
    }

    my $data1 = &data_on_subs_roles_pegs($fig_services,$genome1);
    my $data2 = &data_on_subs_roles_pegs($fig_services,$genome2);

    my $common   = &process_common($fig_services,$genome1,$genome2,$data1,$data2);
    my $in1_not2 = &process_diff($fig_services,$genome1,$genome2,$data1,$data2);
    my $in2_not1 = &process_diff($fig_services,$genome2,$genome1,$data2,$data1);
    return ($common,$in1_not2,$in2_not1);
}

sub process_common {
    my($fig_services,$genome1,$genome2,$data1,$data2) = @_;

    my $common = [];
    my $i1 = 0;
    my $i2 = 0;
    while (($i1 < @$data1) && ($i2 < @$data2))
    {
	if    ( ($data1->[$i1]->[0] lt $data2->[$i2]->[0]) ||
	       (($data1->[$i1]->[0] eq $data2->[$i2]->[0]) && ($data1->[$i1]->[1] lt $data2->[$i2]->[1])))
	{
	    $i1++;
	}
	elsif ( ($data1->[$i1]->[0] gt $data2->[$i2]->[0]) ||
	       (($data1->[$i1]->[0] eq $data2->[$i2]->[0]) && ($data1->[$i1]->[1] gt $data2->[$i2]->[1])))
	{
	    $i2++;
	}
	else
	{
	    my $currS = $data1->[$i1]->[0];
	    my $currR = $data1->[$i1]->[1];
	    my $pegs1 = [];
	    while (($i1 < @$data1) && ($data1->[$i1]->[0] eq $currS) && ($data1->[$i1]->[1] eq $currR))
	    {
		push(@$pegs1,$data1->[$i1]->[2]);
		$i1++;
	    }

	    my $pegs2 = [];
	    while (($i2 < @$data2) && ($data2->[$i2]->[0] eq $currS) && ($data2->[$i2]->[1] eq $currR))
	    {
		push(@$pegs2,$data2->[$i2]->[2]);
		$i2++;
	    }
	    push(@$common,[$currS,$currR,$pegs1,$pegs2]);
	}
    }
    return $common;
}

sub process_diff {
    my($fig_services,$genomeA,$genomeB,$dataA,$dataB) = @_;

    my $in1_not2 = [];
    my $i1 = 0;
    my $i2 = 0;
    while ($i1 < @$dataA)
    {
	if    ( ($i2 == @$dataB) || ($dataA->[$i1]->[0] lt $dataB->[$i2]->[0]) ||
	       (($dataA->[$i1]->[0] eq $dataB->[$i2]->[0]) && ($dataA->[$i1]->[1] lt $dataB->[$i2]->[1])))
	{
	    my $currS = $dataA->[$i1]->[0];
	    my $currR = $dataA->[$i1]->[1];
	    my $pegs1 = [];
	    while (($i1 < @$dataA) && ($dataA->[$i1]->[0] eq $currS) && ($dataA->[$i1]->[1] eq $currR))
	    {
		push(@$pegs1,$dataA->[$i1]->[2]);
		$i1++;
	    }
	    push(@$in1_not2,[$currS,$currR,$pegs1,[$fig_services->seqs_with_role($currR,undef,$genomeB)]]);
	}
	elsif ( ($dataA->[$i1]->[0] gt $dataB->[$i2]->[0]) ||
	       (($dataA->[$i1]->[0] eq $dataB->[$i2]->[0]) && ($dataA->[$i1]->[1] gt $dataB->[$i2]->[1])))
	{
	    $i2++;
	}
	else
	{
	    my $currS = $dataA->[$i1]->[0];
	    my $currR = $dataA->[$i1]->[1];
	    while (($i1 < @$dataA) && ($dataA->[$i1]->[0] eq $currS) && ($dataA->[$i1]->[1] eq $currR))
	    {
		$i1++;
	    }

	    while (($i2 < @$dataB) && ($dataB->[$i2]->[0] eq $currS) && ($dataB->[$i2]->[1] eq $currR))
	    {
		$i2++;
	    }
	}
    }
    return $in1_not2;
}

sub data_on_subs_roles_pegs {
    my($fig_services,$genome) = @_;

    my $sub_data = $fig_services->get_genome_subsystem_data($genome);
    return [sort { ($a->[0] cmp $b->[0]) or ($a->[1] cmp $b->[1]) or ($a->[3] <=> $b->[3]) } map {$_->[2] =~ /\.(\d+)$/;$_->[1] =~ s/\#.+//; [@$_, $1]} @$sub_data];

    # sort by_fig_id inefficient for metagenomes with large numbers of features
#    return [sort { ($a->[0] cmp $b->[0]) or ($a->[1] cmp $b->[1]) or &FIG::by_fig_id($a->[2],$b->[2]) } @$sub_data];
}

1;

