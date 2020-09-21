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

package SameFunc;

use Carp;
use Data::Dumper;
use FIG;

sub same_func {
    my($f1,$f2,$strong) = @_;
    my(@ecs1,@ecs2,$ec,$f1n,$f2n);

    my @test1 = split(/\s*[\/;]\s+/,$f1);
    my @test2 = split(/\s*[\/;]\s+/,$f2);
#    if (@test1 != @test2) { return 0 }

    my $atomic_units1 = &atomic_units(&reduce_to_common($f1));
    my $atomic_units2 = &atomic_units(&reduce_to_common($f2));

    if (&hypo($f1) )
    {
	return &hypo($f2);
    }
    elsif (&hypo($f2))
    {
	return 0;
    }

    if ($strong)
    {
	return &in($atomic_units1,$atomic_units2) && &in($atomic_units2,$atomic_units1);
    }

    if (&in($atomic_units1,$atomic_units2) || &in($atomic_units2,$atomic_units1))
    {
	return 1;
    }

    @ecs1 = ($f1 =~ /\d+\.\d+\.\d+\.\d+/g);
    foreach $ec (@ecs1)
    {
	if ($f2 !~ /\b$ec\b/)
	{
	    return 0;
	}
    }

    @ecs2 = ($f2 =~ /\d+\.\d+\.\d+\.\d+/g);
    foreach $ec (@ecs2)
    {
	if ($f1 !~ /\b$ec\b/)
	{
	    return 0;
	}
    }
    
    if (@ecs1) { return 1; }
    return 0;
}

sub same_func_why {
    my($f1,$f2,$strong) = @_;
    my(@ecs1,@ecs2,$ec,$f1n,$f2n);

#    my @test1 = split(/\s*[\/;]\s+/,$f1);
#    my @test2 = split(/\s*[\/;]\s+/,$f2);
#    if (@test1 != @test2) { return 0, 0 }

    my $atomic_units1 = &atomic_units(&reduce_to_common_why($f1));
    my $atomic_units2 = &atomic_units(&reduce_to_common_why($f2));


    if ($strong)
    {
	    ($rc, $why) = ec_test($f1, $f2);
	    if ($rc != -1) {
		return($rc, $why);
	    }
	    return &in($atomic_units1,$atomic_units2) && &in($atomic_units2,$atomic_units1), 3;
    }

    if (&in($atomic_units1,$atomic_units2) || &in($atomic_units2,$atomic_units1))
    {
	return 1, 4;
    }

    ($rc, $why) = ec_test($f1, $f2, 1);
    if ($rc != -1) {
	return($rc, $why);
    }

#    @ecs1 = ($f1 =~ /\d+\.\d+\.\d+\.\d+/g);
#    foreach $ec (@ecs1)
#    {
#	#print STDERR "EC1", $ec;
#	if ($f2 !~ /\b$ec\b/)
#	{
#	    return 0, 5;
#	}
#    }
#
#
#    @ecs2 = ($f2 =~ /\d+\.\d+\.\d+\.\d+/g);
#    foreach $ec (@ecs2)
#    {
#	#print STDERR "EC2", $ec;
#	if ($f1 !~ /\b$ec\b/)
#	{
#	    return 0, 6;
#	}
#    }
#    
#    if (@ecs1) { return 1, 8; }
    if (&hypo($f1) )
    {
	return &hypo($f2), 1;
    }
    elsif (&hypo($f2))
    {
	#print STDERR Dumper($f1, $f2);
	return 0, 2;
    }
    return 0, 7;
}

sub ec_test {
    my ($f1, $f2, $weak) = @_;
    @ecs1 = ($f1 =~ /\d+\.\d+\.\d+\.\d+/g);
    foreach $ec (@ecs1)
    {
	#print STDERR "EC1", $ec;
	if ($weak) {

		if ($f2 =~ /\b$ec\b/)
		{
		    return 1, 8;
		}
	} else {
		if ($f2 !~ /\b$ec\b/)
		{
		    return 0, 5;
		}
	}
    }

    @ecs2 = ($f2 =~ /\d+\.\d+\.\d+\.\d+/g);
    foreach $ec (@ecs2)
    {
	#print STDERR "EC2", $ec;
	if ($f1 !~ /\b$ec\b/)
	{
	    return 0, 6;
	}
    }
    
    if (@ecs1) { return 1, 8; }
    return (-1, -1);
}


sub atomic_units {
    my($x) = @_;

    my @units = split(/\s+/,$x);
    @units = grep {
	              ($_ !~ /^\d{1,2}/) &&
		      ($_ !~ /^i+$/i)
		  } @units;

    return $x ? [@units] : [];
}

sub in {
    my($atomic_units1,$atomic_units2) = @_;
    my $word;

    if (@$atomic_units1 == 0) { return 0 }  # empty list is not considered in anything
    my %atomic_units2 = map { $_ => 1 } @$atomic_units2;
    foreach $word (@$atomic_units1)
    {
	if (! $atomic_units2{$word}) { return 0 }
    }
    return 1;
}

sub reduce_to_common {
    my($func) = @_;

    $func =  lc($func);
    $func =~ s/ (subunit|chain)$//gi;
    $func =~ s/(integral|protein|processing|enzyme|plasmid)//gi;
    $func =~ s/(probable|putative|precursor|homolog|imported|prime)//gi;
    $func =~ s/(subunit|chain) \S+//gi;
    $func =~ s/ABC transporter/transporter/gi;
    $func =~ s/NADH-quinone oxidoreductase/NADH dehydrogenase I/gi;
    $func =~ s/XAA pro/proline/gi;
    $func =~ s/[\/,\.\(\)\[\]\'\-]/ /g;
    $func =~ s/\(EC\s+\d\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+\)//gi;
    $func =~ s/\bNA(\+?)\b/sodium/gi;
    $func =~ s/\bPI\b/phosphate/gi;
    $func =~ s/\d+ kda? SUBUNIT//gi;
    $func =~ s/^.*transcription.*regulat.*$/regulatory/gi;
    $func =~ s/adenylyltransferase/adenylate/gi;
    $func =~ s/biosynthesis/biosynthetic/gi;
    $func =~ s/cotransporter/transporter/gi;
    $func =~ s/diphosphate/pyrophosphate/gi;
    $func =~ s/exopolysacchride/exopolysaccharide/gi;
    $func =~ s/glycosyl transferase/glycosyltransferase/gi;
    $func =~ s/heavy metal associated/heavy metal binding/gi;
    $func =~ s/major facilitator family/MFS/gi;
    $func =~ s/metallopeptidases/metalloprotease/gi;
    $func =~ s/octaprenyltransferase/polyprenyltransferase/gi;
    $func =~ s/oligopeptide/peptide/gi;
    $func =~ s/oxygen independent/anaerobic/gi;
    $func =~ s/permease/transporter/gi;
    $func =~ s/regulation/regulatory/gi;
    $func =~ s/regulatory/regulator/gi;
    $func =~ s/regulon repressor/transcriptional regulator/gi;
    $func =~ s/sensor histidine kinase/sensory transduction protein kinase/gi;
    $func =~ s/sensory box histidine kinase/sensory transduction protein kinase/gi;
    $func =~ s/single-stranded/single-strand/i;
    $func =~ s/site-specific recombinase, phage integrase family/DNA integration\/recombination\/invertion protein/gi;
    $func =~ s/site-specific recombinase, resolvase family/DNA integration\/recombination\/invertion protein/gi;
    $func =~ s/surfeit locus protein 1/SRF1/gi;
    $func =~ s/synthetase/synthase/gi;
    $func =~ s/tetratricopeptide repeat/TPR repeat/gi;
    $func =~ s/transcriptional activator/transcriptional regulator/gi;
    $func =~ s/translocating/transporting/gi;
    $func =~ s/transporter/transport/gi;
    $func =~ s/two component system histidine kinase/sensory transduction protein kinase/gi;
    $func =~ s/two-component response regulator/dna-binding response regulator/gi;
    $func =~ s/type I secretion/efflux/gi;
    $func =~ s/^\s+//;
    $func =~ s/\s+$//;
    $func =~ s/\s{2,100}/ /g;

    return $func;
}
sub reduce_to_common_why {
    my($func) = @_;
    #print STDERR "Before $func\n";
    if ($func =~ m/ribosomal\s+protein\s+([^p]+)p?/io) {
	return ("ribosomal protein $1");
    }
    $func =  lc($func);
    $func =~ s/\(\w+://;
    $func =~ s/family//gi;
    $func =~ s/ (subunit|chain|family)$//gi;
    $func =~ s/(integral|protein|processing|enzyme|plasmid)//gi;
    $func =~ s/(probable|putative|precursor|homolog|imported|prime)//gi;
    $func =~ s/(subunit|chain) \S+//gi;
    $func =~ s/ABC transporter/transporter/gi;
    $func =~ s/NADH-quinone oxidoreductase/NADH dehydrogenase I/gi;
    $func =~ s/XAA pro/proline/gi;
    $func =~ s/[\/,\.\(\)\[\]\'\-]/ /g;
    $func =~ s/\(EC\s+\d\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+\)//gi;
    $func =~ s/\bNA(\+?)\b/sodium/gi;
    $func =~ s/\bPI\b/phosphate/gi;
    $func =~ s/\d+ kda? SUBUNIT//gi;
    $func =~ s/^.*transcription.*regulat.*$/regulatory/gi;
    $func =~ s/adenylyltransferase/adenylate/gi;
    $func =~ s/biosynthesis/biosynthetic/gi;
    $func =~ s/cotransporter/transporter/gi;
    $func =~ s/diphosphate/pyrophosphate/gi;
    $func =~ s/exopolysacchride/exopolysaccharide/gi;
    $func =~ s/glycosyl transferase/glycosyltransferase/gi;
    $func =~ s/heavy metal associated/heavy metal binding/gi;
    $func =~ s/major facilitator family/MFS/gi;
    $func =~ s/metallopeptidases/metalloprotease/gi;
    $func =~ s/octaprenyltransferase/polyprenyltransferase/gi;
    $func =~ s/oligopeptide/peptide/gi;
    $func =~ s/oxygen independent/anaerobic/gi;
    $func =~ s/permease/transporter/gi;
    $func =~ s/regulation/regulatory/gi;
    $func =~ s/regulatory/regulator/gi;
    $func =~ s/regulon repressor/transcriptional regulator/gi;
    $func =~ s/sensor histidine kinase/sensory transduction protein kinase/gi;
    $func =~ s/sensory box histidine kinase/sensory transduction protein kinase/gi;
    $func =~ s/single-stranded/single-strand/i;
    $func =~ s/site-specific recombinase, phage integrase family/DNA integration\/recombination\/invertion protein/gi;
    $func =~ s/site-specific recombinase, resolvase family/DNA integration\/recombination\/invertion protein/gi;
    $func =~ s/surfeit locus protein 1/SRF1/gi;
    $func =~ s/synthetase/synthase/gi;
    $func =~ s/tetratricopeptide repeat/TPR repeat/gi;
    $func =~ s/transcriptional activator/transcriptional regulator/gi;
    $func =~ s/translocating/transporting/gi;
    $func =~ s/transporter/transport/gi;
    $func =~ s/two component system histidine kinase/sensory transduction protein kinase/gi;
    $func =~ s/two-component response regulator/dna-binding response regulator/gi;
    $func =~ s/type I secretion/efflux/gi;
    $func =~ s/  */ /g;
    $func =~ s/^\s+//;
    $func =~ s/\s+$//;
    $func =~ s/\s{2,100}/ /g;

    #print "after $func\n";
    return $func;
}

sub hypo {
    my $x = (@_ == 1) ? $_[0] : $_[1];

    return &FIG::hypo($x);
}

sub group_funcs {
    my(@functions) = @_;
    my(@groups,$func,$i,%seen,$func1);

    @groups = ();
    while ($func = shift @functions)
    {
	if (! $seen{$func})
	{
	    $group = [$func];
	    $seen{$func} = 1;
	    for ($i=0; ($i < @$group); $i++)
	    {
		foreach $func1 (@functions)
		{
		    if ((! $seen{$func1}) && (&same_func($func1,$group->[$i])))
		    {
			push(@$group,$func1);
			$seen{$func1} = 1;
		    }
		}
	    }
	    push(@groups,$group);
	}
    }
    return @groups;
}

1;
