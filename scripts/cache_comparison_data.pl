########################################################################
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

use FIG_Config;
use FIG;
use Carp;
my $fig = new FIG;

my($ref,@others);
my $usage = "usage: cache_comparison_data Ref Other1 Other2 ... OtherN";

(
 ($ref         = shift @ARGV) &&
 (@others      = @ARGV)
)
    || die $usage;

my($ref_dir,$refG,$ref_typeG) = split(/::/,$ref);
#my $ref_cache_dir = "$ref_dir/GenomeComparisonCache";
my $ref_cache_dir  = &GenomeComparisonCache($refG);

&FIG::verify_dir($ref_cache_dir);

if (! -s "$ref_cache_dir/reference_genome")
{
    my $ref_tab = &load_tbl($fig,$ref_dir,$ref_typeG);
    &write_tab($ref_tab,"$ref_cache_dir/reference_genome");
}

foreach my $other (@others)
{
#   print STDERR "processing other genome: $other\n";
    my($other_dir,$otherG,$other_typeG) = split(/::/,$other);
#   my $other_cache_dir = "$other_dir/GenomeComparisonCache";
    my $other_cache_dir = &GenomeComparisonCache($otherG);

    if (! -s "$other_cache_dir/reference_genome")
    {
	my $other_tab = &load_tbl($fig,$other_dir,$other_typeG);
	&FIG::verify_dir($other_cache_dir);
	&write_tab($other_tab,"$other_cache_dir/reference_genome");
    }
    if ((! &current_request($ref_cache_dir,$otherG)) && (! -s "$ref_cache_dir/$otherG"))
    {
	&record_request($ref_cache_dir,$otherG);
	&record_request($other_cache_dir,$refG);
	&process_sims($ref_dir,$refG,$ref_typeG,$other_dir,$otherG,$other_typeG);
	&clear_request($ref_cache_dir,$refG,$otherG);
	&clear_request($other_cache_dir,$otherG,$refG);
    }
}

sub current_request {
    my($dir,$G) = @_;

    my $age;
    # if (age of request < four hours)
    return (($age = -M "$dir/Requests/$G") && ($age < (4 / 24)));
}

sub record_request {
    my($dir,$G) = @_;

    &FIG::verify_dir("$dir/Requests");
    system "touch $dir/Requests/$G";
}

sub clear_request {
    my($dir,$G) = @_;

    if (-e "$dir/Requests/$G")
    {
	unlink("$dir/Requests/$G");
    }
}


sub process_sims {
    my($dir1,$G1,$typeG1,$dir2,$G2,$typeG2) = @_;

    &process_new_sims($dir1,$G1,$dir2,$G2);
}

sub write_tab {
    my($tab,$file) = @_;
    open(FILE,">$file") || confess "could not open $file";
    foreach my $x (@$tab)
    {
	print FILE join("\t",@$x),"\n";
    }
    close(FILE);
}

sub load_func_hash {
    my(@files) = @_;

    my $funcs = {};
    foreach my $file (@files)
    {
	if (open(FILE,"<$file"))
	{
	    while (defined($_ = <FILE>))
	    {
		chomp;
		my($peg,$func) = split(/\t/,$_);
		$funcs->{$peg} = $func;
	    }
	    close(FILE);
	}
    }
    return $funcs;
}

sub load_tbl {
    my($fig,$dir,$type) = @_;

    my $file = "$dir/Features/peg/tbl";
    my $funcs;
    if (-f "$dir/similarities") {
      # this is a RAST org
      $funcs = &load_func_hash("$dir/assigned_functions",
			       "$dir/proposed_non_ff_functions",
			       "$dir/proposed_functions",
			       "$dir/proposed_user_functions");
    } else {
      # this is a SEED org
      $funcs = &load_func_hash("$dir/proposed_non_ff_functions",
			       "$dir/proposed_functions",
			       "$dir/assigned_functions");
    }

    my %tab;
    if (open(TBL,"<$file"))
    {
	while ($_ = <TBL>)
	{
	    if ($_ =~ /^(\S+)\t(\S+)/)
	    {
		my($peg,$loc) = ($1,$2);
		
		my($contig,$beg,$end) = &FIG::boundaries_of($loc);
		if (($type ne "seed") || (! $fig->is_deleted_fid($peg)))
		{
		    $tab{$peg} = [$peg,$contig,$beg,$end,$beg+$end];
	        }
	    }
	}
        close(TBL);

	my @tabL = map { $tab{$_} } keys(%tab);
    
	my $stab = [];
	my $pegI = 1;
	my $contigI = 0;
	my $contigL  = '';
	foreach my $tuple (sort { ($a->[1] cmp $b->[1]) or ($a->[4] <=> $b->[4]) } @tabL)
	{
	    if ($tuple->[1] ne $contigL) { $contigI++; $contigL = $tuple->[1] }
#           reference_genome: [pegI,peg,contigI,contig,beg,end,function]
	    push(@$stab,   [$pegI,
			    $tuple->[0],
			    $contigI,
			    $tuple->[1],
			    $tuple->[2],
			    $tuple->[3],
			    $funcs->{$tuple->[0]} || ''
			   ]);
	    $pegI++;
	}
	return $stab;
    }
    return undef;
}

sub process_new_sims {
    my($dir1,$G1,$dir2,$G2) = @_;

    my $tmpF = "$FIG_Config::temp/sims.$$";

    system "$FIG_Config::bin/sims_between $dir1/Features/peg/fasta $dir2/Features/peg/fasta P 10 0.1 > $tmpF";
    &process_sims_in_file($dir1,$G1,$dir2,$G2,$tmpF);
    unlink($tmpF);
}

sub process_sims_in_file {
    my($dir1,$G1,$dir2,$G2,$sims) = @_;

    my(%sims1,%sims2);
    open(SIMS,"<$sims") || confess "could not open $sims";
    while (defined($_ = <SIMS>))
    {
	if ($_ =~ /^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)/)
	{
	    my $id1 = $1;
	    my $id2 = $2;
	    my $iden = $3 / 100;
	    my $cov1 = sprintf("%0.3f",(($5 - $4)+1) / $9);
	    my $cov2 = sprintf("%0.3f",(($7 - $6)+1) / $10);
	    push(@{$sims1{$id1}},[$iden,$id2,$cov1,$cov2,$9,$10]);
	    push(@{$sims2{$id2}},[$iden,$id1,$cov2,$cov1,$10,$9]);
	}
    }
    close(SIMS);

    my $id1;
    foreach $id1 (keys(%sims1))
    {
	$sims1{$id1} = [sort { ($b->[0] <=> $a->[0]) or ($b->[2] <=> $a->[2]) or ($a->[1] cmp $b->[1]) } @{$sims1{$id1}}];
    }

    foreach $id1 (keys(%sims2))
    {
	$sims2{$id1} = [sort { ($b->[0] <=> $a->[0]) or ($b->[2] <=> $a->[2]) or ($a->[1] cmp $b->[1]) } @{$sims2{$id1}}];
    }
    &process_one_way($dir1,$G1,$dir2,$G2,\%sims1,\%sims2);
    &process_one_way($dir2,$G2,$dir1,$G1,\%sims2,\%sims1);
}


sub process_one_way {
    my($dir1,$g1,$dir2,$g2,$sims1,$sims2) = @_;

    my $cache_dir1 = &GenomeComparisonCache($g1);
    my $cache_dir2 = &GenomeComparisonCache($g2);

    open(REF2,"<$cache_dir2/reference_genome") 
	|| confess "could not open $cache_dir2/reference_genome";

#   reference_genome: [pegI,peg,contigI,contig,beg,end,function]
    my %refH2;
    while (defined($_ = <REF2>))
    {
	chomp;
	my $flds = [split(/\t/,$_)];
	$refH2{$flds->[1]} = $flds;
    }
    close(REF2);

    open(REF,"<$cache_dir1/reference_genome") 
	|| confess "could not open $cache_dir1/reference_genome";
    open(TO,">$cache_dir1/$g2")
	|| confess "could not open $cache_dir1/$g2";

    while (defined($_ = <REF>))
    {
	chomp;
	my($peg1I,$peg1,$contig1I,$contig1,$beg1,$end1,$func1) = split(/\t/,$_);

        my($s1,$s2);
	if ($s1 = $sims1->{$peg1})
	{
#           sims entry: [iden,id2,cov1,cov2,ln1,ln2]
	    my $peg2 = $s1->[0]->[1];
	    my $ref2 = $refH2{$peg2};
	    if (! $ref2) 
	    {
		print TO join("\t",($peg1I,$peg1,"-",'','','','','')),"\n";
	    }
	    else
	    {
		my($peg2I,undef,$contig2I,$contig2,$beg2,$end2,$func2) = @$ref2;

		my $mouse2;
		if (($s1->[0]->[2] < 0.8) || ($s1->[0]->[3] < 0.8))
		{
		    $mouse2 = join("<br>",("location: $contig2 $beg2 $end2",
					   "length: " . $s1->[0]->[5],
					   "identity: " . $s1->[0]->[0],
					   "coverage[ref]: " . $s1->[0]->[2],
					   "coverage[other]: " . $s1->[0]->[3],
					   "function: " . &html_escape($func2)
					   ));
		}
		else
		{
		    $mouse2 = join("<br>",("location: $contig2 $beg2 $end2",
					   "length: " . $s1->[0]->[5],
					   "identity: " . $s1->[0]->[0],
					   "function: " . &html_escape($func2)
					   ));
		}
		my $type = "->";
		if ($s2 = $sims2->{$peg2})
		{
		    if ($s2->[0]->[1] eq $peg1)
		    {
			$type = "<->";
		    }
		}
		print TO join("\t",($peg1I,$peg1,$type,$contig2I,$peg2I,$peg2,$s1->[0]->[0],$mouse2)),"\n";
	    }
	}
	else
	{
	    print TO join("\t",($peg1I,$peg1,"-",'','','','','')),"\n";
	}
    }
    close(REF);
    close(TO);
}

#-------------------------------------------------------------------------------
#  Escape special characters in text for use as html:
#
#     $html = html_escape( $text )
#
#-------------------------------------------------------------------------------
sub html_escape { local $_ = $_[0]; s/\&/&amp;/g; s/>/&gt;/g; s/</&lt;/g; $_ }

sub GenomeComparisonCache {
    my($genome) = @_;

    my $compDir       = $FIG_Config::GenomeComparisonCache ? $FIG_Config::GenomeComparisonCache : $FIG_Config::temp."/GenomeComparisonCache";
    my $ref_cache_dir = "$compDir/$genome";
    return $ref_cache_dir;
}
