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

package Blast;

use Carp;
use Data::Dumper;

use strict;

sub blastp {
    my($id_seqs,$db,$parms,$keep_ali,$keep_def) = @_;
    my($entry,$sim,$line,$id,$seq,$sims);
    my($id1,$ln1,$id2,$ln2,$score,$exp,$n_match,$n_ident,$b1,$e1,$seq1,$b2,$e2,$seq2);
    my($def2,$ali);

    my $old_eol = $/;
    $/ = "\n";

    if (! $ENV{"BLASTMAT"}) { $ENV{"BLASTMAT"} = "$FIG_Config::blastmat"; }

    if ($ENV{"PATH"} !~ /fig\/bin/)
    {
	$ENV{"PATH"} = "$FIG_Config::bin:" . $ENV{"PATH"};
    }

    my $full_db = (-s $db) ? $db : $ENV{"BLASTDB"} . $db;
    if  (! -s $full_db) 
    {
	print STDERR "could not locate $db\n";
	$/ = $old_eol;
	return undef;
    }

    if ((! -s "$full_db.psq") || (-M $full_db < -M "$full_db.psq"))
    {
	system "$FIG_Config::ext_bin/formatdb -i $full_db -p T";
    }

    open(TMP,">tmp$$.fasta") || die "could not open tmp$$.fasta";
    $sims = {};
    foreach $entry (@$id_seqs)
    {
	($id,$seq) = @$entry;
	$sims->{$id} = [];
	print TMP ">$id\n$seq\n";
#	print STDERR "id=$id ",length($seq),"\n";
    }
    close(TMP);

    my $count = 3;
    while ($count && (system("$FIG_Config::ext_bin/blastall -i tmp$$.fasta -d \"$db\" -p blastp $parms | rationalize_blast > tmp$$.blastout") != 0))
    {
	open(LOG,">>$FIG_Config::data/blast.log") || die "could not open $FIG_Config::data/blast.log";
	print LOG "============\n";
	open(TMP,"<tmp$$.fasta") || die "could not open tmp$$.fasta";
	while (defined($line = <TMP>))
	{
	    print LOG $line;
	}
	close(TMP);
	print LOG "blastall -i tmp$$.fasta -d \"$db\" -p blastp $parms | rationalize_blast > tmp$$.blastout\n";
	close(LOG);

	print STDERR "blast failed\n";
	$count--;
    }

    if ($count == 0) { return $sims }

    open(TMP,"<tmp$$.blastout")
	|| die "could not open tmp$$.blastout";

    while (defined($line = <TMP>))
    {
	chomp $line;
	if ($line =~ /^Query=\t(\S+)\t(.*)\t(\d+)/)
	{
	    $id1 = $1;
	    $ln1 = $3;
	    if (! $sims->{$id1}) { $sims->{$id1} = [] }
	}
	elsif ($line =~ /^>\t(\S+)\t(.*)\t(\d+)/)
	{
	    $id2  = $1;
	    $ln2  = $3;
	    $def2 = ($2 && $keep_def) ? $2 : undef;
	}
	elsif ($line =~ /^HSP/)
	{
#           HSP  score  exp  p_n  p_val  n_match  n_ident  n_sim  n_gap  dir  q1  q2  q_sq  s1  s2  s_sq
	    (undef,$score,$exp,undef,undef,$n_match,$n_ident,undef,undef,undef,$b1,$e1,$seq1,$b2,$e2,$seq2) = split(/\t/,$line);
	    $ali = $keep_ali ? [$seq1,$seq2] : undef;
	    if ($n_match)
	    {
		$sim = [$id1,
			$id2,
			int(($n_ident * 100)/$n_match),
			length($seq1),
			$n_match - $n_ident,
			undef,
			$b1,
			$e1,
			$b2,
			$e2,
			$exp,
			$score,
			$ln1,
			$ln2,
			"blastp",
			$def2,
			$ali
		       ];
		bless($sim,"Sim");
		push(@{$sims->{$id1}},$sim);
	    }
	}
    }
    close(TMP);
    unlink("tmp$$.fasta");
    unlink("tmp$$.blastout");
    $/ = $old_eol;
    return $sims;
}

sub blastn {
    my($id_seqs,$db,$parms,$keep_ali,$keep_def) = @_;
    my($entry,$sim,$line,$id,$seq,$sims);
    my($id1,$ln1,$id2,$ln2,$score,$exp,$n_match,$n_ident,$b1,$e1,$seq1,$b2,$e2,$seq2);
    my($def2,$ali);

    my $old_eol = $/;
    $/ = "\n";

    if (! $ENV{"BLASTMAT"}) { $ENV{"BLASTMAT"} = "$FIG_Config::fig/BLASTMAT"; }

    if ($ENV{"PATH"} !~ /fig\/bin/)
    {
	$ENV{"PATH"} = "$FIG_Config::bin:" . $ENV{"PATH"};
    }

    my $full_db = (-s $db) ? $db : $ENV{"BLASTDB"} . $db;
    if  (! -s $full_db) 
    {
	print STDERR "could not locate $db\n";
	$/ = $old_eol;
	return undef;
    }

    if ((! -s "$full_db.nsq") || (-M $full_db < -M "$full_db.nsq"))
    {
	system "formatdb -i $full_db -p F";
    }

    open(TMP,">tmp$$.fasta") || die "could not open tmp$$.fasta";
    $sims = {};
    foreach $entry (@$id_seqs)
    {
	($id,$seq) = @$entry;
	$sims->{$id} = [];
	print TMP ">$id\n$seq\n";
#	print STDERR "id=$id ",length($seq),"\n";
    }
    close(TMP);

    my $count = 3;
    while ($count && (system("blastall -i tmp$$.fasta -d \"$db\" -p blastn $parms | rationalize_blast 2> /dev/null > tmp$$.blastout") != 0))
    {
	open(LOG,">>$FIG_Config::data/blast.log") || die "could not open $FIG_Config::data/blast.log";
	print LOG "============\n";
	open(TMP,"<tmp$$.fasta") || die "could not open tmp$$.fasta";
	while (defined($line = <TMP>))
	{
	    print LOG $line;
	}
	close(TMP);
	print LOG "blastall -i tmp$$.fasta -d \"$db\" -p blastn $parms | rationalize_blast > tmp$$.blastout\n";
	close(LOG);

	print STDERR "blast failed\n";
	$count--;
    }

    if ($count == 0) { return $sims }

    open(TMP,"<tmp$$.blastout")
	|| die "could not open tmp$$.blastout";

    while (defined($line = <TMP>))
    {
	chomp $line;
	if ($line =~ /^Query=\t(\S+)\t(.*)\t(\d+)/)
	{
	    $id1 = $1;
	    $ln1 = $3;
	    if (! $sims->{$id1}) { $sims->{$id1} = [] }
	}
	elsif ($line =~ /^>\t(\S+)\t(.*)\t(\d+)/)
	{
	    $id2  = $1;
	    $ln2  = $3;
	    $def2 = ($2 && $keep_def) ? $2 : undef;
	}
	elsif ($line =~ /^HSP/)
	{
#           HSP  score  exp  p_n  p_val  n_match  n_ident  n_sim  n_gap  dir  q1  q2  q_sq  s1  s2  s_sq
	    (undef,$score,$exp,undef,undef,$n_match,$n_ident,undef,undef,undef,$b1,$e1,$seq1,$b2,$e2,$seq2) = split(/\t/,$line);
	    $ali = $keep_ali ? [$seq1,$seq2] : undef;
	    if ($n_match)
	    {
		$sim = [$id1,
			$id2,
			int(($n_ident * 100)/$n_match),
			length($seq1),
			$n_match - $n_ident,
			undef,
			$b1,
			$e1,
			$b2,
			$e2,
			$exp,
			$score,
			$ln1,
			$ln2,
			"blastn",
			$def2,
			$ali
		       ];
		bless($sim,"Sim");
		push(@{$sims->{$id1}},$sim);
	    }
	}
    }
    close(TMP);
    unlink("tmp$$.fasta");
    unlink("tmp$$.blastout");
    $/ = $old_eol;
    return $sims;
}

sub blast2p {
    my($id_seqs,$db,$parms,$keep_ali,$keep_def) = @_;
    my($entry,$sim,$line,$id,$seq,$sims);
    my($id1,$ln1,$id2,$ln2,$score,$exp,$n_match,$n_ident,$b1,$e1,$seq1,$b2,$e2,$seq2);
    my($def2,$ali);

    if (! $ENV{"BLASTMAT"}) { $ENV{"BLASTMAT"} = "$FIG_Config::fig/BLASTMAT"; }

    my $full_db = (-s $db) ? $db : $ENV{"BLASTDB"} . $db;
    if  (! -s $full_db) 
    {
	print STDERR "could not locate $db\n";
	return undef;
    }

    if ((! -s "$full_db.psq") || (-M $full_db < -M "$full_db.psq"))
    {
	system "formatdb -i $full_db -p T";
    }

    my($entry1,$entry2) = @$id_seqs;

    if ((ref($entry1) ne "ARRAY") ||
	(ref($entry2) ne "ARRAY"))
    {
	print STDERR "bad input to blast2p\n";
	return ();
    }

    ($id,$seq) = @$entry1;
    open(TMP,">tmp1$$.fasta") || die "could not open tmp1$$.fasta";
    print TMP ">$id\n$seq\n";
    close(TMP);

    ($id,$seq) = @$entry2;
    open(TMP,">tmp2$$.fasta") || die "could not open tmp2$$.fasta";
    print TMP ">$id\n$seq\n";
    close(TMP);

    open(TMP,"bl2seq -i tmp1$$.fasta -j tmp2$$.fasta -p blastp $parms |tee /dev/stderr | rationalize_blast |")
	|| die "could not open pipe for blasting";

    $sims = {};
    while (defined($line = <TMP>))
    {
	print STDERR $line;

	chomp $line;
	if ($line =~ /^Query=\t(\S+)\t(.*)\t(\d+)/)
	{
	    $id1 = $1;
	    $ln1 = $3;
	    if (! $sims->{$id1}) { $sims->{$id1} = [] }
	}
	elsif ($line =~ /^>\t(\S+)\t(.*)\t(\d+)/)
	{
	    $id2  = $1;
	    $ln2  = $3;
	    $def2 = ($2 && $keep_def) ? $2 : undef;
	}
	elsif ($line =~ /^HSP/)
	{
#           HSP  score  exp  p_n  p_val  n_match  n_ident  n_sim  n_gap  dir  q1  q2  q_sq  s1  s2  s_sq
	    (undef,$score,$exp,undef,undef,$n_match,$n_ident,undef,undef,undef,$b1,$e1,$seq1,$b2,$e2,$seq2) = split(/\t/,$line);
	    $ali = $keep_ali ? [$seq1,$seq2] : undef;
	    if ($n_match)
	    {
		$sim = [$id1,
			$id2,
			int(($n_ident * 100)/$n_match),
			length($seq1),
			$n_match - $n_ident,
			undef,
			$b1,
			$e1,
			$b2,
			$e2,
			$exp,
			$score,
			$ln1,
			$ln2,
			"blastp",
			$def2,
			$ali
		       ];
		bless($sim,"Sim");
		push(@{$sims->{$id1}},$sim);
	    }
	}
    }
    close(TMP);
    return $sims;
}
	
1
