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

use Carp;
use FIG_Config;

#&make_sure("formatdb","blastall");

$usage = "usage: sims_between Fasta1 Fasta2 [P|D] IdentityThreshhold CoverageThreshhold";

(
 ($fasta1       = shift @ARGV) &&
 ($fasta2       = shift @ARGV) &&
 ($type         = shift @ARGV) &&
 ($iden_cutoff  = shift @ARGV) &&
 ($cov_cutoff   = shift @ARGV)
)
    || die $usage;
my $length_of = {};
my $length_of1 = &load($fasta1);
my $length_of2 = &load($fasta2);

if ($type =~ /^p$/i)
{
    if ((!-s "$fasta2.pin") || ((-M $fasta2) < (-M "$fasta2.pin")))
    {
	&run("$FIG_Config::ext_bin/formatdb -i $fasta2 -p T");
    }
    
    open(BLAST,"-|", "$FIG_Config::ext_bin/blastall", "-i", $fasta1, "-d", $fasta2, "-e", 
	    "1.0e-5", "-b", "5", "-v", "5", "-FF", "-m", "8", "-p", "blastp") || die "could not open blast: $!";
    my $blast_out;
    while (defined($blast_out = <BLAST>))
    {
	if ($blast_out =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+\t){3}(\S+)\t(\S+)\t(\S+)\t(\S+)\s+(\S+)\s+(\S+)/)
	{
	    $_ = [$1,$2,$3,$5,$6,$7,$8,$9,$length_of1->{$1},$length_of2->{$2},$10];
	    if (($_->[0] ne $_->[1]) &&
		($_->[2] >= $iden_cutoff) &&
                    ((($_->[4] - $_->[3]) / $_->[8]) >= $cov_cutoff) &&
                    ((($_->[6] - $_->[5]) / $_->[9]) >= $cov_cutoff))
	    {
		push(@sims,$_);
	    }
	}
    }
    close(BLAST);
}
else
{
    &run("$FIG_Config::ext_bin/formatdb -i $fasta2 -p F");

    open(BLAST,"$FIG_Config::ext_bin/blastall -i $fasta1 -d $fasta2 -FF -m 8 -p blastn |") || die "could not open blast";
    my $blast_out;
    while (defined($blast_out = <BLAST>))
    {
	if ($blast_out =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+\t){3}(\S+)\t(\S+)\t(\S+)\t(\S+)\s+(\S+)\s+(\S+)/)
	{
	    $_ = [$1,$2,$3,$5,$6,$7,$8,$9,$length_of1->{$1},$length_of2->{$2},$10];
	    if (($_->[0] ne $_->[1]) &&
		($_->[2] >= $iden_cutoff) &&
                    ((($_->[4] - $_->[3]) / $_->[8]) >= $cov_cutoff) &&
                    ((($_->[6] - $_->[5]) / $_->[9]) >= $cov_cutoff))
	    {
		push(@sims,$_);
	    }
	}
    }
    close(BLAST);
}

$last = "";
foreach $sim (@sims)
{
    if ("$sim->[0]\t$sim->[1]" ne $last)
    {
	print join("\t",@$sim), "\n";
	$last = "$sim->[0]\t$sim->[1]";
    }
}


sub run {
    my($cmd) = @_;

#   my @tmp = `date`; chomp @tmp; print STDERR "$tmp[0]: running $cmd\n";
    (system($cmd) == 0) || confess "FAILED: $cmd";
}

sub make_sure {
    my(@progs) = @_;

    my $prog;
    foreach $prog (@progs)
    {
        my @tmp = `which $prog`;
        if ($tmp[0] =~ /^no $prog/)
        {
            print STDERR $tmp[0];
            exit(1);
        }
    }
}

sub load {
    my($fasta) = @_;
    my($id,$seq);
    my $length_of = {};

    $/ = "\n>";
    open(FASTA,"<$fasta") || die "could not open $fasta";
    while (defined($_ = <FASTA>))
    {
        chomp;
        if ($_ =~ /^>?(\S+)[^\n]*\n(.*)/s)
        {
            $id  =  $1;
            $seq =  $2;
            $seq =~ s/\s//gs;
            $length_of->{$id} = length($seq);
        }
    }
    close(FASTA);
    $/ = "\n";
    return $length_of;
}
