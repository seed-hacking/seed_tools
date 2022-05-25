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
use FileHandle;
use Data::Dumper;

my $usage = "usage: pull_fasta_entries [-v] [-chunk filesize basename] FastaFile < entries_desired";

my $invert  = "";
my $trouble = 0;
my $filesize;
my $basename;
my @fasta;
while (@ARGV)
{
    my $arg = shift;
    
    if (-s $arg) {
	push(@fasta, $arg);
    }
    elsif ($arg =~ m/-v/) {
	$invert = $arg;
    }
    elsif ($arg eq '-chunk') {
	$filesize = shift;
	$basename = shift;
	if ($filesize !~ /^\d+$/ or !defined($basename)) {
	    die "Usage: $usage\n";
	}
    }
    else {
	$trouble = 1;
	print STDERR "Invalid arg $arg\n";
    }
}

die "$usage\n" if $trouble;

#my %to_pull = map { chomp; s/^([^|]+\|[^|]+)\|.*$/$1/; s/\s*$//; $_ => 1 } <STDIN>;
my %to_pull;
while (<STDIN>)
{
    if (/(\S+)/)
    {
	$to_pull{$1} = 1;
    }
}

my $next_idx = 1;
my $out_fh;

if (defined($basename))
{
    $out_fh = open_next_output();
}
else
{
    $out_fh = \*STDOUT;
}
    

$/ = "\n>";
for my $fasta (@fasta)
{
    open(FASTA,"<$fasta") || die "could not open $fasta: $!";
    while (defined($_ = <FASTA>))
    {
	chomp;
	if ($_ =~ /^>?(\S+)[^\n]*\n(.*)/s)
	{
	    my $id  =  $1;
	    my $seq =  $2;
	    my $id1 = $id;
	    #$id1 =~ s/^([^|]+\|[^|]+)\|.*$/$1/;
	    $seq =~ s/\s//gs;
	    if ($to_pull{$id1} xor $invert)
	    {
		&display_id_and_seq($id,\$seq, $out_fh);
		if (defined($basename) and tell($out_fh) > $filesize)
		{
		    close($out_fh);
		    $out_fh = open_next_output();
		}
	    }
	}
    }
    close(FASTA);
}

sub open_next_output
{
    my $fh;
    my $fn = sprintf("$basename.%02d", $next_idx++);
    $fh = new FileHandle(">$fn") or die "cannot open $fn: $!\n";
    return $fh;
}

sub display_id_and_seq {
    my( $id, $seq, $fh ) = @_;
    
    if (! defined($fh) )  { $fh = \*STDOUT; }
    
    print $fh ">$id\n";
    &display_seq($seq, $fh);
}

sub display_seq {
    my ( $seq, $fh ) = @_;
    my ( $i, $n, $ln );
    
    if (! defined($fh) )  { $fh = \*STDOUT; }

    $n = length($$seq);
#   confess "zero-length sequence ???" if ( (! defined($n)) || ($n == 0) );
    for ($i=0; ($i < $n); $i += 60)
    {
	if (($i + 60) <= $n)
	{
	    $ln = substr($$seq,$i,60);
	}
	else
	{
	    $ln = substr($$seq,$i,($n-$i));
	}
	print $fh "$ln\n";
    }
}
