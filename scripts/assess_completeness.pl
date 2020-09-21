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
use Carp;
use Data::Dumper;

use FIG;    my $fig = new FIG;

my $usage = "assess_completeness [-f min_fraction] [-s min_size] genome_id ...";

#  assess_completeness
#
#  Reads contigs files for a genome and writes a file called PROBABLY_COMPLETE
#  to the genome directory if the genome is judged to be complete by the
#  criterion described below.
#
#  A simple but remarkably robust heuristic is used.  The program calculates the
#  fraction of the nucleotides in contigs of greater than or equal to 20,000 bp.
#  Emperically, values of 0.7 or better correspond to effectively complete
#  genomes (the vast majority of genes can be found and the data quality tends
#  to be very good).  The threshold can be set by the -f option.
#

my $minfrac =    0.7;
my $minlen  =  20000;
my $minsize = 300000;

my ( $flag, $val );
while ( @ARGV && ( ( $flag, $val ) = $ARGV[0] =~ /^\-(.)(.*)/ ) )
{
    shift;
    if ( ( $flag eq "f" ) && ( $val ||= shift @ARGV ) )
    {
	$minfrac = $val;
    }
    elsif ( ( $flag eq "l" ) && ( $val ||= shift @ARGV ) )
    {
	$minlen = $val;
    }
    elsif ( ( $flag eq "s" ) && ( $val ||= shift @ARGV ) )
    {
	$minsize = $val;
    }
    else
    {
	print STDERR "Problem with flag '$flag'\n";
	die "Usage: $usage\n";
    }
}

@ARGV and ( $minfrac > 0 ) and ( $minsize > 0 ) or die "Usage: $usage\n";

my $orgdir = "$FIG_Config::organisms";
my $genomedir;

my ( $genome, $ttlen, $inbig, $contig, $len );

foreach $genome ( @ARGV )
{
    print STDERR "." if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} == 1));
    
    if (-d $genome)
    { 
	$genomedir = $genome; 
    } else {
	$genomedir = "$orgdir/$genome";   #...Fall back to installed organisms...
    }
    
    if (opendir( GENOME, $genomedir ) )
    {
	print STDERR "Checking genomedir=$genomedir\n" if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} > 1));
    }
    else
    {
        print STDERR "Skipping $genome; could not open directory $genomedir\n";
        next;
    }

    #The next line requires the org to have already been installed.  But
    #add_genome calls assess_completeness as part of the load!
    #
    #next unless ($fig->is_prokaryotic($genome) || $fig->is_eukaryotic($genome));
    #
    #So, intead, go directly to the TAXONOMY file and make the check ourselves:

    my $taxonomy = `cat "$genomedir/TAXONOMY"`;
    if ($taxonomy =~ /^\s*(Archaea|Bacteria|Eukaryota)/) {
	print STDERR "TAXONOMY=$taxonomy" 
	    if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} > 1));
    } else {
	print STDERR "Taxonomy $taxonomy is not an A,B,E genome --- skipping\n"
	    if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} > 1));
	next;
    }
    
    $ttlen = $inbig = 0;
    foreach $contig ( grep { $_ =~ /^contigs\d*$/ } readdir( GENOME ) )
    {
	open( LENGTHS, "fastasize < $genomedir/$contig |" ) || die "Could not pipe-open fastasize < $genomedir/$contig";
	
	while ( defined( $_ = <LENGTHS> ) )
	{
	    chomp;
	    if ( ( $len ) = $_ =~ /\s(\d+)$/ )
	    {
		$ttlen += $len;
		if ($len >= $minlen) { $inbig += $len }
	    }
	}
	close( LENGTHS ) || die "Could not close pipe fastasize < $genomedir/$contig";
    }
    
    if (not $ttlen) 
    {
	print STDERR "Could not get total contig lengths from $genomedir/$contig --- skipping"
	    if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} > 1));
	next;
    }

    print STDERR "$genomedir:\tttlen = $ttlen,\tinbig = $inbig,\tfrac = ", int(0.5 + 10000 * $inbig / $ttlen)/100, "%\n";

    if (($ttlen >= $minsize) && ($inbig >= $minfrac *$ttlen) )
    {
	print STDERR "Write-opening $genomedir/PROBABLY_COMPLETE\n" 
	    if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} > 1));
	open ( COMP, ">$genomedir/PROBABLY_COMPLETE" ) || die "Could not write-open $genomedir/PROBABLY_COMPLETE\n";
	my $frac = sprintf "%.2f", 100* $inbig / $ttlen;
	print COMP "Judged to be complete because $frac% of the nucleotides are in contigs >= $minlen bp\n";
	close COMP;
    }
    else
    {
	if (-e "$genomedir/PROBABLY_COMPLETE")
	{
	    print STDERR "Contigs $genomedir/contigs failed to pass cuts --- deleting existing file $genomedir/PROBABLY_COMPLETE\n"
		if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} > 1));
	    unlink "$genomedir/PROBABLY_COMPLETE"
		|| die "Failed attempt to remove $genomedir/PROBABLY_COMPLETE\n";
	}
    }
}

print STDERR "\n" if $ENV{VERBOSE};
print STDERR "\n" if (defined($ENV{VERBOSE}) && ($ENV{VERBOSE} == 1));

exit(0);
