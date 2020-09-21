# -*- perl -*-
########################################################################
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
########################################################################

use FIG;
use Tracer;

my $fig = new FIG;

# usage: compute_genome_counts [G1 G2 ...]

my($genome);

#  Build list of the genomes to be processed ---------------------------------------

my ($mode, @genomes) = FIG::parse_genome_args(@ARGV);

#  Gather the data on each genome --------------------------------------------

Trace("Gathering genome data.") if T(2);

my $genomedata = "$FIG_Config::temp/tmp$$";
open( REL, ">$genomedata" ) || Confess("Could not write-open $genomedata");

my $genome_md5sum = "$FIG_Config::temp/tmp_md5_$$";
open( MD5, ">$genome_md5sum" ) || Confess("could not write-open $genome_md5sum");

foreach $genome ( @genomes ) {
    my $genome_dir = "$FIG_Config::organisms/$genome";
    
    if ((! (-d $genome_dir)) || (-s "$genome_dir/DELETED")) {
	print STDERR "WARNING: $genome has been deleted\n";
	next;
    }
    
    Trace("Checking $genome_dir.") if T(3);
    
    if ( open( TMP, "<$genome_dir/GENOME" ) && defined( $name = <TMP> ) ) {
	close( TMP );
	chomp $name;
	$name =~ s/^\s*//o;
	$name =~ s/^Candidatus\s*//o;
	$name =~ s/\s*$//o;
    }
    else {
	$name = "Unknown sp.";
    }
    
    $complete     = ( -e "$genome_dir/COMPLETE" )     ? "1" : "0";
    $restrictions = ( -e "$genome_dir/RESTRICTIONS" ) ? "1" : "0";
    
    if (open(TMP,"<$genome_dir/TAXONOMY") && defined($_ = <TMP>) && ($_ =~ /(\S.*\S)/o)) {
	$taxonomy = $1;
	$taxonomy =~ s/^root\;\s*//io;               # Remove root node if present
	$taxonomy =~ s/^\s+//o;                      # Remove leading whitespace if present
	$taxonomy =~ s/\Candidatus\s+//igo ;         # Remove all occurences of "Candidatus" weasel-word
	$taxonomy =~ s/\s+/ /go;                     # Collapse multiple whitespaces
	$taxonomy =~ s/\s+$//o;                      # Remove trailing whitespace if present
	
	if    ($taxonomy =~ /^Eu[ck]ary/io)          { $tax = "Eukaryota"; }
	elsif ($taxonomy =~ /^Bacteri/io)            { $tax = "Bacteria";  }
	elsif ($taxonomy =~ /^Archae/io)             { $tax = "Archaea";   }
	elsif ($taxonomy =~ /^Vir/io)                { $tax = "Virus";     }
	elsif ($taxonomy =~ /^Environ/io)            { $tax = "Environmental Sample"; }
	elsif ($taxonomy =~ /^Un\S* Environ/io)      { $tax = "Environmental Sample"; }
	else                                         { $tax = "Unknown";   }
	
	if ($taxonomy =~ /plasmid/io)                { $tax = "Plasmid";   }  # Override "maindomain" if plasmid
	
	close( TMP );
    }
    else {
	# These were left at previous genome if open failed -- GJO
        $taxonomy = 'Unknown.';
        $tax      = 'Unknown';
    }
    
    
    #  We only allow bacteria, euk, and archaeal "complete" genomes.
    #  (Perhaps we should have left the COMPLETE file.  This might be used
    #  in the future, filtering elsewhere. -- GJO)
    
    if ( $complete && ( ( $tax =~ /^Environ/io ) || ( $tax !~ /^[EBA]/io ) ) ) {
	$complete = 0;
	unlink( "$genome_dir/COMPLETE" );
	print STDERR "NOTE: Clearing \"complete\" status of genome $genome, $name\n";
    }
    
    
    # Import genome statitistics and checksum from COUNTS file
    # (Check to make sure COUNTS is up to date, first)
    
    $sz = $cksum = 0;
    my $countfile = "$genome_dir/COUNTS";
    my $contigs = (<$genome_dir/contigs*>)[0];
    if ((!-e $countfile) || ((-M $countfile) > (-M $contigs))) {
	die qq(File $countfile does not exist or is out of date. Please run 'index_contigs')
    }
    
    open( COUNTS, "<$countfile" )
	|| die qq(Could not read-open $countfile);
    my $counts = <COUNTS>;
    close( COUNTS );
    if ( $counts ) {
	chomp $counts;
	( undef, undef, $sz, $cksum ) = split /\t/o, $counts;
	$sz    ||= 0;
	$cksum ||= 0;
    }
    
    $pegs = &count_features( $genome, "peg" );
    $rnas = &count_features( $genome, "rna" );
    
    print REL join("\t", $genome, $name, $sz, $cksum, $tax, $pegs, $rnas, $complete, $restrictions, $taxonomy), "\n";
    
    if (open(TMP, "<$genome_dir/MD5SUM")) {
	my $md5 = <TMP>;
	chomp($md5);
	print MD5 "$genome\t$md5\n";
    }
    else {
	warn "Could not read-open $genome_dir/MD5SUM";
    }
}

close(MD5);
close(REL);


#  Load the database ---------------------------------------
if (-s $genomedata) {
    $fig->reload_table($mode, "genome",
		       "genome varchar(16) UNIQUE NOT NULL, gname varchar(255), szdna BIGINT, "
		       . "cksum BIGINT, maindomain varchar(20), pegs INTEGER, "
		       . "rnas INTEGER, complete CHAR, restrictions CHAR, taxonomy text, "
		       . "PRIMARY KEY ( genome )",
		       { genome_ix => "genome" },
		       $genomedata, \@genomes);
    
    unlink($genomedata);
}
else {
    print STDERR "WARNING: No genome data to update\n";
}

if (-s $genome_md5sum) {
    $fig->reload_table($mode, "genome_md5sum",
		       "genome varchar(16) UNIQUE NOT NULL, md5sum char(32), PRIMARY KEY ( genome )",
		       { }, $genome_md5sum, \@genomes);
    
    unlink($genome_md5sum);
}
else {
    print STDERR "WARNING: No genome MD5 sums to update\n";
}
Trace("Genome counts complete.") if T(2);
exit(0);



# Subroutines --------------------------------------------------------

sub count_features {
    my($genome,$type) = @_;
    my $genome_dir = "$FIG_Config::organisms/$genome";    
    
    my %deleted;
    my $deleted = qq($genome_dir/Features/$type/deleted.features);
    if (-s $deleted) {
	%deleted = map { ($_ => 1) } $fig->file_read($deleted, qq(*));
    }
    
    my %seen;
    if (open(TBL,"<$genome_dir/Features/$type/tbl")) {
	while (defined($_ = <TBL>)) {
	    if (($_ =~ /^(fig\|\S+)/o) && (not $deleted{$1})) {
		$seen{$1} = 1;
	    }
	}
	close(TBL);
    }
    
    return (scalar keys %seen);
}

1;
