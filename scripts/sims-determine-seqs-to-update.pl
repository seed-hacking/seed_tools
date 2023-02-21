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
use gjoseqlib;
use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o old-nr new-nr target-dir",
	["max-chunk=i" => "Maximum chunk size (number of sequences)", { default => 500_000 }],
	["help|h" => "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 3;

my $old_nr = shift;
my $new_nr = shift;
my $target_dir = shift;

if (! -d $target_dir)
{
    mkdir($target_dir) or die "Cannot mkdir $target_dir: $!";
}

my $old_ids = &load_ids($old_nr);

open(NEW, "<", $new_nr) or die "Cannot open new nr $new_nr: $!";

my $fidx = 1;

my $fn = sprintf("$target_dir/seqs.%03d", $fidx++);
open(SEQS, ">", $fn) or die "Cannot write $fn: $!";
my $cur_count = 0;
while (my($id, $def, $seq) = read_next_fasta_seq(\*NEW))
{
    if (! $old_ids->{$id})
    {
	print_alignment_as_fasta(\*SEQS, [$id, $def, $seq]);
	$cur_count++;
	if ($cur_count >= $opt->max_chunk)
	{
	    close(SEQS);
	    my $fn = sprintf("$target_dir/seqs.%03d", $fidx++);
	    open(SEQS, ">", $fn) or die "Cannot write $fn: $!";
	    $cur_count = 0;
	}
    }
}
close(SEQS);
close(NEW);

sub load_ids {
    my($nr) = @_;

    open(NR, "<", $nr) || die "Could not open $nr: $!";
    my $entries = {};
    while (<NR>)
    {
	if (/^>(\S+)/)
	{
	    $entries->{$1} = 1;
	}
    }
    close(NR);
    return $entries;
}
