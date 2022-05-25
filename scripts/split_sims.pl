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


use FIG;

$usage = "split_sims Dir Prefix  [ File1 File2 ... ]  (Default input is from STDIN)";

(($dir = shift @ARGV) &&
 ($pre = shift @ARGV)
)
    || die "\n\t$usage\n\n";

(-d $dir) || mkdir($dir,0777) || die "Could not create $dir: $!";
open(OUT,">$dir/$pre.1") || die "Could not write-open $dir/$pre.1";
$cnt = 0;
$n = 2;

if (! @ARGV) { push @ARGV, "-"; }

foreach $file (@ARGV)
{
    if ($file =~ /\.gz$/)
    {
	open(IN,"zcat $file |") || die "Could not open zcat-pipe from $file";
    }
    else
    {
	open(IN,"<$file") || die "Could not read-open $file";
    }

    $last = "";
    while (defined($_ = <IN>))
    {
	if ($_ =~ /^(\S+)/)
	{
	    $id = $1;
	    if (($id ne $last) && ($cnt > 3000000))
	    {
		close(OUT);
		open(OUT,">$dir/$pre.$n") || die "Could not write-open $dir/$pre.$n";
		$n++;
		$cnt = 0;
	    }
	    $last = $id;
	    print OUT $_;
	    $cnt++;
	}
    }
    close(IN);
}
close(OUT);
