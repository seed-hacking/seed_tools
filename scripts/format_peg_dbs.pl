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
my $fig = new FIG;

my $path = ( $FIG_Config::blastbin ? $FIG_Config::blastbin : $FIG_Config::ext_bin );
my $genome;
foreach $genome ($fig->genomes("complete"))
{
    my $db = "$FIG_Config::organisms/$genome/Features/peg/fasta";
    next if (! -s $db);
    
    &verify_db($db,"p");
}

sub verify_db {
    my($db,$type) = @_;

    if ($type =~ /^p/i)
    {
		if ((! -s "$db.psq") || (-M "$db.psq" > -M $db))
		{
			system "$path/formatdb -p T -i $db";
		}
    }
    else
    {
		if ((! -s "$db.nsq") || (-M "$db.nsq" > -M $db))
		{
			system "$path/formatdb -p F -i $db";
		}
    }
}	
