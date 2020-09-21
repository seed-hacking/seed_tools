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


###########################################
use strict;

use FIG;
my $fig = new FIG;
use Tracer;

my $usage = "usage: init_maps";

use DBrtns;

&initialize_maps;

undef $fig;

sub initialize_maps{
    my $dbf = $fig->{_dbf};

    $dbf->drop_table( tbl => "ec_map" );
    $dbf->create_table( tbl => "ec_map",
                       flds => "ec varchar(100), map varchar(100)"
                       );
    $dbf->drop_table( tbl => "map_name" );
    $dbf->create_table( tbl => "map_name",
			flds => "map varchar(100) UNIQUE NOT NULL, mapname varchar(200), primary key ( map )"
			);
    
    $dbf->create_index( idx  => "index_ec_map_ec",
                        tbl  => "ec_map",
                        type => "btree",
                        flds => "ec" );
    $dbf->create_index( idx  => "index_ec_map_map",
                        tbl  => "ec_map",
                        type => "btree",
                        flds => "map" );
Trace("Map tables initialized.") if T(2);
}

