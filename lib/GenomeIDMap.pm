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

#
# Module for implementation of genome mapping (for FC computation among others).
#


package GenomeIDMap;

use strict;
use Carp;

sub new
{
    my($class, $fig) = @_;

    my $db = $fig->db_handle;
    my $self = {
	fig => $fig,
	db => $db,
	dbh => $db->{_dbh},
	id_cache => {},
	num_cache => {},
    };

    if (!$db->table_exists('genome_mapping'))
    {
	#
	# Ugh. We're using serials, but have to do it differently between Pg and mysql.
	#

	if ($db->{_dbms} eq "Pg")
	{
	    $db->create_table(tbl => "genome_mapping",
			      flds => qq(genome varchar(32) primary key,
					 gnum smallint  default nextval('genome_mapping_gnum_seq') not null unique
					));
	    $db->{_dbh}->do("create sequence genome_mapping_gnum_seq");
	}
	elsif ($db->{_dbms} eq "mysql")
	{
	    $db->create_table(tbl => "genome_mapping",
			      flds => qq(genome varchar(32) primary key,
					 gnum smallint unique auto_increment
					));
	}
	else
	{
	    confess "We don't know how to set up autonumbers on database other than Postgres and Mysql.";
	}
    }

    return bless($self, $class);
}

sub map_genome_id_to_gnum
{
    my($self, $genome_id) = @_;


    my $id_cache = $self->{id_cache};
    my $num_cache = $self->{num_cache};

    if (exists($id_cache->{$genome_id}))
    {
	return $id_cache->{$genome_id};
    }

    #
    # Look it up.
    #

    my $dbh = $self->{dbh};
    my $res = $dbh->selectcol_arrayref("select gnum from genome_mapping where genome = ?",
				       undef, $genome_id);
    if ($res and @$res == 1)
    {
	my $id = $res->[0];
	$id_cache->{$genome_id} = $id;
	$num_cache->{$id} = $genome_id;
	return $id;
    }

    #
    # It doesn't exist. Insert and retrieve id. In order to get the inserted id,
    # we have to use the DBI routines.
    #

    my $sth = $dbh->prepare("insert into genome_mapping(genome) values(?)");
    $sth->execute($genome_id);
    my $id = $self->{db}->get_inserted_id('genome_mapping', $sth, 'gnum');

    if ($id)
    {
	$id_cache->{$genome_id} = $id;
	$num_cache->{$id} = $genome_id;
	return $id;
    }

    confess "Could not insert into genome_mapping table";
	
}

sub map_gnum_to_gid
{
    my($self, $gnum) = @_;

    my $num_cache = $self->{num_cache};

    if (exists($num_cache->{$gnum}))
    {
	return $num_cache->{$gnum};
    }

    #
    # Look it up.
    #

    my $dbh = $self->{dbh};
    my $res = $dbh->selectcol_arrayref("select genome from genome_mapping where gnum = ?",
				       undef, $gnum);
    if ($res and @$res == 1)
    {
	my $genome_id = $res->[0];
	$self->{id_cache}->{$genome_id} = $gnum;
	$num_cache->{$gnum} = $genome_id;
	return $genome_id;
    }

    confess "No mapping found for gnum $gnum\n";
	
}

sub map_peg_to_nums
{
    my($self, $peg) = @_;

    if ($peg =~ /^fig\|(\d+\.\d+)\.peg\.(\d+)$/)
    {
	my $gnum = $self->map_genome_id_to_gnum($1);
	return ($gnum, $2);
    }
}

sub map_nums_to_peg
{
    my($self, $gnum, $pnum) = @_;

    my $gid = $self->map_gnum_to_gid($gnum);

    if (defined($gid))
    {
	return "fig|$gid.peg.$pnum";
    }
}
    


1;
