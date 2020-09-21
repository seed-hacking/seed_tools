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
# Subsystem construct utilities.
#

package Construct;
use strict;

use Exporter;
use base qw(Exporter);

use vars qw(@EXPORT_OK);

@EXPORT_OK = qw(validate_constructs get_constructs_from_form
		make_html write_constructs_file parse_constructs_file);

#
# Walk a constructs list.
#
# For each require element check to see
# if it is a role reference (if it is a number) or if it
# is a construct name, in which case it should not overlap with
# an active role in the subsystem or any of its aliases.
#
sub validate_constructs
{
    my($cons, $sub, $errs) = @_;

    my %construct_names;

    #
    # Scan once and fill in the list of construct names.
    #

    for my $con_ent (@$cons)
    {
	my($name, $req_list) = @$con_ent;

	#
	# Ensure construct names aren't role names.
	#

	if (defined($sub->get_role_index($name)))
	{
	    push(@$errs, "Construct name <i>$name</i> is also a role name.");
	}
	elsif (defined(my $r = $sub->get_role_from_abbr($name)))
	{
	    push(@$errs, "Construct name <i>$name</i> is a role abbreviation for <i>$r</i>.");
	}
	else
	{
	    $construct_names{$name}++;
	}
    }

    #
    # Scan again and walk the requires lists and validate them.
    #

    for my $con_ent (@$cons)
    {
	my($name, $req_list) = @$con_ent;

	for my $req (@$req_list)
	{
	    my ($type, $req_name) = @$req;

	    #
	    # See if $req_name is numeric and matches a role.
	    #

	    if ($req_name =~ /^\d+$/)
	    {
		my $role = $sub->get_role($req_name - 1);

		if ($role)
		{
		    $req->[0] = 'R';
		    $req->[1] = $sub->get_role_abbr($req_name - 1);
		    $req->[2] = $role;
		}
		else
		{
		    $req->[0] = 'R';
		    $req->[1] = '*' . $req_name;
		    $req->[2] = undef;
		}   
	    }
	    else
	    {
		#
		# otherwise, it can be a construct, role, or role alias
		#

		my $role = $sub->get_role_from_abbr($req_name);

		if (defined($role))
		{
		    my $idx;
		    $idx = $sub->get_role_index($role);
		    if (defined($idx))
		    {
			$req->[0] = 'R';
			$req->[2] = $role;
		    }
		    else
		    {
			#
			# This shouldn't happen.
			#

			warn "construct.cgi: heap big error in subsystem data\n";
			$req->[0] = 'X';
			$req->[2] = undef;
			push(@$errs, "Spreadsheet data may be corrupted: Role $role found from abbr $req_name did not have an index");
		    }
		}
		else
		{
		    #
		    # Not a role or abbreviation. See if it's a construct name.
		    #
		    if (defined($construct_names{$req_name}))
		    {
			$req->[0] = 'C';
		    }
		    else
		    {
			$req->[0] = 'X';
			push(@$errs, "Requirement element <i>$req_name</i> in construct <i>$name</i> is not a role abbreviation or construct name");
		    }
		}   
	    }
	}
    }

    return (@$errs == 0);
}

sub get_constructs_from_form
{
    my($cgi) = @_;
    
    my $row_idx = 0;

    my(@ret);
    
    while (1)
    {
	my $name = $cgi->param("construct_name_$row_idx");
	my $req_str = $cgi->param("construct_req_$row_idx");

	last unless defined($name);

	$name =~ s/^\s+//;
	$name =~ s/\s+$//;

	if ($name ne '')
	{

	    my $reqs = [];
	    for my $req (split(/,/, $req_str))
	    {
		$req =~ s/^\s+//;
		$req =~ s/\s+$//;
		push(@$reqs, [undef, $req]);
	    }
	    push(@ret, [$name, $reqs]);
	}
	$row_idx++;
    }
    return @ret;
}

sub make_html
{
    my($cons, $cgi) = @_;

    my(@table, @ext_cons);

    my $row_idx = 0;

    @ext_cons = @$cons;
    my $n_blanks = 5;
    for  (1..$n_blanks)
    {
	push(@ext_cons, ['', []]);
    }
    
    for my $con_ent (@ext_cons)
    {
	my($name, $req_list) = @$con_ent;
	
	my $row = [];
	
	push(@$row, $cgi->textfield(-name => "construct_name_$row_idx",
				    -override => 1,
				    -value => $name));
	
	my(@name_strings);
	for my $req (@$req_list)
	{
	    my ($type, $req_name) = @$req;
	    
	    #
	    # For now, just use the name from the file.
	    # Eventually, we'll do lookups from role # -> name,
	    # etc. (though no row #s in saved files).
	    #
	    push(@name_strings, $req_name);
	}
	
	push(@$row, $cgi->textfield(-name => "construct_req_$row_idx",
				    -override => 1,
				    -value => join(", ", @name_strings),
				    -size => 80));
	push(@table, $row);
	$row_idx++;
    }

    return HTML::make_table(["Construct Name", "Requires"],
			    \@table,
			    "Constructs");
}

sub write_constructs_file
{
    my($cons, $file) = @_;

    my($fh);

    open($fh, ">$file") or die "Cannot write $file: $!\n";

    for my $con_ent (@$cons)
    {
	my($name, $req_list) = @$con_ent;

	print $fh "$name\n";
	for my $req (@$req_list)
	{
	    my ($type, $req_name, $role_name) = @$req;

	    if ($type eq "R")
	    {
		print $fh "$type $role_name\n";
	    }
	    elsif ($type eq "C")
	    {
		print $fh "$type $req_name\n";
	    }
	}

	print $fh "//\n";
    }
    close($fh);
}


sub parse_constructs_file
{
    my($file, $sub) = @_;

    my($fh);

     if (!open($fh, "<$file"))
    {
	die "Cannot open $file: $!\n";
    }
    
    local($/);

    $/ = "//\n";

    my(@ret);

    while (<$fh>)
    {
	chomp $_;
	next if $_ eq '';
	my($name, @reqs) = split(/\n/, $_);

	my $reqlist = [];

	foreach my $req (@reqs)
	{
	    my($type, $rname) = split(/\s+/, $req, 2);
	    my $abbr;

	    $rname =~ s/\s+$//;

	    if ($type eq "R")
	    {
		my $idx = $sub->get_role_index($rname);
		if (defined($idx))
		{
		    $abbr = $sub->get_role_abbr($idx);
		}
		else
		{
		    $abbr = $rname;
		    $abbr =~ s/,/_/g;
		}
		push(@$reqlist, [$type, $abbr, $rname]);
	    }
	    else
	    {
		push(@$reqlist, [$type, $rname]);
	    }
	}
	push(@ret, [$name, $reqlist]);
    }
    return @ret;
}

1;
