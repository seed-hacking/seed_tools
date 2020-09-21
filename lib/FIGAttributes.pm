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


package FIGAttributes;

#require v5.6.0;

use Exporter;
use base 'Exporter';
@EXPORT = qw(%_FunctionAttributes);

use strict;

no warnings 'redefine';


#
# We use the methods below to process the subroutine attributes.
#
# See
#	http://www.perldoc.com/perl5.8.4/pod/perlsub.html
#	http://www.perldoc.com/perl5.8.4/lib/attributes.html
#	http://www.perldoc.com/perl5.8.4/lib/Attribute/Handlers.html
#
# for details.
#
# We take action on the following attributes:
#
#   Remote  Allow remote access to this function.
#   Scalar  Function returns a scalar.
#   List    Function returns a list.
#	method	This function is actually a method, and expects $self as the first arg.
#
# If an attribute is set for a function in some package, we
# will update a table %Package::_FunctionAttributes which is
# keyed on the (string) for the coderef, and which has value
# of a list of attribute values.
#

sub MODIFY_CODE_ATTRIBUTES
{
    my($pkg, $sub, @attrs) = @_;

    for my $attr (@attrs)
    {
	set_code_attribute($pkg, $sub, $attr);
    }
    return ();
}

sub set_code_attribute
{
    my($package, $coderef, $attr, $data) = @_;

#    print "Set code attr for $package $coderef $attr d=$data\n";
    my $str = "package $package;\n";
    $str .= <<'EOE';
$_FunctionAttributes{$coderef}->{$attr} = $data;
#print "in set_code_attribute: \n", Dumper(\%_FunctionAttributes);
EOE
    eval $str;
    if ($@ ne "")
    {
	warn "EVAL failed: $@\n";
    }
}

1;
