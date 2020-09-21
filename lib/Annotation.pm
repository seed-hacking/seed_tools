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
# Simple object to hold annotations in-memory.
#

=pod

#TITLE OldAnnotationPm

This is a simple object that is used by [[FigPm]] to keep annotations in memory.

=cut

package Annotation;

use strict;

sub new
{
    my($class, $fid, $anno_time, $made_by, @anno) = @_;


    my $self = {
        fid => $fid,
        anno_time => $anno_time,
        made_by => $made_by,
        anno => [@anno],
    };
    
    #
    # See if it's an assignment annotation.
    #

    if ($anno[0] =~ /^Set\s+(\S+)\s+function\s+to$/)
    {
        $self->{is_assignment} = 1;
        $self->{assignment_who} = $1;
        $self->{assignment} = $anno[1];

        # parentheses required for ActivePERL
        $self->{is_master_assignment} = ($1 eq 'master' or $1 eq 'FIG');
    }

    return bless($self, $class);
}

sub fid
{
    my($self) = @_;
    return $self->{fid};
}

sub anno_time
{
    my($self) = @_;
    return $self->{anno_time};
}

sub made_by
{
    my($self) = @_;
    return $self->{made_by};
}

sub is_assignment
{
    my($self) = @_;
    return $self->{is_assignment};
}

sub is_master_assignment
{
    my($self) = @_;
    return $self->{is_master_assignment};
}

sub assignment_who
{
    my($self) = @_;
    return $self->{assignment_who};
}

sub assignment
{
    my($self) = @_;
    return $self->{assignment};
}


sub set_assignment_who
{
    my($self, $val) = @_;
    my $old = $self->{assignment_who};
    $self->{assignment_who} = $val;
    $self->{anno}->[0] = "Set $val function to";
    return $old;
}

sub as_text
{
    my($self) = @_;

    return join("\n", @$self{'fid', 'anno_time', 'made_by'}, @{$self->{anno}});
}

sub anno_text
{
    my($self) = @_;

    return join("\n", @{$self->{anno}});
}

1;
