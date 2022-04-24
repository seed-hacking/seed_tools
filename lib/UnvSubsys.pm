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

package UnvSubsys;

use Subsystem;
use Carp;
use FIG;
use SFXlate;
use Data::Dumper;
use strict;
use Tracer;
use HTML;

=head1 Universal Subsystem Object

=head2 Introduction

The universal subsystem object provides methods used to display useful information
about a subsystem. Its intent is to support both Sprout and SEED.

The soul of a subsystem is its spreadsheet. The columns of the spreadsheet are
the subsystem roles. The rows are the genomes that contain the subsystem. PEGs
are stored in the spreadsheet's cells. Each genome can be identified by genome
ID or row index. Each role can be identified by role name or column index.

It is worth noting that the term I<subset> has two meanings-- a subset of roles
or a subset of genomes. A subset of roles is called a I<column subset>. A
subset of genomes is called a I<row subset>.

The object created contains a great deal of information, so it's worthwhile to
stop for a moment and discuss it. The object will have the following members.

=over 4

=item Roles

A list of 3-tuples, each consisting of a role name, the abbreviated role name,
and a list of URLs for the reactions. Each role is a column in the
subsystem spreadsheet. Indexing into the list by column index yields the
ID, abbreviation, and reactions for the role in that column.

=item RoleIndex

A hash mapping each role to its column index in the spreadsheet.

=item RoleSubsets

A list of 2-tuples, each consisting of a subset name followed by a list of the
column indices for the roles in the subset.

=item Genomes

A list of 2-tuples, each containing the genome ID for the relevant row and the
variant code for that genome. Subsystems can have several variations, and the
variant code indicates which variant of the subsystem that the genome uses.
There is one genome for each row in the spreadsheet. Indexing into the list
yields the ID of the genome for that row and its variant code.

=item GenomeIndex

A hash mapping each genome ID to its row index in the spreadsheet.

=item PegHash

A hash of hashes containing the spreadsheet cells. If C<$pegHash> is the hash
of hashes, C<$row> is a genome index, and C<$col> is a role index, then

    $pegHash->{$row}->{$col}

returns a reference to a list of the IDs for the PEGs in the relevant
spreadsheet cell.

=item ColorHash

A hash mapping each PEG ID to the color that should be used to represent that
PEG in the display.

=item AliasHash

A hash mapping each PEG ID to a list of its aliases.

=item ReactionHash

A hash mapping each role to the list of reactions it catalyzes.

=back

=head2 Public Methods

=head3 new

    my $usub = UnvSubsys->new($ssa, $fig, $active_subsetR, $focus, $show_clusters, $aliases, \@peg_colors);

Construct a new universal subsystem object for a specified subsystem.

=over 4

=item ssa

Name of the subsystem.

=item fig

Data access object, either a FIG object or an SFXlate object.

=item active_subsetR

Name of the active row subset. A row subset names a group of genomes.

=item focus

ID of the genome currently in focus.

=item show_clusters

TRUE if clusters should be painted by color, else FALSE.

=item aliases

TRUE if PEG aliases should be shown, else FALSE.

=item peg_colors (optional)

Reference to a list of 2-tuples, each 2-tuple consisting of a PEG ID and the color
to be assigned to the PEG.

=back

=cut

sub new {
    # Get the parameters.
    my($class, $ssa, $fig, $active_subsetR, $focus, $show_clusters, $aliases, $peg_colors) = @_;
    # Fix the subsystem name. Spaces are replaced by underscores to avoid naming problems
    # in the seed's file system.
    $ssa =~ s/ /_/g;
    # Get the FIG object. At this point in time, we're only getting data from the SEED. On
    # a future pass through the code, we'll get data from a SEED or Sprout.
    if ((ref($fig) eq "FIG") || (ref($fig) eq "FIGV") || (ref($fig) eq 'SFXlate')) {
        # Create a subsystem object. The "get_subsystem" method is provided by both FIG
        # and SFXlate; however, the object returned in each case is different: FIG
        # returns a "Subsystem" object; SFXlate returns a "SproutSubsys" object.
        my $subsystem = $fig->get_subsystem($ssa);
        # Get the key subsystem data. Note that CRs in the notes are converted to LFs
        # in case the notes come from a Mac.
        my $curator = $subsystem->get_curator;
        my $notes = $subsystem->get_notes;
        $notes =~ s/\r/\n/g;
        my @roles = $subsystem->get_roles;
        my $reactions = $subsystem->get_reactions;
        my @genomes = $subsystem->get_genomes;
        my @col_subsets = $subsystem->get_subset_namesC;
        my @diagrams = $subsystem->get_diagrams();
        # Create the data structures for the role list and the role index.
        my $role_info = [];
        my $roleH     = {};
        # Loop through the roles to create the role list. The $i index will move
        # through the columns of the spreadsheet.
        my($i,$j,$subset,$peg);
        for ($i=0; ($i < @roles); $i++)
        {
            # Get the ID of the role in the current column.
            my $role = $roles[$i];
            # Extract its abbreviation.
            my $abbrev = $subsystem->get_role_abbr($i);
            # Get the reactions. Note that its possible no reactions were found. If there were
            # reactions found, however, we convert from reaction IDs to a comma-delimited
            # list of HTML code full of hyperlinks.
            my $react = $reactions ? join(",", map { &HTML::reaction_link($_) } @{$reactions->{$role}}) : [];
            # Form the role name, abbreviation, and reaction list into a 3-tuple and push
            # them onto the main role list.
            push(@$role_info,[$role,$abbrev,$react]);
            # Set the role hash so that we can get back the column index for a given
            # role name.
            $roleH->{$role} = $i;
        }
        # Get the column subsets. A column subset is a list of role IDs, so we need to
        # convert it to a set of column indices. We ignore the special "All" subset
        # that contains everything.
        my $subset_info = [];
        foreach $subset (@col_subsets)
        {
            if ($subset ne 'All')
            {
                push(@$subset_info,[$subset,[map { $roleH->{$_} } $subsystem->get_subsetC_roles($subset)]]);
            }
        }
        # Now we create the genome directories. For each genome we need to be able to
        # hash from the genomeID to the specified role index and we must be able to
        # get the genome's variant code. (The variant code indicates which subsystem
        # variant the genome uses.
        my $genomes_info = [];
        my $genomeH      = {};
        for ($i=0; ($i < @genomes); $i++)
        {
            # Get the genome ID for this row.
            my $genome  = $genomes[$i];
            # Get its variant code.
            my $variant = $subsystem->get_variant_code($i);
            # Form them into a 2-tuple and add the result to the genome list.
            push(@$genomes_info,[$genome,$variant]);
            # Set up the hash to get from the genome ID to the row index.
            $genomeH->{$genome} = $i;
        }

        # Next we gather the data from the actual spreadsheet cells. For an SFXlate
        # object, this is the most expensive part of the operation, since it requires
        # a database call for each cell.
        my $pegH = {};
        # $i is the row index, which cycles through genomes.
        for ($i=0; ($i < @genomes); $i++)
        {
            # $j is the column index, which cycles through roles.
            for ($j=0; ($j < @roles); $j++)
            {
                my @pegs = $subsystem->get_pegs_from_cell($i,$j);
                $pegH->{$i}->{$j} = [@pegs];
            }
        }

        # Get the row subsets. Row subsets are determined by the genome
        # taxonomy, so these are different from the row subsets stored
        # in the subsystem object.
        my $row_subsets    = &row_subsets($fig, $genomeH);
        # Here we try to get a list of the active genomes. The caller
        # gives us hints in the form of the "$focus" and "$active_subsetR"
        # parameters. If, for example, we are using this object to generate
        # a subsystem web page, the focus information would steer us to
        # whatever the user wants to look at on the page.
        my $active_genomes = &active_genomes($fig, $row_subsets,$active_subsetR,$focus,
                                             $genomeH,$genomes_info);

        # Now we generate a table of colors for the various PEGs. If the
        # caller gave us a map of peg IDs to colors, we use that. Otherwise,
        # we allow the option of painting the PEGs by cluster number. (The
        # caller indicates this by setting the "show_clusters" flag in the
        # parameter list.)
        my $colorsH;
        if ($peg_colors)
        {
            # Here the caller gave us a list of peg colors. The list contains
            # 2-tuples, each consisting of a PEG ID followed by the
            # color value. The loop below extracts the pairs and stuffs them
            # into the color hash.
            $colorsH = {};
            foreach $_ (@$peg_colors)
            {
                my($peg,$color) = @$_;
                $colorsH->{$peg} = $color;
            }
        }
        elsif ($show_clusters)
        {
            # Here the user wants us to base the colors on the genome clustering
            # information.
            $colorsH  = &set_colors($fig,$subsystem,$pegH,$active_genomes);
        }
        else
        {
            # Here the user is not interested in coloring the PEGs.
            $colorsH = {};
        }
        # If the user wants to see aliases, compute the alias hash. Aliases
        # will only be computed for PEGs belonging to active (highlighted)
        # genomes, and there is a maximum of one alias per PEG.
        my $aliasesH = $aliases ? &set_aliases($fig,$pegH,$active_genomes) : {};
        # Create and bless the UnvSubsys object.
        my $self = { SSA => $ssa,
		     Roles => $role_info,
                     RoleIndex => $roleH,
                     RoleSubsets => $subset_info,
                     Genomes => $genomes_info,
                     GenomeIndex => $genomeH,
                     GenomeSubsets => $row_subsets,
                     PegHash => $pegH,
                     Colors => $colorsH,
                     Aliases => $aliasesH,
                     Curator => $curator,
                     Notes => $notes,
                     Reactions => $reactions,
                     Diagrams => \@diagrams,
                   };
        bless($self, $class);
        # Return the object.
        return $self;
    }
    else
    {
        # Here the FIG-like object was not recognized, so we return an
        # undefined value.
        return undef;
    }
}

=head3 get_ssa

    my $ssa = $unvsub->get_ssa();

Return the name of the subsystem

=cut

sub get_ssa {
    my($self) = @_;
    return $self->{SSA};
}

=head3 get_ssa_pretty

    my $ssa = $unvsub->get_ssa_pretty();

Return the 'prettyfied' name of the subsystem

=cut

sub get_ssa_pretty{
    my($self) = @_;
    my $ssa = $self->{SSA};
    $ssa =~ s/_/ /g;
    return $ssa;
}


=head3 get_subset_namesR

    my @names = $unvsub->get_subset_namesR();

Return the names of the genome (row) subsets.

=cut

sub get_subset_namesR {
    my($self) = @_;

    return map { $_->[0] } @{$self->{GenomeSubsets}};
}

=head3 get_subsetR

    my @genomes = $unvsub->get_subsetR($set);

Return a list of the genome IDs covered by a row subset.

=over 4

=item set

Name of the row subset whose genomes are desired.

=item RETURN

Returns a list of the IDs for the genomes found in the specified row
set.

=back

=cut

sub get_subsetR {
    # Get the parameters.
    my($self,$set) = @_;
    my($i);
    # Get the list of row subsets.
    my $sets = $self->{GenomeSubsets};
    # Find the row subset with the specified name. The row subset list is a
    # list of 2-tuples, and the first element of each tuple is the set
    # name.
    for ($i=0; ($i < @$sets) && ($sets->[$i]->[0] ne $set); $i++) {}
    if ($i < @$sets)
    {
        # Here we found the named subset. The subset tuple's second element is
        # the list of row indices. We map these to genome IDs before returning
        # them.
        return map { $self->{Genomes}->[$_]->[0] } @{$sets->[$i]->[1]}
    }
    # Here we subset was not found, so we return the undefined value.
    return undef;
}

=head3 get_subsetR

    my @pairs = $unvsub->get_subsetsR();

Return a list of all the row subsets. The subsets are returned in the form
of 2-tuples, each consisting of a subset name followed by a reference to a
list of genome IDs. The genome IDs correspond to the rows in the subset.

=cut

sub get_subsetsR {
    # Get the parameters.
    my($self) = @_;
    # Extract the list of genome subsets. This list is in the form of
    # 2-tuples, but the rows are identified by row index, not genome ID.
    my $sets = $self->{GenomeSubsets};
    # Create the return list.
    my @pairs = ();
    # Loop through the subsets.
    my $pair;
    foreach $pair (@$sets)
    {
        # Convert this subset's member list from row indices to genome IDs
        # and stash the result in the return list.
        my($id,$members) = @$pair;
        push(@pairs,[$id,[map { $self->{Genomes}->[$_]->[0] } @$members]]);
    }
    # Return the list constructed.
    return @pairs;
}

=head3 row_subsets

    my $subsetList = UnvSubsys::row_subsets($fig, \%genomeH);

This method computes the taxonomic row subsets for a subsystem. It takes
as input a hash that maps genome IDs to column indices and a FIG object.
The FIG object provides a list of taxonomic groups of 10 or more complete
genomes. From the list, we extract subsets which have more than 10
genomes in the list of subsystem genomes. If no such subsets exist,
we extract subsets which have at least 1 genome from the list of
subsystem genomes. The subsets are returned as 2-tuples, the first
element being the subset ID and the second being a list of row indices
for the genomes in the subset.

=over 4

=item fig

A FIG-like object for accessing the data store.

=item genomeH

Reference to a hash that maps each genome ID to its row index in the
subsystem spreadsheet.

=item RETURN

Returns a reference to a list of 2-tuples. Each 2-tuple consists of a
subset ID followed by a list of the row indices for the genomes in the
subset.

=back

=cut

sub row_subsets {
    my ($fig, $genomeH) = @_;

    # We need a real FIG object, since SFXlate does not yet support
    # taxonomy trees. The "FIG" method does this for us.
    $fig = $fig->FIG();
    # Initialize the return value.
    my $subsets = [];
    # Get a list of taxonomic groups. This will come back as a list of
    # 2-tuples.
    my $taxonomic_groups = $fig->taxonomic_groups_of_complete(5);
    # Loop through the 2-tuples. We're looking for subsets which
    # contain at least one genome on the subsystem's spreadsheet.
    my($pair,$id,$members);
    foreach $pair (@$taxonomic_groups)
    {
        ($id,$members) = @$pair;
#	warn "Group $id is @$members\n";
#	if ($id eq 'All')
#	{
#	    push(@$members, '372461.6');
#	}

        # Extract the genomes in the member list that participate in this
        # subsystem. To do this, we convert each genome ID to its row
        # index. If no row index exists, the GREP condition discards the
        # member.
        my @mem = grep { defined($_) } map { $genomeH->{$_} } @$members;
        # If there are enough members, save the subset.
        if (@mem > 0)
        {
            push(@$subsets,[$id,[@mem]]);
        }
    }
    # Return the list of row subsets.
    return $subsets;
}

=head3 set_aliases

    my $aliasHash = UnvSubsys::set_aliases($fig, $pegH, $active_genomes);

Return a hash mapping PEG IDs to aliases.

=over 4

=item fig

FIG-like object that can be used to access the data store.

=item pegH

Reference to the spreadsheet hash table. Given a row index I<$row> and a
column index I<$col>,

    $pegH->{$row}->{$col}

will return a reference to a list of PEGs in the specified spreadsheet cell.

=item active_genomes

Reference to a hash whose keys correspond to the spreadsheet row indices
of genomes that should be highlighted.

=item RETURN

Returns a hash that takes as input a PEG ID and returns an alias. Only PEGs
for active genomes will be processed.

=back

=cut

sub set_aliases {
    # Get the parameters.
    my($fig,$pegH,$active_genomes) = @_;
    my($genomeI,$roleI,$pegs,$peg,$roleH);

    # Create the return hash.
    my $aliasesH = {};

    # Loop through each row that corresponds to an active genome.
    # The active genome list contains row indices, and there is
    # one genome per row. Note that the genome ID is never used,
    # only the row index.
    foreach $genomeI (grep { $active_genomes->{$_} } keys(%$pegH))
    {
        # Get the role hash for the specified genome. The role hash
        # maps column indices (representing roles) to lists of PEGs.
        $roleH = $pegH->{$genomeI};
        # Loop through the role (column) indices.
        foreach $roleI (keys(%$roleH))
        {
            # Get the PEG list for this row/column combination.
            $pegs = $roleH->{$roleI};
            # Only proceed if data was found in the cell.
            if (defined $pegs) {
                # Loop through the pegs in the cell.
                foreach $peg (@$pegs)
                {
                    # If we do not already have an alias for this PEG,
                    # compute one.
                    if (! $aliasesH->{$peg})
                    {
                        $aliasesH->{$peg} = scalar &ext_id($fig,$peg);
                    }
                }
            }
        }
    }
    # Return the hash we built.
    return $aliasesH;
}

=head3 set_colors

    my $colorHash = UnvSubsys::set_colors($fig, $sub, \%pegH, \%active_genomes);

Return a hash that maps each PEG in the subsystem spreadsheet to a display
color. Not all PEGs need to be mapped. Those that do not have a color
assigned will generally be displayed in white.

=over 4

=item fig

FIG-like object that can be used to access the data store.

=item sub

Subsystem object for the current subsystem.

=item pegH

Reference to the spreadsheet hash table. Given a row index I<$row> and a
column index I<$col>,

    $pegH->{$row}->{$col}

will return a reference to a list of PEGs in the specified spreadsheet cell.

=item active_genomes

Reference to a hash whose keys correspond to the spreadsheet row indices
of genomes that should be highlighted.

=item RETURN

Returns a hash that takes as input a PEG ID and returns a color. These colors
are used when displaying the PEGs in the subsystem spreadsheet. Only PEGs
for active genomes will be colored.

=back

=cut

sub set_colors {
    # Get the parameters.
    my($fig,$sub,$pegH,$active_genomes) = @_;

    my($genomeI,$roleI,$pegs,$peg,$roleH,%pegs_in_genome);
    # Create the return hash.
    my $colorsH = {};
    # Loop through the active genomes. The keys of "%$pegH" are the row indices
    # for rows that have at least one occupied cell in the spreadsheet. The
    # Grep then reduces this list to those rows that are highlighted.
    foreach $genomeI (grep { $active_genomes->{$_} } keys(%$pegH))
    {
        # We will use the following hash to compile a list of all the PEGs
        # in the spreadsheet for the current genome, that is, the genome
        # represented by row index "$genomeI".
        undef %pegs_in_genome;
        # Get the hash for the current row. This hash maps column indices to
        # lists of pegs.
        $roleH = $pegH->{$genomeI};
        # Loop through the column indices for the specified row.
        foreach $roleI (keys(%$roleH))
        {
            # Here we can finally get a list of the pegs in this spreadsheet
            # cell. We loop through them, marking them in the "%pegs_in_genome"
            # hash.
            $pegs = $roleH->{$roleI};
            foreach $peg (@$pegs)
            {
                $pegs_in_genome{$peg} = 1;
            }
        }
        # Extract the "%pegs_in_genome" hash keys. This gives us a duplicate-free
        # list of all the pegs for the current spreadsheet role.
        my @pegs = keys(%pegs_in_genome);
        my($tuple,$peg,$color);
        # Get a hash that maps the PEG IDs to colors.
        my $colors_for_one_genome = &set_colors_for_genome($fig,$sub, \@pegs);
        # Loop through the hash we got back and assign the colors from that
        # hash to the master hash we're returning to the caller.
        while (($peg,$color) = each %$colors_for_one_genome)
        {
            $colorsH->{$peg} = $colors_for_one_genome->{$peg};
        }
    }
    # Return the color hash.
    return $colorsH;
}

=head3 set_colors_for_genome

    my $colorHash = UnvSubsys::set_colors_for_genome($fig, $sub, \@pegs);

Return a reference to a hash mapping the specified pegs to colors. PEGs that
are physically close to each other will be painted the same color.

=over 4

=item fig

A fig-like object that can be used to access the data store.

=item sub

Subsystem object for the relevant subsystem.

=item pegs

Reference to a list of PEG IDs. All of the peg IDs should be for the
same genome.

=item RETURN

Returns a reference to a hash that maps each PEG ID to a color.

=back

=cut

sub set_colors_for_genome {
    # Get the parameters.
    my($fig, $sub, $pegs) = @_;
    # Default all the PEGs to white.
    my %color_of = map { $_ => '#FFFFFF' } @$pegs;
    # Divide the pegs into clusters.
    my @clusters = $fig->compute_clusters($pegs, $sub);
    # Get a list of useful colors.
    my @colors =  &cool_colors();
    # If we have too many clusters, chop off the big ones at the end. These
    # are least likely to be important.
    if (@clusters > @colors) { splice(@clusters, 0, (@clusters - @colors)) }
    # Loop through the clusters.
    for my $cluster (@clusters) {
        # Get the color for this cluster.
        my $color = shift @colors;
        # Loop through this cluster, putting this color into the color_of
        # entries for each PEG.
        for my $peg (@$cluster) {
            $color_of{$peg} = $color;
        }
    }
    # Return the color map.
    return \%color_of;
}

=head3 cool_colors

    my @colorList = UnvSubsys::cool_colors();

Return a list of web-safe colors.

=cut

sub cool_colors {
 # just an array of "websafe" colors or whatever colors we want to use. Feel free to remove bad colors (hence the lines not being equal length!)
 return (
 '#C0C0C0', '#FF40C0', '#FF8040', '#FF0080', '#FFC040', '#40C0FF', '#40FFC0', '#C08080', '#C0FF00', '#00FF80', '#00C040',
 "#6B8E23", "#483D8B", "#2E8B57", "#008000", "#006400", "#800000", "#00FF00", "#7FFFD4",
 "#87CEEB", "#A9A9A9", "#90EE90", "#D2B48C", "#8DBC8F", "#D2691E", "#87CEFA", "#E9967A", "#FFE4C4", "#FFB6C1",
 "#E0FFFF", "#FFA07A", "#DB7093", "#9370DB", "#008B8B", "#FFDEAD", "#DA70D6", "#DCDCDC", "#FF00FF", "#6A5ACD",
 "#00FA9A", "#228B22", "#1E90FF", "#FA8072", "#CD853F", "#DC143C", "#FF6347", "#98FB98", "#4682B4",
 "#D3D3D3", "#7B68EE", "#2F4F4F", "#FF7F50", "#FF69B4", "#BC8F8F", "#A0522D", "#DEB887", "#00DED1",
 "#6495ED", "#800080", "#FFD700", "#F5DEB3", "#66CDAA", "#FF4500", "#4B0082", "#CD5C5C",
 "#EE82EE", "#7CFC00", "#FFFF00", "#191970", "#FFFFE0", "#DDA0DD", "#00BFFF", "#DAA520", "#008080",
 "#00FF7F", "#9400D3", "#BA55D3", "#D8BFD8", "#8B4513", "#3CB371", "#00008B", "#5F9EA0",
 "#4169E1", "#20B2AA", "#8A2BE2", "#ADFF2F", "#556B2F",
 "#F0FFFF", "#B0E0E6", "#FF1493", "#B8860B", "#FF0000", "#F08080", "#7FFF00", "#8B0000",
 "#40E0D0", "#0000CD", "#48D1CC", "#8B008B", "#696969", "#AFEEEE", "#FF8C00", "#EEE8AA", "#A52A2A",
 "#FFE4B5", "#B0C4DE", "#FAF0E6", "#9ACD32", "#B22222", "#FAFAD2", "#808080", "#0000FF",
 "#000080", "#32CD32", "#FFFACD", "#9932CC", "#FFA500", "#F0E68C", "#E6E6FA", "#F4A460", "#C71585",
 "#BDB76B", "#00FFFF", "#FFDAB9", "#ADD8E6", "#778899",
 );
}

=head3 ext_id

    my $externalID = UnvSubsys::ext_id($fig, $peg);

or

    my @externalIDs = UnvSubsys::ext_id($fig, $peg);

Return a list of non-FIG IDs for the specified feature. In a scalar context, return
a single non-FIG ID for the specified feature.

This method returns IDs that are all of the same type, that is, all UniProt IDs, or
all KEGG IDs, and so forth. To do this, it checks the feature's alias list for IDs
of a given type. If it finds at least one, then all IDs of that type are returned.
Highest priority is given to the UniProt IDs, then SP IDs, GI IDs, and finally
KEGG IDs.

=over 4

=item fig

A FIG-like object for accessing the data store.

=item peg

ID of the feature whose aliases are desired.

=item RETURN

In list context, a list of non-FIG IDs for the feature that are all of the same
type. In scalar context, the first non-FIG ID for the feature of the
highest-priority available type.

=back

=cut

sub ext_id {
    my($fig,$peg) = @_;

    my @tmp;
    my @aliases = $fig->feature_aliases($peg);
    if      ((@tmp = grep { $_ =~ /^uni\|/ } @aliases) > 0)
    {
        @aliases =  @tmp;
    }
    elsif   ((@tmp = grep { $_ =~ /^sp\|/ } @aliases) > 0)
    {
        @aliases = @tmp;
    }
    elsif   ((@tmp = grep { $_ =~ /^gi\|/ } @aliases) > 0)
    {
        @aliases = @tmp;
    }
    elsif   ((@tmp = grep { $_ =~ /^kegg\|/ } @aliases) > 0)
    {
        @aliases = @tmp;
    }
    else
    {
        @aliases = ();
    }

    if (wantarray())
    {
        return @aliases;
    }
    else
    {
        return $aliases[0];
    }
}

=head3 subsystem_curator

    my $name = $unvsub->subsystem_curator();

Return the name of the subsystem curator. The database stores user names as
C<master:>I<name>. This method strips off the C<master:> prefix before it
passes the result back to the caller.

=cut

sub subsystem_curator {
    my($self) = @_;

    my $curator = $self->{Curator};
    $curator =~ s/master://;
    return $curator;
}

=head3 get_roles

    my @roles = $unvsub->get_roles();

Return a list of the roles (columns) for this subsystem. The roles will be
returned in column order, so that if you access the sixth element of the
return list, you'll get the name of the role for the sixth column.

=cut

sub get_roles {
    my($self) = @_;
    # The role index in this object is a list of 3-tuples. The caller only
    # wants the first element of each tuple, which is the role name.
    return map { $_->[0] } @{$self->{Roles}};
}

=head3 get_genome_index

    my $index = $unvsub->get_genome_index($genome);

Return the row index of the specified genome.

=over 4

=item genome

ID of the genome whose row index is desired.

=item RETURN

Returns the index of the row corresponding to the specified genome, or an
undefined value if the genome is not represented in the subsystem
spreadsheet.

=back

=cut

sub get_genome_index {
    my($self,$genome) = @_;

    return $self->{GenomeIndex}->{$genome};
}

=head3 get_genomes

    my @genomes = $unvsub->get_genomes();

Return a list of the genome IDs for the subsystem. The genomes will be
presented in row order. In other words, if you index into the sixth
element of the return list, you will retrieve the ID of the genome for
the sixth row.

=cut

sub get_genomes {
    my($self) = @_;
    # The genome array is a list of 2-tuples. We extract the first
    # element of each tuple, which is the genome ID.
    return map { $_->[0] } @{$self->{Genomes}};
}

=head3 get_variant_code

    my $code = $unvsub->get_variant_code($genome);

Return the variant code for a genome. Each subsystem has several variations.
The variant code indicates which variation of a subsystem is used by a
particular genome.

Genome data is stored in a list of 2-tuples. The first element is the genome
ID; the second is the variant code.

=over 4

=item genome

ID or row index of the genome whose variant code is desired.

=item RETURN

Returns the variant code for the specified genome.

=back

=cut

sub get_variant_code {
    my($self,$genome) = @_;
    # Check to see if we have a row index.
    if ($genome =~ /^\d+$/)
    {
        # Here we have a row index, so use it to find the genome's variant
        # code.
        return $self->{Genomes}->[$genome]->[1];
    }
    else
    {
        # Here we have a genome ID, so we need to convert it to a row index.
        my $genomeI = $self->{GenomeIndex}->{$genome};
        return $self->{Genomes}->[$genomeI]->[1];
    }
}

=head3 get_pegs_from_cell

    my @pegs = $unvsub->get_pegs_from_cell($genome, $role);

Return a list of the features in a specified spreadsheet cell. The cell is specified
by genome ID and role ID.

=over 4

=item genome

ID of the genome relevant to the cell.

=item role

ID of the role relevant to the cell.

=item RETURN

Returns a list of the features in the cell, or an empty list if the cell is empty.

=back

=cut

sub get_pegs_from_cell {
    my($self,$genome,$role) = @_;
    # Convert the genome and role IDs to row and column indices.
    my $genomeI = $self->{GenomeIndex}->{$genome};
    my $roleI   = $self->{RoleIndex}->{$role};
    # Get the pegs from the cell and return them.
    my $pegs    = $self->{PegHash}->{$genomeI}->{$roleI};
    return $pegs ? @$pegs : ();
}

=head3 get_diagrams

    my @list = $sub->get_diagrams();

Return a list of the diagrams associated with this subsystem. Each diagram
is represented in the return list as a 4-tuple C<[diagram_id, diagram_name,
page_link, img_link]> where

=over 4

=item diagram_id

ID code for this diagram.

=item diagram_name

Displayable name of the diagram.

=item page_link

URL of an HTML page containing information about the diagram.

=item img_link

URL of an HTML page containing an image for the diagram.

=back

Note that the URLs are in fact for CGI scripts with parameters that point them
to the correct place. Though Sprout has diagram information in it, it has
no relationship to the diagrams displayed in SEED, so the work is done entirely
on the SEED side.

=cut

sub get_diagrams {
    # Get the parameters.
    my ($self) = @_;
    # Return the diagram list.
    return @{$self->{Diagrams}};
}

sub get_notes {
    my($self) = @_;

    return $self->{Notes};
}

sub get_role_index {
    my($self,$role) = @_;

    return $self->{RoleIndex}->{$role};
}

sub get_role_abbr {
    my($self,$roleI) = @_;

    if ($roleI !~ /^\d+$/)
    {
        $roleI = $self->{RoleIndex}->{$roleI};
    }
    my $roles = $self->{Roles};
    return $roles->[$roleI]->[1];
}

sub get_reactions {
    my($self) = @_;

    return $self->{Reactions};
}

sub get_subset_namesC {
    my($self) = @_;

    return map { $_->[0] } @{$self->{RoleSubsets}};
}

sub get_subsetC_roles {
    my($self,$subset) = @_;
    my($i,$j);

    my $subset_info = $self->{RoleSubsets};
    for ($i=0; ($i < @$subset_info) && ($subset_info->[$i]->[0] ne $subset); $i++) {}
    if ($i < @$subset_info)
    {
        my @roles = ();
        foreach $j (@{$subset_info->[$i]->[1]})
        {
            push(@roles,$self->{Roles}->[$j]->[0]);
        }
        return @roles;
    }
    return undef;
}

sub get_color_of {
    my($self,$peg)  = @_;

    return $self->{Colors}->{$peg};
}

=head3 active_genomes

    my $activeHash = UnvSubsys::active_genomes(\@row_subsets, $active_subsetR, $focus, \%genomeH, \@genomes_info);

Return a hash containing the active genomes for this subsystem display. The
keys of the hash will be the row indices of the genomes to be highlighted on the
display. Each genome ID will map to 1. Thus, if C<< $activeHash->{3} >>
tests TRUE, the fourth row should be highlighted.

The rules for determining the active genomes are as follows. If I<$active_subsetR> is
specified, it is presumed to be the ID of the subset containing the active genomes.
If it is not specified and I<$focus> is specified, then I<$focus> is presumed to be the
ID of the genome currently in focus, and the active genomes will be the ones in the
smallest row subset containing the genome in focus. If neither I<$active_subsetR> nor
I<$focus> are specified, then all genomes are active.

=over 4

=item row_subsets

Reference to a list of 2-tuples. Each tuple consists of a row subset ID followed by
a reference to a list of the row indices for the rows in the identified subset.

=item active_subsetR (optional)

ID of the active subset (if any).

=item focus

ID of the genome currently in focus (if any). If there is no active subset, then
the smallest subset containing this genome will be made active.

=item genomeH

Reference to a hash of genome IDs to row indices. The keys of this hash are
the genomes in this subsystem, which also form the subsystem spreadsheet's
rows.

=item genomes_info

Reference to a list of 2-tuples. The first element of each 2-tuple is the
ID of a genome; the second is the variant code for the subsystem variant
used by the genome. The tuples are ordered by row index, so that the ID
and variant code of the genome in a particular row can be located by indexing
into this parameter using the subsystem spreadsheet row number.

=item RETURN

Returns a reference to a hash that maps the row indices of the active genomes
to 1. This hash can be used to quickly determine whether or not a particular
row is to be highlighted.

=back

=cut

sub active_genomes {
    # Get the parameters.
    my($fig, $row_subsets, $active_subsetR, $focus, $genomeH, $genomes_info) = @_;
    my($i,@bestL);
    # Declare the return variable.
    my $active_genomes = {};
    # Check for an active subset.
    if ($active_subsetR)
    {
        # Search for the active subset in the row subset array.
        for ($i=0; ($i < @$row_subsets) && ($row_subsets->[$i]->[0] ne $active_subsetR); $i++) {}
        if ($i < @$row_subsets)
        {
            # Here we found the active subset, so we extract its list of row indices.
            @bestL = @{$row_subsets->[$i]->[1]};
        }
        else {
            # Here we have the ID of the active subset. First, we search for that ID
            # in the row subset list.
            for ($i=0; ($i < @$row_subsets) && ($row_subsets->[$i]->[0] ne $active_subsetR); $i++) {}
            if ($i < @$row_subsets)
            {
                # Here we found the named subset, so we return its member list.
                @bestL = @{$row_subsets->[$i]->[1]};
            }
            else
            {
                # Here the active subset does not exist. We punt by extracting a
                # list of all the row indices in the spreadsheet.
                $active_subsetR = 'All';
                @bestL = map { $genomeH->{$_} } keys(%$genomeH);
            }
        }
    }
    elsif ($focus)
    {
        # Here we don't have an active row subset, but a particular genome is in
        # focus. We'll look for the smallest row subset containing the genome
        # in focus. First, we need to prime the loop. "$bestN" will be the ID
        # of the best subset found so far; "@bestL" is where we stash the list
        # of IDs in the subset. Our initial selection, then, will be the
        # fake "All" subset, which contains the entire collection of rows.

        if (! $fig->is_complete($focus))
        {
            # Here the gnome in focus is incomplete, so it won't be anywhere
            # in our list. We default to making everything active.
            $active_subsetR = 'All';
            @bestL = map { $genomeH->{$_} } keys(%$genomeH);
        } else {
            my $bestN   = "All";
            @bestL   = map { $genomeH->{$_} } keys(%$genomeH);
            # Next, we get the row index for the genome in focus.
            my $focusIndex = $genomeH->{$focus};
            # Now we loop through all the row subsets.
            my $tuple;
            foreach $tuple (@$row_subsets)
            {
                # Extract the subset ID and its list of row indices. The latter is
                # in "$genomeIs".
                my($id,$genomeIs) = @$tuple;
                # Search for the index of the focus row in the current subset's list
                # of row indices.
                for ($i=0; ($i < @$genomeIs) && ($genomeIs->[$i] != $focusIndex); $i++) {}
                # Now either $i will be the index of the focus row in the subset, or
                # it is going to be equal to the number of rows in the subset.
                if ($i < @$genomeIs)
                {
                    # We've found the focus row in this subset. Select it as the new
                    # best subset if it's smaller than the last one we found.
                    if (@$genomeIs < @bestL)
                    {
                        $bestN  = $id;
                        @bestL = @$genomeIs;
                    }
                }
            }
            # Save the best subset found as the active one.
            $active_subsetR = $bestN;
        }
    } else {
        # Here we have nothing: no active subset, and no focus row. We make
        # all rows active.
        $active_subsetR = 'All';
        @bestL = map { $genomeH->{$_} } keys(%$genomeH);
    }
    # "@bestL" now contains a list of the row indices for the active genomes.
    # We convert it to a hash and return it.
    my %active_genomes = map { $_ => 1 } @bestL;
    return \%active_genomes;
}

1;
