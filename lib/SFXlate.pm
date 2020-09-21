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

package SFXlate;

use strict;

use Carp;

#
# Conditonally import sprout stuff.
#
BEGIN {
    eval {

        require Sprout;
        import Sprout;
        require SproutSubsys;
        import SproutSubsys;
        require CustomAttributes;
        import CustomAttributes;
        require RemoteCustomAttributes;
        import RemoteCustomAttributes;
    };
}

use Data::Dumper;
use FIG_Config;
use SeedDas;
use Tracer;
use FIG;
use FIGRules;
use BasicLocation;
use FullLocation;

=head1 Sprout/FIG API Shim

=head2 Introduction

This object accepts calls using the standard FIG method signatures and
translates them into Sprout calls. The end result is that an SFXlate
object can be used in place of a standard FIG object in legacy code.
provided that no unsupported functions are used.

=cut

#: Constructor SFXlate->new();

=head2 Constructors

=head3 new

    my $sfxlate = SFXlate->new($fig, $sproutDB, $sproutData);

This is the constructor for a full-service SFXlate object.

=over 4

=item fig

This is a legacy parameter that is effectively ignored.

=item sproutDB

Name of the sprout database. If undefined, the B<sproutDB> value from
the FIG configuration module will be used.

=item sproutData

Name of the directory containing the sprout data files. If undefined, the
B<sproutData> value from the FIG configuration module will be used.

=item sproutDBD

Name of the DBD file defining the sprout database.

=back

=cut

sub new {
    my($class, $fig, $sproutDB, $sproutData, $sproutDBD) = @_;
    my $sprout = SFXlate->new_sprout_only($sproutDB, $sproutData, $sproutDBD);
    my $self = {
        sprout => $sprout,
        fig => $fig,
        ca => undef,
    };
    return bless $self, $class;
}

=head3 new_sprout_only

    my $sprout = SFXlate->new_sprout_only($sproutDB, $sproutData, $xmlFile, $noOpen);

This is a special constructor that returns a pure Sprout object using
the FIG configuration defaults. Note that the Sprout database has a completely separate
set of configuration parameters from the SEED database. Thus, C<$FIG_Config::sproutUser>
is the user name for Sprout while C<$FIG_Config::user> is the user name for SEED.

=over 4

=item sproutDB

Name of the sprout database. If undefined, the B<sproutDB> value from
the FIG configuration module will be used.

=item sproutData

Name of the directory containing the sprout data files. If undefined, the
B<sproutData> value from the FIG configuration module will be used.

=item xmlFile

Name of the XML file containing the database definition. If undefined,
the file C<SproutDBD.xml> will be used, either from the main FIG
directory or the Sprout data directory.

=item noOpen

If TRUE, the database will not be opened. The default is FALSE.

=back

=cut

sub new_sprout_only {
    my($class, $sproutDB, $sproutData, $xmlFile, $noOpen) = @_;
    $sproutDB = $FIG_Config::sproutDB if !defined $sproutDB;
    $sproutData = $FIG_Config::sproutData if !defined $sproutData;
    Trace("Using sproutDB=$sproutDB sproutData=$sproutData") if T(1);
    my $openFlag = ($noOpen ? 1 : 0);
    if (! $xmlFile) {
        # Compute the DBD directory.
        my $dbd_dir = (defined($FIG_Config::dbd_dir) ? $FIG_Config::dbd_dir :
                                                      $FIG_Config::fig );
        $xmlFile = "$dbd_dir/SproutDBD.xml";
    }
    my $sprout = new Sprout(dbName => $sproutDB, dbd => $xmlFile,
                            options => { port => $FIG_Config::sproutPort,
                                         dbType => $FIG_Config::sproutDbms,
                                         dataDir => $sproutData,
                                         userData => "$FIG_Config::sproutUser/$FIG_Config::sproutPass",
                                         noDBOpen => $openFlag,
                                         sock => $FIG_Config::sproutSock,
                                         host => ($FIG_Config::sprout_host ||
                                                  $FIG_Config::dbhost),
    });
    return $sprout;
}

=head3 old_sprout_only

    my $sprout = SFXlate->old_sprout_only();

Open a Sprout object for the old NMPDR database. This simply calls L</new_sprout_only>
with the correct parameters required to access the old database.

=cut

sub old_sprout_only {
    # Get the parameters.
    my ($class) = @_;
    # Create and return the sprout object.
    return new_sprout_only($class, $FIG_Config::oldSproutDB, undef, $FIG_Config::oldSproutDBD);
}

=head3 all_usable_subsystems

    my @subs = $fig->all_usable_subsystems();

Return a list of the available subsystems. The FIG version of this method
skips subsystems that are experimental, deleted, or cluster-based;
however, these subsystems are filtered out before the data is loaded into
Sprout, so we simply return the complete list.

=cut

sub all_usable_subsystems {
    # Get the parameters.
    my ($self) = @_;
    # Declare the return variable.
    my @retVal = $self->{sprout}->GetFlat(['Subsystem'], "", [], 'Subsystem(id)');
    # Return the result.
    return @retVal;
}

=head3 new_sprout

    my $sfxlate = SFXlate->new_sprout($fig, $sprout);

This constructor creates an SFXlate object from pre-built FIG and
Sprout objects. Use this constructor if you have a pre-built Sprout
object and wish to pass it in to a module that is expecting a FIG
object.

=over 4

=item fig

This is a legacy parameter that is effectively ignored.

=item sprout

Sprout object to be used for Sprout calls.

=back

=cut

sub new_sprout {
    my($class, $fig, $sprout) = @_;
    my $self = {
        sprout => $sprout,
        ca => undef,
    };

    return bless $self, $class;
}

=head3 get_db

    my $dbObject = SFXLate::get_db($modName);

This method returns an ERDB object for the specified database. The database name
should be C<Sprout> for the Sprout database and C<CustomAttribute> for the custom
attributes database. This method needs to be made more robust; right now it's the
minimum needed for creating database test pages.

=over 4

=item modName

Name of the database to be loaded into the ERDB object.

=item RETURN

Returns an ERDB object for the default database with the specified name.

=back

=cut

sub get_db {
    # Get the parameters.
    my ($modName) = @_;
    # Declare the return variable.
    my $retVal;
    # Process according to the module name passed in.
    if ($modName eq "Sprout") {
        Trace("Selecting Sprout module.") if T(ERDB => 2);
        # Create a sprout object for the default database.
        $retVal = SFXlate->new_sprout_only();
    } elsif ($modName eq "CustomAttributes") {
        Trace("Selecting Custom Attributes.") if T(ERDB => 2);
        $retVal = CustomAttributes->new();
    } else {
        Confess("Invalid module name \"$modName\" specified in get_db.");
    }
    return $retVal;
}

=head2 Public Methods

=head3 minmax

    my ($lo, $hi) = SFXlate::minmax($a, $b);

Return a list containing the incoming pair of arguments sorted
numerically.

=over 4

=item a

First parameter.

=item b

Second parameter.

=item RETURN

Returns the incoming pair of arguments in numeric order.

=back

=cut

sub minmax {
    # Get the parameters.
    my ($a, $b) = @_;
    # Declare the return variables.
    my ($lo, $hi) = ($a <= $b ? ($a, $b) : ($b, $a));
    # Return the results.
    return ($lo, $hi);
}

=head3 location_overlap

    my $overlap = $fig->location_overlap($loc1, $loc2);

or

    my $overlap = SFXlate::location_overlap($loc1, $loc2);

Return the number of overlapping nucleotides between 2 locations.

=over 4

=item loc1

A location string for the first location.

=item loc2

A location string for the second location.

=item RETURN

Returns the number of overlapping base pairs, or zero if there is no
overlap.

=back

=cut

sub location_overlap {
    # Get the parameters.
    shift if UNIVERSAL::isa($_[0], __PACKAGE__);
    my ($loc1, $loc2) = @_;
    # Convert the locations to location objects.
    my $loc1Object = BasicLocation->new($loc1);
    my $loc2Object = BasicLocation->new($loc2);
    # Compute the overlap.
    my $retVal = $loc1Object->Overlap($loc2Object->Left, $loc2Object->Right);
    # Return the result.
    return $retVal;
}


=head3 to_alias

    my @aliases = $sfxlate->to_alias($featureID, $type);

or

    my $alias = $sfxlate->to_alias($featureID, %type);

Return all of a feature's aliases of the specified type. This is a fairly simplistic
method: the alias's type is the string preceding the vertical bar (e.g. C<uni|12345> is
a UniProt alias).

=over 4

=item featureID

ID of the feature whose aliases are desired.

=item type

Type of aliases that are desired.

=item RETURN

In a list context, returns a list of the aliases of the specified type for the specified
feature ID. In a scalar context, returns the first alias of the specified type.

=back

=cut
#: Return Type @;
sub to_alias {
    # Get the parameters.
    my ($self, $featureID, $type) = @_;
    # Get the desired aliases.
    my @retVal = grep { $_ =~ /^$type\|/ } $self->{sprout}->FeatureAliases($featureID);
    # Return the result.
    return @retVal;
}

=head3 genome_list

    my $genomeData = $fig->genome_list();

Return a reference to a list of all the genomes. For each genome there will be a
list reference that contains the genome ID, the genome name, and the genome
domain.

=cut

sub genome_list {
    # Get the parameters.
    my ($self) = @_;
    # Declare the return variable.
    my @retVal = ();
    # Get the genome data.
    my @genomeData = $self->{sprout}->GetAll([qw(Genome)], "", [],
                                             ['Genome(id)',
                                              'Genome(genus)',
                                              'Genome(species)',
                                              'Genome(unique-characterization)',
                                              'Genome(taxonomy)']);
    # Reformat it. We combine the genus, species, and unique-characterization to form the name,
    # and strip off the first taxonomy word to get the domain.
    for my $genomeDatum (@genomeData) {
        my ($id, $genus, $species, $strain, $taxonomy) = @{$genomeDatum};
        my ($domain) = split(/\s*;\s+/, $taxonomy, 2);
        push @retVal, [$id, "$genus $species $strain", $domain];
    }
    # Return the result.
    return \@retVal;
}

=head3 all_compounds

    my @compounds = $sfx->all_compounds();

Return a list containing all of the KEGG compounds.

=cut

sub all_compounds {
    # Get the parameters.
    my ($self) = @_;
    # Get all the compound IDs.
    my @retVal = $self->{sprout}->GetFlat(['Compound'], "", [], 'Compound(id)');
    # Return them to the caller.
    return @retVal;
}

=head3 names_of_compound

    my @names = $sfx->names_of_compound($cid);

Returns a list containing all of the names assigned to the specified KEGG compound. The list
will be ordered with the primary name first.

=over 4

=item cid

ID of the desired compound.

=item RETURN

Returns a list of names for the specified compound.

=back

=cut

sub names_of_compound {
    # Get the parameters.
    my($self, $cid) = @_;
    # Get the specified compound's list of names.
    my @names = $self->{sprout}->GetFlat('HasCompoundName',
                                         "HasCompoundName(from-link) = ?",
                                         [$cid], 'to-link');
    # Put the main label at the front.
    my ($label) = $self->{sprout}->GetEntityValues(Compound => $cid, ['label']);
    my @retVal = grep { $_ ne $label } @names;
    unshift @retVal, $label;
    # Return the result to the caller.
    return @retVal;
}

=head3 comp2react

    my @rids = $sfx->comp2react($cid);

Return a list containing all of the reaction IDs for reactions that take the
specified compound as either a substrate or a product.

=cut

sub comp2react {
    # Get the parameters.
    my ($self, $cid) = @_;
    # Get a list of the reactions connected to this compound.
    my @retVal = $self->{sprout}->GetFlat(['IsAComponentOf'],
                                          'IsAComponentOf(from-link) = ?',
                                          [$cid], 'IsAComponentOf(to-link)');
    return @retVal;
}

=head3 valid_reaction_id

    my $flag = $sfx->valid_reaction_id($rid);

Returns true iff the specified ID is a valid reaction ID.

This will become important as we include non-KEGG reactions

=over 4

=item rid

Reaction ID to test.

=item RETURN

Returns TRUE if the reaction ID is in the data store, else FALSE.

=back

=cut

sub valid_reaction_id {
    # Get the parameters.
    my ($self, $rid) = @_;
    # Check to see if a reaction with the specified ID exists.
    my $retVal = $self->{sprout}->Exists('Reaction', $rid);
    return $retVal;
}

=head3 cas

    my $cas = $sfx->cas($cid);

Return the Chemical Abstract Service (CAS) ID for the compound, if known.

=over 4

=item cid

ID of the compound whose CAS ID is desired.

=item RETURN

Returns the CAS ID of the specified compound, or an empty string if the CAS ID
is not known or does not exist.

=back

=cut

sub cas {
    # Get the parameters.
    my ($self, $cid) = @_;
    # Ask for the CAS ID.
    my ($retVal) = $self->{sprout}->GetFlat(['Compound'], 'Compound(id) = ?',
                                            [$cid], 'Compound(cas-id)');
    # If we didn't find a CAS ID, return an empty string.
    if (! $retVal) {
        $retVal = "";
    }
    return $retVal;
}

=head3 cas_to_cid

    my $cid = $sfx->cas_to_cid($cas);

Return the compound id (cid), given the Chemical Abstract Service (CAS) ID.

=over 4

=item cas

CAS ID of the desired compound.

=item RETURN

Returns the ID of the compound corresponding to the specified CAS ID, or an empty
string if the CAS ID is not in the data store.

=back

=cut

sub cas_to_cid {
    # Get the parameters.
    my ($self, $cas) = @_;
    # Look for the compound.
    my @retVal = $self->{sprout}->GetFlat(['Compound'], "Compound(cas-id) = ?",
                                            [$cas], 'Compound(id)');
    # Return an empty string if the compound was not found or there were too
    # many.
    my $retVal = (@retVal != 1 ? "" : $retVal[0]);
    return $retVal;
}

=head3 all_reactions

    my @rids = $sfx->all_reactions();

Return a list containing all of the KEGG reaction IDs.

=cut

sub all_reactions {
    # Get the parameters.
    my ($self) = @_;
    # Get a list of reaction IDs.
    my @retVal = $self->{sprout}->GetFlat(['Reaction'], "", [], 'Reaction(id)');
    return @retVal;
}

=head3 delete_genomes

    my $stats = $sfx->delete_genomes(\@genomes);

Delete the specified genomes from the database.

=over 4

=item genomes

Reference to a list of the IDs of the genomes to be deleted.

=item RETURN

Returns a statistics object detailing the number of rows deleted from each table.

=back

=cut
#: Return Type $%;
sub delete_genomes {
    # Get the parameters.
    my ($self, $genomes) = @_;
    # Create the statistics object.
    my $retVal = Stats->new('genomeIDs');
    # Loop through the genome IDs, deleting them individually.
    for my $genomeID (@{$genomes}) {
        # Delete the genome.
        my $stats = $self->{sprout}->DeleteGenome($genomeID);
        # Accumulate the statistics.
        $retVal->Accumulate($stats);
        # Denote we've handled another genome.
        $retVal->Add('genomeIDs');
    }
    # Return the result.
    return $retVal;
}

=head3 feature_aliases_bulk

    my $aliasHash = $fig->feature_aliases_bulk(\@fids);

Return a reference to a hash mapping feature IDs to aliases. The aliases
are retrieved using a single query, rather than one at a time.
The FIG version of this method has an additional parameter to suppress
the check for deleted features, but this is not necessary in the Sprout
database, because deleted features are truly deleted instead of being
marked inactive.

=over 4

=item fids

Reference to a list of feature IDs. The aliases returned will all belong to
the specified features.

=item RETURN

Returns a reference to a hash mapping each feature ID to a list of its aliases.

=back

=cut

sub feature_aliases_bulk {
    # Get the parameters.
    my ($self, $fids) = @_;
    # Create a filter for the IDs in the list.
    my @filterMarks = ();
    my @filterParms = ();
    for my $fid (@{$fids}) {
        push @filterMarks, '?';
        push @filterParms, $fid;
    }
    my $filterString = "IsAliasOf(to-link) IN (" . join(", ", @filterMarks) . ")";
    # Get the aliases.
    my @rows = $self->{sprout}->GetAll(['IsAliasOf'], $filterString, \@filterParms, [qw(IsAliasOf(to-link) IsAliasOf(from-link))]);
    # Form them into a hash.
    my %retVal = ();
    for my $row (@rows) {
        push @{$retVal{$row->[0]}}, $row->[1];
    }
    # Return the result.
    return \%retVal;
}

=head3 reversible

    my $flag = $sfx->reversible($rid);

Return TRUE if the specified reaction is reversible. A reversible reaction has no main
direction. The connector is symbolized by C<< <=> >> instead of C<< => >>.

=over 4

=item rid

ID of the ralevant reaction.

=item RETURN

Returns TRUE if the specified reaction is reversible, else FALSE.

=back

=cut

sub reversible {
    # Get the parameters.
    my ($self, $rid) = @_;
    # Assume a reversible reaction unless we prove otherwise.
    my $retVal = 1;
    # Look for the reaction's reversibility flag.
    my ($reversible) = $self->{sprout}->GetFlat(['Reaction'], "Reaction(id) = ?", [$rid],
                                                'Reaction(rev)');
    if (defined $reversible) {
        $retVal = $reversible;
    }
    # Return the result.
    return $retVal;
}

=head3 contig_lengths

    my $contigHash = $fig->contig_lengths($genomeID);

Return a hash reference that maps each contig ID to a length that
indicates the number of base pairs.

=over 4

=item genomeID

ID of the genome whose contigs are desired.

=item RETURN

Returns a hash that maps the IDs of each of the genome's contigs to their lengths in base pairs.

=back

=cut

sub contig_lengths {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    # Declare the return variable.
    my $retVal = {};
    # Get a list of the contigs for this genome.
    my @contigData = $self->{sprout}->GetAll([qw(HasContig IsMadeUpOf)], "HasContig(from-link) = ?", [$genomeID],
                                             [qw(HasContig(to-link) IsMadeUpOf(start-position) IsMadeUpOf(len))]);
    # Loop through the contigs.
    for my $contigRow (@contigData) {
        # Break out the data items in the row.
        my ($contigID, $start, $len) = @{$contigRow};
        Trace("Contig $contigID has start $start and length $len.") if T(3);
        # Compute the position after the end of this sequence. If it's the last sequence,
        # that will be the contig length.
        my $end = $start + $len;
        # If this is the best value so far, save it.
        if (! exists $retVal->{$contigID} || $end > $retVal->{$contigID}) {
            $retVal->{$contigID} = $end;
        }
    }
    # Return the result.
    return $retVal;
}


=head3 reaction2comp

    my @tuples = $fig->reaction2comp($rid, $which);

Return the substrates or products for a reaction.  In any event (i.e.,
whether you ask for substrates or products), you get back a list of
3-tuples.  Each 3-tuple will contain

    [$cid,$stoich,$main]

Stoichiometry indicates how many copies of the compound participate in
the reaction. It is normally numeric, but can be things like "n" or "(n+1)".
$main is 1 iff the compound is considered "main" or "connectable".

=over 4

=item rid

ID of the raction whose compounds are desired.

=item which

TRUE if the products (right side) should be returned, FALSE if the substrates
(left side) should be returned.

=item RETURN

Returns a list of 3-tuples. Each tuple contains the ID of a compound, its
stoichiometry, and a flag that is TRUE if the compound is one of the main
participants in the reaction.

=back

=cut

sub reaction2comp {
    # Get the parameters.
    my ($self, $rid, $which) = @_;
    # Convert the which flag into a set number.
    my $setN = ($which ? 1 : 0);
    # Get the reaction's compounds from the database. We sort them by component ID for
    # compatability with the SEED method.
    my @retVal = $self->{sprout}->GetAll(['IsAComponentOf'],
                                         "IsAComponentOf(to-link) = ? AND IsAComponentOf(product) = ? ORDER BY IsAComponentOf(from-link)",
                                         [$rid, $setN], ['IsAComponentOf(from-link)',
                                                         'IsAComponentOf(stoichiometry)',
                                                         'IsAComponentOf(main)']);
    # Return the result.
    return @retVal;
}

=head3 catalyzed_by

    my @ecs = $sfx->catalyzed_by($rid);

Return the ECs (roles) that are reputed to catalyze the reaction.  Note that we are currently
just returning the ECs that KEGG gives.

=over 4

=item rid

ID of the reaction whose catalyzing roles are desired.

=item RETURN

Returns the EC codes of the roles that catalyze the reaction.

=back

=cut

sub catalyzed_by {
    # Get the parameters.
    my($self, $rid) = @_;
    # Get the roles.
    my @retVal = $self->{sprout}->GetFlat(['Catalyzes', 'Role'], "Catalyzes(to-link) = ?", [$rid],
                                          'Role(EC)');
    # Return them to the caller.
    return @retVal;
}

=head3 catalyzes

    my @ecs = $fig->catalyzes($role);

Returns the reaction IDs of the reactions catalyzed by the specified role (normally an EC).

=over 4

=item role

ID or EC number of the role whose reactions are desired.

=item RETURN

Returns a list containing the IDs of the reactions catalyzed by the role.

=back

=cut

sub catalyzes {
    # Get the parameters.
    my ($self, $role) = @_;
    # Check the type of role specifier.
    my $key = (FIG::is_ec($role) ? "Role(EC)" : "Role(id)");
    # Look for a list of reactions corresponding to the role.
    my @retVal = $self->{sprout}->GetFlat(['Role', 'Catalyzes'], "$key = ?",
                                          [$role], 'Catalyzes(to-link)');
    return @retVal;
}

=head3 displayable_reaction

    my $displayString = $fig->displayable_reaction($rid);

Return the displayable version of a reaction. This is built on the fly from
the B<IsAComponentOf> relationship.

=cut

sub displayable_reaction {
    # Get the parameters.
    my($self, $rid) = @_;
    # Declare the return variable.
    my $retVal = "";
    # Get the reaction's connector type.
    my ($connector) = $self->{sprout}->GetEntityValues('Reaction', $rid, ['Reaction(rev)']);
    if (! defined $connector) {
        # Here the reaction does not exist. We return the ID unmodified.
        $retVal = $rid;
    } else {
        # Determine the connector style. TRUE means the reaction is
        # reversible.
        $connector = ($connector ? "<=>" : "=>");
        # The reaction display consists of the substrate compounds, the
        # connector, and then the product compounds. First, we need the
        # data.
        my @compounds = $self->{sprout}->GetAll(['IsAComponentOf', 'Compound'],
                                                "IsAComponentOf(to-link) = ? ORDER BY IsAComponentOf(product), IsAComponentOf(loc), IsAComponentOf(main) DESC",
                                                [$rid], ['IsAComponentOf(product)',
                                                         'IsAComponentOf(stoichiometry)',
                                                         'Compound(label)']);
        # Each compound is prefixed by an operator. The first substrate
        # has nothing in front of it. The first product has a connector
        # in front of it. Everything else is preceded by a plus sign.
        # We use a saved mode indicator to detect the changes, and a
        # simple array to hold the two types of special operators.
        my $thisMode = -1;
        my @op = ("", " $connector ");
        # Loop through the compounds.
        for my $compoundData (@compounds) {
            # Split out the compound data.
            my ($modeFlag, $stoich, $name) = @{$compoundData};
            # Determine the stoichiometry prefix. If it's 1, there is
            # no prefix. Otherwise, we use the stoichiometry value
            # followed by a space. In other words, two DNAs is "2 DNA",
            # but one DNA is simply "DNA".
            my $prefix = ($stoich == 1 ? "" : "$stoich ");
            # Determine the operator to go in front of this compound.
            if ($thisMode eq $modeFlag) {
                $retVal .= " + ";
            } else {
                # Here the mode has changed, so we pull the appropriate
                # operator out of the @op array. Note that we don't
                # surround it by spaces the way we do with the plus.
                # The empty string is supposed to be empty, and the
                # connector has already had the spaces put in.
                $retVal .= $op[$modeFlag];
                $thisMode = $modeFlag;
            }
            # Add this compound.
            $retVal .= $prefix . $name;
        }
    }
    # Return the result string.
    return $retVal;
}

=head3 all_maps

    my @maps = $fig->all_maps();

Return all of the KEGG maps in the data store. KEGG maps in Sprout are
represented by the B<Diagram> entity.

=cut

sub all_maps {
    # Get the parameters.
    my ($self) = @_;
    # Get all of the diagram IDs.
    my @retVal = $self->{sprout}->GetFlat(['Diagram'], "", [], 'Diagram(id)');
    # Return the list of IDs.
    return @retVal;
}

=head3 ec_to_maps

    my @maps = $fig->ec_to_maps($ec);

Return the set of maps that contain a specific functional role. The role can be
specified by an EC number or a full-blown role ID. Maps in Sprout are stored as
diagrams.

=over 4

=item ec

The EC number or role ID of the role whose maps are desired.

=item RETURN

Returns a list of the IDs for the maps that contain the specified role.

=back

=cut

sub ec_to_maps {
    # Get the parameters.
    my ($self, $ec) = @_;
    # Declare the return value.
    my @retVal;
    # Get the Sprout object.
    my $sprout = $self->{sprout};
    # Determine whether this is an EC number or a role.
    if (FIG::is_ec($ec)) {
        # Here we have an EC number. We determine the roles using the IsIdentifiedByEC
        # relationship.
        @retVal = $sprout->GetFlat([qw(IsIdentifiedByEC RoleOccursIn)], "IsIdentifiedByEC(to-link) = ?",
                                   [$ec], 'RoleOccursIn(to-link)');
    } else {
        # Here we have a role ID, which cuts out a step.
        @retVal = $sprout->GetFlat(['RoleOccursIn'], "RoleOccursIn(from-link) = ?", [$ec],
                                   'RoleOccursIn(to-link)');
    }
    if (T(ECLink => 4)) {
        my $count = @retVal;
        Trace("$count maps returned for $ec.");
    }
    return @retVal;
}

=head3 role_to_maps

This is an alternate name for L</ec_to_maps>.

=cut

sub role_to_maps {
    my ($self, $role) = @_;
    return $self->ec_to_maps($role);
}

=head3 in_pch_pin_with_and_evidence

    my @list = $fig->in_pch_pin_with_and_evidence($peg);

Return a list of the features pinned to the specified feature. Each
element of the list will be a 2-tuple, the first element being the ID of
a pinned feature and the second a flag that is TRUE if the feature has
physically close homologs in diverse organisms and 0 otherwise.

=over 4

=item peg

ID of the relevant feature.

=item RETURN

Returns a list of 2-tuples, each consisting of a pinned feature ID and an indicator of
whether or not the pin is conserved in any diverse organisms.

=back

=cut

sub in_pch_pin_with_and_evidence {
    # Get the parameters.
    my ($self, $peg) = @_;
    # Ask the coupling server for the data.
    my @retVal = FIGRules::NetCouplingData('in_pch_pin_with_and_evidence', id1 => $peg);
    # Return the result.
    return @retVal;
}


=head3 map_to_ecs

    my @ecs = $fig->map_to_ecs($map);

Return the set of functional roles (usually ECs) that are contained in the functionality
depicted by a map.

This method only returns EC numbers. If we need roles that don't have EC numbers, it
will need to be changed.

=over 4

=item map

ID of the KEGG map whose roles are desired.

=item RETURN

Returns a list of EC numbers for the roles in the specified map.

=back

=cut

sub map_to_ecs {
    # Get the parameters.
    my ($self, $map) = @_;
    # Get all the EC numbers for the roles in the specified map.
    my @retVal = $self->{sprout}->GetFlat(['Role', 'RoleOccursIn'], "RoleOccursIn(to-link) = ?",
                                          [$map], 'Role(EC)');
    return @retVal;
}

=head3 map_name

    my $name = $fig->map_name($map);

Return the descriptive name covering the functionality depicted by the specified map.

=over 4

=item map

ID of the map whose description is desired.

=item RETURN

Returns the descriptive name of the map, or an empty string if no description is available.

=back

=cut

sub map_name {
    # Get the parameters.
    my ($self, $map) = @_;
    # Ask for the map name.
    my ($retVal) = $self->{sprout}->GetFlat(['Diagram'], "Diagram(id) = ?", [$map],
                                            'Diagram(name)');
    # If the map was not found, return an empty string.
    if (! defined $retVal) {
        $retVal = "";
    }
    # Return the result.
    return $retVal;
}

=head3 abbrev

    my $abbreviated_name = $sfxlate->abbrev($genome_name);

Abbreviate a genome name to 10 characters or less.

For alignments and such, it is very useful to be able to produce an abbreviation of genus/species.
That's what this does.  Note that multiple genus/species might reduce to the same abbreviation, so
be careful (disambiguate them, if you must).

The abbreviation is formed from the first three letters of the species name followed by the
first three letters of the genus name followed by the first three letters of the species name and
then the next four nonblank characters.

=over 4

=item genome_name

The name to abbreviate.

=item RETURN

An abbreviated version of the specified name.

=back

=cut

sub abbrev {
    my ($self, $genome_name) = @_;
    return FIG::abbrev($genome_name);
}

=head3 add_attribute

    $sfxlate->add_attribute($peg, $key, $value, $url);

Add a new attribute value (Property) to a feature. In the SEED system, attributes can
be added to almost any object. In Sprout, they can only be added to features. In
Sprout, attributes are implemented using I<properties>. A property represents a key/value
pair. If the particular key/value pair coming in is not already in the database, a new
B<Property> record is created to hold it.

=over 4

=item peg

ID of the feature to which the attribute is to be replied.

=item key

Name of the attribute (key).

=item value

Value of the attribute.

=item url

URL or text citation from which the property was obtained.

=back

=cut

sub add_attribute {
    my($self, $peg, $key, $value, $url) = @_;
    $self->{sprout}->AddProperty($peg, $key, $value, $url);
}

=head3 get_system_name

    my $name = $sfxlate->get_system_name;

Returns C<sprout>, indicating that this is object is using the Sprout
database. The same method on a FIG object will return C<seed>.

=cut
#: Return Type $;
sub get_system_name {
    return "sprout";
}

=head3 beg_of

    my $offset = $sfxlate->beg_of($loc);

Return the beginning offset of the specified list of locations.

=over 4

=item loc

Space-delimited list of locations in the Sprout format. In other words,
each location must be of the form I<Contig>C<_>I<beg>C<+>I<len> or
I<Contig>C<_>I<beg>C<->I<len>, where I<Contig> is a contig ID, I<beg>
is the starting offset of the location, and I<len> is the length of
the location.

=item RETURN

Returns a number indicating the starting offset of the first location if it
is for a forward gene, or the starting offset of the last location
if it is for a backward gene.

=back

=cut
#: Return Type $;
sub beg_of {
    my($self, $loc) = @_;

    #
    # Loc is a space-separated list of spans.
    #

    my @spans = split(/\s+/, $loc);

    my $first = $spans[0];
    my $last = $spans[$#spans];

    if ($first =~ /_(\d+)\+\d+$/) {
        return $1;
    } elsif ($last =~ /_(\d+)-\d+$/) {
        return $1;
    } else {
        Confess("Bad beg_of for loc='$loc' spans='@spans'");
    }
}

=head3 end_of

    my $offset = $sfxlate->end_of($loc);

Return the ending offset of the specified list of locations.

=over 4

=item loc

Space-delimited list of locations in the Sprout format. In other words,
each location must be of the form I<Contig>C<_>I<beg>C<+>I<len> or
I<Contig>C<_>I<beg>C<->I<len>, where I<Contig> is a contig ID, I<beg>
is the starting offset of the location, and I<len> is the length of
the location.

=item RETURN

Returns a number indicating the ending offset of the first location if it
is for a backward gene, or the ending offset of the last location
if it is for a forward gene.

=back

=cut
#: Return Type $;
sub end_of {
    my($self, $loc) = @_;

    #
    # Loc is a space-separated list of spans.
    #

    my @spans = split(/\s+/, $loc);

    my $first = $spans[0];
    my $last = $spans[$#spans];

    if ($first =~ /_(\d+)-(\d+)$/) {
        return $1 - $2;
    } elsif ($last =~ /_(\d+)\+(\d+)$/) {
        return $1 + $2;
    } else {
        Confess("Bad end_of for \"$loc\".");
    }
}

=head3 contig_of

    my $contigID = $sfxlate->contig_of($loc);

Return the contig ID from a location string. The location must be of the
form I<Contig>C<_>I<beg>C<+>I<len> or I<Contig>C<_>I<beg>C<->I<len>,
where I<Contig> is a contig ID, I<beg> is the starting offset of the
location, and I<len> is the length of the location.

=cut
#: Return Type $;
sub contig_of {
    my($self, $loc) = @_;

    if ($loc =~ /^(\S+)_\d+[+-]\d+/) {
        return $1;
    } else {
        Confess("Bad contig_of for \"$loc\".");
    }
}

=head3 feature_aliases

    my @aliasList = $sfxlate->feature_aliases($feature);

Return a list of the aliases for the specified feature.

=over 4

=item feature

ID of the feature whose aliases are desired.

=item RETURN

Returns a list of the alias names for the specified feature.

=back

=cut
#: Return Type @;
sub feature_aliases {
    my($self, $feature) = @_;

    return $self->{sprout}->FeatureAliases($feature);
}


=head3 feature_aliases_in_tbl

    my @aliasList = $sfxlate->feature_aliases($feature);

Return a list of the aliases for the specified feature. This method is identical
to L</feature_aliases>. The FIG method with this name gets a smaller set of
aliases than the ones returned by its B<feature_aliases> method; however, the
aliases omitted are ones we don't have in the database.

=over 4

=item feature

ID of the feature whose aliases are desired.

=item RETURN

Returns a list of the alias names for the specified feature.

=back

=cut
#: Return Type @;
sub feature_aliases_in_tbl {
    my($self, $feature) = @_;
    return $self->{sprout}->FeatureAliases($feature);
}

=head3 feature_location

    my @locations = $sfxlate->feature_location($feature);

Return a list of the location descriptors for the specified feature.

=over 4

=item feature

ID of the feature whose locations are desired.

=item RETURN

Returns a list of the descriptors for the specified feature's
locations, in transcription order. In a scalar context, returns
the locations as a single comma-delimited string.

=back

=cut
#: Return Type @;
#: Return Type $;
sub feature_location {
    my($self, $feature) = @_;
    return $self->{sprout}->FeatureLocation($feature);
}

=head3 ftype

    my $ftype = $sfxlate->ftype($feature);

Return the type (peg, rna, etc.) of the specified feature.

=over 4

=item feature

ID of the feature whose type is desired.

=item RETURN

Returns the type of the specified feature.

=back

=cut
#: Return Type $;
sub ftype {
    my($self, $feature) = @_;

    return $self->{sprout}->FType($feature);
}

=head3 contig_ln

    my $length = $sfxlate->contig_ln($genome, $contig);

Return the length of the specified contig.

=over 4

=item genome

ID of the genome to which the contig belongs

=item contig

ID of the contig whose length is desired

=back

=cut
#: Return Type $;
sub contig_ln {
    my($self, $genome, $contig) = @_;
    return $self->{sprout}->ContigLength($contig);
}

=head3 compute_clusters

    my @clusterList = $fig->compute_clusters(\@pegList, $subsystem, $distance);

Partition a list of PEGs into sections that are clustered close together on
the genome. The PEGs must be from a single subsystem row that was recently
retrieved using the C<get_pegs_from_cell> method on the subsystem object
passed in. If this is not the case, the method will still work, but the
PEGs could be mapped into the incorrect clusters.

=over 4

=item pegList

Reference to a list of PEG IDs.

=item subsystem

Subsystem object for the relevant subsystem.

=item distance (optional)

The maximum distance between PEGs that makes them considered close. This
parameter is not used, but is required for compatability with SEED.

=item RETURN

Returns a list of lists. Each sub-list is a cluster of PEGs.

=back

=cut

sub compute_clusters {
    # Get the parameters.
    my ($self, $pegList, $subsystem, $distance) = @_;
    # Compute the result.
    my $retVal = $self->{sprout}->ClusterPEGs($subsystem, $pegList);
    # Return it to the caller.
    return @{$retVal};
}

=head3 maps_to_id

    my $peg = $sfxlate->maps_to_id($id);

The "major synonym" of a feature is one that is selected (and represents the
longest version of a set of essentially identical sequences).  This routine
returns the "major synonym".

=over 4

=item id

Feature ID whose major synonym is desired.

=item RETURN

Returns the major synonym corresponding to the named feature, or the
feature itself if there is no synonym for the named feature.

=back

=cut
#: Return Type $;
sub maps_to_id {
    my ($self, $id) = @_;
    return $self->{sprout}->GetSynonymGroup($id);
}

=head3 cgi_url

    my $url = $sfxlate->cgi_url;

Return the URL of the directory containing the CGI scripts.

=cut
#: Return Type $;
sub cgi_url {
    my($self) = @_;

    return FIG::cgi_url();
}

=head3 function_of

    my $function = $sfxlate->function_of($id, $user);

or

    my @functions = $sfxlate->function_of($id);

In a scalar context, returns the most recently-determined functional
assignment of a specified feature by a particular user. In a list
context, returns a list of 2-tuples, each consisting of a user ID
followed by a functional assighment by that user. In this case,
the list contains all the functional assignments for the feature.

=over 4

=item id

ID of the relevant feature.

=item user

ID of the user whose assignment is desired (scalar context only)

=item RETURN

Returns the most recent functional assignment by the given user in scalar
context, and a list of functional assignments in list context.

=back

=cut
#: Return Type $;
#: Return Type @@;
sub function_of {
    if (wantarray()) {
        my ($self, $id) = @_;
        my %mapping = $self->{sprout}->AllFunctionsOf($id);
        my @retVal = ();
        for my $user (sort keys %mapping) {
            push @retVal, [$user, $mapping{$user}];
        }
        return @retVal;
    } else {
        my($self, $id, $user) = @_;
        return $self->{sprout}->FunctionOf($id, $user);
    }
}



=head3 bbh_list

    my $bbhHash = $sfxlate->bbh_list($genome, \@features);

Return a hash mapping the features in a specified list to their bidirectional best hits
on a specified target genome.

=over 4

=item genomeID

ID of the genome from which the best hits should be taken.

=item featureList

List of the features whose best hits are desired.

=item RETURN

Returns a reference to a hash that maps the IDs of the incoming features to the IDs of
their best hits.

=back

=cut
#: Return Type $%;
sub bbh_list {
    my($self, $genome, $features) = @_;

    return $self->{sprout}->BBHList($genome, $features);
}

=head3 sprout

    my $sprout = $fig->sprout();

Return the embedded Sprout object.

=cut

sub sprout {
    # Get the parameters.
    my ($self) = @_;
    # Return the result.
    return $self->{sprout};
}

=head3 init_das

    my $das = $sfxlate->init_das($url, $dsn);

Create a DAS object for use by the GBrowse facility.

=over 4

=item url

URL of the DAS script.

=item dsn

ID of the relevant Genome. The separating dot may be replaced by an
underscore

=item RETURN

Returns a DAS object for browsing the specified genome using Sprout data.
If this installation is not configured for GBrowse, returns an undefined
value.

=back

=cut
#: Return Type %;
sub init_das {
    my($self, $url, $dsn) = @_;

    my $das_data_dir = "$FIG_Config::var/DAS";

    if (-d $das_data_dir) {
        return new SeedDas($self,$das_data_dir, $url, $dsn);
    } else {
        return undef;
    }
}

=head3 genomes

    my @genomeList = $sfxlate->genomes;

Return a list of the IDs of all the genomes in the system.

=cut
#: Return Type @;
sub genomes {
    my($self) = @_;

    return $self->{sprout}->Genomes();
}

=head3 genus_species

    my $infoString = $sfxlate->genus_species($genome);

Return the genus, species, and unique characterization for the specified
genome.

=over 4

=item genome

ID of the genome whose identifying information is desired.

=item RETURN

Returns the genus and species of the genome, with the unique characterization
(if any). If the genome does not exist, returns an undefined value.

=back

=cut
#: Return Type $;
sub genus_species {
    my($self, $genome) = @_;

    return $self->{sprout}->GenusSpecies($genome);
}

=head3 genome_counts

    my ($arch, $bact, $euk, $vir, $env, $unk) = $fig->genome_counts($complete);

Count the number of genomes in each domain. If I<$complete> is TRUE, only complete
genomes will be included in the counts.

=over 4

=item complete

TRUE if only complete genomes are to be counted, FALSE if all genomes are to be
counted

=item RETURN

A six-element list containing the number of genomes in each of six categories--
Archaea, Bacteria, Eukaryota, Viral, Environmental, and Unknown, respectively.

=back

=cut

sub genome_counts {
    # Get the parameters.
    my ($self, $complete) = @_;
    # Get the Sprout object.
    my $sprout = $self->{sprout};
    # Get the counts.
    my @counts = $sprout->GenomeCounts($complete);
    # Return the counts.
    return @counts;
}

=head3 dna_seq

    my $sequence = $sfxlate->dna_seq($genome, @locations);

Return the sequence represented by a list of locations. The locations
should be in the standard sprout form I<contigID>C<_>I<begin>I<dir>I<end>.

=over 4

=item genome

ID of the relevant genome.

=item location1, location2, ... locationN

List of locations to be included in the DNA sequence.

=item RETURN

Returns a string specifying the DNA nucleotides in the specified locations.

=back

=cut
#: Return Type $;
sub dna_seq {
    my($self, $genome, @locations) = @_;
    return $self->{sprout}->DNASeq(\@locations);
}

=head3 all_contigs

    my @contigs = $sfxlate->all_contigs($genome);

Return a list of the contigs that make up the specified genome.

=over 4

=item genome

ID of the genome whose contigs are desired.

=item RETURN

Returns a list of the IDs for the contigs in the genome.

=back

=cut
#: Return Type @;
sub all_contigs {
    my($self, $genome) = @_;

    return $self->{sprout}->AllContigs($genome);
}

=head3 contigs_of

    my @contig_ids = $fig->contigs_of($genome);

Returns a list of all of the contigs occurring in the designated genome.
This is in fact just an alternate way of calling L</all_contigs>.

=over 4

=item genome

ID of the genome whose contigs are desired.

=item RETURN

Returns a list of the IDs for the contigs occurring in the specified genome.

=back

=cut

sub contigs_of {
    my($self,$genome) = @_;
    return $self->all_contigs($genome);
}

=head3 get_genome_subsystem_count

    my $num_subsytems = $fig->get_genome_subsystem_count($genomeID);

Return the number of subsystems in which a genome participates.

=over 4

=item genomeID

ID of the relevant genome.

=item RETURN

Returns the number of subsystems that have this genome in them.

=back

=cut

sub get_genome_subsystem_count{
    # Get the parameters.
    my ($self, $genomeID) = @_;
    # Get the count of subsystems connected to this genome.
    my $retVal = $self->{sprout}->GetCount([qw(ParticipatesIn)], "ParticipatesIn(from-link) = ?", [$genomeID]);
    # Return the result.
    return $retVal;
}

=head3 all_subsystem_classifications

    my @classifications = $fig->all_subsystem_classifications();

Return a list of all the subsystem classifications. Each element in the
list will contain a main subsystem class and a basic subsystem class.
The resulting list enables us to determine easily what the three-level
subsystem tree would look like.

=cut

sub all_subsystem_classifications {
    # Get the parameters.
    my ($self) = @_;
    # Get a list of all the subsystem classifications. Each element of this
    # list will be a sub-list containing the two classifications as a string
    # with a splitter in it.
    my @classList = $self->{sprout}->GetFlat(['Subsystem'], '', [], 'Subsystem(classification)');
    # Merge out duplicates using a hash.
    my %classHash = map { $_ => 1 } @classList;
    # Convert the hash into the desired return list. Each key is the two class names
    # separated by a splitter string. We use "split" to make this a two-element list reference.
    my @retVal = map { [ split(/$FIG_Config::splitter/, $_)] } keys %classHash;
    # Return the result.
    return @retVal;
}

=head3 genes_in_region

    my @features = $sfxlate->genes_in_region($genome, $contig, $start, $end);

Return a list of the features that overlap the specified region in a contig.

=over 4

=item genome

ID of the genome containing the contig in question.

=item contig

ID of the contig containing the desired region.

=item start

Offset of the first nucleotide in the region.

=item end

Offset of the last nucleotide in the region.

=item RETURN

Returns a list of the IDs of the features that overlap the specified
region.

=back

=cut
#: Return Type @;
sub genes_in_region {
    my($self, $genome, $contig, $start, $end) = @_;
    my @results = $self->{sprout}->GenesInRegion($contig, $start, $end);
    return @results;
}

=head3 get_attribute

    my @values = $sfxlate->get_attribute($feature, $attr);

Return a list of the values for the named attribute of a specified
feature.

=over 4

=item feature

ID of the feature whose attribute values are desired.

=item attr

Name of the attribute whose values are desired.

=item RETURN

Returns a list of values for the named attribute. In most cases the list
will be a singleton, however some attributes (such as C<alias>) can have
multiple values in the list.

=back

=cut
#: Return Type @;
sub get_attribute {
    my($self, $feature, $attr) = @_;

    my $obj = $self->{sprout}->GetEntity("Feature", $feature);
    $obj or return undef;

    return $obj->Value("Feature($attr)");
}

=head3 in_cluster_with

    my @features = $sfxlate->in_cluster_with($peg);

Return a list of the features functionally coupled with the specified
feature. A feature is considered functionally coupled if it tends to reside
on the same chromosome as the specified feature.

=over 4

=item peg

ID of the features whose functional couplings are desired.

=item RETURN

Returns a list of the IDs for the functionally-couplped features.

=back

=cut
#: Return Type @;
sub in_cluster_with {
    my($self, $peg) = @_;
    my %coupleData = $self->{sprout}->CoupledFeatures($peg);
    return sort { FIG::by_fig_id($a,$b); } keys %coupleData;
}

=head3 add_annotation

    my $ok = $sfxlate->add_annotation($featureID, $user, $text);

Add an annotation to a feature. This method inserts new data into the
Sprout database.

=over 4

=item featureID

ID of the feature to annotate.'

=item user

ID of the user making the annotation.

=item text

Text of the annotation.

=item RETURN

Returns 1 if successful, 0 if an error occurred.

=back

=cut
#: Return Type $;
sub add_annotation {
    my ($self, $featureID, $user, $text) = @_;
    Trace("Adding annotation in SFXlate.") if T(Bruce => 4);
    my $timestamp = time;
    my $retVal = $self->{sprout}->Annotate($featureID, $timestamp, $user, $text);
    return $retVal;
}

=head3 boundaries_of

    my ($contig, $beg, $end) = $sfxlate->boundaries_of($locations);

Examine a list of locations and return a location that encompasses the
entire set. The location returned will be a list specifying the relevant
contig and the beginning and ending offsets.

=over 4

=item locations

A reference to a list of the desired locations or a string containing a
comma-delimited list of the locations.

=item RETURN

Returns a three-element list consisting of the contig, beginning offset,
and ending offset of a region containing all of the specified locations.

=back

=cut
#: Return Type @;
sub boundaries_of {
    my ($self, $locations) = @_;
    if (ref($locations) ne "ARRAY") {
        $locations = [ split /\s*,\s*/, $locations ];
    }
    Trace("Boundaries of [" . join(", ", @{$locations})) if T(4);
    my ($contig,$beg,$end) = $self->{sprout}->GetBoundaries(@{$locations});
    Trace("Boundaries are ($beg, $end) in $contig.") if T(3);
    return ($contig,$beg,$end);
}

=head3 coupling_and_evidence

    my @couplings = $sfxlate->coupling_and_evidence($feature_id);

Return a list of the features functionally coupled to the specified feature
along with their scores. Note that the FIG version of this method has
four additional parameters. If provided, these parameters are simply
ignored.

=over 4

=item peg1

ID of the feature whose couplings and evidence are desired.

=item RETURN

Returns a list of 3-tuples. Each 3-tuple consists of a coupling score
followed by the ID of a coupled feature and a list of the evidence for
the coupling to that feature. The evidence format is the same as that
for L</coupling_evidence>.

=back

=cut
#: Return Type @@;
sub coupling_and_evidence {
    my ($self,$peg1) = @_;
    my %featureHash = $self->{sprout}->CoupledFeatures($peg1);
    my @retVal = ();
    for my $peg2 (keys %featureHash) {
        # Only proceed if this is not the bogus reflexive coupling.
        if ($peg2 ne $peg1) {
            my $sc = $featureHash{$peg2};
            my @ev = map { [$_->[0], $_->[1]] } $self->coupling_evidence($peg1,$peg2);
            push @retVal, [$sc,$peg2,\@ev];
        }
    }
    return @retVal;
}

=head3 coupling_evidence

    my @evidence = $sfxlate->coupling_evidence($peg1, $peg2);

Return the evidence for a functional coupling between two features. A coupling
is considered functional if the specified pegs are frequently found together.
The evidence for the coupling is therefore based on finding other genomes in
which similar pegs are clustered on the same chromosome.

=over 4

=item peg1

ID of the first feature of interest.

=item peg2

ID of the second feature of interest.

=item RETURN

Returns a list of 3-tuples. Each tuple consists of a feature similar to the feature
of interest, a feature similar to the functionally coupled feature, and a flag
that is TRUE for a representative piece of evidence and FALSE otherwise.

=back

=cut
#: Return Type @@;
sub coupling_evidence {
    my($self,$peg1,$peg2) = @_;
    return $self->{sprout}->CouplingEvidence($peg1,$peg2);
}

=head3 is_deleted_fid

    my $flag = $sfxlate->is_deleted_fid($fid);

Return TRUE if the specified feature does B<not> exist, FALSE if it does
exist.

=over 4

=item fid

ID of the feature whose existence is to be tested.

=item RETURN

Returns TRUE if the feature does not exist, else FALSE. Note that if TRUE
is returned, there is no guarantee that the feature ever existed, only that
it does not exist now.

=back

=cut
#: Return Type $;
sub is_deleted_fid {
    my ($self, $fid) = @_;
    my $exists = $self->{sprout}->Exists("Feature", $fid);
    return !$exists;
}

=head3 close_enough

    my $flag = $sfxlate->close_enough($locs1, $locs2, $bound);

Return TRUE if the specified locations are within the specified
distance. The locations must be on the same contig, and the midpoints
must be within the specified bound.

=over 4

=item locs1, locs2

Locations to compare. Each location is a 3-tuple consisting of a contig ID,
a starting offset, and an ending offset. Note that the 3-tuple represents
a SEED-style, not a Sprout-style location.

=item bound

Maximum distance between the midpoints of the location.

=item RETURN

Returns TRUE if the two locations are close enough, else FALSE.

=back

=cut
#: Return Type $;
sub close_enough {
    my ($self, $locs1, $locs2, $bound) = @_;
    return FIG::close_enough($locs1, $locs2, $bound);
}

=head3 taxonomy_of

    my $taxonomy = $sfxlate->taxonomy_of($genome);

Return the taxonomy of the specified genome.

=over 4

=item genome

ID of the genome whose taxonomy is desired.

=item RETURN

Returns a complete taxonomy for the organism, with entiries separated
by semi-colons.

=back

=cut
#: Return Type $;
sub taxonomy_of {
    my ($self, $genome) = @_;
    my @retVal = $self->{sprout}->Taxonomy($genome);
    return join "; ", @retVal;
}

=head3 crude_estimate_of_distance

    my $distance = $sfxlate->crude_estimate_of_distance($genome1, $genome2);

Returns a crude taxonomic distance between the two genomes. The distance
will be 0 for genomes with identical taxonomies and 1 for genomes from
different domains.

=cut
#: Return Type $;
sub crude_estimate_of_distance {
    my ($self, $genome1, $genome2) = @_;
    return $self->{sprout}->CrudeDistance($genome1, $genome2);
}

=head3 ec_name

    my $name = $sfxlate->ec_name($ec);

Return the name of the role identified by the specified EC number.

=over 4

=item ec

EC number of the role whose name is desired.

=item RETURN

Returns the name of the role identified by the specified EC number.

=back

=cut
#: Return Type $;
sub ec_name {
    my ($self, $ec) = @_;
    return $self->{sprout}->RoleName($ec);
}

=head3 coupled_to

    my @coupled_to = $fig->coupled_to($peg);

Return a list of functionally coupled PEGs.

=over 4

=item peg

ID of the protein encoding group whose functionally-coupled proteins are desired.

=item RETURN

Returns a list of 2-tuples, each consisting of the ID of a coupled PEG and a score. If
there are no PEGs functionally coupled to the incoming PEG, it will return an empty
list. If the PEG data is not present, it will return C<undef>.

=back

=cut
#: Return Type @@;
sub coupled_to {
    my ($self, $peg) = @_;
    my %couplets = $self->{sprout}->CoupledFeatures($peg);
    Trace(scalar(keys %couplets) . " couplings returned from Sprout.") if T(coupling => 3);
    my @retVal = ();
    for my $otherPeg (sort keys %couplets) {
        push @retVal, [$otherPeg, $couplets{$otherPeg}];
    }
    return @retVal;
}

# Bruce will have to add abstract coupling data later

sub abstract_coupled_to {
    my ($self, $peg) = @_;
    return undef;
}

=head3 coupled_to_batch

    my @couplings = $fig->coupled_to_batch(@pegs);

Return the functional couplings of one or more features. This method essentially
returns the result one would get from calling L</coupled_to> for each individual
feature, but saves some overhead because it only queries the coupling server
once.

=over 4

=item pegs

A list of the relevant feature IDs.

=item RETURN

Returns a list of 3-tuples. Each tuple consists of a feature from the input list, a coupled-to
feature, and the coupling score.

=back

=cut

sub coupled_to_batch {
    # Get the parameters.
    my ($self, @pegs) = @_;
    # Query the coupling server.
    my @retVal = FIGRules::NetCouplingData('coupled_to_batch', id1 => \@pegs);
    # Return the result.
    return @retVal;
}

=head3 genome_info

    my $info = $fig->genome_info();

Return an array reference of information from the genome table.

=over 4

=item RETURN

Returns a reference to a list of lists, one list per genome. Each genome's list entry
contains the genome ID, organism name, number of base pairs, taxonomic domain,
number of PEGs, number of RNAs, and the complete flag, in that order.

=back

=cut

sub genome_info {
    my ($self) = @_;
    # Get the desired data from the genome table.
    my @gData = $self->{sprout}->GetAll(['Genome'], "", [], ['Genome(id)', 'Genome(taxonomy)', 'Genome(dna-size)',
                                                             'Genome(taxonomy)', 'Genome(pegs)', 'Genome(rnas)',
                                                             'Genome(complete)']);
    # We need to convert the taxonomy information to the domain and species name, so we'll
    # loop through the genomes and reconstruct the list.
    my @retVal = ();
    for my $gDatum (@gData) {
        my @taxa = split /\s*;\s*/, $gDatum->[1];
        # The domain is always first.
        my $domain = $taxa[0];
        # The name is the last three pieces joined together.
        my $gname = join(" ", @taxa[-3, -2, -1]);
        # Put the result together.
        push @retVal, [ $gDatum->[0], $gname, $gDatum->[2], $domain, $gDatum->[4], $gDatum->[5],
                        $gDatum->[6] ]
    }
    # Return the result.
    return \@retVal;
}


=head3 feature_annotations

    my @descriptors = $sfxlate->feature_annotations($feature, $rawFlag);

Return the annotations of a feature.

=over 4

=item feature

ID of the feature whose annotations are desired.

=item rawFlag (optional)

If TRUE, the time will be returned as a raw number; otherwise, the time will
be returned in human-readable form.

=item RETURN

Returns a list of annotation descriptors. Each descriptor is a 4-tuple with
the following elements.

* B<featureID> ID of the relevant feature.

* B<timeStamp> time the annotation was made.

* B<user> ID of the user who made the annotation

* B<text> text of the annotation.

=back

=cut
#: Return Type @@;
sub feature_annotations {
    my ($self, $feature, $rawFlag) = @_;
    my @annotations = $self->{sprout}->FeatureAnnotations($feature, $rawFlag);
    # Sprout hands back hashes. We need to convert them to tuples.
    my @retVal = ();
    for my $tupleHash (@annotations) {
        push @retVal, [$tupleHash->{featureID}, $tupleHash->{timeStamp}, $tupleHash->{user}, $tupleHash->{text}];
    }
    return @retVal;
}

=head3 possibly_truncated

    my $flag = $sfxlate->possibly_truncated($fid);

Returns TRUE if the indicated feature occurs near the end of a contig,
else FALSE. This method calls the FIG.pm method passing in the SFXlate
object as the first parameter. As a result, when the FIG method tries to
call another FIG method, it will call the corresponding method on this
object instead.

=over 4

=item fid

ID of the relevant feature.

=item RETURN

Returns TRUE if the feature's location is near the end of the containing
contig, else FALSE.

=back

=cut
#: Return Type $;
sub possibly_truncated {
    my ($self, $fid) = @_;
    return FIG::possibly_truncated($self, $fid);
}

=head3 near_end

    my $flag = $sfxlate->near_end($genome, $contig, $x);

Return TRUE if the offset I<$x> is near either end of the specified contig,
else FALSE.

=over 4

=item genome

ID of the relevant genome.

=item contig

ID of the relevant contig.

=item x

Offset to check for its proximity to either end of the contig.

=item RETURN

Returns TRUE if the specified offset is within 300 positions of either end
of the contig, else FALSE.

=back

=cut
#: Return Type $;
sub near_end {
    my ($self, $genome, $contig, $x) = @_;
    return FIG::near_end($self, $genome, $contig, $x);
}

=head3 genome_of

    my $genomeID = $sprout->genome_of($fid);

Return the ID of the genome containing the specified feature.

=over 4

=item fid

ID of the feature whose genome is desired.

=item RETURN

Returns the ID of the genome that contains the feature.

=back

=cut
#: Return Type $;
sub genome_of {
    my ($self, $fid) = @_;
    my $retVal;
    # FIG.pm returns undefined if the feature ID is undefined. We need to do the same.
    if (defined $fid) {
        $retVal = $self->{sprout}->GenomeOf($fid);
    }
    return $retVal;
}

=head3 get_attributes

    my @properties = $sfxlate->get_attributes(@parms);

Return a list of n-tuples from the attribute server.


=cut
#: Return Type @@;
sub get_attributes {
    my ($self, @parms) = @_;
    my $ca = $self->attribute_object();
    my @retVal;
    if ($ca) {
        push @retVal, $ca->GetAttributes(@parms);
    }
    return @retVal;
}

# Deprecated feature_attributes call. Use get_attributes instead.
sub feature_attributes {
    return get_attributes(@_);
}

=head3 attribute_object

    my $ca = $fig->attribute_object();

Return the attribute object for this Sprout instance. If an attribute
object is already attached, it will be returned immediately; otherwise,
one will be created.

=cut

sub attribute_object {
    # Get the parameters.
    my ($self) = @_;
    # See if we have an attribute object already cached.
    if (! defined $self->{ca}) {
        # We don't, so determine the type of object we need.
        if ($FIG_Config::attrURL) {
            # Here it's a remote attribute object.
            Trace("Remote attribute server $FIG_Config::attrURL chosen.") if T(3);
            $self->{ca} = RemoteCustomAttributes->new($FIG_Config::attrURL);
        } elsif (! $FIG_Config::attrHost) {
            # Here attributes are disabled.
            $self->{ca} = "";
        } else {
            # Here it's a local database.
            Trace("Local attribute database $FIG_Config::attrDbName chosen.") if T(3);
            my $user = ($FIG_Config::arch eq 'win' ? 'self' : scalar(getpwent()));
            # Insure we recover from errors.
            eval {
                $self->{ca} = CustomAttributes->new(user => $user);
            };
            if ($@) {
                Trace("Attribute connection error: $@") if T(0);
            }
        }
    }
    # Return the cached object.
    return $self->{ca};
}

=head3 get_peg_keys

    my @list = $fig->get_peg_keys();

Return a list of the attribute keys that only apply to pegs. This list is
essentially all of the keys in the C<peg> group.

=cut

sub get_peg_keys {
    # Get the parameters.
    my ($self) = @_;
    # Get the attribute object.
    my $ca = $self->attribute_object();
    # Compute the list of peg keys. Note we may not have a valid attribute
    # object, in which case we return nothing.
    my @retVal;
    if ($ca) {
        @retVal = $ca->GetAttributeKeys('peg');
    }
    # Return the result.
    return @retVal;
}

=head3 get_translation

    my $translation = $sfxlate->get_translation($feature);

Return the protein sequence for the specified feature.

=over 4

=item feature

ID of the feature whose protein sequence is desired.

=item RETURN

Returns the protein sequence for the specified feature.

=back

=cut
#: Return Type $;
sub get_translation {
    my ($self, $feature) = @_;
    return $self->{sprout}->FeatureTranslation($feature);
}

=head3 subsystems_for_peg_complete

    my @list = $fig->subsystems_for_peg_complete($peg);

Return information about the subsystems in which the specified feature
participate.

=over 4

=item peg

ID of the relevant feature or a reference to a list of the relevant features.

=item RETURN

Returns a list of 4-tuples. Each 4-tuple will contain a subsystem name, the role the indicated
feature plays in the subsystem, the variant code, and a flag that is TRUE if the role is auxiliary.

=back

=cut

sub subsystems_for_peg_complete {
    # Get the parameters.
    my ($self, $peg) = @_;
    # Insure that we're dealing with a list of pegs.
    my $pegs = (ref($peg) eq 'ARRAY' ? $peg : [$peg]);
    # Get the sprout object.
    my $sprout = $self->{sprout};
    # Declare the return variable.
    my @retVal;
    # Loop through the PEGs.
    for my $fid (@$pegs) {
        # Get the genome ID.
        my $genomeID = $self->genome_of($fid);
        # Get the desired 4-tuples.
        push @retVal, $sprout->GetAll([qw(ContainsFeature IsRoleOf HasSSCell ParticipatesIn OccursInSubsystem)],
                                      "ContainsFeature(to-link) = ? AND ParticipatesIn(from-link) = ? AND " .
                                      "OccursInSubsystem(from-link) = IsRoleOf(from-link)",
                                      [$fid, $genomeID],
                                      [qw(HasSSCell(from-link) IsRoleOf(from-link) ParticipatesIn(variant-code)
                                          OccursInSubsystem(auxiliary))]);
    }
    # Return the result.
    return @retVal;
}

=head3 subsystems_for_pegs_complete

    my %pegHash = $fig->subsystems_for_pegs_complete(\@pegs, $aux_flag);

Return a hash that maps the incoming pegs to the list of subsystems,
roles and variants that the pegs appear in. Each peg is mapped to an
array of 3-tuples. Each 3-tuple contains a subsystem name, a role name,
and a variant code. Thus, for each peg we will get a list of the
subsystems in which it appears along with the relevant role and variant.

=over 4

=item pegs

Reference to a list of feature IDs.

=item aux_flag

TRUE if auxiliary roles are to be included in the results, else FALSE.

=item RETURN

Returns a hash that maps each incoming feature ID to a list of its
subsystem roles. Each role is represented by a 3-tuple consisting of
the subsystem name, the role name, and the relevant variant code.

=back

=cut

sub subsystems_for_pegs_complete {
    # Get the parameters.
    my ($self, $pegs, $aux_flag) = @_;
    # Get the sprout object.
    my $sprout = $self->{sprout};
    # We need to compute the filter string for this query.
    my $filter = "ContainsFeature(to-link) = ? AND ParticipatesIn(from-link) = ? AND " .
                 "OccursInSubsystem(from-link) = IsRoleOf(from-link)";
    # If we do NOT want auxiliary roles, we add an extra filter condition.
    if (! $aux_flag) {
        $filter .= " AND OccursInSubsystem(auxiliary) = 0";
    }
    # Declare the return variable.
    my %retVal;
    # Loop through the PEGs.
    for my $fid (@$pegs) {
        # Get the genome ID.
        my $genomeID = $self->genome_of($fid);
        # Get the desired 3-tuples.
        my @list = $sprout->GetAll([qw(ContainsFeature IsRoleOf HasSSCell
                                       ParticipatesIn OccursInSubsystem)],
                                   $filter, [$fid, $genomeID],
                                   [qw(HasSSCell(from-link) IsRoleOf(from-link)
                                       ParticipatesIn(variant-code))]);
        # Store the list in the return hash.
        $retVal{$fid} = \@list;
    }

    # Return the result.
    return %retVal;
}

=head3 subsystems_for_pegs

    my @list = $fig->subsystems_for_pegs(\@pegs, $noaux);

Return information about the subsystems in which the features in a list
participate. For each incoming feature, this method will return a list of
2-tuples, the first element being the subsystem name and the second being
the feature's role in that subsystem.

=over 4

=item pegs

Reference to a list of features whose subsystem information is desired.

=item noaux (optional)

If TRUE, subsystems in which a feature has an auxiliary role will be omitted
from the results.

=item RETURN

Returns a list that contains feature IDs followed by list references. For example,
if three features were specified in the parameters, the list would have six
elements. The first element will be a feature ID, the second will be a reference
to a list of 2-tuples, the third will be another feature ID, the fourth will be
a reference to another list of 2-tuples, and so on. Each 2-tuple contains the
name of a subsystem followed by a role. The 2-tuples contain the subsystem information
for the preceding feature.

=back

=cut

sub subsystems_for_pegs {
    # Get the parameters.
    my ($self, $pegs, $noaux) = @_;
    # If we're filtering out auxiliary roles, create a filter clause for that.
    my $auxFilter = ($noaux ? "AND OccursInSubsystem(auxiliary) = 0" : "");
    # Build a hash of the result data.
    my %retVal = ();
    # Loop through the pegs.
    for my $peg (@{$pegs}) {
        # Get this peg's subsystem data. It will come back as a series of 2-tuples.
        my @subPairs = $self->{sprout}->GetAll([qw(ContainsFeature IsRoleOf HasSSCell OccursInSubsystem)],
                                               "ContainsFeature(to-link) = ? AND " .
                                               "OccursInSubsystem(from-link) = IsRoleOf(from-link) $auxFilter", [$peg],
                                               [qw(HasSSCell(from-link) IsRoleOf(from-link))]);
        # Add it to the hash we're building.
        $retVal{$peg} = \@subPairs;
    }
    # Return the result. We convert the hash to a list.
    return (%retVal);
}


=head3 families_containing_peg

    my @fams = $fig->families_containing_peg($fid);

Return a list of the names of the FIGfams containing the specified
feature.

=over 4

=item fid

ID of the feature of interest.

=item RETURN

Returns a list of FIGfam IDs for families containing the feature. If the feature
is not in any FIGfam, it returns an empty list. Currently, no feature can be in
more than one FIGfam, but this is not necessarily guaranteed for the future.

=back

=cut

sub families_containing_peg {
    # Get the parameters.
    my ($self, $fid) = @_;
    # Get a lightweight FigFams object.
    my $figfam_data = &FIG::get_figfams_data();
    my $ffs = new FFs($figfam_data, $self);
    my @retVal = $ffs->families_containing_peg($fid);
    # Return the result.
    return @retVal;
}


=head3 merged_related_annotations

    my @annotations = $sfxlate->merged_related_annotations(\@list);

Returns a merged list of the annotations for the features in a list. Each annotation is
represented by a 4-tuple of the form C<($fid, $timestamp, $userID, $annotation)>, where
C<$fid> is the ID of a feature, C<$timestamp> is the time at which the annotation was made,
C<$userID> is the ID of the user who made the annotation, and C<$annotation> is the annotation
text. The list is sorted by the timestamp.

=over 4

=item list

List of the IDs for the features whose annotations are desired.

=item RETURN

Returns a list of annotation descriptions sorted by the annotation time.

=back

=cut
#: Return Type @@;
sub merged_related_annotations {
    my ($self, $list) = @_;
    my @retVal = $self->{sprout}->MergedAnnotations($list);
    return @retVal;
}

=head3 run

    SFXlate::run($cmd);

or

    $fig->run($cmd);

=head3 run

    FIG::run($cmd);

or

    $fig->run($cmd);

Run a command. If the command fails, the error will be traced.

=over 4

=item cmd

Text of the command to run. The C<system> function will be used to
invoke the command.

=back

=cut

sub run {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cmd) = @_;

    if ($ENV{FIG_VERBOSE}) {
        my @tmp = `date`;
        chomp @tmp;
        print STDERR "$tmp[0]: running $cmd\n";
    }
    Trace("Running command: $cmd") if T(3);
    (system($cmd) == 0) || Confess("FAILED: $cmd");
}

=head3 run_gathering_output

    my @lines = SFXlate::run_gathering_output($cmd, @args);

or

    my @lines = $fig->run_gathering_output($cmd, @args);

or

    my $text = $fig->run_gathering_output($cmd, @args);

Run a command, gathering the output. This is similar to the backtick
operator, but it does not invoke the shell. Note that the argument list
must be explicitly passed in the parameter lit.

If the command fails, the error will be traced.

=over 4

=item cmd

Name of the command to run.

=item args

List of the arguments to the command. Each argument must be
passed as a separate element of this list.

=item RETURN

In array context, returns a list of output lines. In a scalar context,
returns the command output as a single string.

=back

=cut

sub run_gathering_output {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    my($cmd, @args) = @_;

    #
    # Run the command in a safe fork-with-pipe/exec.
    #

    my $pid = open(PROC_READ, "-|");

    if ($pid == 0) {
        exec { $cmd } $cmd, @args;
    }
    if ($pid == 0) {
        # This next statement will only execute if the exec function fails. The goofy
        # use of the redundant IF is to avoid a compiler warning.
        Confess("Could not execute $cmd @args: $!");
    }
    if (wantarray) {
        my @out;
        while (<PROC_READ>) {
            push(@out, $_);
        }
        if (!close(PROC_READ)) {
            Confess("FAILED: $cmd @args with error return $?");
        }
        return @out;
    } else {
        my $out = '';

        while (<PROC_READ>) {
            $out .= $_;
        }
        if (!close(PROC_READ)) {
            Confess("FAILED: $cmd @args with error return $?");
        }
        return $out;
    }
}

=head3 assign_function

    my $ok = $sfxlate->assign_function($fid,$user,$function);

Assign a function to the specified feature. In Sprout, an assigned function
is a special type of annotation. The Sprout methods automatically convert
the function text into the correct structured formay.

=over 4

=item fid

ID of the feature to receive the assignment.

=item user

Name of the user making the assignment.

=item function

Text of the functional assignment.

=back

=cut
#: Return Type $;
sub assign_function {
    my ($self, $fid, $user, $function) = @_;
    my $retVal = $self->{sprout}->AssignFunction($fid, $user, $function);
    return $retVal;
}

=head3 neighborhood_of_role

    my @roleList = $sprout->neighborhood_of_role($role);

Returns a list of the roles that occur in the same diagram as the specified role. Because
diagrams and roles are in a many-to-many relationship with each other, the list is
essentially the set of roles from all of the maps that contain the incoming role. Such
roles are considered neighbors because they are used together in cellular subsystems.

=over 4

=item role

ID of the role whose neighbors are desired.

=item RETURN

Returns a list containing the IDs of the roles that are related to the incoming role.

=back

=cut
#: Return Type @;
sub neighborhood_of_role {
    my ($self, $role) = @_;
    my @retVal = $self->{sprout}->RoleNeighbors($role);
    Trace("roles found = " . join(", ", @retVal)) if T(4);
    return @retVal;
}

=head3 org_and_color_of

    my ($orgName, $color) = $sfxlate->org_and_color_of($fid);

Return the name of the organism to which the specified feature belongs.
The organism name is generally the genus and species followed by the
unique characterization.

=over 4

=item fid

ID of the feature whose organism information is desired.

=item RETURN

Returns a 2-tuple. The first element is a string consisting of the genus,
species, and unique characterization of the specified feature's organism.
The second element is an HTML color code based on the domain.

=back

=cut
#: Return Type @;
sub org_and_color_of {
    my ($self, $fid) = @_;
    my $genome = $self->{sprout}->GenomeOf($fid);
    my @taxonomy = $self->{sprout}->Taxonomy($genome);
    my $color = FIG::domain_color($taxonomy[0]);
    my $gs = $taxonomy[$#taxonomy];
    return ($gs, $color);
}

=head3 org_of

    my $orgName = $sfxlate->org_of($fid);

Return the name of the organism to which the specified feature belongs.
The organism name is generally the genus and species followed by the
unique characterization.

=over 4

=item fid

ID of the feature whose organism information is desired.

=item RETURN

Returns the genus, species, and unique characterization of the specified
feature's organism. The result is returned as a single string.

=back

=cut
#: Return Type $;
sub org_of {
    my ($self, $fid) = @_;
    my $genome = $self->{sprout}->GenomeOf($fid);
    my $retVal = $self->{sprout}->GenusSpecies($genome);
    return $retVal;
}

=head3 to_structured_english

    my ($ev_code_list, $subsys_list, $english_string) = $fig->to_structured_english($fig, $peg, $escape_flag);

Create a structured English description of the evidence codes for a PEG,
in either HTML or text format. In addition to the structured text, we
also return the subsystems and evidence codes for the PEG in list form.

=over 4

=item peg

ID of the protein or feature whose evidence is desired.

=item escape_flag

TRUE if the output text should be HTML, else FALSE

=item RETURN

Returns a three-element list. The first element is a reference to a list of evidence codes,
the second is a list of the subsystem containing the peg, and the third is the readable
text description of the evidence.

=back

=cut

sub to_structured_english {
    my ($self, $peg, $escape_flag) = @_;
    return FIGRules::to_structured_english($self, $peg, $escape_flag);
}


=head3 peg_to_subsystems

    my @subsystems = $sfxlate->peg_to_subsystems($fid);

Return a list of the subsystems in which the specified feature participates.
In the Sprout system, a subsystem is connected to features indirectly
via the B<SSCell> object.

=over 4

=item fid

ID of the feature whose subsystems are desired.

=item RETURN

Returns a list of the IDs of the subsystems containing the specified
feature.

=back

=cut
#: Return Type @;
sub peg_to_subsystems {
    my ($self, $fid) = @_;
    my @retVal = $self->{sprout}->SubsystemList($fid);
    return @retVal;
}

=head3 hypo

    my $flag = $sfxlate->hypo($func);

Return TRUE if the specified functional role is hypothetical, else FALSE.
Hypothetical functional roles are identified by key words in the text,
such as I<hypothesis>, I<predicted>, or I<glimmer> (among others).

=over 4

=item func

Text of the functional role whose nature is to be determined.

=item RETURN

Returns TRUE if the role is hypothetical, else FALSE.

=back

=cut
#: Return Type $;
sub hypo {
    my ($self, $func) = @_;
    return FIG::hypo($func);
}

=head3 related_by_func_sim

    my @fids = $sfxlate->related_by_func_sim($fid, $user);

Return a list of the features that have a similar function to the specified
feature as determined by the specified user. This method looks at the
bidirectional best hits of the incoming feature and returns a list of the
ones who have the same functional assignment with respect to the specified
user.

=over 4

=item fid

ID of the feature whose similarities are desired.

=item user

ID of the user whose functional assignments are to be used. The functional
assignments chosen will be the most recent ones by the specified user or
a user trusted by the specified user. If no user is specified, only the
user C<FIG> will be considered.

=item RETURN

Returns a list of the IDs of the desired features.

=back

=cut
#: Return Type @;
sub related_by_func_sim {
    my ($self, $fid, $user) = @_;
    my $function = $self->{sprout}->FunctionOf($fid, $user);
    my @retVal = ();
    if (! $self->hypo($function)) {
        push @retVal, $self->{sprout}->RelatedFeatures($fid, $function, $user);
    }
    return @retVal;
}

=head3 sort_fids_by_taxonomy

    my @sortedFids = $sfxlate->sort_fids_by_taxonomy(@fidList);

Sort the specified list of features according to the taxonomy of the
feature's genome. The intent is to group features belonging to similar
organisms.

=over 4

=item fidList

List of feature IDs.

=item RETURN

Returns a list containing the same feature IDs ordered by their respective
taxonomies.

=back

=cut
#: Return Type @;
sub sort_fids_by_taxonomy {
    my ($self, @fidList) = @_;
    return $self->{sprout}->TaxonomySort(\@fidList);
}

=head3 translatable

    my $flag = $sfxlate->translatable($fid);

Return TRUE if the specified feature has a translation, else FALSE.

=over 4

=item fid

ID of the feature whose translatability is to be determined.

=item RETURN

Returns TRUE if the feature exists in the database, else FALSE.

=back

=cut
#: Return Type $;
sub translatable {
    my ($self, $fid) = @_;
    my $feature = $self->{sprout}->GetEntity('Feature', $fid);
    my $retVal = 0;
    if ($feature) {
        $retVal = 1;
    }
    return $retVal;
}

=head3 peg_links

    my @linkList = $sfxlate->peg_links($fid);

List the links associated with a feature. These are generally HTML
hyperlinks to pages with information about the feature or the feature's
associated protein.

=over 4

=item fid

ID of the feature whose links are desired.

=item RETURN

Returns a list of the HTML links stored for the specified feature, sorted
more or less by the link's target URL.

=back

=cut
#: Return Type @;
sub peg_links {
    my ($self, $fid) = @_;
    my @links = $self->{sprout}->FeatureLinks($fid);
    return sort { $a =~ /\>([^\<]+)\<\/a\>/; my $l1 = $1;
          $b =~ /\>([^\<]+)\<\/a\>/; my $l2 = $1;
          $l1 cmp $l2 } @links;
}

=head3 get_gbrowse_feature_link

    my $url = $sfxlate->get_gbrowse_feature_link;

Compute the URL required to pull up a Gbrowse page for the the
specified feature. In order to do this, we need to pull out
the ID of the feature's Genome, its contig ID, and some rough
starting and stopping offsets.

=over 4

=item feat

ID of the feature whose Gbrowse URL is desired.

=item RETURN

Returns a GET-style URL for the Gbrowse CGI, with parameters
specifying the genome ID, contig ID, starting offset, and
stopping offset.

=back

=cut
#: Return Type $;
sub get_gbrowse_feature_link {
    my($self, $feat) = @_;

    my $genome;

    if ($feat =~ /fig\|(\d+\.\d+)/) {
        $genome = $1;
    } else {
       return undef;
    }

    my $gs = $self->{sprout}->GenusSpecies($genome);
    my $loc = $self->{sprout}->FeatureLocation($feat);

    my($start, $stop, $contig);

    #
    # Eval this code to catch possible badness in the database loads.
    #

    eval {

        $start = $self->beg_of($loc);
        $stop = $self->end_of($loc);
        $contig = $self->contig_of($loc);
    };

    if ($@)
    {
        warn "Error in getting location information for feature $feat\n$@\n";
        return undef;
    }

    my $mid = int(($start + $stop) / 2);

    my $chunk_len = 20000;
    my $max_feature = 40000;
    #
    # Make sure large features show up.
    #
    # However, if the feature is larger than max_feature,
    # show the start.
    #
    my $feat_len = abs($stop - $start);

    if ($feat_len > $chunk_len) {
        if ($feat_len > $max_feature) {
            $chunk_len = $max_feature;
        } else {
            $chunk_len = $feat_len + 100;
        }
    }

    my($show_start, $show_stop);
    if ($chunk_len == $max_feature) {
        $show_start = $start - 300;
    } else {
        $show_start = $mid - int($chunk_len / 2);
    }
    if ($show_start < 1) {
        $show_start = 1;
    }
    $show_stop = $show_start + $chunk_len - 1;

    my $clen = $self->{sprout}->ContigLength($contig);
    if ($show_stop > $clen) {
        $show_stop = $clen;
    }

    my $seg_id = $contig;
    $seg_id =~ s/:/--/g;
    return ("/gbrowse.cgi/GB_$genome?ref=$seg_id&start=$show_start&stop=$show_stop");
}


=head3 subsystems_for_peg

    my @ssList = $sfxlate->subsystems_for_peg($fid);

Return a list of the subsystems in which a specified feature participates.unlike L</peg_to_subsystems>,
this method returns the role the feature plays in the subsystem in addition to the subsystem name.

=over 4

=item fid

ID of the feature whose subsystem list is desired.

=item RETURN

Returns a list of 2-tuples, one per subsystem. Each tuple consists of the
subsystem ID followed by a role the input feature plays in that subsystem.

=back

=cut
#: Return Type @@;
sub subsystems_for_peg {
    # Get the parameters.
    my ($self, $featureID) = @_;
    my $sprout = $self->{sprout};
    # Get the subsystem list.
    my @subsystems = $sprout->GetAll([qw(ContainsFeature HasSSCell IsRoleOf)],
                                    "ContainsFeature(to-link) = ?", [$featureID],
                                    [qw(HasSSCell(from-link) IsRoleOf(from-link))]);
    # Create the return value.
    return @subsystems;
}

=head3 bbhs

    my @bbhList = $sfxlate->bbhs($peg, $cutoff);

Return a list of the specified feature's bidirectional best hits. All the
hits returned will have a score lower than the specified cutoff score.

=over 4

=item peg

ID of the feature whose BBHs are desired.

=item cutoff

Maximum permissible score to be accepted. If omitted, 1e-10 is used.

=item RETURN

Returns a list of 2-tuples. Each tuple will consist of a feature ID followed
by a score. The identified feature will be a bidirectional best hit of the
incoming feature and the score is guaranteed to be no greater than the
cutoff value.

=back

=cut
#: Return Type @@;
sub bbhs {
    my ($self,$peg,$cutoff) = @_;
    if (! defined $cutoff) { $cutoff = 1e-10; }
    my %bbhMap = $self->{sprout}->LowBBHs($peg, $cutoff);
    my @retVal = ();
    for my $peg2 (keys %bbhMap) {
        Trace("Pushing BBH from $peg to $peg2.") if T(4);
        push @retVal, [$peg2, $bbhMap{$peg2}];
    }
    return @retVal;
}

=head3 get_pathogen_groups

    my @groups = $sfxlate->get_pathogen_groups();

Return a list of all the pathogen groups and the genomes in them.

This method is a special case of the Sprout B<GetGroups> function. It
returns a list of 2-tuples, with each tuple containing the name of a
group followed by a list of the genomes in the group, each genome
being represented in the list by its ID. Note that the genome IDs are
a full three layers deep. The code to get the first genome ID for the
third group is C<< $groups[2]->[1]->[0] >>. The C<2> indicates we are
asking for the third group, the C<1> means we're asking for a group member
rather than the group name, and the C<0> means we want the first
group member.  Note that a single genome could be in more than one group,
and that in fact many genomes are not in a group at all.

=cut
#: Return Type @@;
sub get_pathogen_groups {
    my ($self) = @_;
    # Sprout returns a hash of lists. Because there are no parameters,
    # all groups will be returned.
    my %groups = $self->{sprout}->GetGroups();
    # Now we convert the hash to a list of lists.
    my @retVal = ();
    for my $groupName (sort keys %groups) {
        my $groupList = [$groupName, $groups{$groupName}];
        push @retVal, $groupList;
    }
    return @retVal;
}

=head3 subsystem_classification

    my $class = $sfx->subsystem_classification($subsystemName);

Return the classification of a subsystem.

=over 4

=item subsystemName

Name of the subsystem whose classification is desired.

=item RETURN

Returns a reference to a list of the classification elements of the subsystem, or an
empty list if the subsystem is not classified.

=back

=cut

sub subsystem_classification {
    # Get the parameters.
    my ($self, $subsystemName) = @_;
    # Declare the return variable.
    my $retVal;
    # Try to get the subsystem classifications.
    my ($classes) = $self->{sprout}->GetFlat(['Subsystem'], "Subsystem(id) = ?", [$subsystemName],
                                              'Subsystem(classification)');
    if (defined $classes) {
        $retVal = [split($FIG_Config::splitter, $classes)];
    } else {
        $retVal = [];
    }
    # Return the result.
    return $retVal;
}


=head3 usable_subsystem

    my $flag = $fig->usable_subsystem($sub);

Return TRUE if the named subsystem is not experimental or deleted.

=over 4

=item sub

Name of the subsystem to check.

=item RETURN

Returns TRUE if the named subsystem is usable, FALSE if it is experimental or deleted.

=back

=cut

sub usable_subsystem {
    # Get the parameters.
    my($self, $sub) = @_;
    # Declare the return value. We default to not usable.
    my $retVal = 0;
    # Get the subsystem record.
    my $subsys = $self->{sprout}->GetEntity(Subsystem => $sub);
    # Only proceed if we found it.
    if (! defined $subsys) {
        Trace("Subsystem $subsys not found.") if T(3);
    } else {
        # Get the subsystem's classifications.
        my @cats = $subsys->Values(['Subsystem(classification)']);
        # Look for an experimental or deleted marker.
        $retVal = 1;
        for my $cat (@cats) {
            if ($cat =~ /experimental/i || $cat =~ /delete/i) {
                Trace("Subsystem not usable: classification is $cat.") if T(3);
                $retVal = 0;
            }
        }
    }
    return $retVal;
}


=head3 is_genome

    my $flag = $fig->is_genome($genome);

Return TRUE if the specified genome exists, else FALSE.

=over 4

=item genome

ID of the genome to test.

=item RETURN

Returns TRUE if a genome with the specified ID exists in the data store, else FALSE.

=back

=cut

sub is_genome {
    # Get the parameters.
    my ($self, $genome) = @_;
    my $retVal;
    # Only proceed if the genome ID exists. The FIG.pm method allows undefined as a parameter.
    if (defined($genome)) {
        # Test for the genome ID.
        $retVal = $self->{sprout}->Exists('Genome', $genome);
    }
    # Return the test result.
    return $retVal;
}

=head3 function_of_bulk

    my $functionHash = $fig->function_of_bulk(\@fids, $no_del_check);

Return a hash mapping the specified proteins to their master functional assignments.

=over 4

=item fids

Reference to a list of feature IDs.

=item no_del_check

If TRUE, then deleted features B<will not> be removed from the list. The default
is FALSE, which means deleted feature B<will> be removed from the list.

=item RETURN

REturns a reference to a hash mapping feature IDs to their main functional assignments.

=back

=cut

sub function_of_bulk {
    # Get the parameters.
    my ($self, $fids, $no_del_check) = @_;
    # Remove any deleted features from the list according to the value of the
    # no_del_check parameter. Note we copy the list in the process so we don't
    # do any damage to the caller's data.
    my $del_check = ($no_del_check ? 0 : 1);
    my @fids = grep { $del_check || ! $self->is_deleted_fid($_) } @{$fids};
    # Declare the return variable.
    my $retVal = {};
    # Get the underlying Sprout object.
    my $sprout = $self->{sprout};
    # Loop through the features.
    for my $fid (@fids) {
        # Look for a functional assignment.
        my $assignment = $sprout->FunctionOf($fid);
        # If we found one, remember it.
        if ($assignment) {
            $retVal->{$fid} = $assignment;
        }
    }
    # Return the hash.
    return $retVal;
}

=head3 uniprot_aliases_bulk

    my $hash = $fig->uniprot_aliases_bulk(\@fids, $no_del_check);

Return a hash mapping the specified feature IDs to lists of their uniprot
aliases.

=over 4

=item fids

A list of FIG feature IDs.

=item no_del_check

If TRUE, deleted feature IDs B<will not> be removed from the feature ID list
before processing. The default is FALSE, which means deleted feature IDs
B<will> be removed before processing.

=item RETURN

Returns a hash mapping each feature ID to a list of its uniprot aliases.

=back

=cut

sub uniprot_aliases_bulk {
    # Get the parameters.
    my ($self, $fids, $no_del_check) = @_;
    # Remove any deleted features from the list according to the value of the
    # no_del_check parameter. Note we copy the list in the process so we don't
    # do any damage to the caller's data.
    my $del_check = ($no_del_check ? 0 : 1);
    my @fids = grep { $del_check || ! $self->is_deleted_fid($_) } @{$fids};
    # Declare the return variable.
    my $retVal = {};
    # Get the underlying Sprout object.
    my $sprout = $self->{sprout};
    # Loop through the features.
    for my $fid (@fids) {
        # Get this feature's aliases.
        my @aliases = $sprout->GetFlat(['IsAliasOf'], "IsAliasOf(to-link) = ?",
                                      [$fid], 'IsAliasOf(from-link)');
        # Put the uniprot aliases into the hash for the given feature.
        my @unis =  grep { $_ =~ /^uni\|/ } @aliases;
        if (@unis) {
            $retVal->{$fid} = [ sort @unis ];
        }
    }
    # Return the result.
    return $retVal;
}

=head3 by_alias

    my @features = $sfxlate->by_alias($alias);

or

    my $features = $sfxlate->by_alias($alias);

Returns a list of features with the specified alias. The alias is parsed to determine
the type of the alias. A string of digits is a GenBack ID and a string of exactly 6
alphanumerics is a UniProt ID. A built-in FIG.pm method is used to analyze the alias
string and attach the necessary prefix. If the result is a FIG ID then it is returned
unmodified; otherwise, we look for an alias.

=over 4

=item alias

Alias whose feature is desired.

=item RETURN

Returns the ID of the feature with the given alias. In a list context, the feature
ID is returned as a singleton list; in a scalar context, it's returned as a string.

=back

=cut
#: Return Type $;
#: Return Type @;
sub by_alias {
    my ($self, $alias) = @_;
    my @retVal = $self->{sprout}->FeaturesByAlias($alias);
    if (@retVal == 0) {
        return (wantarray ? () : "");
    } else {
        return (wantarray ? @retVal : $retVal[0]);
    }
}

=head3 genome_domain

    my $domain = $sfxlate->genome_domain($genomeID);

Return the domain for a specified genome: Archaea, Bacteria, Eukaryotes, Viruses, or
Environmental Samples.

=cut
#: Return Type $;
sub genome_domain {
    my ($self, $genome) = @_;
    my ($retVal) = $self->{sprout}->Taxonomy($genome);
    return $retVal;
}

=head3 find_role_in_org

    my @table = $fig->find_role_in_org($role, $org, $user, $sims_cutoff);

Find features in a specified organism that probably have the specified
functional role. To do this, we look for features with the specified role
in close genomes, then find features in this genome that are similar. This
will only work if the role in question can be found in a subsystem.

=over 4

=item role

Text of the desired functional role.

=item org

Genome ID for the target organism.

=item user

Name of the user whose annotations are of interest. This parameter is required in
[[FigPm]], but it is ignored here.

=item sims_cutoff

Cutoff value to use for similarities.

=item RETURN

Returns a list of 7-tuples, each containing a p-score, the ID of a feature in
the target organism, its amino acid count, its current functional role,
the ID of the similar feature, its amino acid count,
and its functional role. The desired feature IDs are in the second position
of each tuple; the remaining values are designed to make it easier to interpret
the results.

=back

=cut

sub find_role_in_org {
    # Get the parameters.
    my ($self, $role, $org, $user, $sims_cutoff) = @_;
    # Get the database.
    my $sprout = $self->{sprout};
    # Declare the return variable.
    my @retVal;
    # Find all features with the specified role.
    my @candidatesData = $sprout->GetAll("IsRoleOf ContainsFeature IsInGenome Genome",
                                         "IsRoleOf(from-link) LIKE ?",
                                         [$role],
                                         "IsInGenome(from-link) Genome(taxonomy)");
    Trace(scalar(@candidatesData) . " candidates found for role \"$role\".") if T(3);
    # Get the target genome's taxonomy.
    my @orgTaxonomy = $sprout->Taxonomy($org);
    # For each candidate, compute the taxonomic distance of its genome to ours.
    my %candidates;
    for my $candidate (@candidatesData) {
        # Get this candidate's data.
        my ($peg, $taxonomy) = @$candidate;
        # Split its taxonomy.
        my @taxonomy = split /\s*;\s*/, $taxonomy;
        # Compute the taxonomic distance.
        $candidates{$peg} = FIGRules::CrudeDistanceFormula(\@taxonomy, \@orgTaxonomy);
    }
    # Sort the pegs by distance.
    my @pegs = Tracer::SortByValue(\%candidates);
    # We process at most ten hits.
    if (scalar(@pegs) > 10) {
        splice @pegs, 10, scalar(@pegs) - 10;
    }
    Trace("Processing: " . join(", ", @pegs)) if T(3);
    # Find similarities for these hits. We retrieve a limited number.
    my @possibleSims = $self->sims(\@pegs, $FIG_Config::estimation_sim_limit, $sims_cutoff);
    # Only retain those that are in the target genome. Here we also weed out non-FIG
    # features.
    my @targetSims = grep { $_->id2 =~ /^fig\|$org/ } @possibleSims;
    Trace(scalar(@targetSims) . " target sims found out of " . scalar(@possibleSims) . " possibles.") if T(3);
    # Loop through them, building the output list.
    for my $sim (@targetSims) {
        # Get the two features relevant to this similarity.
        my $ourFeature = $sprout->GetEntity(Feature => $sim->id2);
        # Only proceed if the feature exists.
        if (defined $ourFeature) {
            # This one has to exist because we found it earlier when we were isolating
            # candidates.
            my $otherFeature = $sprout->GetEntity(Feature => $sim->id1);
            # Get both translations.
            my ($ourTran) = $ourFeature->Value('translation');
            my ($otherTran) = $ourFeature->Value('translation');
            # Build the result tuple.
            my $tuple = [$sim->psc,
                         $sim->id2, length $ourTran,
                         $ourFeature->PrimaryValue('assignment'),
                         $sim->id1, length $otherTran,
                         $otherFeature->PrimaryValue('assignment')];
            push @retVal, $tuple;
        }
    }
    # Return the result.
    return @retVal;
}

=head3 is_complete

    my $flag = $fig->is_complete($genome);

Return TRUE if the genome with the specified ID is complete, else FALSE.

=over 4

=item genome

ID of the relevant genome.

=item RETURN

Returns TRUE if there is a complete genome in the database with the specified ID,
else FALSE.

=back

=cut
#: Return Type $;
sub is_complete {
    my ($self, $genome) = @_;
    # return $self->FIG()->is_complete($genome);
    return $self->{sprout}->IsComplete($genome);
}

=head3 all_features_detailed_fast

    my $featureList = $fig->all_features_detailed($genomeID, $min, $max, $contig);

Returns a list of all features in the designated genome, with various useful information
included.

Deleted features are not returned!

=over 4

=item genome

ID of the genome whose features are desired.

=item min (optional)

If specified, the minimum contig location of interest. Features not entirely to the right
of this location are ignored.

=item max (optional)

If specified, the maximum contig location of interest. Features not entirely to the left
of this location are ignore.

=item contig (optional)

If specified, the contig of interest. Features not on this contig are ignored.

=item RETURN

Returns a reference to a list of tuples. Each tuple consists of nine elements: (0) the feature
ID, (1) the feature location (as a comma-delimited list of location specifiers), (2) the feature
aliases (as a comma-delimited list of named aliases), (3) the feature type, (4) the leftmost
index of the feature's leftmost location, (5) the rightmost index of the feature's rightmost location,
(6) the current functional assignment, (7) the user who made the assignment, and (8) the
quality of the assignment (which is usually a space).

=back

=cut

sub all_features_detailed_fast {
    # Get the parameters.
    my ($self, $genome, $min, $max, $contig) = @_;
    # We are going to service this request by reading feature and location records. Our biggest
    # performance hit is going to be the aliases, which have to be read one feature at a time.
    # For the rest, however, we will be getting one value per feature/location pair. Most of
    # the features will be single-location, but we have to beware of the multi-location ones
    # every step of the way. They will require a bit of read-ahead logic.
    # First, however, we create the filter clause. The genome filtering is easy, but the
    # location filtering is complicated by the fact that it's optional.
    my @possibleFilters = ("HasFeature(from-link) = ?", "IsLocatedIn(beg) > ?",
                           "IsLocatedIn(beg) + IsLocatedIn(len) < ?",
                           "IsLocatedIn(to-link) = ?");
    my @possibleParms =   ($genome, $min, $max, $contig);
    my @actualFilters = ();
    my @actualParms = ();
    Trace("All_features_detailed_fast called for $genome.") if T(3);
    for (my $i = 0; $i <= $#possibleFilters; $i++) {
        Trace("Checking detailed_fast parameter $i.") if T(3);
        if (defined $possibleParms[$i]) {
            Trace("detailed_fast parameter $i is $possibleParms[$i].") if T(3);
            push @actualFilters, $possibleFilters[$i];
            push @actualParms, $possibleParms[$i];
        }
    }
    my $actualFilter = join(" AND ", @actualFilters);
    my $actualParmList = \@actualParms;
    # Now we use the filter to make a query.
    my $query = $self->{sprout}->Get(['HasFeature', 'Feature', 'IsLocatedIn'],
                                     "$actualFilter ORDER BY Feature(id)",
                                     $actualParmList);
    # With the query in hand, we declare the return variable.
    my @retVal = ();
    # The following variable will contain the ID of the feature currently being processed. That
    # feature ID is presumed to be the one in the last entry of @retVal.
    my $fid = "";
    # We will fill the return variable using the following loop, which runs through all the query results.
    # The only tricky part is that we may have multiple results for a single feature. We process the first
    # one and ignore the rest.
    while (my $object = $query->Fetch()) {
        # Find out if this is a new feature or a new location for an old feature.
        my $newFid = $object->PrimaryValue('Feature(id)');
        if ($newFid ne $fid) {
            Trace("Processing feature $newFid.") if T(4);
            # Here we have a new feature, so we need to create a new row. We start with the easy stuff.
            my ($type, $assignment, $user, $quality, $locs) = $object->Values(['Feature(feature-type)',
                                                                               'Feature(assignment)',
                                                                               'Feature(assignment-maker)',
                                                                               'Feature(assignment-quality)',
                                                                               'Feature(location-string)']);
            # Next, we get the aliases and convert them to a string.
            my @aliases = $self->{sprout}->GetFlat(['IsAliasOf'], 'IsAliasOf(to-link) = ?', [$newFid], 'IsAliasOf(from-link)');
            my $aliasList = join(",", @aliases);
            # We also need a location object for the locations.
            my $locObject = FullLocation->new($self, $genome, $locs);
            # Get the boundaries.
            my ($boundingLoc) = $locObject->GetBounds();
            # Now we can build the row.
            my @newRow = ($newFid, $locObject->SeedString, $aliasList, $type, $boundingLoc->Left,
                          $boundingLoc->Right, $assignment, $user, $quality);
            # Push it into the return list.
            push @retVal, \@newRow;
            # Remember the feature ID.
            $fid = $newFid;
        }
    }
    # Return the result as a list reference.
    return \@retVal;
}

=head3 all_features

    my @featureIDs = $sfx->all_features($genomeID,$type);

Return a list of the IDs of all the features for a specified genome.

=over 4

=item genomeID

ID of the genome whose features are desired.

=item type (optional)

Type of feature desired (peg, rna, etc.). If omitted, all features will be returned.

=item RETURN

=back

=cut
#: Return Type @;
sub all_features {
    # Get the parameters.
    my ($self, $genomeID, $type) = @_;
    # Form the filter clause.
    my $filter = "HasFeature(from-link) = ?";
    my @parms = ($genomeID);
    if ($type) {
        $filter .= " AND HasFeature(type) = ?";
        push @parms, $type;
    }
    # Ask for the feature IDs.
    my @retVal = $self->{sprout}->GetFlat(['HasFeature'], $filter, \@parms, 'HasFeature(to-link)');
    # Return the result.
    return @retVal;
}

=head3 is_real_feature

    my $flag = $sfxlate->is_real_feature($fid);

Return TRUE if the specified feature is in the database, else FALSE.

=over 4

=item fid

ID of the feature whose existence is in question.

=item RETURN

Returns TRUE if the specified feature exists in the database, else FALSE.

=back

=cut
#: Return Type $;
sub is_real_feature {
    my ($self, $fid) = @_;
    return $self->{sprout}->Exists('Feature', $fid);
}

=head3 active_subsystems

    my $ssHash = $fig->active_subsystems($genome, $allFlag);

Get all the subsystems in which a genome is present. The return value is a hash
which maps each subsystem name to the code for the variant used by the specified
genome.

=over 4

=item genome

ID of the genome whose subsystems are desired.

=item allFlag (optional)

If TRUE, all subsystems are returned, with unknown variants marked by a variant
code of C<-1> and iffy variants marked by a code of C<0>. If FALSE or omitted,
only subsystems in which the variant is definitively known are returned. The
default is FALSE.

=back

=cut

sub active_subsystems {
    # Get the parameters.
    my($self, $genome, $allFlag) = @_;
    # Build the filter from the all-flag.
    my $filter = "ParticipatesIn(from-link) = ?";
    if (! $allFlag) {
        $filter .= " AND ParticipatesIn(variant-code) >= 0";
    }
    # Get the subsystems and codes.
    my @subList = $self->{sprout}->GetAll(['ParticipatesIn'], $filter, [$genome],
                                        ['ParticipatesIn(to-link)', 'ParticipatesIn(variant-code)']);
    # Form them into a hash.
    my %retVal = map { $_->[0] => $_->[1] } @subList;
    # Return the results.
    return \%retVal;
}

=head3 get_subsystem

    my $subsysObject = $sfx->get_subsystem($name);

Return a subsystem object for manipulation of the named subsystem. If the
subsystem does not exist, an undefined value will be returned.

=over 4

=item name

Name of the desired subsystem.

=item RETURN

Returns a blessed object that allows access to subsystem data.

=back

=cut

sub get_subsystem {
    # Get the parameters.
    my ($self, $name) = @_;
    # Declare the return value. If we don't change it, then the caller will know
    # the subsystem was not found.
    my $retVal;
    # Get the database.
    my $sprout = $self->{sprout};
    # Does the subsystem exist?
    if ($sprout->Exists(Subsystem => $name)) {
        # Yes, construct the subsystem object.
        $retVal = SproutSubsys->new($name, $self->{sprout});
    }
    return $retVal;
}

=head3 proteins_in_family

    my @pegs = $fig->in_family($family);

Return a list of the features in a specified protein family.

Sprout does not support protein families, so it calls through to the SEED method.

=over 4

=item family

Name of the protein family whose features are desired.

=item RETURN

Returns a list of the IDs of the features in the specified protein family.

=back

=cut
#: Return Type @;
sub proteins_in_family {
    # Get the parameters.
    my($self, $family) = @_;
    # Get the Sprout object.
    my $sprout = $self->{sprout};
    # Get a list of the PEGs for the specified family. Note that if the family
    # doesn't exist we'll simply get an empty list back.
    my @retVal = $sprout->GetFlat(['IsFamilyForFeature'], "IsFamilyForFeature(from-link) = ?",
                                  [$family], 'IsFamilyForFeature(to-link)');
    return @retVal;
}

=head3 families_for_protein

    my @families = $fig->families_for_protein($peg);

Return a list of all the families containing the specified protein.

=over 4

=item peg

ID of the PEG representing the protein in question.

=item RETURN

Returns a list of the IDs of the families containing the protein.

=back

=cut

sub families_for_protein {
    # Get the parameters.
    my ($self, $peg) = @_;
    # Get the Sprout object.
    my $sprout = $self->{sprout};
    # Read all the families for this protein.
    my @retVal = $sprout->GetFlat(['IsFamilyForFeature'], "IsFamilyForFeature(to-link) = ?",
                                  [$peg], 'IsFamilyForFeature(from-link)');
    return @retVal;
}


=head3 family_function

    my $func = $fig->family_function($family);

Returns the putative function of all of the pegs in a protein family.  Remember, we
are defining "protein family" as a set of homologous proteins that have the
same function.

=over 4

=item family

ID of the relevant protein family.

=item RETURN

Returns the name of the function assigned to the members of the specified family.

=back

=cut

sub family_function {
    # Get the parameters.
    my ($self, $family) = @_;
    # Get the Sprout object.
    my $sprout = $self->{sprout};
    # Find the family's function.
    my ($retVal) = $sprout->GetEntityValues('Family', $family, ['Family(function)']);
    # If it doesn't exist, return a null string.
    if (! defined $retVal) {
        $retVal = "";
    }
    return $retVal;
}

=head3 sz_family

    my $n = $fig->sz_family($family);

Returns the number of proteins in a family.

=over 4

=item family

ID of the relevant protein family.

=item RETURN

Returns the number of proteins in the specified family.

=back

=cut

sub sz_family {
    # Get the parameters.
    my ($self, $family) = @_;
    # Get the Sprout object.
    my $sprout = $self->{sprout};
    # Find the family's size.
    my ($retVal) = $sprout->GetEntityValues('Family', $family, ['Family(size)']);
    # If it doesn't exist, return 0.
    if (! defined $retVal) {
        $retVal = 0;
    }
    return $retVal;
}

=head3 sims

    my @sims = $fig->sims($pegs, $maxN, $maxP, $raw);

Returns a list of similarities for $peg governed by the constraints of the parameters.

=over 4

=item pegs

ID of the feature whose similarities are desired, or a reference to a list of feature IDs,
all of whose similarities are desired.

=item maxN

Maximum number of similarities to be performed.

=item maxP

Maximum allowable similarity score.

=item raw

If equal to C<raw>, then feature IDs are translated from external form to FIG
form.

=item RETURN

Returns a list of [[SimPm]] objects.

=back

=cut

sub sims {
    my ($self, $pegs, $maxN, $maxP, $raw) = @_;
    # Declare the return variable.
    my @retVal = ();
    # Ask the Sprout object for the similarities.
    my $sprout = $self->{sprout};
    my $simList = $sprout->Sims($pegs, $maxN, $maxP, $raw);
    # Check to see if we were successful.
    if (defined $simList) {
        Trace("Similarities returned from network.") if T(3);
        @retVal = @{$simList};
    } else {
        # The network failed us.
        Confess("Unable to retrieve similarities from network server.");
    }
    return @retVal;
}

=head3 get_genome_subsystem_data

    my $roleList = $fig->get_genome_subsystem_data($genomeID);

Return the roles and pegs for a genome's participation in subsystems. The
subsystem name, role ID, and feature ID will be returned for each of
the genome's subsystem-related PEGs. Only subsystems in which the genome
has a definite variant code (1 or greater) will be processed.

=over 4

=item genomeID

ID of the genome whose PEG breakdown is desired.

=item RETURN

Returns a list of 3-tuples. Each tuple consists of a subsystem name, a role ID,
and a feature ID.

=back

=cut

sub get_genome_subsystem_data {
    # Get the parameters.
    my ($self,$genomeID) = @_;
    my $sprout = $self->{sprout};
    # Create the return list.
    my @retVal = $sprout->GetAll("IsGenomeOf HasSSCell IsRoleOf ContainsFeature", 
                                 'IsGenomeOf(from-link) = ?', [$genomeID], 
                                 [qw(HasSSCell(from-link) IsRoleOf(from-link)
                                     ContainsFeature(to-link))]);
    # Return the result.
    return [ sort { FIGRules::FIGCompare($a->[2], $b->[2]) } @retVal ];
}

=head3 get_genome_assignment_data

    my $assignList = $fig->get_genome_subsystem_data($genomeID);

Return the functional assignments and pegs for a genome. The feature ID and assigned
function will be returned for each of the genome's PEGs.

=over 4

=item genomeID

ID of the genome whose PEG breakdown is desired.

=item RETURN

Returns a list of 2-tuples. Each tuple consists of a peg ID and its master
functional assignment.

=back

=cut

sub get_genome_assignment_data {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    # Get the sprout object.
    my $sprout = $self->{sprout};
    # Get the PEGs and annotations for this genome.
    my $featureHash = $sprout->GenomeAssignments($genomeID);
    # Return the result list.
    return [ map { [$_, $featureHash->{$_}] } sort { FIGRules::FIGCompare($a,$b) } keys %{$featureHash} ];
}

=head3 get_genome_stats

    my ($gname,$szdna,$pegs,$rnas,$taxonomy) = $fig->get_genome_stats($genomeID);

Return basic statistics about a genome.

=over 4

=item genomeID

ID of the relevant genome.

=item RETURN

Returns a 5-tuple containing the genome name, number of base pairs, number of PEG
features, number of RNA features, and the taxonomy string.

=back

=cut

sub get_genome_stats {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    my $sprout = $self->{sprout};
    my $genomeData = $sprout->GetEntity('Genome', $genomeID);
    my $name = join(" ", $genomeData->Values(['Genome(genus)', 'Genome(species)',
                                              'Genome(unique-characterization)']));
    my ($taxonomy) = $genomeData->Value('Genome(taxonomy)');
    my $szDNA = $sprout->GenomeLength($genomeID);
    my $numRNA = $sprout->FeatureCount($genomeID,"rna");
    my $numPEG = $sprout->FeatureCount($genomeID,"peg");
    return ($name, $szDNA, $numPEG, $numRNA, $taxonomy);
}

=head3 sort_genomes_by_taxonomy

    my @sortedGenomeIDs = $fig->sort_genomes_by_taxonomy(@genomeIDs);

Taxonomically sort a list of genome IDs.

=over 4

=item genomeIDs

The list of genome IDs to be sorted.

=item RETURN

Returns a list of the same genome IDs sorted by their taxonomy string.

=back

=cut

sub sort_genomes_by_taxonomy {
    # Get the parameters.
    my ($self, @genomeIDs) = @_;
    # Append the taxonomies.
    my @data = map { [$_, $self->taxonomy_of($_)] } @genomeIDs;
    # Sort the IDs by taxonomy.
    my @retVal = map { $_->[0] } sort { $a->[1] cmp $b->[1] } @data;
    # Return the result.
    return @retVal;
}

=head3 get_peg_keys_for_genome

    my @propertyList = $fig->get_peg_keys_for_genome($genomeID);

Return a list of all the properties for features of the specified genome. For each
property, the feature ID, key, value, and evidence will be returned.

=over 4

=item genomeID

ID of the genome whose properties are desired.

=item RETURN

Returns a list of 4-tuples. Each tuple contains a feature ID, a property name (key),
the property value, and the evidence string (which is usually a link).

=back

=cut

sub get_peg_keys_for_genome {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    my $sprout = $self->{sprout};
    # Get the property data.
    my @retVal = $sprout->GetAll(['HasFeature', 'HasProperty', 'Property'],
                                 "HasFeature(from-link) = ?", [$genomeID],
                                 ['HasFeature(to-link)', 'Property(property-name)',
                                  'Property(property-value)', 'HasProperty(evidence)']);
    # Return the result.
    return @retVal;
}

=head3 genome_version

    my $version = $fig->genome_version($genome_id);

Return the version number of the specified genome.

Versions are incremented for major updates.  They are put in as major
updates of the form 1.0, 2.0, ...

Users may do local "editing" of the DNA for a genome, but when they do,
they increment the digits to the right of the decimal.  Two genomes remain
comparable only if the versions match identically.  Hence, minor updating should be
committed only by the person/group responsible for updating that genome.

We can, of course, identify which genes are identical between any two genomes (by matching
the DNA or amino acid sequences).  However, the basic intent of the system is to
support editing by the main group issuing periodic major updates.

The genome version starts with the genome ID itself, so (for example) C<100226.1.1721812195> is a
variant of Streptomyces coelicolor version 1, while C<100226.2> would be an unmodified copy of
Streptomyces coelicolor version 2.

=over 4

=item genome_id

ID of the genome whose version is desired.

=item RETURN

Returns the version number of the specified genome, or C<undef> if the genome is not in
the data store.

=back

=cut

sub genome_version {
    # Get the parameters.
    my($self, $genome) = @_;
    # Extract the version string from the genome's database record.
    my ($retVal) = $self->{sprout}->GetFlat(['Genome'], "Genome(id) = ?", [$genome], 'Genome(version)');
    # Return the result.
    return $retVal;
}

=head3 subsystem_curator

    my $curator = $sfxlate->subsystem_curator($subsystem_name);

Return the curator of a subsystem.

=over 4

=item subsystem_name

Name of the subsystem whose curator is desired.

=item RETURN

Returns the name of the user who is assigned as the subsystem's curator, or C<undef> if the
subsystem is not found.

=back

=cut

sub subsystem_curator {
    my ($self, $subsystem_name) = @_;
    my ($retVal) = $self->{sprout}->GetFlat(['Subsystem'], "Subsystem(id) = ?",
                                            [$subsystem_name], 'Subsystem(curator)');
    return $retVal;
}

=head3 all_subsystems

    my @names = $fig->all_subsystems();

Return a list of all of the subsystems in the data store.

=cut

sub all_subsystems {
    # Get the parameters.
    my ($self) = @_;
    # Ask the Sprout for all subsystems.
    return $self->{sprout}->GetFlat(['Subsystem'], "", [], 'Subsystem(id)');
}

=head3 genome_rnas

    my $num_rnas = $fig->genome_rnas($genome_id);

Return the number of RNA-encoding genes for a genome.

=over 4

=item genome_id

ID of the genome whose RNA count is desired.

=item RETURN

Returns the number of RNAs for the specified genome, or C<undef> if the genome
is not indexed in the database.

=back

=cut

sub genome_rnas {
    # Get the parameters.
    my ($self, $genome_id) = @_;
    # Compute the count.
    return $self->genome_thing_count($genome_id, 'rna');
}

=head3 genome_pegs

    my $num_pegs = $fig->genome_pegs($genome_id);

Return the number of protein-encoding genes for a genome.

=over 4

=item genome_id

ID of the genome whose PEG count is desired.

=item RETURN

Returns the number of PEGs for the specified genome, or C<0> if the genome
is not indexed in the database.

=back

=cut

sub genome_pegs {
    # Get the parameters.
    my ($self, $genome_id) = @_;
    # Compute the count.
    return $self->genome_thing_count($genome_id, 'peg');
}

=head3 genome_thing_count

    my $num_pegs = $fig->genome_thing_count($genome_id, $type);

Return the number of genes for a genome of a specified type.

=over 4

=item genome_id

ID of the genome whose gene count is desired.

=item type

Type of gene to count-- C<peg>, C<rna>, etc.

=item RETURN

Returns the number of genes of the specified type for the specified genome, or C<0> if the genome
is not in the database.

=back

=cut

sub genome_thing_count {
    # Get the parameters.
    my ($self, $genome_id, $type) = @_;
    # Compute the count.
    my $retVal = $self->{sprout}->GetCount(['HasFeature'], "HasFeature(from-link) = ? AND HasFeature(type) = ?", [$genome_id, $type]);
    # Return the result.
    return $retVal;
}

=head3 genome_szdna

    my $szdna = $fig->genome_szdna($genome_id);

Return the number of DNA base-pairs in a genome's contigs.

=over 4

=item genome_id

ID of the genome whose base-pair count is desired.

=item RETURN

Returns the number of base pairs in the specified genome's contigs, or C<undef>
if the genome is not in the database.

=back

=cut

sub genome_szdna {
    # Get the parameters.
    my ($self, $genome_id) = @_;
    # Get the desired value from the genome record.
    my ($retVal) = $self->{sprout}->GetFlat(['Genome'], "Genome(id) = ?", [$genome_id], 'Genome(dna-size)');
    # Return the result.
    return $retVal;
}

=head3 orgid_of_orgname

    my $genomeID = $fig->orgid_of_orgname($genomeName);

Return the ID of the genome corresponding to the specified organism name, or a
null string if the genome is not found.

=over 4

=item genomeName

Name of the organism, consisting of the organism's genus, species, and
unique characterization, separated by spaces.

=item RETURN

Returns the genome ID number for the named organism, or an empty string if
the genome is not found.

=back

=cut

sub orgid_of_orgname {
    # Get the parameters.
    my ($self, $genomeName) = @_;
    # Split the name into its component parts. Note that we assume the genus and species
    # do not contain any internal spaces.
    my ($genus, $species, $subspecies) = split /\s+/, $genomeName, 3;
    # If no subspecies is specified, use an empty string.
    if (! defined($subspecies)) {
        $subspecies = '';
    }
    # Ask for a matching genome ID.
    my ($retVal, $other) = $self->{sprout}->GetFlat(['Genome'], "Genome(genus) = ? AND Genome(species) = ? AND Genome(unique-characterization) LIKE ?",
                                              [$genus, $species, "$subspecies%"], 'Genome(id)');
    # If there's no match or more than one match, return an empty string.
    if (! defined($retVal) || defined($other)) {
        $retVal = '';
    }
    # Return the result.
    return $retVal;
}

=head3 pegs_of

    my @pegs = $fig->pegs_of($genome)

Return a list of all PEGs in the specified genome.

=over 4

=item genome

ID of the relevant genome.

=item RETURN

Returns a list of the IDs of all PEG features in the genome.

=back

=cut

sub pegs_of {
    # Get the parameters.
    my ($self,$genome) = @_;
    # Find the pegs.
    my @retVal = $self->all_features($genome,"peg");
    # Return the result.
    return @retVal;
}

=head3 protein_subsystem_to_roles

    my $roles = $fig->protein_subsystem_to_roles($peg, $subsystem);

Return the roles played by a particular PEG in a particular subsytem. If the protein is not part of the
subsystem, an empty list will be returned.

=over 4

=item peg

ID of the protein whose role is desired.

=item subsystem

Name of the relevant subsystem.

=item RETURN

Returns a reference to a list of the roles performed by the specified PEG in the specified subsystem.

=back

=cut

sub protein_subsystem_to_roles {
    # Get the parameters.
    my ($self, $peg, $subsystem) = @_;
    Trace("Retrieving roles of peg $peg in $subsystem.") if T(3);
    # We find the subsystem spreadsheet cells containing the given PEG, and then extract each cell's role.
    my @roles = $self->{sprout}->GetFlat(['IsRoleOf', 'HasSSCell', 'ContainsFeature'],
                                          "HasSSCell(from-link) = ? AND ContainsFeature(to-link) = ?",
                                          [$subsystem, $peg], 'IsRoleOf(from-link)');
    # Remove duplicates.
    my %retVal = map { $_ => 1 } @roles;
    # Return the resulting list of roles.
    return [sort keys %retVal];
}

=head3 seqs_with_role

    my @pegs = $fig->seqs_with_role($role, $who, $genome);

Return a list of the pegs that implement a particular subsystem role.

=over 4

=item role

Subsystem role whose PEGs are desired.

=item who (optional)

This parameter is provided for compatibility with FIG. It is not used in the Sprout.

=item genome (optional)

If specified, only pegs in the specified genome will be returned.

=item RETURN

Returns a list of PEG IDs.

=back

=cut

sub seqs_with_role {
    # Get the parameters.
    my ($self, $role, $who, $genome) = @_;
    # Start with the base query of features by role.
    my @filters = "IsRoleOf(from-link) = ?";
    my @tables = qw(IsRoleOf ContainsFeature);
    my @parms = $role;
    # Check for genome filtering.
    if (defined $genome) {
        push @filters, "HasFeature(from-link) = ?";
        push @tables, 'HasFeature';
        push @parms, $genome;
    }
    # Exectue the query.
    my @pegs = $self->{sprout}->GetFlat(\@tables, join(" AND ", @filters), \@parms,
                                        'ContainsFeature(to-link)');
    return FIGRules::SortedFids(@pegs);
}

=head3 genus_species_domain

    my ($gs, $domain) = $fig->genus_species_domain($genome_id);

Returns a genome's genus and species (and strain if that has been properly
recorded) in a printable form, along with its domain. This method is similar
to L</genus_species>, except it also returns the domain name (archaea,
bacteria, etc.).

=over 4

=item genome_id

ID of the genome whose species and domain information is desired.

=item RETURN

Returns a two-element list. The first element is the species name and the
second is the domain name.

=back

=cut

sub genus_species_domain {
    # Get the parameters
    my ($self, $genome_id) = @_;
    # Declare the return value.
    my @retVal;
    # Get the genome record.
    my ($genomeData) = $self->{sprout}->GetAll(['Genome'], "Genome(id) = ?", [$genome_id],
                                               ['Genome(genus)', 'Genome(species)',
                                                'Genome(unique-characterization)',
                                                'Genome(taxonomy)']);
    if (! defined $genomeData) {
        # Here the genome does not exist.
        @retVal = ("", "");
    } else {
        # Build the genome name from its three parts.
        my $genomeName = join(" ", $genomeData->[0], $genomeData->[1], $genomeData->[2]);
        # Extract the domain name.
        my ($domain) = split /\s*;\s*/, $genomeData->[3], 1;
        # Put them together.
        @retVal = ($genomeName, $domain);
    }
    # Return the result.
    return @retVal;
}

=head3 wikipedia_link

    my $url = $fig->wikipedia_link($organism_name);

Return the URL of a Wikipedia page for the specified organism,
or C<undef> if no Wikipedia page exists.

=over

=item organism_name

Word or phrase to look for in Wikipedia.

=item RETURN

Returns the Wikipedia URL for the specified organism, or C<undef> if no Wikipedia
page for the organism exists.

=back

=cut

sub wikipedia_link {
    # Get the parameters.
    my ($self, $organism_name) = @_;
    # Return the result.
    return FIGRules::wikipedia_link($organism_name);
}

=head3 orgname_of_orgid

    my $genomeName = $fig->orgname_of_orgid($genomeID);

Return the name of the genome corresponding to the specified organism ID.

=over 4

=item genomeID

ID of the relevant genome.

=item RETURN

Returns the name of the organism, consisting of the organism's genus, species, and
unique characterization, separated by spaces, or a null string if the genome is not
found.

=back

=cut

sub orgname_of_orgid {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    # Get the sprout database object.
    my $sprout = $self->{sprout};
    # Ask for the genome name, storing a null string if it is not find.
    my $retVal = $sprout->GenusSpecies($genomeID) || "";
    # Return the result.
    return $retVal;
}

=head3 search_database

    my $resultList = $fig->search_database($query, $options);

Search the database for objects that match the query string in some way.

=over 4

=item query

A space-delimited list of keywords. The keywords can be genuine Sprout search
keywords, or could be the exact name of a feature, functional role, or
subsystem.

=item options

Reference to a hash of options. Currently, the only option is C<limit>, which puts
an upper limit on the number of results.

=item RETURN

Returns a hash reference. The key C<type> contains the result type and the key C<results> contains
either a single ID or a list of the resulting object IDs. For a type of C<organism>, the result is a genome ID.
For a type of <feature>, the result is a feature ID. For a type of C<subsystem>, the result is a subsystem name.
For a type of C<functional_role>, the result is a list of 2-tuples, each containing a role name followed by a
subsystem ID. For a type of C<proteins>, the result is a list of 3-tuples, containing the peg ID,
its functional role, and the ID of its parent genome. If nothing is found, the return will be an undefined
value.

=back

=cut

sub search_database {
    # Get the parameters.
    my ($self, $query, $options) = @_;
    # Declare the return variable.
    my $retVal;
    # Convert the query string into lower case. We don't use the keyword cleaner here, because we're going
    # to be using exact-match lookups on some things.
    my $lcQuery = lc($query);
    # Convert the query into a subsystem name.
    my $ss_query = $lcQuery;
    $ss_query =~ s/ /_/g;
    # Get the Sprout database.
    my $sprout = $self->{sprout};
    # Check for the exact organism id.
    my ($result) = $sprout->GetFlat(['Genome'], 'Genome(id) = ?', [$query], 'Genome(id)');
    if (defined($result)) {
        # If we found it, we can return an organism.
        Trace("Returning genome for ID $result.") if T(3);
        $retVal = { type => 'organism', result => $result };
    } else {
        # Check for the exact subsystem name.
        ($result) = $sprout->GetFlat(['Subsystem'], 'LOWER(Subsystem(id)) = ?', [$ss_query],
                                     'Subsystem(id)');
        if (defined($result)) {
            # If we found it, we have a subsystem.
            Trace("Returning subsystem for ID $result.") if T(3);
            $retVal = { type => 'subsystem', result => $result };
        } else {
            # Try to find a feature.
            ($result) = $sprout->GetFlat(['Feature'], 'Feature(id) = ?', [$lcQuery], 'Feature(id)');
            if (defined($result)) {
                Trace("Returning feature for ID $result.") if T(3);
                $retVal = { type => 'feature', result => $result };
            } else {
                # Ask for an organism by name.
                my $id = $self->orgid_of_orgname($query);
                if ($id) {
                    # We found one.
                    Trace("Returning genome $id for $query.") if T(3);
                    $retVal = { type => 'organism', result => $id };
                } else {
                    # Nothing for it here but to do a keyword search. Note that even if we're dealing with
                    # an alias, we need to recognize the possibility of multiple features returned,
                    # so this is a whole different kind of search. 
                    my $filter = "";
                    my @parms;
                    # Process the options.
                    if (defined $options) {
                        # If there are any other options, add them before this one, because
                        # the LIMIT clause always goes last in the filter string.
                        if ($options->{limit}) {
                            $filter .= " LIMIT ?";
                            push @parms, $options->{limit};
                        }
                    }
                    my $query = $sprout->Search($query, 0, ['Feature', 'HasFeature'], $filter, \@parms);
                    # Loop through the query, creating a result list.
                    my @results = ();
                    while (my $feature = $query->Fetch()) {
                        push @results, [$feature->Values(['Feature(id)', 'Feature(assignment)', 'HasFeature(from-link)'])];
                    }
                    # Only proceed if there are actual results.
                    if (@results) {
                        Trace("Returning list of " . scalar(@results) . " features from keyword search.") if T(3);
                        $retVal = { type => 'proteins', results => \@results };
                    }
                }
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 taxonomy_list

    my $taxonomyHash = $fig->taxonomy_list()

Return a reference to a hash that maps each genome ID to its taxonomy list. The taxonomy
lists are stored as strings delimited by semicolons.

=cut

sub taxonomy_list {
    # Get the parameters.
    my($self) = @_;
    # Form the taxonomies into a hash.
    my %retVal = map { $_->[0] => $_->[1] } $self->{sprout}->GetAll(['Genome'], "", [], [qw(Genome(id) Genome(taxonomy))]);
    return \%retVal;
}

=head3 adjacent_feature

    my $feature = $fig->adjacent_feature(\%options);

Locate the next or previous feature (optionally filtered by type) in a
contig. The start position for the search can be defined by supplying
genome, contig and position, or by supplying a feature id.  Feature
locations are defined by their midpoint.  If a fid is supplied
with contig and position, the latter are used to resolve ambiguities in
the desired segement of a feature with a complex location.

The following options are supported.

=over 4

=item after

ID of a feature that should preceed the returned feature, or a reference to
a list of the IDs of features that should preceed the returned feature. Note
that this is a local operation, and is only meant to resolve features that
are otherwise tied in location.

=item before

ID of a feature that should follow the returned feature, or a reference to
a list of the IDs of features that should follow the returned feature. This
option is mutually exclusive with C<after>.

=item contig

Name of the contig containing the features. If this option is omitted,
then C<fid> must be specified.

=item exclude

ID of a feature to exclude, or a reference to a list of IDs for features
to be excluded. Note that features listed with the C<after> or C<before>
option are also excluded (and that is most commonly the desired behavior).

=item fid

Alternative to supplying a location.  It is possible to supply a fid and
C<contig> and C<position>, which allows disambiguating the desired segment of
a feature with a complex location.

=item position

The feature midpoint must be >= $position (if C<after> is specified) or
<= position (if C<before> is specified). Note that this can be any multiple
of 1/2.  If the supplied value is negative, the position is taken from the
right end of the contig.

=item type

Type of desired feature, or reference to a list of permissible types (default is any type).

=item RETURN

Returns the ID of the desired feature, or C<undef> if no feature could be found.

=back

=cut

sub adjacent_feature {
    # Get the parameters.
    my ($self, $options) = @_;
    # Declare the return value.
    my $retVal;
    # Insure the option list is a hash.
    if (! $options || ref $options ne 'HASH') {
        Confess("Invalid parameter passed to adjacent_feature.");
    } else {
        # Determine the direction of movement.
        my $direction;
        if (exists $options->{before}) {
            $direction = 'before';
        } elsif (exists $options->{after}) {
            $direction = 'after';
        } else {
            Confess("No direction specified. There must be a \"before\" or \"after\" parameter.");
        }
        Trace("adjacent_feature direction is $direction.") if T(3);
        # Now we need to determine the position on the contig from which we'll be searching.
        my $contig = $options->{contig};
        my $position = $options->{position};
        my $fid = $options->{fid};
        my $genome;
        # Do we have a contig and a position?
        if ($contig && $position) {
            Trace("Contig $contig and position $position specified for adjacent_feature.") if T(3);
            # Yes, compute the genome ID.
            $contig =~ /([^:]+):/;
            $genome = $1;
            # If the position is negative, count it from the end.
            if ($position < 0) {
                $position += $self->contig_ln($genome, $contig) + 1;
            }
        } elsif ($fid) {
            Trace("Feature $fid specified for adjacent_feature.") if T(3);
            # Here we are using a feature ID to create the position. Extract the genome ID from the feature ID.
            my ($genome) = FIGRules::ParseFeatureID($fid);
            # Get a location object for the feature's locations.
            my $loc = FullLocation->new($self, $genome, [$self->feature_location($fid)]);
            # Take the feature's midpoint as the start position.
            $contig = $loc->Contig();
            $position = ($loc->Begin() + $loc->EndPoint()) / 2;
        } else {
            # Here we're stuck.
            Confess("Unable to determine start position. You must specify either a contig ID (contig) and a position (position), or a feature ID (fid).");
        }
        Trace("Contig is $contig and final position is $position.") if T(3);
        # Create a hash of the features in the parameter relevant to the direction. Note that we wouldn't be here
        # if the direction wasn't specified, so we don't need to verify the existence of the parameter. Each feature
        # will be mapped to its contig boundaries.
        my $fidList = $options->{$direction};
        if (ref $fidList ne 'ARRAY') {
            $fidList = [$fidList];
        }
        my %fidMap = ();
        for my $fid (@{$fidList}) {
            my $loc = FullLocation->new($self, $genome, [$self->feature_location($fid)]);
            my ($bounds) = $loc->GetBounds();
            $fidMap{$fid} = $bounds;
        }
        # Create the exclude list.
        my %exclusions = ();
        my $excludes = $options->{exclude};
        if ($excludes) {
            if (ref $excludes eq 'ARRAY') {
                for my $excludeFid (@{$excludes}) {
                    $exclusions{$excludeFid} = 1;
                }
            } else {
                $exclusions{$excludes} = 1;
            }
        }
        # Now get the permissible-type list.
        my @types = ();
        my $type = $options->{type};
        if ($type) {
            if (ref $type ne 'ARRAY') {
                $type = [$type];
            }
            @types = @{$type};
        }
        # If $type is TRUE, then this hash will be used to verify the feature type. Otherwise,
        # all feature types are permitted.
        my %types = map { $_ => 1 } @types;
        # Get the length of the contig.
        my $contigLength = $self->contig_ln($genome, $contig);
        # Compute the minimum and maximum positions.
        my ($minV, $maxV, $step);
        if ($direction eq 'before') {
            $maxV = int($position);
            $minV = $maxV - 9999;
            $step = -10000;
        } else {
            $minV = int($position + 0.5);
            $maxV = $minV + 9999;
            $step = 10000;
        }
        # Get the Sprout object.
        my $sprout = $self->{sprout};
        # Loop until we find something or run out of contig.
        while (! defined $retVal && (FIG::between(0, $minV, $contigLength) || FIG::between(0, $maxV, $contigLength))) {
            # Get the genes in the specified region.
            my @features = $sprout->GeneDataInRegion($contig, $minV, $maxV);
            # Create a hash of the features found and their locations.
            my %featureData = ();
            for my $feature (@features) {
                my $fid = $feature->PrimaryValue('Feature(id)');
                my $floc = FullLocation->new($self, $genome, $feature->PrimaryValue('Feature(location-string)'));
                my ($bloc) = $floc->GetBounds();
                Trace("Feature $fid has bounds " . $bloc->String() . ".") if T(4);
            }
            # Now we look for the best feature in the lot. It will be the one with the closest distance.
            my $bestDistance = $contigLength;
            for my $fid (keys %featureData) {
                # Check this feature's distance.
                my $floc = $featureData{$fid};
                my $distance = abs(($floc->Begin + $floc->EndPoint)/2 - $position);
                Trace("Distance for $fid is $distance. Best is $bestDistance.") if T(4);
                # If it's no worse than the last one, keep checking.
                if ($distance <= $bestDistance) {
                    # Parse the feature ID.
                    my ($thisGenome, $thisType, $thisNum) = FIGRules::ParseFeatureID($fid);
                    # Only proceed if the type is valid and we're not excluded.
                    if (! $exclusions{$fid} && (! $type || $types{$thisType})) {
                        # Check the direction list.
                        my $ok = 1;
                        for my $dirFid (keys %fidMap) {
                            my $dirLoc = $fidMap{$dirFid};
                            if ($direction eq 'before' && $floc->Right >= $dirLoc->Left ||
                              $direction eq 'after' && $floc->Left <= $dirLoc->Right) {
                                $ok = 0;
                            }
                        }
                        # If we passed, save this feature.
                        if ($ok) {
                            $bestDistance = $distance;
                            $retVal = $fid;
                        }
                    }
                }
            }
            # Move the boundaries in case we need to try again.
            $minV += $step;
            $maxV += $step;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 next_feature

This method calls through to L</adjacent_feature>. It creates an C<after> option if none exists.

=cut

sub next_feature {
    my ($self, $options) = @_;
    Trace("next_feature called.") if T(3);
    if (! exists $options->{after}) {
        $options->{after} = [];
    }
    return $self->adjacent_feature($options);
}

=head3 previous_feature

This method calls through to L</adjacent_feature>. If creates a C<before> option if none exists.

=cut

sub previous_feature {
    my ($self, $options) = @_;
    Trace("previous_feature called.") if T(3);
    if (! exists $options->{before}) {
        $options->{before} = [];
    }
    return $self->adjacent_feature($options);
}

=head3 get_representative_genome

    my $rep_id = $fig->get_representative_genome($id);

Return the ID of a specified genome's representative, or an undefined
value if it has no representative. Genomes are divided into sets of
close strains; the representative is the one considered to be the
most accurately annotated.

=over 4

=item id

ID of the relevant genome.

=item RETURN

The ID of the representative for the specified genome, or C<undef> if
the genome is standalone (has no close strains available).

=back

=cut

sub get_representative_genome {
    # Get the parameters.
    my ($self, $id) = @_;
    # Ask for the representative. If none is found, we'll get an undef.
    my ($retVal) = $self->{sprout}->GetFlat("IsRepresentativeOf",
                                            "IsRepresentativeOf(to-link) = ?",
                                            [$id], "from-link");
    # Return the result.
    return $retVal;
}

=head3 get_corresponding_ids

    my @id_list = $fig->get_corresponding_ids($id, $with_type_info);

Return a list of the identifiers that correspond to the given identifier, based on 
the  PIR id correspondence table.

=over 4

=item id

Identifer to look up.

=item with_type_info

Pass a true value here to return tuples [id, source-type, link-information] instead of identifiers.

=item RETURN

A list of identifiers if $with_type_info not true; a  list of tuples 
[id, source-type, link-information] otherwise.

=back

=cut

sub get_corresponding_ids {
    my ($self, $id, $with_type_info) = @_;
    # Get the database.
    my $sprout = $self->{sprout};
    # We need to drive this from the FIG side of things, so our first task is to convert the incoming
    # ID to a list of FIG IDs.
    my @fids;
    if ($id =~ /^fig/) {
        # The incoming ID is a FIG ID.
        @fids = $id;
    } else {
        # We need to get the FIG IDs from the database.
        my @fids = $sprout->GetFlat(['IsAlsoFoundIn'], "IsAlsoFoundIn(alias) = ?", [$id],
                                    "IsAlsoFoundIn(from-link)");
    }
    # The aliases we find will be put in here. Each alias will map to its database name. We prime it with the FIG
    # IDs found so far (usually there's only one).
    my %retVal = map { $_ => 'SEED' } @fids;
    # Only proceed if we have something to look for.
    if (@fids) {
        # Create a filter for the specified FIG IDs.
        my $where = join(" AND ", map { "IsAlsoFoundIn(from-link) = ?" } @fids);
        # Get all the external IDs.
        my @others = $sprout->GetAll(['IsAlsoFoundIn'], $where, \@fids, ['IsAlsoFoundIn(alias)', 'IsAlsoFoundIn(to-link)']);
        # Add them to the hash.
        for my $other (@others) {
            $retVal{$other->[0]} = $other->[1];
        }
    }
    # How we return the result depends on $with_type_info.
    my @retVal;
    if ($with_type_info) {
        @retVal = map { [ $_, $retVal{$_} ] } keys %retVal;
    } else {
        @retVal = keys %retVal;
    }
    return @retVal;
}

=head3 clearinghouse_register_subsystem_id

    my $tax = $fig->clearinghouse_register_subsystem_id($ss_name);

Return a subsystem's short ID. Short IDs are maintained at a special
clearinghouse web site. If the subsystem does not yet have a short ID, a
new one will be assigned by the clearinghouse and returned.

=over 4

=item ss_name

Name of the subsystem whose ID is desired.

=item RETURN

ID of the desired subsystem.

=back

=cut

sub clearinghouse_register_subsystem_id {
    # Get the parameters.
    my ($self, $ss_name) = @_;
    # Return the result.
    return FIGRules::clearinghouse_register_subsystem_id($ss_name);
}

=head3 subsystem_version

    my $version = $fig->subsystem_version($ssn);

Return the version string for a subsystem.

=over 4

=item ssn

Name of the subsystem.

=item RETURN

Returns the subsystem's version string.

=back

=cut

sub subsystem_version {
    # Get the parameters.
    my ($self, $ssn) = @_;
    # Read the version number.
    my ($retVal) = $self->{sprout}->GetEntityValues(Subsystem => $ssn, ['Subsystem(version)']);
    # Return the result.
    return $retVal;
}

=head3 is_locked_fid

    my $flag = $fig->is_locked_fid($fid);

Return 1 if a feature is locked, else 0.

=over 4

=item fid

ID of the relevant feature.

=item RETURN

Returns 1 if the feature is locked, or 0 if it is not locked.

=back

=cut

sub is_locked_fid {
    # Get the parameters.
    my ($self, $fid) = @_;
    # Declare the return variable.
    my ($retVal) = $self->{sprout}->GetEntityValues(Feature => $fid, ['Feature(locked)']);
    # Return the result.
    return $retVal;
}

=head3 peg_in_gendb

    my $flag = $sfxlate->peg_in_gendb($featureID);

Return TRUE if the specified feature is in the GenBank database, else FALSE.

As currently implemented, this function always returns FALSE.

=over 4

=item featureID

ID of the feature whose GenBank status is desired.

=item RETURN

Returns TRUE if the feature is known to be in the GenBank database, else FALSE.

=back

=cut
#: Return Type $;
sub peg_in_gendb {
    # Get the parameters.
    my ($self, $featureID) = @_;
    # Declare the return variable.
    my $retVal = $self->{sprout}->GetEntityValues(Feature => $featureID, ['Feature(in-genbank)']);
    # Return the result.
    return $retVal;
}

=head3 translation_length

    my $len = $sfxlate->translation_length($featureID);

Return the length of the specified feature's translation.

=over 4

=item featureID

ID of the feature whose translation length is desired.

=item RETURN

The number of letters in the translation string. This will often be the number of proteins
in the feature's transcription.

=back

=cut
#: Return Type $;
sub translation_length {
    # Get the parameters.
    my ($self, $featureID) = @_;
    # Declare the return variable.
    my $retVal;
    # Get the translation.
    my ($translation) = $self->{sprout}->GetEntityValues(Feature => $featureID, ['Feature(translation)']);
    # Only proceed if the translation exists.
    if (defined $translation) {
        $retVal = length($translation);
    }
    # Return the result.
    return $retVal;
}

=head3 in_pch_pin_with

    my @pegs = $sfxlate->in_pch_pin_with($peg);

Return the IDs of features that are believed to be pinned to I<$peg> (in the
sense that PCHs occur containing these pegs over significant phylogenetic
distances).

This method is currently implemented by a call through to SEED that
connects to the pin server.

=over 4

=item peg

ID of the feature whose pinned pegs are to be returned.

=item RETURN

Returns a list of the IDs of the features pinned to the specified feature.

=back

=cut

sub in_pch_pin_with {
    my($self, $peg) = @_;
    # Get the pins from the pin server.
    my @pins = FIGRules::NetCouplingData('in_pch_pin_with_and_evidence', id1 => $peg);
    # Extract the peg IDs.
    return map { $_->[0] } @pins;
}

=head3 mapped_prot_ids

    my @mapped = $fig->mapped_prot_ids($fid);

Return a list of the features mapped to the specified feature's protein.
This information is not in the Sprout database, so it is implemented as a call through
to the SEED. The features are returned as 2-tuples containing the feature ID and the length,
sorted by length.

=over 4

=item fid

ID of the feature whose mapped synonyms are desired.

=item RETURN

Returns list of 2-tuples, each consisting of a mapped feature's ID and its length,
sorted by length.

=back

=cut
#: Return Type @@;
sub mapped_prot_ids {
    my ($self, $fid) = @_;
    return $self->FIG()->mapped_prot_ids($fid);
}


=head3 is_aux_role_in_subsystem

    my $flag = $fig->is_aux_role_in_subsystem($subsystem, $role);

Return TRUE if the specified role in the specified subsystem is
auxiliary.

=over 4

=item subsystem

Subsystem relevant to the role.

=item role

Name of the role in question.

=item RETURN

Returns TRUE if the role is auxiliary, else FALSE.

=back

=cut

sub is_aux_role_in_subsystem {
    # Get the parameters.
    my ($self, $subsystem, $role) = @_;
    # Declare the return variable.
    my ($retVal) = $self->{sprout}->GetFlat(['OccursInSubsystem'], "OccursInSubsystem(from-link) = ? AND OccursInSubsystem(to-link) = ?",
                                            [$subsystem, $role], 'OccursInSubsystem(auxiliary)');
    # Return the result.
    return $retVal;
}


=head2 Stubs and Limited Methods

=head3 translate_function

    my $translation = $sfxlate->translate_function($function);

Return a synonym for the specified function.

Function synonyms are not supported by Sprout, so this method invokes the SEED routine.

=over 4

=item function

Name of the function to translate.

=item RETURN

Returns a synonym for the specified function.

=back

=cut

sub translate_function {
    my ($self, $function) = @_;
    return $self->FIG()->translate_function($function);
}

=head3 is_exchangable_subsystem

    my $flag = $fig->is_exchangable_subsystem($ssa);

Return TRUE if the specified subsystem is exchangeable, else FALSE. In
the Sprout, all subsystems are considered exchangeable, so this method
always returns TRUE.

=over 4

=item ssa

ID of the relevant subsystem.

=item RETURN

Returns TRUE.

=back

=cut

sub is_exchangable_subsystem {
    # Get the parameters.
    my ($self, $ssa) = @_;
    # Declare the return variable.
    my $retVal = 1;
    # Return the result.
    return $retVal;
}

=head3 scenario_directory

    my $directoryName = $fig->scenario_directory($organism);

Return the scenario directory of an organism.  If the organism is 'All', return
the directory containing all possible paths through scenarios.

=over 4

=item organism

ID of the desired organism, or 'All'.

=item RETURN

Returns the name of the desired directory.

=back

=cut

sub scenario_directory {
    # Get the parameters.
    my ($self, $organism) = @_;
    # Pick the desired directory.
    my $retVal = "";
    if ($organism eq 'All') {
	$retVal .= "$FIG_Config::global/Models/All/Scenarios";
    }
    else {
	$retVal .= "$FIG_Config::organisms/$organism/Scenarios";
    }
    # Return the result.
    return $retVal;
}

=head3 organism_directory

    my $organism_directory = $fig->organism_directory($genome_id);

Get the directory that contains the organism data. This is just like the
FIGV version.

=over 4

=item genome_id

The id of the organism, e.g. 83333.1.

=item RETURN

A string containing the path to the organism directory.

=back

=cut

sub organism_directory {
    my ($self, $org_id) = @_;
    return "$FIG_Config::organisms/$org_id";
}

=head3 table_exists

    my $flag = $sfxlate->table_exists($table_name);

Return TRUE if the specified table exists, else FALSE. Since the table names in Sprout
are not compatible with the SEED table names, this method currently returns TRUE for
everything. This is because most of the optional capabilities in SEED are available in Sprout;
it's only the required capabilities that are missing.

=over 4

=item table_name

Name of the table of interest.

=item RETURN

Returns TRUE if the table's capability is supported, else FALSE.

=back

=cut
#: Return Type $;
sub table_exists {
    # Get the parameters.
    my ($self, $table_name) = @_;
    return 1;
}

=head3 change_attribute

    my $okFlag = $sfxlate->change_attribute($peg, $key, $value, $url);

Change the value of a property. This method has no effect because Sprout cannot
do updates.

=over 4

=item peg

ID of the feature whose property value is to be changed.

=item key

Name of the property whose value is to be changed.

=item value

Proposed new value for the property.

=item url

URL or citation used as evidence of the new value.

=item RETURN

Returns 1 if successful, 0 if an error occurred.

=back

=cut

sub change_attribute {
    my ($self, $peg, $key, $value, $url) = @_;
    Trace("Attribute updates not supported in Sprout.") if T(0);
    return 0;
}

=head3 delete_attribute

    my $okFlag = $sfxlate->delete_attribute($peg, $key);

Delete the specified property's value for the specified feature. This method
has no effect because Sprout cannot do deletions.

=over 4

=item peg

ID of the feature whose property value is to be deleted.

=item key

Name of the property to be deleted.

=item RETURN

Returns 1 if successful, 0 if an error occurred.

=back

=cut

sub delete_attribute {
    my ($self, $peg, $key) = @_;
    Trace("Attribute deletes not supported in Sprout.") if T(0);
    return 0;
}

=head2 Static Utility Methods

=head3 roles_of_function

    my @roles = $sfxlate->roles_of_function($func);

Returns a list of the functional roles implemented by the specified function. This method
parses the role data out of the function name, and does not require access to the database.

=over 4

=item func

Name of the function whose roles are to be parsed out.

=item RETURN

Returns a list of the roles performed by the specified function.

=back

=cut

sub roles_of_function {
    # Remove the SELF parameter.
    shift @_;
    # Call the FIG method.
    return FIG::roles_of_function(@_);
}

=head3 max

    my $max = FIG::max(@x);

or

    my $max = $sfx->max(@x);

Return the maximum numeric value from a list.

=over 4

=item x1, x2, ... xN

List of numbers to process.

=item RETURN

Returns the numeric value of t/he list entry possessing the highest value. Returns
C<undef> if the list is empty.

=back

=cut
#: Return Type $;
sub max {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    return FIG::max(@_);
}

=head3 min

    my $min = FIG::min(@x);

or

    my $min = $sfx->min(@x);

Return the minimum numeric value from a list.

=over 4

=item x1, x2, ... xN

List of numbers to process.

=item RETURN

Returns the numeric value of the list entry possessing the lowest value. Returns
C<undef> if the list is empty.

=back

=cut
#: Return Type $;
sub min {
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);
    return FIG::min(@_);
}

=head3 FIG

    my $realFIG = $sfx->FIG();

Return a FIG object. If a cached FIG object has not yet been created, one will be.
This is more or less an emergency measure to provide functionality that is in the
SEED system but not in Sprout.

=cut

sub FIG {
    my ($self) = @_;
    my $retVal = $self->{FIG};
    if (! $retVal) {
        # This is the first call, so we create a
        $retVal = FIG->new();
        $self->{FIG} = $retVal;
    }
    return $retVal;
}

=head3 openF

    my $fh = $sfx->openF($fileName);

Open a file for input. The FIG version of this method uses it to cache file handles in
memory and insure only a limited number of files are left open at any given time.
Since almost all the data has been moved from files to the database, however, this is
no longer a problem.

=over 4

=item fileName

Name of the file to open.

=item RETURN

Returns a handle for the open file. If the file does not open, an exception will be thrown.

=back

=cut

sub openF {
    # Get the parameters.
    my ($self, $fileName) = @_;
    # Open the file and create a new handle for it.
    my $retVal = Open(undef, "<$fileName");
    # Return the file handle.
    return $retVal;
}

=head3 clean_spaces

Remove any extra spaces from input fields. This will (currently) remove ^\s, \s$,
and concatenate multiple spaces into one.

    my $input = $sfx->clean_spaces($s);

This method is technically static, in that it does not use the FIG, SFXlate, or
Sprout objects at all; however, it has the signature of an instance method (i.e.
it is invoked using the arrow operator). This method simply does a little dance
with the parameters to get to the FIG space-cleaning method.

=over 4

=item s

String to clean.

=item RETURN

Returns a version of the input string with multiple spaces reduced to a single
space and leading and trailing white space removed.

=back

=cut

sub clean_spaces {
    my ($self, $s) =@_;
    return FIG::clean_spaces($self, $s);
}

=head3 Title

    my $title = $fig->Title();

Return the title of this database. For SEED, this will return SEED, for Sprout
it will return NMPDR, and so forth.

=cut

sub Title {
    return "NMPDR";
}

1;
