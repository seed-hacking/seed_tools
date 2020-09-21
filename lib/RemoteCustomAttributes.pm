#!/usr/bin/perl -w

package RemoteCustomAttributes;

    use strict;
    use Tracer;
    use PageBuilder;
    use Frontier::Client;

=head1 Remote Custom Attribute Manager

=head2 Introduction

This package looks like a limited version of the B<CustomAttributes> object, but
it gets the data from a remote server using the Frontier XML interface. The
methods supported are the four FIG replacement methods of the custom attributes object.

=head2 Public Methods

=head3 new

    my $rca = RemoteCustomAttributes->new($url);

Construct a new RemoteCustomAttributes object served by the script at the
specified URL.

=cut

sub new {
    # Get the parameters.
    my ($class, $url) = @_;
    # Create the rca object.
    my $retVal = {
                  proxy => Frontier::Client->new(url => $url)
                 };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head3 GetAttributes

    my @attributeList = $attrDB->GetAttributes($objectID, $key, @valuePatterns);

In the database, attribute values are sectioned into pieces using a splitter
value specified in the constructor (L</new>). This is not a requirement of
the attribute system as a whole, merely a convenience for the purpose of
these methods. If you are using the static method calls instead of the
object-based calls, the splitter will always be the default value of
double colons (C<::>). If a value has multiple sections, each section
is matched against the correspond criterion in the I<@valuePatterns> list.

This method returns a series of tuples that match the specified criteria. Each tuple
will contain an object ID, a key, and one or more values. The parameters to this
method therefore correspond structurally to the values expected in each tuple.

    my @attributeList = GetAttributes('fig|100226.1.peg.1004', 'structure%', 1, 2);

would return something like

    ['fig}100226.1.peg.1004', 'structure', 1, 2]
    ['fig}100226.1.peg.1004', 'structure1', 1, 2]
    ['fig}100226.1.peg.1004', 'structure2', 1, 2]
    ['fig}100226.1.peg.1004', 'structureA', 1, 2]

Use of C<undef> in any position acts as a wild card (all values). In addition,
the I<$key> and I<@valuePatterns> parameters can contain SQL pattern characters: C<%>, which
matches any sequence of characters, and C<_>, which matches any single character.
(You can use an escape sequence C<\%> or C<\_> to match an actual percent sign or
underscore.)

In addition to values in multiple sections, a single attribute key can have multiple
values, so even

    my @attributeList = GetAttributes($peg, 'virulent');

which has no wildcard in the key or the object ID, may return multiple tuples.

For reasons of backward compatability, we examine the structure of the object ID to
determine the entity type. In that case the only two types allowed are C<Genome> and
C<Feature>. An alternative method is to use a list reference, with the list consisting
of an entity type name and the actual ID. Thus, the above example could equivalently
be written as

    my @attributeList = GetAttributes([Feature => $peg], 'virulent');

The list-reference approach allows us to add attributes to other entity types in
the future. Doing so, however, will require modifying the L</Refresh> method and
updated the database design XML.

The list-reference approach also allows for a more fault-tolerant approach to
getting all objects with a particular attribute.

    my @attributeList = GetAttributes([Feature => undef], 'virulent');

will only return feature attributes, while

    my @attributeList = GetAttributes(undef, 'virulent');

could at some point in the future get you attributes for genomes or even subsystems
as well as features.

=over 4

=item objectID

ID of the genome or feature whose attributes are desired. In general, an ID that
starts with C<fig|> is treated as a feature ID, and an ID that is all digits with a
single period is treated as a genome ID. For other entity types, use a list reference; in
this case the first list element is the entity type and the second is the ID. A value of
C<undef> here will match all objects.

=item key

Attribute key name. Since attributes are stored as fields in the database with a
field name equal to the key name, it is very fast to find a list of all the
matching keys. Each key's values require a separate query, however, which may
be a performance problem if the pattern matches a lot of keys. Wild cards are
acceptable here, and a value of C<undef> will match all attribute keys.

=item valuePatterns

List of the desired attribute values, section by section. If C<undef>
is specified, all values in that section will match.

=item RETURN

Returns a list of tuples. The first element in the tuple is an object ID, the
second is an attribute key, and the remaining elements are the sections of
the attribute value. All of the tuples will match the criteria set forth in
the parameter list.

=back

=cut

sub GetAttributes {
    # Get the parameters.
    my ($self, $objectID, $key, @valuePatterns) = @_;
    Trace("Calling the remote custom attributes server.") if T(3);
    # Call the remote server.
    my $result = $self->{proxy}->call('GetAttributes', $objectID, $key, @valuePatterns);
    Trace("Result type is " . ref($result) . ".") if T(3);
    Trace(scalar(@{$result}) . " rows found.") if T(3) && ref($result) eq 'ARRAY';
    # Unwrap the list and return it.
    my @retVal = @{$result};
    return @retVal;
}

=head3 QueryAttributes

    my @attributeData = $ca->QueryAttributes($filter, $filterParms);

Return the attribute data based on an SQL filter clause. In the filter clause,
the name C<$object> should be used for the object ID, C<$key> should be used for
the key name, C<$subkey> for the subkey value, and C<$value> for the value field.

=over 4

=item filter

Filter clause in the standard ERDB format, except that the field names are C<$object> for
the object ID field, C<$key> for the key name field, C<$subkey> for the subkey field,
and C<$value> for the value field. This abstraction enables us to hide the details of
the database construction from the user.

=item filterParms

Parameters for the filter clause.

=item RETURN

Returns a list of tuples. Each tuple consists of an object ID, a key (with optional subkey), and
one or more attribute values.

=back

=cut

sub QueryAttributes {
    # Get the parameters.
    my ($self, $filter, $filterParms) = @_;
    Trace("Calling the remote custom attributes server.") if T(3);
    # Call the remote server.
    my $result = $self->{proxy}->call('QueryAttributes', $filter, $filterParms);
    Trace("Result type is " . ref($result) . ".") if T(3);
    Trace(scalar(@{$result}) . " rows found.") if T(3) && ref($result) eq 'ARRAY';
    # Unwrap the list and return it.
    my @retVal = @{$result};
    return @retVal;
}

=head3 AddAttribute

    $attrDB->AddAttribute($objectID, $key, @values);

Add an attribute key/value pair to an object. This method cannot add a new key, merely
add a value to an existing key. Use L</StoreAttributeKey> to create a new key.

=over 4

=item objectID

ID of the genome or feature to which the attribute is to be added. In general, an ID that
starts with C<fig|> is treated as a feature ID, and an ID that is all digits and periods
is treated as a genome ID. For IDs of other types, this parameter should be a reference
to a 2-tuple consisting of the entity type name followed by the object ID.

=item key

Attribute key name. This corresponds to the name of a field in the database.

=item values

One or more values to be associated with the key. The values are joined together with
the splitter value before being stored as field values. This enables L</GetAttributes>
to split them apart during retrieval. The splitter value defaults to double colons C<::>.

=back

=cut

sub AddAttribute {
    # Get the parameters.
    my ($self, $objectID, $key, @values) = @_;
    # Don't allow undefs.
    if (! defined($objectID)) {
        Confess("No object ID specified for AddAttribute call.");
    } elsif (! defined($key)) {
        Confess("No attribute key specified for AddAttribute call.");
    } elsif (! @values) {
        Confess("No values specified in AddAttribute call for key $key.");
    } else {
        # Call the remote server.
        my $result = $self->{proxy}->call('AddAttribute', $objectID, $key, @values);
    }
    # Return a one. We do this for backward compatability.
    return 1;
}

=head3 DeleteAttribute

    $attrDB->DeleteAttribute($objectID, $key, @values);

Delete the specified attribute key/value combination from the database.

The first form will connect to the database and release it. The second form
uses the database connection contained in the object.

=over 4

=item objectID

ID of the genome or feature to which the attribute is to be added. In general, an ID that
starts with C<fig|> is treated as a feature ID, and an ID that is all digits and periods
is treated as a genome ID. For IDs of other types, this parameter should be a reference
to a 2-tuple consisting of the entity type name followed by the object ID.

=item key

Attribute key name. This corresponds to the name of a field in the database.

=item values

One or more values to be associated with the key.

=back

=cut

sub DeleteAttribute {
    # Get the parameters.
    my ($self, $objectID, $key, @values) = @_;
    # Don't allow undefs.
    if (! defined($objectID)) {
        Confess("No object ID specified for DeleteAttribute call.");
    } elsif (! defined($key)) {
        Confess("No attribute key specified for DeleteAttribute call.");
    } elsif (! @values) {
        Confess("No values specified in DeleteAttribute call for key $key.");
    } else {
        # Call the remote server.
        my $result = $self->{proxy}->call('DeleteAttribute', $objectID, $key, @values);
    }
    # Return a one. This is for backward compatability.
    return 1;
}

=head3 ChangeAttribute

    $attrDB->ChangeAttribute($objectID, $key, \@oldValues, \@newValues);

Change the value of an attribute key/value pair for an object.

=over 4

=item objectID

ID of the genome or feature to which the attribute is to be changed. In general, an ID that
starts with C<fig|> is treated as a feature ID, and an ID that is all digits and periods
is treated as a genome ID. For IDs of other types, this parameter should be a reference
to a 2-tuple consisting of the entity type name followed by the object ID.

=item key

Attribute key name. This corresponds to the name of a field in the database.

=item oldValues

One or more values identifying the key/value pair to change.

=item newValues

One or more values to be put in place of the old values.

=back

=cut

sub ChangeAttribute {
    # Get the parameters.
    my ($self, $objectID, $key, $oldValues, $newValues) = @_;
    # Don't allow undefs.
    if (! defined($objectID)) {
        Confess("No object ID specified for ChangeAttribute call.");
    } elsif (! defined($key)) {
        Confess("No attribute key specified for ChangeAttribute call.");
    } elsif (! defined($oldValues) || ref $oldValues ne 'ARRAY') {
        Confess("No old values specified in ChangeAttribute call for key $key.");
    } elsif (! defined($newValues) || ref $newValues ne 'ARRAY') {
        Confess("No new values specified in ChangeAttribute call for key $key.");
    } else {
        # Call the remote server.
        my $result = $self->{proxy}->call('ChangeAttribute', $objectID, $key, $oldValues, $newValues);
    }
    # Return a one. We do this for backward compatability.
    return 1;
}

=head3 EraseAttribute

    $attrDB->EraseAttribute($entityName, $key);

Erase all values for the specified attribute key. This does not remove the
key from the database; it merely removes all the values.

=over 4

=item entityName

Name of the entity to which the key belongs. If undefined, all entities will be
examined for the desired key.

=item key

Key to erase.

=back

=cut

sub EraseAttribute {
    # Get the parameters.
    my ($self, $entityName, $key) = @_;
    # Call the remote server.
    my $result = $self->{proxy}->call('EraseAttribute', $entityName, $key);
    # Return a 1, for backward compatability.
    return 1;
}

=head3 GetAttributeKeys

    my @keyList = $attrDB->GetAttributeKeys($entityName);

Return a list of the attribute keys for a particular entity type.

=over 4

=item entityName

Name of the entity whose keys are desired.

=item RETURN

Returns a list of the attribute keys for the specified entity.

=back

=cut

sub GetAttributeKeys {
    # Get the parameters.
    my ($self, $entityName) = @_;
    # Call the remote server.
    my $result = $self->{proxy}->call('GetAttributeKeys', $entityName);
    Trace("Result type is " . ref($result) . ".") if T(3);
    Trace(scalar(@{$result}) . " keys found.") if T(3) && ref($result) eq 'ARRAY';
    # Unwrap the list and return it.
    my @retVal = @{$result};
    return @retVal;
}

=head3 DeleteMatchingAttributes

    my @deleted = $attrDB->DeleteMatchingAttributes($objectID, $key, @values);

Delete all attributes that match the specified criteria. This is equivalent to
calling L</GetAttributes> and then invoking L</DeleteAttribute> for each
row found.

=over 4

=item objectID

ID of object whose attributes are to be deleted. If the attributes for multiple
objects are to be deleted, this parameter can be specified as a list reference. If
attributes are to be deleted for all objects, specify C<undef> or an empty string.
Finally, you can delete attributes for a range of object IDs by putting a percent
sign (C<%>) at the end.

=item key

Attribute key name. A value of C<undef> or an empty string will match all
attribute keys. If the values are to be deletedfor multiple keys, this parameter can be
specified as a list reference. Finally, you can delete attributes for a range of
keys by putting a percent sign (C<%>) at the end.

=item values

List of the desired attribute values, section by section. If C<undef>
or an empty string is specified, all values in that section will match. A
generic match can be requested by placing a percent sign (C<%>) at the end.
In that case, all values that match up to and not including the percent sign
will match. You may also specify a regular expression enclosed
in slashes. All values that match the regular expression will be deleted. For
performance reasons, only values have this extra capability.

=item RETURN

Returns a list of tuples for the attributes that were deleted, in the
same form as L</GetAttributes>.

=back

=cut

sub DeleteMatchingAttributes {
    # Get the parameters.
    my ($self, $objectID, $key, @values) = @_;
    # Call the remote server.
    my $result = $self->{proxy}->call('DeleteMatchingAttributes', $objectID, $key, @values);
    Trace("Result type is " . ref($result) . ".") if T(3);
    Trace(scalar(@{$result}) . " matches found.") if T(3) && ref($result) eq 'ARRAY';
    # Unwrap the list and return it.
    my @retVal = @{$result};
    return @retVal;
}

=head3 GetAttributeData

    my %keys = $attrDB->GetAttributeData($type, @list);

Return attribute data for the selected attributes. The attribute
data is a hash mapping each attribute key name to a n-tuple containing the
data type, the description, and the groups. This is the same format expected in
the L</FieldMenu> and L</ControlForm> methods for the list of attributes to display.

=over 4

=item type

Type of attribute criterion: C<name> for attributes whose names begin with the
specified string, or C<group> for attributes in the specified group.

=item list

List containing the names of the groups or keys for the desired attributes.

=item RETURN

Returns a hash mapping each attribute key name to its data type, description, and
parent groups.

=back

=cut

sub GetAttributeData {
    # Get the parameters.
    my ($self, $type, @list) = @_;
    # Call the remote server.
    my $result = $self->{proxy}->call('GetAttributeData', $type, @list);
    Trace("Result type is " . ref($result) . ".") if T(3);
    Trace(scalar(@{$result}) . " matches found.") if T(3) && ref($result) eq 'ARRAY';
    # Unwrap the list as a hash and return it.
    my %retVal = @{$result};
    return %retVal;
}

1;
