#!/usr/bin/perl -w
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


package FullLocation;

    use strict;
    use Tracer;
    use PageBuilder;
    use BasicLocation;

=head1 Full Location Object

=head2 Introduction

A I<full location object> describes a list of basic locations (segments) in a
particular genome. In addition to an array of basic location objects, it contains
a genome ID and a reference to a FIG-like object. The FIG-like object is always
accessed using the variable I<$fig>. This simplifies the process of determining which
FIG methods must be supported in order to make use of this object's features.

The simplest way to create a full location object is by passing in the genome ID
and a location string. The location string contains a list of basic locations
separated by commas. These are converted into location objects and assembled into
the full location. The full location is considered I<bounded> by the first and last
basic locations in the list. This bounded region has a B<Left>, B<Right>, B<Begin>, and
B<Endpoint>, just like a basic location. So, for example, with a location list of

    RED_100_250, RED_275_325, RED_330_430

the B<Left> and B<Begin> locations are C<RED_100>, while the B<Right> and B<EndPoint>
locations are C<RED_430>. Similarly, with the location list

    BLUE_500_450, BLUE_425_300, BLUE_295_200

the B<Left> and B<EndPoint> locations are C<BLUE_200>, while the B<Right> and B<Begin>
locations are C<BLUE_500>.

A location can be converted to a DNA string using the genome ID and data accessible
through the fig-like object. A location also has a I<translation> that represents the
protein sequence produced DNA. The translation can be computed from the DNA or can be
provided by the client. If it is provided by the client, it will automatically be
extended when the boundaries are moved.

Theoretically, a location can contain basic locations on different contigs and
pointing in different directions. In practice, this does not occur; therefore,
the location will be treated as a set of basic locations in a single direction
on a single contig. If this is not the case, problems will arise with some of
the methods.

=cut

#: Constructor FullLocation->new();

=head2 Public Methods

=head3 new

    my $loc = FullLocation->new($fig, $genomeID, $locList, $translation, $code);

Construct a new FullLocation object.

=over 4

=item fig

A fig-like object that can be used to get the DNA and translation information.

=item genomeID

ID of the genome containing the location.

=item locList

List of locations. This can be a reference to a list of location strings, a
comma-delimited list of location strings, or a reference to a list of basic
location objects.

=item translation (optional)

Protein string representing the DNA inside the boundaries of the location.

=item code (optional)

Translation code table to use for translating DNA triplets to proteins. If
none is specified, the standard code table will be used. If a number is
specified, then the appropriate NCBI code table will be requested from
[[FigPm]].

=back

=cut

sub new {
    # Get the parameters.
    my ($class, $fig, $genomeID, $locList, $translation, $code) = @_;
    # If there's no code, default to the standard translation code.
    if (! defined $code) {
        $code = FIG::standard_genetic_code();
    } elsif (! ref $code) {
        $code = FIG::genetic_code($code);
    }
    # Create the $loc object.
    my $retVal = {
                  fig => $fig,
                  genomeID => $genomeID,
                  translation => $translation,
                  code => $code
                 };
    # The tricky part is the location list. Regardless of its incoming format,
    # we must convert it to a list of location objects.
    my $locType = ref $locList;
    if ($locType eq '') {
        # Here we have a comma-delimited list of locations, so we convert it to
        # an array reference.
        $retVal->{locs} = _ParseLocations($retVal, split /\s*,\s*/, $locList);
    } elsif ($locType eq 'ARRAY') {
        # Here we have an array of location objects or strings, which we can
        # pass directly to the parser.
        $retVal->{locs} = _ParseLocations($retVal, @{$locList});
    } else {
        Confess("Invalid location list parameter of type $locType.");
    }
    # Scan the location array to determine the contig and direction. We choose
    # the most popular contig and direction to represent all the locations.
    my %contigs = ();
    my %dirs = ( '+' => 0, '-' => 0 );
    for my $loc (@{$retVal->{locs}}) {
        $dirs{$loc->Dir}++;
        $contigs{$loc->Contig}++
    }
    # Choose the most popular direction and contig.
    $retVal->{dir} = GetBest(%dirs);
    $retVal->{contig} = GetBest(%contigs);
    if ($dirs{$retVal->{dir}} == 0) {
        # Here the location list was empty.
        Confess("Attempt to create a location for genome $genomeID that has no locations.");
    }
    # Bless and return the object.
    bless $retVal, $class;
    return $retVal;
}

=head3 Locs

    my $locObject = $loc->Locs->[$idx];

Return a reference to the array of location objects.

B<NOTE>: Do not update the locations, as it will mess up the translation. Use the
L</Extend> method to change a full location's begin and/or end points.

=cut
#: Return Type $@;
sub Locs {
    return $_[0]->{locs};
}

=head3 Contig

    my $contigID = $loc->Contig();

Return the ID of the base contig for this location list. The base contig is the
contig used for most of the locations in the list.

=cut
#: Return Type $;
sub Contig {
    # Get the parameters.
    my ($self) = @_;
    # Return the result.
    return $self->{contig};
}

=head3 Dir

    my $dir = $loc->Dir;

Return the base direction for this location list. The base direction is the direction
(C<+> or C<->) used for most of the locations in the list.

=cut
#: Return Type $;
sub Dir {
    # Get the parameters.
    my ($self) = @_;
    # Return the result.
    return $self->{dir};
}

=head3 Length

    my $len = $loc->Length

Return the total length of all the locations.

=cut

sub Length {
    # Get the parameters.
    my ($self) = @_;
    # Declare the result variable.
    my $retVal = 0;
    # Loop through the constituent locations.
    for my $loc (@{$self->{locs}}) {
        $retVal += $loc->Length;
    }
    # Return the total.
    return $retVal;
}

=head3 NextPoint

    my $offset = $loc->NextPoint;

Return the location immediately after the end point of the last location.

=cut
#: Return Type $;
sub NextPoint {
    # Get the parameters.
    my ($self) = @_;
    # Get the last location.
    my (undef, undef, $locN) = $self->GetBounds();
    # Return the point after it.
    return $locN->PointOffset($locN->Length);
}

=head3 PrevPoint

    my $offset = $loc->PrevPoint;

Return the location immediately before the begin point of the first location.

=cut
#: Return Type $;
sub PrevPoint {
    # Get the parameters.
    my ($self) = @_;
    # Get the last location.
    my (undef, $loc1, undef) = $self->GetBounds();
    # Return the point after it.
    return $loc1->PointOffset(-1);
}

=head3 Begin

    my $offset = $loc->Begin;

Return the begin point of the first location.

=cut
#: Return Type $;
sub Begin {
    # Get the parameters.
    my ($self) = @_;
    # Get the first location.
    my (undef, $loc1, undef) = $self->GetBounds();
    # Return its begin point.
    return $loc1->Begin;
}

=head3 EndPoint

    my $offset = $loc->EndPoint;

Return the end point of the last location.

=cut
#: Return Type $;
sub EndPoint {
    # Get the parameters.
    my ($self) = @_;
    # Get the last location.
    my (undef, undef, $locN) = $self->GetBounds();
    # Return its end point.
    return $locN->EndPoint;
}

=head3 SeedString

    my $string = $loc->SeedString;

Return a comma-delimited list of this object's basic locations, in SEED format.

=cut
#: Return Type $;
sub SeedString {
    # Get the parameters.
    my ($self) = @_;
    # Map the location list to SEED strings.
    my @seeds = map { $_->SeedString } @{$self->{locs}};
    # Return the result.
    return join ", ", @seeds;
}

=head3 Adjusted

    my $offset = $loc->Adjusted($oldOffset, $distance);

Adjust the specified offset by the specified distance in the direction of this
location. If this is a forward location, the distance is added; if it is a backward
location, the distance is subtracted.

=over 4

=item oldOffset

Offset to adjust.

=item distance

Distance by which to adjust the offset. This value can be negative.

=item RETURN

Returns a new offset formed by moving the specified distance from the original offset
in this location's direction.

=back

=cut
#: Return Type $;
sub Adjusted {
    # Get the parameters.
    my ($self, $oldOffset, $distance) = @_;
    # Do the adjustment.
    return $oldOffset + ($self->Dir eq '+' ? $distance : -$distance);
}

=head3 GetBest

    my $bestKey = FullLocation::GetBest(%hash);

Return the key of the hash element with the highest positive numeric value.

=over 4

=item hash

A hash mapping keys to numbers.

=item RETURN

Returns the key having the highest value. If the hash is empty, or has no non-negative
values, returns C<undef>.

=back

=cut
#: Return Type $;
sub GetBest {
    # Get the parameters.
    my (%hash) = @_;
    # Declare the return variable and initialize the best-count.
    my ($retVal, $best) = (undef, 0);
    # Search the hash.
    for my $key (keys %hash) {
        my $value = $hash{$key};
        if ($value >= $best) {
            $retVal = $key;
            $best = $value;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 DNA

    my $dnaString = $loc->DNA;

Return the complete DNA string for this location.

=cut
#: Return Type $;
sub DNA {
    # Get the parameters.
    my ($self) = @_;
    my $fig = $self->{fig};
    # Use the FIG object to extract the DNA.
    my $retVal = $fig->dna_seq($self->{genomeID}, $self->Contig, map { $_->SeedString } @{$self->Locs});
    # Return the result.
    return $retVal;
}

=head3 Codon

    my $codon = $loc->Codon($point);

Return the DNA codon at the specified point on this location's contig in this
location's direction.

=over 4

=item point

Offset into the contig of the codon.

=item RETURN

Returns a three-letter DNA codon from the specified point.

=back

=cut
#: Return Type $;
sub Codon {
    # Get the parameters.
    my ($self, $point) = @_;
    # Get the FIG object.
    my $fig = $self->{fig};
    # Compute the codon location.
    my $loc = $self->Contig . "_" . $point . "_" . $self->Adjusted($point,2);
    # Return the DNA.
    return $fig->dna_seq($self->{genomeID}, $loc);
}

=head3 Translation

    my $proteinString = $loc->Translation($code, $fixStart);

Return the protein translation of this location's DNA. The first time a
translation is requested, it will be cached in the object, and returned
unmodified. It is also possible that a translation specified in the
constructor exists, in which case it will be returned. Thus, the
I<$code> and I<$fixStart> parameters only matter on the first call.

=over 4

=item code (optional)

Translation code table, in the form of a hash mapping DNA triples to protein
letters. If omitted, the standard translation table in C<FIG.pm> will be
used.

=item fixStart (optional)

TRUE if the first DNA triple should be given special handling, else FALSE.
If TRUE, then a value of C<TTG> or C<GTG> in the first position will be
translated to C<M> instead of the value specified in the translation code.

=item RETURN

Returns the protein translation for this location.

=back

=cut
#: Return Type $;
sub Translation {
    # Get the parameters.
    my ($self, $code, $fixStart) = @_;
    # Declare the return variable.
    my $retVal;
    # Check for a cached translation.
    if ($self->{translation}) {
        # Return the cahced translation.
        $retVal = $self->{translation};
    } else {
        # Check for a translation code.
        if (! defined $code) {
	    #
	    # Don't default to the standard one, you very carefully
	    # created one in the constructor. Why not use it.
	    #
            $code = $self->{code};
        }
        # Here we have to do some work. Extract our DNA.
        my $dna = $self->DNA;
        # Translate it.
        $retVal = FIG::translate($dna, $code, $fixStart);
        # Chop off the stop codon.
        $retVal =~ s/\*$//;
	
        # Cache the translation and the code table.
	#
	# This probably should be cached based on the code table used. It maybe
	# shouldn't even be cached at all.
        $self->{translation} = $retVal;

	#
	# Why would you cache this; it was probably passed in as an exceptional case.
        #$self->{code} = $code;
    }
    # Return the result.
    return $retVal;
}

=head3 ConstrainPoint

    my $constrainedPoint = $loc->ConstrainPoint($point);

Change a point location value so that it fits inside the base contig. If the point
location is less than 1, it will be set to 1. If it's greater than the length of
the contig, it will be set to the length of the contig.

=over 4

=item point

Location to be constrained.

=item RETURN

Returns the value of the nearest location that is on the base contig of this location.

=back

=cut
#: Return Type $;
sub ConstrainPoint {
    # Get the parameters.
    my ($self, $point) = @_;
    # Declare the return variable.
    my $retVal;
    # Check for a value less than 1.
    if ($point < 1) {
        $retVal = 1;
    } else {
        # Check for a value off the end of the contig.
        my $fig = $self->{fig};
        my $contigEnd = $fig->contig_ln($self->{genomeID}, $self->Contig);
        if ($point > $contigEnd) {
            # Bring the point back onto the contig.
            $retVal = $contigEnd;
        } else {
            # Return the incoming point unmodified.
            $retVal = $point;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 ExtremeCodon

    my $loc->ExtremeCodon($dir);

Return the most extreme codon in the specified direction. This is not always the most
extreme location, since the distance to the appropriate edge of the location must be
a multiple of 3.

=over 4

=item dir

C<first> to get the codon moving away from the beginning of the location, C<last>
to get the codon moving away from the end of the location.

=item RETURN

Returns the edge location of the desired codon. If we are going toward the left, this
is the left point in the codon; if we are going toward the right, this is the right
point in the codon.

=back

=cut

sub ExtremeCodon {
    # Get the parameters.
    my ($self, $dir) = @_;
    my $fig = $self->{fig};
    # The first task is to determine the starting point and direction for the
    # search. We start by converting the direction to the same format as the
    # location direction.
    my $parity = ($dir eq 'first' ? '-' : '+');
    # Get the contig length.
    my $contig_len = $fig->contig_ln($self->{genomeID}, $self->Contig);
    # If we're moving in the opposite direction from the location, we're going to
    # go toward the beginning of the contig; otherwise, we're going toward the
    # end.
    my ($multiplier, $endPoint);
    if ($parity ne $self->Dir) {
        ($multiplier, $endPoint) = (-3, 1);
    } else {
        ($multiplier, $endPoint) = (3, $contig_len);
    }
    # Now we need the start point, which is determined by direction of this method.
    my $beginPoint = ($parity eq '-' ? $self->Begin : $self->EndPoint);
    # Compute the number of positions to move and add it to the begin point.
    my $retVal = int(($endPoint - $beginPoint) / $multiplier) * $multiplier +
                 $beginPoint;
    # Return the codon found.
    return $retVal;
}


=head3 Extend

    my  = $loc->Extend($newBegin, $newEnd, $trimFlag);

Extend this gene to a new begin point and a new end point. If a translation exists,
it will be updated to match the new locations. The I<$trimFlag> indicates whether
or not it is permissible to shrink the location at either end. If an attempt is made
to shrink and I<$trimFlag> is not specified, then a fatal error will occur.

=over 4

=item newBegin

Proposed new beginning offset. If undefined, the begin location will not be changed.

=item newEnd

Proposed new ending offset. If undefined, the ending location will not be changed.

=item trimFlag

If TRUE, the begin and end offsets can shrink the location.

=back

=cut
#: Return Type ;
sub Extend {
    # Get the parameters.
    my ($self, $newBegin, $newEnd, $trimFlag) = @_;
    my $fig = $self->{fig};
    # Get our boundaries.
    my ($boundLoc, $loc1, $locN) = $self->GetBounds;
    # Get the current length of the start location.
    my $len = $boundLoc->Length;
    # Extend the beginning of the bounds.
    if (defined $newBegin) {
        $boundLoc->SetBegin($newBegin);
        my $excess = $boundLoc->Length - $len;
        # Insure this is a real extension.
        if ($excess < 0) {
            # Find out if we can trim.
            if (! $trimFlag) {
                # We can't trim, so it's an error.
                Confess("Invalid begin location $newBegin for location " . $boundLoc->String . ".");
            } elsif ($self->{translation}) {
                # We can trim, and a translation exists, so we lop some characters off
                # the front. Note that we divide by -3 because if we're here, the excess
                # is automatically negative.
                my $proteinExcess = int($excess / -3);
                $self->{translation} = substr $self->{translation}, $proteinExcess;
            }
        } elsif ($excess > 0 && $self->{translation}) {
            # Here we have new stuff to translate. Get its location.
            my $excessLoc = BasicLocation->new($loc1->Contig, $newBegin, $loc1->Dir, $excess);
            # Extract its DNA and translate it.
            my $newDNA = $fig->dna_seq($self->{genomeID}, $excessLoc->SeedString);
            my $newTran = FIG::translate($newDNA, $self->{code});
            # Prefix the new translation to the old one.
            $self->{translation} = $newTran . $self->{translation};
        }
        # We successfully updated the translation (if necessary), so we adjust the
        # start of the first location in the full location's list.
        $loc1->SetBegin($newBegin);
        # Get the new current length of the bounds location.
        $len = $boundLoc->Length;
    }
    # Extend the ending.
    if (defined $newEnd) {
        $boundLoc->SetEnd($newEnd);
        my $excess = $boundLoc->Length - $len;
        # Insure this is a real extension.
        if ($excess < 0) {
            # Find out if we can trim.
            if (! $trimFlag) {
                # We can't trim, so it's an error.
                Confess("Invalid end location $newEnd for location " . $boundLoc->String . ".");
            } elsif ($self->{translation}) {
                # We can trim, and a translation exists, so we lop off some proteins at
                # the end. Note that we divide by 3 and the excess is negative, so the
                # result will be negative. We use it as a negative length in the substr
                # expression to trim end characters.
                my $negativeProteinExcess = int($excess / 3);
                $self->{translation} = substr $self->{translation}, 0, $negativeProteinExcess;
            }
        } elsif ($excess > 0 && $self->{translation}) {
            # Here we have new stuff to translate. Get its location.
            my $excessLoc = BasicLocation->new($locN->Contig, $boundLoc->PointOffset($len), $locN->Dir,
                                               $excess);
            # Extract its DNA and translate it.
            my $newDNA = $fig->dna_seq($self->{genomeID}, $excessLoc->SeedString);
            my $newTran = FIG::translate($newDNA, $self->{code});
            # Append the new translation to the old one.
            $self->{translation} .= $newTran;
        }
        # Here we sucessfully updated the translation and the update is legal, so
        # we can modify the end of the last location.
        $locN->SetEnd($newEnd);
    }
    if ($self->{translation}) {
        # Chop the stop codon off the end of the translation.
        $self->{translation} =~ s/\*$//;
    }
}

=head3 ConstrainCodon

    my $point = $loc->ConstrainCodon($point, $codonPoint);

Constrain the specified point so that it is inside the bounds of this
location's contig and its distance to the specified codon point is a
multiple of 3.

=over 4

=item point

Point index (relative to the contig) to be constrained.

=item codonPoint

Index (relative to the contig) of a point that is the start of a codon.

=item RETURN

Returns the constrained index.

=back

=cut

sub ConstrainCodon {
    # Get the parameters.
    my ($self, $point, $codonPoint) = @_;
    # Declare the return variable.
    my $retVal = $point;
    # Check for too far left.
    if ($retVal < 1) {
        $retVal = ($codonPoint - 1) % 3 + 1;
    } else {
        # Check for too far right.
        my $contigLen = $self->{fig}->contig_ln($self->{genomeID}, $self->Contig);
        if ($retVal > $contigLen) {
            $retVal = $contigLen - ($contigLen - $codonPoint) % 3;
        }
    }
    # Return the result.
    return $retVal;
}


=head3 GetBounds

    my ($boundLoc, $loc1, $locN) = $loc->GetBounds;

Analyze this location and return information about its boundaries. This includes
a location for the bounds, the first location, and the last location. The
bounds essentially define the location as it would be if it were all on a single
contig in the same direction and had no gaps. The first location is the location
object containing the begin point of the bounds, and the last location is the
location object containing the end point of the bounds.

=cut

sub GetBounds {
    # Get the parameters.
    my ($self) = @_;
    # Declare the return variables.
    my ($boundLoc, $loc1, $locN);
    # Get a reference to the location list.
    my $locList = $self->Locs;
    # The most common case is a singleton location list. We handle that first.
    if (@{$locList} == 1) {
        my $bloc = $locList->[0];
        $boundLoc = BasicLocation->new($bloc);
        $loc1 = $bloc;
        $locN = $bloc;
    } else {
        # Here we have a multiple-location list. We search for the leftmost left
        # and rightmost right on the base contig. To do that, we first extract
        # all the eligible locations.
        my $baseContig = $self->Contig;
        my @baseLocs = grep { $_->Contig eq $baseContig } @{$locList};
        # Next we prime the loop with a location popped off the list.
        my $loc0 = pop @baseLocs;
        my ($leftLoc, $rightLoc) = ($loc0, $loc0);
        # Search for the leftmost and rightmost locations.
        for my $loci (@baseLocs) {
            if ($loci->Left < $leftLoc->Left) {
                $leftLoc = $loci;
            }
            if ($loci->Right > $rightLoc->Right) {
                $rightLoc = $loci;
            }
        }
        # Now we have enough information to build the bounding location. First,
        # we get the length.
        my $len = $rightLoc->Right - $leftLoc->Left + 1;
        # Next, we arrange the left and right locations according to the direction.
        if ($self->Dir eq '+') {
            ($loc1, $locN) = ($leftLoc, $rightLoc);
        } else {
            ($loc1, $locN) = ($rightLoc, $leftLoc);
        }
        # Finally, we create the bounding location.
        $boundLoc = BasicLocation->new($baseContig, $loc1->Begin, $self->Dir, $len);
    }
    # Return the results.
    return ($boundLoc, $loc1, $locN);
}

=head2 Codon Search Methods

=head3 UpstreamSearch

    my $loc = $floc->UpstreamSearch($pattern, $limit);

Search upstream from this location for a codon as identified by the
specified pattern, stopping at the end of the contig or when the
specified limit is reached.

=over 4

=item pattern

Codon pattern to search for, expressed as a bar-delimited list of base triplets.
For example, C<taa|tag|tga> would search for a stop codon.

=item limit (optional)

Maximum number of base pairs to search. Must be a multiple of 3.

=item RETURN

Returns a [[BasicLocationPm]] object for the codon found.

=back

=cut

sub UpstreamSearch {
    # Get the parameters.
    my ($self, $pattern, $limit) = @_;
    # Get the FIG object.
    my $fig = $self->{fig};
    # Declare the return variable.
    my $retVal;
    # Locate the starting and ending positions for the search. The search
    # ends immediately upstream of our begin point.
    my $end = $self->PrevPoint;
    my $start = $end - ($self->Dir . 1) * ($limit - 1);
    # Constrain these values to the inside of the contig.
    $start = $self->ConstrainCodon($start, $self->Begin);
    $end = $self->ConstrainCodon($end, $self->Begin + ($self->Dir . 2));
    # Get the DNA to search. Note we convert it automatically to lower case.
    my $dna = lc $fig->dna_seq($self->{genomeID}, $self->Contig . "_${start}_${end}");
    # Insure the pattern is also lower-case.
    $pattern = lc $pattern;
    # Get the location of the last codon in the dna sequence.
    my $i1 = length($dna) - 3;
    Trace("$i1 base pairs in search.") if T(4);
    Trace("Upsearch DNA translation\n" . FIG::translate($dna, $self->{code})) if T(4);
    for (my $i = $i1; $i >= 0 && ! defined($retVal); $i -= 3) {
        # Check for a match.
        if (substr($dna, $i, 3) =~ /$pattern/) {
            # Compute the actual return value. This will also stop the loop.
            $retVal = BasicLocation->new($self->Contig, $i * ($self->Dir . 1) + $start, $self->Dir, 3);
        }
    }
    # Return the result.
    return $retVal;
}


=head3 DownstreamSearch

    my $loc = $floc->DownstreamSearch($pattern, $limit);

Search downstream from this location for a codon as identified by the
specified pattern, stopping at the end of the contig or when the
specified limit is reached.

=over 4

=item pattern

Codon pattern to search for, expressed as a bar-delimited list of base triplets.
For example, C<taa|tag|tga> would search for a stop codon.

=item limit (optional)

Maximum number of base pairs to search. Must be a multiple of 3.

=item RETURN

Returns a [[BasicLocationPm]] object for the codon found.

=back

=cut

sub DownstreamSearch {
    # Get the parameters.
    my ($self, $pattern, $limit) = @_;
    # Get the FIG object.
    my $fig = $self->{fig};
    # Declare the return variable.
    my $retVal;
    # Locate the starting and ending positions for the search. The search starts
    # immediately downstream of our end point.
    my $start = $self->NextPoint;
    my $end = $start + ($self->Dir . 1) * ($limit - 1);
    # Do a down search to find the codon.
    $retVal = $self->DownSearch($pattern, $start, $end);
    # Return the result.
    return $retVal;
}


=head3 InsideSearch

    my $loc = $floc->InsideSearch($pattern);

Search inside this location for a codon as identified by the
specified pattern.

=over 4

=item pattern

Codon pattern to search for, expressed as a bar-delimited list of base triplets.
For example, C<taa|tag|tga> would search for a stop codon.

=item RETURN

Returns a [[BasicLocationPm]] object for the codon found.

=back

=cut

sub InsideSearch {
    # Get the parameters.
    my ($self, $pattern, $limit) = @_;
    # Get the FIG object.
    my $fig = $self->{fig};
    # Declare the return variable.
    my $retVal;
    # Locate the starting and ending positions for the search. The search
    # is entirely inside the location.
    my $start = $self->Begin;
    my $end = $self->EndPoint;
    # Do a down search to find the codon.
    $retVal = $self->DownSearch($pattern, $start, $end);
    # Return the result.
    return $retVal;
}

=head3 DownSearch

    my $loc = $floc->DownSearch($pattern, $start, $end);

Search parallel to this location for a codon as identified by the
specified pattern, starting at the specified start point and stopping
at the specified end point.

=over 4

=item pattern

Codon pattern to search for, expressed as a bar-delimited list of base triplets.
For example, C<taa|tag|tga> would search for a stop codon.

=item $start

Starting position for search.

=item $end

Ending position for search.

=item RETURN

Returns a [[BasicLocationPm]] object for the codon found.

=back

=cut

sub DownSearch {
    # Get the parameters.
    my ($self, $pattern, $start, $end) = @_;
    # Get the FIG object.
    my $fig = $self->{fig};
    # Declare the return variable.
    my $retVal;
    # Insure we're inside the contig.
    my $realStart = $self->ConstrainCodon($start, $start);
    my $realEnd = $self->ConstrainCodon($end, $start + ($self->Dir . 2));
    # Get the DNA to search. Note we convert it automatically to lower case.
    my $dna = lc $fig->dna_seq($self->{genomeID}, $self->Contig . "_${realStart}_${realEnd}");
    # Insure the pattern is also lower-case.
    $pattern = lc $pattern;
    # Get the length of the dna sequence.
    my $i1 = length($dna);
    Trace("$i1 base pairs in search.") if T(4);
    Trace("Downsearch DNA translation\n" . FIG::translate($dna, $self->{code})) if T(4);
    for (my $i = 0; $i < $i1 && ! defined($retVal); $i += 3) {
        # Check for a match.
        if (substr($dna, $i, 3) =~ /$pattern/) {
            # Compute the actual return value. This will also stop the loop.
            $retVal = BasicLocation->new($self->Contig, $i * ($self->Dir . 1) + $realStart, $self->Dir, 3);
        }
    }
    # Return the result.
    return $retVal;
}

=head3 PickGeneBoundaries

    my $rc = $floc->PickGeneBoundaries(-stop => $stopPattern,
        -start => $startPattern, -limit => $limit);

Update this location so that it has a valid start and stop. The basic
algorithm used is to search upstream and downstream for stop codons, then
search between the stop codons for the first start codon.

=over 4

=item stop (optional)

Search pattern for the stop codon, encoded as a bar-delimited list of DNA triplets
(e.g. C<tta|ata|tag>). If omitted, the default is to look for stop codons in the
attached genetic code table.

=item start (optional)

Search pattern for the start codon, encoded as a bar-delimited list of DNA triplets.
If omitted, the default is C<atg|gtg|ttg>.

=item limit (optional)

If a number, then the maximum distance to search when attempting to extend the
location. If a [[BasicLocationPm]] object or a location string, then none of the searches will
go outside the region spanned by the location. If omitted, the default is a scalar value of C<9000>.

=item RETURN

Returns TRUE if successful, FALSE if the process fails. 

=back

=cut

sub PickGeneBoundaries {
    # Get the parameters.
    my ($self, %parms) = @_;
    # Declare the return variable.
    my $retVal;
    # Get the stop pattern.
    my $stopPattern = $parms{-stop} || lc join("|", grep { $self->{code}->{$_} eq '*' } keys %{$self->{code}});
    # Get the start pattern.
    my $startPattern = $parms{-start} || 'atg|gtg|ttg';
    # Compute the limits. We have an upstream limit and a downstream limit.
    my $limit = $parms{-limit};
    my ($upLimit, $downLimit);
    if (! defined $limit) {
        $upLimit = 9000;
        $downLimit = 9000;
    } elsif ($limit =~ /^\d+$/) {
        $upLimit = $limit;
        $downLimit = $limit;
    } else {
        my $limitLoc = (! ref $limit ? new BasicLocation($limit) : $limit);
        $upLimit = abs($limitLoc->Begin - $self->Begin);
        $downLimit = abs($limitLoc->EndPoint - $self->EndPoint);
    }
    # Insure the limits are on codon boundaries.
    $upLimit -= $upLimit % 3;
    $downLimit -= $downLimit % 3;
    # Save the current boundaries.
    my $oldBegin = $self->Begin;
    my $oldEnd = $self->EndPoint;
    # Get the contig boundaries.
    my $contigBegin = $self->ExtremeCodon('first');
    my $contigEnd = $self->ExtremeCodon('last');
    # Get the distance to each extreme.
    my $upExtreme = abs($contigBegin - $oldBegin);
    my $downExtreme = abs($contigEnd - $oldEnd);
    Trace("Up: limit = $upLimit, extreme = $upExtreme. Down: limit = $downLimit, extreme = $downExtreme.") if T(4);
    # Search upstream for a stop.
    my $upLoc = $self->UpstreamSearch($stopPattern, $upLimit);
    # Check to see if we found one.
    my $newBegin;
    if (defined $upLoc) {
        # Here we found a stop, so the new beginning is after the codon found.
        my $newOrfBegin = $upLoc->PointOffset(3);
        $self->Extend($newOrfBegin);
        # Now look for a start between the stop codon found and the end of the location.
        my $startLoc = $self->InsideSearch($startPattern);
        if (! defined $startLoc) {
            Trace("New start not found for " . $self->SeedString() . ".") if T(3);
        } else {
            # We found the start, so it becomes our new beginning.
            $newBegin = $startLoc->Begin;
        }
    } elsif ($upExtreme <= $upLimit) {
        # Here we fell off the end of the contig while searching, so this is ok.
        Trace("Contig beginning used for stop codon after upstream search.") if T(3);
        $newBegin = $contigBegin;
    } else {
        Trace("Upstream stop not found for " . $self->SeedString() . ".") if T(3);
    }
    # Only proceed if we have a new beginning.
    if (defined $newBegin) {
        # Search downstream for a stop codon.
        my $stopLoc = $self->DownstreamSearch($stopPattern, $downLimit);
        # Check the result.
        my $endPoint;
        if (defined $stopLoc) {
            # Here we found the end point.
            $endPoint = $stopLoc->EndPoint;
        } elsif ($downLimit >= $downExtreme) {
            # Here we fell off the end, so we use the end.
            $endPoint = $contigEnd;
            Trace("Contig end used for stop codon after downstream search.") if T(3);
        } else {
            Trace("Downsream stop not found for " . $self->SeedString() . ".") if T(3);
        }
        if (defined $endPoint) {
            # Extend the location to the start and stop found. This will almost certainly
            # require trimming in the start-codon direction.
            $self->Extend($newBegin, $endPoint, 'trim');
            # Denote we've succeeded.
            $retVal = 1;
        }
    }
    # If we failed, restore the location.
    if (! $retVal) {
        $self->Extend($oldBegin, $oldEnd, 'trim');
    }
    # Return the success indicator.
    return $retVal;
}


=head2 Internal Methods

=head3 ParseLocations

Parse an array into a list of basic location objects. The array can contain
basic location objects or location strings. The first parameter must be the
full location object being constructed.

This is a static method.

=cut

sub _ParseLocations {
    # Get the location list.
    my ($parent, @locs) = @_;
    # Create the return array.
    my @retVal = ();
    # Create a location index counter.
    my $idx = 0;
    # Loop through the locations.
    for my $loc (@locs) {
        # Create a variable to hold the location object created.
        my $locObject;
        # Check to see if this is a string or a location object.
        if (ref $loc eq '') {
            # It's a string, so parse it into a location object.
            $locObject = BasicLocation->new($loc, $parent, $idx);
        } elsif (UNIVERSAL::isa($loc, "BasicLocation")) {
            # It's a location object, so copy it and set the parent and
            # index.
            $locObject = BasicLocation->new($loc);
            $locObject->Attach($parent, $idx);
        } else {
            # Here we have an error.
            my $type = ref $loc;
            Confess("Invalid location object of type $type found at index $idx.");
        }
        # Add the location to the list.
        push @retVal, $locObject;
        $idx++;
    }
    # Return a reference to the location list.
    return \@retVal;
}

=head3 FindPattern

Locate the index of a specified pattern in a DNA string.

This is a static method.

=over 4

=item dna

DNA string to search.

=item pattern

Pattern for which to search (see L</Search>).

=item RETURN

Returns the index of the specified pattern, or C<undef> if the pattern is not found.

=back

=cut
#: Return Type $;
sub _FindPattern {
    # Get the parameters.
    my ($dna, $pattern) = @_;
    # Declare the return variable.
    my $retVal;
    # Insure the pattern is lower case.
    my $realPattern = lc $pattern;
    # Start at the beginning of the string. We will chop stuff off the string
    # as we search through it. This smount chopped must then be added to the offset
    # of the found string in order to get the return value.
    my $pos = 0;
    # Do the search. We search for stop codons, and stop at the first one
    # that lies on a codon boundary.
    while (!defined($retVal) && $dna =~ m/$realPattern/g) {
        # We have a match. Get its location. Note that the "pos" function returns
        # the point where the search left off, so we have to back off a bit to
        # get the starting point of the codon.
        my $newPos = (pos $dna) - 3;
        # See if this is on a codon boundary.
        my $mod = $newPos % 3;
        if ($mod == 0) {
            # It is, so return the offset from the original start of the string.
            $retVal = $pos + $newPos;
        } else {
            # Here we need to keep searching. First, however, we move the search
            # position forward to the next codon boundary. This avoids useless checking
            # of the next byte or two and insures that we don't find the same value
            # again.
            $newPos += 3 - $mod;
            $pos += $newPos;
            $dna = substr($dna, $newPos);
        }
    }
    # Return the result.
    return $retVal;
}


1;

