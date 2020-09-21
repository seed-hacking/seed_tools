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

package SeedDas;
use Data::Dumper;
#$Data::Dumper::Deparse = 1;

use strict;

use Text::ParseWords;

use FIG_Config;
use SproutSearch;

#
# DAS Support
#

sub new
{
    my($class, $fig, $data_dir, $url, $dsn) = @_;

    my $genome;

    if ($dsn =~ /^(\d+)[_.](\d+)$/)
    {
        $genome = "$1.$2";
    }
    else
    {
        warn "Genome not found for dsn $dsn\n";
    }

    my $self = {
        fig => $fig,
        dir => $data_dir,
        url => $url,
        dsn => $dsn,
        genome => $genome,
        peg_index => {},
        cache_peg_name => {},
        cache_peg_func => {},
    };

    bless $self, $class;
    
    $self->load_track_definitions();
    $self->load_track_bindings();

    # warn Dumper($self);

    return $self;
}


=pod

=head1 load_tracks

Load the track definitions (binding between track name and code that
implements the track), and the track bindings (binding between a
genome and the tracks that are available on that genome).


Track definitions are found in FIG/Data/Global/DAS/Track.defs.

Track bindings are found in FIG/Data/Global/DAS/Default.bindings (bindings
for all genomes), with additional per-genome bindings in
DAS/<GenomeID>.bindings.
    
After loading a binding, we perform variable expansion on the fields in the
track definition that maps to this binding, and store a copy in the binding.

=cut

sub load_track_definitions
{
    my($self) = @_;

    my $trackdefs = "$self->{dir}/Track.defs";
    my $tfh;

    if (!open($tfh, "<$trackdefs"))
    {
        warn "Could not open $trackdefs: $!\n";
        return undef;
    }

    local $/ = "//\n";

    while (<$tfh>)
    {
        chomp;

        next if $_ eq "";

        my($name, @defs) = split(/\n/, $_);

        my $track = {};

        for my $defline (@defs)
        {
            my($k, $v) = split(/\t/, $defline, 2);
            $track->{$k} = $v;
        }

        $self->{tracks}->{$name} = $track;
        push(@{$self->{track_list}}, $track);

        #
        # map perl binding to a code ref.
        #

        if (my $func = $track->{perl_binding})
        {
            my $ref = \&$func;

            if (defined(&$ref))
            {
                $track->{perl_binding_ref} = $ref;
            }
            else
            {
                warn "No function found for perl binding $func\n";
            }
        }
    }
}
    
sub load_track_bindings
{
    my($self) = @_;

    #
    # Load defaults.
    #

    my $bindings = "$self->{dir}/Default.bindings";
    my $bfh;

    if (!open($bfh, "<$bindings"))
    {
        warn "Could not open $bindings: $!\n";
        return undef;
    }

    local $/ = "//\n";

    while (<$bfh>)
    {
        chomp;

        next if $_ eq "";

        my($name, @rest) = split(/\n/, $_);

        #
        # Read @rest into a binding hash.
        #
        # The associated track is binding_hash{track} if present,
        # otherwise it is assumed to be the same as the name of the binding.
        #

        my $binding = {};
        for my $ent (@rest)
        {
            my($k, $v) = split(/\t/, $ent, 2);
            $binding->{$k} = $v;
        }

        my $track_name = $binding->{track};
        if (!$track_name)
        {
            $track_name = $name;
        }
        
        $self->{bindings}->{$name} = [$track_name, undef, $name, $binding];
    }

    close($bfh);

    #
    # Now load the genome-specific binding information, if any.
    #

    if (open($bfh, "$self->{dir}/$self->{genome}.bindings"))
    {
        while (<$bfh>)
        {
            chomp;
            
            next if $_ eq "";
            
            my($name, @rest) = split(/\n/, $_);
            
            #
            # Read @rest into a binding hash.
            #
            # The associated track is binding_hash{track} if present,
            # otherwise it is assumed to be the same as the name of the binding.
            #
            
            my $binding_info = $self->{bindings}->{$name};
            my $binding;
            
            #
            # See if this is a new binding. If it is, create a new hash.
            # If not, use the existing one.
            #
            if (!defined($binding_info))
            {
                $binding = {};
            }
            else
            {
                $binding = $binding_info->[3];
            }
            
            for my $ent (@rest)
            {
                my($k, $v) = split(/\t/, $ent, 2);
                $binding->{$k} = $v;
            }
            
            #
            # See if we are setting or overriding the track.
            #
            
            my $track_name;
            if (!defined($binding_info))
            {
                $track_name = $binding->{track};
                if (!$track_name)
                {
                    $track_name = $name;
                }
            }
            else
            {
                $track_name = $binding_info->[0];
                if (defined($binding->{track}))
                {
                    $track_name = $binding->{track};
                }
            }
            
            if (!$track_name)
            {
                warn "Track not found for binding '$name'\n";
                next;
            }
            
            $self->{bindings}->{$name} = [$track_name, undef, $name, $binding];
        }
        close($bfh);
    }

    #
    # We have now loaded all the binding information.
    #
    # Walk them, filling out the tracks based on the bindings.
    #

    while (my($bname, $binding_info) = each(%{$self->{bindings}}))
    {
        my($track_name, undef, $name, $binding) = @$binding_info;

        my $orig_track = $self->{tracks}->{$track_name};

        my $track = {};

        for my $k (keys(%$orig_track))
        {
            #
            # Do magic variable replacement on values.
            #

            my $v = $orig_track->{$k};
            $v =~ s,\$\{(\w+)\},defined($binding->{$1}) ? $binding->{$1} : $1,eg;
            $track->{$k} = $v;
        }
        $binding_info->[1] = $track;

        #
        # If we have any perl function args, split them up here using shellwords.
        #
        if (defined($track->{perl_binding_args}))
        {
            $track->{perl_binding_args} = [shellwords($track->{perl_binding_args})];
        }
        else
        {
            $track->{perl_binding_args} = [];
        }

        #
        # The track type is updated to be the specific track type we've
        # determined here.
        #

        $track->{type} = $name;
    }
}

sub types_header
{
    my($self) = @_;

    return <<END;
<?xml version="1.0" standalone="yes"?>
<!DOCTYPE DASTYPES SYSTEM "http://www.biodas.org/dtd/dastypes.dtd">
<DASTYPES>
<GFF version="1.2" summary="yes" href="$self->{url}">
<SEGMENT>
END
    ;

}

sub types_footer
{
    my($self) = @_;

    return <<END;
</SEGMENT>
</GFF>
</DASTYPES>
END
}

=pod

=head1 all_types

Return XML for the list of all types we support.

This comes from the list of bindings for the current genome (dsn).

=cut

sub all_types
{
    my($self) = @_;

    my $out = "";

    #
    # First entry in the bindings list entry is the track.
    #

    for my $ent (values(%{$self->{bindings}}))
    {
        my($track_name, $track, $name, $other_info) = @$ent;
        $out .= $self->get_track_xml($track);
    }
    return $out;
}

sub get_track_xml
{
    my($self, $track) = @_;

    return sprintf(qq(<TYPE id="%s" category="%s" method="%s" source="%s"/>\n),
                   $track->{type},
                   $track->{category},
                   $track->{method},
                   $track->{source});
}


#
# Features.
#

sub features_header
{
    my($self) = @_;

    return <<END
<?xml version="1.0" standalone="yes"?>
<!DOCTYPE DASGFF SYSTEM "http://www.biodas.org/dtd/dasgff.dtd">
<DASGFF>
<GFF version="1.01" href="$self->{url}">
END

}

sub features_footer
{
    my($self) = @_;
    return <<END
</GFF>
</DASGFF>
END

}

sub features
{
    my($self, $segments, $types, $feature_ids)  = @_;

    #
    # Determine if this is a normal get-the-features-from-a-segment
    # call, or if it was likely performed as a result of a gbrowse search.
    #
    # We detect the latter by noticing if @$segments is empty and @$feature_ids is not.
    #

    if (@$segments == () and @$feature_ids > 0)
    {
        return $self->features_search($types, $feature_ids);
    }


    my %types;

    #
    # Build a hash of the types we're looking for.
    #

    map { $types{$_}++ } @$types;

    my $fig = $self->{fig};

    #warn "features: segments=", Dumper($segments);
    #warn "          types=", Dumper($types);
    #warn "          feature_ids=", Dumper($feature_ids);

    #
    # Determine what tracks we need to extract by matching against
    # the list of types requested. Also determine, for the tracks
    # we are keeping, if we need to map segments to features.
    #

    my $need_features = 0;
    
    my @bindings = values(%{$self->{bindings}});
    my @all_tracks = map { $_->[1] } @bindings;
    my @tracks = ();

    # warn "Types: ",  Dumper(\%types);

    #
    # We keep track of all the tracks of method Related so that
    # they can be processed first. This is a little ugly, but
    # there is some cross-track processing that needs to be supported.
    #

    my @related_tracks;
    for my $track (@all_tracks)
    {
        # warn "Testing track ", Dumper($track);
        if ($types{$track->{type}})
        {
            # warn "Adding track $track->{type}\n";

            if ($track->{method} eq "Related")
            {
                push(@related_tracks, $track);
            }
            else
            {
                push(@tracks, $track);
            }
            $need_features++ if $track->{need_pegs};
        }
    }

    #
    # Stuff the related tracks onto the front of the track list.
    #

    unshift(@tracks, @related_tracks);

    #
    # Walk the segments. For each segment, we emit the segment header,
    # and determine if we have any active tracks. If we do,
    # invoke that code to generate features. Then emit the segment footer.
    #
    # At the start of each segment we determine if we have any peg mapping to
    # determine, and do it once.
    #

    # warn "have segs ",  Dumper($segments);
    my $output = [];
    for my $segment (@$segments)
    {
        #
        # Each segment is a list [$reference, $class, $start, $end]
        #
        # $reference is the contig.
        # $class is the class of this segment.
        # $start and $end are the starting and ending points on the contig.
        #

        my($seg_id, $class, $start, $end) = @$segment;

        my $contig = $seg_id;
        $contig =~ s/--/:/g;

        #
        # Determine the start and stop if it was not passed in
        # (an implicit reference to the entire contig).
        #
        if (!defined($start))
        {
            $start = 1;
            $end = $fig->contig_ln($self->{genome}, $contig);
            
            $segment->[2] = $start;
            $segment->[3] = $end;
        }

        push(@$output, qq(<SEGMENT id="$seg_id" start="$start" stop="$end" version="1.0">));
        
        my($pegs, $peg_start, $peg_end);

        if ($need_features)
        {
            ($pegs, $peg_start, $peg_end) = $self->{fig}->genes_in_region($self->{genome},
                                                                          $contig,
                                                                          $start,
                                                                          $end);

            # warn "Checking segment. start=$start end=$end contig=$contig\n";
            # warn "   pegs $pegs $peg_start $peg_end\n";

            # warn "$fig have pegs for $contig $start $end ", Dumper($pegs);
            #
            # We assign a peg index to each peg, and store it in $self->{peg_index}.
            #
            # This is used for coloring related pegs when enabled.
            #

            my $idx = 1;
            for my $peg (@$pegs)
            {
                $self->{peg_index}->{$peg} = $idx++;
            }
        }

        for my $track (@tracks)
        {
            my $func = $track->{perl_binding_ref};

            if ($func)
            {
                my $track_output = $self->$func($track, $segment, $pegs);
                push(@$output, @$track_output);
            }
            else
            {
                die "No function found for ", Dumper($track);
            }
        }
        push(@$output, qq(</SEGMENT>\n));
    }
    return $output;
}

#
# This code is a close replica  of the features code. In fact, at some point,
# it might be possible to factor the common code into two chunks, one that determines
# a set of segments and a set of pegs and tracks in those segments to show.
#
# Until then, we'll keep them separate.
#

sub features_search
{
    my($self, $types, $feature_ids) = @_;
    
    #
    # Create a search object and perform the search.
    #

    my $fig = $self->{fig};
    my $search = new SproutSearch($fig);

    my $search_str = join(" ", @$feature_ids);
    #
    # Elide the trailing * that might be on there, gbrowse gives it meaning..
    #
    $search_str =~ s/\*$//;

    my($features, $genomes, $subsystems) = $search->search($search_str);

    # warn "Search returns features ", Dumper($features);

    #
    # Walk the list of features, determining which contigs the features are
    # on. This will create our list of segments.
    #

    my(%contig_pegs, %contig_min, %contig_max);

    for my $feature (@$features)
    {
        my($peg, $anno, $alias) = @$feature;

        my $loc = $fig->feature_location($peg);
        my $contig = $fig->contig_of($loc);

        my $beg = $fig->beg_of($loc);
        my $end = $fig->end_of($loc);

        my $min = \$contig_min{$contig};
        my $max = \$contig_max{$contig};

        $$min = $beg if !defined($$min) or $beg < $$min;
        $$max = $end if !defined($$max) or $end > $$max;

        push(@{$contig_pegs{$contig}}, [$peg, $loc, $beg, $end]);
    }


    #
    # Now we can create our segment list from the contigs.
    #

    my $segments = [];

    for my $contig (keys(%contig_pegs))
    {
        my $pegs = $contig_pegs{$contig};
        my $min = $contig_min{$contig};
        my $max = $contig_max{$contig};

        #
        # The undef here makes this code Sprout-specific.
        #
        my $len = $fig->contig_ln(undef, $contig);

        # warn "contig $contig: min=$min max=$max\n";
        
        #
        # Ensure we have at least a 10k wide context.
        #
        if ($max - $min < 10000)
        {
            my $pad = int((10000 - ($max - $min)) / 2);

            $min -= $pad;
            $min = 0 if $min < 0;

            $max += $pad;

            $max = $len if $max > $len;

            # warn "pad=$pad min=$min max=$max\n";
        }

        # Each segment is a list [$reference, $class, $start, $end]

        my $seg_id = $contig;
        $seg_id =~ s/:/--/g;
        push(@$segments, [$seg_id, "CDS", $min, $max, $pegs]);
    }

    if (!$types or @$types == 0)
    {
        $types = ['CDS:curated'];
    }
    else
    {
        return $self->emit_segments_only_for_match($segments);
    }
    my %types;

    #
    # Build a hash of the types we're looking for.
    #

    map { $types{$_}++ } @$types;

    $fig = $self->{fig};

    # warn "features: segments=", Dumper($segments);
    # warn "          types=", Dumper($types);
    # warn "          feature_ids=", Dumper($feature_ids);

    #
    # Determine what tracks we need to extract by matching against
    # the list of types requested. Also determine, for the tracks
    # we are keeping, if we need to map segments to features.
    #

    my $need_features = 0;
    
    my @bindings = values(%{$self->{bindings}});
    my @all_tracks = map { $_->[1] } @bindings;
    my @tracks = ();

    # warn "Types: ",  Dumper(\%types);

    #
    # We keep track of all the tracks of method Related so that
    # they can be processed first. This is a little ugly, but
    # there is some cross-track processing that needs to be supported.
    #

    my @related_tracks;
    for my $track (@all_tracks)
    {
        # warn "Testing track ", Dumper($track);
        if ($types{$track->{type}})
        {
            # warn "Adding track $track->{type}\n";

            if ($track->{method} eq "Related")
            {
                push(@related_tracks, $track);
            }
            else
            {
                push(@tracks, $track);
            }
            $need_features++ if $track->{need_pegs};
        }
    }

    #
    # Stuff the related tracks onto the front of the track list.
    #

    unshift(@tracks, @related_tracks);

    #
    # Walk the segments. For each segment, we emit the segment header,
    # and determine if we have any active tracks. If we do,
    # invoke that code to generate features. Then emit the segment footer.
    #
    # At the start of each segment we determine if we have any peg mapping to
    # determine, and do it once.
    #

    # warn "have segs ",  Dumper($segments);
    my $output = [];
    for my $segment (@$segments)
    {
        #
        # Each segment is a list [$reference, $class, $start, $end, $peg_list]
        #
        # $reference is the contig.
        # $class is the class of this segment.
        # $start and $end are the starting and ending points on the contig.
        #

        my($seg_id, $class, $start, $end, $peg_list) = @$segment;

        my $pegs = [map { $_->[0] } @$peg_list];

        # warn "segment $seg_id has pegs ", Dumper($pegs);

        my $contig = $seg_id;
        $contig =~ s/--/:/g;

        push(@$output, qq(<SEGMENT id="$seg_id" start="$start" stop="$end" version="1.0">));
        
        my $idx = 1;
        for my $peg (@$pegs)
        {
            $self->{peg_index}->{$peg} = $idx++;
        }

        for my $track (@tracks)
        {
            my $func = $track->{perl_binding_ref};

            if ($func)
            {
                my $track_output = $self->$func($track, $segment, $pegs);
                push(@$output, @$track_output);
            }
            else
            {
                die "No function found for ", Dumper($track);
            }
        }
        push(@$output, qq(</SEGMENT>));
    }
    return $output;
}

sub emit_segments_only_for_match
{
    my($self, $segments) = @_;

    my $output = [];
    
    for my $segment (@$segments)
    {
        #
        # Each segment is a list [$reference, $class, $start, $end, $peg_list]
        #
        # $reference is the contig.
        # $class is the class of this segment.
        # $start and $end are the starting and ending points on the contig.
        #

        my($seg_id, $class, $start, $end, $peg_list) = @$segment;

        for my $peg_ent (@$peg_list)
        {
            my($peg, $loc, $peg_start, $peg_end) = @$peg_ent;

            my $contig = $seg_id;
            $contig =~ s/--/:/g;

            my $label = $self->get_peg_name($peg);

            if ($peg_start > $peg_end)
            {
                ($peg_start, $peg_end) = ($peg_end, $peg_start);
            }

            push(@$output, qq(<SEGMENT id="$seg_id" start="$peg_start" stop="$peg_end" Note="$label" version="1.0">));
            push(@$output, qq(</SEGMENT>));
        }
    }
    return $output;
}


#
# Track binding functions
#

=pod

=head1 track_cds

Generate the CDS features for a segment

=cut

sub track_cds
{
    my($self, $track, $segment, $pegs) = @_;

    #
    # $track is the track configuration hash.
    # $segments is a list of segments, including the SEED features on that segment.
    # @$pegs is a list of pegs in that segment.
    #

    my $fig = $self->{fig};
    my $type = $track->{type};
    my $method = $track->{method};
    my $category = $track->{category};

    # warn "track_cds: track=", Dumper($track);
    # warn "track_cds: segment=", Dumper($segment);
    # warn "track_cds: pegs=", Dumper($pegs);

    my($feature_type) = @{$track->{perl_binding_args}};

    my $features = [];

    my($contig, $class, $start, $end) = @$segment;

    foreach my $peg (@$pegs)
    {
        next unless $peg =~ /\.$feature_type\./;
        my $loc = $self->get_feature_location($peg);
        my $start = $fig->beg_of($loc);
        my $stop = $fig->end_of($loc);
        my $orientation;

        # warn "Checking peg $peg\n";

        if ($start < $stop)
        {
            $orientation = '+';
        }
        else
        {
            $orientation = '-';
            my $t = $stop;
            $stop = $start;
            $start = $t;
        }
        
        my $len = $stop - $start + 1;
        my $phase = $len % 3;

        #
        # Determine function and label.
        #

        my($func, $name);

        $name = $self->get_peg_name($peg);
        $func = $self->get_peg_func($peg);

        my $score;
        if (defined($self->{peg_in_related_track}->{$peg}))
        {
           $score = $self->{peg_index}->{$peg};
        }
        else
        {
            $score = 0;
        }

        #
        # oh, don't look too closely here.
        #
        # we want to use relative urls here, so that proxying does the right thing.
        #
        # These links are exposed from gbrowse, where the urls look like
        # http://www.nmpdr.org/nmpdr_v3/FIG/gbrowse.cgi/GB_195099.1?ref=195099.1--CP000025;start=925978;stop=966477;nav4=1;plugin=
        #
        # the brwoser doesn't know that /GB_195099.1 isn't a real part of the path, so
        # the usual relative links don't work. we put in .. to make that work.
        # 

        my $link_url;
        if ($feature_type eq "peg")
        {
            $link_url = "../protein.cgi?prot=$peg&SPROUT=1";
        }
        else
        {
            $link_url = "../feature.cgi?prot=$peg;feature=$peg&SPROUT=1";
        }
        my $link = qq(<LINK href="$link_url">$peg</LINK>);
        
        my $feature = <<END;
    <FEATURE id="Gene:$peg" label="$name">
        <TYPE id="$type" category="$category" reference="no" subparts="no">$type</TYPE>
        <METHOD id="$method">$method</METHOD>
        <START>$start</START>
        <END>$stop</END>
        <SCORE>$score</SCORE>
        <ORIENTATION>$orientation</ORIENTATION>
        <NOTE>$func</NOTE>
        <PHASE>$phase</PHASE>
        $link
    </FEATURE>
END
        push(@$features, $feature);
    }
    return $features;
}


#
# Return a track containing features that have a given attribute.
#
# If the feature has the attribute, the score contains the value of that attribute.
#
sub track_attribute
{
    my($self, $track, $segment, $pegs) = @_;

    #
    # $track is the track configuration hash.
    # $segments is a list of segments, including the SEED features on that segment.
    # @$pegs is a list of pegs in that segment.
    #

    my $fig = $self->{fig};
    my $type = $track->{type};
    my $method = $track->{method};
    my $category = $track->{category};

    my($attribute) = @{$track->{perl_binding_args}};

    my $features = [];

    my($contig, $class, $start, $end) = @$segment;

    foreach my $peg (@$pegs)
    {
        # Old mechansim.
        #
        # my $val = $fig->get_attribute($peg, $attribute);
        # next unless defined($val);
        # next unless $val > 0;
        #
        # Now we use feature_attributes, which returns the attributes for a feature
        # as a list of triples [att-name, value, associated-url].
        # RAE: now a list of fourples: [peg, att-name, value, associated-url]
        #

        my @atts = $fig->feature_attributes($peg);

        #
        # Search @atts for the one we're looking for with a true value.
        #
        
        @atts = grep { $_->[0] eq $attribute and $_->[1] != 0 } @atts;

        next unless @atts > 0;
        # RAE: I changed $fig->feature_attributes to return peg, tag, value, url
        #my $val = $atts[0]->[1];
        #my $url = $atts[0]->[2];
        
        my $val = $atts[0]->[2];
        my $url = $atts[0]->[3];

        my $link;

        if ($url)
        {
            $link = qq(<LINK href="$url">$peg</LINK>);
        }
        
        #
        # Computed val, finish up the rest of the feature computation.
        #

        my $loc = $self->get_feature_location($peg);
        my $start = $fig->beg_of($loc);
        my $stop = $fig->end_of($loc);
        my $orientation;

        $orientation = ".";
        if ($start > $stop)
        {
            my $t = $stop;
            $stop = $start;
            $start = $t;
        }
        
        my $len = $stop - $start + 1;
        my $phase = $len % 3;

        my $name = $self->get_peg_name($peg);
        
        my $feature = <<END;
    <FEATURE id="Gene:$peg" label="$name">
        <TYPE id="$type" category="$category" reference="no" subparts="no">$type</TYPE>
        <METHOD id="$method">$method</METHOD>
        <START>$start</START>
        <END>$stop</END>
        <SCORE>$val</SCORE>
        <ORIENTATION>$orientation</ORIENTATION>
        <PHASE>$phase</PHASE>
        $link
    </FEATURE>
END
        push(@$features, $feature);
    }
    return $features;
}

#
# Return a track containing pegs grouped into functional clusters.
#
# Two pegs with the same score are in the same cluster.
#
sub track_functional_cluster
{
    my($self, $track, $segment, $pegs) = @_;

    #
    # $track is the track configuration hash.
    # $segments is a list of segments, including the SEED features on that segment.
    # @$pegs is a list of pegs in that segment.
    #

    my $fig = $self->{fig};
    my $type = $track->{type};
    my $method = $track->{method};
    my $category = $track->{category};

    my $features = [];

    my($contig, $class, $start, $end) = @$segment;

    my(%pool, %assigned);

    map { $pool{$_}++ } @$pegs;
    
    #
    # Compute the clusters for these pegs.
    #
    
    my %clusters;
    
    for my $peg (@$pegs)
    {
        my @cluster = grep { defined($pool{$_}) }  $fig->in_cluster_with($peg);
        $clusters{$peg} = [$peg, int(@cluster), [@cluster]];
    }
    
    
    my $score = 1;
    while (%pool > 0)
    {
        #
        # Scan the remaining pegs, and find the largest cluster.
        #
        
        my @peg_clusters;
        
        #
        # Need to remove pegs that have already been assigned
        # from the original clusters.
        #
        
        for my $c (@clusters{keys(%pool)})
        {
            my($peg, $len, $cluster) = @$c;
            my @new_cluster = grep { defined($pool{$_}) } @$cluster;
            push(@peg_clusters, [$peg, int(@new_cluster), [@new_cluster]]);
        }
        
        # print "peg cluster len=", int(@peg_clusters), "\n";
        
        my $winner = (sort { $b->[1] <=> $a->[1] } @peg_clusters)[0];
        # print "Winner=", Dumper($winner);
        
        my($peg, $count, $cluster) = @$winner;
        
        delete $pool{$peg};
        
        # print "Cluster is @$cluster\n";
        
        for my $c (@$cluster)
        {
            if ($pool{$c})
            {
                delete $pool{$c};
                $assigned{$c} = $score;
            }
        }
        $assigned{$peg} = $score;
        $score++;
    }

    foreach my $peg (@$pegs)
    {
        my $val = $assigned{$peg};
        next unless defined($val);
        next unless $val > 0;

        my $loc = $self->get_feature_location($peg);
        my $start = $fig->beg_of($loc);
        my $stop = $fig->end_of($loc);
        my $orientation;

        if ($start < $stop)
        {
            $orientation = '+';
        }
        else
        {
            $orientation = '-';
            my $t = $stop;
            $stop = $start;
            $start = $t;
        }
        
        my $len = $stop - $start + 1;
        my $phase = $len % 3;

        my $name = $self->get_peg_name($peg);
        
        my $feature = <<END;
    <FEATURE id="Gene:$peg" label="$name">
        <TYPE id="$type" category="$category" reference="no" subparts="no">$type</TYPE>
        <METHOD id="$method">$method</METHOD>
        <START>$start</START>
        <END>$stop</END>
        <SCORE>$val</SCORE>
        <ORIENTATION>$orientation</ORIENTATION>
        <PHASE>$phase</PHASE>
    </FEATURE>
END
        push(@$features, $feature);
    }
    return $features;
}


sub track_contig
{
    my($self, $track, $segment, $pegs) = @_;

    #
    # $track is the track configuration hash.
    # $segments is a list of segments, including the SEED features on that segment.
    # @$pegs is a list of pegs in that segment.
    #

    my $fig = $self->{fig};
    my $type = $track->{type};
    my $method = $track->{method};
    my $category = $track->{category};

    my $features = [];

    my($contig, $class, $start, $end) = @$segment;

    my $feature = <<EOF;
    <FEATURE id="sequence:$contig/1" label="$contig">
        <TYPE id="$type"  category="$category" reference="no" subparts="no">$type</TYPE>
        <METHOD id="Component">Component</METHOD>
        <START>$start</START>
        <END>$end</END>
        <SCORE>0</SCORE>
        <ORIENTATION>+</ORIENTATION>
        <PHASE>0</PHASE>
        <GROUP id="Sequence:$contig" type="Sequence" />
    </FEATURE>
EOF

    push(@$features, $feature);

    return $features;
}


#
# Map a set of aliases to a canonical name, based on a hierarchy
# of naming preferences.
#
sub anti_alias() {
    my(@aliases) = @_;

    my(%names);

    for $_ (@aliases)
    {
        #
        # Match the various forms of names, and mark in %matches.
        #

        if (/^gi\|/) {
            $names{"gi"} = $_;
        }
        elsif (/^uni\|/)
        {
            $names{"uni"} = $_;
        }
        elsif (/^[a-z]{2,3}[A-Z][0-9]?/)
        {
            $names{"gene_name"} = $_;
        }
        elsif (/^[A-Za-z]{1,2}\d+$/)
        {
            $names{seq} = $_;
        }
        elsif (/^[A-Z]{2}_/) {
            $names{"ref"} = $_;
        }
        elsif (/^fig\|/)
        {
            $names{fig} = $_;
        }

    }

    # warn "generated map ", join(" ", %names), "\n";

    my(@prec) = qw(gene_name seq uni ref fig kegg);
 
    my @matches = grep { $_ ne "" } @names{@prec};

    # warn "mapped to @matches\n";
    return $matches[0];
}

=pod

4. Bruce and Bob need to help implement a "track" for showing each
       genome related to the given one.  The logic for this track is based
       on a subroutine that takes as input

          the given genome X
          the given contig
          the begin (leftmost coord) of the map B
          the end of the map E
          a second genome Y

       The routine must return a set of features with coordinates in the
       range begin-end (i.e., coordinates in the given genome's contig)

   This routine should have roughly this logic:

        a. Let S1 be the features in X in the region B-E
        b. Let O1 be the BBHs of those features in Y
        c. Let O1r be the largest subset of features in O1 that occur
           within a single contig and within a distance of E-B.
        d. pick the bounding pairs (of S1 entries with O1r entries).
           Align the midpoints of each bounding pair.  Now find all
           genes in Y that fall within the appropriate range around
           the aligned midpoints (remember the genomes may have
           reversed directions, so pick the appropriate coords).
        e. adjust all the coordinates to fall in the given genomes
           range.
        f. return the features.

=cut

sub track_related
{
    my($self, $track, $segment, $pegs) = @_;

    # warn "Track_related ", Dumper(\@_);
    # warn "Segment: ", Dumper($segment);
    # warn "pegs ", Dumper($pegs);
    #
    # $track is the track configuration hash.
    # $segments is a list of segments, including the SEED features on that segment.
    # @$pegs is a list of pegs in that segment.
    #

    my $fig = $self->{fig};
    my $type = $track->{type};
    my $method = $track->{method};
    my $category = $track->{category};
    
    my($y) = @{$track->{perl_binding_args}};

    my $features = [];
    
    my($contig, $class, $start, $end) = @$segment;

    my $s1 = $pegs;

    my %contig;

    #
    # Remember the index of the related peg.
    #

    my %related_index;
    my %related_peg_to_original_peg;

    my %org_peg_info;

    my($my_peg_min, $my_peg_max);

    for my $peg (@$s1)
    {
        # print "Handling $peg\n";

        my $peg_bbhs = $fig->bbh_list($y, [$peg]);

        # warn "bbhs for $peg: ", Dumper ($peg_bbhs);

        #
        # Is it one peg or a list?
        #
        my @o1;
        my $tmp = $peg_bbhs->{$peg};
        if (ref($tmp) eq "ARRAY")
        {
            @o1 = @$tmp;
        }
        elsif ($tmp ne "")
        {
            @o1 = ($tmp);
        }           

        for my $o (@o1)
        {
            $related_peg_to_original_peg{$o} = $peg;
            
            my $loc = $self->get_feature_location($o);
            # warn "o=$o loc=$loc\n";

            if ($loc eq '')
            {
                #
                # BUG.
                #
                # We hit a feature not in the db. it may have been deleted
                # (the bug as of 2005-0830 is that we renamed PEGs and
                # didnt recompute bbhs).
                #

                next;
                
            }

            my $b = $fig->beg_of($loc);
            my $e = $fig->end_of($loc);
            my $contig = $fig->contig_of($loc);

            #
            # Determine min and max and stash those for later use.
            #

            my($min, $max);
            if ($b < $e)
            {
                $min = $b;
                $max = $e;
            }
            else
            {
                $min = $e;
                $max = $b;
            }

            push(@{$contig{$contig}}, [$o, $b, $e, $min, $max]);

            $related_index{$o} = $self->{peg_index}->{$peg};
        }

        #
        # While we're scanning, compute min/max extent of our pegs.
        #

        my $ploc = $self->get_feature_location($peg);

        #
        # Work around the bug here too.
        #
        next if $ploc eq '';
        
        my $b = $fig->beg_of($ploc);
        my $e = $fig->end_of($ploc);

        my($max, $min);
        if ($b < $e)
        {
            $min = $b;
            $max = $e;
        }
        else
        {
            $min = $e;
            $max = $b;
        }

        $org_peg_info{$peg} = [$b, $e, $min, $max];
        
        if (!defined($my_peg_min))
        {
            $my_peg_min = $min;
            $my_peg_max = $max;
        }
        else
        {
            $my_peg_min = $min if $min < $my_peg_min;
            $my_peg_max = $max if $max > $my_peg_max;
        }
             
    }

    # warn "Contig info ", Dumper(\%contig);

    #
    # Pick the contig with the largest set of pegs.
    #

    #
    # Flatten contig map, generating an array @x
    # of pairs [contig-name, list of [peg, beg, end, min, max] tuples]
    #
    my @x = map { [$_, $contig{$_}] } keys(%contig);

    #
    # If we found no bbhs, just return.
    #
    return [] if @x == 0;

    #
    # Sort by number of hits per contig.
    #

    @x = sort { @{$b->[1]} <=> @{$a->[1]} } @x;

    #
    # And select the largest.
    #

    my($o1r_contig, $o1r_pegs) = @{$x[0]};

    # print Dumper($o1r_contig, $o1r_pegs);

    #
    # We now have a candidate contig and a set of pegs on that contig.
    # The peg set is a list of tuples [peg, start, stop, min, max].
    # start > stop if on the - strand.
    #

    #
    # Sort the candidates by lowest-position.
    #

    my @csorted = sort { $a->[3] <=> $b->[3] } @$o1r_pegs;

    # print Dumper(\@csorted);

    #
    # Now scan them in order, creating an array that contains
    # for each element, the number of pegs in a span-wide section that
    # overlap.
    #

    my @overlaps;

    #
    # Remember, $start and $end bound the contig we're searching from.
    #

    my $search_size = $end - $start + 1;
    
    for my $base_idx (0..@csorted - 1)
    {
        my $dist = 0;
        my $base = $csorted[$base_idx];
        my $reach = $base->[4] + $search_size;
        my $end_idx;

        # warn "Base: $base_idx reach=$reach base=@$base\n";
        
        for my $extent_idx ($base_idx + 1.. @csorted - 1)
        {
            my $extent = $csorted[$extent_idx];

            # warn "Look: $extent_idx extent=@$extent\n";

            #
            # If this peg is out of range, we're done.
            #
            if ($extent->[3] > $reach)
            {
                last;
            }
            $end_idx = $extent_idx;
        }
        if ($end_idx)
        {
            # print "Reached end $end_idx\n";
            push(@overlaps, [$base_idx, $end_idx, $end_idx - $base_idx]);

            #
            # Special case: If we cover all of the pegs, we don't have to search any more
            #

            if ($base_idx == 0 and $end_idx == (@csorted - 1))
            {
                last;
            }
        }
    }

    #
    # Overlaps is now a list of tuples (starting peg, ending peg, difference);
    # Sort on the number of pegs there. The maximal size is the set
    # we will map to.
    #

    my @sorted_overlaps = sort { $b->[2] <=> $a->[2] } @overlaps;

    return [] if @sorted_overlaps == 0;

    my ($start_peg, $end_peg) = @{$sorted_overlaps[0]};
    
    # print "got $start_peg $end_peg ", Dumper(\@sorted_overlaps);

    #
    # We now have the set of target pegs.
    # Compute the midpoint of this set of pegs, and the midpoint
    # of the related pegs. This defines the translation
    # between the contigs.
    #

    my @related_pegs = @csorted[$start_peg .. $end_peg];

    # print Dumper(\@related_pegs);

    #
    # Scan our pegs, searching for the min and max points.
    #

    my($rel_min, $rel_max, $org_min, $org_max);
    $rel_min = $org_min = 1e10;
    $org_max = $rel_max = -1;

    my $flip;
    
    for my $rel_info (@related_pegs)
    {
        my($rel_peg, $rel_b, $rel_e, $rel_pmin, $rel_pmax) = @$rel_info;
        my $org_peg = $related_peg_to_original_peg{$rel_peg};
        
        my($org_b, $org_e, $org_pmin, $org_pmax) = @{$org_peg_info{$org_peg}};

        print "Rel peg: @$rel_info\n";
        print "Org peg: $org_peg  @{$org_peg_info{$org_peg}}\n";

        if ($rel_pmin < $rel_min)
        {
            $rel_min = $rel_pmin;
            $org_min = $org_pmin;
        }

        if ($rel_pmax > $rel_max)
        {
            $rel_max = $rel_pmax;
            $org_max = $org_pmax;
        }
        
        #$org_min = $org_pmin if $org_pmin < $org_min;
        #$org_max = $org_pmax if $org_pmax > $org_max;

        $self->{peg_in_related_track}->{$org_peg}++;

        #
        # Check one peg to see if we need to flip.
        #

        if (!defined($flip))
        {
            my $rel_strand = $rel_b > $rel_e;
            my $org_strand = $org_b > $org_e;
            $flip = ($rel_strand == $org_strand) ? 1 : -1;
        }
    }

    my $my_center = int(($org_min + $org_max) / 2);
    my $rel_center = int(($rel_min + $rel_max) / 2);

    #
    # The centers induce the mapping of peg locations.
    #
    # Given a location in related-org coords L,
    # the corresponding location in  my-org coords
    # is (L - $rel_center) + $my_center
    # or L + ($my_center - $rel_center)
    #
    # However, in order to support flip, we need this:
    #
    # $org = $my_center + flip * (L - $rel_center)
    #
    # To map backwards:
    # 
    # L = $flip * ($org - $my_center) + $rel_center
    #

    my $offset = int($my_center - $rel_center);


    #
    #  Map from related coords to orig coords
    #
    my $map = sub {
        my($val) = @_;
        return $my_center + $flip * ($val - $rel_center);
    };

    #
    # Map from org coords to related coords
    #
    my $revmap  = sub {
        my($val) = @_;
        return ($val - $my_center) + $rel_center;
    };
    
    #
    # Now walk the pegs in the area corresponding to the chosen contig,
    # returning translated features
    #


    # warn "my_center=$my_center rel_center=$rel_center offset=$offset\n";
    my ($rel_pegs, $rel_start, $rel_end) = $self->{fig}->genes_in_region($y,
                                                                         $o1r_contig,
                                                                         &$revmap($start),
                                                                         &$revmap($end));

    # warn Dumper($rel_pegs);
    # warn Dumper(\%related_index);
    # warn Dumper($self);
    for my $peg (@$rel_pegs)
    {
        next if $fig->ftype($peg) ne "peg";
        my $loc = $self->get_feature_location($peg);

        my $pstart = $fig->beg_of($loc);
        my $pstop = $fig->end_of($loc);

        my($orientation, $max, $min);

        if ($pstart < $pstop)
        {
            $max = $pstop;
            $min = $pstart;
            $orientation = "+";
        }
        else
        {
            $max = $pstart;
            $min = $pstop;
            $orientation = "-";
        }

        my $len = $max - $min + 1;
        my $name = $self->get_peg_name($peg);


        if ($flip == -1)
        {
            $orientation = $orientation eq "-" ? "+" : "-";
        }
        $pstart = &$map($pstart);
        $pstop = &$map($pstop);
        
        my $score = $related_index{$peg};
        $score = 0 unless defined($score);
        my $phase = $len % 3;

        #
        # This computes links that go to the SEED protein page.
        #

        my $protein_url = "../protein.cgi?prot=$peg&SPROUT=1";

        #
        # Compute a link back to gbrowse, but on *this* organism instead
        # of the reference. We keep the same zoom factor by mapping
        # the endpoints of the reference contig area to this org's contig.
        #

        {
            my $genome;
            $genome = $1 if $peg =~ /fig\|(\d+\.\d+)/;
                
            my $lstart = &$revmap($start);
            my $lend = &$revmap($end);
            my $seg = $o1r_contig;
            $seg =~ s/:/--/g;
            $protein_url = "../gbrowse.cgi/GB_$genome?ref=$seg&start=$lstart&stop=$lend";
        }
        
        my $link = qq(<LINK href="$protein_url">$name</LINK>);

        my $feature = <<END;
    <FEATURE id="Gene:$peg" label="$name">
        <TYPE id="$type" category="$category" reference="no" subparts="no">$type</TYPE>
        <METHOD id="$method">$method</METHOD>
        <START>$pstart</START>
        <END>$pstop</END>
        <SCORE>$score</SCORE>
        <ORIENTATION>$orientation</ORIENTATION>
        <PHASE>$phase</PHASE>
        $link
    </FEATURE>
END
        push(@$features, $feature);
            
    }
    
    return $features;
}

#
# Some utility routines.
#

sub get_peg_func
{
    my($self, $peg) = @_;

    if (my $cached = $self->{cache_peg_func}->{$peg})
    {
        return $cached;
    }
    
    my $func = $self->{fig}->function_of($peg);
    $self->{cache_peg_func}->{$peg} = $func;
    return $func;
}


sub get_peg_name
{
    my($self, $peg) = @_;

    if (my $cached = $self->{cache_peg_name}->{$peg})
    {
        return $cached;
    }
    
    my @aliases = $self->{fig}->feature_aliases($peg);
    my $name = &anti_alias(@aliases);
    if (!$name) {
        $name = $peg;
        $name =~ s,^fig\|\d+\.\d+\.peg\.,,;
    }
    $self->{cache_peg_name}->{$peg} = $name;

    return $name;
}

sub get_feature_location
{
    my($self, $peg) = @_;

    if (my $cached = $self->{cache_feature_location}->{$peg})
    {
        return $cached;
    }
    
    my $loc = $self->{fig}->feature_location($peg);
    $self->{cache_feature_location}->{$peg} = $loc;
    return $loc;
}

1;
