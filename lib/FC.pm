
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
package FC;

=head1 Functional Coupling Methods

=cut

use strict;
use ERDB;
use Data::Dumper;
use Carp;
use SeedUtils;
use PartitionSeqs;
use Tracer;
use FIG;
no warnings 'redefine';

sub co_occurring_FIGfams {
    my($db,$ff) = @_;

    my $query = $db->Get("Family HasMember Feature IsInPair Pairing Determines PairSet","Family(id) = ?",[$ff]);

    my(%hits,%pairsets,$row,$pairset,$sc);
    while ($row = $query->Fetch()) 
    {
	($pairset,$sc) = $row->Values(['PairSet(id)','PairSet(score)']);
	$pairsets{$pairset} = $sc;
    }

    while (($pairset,$sc) = each(%pairsets))
    {
	$query = $db->Get("IsDeterminedBy Pairing IsPairOf Feature IsMemberOf Family","IsDeterminedBy(from-link) = ?",[$pairset]);
	while ($row = $query->Fetch()) 
	{
	    my($ff1,$fam_func) = $row->Values(['IsMemberOf(to-link)',
                                               'Family(family-function)'
                                              ]);
	    if ($ff1 ne $ff)
	    {
		if ((! $hits{$ff1}) || ($hits{$ff1}->[0] < $sc))
		{
		    $hits{$ff1} = [$sc,$fam_func];
		}
	    }
	}
    }
    return map { [$_,$hits{$_}] } sort { $hits{$b} <=> $hits{$a} } keys(%hits);
}

sub co_occurs {
    my($db,$peg) = @_;

    my @tuples;
    my $query = $db->Get("IsInPair Pairing Determines PairSet",
			 "IsInPair(from-link) = ?", [$peg]);
    while (my $row = $query->Fetch()) 
    {
	my ($pair, $pairset, $score) = $row->Values(['Determines(from-link)','PairSet(id)','PairSet(score)']);
	my ($fid) = grep { $_ ne $peg } split(/:/, $pair);
	push @tuples, [$score,$fid,$pairset];
    }
    return @tuples;
}

sub co_occurs_batch {
    my($db,$pegs) = @_;

    my @tuples;
    my $query = $db->Get("IsInPair Pairing Determines PairSet",
			 "IsInPair(from-link) IN (" . join(", ", map { "'$_'" } @$pegs) . ")", []);
    while (my $row = $query->Fetch()) 
    {
	my ($peg1,$pair,$score) = $row->Values(['IsInPair(from-link)','Determines(from-link)','PairSet(score)']);
	my($peg2) = grep { $_ ne $peg1 } split(/:/,$pair);
	push(@tuples, [$peg1,$peg2,$score]);
    }
    return wantarray ? @tuples : \@tuples;
}

sub coupled_to {
    my($db,$peg) = @_;
    return map { [$_->[1],$_->[0]] } &FC::co_occurs($db,$peg);
}

sub coupled_to_batch {
    my($db,@pegs) = @_;

    return &FC::co_occurs_batch($db,\@pegs);
}

sub in_co_occurrence_cluster {
    my($db,$peg) = @_;

    my $query = $db->Get("OccursIn Cluster IsOccurrenceOf",
			 "OccursIn(from-link) = ?", [$peg]);
    my @pegs = ();    
    while (my $row = $query->Fetch())
    {
	my($peg1) = $row->Values(['IsOccurrenceOf(to-link)']);
	push(@pegs,$peg1);
    }
    return (@pegs > 0) ? \@pegs : undef;
}    

sub co_occurrence_evidence {
    my($db,$peg1,$peg2) = @_;
    my ($flipped, $pairing_id);
    if ($peg1 gt $peg2)
    {
	$flipped = 1;
	$pairing_id = "$peg2:$peg1";
    } else {
	$flipped = 0;
	$pairing_id = "$peg1:$peg2";
    } 
    my $query = $db->Get("Pairing Determines PairSet IsDeterminedBy",
			 "Pairing(id) = ?", [$pairing_id]);
    my $tuples;
    while (my $row = $query->Fetch()) 
    {
	my($pair_id,$inverted,$original_inverted) = $row->Values(['IsDeterminedBy(to-link)','IsDeterminedBy(inverted)','Determines(inverted)']);
	my($peg3,$peg4) = split(":",$pair_id);
	# We need to know if ($peg1, $peg2) is oriented the same as the prototypes of this
	# pairset.
	my $notOriented = ($flipped == $original_inverted ? 0 : 1);
	# If ($peg3, $peg4) is NOT oriented in the same relationship, we flip it.
	if ($notOriented != $inverted) { ($peg3,$peg4) = ($peg4,$peg3) }
	# Only include this result if it's NOT our original.
	if ($peg3 ne $peg1)
	{
	    push(@$tuples,[$peg3,$peg4]);
	}
    }
    return $tuples;
}

sub coupling_evidence {
    my($db,$peg1,$peg2) = @_;

    my $tuples = &co_occurrence_evidence($db,$peg1,$peg2);
    if (! $tuples) { return () }

    my $i;
    my %seen;
    for ($i=0; ($i < @$tuples); $i++)
    {
	my($peg3,$peg4) = @{$tuples->[$i]};
	$peg3 =~ /^fig\|(\d+\.\d+)/;
	my $genome = $1;
	my $otu    = $db->OTU($genome);
	if (! $seen{$otu})
	{
	    $seen{$otu} = 1;
	    push(@{$tuples->[$i]},1);
	}
	else
	{
	    push(@{$tuples->[$i]},0);
	}
    }
    return wantarray ? @{$tuples} : $tuples;
}

sub co_occurrence_and_evidence {
    my($db,$peg) = @_;

    my @tuples;
    my $query = $db->Get("IsInPair Pairing Determines PairSet IsDeterminedBy",
			 "IsInPair(from-link) = ?", [$peg]);
    my @tuoples = ();
    while (my $row = $query->Fetch()) 
    {
	push(@tuples,[$row->Values(['Determines(from-link)',
				    'PairSet(score)',
                                    'IsDeterminedBy(inverted)',
				    'IsDeterminedBy(to-link)'])]);
    }
    @tuples = sort { $a->[0] <=> $b->[0] } @tuples;
    my @co_occurs = ();

    my $x = shift @tuples;
    while ($x)
    {
	my $curr = $x->[0];
	my $sc   = $x->[1];
	my @set = ();
	while ($x && ($x->[0] eq $curr))
	{
	    my($peg3,$peg4) = split(/:/,$x->[3]);
	    push(@set,($x->[2]) ? [$peg3,$peg4] : [$peg4,$peg3]);
	    $x = shift @tuples;
	}
	my $i;
	for ($i=0; ($i < @set) && ($peg ne $set[$i]->[0]); $i++) {}
	if ($i == @set)
	{
	    @set = map { [$_->[1],$_->[0]] } @set;
	}
	my($peg2) = grep { $_ ne $peg } split(/:/,$curr);
	push(@co_occurs,[$sc,$peg2,[grep { $_->[0] ne $peg } @set]]);
    }
    return sort { $b->[0] <=> $a->[0] } @co_occurs; 
}

sub co_occurrence_set {
    my($db,$set) = @_;

    my $pairs;
    my $sc;

    my $query = $db->Get("PairSet",
			 "PairSet(id) = ?", [$set]);
    my $row = $query->Fetch();
    if ($row)
    {
	($sc) = $row->Values(['PairSet(score)']);
	$pairs = [];
	my $query = $db->Get("PairSet IsDeterminedBy",
			     "PairSet(id) = ?", [$set]);
	while (my $row = $query->Fetch())
	{
	    my($pairing,$inverted) = $row->Values(['IsDeterminedBy(to-link)','IsDeterminedBy(inverted)']);
	    my($peg1,$peg2) = split(/:/,$pairing);
	    push(@$pairs,$inverted ? [$peg2,$peg1] : [$peg1,$peg2]);
	}
    }
    return ($sc,$pairs);
}

sub all_co_occurrence_pair_sets {
    my($db) = @_;

    my $query = $db->Get("PairSet","",[]);
    my @PairSets = ();
    while (my $row = $query->Fetch())
    {
	my($pair_set) = $row->Values(['PairSet(id)']);
	push(@PairSets,$pair_set);
    }
    return @PairSets;
}

sub all_co_occurrence_clusters {
    my($db) = @_;

    my $query = $db->Get("Cluster","",[]);
    my @Clusters = ();
    while (my $row = $query->Fetch())
    {
	my($id) = $row->Values(['Cluster(id)']);
	push(@Clusters,$id);
    }
    return @Clusters;
}

sub largest_co_occurrence_clusters {
    my($db,$peg) = @_;

    my %pegs = map { $_->id2 => $_->psc } sims($peg,1000,1.0e-30,'fig');

    my @pegs = sort { $pegs{$a} <=> $pegs{$b} } grep { $_ ne $peg } keys(%pegs);
    unshift @pegs,$peg;
    $pegs{$peg} = 0;
    my @clusters = ();
    foreach my $peg1 (@pegs)
    {
	my $cluster = &in_co_occurrence_cluster($db,$peg1);
	if (defined($cluster))
	{
	    my $tuple = [$peg1,$pegs{$peg1}, [grep { $_ ne $peg1 } @$cluster]];
	    push(@clusters,$tuple);
	}
    }
    return sort { @{$b->[2]} <=> @{$a->[2]} } @clusters;
}

sub co_occurrence_cluster {
    my($db,$cluster) = @_;

    my $query = $db->Get("Cluster IsOccurrenceOf",
			 "Cluster(id) = ?", [$cluster]);
    my @pegs = ();    
    while (my $row = $query->Fetch())
    {
	my($peg) = $row->Values(['IsOccurrenceOf(to-link)']);
	push(@pegs,$peg);
    }
    return @pegs;
}

#sub co_occurrence_evidence {
#    my($db,$peg1,$peg2) = @_;
#    
#    my $key = ($peg1 lt $peg2) ? "$peg1:$peg2" : "$peg2:$peg1";
#    my $query = $db->Get("Determines PairSet",
#			 "Determines(from-link) = ?", [$key]);
#    my($sc,$pairs);
#
#    if ($query)
#    {
#	my $row = $query->Fetch();
#	if ($row)
#	{
#	    my($set) = $row->Values(['PairSet(id)']);
#	    ($sc,$pairs) = &co_occurrence_set($db,$set);
#	    my $i;
#	    for ($i=0; ($i < @$pairs) && ($pairs->[$i]->[0] ne $peg1); $i++) {}
#	    if ($i == @$pairs)
#	    {
#		$pairs = [map { [$_->[1],$_->[0]] } @$pairs];
#		for ($i=0; ($i < @$pairs) && ($pairs->[$i]->[0] ne $peg1); $i++) {}
#		if ($i == @$pairs)
#		{
#		    print STDERR &Dumper($pairs,$i,$peg1,$peg2);
#		    die "Something is screwed up";
#		}
#	    }
#	    splice(@$pairs,$i,1);
#	}
#	return $pairs;
#    }
#    return [];
#}
#
#sub is_co_occurrence_pair {
#    my($db,$peg1,$peg2) = @_;
#
#    my $key = ($peg1 lt $peg2) ? "$peg1:$peg2" : "$peg2:$peg1";
#    my $query = $db->Get("Pairing",
#			 "Pairing(id) = ?",[$key]);
#    return ($query && $query->Fetch()) ? 1 : 0;
#}

sub in_pair_set {
    my($db,$peg1,$peg2) = @_;

    my($key,$inv,$set);
    $key = ($peg1 lt $peg2) ? "$peg1:$peg2" : "$peg2:$peg1";

    my $query = $db->Get("Pairing Determines",
			 "Pairing(id) = ?",[$key]);
    if ($query)
    {
	my $row = $query->Fetch();
	if ($row)
	{
	    my($pair,$inverted) = $row->Values(['Determines(to-link)','Determines(inverted)']);
	    if ($peg1 gt $peg2)
	    {
		$inverted = ! $inverted;
	    }
	    return ($pair,$inverted);
	}
    }
    return undef;
}

# This routine inserts a Pairing for $peg1 and $peg2
sub insert_pair {
    my($db,$peg1,$peg2,$set) = @_;
    # Compute the key and the inversion flag.
    my ($invertFlag, $key);
    if ($peg1 lt $peg2) {
	$key = "$peg1:$peg2";
	$invertFlag = 0;
    } else {
	$key = "$peg2:$peg1";
	$invertFlag = 1;
    }
    # Create the pairing.

    my $query = $db->Get("Pairing",
			 "Pairing(id) = ?",[$key]);
    if ((! $query) || (! $query->Fetch()))
    {
	$db->InsertObject('Pairing', id => $key);
	# Connect it to its constituent features.
	$db->InsertObject('IsInPair', from_link => $peg1, to_link => $key);
	$db->InsertObject('IsInPair', from_link => $peg2, to_link => $key);
    }

    $query = $db->Get("Determines",
		      "Determines(from-link) = ?",[$key]);
    if ($query && $query->Fetch())
    {
	return 0;   # it is already connected to another PairSet
    }

    # Insert it into the pair set.
    $db->InsertObject('IsDeterminedBy', from_link => $set, to_link => $key,
		      inverted => $invertFlag);
    Trace("inserted $peg1:$peg2 into $set") if T(3);
    return 1;
}

# This routine creates a new PairSet, returning the ID
sub new_set {
    my($db) = @_;
    # Create the new pair set.
    my $retVal = $db->InsertNew('PairSet', score => 0);
    # Return its ID.
    return $retVal;
}
    

# This routine deletes a PairSet and relationships to Pairings
# (but not the Pairings).

sub delete_PairSet {
    my($db,$pair_set) = @_;
    # Disconnect the pair set from all of its pairings.
    $db->Disconnect('IsDeterminedBy', PairSet => $pair_set);
    # Delete the pair set. Because we've already disconnected, the pairings
    # are safe.
    $db->Delete(PairSet => $pair_set, onlyRoot => 1);
    Trace("deleted PairSet $pair_set") if T(3);
}

# This routine deletes the relationship between a Pairing and a PairSet
# (but does not delete either entity).

sub delete_pair_from_PairSet {
    my($db,$pair,$PairSet) = @_;
    # Delete the relationship connecting this set to this pair.
    $db->DeleteRow('IsDeterminedBy', $PairSet, $pair);
}

# This routine updates the "score" field in a PairSet
sub set_pair_set_sc {
    my($db,$pair_set,$sc) = @_;
    # Update the score in place.
    $db->UpdateEntity(PairSet => $pair_set, score => $sc);
    Trace("reset score of PairSet to $sc") if T(3);
}

sub check_and_rescore {
    my($db,$pair_set) = @_;
    
    my $sc = &rescore($db,$pair_set);
    if ($sc < 5)
    {
	Trace("Rescored PairSet $pair_set to $sc, so we are deleting it") if T(3);
	&delete_PairSet($db,$pair_set);
	return 0;
    }
    return 1;
}

sub rescore {
    my($db,$pair_set) = @_;

    my $fig = new FIG;
    my(undef,$pairs) = &co_occurrence_set($db,$pair_set);
    defined($pairs) || confess $pair_set;

    my $sz = @$pairs;
    my %genome_sets;
    foreach my $pair (@$pairs)
    {
	$genome_sets{$fig->get_representative_genome(&FIG::genome_of($pair->[0]))} = 1;
    }
    my $sz_reps = keys(%genome_sets);
    Trace("rescoring: $sz pairs in $pair_set, score=$sz_reps") if T(3);
    my $sc = keys(%genome_sets);
    &set_pair_set_sc($db,$pair_set,$sc);
    return $sc;
}

=head2 PairSet Update Methods

This section contains methods used to update the pair-set data in the Sapling
database.

=cut

=head3 extend_pairs_and_pair_sets

    FC::extend_pairs_and_pair_sets($db, $pchF, $stats);

This method will read raw PCH data for a single genome and update the
pairs and pair sets in the Sapling database.

=over 4

=item db

L<Sapling> object for accessing the database.

=item pchF

The name of a raw PCH file. The file must be tab-delimited, each line containing
four feature IDs. The first two features are physically close, and the second
pair of features are physically close homologs at a different location. The
first and third features are similar and the second and fourth features are
similar; the set of four is considered evidence that the features involved have
related functions.

=item stats

A L<Stats> object used to track information about what happened during the
update.

=back

=cut

sub extend_pairs_and_pair_sets {
    my($db,$pchF,$stats) = @_;

    Open(\*PCHS,"sort $pchF |");

    my $line = <PCHS>;
    while ($line && ($line =~ /^((\S+)\t(\S+))\t(\S+)\t(\S+)/))
    {
	Trace($stats->Ask('PairsIn') . " pairs read.") if $stats->Check(PairsIn => 500) && T(3);
	my $curr = $1;
	my $set = [];
	while ($line && ($line =~ /^((\S+)\t(\S+))\t(\S+)\t(\S+)/) && ($1 eq $curr))
	{
	    Trace($stats->Ask('PairsIn') . " pairs read.") if $stats->Check(PairsIn => 500) && T(3);
	    push(@$set,[$2,$3,$4,$5]);
	    $line = <PCHS>;
	}
	&add_pch_set($db,$set,$stats);
    }
    close(PCHS);
}

sub add_pch_set {
    my($db,$set,$stats) = @_;

    my $fig = new FIG;
    my @pairs = ([$set->[0]->[0],$set->[0]->[1]],map { [$_->[2],$_->[3]] } @$set);

    my($pair,@unplaced,%in,%members);
    my $ok = 1;
    foreach $pair (@pairs)
    {
	my($pair_set,$inverted) = &in_pair_set($db,@$pair);
	$inverted = $inverted ? 1 : 0;   # force to {0,1} for boolean
	if (! defined($pair_set))
	{
	    push(@unplaced,$pair);
	    $stats->Add(unplaced => 1);
	}
	else
	{
	    my($keyPI) = "$pair_set:$inverted";
	    $in{"$keyPI"}++;
	    push(@{$members{$keyPI}},join(",",@$pair));

	    my $opposite = $inverted ? 0 : 1;
	    if ($in{"$pair_set:$opposite"})
	    {
		# We have a case in which a set contains both peg1:peg2 and peg2:peg1
		$stats->Add(circular => 1);
		$ok = 0;
	    }
	}
    }

    if ($ok) {

	my @in_sets = map { [split(/:/,$_)] } sort { $in{$b} <=> $in{$a} } keys(%in);
	$stats->Add(InSet => scalar(@in_sets));
    
	if (@in_sets > 1)
	{
	    my $i;
	    for ($i=1; ($i < @in_sets); $i++)
	    {
		my(undef,$in_old) = &co_occurrence_set($db,$in_sets[$i]->[0]);
		if ($in_sets[$i]->[1] ne $in_sets[0]->[1])
		{
		    $in_old = [ map { [$_->[1],$_->[0]] } @$in_old];
		}
		foreach $pair (@$in_old)
		{
		    $stats->Add(PairInserted => 1);
		    &insert_pair($db,@$pair,$in_sets[0]->[0]);
		}
		$stats->Add(SetDeleted => 1);
		&delete_PairSet($db,$in_sets[$i]->[0]);
	    }
	    $stats->Add(SetRescored => 1);
	    &check_and_rescore($db,$in_sets[0]->[0]);
	}
    
	if (@in_sets == 1)
	{
	    $stats->Add(Unplaced => scalar(@unplaced));
	    if (@unplaced > 0)
	    {
		foreach $pair (@unplaced)
		{
		    my($peg3,$peg4) = $in_sets[0]->[1] ? ($pair->[1],$pair->[0]) : @$pair;
		    $stats->Add(PairInserted => 1);
		    &insert_pair($db,$peg3,$peg4,$in_sets[0]->[0]);
		}
		$stats->Add(SetRescored => 1);
		&check_and_rescore($db,$in_sets[0]->[0]);
	    }
	}
	elsif (@unplaced >= 5)
	{
    
	    my %genome_sets;
	    foreach my $pair (@unplaced)
	    {
		$genome_sets{$fig->get_representative_genome(&FIG::genome_of($pair->[0]))} = 1;
	    }
	    my $sc = keys(%genome_sets);
	    if ($sc >= 5)
	    {
		my $new_set = &new_set($db);
		$stats->Add(NewSet => 1);
		
		if (! $new_set)
		{
		    Confess("Failed to acquire a new PairSet.");
		}
		foreach my $pair (@unplaced)
		{
		    $stats->Add(PairInserted => 1);
		    my $rc = &insert_pair($db,@$pair,$new_set);
		    if (! $rc) {
			$stats->Add(PairInsertFailed => 1);
		    }
		}
		$stats->Add(SetRescored => 1);
		&check_and_rescore($db,$new_set);
	    }
	}
    }
}

sub cleanup_pair_sets {
    my($db) = @_;

    my $fig = new FIG;
    foreach my $pair_set (&all_co_occurrence_pair_sets)
    {
	my(undef,$set) = &co_occurrence_set($db,$pair_set);
	next if (! $set);
	my @partitioned = sort { @$b <=> @$a } &partition_pair_set($fig,$db,$set);

	my $N = @partitioned;

	my $keep = shift @partitioned;

	if ((! $keep) || (@$keep < 5))
	{
	    &delete_PairSet($db,$pair_set);
	}
	else
	{
	    my $sz = @$keep;
	    my %to_keep = map { join(":",@$_) => 1 } @$keep;
	    foreach my $pair (@$set)
	    {
		my $key = join(":",@$pair);
		if (! $to_keep{ $key } )
		{
		    &delete_pair_from_PairSet($db,$pair,$pair_set);
		}
	    }
	    if (@$keep < @$set) { &check_and_rescore($db,$pair_set) }

	    foreach my $new_set (@partitioned)
	    {
		my $sz = @$new_set;
		if (@$new_set >= 5)
		{
		    my $next_set = &new_set($db);
		    foreach my $pair1 (@$new_set)
		    {
			my $key = join(":",@$pair1);
			Trace("\t\tadding $key") if T(3);
			&insert_pair($db,@$pair1,$next_set);
		    }
		    &check_and_rescore($db,$next_set);
		}
	    }
	}
    }
}

sub partition_pair_set {
    my($fig,$db,$set) = @_;

    my @split = ();
    my @sets = &part1($fig,$set,0);
    foreach my $part1 (@sets)
    {
	foreach my $part2 (&part1($fig,$part1,1))
	{
	    push(@split,$part2);
	}
    }
    return sort { @$b <=> @$a } @split;
}

sub part1 {
    my($fig,$set,$index) = @_;

    my %entries = map { $_->[$index] => $_ } @$set;
    my $ids     = [grep { $fig->is_real_feature($_) } map { $_->[$index] } @$set];
    my @partition = &PartitionSeqs::partition({ pegs => $ids, use => 'blast', identity => 0.4, coverage => 0.7 });
    return map { [map { $entries{$_} } @$_] } @partition;
}

1;
