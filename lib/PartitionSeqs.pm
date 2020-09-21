package PartitionSeqs;

use strict;
use gjoseqlib;
use FIG;
use Data::Dumper;
use FIG_Config;


# This module can be used to take a set of similar pegs and partition them,
# putting those that are quite similar and similar over most of the length
# of each sequence into the same partition.
#
#  \@sets = PartitionSeqs:partition(  \%options )
#
#  Options:
#
#      seqs     =>  [[id1,comment1,seq1],[id2,comment2,seq2]...]
#      pegs     =>  [peg1,peg2,...]
#      use      =>  blast|sims            # D = sims
#      expect   => maximum e-value for blast  
#      coverage =>  $min_coverage         # D = 0.70
#      identity =>  $min_identity         # D = 0.60
#
#  OUTPUT: ([id1,id2,...],[id3,id4,...]...)
#

sub partition
{
    my ($options )     = @_;
    my $min_cover      = $options->{ coverage } ||=   0.70;      # Minimum fraction of reference covered
    my $max_exp        = $options->{ expect }   ||=   1.0e-10;   # Maximum e-value for blast
    my $min_identity   = $options->{ expect }   ||=   0.3;       # Maximum e-value for blast
    my $use            = $options->{ use }      ||=   'sims';    # use sims by default

    my($pegs,$connections);
    $connections = {};
    my $fig = new FIG;

    if (($use eq 'sims') && ($pegs = $options->{pegs}))
    {
	$connections = &connect_by_sims($fig,$pegs,$max_exp,$min_identity,$min_cover);
    }
    else
    {
	my @seqs;
	if ($pegs = $options->{pegs})
	{
	    @seqs = grep { length($_->[2]) > 10 } map { [$_,'',$fig->get_translation($_)] } @$pegs;
	}
	else
	{
	    @seqs = @{$options->{seqs}};
	}
	(@seqs > 0) || return ();

	my %length = map { $_->[0] => length($_->[2]) } @seqs;

	my $blastF = "$FIG_Config::temp/tmp.$$.fasta";
	&gjoseqlib::print_alignment_as_fasta($blastF,\@seqs);
	system "$FIG_Config::ext_bin/formatdb -p T -i $blastF";
	open(BLAST,"$FIG_Config::ext_bin/blastall -i $blastF -d $blastF -m 8 -e $max_exp -p blastp |")
	    || die "blast failed";
	while (defined($_ = <BLAST>))
	{
	    my($id1,$id2,$iden,undef,undef,undef,$b1,$e1,$b2,$e2) = split(/\t/,$_);
	    if (($id1 ne $id2) && ($min_identity <= $iden))
	    {
		my $ln1 = $length{$id1};
		my $ln2 = $length{$id2};
		if (((($e1+1-$b1)/$ln1) >= $min_cover) && ((($e2+1-$b2)/$ln2) >= $min_cover))
		{
		    push(@{$connections->{$id1}},$id2);
		    push(@{$connections->{$id2}},$id1);
		}
	    }
	}
	close(BLAST);
	unlink($blastF);
    }
    my @ans = &cluster($connections);
    return wantarray ? @ans : \@ans;
}

sub connect_by_sims {
    my($fig,$pegs,$max_exp,$min_identity,$min_cover) = @_;

    my $conn = {};
    my %pegsH = map { $_ => 1 } @$pegs;
    foreach my $peg (@$pegs)
    {
	foreach my $sim ($fig->sims($peg,10000,$max_exp,'fig'))
	{
	    if (($pegsH{$sim->id2}) &&
		($min_identity <= $sim->iden) &&
		((($sim->e1 + 1 - $sim->b1) / $sim->ln1) >= $min_cover) &&
		((($sim->e2 + 1 - $sim->b2) / $sim->ln2) >= $min_cover))
	    {
		push(@{$conn->{$peg}},$sim->id2);
		push(@{$conn->{$sim->id2}},$peg);
	    }
	}
    }
    return $conn;
}

sub cluster {
    my($conn) = @_;

    my %seen;
    my @clusters = ();
    my @pegs = keys(%$conn);
    my($peg);
    while ($peg = shift @pegs)
    {
	if (! $seen{$peg})
	{
	    my $cluster = [$peg];
	    $seen{$peg} = 1;
	    my $i;
	    for ($i=0; ($i < @$cluster); $i++)
	    {
		my $x = $conn->{$cluster->[$i]};
		foreach my $y (@$x)
		{
		    if (! $seen{$y}) 
		    {
			push(@$cluster,$y);
			$seen{$y} = 1;
		    }
		}
	    }
	    push(@clusters,$cluster);
	}
    }
    return @clusters;
}

1;
