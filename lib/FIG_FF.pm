use FF;

=head3 get_figfams_data

usage: $dir = $fig->get_figfams_data($mydir)
usage: $dir = &FIG::get_figfams_data($mydir)

Returns the Figfams data directory to use. If $mydir is passed, use that value. Otherwise
see if $FIG_Config::FigfamsData is defined, and use that. Otherwise default to
$FIG_Config::data/FigfamsData.

=cut

sub get_figfams_data
{
    shift if UNIVERSAL::isa($_[0],__PACKAGE__);

    my ($dir) = @_;

    if (defined($dir))
    {
	return $dir;
    }
    if ($FIG_Config::FigfamsData ne '')
    {
	return $FIG_Config::FigfamsData;
    }
    return "$FIG_Config::data/FigfamsData";
}

=head3 dsims

usage: @sims = $fig->dsims($id,$seq,$maxN,$min_nbsc)

Returns a list of similarities for $seq against PEGs from FIGfams such that

    there will be at most $maxN similarities, and

    each similarity will have a normalized bit-score >= $min_nbsc

The "dsims" or "dynamic sims" are not precomputed.  They are computed using a heuristic which
is much faster than blast, but misses some similarities.  Essentially, you have an "index" or
representative sequences, a quick blast is done against it, and if there are any hits these are
used to indicate which sub-databases to blast against.  This implies that the p-scores are fairly
meaningless; use the normalized bit-scores ($sim->nbsc)

=cut

sub dsims {
    my($self,$id,$seq,$maxN,$min_nbsc,$figfams_data,$blast_parms) = @_;
    my($sim,$partition,%hits);

    $figfams_data = $self->get_figfams_data($figfams_data);

    my $reps_db = "$figfams_data/repdb";

    (-s $reps_db) || return ();

    my @index1 = &blastitP('query',$seq,$reps_db,1.0e-5);
    my %indexH;
    foreach $_ (@index1)
    {
	my $id2 = $_->id2;
	if ($id2 =~ /^(FIG\d+)/)
	{
	    my $nbsc = $_->nbsc;
	    if ((! $indexH{$id2}) || ($indexH{$id2} < $nbsc))
	    {
		$indexH{$id2} = $nbsc;
	    }
	}
    }
    my @index = sort {$indexH{$b} <=> $indexH{$a} } keys(%indexH);
    if (@index > 20) { $#index = 19 }

#   print STDERR "index contains ",scalar @index, " entries\n";
    foreach my $id2 (@index)
    {
	if ($id2 =~ /^(FIG\d+)/)
	{
	    my $fam_id = $1;
#	    print STDERR "\ttrying $fam_id\n";
	    my $fam_dir = &FF::fam_dir($figfams_data,$fam_id);
	    if ((-s "$fam_dir/blast.partition") && open(PARTITION,"<$fam_dir/blast.partition"))
	    {
		if (defined($_ = <PARTITION>) && ($_ =~ /^(\d+)/))
		{
		    my $partition = $1;
		    my $partitionF = "$figfams_data/Partitions/" . ($partition % 1000) . "/$partition/fasta";
		    foreach $sim (&blastitP('query',$seq,$partitionF,1.0e-3))
		    {
			if ($sim->nbsc >= $min_nbsc)
			{
			    my $hit = $sim->id2;
			    if ((! $hits{$hit}) || ($sim->nbsc > $hits{$hit}->nbsc))
			    {
				$sim->[0] = $id;
				$hits{$hit} = $sim;
			    }
			}
		    }
		}
		close(PARTITION);
	    }
	}
    }

    return sort { $b->nbsc <=> $a->nbsc } map { $hits{$_} } keys(%hits);
}

