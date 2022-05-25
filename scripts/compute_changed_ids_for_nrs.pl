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


$usage = "usage: compute_changed_ids_for_nrs OldNR OldSyn NewNR AddedIDs ChangedIDs DeletedIDs";

(  ($old_nr  = shift @ARGV)  
&& ($old_syn = shift @ARGV)
&& ($new_nr  = shift @ARGV)  
&& ($added   = shift @ARGV) && open(ADDED,   ">$added") 
&& ($changed = shift @ARGV) && open(CHANGED, ">$changed")
&& ($deleted = shift @ARGV) && open(DELETED, ">$deleted")
)
    || die $usage;

#$old = &load_ids($old_nr);
$new = &load_ids($new_nr);
#$alt_peg = &load_equiv($old_syn);

#while (($key,undef) = each(%$old))

open(OLD, "<", $old_nr) or die "cannot open old nr $old_nr: $!";

while (<OLD>)
{
    next unless /^>(\S+)/;
    my $key = $1;
    
    if (! $new->{$key})
    {
	#
	# The way we build NR these days we will not have any IDs in the
	# NR that would be on the RHS of a peg.synonyms id list.
	#
	print DELETED  "$key\n";
# 	if ($alts = $alt_peg->{$key})
# 	{
# 	    for ($i=0; ($i < @$alts) && (! $new->{$alts->[$i]}); $i++) {}
	    
# 	    if ($i == @$alts)
# 	    {
# 		print DELETED  "$key\n";
# 	    }
# 	    else
# 	    {
# 		delete $new->{$alts->[$i]};
# 		print CHANGED "$key\t$alts->[$i]\n";
# #		print DELETED "$key\n";
# #		print ADDED   "$alts->[$i]\n";
# 	    }
# 	}
# 	else
# 	{
# 	    print DELETED  "$key\n";
# 	}
    }
    else
    {
	delete $new->{$key};
    }
}

while (($key,undef) = each %$new)
{
    print ADDED "$key\n";
}

sub load_ids {
    my($nr) = @_;
    my($entries,$x);

    print "Loading ids\n";
    open(NR, "<", $nr) || die "could not open $nr";
    my $entries = {};
    while (defined($x = <NR>))
    {
	if ($x =~ /^>(\S+)/)
	{
	    $entries->{$1} = 1;
	}
    }
    close(NR);
    print "done\n";
    return $entries;
}

sub load_equiv {
    my($file) = @_;
    my($main,$alt,$main_id,$main_ln,@alt);

    my $alt_peg = {};
    open(TMP,"<$file") || die "could not open $file";
    while (defined($x = <TMP>))
    {
	chop;
	($main,$alt) = split(/\t/,$x);
	($main_id,$main_ln) = split(/,/,$main);
	@alt = map { $_ =~ /^(\S+),(\d+)$/; ($2 == $main_ln) ? $1 : () } split(/;/,$alt);
	if (@alt > 0)
	{
	    $alt_peg->{$main_id} = [@alt];
	}
    }
    close(TMP);
    return $alt_peg;
}
