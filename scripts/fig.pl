# -*- perl -*-
########################################################################
# Copyright (c) 2003-2008 University of Chicago and Fellowship
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
########################################################################

use Carp;
use Data::Dumper;
use Time::HiRes qw(gettimeofday);
use Time::Local;

use FIG;
use FIGV;
eval { require FFs; };
eval { require SeedSims; };
eval { require ERDB; };
eval { require FC; };

my($user);

# usage: fig [-echo] [-time] [-url seed-url] [-orgdir RAST-genome-dir] [command]

$echo       = 0;
$time_cmds  = 0;
$orgdir = '';
while ((@ARGV > 0) && ($ARGV[0] =~ /^-/))
{
    $arg = shift @ARGV;
    if ($arg =~ /^-time/i) { $time_cmds = 1 }
    if ($arg =~ /^-echo/i) { $echo      = 1 }
    if ($arg =~ /^-url/i) 
    { 
	$url = shift(@ARGV);
    }
    if ($arg =~ /^-orgdir/)
    {
	$orgdir = shift @ARGV;
    }
}

if ($url)
{
    $fig = new XmlrpcFIGnet;
}
elsif ($orgdir ne '')
{
    $fig = new FIGV($orgdir);
}
else
{
    $fig = new FIG;
}

my($t1,$t2);
if (@ARGV > 0)  { $req = join( " ", @ARGV ); }
while ( (defined($req) && $req) || ((@ARGV == 0) && ($req = &get_req)) )
{
    if ($time_cmds)
    {
	$t1 = gettimeofday;
    }

    if ($req =~ /^\s*h\s*$/ || $req =~ /^\s*help\s*$/)
    {
     &help;
    }
    elsif ($req =~ /^\s*DB\s*$/)
    {
	print "\n$FIG_Config::db\n\n";
    }
    elsif ($req =~ /^\s*translate\s+(\S+)(\s+\S+)?\s*$/)
    {
	$dna    = $1;
	$start  = $2 ? 1 : 0;
	print &FIG::translate($dna,undef,$start),"\n";
    }
    elsif ($req =~ /^\s*reverse_comp\s+(\S+)\s*$/)
    {
	$dna    = $1;
	print &FIG::reverse_comp($dna),"\n";
    }
    elsif ($req =~ /^\s*change_funcrole\s+\'(\S[^\']*\S)\'\s+\'(\S[^\']*\S)\'\s+(\S+)/)
    {
	my $role = $1;
	my $newname = $2;
        my $seed_user = $3;
        $fig->change_funcrole($role,$newname,$seed_user,0);
        print "changed $role\nto $newname\n\n";
    }
    elsif ($req =~ /^\s*standard_genetic_code\s*$/)
    {
	print &Dumper(&FIG::standard_genetic_code),"\n";
    }
    elsif ($req =~ /^\s*cgi_url\s*$/)
    {
	print &FIG::cgi_url,"\n";
    }
    elsif ($req =~ /^\s*temp_url\s*$/)
    {
	print &FIG::temp_url,"\n";
    }
    elsif ($req =~ /^\s*max\s+(((\d+)\s+)*\d+)\s*$/)
    {
	print &FIG::max(split(/\s+/,$1)),"\n\n";
    }
    elsif ($req =~ /^\s*fids_with_link_to\s+(.*)$/)
    {
        print &Dumper([$fig->fids_with_link_to($1)]);
    }
    elsif ($req =~ /^\s*fid_links\s+(fig\|\d+\.\d+\.[^.]+\.\d+)\s*$/)
    {
	print &Dumper([$fig->fid_links($1)]);
    }
    elsif ($req =~ /^\s*add_fid_link\s+(fig\|\d+\.\d+\.[^.]+\.\d+)\s+(\S.*\S)\s*$/)
    {
	$fig->add_fid_links([$1,$2]);
	print "added link to $1\n";
    }
    elsif ($req =~ /^\s*add_fid_links\s+(\S+)\s*$/)
    {
	$file = $1;
	my @links = map { ($_ =~ /^(\S+)\t(\S.*\S)\s*$/) ? [$1,$2] : () } `cat $file`;
	$fig->add_fid_links(@links);
	print "added links\n";
    }
    elsif ($req =~ /^\s*delete_fid_link\s+(fig\|\d+\.\d+\.[^.]+\.\d+)\s+(\S+)\s*$/)
    {
	$fig->delete_fid_link($1,$2);
	print "deleted link for $1\n";
    }
    elsif ($req =~ /^\s*delete_all_fid_links\s+(fig\|\d+\.\d+\.[^.]+\.\d+)\s*$/)
    {
	$fig->delete_all_fid_links($1);
	print "deleted all links for $1\n";
    }
    elsif ($req =~ /^\s*change_location_of_feature\s+(fig\|\S+)\s+(\S+)(\s+(\S+))?\s*$/)
    {
	$fid = $1;
	$loc = $2;
	$tran = $4 ? $4 : "";
	print $fig->change_location_of_feature($fid,$loc,$tran);
    }
    elsif ($req =~ /^\s*add_features\s+(\S+)\s+(\S+)\s*$/)
    {
	$user = $1;
	$file = $2;
	open(FILE,"<$file") || die "could not open $file";
	$n = 0;
	while (defined($_ = <FILE>))
	{
	    chop;
	    ($genome,$contig,$beg,$end,$func,$type,$alias) = split(/\t/,$_);
	    if ($type eq 'peg') { die 'this does not handle pegs properly (the sequence needs to be translated)' }
	    $loc = join("_",($contig,$beg,$end));
	    $seq = $fig->dna_seq($genome,$loc);
	    $alias = $alias ? $alias : "";
	    if ($fid = $fig->add_feature($user,$genome,$type,$loc,$alias,$seq))
	    {
		$n++;
		#  assign_function addes the annotation
		$fig->assign_function($fid,$user,$func);
	    }
	}
	close(FILE);
	print "added $n features\n\n";
    }
    elsif ($req =~ /^\s*add_feature\s+(\S+)\s+(\d+\.\d+)\s+(\S+)\s+(\S+)((\s+(\S+)){0,2})\s*$/)
    {
	$user   = $1;
	$genome = $2;
	$type   = $3;
	$loc    = $4;
	$rest   = $5;
	$aliases = ($rest =~ /aliases=(\S+)/) ? $1 : "";
	$tran    = ($rest =~ /tran=(\S+)/) ? $1 : "";

	if ($fid = $fig->add_feature($user,$genome,$type,$loc,$aliases,$tran))
	{
	    print "$fid\n";
	}
	else
	{
	    print "add_feature failed\n";
	}
    }
    elsif ($req =~ /^\s*delete_feature\s+(\S+)\s+(fig\|\S+)\s*$/)
    {
	$user = $1;
	$fid = $2;
	$fig->delete_feature($user,$fid);
    }
    elsif ($req =~ /^\s*co_occurs\s+(fig\|\S+)\s*$/)
    {
	$fid = $1;
	my $db = ERDB::GetDatabase('Sapling');
	my @fc = &FC::co_occurs($db,$fid);
	print &Dumper(\@fc);
    }
    elsif ($req =~ /^\s*co_occurrence_and_evidence\s+(fig\|\S+)\s*$/)
    {
	$fid = $1;
	my $db = ERDB::GetDatabase('Sapling');
	my @fc = &FC::co_occurrence_and_evidence($db,$fid);
	print &Dumper(\@fc);
    }
    elsif ($req =~ /^\s*co_occurrence_evidence\s+(fig\|\S+)\s+(fig\|\S+)\s*$/)
    {
	my($peg1,$peg2) = ($1,$2);
	my $db = ERDB::GetDatabase('Sapling');
	my $pairs = &FC::co_occurrence_evidence($db,$peg1,$peg2);
	print &Dumper($pairs);
    }
    elsif ($req =~ /^\s*is_co_occurrence_pair\s+(\S+)\s+(\S+)\s*$/)
    {
	my $peg1 = $1;
	my $peg2 = $2;
	my $db  = ERDB::GetDatabase('Sapling');
	my $rc = &FC::is_co_occurrence_pair($db,$peg1,$peg2) ? "yes" : "no";
	print "$rc\n\n";
    }
    elsif ($req =~ /^\s*co_occurrence_set\s+(\S+)\s*$/)
    {
	my $set = $1;
	my $db  = ERDB::GetDatabase('Sapling');
	my($sc,$pair_set) = &FC::co_occurrence_set($db,$set);
	if (! defined($pair_set))
	{
	    print "PairSet $set does not exist\n";
	}
	else
	{
	    print "PairSet $set [score=$sc]\n";
	    foreach $pair (@$pair_set)
	    {
		my($peg1,$peg2) = @$pair;
		my $func1 = $fig->function_of($peg1);
		my $func2 = $fig->function_of($peg2);
		print "\t$peg1\t$func1\t$peg2\t$func2\n";
	    }
	}
	print "\n";
    }
    elsif ($req =~ /^\s*co_occurrence_cluster\s+(\S+)\s*$/)
    {
	my $cluster = $1;
	my $db  = ERDB::GetDatabase('Sapling');
	my @cluster  = &FC::co_occurrence_cluster($db,$cluster);
	print &Dumper(\@cluster);
    }
    elsif ($req =~ /^\s*extend_pairs_and_pair_sets\s+(\S+)\s*$/)
    {
	my $pchF = $1;
	my $db  = ERDB::GetDatabase('Sapling');
	&FC::extend_pairs_and_pair_sets($db,$pchF);
	print "extended Pairings and PairSets using PCHs in $pchF\n";
    }
    elsif ($req =~ /^\s*in_pair_set\s+(\S+)\s+(\S+)\s*$/)
    {
	my $peg1 = $1;
	my $peg2 = $2;
	my $db  = ERDB::GetDatabase('Sapling');
	my($pair_set,$inverted)  = &FC::in_pair_set($db,$peg1,$peg2);
	if ($pair_set)
	{
	    print "pair=$pair_set inverted=$inverted";
	}
	else
	{
	    print "Either pairing does not exist, or it does not connect to a PairSet\n";
	}
    }
    elsif ($req =~ /^\s*all_co_occurrence_pair_sets\s*$/)
    {
	my $db  = ERDB::GetDatabase('Sapling');
	my(@sets) = &FC::all_co_occurrence_pair_sets($db);
	print join(",",@sets),"\n\n";
    }
    elsif ($req =~ /^\s*all_co_occurrence_clusters\s*$/)
    {
	my $db  = ERDB::GetDatabase('Sapling');
	my(@sets) = &FC::all_co_occurrence_clusters($db);
	print join(",",@sets),"\n\n";
    }
    elsif ($req =~ /^\s*cleanup_pair_sets\s*$/)
    {
	my $db  = ERDB::GetDatabase('Sapling');
	&FC::cleanup_pair_sets($db);
	print "Cleaned up PairSets\n\n";
    }
    elsif ($req =~ /^\s*delete_features\s+(\S+)\s+(\S+)\s*$/)
    {
	$user = $1;
	$file = $2;
	my $n = 0;
	if (open(FILE,"<$file"))
	{
	    while (defined($_ = <FILE>))
	    {
		if ($_ =~ /(fig\|\d+\.\d+\.[a-zA-Z_][a-zA-Z_0-9]*\.\d+)/)
		{
		    $fid = $1;
		    $fig->delete_feature($user,$fid);
		    $n++;
		}
	    }
	    close(FILE);
	}
	else
	{
	    print "Cannot open $file: $!\n";
	}
	print "deleted $n features\n";
    }
    elsif ($req =~ /^\s*undelete_feature\s+(\S+)\s+(fig\|\S+)\s*$/)
    {
	$user = $1;
	$fid = $2;
	$fig->undelete_feature($user,$fid);
    }
    elsif ($req =~ /^\s*who_set_function\s+(fig\|\d+\.\d+\.peg\.\d+)\s*$/)
    {
	$peg = $1;
	$who = $fig->who_set_function($peg);
	if ($who)
	{
	    print "$who\n\n";
	}
	else
	{
	    print "We have no record of who set the function\n";
	}
    }
    elsif ($req =~ /^\s*insert_dynamic_sims\s+(\S+)\s*$/)
    {
	$file = $1;
	$rc = 0;
	my @input = map { chomp; [split(/\t/,$_)] } `cat $file`;
	my $rc = $fig->insert_dynamic_sims(\@input);
	print "$rc\n";
    }
    elsif ($req =~ /^\s*recast_ids (\S+)\s+(\S.*\S)/)
    {
	my $pat = $1;
	my @ids = split(/\s+/,$2);
	print join(",",$fig->recast_ids($pat,\@ids)),"\n\n";
    }
    elsif ($req =~ /^\s*add_genome\s+/)
    {
	add_genome_block:
    {
	# I was to to add some options here
	chomp($req);
	my ($genomeF, $force, $skipnr, $dont_mark_complete) = ('', 0, 0, 0);
	my @pieces=split /\s+/, $req;
	shift @pieces;
	my $user = shift @pieces;
	
	if (@pieces == 0)
	{
	    print STDERR "Usage:     add_genome                      User GenomeDir [-force] [-skipnr] [-dont_mark_complete] (-force will ignore verify dir, -skipnr will not add proteins to nr, -dont_mark_complete will keep PROBABLY_COMPLETE from being copied to COMPLETE if assess_completeness ran)\n";
	    last add_genome_block;
	}
		
	foreach  my $test (@pieces) {
	    if ($test eq "-force") {$force=1}
	    elsif ($test eq "-skipnr") {$skipnr=1}
	    elsif ($test eq "-dont_mark_complete") {$dont_mark_complete=1}
	    elsif (-d $test) {$genomeF = $test}
	    else {print STDERR "fig add_genomes: don't understand what $test is\n"}
	}
	if ($genomeF) 
	{ 
	  $genomeF =~ /(\d+\.\d+)/;
	  my $genome=$1;

	  if (-d "$FIG_Config::organisms/$genome")
	  {
	      print "$genome is already installed in this seed ($FIG_Config::organisms/$genome), not adding.\n";
	  }
	  else
	  {
	      my $ok = $fig->add_genome($user,$genomeF, $force, $skipnr, $dont_mark_complete);

	      if ($ok)
	      {
		  print "Added $genome\n";
	      }
	      else
	      {
		  print "Failed to add $genome from $genomeF\n";
	      }
	  }
        }
	else
	{
	 print STDERR "Couldn't find a genome directory in $_\n";
        }
    } }
    elsif ($req =~ /^\s*log_corr\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S.*\S)/)
    {
	$user = $1;
	$genome = $2;
	$corr = $3;
	$msg  = $4;
	$fig->log_corr($user,$genome,$corr,$msg);
    }
    elsif ($req =~ /^\s*replace_features_with\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(\s+(\S+))?/)
    {
	my($old_fid,$user,$genome,$type,$location,$translation) = ($1,$2,$3,$4,$5,$7);
	my $new_fid = $fig->replace_features_with(old_fids => $old_fid,
						  user => $user,
						  genome =>$genome,
						  type => $type,
						  location => $location,
						  translation => $translation);
	print "NEW FID=$new_fid\n\n";
    }
    elsif ($req =~ /^\s*replace_genome\s+/)
    {
        # I was to to add some options here
	chomp($req);
	my ($genomeF, $force, $skipnr)=('', 0, 0);
	my @pieces=split /\s+/, $req;
	shift @pieces;
	my $user    = shift @pieces;
	my $oldG    = shift @pieces;
	my $genomeF = shift @pieces;
	my $mapping = shift @pieces;
	
	if (-s $mapping)
	{
	    $err = 0;
	}
	else
	{
	    print STDERR "Missing a mapping file for the PEG IDs\n";
	    $err = 1;
	}

	foreach  my $test (@pieces) {
	 if ($test eq "-force") {$force=1}
	 elsif ($test eq "-skipnr") {$skipnr=1}
	 else {$err=1; print STDERR "fig replace_genome: don't understand what $test is\n"}
	}

	if (! $err)
	{
	    if ($genomeF) 
	    { 
		$genomeF =~ /(\d+\.\d+)/;
		my $genome=$1;
		if ((! -d "$FIG_Config::organisms/$genome") && ($fig->replace_genome($user,$oldG,$genomeF,$mapping, $force,$skipnr)))
		{
		    print "Replaced $oldG with $genome using mapping file $mapping\n";
		}
		else
		{
		    print "Failed to replace $oldG with $genome from $genomeF\n";
		}
	    }
	    else
	    {
		print STDERR "Couldn't find a genome directory in $_\n";
	    }
	}
    }
    elsif ($req =~ /^\s*delete_genomes\s+(\d.*\d)\s*$/)
    {
	$genomes = $1;
	@genomes = split(/\s+/,$genomes);
	if (@genomes > 0)
	{
	    @bad_genomes = grep { ! -d "$FIG_Config::organisms/$_" } @genomes;
	    if (@bad_genomes > 0)
	    {
		print "Bad Genomes: ",join(" ",@bad_genomes),"\n";
	    }
	    else
	    {
		$fig->delete_genomes(\@genomes);
		print "Completed deleting ",join(" ",@genomes),"\n";
	    }
	}
    }
    elsif ($req =~ /^\s*mark_deleted_genomes\s+(\S+)\s+(\d.*\d)\s*$/)
    {
	$user    = $1;
	$genomes = $2;
	@genomes = split(/\s+/,$genomes);
	if (@genomes > 0)
	{
	    @bad_genomes = grep { ! $fig->is_genome($_) } @genomes;
	    if (@bad_genomes > 0)
	    {
		print "Bad Genomes: ",join(" ",@bad_genomes),"\n";
	    }
	    else
	    {
		$fig->mark_deleted_genomes($user,\@genomes);
		print "Completed marking deleted genomes: ",join(" ",@genomes),"\n";
	    }
	}
    }
    elsif ($req =~ /^\s*unmark_deleted_genomes\s+(\S+)\s+(\d.*\d)\s*$/)
    {
	$user = $1;
	$genomes = $2;
	@genomes = split(/\s+/,$genomes);
	if (@genomes > 0)
	{
	    @bad_genomes = grep { $fig->is_genome($_) } @genomes;
	    if (@bad_genomes > 0)
	    {
		print "Bad Genomes: ",join(" ",@bad_genomes),"\n";
	    }
	    else
	    {
		$fig->unmark_deleted_genomes($user,\@genomes);
		print "Completed unmarking deleted genomes: ",join(" ",@genomes),"\n";
	    }
	}
    }
    elsif ($req =~ /^\s*lock_fid\s+(\S+)\s+(fig\|\d+\.\d+\.[a-zA-Z_]+\.\d+)\s*$/)
    {
	$user = $1;
	$fid = $2;
	$fig->lock_fid($user,$fid);
	print "$user locked $fid\n";
    }
    elsif ($req =~ /^\s*unlock_fid\s+(\S+)\s+(fig\|\d+\.\d+\.[a-zA-Z_]+\.\d+)\s*$/)
    {
	$user = $1;
	$fid = $2;
	$fig->unlock_fid($user,$fid);
	print "$user unlocked $fid\n";
    }
    elsif ($req =~ /^\s*min\s+(((\d+)\s+)*\d+)\s*$/)
    {
	print &FIG::min(split(/\s+/,$1)),"\n\n";
    }
    elsif ($req =~ /^\s*between\s+(\d+)\s+(\d+)\s+(\d+)\s*$/)
    {
	print &FIG::between($1,$2,$3),"\n\n";
    }
    elsif ($req =~ /^\s*is_deleted_fid\s+(\S+)\s*$/)
    {
	$fid = $1;
	print &Dumper($fig->is_deleted_fid($fid));
    }
    elsif ($req =~ /^\s*feature_attributes\s+(\S+)\s*$/)
    {
	$fid = $1;
	print "Attributes of $fid\n";
	foreach $x ($fig->get_attributes($fid))
	{
		print join "\t", @$x, "\n";
	}
    }
    elsif ($req =~ /^\s*get_key_value\s+(.*)\s*$/)
    {
        my $want=$1;
	my ($key, $val)=('','');
	if ($want =~ /key\=(\S+)/i) {$key=$1}
	if ($want =~ /value\=(\S+)/i) {$value=$1}
	print "Attributes with key=$key and value=$value\n";
	foreach $x ($fig->get_key_value($key, $value))
	{
	                print $x, ", ";
        }
    }
    elsif ($req =~ /^\s*erase_attribute_entirely\s+(\S+)\s*$/)
    {
        my $key=$1;
	my $deleted=$fig->erase_attribute_entirely($key);
	print "$deleted files were removed\n";
    }
    elsif ($req =~ /^\s*find_by_attribute\s+(\S+)\s*$/)
    {
    	my $searchTerm=$1;
	my @ret=$fig->find_by_attribute($searchTerm) if ($searchTerm);
	if (@ret) {
	    foreach my $tuple (@ret) {
		print "@$tuple[0]\t@$tuple[1]\t@$tuple[2]\t\n";
	    }
	}
    }

    elsif ($req =~ /\s*add_attribute\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S*)/)
    {
    	my ($fid, $key, $val, $url)=($1, $2, $3, $4);
	my $ret=$fig->add_attribute($fid, $key, $val, $url) if ($fid && $key && $val);
	if ($ret) {print "Added\n\tpeg: $fid\n\tkey: $key\n\tval: $val\n\turl: $url\n"}
	else {print STDERR "Some error adding\n\tpeg: $fid\n\tkey: $key\n\tval: $val\n\turl: $url\nDid not get a value returned from FIG.pm\n"}
    }
    elsif ($req =~ /\s*get_keys\s*/)
    {
      my $keys=$fig->get_keys;
      foreach my $type (keys %$keys)
      {
        foreach my $label (keys %{$keys->{$type}})
        {
          foreach my $peg (@{$keys->{$type}->{$label}})
          {
            print "$type\t$label\t$peg\n";
          }
        }
      }
    } 
    elsif ($req =~ /^\s*get_all_attributes/)
    {
     foreach my $res ($fig->get_attributes)
     {
      print join("\t", @$res), "\n";
     }
    }
    elsif ($req =~ /^\s*get_all_keys/)
    {
     foreach my $res ($fig->get_all_keys)
     {
      print join("\t", @$res), "\n";
     }
    }
    elsif ($req =~ /^\s*merge_related_annotations\s+(\S.*\S)\s*$/)
    {
	$pegs = [split(/\s+/,$1)];
	my @tmp = $fig->merged_related_annotations($pegs);
	print &Dumper(\@tmp);
	print "\n";
    }
    elsif ($req =~ /^\s*related_by_func_sim\s+(\S+)(\s+(\S+))?\s*$/)
    {
	$peg = $1;
	$user = $3 ? $3 : "";
	$func = $fig->function_of($peg,$user);
	print "PEGs related to $peg: $func\n\n";
	foreach $peg ($fig->related_by_func_sim($peg,$user))
	{
	    $func = $fig->function_of($peg1,$user);
	    print "$peg: $func\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*genomes(\s+(\S.*\S))?\s*$/)
    {
	$pat = $2;
	$domain = "";
	if ($pat =~ /^(\S.*\S)\s+([abepv])\s*$/i)
	{
	    $domain = $2;
	    $pat = $1;
	}
	
	my @genomes = ($pat && ($pat =~ /^complete/i)) ? $fig->genomes("complete") : $fig->genomes;
	foreach $genome (@genomes)
	{
	    next if (($domain eq 'a') && (! $fig->is_archaeal($genome)));
	    next if (($domain eq 'b') && (! $fig->is_bacterial($genome)));
	    next if (($domain eq 'p') && (! $fig->is_prokaryotic($genome)));
	    next if (($domain eq 'e') && (! $fig->is_eukaryotic($genome)));
	    next if (($domain eq 'v') && (! $fig->is_viral($genome)));

	    my $genus_species = $fig->genus_species($genome);
	    if ((! $pat) || ($pat =~ /^complete/i) || ($genus_species =~ /$pat/))
	    {
		print join("\t",(&padded($genome,15),$genus_species)),"\n";
	    }
	}
	print "\n";
    }
    elsif ($req =~ /^\s*restricted\s*$/)
    {
	foreach $genome ($fig->genomes)
	{
	    my $genus_species = $fig->genus_species($genome);
	    if (-s "$FIG_Config::organisms/$genome/RESTRICTIONS")
	    {
		print join("\t",(&padded($genome,15),$genus_species)),"\n";
	    }
	}
	print "\n";
    }
    elsif ($req =~ /^\s*genus_species\s+(\d+\.\d+)\s*$/)
    {
	print $fig->genus_species($1),"\n\n";
    }
    elsif ($req =~ /^\s*genome_counts(\s+complete)?\s*$/i)
    {
	($a,$b,$e,$v) = $fig->genome_counts($1);
	print "archaea=$a bacteria=$b eukaryotes=$e viral=$v\n";
    }
    elsif ($req =~ /^\s*file2N\s+(\S+)\s*$/)
    {
	print $fig->file2N($1),"\n\n";
    }
    elsif ($req =~ /^\s*genome_version\s+(\S+)\s*$/)
    {
	$_ = defined($_ = $fig->genome_version($1)) ? $_ : "undefined";
	print "$_\n\n";
    }
    elsif ($req =~ /^\s*crude_estimate_of_distance\s+(\S+)\s+(\S+)\s*$/)
    {
	print $fig->crude_estimate_of_distance($1,$2),"\n\n";
    }
    elsif ($req =~ /^\s*org_of\s+(\S+)\s*$/)
    {
	print $fig->org_of($1),"\n\n";
    }
    elsif ($req =~ /^\s*taxonomy_of\s+(\d+\.\d+)\s*$/)
    {
	print $fig->taxonomy_of($1),"\n\n";
    }
    elsif ($req =~ /^\s*abstract_coupled_to\s+(fig\|\d+\.\d+\.peg\.\d+)/)
    {
	$peg = $1;
	@tmp = $fig->abstract_coupled_to($peg);
	print &Dumper(\@tmp),"\n\n";
    }
    elsif ($req =~ /^\s*coupled_to\s+(fig\|\d+\.\d+\.peg\.\d+)/)
    {
	$peg = $1;
	@tmp = $fig->coupled_to($peg);
	print &Dumper(\@tmp),"\n\n";
    }
    elsif ($req =~ /^\s*coupled_to_batch\s+(\S.*\S)/)
    {
	my @pegs = split(/\s+/,$1);
	my @tmp = $fig->coupled_to_batch(@pegs);
	print &Dumper(\@tmp),"\n\n";
    }
    elsif ($req =~ /^\s*coupling_evidence\s+(fig\|\d+\.\d+\.peg\.\d+)\s+(fig\|\d+\.\d+\.peg\.\d+)\s*$/)
    {
	$peg1 = $1;
	$peg2 = $2;
	@tmp = $fig->coupling_evidence($peg1,$peg2);
	print &Dumper(\@tmp),"\n\n";
    }
    elsif ($req =~ /^\s*coupling_and_evidence\s+(fig\|\d+\.\d+\.peg\.\d+)\s+(\d+)\s+(\S+)\s+(\S+)/)
    {
	$fid = $1;
	$bound = $2;
	$sim_cutoff = $3;
	$coup_cutoff = $4;
	@tmp = $fig->coupling_and_evidence($fid,$bound,$sim_cutoff,$coup_cutoff);
	print &Dumper(\@tmp),"\n\n";
    }
    elsif ($req =~ /^\s*fast_coupling\s+(fig\|\d+\.\d+\.peg\.\d+)\s+(\d+)\s+(\S+)/)
    {
	$fid = $1;
	$bound = $2;
	$coup_cutoff = $3;
	@tmp = $fig->fast_coupling($fid,$bound,$coup_cutoff);
	print &Dumper(\@tmp),"\n\n";
    }
    elsif ($req =~ /^\s*assignments_made\s+(\S+)\s+(\d{1,2}\/\d{1,2}\/\d{4})(\s+(\S.*\S))?\s*$/)
    {
	$who = $1;
	$date = $2;
	if ($3)
	{
	    $genomes = [split(/\s+/,$4)];
	}
	else
	{
	    $genomes = [$fig->genomes];
	}
	    
	foreach $assignment ($fig->assignments_made($genomes,$who,$date))
	{
	    print join("\t",@$assignment),"\n";
	}
    }
    elsif ($req =~ /^\s*auto_assign\s+(fig\|\d+\.\d+\.peg\.\d+)(\s+(\S+))?\s*$/)
    {
	$fid = $1;
	$seq = $3;
	print &FIG::auto_assign($fid,$seq),"\n\n";
    }
    elsif ($req =~ /^\s*auto_assignG\s+(\S+)\s+(\S.*\S)\s*$/)
    {
	$user = $1;
	@genomes = split(/\s+/,$2);
	foreach $genome (@genomes)
	{
	    if (-s "$FIG_Config::organisms/$genome/Features/peg/tbl")
	    {
		print STDERR "making assignments for $genome\n";
		system "cut -f1 $FIG_Config::organisms/$genome/Features/peg/tbl > tmp$$ ; fig auto_assignF $user tmp$$";
		unlink("tmp$$");
	    }
	}
	print "\n";
    }
    elsif (($req =~ /^\s*auto_assignF\s+(\S+)\s+(\S+)\s*$/) && 
	   ($file = $2) && (-s $file) &&
	   ($user = $1))
    {
        $n = 0; 
        foreach $assignment (`cut -f1 $file | auto_assign | make_calls`)
        {
            if ($assignment =~ /^(\S+)\t(\S[^\t]+\S)(\t(\S))?\s*$/)
            {
		$prot = $1;
		$function = $2;
		$conf = $4 ? $4 : "";
		$n++;
		#  Everyone is master, and assign_function writes annotation.
		$fig->assign_function( $prot, $user, $function, $conf );
	    }
	}
	print STDERR "$n assignments made (all attributed to $user)\n";
    }
    elsif ($req =~ /^\s*genes_in_region\s+(\d+\.\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s*$/)
    {
	print &Dumper($fig->genes_in_region($1,$2,$3,$4)),"\n\n";
    }
    elsif ($req =~ /^\s*all_features\s+(\d+\.\d+)\s+(\S+)\s*$/)
    {
	@ids =  sort { &FIG::by_fig_id($a,$b) } $fig->all_features($1,$2);
	while (@ids > 0)
	{
	    @tmp = splice(@ids,0,&FIG::min(3,scalar @ids));
	    print join("\t",@tmp),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*function_of\s+(\S+)\s*(\S+)?\s*$/)
    {
	if ($2)
	{
	    $extra = $2;
	    $all   = ($extra =~ /all/i);
	    $tran  = ($extra =~ /tran/i);
	}
	if ($all)
	{
	    foreach $x ($fig->function_of($1))
	    {
		($who,$func) = @$x;
		if ($tran) { $func = $fig->translate_function($func) }
		print &padded($who,30)," $func\n";
	    }
	}
	else
	{
	    $func = $fig->function_of($1);
	    if ($tran) { $func = $fig->translate_function($func) }
	    print "$func\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*genome_of\s+(\S+)\s*$/)
    {
	print $fig->genome_of($1),"\n\n";
    }
    elsif ($req =~ /^\s*in_cluster_with\s+(\S+)\s*$/)
    {
	$peg = $1;
	@pegs = $fig->in_cluster_with($peg);
	$func = $fig->function_of($peg);
	print "$peg\t$func\n\n";
	foreach $peg (@pegs)
	{
	    $func = $fig->function_of($peg);
	    print "$peg\t$func\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*feature_location\s+(\S+)\s*$/)
    {
	$fid = $1;
	if ($loc = $fig->feature_location($fid))
	{
	    print "$loc\n";
	}
	else
	{
	    print "no location\n";
	}
    }
    elsif ($req =~ /^\s*translation_length\s+(\S+)\s*$/)
    {
	$fid = $1;
	if ($len = $fig->translation_length($fid))
	{
	    print "$len\n";
	}
	else
	{
	    print "no length for $fid\n";
	}
    }
    elsif ($req =~ /^\s*translatable\s+(\S+)\s*$/)
    {
	print $fig->translatable($1),"\n\n";
    }
    elsif ($req =~ /^\s*mapped_prot_ids\s+(\S+)\s*$/)
    {
	foreach $x ($fig->mapped_prot_ids($1))
	{
	    ($id,$len) = @$x;
	    print &padded($id,20),"\t$len\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*get_corresponding_ids\s+(\S+)\s*$/)
    {
	foreach $x ($fig->get_corresponding_ids($1, 1))
	{
	    ($id, $type, $link) = @$x;
	    print &padded($id,20),"\t$type\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*boundaries_of\s+(\S+)\s*$/)
    {
	print join("_",$fig->boundaries_of($1)),"\n\n";
    }
    elsif ($req =~ /^\s*contig_ln\s+(\S+)\s+(\S+)\s*$/)
    {
	print $fig->contig_ln($1,$2),"\n\n";
    }
    elsif ($req =~ /^\s*all_families\s*$/)
    {
	my $ff = new FFs($fig->get_figfams_data());
	my @all = $ff->all_families;
	foreach $fam (@all)
	{
#	    $famO = new FigFam($fig,$fam);
#	    $func = $famO->family_function;
#	    print "$fam\t$func\n";
	    print "$fam\n";
	}
    }
    elsif ($req =~ /^\s*all_contigs\s+(\S+)\s*$/)
    {
	$genome = $1;
	foreach $contig ($fig->all_contigs($genome))
	{
	    $ln = $fig->contig_ln($genome,$contig);
	    print "$contig\t$ln\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*n_contigs\s+(\S+)\s*$/)
    {
        my $genome=$1;
        # this is the old style way of counting contigs
        my $n = scalar $fig->all_contigs($genome);
        print "$n\n";
    }
    elsif ($req =~ /^\s*number_of_contigs\s+(\S+)\s*$/)
    {
        my $genome=$1;
        # this is the old style way of counting contigs
        my $n = $fig->number_of_contigs($genome);
        print "$n\n";
    }
    elsif ($req =~ /^\s*family_function\s+(\S+)\s*$/)
    {
	$family     = $1;
	$func       = $fig->family_function($family);
	print "family $family: $func\n\n";
    }
    elsif ($req =~ /^\s*in_family\s+(\S+)\s*$/)
    {
	$id         = $1;
	@fam        = $fig->in_family($id);
	print &Dumper(\@fam),"\n";
    }
    elsif ($req =~ /^\s*families_for_protein\s+(\S+)\s*$/)
    {
	$id         = $1;
	@fam        = $fig->families_for_protein($id);
	print &Dumper(\@fam),"\n";
    }
    elsif ($req =~ /^\s*proteins_in_family\s+(\S+)\s*$/)
    {
	$fam        = $1;
	@prot       = $fig->proteins_in_family($fam);
	print &Dumper(\@prot),"\n";
    }
    elsif ($req =~ /^\s*external_ids_in_family\s+(\S+)\s*$/)
    {
        my $fam     = $1;
	my @prot    = $fig->ext_ids_in_family($fam);
	print &Dumper(\@prot),"\n";
    }
    elsif ($req =~ /^\s*external_family_for_id\s+(\S+)\s*$/)
    { 
        my $id      = $1;
	my @prot    = $fig->ext_in_family($id);
	print &Dumper(\@prot),"\n";
    }
    elsif ($req =~ /^\s*cid_to_prots\s+(\S+)\s*$/)
    {
	$cid        = $1;
	@prot       = $fig->cid_to_prots($cid);
	print &Dumper(\@prot),"\n";
    }
    elsif ($req =~ /^\s*all_protein_families\s*$/)
    {
	print "Fam.\tSize\n\n";
	foreach $family ($fig->all_protein_families)
	{
	    print "$family\t",$fig->sz_family($family),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*sz_family\s+(\d+)\s*$/)
    {
	$family     = $1;
	$sz         = $fig->sz_family($1);
	print "family $family contains $sz ids\n\n";
    }
    elsif ($req =~ /^\s*ids_in_family\s+(\S+)\s*$/)
    {
	$family     = $1;
	$func       = $fig->family_function($family);
	print "family $family: $func\n\n";
	@ids = $fig->ids_in_family($family);
	while (@ids > 0)
	{
	    @tmp = splice(@ids,0,&FIG::min(3,scalar @ids));
	    print join("\t",@tmp),"\n";
	}
	print "\n";

    }
    elsif ($req =~ /^\s*prot_to_cid\s+(\S+)\s*$/)
    {
	$prot       = $1;
	print $fig->prot_to_cid($prot);
	print "\n\n";

    }
    elsif ($req =~ /^\s*cid_to_prots\s+(\d+)\s*$/)
    {
	$cid        = $1;
	@prots = $fig->cid_to_prots($cid);
	print &Dumper(\@prots),"\n\n";

    }
    elsif ($req =~ /^\s*differentiate_families\s+(\S+)\s+(\S+)/)
    {
      my ($family_id1, $family_id2)=($1, $2);
      my ($fam1, $fam2)=([$fig->ids_in_family($family_id1)], [$fig->ids_in_family($family_id2)]);
      print join("\n", "Proteins that are in $family_id1 and not in $family_id2",
       map {$fig->cid_to_prots($_)} @{&set_utilities::set_diff($fam1, $fam2)});
    }
    elsif ($req =~ /^\s*number_of_families\s+(\S+)/)
    {
        my $source=$1;
        print "There are  ", $fig->number_of_families($source), " families in $source\n";
    }
    elsif ($req =~ /^\s*number_of_proteins_in_families\s+(\S+)/)
    {
        my $source=$1;
        print "There are  ", $fig->number_of_proteins_in_families($source), " proteins in $source\n";
    }
    elsif ($req =~ /^\s*number_of_distinct_proteins_in_families\s+(\S+)/)
    {
        my $source=$1;
        print "There are  ", $fig->number_of_proteins_in_families($source, 1), " proteins in $source\n";
    }
    elsif ($req =~ /^\s*in_sets\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
    {
	$id         = $1;
	$relation   = $2;
	$set_name   = $3;
	print "$id is in set(s) ",join(",",$fig->in_sets($id,$relation,$set_name)),"\n\n";
    }
    elsif ($req =~ /^\s*in_pch_pin_with\s+(\S+)\s*$/)
    {
	$id         = $1;
	print "$id is in set(s) ",join(",",$fig->in_pch_pin_with($id)),"\n\n";
    }
    elsif ($req =~ /^\s*exportable_subsystem\s+(\S.*\S+)\s*$/)
    {
	my $ssa = $1;
	my($spreadsheet,$notes) = $fig->exportable_subsystem($ssa);
	print join("",@$spreadsheet),join("",@$notes),"\n";
    }
    elsif ($req =~ /^\s*is_cluster_based_subsystem\s+(\S.*\S+)\s*$/)
    {
	my $ssa = $1;
	print $fig->is_cluster_based_subsystem($ssa),"\n";
    }
    elsif ($req =~ /^\s*is_exchangable_subsystem\s+(\S.*\S+)\s*$/)
    {
	my $ssa = $1;
	print &FIG::is_exchangable_subsystem($ssa),"\n";
    }
    elsif ($req =~ /^\s*all_exchangable_subsystems\s*$/)
    {
	print join("\n",&FIG::all_exchangable_subsystems),"\n";
    }
    elsif ($req =~ /^\s*export_set\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
    {
	$relation   = $1;
	$set_name   = $2;
	$file       = $3;
	print $fig->export_set($relation,$set_name,$file),"\n\n";
    }
    elsif ($req =~ /^\s*export_chromosomal_clusters\s*$/)
    {
	print $fig->export_set("chromosomal_clusters","cluster_id","$FIG_Config::global/chromosomal_clusters"),"\n\n";
    }
    elsif ($req =~ /^\s*export_pch_pins\s*$/)
    {
	print $fig->export_set("pch_pins","pin","$FIG_Config::global/pch_pins"),"\n\n";
    }
    elsif ($req =~ /^\s*add_chromosomal_clusters\s+(\S+)\s*$/)
    {
	$file = $1;
	print $fig->add_chromosomal_clusters($file),"\n\n";
    }
    elsif ($req =~ /^\s*add_pch_pins\s+(\S+)\s*$/)
    {
	$file = $1;
	print $fig->add_pch_pins($file),"\n\n";
    }
    elsif ($req =~ /^\s*all_sets\s+(\S+)\s+(\S+)\s*$/)
    {
	$relation   = $1;
	$set_name   = $2;
	print "Set\tSize\n\n";
	foreach $set ($fig->all_sets($relation,$set_name))
	{
	    print "$set\t",$fig->sz_set($set,$relation,$set_name),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*sz_set\s+(\d+)\s+(\S+)\s+(\S+)\s*$/)
    {
	$set        = $1;
	$relation   = $2;
	$set_name   = $3;
	$sz         = $fig->sz_set($1,$2,$3);
	print "$set_name $set contains $sz ids\n\n";
    }
    elsif ($req =~ /^\s*next_set\s+(\S+)\s+(\S+)\s*$/)
    {
	$relation   = $1;
	$set_name   = $2;
	print $fig->next_set($relation,$set_name)," is the next $set_name for $set\n\n";
    }
    elsif ($req =~ /^\s*ids_in_set\s+(\d+)\s+(\S+)\s+(\S+)\s*$/)
    {
	$set        = $1;
	$relation   = $2;
	$set_name   = $3;
	print "$set_name $set\n\n";
	@ids = $fig->ids_in_set($set,$relation,$set_name);
	while (@ids > 0)
	{
	    @tmp = splice(@ids,0,&FIG::min(3,scalar @ids));
	    print join("\t",@tmp),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*seqs_with_roles_in_genomes\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
    {
	$genomesF = $1;
	$rolesF   = $2;
	$who = $3;
	$genomes = [map { chop; $_ } `cat $genomesF`];
	$roles   = [map { chop; $_ } `cat $rolesF`];
	print &Dumper($fig->seqs_with_roles_in_genomes($genomes,$roles,$who));
	print "\n";
    }
    elsif ($req =~ /^\s*seqs_with_role\s+(\S+)\s+(\S+)\s+(\S.*\S)\s*$/)
    {
	$genome = $1;
	$role   = $3;
	$who    = $2;
	foreach $peg ($fig->seqs_with_role($role,$who,$genome))
	{
	    $func = $fig->function_of($peg,$who);
	    print "$peg\t$func\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*load_all(\s+(\d+|list))?\s*$/)
    {
	if (&FIG::load_all($2))
	{
	    print "imported everything\n";
	}
    }
    elsif ($req =~ /^\s*load_table\s*(\S+)\s*(\S+)/)
    {
        my $dbf = $fig->{_dbf};
	$file    = $1;
	$tbl     = $2;
	my $result = $dbf->load_table( tbl => $tbl, file => $file );
	print STDERR "Trying to load table '$tbl' and file '$file' got result -->$result<--\n";
    }
    elsif ($req =~ /^\s*dsims\s+(\S+)(\s+(\S+)\s+(\S+)\s+(\S+))?\s*$/)
    {
	$file    = $1;
	$maxN    = defined($3) ? $3 : 500;
	$maxP    = defined($4) ? $4 : 1,0e-5;
	$select  = defined($5) ? $5 : "raw";
	if (open(TMP,"<$file"))
	{
	    $_ = <TMP>;
	    while (defined($_) && ($_ =~ /^>(\S+)/))
	    {
		$id = $1;
		@seqs = ();
		while (defined($_ = <TMP>) && ($_ !~ /^>/))
		{
		    push(@seqs,$_);
		}
		$seq = join("",@seqs);
		$seq =~ s/\s//gs;
		if (@sims = $fig->dsims($id,$seq,$maxN,$maxP,$select))
		{
		    print "$id\t",length($seq),"\n\n";
		    foreach $sim (@sims)
		    {
			$func = $fig->function_of($sim->id2);
			print join("\t",($sim->id2(),
					 $sim->ln2,
					 $sim->psc,
					 $sim->b1,
					 $sim->e1,
					 $sim->b2,
					 $sim->e2,
					 $sim->iden,
					 $sim->bsc,
					 $func)),"\n";
		    }
		    print "\n";
		}
	    }
	    close(TMP);
	}
	print "\n";
    }
    elsif ($req =~ /^\s*access_sims\s+(\S+)?\s*$/)
    {
	$opt = $1 ? $1 : "";
	foreach $peg ("fig|562.1.peg.14","fig|630.1.peg.589","fig|731.2.peg.1626","fig|75985.1.peg.2319","fig|294.1.peg.4630","fig|305.1.peg.2718","fig|292.1.peg.3262","fig|519.1.peg.4061","fig|633.1.peg.2431","fig|384.1.peg.2479")
	{
	    @sims = $fig->sims($peg,600,1.0e-3,"raw");
	    if ($opt =~ /expand/)
	    {
		foreach $sim (@sims)
		{
		    @map = $fig->mapped_prot_ids($sim->id2);
		    foreach $to (@map)
		    {
			if ($opt =~ /func/)
			{
			    @funcs = $fig->function_of($to->[0]);
			}
		    }
		}
	    }
	    elsif ($opt =~ /func/)
	    {
		foreach $sim (@sims)
		{
		    @funcs = $fig->function_of($sim->id2);
		}
	    }
	}
    }
    elsif ($req =~ /^\s*localtime\s+(\S+)\s*$/)
    {
	print "localtime: ",$fig->epoch_to_readable($1),"\n\n";
    }
    elsif ($req =~ /^\s*peg_to_roles_in_subsystems\s+(\S+)\s*$/)
    {
	$peg = $1;
	my @in = $fig->peg_to_roles_in_subsystems($peg);
	print &Dumper(\@in);
	print "\n";
    }
    elsif ($req =~ /^\s*bbhs\s+(\S+)\s+(\S+)\s*$/)
    {
	$id      = $1;
	$cutoff = $2;
	foreach $_ ($fig->bbhs($id,$cutoff))
	{
	    print join("\t",@$_,$fig->org_of($_->[0])),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*partition_coupling_evidence\s+(\d+)\s+(\d+)\s*$/)
    {
	$p1 = $1;
	$p2 = $2;
	my($score,$pairs) = &SeedSims::partition_coupling_evidence($fig,$p1,$p2);
	print "score=$score\n";
	foreach $pair (@$pairs)
	{
	    my($x,$y) = @$pair;
	    print join("\t",($x->[2],scalar $fig->function_of($x->[2]))),"\n";
	    print join("\t",($y->[2],scalar $fig->function_of($y->[2]))),"\n";
	    print "\n";
	}
    }
    elsif ($req =~ /^\s*seed_sims\s+(\S+)(\s+(\d+))?\s*$/)
    {
	$peg      = $1;
	$maxsc    = $3 ? $3 : 500;
	$func = $fig->function_of($peg);
	my @sims = &SeedSims::seed_sims($peg,"-b $maxsc");
#	print STDERR "got sims\n"; @sims = ();
	if (@sims == 0)
	{
	    print "no similarities for $peg\n";
	}
	else
	{
	    print &padded($peg,20),"\t",$sims[0]->ln1."\t$func\n\n";
	    foreach $sim (@sims)
	    {
		$func = $fig->function_of($sim->id2);
		if ($tran) { $func = $fig->translate_function($func) }
	    
		$org  = $fig->org_of($sim->id2);
		print join("\t",(&padded($sim->id2,20),
				 $sim->ln2,
				 $sim->psc,
				 $sim->b1,
				 $sim->e1,
				 $sim->b2,
				 $sim->e2,
				 $sim->bsc,
				 $sim->iden,
				 $org,
				 $func)),"\n";
	    }
	    print "\n";
	}
    }
    elsif ($req =~ /^\s*md5_sims\s+(\S+)\s*$/)
    {
	$md5      = $1;
	$func = $fig->function_of($md5);
	my @sims = &SeedSims::md5_sims($md5);
	print $md5,"\t",$sims[0]->ln1."\t$func\n\n";
	foreach $sim (@sims)
	{
	    $org  = $fig->org_of($sim->id2);
	    print join("\t",($sim->id2,,
			     $sim->ln2,
			     $sim->psc,
			     $sim->b1,
			     $sim->e1,
			     $sim->b2,
			     $sim->e2,
			     $sim->bsc,
			     $sim->iden,
			    )),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*sims\s+(\S+)(\s+(\S+)\s+(\S+)\s+(\S+)(\s+(\d))?)?\s*$/)
    {
	$id      = $1;
	$maxN    = defined($3) ? $3 : 500;
	$maxP    = defined($4) ? $4 : 1.0e-5;
	$select  = defined($5) ? $5 : "raw";
	$tran = defined($7) ? $7 : 0;

	$func = $fig->function_of($id);
	if (@sims = $fig->sims($id,$maxN,$maxP,$select))
	{
	    print &padded($id,20),"\t",$sims[0]->ln1."\t$func\n\n";
	    foreach $sim (@sims)
	    {
		$func = $fig->function_of($sim->id2);
		if ($tran) { $func = $fig->translate_function($func) }

		$org  = $fig->org_of($sim->id2);
		print join("\t",(&padded($sim->id2,20),
				 $sim->ln2,
				 $sim->psc,
				 $sim->b1,
				 $sim->e1,
				 $sim->b2,
				 $sim->e2,
				 $sim->bsc,
				 $sim->iden,
				 $org,
				 $func)),"\n";
	    }
	}
	print "\n";
    }
    elsif ($req =~ /^\s*get_dna_seq\s+(\S+)\s*$/)
    {
	$id   = $1;
	$seq  = $fig->get_dna_seq($id);
	print   $fig->display_id_and_seq($id, \$seq);
    }
    elsif ($req =~ /^\s*get_representative_genome\s+(\S+)\s*$/)
    {
	$id   = $1;
	print $fig->get_representative_genome($id);
	print "\n\n";
    }
    elsif ($req =~ /^\s*get_translation\s+(\S+)\s*$/)
    {
	$id   = $1;
	$seq  = $fig->get_translation($id);
	if ($seq)
	{
	    print   $fig->display_id_and_seq($id, \$seq);
	}
	else
	{
	    print "NO TRANSLATION AVAILABLE -- maybe a deleted peg?\n";
	}
    }
    elsif ($req =~ /^\s*get_translations\s+(\S+)(\s+tab)?\s*$/i)
    {
	$file = $1;
	$tab = $2 ? 1 : 0;
	if (open(TMP,"<$file"))
	{
	    while (defined($_ = <TMP>))
	    {
		foreach $id ($_ =~ /\S+/g)
		{
		    if ($seq = $fig->get_translation($id))
		    {
			if ($tab)
			{
			    print "$id\t$seq\n";			
                        }
			else
			{
			    &FIG::display_id_and_seq($id,\$seq);
			}
		    }
		}
	    }
	    close(TMP);
	}
	else
	{
	    print STDERR "could not open $file\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*blast\s+(\S+)\s*$/)
    {
	$id     = $1;
	$func = $fig->function_of($id);

	if ($sims = &Blast::blastp([[$id,$fig->get_translation($id)]],"$FIG_Config::global/nr","-e 1.0e-5",0,0))
	{
	    $sims = $sims->{$id};

	    print &padded($id,20),"\t",$sims->[0]->ln1."\t$func\n\n";
	    foreach $sim (@$sims)
	    {
		$func = $fig->function_of($sim->id2);
		$org  = $fig->org_of($sim->id2);
		print join("\t",(&padded($sim->id2,20),
				 $sim->ln2,
				 $sim->psc,
				 $sim->b1,
				 $sim->e1,
				 $sim->b2,
				 $sim->e2,
				 $org,
				 $func)),"\n";
	    }
	}
	print "\n";
    }
    elsif ($req =~ /^\s*neighborhood_of_role\s+(\S.*\S)\s*$/)
    {
	$role     = $1;
	@tmp      = $fig->neighborhood_of_role($role);
	print     "$role\n\n",join("\n",@tmp),"\n";
	print "\n";
    }
    elsif ($req =~ /^\s*feature_annotations\s+(\S+)\s*$/)
    {
	$fid = $1;
	print "Annotations of $fid\n\n";
	@annotations = $fig->feature_annotations($fid);
	foreach $x (@annotations)
	{
	    (undef,$ts,$who,$annotation) = @$x;
	    print "$ts\t$who\n$annotation\n============\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*dna_seq\s+(\S+)\s+(\S+)\s*$/)
    {
	$seq = $fig->dna_seq($1,split(/,/,$2));
	print "$seq\n\n";
    }
    elsif ($req =~ /^\s*roles_of_function\s+(\S.*\S)\s*$/)
    {
	$func     = $1;
	@tmp      = &FIG::roles_of_function($func);
	print     "$func\n\n",join("\n",@tmp),"\n";
	print "\n";
    }
    elsif ($req =~ /^\s*assign_function\s+(\S+)\s+(\S+)(\s+conf=(\S))?\s+(\S.*\S)\s*$/)
    {
	$peg      = $1;
	$user     = $2;
	$conf     = $3 ? $4 : "";
	$func     = $5;
	print $fig->assign_function($peg,$user,$func,$conf),"\n\n";
    }
    elsif ($req =~ /^\s*install_assignment_set\s+(\S+)\s+(\S+)\s*$/)
    {
	$user     = $1;
	$file     = $2;
	if (-s $file)
	{
	    &FIG::verify_dir("$FIG_Config::data/Assignments/$user");
	    my $fileA = "$FIG_Config::data/Assignments/$user/" . &FIG::epoch_to_readable(time) . ":$user:imported";
	    system "cp $file $fileA; chmod 777 $fileA";
	}
    }
    elsif ($req =~ /^\s*assign_functionF\s+(\S+)\s+(\S+)(\s+no_annotations)?\s*$/)
    {
	$user     = $1;
	$file     = $2;
	$no_annotations = $3;
	if (open(TMP,"<$file"))
	{
	    while (defined($_ = <TMP>))
	    {
		chop;
		($peg,$func,$conf) = split(/\t/,$_);
		if (! $conf) { $conf = "" }

	        $funcO = $fig->function_of( $peg );
		if ($funcO ne $func)
		{
		    #  assign_function writes the annotation
		    $fig->assign_function( $peg, $user, $func, $conf );
		}
	    }
	    close(TMP);
	}
    }
    elsif ($req =~ /^\s*feature_location\s+(\S+)\s*$/)
    {
	$id     = $1;
	$loc = $fig->feature_location($id);
	print "$loc\n";
	print "\n";
    }
    elsif ($req =~ /^\s*all_compounds\s*$/)
    {
	foreach $cid ($fig->all_compounds)
	{
	    @names = $fig->names_of_compound($cid);
	    print "$cid\t$names[0]\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*names_of_compound\s+(\S+)\s*$/)
    {
	$cid = $1;
	@names = $fig->names_of_compound($cid);
	foreach $name (@names)
	{
	    print "$cid\t$name\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*comp2react\s+(\S+)\s*$/)
    {
	$cid = $1;
	foreach $rid ($fig->comp2react($cid))
	{
	    $x = $fig->displayable_reaction($rid);
	    print "$x\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*reversible\s+(\S+)\s*$/)
    {
	$rid = $1;
	$rev = $fig->reversible($rid) ? "reversible" : "not reversible";;
	$displayable = $fig->displayable_reaction($rid);
	print "$rev\n$displayable\n\n";
    }
    elsif ($req =~ /^\s*catalyzed_by\s+(\S+)\s*$/)
    {
	$rid = $1;
	@ecs = $fig->catalyzed_by($rid);
	$displayable = $fig->displayable_reaction($rid);
	print "$displayable\n";
	foreach $ec (@ecs)
	{
	    $x = $fig->expand_ec($ec);
	    print "$x\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*catalyzes\s+(\S+)\s*$/)
    {
	$ec = $1;
	@rids = $fig->catalyzes($ec);
	print $fig->expand_ec($ec),"\n";
	foreach $rid (@rids)
	{
	    $displayable = $fig->displayable_reaction($rid);
	    print "$displayable\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*reaction2comp\s+(\S+)\s+(\d)\s*$/)
    {
	$rid = $1;
	$which = $2;
	@tmp = $fig->reaction2comp($rid,$which);
	print &Dumper(\@tmp);
	$displayable = $fig->displayable_reaction($rid);
	print "$displayable\n\n";
    }
    elsif ($req =~ /^\s*cas\s+(\S+)\s*$/)
    {
	$cid = $1;
	$cas = $fig->cas($cid);
	@names = $fig->names_of_compound($cid);
	print "$cid\t$cas\t$names[0]\n\n";
    }
    elsif ($req =~ /^\s*cas_to_cid\s+(\S+)\s*$/)
    {
	$cas = $1;
	$cid = $fig->cas_to_cid($cas);
	@names = $fig->names_of_compound($cid);
	print "$cid\t$cas\t$names[0]\n\n";
    }
    elsif ($req =~ /^\s*all_reactions\s*$/)
    {
	foreach $rid ($fig->all_reactions)
	{
	    print $fig->displayable_reaction($rid),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*ec_to_maps\s+(\S+)\s*$/)
    {
	print join("\n",$fig->ec_to_maps($1)),"\n\n";
    }
    elsif ($req =~ /^\s*map_to_ecs\s+(\S+)\s*$/)
    {
	print join("\n",$fig->map_to_ecs($1)),"\n\n";
    }
    elsif ($req =~ /^\s*map_name\s+(\S+)\s*$/)
    {
	print $fig->map_name($1),"\n\n";
    }
    elsif ($req =~ /^\s*add_annotation\s+(\S+)\s+(\S+)\s*$/)
    {
	$fid     = $1;
	$user    = $2;
	print "Type in annotation (end with \".\" at start of line)\n>> ";
	@ann = ();
	while (defined($_ = <STDIN>) && ($_ !~ /^\./))
	{
	    push(@ann,$_);
	    print ">> ";
	}
	print $fig->add_annotation($fid,$user,join("",@ann)),"\n\n";
    }
    elsif ($req =~ /^\s*add_annotations\s+(\S+)\s*$/)
    {
	my $file = $1;

	my $dbh = $fig->db_handle()->{_dbh};
	local $dbh{RaiseError} = 1;
	my($count, $errors) = $fig->add_annotation_batch($file);
	
	print "Added $count annotations\n";
	if (@$errors > 0)
	{
	    print Dumper($errors);
	}
    }
    elsif ($req =~ /^\s*h\s+(\S+)/)
    {
	if (-s "Help/$1")
	{
	    print "Sorry, we have not implemented the detailed help feature yet\n";
#	    @tmp = `cat Help/$1`;
#	    print "\n",join("",@tmp),"\n";
	}
	else
	{
	    print "sorry, no help for $1\n\n";
	}
    }
    elsif ($req =~ /^\s*verify_dir\s+(\S+)\s*$/)
    {
	$dir = $1;
	&FIG::verify_dir($dir);
	print "$dir is now a directory\n";
    }
    elsif ($req =~ /^\s*ec_name\s+(\S+)\s*$/)
    {
	$ec= $1;
        print $fig->ec_name($ec),"\n";
    }
    elsif ($req =~ /^\s*all_roles\s*$/)
    {
	my @tmp = $fig->all_roles;
	print &Dumper(\@tmp);
    }

    elsif ($req =~ /^\s*expand_ec\s+(\S+)\s*$/)
    {
	$ec = $1;
        print $fig->expand_ec($ec),"\n";
    }
    elsif ($req =~ /^\s*clean_tmp\s*$/)
    {
	&FIG::clean_tmp;
	print "Cleaned $FIG_Config::temp\n";
    }
    elsif ($req =~ /^\s*org_of\s+(\S+)\s*$/)
    {
        $prot_id = $1;
	print $fig->org_of($prot_id),"\n";
    }
    elsif ($req =~ /^\s*abbrev\s+(\S.*\S)\s*$/)
    {
        $genome_name = $1;
	print &FIG::abbrev($genome_name),"\n";;
    }
    elsif ($req =~ /^\s*ftype\s+(\S+)\s*$/)
    {
        $feature_id = $1;
	print &FIG::ftype($feature_id),"\n";;
    }
    elsif ($req =~ /^\s*genome_of\s+(\S+)\s*$/)
    {
        $feature_id = $1;
	print &FIG::genome_of($feature_id),"\n";
    }
    elsif ($req =~ /^\s*by_fig_id\s+(\S+)\s+(\S+)\s*$/)
    {
        $feature_id_a = $1;
        $feature_id_b = $2;
	print &FIG::by_fig_id($feature_id_a,$feature_id_b),"\n";
    }
    elsif ($req =~ /^\s*close_genes\s+(\S+)\s+(\S+)\s*$/)
    {
        $feature_id = $1;
        $dist = $2;
	print join(",",$fig->close_genes($feature_id,$dist)),"\n";
    }
    elsif ($req =~ /^\s*pegs_of\s+(\S+)\s*$/)
    {
        $genome = $1;
	print join(",",$fig->pegs_of($genome)),"\n";
    }
    elsif ($req =~ /^\s*rnas_of\s+(\S+)\s*$/)
    {
        $genome = $1;
	my @tmp = $fig->rnas_of($genome);
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*subsystems_for_genome\s+(\S+.*)$/)
    {
        my ($genome, $zero) = split /\s+/, $1;
	my @tmp = $fig->subsystems_for_genome($genome, $zero);
	print &Dumper(@tmp);
    }
    elsif ($req =~ /^\s*subsystems_for_peg\s+(\S+)\s*$/)
    {
        my $peg = $1;
	my @tmp = $fig->subsystems_for_peg($peg);
#	print "Getting ss for '$peg'\n";
	print &Dumper(@tmp);
    }
    elsif ($req =~ /\s*subsystem_genomes\s+(\S+)\s*$/)
    {
        my $ss=$1;
        my @tmp = $fig->subsystem_genomes($ss);
        print &Dumper(@tmp);
    }
    elsif ($req =~ /\s*subsystem_classification\s+(\S+)\s*$/)
    {
        my $ss=$1;
        my $tmp = $fig->subsystem_classification($ss);
	if ($tmp)
	{
	    my($a, $b) = @$tmp;
	    print qq("$a" "$b"\n);
	}
	else
	{
	    print "No classification for $ss\n";
	}
    }
    elsif ($req =~ /\s*all_subsystem_classifications\s*$/)
    {
	my @tmp = $fig->all_subsystem_classifications();
	print Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*get_titles_by_gi\s+(\S+)\s*$/)
    {
	print &Dumper([$fig->get_titles_by_gi($1)]);
    }
    elsif ($req =~ /^\s*get_titles_by_peg\s+(\S+)\s*$/)
    {
	print &Dumper([$fig->get_titles_by_peg($1)]);
    }
    elsif ($req =~ /^\s*get_title_by_pmid\s+(\S+)\s*$/)
    {
	print &Dumper($fig->get_title_by_pmid($1));
    }
    elsif ($req =~ /^\s*feature_aliases\s+(\S+)\s*$/)
    {
        $feature_id = $1;
	print join(",",$fig->feature_aliases($feature_id)),"\n";
    }
    elsif ($req =~ /^\s*to_alias\s+(\S+)\s+(\S+)\s*$/)
    {
	print &Dumper($fig->to_alias($1,$2));
    }
    elsif ($req =~ /^\s*by_alias\s+(\S+)\s*$/)
    {
        $alias = $1;
	print &Dumper($fig->by_alias($alias));
    }
    elsif ($req =~ /^\s*possibly_truncated\s+(\S+)\s*$/)
    {
        $feature_id = $1;
	print $fig->possibly_truncated($feature_id),"\n";
    }
    elsif ($req =~ /^\s*is_real_feature\s+(\S+)\s*$/)
    {
        $feature_id = $1;
	print $fig->is_real_feature($feature_id),"\n";
    }
    elsif ($req =~ /^\s*is_locked_fid\s+(\S+)\s*$/)
    {
        $feature_id = $1;
	print $fig->is_locked_fid($feature_id),"\n";
    }
    elsif ($req =~ /^\s*translatable\s+(\S+)\s*$/)
    {
        $prot_id = $1;
	print &FIG::translatable($prot_id),"\n";
    }
    elsif ($req =~ /^\s*maps_to_id\s+(\S+)\s*$/)
    {
        $id = $1;
	print $fig->maps_to_id($id),"\n";
    }
    elsif ($req =~ /^\s*translated_function_of\s+(\S+)\s+(\S+)\s*$/)
    {
        $peg = $1;
        $user = $2;
	print $fig->translated_function_of($peg,$user),"\n";
    }
    elsif ($req =~ /^\s*translate_function\s+(\S.*\S)\s*$/)
    {
        $func = $1;
	print $fig->translate_function($func),"\n";
    }
    elsif ($req =~ /^\s*hypo\s+(\S.*\S)\s*$/)
    {
	print &FIG::hypo($1),"\n";
    }
    elsif ($req =~ /^\s*blastit\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/)
    {
        $peg = $1;
        $seq = $2;
        $db = $3;
        $maxP = $4;
	@tmp = &FIG::blastit($peg,$seq,$db,$maxP);
	print &Dumper(\@tmp);
    }    
    elsif ($req =~ /^\s*related_by_func_sim\s+(\S+)\s+(\S+)\s*$/)
    {
        $peg = $1;
        $user = $2;
        print join(",",$fig->related_by_func_sim($peg,$user)),"\n";
    }
    elsif ($req =~ /^\s*epoch_to_readable(\s+(\d+))?\s*$/)
    {
	if (defined($2))
	{
	    print &FIG::epoch_to_readable($2);
	}
	else
	{
	    print &FIG::epoch_to_readable(time);
	}
    }
    elsif ($req =~ /^\s*epoch_date\s+(\d{1,2}\/\d{1,2}\/\d{4})\s*$/)
    {
	print $fig->parse_date($1),"\n\n";
    }
    elsif ($req =~ /^\s*search_index\s+(\S(.*\S)?)\s*$/)
    {
        $pattern = $1;
	my @tmp = $fig->search_index($pattern);
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*init_local_sim_computation\s*$/)
    {
	print "Creating sim askfor pool... (this may take a few minutes)\n";
	my $id = $fig->create_sim_askfor_pool();
	print "Created sim askfor pool id=$id\n";
	print "Run compute_sim_chunk to kick off the actual computation\n";
    }
    elsif ($req =~ /^\s*sz_family\s+(\S+)\s*$/)
    {
        $family = $1;
	print $fig->sz_family($family),"\n";
    }
    elsif ($req =~ /^\names_of_compound\s+(\S+)\s*$/)
    {
	my @tmp = $fig->names_of_compound($1);
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*reversible\s+(\S+)\s*$/)
    {
        $rid = $1;
	print $fig->reversible($rid),"\n";
    }
    elsif ($req =~ /^\s*reaction2comp\s+(\S+)\s+(\S+)\s*$/)
    {
        $rid = $1;
        $which = $2;
	my @tmp = $fig->reaction2comp($rid,$which);
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*displayable_reaction\s+(\S+)\s*$/)
    {
        $rid = $1;
	print $fig->displayable_reaction($rid),"\n";
    }
    elsif ($req =~ /^\s*all_maps\s*$/)
    {
       	my @tmp = $fig->all_maps;
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*neighborhood_of_role\s+(\S+)\s*$/)
    {
        $role = $1;
	my @tmp = $fig->neighborhood_of_role($role);
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*roles_of_function\s+(\S+)\s*$/)
    {
        $func = $1;
	my @tmp = $fig->roles_of_function($func);
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*largest_clusters\s+(\S+)\s+(\S+)\s*$/)
    {
        $file = $1;
        $user = $2;
	my @tmp = `cat $file`;
	chop @tmp;;
	my @tmp1 = $fig->largest_clusters(\@tmp,$user);
	print &Dumper(\@tmp1);
    }
    elsif ($req =~ /^\s*largest_co_occurrence_clusters\s+(\S+)\s*$/)
    {
        my $peg = $1;
	my $db = ERDB::GetDatabase('Sapling');
	my @tmp1 = &FC::largest_co_occurrence_clusters($db,$peg);
	print &Dumper(\@tmp1);
    }
    elsif ($req =~ /^\s*in_co_occurrence_cluster\s+(\S+)\s*$/)
    {
        my $peg = $1;
	my $db = ERDB::GetDatabase('Sapling');
	my @tmp1 = &FC::in_co_occurrence_cluster($db,$peg);
	print &Dumper(\@tmp1);
    }
    elsif ($req =~ /^\s*co_occurring_FIGfams\s+(\S+)\s*$/)
    {
        my $ff = $1;
	my $db = ERDB::GetDatabase('Sapling');
	my @tmp1 = &FC::co_occurring_FIGfams($db,$ff);
	print &Dumper(\@tmp1);
    }
    elsif ($req =~ /^\s*unique_functions\s+(\S+)\s+(\S+)\s*$/)
    {
        $pegs = $1;
        $user = $2;
	my @tmp = $fig->unique_functions([split(/,/,$pegs)],$user);
	print &Dumper(\@tmp);
    }
    elsif ($req =~ /^\s*candidates_for_role\s+(\d+\.\d+)\s+(\S+)\s+(\S.*\S)\s*$/)
    {
	$genome = $1;
	$cutoff = $2;
	$role   = $3;
	my @tmp = $fig->candidates_for_role($role,$genome,$cutoff);
	foreach $peg (@tmp)
	{
	    print "$peg\t",scalar $fig->function_of($peg),"\n";
	}
	print "\n";
    }
    elsif ($req =~ /^\s*external_calls\s+(\S.*\S)\s*$/)
    {
	my @pegs = split(/[, \t]+/,$1);
	print join("\n",map { join("\t",@$_) } $fig->external_calls(\@pegs)),"\n";
	print "\n";
    }
    elsif ($req =~ /^\s*same_func\s+\'([^\']+)\'\s+\'([^\']+)\'\s*$/)
    {
	print $fig->same_func($1,$2),"\n";
	print "\n";
    }
    elsif ($req =~ /^\s*screwed_up\s+(fig\|\d+\.\d+\.peg\.\d+)/)
    {
	print &Dumper($fig->screwed_up($1));
	print "\n";
    }
    elsif ($req =~ /^\s*best_bbh_candidates\s+(\d+\.\d+)\s+(\S+)\s+(\d+)\s+(\S.*\S)\s*$/)
    {
	$genome = $1;
	$cutoff = $2;
	$requested = $3;
	$known     = [split(/[, \t]+/,$4)];
	
	my @tmp = $fig->best_bbh_candidates($genome,$cutoff,$requested,$known);
	foreach $_ (@tmp)
	{
	    print &Dumper($_);
	}
	print "\n";
    }
    elsif ($req =~ /^\s*extract_seq\s+(\S+)\s+(\S+)\s*$/)
    {
        $contigsF = $1;
        $loc = $2;
	$contigs  = &load_contigs($contigsF);
	print &Dumper(&FIG::extract_seq($contigs,$loc));
    }
    elsif ($req =~ /^\s*is_bacterial\s+(\S+)\s*$/)
    {
        $genome = $1;
        print $fig->is_bacterial($genome),"\n";;
    }
    elsif ($req =~ /^\s*is_complete\s+(\S+)\s*$/)
    {
        $genome = $1;
        print $fig->is_complete($genome),"\n";;
    }
    elsif ($req =~ /^\s*is_archaeal\s+(\S+)\s*$/)
    {
        $genome = $1;
        print $fig->is_archaeal($genome),"\n";
    }
    elsif ($req =~ /^\s*is_aux_role_in_subsystem\s+\'([^\']+)\'\s+\'([^\']+)\'\s*$/)
    {
        print $fig->is_aux_role_in_subsystem($1,$2),"\n";
    }
    elsif ($req =~ /^\s*is_prokaryotic\s+(\S+)\s*$/)
    {
        $genome = $1;
        print $fig->is_prokaryotic($genome),"\n";
    }
    elsif ($req =~ /^\s*is_eukaryotic\s+(\S+)\s*$/)
    {
        $genome = $1;
        print $fig->is_eukaryotic($genome),"\n";
    }
    elsif ($req =~ /^\s*sort_genomes_by_taxonomy\s+(\S+(\s+\S+)*)\s*$/)
    {
        @list_of_genomes = split(/\s+/,$1);
        print join(",",$fig->sort_genomes_by_taxonomy(@list_of_genomes)),"\n";
    }
    elsif ($req =~ /^\s*sort_fids_by_taxonomy\s+(\S+(\s+\S+)*)\s*$/)
    {
        @fids = split(/\s+/,$1);
        print join(",",$fig->sort_fids_by_taxonomy(@fids)),"\n";
    }
    elsif ($req =~ /^\s*build_tree_of_complete\s+(\S+)\s*$/)
    {
        $min_for_label = $1;
        print &Dumper($fig->build_tree_of_complete($min_for_label));
    }
    elsif ($req =~ /^\s*taxonomic_groups_of_complete\s+(\S+)\s*$/)
    {
        $min_for_labels = $1;
        print &Dumper($fig->taxonomic_groups_of_complete($min_for_labels));
    }
    elsif ($req =~ /^\s*update_subsys_conflicts\s*$/)
    {
	&FIG::run("pegs_in_conflict | peg_to_subsystems > $FIG_Config::global/conflicted.pegs");
	print "Updated record of subsystem conflicts\n\n";
    }
    elsif ($req =~ /^\s*pegs_in_subsystems\s+(\S+)/)
    {
        my $genome=$1;
	print map  {my $peg=$_; map {"$peg\t".$_->[0]."\n"} $fig->subsystems_for_peg($peg)} $fig->pegs_of($genome);
    }
    elsif ($req =~ /^\s*pegs_in_proteinfams\s+(\S+)/)
    {
        my $genome=$1;
	print map  {my $peg=$_; map {"$peg\t".$_."\n"} $fig->families_for_protein($peg)} $fig->pegs_of($genome);
    }

    elsif ($req =~ /^\s*pegs_in_proteinfams_by_homology/)
    {
	chomp $req; s/^\s+//;
	my ($request, $genome, $maxN, $maxP)=(split /\s+/, $req);
	$maxN or ($maxN=5);
	$maxP or ($maxP=1e-20);
	foreach my $peg ($fig->pegs_of($genome)) {
	  # I prefer this way over a another map layer since it appears to use less memory
	  print 
	     map{
	       my $sims=$_;
	       map {"$peg\t".$sims->[1]."\t".$_."\n"} $fig->families_for_protein($sims->[1]);
	     } $fig->sims($peg, $maxN, $maxP, 'fig', 'figx'); # figx should force maximum expansion 
	}
    }
    elsif ($req =~ /^\s*md5_of_peg\s+(\S+)/)
    {
	my $peg = $1;
	print $fig->md5_of_peg($peg),"\n";
    }
    elsif ($req =~ /^\s*pegs_with_md5\s+(\S+)/)
    {
	my $md5 = $1;
	print join(",",$fig->pegs_with_md5($md5)),"\n\n";
    }
    elsif ($req =~ /^\s*is_NMPDR_genome\s+(\S+)\s*$/)
    {
	$is = $fig->is_NMPDR_genome($1);
	print "is=$is\n";
    }
    elsif ($req =~ /^\s*fr_to_go\s+(\S.*\S)\s*$/)
    {
	$role = $1;
	my @roles = $fig->fr_to_go($role);
	print &Dumper(\@roles);
    }
    elsif ($req =~ /^\s*add_dlit\s+\'([ DRGN])\'\s+(fig\|\d+\.\d+\.peg\.\d+)\s+(\d+)\s+(\S+)\s+\'([A-Z]+(,[A-Z]+)*)?\'\s+([01])\s*$/)
    {
	my($status,$peg,$pubmed,$curator,$go,$override) = ($1,$2,$3,$4,$5,$7);
	my $rc = $fig->add_dlit(-status    => $status,
				-peg       => $peg,
				-pubmed    => $pubmed,
				-curator   => $curator,
				-go        => $go,
				-override  => $override);
	print "RC=$rc\n";
    }
    elsif ($req =~ /^\s*add_dlits\s+(\S+)(.*)\s*$/)
    {
	my $file = $1;
	my $args = $2;
	my $status   = ($args =~ /-s\s*\'([ NGR])\'/)   ? $1 : "D";
	my $go       = ($args =~ /-g\s*\'([A-Za-z])\'/) ? $1 : "";
	my $curator  = ($args =~ /-c\s*([A-Za-z]+)/)    ? $1 : "";
	my $override = ($args =~ /-o\s*([01])/)         ? $1 : 0;

	open(DLITS,"<$file") || die "could not open $file";
	while (defined(my $x = <DLITS>))
	{
	    chomp $x;
	    my $rc = 0;
	    my @flds = split(/\t/,$x);
	    if (&FIG::between(2,scalar @flds,4) && 
		(($flds[0] =~ /^fig\|\d+\.\d+\.peg\.\d+/) || (length($flds[0]) == 32)) && ($flds[1] =~ /^\d+$/))
	    {
		my $curator1 = ((@flds >= 3) && ($flds[2] !~ /\s/)) ? $flds[2] : $curator;
		my $go1      = $go;
		my $peg;
		if ($flds[0] =~ /^fig\|/) 
		{
		    $peg      = $flds[0];
		}
		else
		{
		    my @pegs = $fig->pegs_with_md5($flds[0]);
		    $peg = $pegs[0];
		}
		my $pubmed   = $flds[1];
		if ((@flds == 3) && ($flds[2] =~ /\s/))
		{
		    $fig->add_title($pubmed,$flds[2]);
		}
		elsif ((@flds == 4) && ($flds[3] =~ /\s/))
		{
		    $fig->add_title($pubmed,$flds[3]);
		}
		$rc = $fig->add_dlit(-status    => $status,
				     -peg       => $peg,
				     -pubmed    => $pubmed,
				     -curator   => $curator1,
				     -go        => $go1,
				     -override  => $override);
	    }
	    if (&FIG::between(3,scalar @flds,5) && 
		(($flds[1] =~ /^fig\|\d+\.\d+\.peg\.\d+/) || (length($flds[1]) == 32)) && ($flds[2] =~ /^\d+$/))
	    {
		my $curator1 = ((@flds >= 4) && ($flds[3] !~ /\s/)) ? $flds[3] : $curator;
		my $go1      = $go;
		my $peg      = $flds[1];
		my $pubmed   = $flds[2];
		my $status1  = $flds[0];

		my $peg;
		if ($flds[1] =~ /^fig\|/) 
		{
		    $peg      = $flds[1];
		}
		else
		{
		    my @pegs = $fig->pegs_with_md5($flds[1]);
		    $peg = $pegs[0];
		}

		if ((@flds == 4) && ($flds[3] =~ /\s/))
		{
		    $fig->add_title($pubmed,$flds[3]);
		}
		elsif ((@flds == 5) && ($flds[4] =~ /\s/))
		{
		    $fig->add_title($pubmed,$flds[4]);
		}
		$rc = $fig->add_dlit(-status    => $status1,
				     -peg       => $peg,
				     -pubmed    => $pubmed,
				     -curator   => $curator1,
				     -go        => $go1,
				     -override  => $override);
	    }
	    elsif (@flds == 5)
	    {
		my $peg;
		my($status,$peg_or_hash,$pubmed,$curator,$go) = @flds;
		if ($peg_or_hash !~ /^fig\|/)
		{
		    my @pegs = $fig->pegs_with_md5($peg_or_hash);
		    $peg = $pegs[0];
		}
		else
		{
		    $peg = $peg_or_hash;
		}

		$rc = $fig->add_dlit(-status    => $status,
				     -peg       => $peg,
				     -pubmed    => $pubmed,
				     -curator   => $curator,
				     -go        => $go,
				     -override  => $override);
	    }
	    if (! $rc) { print STDERR "$x\n" }
	}
	print "loaded dlits from $file\n\n";
    }
    elsif ($req =~ /^\s*go_number_to_term\s+(\S+)\s*$/)
    {
	print $fig->go_number_to_term($1),"\n";
    }
    elsif ($req =~ /^\s*go_number_to_info\s+(\S+)\s*$/)
    {
	print &Dumper($fig->go_number_to_info($1));
    }
    
    elsif ($req =~ /^\s*function_to_subsystems\s+(\S.*\S)/) {
	    my $fn = $1;
	    $fn =~ s/^\'//; $fn =~ s/^\'//;
	    $fn =~ s/^\"//; $fn =~ s/^\"//;
	    foreach $sub (sort $fig->function_to_subsystems($fn)) 
	    {
		    print "$sub\n";
	    }
	    print "\n\n";
    }
    elsif ($req =~ /^\s*role_to_subsystems\s+(\S.*\S)/)
    {
	foreach $sub (sort $fig->role_to_subsystems($1))
	{
	    print "$sub\n";
	}
	print "\n\n";
    }
    elsif ($req =~ /^\s*prots_for_role\s+(\S.*\S)/)
    {
	my @pegs = $fig->prots_for_role($1);
	print join(",",@pegs),"\n\n";
    }
    elsif ($req =~ /^\s*role_to_pegs\s+(\S.*\S)/)
    {
	my @pegs = $fig->role_to_pegs($1);
	print join(",",@pegs),"\n\n";
    }
    elsif ($req =~ /^\s*variant_code\s+\'(.*)\'\s+(\d+\.\d+)\s*/)
    {
	print $fig->variant_code($1,$2),"\n";
    }
    elsif ($req =~ /^\s*pegs_in_subsystems_by_homology/)
    {
	chomp $req; s/^\s+//;
	my ($request, $genome, $maxN, $maxP)=(split /\s+/, $req);
	$maxN or ($maxN=5);
	$maxP or ($maxP=1e-20);
	print "Looking for pegs from $genome using a cutoff of $maxP and looking through $maxN sims\n";
	foreach my $peg ($fig->pegs_of($genome)) {
	   print 
	     map{
	       my $sims=$_;
	       map {"$peg\t".$sims->[1]."\t".$_->[0]."\n"} $fig->subsystems_for_peg($sims->[1]);
	     } $fig->sims($peg, $maxN, $maxP, 'fig', 'figx'); # figx should force maximum expansion
        }
    }
    elsif ($req =~ /^\s*change_func_roles\s+(\S+)/)
    {
	my $file = $1;
	open (NEW_ROLES, "<$file") or die "cannot open $file\n";
	
	while (my $line = <NEW_ROLES>)
	{
	    chomp ($line);
	    my ($seeduser, $synFlag, $oldrole, $newrole) = split (/\t/, $line);
	    $fig->change_funcrole($oldrole,$newrole,$seeduser, $synFlag);
	    print STDERR "changed functional role $oldrole to $newrole\n";
	}
	close NEW_ROLES;
    }
    elsif ($req eq "phages") {
    	eval {
	    require Phage;
	    my $p=new Phage;
	    map {print "$_\t", $fig->genus_species($_), "\n"} $p->phages;
	};
    }
    else
    {
	print "invalid command\n";
    }
    print "\n";
    $req = "";

    if ($time_cmds)
    {
	$t2 = gettimeofday;
	print $t2-$t1," seconds to execute command\n\n";
    }
}
sub padded {
    my($x,$n) = @_;

    if (length($x) < $n)
    {
	return $x . (" " x ($n - length($x)));
    }
    return $x;
}

sub get_req {
    my($x);

    print "?? ";
    $x = <STDIN>;
    while (defined($x) && ($x =~ /^h$/i) )
    { 
	&help;
	print "?? ";
	$x = <STDIN>;
    }
    
    if ((! defined($x)) || ($x =~ /^\s*[qQxX]/))
    {
	return "";
    }
    else
    {
        if ($echo)
	{
	    print ">> $x\n";
	}
	return $x;
    }
}

sub load_contigs {
    my($file) = @_;
    my($id,$seq);

    my($contigs) = {};
    open(CONTIGS,"<$file") || die "could not open $file";
    $/ = "\n>";
    while (defined($_ = <CONTIGS>))
    {
	chomp;
	if ($_ =~ /^>?(\S+)[^\n]*\n(.*)/s)
	{
	    $id  =  $1;
	    $seq =  $2;
	    $seq =~ s/\s//g;
	    $contigs->{$id} = $seq;
	}
    }
    close(CONTIGS);
    $/ = "\n";
    return $contigs;
};


sub help {
    print <<END;
    DB                              Current DB
    abbrev                          genome_name
    abstract_coupled_to             PEG
    access_sims                     [expand|func|expandfunc]       (for timeing - access 10 distinct sims)
    add_annotation                  FID User [ prompted for annotation; terminate with "." at start of line ]
    add_annotations                 File
    add_attribute		    FID  key  value    url
    add_chromosomal_clusters        File
    add_dlit                        'status' PEG PUBMED CURATOR 'GO' OVERRIDE
    add_dlits                       File [-s 'Status'] [-c Curator] [-g 'GO'] [-o OVERRIDE]
    add_feature                     User Genome Type Location [aliases=Aliases] [tran=Translation]
    add_features                    User File [Genome,Contig,Begin,End,Function,feature-type,alias] # tab-separated
    add_fid_link                    FID Link (e.g., <a href=http//...>link to somewhere</a>)
    add_fid_links                   File of [PEG,Link] pairs (tab-separated)
    add_genome                      User GenomeDir [-force] [-skipnr] [-dont_mark_complete] (-force will ignore verify dir, -skipnr will not add proteins to nr, -dont_mark_complete will keep PROBABLY_COMPLETE from being copied to COMPLETE if assess_completeness ran)
    add_pch_pins                    File
    all_co_occurrence_clusters      
    all_co_occurrence_pair_sets
    all_compounds               
    all_contigs                     GenomeID
    all_exchangable_subsystems
    all_families
    all_features                    GenomeID Type
    all_maps                 
    all_protein_families
    all_reactions
    all_roles
    all_sets                        Relation SetName
    assign_function                 PEG User [conf=X] Function
    assign_functionF                User File
    assignments_made                who date G1 G2 ...  [ none => all ]
    auto_assign                     PEG [Seq]
    auto_assignF                    User FileOfPEGs
    auto_assignG                    User GenomeID1 GenomeID2 ...
    bbhs                            PEG Cutoff
    best_bbh_candidates             Genome Cutoff Requested Known
    between                         n1 n2 n3
    blast                           PEG  [against nr]
    blastit                         PEG seq db maxP
    boundaries_of                   Loc
    build_tree_of_complete          min_for_label   
    by_alias                        alias
    by_fig_id                       FID1 FID2
    candidates_for_role             Genome Cutoff Role
    cas                             CID
    cas_to_cid                      cas
    catalyzed_by                    RID
    catalyzes                       role
    cgi_url
    change_location_of_feature      FID Location [Translation]
    cid_to_prots                    CID
    clean_tmp
    cleanup_pair_sets
    close_genes                     FID distance
    co_occurs                       PEG
    co_occurrence_cluster           CLUSTER
    co_occurrence_and_evidence      PEG1
    co_occurrence_evidence          PEG1 PEG2
    co_occurrence_set               SET
    co_occurring_FIGfams            FIGFAM
    comp2react                      CID
    contig_ln                       GenomeID Contig
    coupled_to                      PEG
    coupling_and_evidence           FID Bound SimCutoff CouplingCutoff
    coupling_evidence               PEG1 PEG2
    crude_estimate_of_distance      GenomeID1 GenomeID2
    delete_all_fid_links            FID
    delete_feature                  User FID
    delete_features                 User File
    delete_genomes                  G1 G2 G3 ...Gn
    delete_fid_link                 FID URL
    differentiate_families	    FAM1 FAM2
    displayable_reaction            RID
    dna_seq                         GenomeID Loc
    dsims                           FastaFile [MaxN MaxPsc Select]   *** select defaults to raw
    ec_to_maps                      EC
    ec_name                         EC
    epoch_date                      MM/DD/YYYY
    expand_ec                       EC
    epoch_to_readable               [time]   [gives readable version of time, or of current time if not specified]
    export_chromosomal_clusters
    export_pch_pins
    export_set                      Relation SetName File
    exportable_subsystem            Subsystem
    extend_pairs_and_pair_sets      PCH-file
    external_family_for_id          ID (ID should be e.g. kegg|YPS:YPTB1307)
    external_ids_in_family          Family (family should be e.g. kegg|K02035)
    external_calls                  PEG1 PEG2 PEG3...
    extract_seq                     ContigsFile loc 
    families_for_protein            PEG
    family_function                 Family
    fast_coupling                   FID Bound CouplingCutoff
    feature_aliases                 FID
    feature_annotations             FID
    feature_attributes              FID
    function_to_subsystems          Func
    get_key_value		    key= value=  (key or value are optional)
    get_all_keys
    get_all_attributes
    erase_attribute_entirely        Key
    feature_location                FID	
    fid_links                       FID
    fids_with_link_to               text string
    file2N                          File
    find_by_attribute               searchString (substring of tag or value)
    ftype                           FID 
    function_of                     ID [all] [trans]
    genes_in_region                 GenomeID Contig Beg End 
    genome_of                       PEG
    genomes                         [complete|Pat]
    genome_counts                   [complete]
    genome_version                  GenomeID
    genus_species                   GenomeID
    get_dna_seq                     ID
    get_representative_genome	    GenomeID
    get_title_by_pmid               PMID
    get_titles_by_gi                GI
    get_titles_by_peg               PEG
    get_translation                 ID
    get_translations                File [tab]
    go_number_to_info               ID
    go_number_to_term               ID
    h                               [command]  *** h h for help on how to use help ***
    hypo                            Function
    ids_in_family                   Family
    ids_in_set                      WhichSet Relation SetName
    in_cluster_with                 PEG
    in_co_occurrence_cluster        PEG
    in_family                       FID
    in_pair_set                     PEG1 PEG2
    in_pch_pin_with                 FID
    in_sets                         Which Relation SetName
    init_local_sim_computation
    insert_dynamic_sims             File
    install_assignment_set          User File
    is_archaeal                     GenomeID
    is_aux_role_in_subsystem        'Subsystem' 'Role'
    is_bacterial                    GenomeID
    is_cluster_based_subsystem      Subsystem
    is_co_occurrence_pair           PEG1 PEG2
    is_complete                     GenomeID
    is_deleted_fid                  FID
    is_eukaryotic                   GenomeID
    is_NMPDR_genome                 GenomeID
    is_prokaryotic                  GenomeID
    is_exchangable_subsystem        Subsystem
    is_locked_fid                   FID
    is_real_feature                 FID   
    largest_clusters                FileOfRoles user 
    load_all
    load_table			    file  tbl
    localtime                       TimeStamp
    lock_fid                        User FID 
    log_corr                        User Genome CorrespondenceFile AttachedMsg
    map_to_ecs                      map
    map_name                        map
    mapped_prot_ids                 ID
    maps_to_id                      ID 
    mark_deleted_genomes            User G1 G2 G3 ...Gn
    max                             n1 n2 n3 ...
    md5_of_peg                      PEG
    md5_sims                        HASH
    merged_related_annotations      PEG1 PEG2 ... PEGn
    min                             n1 n2 n3 ...
    names_of_compound
    neighborhood_of_role            role        
    n_contigs                       GenomeID (Number of contigs using old-style getting all the data)
    number_of_contigs               GenomeID (Number of contigs using SQL Count)
    number_of_families                       fig, pfam, pir, etc 
    number_of_proteins_in_families           fig, pfam, pir, etc 
    number_of_distinct_proteins_in_families  fig, pfam, pir, etc 
    org_of                          prot_id
    partition_coupling_evidence     Partition1 Partition2
    pegs_of                         GenomeID
    pegs_to_roles_in_subsystems     PEG
    pegs_in_subsystems              GenomeID
    pegs_in_subsystems_by_homology  GenomeID maxN maxP (default maxN=5 maxP=1e-20)
    pegs_in_proteinfams             GenomeID
    pegs_in_proteinfams_by_homology GenomeID maxN maxP (default maxN=5 maxP=1e-20)
    pegs_with_md5                   Hash
    phages
    prots_for_role                  Role
    proteins_in_family              Family
    possibly_truncated              FID
    prot_to_cid                     Prot
    reaction2comp                   RID
    recast_ids                      Pat Id1 Id2 ....
    related_by_func_sim             PEG user
    replace_features_with           OldFID User Genome Type Location [Translation]
    replace_genome                  User OldGenome NewGenome MappingFile [-force] [-skipnr]
    reversible                      RID  
    rnas_of                         GenomeID
    role_to_pegs                    Role
    role_to_subsystems              Role
    roles_of_function               function
    same_func                       'function1'   'function2'
    screwed_up                      PEG
    search_index                    pattern
    seed_sims                       PEG
    seqs_with_role                  Genome who role
    seqs_with_roles_in_genomes      genomes roles made_by
    sims                            PEG maxN maxP select
    sort_fids_by_taxonomy           FID1 FID2 FID3 ...
    sort_genomes_by_taxonomy        GenomeID1 GenomeID2 GenomeID3 ...
    subsystems_for_genome           GenomeID
    subsystems_for_peg              FID
    subsystem_genomes               Subsytem
    sz_family                       family
    taxonomic_groups_of_complete    min_for_labels
    taxonomy_of                     GenomeID
    to_alias                        PEG Type
    translatable                    prot_id
    translate_function              PEG user 
    translated_function_of          PEG user
    translation_length              ID
    undelete_feature                User FID
    unique_functions                PEGs user  [make the pegs comma separated]
    unlock_fid                      User FID
    unmark_deleted_genomes          User G1 G2 G3 ...Gn
    update_subsys_conflicts
    variant_code                    'Subsystem' GENOME
    verify_dir                      dir
    who_set_function                PEG
    change_func_roles               FILE (two column file -> old_func new_func)
    
END
}





