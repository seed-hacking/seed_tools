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

package SproutSearch;
use strict;
use FIG_Config;
use File::Basename;
use Data::Dumper;
use Tracer;

my $uselc=1;
eval {
 require List::Compare;
};
undef $uselc if ($@);

# my $genome_list = qq(
# 93061.1		Staphylococcus aureus NCTC 8325
# 1314.1		Streptococcus pyogenes M5
# 192222.1	Campylobacter jejuni subsp. jejuni NCTC 11168
# 169963.1	Listeria monocytogenes EGD-e
# 265669.1	Listeria monocytogenes str. 4b F2365
# 159288.1	Staphylococcus aureus EMRSA-16 (Str. 252)
# 159289.1	Staphylococcus aureus MSSA (Str. 476)
# 196620.1	Staphylococcus aureus subsp. aureus MW2
# 158878.1	Staphylococcus aureus subsp. aureus Mu50
# 158879.1	Staphylococcus aureus subsp. aureus N315
# 216600.1	Streptococcus pneumoniae 23F
# 171101.1	Streptococcus pneumoniae R6
# 170187.1	Streptococcus pneumoniae TIGR4
# 160490.1	Streptococcus pyogenes M1 GAS
# 198466.1	Streptococcus pyogenes MGAS315
# 186103.1	Streptococcus pyogenes MGAS8232
# 193567.1	Streptococcus pyogenes SSI-1
# 243277.1	Vibrio cholerae O1 biovar eltor str. N16961
# 223926.1	Vibrio parahaemolyticus RIMD 2210633
# 216895.1	Vibrio vulnificus CMCP6
# 196600.1	Vibrio vulnificus YJ016
# );

# TSetup("3 SproutSearch ", "ERROR");


# my $genera = ["Staphylococcus",
# 	      "Campylobacter",
# 	      "Listeria",
# 	      "Vibrio",
# 	      "Streptococcus"];

my($genome_list, $genera, %genus_orgs, %genera);

# for my $a (split(/\n/, $genome_list))
# {

#     if ($a =~ /^(\d+\.\d+)\s+(\S+)/)
#     {
# 	if ($genera{$2})
# 	{
# 	    $genus_orgs{$2}->{$1}++;
# 	}
#     }
# }

sub configure_groups
{
    my($sprout) = @_;

    my %groups = $sprout->GetGroups();

    $genera = [];
    for my $group (keys %groups)
    {
	push(@$genera, $group);
	$genera{$group}++;
	    
	my $org_list = $groups{$group};

	map { $genus_orgs{$group}->{$_}++ } @$org_list;
    }
}

sub new
{
    my($class, $fig, $genus_defs) = @_;

    my $self = {
	genus_defs => $genus_defs,
	genome_filters => undef,
	fig => $fig,
	options => {
	    regexp => 0,
	},
	index_dir => "$FIG_Config::sproutData/Indexes",
    };

    if (!$genera)
    {
	configure_groups($fig->{sprout});
    }
	

    return bless($self, $class);
}

sub option
{
    my($self, $name, $val) = @_;

    if (defined($val))
    {
	my $old = $self->{options}->{name};
	$self->{options}->{$name} = $val;
	return $old;
    }
    else
    {
	return $self->{options}->{$name};
    }
}

sub add_genome_filter
{
    my($self, @genomes) = @_;

    map { $self->{genome_filters}->{$_}++ } @genomes;
}

sub add_genus_filters
{
    my($self, @genera) = @_;

    for my $genus (@genera)
    {
	$self->add_genome_filter(@{$self->{genus_defs}->{$genus}});
    }
}


#
# A sprout search returns a set of objects that somehow match the search.
#
# These objects can be features, genomes, or subsystems.
# 

sub search
{
    my($self, $search_string) = @_;

    my $features = [];
    my $genomes = [];
    my $subsystems = [];
    #
    # First determine if this is a feature id.
    #

    $search_string =~ s/^\s+//;
    $search_string =~ s/\s+$//;

    if ($search_string =~ /^fig\|/)
    {
#	my @annos = $self->{fig}->feature_annotations($search_string);
	return ([$search_string], [], [], []);
    }

    #
    # Compute filtering stuff.
    #
    my @genomes = $self->{genome_filters} ? @{$self->{genome_filters}} : ();
    my @filter_genus = $self->{genus_filters} ? @{$self->{genus_filters}} : ();

    #
    # If neither genus nor genome is set, set to include all genera.
    #

    for my $genus (@filter_genus)
    {
	if ($genus_orgs{$genus})
	{
	    push(@genomes, keys %{$genus_orgs{$genus}});
	}
    }

    my %genomes;
    map { $genomes{$_}++} @genomes;

    #
    # Now, any returned peg must be in %genomes to be displayed.
    #

    my @words = $search_string =~ /\S+/g;

    my @result_sets;
    my $stats;

    for my $word (@words)
    {
	my $set = $self->search_word($word);
	push(@result_sets, $set);
	push(@$stats, [$word, int(@$set)]);
    }
    open(O, ">/tmp/oo");
    print O Dumper(@result_sets);
    close(O);

    my $out;
    if (@result_sets > 1)
    {
	if ($uselc)
	{
	    my $lc = new List::Compare(@result_sets);
	    $out = [$lc->get_intersection()];
	}
	else
	{
	    # here's a quick intesector
	    # @results_sets has references to arrays
	    my %result_count;
	    my $array_count;
	    foreach my $arr (@result_sets)
	    {
		$array_count++;
		foreach (@$arr)
		{
		    $result_count{$_}++
		    }
	    } # can we do this with maps. Probably.

	    $out = [];
	    foreach (keys %result_count)
	    {
		# this requires that the elment is in all arrays
		push @$out, $_ if ($result_count{$_} == $array_count);
	    } 
        }
    }
    else
    {
	#warn "result sets else: ", Dumper(@result_sets);
	$out = $result_sets[0];
    }

    #warn "Returning out=$out\n";
    return ($out, $genomes, $subsystems, $stats, [@words]);
}

sub search_phrase
{
    my($self, $search_string) = @_;

    my $features = [];
    my $genomes = [];
    my $subsystems = [];

    #
    # First determine if this is a feature id.
    #

    $search_string =~ s/^\s+//;
    $search_string =~ s/\s+$//;

    if ($search_string =~ /^fig\|/)
    {
	my @annos = $self->{fig}->feature_annotations($search_string);
	return ([$search_string], [], [], []);
    }

    #
    # Compute filtering stuff.
    #
    my @genomes = $self->{genome_filters} ? @{$self->{genome_filters}} : ();
    my @filter_genus = $self->{genus_filters} ? @{$self->{genus_filters}} : ();

    #
    # If neither genus nor genome is set, set to include all genera.
    #

    for my $genus (@filter_genus)
    {
	if ($genus_orgs{$genus})
	{
	    push(@genomes, keys %{$genus_orgs{$genus}});
	}
    }

    my %genomes;
    map { $genomes{$_}++} @genomes;

    #
    # Now, any returned peg must be in %genomes to be displayed.
    #

    my $phrase = $search_string;
    $phrase =~ s/\s+/;/g;
    $phrase =~ s/\./\\./g;
    
    my @result_sets;
    my $stats;

    my $was_re = $self->option("regexp", 1);

    my $set = $self->search_word($phrase);

    $self->option("regexp", $was_re);
    
    push(@result_sets, $set);
    push(@$stats, [$phrase, int(@$set)]);

    my $out = $set;
	
    return ($out, $genomes, $subsystems, $stats, [$search_string]);
}

sub search_word
{
    my($self, $word) = @_;

    #
    # Do a glimpse search on word, returning a list of ids.
    #

    #
    # First check if it's in the inverted index.
    #

    if (1 or $ENV{USE_DB})
    {
	my $stmt = $self->{inverted_stmt};
	
	if (!$stmt)
	{
	    my $dbh = $self->{fig}->{fig}->db_handle()->{_dbh};
	    $stmt = $dbh->prepare("select distinct peg from sprout_search_terms where word = ?");
	    $self->{inverted_stmt} = $stmt;
	}
	
	$stmt->execute($word);
	my $res = $stmt->fetchall_arrayref();
	
	if ($res and @$res > 0)
	{
	    my $n = @$res;
	    #warn "inverted index found $n hits for $word\n";
	    my $ret = [map { $_->[0] } @$res];
	    return $ret;
	}
    }
	

    my @glimpse_args = ('-w', '-h', '-y', '-i', '-H', $self->{index_dir});
    
    if ($self->option("regexp"))
    {
    }
    else
    {
	push(@glimpse_args, '-k');
    }

    push(@glimpse_args, $word);

    Trace("Glimpse: $FIG_Config::ext_bin/glimpse @glimpse_args") if T(3);

    open(GL, "-|", "$FIG_Config::ext_bin/glimpse", @glimpse_args) or die "Cannot open glimpse:$ !";

    my %set;

    #
    # See if we appear to be searching for an EC # (or an IP address :-)
    # If we are, require that word to be present in the output as a word (not a substring)
    #

    my $re = '^(fig\|[^\t]*)\t';

    if ($word =~ /^\d+\.\d+\.\d+\.\d+$/)
    {
	my $wordre = $word;
	$wordre =~ s/\./\\\./g;
	    
	$re .= ".*\\b$wordre\\b";
    }
    # print "word=$word re=$re\n";

    while (<GL>)
    {
	chomp;
	if (/$re/)
	{
	    # warn "$_\n";
	    $set{$1}++;
	}

    }
    close(GL);

    return [keys(%set)];
}

sub search_old
{
    my($self, $search_string) = @_;

    #
    # First determine if this is a feature id.
    #

    if ($search_string =~ /^fig\|/)
    {
	my @annos = $self->{fig}->feature_annotations($search_string);
	return ([[$search_string, $annos[0]->{text}, undef]], [], []);
    }

    #
    # Is it an alias?
    #

    if (my @feats = $self->{fig}->{sprout}->FeaturesByAlias($search_string))
    {
	my $featret = [];
	for my $feat (@feats)
	{
	    my @annos = $self->{fig}->feature_annotations($feat);
	    # warn "$feat ", Dumper(@annos);
	    push(@$featret, [$feat, $annos[0]->[3], undef]);
	}
	return ($featret, undef, undef);
    }

    #
    # Compute filtering stuff.
    #
    my @genomes = $self->{genome_filters} ? @{$self->{genome_filters}} : ();
    my @filter_genus = $self->{genus_filters} ? @{$self->{genus_filters}} : ();

    #
    # If neither genus nor genome is set, set to include all genera.
    #

    for my $genus (@filter_genus)
    {
	if ($genus_orgs{$genus})
	{
	    push(@genomes, keys %{$genus_orgs{$genus}});
	}
    }

    my %genomes;
    map { $genomes{$_}++} @genomes;

    #
    # Now, any returned peg must be in %genomes to be displayed.
    #
    
    my $index_dir = $self->{index_dir};
    
    my @glimpse_args = ('-y', '-i', '-H', $index_dir);
    
    if ($self->option("regexp"))
    {
    }
    else
    {
	push(@glimpse_args, '-k');
    }
    
    
    push(@glimpse_args, $search_string);
    
    warn "args: @glimpse_args\n";
    
    open(GL, "-|", "$FIG_Config::ext_bin/glimpse", @glimpse_args) or die "Cannot open glimpse:$ !";
    
    my (@annos, @alias, @org, @path);

    #
    # The general scheme here is that we match hits from the glimpse output based
    # on the table that they hit. 
    #

    #
    # Output lists.
    #
    # Features are of the form [$fid, $annotation, $alias]
    #
    # Genomes are of the form [$gid, $name].
    #

    my $features = [];
    my $genomes = [];
    my $subsystems = [];

    while (<GL>)
    {
	chomp;
	s/\r//;

	#
	# Detect the filename part of the output.
	#
	if (/(^[^:]+):\s+(.*)$/)
	{
	    my $file = $1;
	    my $rest = $2;

	    my $table = basename($file);
	    $table =~ s/\.dtx$//;

	    if ($table eq "Annotation")
	    {
		my($key, $time, $anno) = split(/\t/, $rest, 3);
		my $peg;

		if ($key !~ /^fig\|/)
		{
		    #
		    # Need to find the fig id that this annotation is the target of.
		    #

		    my $ret = $self->{fig}->{sprout}->Get(['IsTargetOfAnnotation'],
						  'IsTargetOfAnnotation(to-link) = ?',
						  [$key]);
		    my $data = $ret->Fetch();
		    my @pegs = $data->Values(['IsTargetOfAnnotation(from-link)']);
		    $peg = $pegs[0];
		}
		else
		{
		    $peg = $key;
		    $peg =~ s/:.*$//;
		}
		$anno =~ s/\\n/\n/g;

		if ($self->feature_survives_filter($peg))
		{
		    push(@$features, [$peg, $anno, undef]);
		}
	    }
	    elsif ($table eq "ComesFrom")
	    {
	    }
	    elsif ($table eq "FeatureAlias")
	    {
		my($peg, $alias) = split(/\t/, $rest);
		if ($self->feature_survives_filter($peg))
		{
		    push(@$features, [$peg, undef, $alias]);
		}
	    }
	    elsif ($table eq "Genome")
	    {
		my($org, $whatisthis, $genus, $species, $tax) = split(/\t/, $rest);
		if ($self->genome_survives_filter($org))
		{
		    push(@$genomes, [$org, $genus, $species]);
		}
	    }
	}
    }

    return ($features, $genomes, $subsystems);
}

sub genome_survives_filter
{
    my($self, $gid) = @_;

    if (defined($self->{genome_filters}))
    {
	return $self->{genome_filters}->{$gid};
    }
    else
    {
	return 1;
    }
}

sub feature_survives_filter
{
    my($self, $feature) = @_;

    #
    # Right now, just filter for fig id's.
    #

    if ($feature =~ /^fig\|(\d+\.\d+)\./)
    {
	my $genome = $1;
	return $self->genome_survives_filter($genome);
    }
    else
    {
	return 0;
    }
}
1;
