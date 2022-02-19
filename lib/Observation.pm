package Observation;

#use lib '/vol/ontologies';
use DBMaster;
use Data::Dumper;

require Exporter;
@EXPORT_OK = qw(get_objects get_sims_objects); 

use WebColors;
use WebConfig;

use FIG_Config;
use LWP::Simple;
#use strict;
#use warnings;
use HTML;
use FFs;

1;

=head1 NAME

Observation -- A presentation layer for observations in SEED.

=head1 DESCRIPTION

The SEED environment contains various sources of information for sequence features. The purpose of this library is to provide a 
single interface to this data.

The data can be used to display information for a given sequence feature (protein or other, but primarily information is computed for proteins). 

=cut

=head1 BACKGROUND

=head2 Data incorporated in the Observations 

As the goal of this library is to provide an integrated view, we combine diverse sources of evidence.

=head3 SEED core evidence

The core SEED data structures provided by FIG.pm. These are Similarities, BBHs and PCHs.

=head3 Attribute based Evidence

We use the SEED attribute infrastructure to store information computed by a variety of computational procedures.

These are e.g. InterPro hits via InterProScan (ipr), NCBI Conserved Domain Database Hits via PSSM(cdd), 
PFAM hits via HMM(pfam), SignalP results(signalp), and various others.

=head1 METHODS

The public methods this package provides are listed below:


=head3 context()

Returns close or diverse for purposes of displaying genomic context

=cut

sub context {
  my ($self) = @_;

  return $self->{context};
}

=head3 rows()

each row in a displayed table

=cut

sub rows {
  my ($self) = @_;

  return $self->{rows};
}

=head3 acc()

A valid accession or remote ID (in the style of a db_xref) or a valid local ID (FID) in case this is supported.

=cut

sub acc {
  my ($self) = @_;
  return $self->{acc};
}

=head3 query()

The query id

=cut

sub query {
    my ($self) = @_;
    return $self->{query};
}


=head3 class()

The class of evidence (required). This is usually simply the name of the tool or the name of the SEED data structure.
B<Please note> the connection of class and display_method and URL.
    
Current valid classes are:

=over 9

=item IDENTICAL (seq)

=item SIM (seq)

=item BBH (seq)

=item PCH (fc)

=item FIGFAM (seq)

=item IPR (dom)

=item CDD (dom)

=item PFAM (dom)

=item SIGNALP_CELLO_TMPRED (loc)

=item PDB (seq)

=item TMHMM (loc)

=item HMMTOP (loc)

=back

=cut

sub class {
  my ($self) = @_;

  return $self->{class};
}

=head3 type()

The type of evidence (required).

Where type is one of the following:

=over 8

=item seq=Sequence similarity

=item dom=domain based match

=item loc=Localization of the feature

=item fc=Functional coupling.

=back

=cut

sub type {
  my ($self) = @_;

  return $self->{type};
}

=head3 start()

Start of hit in query sequence.

=cut

sub start {
  my ($self) = @_;

  return $self->{start};
}

=head3 end()

End of the hit in query sequence.

=cut

sub stop {
  my ($self) = @_;

  return $self->{stop};
}

=head3 start()

Start of hit in query sequence.

=cut

sub qstart {
    my ($self) = @_;

    return $self->{qstart};
}

=head3 qstop()

End of the hit in query sequence.

=cut

sub qstop {
    my ($self) = @_;

    return $self->{qstop};
}

=head3 hstart()

Start of hit in hit sequence.

=cut

sub hstart {
    my ($self) = @_;

    return $self->{hstart};
}

=head3 end()

End of the hit in hit sequence.

=cut

sub hstop {
    my ($self) = @_;

    return $self->{hstop};
}

=head3 qlength()

length of the query sequence in similarities

=cut

sub qlength {
    my ($self) = @_;

    return $self->{qlength};
}

=head3 hlength()

length of the hit sequence in similarities

=cut

sub hlength {
    my ($self) = @_;

    return $self->{hlength};
}

=head3 evalue()

E-value or P-Value if present.

=cut

sub evalue {
  my ($self) = @_;

  return $self->{evalue};
}

=head3 score()

Score if present. 

=cut

sub score {
  my ($self) = @_;
  return $self->{score};
}

=head3 display()

will be different for each type

=cut 

sub display {
  
  die "Abstract Method Called\n";

}

=head3 display_table()

will be different for each type

=cut 

sub display_table {
  
  die "Abstract Table Method Called\n";

}

=head3 get_objects()

This is the B<REAL WORKHORSE> method of this Package.

=cut

sub get_objects {
    my ($self,$fid,$fig,$parameters,$scope) = @_;

    my $objects = [];
    my @matched_datasets=();

    # call function that fetches attribute based observations
    # returns an array of arrays of hashes
 
    if($scope){
	get_cluster_observations($fid,\@matched_datasets,$scope);
    }
    else{
	my %domain_classes;
	my @attributes = $fig->get_attributes($fid);
	#$domain_classes{'CDD'} = 1;
	$domain_classes{'PFAM'} = 1;
	get_identical_proteins($fid,\@matched_datasets,$fig);
	get_attribute_based_domain_observations($fid,\%domain_classes,\@matched_datasets,\@attributes,$fig);
	get_sims_observations($fid,\@matched_datasets,$fig,$parameters);
	get_functional_coupling($fid,\@matched_datasets,$fig);
	get_attribute_based_location_observations($fid,\@matched_datasets,\@attributes,$fig);
	get_pdb_observations($fid,\@matched_datasets,\@attributes,$fig);
    }
 
    foreach my $dataset (@matched_datasets) {
	my $object;
	if($dataset->{'type'} eq "dom"){
	    $object = Observation::Domain->new($dataset);
	}
	elsif($dataset->{'class'} eq "PCH"){
            $object = Observation::FC->new($dataset);
        }
	elsif ($dataset->{'class'} eq "IDENTICAL"){
	    $object = Observation::Identical->new($dataset);
	}
	elsif ($dataset->{'class'} eq "SIGNALP_CELLO_TMPRED"){
	    $object = Observation::Location->new($dataset);
	}
	elsif ($dataset->{'class'} eq "SIM"){
            $object = Observation::Sims->new($dataset);
        }
	elsif ($dataset->{'class'} eq "CLUSTER"){
            $object = Observation::Cluster->new($dataset);
        }
	elsif ($dataset->{'class'} eq "PDB"){
            $object = Observation::PDB->new($dataset);
        }
	
	push (@$objects, $object);
    }
    
    return $objects;

}

=head3 get_attributes
    provides layer of abstraction between tools and underlying access method to Attribute Server
=cut 

sub get_attributes{
    my ($self,$fig,$search_set,$search_term,$value_array_ref) = @_;
    my @attributes = $fig->get_attributes($search_set,$search_term,@$value_array_ref);
    return @attributes;
}

=head3 get_sims_objects()

This is the B<REAL WORKHORSE> method of this Package.

=cut

sub get_sims_objects {
    my ($self,$fid,$fig,$parameters) = @_;

    my $objects = [];
    my @matched_datasets=();

    # call function that fetches attribute based observations
    # returns an array of arrays of hashes
    get_sims_observations($fid,\@matched_datasets,$fig,$parameters);
 
    foreach my $dataset (@matched_datasets) {
	my $object;
	if ($dataset->{'class'} eq "SIM"){
            $object = Observation::Sims->new($dataset);
        }
	push (@$objects, $object);
    }
    return $objects;
}


=head3 display_housekeeping
This method returns the housekeeping data for a given peg in a table format

=cut
sub display_housekeeping {
    my ($self,$fid,$fig) = @_;
    my $content = [];
    my $row = [];

    my $org_name = "Data not available";
    if ( $fig->org_of($fid)){
	$org_name = $fig->org_of($fid);
    }
    my $org_id = $fig->genome_of($fid);
    my $function = $fig->function_of($fid);
    #my $taxonomy = $fig->taxonomy_of($org_id);
    my $length = $fig->translation_length($fid);
    
    push (@$row, $org_name);
    push (@$row, $fid);
    push (@$row, $length);
    push (@$row, $function);

    # initialize the table for commentary and annotations
    #$content .= qq(<b>My Sequence Data</b><br><table border="0">);
    #$content .= qq(<tr width=15%><td >FIG ID</td><td>$fid</td></tr>\n);
    #$content .= qq(<tr width=15%><td >Organism Name</td><td>$org_name</td></tr>\n);
    #$content .= qq(<tr><td width=15%>Taxonomy</td><td>$taxonomy</td></tr>\n);
    #$content .= qq(<tr width=15%><td>Function</td><td>$function</td></tr>\n);
    #$content .= qq(<tr width=15%><td>Sequence Length</td><td>$length aa</td></tr>\n);
    #$content .= qq(</table><p>\n);
    
    push(@$content, $row);

    return ($content);
}

=head3 get_sims_summary
This method uses as input the similarities of a peg and creates a tree view of their taxonomy

=cut

sub get_sims_summary {
    my ($observation, $dataset, $fig) = @_;
    my %families;
    my $taxes = $fig->taxonomy_list();
    
    foreach my $thing (@$dataset) {
	my ($id, $evalue);
	if ($thing =~ /fig\|/){
	    $id = $thing;
	    $evalue = -1;
	}
	else{
	    next if ($thing->class ne "SIM");
	    $id      = $thing->acc;
	    $evalue  = $thing->evalue;
	}
        next if ($id !~ /fig\|/);
	next if ($fig->is_deleted_fid($id));

        my $genome = $fig->genome_of($id);
	#my ($genome1) = ($genome) =~ /(.*)\./;
	my $taxonomy = $taxes->{$genome};
        my $parent_tax = "Root";
	my @currLineage = ($parent_tax);
	push (@{$families{figs}{$parent_tax}}, $id);
	my $level = 2;

        foreach my $tax (split(/\; /, $taxonomy),$id){
	  next if ($tax eq $parent_tax);
	  push (@{$families{children}{$parent_tax}}, $tax) if ($tax ne $parent_tax);
	  push (@{$families{figs}{$tax}}, $id) if ($tax ne $parent_tax);
	  $families{level}{$tax} = $level;
	  push (@currLineage, $tax);
	  $families{parent}{$tax} = $parent_tax;
	  $families{lineage}{$tax} = join(";", @currLineage);
	  if (defined ($families{evalue}{$tax})){
	    if ($evalue < $families{evalue}{$tax}){
	      $families{evalue}{$tax} = $evalue;
	      $families{color}{$tax} = &get_taxcolor($evalue);
	    }
	  }
	  else{
	    $families{evalue}{$tax} = $evalue;
	    $families{color}{$tax} = &get_taxcolor($evalue);
	  }
	  
	  $parent_tax = $tax;
	  $level++;
        }
    }

    foreach my $key (keys %{$families{children}}){
        $families{count}{$key} = @{$families{children}{$key}};

        my %saw;
        my @out = grep(!$saw{$_}++, @{$families{children}{$key}});
        $families{children}{$key} = \@out;
    }

    return \%families;
}

=head1 Internal Methods 

These methods are not meant to be used outside of this package. 

B<Please do not use them outside of this package!>

=cut

sub get_taxcolor{
    my ($evalue) = @_;
    my $color;
    if ($evalue == -1){            $color = "black";      }
    elsif (($evalue <= 1e-170) && ($evalue >= 0)){        $color = "#FF2000";    }
    elsif (($evalue <= 1e-120) && ($evalue > 1e-170)){        $color = "#FF3300";    }
    elsif (($evalue <= 1e-90) && ($evalue > 1e-120)){        $color = "#FF6600";    }
    elsif (($evalue <= 1e-70) && ($evalue > 1e-90)){        $color = "#FF9900";    }
    elsif (($evalue <= 1e-40) && ($evalue > 1e-70)){        $color = "#FFCC00";    }
    elsif (($evalue <= 1e-20) && ($evalue > 1e-40)){        $color = "#FFFF00";    }
    elsif (($evalue <= 1e-5) && ($evalue > 1e-20)){        $color = "#CCFF00";    }
    elsif (($evalue <= 1) && ($evalue > 1e-5)){        $color = "#66FF00";    }
    elsif (($evalue <= 10) && ($evalue > 1)){        $color = "#00FF00";    }
    else{        $color = "#6666FF";    }
    return ($color);
}


sub get_attribute_based_domain_observations{

    # we read a FIG ID and a reference to an array (of arrays of hashes, see above)
    my ($fid,$domain_classes,$datasets_ref,$attributes_ref,$fig) = (@_);
    my $seen = {};
    foreach my $attr_ref (@$attributes_ref) {
	my $key = @$attr_ref[1];
	my @parts = split("::",$key);
	my $class = $parts[0];
	my $name = $parts[1];
	next if ($seen->{$name});
	$seen->{$name}++;
	#next if (($class eq "PFAM") && ($name !~ /interpro/));

	if($domain_classes->{$parts[0]}){
	    my $val = @$attr_ref[2];
	    if($val =~/^(\d+\.\d+|0\.0);(\d+)-(\d+)/){
		my $raw_evalue = $1;
		my $from = $2;
		my $to = $3;
		my $evalue;
		if(($raw_evalue =~/(\d+)\.(\d+)/) && ($class ne "PFAM")){
		    my $part2 = 1000 - $1;
		    my $part1 = $2/100;
		    $evalue = $part1."e-".$part2;
		}
		elsif(($raw_evalue =~/(\d+)\.(\d+)/) && ($class eq "PFAM")){
		    #$evalue=$raw_evalue;
		    my $part2 = 1000 - $1;
                    my $part1 = $2/100;
                    $evalue = $part1."e-".$part2;

		}
		else{
		    $evalue = "0.0";
		}
		
		my $dataset = {'class' => $class,
			       'acc' => $key,
			       'type' => "dom" ,
			       'evalue' => $evalue,
			       'start' => $from,
			       'stop' => $to,
			       'fig_id' => $fid,
			       'score' => $raw_evalue
			       };
		
		push (@{$datasets_ref} ,$dataset);
	    }
	}
    }
}

sub get_attribute_based_location_observations{

    my ($fid,$datasets_ref, $attributes_ref,$fig) = (@_);
    #my $fig = new FIG;
    
    my $location_attributes = ['SignalP','CELLO','TMPRED','Phobius'];
    
    my $dataset = {'type' => "loc", 
		   'class' => 'SIGNALP_CELLO_TMPRED',
		   'fig_id' => $fid
		   };

    foreach my $attr_ref (@$attributes_ref){
	my $key = @$attr_ref[1];
	next if (($key !~ /SignalP/) && ($key !~ /CELLO/) && ($key !~ /TMPRED/)  && ($key !~/Phobius/) );
	my @parts = split("::",$key);
	my $sub_class = $parts[0];
	my $sub_key = $parts[1];
	my $value = @$attr_ref[2];
	if($sub_class eq "SignalP"){
	    if($sub_key eq "cleavage_site"){
		my @value_parts = split(";",$value);
		$dataset->{'cleavage_prob'} = $value_parts[0];
		$dataset->{'cleavage_loc'} = $value_parts[1];
	    }
	    elsif($sub_key eq "signal_peptide"){
		$dataset->{'signal_peptide_score'} = $value;
	    }
	}
	
	elsif($sub_class eq "CELLO"){
	    $dataset->{'cello_location'} = $sub_key;
	    $dataset->{'cello_score'} = $value;
	}
	
	elsif($sub_class eq "Phobius"){
	    if($sub_key eq "transmembrane"){
		$dataset->{'phobius_tm_locations'} = $value;
	    }
	    elsif($sub_key eq "signal"){
		$dataset->{'phobius_signal_location'} = $value;
	    }
	}
	
	elsif($sub_class eq "TMPRED"){
	    my @value_parts = split(/\;/,$value);
	    $dataset->{'tmpred_score'} = $value_parts[0];
	    $dataset->{'tmpred_locations'} = $value_parts[1];
	}
    }
    
    push (@{$datasets_ref} ,$dataset);
    
}

=head3 get_pdb_observations() (internal)

This methods sets the type and class for pdb observations

=cut

sub get_pdb_observations{
    my ($fid,$datasets_ref, $attributes_ref,$fig) = (@_);
    
    #my $fig = new FIG;
    
    foreach my $attr_ref (@$attributes_ref){
	my $key = @$attr_ref[1];
	next if ( ($key !~ /PDB/));
	my($key1,$key2) =split("::",$key);
	my $value = @$attr_ref[2];
	my ($evalue,$location) = split(";",$value);
	
	if($evalue =~/(\d+)\.(\d+)/){
	    my $part2 = 1000 - $1;
	    my $part1 = $2/100;
	    $evalue = $part1."e-".$part2;
	} 

	my($start,$stop) =split("-",$location);

	my $url = @$attr_ref[3];
	my $dataset = {'class' => 'PDB',
		       'type' => 'seq' ,
		       'acc' => $key2,
		       'evalue' => $evalue,
                       'start' => $start,
                       'stop' => $stop,
		       'fig_id' => $fid
		       };

	push (@{$datasets_ref} ,$dataset);
    }
}

=head3 get_cluster_observations() (internal)

This methods sets the type and class for cluster observations

=cut

sub get_cluster_observations{
    my ($fid,$datasets_ref,$scope) = (@_);

    my $dataset = {'class' => 'CLUSTER',
		   'type' => 'fc',
		   'context' => $scope,
		   'fig_id' => $fid
		   };
    push (@{$datasets_ref} ,$dataset);
}


=head3 get_sims_observations() (internal)

This methods retrieves sims fills the internal data structures.

=cut

sub get_sims_observations{
    my ($fid,$datasets_ref,$fig,$parameters) = (@_);

    my ($max_sims, $max_expand, $max_eval, $sim_order, $db_filter, $sim_filters);
    if ( (defined $parameters->{flag}) && ($parameters->{flag})){
      $max_sims = $parameters->{max_sims};
      $max_expand = $parameters->{max_expand};
      $max_eval = $parameters->{max_eval};
      $db_filter = $parameters->{db_filter};
      $sim_filters->{ sort_by } = $parameters->{sim_order};
      #$sim_order = $parameters->{sim_order};
      $group_by_genome = 1 if (defined ($parameters->{group_genome}));
    }
    elsif ( (defined $parameters->{sims_db}) && ($parameters->{sims_db} eq 'all')){
      $max_sims = 50;
      $max_expand = 5;
      $max_eval = 1e-5;
      $db_filter = "all";
      $sim_filters->{ sort_by } = 'id';
    }
    else{
      $max_sims = 50;
      $max_expand = 5;
      $max_eval = 1e-5;
      $db_filter = "figx";
      $sim_filters->{ sort_by } = 'id';
      #$sim_order = "id";
    }
    
    my($id, $genome, @genomes, %sims);
#    my @tmp= $fig->sims($fid,$max_sims,$max_eval,$db_filter,$max_expand,$sim_filters);
    my @tmp= $fig->sims($fid,1000000,$max_eval,$db_filter,$max_expand,$sim_filters);
    @tmp = grep { !($_->id2 =~ /^fig\|/ and $fig->is_deleted_fid($_->id2)) } @tmp;
    my ($dataset);

    if ($group_by_genome){
      #  Collect all sims from genome with the first occurance of the genome:
      foreach $sim ( @tmp ){
	$id = $sim->id2;
	$genome = ($id =~ /^fig\|(\d+\.\d+)\.peg\.\d+/) ? $1 : $id;
	if (! defined( $sims{ $genome } ) ) { push @genomes, $genome }
	push @{ $sims{ $genome } }, $sim;
      }
      @tmp = map { @{ $sims{$_} } } @genomes;
    }
    
    my $seen_sims={};
    my $count=1;
    foreach my $sim (@tmp){

	my $hit = $sim->[1];
	next if ($seen_sims->{$hit});
	next if ($hit =~ /nmpdr\||gnl\|md5\|/);
	$seen_sims->{$hit}++;

	last if ($count>$max_sims);
	$count++;

	my $percent = $sim->[2];
	my $evalue = $sim->[10];
	my $qfrom = $sim->[6];
	my $qto = $sim->[7];
	my $hfrom = $sim->[8];
	my $hto = $sim->[9];
	my $qlength = $sim->[12];
	my $hlength = $sim->[13];
	my $db = get_database($hit);
	my $func = $fig->function_of($hit);
	my $organism;
	if ($fig->org_of($hit)){
	    $organism = $fig->org_of($hit);
	}
	else{
	    $organism = "Data not available";
	}

	$dataset = {'class' => 'SIM',
		    'query' => $sim->[0],
		    'acc' => $hit,
		    'identity' => $percent,
		    'type' => 'seq',
		    'evalue' => $evalue,
		    'qstart' => $qfrom,
		    'qstop' => $qto,
		    'hstart' => $hfrom,
                    'hstop' => $hto,
		    'database' => $db,
		    'organism' => $organism,
		    'function' => $func,
		    'qlength' => $qlength,
		    'hlength' => $hlength,
		    'fig_id' => $fid
		    };

	push (@{$datasets_ref} ,$dataset);
    }
}

=head3 get_database (internal)
This method gets the database association from the sequence id

=cut

sub get_database{
    my ($id) = (@_);
    
    my ($db);
    if ($id =~ /^fig\|/)              { $db = "SEED" }
    elsif ($id =~ /^gi\|/)            { $db = "NCBI" }
    elsif ($id =~ /^gb\|/)            { $db = "GenBank" }
    elsif ($id =~ /^^[NXYZA]P_/)      { $db = "RefSeq" }
    elsif ($id =~ /^ref\|/)           { $db = "RefSeq" }
    elsif ($id =~ /^sp\|/)            { $db = "SwissProt" }
    elsif ($id =~ /^uni\|/)           { $db = "UniProt" }
    elsif ($id =~ /^tigr\|/)          { $db = "TIGR" }
    elsif ($id =~ /^pir\|/)           { $db = "PIR" }
    elsif (($id =~ /^kegg\|/) || ($id =~ /Spy/))    { $db = "KEGG" }
    elsif ($id =~ /^tr\|/)                          { $db = "TrEMBL" }
    elsif ($id =~ /^eric\|/)          { $db = "ASAP" }
    elsif ($id =~ /^img\|/)           { $db = "JGI" }
    elsif ($id =~ /^pdb\|/)           { $db = "PDB" }
    elsif ($id =~ /^img\|/)           { $db = "IMG" }
    elsif ($id =~ /^cmr\|/)           { $db = "CMR" }
    elsif ($id =~ /^dbj\|/)           { $db = "DBJ" }

    return ($db);

}


=head3 get_identical_proteins() (internal)

This methods retrieves sims fills the internal data structures.

=cut

sub get_identical_proteins{

    my ($fid,$datasets_ref,$fig) = (@_);
    #my $fig = new FIG;
    my $funcs_ref;

    my @maps_to = grep { $_ ne $fid and $_ !~ /^xxx/ } map { $_->[0] } $fig->mapped_prot_ids($fid);
    foreach my $id (@maps_to) {
        my ($tmp, $who);
	if (($id ne $fid) && ($tmp = $fig->function_of($id))) {
	    $who = &get_database($id);
            push(@$funcs_ref, [$id,$who,$tmp]);
        }
    }

    my $dataset = {'class' => 'IDENTICAL',
		   'type' => 'seq',
		   'fig_id' => $fid,
		   'rows' => $funcs_ref   
		   };
    
    push (@{$datasets_ref} ,$dataset);
    

}

=head3 get_functional_coupling() (internal)

This methods retrieves the functional coupling of a protein given a peg ID

=cut

sub get_functional_coupling{

    my ($fid,$datasets_ref,$fig) = (@_);
    #my $fig = new FIG;
    my @funcs = ();

    # initialize some variables
    my($sc,$neigh);

    # set default parameters for coupling and evidence
    my ($bound,$sim_cutoff,$coupling_cutoff) = (5000, 1.0e-10, 4);

    # get the fc data
    my @fc_data = $fig->coupling_and_evidence($fid,$bound,$sim_cutoff,$coupling_cutoff);

    # retrieve data
    my @rows = map { ($sc,$neigh) = @$_;
		     [$sc,$neigh,scalar $fig->function_of($neigh)]
		  } @fc_data;
		     
    my $dataset = {'class' => 'PCH',
		   'type' => 'fc',
		   'fig_id' => $fid,
		   'rows' => \@rows
		   };
    
    push (@{$datasets_ref} ,$dataset);

}

=head3 new (internal)

Instantiate a new object.

=cut

sub new {
  my ($class,$dataset) = @_;
  
  my $self = { class => $dataset->{'class'},
	       type => $dataset->{'type'},
	       fig_id => $dataset->{'fig_id'},
	       score => $dataset->{'score'},
	   };
  
  bless($self,$class);
  
  return $self;
}

=head3 identity (internal)

Returns the % identity of the similar sequence

=cut

sub identity {
    my ($self) = @_;

    return $self->{identity};
}

=head3 fig_id (internal)

=cut

sub fig_id {
  my ($self) = @_;
  return $self->{fig_id};
}

=head3 feature_id (internal)


=cut

sub feature_id {
  my ($self) = @_;

  return $self->{feature_id};
}

=head3 id (internal)

Returns the ID  of the identical sequence

=cut

sub id {
    my ($self) = @_;

    return $self->{id};
}

=head3 organism (internal)

Returns the organism  of the identical sequence

=cut

sub organism {
    my ($self) = @_;

    return $self->{organism};
}

=head3 function (internal)

Returns the function of the identical sequence

=cut

sub function {
    my ($self) = @_;

    return $self->{function};
}

=head3 database (internal)

Returns the database of the identical sequence

=cut

sub database {
    my ($self) = @_;

    return $self->{database};
}

############################################################
############################################################
package Observation::PDB;

use base qw(Observation);

sub new {
    
    my ($class,$dataset) = @_; 
    my $self = $class->SUPER::new($dataset);
    $self->{acc} = $dataset->{'acc'};
    $self->{evalue} = $dataset->{'evalue'};
    $self->{start} = $dataset->{'start'};
    $self->{stop} = $dataset->{'stop'};
    bless($self,$class);
    return $self;
}

=head3 display()

displays data stored in best_PDB attribute and in Ontology server for given PDB id

=cut

sub display{
    my ($self,$gd,$fig) = @_;

    my $fid = $self->fig_id;
    my $dbmaster = DBMaster->new(-database =>'Ontology',
				 -host     => $WebConfig::DBHOST,
				 -user     => $WebConfig::DBUSER,
				 -password => $WebConfig::DBPWD);
    
    my $acc = $self->acc;
  
    my ($pdb_description,$pdb_source,$pdb_ligand);
    my $pdb_objs = $dbmaster->pdb->get_objects( { 'id' => $acc } );
    if(!scalar(@$pdb_objs)){
	$pdb_description = "not available";
	$pdb_source = "not available";
	$pdb_ligand = "not available";
    }
    else{
	my $pdb_obj = $pdb_objs->[0];
	$pdb_description = $pdb_obj->description;
	$pdb_source = $pdb_obj->source;
	$pdb_ligand = $pdb_obj->ligand;
    }

    my $lines = [];
    my $line_data = [];
    my $line_config = { 'title' => "PDB hit for $fid",
			'hover_title' => 'PDB',
			'short_title' => "best PDB",
			'basepair_offset' => '1' };

    #my $fig = new FIG;
    my $seq = $fig->get_translation($fid);
    my $fid_stop = length($seq);

    my $fid_element_hash = {
	"title" => $fid,
	"start" => '1',
	"end" =>  $fid_stop,
	"color"=> '1',
	"zlayer" => '1'
	};
    
    push(@$line_data,$fid_element_hash);
    
    my $links_list = [];
    my $descriptions = [];

    my $name;
    $name = {"title" => 'id',
	     "value" => $acc};
    push(@$descriptions,$name);

    my $description;
    $description = {"title" => 'pdb description',
		    "value" => $pdb_description};
    push(@$descriptions,$description);
    
    my $score;
    $score = {"title" => "score",
	      "value" => $self->evalue};
    push(@$descriptions,$score);

    my $start_stop;
    my $start_stop_value = $self->start."_".$self->stop; 
    $start_stop = {"title" => "start-stop",
		   "value" => $start_stop_value};
    push(@$descriptions,$start_stop);
        
    my $source;
    $source = {"title" => "source",
	      "value" => $pdb_source};
    push(@$descriptions,$source);

    my $ligand;
    $ligand = {"title" => "pdb ligand",
	       "value" => $pdb_ligand};
    push(@$descriptions,$ligand);
 
    my $link;
    my $link_url ="http://www.rcsb.org/pdb/explore/explore.do?structureId=".$acc;
    
    $link = {"link_title" => $acc,
	     "link" => $link_url};
    push(@$links_list,$link);
    
    my $pdb_element_hash = {
	"title" => "PDB homology",
	"start" => $self->start,
	"end" =>  $self->stop,
	"color"=> '6',
	"zlayer" => '3',
	"links_list" => $links_list,
	"description" => $descriptions};
    
    push(@$line_data,$pdb_element_hash);
    $gd->add_line($line_data, $line_config);

    return $gd;
}

1;

############################################################
############################################################
package Observation::Identical;

use base qw(Observation);

sub new {
    
    my ($class,$dataset) = @_; 
    my $self = $class->SUPER::new($dataset);
    $self->{rows} = $dataset->{'rows'};
    
    bless($self,$class);
    return $self;
}

=head3 display_table()

If available use the function specified here to display the "raw" observation.
This code will display a table for the identical protein


B<Please note> that URL linked to in display_method() is an external component and needs to added to the code for every class of evi
dence.

=cut


sub display_table{
    my ($self,$fig) = @_;
    
    #my $fig = new FIG;
    my $fid = $self->fig_id;
    my $rows = $self->rows;
    my $cgi = new CGI;
    my $all_domains = [];
    my $count_identical = 0;
    my $content;
    foreach my $row (@$rows) {
	my $id = $row->[0];
	my $who = $row->[1];
	my $assignment = $row->[2];
	my $organism = "Data not available";
	if ($fig->org_of($id)){
	    $organism = $fig->org_of($id);
	}
        my $single_domain = [];
        push(@$single_domain,$who);
	push(@$single_domain,$self->get_url_for_id($id));
        push(@$single_domain,$organism);
	push(@$single_domain,$assignment);
        push(@$all_domains,$single_domain);
	$count_identical++;
    }

    if ($count_identical >0){
        $content = $all_domains;
    }
    else{
        $content = "<p>This PEG does not have any essentially identical proteins</p>";
    }
    return ($content);
}

sub get_url_for_id {
  my ($self, $id) = @_;

  my $copy = $id;
  if ($copy =~ s/^kegg\|//) {
    return "<a href='http://www.genome.jp/dbget-bin/www_bget?$copy'>$id</a>";
  }
  elsif ($copy =~ s/^sp\|//) {
    return "<a href='http://www.uniprot.org/entry/$copy'>$id</a>";
  }
  elsif ($copy =~ s/^tr\|//) {
    return "<a href='http://www.uniprot.org/entry/$copy'>$id</a>";
  }
  elsif ($copy =~ s/^uni\|//) {
    return "<a href='http://www.uniprot.org/entry/$copy'>$id</a>";
  }
  elsif ($copy =~ s/^gi\|//) {
    return "<a href='http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$copy'>$id</a>";
  }
  elsif ($copy =~ s/^ref\|//) {
    return "<a href='http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$copy'>$id</a>";
  }
  elsif ($copy =~ s/^gb\|//) {
    return "<a href='http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id=$copy'>$id</a>";
  }
  elsif ($copy =~ s/^cmr\|// or $copy =~ s/^tigrcmr\|//) {
    return "<a href='http://cmr.tigr.org/tigr-scripts/CMR/shared/GenePage.cgi?locus=$copy'>$id</a>";
  }
  elsif ($copy =~ /^fig\|/) {
    return "<a href='http://seed-viewer.theseed.org/linkin.cgi?id=$id'>$id</a>";
  }
  elsif ($copy =~ s/^img\|//) {
    return "<a href='http://img.jgi.doe.gov/cgi-bin/pub/main.cgi?section=GeneDetail&page=geneDetail&gene_oid=$copy'>$id</a>";
  }
  else {
    return $id;
  }

}



1;

#########################################
#########################################
package Observation::FC;
1;

use base qw(Observation);

sub new {
    
    my ($class,$dataset) = @_; 
    my $self = $class->SUPER::new($dataset);
    $self->{rows} = $dataset->{'rows'};
    
    bless($self,$class);
    return $self;
}

=head3 display_table()

If available use the function specified here to display the "raw" observation.
This code will display a table for the identical protein


B<Please note> that URL linked to in display_method() is an external component and needs to added to the code for every class of evi
dence.

=cut

sub display_table {

    my ($self,$dataset,$fig) = @_;
    my $fid = $self->fig_id;
    my $rows = $self->rows;
    my $cgi = new CGI;
    my $functional_data = [];
    my $count = 0;
    my $content;

    foreach my $row (@$rows) {
	my $single_domain = [];
	$count++;

	# construct the score link
	my $score = $row->[0];
	my $toid = $row->[1];
	my $link = $cgi->url(-relative => 1) . "?page=Annotation&feature=$fid";
	my $sc_link = "<a href='$link'>$score</a>";

	push(@$single_domain,$sc_link);
	push(@$single_domain,$row->[1]);
	push(@$single_domain,$row->[2]);
	push(@$functional_data,$single_domain);
    }

    if ($count >0){
	$content = $functional_data;
    }
    else
    {
	$content = "<p>This PEG does not have any functional coupling</p>";
    }
    return ($content);
}


#########################################
#########################################
package Observation::Domain;

use base qw(Observation);

sub new {
    
    my ($class,$dataset) = @_; 
    my $self = $class->SUPER::new($dataset);
    $self->{evalue} = $dataset->{'evalue'};
    $self->{acc} = $dataset->{'acc'};
    $self->{start} = $dataset->{'start'};
    $self->{stop} = $dataset->{'stop'};
    
    bless($self,$class);
    return $self;
}

sub display {
    my ($thing,$gd) = @_;
    my $lines = [];
#    my $line_config = { 'title' => $thing->acc,
#			'short_title' => $thing->type,
#			'basepair_offset' => '1' };
    my $color = "4";
    
    my $line_data = [];
    my $links_list = [];
    my $descriptions = [];

    my $db_and_id = $thing->acc;
    my ($db,$id) = split("::",$db_and_id);

    my $dbmaster = DBMaster->new(-database =>'Ontology',
				-host     => $WebConfig::DBHOST,
				-user     => $WebConfig::DBUSER,
				-password => $WebConfig::DBPWD);
    
    my ($name_title,$name_value,$description_title,$description_value);

    if($db =~ /PFAM/){
	my $new_id;
	if ($id =~ /_/){
	    ($new_id) = ($id) =~ /(.*?)_/;
	}
	else{
	    $new_id = $id;
	}

        my $pfam_objs = $dbmaster->pfam->get_objects( { 'id' => $new_id } );
        if(!scalar(@$pfam_objs)){
            $name_title = "name";
            $name_value = "not available";
            $description_title = "description";
            $description_value = "not available";
        }
        else{
            my $pfam_obj = $pfam_objs->[0];
	    $name_title = "name";
            $name_value = $pfam_obj->term;	    
            #$description_title = "description";
            #$description_value = $pfam_obj->description;
        }
    }
    
    my $short_title = $thing->acc;
    $short_title =~ s/::/ - /ig;
    my $new_short_title=$short_title;
    if ($short_title =~ /interpro/){
	($new_short_title) = ($short_title) =~ /(.*?)_/;
    }
    my $line_config = { 'title' => $name_value,
			'hover_title', => 'Domain',
                        'short_title' => $new_short_title,
                        'basepair_offset' => '1' };
    
    my $name;
    my ($new_id) = ($id) =~ /(.*?)_/;
    $name = {"title" => $db,
	     "value" => $new_id};
    push(@$descriptions,$name);

#    my $description;
#    $description = {"title" => $description_title,
#		    "value" => $description_value};
#    push(@$descriptions,$description);
    
    my $score;
    $score = {"title" => "score",
	      "value" => $thing->evalue};
    push(@$descriptions,$score);
    
    my $location;
    $location = {"title" => "location",
		 "value" => $thing->start . " - " . $thing->stop};
    push(@$descriptions,$location);

    my $link_id;
    if ($thing->acc =~/::(.*)/){
	$link_id = $1;
    }
    
    my $link;
    my $link_url;
#    if ($thing->class eq "CDD"){$link_url = "http://0-www.ncbi.nlm.nih.gov.library.vu.edu.au:80/Structure/cdd/cddsrv.cgi?uid=$link_id"}
    if($thing->class eq "PFAM"){$link_url = "http://pfam.sanger.ac.uk/family?acc=$link_id"}
    else{$link_url = "NO_URL"}
    
    $link = {"link_title" => $thing->acc,
	     "link" => $link_url};
    push(@$links_list,$link);
    
    my $element_hash = {
	"title" => $name_value,
	"start" => $thing->start,
	"end" =>  $thing->stop,
	"color"=> $color,
	"zlayer" => '2',
	"links_list" => $links_list,
	"description" => $descriptions};
    
    push(@$line_data,$element_hash);
    $gd->add_line($line_data, $line_config);
    
    return $gd;

}

sub display_table {
    my ($self,$dataset) = @_;
    my $cgi = new CGI;
    my $data = [];
    my $count = 0;
    my $content;
    my $seen = {};

    foreach my $thing (@$dataset) {
	next if ($thing->type !~ /dom/);
        my $single_domain = [];
	$count++;

	my $db_and_id = $thing->acc;
	my ($db,$id) = split("::",$db_and_id);

	my $dbmaster = DBMaster->new(-database =>'Ontology',
				-host     => $WebConfig::DBHOST,
				-user     => $WebConfig::DBUSER,
				-password => $WebConfig::DBPWD);

	my ($name_title,$name_value,$description_title,$description_value);

	my $new_id;
	if($db =~ /PFAM/){
	    if ($id =~ /_/){
		($new_id) = ($id) =~ /(.*?)_/;
	    }
	    else{
		$new_id = $id;
	    }

	    next if ($seen->{$new_id});
	    $seen->{$new_id}=1;

	    my $pfam_objs = $dbmaster->pfam->get_objects( { 'id' => $new_id } );
#	    print STDERR "VALUES: " . $pfam_objs . "\n";
	    if(!scalar(@$pfam_objs)){
		$name_title = "name";
		$name_value = "not available";
		$description_title = "description";
		$description_value = "not available";
	    }
	    else{
		my $pfam_obj = $pfam_objs->[0];
		$name_title = "name";
		$name_value = $pfam_obj->term;
		#$description_title = "description";
		#$description_value = $pfam_obj->description;
	    }
	}

	my $location =  $thing->start . " - " . $thing->stop;
	
        push(@$single_domain,$db);
	push(@$single_domain,$new_id);
	push(@$single_domain,$name_value);
	push(@$single_domain,$location);
	push(@$single_domain,$thing->evalue);
        push(@$single_domain,$description_value);
        push(@$data,$single_domain);
    }

    if ($count >0){
        $content = $data;
    }
    else
    {
        $content = "<p>This PEG does not have any similarities to domains</p>";
    }
}

 
#########################################
#########################################
package Observation::Location;

use base qw(Observation);

sub new {

    my ($class,$dataset) = @_; 
    my $self = $class->SUPER::new($dataset);
    $self->{cleavage_prob} = $dataset->{'cleavage_prob'};
    $self->{cleavage_loc} = $dataset->{'cleavage_loc'};
    $self->{signal_peptide_score} = $dataset->{'signal_peptide_score'};
    $self->{cello_location} = $dataset->{'cello_location'};
    $self->{cello_score} = $dataset->{'cello_score'};	
    $self->{tmpred_score} = $dataset->{'tmpred_score'};
    $self->{tmpred_locations} = $dataset->{'tmpred_locations'};	
    $self->{phobius_signal_location} = $dataset->{'phobius_signal_location'};	
    $self->{phobius_tm_locations} = $dataset->{'phobius_tm_locations'};	
    
    bless($self,$class);
    return $self;
}

sub display_cello {
    my ($thing) = @_;
    my $html;
    my $cello_location = $thing->cello_location;
    my $cello_score = $thing->cello_score;
    if($cello_location){
	$html .= "<p><font type=verdana size=-2>Subcellular location  prediction: $cello_location, score: $cello_score</font> </p>";   
	#$html .= "<p>CELLO score: $cello_score </p>";   
    }
    return ($html);
}

sub display {
    my ($thing,$gd,$fig) = @_;
    
    my $fid = $thing->fig_id;
    #my $fig= new FIG;
    my $length = length($fig->get_translation($fid));
    
    my $cleavage_prob;
    if($thing->cleavage_prob){$cleavage_prob = $thing->cleavage_prob;}
    my ($cleavage_loc_begin,$cleavage_loc_end) = split("-",$thing->cleavage_loc);
    my $signal_peptide_score = $thing->signal_peptide_score;
    my $cello_location = $thing->cello_location;
    my $cello_score = $thing->cello_score;
    my $tmpred_score = $thing->tmpred_score;
    my @tmpred_locations = split(",",$thing->tmpred_locations);
    
    my $phobius_signal_location = $thing->phobius_signal_location;
    my @phobius_tm_locations = split(",",$thing->phobius_tm_locations);
 
    my $lines = [];
    
    #color is 
    my $color = "6";



#    if($cello_location){
#	my $cello_descriptions = [];
#	my $line_data =[];
#
#	my $line_config = { 'title' => 'Localization Evidence',
#			    'short_title' => 'CELLO',
#                            'hover_title' => 'Localization',
#			    'basepair_offset' => '1' };
#
#	my $description_cello_location = {"title" => 'Best Cello Location',
#					  "value" => $cello_location};
#	
#	push(@$cello_descriptions,$description_cello_location);
#	
#	my $description_cello_score = {"title" => 'Cello Score',
#				       "value" => $cello_score};
#    
#	push(@$cello_descriptions,$description_cello_score);
#    
#	my $element_hash = {
#	    "title" => "CELLO",
#	    "color"=> $color,
#	    "start" => "1",
#	    "end" =>  $length + 1,
#	    "zlayer" => '1',
#	    "description" => $cello_descriptions};
#	
#	push(@$line_data,$element_hash);
#	$gd->add_line($line_data, $line_config);
#    }
#    
#    $color = "2";
#    if($tmpred_score){
#	my $line_data =[];
#	my $line_config = { 'title' => 'Localization Evidence',
#			    'short_title' => 'Transmembrane',
#			    'basepair_offset' => '1' };
#
#	foreach my $tmpred (@tmpred_locations){
#	    my $descriptions = [];
#	    my ($begin,$end) =split("-",$tmpred);
#	    my $description_tmpred_score = {"title" => 'TMPRED score',
#			     "value" => $tmpred_score};
#	
#	    push(@$descriptions,$description_tmpred_score);
#	    
#	    my $element_hash = {
#	    "title" => "transmembrane location",
#	    "start" => $begin + 1,
#	    "end" =>  $end + 1,
#	    "color"=> $color,
#	    "zlayer" => '5',
#	    "type" => 'box',
#	    "description" => $descriptions};
#	    
#	    push(@$line_data,$element_hash);
#
#	}
#	$gd->add_line($line_data, $line_config);
#    }


    if((scalar(@phobius_tm_locations) > 0) || $phobius_signal_location){
	my $line_data =[];
	my $line_config = { 'title' => 'Localization Evidence, Transmembrane and Signal Peptide',
			    'short_title' => 'TM and SP',
                            'hover_title' => 'Localization',
			    'basepair_offset' => '1' };

	foreach my $tm_loc (@phobius_tm_locations){
	    my $descriptions = [];
	    my $description_phobius_tm_locations = {"title" => 'transmembrane location',
			     "value" => $tm_loc};
	    push(@$descriptions,$description_phobius_tm_locations);
	    
	    my ($begin,$end) =split("-",$tm_loc);
	    
	    my $element_hash = {
	    "title" => "Phobius",
	    "start" => $begin + 1,
	    "end" =>  $end + 1,
	    "color"=> '6',
	    "zlayer" => '4',
	    "type" => 'bigbox',
	    "description" => $descriptions};
	    
	    push(@$line_data,$element_hash);

	}
	
	if($phobius_signal_location){
	    my $descriptions = [];
	    my $description_phobius_signal_location = {"title" => 'Phobius Signal Location',
			     "value" => $phobius_signal_location};
	    push(@$descriptions,$description_phobius_signal_location);
	    

	    my ($begin,$end) =split("-",$phobius_signal_location);
	    my $element_hash = {
	    "title" => "phobius signal locations",
	    "start" => $begin + 1,
	    "end" =>  $end + 1,
	    "color"=> '1',
	    "zlayer" => '5',
	    "type" => 'box',
	    "description" => $descriptions};
	    push(@$line_data,$element_hash);
	}
	
	$gd->add_line($line_data, $line_config);
    }


#    $color = "1";
#    if($signal_peptide_score){
#	my $line_data = [];
#	my $descriptions = [];
#
#	my $line_config = { 'title' => 'Localization Evidence',
#			    'short_title' => 'SignalP',
#                            'hover_title' => 'Localization',
#			    'basepair_offset' => '1' };
#
#	my $description_signal_peptide_score = {"title" => 'signal peptide score',
#						"value" => $signal_peptide_score};
#	
#	push(@$descriptions,$description_signal_peptide_score);
#
#	my $description_cleavage_prob = {"title" => 'cleavage site probability',
#					 "value" => $cleavage_prob};
#	
#	push(@$descriptions,$description_cleavage_prob);
#	    
#	my $element_hash = {
#	    "title" => "SignalP",
#	    "start" => $cleavage_loc_begin - 2,
#	    "end" =>  $cleavage_loc_end + 1,
#	    "type" => 'bigbox',
#	    "color"=> $color,
#	    "zlayer" => '10',
#	    "description" => $descriptions};
#	
#	push(@$line_data,$element_hash);
#	$gd->add_line($line_data, $line_config);
#    }


    return ($gd);

}

sub cleavage_loc {
  my ($self) = @_;

  return $self->{cleavage_loc};
}

sub cleavage_prob {
  my ($self) = @_;

  return $self->{cleavage_prob};
}

sub signal_peptide_score {
  my ($self) = @_;

  return $self->{signal_peptide_score};
}

sub tmpred_score {
  my ($self) = @_;

  return $self->{tmpred_score};
}

sub tmpred_locations {
  my ($self) = @_;

  return $self->{tmpred_locations};
}

sub cello_location {
  my ($self) = @_;

  return $self->{cello_location};
}

sub cello_score {
  my ($self) = @_;

  return $self->{cello_score};
}

sub phobius_signal_location {
  my ($self) = @_;
  return $self->{phobius_signal_location};
}

sub phobius_tm_locations {
  my ($self) = @_;
  return $self->{phobius_tm_locations};
}



#########################################
#########################################
package Observation::Sims;

use base qw(Observation);

sub new {

    my ($class,$dataset) = @_;
    my $self = $class->SUPER::new($dataset);
    $self->{identity} = $dataset->{'identity'};
    $self->{acc} = $dataset->{'acc'};
    $self->{query} = $dataset->{'query'};
    $self->{evalue} = $dataset->{'evalue'};
    $self->{qstart} = $dataset->{'qstart'};
    $self->{qstop} = $dataset->{'qstop'};
    $self->{hstart} = $dataset->{'hstart'};
    $self->{hstop} = $dataset->{'hstop'};
    $self->{database} = $dataset->{'database'};
    $self->{organism} = $dataset->{'organism'};
    $self->{function} = $dataset->{'function'};
    $self->{qlength} = $dataset->{'qlength'};
    $self->{hlength} = $dataset->{'hlength'};

    bless($self,$class);
    return $self;
}

=head3 display()

If available use the function specified here to display a graphical observation.
This code will display a graphical view of the similarities using the genome drawer object

=cut

sub display {
    my ($self,$gd,$thing,$fig,$base_start,$in_subs,$cgi) = @_;

    # declare variables
    my $window_size = $gd->window_size;
    my $peg = $thing->acc;
    my $query_id = $thing->query;
    my $organism = $thing->organism;
    my $abbrev_name = $fig->abbrev($organism);
    if (!$organism){
      $organism = $peg;
      $abbrev_name = $peg;
    }
    my $genome = $fig->genome_of($peg);
    my ($org_tax) = ($genome) =~ /(.*)\./;
    my $function = $thing->function;
    my $query_start = $thing->qstart;
    my $query_stop = $thing->qstop;
    my $hit_start = $thing->hstart;
    my $hit_stop = $thing->hstop;
    my $ln_query = $thing->qlength;
    my $ln_hit = $thing->hlength;
    my $query_color = match_color($query_start, $query_stop, abs($query_stop-$query_start)+1, 1);
    my $hit_color = match_color($hit_start, $hit_stop, abs($query_stop-$query_start)+1, 1);
    
    my $tax_link = "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" . $org_tax;

    # hit sequence title
    my $line_config = { 'title' => "$organism [$org_tax]",
			'short_title' => "$abbrev_name",
			'title_link' => '$tax_link',
			'basepair_offset' => '0',
			'no_middle_line' => '1'
			};
	    
    # query sequence title
    my $replace_id = $peg;
    $replace_id =~ s/\|/_/ig;
    my $anchor_name = "anchor_". $replace_id;
    my $query_config = { 'title' => "Query",
			 'short_title' => "Query",
			 'title_link' => "changeSimsLocation('$replace_id', 1)",
			 'basepair_offset' => '0',
			 'no_middle_line' => '1'
			 };
    my $line_data = [];
    my $query_data = [];
	    
    my $element_hash;
    my $hit_links_list = [];
    my $hit_descriptions = [];
    my $query_descriptions = [];
    
    # get sequence information
    # evidence link
    my $evidence_link;
    if ($peg =~ /^fig\|/){
      $evidence_link = "?page=Annotation&feature=".$peg;
    }
    else{
      my $db = &Observation::get_database($peg);
      my ($link_id) = ($peg) =~ /\|(.*)/;
      $evidence_link = &HTML::alias_url($link_id, $db);
      #print STDERR "LINK: $db    $evidence_link";
    }
    my $link = {"link_title" => $peg,
		"link" => $evidence_link};
    push(@$hit_links_list,$link) if ($evidence_link);
    
    # subsystem link
    my $subs = $in_subs->{$peg} if (defined $in_subs->{$peg});
    my @subsystems;
    foreach my $array (@$subs){
	my $subsystem = $$array[0];
	push(@subsystems,$subsystem);
	my $link = {"link" => "?page=Subsystems&subsystem=$subsystem",
		    "link_title" => $subsystem};
	push(@$hit_links_list,$link);
    }
    
    # blast alignment
    $link = {"link_title" => "view blast alignment",
	     "link" => "$FIG_Config::cgi_url/seedviewer.cgi?page=ToolResult&tool=bl2seq&peg1=$query_id&peg2=$peg"};
    push (@$hit_links_list,$link) if ($peg =~ /^fig\|/);
    
    # description data
    my $description_function;
    $description_function = {"title" => "function",
			     "value" => $function};
    push(@$hit_descriptions,$description_function);
    
    # subsystem description
    my $ss_string = join (",", @subsystems);
    $ss_string =~ s/_/ /ig;
    my $description_ss = {"title" => "subsystems",
			  "value" => $ss_string};
    push(@$hit_descriptions,$description_ss);
    
    # location description
    # hit
    my $description_loc;
    $description_loc = {"title" => "Hit Location",
			"value" => $hit_start . " - " . $hit_stop};
    push(@$hit_descriptions, $description_loc);
    
    $description_loc = {"title" => "Sequence Length",
			"value" => $ln_hit};
    push(@$hit_descriptions, $description_loc);

    # query
    $description_loc = {"title" => "Hit Location",
			"value" => $query_start . " - " . $query_stop};
    push(@$query_descriptions, $description_loc);
    
    $description_loc = {"title" => "Sequence Length",
			"value" => $ln_query};
    push(@$query_descriptions, $description_loc);

    
    
    # evalue score description
    my $evalue = $thing->evalue;
    while ($evalue =~ /-0/)
    {
	my ($chunk1, $chunk2) = split(/-/, $evalue);
	$chunk2 = substr($chunk2,1);
	$evalue = $chunk1 . "-" . $chunk2;
    }
    
    my $color = &color($evalue);
    my $description_eval = {"title" => "E-Value",
			    "value" => $evalue};
    push(@$hit_descriptions, $description_eval);
    push(@$query_descriptions, $description_eval);

    my $identity = $self->identity;
    my $description_identity = {"title" => "Identity",
				"value" => $identity};
    push(@$hit_descriptions, $description_identity);
    push(@$query_descriptions, $description_identity);
    

    my $number = $base_start + ($query_start-$hit_start);
    #print STDERR "START: $number";
    $element_hash = {
	"title" => $query_id,
	"start" => $base_start,
	"end" => $base_start+$ln_query,
	"type"=> 'box',
	"color"=> $color,
	"zlayer" => "2",
	"links_list" => $query_links_list,
	"description" => $query_descriptions
	};
    push(@$query_data,$element_hash);
    
    $element_hash = {
	"title" => $query_id . ': HIT AREA',
	"start" => $base_start + $query_start,
	"end" =>  $base_start + $query_stop,
	"type"=> 'smallbox',
	"color"=> $query_color,
	"zlayer" => "3",
	"links_list" => $query_links_list,
	"description" => $query_descriptions
	};
    push(@$query_data,$element_hash);
    
    $gd->add_line($query_data, $query_config);
    
    
    $element_hash = {
		"title" => $peg,
		"start" => $base_start + ($query_start-$hit_start),
		"end" => $base_start + (($query_start-$hit_start)+$ln_hit),
		"type"=> 'box',
		"color"=> $color,
		"zlayer" => "2",
		"links_list" => $hit_links_list,
		"description" => $hit_descriptions
		};
    push(@$line_data,$element_hash);
    
    $element_hash = {
	"title" => $peg . ': HIT AREA',
	"start" => $base_start + $query_start,
	"end" =>  $base_start + $query_stop,
	"type"=> 'smallbox',
	"color"=> $hit_color,
	"zlayer" => "3",
	"links_list" => $hit_links_list,
	"description" => $hit_descriptions
	};
    push(@$line_data,$element_hash);
    
    $gd->add_line($line_data, $line_config);
    
    my $breaker = [];
    my $breaker_hash = {};
    my $breaker_config = { 'no_middle_line' => "1" };
    
    push (@$breaker, $breaker_hash);
    $gd->add_line($breaker, $breaker_config);
    
    return ($gd);
}

=head3 display_domain_composition()

If available use the function specified here to display a graphical observation of the CDD(later Pfam or selected) domains that occur in the set of similar proteins

=cut

sub display_domain_composition {
    my ($self,$gd,$fig) = @_;
    
    #$fig = new FIG;
    my $peg = $self->acc;

    my $line_data = [];
    my $links_list = [];
    my $descriptions = [];
    
    my @domain_query_results =$fig->get_attributes($peg,"CDD");
    #my @domain_query_results = ();
    foreach $dqr (@domain_query_results){
	my $key = @$dqr[1];
	my @parts = split("::",$key);
	my $db = $parts[0];
	my $id = $parts[1];
	my $val = @$dqr[2];
	my $from;
	my $to;
	my $evalue;
	
	if($val =~/^(\d+\.\d+|0\.0);(\d+)-(\d+)/){
	    my $raw_evalue = $1;
	    $from = $2;
	    $to = $3;
	    if($raw_evalue =~/(\d+)\.(\d+)/){
		my $part2 = 1000 - $1;
		my $part1 = $2/100;
		$evalue = $part1."e-".$part2;
	    }
	    else{
		$evalue = "0.0";
	    }
	}

	my $dbmaster = DBMaster->new(-database =>'Ontology',
				-host     => $WebConfig::DBHOST,
				-user     => $WebConfig::DBUSER,
				-password => $WebConfig::DBPWD);
	my ($name_value,$description_value);

	if($db eq "CDD"){
	    my $cdd_objs = $dbmaster->cdd->get_objects( { 'id' => $id } );
	    if(!scalar(@$cdd_objs)){
		$name_title = "name";
		$name_value = "not available";
		$description_title = "description";
		$description_value = "not available";
	    }
	    else{
		my $cdd_obj = $cdd_objs->[0];
		$name_value = $cdd_obj->term;
		$description_value = $cdd_obj->description;
	    }
	}

	my $domain_name;
	$domain_name = {"title" => "name",
			"value" => $name_value};
	push(@$descriptions,$domain_name);

	my $description;
	$description = {"title" => "description",
			"value" => $description_value};
	push(@$descriptions,$description);
    
	my $score;
	$score = {"title" => "score",
		  "value" => $evalue};
	push(@$descriptions,$score);

	my $link_id = $id;
	my $link;
	my $link_url;
	if ($db eq "CDD"){$link_url = "http://0-www.ncbi.nlm.nih.gov.library.vu.edu.au:80/Structure/cdd/cddsrv.cgi?uid=$link_id"}
	elsif($db eq "PFAM"){$link_url = "http://pfam.sanger.ac.uk/family?acc=$link_id"}
	else{$link_url = "NO_URL"}
    
	$link = {"link_title" => $name_value,
		 "link" => $link_url};
	push(@$links_list,$link);
    
	my $domain_element_hash = {
	    "title" => $peg,
	    "start" => $from,
	    "end" =>  $to,
	    "type"=> 'box',
	    "zlayer" => '4',
	    "links_list" => $links_list,
	    "description" => $descriptions
	    };

	push(@$line_data,$domain_element_hash);

	#just one CDD domain for now, later will add option for multiple domains from selected DB
	last;
    }

    my $line_config = { 'title' => $peg,
			'hover_title' => 'Domain',
			'short_title' => $peg,
			'basepair_offset' => '1' };

    $gd->add_line($line_data, $line_config);
    
    return ($gd);
    
}

=head3 display_table()

If available use the function specified here to display the "raw" observation.
This code will display a table for the similarities protein

B<Please note> that URL linked to in display_method() is an external component and needs to added to the code for every class of evidence.

=cut

sub display_table {
    my ($self,$dataset, $show_columns, $query_fid, $fig, $application, $cgi) = @_;
    my ($count, $data, $content, %box_column, $subsystems_column, $evidence_column, %e_identical, $function_color, @ids);
    
    my $scroll_list;
    foreach my $col (@$show_columns){
	push (@$scroll_list, $col->{key});
    }

    push (@ids, $query_fid);
    foreach my $thing (@$dataset) {
	next if ($thing->class ne "SIM");
	push (@ids, $thing->acc);
    }

    $lineages = $fig->taxonomy_list() if (grep /lineage/, @$scroll_list);
    my @attributes = $fig->get_attributes(\@ids) if ( (grep /evidence/, @$scroll_list) || (grep /(pfam|mw)/, @$scroll_list) );

    # get the column for the subsystems
    $subsystems_column = &get_subsystems_column(\@ids,$fig,$cgi,'hash');

    # get the column for the evidence codes
    $evidence_column = &get_evidence_column(\@ids, \@attributes, $fig, $cgi, 'hash');

    # get the column for pfam_domain
    $pfam_column = &get_attrb_column(\@ids, \@attributes, $fig, $cgi, 'pfam', 'PFAM', 'hash') if (grep /^pfam$/, @$scroll_list);

    # get the column for molecular weight
    $mw_column = &get_attrb_column(\@ids, \@attributes, $fig, $cgi, 'mw', 'molecular_weight', 'hash') if (grep /^mw$/, @$scroll_list);

    # get the column for organism's habitat
    my $habitat_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'habitat', 'Habitat', 'hash') if (grep /^habitat$/, @$scroll_list);

    # get the column for organism's temperature optimum
    my $temperature_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'temperature', 'Optimal_Temperature', 'hash') if (grep /^temperature$/, @$scroll_list);
    
    # get the column for organism's temperature range
    my $temperature_range_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'temp_range', 'Temperature_Range', 'hash') if (grep /^temp_range$/, @$scroll_list);

    # get the column for organism's oxygen requirement
    my $oxygen_req_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'oxygen', 'Oxygen_Requirement', 'hash') if (grep /^oxygen$/, @$scroll_list);

    # get the column for organism's pathogenicity
    my $pathogenic_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'pathogenic', 'Pathogenic', 'hash') if (grep /^pathogenic$/, @$scroll_list);

    # get the column for organism's pathogenicity host
    my $pathogenic_in_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'pathogenic_in', 'Pathogenic_In', 'hash') if (grep /^pathogenic_in$/, @$scroll_list);

    # get the column for organism's salinity
    my $salinity_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'salinity', 'Salinity', 'hash') if (grep /^salinity$/, @$scroll_list);

    # get the column for organism's motility
    my $motility_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'motility', 'Motility', 'hash') if (grep /^motility$/, @$scroll_list);

    # get the column for organism's gram stain
    my $gram_stain_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'gram_stain', 'Gram_Stain', 'hash') if (grep /^gram_stain$/, @$scroll_list);

    # get the column for organism's endospores
    my $endospores_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'endospores', 'Endospores', 'hash') if (grep /^endospores$/, @$scroll_list);

    # get the column for organism's shape
    my $shape_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'shape', 'Shape', 'hash') if (grep /^shape$/, @$scroll_list);

    # get the column for organism's disease
    my $disease_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'disease', 'Disease', 'hash') if (grep /^disease$/, @$scroll_list);

    # get the column for organism's disease
    my $gc_content_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'gc_content', 'GC_Content', 'hash') if (grep /^gc_content$/, @$scroll_list);

    # get the column for transmembrane domains
    my $transmembrane_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'transmembrane', 'Phobius::transmembrane', 'hash') if (grep /^transmembrane$/, @$scroll_list);

    # get the column for similar to human
    my $similar_to_human_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'similar_to_human', 'similar_to_human', 'hash') if (grep /^similar_to_human$/, @$scroll_list);

    # get the column for signal peptide
    my $signal_peptide_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'signal_peptide', 'Phobius::signal', 'hash') if (grep /^signal_peptide$/, @$scroll_list);

    # get the column for transmembrane domains
    my $isoelectric_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'isoelectric', 'isoelectric_point', 'hash') if (grep /^isoelectric$/, @$scroll_list);

    # get the column for conserved neighborhood
    my $cons_neigh_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'conserved_neighborhood', undef, 'hash') if (grep /^conserved_neighborhood$/, @$scroll_list);

    # get the column for cellular location
    my $cell_location_column = &get_attrb_column(\@ids, undef, $fig, $cgi, 'cellular_location', 'PSORT::', 'hash') if (grep /^isoelectric$/, @$scroll_list);

    # get the aliases
    my $alias_col;
    if ( (grep /asap_id/, @$scroll_list) || (grep /ncbi_id/, @$scroll_list) ||
         (grep /refseq_id/, @$scroll_list) || (grep /swissprot_id/, @$scroll_list) ||
         (grep /uniprot_id/, @$scroll_list) || (grep /tigr_id/, @$scroll_list) ||
         (grep /kegg_id/, @$scroll_list) || (grep /pir_id/, @$scroll_list) ||
         (grep /trembl_id/, @$scroll_list) || (grep /jgi_id/, @$scroll_list) ) {
	$alias_col = &get_db_aliases(\@ids,$fig,'all',$cgi,'hash');
    }

    # get the colors for the function cell
    my $functions = $fig->function_of_bulk(\@ids,1);
    $functional_color = &get_function_color_cell($functions, $fig);
    my $query_function = $fig->function_of($query_fid);

    my %e_identical = &get_essentially_identical($query_fid,$dataset,$fig);

    my $figfam_data = &FIG::get_figfams_data();
    my $figfams = new FFs($figfam_data);
    my $same_genome_flag = 0;

    my $func_color_offset=0;
    unshift(@$dataset, $query_fid);
    for (my $thing_count=0;$thing_count<scalar @$dataset;$thing_count++){
#    foreach my $thing ( @$dataset){
	my $thing = $dataset->[$thing_count];
	my $next_thing = $dataset->[$thing_count+1] if (defined $dataset->[$thing_count+1]);
	my ($id, $taxid, $iden, $ln1,$ln2,$b1,$b2,$e1,$e2,$d1,$d2,$color1,$color2,$reg1,$reg2, $next_org);
	if ($thing eq $query_fid){
	    $id = $thing;
	    $taxid   = $fig->genome_of($id);
	    $organism = $fig->genus_species($taxid);
	    $current_function = $fig->function_of($id);
	}
	else{
	    next if ($thing->class ne "SIM");

	    $id      = $thing->acc;
	    $evalue  = $thing->evalue;
	    $taxid   = $fig->genome_of($id);
	    $iden    = $thing->identity;
	    $organism= $thing->organism;
	    $ln1     = $thing->qlength;
	    if ($ln1 < 1) { $ln1 = 1; }
	    $ln2     = $thing->hlength;
	    if ($ln2 < 1) { $ln2 = 1; }
	    $b1      = $thing->qstart;
	    $e1      = $thing->qstop;
	    $b2      = $thing->hstart;
	    $e2      = $thing->hstop;
	    $d1      = abs($e1 - $b1) + 1;
	    $d2      = abs($e2 - $b2) + 1;
	    $color1  = match_color( $b1, $e1, $ln1 );
	    $color2  = match_color( $b2, $e2, $ln2 );
	    $reg1    = {'data'=> "$b1-$e1 (<b>$d1/$ln1</b>)", 'highlight' => $color1};
	    $reg2    = {'data'=> "$b2-$e2 (<b>$d2/$ln2</b>)", 'highlight' => $color2};
	    $current_function = $thing->function;
	    $next_org = $next_thing->organism if (defined $next_thing);
	}

	next if ($id =~ /nmpdr\||gnl\|md5\|/);

        my $single_domain = [];
        $count++;

	# organisms cell
	my ($org, $org_color) = $fig->org_and_color_of($id);
	
	my $org_cell;
	if ( ($next_org ne $organism) && ($same_genome_flag == 0) ){
	    $org_cell = { 'data' =>  $organism, 'highlight' => $org_color};
	}
	elsif ($next_org eq $organism){
	    $org_cell = { 'data' =>  "<b>" . $organism . "</b>", 'highlight' => $org_color};
	    $same_genome_flag = 1;
	}
	elsif ($same_genome_flag == 1){
	    $org_cell = { 'data' =>  "<b>" . $organism . "</b>", 'highlight' => $org_color};
	    $same_genome_flag = 0;
	}

	# checkbox cell
	my ($box_cell,$tax, $radio_cell);
	my $field_name = "tables_" . $id;
	my $pair_name = "visual_" . $id;
	my $cell_name = "cell_". $id;
	my $replace_id = $id;
	$replace_id =~ s/\|/_/ig;
	my $white = '#ffffff';
	$white = '#999966' if ($id eq $query_fid);
	$org_color = '#999966' if ($id eq $query_fid);
	my $anchor_name = "anchor_". $replace_id;
	my $checked = ""; 
	#$checked = "checked" if ($id eq $query_fid);
#	if ($id =~ /^fig\|/){
	  my $box = qq~<a name="$anchor_name"></a><input type="checkbox" name="seq" value="$id" id="$field_name" onClick="VisualCheckPair('$field_name', '$pair_name','$cell_name');" $checked>~;
	  $box_cell = { 'data'=>$box, 'highlight'=>$org_color};
	  $tax = $fig->genome_of($id) if ($id =~ /^fig\|/);
#	}
#	else{
#	  my $box = qq(<a name="$anchor_name"></a>);
#	  $box_cell = { 'data'=>$box, 'highlight'=>$org_color};
#	}

	# create the radio cell for any sequence, not just fig ids
	my $radio = qq(<input type="radio" name="function_select" value="$current_function" id="$field_name" onClick="clearText('new_text_function')">);
	$radio_cell = { 'data'=>$radio, 'highlight'=>$white};

	# get the linked fig id
	my $anchor_link = "graph_" . $replace_id;

	my $fig_data;
	if ($id =~ /^fig\|/)
	{
	    $fig_data =  "<table><tr><td><a href='?page=Annotation&feature=$id'>$id</a></td>" . "&nbsp;" x 2;
	}
	else
	{
	    my $url_link = &HTML::set_prot_links($cgi,$id);
	    $fig_data = "<table><tr><td>$url_link</td>". "&nbsp;" x 2;
	}
	$fig_data .= qq(<td><img height='10px' width='20px' src='$FIG_Config::cgi_url/Html/anchor_alignment.png' alt='View Graphic View of Alignment' onClick='changeSimsLocation("$anchor_link", 0)'/></td></tr></table>);
	my $fig_col = {'data'=> $fig_data,
		       'highlight'=>$white};

	$replace_id = $peg;
	$replace_id =~ s/\|/_/ig;
	$anchor_name = "anchor_". $replace_id;
	my $query_config = { 'title' => "Query",
			     'short_title' => "Query",
			     'title_link' => "changeSimsLocation('$replace_id')",
			     'basepair_offset' => '0'
			     };

	# function cell
	my $function_cell_colors = {0=>"#ffffff", 1=>"#eeccaa", 2=>"#ffaaaa",
				    3=>"#ffcc66", 4=>"#ffff00", 5=>"#aaffaa",
				    6=>"#bbbbff", 7=>"#ffaaff", 8=>"#dddddd"};
	
	my $function_color;
	if ( (defined($functional_color->{$query_function})) && ($functional_color->{$query_function} == 1) ){
	    $function_color = $function_cell_colors->{ $functional_color->{$current_function} - $func_color_offset};
	}
	else{
	    $function_color = $function_cell_colors->{ $functional_color->{$current_function}};
	}
	my $function_cell;
	if ($current_function){
	  if ($current_function eq $query_function){
	    $function_cell = {'data'=>$current_function, 'highlight'=>$function_cell_colors->{0}};
	    $func_color_offset=1;
	  }
	  else{
	      $function_cell = {'data'=>$current_function,'highlight' => $function_color};
	  }
	}
	else{
	  $function_cell = {'data'=>$current_function,'highlight' => "#dddddd"};
	}
       
	if ($id eq $query_fid){
	    push (@$single_domain, $box_cell, {'data'=>qq~<i>Query Sequence: </i>~  . qq~<b>$id</b>~ , 'highlight'=>$white}, {'data'=> 'n/a', 'highlight'=>$white},
		  {'data'=>'n/a', 'highlight'=>$white}, {'data'=>'n/a', 'highlight'=>$white}, {'data'=>'n/a', 'highlight'=>$white},
		  {'data' =>  $organism, 'highlight'=> $white}, {'data'=>$current_function, 'highlight'=>$white},
		  {'data'=>$subsystems_column->{$id},'highlight'=>$white},
		  {'data'=>$evidence_column->{$id},'highlight'=>$white});  # permanent columns   
	}
	else{
	    push (@$single_domain, $box_cell, $fig_col, {'data'=> $evalue, 'highlight'=>"#ffffff"},
		  {'data'=>"$iden\%", 'highlight'=>"#ffffff"}, $reg1, $reg2, $org_cell, $function_cell,
		  {'data'=>$subsystems_column->{$id},'highlight'=>"#ffffff"},
		  {'data'=>$evidence_column->{$id},'highlight'=>"#ffffff"});  # permanent columns
	    
	}

	if ( ( $application->session->user) ){
	    my $user = $application->session->user;
	    if ($user && $user->has_right(undef, 'annotate', 'genome')) {
		push (@$single_domain,$radio_cell);
	    }
	}

	my ($ff) = $figfams->families_containing_peg($id);

	foreach my $col (@$scroll_list){
	    if ($id eq $query_fid) { $highlight_color = "#999966"; }
	    else { $highlight_color = "#ffffff"; }

	    if ($col =~ /pfam/)                       {push(@$single_domain,{'data'=>$pfam_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /mw/)                         {push(@$single_domain,{'data'=>$mw_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /habitat/)                    {push(@$single_domain,{'data'=>$habitat_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /temperature/)                {push(@$single_domain,{'data'=>$temperature_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /temp_range/)                 {push(@$single_domain,{'data'=>$temperature_range_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /oxygen/)                     {push(@$single_domain,{'data'=>$oxygen_req_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /^pathogenic$/)               {push(@$single_domain,{'data'=>$pathogenic_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /^pathogenic_in$/)            {push(@$single_domain,{'data'=>$pathogenic_in_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /salinity/)                   {push(@$single_domain,{'data'=>$salinity_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /motility/)                   {push(@$single_domain,{'data'=>$motility_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /gram_stain/)                 {push(@$single_domain,{'data'=>$gram_stain_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /endospores/)                 {push(@$single_domain,{'data'=>$endospores_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /shape/)                      {push(@$single_domain,{'data'=>$shape_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /disease/)                    {push(@$single_domain,{'data'=>$disease_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /gc_content/)                 {push(@$single_domain,{'data'=>$gc_content_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /transmembrane/)              {push(@$single_domain,{'data'=>$transmembrane_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /signal_peptide/)             {push(@$single_domain,{'data'=>$signal_peptide_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /isoelectric/)                {push(@$single_domain,{'data'=>$isoelectric_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /conerved_neighborhood/)     {push(@$single_domain,{'data'=>$cons_neigh_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /cellular_location/)          {push(@$single_domain,{'data'=>$cell_location_column->{$id},'highlight'=>$highlight_color});}
	    elsif ($col =~ /ncbi_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"NCBI"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /refseq_id/)                  {push(@$single_domain,{'data'=>$alias_col->{$id}->{"RefSeq"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /swissprot_id/)               {push(@$single_domain,{'data'=>$alias_col->{$id}->{"SwissProt"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /uniprot_id/)                 {push(@$single_domain,{'data'=>$alias_col->{$id}->{"UniProt"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /tigr_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"TIGR"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /pir_id/)                     {push(@$single_domain,{'data'=>$alias_col->{$id}->{"PIR"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /kegg_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"KEGG"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /trembl_id/)                  {push(@$single_domain,{'data'=>$alias_col->{$id}->{"TrEMBL"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /asap_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"ASAP"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /jgi_id/)                     {push(@$single_domain,{'data'=>$alias_col->{$id}->{"JGI"},'highlight'=>$highlight_color});}
	    elsif ($col =~ /lineage/)                   {push(@$single_domain,{'data'=>$lineages->{$tax},'highlight'=>$highlight_color});}
	    elsif ($col =~ /figfam/)                     {push(@$single_domain,{'data'=>"<a href='?page=FigFamViewer&figfam=" . $ff . "' target='_new'>" . $ff . "</a>",'highlight'=>$highlight_color});}
	}
        push(@$data,$single_domain);
    }
    if ($count >0 ){
	$content = $data;
    }
    else{
        $content = "<p>This PEG does not have any similarities</p>";
    }
    shift(@$dataset);
    return ($content);
}


=head3 display_figfam_table()

If available use the function specified here to display the "raw" observation.
This code will display a table for the similarities protein

B<Please note> that URL linked to in display_method() is an external component and needs to added to the code for every class of evidence.

=cut

sub display_figfam_table {
  my ($self,$ids, $show_columns, $fig, $application, $cgi) = @_;
  my ($count, $data, $content, %box_column, $subsystems_column, $evidence_column, %e_identical, $function_color, @ids);
  
  my $scroll_list;
  foreach my $col (@$show_columns){
    push (@$scroll_list, $col->{key});
  }
  
  $lineages = $fig->taxonomy_list() if (grep /lineage/, @$scroll_list);
  my @attributes = $fig->get_attributes($ids) if ( (grep /evidence/, @$scroll_list) || (grep /(pfam|mw)/, @$scroll_list) );
  
  # get the column for the subsystems
  $subsystems_column = &get_subsystems_column($ids,$fig,$cgi,'hash');
  
  # get the column for the evidence codes
  $evidence_column = &get_evidence_column($ids, \@attributes, $fig, $cgi, 'hash') if (grep /^evidence$/, @$scroll_list);
  
  # get the column for pfam_domain
  $pfam_column = &get_attrb_column($ids, \@attributes, $fig, $cgi, 'pfam', 'PFAM', 'hash') if (grep /^pfam$/, @$scroll_list);
  
  # get the column for molecular weight
  $mw_column = &get_attrb_column($ids, \@attributes, $fig, $cgi, 'mw', 'molecular_weight', 'hash') if (grep /^mw$/, @$scroll_list);
  
  # get the column for organism's habitat
  my $habitat_column = &get_attrb_column($ids, undef, $fig, $cgi, 'habitat', 'Habitat', 'hash') if (grep /^habitat$/, @$scroll_list);
  
  # get the column for organism's temperature optimum
  my $temperature_column = &get_attrb_column($ids, undef, $fig, $cgi, 'temperature', 'Optimal_Temperature', 'hash') if (grep /^temperature$/, @$scroll_list);
  
  # get the column for organism's temperature range
  my $temperature_range_column = &get_attrb_column($ids, undef, $fig, $cgi, 'temp_range', 'Temperature_Range', 'hash') if (grep /^temp_range$/, @$scroll_list);
  
  # get the column for organism's oxygen requirement
  my $oxygen_req_column = &get_attrb_column($ids, undef, $fig, $cgi, 'oxygen', 'Oxygen_Requirement', 'hash') if (grep /^oxygen$/, @$scroll_list);
  
  # get the column for organism's pathogenicity
  my $pathogenic_column = &get_attrb_column($ids, undef, $fig, $cgi, 'pathogenic', 'Pathogenic', 'hash') if (grep /^pathogenic$/, @$scroll_list);
  
  # get the column for organism's pathogenicity host
  my $pathogenic_in_column = &get_attrb_column($ids, undef, $fig, $cgi, 'pathogenic_in', 'Pathogenic_In', 'hash') if (grep /^pathogenic_in$/, @$scroll_list);
  
  # get the column for organism's salinity
  my $salinity_column = &get_attrb_column($ids, undef, $fig, $cgi, 'salinity', 'Salinity', 'hash') if (grep /^salinity$/, @$scroll_list);
  
  # get the column for organism's motility
  my $motility_column = &get_attrb_column($ids, undef, $fig, $cgi, 'motility', 'Motility', 'hash') if (grep /^motility$/, @$scroll_list);
  
  # get the column for organism's gram stain
  my $gram_stain_column = &get_attrb_column($ids, undef, $fig, $cgi, 'gram_stain', 'Gram_Stain', 'hash') if (grep /^gram_stain$/, @$scroll_list);

  # get the column for organism's endospores
  my $endospores_column = &get_attrb_column($ids, undef, $fig, $cgi, 'endospores', 'Endospores', 'hash') if (grep /^endospores$/, @$scroll_list);

  # get the column for organism's shape
  my $shape_column = &get_attrb_column($ids, undef, $fig, $cgi, 'shape', 'Shape', 'hash') if (grep /^shape$/, @$scroll_list);
  
  # get the column for organism's disease
  my $disease_column = &get_attrb_column($ids, undef, $fig, $cgi, 'disease', 'Disease', 'hash') if (grep /^disease$/, @$scroll_list);
  
  # get the column for organism's disease
  my $gc_content_column = &get_attrb_column($ids, undef, $fig, $cgi, 'gc_content', 'GC_Content', 'hash') if (grep /^gc_content$/, @$scroll_list);
  
  # get the column for transmembrane domains
  my $transmembrane_column = &get_attrb_column($ids, undef, $fig, $cgi, 'transmembrane', 'Phobius::transmembrane', 'hash') if (grep /^transmembrane$/, @$scroll_list);
  
  # get the column for similar to human
  my $similar_to_human_column = &get_attrb_column($ids, undef, $fig, $cgi, 'similar_to_human', 'similar_to_human', 'hash') if (grep /^similar_to_human$/, @$scroll_list);
  
  # get the column for signal peptide
  my $signal_peptide_column = &get_attrb_column($ids, undef, $fig, $cgi, 'signal_peptide', 'Phobius::signal', 'hash') if (grep /^signal_peptide$/, @$scroll_list);
  
  # get the column for transmembrane domains
  my $isoelectric_column = &get_attrb_column($ids, undef, $fig, $cgi, 'isoelectric', 'isoelectric_point', 'hash') if (grep /^isoelectric$/, @$scroll_list);
  
  # get the column for conserved neighborhood
  my $cons_neigh_column = &get_attrb_column($ids, undef, $fig, $cgi, 'conserved_neighborhood', undef, 'hash') if (grep /^conserved_neighborhood$/, @$scroll_list);

  # get the column for cellular location
  my $cell_location_column = &get_attrb_column($ids, undef, $fig, $cgi, 'cellular_location', 'PSORT::', 'hash') if (grep /^isoelectric$/, @$scroll_list);
  
  # get the aliases
  my $alias_col;
  if ( (grep /asap_id/, @$scroll_list) || (grep /ncbi_id/, @$scroll_list) ||
       (grep /refseq_id/, @$scroll_list) || (grep /swissprot_id/, @$scroll_list) ||
       (grep /uniprot_id/, @$scroll_list) || (grep /tigr_id/, @$scroll_list) ||
       (grep /kegg_id/, @$scroll_list) || (grep /pir_id/, @$scroll_list) ||
       (grep /trembl_id/, @$scroll_list) || (grep /jgi_id/, @$scroll_list) ) {
    $alias_col = &get_db_aliases($ids,$fig,'all',$cgi,'hash');
  }
  
  foreach my $id ( @$ids){
    my $current_function = $fig->function_of($id);
    my $organism = $fig->org_of($id);
    my $single_domain = [];
    
    # organisms cell comehere2
    my ($org, $org_color) = $fig->org_and_color_of($id);	
    my $org_cell = { 'data' =>  $organism, 'highlight' => $org_color};
    
    # get the linked fig id
    my $fig_data;
    if ($id =~ /^fig\|/)
    {
	$fig_data =  "<a href='?page=Annotation&feature=$id'>$id</a>";
    }
    else
    {
	my $url_link = &HTML::set_prot_links($cgi,$id);
	$fig_data = "<table><tr><td>$url_link</td>". "&nbsp;" x 2;
    }

    my $fig_col = {'data'=> $fig_data,
		   'highlight'=>"#ffffff"};

    # get sequence length
    my $length_col = {'data'=> $fig->translation_length($id),
		      'highlight'=>"#ffffff"};

    # function cell
    $function_cell = {'data'=>$current_function, 'highlight'=> "#ffffff"};
    
    # insert data
    push (@$single_domain, $fig_col, $length_col, $org_cell, {'data'=>$subsystems_column->{$id},'highlight'=>"#ffffff"}, $function_cell);
    
    foreach my $col (@$scroll_list){
      my $highlight_color = "#ffffff";
      
      if ($col =~ /evidence/)                   {push(@$single_domain,{'data'=>$evidence_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /pfam/)                       {push(@$single_domain,{'data'=>$pfam_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /mw/)                         {push(@$single_domain,{'data'=>$mw_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /habitat/)                    {push(@$single_domain,{'data'=>$habitat_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /temperature/)                {push(@$single_domain,{'data'=>$temperature_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /temp_range/)                 {push(@$single_domain,{'data'=>$temperature_range_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /oxygen/)                     {push(@$single_domain,{'data'=>$oxygen_req_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /^pathogenic$/)               {push(@$single_domain,{'data'=>$pathogenic_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /^pathogenic_in$/)            {push(@$single_domain,{'data'=>$pathogenic_in_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /salinity/)                   {push(@$single_domain,{'data'=>$salinity_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /motility/)                   {push(@$single_domain,{'data'=>$motility_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /gram_stain/)                 {push(@$single_domain,{'data'=>$gram_stain_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /endospores/)                 {push(@$single_domain,{'data'=>$endospores_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /shape/)                      {push(@$single_domain,{'data'=>$shape_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /disease/)                    {push(@$single_domain,{'data'=>$disease_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /gc_content/)                 {push(@$single_domain,{'data'=>$gc_content_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /transmembrane/)              {push(@$single_domain,{'data'=>$transmembrane_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /signal_peptide/)             {push(@$single_domain,{'data'=>$signal_peptide_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /isoelectric/)                {push(@$single_domain,{'data'=>$isoelectric_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /conerved_neighborhood/)     {push(@$single_domain,{'data'=>$cons_neigh_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /cellular_location/)          {push(@$single_domain,{'data'=>$cell_location_column->{$id},'highlight'=>$highlight_color});}
      elsif ($col =~ /ncbi_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"NCBI"},'highlight'=>$highlight_color});}
      elsif ($col =~ /refseq_id/)                  {push(@$single_domain,{'data'=>$alias_col->{$id}->{"RefSeq"},'highlight'=>$highlight_color});}
      elsif ($col =~ /swissprot_id/)               {push(@$single_domain,{'data'=>$alias_col->{$id}->{"SwissProt"},'highlight'=>$highlight_color});}
      elsif ($col =~ /uniprot_id/)                 {push(@$single_domain,{'data'=>$alias_col->{$id}->{"UniProt"},'highlight'=>$highlight_color});}
      elsif ($col =~ /tigr_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"TIGR"},'highlight'=>$highlight_color});}
      elsif ($col =~ /pir_id/)                     {push(@$single_domain,{'data'=>$alias_col->{$id}->{"PIR"},'highlight'=>$highlight_color});}
      elsif ($col =~ /kegg_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"KEGG"},'highlight'=>$highlight_color});}
      elsif ($col =~ /trembl_id/)                  {push(@$single_domain,{'data'=>$alias_col->{$id}->{"TrEMBL"},'highlight'=>$highlight_color});}
      elsif ($col =~ /asap_id/)                    {push(@$single_domain,{'data'=>$alias_col->{$id}->{"ASAP"},'highlight'=>$highlight_color});}
      elsif ($col =~ /jgi_id/)                     {push(@$single_domain,{'data'=>$alias_col->{$id}->{"JGI"},'highlight'=>$highlight_color});}
      elsif ($col =~ /lineage/)                   {push(@$single_domain,{'data'=>$lineages->{$tax},'highlight'=>$highlight_color});}
      elsif ($col =~ /figfam/)                     {push(@$single_domain,{'data'=>"<a href='?page=FigFamViewer&figfam=" . $ff . "' target='_new'>" . $ff . "</a>",'highlight'=>$highlight_color});}
    }
    push(@$data,$single_domain);
  }
 
  $content = $data;
  return ($content);
}

sub get_box_column{
    my ($ids) = @_;
    my %column;
    foreach my $id (@$ids){
        my $field_name = "tables_" . $id;
        my $pair_name = "visual_" . $id;
	my $cell_name = "cell_" . $id;
	$column{$id} = qq(<input type=checkbox name=seq value="$id" id="$field_name" onClick="VisualCheckPair('$field_name', '$pair_name', '$cell_name');">);
    }
    return (%column);
}

sub get_figfam_column{
    my ($ids, $fig, $cgi) = @_;
    my $column;

    my $figfam_data = &FIG::get_figfams_data();
    my $figfams = new FFs($figfam_data);

    foreach my $id (@$ids){
	my ($ff);
	if ($id =~ /\.peg\./){
	    ($ff) =  $figfams->families_containing_peg($id);
	}
	if ($ff){
	    push (@$column, "<a href='?page=FigFamViewer&figfam=" . $ff . "' target='_new'>" . $ff . "</a>");
	}
	else{
	    push (@$column, " ");
	}
    }

    return $column;
}

sub get_subsystems_column{
    my ($ids,$fig,$cgi,$returnType) = @_;

    my %in_subs  = $fig->subsystems_for_pegs($ids,1);
    my ($column, $ss);
    foreach my $id (@$ids){
	my @in_sub = @{$in_subs{$id}} if (defined $in_subs{$id});
	my @subsystems;
        if (scalar(@in_sub)) {
	    foreach my $array (@in_sub){
		my $ss_name = $array->[0];
		$ss_name =~ s/_/ /ig;
		push (@subsystems, "-" . $ss_name);
	    }
            my $in_sub_line = join ("<br>", @subsystems);
	    $ss->{$id} = $in_sub_line;
        } else {
            $ss->{$id} = "None added";
        }
	push (@$column, $ss->{$id});
    }

    if ($returnType eq 'hash') { return $ss; }
    elsif ($returnType eq 'array') { return $column; }
}

sub get_lineage_column{
    my ($ids, $fig, $cgi) = @_;

    my $lineages = $fig->taxonomy_list();

    foreach my $id (@$ids){
	my $genome = $fig->genome_of($id);
	if ($lineages->{$genome}){
#	    push (@$column, qq~<table style='border-style:hidden;'><tr><td style='background-color: #ffffff;'>~ . $lineages->{$genome} . qq~</td></tr</table>~);
	    push (@$column, $lineages->{$genome});
	}
	else{
	    push (@$column, " ");
	}
    }
    return $column;
}

sub match_color {
    my ( $b, $e, $n , $rgb) = @_;
    my ( $l, $r ) = ( $e > $b ) ? ( $b, $e ) : ( $e, $b );
    my $hue = 5/6 * 0.5*($l+$r)/$n - 1/12;
    my $cov = ( $r - $l + 1 ) / $n;
    my $sat = 1 - 10 * $cov / 9;
    my $br  = 1;
    if ($rgb){
	return html2rgb( rgb2html( hsb2rgb( $hue, $sat, $br ) ) );
    }
    else{
	rgb2html( hsb2rgb( $hue, $sat, $br ) );
    }
}

sub hsb2rgb {
    my ( $h, $s, $br ) = @_;
    $h = 6 * ($h - floor($h));
    if ( $s  > 1 ) { $s  = 1 } elsif ( $s  < 0 ) { $s  = 0 }
    if ( $br > 1 ) { $br = 1 } elsif ( $br < 0 ) { $br = 0 }
    my ( $r, $g, $b ) = ( $h <= 3 ) ? ( ( $h <= 1 ) ? ( 1,      $h,     0      )
                                      : ( $h <= 2 ) ? ( 2 - $h, 1,      0      )
                                      :               ( 0,      1,      $h - 2 )
                                      )
                                    : ( ( $h <= 4 ) ? ( 0,      4 - $h, 1      )
                                      : ( $h <= 5 ) ? ( $h - 4, 0,      1      )
                                      :               ( 1,      0,      6 - $h )
                                      );
    ( ( $r * $s + 1 - $s ) * $br,
      ( $g * $s + 1 - $s ) * $br,
      ( $b * $s + 1 - $s ) * $br
    )
}

sub html2rgb {
    my ($hex) = @_;
    my ($r,$g,$b) = ($hex) =~ /^\#(\w\w)(\w\w)(\w\w)/;
    my $code = { 'A'=>10, 'B'=>11, 'C'=>12, 'D'=>13, 'E'=>14, 'F'=>15, 
		 1=>1, 2=>2, 3=>3, 4=>4, 5=>5, 6=>6, 7=>7, 8=>8, 9=>9};

    my @R = split(//, $r);
    my @G = split(//, $g);
    my @B = split(//, $b);

    my $red = ($code->{uc($R[0])}*16)+$code->{uc($R[1])};
    my $green = ($code->{uc($G[0])}*16)+$code->{uc($G[1])};
    my $blue = ($code->{uc($B[0])}*16)+$code->{uc($B[1])};
    
    my $rgb = [$red, $green, $blue];
    return $rgb;

}

sub rgb2html {
    my ( $r, $g, $b ) = @_;
    if ( $r > 1 ) { $r = 1 } elsif ( $r < 0 ) { $r = 0 }
    if ( $g > 1 ) { $g = 1 } elsif ( $g < 0 ) { $g = 0 }
    if ( $b > 1 ) { $b = 1 } elsif ( $b < 0 ) { $b = 0 }
    sprintf("#%02x%02x%02x", int(255.999*$r), int(255.999*$g), int(255.999*$b) )
}

sub floor {
    my $x = $_[0];
    defined( $x ) || return undef;
    ( $x >= 0 ) || ( int($x) == $x ) ? int( $x ) : -1 - int( - $x )
}

sub get_function_color_cell{
  my ($functions, $fig) = @_;

  # figure out the quantity of each function
  my %hash;
  foreach my $key (keys %$functions){
    my $func = $functions->{$key};
    $hash{$func}++;
  }

  my %func_colors;
  my $count = 1;
  foreach my $key (sort {$hash{$b}<=>$hash{$a}} keys %hash){
    $func_colors{$key}=$count;
    $count++;
  }

  return \%func_colors;
}

sub get_essentially_identical{
    my ($fid,$dataset,$fig) = @_;
    #my $fig = new FIG;
    
    my %id_list;
    #my @maps_to = grep { $_ ne $fid and $_ !~ /^xxx/ } map { $_->[0] } $fig->mapped_prot_ids($fid);

    foreach my $thing (@$dataset){
	if($thing->class eq "IDENTICAL"){
	    my $rows = $thing->rows;
	    my $count_identical = 0;
	    foreach my $row (@$rows) {
		my $id = $row->[0];
		if (($id ne $fid) && ($fig->function_of($id))) {
		    $id_list{$id} = 1;
		}
	    }
	}
    }

#    foreach my $id (@maps_to) {
#        if (($id ne $fid) && ($fig->function_of($id))) {
#	    $id_list{$id} = 1;
#        }
#    }
    return(%id_list);
}


sub get_evidence_column{
    my ($ids,$attributes,$fig,$cgi,$returnType) = @_;
    my ($column, $code_attributes);

    if (! defined $attributes) {
	my @attributes_array = $fig->get_attributes($ids);
	$attributes = \@attributes_array;
    }

    my @codes = grep { $_->[1] =~ /^evidence_code/i } @$attributes;
    foreach my $key (@codes){
        push (@{$code_attributes->{$key->[0]}}, $key);
    }

    foreach my $id (@$ids){
	# add evidence code with tool tip
        my $ev_codes=" &nbsp; ";
	
	my @codes = @{$code_attributes->{$id}} if (@{$code_attributes->{$id}});
	my @ev_codes = ();
	foreach my $code (@codes) {
	    my $pretty_code = $code->[2];
	    if ($pretty_code =~ /;/) {
		my ($cd, $ss) = split(";", $code->[2]);
		if ($cd =~ /ilit|dlit/){
		    my ($type,$pubmed_id) = ($cd) =~ /(.*?)\((.*)\)/;
		    my $publink = &HTML::alias_url($pubmed_id,'PMID');
		    $cd = $type . "(<a href='" . $publink . "'>" . $pubmed_id . "</a>)";
		}
		$ss =~ s/_/ /g;
		$pretty_code = $cd;# . " in " . $ss;
	    }
	    push(@ev_codes, $pretty_code);
	}
	
        if (scalar(@ev_codes) && $ev_codes[0]) {
            my $ev_code_help=join("<br />", map {&HTML::evidence_codes_explain($_)} @ev_codes);
            $ev_codes = $cgi->a(
                                {
                                    id=>"evidence_codes", onMouseover=>"javascript:if(!this.tooltip) this.tooltip=new Popup_Tooltip(this, 'Evidence Codes', '$ev_code_help', ''); this.tooltip.addHandler(); return false;"}, join("<br />", @ev_codes));
        }
	
	if ($returnType eq 'hash') { $column->{$id}=$ev_codes; }
	elsif ($returnType eq 'array') { push (@$column, $ev_codes); }
    }
    return $column;
}

sub get_attrb_column{
    my ($ids, $attributes, $fig, $cgi, $colName, $attrbName, $returnType) = @_;

    my ($column, %code_attributes, %attribute_locations);
    my $dbmaster = DBMaster->new(-database =>'Ontology',
				 -host     => $WebConfig::DBHOST,
				 -user     => $WebConfig::DBUSER,
				 -password => $WebConfig::DBPWD);

    if ($colName eq "pfam"){    
	if (! defined $attributes) {
	    my @attributes_array = $fig->get_attributes($ids);
	    $attributes = \@attributes_array;
	}

	my @codes = grep { $_->[1] =~ /^$attrbName/i } @$attributes;
	foreach my $key (@codes){
	    my $name = $key->[1];
	    if ($name =~ /_/){
		($name) = ($key->[1]) =~ /(.*?)_/;
	    }
	    push (@{$code_attributes{$key->[0]}}, $name);
	    push (@{$attribute_location{$key->[0]}{$name}}, $key->[2]);
	}
	
	foreach my $id (@$ids){
	    # add pfam code
	    my $pfam_codes=" &nbsp; ";
	    my @pfam_codes = "";
	    my %description_codes;
	    
	    if ($id =~ /^fig\|\d+\.\d+\.peg\.\d+$/) {
		my @ncodes = @{$code_attributes{$id}} if (@{$code_attributes{$id}});
		@pfam_codes = ();
		
		# get only unique values
		my %saw;
		foreach my $key (@ncodes) {$saw{$key}=1;}
		@ncodes = keys %saw;
		
		foreach my $code (@ncodes) {
		    my @parts = split("::",$code);
		    my $pfam_link = "<a href=http://pfam.sanger.ac.uk/family?acc=" . $parts[1] . ">$parts[1]</a>";
		    
#		    # get the locations for the domain
#		    my @locs;
#		    foreach my $part (@{$attribute_location{$id}{$code}}){
#			my ($loc) = ($part) =~ /\;(.*)/;
#			push (@locs,$loc);
#		    }
#		    my %locsaw;
#		    foreach my $key (@locs) {$locsaw{$key}=1;}
#		    @locs = keys %locsaw;
#		    
#		    my $locations = join (", ", @locs);
#		    
		    if (defined ($description_codes{$parts[1]})){
			push(@pfam_codes, "$parts[1]");
		    }
		    else {
			my $description = $dbmaster->pfam->get_objects( { 'id' => $parts[1] } );
			$description_codes{$parts[1]} = $description->[0]->{term};
		        push(@pfam_codes, "$pfam_link");
		    }
		}
		
		if ($returnType eq 'hash') { $column->{$id} = join("<br><br>", @pfam_codes); }
	        elsif ($returnType eq 'array') { push (@$column, join("<br><br>", @pfam_codes)); }
	    }
	}
    }
    elsif ($colName eq 'cellular_location'){
	if (! defined $attributes) {
            my @attributes_array = $fig->get_attributes($ids);
            $attributes = \@attributes_array;
        }

	my @codes = grep { $_->[1] =~ /^$attrbName/i } @$attributes;
        foreach my $key (@codes){
	    my ($loc) = ($key->[1]) =~ /::(.*)/;
	    my ($new_loc, @all);
	    @all = split (//, $loc);
	    my $count = 0;
	    foreach my $i (@all){
		if ( ($i eq uc($i)) && ($count > 0) ){
		    $new_loc .= " " . $i;
		}
		else{
		    $new_loc .= $i;
		}
		$count++;
	    }
	    push (@{$code_attributes{$key->[0]}}, [$new_loc, $key->[2]]);
	}
	
	foreach my $id (@$ids){
            my (@values, $entry);
            #@values = (" ");
            if (@{$code_attributes{$id}}){
                my @ncodes = @{$code_attributes{$id}};
                foreach my $code (@ncodes){
                    push (@values, $code->[0] . ", " . $code->[1]);
                }
            }
            else{
                @values = ("Not available");
            }

	    if ($returnType eq 'hash') { $column->{$id} = join ("<BR>", @values); }
            elsif ($returnType eq 'array') { push (@$column, join ("<BR>", @values)); }
	}
    }
    elsif ( ($colName eq 'mw') || ($colName eq 'transmembrane') || ($colName eq 'similar_to_human') ||
	    ($colName eq 'signal_peptide') || ($colName eq 'isoelectric') ){
	if (! defined $attributes) {
	    my @attributes_array = $fig->get_attributes($ids);
	    $attributes = \@attributes_array;
	}

	my @codes = grep { $_->[1] =~ /^$attrbName/i } @$attributes;
        foreach my $key (@codes){
	    push (@{$code_attributes{$key->[0]}}, $key->[2]);
	}
	
	foreach my $id (@$ids){
	    my (@values, $entry);
	    #@values = (" ");
	    if (@{$code_attributes{$id}}){
		my @ncodes = @{$code_attributes{$id}};
		foreach my $code (@ncodes){
		    push (@values, $code);
		}
	    }
	    else{
		@values = ("Not available");
	    }
	    
	    if ($returnType eq 'hash') { $column->{$id} = join ("<BR>", @values); }
	    elsif ($returnType eq 'array') { push (@$column, join ("<BR>", @values)); }
	}
    }
    elsif ( ($colName eq 'habitat') || ($colName eq 'temperature') || ($colName eq 'temp_range') ||
	    ($colName eq 'oxygen') || ($colName eq 'pathogenic') || ($colName eq 'pathogenic_in') || 
	    ($colName eq 'salinity') || ($colName eq 'motility') || ($colName eq 'gram_stain') ||
	    ($colName eq 'endospores') || ($colName eq 'shape') || ($colName eq 'disease') || 
	    ($colName eq 'gc_content') ) {
	if (! defined $attributes) {
	    my @attributes_array = $fig->get_attributes(undef,$attrbName);
	    $attributes = \@attributes_array;
	}

	my $genomes_with_phenotype;
	foreach my $attribute (@$attributes){
	    my $genome = $attribute->[0];
	    $genomes_with_phenotype->{$genome} = $attribute->[2];
	}
	
	foreach my $id (@$ids){
	    my $genome = $fig->genome_of($id);
	    my @values = (' ');
	    if (defined $genomes_with_phenotype->{$genome}){
		push (@values, $genomes_with_phenotype->{$genome});
            }
	    if ($returnType eq 'hash') { $column->{$id} = join ("<BR>", @values); }
	    elsif ($returnType eq 'array') { push (@$column, join ("<BR>", @values)); }
	}
    }

    return $column;
}

sub get_aclh_aliases {
    my ($ids,$fig,$db,$cgi,$returnType) = @_;
    my $db_array;

    my $id_line = join (",", @$ids);
    my $aclh_url = "http://clearinghouse.nmpdr.org/aclh.cgi?page=SearchResults&raw_dump=1&query=" . $id_line;

    
}

sub get_id_aliases {
    my ($id, $fig) = @_;
    my $aliases = {};
    
    my $org = $fig->org_of($id);
    my $url = "http://clearinghouse.nmpdr.org/aclh.cgi?page=SearchResults&raw_dump=1&query=$id";
    if ( my $form = &LWP::Simple::get($url) ) {
	my ($block) = ($form) =~ /<pre>(.*)<\/pre>/s;
	foreach my $line (split /\n/, $block){
	    my @values = split /\t/, $line;
	    next if ($values[3] eq "Expert");
	    if (($values[1] =~ /$org/) || ($org =~ /$values[1]/) && (! defined $aliases->{$values[4]}) ){
		$aliases->{$values[4]} = $values[0];
	    }
	}
    }

    return $aliases;
}

sub get_db_aliases {
    my ($ids,$fig,$db,$cgi,$returnType) = @_;
    my $db_array;
    my $all_aliases = $fig->feature_aliases_bulk($ids);
    foreach my $id (@$ids){
#	my @all_aliases = grep { $_ ne $id and $_ !~ /^xxx/ } map { $_->[0] } $fig->mapped_prot_ids($id);
	my $id_org = $fig->org_of($id);
	
	foreach my $alias (@{$$all_aliases{$id}}){
#	foreach my $alias (@all_aliases){
	    my $id_db = &Observation::get_database($alias);
	    next if ( ($id_db ne $db) && ($db ne 'all') );
	    next if ($aliases->{$id}->{$db});
	    my $alias_org = $fig->org_of($alias);
#	    if (($id ne $peg) && ( ($alias_org =~ /$id_org/) || ($id_org =~ /$alias_org/)) ) {
		#push(@funcs, [$id,$id_db,$tmp]);
		$aliases->{$id}->{$id_db} = &HTML::set_prot_links($cgi,$alias);
#	    }
	}
	if (!defined( $aliases->{$id}->{$db})){
	    $aliases->{$id}->{$db} = " ";
	}
	#push (@$db_array, {'data'=>  $aliases->{$id}->{$db},'highlight'=>"#ffffff"});
	push (@$db_array, $aliases->{$id}->{$db});
    }
	
    if ($returnType eq 'hash') { return $aliases; }
    elsif ($returnType eq 'array') { return $db_array; }
}



sub html_enc { $_ = $_[0]; s/\&/&amp;/g; s/\>/&gt;/g; s/\</&lt;/g; $_ }

sub color {
    my ($evalue) = @_;
    my $palette = WebColors::get_palette('vitamins');
    my $color;
    if ($evalue <= 1e-170){        $color = $palette->[0];    }
    elsif (($evalue <= 1e-120) && ($evalue > 1e-170)){        $color = $palette->[1];    }
    elsif (($evalue <= 1e-90) && ($evalue > 1e-120)){        $color = $palette->[2];    }
    elsif (($evalue <= 1e-70) && ($evalue > 1e-90)){        $color = $palette->[3];    }
    elsif (($evalue <= 1e-40) && ($evalue > 1e-70)){        $color = $palette->[4];    }
    elsif (($evalue <= 1e-20) && ($evalue > 1e-40)){        $color = $palette->[5];    }
    elsif (($evalue <= 1e-5) && ($evalue > 1e-20)){        $color = $palette->[6];    }
    elsif (($evalue <= 1) && ($evalue > 1e-5)){        $color = $palette->[7];    }
    elsif (($evalue <= 10) && ($evalue > 1)){        $color = $palette->[8];    }
    else{        $color = $palette->[9];    }
    return ($color);
}


############################
package Observation::Cluster;

use base qw(Observation);

sub new {

    my ($class,$dataset) = @_;
    my $self = $class->SUPER::new($dataset);
    $self->{context} = $dataset->{'context'}; 
    bless($self,$class);
    return $self;
}

sub display {
    my ($self,$gd,$selected_taxonomies,$taxes,$sims_array,$fig) = @_;
    
    $taxes = $fig->taxonomy_list();

    my $fid = $self->fig_id;
    my $compare_or_coupling = $self->context;
    my $gd_window_size = $gd->window_size;
    my $range = $gd_window_size;
    my $all_regions = [];
    my $gene_associations={};

    #get the organism genome
    my $target_genome = $fig->genome_of($fid);
    $gene_associations->{$fid}->{"organism"} = $target_genome;
    $gene_associations->{$fid}->{"main_gene"} = $fid;
    $gene_associations->{$fid}->{"reverse_flag"} = 0;

    # get location of the gene
    my $data = $fig->feature_location($fid);
    my ($contig, $beg, $end);
    my %reverse_flag;

    if ($data =~ /(.*)_(\d+)_(\d+)$/){
	$contig = $1;
	$beg = $2;
	$end = $3;
    }

    my $offset;
    my ($region_start, $region_end);
    if ($beg < $end)
    {
	$region_start = $beg - ($range);
	$region_end = $end+ ($range);
	$offset = ($2+(($3-$2)/2))-($gd_window_size/2);
    }
    else
    {
	$region_start = $end-($range);
	$region_end = $beg+($range);
	$offset = ($3+(($2-$3)/2))-($gd_window_size/2);
	$reverse_flag{$target_genome} = $fid;
	$gene_associations->{$fid}->{"reverse_flag"} = 1;
    }

    # call genes in region
    my ($target_gene_features, $reg_beg, $reg_end) = $fig->genes_in_region($target_genome, $contig, $region_start, $region_end);
    #foreach my $feat (@$target_gene_features){
    #   push (@$all_regions, $feat) if ($feat =~ /peg/);
    #}
    push(@$all_regions,$target_gene_features);
    my (@start_array_region);
    push (@start_array_region, $offset);

    my %all_genes;
    my %all_genomes; 
    foreach my $feature (@$target_gene_features){ 
	#if ($feature =~ /peg/){
	    $all_genes{$feature} = $fid; $gene_associations->{$feature}->{"main_gene"}=$fid;
	#}
    }

    my @selected_sims;

    if ($compare_or_coupling eq "sims"){
	# get the selected boxes
	my @selected_taxonomy = @$selected_taxonomies;

	# get the similarities and store only the ones that match the lineages selected
	if (@selected_taxonomy > 0){
	    foreach my $sim (@$sims_array){
		next if ($sim->class ne "SIM");
		next if ($sim->acc !~ /fig\|/);

		#my $genome = $fig->genome_of($sim->[1]);
		my $genome = $fig->genome_of($sim->acc);
		#my ($genome1) = ($genome) =~ /(.*)\./;
		my $lineage = $taxes->{$genome};
		#my $lineage = $fig->taxonomy_of($fig->genome_of($genome));
		foreach my $taxon(@selected_taxonomy){
		    if ($lineage =~ /$taxon/){
			#push (@selected_sims, $sim->[1]);
			push (@selected_sims, $sim->acc);
		    }
		}
	    }
	}
	else{
	    my $simcount = 0;
	    foreach my $sim (@$sims_array){
		next if ($sim->class ne "SIM");
		next if ($sim->acc !~ /fig\|/);

		push (@selected_sims, $sim->acc);
		$simcount++;
		last if ($simcount > 4);
	    }
	}

	my %saw;
	@selected_sims = grep(!$saw{$_}++, @selected_sims);

	# get the gene context for the sorted matches
	foreach my $sim_fid(@selected_sims){
	    #get the organism genome
	    my $sim_genome = $fig->genome_of($sim_fid);
	    $gene_associations->{$sim_fid}->{"organism"} = $sim_genome;
	    $gene_associations->{$sim_fid}->{"main_gene"} = $sim_fid;
	    $gene_associations->{$sim_fid}->{"reverse_flag"} = 0;

	    # get location of the gene
	    my $data = $fig->feature_location($sim_fid);
	    my ($contig, $beg, $end);

	    if ($data =~ /(.*)_(\d+)_(\d+)$/){
		$contig = $1;
		$beg = $2;
		$end = $3;
	    }
	    
	    my $offset;
	    my ($region_start, $region_end);
	    if ($beg < $end)
	    {
		$region_start = $beg - ($range/2);
		$region_end = $end+($range/2);
		$offset = ($beg+(($end-$beg)/2))-($gd_window_size/2);
	    }
	    else
	    {
		$region_start = $end-($range/2);
		$region_end = $beg+($range/2);
		$offset = ($end+(($beg-$end)/2))-($gd_window_size/2);
		$reverse_flag{$sim_genome} = $sim_fid;
		$gene_associations->{$sim_fid}->{"reverse_flag"} = 1;
	    }

	    # call genes in region
	    my ($sim_gene_features, $reg_beg, $reg_end) = $fig->genes_in_region($sim_genome, $contig, $region_start, $region_end);
	    push(@$all_regions,$sim_gene_features);
	    push (@start_array_region, $offset);
	    foreach my $feature (@$sim_gene_features){ $all_genes{$feature} = $sim_fid;$gene_associations->{$feature}->{"main_gene"}=$sim_fid;}
	    $all_genomes{$sim_genome} = 1;
	}

    }

    #print STDERR "START CLUSTER OF GENES IN COMP REGION: " . `date`;
    # cluster the genes
    my @all_pegs = keys %all_genes;
    my $color_sets = &cluster_genes($fig,\@all_pegs,$fid);
    #print STDERR "END CLUSTER OF GENES IN COMP REGION: ". `date`;
    my %in_subs  = $fig->subsystems_for_pegs(\@all_pegs,1);

    foreach my $region (@$all_regions){
	my $sample_peg = @$region[0];
	my $region_genome = $fig->genome_of($sample_peg);
	my $region_gs = $fig->genus_species($region_genome);
	my $abbrev_name = $fig->abbrev($region_gs);
	#my ($genome1) = ($region_genome) =~ /(.*?)\./;
	my $lineage = $taxes->{$region_genome};
	#my $lineage = $fig->taxonomy_of($region_genome);
	#$region_gs .= "Lineage:$lineage";
        my $line_config = { 'title' => $region_gs,
			    'short_title' => $abbrev_name,
			    'basepair_offset' => '0'
			    };
	
	my $offsetting = shift @start_array_region;

	my $second_line_config = { 'title' => "$lineage",
				   'short_title' => "",
				   'basepair_offset' => '0',
				   'no_middle_line' => '1'
				   };
                
	my $line_data = [];
	my $second_line_data = [];

	# initialize variables to check for overlap in genes
	my ($prev_start, $prev_stop, $prev_fig, $second_line_flag);
	my $major_line_flag = 0;
	my $prev_second_flag = 0;

	foreach my $fid1 (@$region){
	    $second_line_flag = 0;
	    my $element_hash;
	    my $links_list = [];
	    my $descriptions = [];

	    my $color = $color_sets->{$fid1};

	    # get subsystem information
	    my $function = $fig->function_of($fid1);
	    my $url_link = "?page=Annotation&feature=".$fid1;

	    my $link;
            $link = {"link_title" => $fid1,
                     "link" => $url_link};
            push(@$links_list,$link);

	    my @subs = @{$in_subs{$fid1}} if (defined $in_subs{$fid1});
	    my @subsystems;
	    foreach my $array (@subs){
		my $subsystem = $$array[0];
		my $ss = $subsystem;
		$ss =~ s/_/ /ig;
		push (@subsystems, $ss);
		my $link;
		$link = {"link" => "?page=Subsystems&subsystem=$subsystem",
			 "link_title" => $ss};
		push(@$links_list,$link);
	    }

	    if ($fid1 eq $fid){
		my $link;
		$link = {"link_title" => "Annotate this sequence",
			 "link" => "$FIG_Config::cgi_url/seedviewer.cgi?page=Commentary"};
		push (@$links_list,$link);
	    }

            my $description_function;
            $description_function = {"title" => "function",
                                     "value" => $function};
            push(@$descriptions,$description_function);

            my $description_ss;
            my $ss_string = join (", ", @subsystems);
            $description_ss = {"title" => "subsystems",
                               "value" => $ss_string};
            push(@$descriptions,$description_ss);


	    my $fid_location = $fig->feature_location($fid1);
	    if($fid_location =~/(.*)_(\d+)_(\d+)$/){
		my($start,$stop);
		$start = $2 - $offsetting;
		$stop = $3 - $offsetting;
		
		if ( (($prev_start) && ($prev_stop) ) &&
		     ( ($start < $prev_start) || ($start < $prev_stop) ||
		       ($stop < $prev_start) || ($stop < $prev_stop) )){
		    if (($second_line_flag == 0) && ($prev_second_flag == 0)) { 
			$second_line_flag = 1; 
			$major_line_flag = 1;
		    }
		}
		$prev_start = $start;
		$prev_stop = $stop;
		$prev_fig = $fid1;
		
		if ((defined($reverse_flag{$region_genome})) && ($reverse_flag{$region_genome} eq $all_gnes{$fid1})){
		    $start = $gd_window_size - $start;
		    $stop = $gd_window_size - $stop;
		}
		
		my $title = $fid1;
		if ($fid1 eq $fid){
		    $title = "My query gene: $fid1";
		}

		$element_hash = {
		    "title" => $title,
		    "start" => $start,
		    "end" =>  $stop,
		    "type"=> 'arrow',
		    "color"=> $color,
		    "zlayer" => "2",
		    "links_list" => $links_list,
		    "description" => $descriptions
		};

		# if there is an overlap, put into second line
		if ($second_line_flag == 1){ push(@$second_line_data,$element_hash); $prev_second_flag = 1;}
		else{ push(@$line_data,$element_hash); $prev_second_flag = 0;}

		if ($fid1 eq $fid){
		    $element_hash = {
			"title" => 'Query',
			"start" => $start,
			"end" =>  $stop,
			"type"=> 'bigbox',
			"color"=> $color,
			"zlayer" => "1"
			};
		    
		    # if there is an overlap, put into second line
		    if ($second_line_flag == 1){ push(@$second_line_data,$element_hash); $prev_second_flag = 1;}
		    else{ push(@$line_data,$element_hash); $prev_second_flag = 0;}
		}
	    }
	}
	$gd->add_line($line_data, $line_config);
	$gd->add_line($second_line_data, $second_line_config); # if ($major_line_flag == 1);
    }
    return ($gd, \@selected_sims);
}

sub cluster_genes {
    my($fig,$all_pegs,$peg) = @_;
    my(%seen,$i,$j,$k,$x,$cluster,$conn,$pegI,$red_set);

    my @color_sets = ();

    $conn = &get_connections_by_similarity($fig,$all_pegs);

    for ($i=0; ($i < @$all_pegs); $i++) {
        if ($all_pegs->[$i] eq $peg) { $pegI = $i }
        if (! $seen{$i}) {
            $cluster = [$i];
            $seen{$i} = 1;
            for ($j=0; ($j < @$cluster); $j++) {
                $x = $conn->{$cluster->[$j]};
                foreach $k (@$x) {
                    if (! $seen{$k}) {
                        push(@$cluster,$k);
                        $seen{$k} = 1;
                    }
                }
            }

            if ((@$cluster > 1) || ($cluster->[0] eq $pegI)) {
                push(@color_sets,$cluster);
            }
        }
    }
    for ($i=0; ($i < @color_sets) && (! &in($pegI,$color_sets[$i])); $i++) {}
    $red_set = $color_sets[$i];
    splice(@color_sets,$i,1);
    @color_sets = sort { @$b <=> @$a } @color_sets;
    unshift(@color_sets,$red_set);

    my $color_sets = {};
    for ($i=0; ($i < @color_sets); $i++) {
        foreach $x (@{$color_sets[$i]}) {
            $color_sets->{$all_pegs->[$x]} = $i;
        }
    }
    return $color_sets;
}

sub get_connections_by_similarity {
    my($fig,$all_pegs) = @_;
    my($i,$j,$tmp,$peg,%pos_of);
    my($sim,%conn,$x,$y);

    for ($i=0; ($i < @$all_pegs); $i++) {
        $tmp = $fig->maps_to_id($all_pegs->[$i]);
        push(@{$pos_of{$tmp}},$i);
        if ($tmp ne $all_pegs->[$i]) {
            push(@{$pos_of{$all_pegs->[$i]}},$i);
        }
    }

    foreach $y (keys(%pos_of)) {
	$x = $pos_of{$y};
        for ($i=0; ($i < @$x); $i++) {
            for ($j=$i+1; ($j < @$x); $j++) {
                push(@{$conn{$x->[$i]}},$x->[$j]);
                push(@{$conn{$x->[$j]}},$x->[$i]);
            }
        }
    }

    for ($i=0; ($i < @$all_pegs); $i++) {
        foreach $sim ($fig->sims($all_pegs->[$i],500,10,"raw")) {
            if (defined($x = $pos_of{$sim->id2})) {
                foreach $y (@$x) {
                    push(@{$conn{$i}},$y);
                }
            }
        }
    }
    return \%conn;
}

sub in {
    my($x,$xL) = @_;
    my($i);

    for ($i=0; ($i < @$xL) && ($x != $xL->[$i]); $i++) {}
    return ($i < @$xL);
}

#############################################
#############################################
package Observation::Commentary;

use base qw(Observation);

=head3 display_protein_commentary()

=cut

sub display_protein_commentary {
    my ($self,$dataset,$mypeg,$fig) = @_;

    my $all_rows = [];
    my $content;
    #my $fig = new FIG;
    my $cgi = new CGI;
    my $count = 0;
    my $peg_array = [];
    my ($evidence_column, $subsystems_column,  %e_identical);

    if (@$dataset != 1){
	foreach my $thing (@$dataset){
	    if ($thing->class eq "SIM"){
		push (@$peg_array, $thing->acc);
	    }
	}
	# get the column for the evidence codes
	$evidence_column = &Observation::Sims::get_evidence_column($peg_array, undef, $fig, $cgi, 'hash');

	# get the column for the subsystems
        $subsystems_column = &Observation::Sims::get_subsystems_column($peg_array,$fig, $cgi, 'array');

	# get essentially identical seqs
        %e_identical = &Observation::Sims::get_essentially_identical($mypeg,$dataset,$fig);
    }
    else{
	push (@$peg_array, @$dataset);
    }
    
    my $selected_sims = [];
    foreach my $id (@$peg_array){
	last if ($count > 10);
	my $row_data = [];
	my ($set, $org, $ss, $ev, $function, $function_cell, $id_cell);
	if ($fig->org_of($id)){
	    $org = $fig->org_of($id);
	}
	else{
	    $org = "Data not available";
	}
	$function = $fig->function_of($id);
	if ($mypeg ne $id){
	    $function_cell = "<input type='radio' name='function' id='$id' value='$function' onClick=\"clearText('setAnnotation');\">&nbsp;&nbsp;$function";
	    $id_cell .= "<a href='?page=Annotation&feature=$id'>$id</a>"; # &HTML::set_prot_links($cgi,$id);
	    if (defined($e_identical{$id})) { $id_cell .= "*";}
	}
	else{
	    $function_cell = "&nbsp;&nbsp;$function";
	    $id_cell = "<input type='checkbox' name='peg' id='peg$count' value='$id' checked='true'>";
	    $id_cell .= "<a href='?page=Annotation&feature=$id'>$id</a>"; # &HTML::set_prot_links($cgi,$id);
	}

        push(@$row_data,$id_cell);
        push(@$row_data,$org);
	push(@$row_data, $subsystems_column->{$id}) if ($mypeg ne $id);
	push(@$row_data, $evidence_column->{$id}) if ($mypeg ne $id);
	push(@$row_data, $fig->translation_length($id));
        push(@$row_data,$function_cell);
        push(@$all_rows,$row_data);
	push (@$selected_sims, $id);
        $count++;	
    }

    if ($count >0){
        $content = $all_rows;
    }
    else{
        $content = "<p>This PEG does not have enough similarities to change the commentary</p>";
    }
    return ($content,$selected_sims);
}

sub display_protein_history {
    my ($self, $id,$fig) = @_;
    my $all_rows = [];
    my $content;

    my $cgi = new CGI;
    my $count = 0;
    foreach my $feat ($fig->feature_annotations($id)){
	my $row = [];
	my $col1 = $feat->[2];
	my $col2 = $feat->[1];
	#my $text = "<pre>" . $feat->[3] . "<\pre>";
	my $text = $feat->[3];

	push (@$row, $col1);
	push (@$row, $col2);
	push (@$row, $text);
	push (@$all_rows, $row);
	$count++;
    }
    if ($count > 0){
	$content = $all_rows;
    }
    else {
	$content = "There is no history for this PEG";
    }

    return($content);
}

