# _*_ Perl _*_
#
# 6/4/07
# Scenario.pm contains the package/class for creating scenario data instances
#  and completing analysis based on that information
#
package Scenario;

no warnings 'redefine';
use strict;
use FIGV;
use Subsystem;
use model;

my $fig;

#**************
# Constructors for Instances of Scenarios
#
#**************

sub new
{
    my ($class, $dir, $scenario_id, $scenario_name, $genome_id) = @_;
    my %self = ("dir" => $dir, "id" =>$scenario_id,
		"scenario_name" => $scenario_name,
		"substrates" => {}, "substrates_cofactors" => {},
		"products" => {}, "products_cofactors" => {},
		"reactions" => {}, "genome_id" => $genome_id
		); 
    bless \%self, $class;
}

#**************
# Static Functions
# 
#**************

#You must pass a FIG object which will help it locate the scenarios.
# needed to have rast/normal seed compatablity
sub set_fig
{
    my ($class, $figobj) = @_;
    $fig = $figobj;
    model::set_fig($fig);
}


# Arguments : Class , Genome_id
# Returns : Array of Scenario instances
# This function will create scenario instances for all scenarios of a given genome
# Example my @scenarios = @{Scenario->get_genome_scenarios("83333.1")};
sub get_genome_scenarios
{
	my ($class, $genome_id, $load_reactions) = @_;
	my $scenarios = create_scenario_objects($class, $genome_id, undef, $load_reactions);
	return $scenarios;
}

# Arguments : Class , Genome_id , Subsystem_name
# Returns : Array of Scenario instances
# This function will create scenario instances for scenarios in the given subsystem for the given genome
# Example my @scenarios = @{Scenario->get_genome_scenarios_for_subsystem("83333.1", "Glycolysis_and_Gluconeogenesis")};
sub get_genome_scenarios_for_subsystem
{
	my ($class, $genome_id, $ss_name, $load_reactions) = @_;
	my $ss_hash_ref = model::load_superset_file(); 
	my $scenarios = $class->create_scenario_objects($genome_id, "$ss_hash_ref->{$ss_name}/$ss_name", $load_reactions);
	return $scenarios;
}

# Arguments : Class , Genome_id , Subsystem_name , Scenario_name
# Returns : Array of Scenario instances
# This function will create scenario instances for the given scenario and genome
# Example my @scenarios = @{Scenario->get_genome_scenarios_for_scenario("83333.1", "Glycolysis_and_Gluconeogenesis", "Glucose_to_Pyruvate")};
sub get_genome_scenarios_for_scenario
{
	my ($class, $genome_id, $ss_name, $scen_name, $load_reactions) = @_;
	my $ss_hash_ref = model::load_superset_file(); 
	my $scenarios = $class->create_scenario_objects($genome_id, "$ss_hash_ref->{$ss_name}/$ss_name/$scen_name", $load_reactions);
	return $scenarios;
}

# Arguments : Class , Genome_id , Subsystem_name , Scenario_name, Path
# Returns : Scenario instance
# This function will create a scenario instance for the given path and genome
# Example my $scenario = Scenario->get_genome_scenarios_for_scenario("83333.1", "Glycolysis_and_Gluconeogenesis", "Glucose_to_Pyruvate", "path_1");
sub get_genome_scenarios_for_path
{
	my ($class, $genome_id, $ss_name, $scen_name, $path, $load_reactions) = @_;
	my $ss_hash_ref = model::load_superset_file(); 
	my $scenarios = $class->create_scenario_objects($genome_id, "$ss_hash_ref->{$ss_name}/$ss_name/$scen_name/$path", $load_reactions);
	return $scenarios->[0];
}

sub create_scenario_objects
{
    my ($class, $genome_id, $path_description, $load_reactions) = @_;

    my $scenario_dir = $fig->scenario_directory($genome_id);
    if (! -d $scenario_dir) {
	print STDERR "No Scenarios computed for $genome_id\n";
	return [];
    }
    my @scenario_obj;
    my @scenario_paths_dirs;

    if ($genome_id eq 'All' || (-d "$scenario_dir/PathInfo")) {
	my $dir_to_parse =  "/$genome_id/PathInfo/$path_description";
	@scenario_paths_dirs = @{model::parse_assembly_scenarios([$dir_to_parse])};
    }
    elsif (-e "$scenario_dir/path_info.txt") {
	open (PI, "$scenario_dir/path_info.txt");
	while (my $line = <PI>) {
	    next if (defined $path_description && $line !~ /$path_description/);
	    my @path_info = split "/", $line;
	    # get rid of path_info at end
	    pop @path_info;
	    push @scenario_paths_dirs, \@path_info;
	}
	close PI;
	# from here we'll look up path info in the "All" Scenarios structure
	$scenario_dir = $fig->scenario_directory("All");
    }
    else {
	print STDERR "No path info for $genome_id\n";
	return [];
    }

    foreach my $scenario_path (@scenario_paths_dirs)
    {
	my $path = "$scenario_dir/" . join "/" , @$scenario_path;
	my $scenario_id = "$scenario_path->[-3]/$scenario_path->[-2]/$scenario_path->[-1]";
	my $scenario_name = "$scenario_path->[-3]/$scenario_path->[-2]";
	my $scenario = $class->new($path,$scenario_id,$scenario_name,$genome_id);
	$scenario->load_information();
	$scenario->read_reactions() if($load_reactions);
	push @scenario_obj, $scenario;
    }
    return \@scenario_obj;
}

sub same_scenario
{
	my ($class, $scenario1, $scenario2) = @_;
    return ($scenario1->get_scenario_name() eq $scenario2->get_scenario_name());
    #return (split '/', $scenario1->get_id())[1] eq (split '/', $scenario2->get_id())[1];
}


#***************
# Data Manipulation (for instances of Scenarios only)
#
#***************

sub load_information
{
    my $self = shift;
    $self->read_data("/inputs", "substrates");
    $self->read_data("/outputs", "products");
    $self->read_path_info;
}

sub read_data
{
	my ($self, $file, $store_loc) = @_;
	
	open(MAINDATA, $self->{dir} . $file . "_main") or
	die("Failed to open " . $self->{dir} . $file . "_main" . " in " . $self->{id});
	my %main_cpds;
    while(<MAINDATA>)
    {
		chomp;
		my ($cid, $main) = split "\t";		
		if ($main)
		{
			$main_cpds{$cid} = 1;
		}
    }
    close(MAINDATA);
	open(DATA, $self->{dir} . $file) or
	die("Failed to open " . $self->{dir} . $file . " in " . $self->{id});
	while(<DATA>)
	{
		chomp;
		my ($cid, $stoich, $name) = split "\t";		
		if (exists $main_cpds{$cid})
		{
			$self->{$store_loc}->{$cid}->{stoich} = $stoich;
			$self->{$store_loc}->{$cid}->{name} = $name;
		} else {
			$self->{$store_loc."_cofactors"}->{$cid}->{stoich} = $stoich;
			$self->{$store_loc."_cofactors"}->{$cid}->{name} = $name;
		}
	}
	close(DATA);
}

sub read_reactions
{
	my $self = shift;
	open(REACTIONS,$self->{dir}."/reactions") or
	die("Failed to open ".$self->{dir}."/reactions file.\n");

    while(<REACTIONS>)
    {
	my $reaction_name;
	my %substrates;
	my %products;
	my $direction = 0; #by default reactions go left to right
	
	chomp;
	my @line = split "\t";
	chomp @line;
	#section 0 contains the reaction name
	if($line[0] =~ /^R\d{5}/)
	{
	    $reaction_name = $line[0];
	}
	else { next; } #This isn't a reaction, skip this line
	#section 1 contains compound and stoich information
	my @temp = split "=", $line[1];
	chomp @temp;
	my @substrates = split /\+/, $temp[0];
	my @products = split /\+/, $temp[1];
	chomp @substrates;
	chomp @products;
	foreach (@substrates)
	{
	    my @temp1 = split "\ ";
	    $substrates{$temp1[1]} = $temp1[0];
	}
	foreach (@products)
	{
	    my @temp1 = split "\ ";
	    $products{$temp1[1]} = $temp1[0];
	}
	#section 4 contains direction information
	my @temp2 = split "\ ", $line[4];
	if($temp2[0] eq "-Inf")
	{
	    $direction = 1; #Reversible direction
	}
	#Store the reaction..
	# Array: ($direction, %substrates, %products, genome_pegs)
	my @temp3 = ($direction,\%substrates,\%products,[]);
	$self->{reactions}->{$reaction_name} = \@temp3;
    }
    close(REACTIONS);
}

sub read_path_info
{
	my $self = shift;
    my $file = "/path_info";
    my $store_loc = "path_info";
    open(DATA,$self->{dir}.$file) or
	die("Failed to open ".$self->{dir}.$file." in ".$self->{id});
    while(<DATA>)
    {
	chomp;
	push @{$self->{$store_loc}}, $_;
    }
    close(DATA);
}

sub load_reaction_pegs
{
    my $self = shift;
    my $fig = shift;
    my @temp =  split "/", $self->{id};
    chomp @temp;
    my $reactions_for_ss = model::get_reactions_for_genome_in_subsystem_fast($self->{genome_id});
    foreach my $reaction (keys %$reactions_for_ss)
    {
	if(defined $self->{reactions}->{$reaction})
	{
	    push @{$self->{reactions}->{$reaction}->[3]} , @{$reactions_for_ss->{$reaction}};
	}
    }
}

#*********************
# Accessor Subroutines
#*********************

#main compounds only
sub get_substrates
{
    my $self = shift;
    return $self->{substrates};
}
#main compounds only
sub get_products
{
    my $self = shift;
    return $self->{products};
}
#reactions in the path only
sub get_path_info
{
    my $self = shift;
    return $self->{path_info};
}

sub get_substrate_cofactors
{
    my $self = shift;
    return $self->{substrates_cofactors};
}

sub get_product_cofactors
{
    my $self = shift;
    return $self->{products_cofactors};
}

#main and cofactor compounds
sub get_substrates_all
{
    my $self = shift;
    my $mains = $self->{substrates};
    my $cofactors = $self->{substrates_cofactors};
	map {$mains->{$_} = $cofactors->{$_}} keys %$cofactors;
    return $mains;
}
#main and cofactor compounds
sub get_products_all
{
    my $self = shift;
    my $mains = $self->{products};
    my $cofactors = $self->{products_cofactors};
	map {$mains->{$_} = $cofactors->{$_}} keys %$cofactors;
    return $mains;
}

sub get_id
{
    my $self = shift;
    return $self->{id};
}

sub get_scenario_name
{
    my $self = shift;
    return $self->{scenario_name};
}


sub get_scenario_name_only
{
    my $self = shift;
    my @parts = split "/" , $self->{scenario_name};
    chomp @parts;
    return $parts[1];
}

sub get_path_name
{
    my $self = shift;
    my @parts = split "/" , $self->{id};
    chomp @parts;
    return $parts[2];
}

sub get_subsystem_name
{
    my $self = shift;
    my @parts = split "/" , $self->{scenario_name};
    chomp @parts;
    return $parts[0];
}

sub get_reaction_ids
{
    my $self = shift;
    $self->read_reactions() unless scalar keys %{$self->{reactions}} > 0;
    my @names = keys %{$self->{reactions}};
    return \@names;
}

sub get_reaction_pegs
{
    my $self = shift;
    my $reaction_id = shift;
    return $self->{reactions}->{$reaction_id}->[3];
}


# Returns false if left-to-right, and true if reversible
sub get_reaction_reversibility
{
    my $self = shift;
    my $reaction_id = shift;
    return $self->{reactions}->{$reaction_id}->[0];
}

# Returns a hash of keys=compound_ids and values=stoich
sub get_reaction_substrates
{
    my $self = shift;
    my $reaction_id = shift;
    return $self->{reactions}->{$reaction_id}->[1];
}

# Returns a hash of keys=compound_ids and values=stoich
sub get_reaction_products
{
    my $self = shift;
    my $reaction_id = shift;
    return $self->{reactions}->{$reaction_id}->[2];
}

1;
