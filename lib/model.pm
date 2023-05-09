## _*_ Perl _*_ ##
#
# model.pm
#
# Kevin Formsma
# Hope College
# Created: 6/1/2006
#
##################

### BrainStorm/Planning ###
#
# Functions 
#   -Scenario access methods for subsystem
#   -Access current genome models
#   -Interface to/Integration of find_reaction_paths.cgi
#   -Scenario Relationships:  
#       1. Automanaged - map togather based on starting/ending compounds
#       2. User Defined - through some type of map or XML 
#   -Generate files for Scenarios on a organism
#       1. Path Picking - All or Select?
#   -Filter written Scenarios based on Relationships
#   -Report valid Scenarios for a genome
#
#
#  Basic Process for a Control Script
#    1. Input Genomes for model creation and Subsystems
#    2. Define selected Scenario relationships
#    3. Create Scenario Paths for each genome and subsystem
#    4. Filter invalid Scenarios based on relationships
#    5. Combine Subsystem models for each genome into genome models
#    6. Report on the model results for each genome
#

package model;

use strict;
no warnings 'redefine';
use Scenario;
use FIG;
use Subsystem;
use File::Path;
use Data::Dumper;

my $fig = eval("FIG->new()"); ## So that this will compile properly on non-FIG systems!

#global variables that make process_paths work
#These need to be cleared and reloaded frequently for process path and flux writing
my (%reactions_to_substrate_arrays, %reactions_to_product_arrays, %all_compounds_to_main);
my %all_reactions;
my %scenario_cycles;
my @all_outputs_lists;
my %all_inputs;
my %all_outputs;

#this variable is used to set the loop count for single scenario and then assembly runs, 100/25 by default.
my $loop_max = 100;
my $loop_max_assembly = 25;

#These store what supersets/Subsystems we are using
my %superset_to_ss;
my %ss_to_superset;

#Flip this bit to enable debugging
my $debug = int($ENV{HOPE_DEBUG});
$debug = 0;

sub set_fig
{
    my($newfig) = @_;
    $fig = $newfig;
}

sub new {
    my $type = shift;
    my $self = {};
    return bless $self, $type;
}

sub get_ss_scenarios 
{
    my ($ss_name) = @_;
    # This is the scenario data structure storage
    # %scenario_data{scenario_name}
    #
    my %scenario_data;

    my $subsystem = $fig->get_subsystem($ss_name);
    if (!$subsystem)
    {
	warn "Cannot open subsystem $subsystem\n";
	return;
    }
    
    my @scenario_names = $subsystem->get_hope_scenario_names;
    foreach my $name (@scenario_names)
    {
	$scenario_data{$name} = &get_scenario($subsystem,$name);
    }

    return \%scenario_data;
}


sub get_scenario
{
    my($subsystem,$name) = @_;

    my %data;

    my @inputs = $subsystem->get_hope_input_compounds($name);
    my @outputs = $subsystem->get_hope_output_compounds($name);
    my @map_ids = $subsystem->get_hope_map_ids($name);
    my @additional_reactions = $subsystem->get_hope_additional_reactions($name);
    my @ignore_reactions = $subsystem->get_hope_ignore_reactions($name);
    
    $data{inputs} = \@inputs;
    $data{outputs} = \@outputs;
    $data{map_ids} = \@map_ids;
    $data{additional_reactions} = \@additional_reactions;
    $data{ignore_reactions} = \@ignore_reactions;

    return \%data;

}

sub process_init
{
    my ($ss_name,$scenario_name,$genome,$assembly) = @_;
    my (%sc_inputs, %sc_outputs);

    if ($genome eq "")
    {
	$genome = "All";
    }
    print STDERR "\nSubsystem : ".$ss_name." Scenario: $scenario_name \n" if $debug;
    my $subsystem = $fig->get_subsystem($ss_name);
    my $scenario_data = &get_scenario($subsystem,$scenario_name);    

    #load the other arrays
    my %ss_reactions;    

    if ($genome eq "All")
    {
	my %all_reactions = $subsystem->get_hope_reactions;
	foreach my $role (keys %all_reactions)
	{
	    map { $ss_reactions{$_} = 1 } @{$all_reactions{$role}};
	}
    }
    else
    {
	my $reactions_for_genome = get_reactions_for_genome_in_subsystem($genome,$ss_name);
	map { $ss_reactions{$_} = 1 } keys %$reactions_for_genome if keys %$reactions_for_genome;
    }

    map { $sc_inputs{$_} = 1 } @{$scenario_data->{inputs}};
    
    foreach my $list (@{$scenario_data->{outputs}})
    {
	map { $sc_outputs{$_} = 1 } @$list;
	push @all_outputs_lists, $list;
    }

    map { $scenario_cycles{$_} = 1 if defined $sc_outputs{$_} } keys %sc_inputs;
    map { $all_inputs{$_} = 1 } keys %sc_inputs;
    map { $all_outputs{$_} = 1 } keys %sc_outputs;

    my @hope_additional_reactions = @{$scenario_data->{additional_reactions}};
    my @hope_ignore_reactions = @{$scenario_data->{ignore_reactions}};
    my %sc_reactions;
    map { $sc_reactions{$_} = 1 } keys %ss_reactions;
    
    # flag additional reactions so we won't check if they are in a map
    foreach my $rid (@hope_additional_reactions)
    {
	$sc_reactions{$rid} = 2;
    }
    
    foreach my $rid (@hope_ignore_reactions)
    {
	delete $sc_reactions{$rid};
    }

    #for now we do this outside of the if statment, but that might need to change
    &load_substrate_and_product_arrays(\%sc_reactions,$scenario_data->{map_ids}, \%sc_inputs);

    

    
}

sub execute_paths
{
    my ($assembly_paths,$find_first,$input_path, $output_path) = @_;

    my $num_paths = scalar @{$assembly_paths};

    my (%substrates_to_reactions, %products_to_reactions,
	%reactions_to_substrates, %reactions_to_products);

    if($num_paths  ==0 ) #Load the reactions normally if we are creating scenarios
    {
	create_reactions(\%substrates_to_reactions,\%products_to_reactions,
				  \%reactions_to_substrates,
				  \%reactions_to_products);
    }
    #create assembly reactions that are closed 'paths' from scenarios or assemblies
    if($num_paths  > 0)
    {
	create_assembly_reactions(\%substrates_to_reactions,\%products_to_reactions,
				  \%reactions_to_substrates,
				  \%reactions_to_products,$assembly_paths);
    }

    #This deals with user specifed input/output paths, and uses these to generate
    #what our input and output compounds should be for an assembly. 
    #This is only used for creating assemblies, has no effect on scenario creation

    if(scalar @{$input_path} && scalar @{$output_path})
    {

	%all_inputs = ();
	%all_outputs = ();
	#we should have loaded these paths above, so lets just pull out the information we need
	foreach my $path (@$input_path)
	{
	    my $input_rxn = "$path->[-3]/$path->[-2]/$path->[-1]_R";	    
	    my @user_in = @{$reactions_to_substrates{$input_rxn}};
	    map{ $all_inputs{$_} = 1 } @user_in;
	}
	foreach my $path (@$output_path)
	{
	    my $output_rxn = "$path->[-3]/$path->[-2]/$path->[-1]_R";
	    my @user_out = @{$reactions_to_products{$output_rxn}};
	    map{ $all_outputs{$_} = 1 } @user_out;
	}
	#map { $scenario_cycles{$_} = 1 if defined $all_outputs{$_} } keys %all_inputs;
    }

    print STDERR "Inputs :\n" if $debug;
    print STDERR map { $_."\n" } keys %all_inputs if $debug;
    print STDERR "Outputs:\n" if $debug;
    print STDERR map { $_."\n" } keys %all_outputs if $debug;

    #filter the input/outputs lists, removing the intersection unless something
    # is a known cycle
    foreach my $input (keys %all_inputs)
    {
	if(defined $all_outputs{$input} && ! defined $scenario_cycles{$input})
	{
	    print STDERR "Deleting $input from input and output lists\n" if $debug;
	    delete $all_inputs{$input};
	    delete $all_outputs{$input};
	}
    }

    my $create_assembly = 0;
    $create_assembly = 1 if(scalar @{$assembly_paths} !=0);
    return process_paths(\%all_inputs, \%all_outputs, \@all_outputs_lists,
			 \%reactions_to_substrates, \%reactions_to_products,
			 \%substrates_to_reactions,\%products_to_reactions,$create_assembly,$find_first);
}


sub create_reactions
{
    my ($substrates_to_reactions, $products_to_reactions,
	$reactions_to_substrates, $reactions_to_products) = @_;
    print STDERR "building SS reactions\n" if $debug;
    # use subsystem reactions
    foreach my $drxn (map { ($_."_L", $_."_R") } keys %all_reactions)
    {
	foreach my $substrArr ($reactions_to_substrate_arrays{$drxn})
	{
	    foreach my $cinfo (@$substrArr)
	    {
		my $cpd = $cinfo->[0];
		my $main = $cinfo->[2] || defined $all_inputs{$cpd}; # main in this reaction
		
		if ($main)
		{
		    push(@{$reactions_to_substrates->{$drxn}}, $cpd);
		    push(@{$substrates_to_reactions->{$cpd}}, $drxn);
		}
	    }
	}
	
	foreach my $prodArr ($reactions_to_product_arrays{$drxn})
	{
	    foreach my $cinfo (@$prodArr)
	    {
		my $cpd = $cinfo->[0];
		my $main = $cinfo->[2] || defined $all_outputs{$cpd}; # main in this reaction
		
		if ($main)
		{
		    push(@{$reactions_to_products->{$drxn}}, $cpd);
		    push(@{$products_to_reactions->{$cpd}}, $drxn);
		}
	    }
	}
    }
}

sub create_assembly_reactions
{
    my ($substrates_to_reactions, $products_to_reactions,
	$reactions_to_substrates, $reactions_to_products,$assembly_paths) = @_;
    my %intersection;
	
    foreach my $path (@$assembly_paths)
    {
	my $genome = shift @$path;
	my $paths_dir = get_scenario_directory($genome) . "/" . join "/" , @$path;
	
	my $drxn = "$path->[-3]/$path->[-2]/$path->[-1]_R";
	
	$all_reactions{"$path->[-3]/$path->[-2]/$path->[-1]"} = "R";
	
	print STDERR "Making reaction: $drxn\n" if $debug;
	
	open (M_INPUTS , $paths_dir."/inputs") or die ("Failed to open $paths_dir"."/inputs");
	my @substrArr;
	    
	print STDERR "Gathering Inputs:\n" if $debug;
	
	while (<M_INPUTS>)
	{
	    my ($cpd, $stoich) = split "\t" , $_;
	
	    #We are going to assume everything is a main...
	    if(!defined $all_compounds_to_main{$cpd})
	    {
		$all_compounds_to_main{$cpd} = 1;
	    }
    
	    if ($all_compounds_to_main{$cpd})
	    {
		push(@{$reactions_to_substrates->{$drxn}}, $cpd);
		push(@{$substrates_to_reactions->{$cpd}}, $drxn);
	    }
	    
	    if($all_compounds_to_main{$cpd} !=0 || !defined $all_compounds_to_main{$cpd})
	    {
		$all_inputs{$cpd} = 1;
	    }
	    push @substrArr, [$cpd, $stoich, $all_compounds_to_main{$cpd}];
	    
	    my @names = $fig->names_of_compound($cpd);
	    print STDERR "\t$stoich\t$cpd\t$names[0]\t$all_compounds_to_main{$cpd}\n" if $debug;
	}
	
	$reactions_to_substrate_arrays{$drxn} = \@substrArr;
	
	close M_INPUTS;
	
	open (M_OUTPUTS, $paths_dir."/outputs") or die("Failed to open $paths_dir"."/outputs");
	my @prodArr;
	
	print STDERR "Gathering outputs:\n" if $debug;
	
	while (<M_OUTPUTS>)
	{
	    my ($cpd, $stoich) = split "\t", $_;
	    print STDERR "Found $stoich $cpd\n" if $debug;
	    
	    #We are going to assume everything is a main...
	    if(!defined $all_compounds_to_main{$cpd})
	    {
		$all_compounds_to_main{$cpd} = 1;
	    }
	    
	    if ($all_compounds_to_main{$cpd})
	    {
		push(@{$reactions_to_products->{$drxn}}, $cpd);
		push(@{$products_to_reactions->{$cpd}}, $drxn);
	    }
	    
	    
	    #This adds cycles from 'assemblys' because they weren't added earlier
	    foreach my $ele (@{$reactions_to_substrates->{$drxn}})
	    {
		if($ele eq $cpd)
		{
		    $scenario_cycles{$cpd} = 1;
		}
	    }
	    
	    if($all_compounds_to_main{$cpd} !=0 || !defined $all_compounds_to_main{$cpd})
	    {
		$all_outputs{$cpd} = 1;
	    }
	    push @prodArr, [$cpd,$stoich,$all_compounds_to_main{$cpd}];
	    my @names = $fig->names_of_compound($cpd);
	    print STDERR "\t$cpd\t$names[0]\t$all_compounds_to_main{$cpd}\n" if $debug;
	    }
	
	$reactions_to_product_arrays{$drxn} = \@prodArr;
	
	close M_OUTPUTS;
	
    }
}


sub load_substrate_and_product_arrays
{
    my ($reactions, $map_ids, $sc_inputs) = @_;

    # determine whether the reaction is in one of the maps, and get directionality accordingly
    my (%reactions_in_maps, %reactions_not_in_maps, %reactions_not_in_any_map);

    foreach my $rxn (keys %$reactions)
    {
	my $direction;
	
	if($fig->valid_reaction_id($rxn))
	{
	    # get an array of triplets. The triplets are [reaction id][map id]
	    # [left to right - R, right to left - L, or both - B]
	    my @triplets = $fig->reaction_direction($rxn);

	    foreach my $trip (@triplets)
	    {
		foreach my $map_id (@$map_ids)
		{
		    if (@{$trip}[1] eq $map_id)
		    {
			my $this_direction = @{$trip}[2];

			# bidirectional in one map overrules unidirectional in another.
			# opposite directions in two maps becomes bidirectional
			if (! defined $direction)
			{
			    $direction = $this_direction;
			}
			elsif ($direction ne "B" && ($this_direction eq "B" || 
						     $this_direction ne $direction))
			{
			    $direction = "B";
			}

			$reactions_in_maps{$rxn} = 1;
		    }
		}
	    }
	    
	    
	    if(! $reactions_in_maps{$rxn})
	    {
		my $found_in_other_map = 0;
		
		#reaction not in scenario map ids, try to get directionality from other maps
		foreach my $trip (@triplets)
		{
		    my $this_direction =  @{$trip}[2];
		    
		    #bidreactional in one map overrules unidirectional in another
		    #opposite directions in two maps becomes bidirectional
		    if(! defined $direction)
		    {
			$direction = $this_direction;
		    }
		    elsif($direction ne "B" && ($this_direction eq "B" || $this_direction ne $direction))
		    {
			$direction = "B";
		    }
		    
		    $found_in_other_map = 1;
		}
		if (!$found_in_other_map)
		{
		    # reaction not in any map, get directionality without reference to map
		    if($fig->reversible($rxn) eq "1")
		    {
			$direction = "B";
		    }
		    else
		    {
			$direction = "R";
		    }

		    $reactions_not_in_any_map{$rxn} = 1;
		}
	    }

	    if (! defined $all_reactions{$rxn} || $direction eq "B")
	    {
		$all_reactions{$rxn} = $direction;
	    }
	    elsif ($all_reactions{$rxn} ne $direction)
	    {
		$all_reactions{$rxn} = "B";
	    }

	    my (@substrArr, @prodArr);

	    if ($direction eq "L")
	    {
		@substrArr = $fig->reaction2comp($rxn, 1, $map_ids);
		@prodArr = $fig->reaction2comp($rxn, 0, $map_ids);
		$reactions_to_substrate_arrays{$rxn . "_L"} = \@substrArr;
		$reactions_to_product_arrays{$rxn . "_L"} = \@prodArr;

	    }
	    else
	    {
		if ($reactions_in_maps{$rxn})
		{
		    @substrArr = $fig->reaction2comp($rxn, 0, $map_ids);
		    @prodArr = $fig->reaction2comp($rxn, 1, $map_ids);
		}
		else
		{
		    @substrArr = $fig->reaction2comp($rxn, 0);
		    @prodArr = $fig->reaction2comp($rxn, 1);
		}
		
		$reactions_to_substrate_arrays{$rxn."_R"} = \@substrArr;
		$reactions_to_product_arrays{$rxn."_R"} = \@prodArr;

		if ($direction eq "B")
		{
		    $reactions_to_substrate_arrays{$rxn."_L"} = \@prodArr;
		    $reactions_to_product_arrays{$rxn."_L"} = \@substrArr;
		}
	    }

	    print STDERR "\nFor $rxn, found substrates:\n" unless !$debug;
	    map { print STDERR "\t$_->[0]\t$_->[1]\t$_->[2]\n" unless !$debug } @substrArr;
	    print STDERR "For $rxn, found products:\n" unless !$debug;
	    map { print STDERR "\t$_->[0]\t$_->[1]\t$_->[2]\n" unless !$debug } @prodArr;

	    # load "main" designation based on reactions that are in the maps
	    if ($reactions_in_maps{$rxn})
	    {
		foreach my $cinfo ((@substrArr, @prodArr))
		{
		    my $cpd = $cinfo->[0];
		    my $main = $cinfo->[2];

		    if (defined $all_compounds_to_main{$cpd})
		    {
			# only main if it's main in all reactions
			$all_compounds_to_main{$cpd} &= $main;
		    }
		    else
		    {
			$all_compounds_to_main{$cpd} = $main;
		    }
		}
	    }
	    else
	    {
		# save subs and prods for processing at end
		$reactions_not_in_maps{$rxn} = [ (@substrArr, @prodArr) ];
	    }
	}
    }

    # now load "main" designation based on reactions not in the maps - but don't overrule
    # what's already been loaded
    my %additional_compounds_to_main;

    foreach my $rxn (keys %reactions_not_in_maps)
    {
	foreach my $cinfo (@{$reactions_not_in_maps{$rxn}})
	{
	    my $cpd = $cinfo->[0];
	    my $main = $cinfo->[2];

	    if (defined $additional_compounds_to_main{$cpd})
	    {
		# main if it's main in any reactions not in map
		$additional_compounds_to_main{$cpd} |= $main;
	    }
	    else
	    {
		$additional_compounds_to_main{$cpd} = $main;
	    }
	}
    }

    foreach my $cpd (keys %additional_compounds_to_main)
    {
	if (! defined $all_compounds_to_main{$cpd})
	{
	    $all_compounds_to_main{$cpd} = $additional_compounds_to_main{$cpd};
	}
    }

    # the reactions that aren't in any map at all won't have any main compounds.
    # mark those compounds that are in all_compounds_to_main as main.
    foreach my $rxn (keys %reactions_not_in_any_map)
    {
	print STDERR "Checking reaction not in any map: $rxn\n" if $debug;
	foreach my $cpd_array ($reactions_to_substrate_arrays{$rxn."_L"},
			       $reactions_to_substrate_arrays{$rxn."_R"},
			       $reactions_to_product_arrays{$rxn."_L"},
			       $reactions_to_product_arrays{$rxn."_R"})
	{
	    if (defined $cpd_array)
	    {
		my $at_least_one_is_main = 0;

		foreach my $cinfo (@{$cpd_array})
		{
		    my $cpd = $cinfo->[0];
		    if ($all_compounds_to_main{$cpd} || exists $sc_inputs->{$cpd})
		    {
			print STDERR "\t Setting $cpd to main\n" if $debug;
			$cinfo->[2] = 1;
			$at_least_one_is_main = 1;
		    }
		}

		if ($at_least_one_is_main == 0)
		{
		    foreach my $cinfo (@{$cpd_array})
		    {
			my $cpd = $cinfo->[0];
			print STDERR "\t Setting $cpd to provisional main\n" if $debug;
			$cinfo->[2] = 2;
		    }
		}
	    }
	}
    }
}

sub process_paths
{
    my ($input_cpds, $output_cpds, $outputs_lists, $reactions_to_substrates, $reactions_to_products, $substrates_to_reactions, $products_to_reactions, $create_assembly, $find_first ) = @_;

    my (%path_inputs, %path_outputs);
    map { $path_inputs{$_} = 0 } keys %$input_cpds;
    map { $path_outputs{$_} = 0 } keys %$output_cpds;

    my %data_results = ("infinite" => 0);

    # %compounds_to_tokens maps from compound ids to tokens placed on those compounds
    # the tokens are organized in hashes mapping from token id to number of tokens with that id
    my %compounds_to_tokens;
    map { $compounds_to_tokens{$_} = {} } keys %all_compounds_to_main;

    # %tokens maps from token_ids to the token data structures
    my %tokens;
    my $token_id_counter = 1;

    # %compounds_borrowed_to_tokens maps from compounds to lists of token ids that borrowed
    # compound in order to run a reaction
    my %compounds_borrowed_to_tokens;

    print STDERR "\nIn process_paths, path_inputs are @{[ keys %path_inputs ]}, path_outputs are @{[ keys %path_outputs] }, scenario_cycles are @{[ keys %scenario_cycles ]} \n\n" unless !$debug;

    my $initial_pass = 1;
    my $done = 0;
    my $loop_counter = 1;
    my $infinite_loop_check = $create_assembly ? $loop_max_assembly : $loop_max;
    
    # we may get to the point where we need to add some more path inputs into the mix
    # to push stalled tokens
    my $add_path_inputs = 0;

    while(!$done)
    {		
	if ($initial_pass || $add_path_inputs)
	{
	    foreach my $cpd (keys %path_inputs)
	    {
		# place a token on each path input
		my $new_token_id = $token_id_counter++;
		my %new_token;
		$new_token{visited_reactions} = {};
		$new_token{visited_compounds} = { $cpd => 0 }; # 0 means supplied from "outside"
		$new_token{token_path_inputs} = { $cpd => 1 }; # 1 means one was supplied
		$new_token{initial_pass} = $initial_pass;
		$compounds_to_tokens{$cpd}->{$new_token_id}++;
		$tokens{$new_token_id} = \%new_token;
			
		print STDERR "\t\tCreated new token '$new_token_id' for path input $cpd\n" unless !$debug;
	    }

	    $initial_pass = 0;
	    $add_path_inputs = 0;
	}

	# Find the reactions that can run.  
	my %reactions_to_try;

	foreach my $cpd (keys %compounds_to_tokens)
	{
	    foreach my $token_id (keys %{$compounds_to_tokens{$cpd}})
	    {
		if ($compounds_to_tokens{$cpd}->{$token_id} > 0 &&
		    ! $tokens{$token_id}->{done})
		{
		    # this compound has tokens that aren't done
		    map { $reactions_to_try{$_} = 1 } @{$substrates_to_reactions->{$cpd}};
		    last;
		}
	    }
	}

	print STDERR "\n\tIn loop, trying reactions @{[keys %reactions_to_try]}\n\n" unless !$debug;

	# Map the reactions that can run to the tokens they can use
	my %reactions_to_tokens_available;
	# Keep track of the main substrates defined by each reaction
	my %reactions_to_main_substrates;

	# count up the total number of tokens needed for each compound to run
	# every reaction that is ready to go
	my %reactions_to_tokens_needed;

      rxn: foreach my $reaction (keys %reactions_to_try)
	{	
	    print STDERR "\tChecking reaction $reaction\n" unless !$debug;

	    my @substrArr = @{$reactions_to_substrate_arrays{$reaction}};
	    my @prodArr = @{$reactions_to_product_arrays{$reaction}};
	    my %main_substrates;

	    # Determine if this reaction has necessary inputs to run.
	    # There must be tokens available for at least one main substrate that isn't
	    # a path output (unless it's an initial scenario cycled compound).
	    # Also, any path input must have a token.
	    # Count number of tokens needed for each main substrate.
	    my %tokens_available;
	    my %tokens_needed;
	    my $reaction_can_run = 0;

	    foreach my $substr (@substrArr)
	    {
		my $cpd = @{$substr}[0];
		my $stoich = @{$substr}[1];
		my $is_it_main = @{$substr}[2];
		my $main = $is_it_main == 1 || $all_compounds_to_main{$cpd}; #main either way

		if ($is_it_main == 2)
		{
		    # provisionally main - need to see if there is a token that has visited
		    # this compound
		    
		    foreach my $token_id (keys %{$compounds_to_tokens{$cpd}})
		    {
			my %visited_compounds = %{$tokens{$token_id}->{visited_compounds}};
			$main = 1 if exists $visited_compounds{$cpd};
		    }
		}
		
		$main_substrates{$cpd} = 1 if $main;
		print STDERR "\t\tSubstrate: $cpd\tstoich: $stoich\tmain: $main\n" unless !$debug;

		if (! $main)
		{
		    # on any pass we can take in non-main compounds
		    $tokens_needed{$cpd} = $stoich;
		}
		else
		{
		    # if tokens are available for compound, check their history for the
		    # main compounds produced by the reaction so we don't loop back over
		    # previous main compounds.  Also, don't use tokens on scenario inputs
		    # moved to in steps other than initial token creation.
		    my %ok_tokens;
		    my $num_ok_tokens = 0;

		    foreach my $token_id (keys %{$compounds_to_tokens{$cpd}})
		    {
			next if $tokens{$token_id}->{done} || 
			    $compounds_to_tokens{$cpd}->{$token_id} == 0;

			# now check that we aren't running an already visited reaction in reverse
			my %visited_reactions = %{$tokens{$token_id}->{visited_reactions}};

			if (($reaction =~ /(.*)_R/ && defined $visited_reactions{$1."_L"}) ||
			    ($reaction =~ /(.*)_L/ && defined $visited_reactions{$1."_R"}))
			{
			    print STDERR "\t\tToken '$token_id' has run the reverse reaction already\n" unless !$debug;
			    next;
			}

			my %visited_compounds = %{$tokens{$token_id}->{visited_compounds}};

			print STDERR "\t\tToken '$token_id' has visited @{[map { ($_, $visited_compounds{$_} ) } sort { $visited_compounds{$a} <=> $visited_compounds{$b} } keys %visited_compounds]}\n" unless !$debug;

			# need to find at least one prod that hasn't been visited yet
			# or was visited in a loop cycle not before the loop cycle in which
			# it visited the substrate,
			# or is a path output
			my $prods_are_ok = 0;

			# check each main product
			foreach my $prod (@{$reactions_to_products->{$reaction}})
			{
			    if (! defined $visited_compounds{$prod} || 
				$visited_compounds{$prod} >= $visited_compounds{$cpd} ||
				defined $path_outputs{$prod} ||
				$compounds_borrowed_to_tokens{$prod}->{$token_id} > 0)
			    {
				print STDERR "\t\tToken can visit $prod\n" unless !$debug;
				$prods_are_ok = 1;
				last;
			    }
			}

			if ($prods_are_ok)
			{
			    print STDERR "\t\tToken is OK\n" unless !$debug;
			    $ok_tokens{$token_id} = $compounds_to_tokens{$cpd}->{$token_id};
			    $num_ok_tokens += $compounds_to_tokens{$cpd}->{$token_id};
			}
		    }

		    map { $tokens_available{$_}->{$cpd} = $ok_tokens{$_} } keys %ok_tokens;
		    $tokens_needed{$cpd} = $stoich;

		    if ($main && $num_ok_tokens >= 1)
		    {
			if (! defined $path_outputs{$cpd} || defined $scenario_cycles{$cpd})
			{
			    print STDERR "\t\tgot at least one token on main compound\n" unless !$debug;
			    $reaction_can_run = 1;
			}
		    }
		    elsif (defined $path_inputs{$cpd})
		    {
			print STDERR "\t\tno tokens available for path input: $cpd\n" unless !$debug;
			next rxn;
		    }
		    elsif (defined $path_outputs{$cpd} && ! defined $scenario_cycles{$cpd})
		    {
			print STDERR "\t\tno tokens available for path output: $cpd\n" unless !$debug;
			next rxn;
		    }
		}
	    }

	    if ($reaction_can_run)
	    {
		$reactions_to_tokens_available{$reaction} = \%tokens_available;
		$reactions_to_tokens_needed{$reaction} = \%tokens_needed;
		$reactions_to_main_substrates{$reaction} = \%main_substrates;
		print STDERR "\tReaction $reaction can run\n" unless !$debug;
	    }
	}

	# keep track of tokens used that will be used to run rxns.  Clone tokens if necessary.
	my %reactions_to_tokens_to_use;
	my %tokens_to_use_to_reactions;
	my %copy_of_compounds_to_tokens; # for determining which reaction uses which tokens

	foreach my $cpd (keys %compounds_to_tokens)
	{
	    my %cpd_token_ids = %{$compounds_to_tokens{$cpd}};
	    my %new_cpd_token_ids;
	    map { $new_cpd_token_ids{$_} = $cpd_token_ids{$_} } keys %cpd_token_ids;
	    $copy_of_compounds_to_tokens{$cpd} = \%new_cpd_token_ids;
	}

	foreach my $reaction (keys %reactions_to_tokens_available)
	{
	    print STDERR "\tPreparing to run reaction $reaction\n" unless !$debug;

	    # assemble tokens to run reaction, cloning ones that were used by
	    # other reactions during this cycle if necessary
	    my %tokens_available = %{$reactions_to_tokens_available{$reaction}};
	    my %tokens_needed = %{$reactions_to_tokens_needed{$reaction}};
	    my @final_tokens_to_use; # list of maps from substrates to tokens to use
	    my %clone_history; # map from token ids to ids of their new clones

	    # check to see if any available tokens are already commited to the reverse reaction
	    my @token_ids = keys %tokens_available;

	    foreach my $token_id (keys %tokens_available)
	    {
		if (($reaction =~ /(.*)_R/ && defined $tokens_to_use_to_reactions{$token_id}->{$1."_L"}) ||
		    ($reaction =~ /(.*)_L/ && defined $tokens_to_use_to_reactions{$token_id}->{$1."_R"}))
		{
		    # clone the token
		    my $new_token_id = $token_id_counter++;
		    &clone_token($token_id, $new_token_id, \%tokens, \%compounds_to_tokens, 
				 \%compounds_borrowed_to_tokens);
		    
		    foreach my $icpd (keys %compounds_to_tokens)
		    {
			if ($compounds_to_tokens{$icpd}->{$token_id} > 0)
			{
			    $copy_of_compounds_to_tokens{$icpd}->{$new_token_id} = 
				$compounds_to_tokens{$icpd}->{$token_id};
			}
		    }

		    print STDERR "\t\tCloned token '$token_id', new token is '$new_token_id'\n" unless !$debug;
		    
		    $tokens_available{$new_token_id} = $tokens_available{$token_id};
		    delete $tokens_available{$token_id};
		}
	    }

	    # first, assemble tokens that have all the main compounds they need to run
	    @token_ids = keys %tokens_available;

	    foreach my $token_id (keys %tokens_available)
	    {
		my $has_all_main_cpds = 1;
		my $need_to_clone = 0;

		foreach my $cpd (keys %tokens_needed)
		{
		    if ($reactions_to_main_substrates{$reaction}->{$cpd})
		    {
			if ($compounds_to_tokens{$cpd}->{$token_id} >= $tokens_needed{$cpd})
			{
			    if ($copy_of_compounds_to_tokens{$cpd}->{$token_id} < $tokens_needed{$cpd})
			    {
				$need_to_clone = 1;
			    }
			}
			else
			{
			    $has_all_main_cpds = 0;
			    last;
			}
		    }		    
		}

		if ($has_all_main_cpds)
		{
		    print STDERR "\t\ttoken '$token_id' has all main compounds\n" unless !$debug;

		    delete $tokens_available{$token_id};

		    if ($need_to_clone)
		    {
			my $new_token_id = $token_id_counter++;
			&clone_token($token_id, $new_token_id, \%tokens, 
				     \%compounds_to_tokens, 
				     \%compounds_borrowed_to_tokens);
			
			foreach my $icpd (keys %compounds_to_tokens)
			{
			    if ($compounds_to_tokens{$icpd}->{$token_id} > 0)
			    {
				$copy_of_compounds_to_tokens{$icpd}->{$new_token_id} = 
				    $compounds_to_tokens{$icpd}->{$token_id};
			    }
			}

			print STDERR "\t\tCloned token '$token_id', new token is '$new_token_id'\n" unless !$debug;
			
			$clone_history{$token_id} = $new_token_id;
			$token_id = $new_token_id;
		    }
		    
		    # now assemble the compound to tokens map
		    my %cpd_to_tokens;

		    foreach my $cpd (keys %tokens_needed)
		    {
			if ($reactions_to_main_substrates{$reaction}->{$cpd})
			{
			    for (my $i = 0; $i < $tokens_needed{$cpd}; $i++)
			    {
				push @{$cpd_to_tokens{$cpd}}, $token_id;
				$copy_of_compounds_to_tokens{$cpd}->{$token_id}--;
			    }
			}
		    }

		    push @final_tokens_to_use, \%cpd_to_tokens;
		}
	    }

	    if (scalar keys %tokens_available > 0)
	    {
		# try to merge left over available tokens in all combinations that fulfill 
		# needed main substrates.  Create map from main substrates to tokens.
		my %tokens_available_for_cpds;
		my %tokens_to_use_for_cpds;

		foreach my $token_id (keys %tokens_available)
		{
		    foreach my $cpd (keys %{$tokens_available{$token_id}})
		    {
			if ($reactions_to_main_substrates{$reaction}->{$cpd})
			{
			    $tokens_available_for_cpds{$cpd}->{$token_id} =
				$tokens_available{$token_id}->{$cpd};
			}
		    }
		}

		foreach my $cpd (keys %{$reactions_to_main_substrates{$reaction}})
		{
		    my %available_token_ids;

		    foreach my $token_id (keys %{$tokens_available_for_cpds{$cpd}})
		    {
			my $updated_token_id = $token_id;

			# in case the token has already been cloned for another cpd in this reaction
			while ($clone_history{$updated_token_id})
			{
			    $updated_token_id = $clone_history{$updated_token_id};
			}

			$available_token_ids{$updated_token_id} = 1;

			if ($token_id != $updated_token_id)
			{
			    $tokens_available_for_cpds{$cpd}->{$updated_token_id} = 
				$tokens_available_for_cpds{$cpd}->{$token_id};
			    delete $tokens_available_for_cpds{$cpd}->{$token_id};
			}
		    }

		    my $num_tokens_needed = $tokens_needed{$cpd};
		    my $num_available_tokens = 0;

		    foreach my $token_id (keys %available_token_ids)
		    {
			$num_available_tokens += $tokens_available_for_cpds{$cpd}->{$token_id};
		    }

		    # if not enough tokens are available to fill out sets, create new ones
		    my $num_short_of_full;

		    if ($num_available_tokens == 0)
		    {
			$num_short_of_full = $num_tokens_needed;
		    }
		    else
		    {
			$num_short_of_full = ($num_tokens_needed - ($num_available_tokens % $num_tokens_needed)) % $num_tokens_needed;
		    }

		    print STDERR "\t\tFor $cpd, $num_tokens_needed tokens are needed, $num_available_tokens are available. Need to create $num_short_of_full to fill out sets\n" unless !$debug;

		    for (my $i = 0; $i < $num_short_of_full; $i++) 
		    {
			my $new_token_id = $token_id_counter++;
			my %new_token;
			print STDERR "\t\tCreating new token '$new_token_id' for $cpd\n" unless !$debug;

			$new_token{visited_reactions} = {};
			$new_token{visited_compounds} = {}; 
			$compounds_to_tokens{$cpd}->{$new_token_id}++;
			$copy_of_compounds_to_tokens{$cpd}->{$new_token_id}++;
			$tokens{$new_token_id} = \%new_token;
			$available_token_ids{$new_token_id} = 1;
			$tokens_available_for_cpds{$cpd}->{$new_token_id}++;

			# if it's not a path input, remember that we've "borrowed" it and
			# will need to pay it back.
			if (! defined $path_inputs{$cpd})
			{
			    $new_token{visited_compounds}->{$cpd} = $loop_counter;
			    $compounds_borrowed_to_tokens{$cpd}->{$new_token_id}++;
			}
			else
			{
			    $new_token{token_path_inputs} = { $cpd => 1 }; 
			    $new_token{visited_compounds}->{$cpd} = 0;
			}
		    }

		    # for main compounds, there may be more tokens available than needed,
		    # so we may assemble multiple token sets.
		    my %tokens_not_yet_used;

		    foreach my $token_id (keys %available_token_ids)
		    {
			if ($copy_of_compounds_to_tokens{$cpd}->{$token_id} > 0)
			{
			    $tokens_not_yet_used{$token_id} = $tokens_available_for_cpds{$cpd}->{$token_id};
			}
		    }

		    my @token_sets_for_cpd;
		    
		    print STDERR "\t\tNeed $num_tokens_needed tokens for $cpd, '@{[ sort { $a <=> $b } keys %available_token_ids ]}' are usable, '@{[ sort { $a <=> $b } keys %tokens_not_yet_used ]}' are not yet used\n" unless !$debug;

		    my @token_set;
		    my $num_tokens_still_needed = $num_tokens_needed;

		    foreach my $token_id (sort { $a <=> $b } keys %available_token_ids)
		    {
			# need to clone the token if it is all used up
			if ($copy_of_compounds_to_tokens{$cpd}->{$token_id} == 0)
			{
			    my $new_token_id = $token_id_counter++;
			    &clone_token($token_id, $new_token_id, \%tokens, 
					 \%compounds_to_tokens, 
					 \%compounds_borrowed_to_tokens);
			    
			    foreach my $icpd (keys %compounds_to_tokens)
			    {
				if ($compounds_to_tokens{$icpd}->{$token_id} > 0)
				{
				    $copy_of_compounds_to_tokens{$icpd}->{$new_token_id} = 
					$compounds_to_tokens{$icpd}->{$token_id};
				}
			    }

			    print STDERR "\t\tCloned token '$token_id' for $cpd, new token is '$new_token_id'\n" unless !$debug;
			    
			    $clone_history{$token_id} = $new_token_id;
			    $token_id = $new_token_id;
			}

			while ($copy_of_compounds_to_tokens{$cpd}->{$token_id} > 0)
			{
			    push @token_set, $token_id;
			    $copy_of_compounds_to_tokens{$cpd}->{$token_id}--;
			    $num_tokens_still_needed--;

			    if ($num_tokens_still_needed == 0)
			    {
				print STDERR "\t\tPushing token set '@token_set' for $cpd\n" unless !$debug;
				my @copy_of_token_set = @token_set;
				push @token_sets_for_cpd, \@copy_of_token_set;
				@token_set = ();
				$num_tokens_still_needed = $num_tokens_needed;
			    }
			}
		    }

		    $tokens_to_use_for_cpds{$cpd} = \@token_sets_for_cpd;
		}

		# in case a token had to be cloned for this reaction after another
		# compound already determined to use it, check history and use new token id
		foreach my $cpd (keys %tokens_to_use_for_cpds)
		{
		    foreach my $token_set (@{$tokens_to_use_for_cpds{$cpd}})
		    {
			for (my $i = 0; $i < scalar @$token_set; $i++)
			{
			    my $updated_token_id = $token_set->[$i];

			    while ($clone_history{$updated_token_id})
			    {
				$updated_token_id = $clone_history{$updated_token_id};
			    }

			    if ($token_set->[$i] != $updated_token_id)
			    {
				print STDERR "\t\tReplacing '$token_set->[$i]' with '$updated_token_id' for $cpd\n" unless !$debug;
				splice @$token_set, $i, 1, ($updated_token_id);
			    }
			}
		    }
		}

		# we may have multiple sets of tokens to use for some compounds, so
		# make sure we're prepared for each combination by cloning sets as necessary
		my %cpd_to_token_set_index;
		my $num_combinations = 1;
		my %used_token_sets_for_cpd;

		foreach my $cpd (keys %tokens_to_use_for_cpds)
		{
		    $cpd_to_token_set_index{$cpd} = 0;
		    $num_combinations *= scalar @{$tokens_to_use_for_cpds{$cpd}};
		}
		
		for (my $i = 0; $i < $num_combinations; $i++)
		{
		    my %combination_cpds_to_tokens;

		    print STDERR "\t\tPreparing combination ", $i+1, " out of $num_combinations, indices are @{[ map { $cpd_to_token_set_index{$_} } sort keys %cpd_to_token_set_index ]}\n" unless !$debug;

		    # keep track of who was cloned to what for this combination in case
		    # different compounds are using the same tokens
		    my %clone_history_this_combination;

		    foreach my $cpd (sort keys %cpd_to_token_set_index)
		    {
			my $token_set_index = $cpd_to_token_set_index{$cpd};
			my @token_set = @{$tokens_to_use_for_cpds{$cpd}->[$token_set_index]};
			my @new_token_set;

			# don't clone if this is the first time using this token set
			# should be able to do this mathematically
			if (! defined $used_token_sets_for_cpd{$cpd}->{$token_set_index})
			{
			    @new_token_set = @token_set;
			    $used_token_sets_for_cpd{$cpd}->{$token_set_index} = 1;
			}
			else
			{		    
			    foreach my $token_id (@token_set)
			    {
				my $new_token_id;

				if ($clone_history_this_combination{$token_id})
				{
				    $new_token_id = $clone_history_this_combination{$token_id};
				}
				else
				{
				    $new_token_id = $token_id_counter++;
				    &clone_token($token_id, $new_token_id, \%tokens, 
						 \%compounds_to_tokens, 
						 \%compounds_borrowed_to_tokens);
				    
				    print STDERR "\t\t\tCloned token '$token_id' for $cpd, new token is '$new_token_id'\n" unless !$debug;
				    
				    $clone_history_this_combination{$token_id} = $new_token_id;
				}
				
				push @new_token_set, $new_token_id;
			    }
			}
			
			push @{$combination_cpds_to_tokens{$cpd}}, @new_token_set;
		    }

		    push @final_tokens_to_use, \%combination_cpds_to_tokens;

		    # move the cpd to token set indices for the next combination
		    foreach my $cpd (sort keys %cpd_to_token_set_index)
		    {
			if ($cpd_to_token_set_index{$cpd} < scalar @{$tokens_to_use_for_cpds{$cpd}} - 1)
			{
			    $cpd_to_token_set_index{$cpd}++;
			    last;
			}
			else
			{
			    $cpd_to_token_set_index{$cpd} = 0;
			}
		    }
		}
	    }

	    # last step, need to create place-holder tokens for the non-main substrates
	    # and record the tokens being used for this reaction
	    foreach my $combination (@final_tokens_to_use)
	    {
		foreach my $cpd (keys %tokens_needed)
		{
		    if (! defined $reactions_to_main_substrates{$reaction}->{$cpd})
		    {
			my $new_token_id = $token_id_counter++;
			my %new_token;
			$new_token{visited_reactions} = {};
			$new_token{visited_compounds} = {}; 
			$compounds_to_tokens{$cpd}->{$new_token_id} += $tokens_needed{$cpd};
			$copy_of_compounds_to_tokens{$cpd}->{$new_token_id} += $tokens_needed{$cpd};
			$new_token{token_path_inputs} = { $cpd => $tokens_needed{$cpd} }; 
			$tokens{$new_token_id} = \%new_token;

			for (my $i = 0; $i < $tokens_needed{$cpd}; $i++)
			{
			    push @{$combination->{$cpd}}, $new_token_id;
			}
		    }
		}

		foreach my $cpd (keys %$combination)
		{
		    map { $tokens_to_use_to_reactions{$_}->{$reaction} = 1 } @{$combination->{$cpd}};
		}
	    }

	    $reactions_to_tokens_to_use{$reaction} = \@final_tokens_to_use;
	}
	
	# keep track of tokens merged during this round
	my %token_merge_history;

	foreach my $reaction (keys %reactions_to_tokens_to_use)
	{
	    # we may have multiple sets of tokens to use 
	    my @final_tokens_to_use = @{$reactions_to_tokens_to_use{$reaction}};
	    my $num_combinations = scalar @final_tokens_to_use;

	    # since we may be running this reaction several times with some of the same
	    # tokens, don't process token merge history until all sets have been run
	    my %token_merge_history_this_reaction;
	    
	    for (my $i = 0; $i < $num_combinations; $i++)
	    {
		print STDERR "\tRunning reaction $reaction (", $i+1, " out of $num_combinations)\n" unless !$debug;

		# assemble list of tokens to use for this combination
		my %tokens_to_use;
		# find the most recently visited compound's step
		my $most_recent_step = 0;

		# loop through substrates
		foreach my $cpd (keys %{$final_tokens_to_use[$i]})
		{
		    # remove the token ids from the real compounds_to_tokens map
		    my @cpd_token_ids = @{$final_tokens_to_use[$i]->{$cpd}};
		    print STDERR "\t\tFound tokens '@cpd_token_ids' for $cpd\n" unless !$debug;

		    foreach my $token_id (@cpd_token_ids)
		    {
			# check if token_id has been merged into a new token by a previous reaction
			while (defined $token_merge_history{$token_id})
			{
			    $token_id = $token_merge_history{$token_id};
			}

			$compounds_to_tokens{$cpd}->{$token_id}--;

			if ($compounds_to_tokens{$cpd}->{$token_id} == 0)
			{
			    delete $compounds_to_tokens{$cpd}->{$token_id};
			}

			$tokens_to_use{$token_id} = 1;

			if ($tokens{$token_id}->{visited_compounds}->{$cpd} > $most_recent_step)
			{
			    $most_recent_step = $tokens{$token_id}->{visited_compounds}->{$cpd};
			}			
		    }
		}

		# process list of unique token ids
		my @tokens_to_use = sort { $a <=> $b } keys %tokens_to_use;
		my $go_forward_token_id;

		if (scalar @tokens_to_use == 1)
		{
		    $go_forward_token_id = shift @tokens_to_use;
		    print STDERR "\t\tGoing forward with '$go_forward_token_id'\n" unless !$debug;
		}
		else
		{
		    $go_forward_token_id = $token_id_counter++;
		    $tokens{$go_forward_token_id} = {};
		    print STDERR "\t\tRemember to merge tokens '@tokens_to_use' into '$go_forward_token_id'\n" unless !$debug;

		    # record the need to merge - we'll do it after processing all sets of substrates
		    foreach my $token_id (@tokens_to_use)
		    {
			push @{$token_merge_history_this_reaction{$token_id}}, $go_forward_token_id;
		    }
		}

		my $go_forward_token = $tokens{$go_forward_token_id};
		my @prodArr = @{$reactions_to_product_arrays{$reaction}};

		# add current reaction and products to accumulated token history.
		# reaction is mapped to loop counter to maintain history of order of execution
		$go_forward_token->{visited_reactions}->{$reaction} = $loop_counter;

		foreach my $prod (@prodArr)
		{
		    my $cpd = @{$prod}[0];
		    my $stoich = @{$prod}[1];
		    my $main = @{$prod}[2];

		    # keep track of path outputs we've seen
		    if (defined $path_outputs{$cpd})
		    {
			$path_outputs{$cpd} += $stoich;
		    }

		    if ($main)
		    {
			$go_forward_token->{visited_compounds}->{$cpd} = $most_recent_step + 1;
		    }

		    # push tokens
		    $compounds_to_tokens{$cpd}->{$go_forward_token_id} += $stoich;
		}
	    }

	    # now process this reaction's token merge history
	    foreach my $token_id (sort { $a <=> $b } keys %token_merge_history_this_reaction)
	    {
		# find the unique set of up to date merge ids
		my @merge_ids = @{$token_merge_history_this_reaction{$token_id}};
		my %updated_merge_ids;

		foreach my $token_id (@merge_ids)
		{
		    while (defined $token_merge_history{$token_id})
		    {
			$token_id = $token_merge_history{$token_id};
		    }

		    $updated_merge_ids{$token_id} = 1;
		}

		print STDERR "\t\tupdated merge id list for '$token_id': '@{[ keys %updated_merge_ids ]}'\n" unless !$debug;

		my $wrap_up_token_id = $token_id_counter++;
		$tokens{$wrap_up_token_id} = {};
		my $wrap_up_token = $tokens{$wrap_up_token_id};
		my @tokens_to_merge;
		push @tokens_to_merge, keys %updated_merge_ids, $token_id;

		# merge from oldest to youngest to update visited_compounds history
		foreach my $itoken_id (sort { $a <=> $b } @tokens_to_merge)
		{
		    print STDERR "\t\t\tmerging '$itoken_id' into '$wrap_up_token_id'\n" unless !$debug;

		    my $itoken = $tokens{$itoken_id};
		    map { $wrap_up_token->{visited_reactions}->{$_} = $itoken->{visited_reactions}->{$_} } keys %{$itoken->{visited_reactions}};
		    map { $wrap_up_token->{visited_compounds}->{$_} = $itoken->{visited_compounds}->{$_} } keys %{$itoken->{visited_compounds}};
		    map { $wrap_up_token->{token_path_inputs}->{$_} += $itoken->{token_path_inputs}->{$_} } keys %{$itoken->{token_path_inputs}};
		    $wrap_up_token->{initial_pass} |= $itoken->{initial_pass};
		    $token_merge_history{$itoken_id} = $wrap_up_token_id;

		    # tokens might be spread across multiple compounds; change them all
		    # to new id
		    foreach my $cpd (keys %compounds_to_tokens)
		    {
			if ($compounds_to_tokens{$cpd}->{$itoken_id} > 0)
			{
			    $compounds_to_tokens{$cpd}->{$wrap_up_token_id} +=
				$compounds_to_tokens{$cpd}->{$itoken_id};
			}
		    }

		    # tokens might have borrowed compounds; change them all to new id
		    foreach my $cpd (keys %compounds_borrowed_to_tokens)
		    {
			if ($compounds_borrowed_to_tokens{$cpd}->{$itoken_id} > 0)
			{
			    $compounds_borrowed_to_tokens{$cpd}->{$wrap_up_token_id} +=
				$compounds_borrowed_to_tokens{$cpd}->{$itoken_id};
			}
		    }
		}
	    }		
	}

	# now delete the tokens that were used merged in these reactions
	foreach my $token_id (keys %token_merge_history)
	{
	    foreach my $icpd (keys %compounds_to_tokens)
	    {
		delete $compounds_to_tokens{$icpd}->{$token_id};
	    }
		
	    foreach my $icpd (keys %compounds_borrowed_to_tokens)
	    {
		delete $compounds_borrowed_to_tokens{$icpd}->{$token_id};
	    }
		
	    print STDERR "\t\tDeleting token '$token_id'\n" unless !$debug;
	    delete $tokens{$token_id};
	}

	print STDERR "\nBalancing tokens\n" unless !$debug;

	foreach my $token_id (keys %tokens)
	{
	    &balance_borrowing_and_giving($token_id, \%compounds_to_tokens,
					  \%compounds_borrowed_to_tokens);
	}
	
	&print_token_status([sort { $a <=> $b } keys %tokens], \%tokens, \%compounds_to_tokens, \%compounds_borrowed_to_tokens, $fig);

	print STDERR "\nChecking for done\n" unless !$debug;

	print STDERR "\n\ntoken ids: @{[ sort { $a <=> $b } map { $_ if ! defined $tokens{$_}->{done} } keys %tokens ]}\n" unless !$debug;

	# we're done when all the main compounds in initial-pass tokens
	# have reached path outputs and repaid their borrowed tokens,
	# or have reached a dead end.  
	# Check if we're done pushing and borrowing tokens first.
	my %not_done_tokens;

	foreach my $token_id (keys %tokens)
	{
	    if (! defined $tokens{$token_id}->{done})
	    {
		if (! &check_token_for_done($token_id, \%compounds_to_tokens, 
					    \%compounds_borrowed_to_tokens,
					    \%all_compounds_to_main, \%path_outputs, 
					    \%scenario_cycles, \%tokens, $fig, $outputs_lists))
		{
		    $not_done_tokens{$token_id} = 1;
		}
	    }
	}

	print STDERR "\nChecking if we can pay back borrowed compounds from other tokens\n" unless !$debug;

	foreach my $bcpd (keys %compounds_borrowed_to_tokens)
	{
	    next if $path_outputs{$bcpd}; # tokens must manage their own path outputs

	    if (scalar keys %{$compounds_borrowed_to_tokens{$bcpd}} > 0 &&
		scalar keys %{$compounds_to_tokens{$bcpd}} > 0)
	    {
		my %borrowers_to_givers;

		foreach my $borrower_id (keys %{$compounds_borrowed_to_tokens{$bcpd}})
		{
		    next if defined $tokens{$borrower_id}->{done}; # don't repay deadenders

		    my $num_needed = $compounds_borrowed_to_tokens{$bcpd}->{$borrower_id};

		  giver: foreach my $giver_id (keys %{$compounds_to_tokens{$bcpd}})
		  {
		      next if defined $tokens{$giver_id}->{done};
		    
		      my $num_to_give = $compounds_to_tokens{$bcpd}->{$giver_id};

		      print STDERR "\tToken '$giver_id' has $num_to_give $bcpd to give to '$borrower_id', which needs $num_needed\n" unless !$debug;

		      my $borrower = $tokens{$borrower_id};
		      my $giver = $tokens{$giver_id};

		      # check whether the giver and borrower have conflicting histories
		      foreach my $visited_reaction (keys %{$giver->{visited_reactions}})
		      {
			  if (($visited_reaction =~ /(.*)_R/ && 
			       defined $borrower->{visited_reactions}->{$1."_L"}) ||
			      ($visited_reaction =~ /(.*)_L/ && 
			       defined $borrower->{visited_reactions}->{$1."_R"}))
			  {
			      print STDERR "\t\tConflict on $visited_reaction\n" unless !$debug;
			      next giver;
			  }
			  
		      }

		      push @{$borrowers_to_givers{$borrower_id}}, $giver_id;
		  }
		}

		# we have a list of potential givers for each borrower for this compound.  
		# Now figure out who the lucky givers will be.
		my %givers_to_borrowers;

		foreach my $borrower_id (keys %borrowers_to_givers)
		{
		    my %lucky_givers;
		    my $num_needed = $compounds_borrowed_to_tokens{$bcpd}->{$borrower_id};

		    print STDERR "\tCollecting lucky givers for '$borrower_id' for $bcpd, need $num_needed\n" unless !$debug;

		    # check potential givers starting with those with the most
		    my @potential_givers = reverse sort { $compounds_to_tokens{$bcpd}->{$a} <=> $compounds_to_tokens{$bcpd}->{$b} } @{$borrowers_to_givers{$borrower_id}};
			    
		    foreach my $giver_id (@potential_givers)
		    {
			if (! defined $lucky_givers{$giver_id})
			{
			    my $num_to_give = $compounds_to_tokens{$bcpd}->{$giver_id};
			    print STDERR "\t\t'$giver_id' has $num_to_give to give\n" unless !$debug;
			    $lucky_givers{$giver_id} = 1;
			    $num_needed -= $num_to_give;
			    last if $num_needed <= 0;
			}
		    }

		    foreach my $giver_id (keys %lucky_givers)
		    {
			push @{$givers_to_borrowers{$giver_id}}, $borrower_id;
		    }
		}


		foreach my $orig_giver_id (keys %givers_to_borrowers)
		{
		    my @borrowers_list = @{$givers_to_borrowers{$orig_giver_id}};
		    my @givers_list = ($orig_giver_id);

		    # clone enough givers so that every borrower gets one.  Last borrower
		    # gets the orgiinal giver.

		    while (scalar @borrowers_list > scalar @givers_list)
		    {
			my $new_giver_id = $token_id_counter++;
			&clone_token($orig_giver_id, $new_giver_id, \%tokens, 
				     \%compounds_to_tokens, \%compounds_borrowed_to_tokens);
			push @givers_list, $new_giver_id;
		    }
		    
		    for (my $k = 0; $k < scalar @givers_list; $k++)
		    {
			my $giver_id = $givers_list[$k];
			my $giver = $tokens{$giver_id};
			my $borrower_id = $borrowers_list[$k];
			my $borrower = $tokens{$borrower_id};

			print STDERR "\n\tMerging '$giver_id' into '$borrower_id'\n" unless !$debug;

			# bump the borrower's visited reactions and compounds counters forward,
			# then merge the giver's visited reactions and compounds, unless
			# the borrower has already visited them
			map { $borrower->{visited_reactions}->{$_} += $loop_counter } keys %{$borrower->{visited_reactions}};
			map { $borrower->{visited_compounds}->{$_} += $loop_counter } keys %{$borrower->{visited_compounds}};
			map { $borrower->{visited_reactions}->{$_} = $giver->{visited_reactions}->{$_} unless defined $borrower->{visited_reactions}->{$_} } keys %{$giver->{visited_reactions}};
			map { $borrower->{visited_compounds}->{$_} = $giver->{visited_compounds}->{$_} unless defined $borrower->{visited_compounds}->{$_} } keys %{$giver->{visited_compounds}};
			map { $borrower->{token_path_inputs}->{$_} += $giver->{token_path_inputs}->{$_} } keys %{$giver->{token_path_inputs}};
			$borrower->{initial_pass} |= $giver->{initial_pass};

			foreach my $icpd (keys %compounds_to_tokens)
			{
			    if ($compounds_to_tokens{$icpd}->{$giver_id} > 0)
			    {
				$compounds_to_tokens{$icpd}->{$borrower_id} += $compounds_to_tokens{$icpd}->{$giver_id};
				delete $compounds_to_tokens{$icpd}->{$giver_id};
			    }
			}
			
			foreach my $icpd (keys %compounds_borrowed_to_tokens)
			{
			    if ($compounds_borrowed_to_tokens{$icpd}->{$giver_id} > 0)
			    {
				$compounds_borrowed_to_tokens{$icpd}->{$borrower_id} += $compounds_borrowed_to_tokens{$icpd}->{$giver_id};
				delete $compounds_borrowed_to_tokens{$icpd}->{$giver_id};
			    }
			}

			delete $tokens{$giver_id};
			delete $not_done_tokens{$giver_id};

			&balance_borrowing_and_giving($borrower_id, \%compounds_to_tokens,
						      \%compounds_borrowed_to_tokens);
			
			&print_token_status([$borrower_id], \%tokens, \%compounds_to_tokens,
					    \%compounds_borrowed_to_tokens, $fig);

			if (&check_token_for_done($borrower_id, \%compounds_to_tokens, 
						  \%compounds_borrowed_to_tokens,
						  \%all_compounds_to_main, \%path_outputs, 
						  \%scenario_cycles,\%tokens, $fig, $outputs_lists))
			{
			    delete $not_done_tokens{$borrower_id};
			}
		    }
		}
	    }
	}

	# Now check if we've reached a dead end, either a compound we can't push or a
	# borrowed compound we can't repay.  Also determine whether there is a reaction
	# to run that can move an initial pass token forward.

	print STDERR "\nChecking for dead ends\n" unless !$debug;

	my $found_reaction_for_initial_pass_token = 0;

      check: foreach my $token_id (keys %not_done_tokens)
      {
	  # determine which compounds the token is sitting on, and whether a
	  # reaction can proceed from those compounds that isn't a loop
	  # back to compounds already visited
	  my %visited_compounds = %{$tokens{$token_id}->{visited_compounds}};
	  my %visited_reactions = %{$tokens{$token_id}->{visited_reactions}};
	  my $dead_end_cpd;

	  print STDERR "\tChecking if '$token_id' can run\n" unless !$debug;

	  my $found_reaction_for_token = 0;

	substrate: foreach my $cpd (keys %compounds_to_tokens)
	  {
	      if ($compounds_to_tokens{$cpd}->{$token_id} > 0 &&
		  ! defined $path_outputs{$cpd})
	      {
		  print STDERR "\t\tChecking substrate $cpd (main: $all_compounds_to_main{$cpd})\n" unless !$debug;

		  foreach my $reaction (@{$substrates_to_reactions->{$cpd}})
		  {
		      next if ($reaction =~ /(.*)_R/ && defined $visited_reactions{$1."_L"}) ||
			  ($reaction =~ /(.*)_L/ && defined $visited_reactions{$1."_R"});

		      print STDERR "\t\t\tChecking reaction $reaction\n" unless !$debug;

		      my $prods_are_ok = 0;
		  
		      foreach my $prod (@{$reactions_to_products->{$reaction}})
		      {
			  print STDERR "\t\t\t\tChecking product $prod\n" unless !$debug;

			  if (! defined $visited_compounds{$prod} || 
			      $visited_compounds{$prod} >= $visited_compounds{$cpd} ||
			      defined $path_outputs{$prod} ||
			      $compounds_borrowed_to_tokens{$prod}->{$token_id} > 0)
			  {
			      $prods_are_ok = $prod;
			      last;
			  }
		      }

		      if ($prods_are_ok)
		      {
			  print STDERR "\tToken '$token_id' can run $reaction on $cpd to produce $prods_are_ok\n" unless !$debug;
			  $found_reaction_for_token = 1;

			  if ($tokens{$token_id}->{initial_pass})
			  {
			      $found_reaction_for_initial_pass_token = 1;
			  }

			  next substrate;
		      }
		  }

		  # didn't find a reaction for this substrate
		  $dead_end_cpd = $cpd if $all_compounds_to_main{$cpd} && 
		      ! defined $path_outputs{$cpd};
	      }
	  }

	  if ($dead_end_cpd)
	  {
	      print STDERR "\tToken '$token_id' has reached a dead end on $dead_end_cpd\n" unless !$debug;
	      $tokens{$token_id}->{done} = "dead end on $dead_end_cpd";
	  }
	  elsif (! $found_reaction_for_token)
	  {
	      # didn't find any reaction to run.
	      # check to see if there are borrowed compounds to repay
	    product: foreach my $cpd (keys %compounds_borrowed_to_tokens)
	    {
		if ($compounds_borrowed_to_tokens{$cpd}->{$token_id} > 0)
		{
		    print STDERR "\t\tChecking product $cpd\n" unless !$debug;

		  rxn: foreach my $reaction (@{$products_to_reactions->{$cpd}})
		  {
		      next if ($reaction =~ /(.*)_R/ && defined $visited_reactions{$1."_L"})
			  || ($reaction =~ /(.*)_L/ && defined $visited_reactions{$1."_R"});

		      print STDERR "\t\t\tChecking reaction $reaction\n" unless !$debug;

		      my $substrates_are_ok = 0;
		      
		      foreach my $sub (@{$reactions_to_substrates->{$reaction}})
		      {
			  print STDERR "\t\t\t\tChecking substrate $sub\n" unless !$debug;

			  if (defined $path_outputs{$sub})
			  {
			      # don't run reactions that use up outputs
			      next rxn;
			  }

			  if (! defined $visited_compounds{$sub} || 
			      $visited_compounds{$sub} <= $visited_compounds{$cpd} ||
			      defined $scenario_cycles{$sub})
			  {
			      $substrates_are_ok = $sub;
			  }
		      }

		      if ($substrates_are_ok)
		      {
			  print STDERR "\tToken '$token_id' can wait for $reaction on $substrates_are_ok to produce borrowed compound $cpd\n" unless !$debug;
			  $found_reaction_for_token = 1;

			  if ($tokens{$token_id}->{initial_pass})
			  {
			      $found_reaction_for_initial_pass_token = 1;
			  }

			  last product;
		      }
		  }

		    # didn't find a reaction for this product
		    $dead_end_cpd = $cpd if $all_compounds_to_main{$cpd} && 
			! defined $path_outputs{$cpd};
		}
	    }

	      # didn't find any reaction to run.
	      if (! $found_reaction_for_token)
	      {
		  if ($dead_end_cpd)
		  {
		      print STDERR "\tToken '$token_id' has reached a dead end on borrowed compound $dead_end_cpd\n" unless !$debug;
		      $tokens{$token_id}->{done} = "dead end on borrowed compound $dead_end_cpd";
		  }
		  else
		  {
		      # nothing to push, borrow or do
		      print STDERR "\tToken '$token_id' has reached a dead end\n" unless !$debug;
		      $tokens{$token_id}->{done} = "dead end";
		  }
	      }
	  }
      }

	if($find_first)
	{
	    foreach my $token_id (sort { $tokens{$a}->{done} <=> $tokens{$b}->{done} } keys %tokens)
	    {
		if($tokens{$token_id}->{done} == 1)
		{
		    $done = 1;
		}
	    }
	}
	# is there an initial pass token that can make progress?
	if ($found_reaction_for_initial_pass_token)
	{
	    if (scalar keys %reactions_to_tokens_to_use == 0)
	    {
		# Couldn't run any reactions this time around.
		# Push more tokens through from the beginning of the path to
		# supply more substrates.
		$add_path_inputs = 1;
		print STDERR "\nSupplying more path inputs to push stalled tokens\n" unless !$debug;
	    }
	}
	else
	{
	    $done = 1;
	}

	$loop_counter++;

	if ($loop_counter >= $infinite_loop_check)
	{
	    $data_results{"infinite"} = 1;
	    print STDERR "Encountered an infinite loop\n" unless !$debug;
	    $done = 1;
	}
    }

    # reverse %compounds_to_tokens, since all tokens should be at path outputs now
    my %tokens_to_compounds;

    foreach my $cpd (keys %compounds_to_tokens)
    {
	foreach my $token_id (keys %{$compounds_to_tokens{$cpd}})
	{
      	    my $num_tokens = $compounds_to_tokens{$cpd}->{$token_id};
	    $tokens_to_compounds{$token_id}->{$cpd} = $num_tokens if $num_tokens > 0;
	}
    }    

    print STDERR "\n\ntoken ids: @{[ sort { $a <=> $b } map { $_ if ! defined $tokens{$_}->{done} } keys %tokens ]}\n" unless !$debug;

    my $path_counter = 1;

    foreach my $token_id (sort { $tokens{$a}->{done} <=> $tokens{$b}->{done} } keys %tokens)
    {
	my $token = $tokens{$token_id};
	my %visited_reactions = %{$token->{visited_reactions}};
	my @path = sort { $visited_reactions{$a} <=> $visited_reactions{$b} }
	    keys %visited_reactions;
	my %visited_compounds = %{$token->{visited_compounds}};
	my @compounds = sort { $visited_compounds{$a} <=> $visited_compounds{$b} }
	    keys %visited_compounds;

	print STDERR "Adding token id: $token_id\n" unless !$debug;


	#each key in data_results is a token which points to an array
	# [0]=initial pass [1]=0/1 if its done [2]=reaction path [3]=compounds 
	# [4]= html string of inputs
	# [5]= html string of outputs
	# [6]= html string of borrowed compounds
	# [7] = array of path input compounds
	# [8] = array of path output compounds


	$data_results{$token_id} = [$token->{initial_pass},$token->{done},\@path,\@compounds,[],[],[],$token->{token_path_inputs},$tokens_to_compounds{$token_id}];
	
	foreach my $input (keys %{$token->{token_path_inputs}})
	{
	    my $input_stoich = $token->{token_path_inputs}->{$input};
	    my $output_stoich = $tokens_to_compounds{$token_id}->{$input};

	    # don't balance scenario cycled compounds until final assembly
	    if ($scenario_cycles{$input} && ! $create_assembly)
	    {
		my @names = $fig->names_of_compound($input);
		push @{$data_results{$token_id}->[4]}, "\t\t$input_stoich\t$input $names[0]\n";
		next;
	    }

	    if ($input_stoich > $output_stoich)
	    {
		delete $tokens_to_compounds{$token_id}->{$input};
		$input_stoich -= $output_stoich;
		$token->{token_path_inputs}->{$input} -= $output_stoich;
		my @names = $fig->names_of_compound($input);
		push @{$data_results{$token_id}->[4]}, "\t\t$input_stoich\t$input $names[0]\n";
	    }
	    elsif ($output_stoich > $input_stoich)
	    {
		delete $token->{token_path_inputs}->{$input};
		$tokens_to_compounds{$token_id}->{$input} -= $input_stoich;
	    }
	    else
	    {
		delete $token->{token_path_inputs}->{$input};
		delete $tokens_to_compounds{$token_id}->{$input};
	    }
	}

	foreach my $output (keys %{$tokens_to_compounds{$token_id}})
	{
	    my @names = $fig->names_of_compound($output);
	    push @{$data_results{$token_id}->[5]},"\t\t$tokens_to_compounds{$token_id}->{$output}\t$output $names[0]\n";
	}
			
	if ($token->{done} != 1)
	{
	    foreach my $cpd (sort keys %compounds_borrowed_to_tokens)
	    {
		my $num = $compounds_borrowed_to_tokens{$cpd}->{$token_id};
		my @names = $fig->names_of_compound($cpd);
		push @{$data_results{$token_id}->[6]},"\t\t$num  $cpd\t$names[0]\n" if ($num > 0);
	    }
	}
    }

    
    

    return \%data_results;
}

sub balance_borrowing_and_giving
{
    my ($token_id, $compounds_to_tokens, $compounds_borrowed_to_tokens) = @_;

    my %merged_compounds;
    map { $merged_compounds{$_} = 1 } keys %{$compounds_to_tokens};
    map { $merged_compounds{$_} = 1 } keys %{$compounds_borrowed_to_tokens};

    foreach my $icpd (keys %merged_compounds)
    {
	my $inum_to_give = $compounds_to_tokens->{$icpd}->{$token_id};

	if ($inum_to_give > 0)
	{
	    my $inum_needed = $compounds_borrowed_to_tokens->{$icpd}->{$token_id};

	    if ($inum_to_give == $inum_needed)
	    {
		delete $compounds_borrowed_to_tokens->{$icpd}->{$token_id};
		delete $compounds_to_tokens->{$icpd}->{$token_id};
	    }
	    elsif ($inum_to_give > $inum_needed)
	    {
		delete $compounds_borrowed_to_tokens->{$icpd}->{$token_id};
		$compounds_to_tokens->{$icpd}->{$token_id} -= $inum_needed;
	    }
	    else
	    {
		$compounds_borrowed_to_tokens->{$icpd}->{$token_id} -= $inum_to_give;
		delete $compounds_to_tokens->{$icpd}->{$token_id};
	    }
	}
    }
}

sub print_token_status
{
    my ($token_id_list, $tokens, $compounds_to_tokens, $compounds_borrowed_to_tokens, $fig) = @_;

    print STDERR "\nToken status:\n" unless !$debug;

    foreach my $token_id (@$token_id_list)
    {
	next if defined $tokens->{$token_id}->{done};

	print STDERR "\n\ttoken: '$token_id', initial: $tokens->{$token_id}->{initial_pass}\n" unless !$debug;

	foreach my $cpd (sort keys %{$tokens->{$token_id}->{token_path_inputs}})
	{
	    my $num = $tokens->{$token_id}->{token_path_inputs}->{$cpd};
	    my @names = $fig->names_of_compound($cpd);
	    print STDERR "\t\tInput: $num  $cpd\t$names[0]\n" unless !$debug;
	}

	foreach my $cpd (sort keys %$compounds_to_tokens)
	{
	    my $num = $compounds_to_tokens->{$cpd}->{$token_id};
	    my @names = $fig->names_of_compound($cpd);
	    print STDERR "\t\tStatus: $num  $cpd\t$names[0]\n" if ($num > 0 && $debug);
	}

	foreach my $cpd (sort keys %$compounds_borrowed_to_tokens)
	{
	    my $num = $compounds_borrowed_to_tokens->{$cpd}->{$token_id};
	    my @names = $fig->names_of_compound($cpd);
	    print STDERR "\t\tBorrowed: $num  $cpd\t$names[0]\n" if ($num > 0 && $debug);
	}
	
	my %visited_compounds = %{$tokens->{$token_id}->{visited_compounds}};
	print STDERR "\t\tvisited_compounds: @{[map { ($_, $visited_compounds{$_} ) } sort { $visited_compounds{$a} <=> $visited_compounds{$b} } keys %visited_compounds]}\n" unless !$debug;

	my %visited_reactions = %{$tokens->{$token_id}->{visited_reactions}};
	print STDERR "\t\tvisited_reactions: @{[map { ($_, $visited_reactions{$_} ) } sort { $visited_reactions{$a} <=> $visited_reactions{$b} } keys %visited_reactions]}\n" unless !$debug;

    }

    print STDERR "\n" unless !$debug;
}

sub clone_token
{
    my ($clone_id, $new_token_id, $tokens, $compounds_to_tokens, $compounds_borrowed_to_tokens) = @_;
    my (%new_token, %new_visited_reactions, %new_visited_compounds, %new_token_path_inputs);

    $tokens->{$new_token_id} = \%new_token;

    my $clone_token = $tokens->{$clone_id};
    map {$new_visited_reactions{$_} = $clone_token->{visited_reactions}->{$_}} keys %{$clone_token->{visited_reactions}};
    map {$new_visited_compounds{$_} = $clone_token->{visited_compounds}->{$_}} keys %{$clone_token->{visited_compounds}};
    map {$new_token_path_inputs{$_} = $clone_token->{token_path_inputs}->{$_}} keys %{$clone_token->{token_path_inputs}};

    $new_token{visited_reactions} = \%new_visited_reactions;
    $new_token{visited_compounds} = \%new_visited_compounds;
    $new_token{token_path_inputs} =  \%new_token_path_inputs;
    $new_token{initial_pass} = $clone_token->{initial_pass};

    # tokens might be spread across multiple compounds
    foreach my $icpd (keys %$compounds_to_tokens)
    {
	if ($compounds_to_tokens->{$icpd}->{$clone_id} > 0)
	{
	    $compounds_to_tokens->{$icpd}->{$new_token_id} = 
		$compounds_to_tokens->{$icpd}->{$clone_id};
	}
    }

    # tokens might have borrowed compounds
    foreach my $icpd (keys %$compounds_borrowed_to_tokens)
    {
	if ($compounds_borrowed_to_tokens->{$icpd}->{$clone_id} > 0)
	{
	    $compounds_borrowed_to_tokens->{$icpd}->{$new_token_id} = 
		$compounds_borrowed_to_tokens->{$icpd}->{$clone_id};
	}
    }

    return \%new_token;
}

sub check_token_for_done
{
    my ($token_id, $compounds_to_tokens, $compounds_borrowed_to_tokens, $all_compounds_to_main, 
	$path_outputs, $scenario_cycles, $tokens, $fig, $outputs_lists) = @_;

    my $token_is_done_pushing = 1;
    my $token_is_done_borrowing = 1;

    # first determine if there is a main compound that isn't a path output
    foreach my $cpd (keys %$compounds_to_tokens)
    {
	if ($compounds_to_tokens->{$cpd}->{$token_id} > 0)
	{
	    # also check if scenario cycle compounds need to be pushed
	    if ($all_compounds_to_main->{$cpd} && 
		(! defined $path_outputs->{$cpd} || 
		 ($scenario_cycles->{$cpd} && 
		  $tokens->{$token_id}->{visited_compounds}->{$cpd} == 0)))
	    {
		my @names = $fig->names_of_compound($cpd);
		print STDERR "\ttoken '$token_id' needs to push $cpd $names[0]\n" unless !$debug;
		$token_is_done_pushing = 0;
	    }
	}
    }

    # now determine if one of the output lists has been satisfied
    if ($token_is_done_pushing)
    {
	my $found_a_list = 0;

	foreach my $cpd_list (@$outputs_lists)
	{
	    print STDERR "\t\tchecking outputs_list @$cpd_list\n" unless !$debug;

	    $found_a_list = 1;

	    foreach my $cpd (@$cpd_list)
	    {
		if (! $compounds_to_tokens->{$cpd}->{$token_id} > 0)
		{
		    $found_a_list = 0;
		    last;
		}
	    }

	    last if $found_a_list;
	}

	if (! $found_a_list)
	{
	    print STDERR "\ttoken '$token_id' hasn't satisfied output compound list\n" unless !$debug;
	    $token_is_done_pushing = 0;
	}
    }

    if ($token_is_done_pushing)
    {
	foreach my $cpd (keys %$compounds_borrowed_to_tokens)
	{
	    if ($compounds_borrowed_to_tokens->{$cpd}->{$token_id} > 0)
	    {
		# I don't know why I put in this "if" statement, so I'm going to do the same thing
		# in both cases, but print a distinguishing debug statement in case I need to find
		# the occurrences later
		if ($all_compounds_to_main{$cpd})
		{
		    my @names = $fig->names_of_compound($cpd);
		    print STDERR "\ttoken '$token_id' has borrowed $cpd $names[0] (main)\n" unless !$debug;
		    $token_is_done_borrowing = 0;
		}
		else
		{
		    my @names = $fig->names_of_compound($cpd);
		    print STDERR "\ttoken '$token_id' has borrowed $cpd $names[0] (not main)\n" unless !$debug;
		    $token_is_done_borrowing = 0;
		}
	    }
	}
    }

    if (! $token_is_done_pushing || ! $token_is_done_borrowing)
    {
	return 0;
    }
    else
    {
	$tokens->{$token_id}->{done} = 1;
	print STDERR "\tToken '$token_id' is done\n" unless !$debug;
	return 1;
    }
}


sub write_fluxanalyzer_files
{
    my ($dir, $path_inputs, $path_outputs, $path_array,
	$all_reactions,$reactions_to_substrate_arrays,$reactions_to_product_arrays,
	$cidToName) = @_;

    my $x_pos = 10;
    my $y_pos = 20;
   
    #Write the inputs/outputs to a seperate file, along with $stoich and if its main
    open(A_INPUT, ">$dir/inputs_main");
    
    print A_INPUT map {"$_\t$all_compounds_to_main{$_}\n"} keys %$path_inputs;
    
    close(A_INPUT);
    
    open(A_OUTPUT, ">$dir/outputs_main");
    
    print A_OUTPUT map {"$_\t$all_compounds_to_main{$_}\n"} keys %$path_outputs;
    
    close(A_OUTPUT);



    open(REACTIONS,  ">$dir/reactions");
    open(INPUTS, ">$dir/inputs");
    open(OUTPUTS, ">$dir/outputs");
    open(PATH, ">$dir/path_info");

    my @inputs = keys %$path_inputs;

    foreach my $elem(@{$path_array}){
	print PATH $elem ."\n";
    }
    close(PATH);

    foreach my $cpd (@inputs)
    {
	my $toPrint = $cpd."up\t = 1 $cpd \t| \t";
	$toPrint .= $path_inputs->{$cpd};
	$toPrint .= " \t0 100  0 \t$x_pos $y_pos 1 1\t0.01\n";
	print REACTIONS $toPrint;
	$y_pos += 20;

	if ($y_pos == 300)
	{
	    $x_pos += 60;
	    $y_pos = 20;
	}
    }

    $x_pos += 60;
    $y_pos = 20;

    my (@display_array);
    foreach my $rxn (keys %$all_reactions)
    {
	my $direction = $all_reactions->{$rxn};

	my (@substrate_array, @product_array);

	push @display_array, $rxn;

	if ($direction eq "L")
	{
	    @product_array = @{$reactions_to_substrate_arrays->{$rxn."_L"}};
	    @substrate_array = @{$reactions_to_product_arrays->{$rxn."_L"}};
	}
	else
	{
	    @substrate_array = @{$reactions_to_substrate_arrays->{$rxn."_R"}};
	    @product_array = @{$reactions_to_product_arrays->{$rxn."_R"}};
	}
	
	foreach my $subTuple (@substrate_array)
	{
	    my @temp = $fig->names_of_compound($subTuple->[0]);
	    $cidToName->{$subTuple->[0]} = $temp [0] if ! defined $cidToName->{$subTuple->[0]};
	}

	foreach my $prodTuple(@product_array)
	{
	    my @temp = $fig->names_of_compound($prodTuple->[0]);
	    $cidToName->{$prodTuple->[0]} = $temp[0] if ! defined $cidToName->{$prodTuple->[0]};
	}

	if($direction eq "R" || $direction eq "B") 
	{
	    #write data in a strign for copying to file later
	    my $toFile = '';
	    $toFile.= $rxn."\t";
	    
	    #add all the substrates
	    foreach my $curSub(@substrate_array){
		$toFile .= $curSub -> [1].' '. $curSub -> [0].' + ';
	    }

	    ##chop off the +
	    chop($toFile);
	    chop($toFile);

	    $toFile.='= ';

	    #add all the products 
	    foreach my $curProd(@product_array){
		$toFile .= $curProd -> [1].' '. $curProd -> [0].' + ';
	    }

	    #chop off the plus
	    chop($toFile);
	    chop($toFile);	

	    $toFile.="\t|\t#\t";

	    if($direction eq "B"){
		$toFile.="-Inf";
	    }
	    else{
		$toFile.="0";
	    }

	    $toFile.=" Inf 0\t$x_pos $y_pos 1 1 \t0.01\n";
	    
	    print REACTIONS $toFile;
	}
	elsif($direction eq "L")
	{
	    #write data in a strign for copying to file later
	    my $toFile = '';
	    $toFile.= $rxn."\t";
	    
	    #add all the substrates
	    foreach my $curProd(@product_array){
		$toFile .= $curProd -> [1].' '. $curProd -> [0]." + ";
	    }
	    
	    ##chop off the +
	    chop($toFile);
	    chop($toFile);
	    
	    $toFile.="= ";
	    
	    #add all the products 
	    foreach my $curSubstrate(@substrate_array){
		$toFile .= $curSubstrate -> [1].' '. $curSubstrate -> [0]." + ";
	    }
	    
	    #chop off the plus
	    chop($toFile);
	    chop($toFile);	
	    
	    $toFile.="\t|\t#\t";
	    
	    $toFile.="0";
	    
	    $toFile.=" Inf 0\t$x_pos $y_pos 1 1 \t0.01\n";
	    
	    print REACTIONS $toFile;
	}

	$y_pos += 20;	

	if ($y_pos == 300)
	{
	    $x_pos += 60;
	    $y_pos = 20;
	}
    }
    
    my @outputs = keys %$path_outputs;

    $x_pos += 60;
    $y_pos = 20;

    foreach my $cpd (@inputs)
    {
	print INPUTS $cpd, "\t", $path_inputs->{$cpd}, "\t", $cidToName->{$cpd}, "\n";
    }

    foreach my $cpd (@outputs)
    {
	print OUTPUTS $cpd, "\t", $path_outputs->{$cpd}, "\t", $cidToName->{$cpd}, "\n";
	my $toPrint = $cpd."ex\t 1 $cpd = \t| \t# \t0 100  0 \t$x_pos $y_pos 1 1\t0.01\n";
	print REACTIONS $toPrint;
	$y_pos += 20;

	if ($y_pos == 300)
	{
	    $x_pos += 60;
	    $y_pos = 20;
	}
    }

    $x_pos += 60;
    $y_pos = 20;

    #print the macromolecule_synthesis and assembly file
    open(MACRO_SYTH,">$dir/macromolecule_synthesis");
    open(ASSEM,">$dir/assembly");
    my $toPrint = "M1 = ";

    foreach my $cpd (keys %$path_outputs)
    {
	$toPrint.="$path_outputs->{$cpd} $cpd + ";
	print ASSEM "$cpd\tM1\t-100 -100 1\n";
	$y_pos += 25;
    }

    chop $toPrint;
    chop $toPrint;
    chop $toPrint;
    print MACRO_SYTH $toPrint;
    close(MACRO_SYTH);
    close(ASSEM);

    #Print the metabolites for these subsystems.
    open(METABOLITES,">$dir/metabolites");

    foreach my $cid (keys %$cidToName)
    {
	my $name = $cidToName->{$cid};
	$name =~ s/\s/-/g;
	print METABOLITES $cid."\t".$name."\t0.001\t0\n";
    }
    close(METABOLITES);

    #Print the macromolucules file
    open(MACRO,">$dir/macromolecules");
    print MACRO "M1 \tM1 \t1 \t-100 -100  1 1\n";
    close(MACRO);

    $x_pos += 60;
    $y_pos = 20;

    print REACTIONS "mue\t\t\t|\t#\t0  100  0\t$x_pos $y_pos 1 1\t0.01\n";

    #close reaction equation file
    close(REACTIONS);

    # FluxAnalyzer requires this file
    open(APP, ">$dir/app_para.m");
    print APP "epsilon=1e-10;\nbasic_color=[0.7         0.7         0.7];\ncr_color=[0.5         0.5           1];\nbr_color=[1         0.2         0.2];\nnbr_color=[0.2           1         0.2];\ntext_color=[0  0  0];\nmacro_synth_color=[0  0  1];\nmacro_color=[0.6         0.6           1];\nbox_reaction_width=[0.12];\nbox_reaction_height=[0.06];\nbox_macro_width=[0.08];\nbox_macro_height=[0.06];\nfontsize_reaction=[11];\nfontsize_macro=[11];\nfluxmaps={'Fluxmap','dummy.pcx'};\n";
    close(APP);
}

sub write_final_fluxanalyzer_files
{
    my ($dir, $path_inputs, $path_outputs, $all_reactions, $transport_reactions,
	$reactions_to_substrate_arrays,$reactions_to_product_arrays,
	$cidToName,$bioMass,$minSubstrates) = @_;

    #Write the inputs/outputs to a seperate file, along with $stoich and if its main
    open(A_INPUT, ">$dir/inputs_main");
    
    print A_INPUT map {"$_\t$all_compounds_to_main{$_}\n"} keys %$path_inputs;
    
    close(A_INPUT);
    
    open(A_OUTPUT, ">$dir/outputs_main");
    
    print A_OUTPUT map {"$_\t$all_compounds_to_main{$_}\n"} keys %$path_outputs;
    
    close(A_OUTPUT);

    my %open_transports;

    foreach my $cpd (keys %$minSubstrates)
    {
	map { $open_transports{$_} = 1 } @{$minSubstrates->{$cpd}};
    }

    open(REACTIONS,  ">$dir/reactions");
    open(INPUTS, ">$dir/inputs");
    open(OUTPUTS, ">$dir/outputs");

    my @inputs = keys %$path_inputs;

    my $x_pos = 10;
    my $y_pos = 30;
   
    foreach my $cpd (@inputs)
    {
	my $toPrint = $cpd."up\t = 1 $cpd \t| \t";

	if (defined $minSubstrates->{$cpd})
	{
	    $toPrint .= "#";
	}
	else
	{
	    $toPrint .= "0";
	}

	$toPrint .= " \t0 Inf  0 \t$x_pos $y_pos 1 1\t0.01\n";
	print REACTIONS $toPrint;
	$y_pos += 30;

	if ($y_pos > 600)
	{
	    $x_pos += 60;
	    $y_pos = 30;
	}
    }

    $x_pos += 60;
    $y_pos = 30;

    foreach my $rxn (keys %$all_reactions)
    {
	my $direction = $all_reactions->{$rxn};

	my (@substrate_array, @product_array);

	if ($direction eq "L")
	{
	    @product_array = @{$reactions_to_substrate_arrays->{$rxn."_L"}};
	    @substrate_array = @{$reactions_to_product_arrays->{$rxn."_L"}};
	}
	else
	{
	    @substrate_array = @{$reactions_to_substrate_arrays->{$rxn."_R"}};
	    @product_array = @{$reactions_to_product_arrays->{$rxn."_R"}};
	}
	
	foreach my $subTuple (@substrate_array)
	{
	    my @temp = $fig->names_of_compound($subTuple->[0]);
	    $cidToName->{$subTuple->[0]} = $temp [0] if ! defined $cidToName->{$subTuple->[0]};
	}

	foreach my $prodTuple(@product_array)
	{
	    my @temp = $fig->names_of_compound($prodTuple->[0]);
	    $cidToName->{$prodTuple->[0]} = $temp[0] if ! defined $cidToName->{$prodTuple->[0]};
	}

	if($direction eq "R" || $direction eq "B") 
	{
	    #write data in a strign for copying to file later
	    my $toFile = '';
	    $toFile.= $rxn."\t";
	    
	    #add all the substrates
	    foreach my $curSub(@substrate_array){
		$toFile .= $curSub -> [1].' '. $curSub -> [0].' + ';
	    }

	    ##chop off the +
	    chop($toFile);
	    chop($toFile);

	    $toFile.='= ';

	    #add all the products 
	    my $found_prod = 0;
	    foreach my $curProd(@product_array){
		$toFile .= $curProd -> [1].' '. $curProd -> [0].' + ';
		$found_prod = 1;
	    }

	    if ($found_prod)
	    {
		#chop off the plus
		chop($toFile);
		chop($toFile);	
	    }

	    if (defined $transport_reactions->{$rxn})
	    {
		if (defined $open_transports{$rxn})
		{
		    $toFile.="\t|\t#\t";
		}
		else
		{
		    $toFile.="\t|\t0\t";
		}

		if($direction eq "B"){
		    $toFile.="-Inf";
		}
		else{
		    $toFile.="0";
		}

		$toFile.=" Inf 0\t$x_pos $y_pos 1 1 \t0.01\n";
		$y_pos += 30;	
	    }
	    elsif ($rxn =~ /sink/)
	    {
		$toFile.="\t|\t0\t";

		if($direction eq "B"){
		    $toFile.="-0.00001";
		}
		else{
		    $toFile.="0";
		}

		$toFile.=" 0.00001 0\t$x_pos $y_pos 1 1 \t0.01\n";
		$y_pos += 30;	
	    }
	    else
	    {
		$toFile.="\t|\t#\t";

		if($direction eq "B"){
		    $toFile.="-Inf";
		}
		else{
		    $toFile.="0";
		}

		$toFile.=" Inf 0\t-10 -10 1 1 \t0.01\n";
	    }		
	    
	    print REACTIONS $toFile;
	}
	elsif($direction eq "L")
	{
	    #write data in a strign for copying to file later
	    my $toFile = '';
	    $toFile.= $rxn."\t";
	    
	    #add all the substrates
	    foreach my $curProd(@product_array){
		$toFile .= $curProd -> [1].' '. $curProd -> [0]." + ";
	    }
	    
	    ##chop off the +
	    chop($toFile);
	    chop($toFile);
	    
	    $toFile.="= ";
	    
	    #add all the products 
	    foreach my $curSubstrate(@substrate_array){
		$toFile .= $curSubstrate -> [1].' '. $curSubstrate -> [0]." + ";
	    }
	    
	    #chop off the plus
	    chop($toFile);
	    chop($toFile);	
	    
	    
	    $toFile.="\t|\t#\t";
	    $toFile.="0";
	    $toFile.=" Inf 0\t-10 -10 1 1 \t0.01\n";
	    
	    print REACTIONS $toFile;
	}

	if ($y_pos > 600)
	{
	    $x_pos += 60;
	    $y_pos = 30;
	}
    }
    
    my @outputs = keys %$path_outputs;

    $x_pos += 60;
    $y_pos = 30;

    foreach my $cpd (@inputs)
    {
	print INPUTS $cpd, "\t", $path_inputs->{$cpd}, "\t", $cidToName->{$cpd}, "\n";
    }

    foreach my $cpd (@outputs)
    {
	print OUTPUTS $cpd, "\t", $path_outputs->{$cpd}, "\t", $cidToName->{$cpd}, "\n";
	my $toPrint = $cpd."ex\t 1 $cpd = \t| \t# \t0 Inf  0 \t$x_pos $y_pos 1 1\t0.01\n";
	print REACTIONS $toPrint;
	$y_pos += 30;

	if ($y_pos > 600)
	{
	    $x_pos += 60;
	    $y_pos = 30;
	}
    }

    $x_pos += 60;
    $y_pos = 30;

    #print the macromolecule_synthesis and assembly file
    open(MACRO,">$dir/macromolecules");
    open(MACRO_SYTH,">$dir/macromolecule_synthesis");
    open(ASSEM,">$dir/assembly");

    my $toPrint = "M1 = ";

    foreach my $cpd (keys %$bioMass)
    {
	$toPrint .= "$bioMass->{$cpd} $cpd + ";
	print ASSEM "$cpd\tM1\t-100 -100 1\n";
    }

    chop $toPrint;
    chop $toPrint;
    chop $toPrint;

    print MACRO "M1 \tM1 \t1 \t-100 -100  1 1\n";
    print MACRO_SYTH $toPrint, "\n";

    close(MACRO);
    close(MACRO_SYTH);
    close(ASSEM);

    #Print the metabolites for these subsystems.
    open(METABOLITES,">$dir/metabolites");

    foreach my $cid (keys %$cidToName)
    {
	my $name = $cidToName->{$cid};
	$name =~ s/\s/-/g;
	print METABOLITES $cid."\t".$name."\t0.001\t0\n";
    }
    close(METABOLITES);

    $x_pos += 60;
    $y_pos = 30;

    print REACTIONS "mue\t\t\t|\t#\t0  100  0\t$x_pos $y_pos 1 1\t0.01\n";

    #close reaction equation file
    close(REACTIONS);

    # FluxAnalyzer requires this file
    open(APP, ">$dir/app_para.m");
    print APP "epsilon=1e-10;\nbasic_color=[0.7         0.7         0.7];\ncr_color=[0.5         0.5           1];\nbr_color=[1         0.2         0.2];\nnbr_color=[0.2           1         0.2];\ntext_color=[0  0  0];\nmacro_synth_color=[0  0  1];\nmacro_color=[0.6         0.6           1];\nbox_reaction_width=[0.06];\nbox_reaction_height=[0.03];\nbox_macro_width=[0.06];\nbox_macro_height=[0.03];\nfontsize_reaction=[11];\nfontsize_macro=[11];\nfluxmaps={'Fluxmap','dummy_medium.pcx'};\n";
    close(APP);
}


sub clear_arrays
{
    undef %reactions_to_substrate_arrays;
    undef %reactions_to_product_arrays;
    undef %all_compounds_to_main;
    undef %all_reactions;
    undef %scenario_cycles;
    undef @all_outputs_lists;
    undef %all_inputs;
    undef %all_outputs;
    
    %reactions_to_substrate_arrays = ();
    %reactions_to_product_arrays = ();
    %all_compounds_to_main = ();
    %all_reactions = ();
    %scenario_cycles = ();
    @all_outputs_lists = ();
    %all_inputs = ();
    %all_outputs = ();
}

sub load_superset_file
{
    # only need to load once during execution context
    return \%ss_to_superset if %ss_to_superset;

    %superset_to_ss = ();
    %ss_to_superset = ();    
    
    if (open(FILE,"<$FIG_Config::global/Models/hope_supersets.txt"))
    {
	while(<FILE>)
	{
	    my @line = split(/\t/,$_);
	    map { s/ /_/g } @line;
	    map { s/\"//g } @line;
	    map { chomp } @line;
	    $superset_to_ss{$line[0]} = [] if !defined $superset_to_ss{$line[0]};
	    $ss_to_superset{$line[1]} = $line[0];
	    push(@{$superset_to_ss{$line[0]}},$line[1]);
	}
	close FILE;
    }

    return \%ss_to_superset
}

# This function runs a given scenario that is defined the specified subsystem
# It returns the data as it is found by process_paths in a hash reference
# This is used internally by model.pm, and shouldn't be called externally

sub load_scenario
{
    my ($genome,$ssa,$scenario) = @_;

    #load up the arrays with the info we need
    process_init($ssa,$scenario,$genome,0);

    #assume all path inputs and outputs are main
    map { $all_compounds_to_main{$_} = 1 } keys %all_inputs;
    map { $all_compounds_to_main{$_} = 1 } keys %all_outputs;

}

sub internal_scenario
{
    my ($genome,$ssa,$scenario,$find_first) = @_;

    print STDERR "\nIn internal_scenario with '$genome', '$ssa', '$scenario', '$find_first'\n" if $debug;

    load_scenario($genome,$ssa,$scenario);
    
    return execute_paths([],$find_first,[],[]);
}

sub run_scenario
{
    my($genome,$superset,$subsystem,$scenario,$find_first) = @_;
    my $scenario_dir = get_scenario_directory($genome) . "/PathInfo/$superset/$subsystem/$scenario";

    system("rm", "-rf", $scenario_dir);
    &FIG::verify_dir($scenario_dir);
    
    #make sure the arrays are empty to start out
    &clear_arrays;
    write_scenario(internal_scenario($genome,$subsystem,$scenario,$find_first),$scenario_dir);
}

sub compare_scenario
{
    my($genome,$superset,$ss_name,$scenario_name,$dont_copy) = @_;
    my @genome_paths;
    my $scenario_dir_all = get_scenario_directory('All') . "/PathInfo/$superset/$ss_name/$scenario_name";
    my $subsystem = $fig->get_subsystem($ss_name);
    my @additional_reactions = $subsystem->get_hope_additional_reactions($scenario_name);
    my %additional_reactions;
    map { $additional_reactions{$_} = 1 } @additional_reactions;

    my %ss_reactions;    

    if ($genome eq "All")
    {
	my %all_reactions = $subsystem->get_hope_reactions;
	foreach my $role (keys %all_reactions)
	{
	    map { $ss_reactions{$_} = 1 } @{$all_reactions{$role}};
	}
    }
    else
    {
	my $reactions_for_genome = get_reactions_for_genome_in_subsystem($genome,$ss_name);
	map { $ss_reactions{$_} = 1 } keys %$reactions_for_genome if keys %$reactions_for_genome;
    }

    # first find paths in the All directory that should be valid for the genome
    # based on the reactions associated with it in the subsystems
    opendir (DIR_ALL,$scenario_dir_all) or return [];
    my @sub_dirs  = readdir DIR_ALL;
    close DIR_ALL;

    my %paths_all;

    for my $path (@sub_dirs)
    {
	next if $path !~ /path/; # skip . and ..
	my $match = 1;
	open (PATH, "$scenario_dir_all/$path/path_info");
	my @reactions = <PATH>;
	close PATH;

	my $reaction_string = "";

	foreach my $reaction (sort @reactions)
	{
	    if ($reaction =~ /(R\d\d\d\d\d)/)
	    {
		$reaction_string .= $reaction;

		if (! exists($ss_reactions{$1}) && ! exists($additional_reactions{$1}))
		{
		    $match = 0;
		}
	    }
	}
	
	$paths_all{$reaction_string} = $match;
	push @genome_paths, $path if $match;
    }

    # now check all paths found for this particular organism, and make sure they
    # are in the appropriate All subdirectory
    my $scenario_dir_genome = get_scenario_directory($genome). "/PathInfo/$superset/$ss_name/$scenario_name";

    opendir (DIR_GENOME,$scenario_dir_genome) or return \@genome_paths;
    @sub_dirs  = readdir DIR_GENOME;
    close DIR_GENOME;

    my $path_counter = scalar keys %paths_all; 

    for my $path (@sub_dirs)
    {
	next if $path !~ /path/; # skip . and ..
	my $match = 1;
	open (PATH, "$scenario_dir_genome/$path/path_info");
	my @reactions = <PATH>;
	close PATH;

	my $reaction_string = "";

	foreach my $reaction (sort @reactions)
	{
	    if ($reaction =~ /(R\d\d\d\d\d)/)
	    {
		$reaction_string .= $reaction;
	    }
	}
	
	if (! exists($paths_all{$reaction_string}))
	{
	    if ($dont_copy)
	    {
		print STDERR "$scenario_dir_genome/$path not found in All\n";
	    }
	    else
	    {
		$path_counter++;
		my $new_path_name = "path_".$path_counter;
		my $temp_sdg = $scenario_dir_genome;
		$temp_sdg =~ s/\(/\\\(/g;
		$temp_sdg =~ s/\)/\\\)/g;
		my $temp_sda = $scenario_dir_all;
		$temp_sda =~ s/\(/\\\(/g;
		$temp_sda =~ s/\)/\\\)/g;
		`cp -R $temp_sdg/$path $temp_sda/$new_path_name`;
		push @genome_paths, $new_path_name;
		print STDERR "Copied $temp_sdg/$path to $temp_sda/$new_path_name\n";
	    }
	}

	unless ($dont_copy)
	{
	    # remove genome-specific paths
	    rmtree("$scenario_dir_genome/$path");
	}
    }

    unless ($dont_copy)
    {
	# copy each genome-specific path, with the same name as the path in the "All" directory
	foreach my $path (@genome_paths)
	{
	    my $temp_sdg = $scenario_dir_genome;
	    $temp_sdg =~ s/\(/\\\(/g;
	    $temp_sdg =~ s/\)/\\\)/g;
	    my $temp_sda = $scenario_dir_all;
	    $temp_sda =~ s/\(/\\\(/g;
	    $temp_sda =~ s/\)/\\\)/g;
	    `cp -r $temp_sda/$path $temp_sdg`;
	}	
    }

    return \@genome_paths;
}

sub write_scenario
{
    my($scenario_data,$scenario_dir) = @_;
    delete $scenario_data->{"infinite"};
    my $path_count = 1;
    my @list_of_done_tokens=();

    print STDERR "Paths: ", keys %{$scenario_data}, "\n" if $debug;

    foreach my $try_path (keys %{$scenario_data})
    {
	print STDERR "\t Checking $try_path\n" if $debug;
	if($scenario_data->{$try_path}->[1] != 1) 
	{
	    print STDERR "\t Token $try_path is not complete\n" if $debug;
	    next;
	}
	print STDERR "These are the contents of the list: " , @list_of_done_tokens, "\n" if $debug;
	if(scalar @list_of_done_tokens == 0){
	    push @list_of_done_tokens, $try_path;
	    next;
	}	
	
	# Check this token's ($try_path) values against the values of the keys($elem) stored in 
	# list_of_done_tokens. If they match, don't add it to the finished token list, if it doesn't
	# match, add it.

	print STDERR "These are the contents of the list: " , @list_of_done_tokens, "\n" if $debug;

	my $found_match = 0;

	foreach my $elem(@list_of_done_tokens){
	    my @done_reactions = @{$scenario_data->{$elem}->[2]};
	    
	    print STDERR "This is the path we're trying: " , $try_path, "\t","This is the path already in the array: " ,  $elem , "\nThis is the size of the array: " . @list_of_done_tokens ."\n" if $debug;

	    # if the list of reactions match, they represent the same path

	    my @path_reactions =  @{$scenario_data->{$try_path}->[2]};
	    my (%diff_reactions_1, %diff_reactions_2);
	    
	    map {$diff_reactions_1{$_} = 1} @path_reactions;
	    map {delete $diff_reactions_1{$_}} @done_reactions;
	    map {$diff_reactions_2{$_} = 1} @done_reactions;
	    map {delete $diff_reactions_2{$_}} @path_reactions;

	    if (scalar keys %diff_reactions_1 == 0 && scalar keys %diff_reactions_2 == 0)
	    {
		print STDERR "They match.\n" if $debug;
		$found_match = 1;
		last;
	    }		
	}

	if (! $found_match)
	{
	    push @list_of_done_tokens, $try_path;
	    print STDERR $try_path, " Added to the array\n" if $debug;

	}
    }

    foreach my $path (@list_of_done_tokens){
	if($scenario_data->{$path}->[1] != 1) 
	{
	    next;
	}
	if(@list_of_done_tokens == 0){
	    push @list_of_done_tokens, $path;
	}	
		
	#create input/output info
	my $input_hash = $scenario_data->{$path}->[7];
	my $output_hash = $scenario_data->{$path}->[8];
	my $reaction_path = $scenario_data->{$path}->[2];
	my @reaction_array = @$reaction_path;
	
 	print STDERR "\nInputs:\n" if $debug;
	print STDERR map{"$_ => $input_hash->{$_}" } keys %$input_hash, "\n" if $debug;
	print STDERR "\nOutputs:\n" if $debug;
	print STDERR map{"$_ => $output_hash->{$_}"} keys %$output_hash, "\n" if $debug;

	# divide stoichiometry by greatest common denominator
	my ($min_stoich, @all_stoichs);
	
	map { push @all_stoichs, $input_hash->{$_}; $min_stoich = $input_hash->{$_} if $input_hash->{$_} < $min_stoich || $min_stoich == 0 } keys %{$input_hash};
	map { push @all_stoichs, $output_hash->{$_}; $min_stoich = $output_hash->{$_} if $output_hash->{$_} < $min_stoich || $min_stoich == 0 } keys %{$output_hash};
	
	my ($gcd, @gcd_candidates);
	
      outer: for ($gcd = $min_stoich; $gcd > 1; $gcd--)
      {
	  foreach my $stoich (@all_stoichs)
	  {
	      next outer if $stoich % $gcd != 0;
	  }
	  
	  last; # found a gcd
      }
	
	map { $input_hash->{$_} /= $gcd } keys %{$input_hash};
	map { $output_hash->{$_} /= $gcd } keys %{$output_hash};
	
	mkdir "$scenario_dir/path_$path_count";
	
	open(FILE, ">$scenario_dir/path_$path_count/path_info");
	foreach my $elem(@reaction_array){
	    print FILE scalar @reaction_array, "\t", $elem, "\n";
	}
	close(FILE);
	
	&write_fluxanalyzer_files("$scenario_dir/path_$path_count",$input_hash,
				  $output_hash, \@reaction_array,\%all_reactions,
				  \%reactions_to_substrate_arrays,
				  \%reactions_to_product_arrays,
				  {});
	$path_count++;
    }
    print STDERR "no paths for $scenario_dir\n" if $path_count == 1 && $debug;
    @list_of_done_tokens=();
    undef %{$scenario_data};

}


sub load_subsystem
{
    my ($genome,$ss_name) = @_;
    my $subsystem = $fig->get_subsystem($ss_name);
    my @ss_scenarios = $subsystem->get_hope_scenario_names;
    foreach my $name (@ss_scenarios)
    {
	load_scenario($genome,$ss_name,$name);
    }
}

sub internal_subsystem
{
    my ($genome,$ss_name,$find_first) = @_;
    my %scenario_to_paths;
    
    my $subsystem = $fig->get_subsystem($ss_name);
    my @ss_scenarios = $subsystem->get_hope_scenario_names;
    
    foreach my $name (@ss_scenarios)
    {
	$scenario_to_paths{$name} = internal_scenario($genome,$ss_name,$name,$find_first);
    }
    
    return \%scenario_to_paths;
}

sub run_subsystem
{
    my ($genome,$superset,$subsystem,$find_first) = @_;
    
    my $subsystem_obj = $fig->get_subsystem($subsystem);

    if (!$subsystem_obj)
    {
	warn "Cannot open subsystem $subsystem\n";
	return;
    }
    
    my @ss_scenarios = $subsystem_obj->get_hope_scenario_names;

    my $dir = get_scenario_directory($genome) . "/PathInfo/$superset/$subsystem";
    system("rm", "-rf",  $dir);
    &FIG::verify_dir($dir);

    foreach my $name (@ss_scenarios)
    {
	print STDERR "\t Running $name\n";
	run_scenario($genome,$superset,$subsystem,$name,$find_first);
    }
}

sub compare_subsystem
{
    my ($genome,$superset,$subsystem,$dont_copy) = @_;
    my %genome_scenarios;

    my $subsystem_obj = $fig->get_subsystem($subsystem);

    if (!$subsystem_obj)
    {
	warn "Cannot open subsystem $subsystem\n";
	return;
    }

    my @ss_scenarios = $subsystem_obj->get_hope_scenario_names;
    
    foreach my $name (@ss_scenarios)
    {
	$genome_scenarios{$name} = compare_scenario($genome,$superset,$subsystem,$name,$dont_copy);
    }

    return \%genome_scenarios;
}

sub load_superset
{
    my($genome, $superset_name) = @_;

    my @subsystems = @{$superset_to_ss{$superset_name}};
    foreach my $ss_name (@subsystems)
    {
	load_subsystem($genome,$ss_name);
    }
}

sub internal_superset
{
    my($genome, $superset_name,$find_first) = @_;
    
    my @subsystems = @{$superset_to_ss{$superset_name}};

    my %supersets_data;
    
    foreach my $ss_name (@subsystems)
    {
	$supersets_data{$ss_name} = internal_subsystem($genome,$ss_name,$find_first);
    }

    return \%supersets_data;
}

sub run_superset
{
    my($genome, $superset_name,$find_first) = @_;
    
    my @subsystems = @{$superset_to_ss{$superset_name}};
    

    my $dir = get_scenario_directory($genome) . "/PathInfo/$superset_name";
    system("rm", "-rf", $dir);
    &FIG::verify_dir($dir);

    foreach my $ss_name (@subsystems)
    {
	print STDERR "Running Scenarios for $genome for subsystem $ss_name\n";
	run_subsystem($genome,$superset_name,$ss_name,$find_first);
    }
}

sub compare_superset
{
    my($genome, $superset_name, $dont_copy) = @_;
    
    my @subsystems = @{$superset_to_ss{$superset_name}};
    my %genome_subsystems;

    foreach my $ss_name (@subsystems)
    {
	print STDERR "Comparing Scenarios for $genome in subsystem $ss_name\n";
	$genome_subsystems{$ss_name} = compare_subsystem($genome,$superset_name,$ss_name,$dont_copy);
    }

    return \%genome_subsystems;
}


sub load_supersets
{
    my($genome) = @_;
    &load_superset_file;

    foreach my $superset (keys %superset_to_ss)
    {
	load_superset($genome,$superset);
    }
    
    return (\%all_reactions,\%reactions_to_substrate_arrays,\%reactions_to_product_arrays);

}

sub run_supersets
{
    my($genome,$find_first) = @_;
    &load_superset_file;

    my $dir = get_scenario_directory($genome) . "/PathInfo";
    system("rm", "-rf", $dir);
    &FIG::verify_dir($dir);

    foreach my $superset (keys %superset_to_ss)
    {
	run_superset($genome,$superset,$find_first);
    }

}

sub compare_supersets
{
    my($genome, $dont_copy) = @_;
    my %genome_supersets;

    &load_superset_file;

    foreach my $superset (keys %superset_to_ss)
    {
	$genome_supersets{$superset} = compare_superset($genome,$superset,$dont_copy);
    }
    
    return \%genome_supersets;
}

sub run_genome_report
{
    my ($genome) = @_;
    my @string_out;
    push @string_out,"Genome $genome\n";
    #get all the subsystems this genome is involved in
    foreach my $superset (keys %superset_to_ss)
    {
	foreach my $name (@{$superset_to_ss{$superset}})
	{
	    push @string_out, @{print_ss_report($name,internal_subsystem($genome,$name,0))};
	}
    }
    return \@string_out;
}

sub print_ss_report
{
    my ($ss_name,$scenario_to_paths) = @_;
    my @output = ();
    my %scenario_path_count;

    push(@output,"\tSubsystem $ss_name\n");
    
    foreach my $scenario (keys %$scenario_to_paths)
    {
	if($scenario_to_paths->{$scenario}->{"infinite"})
	{
	    push(@output,"\tWarning: Possible Infinite loop\n");
	}
	delete $scenario_to_paths->{$scenario}->{"infinite"};
	foreach my $token (keys %{$scenario_to_paths->{$scenario}})
	{
	    $scenario_path_count{$scenario}++ if ($scenario_to_paths->{$scenario}->{$token}->[1]);
	}
    }

    push @output , map { "\t\t$_ has $scenario_path_count{$_} path(s).\n" } keys %scenario_path_count;
   
    return \@output;
}


sub internal_assembly
{
    #lets get the genome, and a array reference to the paths we want to build togather
    my ($paths,$input_path,$output_path,$one_path) = @_;


    print STDERR $paths."\n" if $debug;
    clear_arrays();

    #This gets us an array of arrays, each subarray holds
    # [0] = genome   [1] = Scenarios [2] = superset   [3] = subsystem   [4] = scenario   [5] = path
    # OR [0] = genome [1] = assembly [2] = path_name
    my @assembly_scenarios = @{parse_assembly_scenarios($paths)};

    #split and the input/output path for later as well
    #these should only return one path array...so just grab that one
    my @input_arr;
    my @output_arr;
    
    if($input_path != undef  && $output_path != undef )
    {
	@input_arr = @{parse_assembly_scenarios($input_path)}; 
	@output_arr = @{parse_assembly_scenarios($output_path)};
    }

    #load all the kegg information for each 'scenario' from the paths we have selected
    foreach my $scenario (@assembly_scenarios)
    {
	if(scalar @$scenario > 5) #this is a normal scenario path
	{
	    print STDERR "Checking $scenario->[3] $scenario->[4] $scenario->[5] \n" if $debug; 
	    process_init($scenario->[3],$scenario->[4],$scenario->[0],1);
	}
	else #This is a assembly path
	{
	    #read in the input/output compounds and mark main's correctly
	    my $genome = shift @$scenario;
	    my $path = get_scenario_directory($genome) . "/" . join "/" , @$scenario;
	    
	    open(M_IN,"$path/inputs_main");
	    while(<M_IN>)
	    {
		my @line = split(/\t/,$_);
		$all_compounds_to_main{$line[0]} = $line[1];
	    }
	    close(M_IN);

	    open(M_OUT,"$path/outputs_main");
	    while(<M_OUT>)
	    {
		my @line = split(/\t/,$_);
		$all_compounds_to_main{$line[0]} = $line[1];
	    }
	    close(M_OUT);
	}
    }
    #assume all path inputs and outputs are main
    map { $all_compounds_to_main{$_} = 1 } keys %all_inputs;
    map { $all_compounds_to_main{$_} = 1 } keys %all_outputs;

    print STDERR "Inputs: " if $debug;
    print STDERR map { $_."\n" } keys %all_inputs if $debug;
    print STDERR "Outputs: " if $debug;
    print STDERR map { $_."\n" } keys %all_outputs if $debug;
    
    #run process paths
    return execute_paths(\@assembly_scenarios,$one_path,\@input_arr,\@output_arr);
}

sub run_assembly
{
    my ($paths,$genome,$write_name,$one_path) = @_;
    print STDERR "\nThis is the passed information.\n";
    print STDERR $paths , "\n" , @$paths, "\n";
    print STDERR $genome . "\n";
    print STDERR $write_name . "\n";
    my $dir = get_scenario_directory($genome) . "/Assemblies/$write_name";
    system("rm", "-rf", $dir);
    &FIG::verify_dir($dir);
    write_scenario(internal_assembly($paths,[],[],$one_path),$dir);
}

sub expand_paths
{
    my ($paths) = @_;
    
    my @final_paths;
    
    foreach my $path (@$paths)
    {
	if($path eq "" || $path eq "//")
	{
	    next;
	}
	my $length = 6;
	if($path =~ /Assemblies/)
	{
	    $length = 4;
	}
	my @parts = split "/", $path;
	shift @parts; # get ride of the first blank entry from /$genome
	$length =$length - scalar @parts;
	my $genome = shift @parts;
	$path = join "/" , @parts;
	my @temp = @{expand_recursive($genome, $path, $length)};
	push @final_paths , @temp if scalar @temp > 0;
    }


    return \@final_paths;
}

sub expand_recursive
{
    my ($genome,$path,$count) = @_;
    my @sub_dirs;
    if($path =~ m\Assemblies$\ || $path =~ m\Analysis$\ || $path =~ m\Curation$\)
    {
	return [];
    }
    if($count !=0)
    {
        #read this path, and pull out all the sub-directories.
	my $model_dir = get_scenario_directory($genome) . "/$path/";
	return [] if (! -d $model_dir);
	opendir (DIR, $model_dir) or die("$model_dir");
	@sub_dirs  = readdir DIR;
	close DIR;
    }   
    else
    {
	
	return [$path];
    }
    $count--;
    my @to_return;
    foreach my $sub_path (@sub_dirs)
    {
	next if $sub_path =~ /^\.$/ || $sub_path =~ /^\.\.$/;
	push @to_return , @{expand_recursive($genome, "$path/$sub_path",$count)};
    }
    return \@to_return;

}

sub parse_assembly_scenarios
{
    my ($paths) = @_;

    $paths = expand_paths($paths);
    my @array_of_path_arrays;

    foreach my $path (@$paths)
    {
	my @parts = split "/", $path;
	#shift @parts;
	push(@array_of_path_arrays, \@parts);
    }

    return \@array_of_path_arrays;
}


sub write_selected_scenarios
{
    my($checked,$genome,$ssa,$sc_name) = @_;
    my (@tempArray);

    #Load this scenario again with all of its rxns and cpds

    model::clear_arrays();
    
    model::process_init($ssa,$sc_name,$genome,0);
    
    #assume all path inputs and outputs are main
    map { $all_compounds_to_main{$_} = 1 } keys %all_inputs;
    map { $all_compounds_to_main{$_} = 1 } keys %all_outputs;

    model::create_reactions({},{},{},{});

    ##End of scenario loading


    #setup the filesystem to store the scenario/paths
    my $superset = $ss_to_superset{$ssa};
    my $base_dir = get_scenario_directory($genome) . "/PathInfo/$superset/$ssa/$sc_name/";
    system("rm", "-rf", $base_dir);
    &FIG::verify_dir($base_dir);

    #for the selected paths, lets gather their cpd from the checkbox and write 
    #the path that we need
    foreach my $path (@$checked)
    {
	my (%input_hash, %output_hash, $path_name);
	#process the strings to get the information from the parameters
	my @items = split(";", $path);
	$path_name = $items[0];
	#next we have the input compounds ids/stoich/main
	map { if ($_ =~ /(.*):(.*):(.*)/)
	      { $input_hash{$1}+= $2 } } split ",", $items[1];
	#the third part has the output compounds ids/stoich/main
	map { if ($_ =~ /(.*):(.*):(.*)/)
	      { $output_hash{$1} += $2 } } split ",",$items[2];
	#the fourth part has the strings of reactions visited
	@tempArray = split("#" , $items[3]);
	
	
 	print STDERR "\nInputs:\n" if $debug;
	print STDERR map{"$_ => $input_hash{$_}" } keys %input_hash, "\n" if $debug;
	print STDERR "\nOutputs:\n" if $debug;
	print STDERR map{"$_ => $output_hash{$_}"} keys %output_hash, "\n" if $debug;
	
	
	# divide stoichiometry by greatest common denominator
	my ($min_stoich, @all_stoichs);
	
	map { push @all_stoichs, $input_hash{$_}; $min_stoich = $input_hash{$_} if $input_hash{$_} < $min_stoich || $min_stoich == 0 } keys %input_hash;
	map { push @all_stoichs, $output_hash{$_}; $min_stoich = $output_hash{$_} if $output_hash{$_} < $min_stoich || $min_stoich == 0 } keys %output_hash;
	
	my ($gcd, @gcd_candidates);
	
      outer: for ($gcd = $min_stoich; $gcd > 1; $gcd--)
      {
	  foreach my $stoich (@all_stoichs)
	  {
	      next outer if $stoich % $gcd != 0;
	  }
  
	  last; # found a gcd
      }
	
	map { $input_hash{$_} /= $gcd } keys %input_hash;
	map { $output_hash{$_} /= $gcd } keys %output_hash;
	
	mkdir "$base_dir/$path_name";
	&write_fluxanalyzer_files("$base_dir/$path_name",\%input_hash,
				  \%output_hash,\@tempArray,\%all_reactions,
				  \%reactions_to_substrate_arrays,
				  \%reactions_to_product_arrays,
				  {});
    }
    return $base_dir;
}

#This write function assumes that we have just run a assembly (and we haven't cleared the arrays)
# becuase the write_fluxanalyzer_files function is depended on those global arrays for the reactions.

sub write_assembly
{
    my($input,$genome,$name) = @_;

    my $paths = $input->[1];
    my $file_paths = $input->[0];
    print STDERR "Paths: @$paths \n File_Dirs: @$file_paths \n" if $debug;
    my @tempArray;
    
    chomp $genome;
    chomp $name;
    #setup the filesystem to store the assembly
    my $base_dir = get_scenario_directory($genome) . "/Assemblies/$name";
    system("rm", "-rf", $base_dir);
    &FIG::verify_dir($base_dir);

    
    ##Here we want to reload the cpd and rxn info so we can write it later
    clear_arrays();

    #This gets us an array of arrays, each subarray holds
    # [0] = genome [1] = Scenarios  [2] = superset   [3] = subsystem   [4] = scenario   [5] = path
    my @assembly_scenarios = @{parse_assembly_scenarios($file_paths)};

    #load all the kegg information for each 'scenario' from the paths we have selected
    foreach my $scenario (@assembly_scenarios)
    {
	print STDERR "Checking $scenario->[3] $scenario->[4] $scenario->[5] \n" if $debug; 
	process_init($scenario->[3],$scenario->[4],$scenario->[0],1);
    }
    #assume all path inputs and outputs are main
    map { $all_compounds_to_main{$_} = 1 } keys %all_inputs;
    map { $all_compounds_to_main{$_} = 1 } keys %all_outputs;

    create_assembly_reactions({},{},{},{});

    ##End of rxn,cpd loading
    
    #for the selected paths, lets gather their cpd from the checkbox and write 
    #the path that we need
    foreach my $path (@$paths)
    {
	my (%input_hash, %output_hash, $path_name);
	#process the strings to get the information from the parameters
	my @items = split(";", $path);
	$path_name = $items[0];
	#next we have the input compounds ids/stoich/main
	map { if ($_ =~ /(.*):(.*):(.*)/)
	      { $input_hash{$1}+= $2 } } split ",", $items[1];
	#the third part has the output compounds ids/stoich/main
	map { if ($_ =~ /(.*):(.*):(.*)/)
	      { $output_hash{$1} += $2 } } split ",",$items[2];
	#the fourth part has the strings of reactions visited
	@tempArray = split("#" , $items[3]);
		
 	print STDERR "\nInputs:\n" if $debug;
	print STDERR map{"$_ => $input_hash{$_}" } keys %input_hash, "\n" if $debug;
	print STDERR "\nOutputs:\n" if $debug;
	print STDERR map{"$_ => $output_hash{$_}"} keys %output_hash, "\n" if $debug;

	# divide stoichiometry by greatest common denominator
	my ($min_stoich, @all_stoichs);
	
	map { push @all_stoichs, $input_hash{$_}; $min_stoich = $input_hash{$_} if $input_hash{$_} < $min_stoich || $min_stoich == 0 } keys %input_hash;
	map { push @all_stoichs, $output_hash{$_}; $min_stoich = $output_hash{$_} if $output_hash{$_} < $min_stoich || $min_stoich == 0 } keys %output_hash;
	
	my ($gcd, @gcd_candidates);
	
      outer: for ($gcd = $min_stoich; $gcd > 1; $gcd--)
      {
	  foreach my $stoich (@all_stoichs)
	  {
	      next outer if $stoich % $gcd != 0;
	  }
  
	  last; # found a gcd
      }
	
	map { $input_hash{$_} /= $gcd } keys %input_hash;
	map { $output_hash{$_} /= $gcd } keys %output_hash;
	
	print STDERR "Making directory $base_dir" if $debug;

	
	mkdir "$base_dir/$path_name";
	&write_fluxanalyzer_files("$base_dir/$path_name",\%input_hash,\%output_hash, 
				  \@tempArray,\%all_reactions,
				  \%reactions_to_substrate_arrays,
				  \%reactions_to_product_arrays,
				  {});
    }
}

sub show_path_results
{
    my ($data_results,$html,$cgi) = @_;
    
    print STDERR "Starting Results Display\n";
 
   
    #Display infinite loop warning if the indicator is on
    if($data_results->{"infinite"})
    {
 	push(@$html, "<h3>Warning:  Looks like an infinite loop</h3>");
    }
    
    #Delete the infinite loop indicator, so we don't need to have a if statment to check for it
    delete $data_results->{"infinite"};
    
    my $path_counter = 1;
    my $reactionPath;
    foreach my $token_id (sort { $data_results->{$a}->[1] <=> $data_results->{$b}->[1] }keys %$data_results)
    {
	print STDERR "Token id : $token_id\n";

	if(!($token_id =~ /^\d/))
	{
	    next;
	}

	my @path = @{$data_results->{$token_id}->[2]};
	my @compounds = @{$data_results->{$token_id}->[3]};
	
	push(@$html, "<pre>Token: $token_id\tInitial Pass: $data_results->{$token_id}->[0]\tDone:$data_results->{$token_id}->[1]\n\tReactions: @path\n\tVisted Compounds: @compounds\n\tPath Inputs\n@{$data_results->{$token_id}->[4]}\n\tPath Outputs\n@{$data_results->{$token_id}->[5]}\n\tBorrowed\n@{$data_results->{$token_id}->[6]}\n</pre>");

	
	if ($data_results->{$token_id}->[1] == 1)
	{
	    my $path_name = "path_".$path_counter++;
	    my @tempArray;
	    foreach my $elem(@{$data_results->{$token_id}->[2]}){
		push @tempArray, $elem;
	    }
	    $reactionPath = join "#",  @tempArray;
	    my $checkbox=$cgi->checkbox(-name=>"$path_name", -label=>'', 
					-value=>"$path_name;@{[ join ',', map { $_ . ':' . $data_results->{$token_id}->[7]->{$_} . ':' . $all_compounds_to_main{$_} }  keys %{$data_results->{$token_id}->[7]} ]};@{[ join ',', map { $_ . ':' . $data_results->{$token_id}->[8]->{$_} . ':' . $all_compounds_to_main{$_} } keys %{$data_results->{$token_id}->[8]} ]};$reactionPath");
	    push @$html, $checkbox, "&nbsp;$path_name", $cgi->br;
	}
	
	push @$html, "<hr>";
	
	
    }

  
    push @$html, $cgi->hidden(-name=>'reaction_info',
			      -value=>$reactionPath);

}

sub show_path_results_two
{
    my ($data_results,$html,$cgi) = @_;
    
    print STDERR "Starting Results Display\n";
 
   
    #Display infinite loop warning if the indicator is on
    if($data_results->{"infinite"})
    {
 	push(@$html, "<h3>Warning:  Looks like an infinite loop</h3>");
    }
    
    #Delete the infinite loop indicator, so we don't need to have a if statment to check for it
    delete $data_results->{"infinite"};
    
    my $path_counter = 1;
    my $reactionPath;
    foreach my $token_id (sort { $data_results->{$a}->[1] <=> $data_results->{$b}->[1] }keys %$data_results)
    {
	print STDERR "Token id : $token_id\n";

	if(!($token_id =~ /^\d/))
	{
	    next;
	}

	my @path = @{$data_results->{$token_id}->[2]};
	my @compounds = @{$data_results->{$token_id}->[3]};
	
	push(@$html, "<pre>Token: $token_id\tInitial Pass: $data_results->{$token_id}->[0]\tDone:$data_results->{$token_id}->[1]\n\tReactions: @path\n\tVisted Compounds: @compounds\n\tPath Inputs\n@{$data_results->{$token_id}->[4]}\n\tPath Outputs\n@{$data_results->{$token_id}->[5]}\n\tBorrowed\n@{$data_results->{$token_id}->[6]}\n</pre>");

	
	if ($data_results->{$token_id}->[1] == 1)
	{
	    my $path_name = "path_".$path_counter++;
	    my @tempArray;
	    foreach my $elem(@{$data_results->{$token_id}->[2]}){
		push @tempArray, $elem;
	    }
	    $reactionPath = join "#", @tempArray;
	    my $checkbox=$cgi->checkbox(-name=>"checked", -label=>'', 
					-value=>"$path_name;@{[ join ',', map { $_ . ':' . $data_results->{$token_id}->[7]->{$_} . ':' . $all_compounds_to_main{$_} }  keys %{$data_results->{$token_id}->[7]} ]};@{[ join ',', map { $_ . ':' . $data_results->{$token_id}->[8]->{$_} . ':' . $all_compounds_to_main{$_} } keys %{$data_results->{$token_id}->[8]} ]};$reactionPath");
	    push @$html, $checkbox, "&nbsp;$path_name", $cgi->br;
	}
	
	push @$html, "<hr>";
    }

}


sub set_loop_max
{
    my ($number) = @_;
    $loop_max = $number;
}

sub set_loop_max_assembly
{
    my ($number) = @_;
    $loop_max_assembly = $number;
}

sub analyze_scenario_connections
{
    my ($genome_id) = @_;
    my $scenario_dir = get_scenario_directory($genome_id);
    my $pathinfo_dir = $scenario_dir  . "/PathInfo";
    my @paths = `find $pathinfo_dir -type d -name "path_*"`;
    my %inputs;
    my %outputs;

    foreach my $dir (@paths)
    {
	chomp $dir;
	my ($cat, $subsys, $scenario);

	if ($dir =~ (/\/PathInfo\/(.*)\/(.*)\/(.*)\//)){
	    $cat = $1;
	    $subsys = $2;
	    $scenario = $3;
	} 
	else
	{
	    next;
	}
    
	open(M_INPUTS,$dir."/inputs_main") or die("Failed to open $dir/inputs_main");

	while (<M_INPUTS>)
	{   
	    chomp;
	    my ($cpd, $main) = split "\t" , $_;
	    my $info = join "\t", $subsys, $scenario;
	    $inputs{$cpd}->{$info} = 1 if $main eq "1";
	}
	close M_INPUTS;

	open(M_OUTPUTS, $dir."/outputs_main");

	while (<M_OUTPUTS>)
	{
	    chomp;
	    my ($cpd, $main) = split "\t" , $_;
	    my $info = join "\t", $subsys, $scenario;
	    $outputs{$cpd}->{$info} = 1 if $main eq "1";
	}
	close M_OUTPUTS;
    }

    my $analysis_dir = $scenario_dir . "/Analysis";
    mkdir $analysis_dir;

    open (IN_CONN, ">$analysis_dir/inputs_to_scenarios");
    foreach my $cpd (sort keys %inputs)
    {
	map { print IN_CONN "$cpd\t$_\n"; } keys %{$inputs{$cpd}};
    }
    close IN_CONN;

    open (OUT_CONN, ">$analysis_dir/outputs_to_scenarios");
    foreach my $cpd (sort keys %outputs)
    {
	map { print OUT_CONN "$cpd\t$_\n"; } keys %{$outputs{$cpd}};
    }
    close OUT_CONN;
}

sub predict_pegs_used
{
    my ($genome_id) = @_;
    my @scenarios = @{Scenario->get_genome_scenarios("All",1)};
    unless(scalar(@scenarios))
    {
	return undef;
    }
    my %reaction_to_pegs;
    
    my @ss_names;
    &load_superset_file;
    foreach (keys %superset_to_ss)
    {
	foreach my $subsys (@{$superset_to_ss{$_}})
	{
	    push @ss_names, $subsys;
	}
    
    }
    foreach my $subsystem_name(@ss_names)
    {
	my $subsystem = $fig->get_subsystem($subsystem_name);
	next if(!defined $subsystem);
	my $reactions_for_ss = get_reactions_for_genome_in_subsystem($genome_id,$subsystem_name);
	next if(! keys %$reactions_for_ss);
	foreach my $reaction (keys %$reactions_for_ss)
	{
	    if(defined $reaction_to_pegs{$reaction})
	    {
		push @{$reaction_to_pegs{$reaction}} , @{$reactions_for_ss->{$reaction}};
	    }	
	    else
	    {
		$reaction_to_pegs{$reaction} = $reactions_for_ss->{$reaction};
	    }
	}
    }

    my %peg_to_scenario;
    foreach my $scenario (@scenarios)
    {
	my @scenario_reactions = @{$scenario->get_reaction_ids};
	my $path_valid = 1;
	my %pegs;
	foreach my $reaction (@scenario_reactions)
	{
	    if(!defined $reaction_to_pegs{$reaction})
	    {
		$path_valid = 0;
		last;
	    }
	    else
	    {
		map {$pegs{$_} = 1} @{$reaction_to_pegs{$reaction}};
	    }
	}
	if($path_valid)
	{
	    foreach my $peg (keys %pegs)
	    {
		if(!defined $peg_to_scenario{$peg})
		{
		    $peg_to_scenario{$peg} = [$scenario->get_id()];
		}
		else
		{
		    push @{$peg_to_scenario{$peg}} , $scenario->get_id();
		}
	    }
	}
    }
    
    return \%peg_to_scenario;
}

sub analyze_network_for_genome
{
    my ($genome_id, $transports, $biomass) = @_;

    my ($node_type, $node_connections, %compound_ids, $path_info, $reaction_presence);
    my $unreachable = 0;
    chomp $genome_id;

    Scenario->set_fig($fig);
    print STDERR "Running check_model_hope on genome $genome_id\n" if($debug);
    print STDERR "\tStage 1 - Loading Scenarios, transports and biomass\n" if($debug);
    my @scenarios = @{Scenario->get_genome_scenarios($genome_id,0)};

    my %scenarios_used;
    my %scenarios_id;

    foreach (@scenarios)
    {
	print STDERR "\t\tLoading ", $_->get_scenario_name(), "\n" if $debug;
	$scenarios_used{$_->get_scenario_name()} = 0;
	$scenarios_id{$_->get_id()} = 0;
    }

    my %biomass_lib_cpds;
    foreach my $biomass_cpd (@$biomass)
    {
	$biomass_lib_cpds{$biomass_cpd} = 1;
    }

    my %transportable_cpds = map { $_ => 1 } @$transports;
    
    print STDERR "\tStage 1 - Complete\n" if($debug);
    print STDERR "\tStage 2 - Run Analysis Tools\n" if($debug);
    
#   Total list of what compounds we reach.
    my %list_generated;
    
#   Here we are running all the scenarios we have inputs for
#   to analyze what biomass we can reach. We are also identifying 
#   any Scenarios that are not being utilized.
    
    # current_inputs contains the compounds that we have "reached" so far in 
    # traversing the network; we preload transported compounds, as well as
    # some cofactors that are "main" inputs to scenarios
    my %current_inputs;
    $current_inputs{"C00342"}->{"[]cofactor (thioredoxin)[C00342]"} = 1;


    foreach (keys %transportable_cpds)
    {
	# keep track of the paths that lead to each compound
	my @names = $fig->names_of_compound($_);
	my $name = @names > 0 ? $names[0] : $_;
	$current_inputs{$_}->{"[]transport of ${name}[$_]"} = 1;
    }
    
    my $iter_num = 1;

    while(1)
    {
	my ($temp_in,$break) = generate_new_inputs(\%current_inputs);
	
	foreach my $cpd (keys %$temp_in)
	{
	    map { $current_inputs{$cpd}->{$_} = 1 }  keys %{$temp_in->{$cpd}};
	}
	
	map { $list_generated{$_} = 1 } keys %$temp_in;
	
	if($break)
	{
	    last;
	}

	print STDERR "Rechecking organism scenarios with inputs\n" if($debug);
	$iter_num++;
    }
    print STDERR "Finished checking organism scenarios\n" if($debug);
    
    my %biomass_reached;
    
    foreach my $found (keys %list_generated)
    {
	if(defined $biomass_lib_cpds{$found})
	{
	    $biomass_reached{$found} = $current_inputs{$found};
	}
    }

    foreach (keys %transportable_cpds)
    {
	push (@{$node_type->{'transport'}}, $_);
	$compound_ids{$_} = 1;
    }
    
    print STDERR "\tStage 2 - Complete!\n" if($debug);
    print STDERR "\tStage 3 - Loading All Scenarios\n" if($debug);

    @scenarios = @{Scenario->get_genome_scenarios("All",0)};
    
    my %old_scenarios_used = %scenarios_used;
    my %old_scenarios_id = %scenarios_id;
    
    my %old_biomass_reached;
    foreach my $cid (keys %biomass_reached)
    {
	my %inner_hash_copy = %{$biomass_reached{$cid}};
	$old_biomass_reached{$cid} = \%inner_hash_copy
    }
    
    my %old_current_inputs;
    foreach my $cid (keys %current_inputs)
    {
	my %inner_hash_copy = %{$current_inputs{$cid}};
	$old_current_inputs{$cid} = \%inner_hash_copy;
    }
    
    
    %list_generated = ();
    
    foreach (@scenarios)
    {
	$scenarios_used{$_->get_scenario_name()} = 0 if ! exists $scenarios_used{$_->get_scenario_name()};
	$scenarios_id{$_->get_id()} = 0 if ! exists $scenarios_id{$_->get_id()};
    }
    
    while(1)
    {
	print STDERR "Rechecking all scenarios with inputs\n" if($debug);
	$iter_num++;

	my ($temp_in,$break) = generate_new_inputs(\%current_inputs);
	
	foreach my $cpd (keys %$temp_in)
	{
	    map { $current_inputs{$cpd}->{$_} = 1 } keys %{$temp_in->{$cpd}};
	}
	
	map { $list_generated{$_} = 1 } keys %$temp_in;
	
	if($break)
	{
	    last;
	}
    }
    print STDERR "Finished checking all scenarios\n" if($debug);
    
    %biomass_reached = ();
    
    foreach my $found (keys %list_generated)
    {
	if(defined $biomass_lib_cpds{$found})
	{
	    $biomass_reached{$found} = $current_inputs{$found};
	}
    }
    
    my (%cpds_to_print);
    my %org_cpds_to_check;
    map { $org_cpds_to_check{$_} = 1 } keys %old_biomass_reached;
    foreach my $cpd (sort keys %biomass_reached)
    {
	next if exists $old_biomass_reached{$cpd};
	
	$cpds_to_print{$cpd} = 1;
	my $done = 0;
	
	my @cpd_list = ($cpd);
	while (! $done)
	{
	    my @new_substrates = collect_substrates(\%current_inputs, @cpd_list);
	    @cpd_list = ();
	    foreach my $substrate (@new_substrates)
	    {
	    if (exists $old_current_inputs{$substrate})
	    {
		$org_cpds_to_check{$substrate} = 1;
		next;
	    }
	    push (@cpd_list, $substrate) if ! exists $cpds_to_print{$substrate};
	    $cpds_to_print{$substrate} = 1;
	    }
	    $done = 1 if scalar @new_substrates == 0;
	}
    }

    foreach my $cpd (keys %cpds_to_print)
    {
	foreach my $path (keys %{$current_inputs{$cpd}})
	{
	    my ($substrates,$scenario,$products) = $path =~ /\[(.*)\](.*)\[(.*)\]/;
	    
	    foreach my $substrate (split ",", $substrates)
	    { 
		foreach my $product (split ",", $products)
		{
		    print STDERR "cpd: $cpd  scenario: $scenario  path: $path  substrate: $substrate  product: $product\n" if ($debug);
		    $node_connections->{$substrate}->{$scenario} = 1;
		    $compound_ids{$substrate} = 1;
		    $node_connections->{$scenario}->{$product} = 1;
		    $compound_ids{$product} = 1;
		}
	    }
	}
    }
    
    %cpds_to_print = ();
    
    foreach my $cpd (sort keys %org_cpds_to_check)
    {
	$cpds_to_print{$cpd} = 1;
	my $done = 0;
	
	my @cpd_list = ($cpd);
	while (! $done)
	{
	    my @new_substrates = collect_substrates(\%old_current_inputs, @cpd_list);
	    @cpd_list = ();
	    foreach my $substrate (@new_substrates)
	    {
		push (@cpd_list, $substrate) if ! exists $cpds_to_print{$substrate};
		$cpds_to_print{$substrate} = 1;
	    }
	    $done = 1 if scalar @new_substrates == 0;
	}
    }

    foreach my $cpd (keys %cpds_to_print)
    {
	foreach my $path (keys %{$old_current_inputs{$cpd}})
	{
	    my ($substrates,$scenario,$products) = $path =~ /\[(.*)\](.*)\[(.*)\]/;
	    
	    foreach my $substrate (split ",", $substrates)
	    { 
		foreach my $product (split ",", $products)
		{
		    $node_connections->{$substrate}->{$scenario} = 1;
		    $compound_ids{$substrate} = 1;
		    $node_connections->{$scenario}->{$product} = 1;
		    $compound_ids{$product} = 1;
		}
	    }
	}
    }

    foreach (sort keys %biomass_lib_cpds)
    {
	if ($old_biomass_reached{$_})
	{
	    push (@{$node_type->{'biomass reached'}}, $_);
	    $compound_ids{$_} = 1;
	}
	elsif ($biomass_reached{$_})
	{
	    push (@{$node_type->{'biomass reached all'}}, $_);
	    $compound_ids{$_} = 1;
	}
    }
    
    my %hope_ss_info_cache;

    foreach my $ss_and_scen (sort keys %scenarios_used)
    {
	if($scenarios_used{$ss_and_scen} == 1)
	{
	    print STDERR "$ss_and_scen was used\n";
	    my ($ss_name, $scen_name) = split "/", $ss_and_scen;

	    if (! exists $hope_ss_info_cache{$ss_name}) {
		$hope_ss_info_cache{$ss_name}->{"kegg_info"} = $fig->get_scenario_info($scen_name);
		$hope_ss_info_cache{$ss_name}->{"org_reactions"} = get_reactions_for_genome_in_subsystem_fast($genome_id, $ss_name);
	    }

	    my $org_ss_reactions = $hope_ss_info_cache{$ss_name}->{"org_reactions"};
	    my @scenarios; # path info for the scenarios

	    if (exists ($old_scenarios_used{$ss_and_scen}))
	    {
		push (@{$node_type->{'scenario present'}}, $ss_and_scen);
		@scenarios = @{Scenario->get_genome_scenarios_for_scenario($genome_id, $ss_name, $scen_name, 0)};
	    }
	    else
	    {
		# determine how many reactions the organism is missing
		@scenarios = @{Scenario->get_genome_scenarios_for_scenario("All", $ss_name, $scen_name, 0)};
		my @percentiles = ("0", "25", "50", "75");
		my $percentile = 3; # assume the greatest number of missing reactions

		foreach my $scenario (@scenarios)
		{
		    my @reactions = @{$scenario->get_path_info()};
		    map {s/_\w//} @reactions;
		    my $num_missing_rxns = 0;
		    map {$num_missing_rxns++ unless exists $org_ss_reactions->{$_}} @reactions;
		    my $sc_percentile = int( ($num_missing_rxns / scalar @reactions) * 4 );
		    $percentile = $sc_percentile < $percentile ? $sc_percentile : $percentile;
		}

		push (@{$node_type->{"scenario present all $percentiles[$percentile]"}}, $ss_and_scen);
	    }

	    # get the path info for scenarios
	    my @scenario_path_info;

	    foreach my $scenario (@scenarios)
	    {
		my $map_ids = $hope_ss_info_cache{$ss_name}->{"kegg_info"}->{"scenarios"}->{$scen_name}->{"map_ids"};
		my $reactions = $scenario->get_path_info();
		my @rxns_and_cpds;
		
		foreach my $reaction (@$reactions)
		{
		    my ($rxn_id, $dir) = split "_", $reaction;
		    my ($subside, $prodside) = $dir eq "R" ? (0,1) : (1,0);
		    my (@main_subs, @main_prods, @reaction_info);
		    my @subs_info = $fig->reaction2comp($rxn_id,$subside,$map_ids);
		    map { if ($_->[2]) { push @main_subs, $_->[0]; $compound_ids{$_->[0]} = 1 } } @subs_info;
		    if (@main_subs == 0) {
			# no main info
			map { push @main_subs, $_->[0]; $compound_ids{$_->[0]} = 1 } @subs_info;
		    }
		    my @prods_info = $fig->reaction2comp($rxn_id,$prodside,$map_ids);
		    map { if ($_->[2]) { push @main_prods, $_->[0]; $compound_ids{$_->[0]} = 1 } } @prods_info;
		    if (@main_prods == 0) {
			# no main info
			map { push @main_prods, $_->[0]; $compound_ids{$_->[0]} = 1 } @prods_info;
		    }

		    push @reaction_info, @main_subs, $rxn_id, @main_prods;
		    push @rxns_and_cpds, (join " ", @reaction_info);
		    $reaction_presence->{$rxn_id} = exists $org_ss_reactions->{$rxn_id} ? "reaction present" : "reaction absent";
		}
		
		push @scenario_path_info, (join ";", @rxns_and_cpds);
	    }

	    $path_info->{$ss_name."/".$scen_name} = join "&", @scenario_path_info;
	}
    }
    
    foreach my $cpd (keys %current_inputs)
    {
	# compounds present with All scenarios that aren't biomass
	if (! exists($old_current_inputs{$cpd}) && ! exists ($biomass_reached{$cpd}))
	{
	    push (@{$node_type->{'compound present all'}}, $cpd);
	    $compound_ids{$cpd} = 1;
	}
    }

    print STDERR "\tStage 3 - Complete!\n" if($debug);
    
    return ($node_type, $node_connections, \%compound_ids, $path_info, $reaction_presence);

    sub generate_new_inputs
    {
	my ($inputs) = @_;
	my %new_inputs;
	my $break = 1;
	
	my %used_scenario_paths_this_time;
	
	foreach my $scenario (@scenarios)
	{
	    if($scenarios_id{$scenario->get_id()} == 1) # || $scenarios_used{$scenario->get_scenario_name()} != 0)
	    {
		next;
	    }
	    
	    my $sids = join ",", sort keys %{$scenario->get_substrates};
	    my $pids = join ",", sort keys %{$scenario->get_products};
	    my $cid_scenario_name = "[".$sids."]".$scenario->get_scenario_name()."[".$pids."]";
	    
	    my @matched;
	    
	    foreach my $substrate (keys %{$scenario->get_substrates()})
	    {
		if(exists $inputs->{$substrate})
		{
		    push @matched, $substrate;
		}
	    }
	    
	    if (scalar @matched != scalar(keys %{$scenario->get_substrates()}))
	    {
		print STDERR "($iter_num) Scenario ".$cid_scenario_name." did not run.  matched: @matched\n" if($debug);
		next;
	    }
	    else
	    {
		print STDERR "($iter_num) Scenario ".$cid_scenario_name." ran. matched: @matched\n" if($debug);
		$break = 0;   
		$scenarios_id{$scenario->get_id()} = 1;
		$scenarios_used{$scenario->get_scenario_name()} = 1;
		
		# local determination that we've used an equivalent path already in this
		# time through the subroutine, so no need to process further
		next if (exists $used_scenario_paths_this_time{$cid_scenario_name});
		$used_scenario_paths_this_time{$cid_scenario_name} = 1;		 
		
		foreach my $product (keys %{$scenario->get_products()})
		{
		    $new_inputs{$product}->{$cid_scenario_name} = 1;
		}
	    }
	}
	return (\%new_inputs,$break);
    }

    sub collect_substrates
    {
	my $input_hash = shift @_;
	my @cpds = @_;
	my %substrate_hash;
	
	foreach my $cpd (@cpds)
	{
	    foreach my $path (keys %{$input_hash->{$cpd}})
	    {
		my ($substrates,$scenario,$products) = $path =~ /\[(.*)\](.*)\[(.*)\]/;
		
		foreach my $substrate (split ",", $substrates)
		{ 
		    $substrate_hash{$substrate} = 1;
		}
	    }
	}
	return keys %substrate_hash;
	
	
    }
}

sub get_reactions_for_genome_in_subsystem
{
    my ($genome_id, $ss_name) = @_;
    my $subsystem = $fig->get_subsystem($ss_name);
    my %genome_reactions = ();
    my %hope_reactions = $subsystem->get_hope_reactions;

    foreach my $role (keys %hope_reactions)
    {
	my @peg_list = $fig->seqs_with_role($role, "", $genome_id);

	if (scalar @peg_list > 0) 
	{
	    foreach my $reaction (@{$hope_reactions{$role}}) 
	    {
		push @{$genome_reactions{$reaction}}, @peg_list;
	    }
	}	
    }

    return \%genome_reactions;
}

sub get_reactions_for_genome_in_subsystem_fast
{
	my ($genome_id, $ss_name) = @_;
	my $ss_dir = Subsystem::get_dir_from_name($ss_name);
	my (%hope_reactions, %genome_reactions);
	
	# this is taken directly from FIGRules::GetHopeReactions with a few modifications
	# this is called when a subsystem object is constructed normally,
	# but here all we need are the reactions so this is quite a bit faster.
	# returns the same reactions as get_reactions_for_genome_in_subsystem
	my $hopeFileName = "$ss_dir/hope_reactions";
	if (-f $hopeFileName)
	{
		# Open the hope reaction file.
		open(RXNFILE, "<$hopeFileName");
		# Loop through the hope reaction file.
		while (<RXNFILE>)
		{
			# Parse out the role and the reaction list.
			if ($_ =~ /^(\S.*\S)\t(\S+)/)
			{
				my ($role, $reactions) = ($1, $2);
				# FIGRules checks to see if the role is in this subsystem but we can't. Watch for this!
				push( @{$hope_reactions{$role}}, split(/,\s*/,$reactions) );
			}
		}
		close(RXNFILE);
	}
	
	foreach my $role (keys %hope_reactions)
	{
	    my @peg_list = $fig->seqs_with_role($role, "", $genome_id);

	    if (scalar @peg_list > 0) 
	    {
		foreach my $reaction (@{$hope_reactions{$role}}) 
		{
		    push @{$genome_reactions{$reaction}}, @peg_list;
		}
	    }	
	}
	
	return \%genome_reactions;
}

sub get_scenario_directory
{
    my ($genome) = @_;
    # want to match "All" and "All.bk", etc.
    return $genome =~ "^All" ? $FIG_Config::global."/Models/$genome/Scenarios" : $fig->scenario_directory($genome);

}

# return three values:
# 1) reference to a hash of invalid roles, i.e., roles no longer in the subsystem.
# 2) reference to a hash of invalid reactions, i.e., reactions no longer in KEGG.
# 3) flag set to main classification value if the subsystem classification has changed 
#    to Experimental or none.
sub check_hope_info
{
    my ($ss) = @_;
    my (%invalid_roles, %invalid_reactions, $invalid_classification, %invalid_kegg_maps, %reaction_ec_mismatch);
    my $sub = $fig->get_subsystem($ss);

    if (! defined $sub)
    {
	return ( {}, {}, "Subsystem no longer exists", {}, {} );
    }

    my $classification = $sub->get_classification->[0];
    $invalid_classification = $classification if ($classification eq '' || $classification =~ /Experimental/);
    my %hope_reactions = $sub->get_hope_reactions;
    my @roles = $sub->get_roles;
    my %roles;
    my %reactions;
    
    map { $roles{$_} = 1 } @roles;

    foreach my $role (keys %hope_reactions)
    {
	if (! exists($roles{$role}))
	{
	    $invalid_roles{$role} = $hope_reactions{$role};
	}
	
	map {push @{$reactions{$_}}, $role} @{$hope_reactions{$role}};
    }

    foreach my $name ($sub->get_hope_scenario_names)
    {
	my @additional_reactions = $sub->get_hope_additional_reactions($name);
	map { push @{$reactions{$_}}, "ADDITIONAL_REACTIONS" } @additional_reactions;

	foreach my $id ($sub->get_hope_map_ids($name))
	{
	    my $map_name = $fig->map_name("map".$id);
	    $invalid_kegg_maps{$id} = 1 if $map_name eq "";
	}
    }


    foreach my $reaction (keys %reactions)
    {
	unless ($reaction eq "not_in_KEGG" || $fig->valid_reaction_id($reaction))
	{
	    $invalid_reactions{$reaction} = $reactions{$reaction};
	}
	else
	{
	    my %kegg_ecs;
	    @kegg_ecs{$fig->catalyzed_by($reaction)} = ();
	    foreach my $role ( @{$reactions{$reaction}})
	    {
		if ($role =~ /(\d+\.\d+\.\d+\.\d+)/ && ! exists($kegg_ecs{$1})) {
		    push @{$reaction_ec_mismatch{$reaction}}, $role;
		}
	    }
	}
    }

    return (\%invalid_roles, \%invalid_reactions, $invalid_classification, \%invalid_kegg_maps, \%reaction_ec_mismatch);

}

1;
