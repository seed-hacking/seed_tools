# -*- perl -*-
#
# Copyright (c) 2003-2011 University of Chicago and Fellowship
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
#!/usr/bin/perl -w

=head1 Make Indexes

This script builds the Glimpse indexes for the SEED. Two indexes are created--
one for features and one for roles. For each index, we create a sequential
file and then feed those files to the Glimpse indexing utility.

The currently-supported command-line options are as follows.

=over 4

=item user

Name suffix to be used for log files. If omitted, the PID is used.

=item trace

Numeric trace level. A higher trace level causes more messages to appear. The
default trace level is 2. Tracing will be directly to the standard output
as well as to a C<trace>I<User>C<.log> file in the FIG temporary directory,
where I<User> is the value of the B<user> option above.

=item sql

If specified, turns on tracing of SQL activity.

=item background

Save the standard and error output to files. The files will be created
in the FIG temporary directory and will be named C<err>I<User>C<.log> and
C<out>I<User>C<.log>, respectively, where I<User> is the value of the
B<user> option above.

=item h

Display this command's parameters and options.

=item phone

Phone number to message when the script is complete.

=item dir

Directory in which to put the indexes. The default is C<Indexes> in the Data directory.

=back

=cut

use strict;
use FIG;
use Tracer;
my $fig = new FIG;
# Get the command-line options and parameters.
my ($options, @parameters) = StandardSetup([qw(FIG CustomAttributes) ],
                                           {
                                              dir => ["$FIG_Config::data/Indexes", "directory in which to put the indexes"],
                                              phone => ["", "phone number (international format) to call when load finishes"],
                                           },
                                           "",
                                           @ARGV);
# Set a variable to contain return type information.
my $rtype;
# Insure we catch errors.
eval {
    # Insure we can get to Glimpse.
    FIG::augment_path($FIG_Config::bin);
    # Compute the directory name.
    my $dirName = $options->{dir};
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #...Generate aliases to be added
    #-----------------------------------------------------------------------
    my $aliasFile = "$FIG_Config::temp/aliases$$";
    Open(\*SYN, "<$FIG_Config::global/peg.synonyms");

    # Open(\*ALIASES, "| sort +0n -1 +1 -2 -T $FIG_Config::temp > $aliasFile");
    # sort: invalid option -- 1
    # Try `sort --help' for more information.
    Open(\*ALIASES, "| sort -k 1,1n -k 2,2n -T $FIG_Config::temp > $aliasFile");

    Trace("Loading synonyms.") if T(2);
    while (defined($_ = <SYN>)) {
        # Each line represents a list of synonymous features. We separate the list
        # into FIG features and non-FIG features.
        chomp;
        my($x,$y) = split(/\t/,$_);
        my @ids = map { $_ =~ /^([^,]+),/; $1 } ($x,split(/;/,$y));
        my @fig = ();
        my(@nonfig) = ();
        foreach $_ (@ids) {
            if ($_ =~ /^fig\|/) {
                push(@fig,$_);
            } else {
                push(@nonfig,$_);
            }
        }
        # If we have both FIG and non-FIG features, we establish all the non-FIG
        # synonyms as aliases for the FIG features. Note that only PEGs are
        # used in this process, and we save some space by storing only the
        # organism ID and the PEG number in the alias file.
        if ((@fig > 0) && (@nonfig > 0)) {
            my $non = join(" ",@nonfig);
            foreach my $peg (@fig) {
                if ($peg =~ /^fig\|(\d+\.\d+)\.peg\.(\d+)/) {
                    print ALIASES "$1\t$2\t$non\n";
                }
            }
        }
    }
    close(SYN);
    close(ALIASES);

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #...Remove old indices and re-create directory structure...
    #----------------------------------------------------------------------- 
    Trace("Removing old indices...") if T(2);
    &FIG::run("rm -rf $dirName/*index");
    
    Trace("Creating new index files.") if T(2);
    (-d "$dirName")
        || mkdir("$dirName",0777)
        || Confess("Could not make $dirName");

    my $peg_index_idx = "00";
    my $max_peg_index_size = 1_500_000_000;

    Open(\*PEG, ">$dirName/peg.index.$peg_index_idx");
    $peg_index_idx++;
    
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #...Load organism features and assignments...
    #-----------------------------------------------------------------------
    # Prime the process by getting the first line in the alias file we
    # created. The alias file is sorted by organism ID, just like the
    # main loop.
    Open(\*ALIASES, "<$aliasFile");
    my $syn_based = <ALIASES>;
    Trace("Main organism loop starting.") if T(2);
    # Loop through the organisms in order.
    foreach my $org (sort {$a <=> $b} $fig->genomes) {
        Trace("Loading features and assignments for $org.") if T(3);
        Trace("Processing synonyms and aliases.") if T(3);
        # This hash will contain a list of all the non-FIG synonyms for each PEG
        # in this organism that has them. The key will be the PEG number.
        my %syn_aliases;
        while ($syn_based && ($syn_based =~ /^(\S+)\t(\S+)\t(\S.*\S)/) && ($1 == $org)) {
            $syn_aliases{$2} = $3;
            $syn_based = <ALIASES>;
        }
        # Now we start collecting data on the PEGs. This hash will be keys on
        # FIG ID, and each entry will contain data we want to send to Glimpse.
        my %by_peg = ();
        # First, we process the feature tables. The feature's location and all its aliases
        # are crammed into %by_peg. Note that because only PEGs have synonym aliases, we do
        # not pass the synonym data into the loads for PI and PP features.
        Trace("Loading feature tables.") if T(3);
        &load_tbl("$FIG_Config::organisms/$org/Features/peg/tbl",\%syn_aliases, \%by_peg);
        &load_tbl("$FIG_Config::organisms/$org/Features/pp/tbl",{ }, \%by_peg);
        &load_tbl("$FIG_Config::organisms/$org/Features/pi/tbl",{ }, \%by_peg);
        &load_tbl("$FIG_Config::organisms/$org/Features/rna/tbl",{ }, \%by_peg);
        # Now %by_peg has a list of aliases and stuff for each feature possessed by this
        # organism. Next, we add the assignments.
        Trace("Loading assignments.") if T(3);
        # First we load the master assignments.
        &load_assign("$FIG_Config::organisms/$org/assigned_functions","master",\%by_peg);
        # Now we load the assignments made by any individual users. These are stored in
        # the UserModels subdirectory. Most organisms won't have this directory, so we
        # start with an IF.
        if (opendir(TMP, "$FIG_Config::organisms/$org/UserModels")) {
            my @users = grep { $_ =~ /^[A-Za-z]/ } readdir(TMP);
            closedir(TMP);
            foreach my $user (@users) {
                &load_assign("$FIG_Config::organisms/$org/UserModels/$user/assigned_functions",$user);
            }
        }
        # Now we get all the attributes for the organism's feature and put them into a hash by
        # feature ID.
        Trace("Building attribute cache.") if T(3);
        #my @attributeList = $fig->get_attributes("fig|$org.%");
	# Disable attribute indexing.
	my @attributeList = ();
        my $count = @attributeList;
        # Pour each attribute's data into the %by_peg has.
        Trace("Processing attributes. $count found for $org.") if T(3);
        for my $tuple (@attributeList) {
            # Get this attribute's feature ID.
            my $fid = shift @{$tuple};
            # Only proceed if the ID is in the PEG hash. If it isn't, then it's probably a feature
            # from a different SEED.
            if (exists $by_peg{$fid}) {
                # Push the attribute's data into the peg's data list. We push in the whole
                # tuple, because the SHIFT operation above removed the feature ID from it.
                # First, however, we need to make a string out of it. It starts with the key.
                my $aKey = shift @$tuple;
                # Convert the key pieces to words.
                $aKey =~ s/_/ /g;
                # Join the key with its values.
                my $aString = "attribute: $aKey = " . join(" ", @{$tuple});
                # Push the result into the hash.
                push @{$by_peg{$fid}}, $aString;
            } else {
                Trace("No entry found for attribute peg $fid.") if T(3);
            }
        }
        # Clean the attribute list from memory. It can get pretty big.
        undef @attributeList;
        # Now we're about write the peg data to the output file. Get the taxon ID without it's
        # strain number and the full genus-species name.
        my $tax = $org;
        $tax =~ s/\.\d+$//;
        my $full = $fig->genus_species($org);
        Trace("Writing peg index data.") if T(3);

        foreach my $peg (sort { &FIG::by_fig_id($a,$b) } keys(%by_peg)) {
            Trace("Processing PEG $peg.") if T(4);
            my $x = $by_peg{$peg};
            print PEG "$peg\t$full $tax\t",join("\t",@$x),"\n";
	    if (tell(PEG) > $max_peg_index_size)
	    {
		close(PEG);
		Open(\*PEG, ">$dirName/peg.index.$peg_index_idx");
		$peg_index_idx++;
	    }
        }
        Trace("Organism $org complete.") if T(3);
    }
    Trace("PEG index cleanup.") if T(2);
    # Close the PEG file. It's now ready to load into Glimpse.
    close(PEG);
    # Clean up the alias file we built.
    close(ALIASES);
    unlink($aliasFile);
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # ...Load EC numbers and roles...
    #
    # Many things were broken here by changes in configuration to some
    # SEEDs, but not others.
    #----------------------------------------------------------------------- 
    Trace("Loading EC numbers and roles.") if T(2);
    my @roles;
    # Loop through the enzyme file. Note that enzyme data is multi-line, with each
    # logical record separated by a line with two slashes on it.
    # my $kegg = $FIG_Config::kegg || "$FIG_Config::data/KEGG";
    my $kegg = $FIG_Config::kegg || "$FIG_Config::data/KEGG";
    my ( $enzyme ) = grep { -f $_ }
                     map  { "$kegg/$_" }
                     qw( ligand/enzyme/enzyme enzyme );
    if ( $enzyme && Open(\*ENZ, "<$enzyme" ) )
    {
        my $old_EOR = $/;
        $/ = "\n///\n";
        # Loop through the KEGG file entries.
        while (defined($_ = <ENZ>)) {
            # Check for an entry line followed by a name line and parse out the pieces.
            # $1 will be an EC number, $2 will be a list of names delimited by new-lines.
            if ($_ =~ /ENTRY\s+EC\s+(\d+\.\d+\.\d+\.\d+).*\nNAME((\s+\S[^\n]+\n)+)/s) {
                my $ec = $1;
                my $names = $2;
                # Now loop through names in the record. For each name, we trim off the
                # excess spaces and output a statement that the EC number has that name.
                foreach my $name (map { $_ =~ /^\s+(\S.*\S)/; $1 } split(/\n/,$names)) {
                    push(@roles,"$ec - $name");
                }
            }
        }
        close(ENZ);
        $/ = $old_EOR;
        Trace("EC numbers and roles loaded.") if T(2);
    } else {
        Trace("EC number and role skipped.") if T(2);
    }

    #open(REAC,"<$FIG_Config::data/KEGG/reaction")
    #    || die "could not open $FIG_Config::data/KEGG/reaction";
    #$/ = "\n///\n";
    #while (defined($_ = <REAC>))
    #{
    #    if ($_ =~ /NAME\s+(\S[^\n]+\S).*ENZYME\s+(\d+\.\d+\.\d+\.\d+)/s)
    #    {
    #   push(@roles,"$2 - $1");
    #    }
    #}
    #close(REAC);

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #...Load functional neighborhoods...
    #
    # Many things were broken here by changes in configuration to some
    # SEEDs, but not others.
    #-----------------------------------------------------------------------
    Trace("Loading functional neighborhoods.") if T(2);
    # The neighborhood file, like the KEGG file, has multi-line records
    # separated by double slashes.
    my $neighbor = "FIG_Config::global/role.neighborhoods";
    if ( -f $neighbor && Open( \*NEIGH, "<$neighbor" ) )
    {
        my $old_EOR = $/;
        $/ = "\n//\n";
        # Each functional neighborhood that is not an EC number is added to the role list.
        while (defined($_ = <NEIGH>)) {
            if ($_ =~ /^(\S[^\n]+\S)\n\n/s) {
                my $role = $1;
                if ($role !~ /^\d+\.\d+\.\d+\.\d+$/) {
                    push(@roles,$role);
                }
            }
        }
        $/ = $old_EOR;
        close(NEIGH);
        Trace("Neighborhoods loaded.") if T(2);
    }

    # Now we write out the roles.
    # Open the output file for the role index.
    Open( \*ROLE, ">$dirName/role.index" );
    foreach ( sort { lc $a cmp lc $b } @roles ) { print ROLE "$_\n" }
    close(ROLE);
    Trace("Role index created.") if T(2);
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #...Run GLIMPSE indexing utility...
    #-----------------------------------------------------------------------
    Trace("Running GLIMPSE indexing utility --- this may take some time...") if T(2);
    
    Open(\*EXCLUDE, ">$dirName/.glimpse_exclude");
    print EXCLUDE "$dirName/.glimpse*\n";
    close(EXCLUDE);
    
    &FIG::run("glimpseindex -o -n 100 -E -H $dirName -M 40 $dirName");
    &FIG::run("chmod -R 777 $dirName");
    
    Trace("GLIMPSE index completed.") if T(2);
    Trace("Index Make complete.") if T(2);
};
if ($@) {
    Trace("Script failed with error: $@") if T(0);
    $rtype = "error";
} else {
    Trace("Script complete.") if T(2);
    $rtype = "no error";
}
if ($options->{phone}) {
    my $msgID = Tracer::SendSMS($options->{phone}, "Make Indexes terminated with $rtype.");
    if ($msgID) {
        Trace("Phone message sent with ID $msgID.") if T(2);
    } else {
        Trace("Phone message not sent.") if T(2);
    }
}

## Read all of a feature's location and alias data into the $by_peg hash. The
## features are read from $file, and additional aliases may be in the $sin hash.
sub load_tbl { 
    my($file,$syn,$by_peg) = @_;
    my($id,$alias,@aliases);

    if (open(TMP,"<$file")) {
        while (defined($_ = <TMP>)) {
            chomp;
            ($id,undef,@aliases) = split(/\t/,$_);
            my $aliases = join(" ",@aliases);
            $id =~ /(\d+)$/;
            if ($_ = $syn->{$1}) {
                $aliases .= " " . $_;
            }
            push(@{$by_peg->{$id}},"aliases:" . $aliases);
        }
        close(TMP);
    }
}

# Read the assignments for the user $user into $by_peg from the file $file.
sub load_assign {
    my($file,$user,$by_peg) = @_;
    my $old_eol;
    ($old_eol, $/) = ($/, "\n");
    my %num_occurrences;
    my %tmp;

    if (open(TMP,"<$file")) {
        while (defined($_ = <TMP>)) {
            chomp;
            my ($id,$func) = split(/\t/,$_);
            
            if ((not &FIG::hypo($func)) && ($num_occurrences{$func} < 1000)) {
                $tmp{$id} = "function:$func\#$user";
            }
            ++$num_occurrences{$func};
        }
        close(TMP);
    
        foreach my $id (keys(%tmp)) {
            push(@{$by_peg->{$id}},$tmp{$id});
        }
    }
    
    $/ = $old_eol;
    return;
}
