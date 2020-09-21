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

package Subsystem;

use Carp;

use POSIX;
use DirHandle;
use Data::Dumper;
use File::Copy;
use File::Spec;
use IPC::Open2;
use FileHandle;
use Tracer;

use strict;

our $SS_version_checked = 0;

my $notes_separator = "###############################";
my @section_order = qw(description notes literature variants);
my %defined_sections = map { $_ => 1 } @section_order;

=head1 Subsystem Manipulation

Any manipulation of subsystem data should happen through this interface.
This allows us to assure ourselves that the relational tables that
mirror and index the subsystem data are kept up to date with the
canonical version of the subsystem information in the flat-files
kept in $FIG_Config::data/Subsystems.

=head2 Objects.

We define the following perl objects:

Subsystem: represents a subsystem. It can be read from disk and
written to disk, and manipulated via its methods when in memory.

If we were completely on the OO side of the world, we would also
define the following set of objects. However, we are not, so they are
only objects in a conceptual sense. They are implemented using the
basic perl datatypes.

Role: represents a single role. A role has a name and an abbreviation.

RoleSubset: represents a subset of available roles. A subset has a
name and a list of role names that comprise the subset.

=head2 Thoughts on locking

It is currently dangerous for multiple users to modify spreadsheets at once.
It will likely remain dangerous while the subsystem backend is fairly
stateless, as it is with the CGI mechanism.

We'd like to make this a little safer. One mechanism might be to allow
a user to open a subsystem for modification, and others for readonly access.
For this to work we have to be able to tell which users is allowed; current
implementation uses the curator of the subsystem for this purpose.

NB: This module does not currently attempt to handle locking or exclusion.
It is up to the caller (user application, CGI script, etc) to do so.
It does attempt to use locking internally where appropriate.

=head2 Data structures

We maintain the following data structures (all members of %$self).

=over 4

=item dir

Directory in which the subsystem is stored.

=item notes

The current notes contents for the subsystem

=item version

Current subsystem version.

=item exchangable

1 if subsystem is exchangable, 0 otherwise.

=item roles

List of role names.

=item role_index

hash that maps from role name to index

=item role_abbrs

list of role abbreviations

=item abbr

hash mapping from role abbreviation to role name

=item col_subsets

list of column subset names

=item col_subset_members

hash that maps from column subset name to subset members

=item col_active_subset

currently-active column subset

=item row_active_subset

currently-active row subset

=item genome

List  of genome IDs.

=item variant_code

List of variant codes.

=item genome_index

Hash mapping from genome ID to genome index.

=item spreadsheet

Spreadsheet data. Structured as a list of rows, each of  which
is a list of entries. An entry is a list of PEG numbers.

=item spreadsheet_inv

Inverted structure of spreadsheet - list of columns, each of which is a list
of rows.

=back

=head2 Public Methods

=head3 new

    my $sub = Subsystem->new($subName, $fig, $createFlag);

Load the subsystem. If it does not exist, and $createFlag is true, create
a new, empty subsystem.

=over 4

=item subName

Name of the desired subsystem.

=item fig

FIG object for accessing the SEED data store.

=item createFlag

TRUE if an empty subsystem should be created with the given name, else FALSE. If a
subsystem with the name already exists, this parameter has no effect.

=back

=cut

sub new
{
    my($class, $name, $fig, $create) = @_;

    my $ssa_dir = get_dir_from_name($name);
    #
    # For loading, the subsystem directory must already exist.
    #

    if (! -d $ssa_dir and not $create)
    {
        #           warn "Subsystem $name does not exist\n";
        return undef;
    }

    # RAE: Please do this:
    $name =~ s/^\s+//; $name =~ s/\s+$//;
    $name =~ s/ /_/g;

    my $self = {
        dir => $ssa_dir,
        name => $name,
        fig => $fig,
    };

    bless($self, $class);

    #
    # Check to see if the database we're running against has a variant column.
    #
    $self->detect_db_version();

    if ($create)
    {
        $self->create_subsystem();
    }
    else
    {
        $self->load();
    }

    return $self;
}

=head3 all_functions

    my $pegRoles = $sub->all_functions();

Return a hash of all the features in the subsystem. The hash maps each
feature ID to its functional assignment.

=cut

sub all_functions {
    my ($self) = @_;
    return $self->{fig}->function_of_bulk([keys %{$self->{peg_roles}}]);
}

sub detect_db_version
{
    my($self) = @_;
    my $db = $self->{fig}->db_handle();
    my $dbh = $db->{_dbh};
    local $dbh->{RaiseError} = 1;
    local $dbh->{PrintError} = 0;

    if (!$SS_version_checked)
    {
	eval {
	    my $x = $db->SQL("select variant from subsystem_index where subsystem = '' limit 1");
	};
	
	#
	# If this failed, it's an old database.
	#
	if ($@ =~ /variant/)
	{
	    warn "Please rerun index_subsystems: current table does not have a variant column\n";
	    $self->{old_database} = 1;
	}
	$SS_version_checked = 1;
    }
}


sub new_from_dir
{
    my($class, $dir, $fig) = @_;

    my $ssa_dir = $dir;
    my $name = $dir;
    $name =~ s,.*/,,;

    #
    # For loading, the subsystem directory must already exist.
    #

    my $self = {
        dir => $ssa_dir,
        name => $name,
        fig => $fig,
    };

    bless($self, $class);

    $self->load();

    return $self;
}

=head3 create_subsystem

Create a new subsystem. This creates the subsystem directory in the
correct place ($FIG_Config::data/Subsystems), and populates it with
the correct initial data.

=cut

sub create_subsystem
{
    my($self) = @_;

    my $dir = $self->{dir};
    my $fig = $self->{fig};

    if (-d $dir)
    {
        warn "Not creating: Subsystem directory $dir already exists";
        return;
    }

    $fig->verify_dir($dir);

    #
    # Initialize empty data structures.
    #

    $self->{genome} = [];
    $self->{genome_index} = {};
    $self->{variant_code} = [];

    $self->{abbr} = {};
    $self->{role_index} = {};
    $self->{roles} = [];
    $self->{role_abbrs} = [];

    $self->{spreadsheet} = [];
    $self->{spreadsheet_inv} = [];

    # added by DB for empty cell annotation
    $self->{emptycells} = {};

    $self->{col_subsets} = [];
    $self->{col_subset_members} = {};

    $self->{row_subsets} = [];
    $self->{row_subset_members} = {};
    $self->load_row_subsets();

    $self->{row_active_subset} = "All";
    $self->{col_active_subset} = "All";

    $self->{version} = 0;
    $self->{exchangable} = 0;
    $self->{classification} = [];

    $self->write_subsystem();
}

=head3 get_diagrams

    my @list = $sub->get_diagrams();

Return a list of the diagrams associated with this subsystem. Each diagram
is represented in the return list as a 4-tuple C<[diagram_id, diagram_name,
page_link, img_link]> where

=over 4

=item diagram_id

ID code for this diagram.

=item diagram_name

Displayable name of the diagram.

=item page_link

URL of an HTML page containing information about the diagram.

=item img_link

URL of an HTML page containing an image for the diagram.

=back

Note that the URLs are in fact for CGI scripts with parameters that point them
to the correct place.

=cut

sub get_diagrams {
    # Get the parameters.
    my ($self) = @_;
    # Look for all the sub-directories in the subsystem's diagram directory. Each
    # one of these will be a diagram ID. If the diagram directory doesn't exist,
    # we'll get an empty list.
    my @ids = GetDiagramIDs($self->{dir});
    # Declare the return variable.
    my @ret;
    # Loop through the diagram IDs found (if any).
    for my $id (@ids) {
        # Get the name and URLs for this diagram.
        my(@diag) = $self->get_diagram($id);
        # If we found a name and URLs, then this diagram must be added to the
        # return list.
        if (@diag) {
            push(@ret, [$id, @diag]);
        }
    }
    # Return the list of diagrams.
    return @ret;
}

=head3 get_diagram

    my ($name, $pageURL, $imgURL) = $sub->get_diagram($id);

Get the information (if any) for the specified diagram. The diagram corresponds
to a subdirectory of the subsystem's C<diagrams> directory. For example, if the
diagram ID is C<d03>, the diagram's subdirectory would be C<$dir/diagrams/d03>,
where I<$dir> is the subsystem directory. The diagram's name is extracted from
a tiny file containing the name, and then the links are computed using the
subsystem name and the diagram ID. The parameters are as follows.

=over 4

=item id

ID code for the desired diagram.

=item RETURN

Returns a three-element list. The first element is the diagram name, the second
a URL for displaying information about the diagram, and the third a URL for
displaying the diagram image.

=back

=cut

sub get_diagram
{
    my($self, $id) = @_;

    my $name = GetDiagramName($self->{dir}, $id);
    my ($link, $img_link) = ComputeDiagramURLs($self, $self->{name}, $id);

    return($name, $link, $img_link);
}

sub get_diagram_html_file
{
    my($self, $id) = @_;

    my $ddir = "$self->{dir}/diagrams/$id";

    return unless -d $ddir;

    my $html = "$ddir/diagram.html";

    if (-f $html)
    {
        return $html;
    }
    else
    {
        return undef;
    }
}

sub is_new_diagram {
  my ($self, $id) = @_;

  my $image_map = $self->get_diagram_html_file($id);
  if ($image_map) {

    open(IN, "$image_map") or die "Unable to open file $image_map.";
    my $header = <IN>;
    close(IN);

    if ($header =~ /\<map name=\"GraffleExport\"\>/) {
      return 1;
    }
  }

  return undef;
}

sub get_link_for_new_diagram {
  my ($self, $id) = @_;

  my $ss_name = $self->{name};
  return "./diagram.cgi?subsystem_name=$ss_name&diagram=$id";
}


sub open_diagram_image
{
    my($self, $id) = @_;

    my $img_base = "$self->{dir}/diagrams/$id/diagram";

    my @types = ([".png", "image/png"],
                 [".gif", "image/gif"],
                 [".jpg", "image/jpeg"]);

    for my $tent (@types)
    {
        my($ext, $type) = @$tent;

        my $file = "$img_base$ext";

        if (open(my $fh, "<$file"))
        {
            return($type, $fh);
        }
    }

    return undef;
}

sub delete_diagram
{
    my($self, $id) = @_;

    my $dir = "$self->{dir}/diagrams/$id";

    if (-d $dir)
    {
        system("rm", "-r", $dir);
    }
}

sub rename_diagram
{
    my($self, $id, $new_name) = @_;

    my $dir = "$self->{dir}/diagrams/$id";

    if (-d $dir)
    {
        open(F, ">$dir/NAME");
        $new_name =~ s/\n.*$//s;
        print F "$new_name\n";
        close(F);
    }
}

sub create_new_diagram
{
    my($self, $fh, $html_fh, $name, $id, $overwrite)  = @_;

    #
    # Get a new id.
    #

    my $dir = "$self->{dir}/diagrams";
    my $old_dir = "$self->{dir}/old_diagrams";

    Tracer::Insure($dir);
    Tracer::Insure($old_dir);

    my $path;

    if (defined($id))
    {
        #
        # Ensure this id doesn't already exist.
        #

        $path = "$dir/$id";

        if (-d $path)
        {
	    if (!$overwrite)
	    {
		confess "Diagram id $id already exists in subsystem $self->{name}";
	    }
	    else
	    {
		my $opath = "$old_dir/$id." . time;
		rename($path, $opath);
	    }
        }

    }
    else
    {
        $id = "d01";

        while (1)
        {
            $path = "$dir/$id";
            last unless -e $path;
            $id++;
        }
    }

    Tracer::Insure($path);

    if ($name)
    {
        open(F, ">$path/NAME");
        $name =~ s/\n.*$//s;
        print F "$name\n";
        close(F);
    }

    #
    # Write the file if we have one.
    #

    if ($fh)
    {
        my($ext, $buf);

        if (read($fh, $buf, 4096))
        {
            my($ext) = $self->classify_image_type($buf);
            open(D, ">$path/diagram$ext");
            print D $buf;

            while (read($fh, $buf, 4096))
            {
                print D $buf;
            }
            close(D);
        }
        close($fh);
    }

    #
    # And write the HTML file if we have one.
    #
    if ($html_fh)
    {
        my $buf;
        open(D, ">$path/diagram.html");

        while (read($html_fh, $buf, 4096))
        {
            print D $buf;
        }
        close(D);
        close($html_fh);
    }

    return $id;
}

sub create_new_illustration
{
    my( $self, $fh, $name, $id, $overwrite )  = @_;

    #
    # Get a new id.
    #

    my $dir = "$self->{dir}/diagrams";
    my $old_dir = "$self->{dir}/old_diagrams";

    Tracer::Insure($dir);
    Tracer::Insure($old_dir);

    my $path;

    if (defined($id))
    {
        #
        # Ensure this id doesn't already exist.
        #

        $path = "$dir/$id";

        if (-d $path)
        {
	    if (!$overwrite)
	    {
		confess "Diagram id $id already exists in subsystem $self->{name}";
	    }
	    else
	    {
		my $opath = "$old_dir/$id." . time;
		rename($path, $opath);
	    }
        }

    }
    else
    {
        $id = "d01";

        while (1)
        {
            $path = "$dir/$id";
            last unless -e $path;
            $id++;
        }
    }

    Tracer::Insure($path);

    if ($name)
    {
        open(F, ">$path/NAME");
        $name =~ s/\n.*$//s;
        print F "$name\n";
        close(F);
    }

    #
    # Write the file if we have one.
    #

    if ($fh)
    {
        my($ext, $buf);

        if (read($fh, $buf, 4096))
        {
            my($ext) = $self->classify_image_type($buf);
            open(D, ">$path/diagram$ext");
            print D $buf;

            while (read($fh, $buf, 4096))
            {
                print D $buf;
            }
            close(D);
        }
        close($fh);
    }

    return $id;
}

sub upload_new_image
{
    my($self, $id, $fh) = @_;

    if (!$fh)
    {
        warn "Subsystem::upload_new_image aborting: fh is undef\n";
        return;
    }


    my $dir = "$self->{dir}/diagrams/$id";

    if (not -d $dir)
    {
        warn "Subsystem::upload_new_image aborting: $dir does not exist\n";
        return;
    }

    #
    # remove any old diagram images.
    #

    for my $path (<$dir/diagram.{png,gif,jpg}>)
    {
        unlink($path);
    }

    my($ext, $buf);

    if (read($fh, $buf, 4096))
    {
        my($ext) = $self->classify_image_type($buf);

        if (!open(D, ">$dir/diagram$ext"))
        {
            warn "Subsystem::upload_new_image open failed for $dir/diagram$ext: $!\n";
            close($fh);
            return;
        }

        warn "Subsystem::upload_new_image classified new image as $ext\n";
        print D $buf;

        while (read($fh, $buf, 4096))
        {
            print D $buf;
        }
        close(D);
    }
    else
    {
        warn "Subsystem::upload_new_image read failed for $fh: $!\n";
    }

    warn "Subsystem::upload_new_image complete: " . `/bin/ls -l '$dir'`;

    close($fh);
}

sub upload_new_html
{
    my($self, $id, $fh) = @_;

    if (!$fh)
    {
        warn "Subsystem::upload_new_html aborting: fh is undef\n";
        return;
    }

    my $dir = "$self->{dir}/diagrams/$id";

    if (not -d $dir)
    {
        warn "Subsystem::upload_new_html aborting: $dir does not exist\n";
        return;
    }

    my($buf);

    if (!open(D, ">$dir/diagram.html"))
    {
        warn "Subsystem::upload_new_html open failed for $dir/diagram.html: $!\n";
        return;
    }

    my $rc;
    while ($rc = read($fh, $buf, 4096))
    {
        print D $buf;
    }
    if (!defined($rc))
    {
        warn "Subsystem::upload_new_html read failed for $fh: $!\n";
    }

    warn "Subsystem::upload_new_html complete: " . `/bin/ls -l '$dir'`;

    close(D);
    close($fh);
}

sub classify_image_type
{
    my($self, $buf) = @_;

    my $ext;

    #
    # Determine file type, for PNG / JPG / GIF. If we could be assured
    # the ImageMagick identify app worked properly, we'd use that instead.
    #
    # Maybe later.
    #

    if (substr($buf, 0, 8) eq "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a")
    {
        $ext = ".png";
    }
    elsif (substr($buf, 0, 3) eq "GIF")
    {
        $ext = ".gif";
    }
    elsif (substr($buf, 0, 2) eq "\xff\xd8" and substr($buf, 6, 4) eq "JFIF")
    {
        $ext = ".jpg";
    }
    else
    {
        warn "Unknown file type in new diagram\n";
        $ext = ".png";
    }

    return $ext;
}


#
# Synchronize the database index for this subsystem to the
# subsystem data.
#
# We assume the table already exists.
#

sub db_sync
{
    my($self, $skip_delete) = @_;

    if ($self->{empty_ss})
    {
	warn "Not synching empty subsystem $self->{name}\n";
	return;
    }

    my $fig = $self->{fig};
    my $rdbH = $fig->db_handle();

    if (!$skip_delete)
    {
        $self->delete_indices();
	$self->delete_aux();
    }

    my @aux_roles = map { $self->get_subsetC_roles($_) }
                    grep { $_ =~ /^(AUX|auxiliary)/ }
                    $self->get_subset_namesC;
    my %aux_roles;
    $aux_roles{$_} = 1 foreach @aux_roles;

    my $nameQ = quotemeta $self->{name};
    foreach my $role (@aux_roles)
    {
	my $roleQ = quotemeta $role;
	$rdbH->SQL("INSERT INTO aux_roles ( subsystem, role) VALUES ('$nameQ','$roleQ')");
    }

    my $tmp = "$FIG_Config::temp/ixsub.$$";
    open(TMP, ">$tmp") or die "Cannot open tmpfile $tmp: $!\n";

    my $tmp_roles_in_genome = "$FIG_Config::temp/roles_in_genome.$$";
    open(TMPROLES, ">", $tmp_roles_in_genome) or die "Cannot open tmpfile $tmp_roles_in_genome: $!";

    my $tmp_variant = "$FIG_Config::temp/variant.$$";
    open(TMPVARIANT, ">", $tmp_variant) or die "Cannot open tmpfile $tmp_variant: $!";

    my $tmp_nonaux_roles = "$FIG_Config::temp/nonaux.$$";
    open(TMPNONAUX, ">", $tmp_nonaux_roles) or die "Cannot open tmpfile $tmp_nonaux_roles: $!";

    #
    # We run thru all the cells, writing an entry in the database for the peg/subsystem/role.
    #

    my @roles = $self->get_roles();
    for my $role (grep { ! $aux_roles{$_} } @roles)
    {
	print TMPNONAUX "$self->{name}\t$role\n";
    }
    for my $genome ($self->get_genomes())
    {
	my @roles_in_genome;
	my $gidx = $self->get_genome_index($genome);
	my $variant = $self->get_variant_code($gidx);
	print TMPVARIANT "$self->{name}\t$genome\t$variant\n";

	# print "Index $genome variant=$variant\n";
	my $row = $self->get_row($gidx);

	for my $i (0..$#$row)
	{
	    my $cell = $row->[$i];
	    my $role = $roles[$i];
	    if ($cell)
	    {
		if (@$cell && $variant ne '0' && $variant !~ /\*/ && !$aux_roles{$role})
		{
		    print TMPROLES "$self->{name}\t$genome\t$role\n";
		}
		    
		for my $peg (@$cell)
		{
		    # $sth->execute($peg, $self->{name}, $role);
		    if ($self->{old_database})
		    {
			print TMP "$peg\t$self->{name}\t$role\n";
		    }
		    else
		    {
			print TMP "$peg\t$self->{name}\t$role\t$variant\n";
		    }
		}
	    }
	}
    }
    close(TMPROLES);
    close(TMPVARIANT);
    close(TMPNONAUX);
    close(TMP);
    $rdbH->load_table(file => $tmp,
                      tbl => 'subsystem_index');
    $rdbH->load_table(file => $tmp_nonaux_roles,
		      tbl => 'subsystem_nonaux_role');
    $rdbH->load_table(file => $tmp_roles_in_genome,
		      tbl => 'subsystem_genome_role');
    $rdbH->load_table(file => $tmp_variant,
		      tbl => 'subsystem_genome_variant');

    #
    # Update the metadata table with this subsystem's info.
    #
    my $classification = $self->get_classification;
    my @info = (join("\t", @$classification),
		$classification->[0],
		$classification->[1],
		$self->get_curator,
		$self->get_created,
		$self->get_last_updated,
		$self->get_version,
		$self->{exchangable});

    my $dbh = $rdbH->{_dbh};
    my $n;
    eval {
	local($dbh->{RaiseError})=1;
	$n = $dbh->do(qq(UPDATE subsystem_metadata
		    SET classification = ?, class_1 = ?, class_2 = ?,
		    	curator = ?, creation_date = ?, last_update = ?,
		    	version = ?, exchangable = ?
		    WHERE subsystem = ?),
		 undef, @info, $self->get_name());
    };
    if ($@ or $n == 0)
    {
	#
	# Row probably missing.
	#
	warn "error on update: $@" if $@;
	eval {
	    $dbh->do(qq(INSERT INTO subsystem_metadata(subsystem, classification,class_1, class_2 ,
						       curator, creation_date, last_update,
						       version, exchangable)
			VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)),
		     undef, $self->get_name(), @info);
	};
    }
    $fig->mark_subsystems_modified();
}

#
# Delete this subsystem's entries from the database index.
#
sub delete_indices
{
    my($self) = @_;

    my $rdbH = $self->{fig}->db_handle();

    $rdbH->SQL("DELETE FROM subsystem_index where subsystem = ?", undef, $self->{name});
    $rdbH->SQL("DELETE FROM subsystem_metadata WHERE subsystem = ?", undef, $self->{name});
    $rdbH->SQL("DELETE FROM subsystem_genome_variant WHERE subsystem = ?", undef, $self->{name});
    $rdbH->SQL("DELETE FROM subsystem_genome_role WHERE subsystem = ?", undef, $self->{name});
    $rdbH->SQL("DELETE FROM subsystem_nonaux_role WHERE subsystem = ?", undef, $self->{name});
}

sub delete_aux
{
    my($self) = @_;

    my $rdbH = $self->{fig}->db_handle();

    $rdbH->SQL("DELETE FROM aux_roles where subsystem = ?", undef, $self->{name});
}

sub is_aux_role {
    my($self,$role) = @_;

    my $rdbH = $self->{fig}->db_handle();
    if (! $rdbH->table_exists('aux_roles')) { return 0 }
    my $nameQ = quotemeta $self->{name};
    my $roleQ = quotemeta $role;
    my $q = "SELECT subsystem FROM aux_roles WHERE subsystem = '$nameQ' AND role = '$roleQ'";

    my $relational_db_response;
    return (($relational_db_response = $rdbH->SQL($q)) && (@$relational_db_response > 0));
}


sub load
{
    my($self) = @_;

    #
    # Load the subsystem.
    #
    
    my $ssa;
    if  (!open($ssa,"<$self->{dir}/spreadsheet"))
    {
	Trace("Spreadsheet does not exist in subsystem $self->{name}") if T(1);
	$self->{empty_ss}++;
        return;
    }

    local $/ = "//\n";

    my $roles = <$ssa>;
    if ($roles)
    {
        $roles =~ s,$/$,,;
        #
        # Split on newline, filter for non-empty lines.
        #
        my @roles = split("\n", $roles);

        @roles = grep { $_ ne "" } @roles;

        $self->load_roles(@roles);
    }

    my $subsets = <$ssa>;
    if ($subsets)
    {
        $subsets =~ s,$/$,,;
        $self->load_subsets($subsets);
    }

    $/ = "\n";

    $self->load_row_subsets();
        
    $self->load_genomes($ssa);

    # now load the empty cell information
    $self->load_emptycells();

    #
    # Now load the rest of the info.
    #

    $self->load_reactions();
    $self->load_hope_kegg_info();
    $self->load_hope_reactions();
    $self->load_hope_reaction_notes();
    $self->load_hope_reaction_links();
    $self->load_hope_curation_notes();
    $self->load_notes();
    $self->load_classification();
    $self->load_version();
    $self->load_exchangable();
    $self->load_curation();

    return 1;
}

sub load_emptycells
{
    my($self) = @_;

    my $absencehash = {};
    if  (open(ECS,"<$self->{dir}/emptycells"))
    {
        while (defined($_ = <ECS>))
        {
	    chomp;
	  my ( $frabbr, $genome, $value ) = split( "\t", $_ );
	  $absencehash->{ $frabbr }->{ $genome } = $value;
        }
        close(ECS);
    }
    $self->{emptycells} = $absencehash;
}


sub load_notes
{
    my($self) = @_;

    #$self->{notes} = &FIG::file_read(File::Spec->catfile($self->{dir}, "notes"));

    my $nh = new FileHandle(File::Spec->catfile($self->{dir}, "notes"));

    $nh or return;

    my $section = "NOTES";
    my $text;

    my $all;
    $_ = <$nh>;
    while (defined($_))
    {
	$all .= $_;
	if (/^$notes_separator$/)
	{
	    #
	    # Next line after the separator is the new section name, a single word. If the
	    # line doesn't match that, push it (and the current line) into the text and continue.
	    #
	    my $new_section = <$nh>;
	    if ($new_section =~ /^(\w+)\s*$/)
	    {
		$new_section = $1;
		$self->handle_note_section($section, $text) if defined($text);
		$section = $new_section;
		$text = '';
		$_ = <$nh>;
	    }
	    else
	    {
		$text .= $_;
		$_ = $section;
	    }
	}
	else
	{
	    $text .= $_;
	    $_ = <$nh>;
	}
    }
    $self->handle_note_section($section, $text);
    $self->{raw_notes} = $all;
}

sub handle_note_section
{
    my($self, $section, $text) = @_;

    # print "Got section $section text=$text\n";

    my $sname = lc($section);
    if (defined($defined_sections{$sname}))
    {
	$self->{$sname} = $text;
    }
    else
    {
	$self->{other_sections}->{$section} = $text;
	warn "Loaded unknown section name $section\n";
    }
}

sub load_hope_kegg_info
{
    my($self) =@_;

    $self->{hope_scenarios} = {};

    if (open(HOPE_KEGG,"<$self->{dir}/hope_kegg_info"))
    {
        my @lines = <HOPE_KEGG>;

	for (my $i = 0; $i < scalar @lines; $i += 6)
	{
	    if (defined $lines[$i])
	    {
		chomp $lines[$i];
		$self->add_hope_scenario($lines[$i]);
	    }
	    if (defined $lines[$i+1])
	    {
		chomp $lines[$i+1];
		$self->set_hope_input_compounds($lines[$i], $lines[$i+1]);
	    }
	    if (defined $lines[$i+2])
	    {
		chomp $lines[$i+2];
		$self->set_hope_output_compounds($lines[$i], $lines[$i+2]);
	    }
	    if (defined $lines[$i+3])
	    {
		chomp $lines[$i+3];
		$self->set_hope_map_ids($lines[$i], $lines[$i+3]);
	    }
	    if (defined $lines[$i+4])
	    {
		chomp $lines[$i+4];
		$self->set_hope_additional_reactions($lines[$i], $lines[$i+4]);
	    }
	    if (defined $lines[$i+5])
	    {
		chomp $lines[$i+5];
		$self->set_hope_ignore_reactions($lines[$i], $lines[$i+5]);
	    }
	}

        close(HOPE_KEGG);
    }
}

sub load_hope_curation_notes
{
    my($self) = @_;

    $self->{hope_curation_notes} = &FIG::file_read(File::Spec->catfile($self->{dir}, "hope_curation_notes"));
}

sub load_reactions
{
    my($self) = @_;

    my $reactions = undef;
    if (open(REACT,"<$self->{dir}/reactions"))
    {
        while (defined($_ = <REACT>))
        {
            if ($_ =~ /^(\S.*\S)\t(\S+)/)
            {
                push(@{$reactions->{$1}},split(/,\s*/,$2));
            }
        }
        close(REACT);
    }

    $self->{reactions} = $reactions;
}
sub load_hope_reactions
{
    my($self) = @_;

    my $hope_reactions = FIGRules::GetHopeReactions($self, $self->{dir});

    $self->{hope_reactions} = $hope_reactions;
}

sub load_hope_reaction_notes
{
    my($self) = @_;

    my $hope_reaction_notes = {};
    if (open(HOPE_REACTION_NOTES,"<$self->{dir}/hope_reaction_notes"))
    {
        while (defined($_ = <HOPE_REACTION_NOTES>))
        {
            if ($_ =~ /^(\S.*\S)\t(.+)$/)
            {
                $hope_reaction_notes->{$1} = $2;
            }
        }
        close(HOPE_REACTION_NOTES);
    }

    $self->{hope_reaction_notes} = $hope_reaction_notes;
}

sub load_hope_reaction_links
{
    my($self) = @_;

    my $hope_reaction_links = {};
    if (open(HOPE_REACTION_LINKS,"<$self->{dir}/hope_reaction_links"))
    {
        while (defined($_ = <HOPE_REACTION_LINKS>))
        {
            if ($_ =~ /^(\S.*\S)\t(.+)$/)
            {
                $hope_reaction_links->{$1} = $2;
            }
        }
        close(HOPE_REACTION_LINKS);
    }

    $self->{hope_reaction_links} = $hope_reaction_links;
}

sub load_classification
{
    my($self) = @_;

    my $class = &FIG::file_read(File::Spec->catfile($self->{dir}, "CLASSIFICATION"));
    my @tmp = grep { $_ =~ /\S/ } split(/\n/,$class);
    $class = join("\n",@tmp);
    if ($class) {$self->{classification} = [split /\t/, $class]} else {$self->{classification} = ['', '', '']}
}

sub load_curation
{
    my($self) = @_;

#    my @l = &FIG::file_head(File::Spec->catfile($self->{dir}, "curation.log"), 1);
#
#    $_ = $l[0];
#    chomp;

    if  (open(LOG,"<$self->{dir}/curation.log"))
    {
	my $last = 0;
        while (defined($_ = <LOG>))
        {
            if (/^(\d+)\t(\S+)\s+started/)
            {
                $self->{curator} = $2;
		$self->{created} = $1;
            }
	    if ((/^(\d+)/) && ($1 > $last))
	    {
		$last = $1;
	    }
        }
        close(LOG);
	if ($last) { $self->{last_updated} = $last; }
    }

    # 
    # If a file OWNER exists in the subsys dir, use that instead.
    #
}

sub load_version
{
    my($self) = @_;

    my @l = &FIG::file_head(File::Spec->catfile($self->{dir}, "VERSION"), 1);
    my $l = $l[0];
    chomp $l;
    $self->{version} = $l;
}

sub load_exchangable
{
    my($self) = @_;

    my $file = File::Spec->catfile($self->{dir}, "EXCHANGABLE");

    if (-f $file)
    {
        my($l, @l);

        @l = &FIG::file_head($file, 1);
        $l = $l[0];
        chomp $l;
        $self->{exchangable} = $l;
    }
    else
    {
        $self->{exchangable} = 0;
    }
}


sub load_roles
{
    my($self, @roles) = @_;

    $self->{abbr} = {};
    $self->{role_index} = {};
    $self->{roles} = [];
    $self->{role_abbrs} = [];

    my $i = 0;
    for my $role (@roles)
    {
        my($abbr, $name) = split(/\t/, $role);
        $abbr =~ s/^\s+//;
        $abbr =~ s/\s+$//;
        $name =~ s/^\s+//;
        $name =~ s/\s+$//;
        # print "Role $i: abbr=$abbr name=$name\n";

        $self->{abbr}->{$abbr} = $name;
        $self->{role_index}->{$name} = $i;
        $self->{roles}->[$i] = $name;
        $self->{role_abbrs}->[$i] = $abbr;
        $i++;
    }
}

sub load_subsets
{
    my($self, $subsets) = @_;

    #
    # Column and row subsets.
    #
    my($subsetsC, $subsetsR) = split(/\n\n/, $subsets);

    #
    # Handle column subsets.
    #

    my @subsetsC = split(/\n/, $subsetsC);

    #
    # Determine active subset.
    #

    my $active_subsetC;
    if (@subsetsC > 0)
    {
        $active_subsetC = pop(@subsetsC);
    }
    else
    {
        $active_subsetC = 'All';
    }

    $self->{col_active_subset} = $active_subsetC;

    $self->{col_subsets} = [];
    $self->{col_subset_members} = {};

    for my $subset (@subsetsC)
    {
        my($name, @members) = split(/\s+/, $subset);

        #
        # File format has members 1-based.
        #

        @members = map { $_ - 1 } @members;

        push(@{$self->{col_subsets}}, $name);

        #
        # Map role members from name to index if necessary.
        #
        # Is it really necessary? ssa2 code was looking up in %pos for this.
        #
        @members = map {
            if (my $new = $self->{role_index}->{$_})
            {
                $new;
            }
            else
            {
                $_;
            }
        } @members;

        @{$self->{col_subset_members}->{$name}} = @members;
    }

    #
    # Now the row subsets.
    #

    chomp($subsetsR);

    if ($subsetsR =~ /(\S+.*\S+)/)
    {
        $self->{row_active_subset} = $1;
    }
    else
    {
        $self->{row_active_subset} = 'All';
    }
    $self->{row_subsets} = [];
}

sub load_genomes
{
    my($self, $fh) = @_;
    my(%seen);

    $self->{spreadsheet} = [];
    $self->{spreadsheet_inv} = [];
    $self->{genome} = [];
    $self->{genome_index} = {};
    $self->{variant_code} = [];
    $self->{peg_roles} = {};

    # cache the genome list to reduce the number of sql queries
    my %all_genomes = map { $_ => 1 } $self->{fig}->genomes;

    my $nr = @{$self->{roles}};

    my $i = 0;
    while (<$fh>)
    {
        next if ($_ =~ /^\/\//);
        chomp;

        my($genome, $variant_code, @row) = split(/\t/, $_, $nr + 2);
        $variant_code =~ s/ //g;
        next if ($seen{$genome} || (($genome =~ /^(\d+\.\d+)/) && (! $all_genomes{$1})));
        $seen{$genome}++;

	$genome =~ /^(\d+\.\d+)/;
	my $just_genome = $1;

        my $j = 0;

        $self->{genome}->[$i] = $genome;
        $self->{genome_index}->{$genome} = $i;
        $self->{variant_code}->[$i] = $variant_code;

        my $thislen = @row;

#       if ($thislen != $nr)
#       {
#           warn "Genome $genome has wrong column count ($thislen != $nr)\n";
#           warn "<$_> $genome $variant_code '", join(":", @row), "'\n";
#       }

        for my $j (0..$nr - 1)
        {
            my $entry = $row[$j];
# OLD       my $e2 = [map("fig|$genome.peg.$_", split(/,/, $entry))];
            my $e2 = [map { ($_ =~ /^[a-zA-Z]+\.\d+$/) ? "fig|$just_genome.$_" : "fig|$just_genome.peg.$_" }
                      split(/,/, $entry)
		     ];
            $self->{spreadsheet}->[$i]->[$j] = $e2;
            $self->{spreadsheet_inv}->[$j]->[$i] = $e2;
            for my $fidj (@{$e2}) {
                push @{$self->{peg_roles}->{$fidj}}, $j;
            }
            $j++;
        }
        $i++;

    }
}

sub add_virtual_genome
{
    my($self, $name, $genome, $variant, $bindings) = @_;

    my $gidx = @{$self->{genome}};

    $self->{genome}->[$gidx] = $genome;
    $self->{genome_index}->{$genome} = $gidx;
    $self->{variant_code}->[$gidx] = $variant;

    for my $role (keys %$bindings)
    {
	my $role_idx = $self->get_role_index($role);

	my $plist = $bindings->{$role};

#	warn "Role $role maps to $role_idx with pegs @$plist\n";

	$self->{spreadsheet}->[$gidx]->[$role_idx] = $plist;
	$self->{spreadsheet_inv}->[$role_idx]->[$gidx] = $plist;
	for my $peg (@$plist)
	{
#	    warn "Peg $peg => $role_idx\n";
	    push(@{$self->{peg_roles}->{$peg}}, $role_idx);
	}
    }
#    warn "After add virtual: \n",  Dumper($self), "\n";
}

=head3 in_genome($genome,$peg)

    if ($sub->in_genome($genome,$peg))
    {
         process a PEG from the genome or region of a genome
    }	     

Return a boolean: "true" -> PEG falls within the genome or region of a genome

=over 4

=item genome

either \d+\.\d+ (a typical genome ID) OR
       \d+\.\d+:Contig_Beg_End

where Contig is a string of non-whitespace characters
      Beg    is an integer
      End    is an integer

=item peg

ID of the peg we are checking

=item RETURN

Returns a boolean

=back

=cut

sub in_genome {
    my ($self, $genome,$fid) = @_;

    if ($genome =~ /^(\d+\.\d+)(:(\S+)_(\d+)_(\d+))?$/)
    {
	my $just_genome = $1;
	my($contig,$beg,$end) = $2 ? ($3,$4,$5) : (undef,undef,undef);
	my $fidG = &FIG::genome_of($fid);
	if (! $contig) { return ($just_genome eq $fidG) }
	my $fig = $self->{fig};
	my $loc = $fig->feature_location($fid);
	my($contig1,$beg1,$end1) = $fig->boundaries_of($loc);
	return (($contig1 eq $contig) && 
		&FIG::between($beg,$beg1,$end) && 
		&FIG::between($beg,$end1,$end));
    }
    else
    {
	return 0;
    }
}

=head3 get_peg_roles

    my @cols = $sub->get_peg_roles($peg);

Return the column numbers in which the specified PEG appears.

=over 4

=item peg

ID of the feature whose roles are desired.

=item RETURN

Returns a list of the column numbers in which the peg appears, or an empty
list if it is not found.

=back

=cut

sub get_peg_roles {
    # Get the parameters.
    my ($self, $peg) = @_;
    # Declare the return variable.
    my @retVal;
    # Find this peg's roles.
    if (exists $self->{peg_roles}->{$peg}) {
        @retVal = @{$self->{peg_roles}->{$peg}};
    }
    # Return the result.
    return @retVal;
}

=head3 get_all_pegs

    my @pegs = $sub->get_all_pegs();

Return all pegs appearing in the subsystem.

=cut

sub get_all_pegs {
    my ($self) = @_;
    return keys %{$self->{peg_roles}};
}

=head3 write_subsystem

Write the subsystem to the disk.  Updates on-disk data with notes,
etc. Perform backups when necessary.

=cut

sub write_subsystem
{
    my($self, $force_backup) = @_;

    my $dir = $self->{dir};
    my $fig = $self->{fig};

    #
    # We first move the existing spreadsheet and notes files (if present)
    # to spreadsheet~ and notes~, and current state.
    #

    my $ss_file = "$dir/spreadsheet";
    my $ss_bak = "$dir/spreadsheet~";
    my $notes_file = "$dir/notes";
    my $notes_bak = "$dir/notes~";
    my $reactions_file = "$dir/reactions";
    my $reactions_bak = "$dir/reactions~";
    my $hope_kegg_info_file = "$dir/hope_kegg_info";
    my $hope_kegg_info_bak = "$dir/hope_kegg_info~";
    my $hope_reactions_file = "$dir/hope_reactions";
    my $hope_reactions_bak = "$dir/hope_reactions~";
    my $hope_reaction_notes_file = "$dir/hope_reaction_notes";
    my $hope_reaction_notes_bak = "$dir/hope_reaction_notes~";
    my $hope_reaction_links_file = "$dir/hope_reaction_links";
    my $hope_reaction_links_bak = "$dir/hope_reaction_links~";
    my $hope_curation_notes_file = "$dir/hope_curation_notes";
    my $hope_curation_notes_bak = "$dir/hope_curation_notes~";
    my $emptycells_file = "$dir/emptycells";
    my $emptycells_bak = "$dir/emptycells~";
    my $classification_file = "$dir/CLASSIFICATION";

    if (-f $ss_file)
    {
        rename($ss_file, $ss_bak);
    }

    if (-f $notes_file)
    {
        rename($notes_file, $notes_bak);
    }

    if (-f $reactions_file)
    {
        rename($reactions_file, $reactions_bak) or warn "rename $reactions_file $reactions_bak failed $!";
#       print STDERR "wrote $reactions_bak\n";
    }

    if( -f $hope_kegg_info_file)
    {
    	rename($hope_kegg_info_file, $hope_kegg_info_bak) or warn "rename $hope_kegg_info_file $hope_kegg_info_bak failed $!";
    }

    if (-f $hope_reactions_file)
    {
        rename($hope_reactions_file, $hope_reactions_bak) or warn "rename $hope_reactions_file $hope_reactions_bak failed $!";
#       print STDERR "wrote $hope_reactions_bak\n";
    }

    if (-f $hope_reaction_notes_file)
    {
        rename($hope_reaction_notes_file, $hope_reaction_notes_bak) or warn "rename $hope_reaction_notes_file $hope_reaction_notes_bak failed $!";
#       print STDERR "wrote $hope_reaction_notes_bak\n";
    }

    if (-f $hope_reaction_links_file)
    {
        rename($hope_reaction_links_file, $hope_reaction_links_bak) or warn "rename $hope_reaction_links_file $hope_reaction_links_bak failed $!";
#       print STDERR "wrote $hope_reaction_links_bak\n";
    }

    if (-f $hope_curation_notes_file)
    {
        rename($hope_curation_notes_file, $hope_curation_notes_bak) or warn "rename $hope_curation_notes_file $hope_curation_notes_bak failed $!";
#       print STDERR "wrote $hope_curation_notes_bak\n";
    }

    if (-f $emptycells_file)
    {
        rename($emptycells_file, $emptycells_bak) or warn "rename $emptycells_file $emptycells_bak failed $!";
#       print STDERR "wrote $hope_curation_notes_bak\n";
    }

    #
    # Eval this whole chunk, so that if we get any fatal errors, we can
    # roll back to the old saved data.
    #

    eval {
        my $fh;
        open($fh, ">$ss_file") or die "Cannot open $ss_file for writing: $!\n";
        $self->write_spreadsheet($fh);
        close($fh);
        chmod(0777,$ss_file);

        open($fh, ">$notes_file") or die "Cannot open $notes_file for writing: $!\n";
	for my $section (@section_order)
	{
	    print $fh "$notes_separator\n";
	    print $fh uc($section) . "\n";
	    my $textsection = $self->{ $section };
	    $textsection =~ s/\r//g;
	    chomp $textsection;	   
	    print $fh $textsection."\n";
	}
	for my $osection (keys %{$self->{other_sections}})
	{
	    print $fh "$notes_separator\n";
	    print $fh uc($osection) . "\n";
	    my $textosection = $self->{other_sections}->{$osection};
	    chomp $textosection;
	    $textosection =~ s/\r//g;
	    print $fh $textosection."\n";
#	    print $fh $self->{other_sections}->{$osection};
	}
        # print $fh "$self->{notes}";

        close($fh);
        chmod(0777,$notes_file);

        open($fh, ">$reactions_file") or die "Cannot open $reactions_file for writing: $!\n";
        my $reactions = $self->{reactions};
        foreach $_ (sort keys(%$reactions))
        {
            print $fh "$_\t" . join(",", @{$reactions->{$_}}), "\n";
        }
        close($fh);
        chmod(0777,$reactions_file);

	open($fh, ">$hope_kegg_info_file") or die "Cannot open $hope_kegg_info_file for writing: $!\n";
	foreach my $scenario_name (keys %{$self->{hope_scenarios}})
	{
	    print $fh $scenario_name, "\n";
	    my $scenario = $self->{hope_scenarios}->{$scenario_name};
	    my $input_compounds = $scenario->{input_compounds};
	    my $temp = join ",", @$input_compounds;
	    print $fh $temp , "\n";
	    my $output_compounds = $scenario->{output_compounds};
	    my @output_compounds_lists;
	    foreach my $cpd_list (@$output_compounds)
	    {
		if (scalar @$cpd_list > 1)
		{
		    push @output_compounds_lists, "(".join(",",@$cpd_list).")";
		}
		else
		{
		    push @output_compounds_lists, @$cpd_list;
		}
	    }
	    $temp = join ",", @output_compounds_lists;
	    print $fh $temp , "\n";
	    my $map_ids = $scenario->{map_ids};
	    $temp = join ",", @$map_ids;
	    print $fh $temp , "\n";
	    my $additional_reactions = $scenario->{additional_reactions};
	    $temp = join ",", @$additional_reactions;
	    print $fh $temp , "\n";
	    my $ignore_reactions = $scenario->{ignore_reactions};
	    $temp = join ",", @$ignore_reactions;
	    print $fh $temp , "\n";
	}

        close($fh);
        chmod(0777,$hope_kegg_info_file);

        open($fh, ">$hope_reactions_file") or die "Cannot open $hope_reactions_file for writing: $!\n";
        my $hope_reactions = $self->{hope_reactions};
        foreach $_ (sort keys(%$hope_reactions))
        {
            print $fh "$_\t" . join(",", @{$hope_reactions->{$_}}), "\n";
        }
        close($fh);
        chmod(0777,$hope_reactions_file);

        open($fh, ">$hope_reaction_notes_file") or die "Cannot open $hope_reaction_notes_file for writing: $!\n";
        my $hope_reaction_notes = $self->{hope_reaction_notes};
        foreach $_ (sort keys(%$hope_reaction_notes))
        {
            print $fh "$_\t" . $hope_reaction_notes->{$_}, "\n";
        }
        close($fh);
        chmod(0777,$hope_reaction_notes_file);

        open($fh, ">$hope_reaction_links_file") or die "Cannot open $hope_reaction_links_file for writing: $!\n";
        my $hope_reaction_links = $self->{hope_reaction_links};
        foreach $_ (sort keys(%$hope_reaction_links))
        {
            print $fh "$_\t" . $hope_reaction_links->{$_}, "\n";
        }
        close($fh);
        chmod(0777,$hope_reaction_links_file);

        open($fh, ">$hope_curation_notes_file") or die "Cannot open $hope_curation_notes_file for writing: $!\n";
        print $fh "$self->{hope_curation_notes}";
        close($fh);
        chmod(0777,$hope_curation_notes_file);

        open($fh, ">$emptycells_file") or die "Cannot open $emptycells_file for writing: $!\n";
	my $gahash = $self->{emptycells};
	foreach my $k1 ( keys %$gahash ) {
	  foreach my $k2 ( keys %{ $gahash->{ $k1 } } ) {
	    print $fh $k1."\t".$k2."\t".$gahash->{ $k1 }->{ $k2 }."\n";
	  }
	}
        close($fh);
        chmod(0777,$emptycells_file);

        open($fh, ">$classification_file") or die "Can not open $classification_file for writing: $!\n";
        print $fh join "\t", (@{$self->{classification}}), "\n";
        close($fh);
        chmod(0777,$classification_file);

        $self->update_curation_log();

        #
        # Write out the piddly stuff.
        #

        open($fh, ">$dir/EXCHANGABLE") or die "Cannot write $dir/EXCHANGABLE: $!\n";
        print $fh "$self->{exchangable}\n";
        close($fh);
        chmod(0777,"EXCHANGABLE");

        #
        # Process backup files. This is the smae process that determines when the
        # version number should be bumped, so write the version file afterward.
        #

        $self->update_backups($force_backup);

        if ($self->{version} < 100) { $self->{version} += 100 }
        open($fh, ">$dir/VERSION") or die "Cannot write $dir/VERSION: $!\n";
        print $fh "$self->{version}\n";
        close($fh);
        chmod(0777,"VERSION");
    };

    if ($@ ne "")
    {
        warn "Spreadsheet write failed, reverting to backup. Error was\n$@\n";
    }

}

sub update_curation_log
{
    my($self) = @_;

    my $fh;
    my $file = "$self->{dir}/curation.log";

    my $now = time;
    my $user = $self->{fig}->get_user();

    if (-f $file)
    {
        open($fh, ">>$file") or die "Cannot open $file for writing: $!\n";
    }
    else
    {
        open($fh, ">$file") or die "Cannot open $file for writing: $!\n";
        print $fh "$now\t$user\tstarted\n";
    }
    print $fh "$now\t$user\tupdated\n";
    close($fh);
}

sub update_backups
{
    my($self, $force_backup) = @_;

    my $dir = $self->{dir};
    my $fig = $self->{fig};

    my $ss_file = "$dir/spreadsheet";
    my $ss_bak = "$dir/spreadsheet~";
    my $notes_file = "$dir/notes";
    my $notes_bak = "$dir/notes~";
    my $reactions_file = "$dir/reactions";
    my $reactions_bak = "$dir/reactions~";
    my $hope_reactions_file = "$dir/hope_reactions";
    my $hope_reactions_bak = "$dir/hope_reactions~";
    my $hope_reaction_notes_file = "$dir/hope_reaction_notes";
    my $hope_reaction_notes_bak = "$dir/hope_reaction_notes~";
    my $hope_reaction_links_file = "$dir/hope_reaction_links";
    my $hope_reaction_links_bak = "$dir/hope_reaction_links~";
    my $hope_curation_notes_file = "$dir/hope_curation_notes";
    my $hope_curation_notes_bak = "$dir/hope_curation_notes~";
    my $emptycells_file = "$dir/emptycells";
    my $emptycells_bak = "$dir/emptycells~";
    my $hope_kegg_info_file = "$dir/hope_kegg_info"; 
    my $hope_kegg_info_bak = "$dir/hope_kegg_info~";

    my $ss_diff = abs((-s $ss_file) - (-s $ss_bak));
    my $notes_diff = abs((-s $notes_file) - (-s $notes_bak));
    my $reactions_diff = (system("cmp",  "-s",  $reactions_file, $reactions_bak) != 0);
    my $hope_reactions_diff = (system("cmp",  "-s",  $hope_reactions_file, $hope_reactions_bak) != 0);
    my $hope_reaction_notes_diff = (system("cmp",  "-s",  $hope_reaction_notes_file, $hope_reaction_notes_bak) != 0);
    my $hope_reaction_links_diff = (system("cmp",  "-s",  $hope_reaction_links_file, $hope_reaction_links_bak) != 0);
    my $hope_curation_notes_diff = (system("cmp",  "-s",  $hope_curation_notes_file, $hope_curation_notes_bak) != 0);
    my $emptycells_diff = (system("cmp",  "-s",  $emptycells_file, $emptycells_bak) != 0);
    my $hope_kegg_info_diff = (system("cmp",  "-s",  $hope_kegg_info_file, $hope_kegg_info_bak) != 0);

    if ($force_backup or ($ss_diff > 10) or ($notes_diff > 10) or $reactions_diff or $hope_reactions_diff or $hope_reaction_notes_diff or $hope_reaction_links_diff or $hope_curation_notes_diff or $hope_kegg_info_diff)
    {
        $self->make_backup();
    }
}

sub make_backup
{
    my($self) = @_;

    my $dir = $self->{dir};
    my $bak = "$dir/Backup";

    $self->{fig}->verify_dir($bak);

    my $ts = time;

    rename("$dir/spreadsheet~", "$bak/spreadsheet.$ts");
    rename("$dir/notes~", "$bak/notes.$ts");
    rename("$dir/reactions~", "$bak/reactions.$ts");
    rename("$dir/hope_reactions~", "$bak/hope_reactions.$ts");
    rename("$dir/hope_reaction_notes~", "$bak/hope_reaction_notes.$ts");
    rename("$dir/hope_reaction_links~", "$bak/hope_reaction_links.$ts");
    rename("$dir/hope_curation_notes~", "$bak/hope_curation_notes.$ts");
    rename("$dir/emptycells~", "$bak/emptycells.$ts");
    rename("$dir/hope_kegg_info~", "$bak/hope_kegg_info.$ts");
    $self->{version}++;
}



=head3 write_spreadsheet

    $sub->write_spreadsheet($fh);

Write the spreadsheet for this subsystem to filehandle $fh.

=cut

sub write_spreadsheet
{
    my($self, $fh) = @_;

    $self->_write_roles($fh);
    print $fh "//\n";

    $self->_write_subsets($fh);
    print $fh "//\n";

    $self->_write_spreadsheet($fh);
}

sub _write_roles
{
    my($self, $fh) = @_;

    my(@roles, @abbrs);

    @roles = $self->get_roles();
    @abbrs = $self->get_abbrs();

    #
    # Check abbreviations for validity. We disallow spaces, commas, and colons,
    # and enforce uniqueness.
    #

    my %abbrs;
    map { s/[\s,:]*//g; } @abbrs;
    map { $abbrs{$_}++ } @abbrs;

    for (my $i = 0; $i < @abbrs; $i++)
    {
	my $a = $abbrs[$i];
	if ($abbrs{$a} > 1)
	{
	    #
	    # abbrev is not unique
	    #
	    $a = "${a}_" . ($i + 1);
	    $abbrs[$i] = $a;
	}
    }

    while (@roles)
    {
        my $role = shift(@roles);
        my $abbr = shift(@abbrs);

        print $fh "$abbr\t$role\n";
    }
}

sub _write_subsets
{
    my($self, $fh) = @_;

    for my $sub ($self->get_subset_namesC())
    {
        next if ($sub eq "All");
        my @members= $self->get_subsetC($sub);

        #
        # member list on disk is 1-based
        #

        @members = map { $_ + 1 } @members;

        print $fh join("\t", $sub, @members), "\n";
    }
    my $active_row_subset = $self->{row_active_subset};
    my $active_col_subset = $self->{col_active_subset};

    print $fh "$active_col_subset\n";

    #
    # separator
    #

    print $fh "\n";

    #
    # genome subsets.
    #

    print $fh "$active_row_subset\n";
}

sub _write_spreadsheet
{
    my($self, $fh) = @_;

    my(@genomes);

    @genomes= $self->get_genomes();

    for (my $i = 0; $i < @genomes; $i++)
    {
        my $genome = $genomes[$i];
        my $vc = $self->get_variant_code($i);

        my $row = $self->get_row($i);

        if ($vc eq "")
        {
            $vc = "0";
        }

	#
	# Validate genome before writing.
	#

	if ($genome !~ /^\d+\.\d+/)
	{
	    next;
	}

	# Not sure if we want this case or not.
#	if (!$self->{fig}->is_genome($genome))
#	{
#	    next;
#	}
	
        print $fh "$genome\t$vc";

        for my $entry (@$row)
        {
            my(@p);

            for my $peg (@$entry)
            {
		$genome =~ /^(\d+\.\d+)/;
		my $just_genome = $1;
                if ($peg =~ /fig\|$just_genome\.(([a-zA-Z]+)\.(\d+))$/)
                {
                    push(@p, ($2 eq "peg") ? $3 : $1);
                }
                else
                {
                    warn "Bad peg $peg in cell for $genome";
                }
            }
            print $fh "\t", join(",", @p);
        }
        print $fh "\n";
    }
}

=head3 get_genomes

    my @genomeList = $sub->get_genomes();

Return a list of the genome IDs for this subsystem. Each genome corresponds to a row
in the subsystem spreadsheet. Indexing into this list returns the ID of the genome
in the specified row.

=cut

sub get_genomes
{
    my($self) = @_;

    my $glist = $self->{genome};

    return ref($glist) ? @$glist : ();
}

=head3 get_variant_codes

    my @codes = $sub->get_variant_codes();

Return a list of the variant codes for each genome, in row index order. The variant
code indicates which variation of the subsystem is used by the given genome.

=cut

sub get_variant_codes
{
    my($self) = @_;

    my $glist = $self->{variant_code};

    return @$glist;
}

=head3 get_variant_code

    my $code = $sub->get_variant_code($gidx);

Return the variant code for the specified genome. Each subsystem has multiple
variants which involve slightly different chemical reactions, and each variant
has an associated variant code. When a genome is connected to the spreadsheet,
the subsystem variant used by the genome must be specified.

=over 4

=item gidx

Row index for the genome whose variant code is desired.

=item RETURN

Returns the variant code for the specified genome.

=back

=cut

sub get_variant_code
{
    my($self, $gidx) = @_;
    my $c = $self->{variant_code}->[$gidx];
    $c =~ s/ //g;
    return $c;
}

sub set_variant_code
{
    my($self, $gidx, $val) = @_;
    $self->{variant_code}->[$gidx] = $val;
    #
    # Update the index for all the pegs in this row.
    # (only if we have a new database)
    #

    if ($self->{old_database})
    {
	return;
    }

    my $rdbH = $self->{fig}->db_handle();
    my $dbh = $rdbH->{_dbh};
    my $cells = $self->get_row($gidx);
    my $sub_name = $self->{name};

    my $sth = $dbh->prepare(qq(UPDATE subsystem_index
			       SET variant = ?
			       WHERE (subsystem = ? AND
				      role = ? AND
				      protein = ?)
			      ));
    for my $i (0 .. $#$cells)
    {
	my $cell = $cells->[$i];
	my $role = $self->get_role($i);

	for my $peg (@$cell)
	{
	    $sth->execute($val, $sub_name, $role, $peg);
	    #warn "Update variant $sub_name $role $peg v='$val'\n";
	}
    }

    return;
}

sub get_variant_code_for_genome
{
    my($self, $genome) = @_;
    my $index = $self->{genome_index}->{$genome};
    if (defined $index) {
     return $self->{variant_code}->[$index];
    }
    else {
     return undef;
    }
}

=head3 get_roles

    my @roles = $sub->get_roles();

Return a list of the subsystem's roles. Each role corresponds to a column
in the subsystem spreadsheet. The list entry at a specified position in
the list will contain the ID of that column's role.

=cut

sub get_roles
{
    my($self) = @_;

    my $rlist = $self->{roles};

    return ref($rlist) ? @$rlist : ();
}

sub get_abbrs
{
    my($self) = @_;

    my $rlist = $self->{role_abbrs};

    return ref($rlist) ? @$rlist : ();
}

=head3 get_abbr_for_role

    my $abbr = $sub->get_abbr_for_role($name);

Return the abbreviation for the given role name.

=cut

sub get_abbr_for_role
{
    my($self, $name) = @_;
    my $idx = $self->{role_index}->{$name};
    if (defined($idx))
    {
	return $self->{role_abbrs}->[$idx];
    }
    else
    {
	return undef;
    }
}


=head3 set_abbr_for_role

    my $abbr = $sub->set_abbr_for_role($name, $abbr);

Set the abbreviation for the given role name.
Return the abbreviation for the given role name.

=cut

sub set_abbr_for_role
{
    my($self, $name, $abbr) = @_;
    my $idx = $self->{role_index}->{$name};
    if (defined($idx))
    {
	$self->{role_abbrs}->[$idx] = $abbr;
	return  $self->{role_abbrs}->[$idx];
    }
    else
    {
	return undef;
    }
}




=head3 get_roles_for_genome

    my $abbr = $sub->get_roles_for_genome($genome_id);

Return the list of roles for which the given genome has nonempty cells.

=cut

sub get_roles_for_genome
{
    my($self, $genome) = @_;

    my $gidx = $self->{genome_index}->{$genome};
    return undef unless defined($gidx);

    my $row = $self->{spreadsheet}->[$gidx];

    my @out;
    for my $ridx (0 .. $#$row)
    {
	my $cell = $row->[$ridx];
	if ($cell && @$cell > 0)
	{
	    push(@out, $self->{roles}->[$ridx]);
	}
    }
    return @out;
}

sub roles_with_abbreviations
{
    my($self) = @_;

    my @ret;
    return @ret unless (defined $self->{roles}); # return an empty list if we have no roles!

    for my $i (0..@{$self->{roles}} - 1)
    {
        push(@ret, [$self->{role_abbrs}->[$i], $self->{roles}->[$i]]);
    }
    return @ret;
}


sub get_sorted_rows
{
    my($self, $sort_order) = @_;

    my $fig = $self->{fig};

    my @rows;
    for (my $i = 0; $i < @{$self->{genome}}; $i++)
    {
        my $gid = $self->{genome}->[$i];
	$gid =~ s/:.*$//;

        my $gs = $fig->genus_species($gid);

        my $q = quotemeta($gid);
        my $cells = [];
        for my $c (@{$self->{spreadsheet}->[$i]})
        {
            push(@$cells, [map { s/^fig\|$q\.peg\.//; $_ } @$c]);
        }

        push(@rows, [$self->{genome}->[$i], $gs, $self->{variant_code}->[$i], $cells]);
    }

    if ($sort_order eq "by_phylo")
    {
        return(map  { $_->[0] }
               sort { ($a->[1] cmp $b->[1]) or ($a->[0]->[1] cmp $b->[0]->[1]) }
               map  { [$_, $fig->taxonomy_of($_->[0]) ] } @rows);
    }
    elsif ($sort_order eq "alphabetic")
    {
        return sort { ($a->[1] cmp $b->[1]) or ($a->[0] <=> $b->[0]) } @rows;
    }
    elsif ($sort_order eq "by_tax_id")
    {
        return sort { $a->[0] <=> $b->[0] } @rows;
    }
    else
    {
        return @rows;
    }
}


sub get_row
{
    my($self, $row) = @_;

    return $self->{spreadsheet}->[$row];
}

sub get_col
{
    my($self, $col) = @_;

    return $self->{spreadsheet_inv}->[$col];
}

sub get_cell
{
    my($self, $row, $col) = @_;

    my $cell = $self->{spreadsheet}->[$row]->[$col];
    if (! defined($cell))
    {
        $cell = $self->{spreadsheet}->[$row]->[$col] = [];
    }
    return $cell;
}

=head3 get_genome_index

    my $idx = $sub->get_genome_index($genome);

Return the row index for the genome with the specified ID.

=over 4

=item genome

ID of the genome whose row index is desired.

=item RETURN

Returns the row index for the genome with the specified ID, or an undefined
value if the genome does not participate in the subsystem.

=back

=cut

sub get_genome_index
{
    my($self, $genome) = @_;

    return $self->{genome_index}->{$genome};
}

sub get_genome
{
    my($self, $gidx) = @_;

    return $self->{genome}->[$gidx];
}

=head3 get_role_index

    my $idx = $sub->get_role_index($role);

Return the column index for the role with the specified ID.

=over 4

=item role

ID (full name) of the role whose column index is desired.

=item RETURN

Returns the column index for the role with the specified name.

=back

=cut

sub get_role_index
{
    my($self, $role) = @_;

    return $self->{role_index}->{$role};
}

sub get_role
{
    my($self, $ridx) = @_;

    return $self->{roles}->[$ridx];
}

=head3 get_role_abbr

    my $abbr = $sub->get_role_abbr($ridx);

Return the abbreviation for the role in the specified column. The abbreviation
is a shortened identifier that is not necessarily unique, but is more likely to
fit in a column heading.

=over 4

=item ridx

Column index for the role whose abbreviation is desired.

=item RETURN

Returns an abbreviated name for the role corresponding to the indexed column.

=back

=cut

sub get_role_abbr
{
    my($self, $ridx) = @_;

    return $self->{role_abbrs}->[$ridx];
}

sub get_role_from_abbr
{
    my($self, $abbr) = @_;

    return $self->{abbr}->{$abbr};
}

=head3 set_pegs_in_cell

    $sub->set_pegs_in_cell($genome, $role, $peg_list);

Set the cell for the given genome and role to $peg_list.

=cut

sub set_pegs_in_cell
{
    my($self, $genome, $role, $peg_list) = @_;
    my($row, $col);

    #
    # If row isn't numeric, look it up in the genomes list.
    #

    if ($genome !~ /^\d+$/)
    {
        $row = $self->{genome_index}->{$genome};
    }
    else
    {
        $row = $genome
    }

    if (! defined($row))
    {
#        print &Dumper($self->{genome_index});
        confess "Cannot find row for $genome\n";
        return undef;
    }

    #
    # If col isn't numeric, look it up in the roles and role abbreviations.
    #

    my $role_name;
    if ($role !~ /^\d+$/)
    {
        #
        # See if it's an abbr
        #

        my $a = $self->{abbr}->{$role};
        $role = $a if $a;

        $col = $self->{role_index}->{$role};
	$role_name = $role;
    }
    else
    {
        $col = $role;
	$role_name = $self->get_role($col);
    }

    if (! defined($col))
    {
        print STDERR &Dumper($self->{role_index});
        confess "Cannot find col for $role\n";
        return undef;
    }
    my $cell = $self->get_cell($row, $col);


    if (defined($cell))
    {
	my $sub_name = $self->{name};
	my $peg;
	my $rdbH = $self->{fig}->db_handle();
	my $dbh = $rdbH->{_dbh};

	my $variant = $self->get_variant_code($row);

	if (@$cell > 0)
	{
	    my $sth = $dbh->prepare(qq(DELETE FROM subsystem_index
				       WHERE (subsystem = ? AND
					      role = ? AND
					      protein = ?)
				       ));
	    foreach $peg (@$cell)
	    {
		$sth->execute($sub_name, $role_name, $peg);
		# warn "Deleting $sub_name $role $peg\n";
	    }
	}

	@$cell = @$peg_list;

	if ($self->{old_database})
	{
	    my $sth = $rdbH->{_dbh}->prepare(qq(INSERT INTO subsystem_index (protein,subsystem,role)
						VALUES (?, ?, ?)));
	    foreach $peg (@$cell)
	    {
		$sth->execute($peg, $sub_name, $role_name);
		# warn "Add old $peg $sub_name $role\n";
	    }
	}
	else
	{
	    my $sth = $rdbH->{_dbh}->prepare(qq(INSERT INTO subsystem_index (protein,subsystem,role,variant)
						VALUES (?, ?, ?, ?)));
	    foreach $peg (@$cell)
	    {
		$sth->execute($peg, $sub_name, $role_name, $variant);
		#warn "Add new $peg $sub_name $role_name v='$variant'\n";
	    }
	}
    }
    else
    {
        warn "set_pegs_in_cell: Could not find cell!";
    }
}

=head3 get_pegs_from_cell

    my @pegs = $sub->get_pegs_from_cell($rowstr, $colstr);

Return a list of the peg IDs for the features in the specified spreadsheet cell.

=over 4

=item rowstr

Genome row, specified either as a row index or a genome ID.

=item colstr

Role column, specified either as a column index, a role name, or a role
abbreviation.

=item RETURN

Returns a list of PEG IDs. The PEGs in the list belong to the genome in the
specified row and perform the role in the specified column. If the indicated
row and column does not exist, returns an empty list.

=back

=cut

sub get_pegs_from_cell
{
    my($self, $rowstr, $colstr) = @_;
    my($row, $col);

    #
    # If row isn't numeric, look it up in the genomes list.
    #

    if ($rowstr !~ /^\d+$/)
    {
        $row = $self->{genome_index}->{$rowstr};
    }
    else
    {
        $row = $rowstr;
    }

    if (! defined($row))
    {
        print STDERR &Dumper($self->{genome_index});
        confess "Cannot find row for $rowstr\n";
        return undef;
    }

    #
    # If col isn't numeric, look it up in the roles and role abbreviations.
    #

    if ($colstr !~ /^\d+$/)
    {
        #
        # See if it's an abbr
        #

        my $a = $self->{abbr}->{$colstr};
        $colstr = $a if $a;

        $col = $self->{role_index}->{$colstr};
    }
    else
    {
        $col = $colstr;
    }

    if (! defined($col))
    {
        warn "Cannot find col for $colstr\n";
        return undef;
    }
    my $cell = $self->get_cell($row, $col);

    if ($cell)
    {
        return @$cell;
    }
    else
    {
        return undef;
    }
}

#
# Subset support
#

sub get_active_subsetC
{
    my($self) = @_;

    return $self->{col_active_subset};
}

sub get_active_subsetR
{
    my($self) = @_;

    return $self->{row_active_subset};
}

sub set_active_subsetC
{
    my($self, $subset) = @_;

    $self->{col_active_subset} = $subset;
}


sub set_active_subsetR
{
    my($self, $subset) = @_;

    $self->{row_active_subset} = $subset;
}

# This method is deprecated. Use get_subset_namesC.
sub get_subset_names
{
    my($self) = @_;

    return $self->get_subset_namesC;
}

=head3 get_subset_namesC

    my @subsetNames = $sub->get_subset_namesC();

Return a list of the names for all the column (role) subsets. Given a subset
name, you can use the L</get_subsetC_roles> method to get the roles in the
subset.

=cut

sub get_subset_namesC
{
    my($self) = @_;
    my $col_subsets = $self->{col_subsets};
    return ("All",defined($col_subsets) ? @$col_subsets : ());
#   return ("All",@{$self->{col_subsets}});
}

sub get_subset_namesR
{
    my($self) = @_;

    return ("All",@{$self->{row_subsets}});
}

=head3 get_subsetC_roles

    my @roles = $sub->get_subsetC_roles($subname);

Return the names of the roles contained in the specified role (column) subset.

=over 4

=item subname

Name of the role subset whose roles are desired.

=item RETURN

Returns a list of the role names for the columns in the named subset.

=back

=cut

sub get_subsetC_roles
{
    my($self, $subname) = @_;
    return map { $self->get_role($_) } $self->get_subsetC($subname);
}

sub get_subsetC
{
    my($self, $subname) = @_;
    if ($subname eq "All") { return map { $self->get_role_index($_) } $self->get_roles }

    if (!defined($self->{col_subset_members}->{$subname}))
    {
        $self->{col_subset_members}->{$subname} = [];
    }

    return @{$self->{col_subset_members}->{$subname}};
}

sub get_subset
{
    my($self, $subname) = @_;
    return $self->get_subsetC($subname);
}

=head3 get_subsetR

    my @genomes = $sub->get_subsetR($subName);

Return the genomes in the row subset indicated by the specified subset name.

=over 4

=item subName

Name of the desired row subset, or C<All> to get all of the rows.

=item RETURN

Returns a list of genome IDs corresponding to the named subset.

=back

=cut

sub get_subsetR {
    my($self, $subname) = @_;
    my($pair,$id,$members,$genome);

    if ($subname eq "All") { return $self->get_genomes }
    my %genomes = map { $_ => 1 } $self->get_genomes;

    return grep { $genomes{$_} } @{$self->{row_subset_members}->{$subname}};
}

sub load_row_subsets {
    my($self) = @_;
    my($id,$members,$pair);

    my $taxonomic_groups = $self->{fig}->taxonomic_groups_of_complete(10);
    foreach $pair (@$taxonomic_groups)
    {
        ($id,$members) = @$pair;
        if ($id ne "All")
        {
            push(@{$self->{row_subsets}},$id);
        }
        $self->{row_subset_members}->{$id} = $members;
    }
}

=head3 load_row_subsets_by_kv

Load a row subset based on a key/value pair. This will take a single key/value pair and only show that subset

It is just a modification of load_row_subsets to deal with kv pairs

This takes a required argument: the key that the genome must have, and a second optional argument, the value that key must hold.

=cut

sub load_row_subsets_by_kv {
 my ($self, $key, $want) = @_;
 my($id,$members,$pair);
 my $keep;
 #
 # First do a single call to retrieve all the values for the subset key.
 #
 my @attr_values = ();
 if ( $key eq 'hundred_hundred' ) {
   @attr_values = $self->{fig}->get_attributes($self->{genome}, 'collection', $key);
 }
 else {
   @attr_values = $self->{fig}->get_attributes($self->{genome}, $key);
 }

 my %amap;
 map {  push(@{$amap{$_->[0]}}, [@$_]); } @attr_values;

 foreach my $genome (@{$self->{genome}}) {
  #my @results=$self->{fig}->get_attributes($genome, $key);
  my $results = $amap{$genome};
  next unless $results;
  foreach my $res (@$results) {
   my ($gotid, $gottag, $value, $url)=@$res;
   next if ($value && $want && $value ne $want);
   next if ($gotid ne $genome);
   push @$keep, $genome;
   last;
  }
 }
 $self->{row_subset_members}->{$key}=$keep;
}

=head3 set_subsetC

    $sub->set_subsetC($name, $members);

Create a subset with the given name and members.

$members is a list of role names.

=cut

sub set_subsetC
{
    my($self, $subname, $list) = @_;

    my $nl = [map { $self->get_role_index($_) } @$list];

    $self->_set_subset($subname, $nl);
}

sub set_subset
{
    my($self, $subname, $list) = @_;

    $self->set_subsetsC($subname,$list);
}

=head3 _set_subset

Create a subset with the given name and members.

Internal version  - here, members is a list of role indices.

=cut

sub _set_subset
{
    my($self, $subname, $list) = @_;
    $self->{col_subset_members}->{$subname} = $list;
    my($i,$x);
    $x = $self->{col_subsets};
    for ($i=0; ($i < @$x) && ($x->[$i] ne $subname); $i++) {}
    if ($i == @$x)
    {
        push(@$x,$subname);
    }
}

sub delete_subsetC
{
    my($self, $subname) = @_;
    my($i,$x);

    $x = $self->{col_subsets};
    for ($i=0; ($i < @$x) && ($x->[$i] ne $subname); $i++) {}
    if ($i < @$x)
    {
        splice(@$x,$i,1);
    }
    delete $self->{col_subset_members}->{$subname};
}

#
# Role manipulation.
#


=head3 set_roles

    $sub->set_roles($role_list);

Set the list of roles. C<$role_list> is a list of tuples C<[$role_name, $abbreviation]>.

If a role already exists, it is used. If it does not exist, it is created empty.

=cut

sub set_roles
{
    my($self, $roles) = @_;

    #
    # We do this by first creating a new spreadsheet.
    #
    # It is easiest to do this by manipulating the inverted spreadsheet
    # (role-major), and then creating the non-inverted spreadsheet from it.
    #

    my $oldss = $self->{spreadsheet};
    my $oldssinv = $self->{spreadsheet_inv};

    my $ss = [];
    my $ssinv = [];

    my $g = $self->{genome};
    my $ng = @$g;

    my $old_roles = $self->{role_index};

    my @role_index_conversion;
    my @old_role_list = @{$self->{roles}};

    #
    # Since we're setting up completely new roles, wipe the
    # existing state.
    #

    $self->{abbr} = {};
    $self->{role_index} = {};
    $self->{roles} = [];
    $self->{role_abbrs} = [];

    #
    # Initialize %defunct_roles with the list of all roles.
    # Remove entries as we walk the list of new roles below.
    # Any that are remaining need to be pulled from the index.
    #

    my %defunct_roles = map { $_ => 1 } @old_role_list;

    # warn "Defunct at start: ", Dumper(\%defunct_roles);
    for (my $idx = 0; $idx < @$roles; $idx++)
    {
        my $role = $roles->[$idx]->[0];
        my $abbr = $roles->[$idx]->[1];

        my $old_idx = $old_roles->{$role};

	if (defined($old_idx))
	{
	    # warn "Found old idx $old_idx for $role $idx\n";
	    # warn $oldssinv->[$old_idx];
	    $ssinv->[$idx] = $oldssinv->[$old_idx];

	    $role_index_conversion[$old_idx] = $idx;

	    #
	    # We're keeping it, so it's not defunct anymore.
	    #
	    delete $defunct_roles{$role};
	}
	else
	{
	    # warn "Did not find old role for $role $idx\n";
	    # warn Dumper($old_roles);
	    my $l = [];
	    for (my $j = 0; $j < $ng; $j++)
	    {
		$l->[$j] = [];
	    }

	    $ssinv->[$idx] = $l;
	}


	#
	# While we're here, update the new role and abbrev indexes
	#
	$self->{role_index}->{$role} = $idx;
	$self->{abbr}->{$abbr} = $role;
	$self->{roles}->[$idx] = $role;
	$self->{role_abbrs}->[$idx] = $abbr;
    }

    #
    # Now we delete the pegs showing up for the list of defunct roles.
    #
    # warn "Defunct at finish: ", Dumper(\%defunct_roles);

    my $rdbH = $self->{fig}->db_handle();
    my $dbh = $rdbH->{_dbh};
    my $sub_name = $self->{name};

    my $sth = $dbh->prepare(qq(DELETE FROM subsystem_index
			       WHERE (subsystem = ? AND
				      role = ? AND
				      protein = ?)
			      ));


    for my $defunct_role (keys(%defunct_roles))
    {
	my $defunct_role_idx = $old_roles->{$defunct_role};
	my $col = $oldssinv->[$defunct_role_idx];
	# warn "Handle defunct role $defunct_role idx=$defunct_role_idx\n", Dumper($col);

	for my $cell (@$col)
	{
	    for my $peg (@$cell)
	    {
		$sth->execute($sub_name, $defunct_role, $peg);
		warn "Deleting $sub_name $defunct_role $peg\n";
	    }
	}
    }


    #
    # Now create the uninverted spreadsheet.
    #

    for (my $gidx = 0; $gidx < $ng; $gidx++)
    {
        my $row = [];
        $ss->[$gidx] = $row;
        for (my $ridx = 0; $ridx < @$roles; $ridx++)
        {
            $row->[$ridx] = $ssinv->[$ridx]->[$gidx];
        }
    }

    $self->{spreadsheet} = $ss;
    $self->{spreadsheet_inv} = $ssinv;

    #
    # Fix up the subsets.
    #


    for my $subset (grep { $_ ne "All" } $self->get_subset_names())
    {
        my $n = [];
        for my $idx ($self->get_subset($subset))
        {
            my $new = $role_index_conversion[$idx];
            if (defined($new))
            {
                push(@$n, $new);
            }
        }
        $self->_set_subset($subset, $n);
    }

}



=head3 add_role($role, $abbr)

Add the given role to the spreadsheet.

This causes a new column to be added, with empty values in each cell.

We do nothing if the role is already present.

Return the index of the new role.

=cut

sub add_role
{
    my($self, $role, $abbr) = @_;

    if (defined($self->get_role_index($role)))
    {
        warn "Role $role already present\n";
        return undef;
    }

    #
    # Add to the roles list. It goes at the end.
    #

    if (!defined $self->{roles}) {$self->{roles}=[]}

    my $idx = @{$self->{roles}};
    $self->{roles}->[$idx] = $role;
    $self->{role_abbrs}->[$idx] = $abbr;
    $self->{role_index}->{$role} = $idx;
    $self->{abbr}->{$abbr} = $role;

    #
    # Update the spreadsheet.
    # On the standard one, we have to go through all the rows adding
    # a columnt to each.
    #
    # On the inverted one, we add a column with [] in each entry.
    #

    if (!defined $self->{genome}) {$self->{genome}=[]}

    my $ng = @{$self->{genome}};
    my $newcol = [];

    for (my $i = 0; $i < $ng; $i++)
    {
        my $cell = [];
        # print "nr: Adding cell $cell for gidx=$i ridx=$idx\n";
        $self->{spreadsheet}->[$i]->[$idx] = $cell;
        $newcol->[$i] = $cell;
    }

    $self->{spreadsheet_inv}->[$idx] = $newcol;

    return $idx;
}


=head3 change_role( $oldrole, $newrole )

Change just the function of a role

=cut

sub change_role
{
    my( $self, $oldrole, $newrole) = @_;

    my $oldindex = $self->get_role_index( $oldrole );
    unless ( defined( $oldindex ) ) {
      return ( 0, "The role $oldrole does not exist in this subsystem.<BR>\n" );
    }

    my $abbr = $self->{role_abbrs}->[$oldindex];

    $self->{roles}->[$oldindex] = $newrole;
    delete $self->{role_index}->{$oldrole};
    $self->{role_index}->{$newrole} = $oldindex;
    $self->{abbr}->{$abbr} = $newrole;

    for my $key (qw(reactions hope_reactions hope_reaction_notes hope_reaction_links))
    {
	my $val = delete $self->{$key}->{$oldrole};
	if ($val)
	{
	    $self->{$key}->{$newrole} = $val;
	}
    }
					 
    return ( 1 );
}

=head3 change_role_in_column_of_file($oldrole, $newrole, $colnum, $file, $backup_file)

Copy $file to $backup_file, then rewrite $file mapping any occurences of $oldrole with $newrole
in column $colnum. Assume file is tab-separated.

Used for remapping the various reactions files.

$file and $backup_file are paths relative to the subsystem directory.

=cut

sub change_role_in_column_of_file
{
    my($self, $oldrole, $newrole, $colnum, $file, $backup_file) = @_;

    my $dir = $self->{dir};

    if ($file !~ m,^/,)
    {
	$file = "$dir/$file";
    }

    if ($backup_file !~ m,^/,)
    {
	$backup_file = "$dir/$backup_file";
    }
    
    my $rc = system("cp", $file, $backup_file);
    if ($rc != 0)
    {
	die "Subsystem::change_role_in_column_of_file: cp $file $backup_file failed with rc $rc";
    }

    my($from, $to);
    open($from, "<", $backup_file) or die "cannot open $backup_file for reading: $!";
    open($to, ">", $file) or die "cannot open $file for writing: $!";

    while (<$from>)
    {
	chomp;
	my @a = split(/\t/);
	if ($colnum < @a and $a[$colnum] eq $oldrole)
	{
	    $a[$colnum] = $newrole;
	}
	print $to join("\t", @a), "\n";
    }

    close($from);
    close($to);
}

=head3 remove_role

Remove the role from the spreadsheet.

We do nothing if the role is not present.

=cut

sub remove_role
{
    my($self, $role) = @_;

    my $idx = $self->get_role_index($role);
    if (!defined($idx))
    {
        warn "Role $role not present\n";
        return undef;
    }

    #
    # Update the index. Again, do this before removing roles
    # so we have full data to work with.
    # We walk the role's column of the spreadsheet removing pegs from the index.
    #

    my $rdbH = $self->{fig}->db_handle();
    my $dbh = $rdbH->{_dbh};
    my $sub_name = $self->{name};

    my $sth = $dbh->prepare(qq(DELETE FROM subsystem_index
			       WHERE (subsystem = ? AND
				      role = ? AND
				      protein = ?)
			      ));
    my $col = $self->get_col($idx);
    for my $cell (@$col)
    {
	 for my $peg (@$cell)
	 {
	    $sth->execute($sub_name, $role, $peg);
	    warn "Deleting $sub_name $role $peg\n";
	 }
    }

    #
    # Remove from the roles list.
    #

    my $abbr = $self->{role_abbrs}->[$idx];

    splice(@{$self->{roles}}, $idx, 1);
    splice(@{$self->{role_abbrs}}, $idx, 1);
    delete $self->{role_index}->{$role};
    delete $self->{abbr}->{$abbr};


    #
    # Update the spreadsheet.
    # On the standard one, we have to go through all the rows removing
    # the column from each.
    #
    # On the inverted one, we just remove the column.
    #

    my $ng = @{$self->{genome}};
    my $newcol = [];

    for (my $i = 0; $i < $ng; $i++)
    {
        splice(@{$self->{spreadsheet}->[$i]}, $idx, 1);
    }

    splice(@{$self->{spreadsheet_inv}}, $idx, 1);

    #
    # We need to rewrite the subsets. if $idx was present in one, it is
    # removed. Any index >$idx is decremented.
    #

    for my $subset ($self->get_subset_names())
    {
        my @n;

        for my $sidx ($self->get_subset($subset))
        {
            if ($sidx < $idx)
            {
                push(@n, $sidx);
            }
            elsif ($sidx > $idx)
            {
                push(@n, $sidx - 1);
            }
        }

        $self->_set_subset($subset, [@n]);
    }
}

=head3 add_genome($genome, $abbr)

Add the given genome to the spreadsheet.

This causes a new row to be added, with empty values in each cell.

We do nothing if the genome is already present.

Return the index of the new genome.

=cut

sub add_genome
{
    my($self, $genome) = @_;

    my $idx = $self->get_genome_index($genome);
    if ( defined( $idx ) ) {
      warn "Genome $genome already present\n";
      return $idx;
    }

    #
    # Add to the genomes list. It goes at the end.
    #

    if (!defined $self->{genome}) {$self->{genome} = []}

    $idx = @{$self->{genome}};
    $self->{variant_code}->[$idx] = 0;
    $self->{genome}->[$idx] = $genome;
    $self->{genome_index}->{$genome} = $idx;

    #
    # Update the spreadsheet.
    # On the inverted one, we have to go through all the columns adding
    # a row to each.
    #
    # On the regular one, we add a row with [] in each entry.
    #

    if (!defined $self->{roles}) {$self->{roles} = []}
    my $nr = @{$self->{roles}};
    my $newrow = [];

    for my $i (0.. $nr - 1)
    {
        my $cell = [];
        # print "ng: Adding cell $cell for gidx=$idx ridx=$i\n";
        $self->{spreadsheet_inv}->[$i]->[$idx] = $cell;
        $newrow->[$i] = $cell;
    }

    $self->{spreadsheet}->[$idx] = $newrow;

    return $idx;
}

=head3 remove_genome

Remove the genome from the spreadsheet.

We do nothing if the genome is not present.

=cut

sub remove_genome
{
    my($self, $genome) = @_;

    my $idx = $self->get_genome_index($genome);
    if (!defined($idx))
    {
        warn "Genome $genome not present\n";
        return undef;
    }

    #
    # Remove from database index (before we delete stuff from here,
    # so we have access to th e data structures).
    #

    my $rdbH = $self->{fig}->db_handle();
    my $dbh = $rdbH->{_dbh};
    my $cells = $self->get_row($idx);
    my $sub_name = $self->{name};

    my $sth = $dbh->prepare(qq(DELETE FROM subsystem_index
			       WHERE (subsystem = ? AND
				      role = ? AND
				      protein = ?)
			      ));
    for my $i (0 .. $#$cells)
    {
	my $cell = $cells->[$i];
	my $role = $self->get_role($i);

	for my $peg (@$cell)
	{
	    $sth->execute($sub_name, $role, $peg);
	    warn "Deleting $sub_name $role $peg\n";
	}
    }

    #
    # Remove from the genomes list.
    #

    splice(@{$self->{genome}}, $idx, 1);

    my $genome1;
    foreach $genome1 (@{$self->{genome}})
    {
        if ($self->{genome_index}->{$genome1} > $idx)
        {
            $self->{genome_index}->{$genome1}--;
        }
    }
    splice(@{$self->{variant_code}}, $idx, 1);

    delete $self->{genome_index}->{$genome};

    #
    # Update the spreadsheet.
    # On the inverted one, we have to go through all the columns removing
    # the row from each.
    #
    # On the standard one, we just remove the row.
    #

    my $nr = @{$self->{roles}};

    for my $i (0 .. $nr - 1)
    {
        splice(@{$self->{spreadsheet_inv}->[$i]}, $idx, 1);
    }

    splice(@{$self->{spreadsheet}}, $idx, 1);

}

sub get_name
{
    my($self) = @_;
    my $name = $self->{name};
    $name =~ s/ /_/g;
    return $name;
}

sub get_dir
{
    my($self) = @_;
    return $self->{dir};
}


sub get_version
{
    my($self) = @_;
    return $self->{version};
}

=head3 get_notes

    my $text = $sub->get_notes();

Return the descriptive notes for this subsystem.

=cut

sub get_notes
{
    my($self) = @_;

    return $self->{notes};
}

=head3 get_description

    my $text = $sub->get_description();

Return the description for this subsystem.

=cut

sub get_description
{
    my($self) = @_;

    return $self->{description};
}

=head3 get_variants

    my $text = $sub->get_variants();

Return the variants for this subsystem.

=cut

sub get_variants
{
    my($self) = @_;

    my $text = $self->{variants};
    my %vars;

    my @lines = split( "\n", $text );
    foreach ( @lines ) {
      my ( $v, $d ) = split( "\t", $_ );
      $vars{ $v } = $d;
    }

    return \%vars;
}

=head3 get_literature

    my $text = $sub->get_literature();

Return the literature for this subsystem.

=cut

sub get_literature
{
    my($self) = @_;

    return $self->{literature};
}

=head3 get_reactions

    my $reactHash = $sub->get_reactions();

Return a reference to a hash that maps each role ID to a list of the reactions
catalyzed by the role.

=cut

sub get_reactions
{
    my($self) = @_;

    return $self->{reactions};
}

sub set_reaction {
    my($self,$role,$rstring) = @_;

    $self->{reactions}->{$role} = [split(/,\s*/,$rstring)];
}

sub get_hope_scenario_names
{
    my($self) = @_;

    return sort keys %{$self->{hope_scenarios}};
}

sub change_hope_scenario_name
{
    my($self, $old_name, $new_name) = @_;
    my $hope_scenarios = $self->{hope_scenarios};
    my $scenario = $hope_scenarios->{$old_name};
    delete $hope_scenarios->{$old_name};
    $hope_scenarios->{$new_name} = $scenario;
}

sub add_hope_scenario
{
    my($self, $new_name) = @_;
    $new_name =~ s/\// /g;
    $self->{hope_scenarios}->{$new_name} = { input_compounds => [], output_compounds => [], map_ids => [], additional_reactions => [], ignore_reactions => [] };
}

sub delete_hope_scenario
{
    my($self, $scenario_name) = @_;

    delete $self->{hope_scenarios}->{$scenario_name};
}

sub get_hope_input_compounds
{
	my($self, $scenario_name) = @_;

	return @{$self->{hope_scenarios}->{$scenario_name}->{input_compounds}};
}

sub set_hope_input_compounds
{
	my($self,$scenario_name,$compounds) = @_;
	$compounds =~ s/^\s+//g;
	$compounds =~ s/\s+$//g;
	$compounds =~ s/,\s+/,/g;
	$compounds =~ s/\s+,/,/g;
	$compounds =~ s/\s+/,/g;
	$self->{hope_scenarios}->{$scenario_name}->{input_compounds} = [split(/,/,$compounds)];
}

sub get_hope_output_compounds
{
	my($self,$scenario_name) = @_;

	return @{$self->{hope_scenarios}->{$scenario_name}->{output_compounds}};
}

sub set_hope_output_compounds
{
	my($self,$scenario_name,$compounds) = @_;
	$compounds =~ s/^\s+//g;
	$compounds =~ s/\s+$//g;
	$compounds =~ s/,\s+/,/g;
	$compounds =~ s/\s+,/,/g;
	$compounds =~ s/\s+/,/g;

	# allow one level of nesting with parentheses
	my @output_compounds_lists;
	my @inner_list;

	foreach my $cpd (split(/,/,$compounds))
	{
	    if ($cpd =~ /\(/)
	    {
		$cpd =~ s/\(//g;
		push @inner_list, $cpd;
	    }
	    elsif (scalar @inner_list > 0 && $cpd =~ /\)/)
	    {
		$cpd =~ s/\)//g;
		push @inner_list, $cpd;
		my @new_inner_list = @inner_list;
		push @output_compounds_lists, \@new_inner_list;
		@inner_list = ();
	    }
	    elsif (scalar @inner_list > 0)
	    {
		push @inner_list, $cpd;
	    }
	    else
	    {
		push @output_compounds_lists, [$cpd];
	    }
	}

	$self->{hope_scenarios}->{$scenario_name}->{output_compounds} = \@output_compounds_lists;
}

sub get_hope_map_ids
{
	my($self,$scenario_name) = @_;

	return @{$self->{hope_scenarios}->{$scenario_name}->{map_ids}};
}

sub set_hope_map_ids
{
	my($self,$scenario_name,$ids) = @_;
	$ids =~ s/^\s+//g;
	$ids =~ s/\s+$//g;
	$ids =~ s/,\s+/,/g;
	$ids =~ s/\s+,/,/g;
	$ids =~ s/\s+/,/g;
	$self->{hope_scenarios}->{$scenario_name}->{map_ids} = [split(/,/,$ids)];
}

sub get_hope_additional_reactions
{
	my($self,$scenario_name) = @_;

	return @{$self->{hope_scenarios}->{$scenario_name}->{additional_reactions}};
}

sub set_hope_additional_reactions
{
	my($self,$scenario_name,$rids) = @_;
	$rids =~ s/^\s+//g;
	$rids =~ s/\s+$//g;
	$rids =~ s/,\s+/,/g;
	$rids =~ s/\s+,/,/g;
	$rids =~ s/\s+/,/g;
	$self->{hope_scenarios}->{$scenario_name}->{additional_reactions} = [split(/,/,$rids)];
}

sub get_hope_ignore_reactions
{
	my($self,$scenario_name) = @_;

	return @{$self->{hope_scenarios}->{$scenario_name}->{ignore_reactions}};
}

sub set_hope_ignore_reactions
{
	my($self,$scenario_name,$rids) = @_;
	$rids =~ s/^\s+//g;
	$rids =~ s/\s+$//g;
	$rids =~ s/,\s+/,/g;
	$rids =~ s/\s+,/,/g;
	$rids =~ s/\s+/,/g;
	$self->{hope_scenarios}->{$scenario_name}->{ignore_reactions} = [split(/,/,$rids)];
}

sub get_hope_reactions
{
    my($self) = @_;

    return %{$self->{hope_reactions}};
}

sub get_emptycells
{
    my($self) = @_;

    return $self->{emptycells};
}

sub set_hope_reaction {
    my($self,$role,$rids) = @_;
    $rids =~ s/,\s+/,/g;
    $rids =~ s/\s+/,/g;
    $self->{hope_reactions}->{$role} = [split(/,/,$rids)];
}

sub get_hope_reaction_notes
{
    my($self) = @_;

    return %{$self->{hope_reaction_notes}};
}

sub set_hope_reaction_note {
    my($self,$role,$rstring) = @_;

    $self->{hope_reaction_notes}->{$role} = $rstring;
}

sub get_hope_reaction_links
{
    my($self) = @_;

    return %{$self->{hope_reaction_links}};
}

sub set_hope_reaction_link {
    my($self,$role,$rstring) = @_;

    $self->{hope_reaction_links}->{$role} = $rstring;
}

sub get_hope_curation_notes
{
    my($self) = @_;

    return $self->{hope_curation_notes};
}

sub set_hope_curation_notes
{
    my($self, $hope_curation_notes) = @_;

    $self->{hope_curation_notes} = $hope_curation_notes;
}

sub set_emptycells
{
    my($self, $emptycells) = @_;

    $self->{emptycells} = $emptycells;
}

sub set_notes
{
    my($self, $notes) = @_;

    $self->{notes} = $notes;
}

sub set_description
{
    my($self, $desc) = @_;

    $self->{description} = $desc;
}

sub set_variants
{
    my($self, $var) = @_;

    my $text = '';
    foreach my $k ( sort keys %$var ) {
      $text .= "$k\t".$var->{ $k }."\n";
    }

    $self->{variants} = $text;
}

sub set_literature
{
    my($self, $lit) = @_;

    $self->{literature} = $lit;
}

sub get_classification
{
    my($self) = @_;

    return $self->{classification};
}

sub set_classification
{
    my($self, $classification) = @_;

    $self->{classification}=$classification;
}


=head3 get_curator

    my $userName = $sub->get_curator();

Return the name of this subsystem's official curator.

=cut

sub get_curator
{
    my($self) = @_;
    return $self->{curator};
}

=head3 set_curator

    my $userName = $sub->set_curator("curator");

Sets the name of this subsystem's official curator.
Return the name of this subsystem's official curator.

=cut

sub set_curator
{
    my($self, $curator) = @_;
    $self->{fig}->set_user($curator);
    $self->{curator} = $curator;
    return $self->{curator};
}



sub get_created
{
    my($self) = @_;
    return $self->{created};
}

sub get_last_updated
{
    my($self) = @_;
    return $self->{last_updated};
}

sub get_checkvariant_definitions {
  my ( $self ) = @_;

  my $cvd_string = '';
  if ( open( CVD, "<$self->{dir}/checkvariant_definitions" ) ) {
    while ( defined( $_ = <CVD> ) ) {
      $cvd_string .= $_;
    }
    close( CVD );
  }
  return $cvd_string;
}

sub get_checkvariant_rules {
  my ( $self ) = @_;

  my $cvr_string = '';
  if ( open( CVR, "<$self->{dir}/checkvariant_rules" ) ) {
    while ( defined( $_ = <CVR> ) ) {
      $cvr_string .= $_;
    }
    close( CVR );
  }
  return $cvr_string;
}

sub save_checkvariant_definitions {
  my ( $self, $cvd_string ) = @_;
  open( CVD, ">$self->{dir}/checkvariant_definitions" ) || die "could not open $self->{dir}/checkvariant_definitions file: $!";
  print CVD "$cvd_string\n";
  close( CVD );
}

sub save_checkvariant_rules {
  my ( $self, $cvs_string ) = @_;

  open( CVS, ">$self->{dir}/checkvariant_rules" ) || die "could not open checkvariant_rules file";
  print CVS "$cvs_string\n";
  close( CVS );
}

#
# Subsystem copying logic
#

=head3 add_to_subsystem($subsystem_name, $columns, $notes_flag)

Merge the given columns from $subsystem_name into this subsystem. Append the
notes from the subsystem if $notes_flag is true.

=cut

sub add_to_subsystem
{
    my($self, $subsystem_name, $cols, $add_notes) = @_;

    my $ss = $self->{fig}->get_subsystem($subsystem_name);

    if (!$ss)
    {
        warn "Cannot open subsystem '$subsystem_name' to copy from";
        return;
    }

    #
    # Merge the data from the other subsystem.
    #
    # First we assure ourselves that we have the appropriate roles. While
    # we do this, build the list of row indices (in this subsystem) that
    # map to the roles we are adding.
    #

    #
    # local_roles[$his_role] = $my_role (map from other role idx to local role idx)
    #

    my @local_roles;

    #
    # his_roles = list of role indices corresponding to the remote roles.
    #
    if ($cols->[0] eq "all")
    {
        $cols = [$ss->get_roles];
    }

    my @his_roles;

    for my $his_role (@$cols)
    {
        my $idx = $self->get_role_index($his_role);
        my $his_idx = $ss->get_role_index($his_role);

        if (!defined($his_idx))
        {
            confess "Cannot map his role $his_role\n";
        }
        push(@his_roles, $his_idx);

        if (!defined($idx))
        {
            my $his_abbr = $ss->get_role_abbr($his_idx);

            $idx = $self->add_role($his_role, $his_abbr);
#           print "Adding missing role $his_role idx=$idx\n";
        }
        else
        {
#           print "Found existing role $his_role idx=$idx\n";
        }


        $local_roles[$his_idx] = $idx;
    }

    #
    # Similar scan to ensure that we have rows for the genomes
    # that are in the other subsystem.
    #

    my @local_genomes;

    my @his_genomes = $ss->get_genomes();

    for my $his_idx (0..@his_genomes - 1)
    {
        my $genome = $his_genomes[$his_idx];


        my $my_idx = $self->get_genome_index($genome);

        if (!defined($my_idx))
        {
            #
            # Not there, need to add.
            #

            $my_idx = $self->add_genome($genome);
#           print "Adding missing genome $genome idx=$my_idx\n";
        }
        else
        {
#           print "Found existing genome $genome idx=$my_idx\n";
        }

        $local_genomes[$his_idx] = $my_idx;
    }


    #
    # Now that we have our local roles set up to receive the data,
    # process the incoming roles one at a time.
    #


    for my $his_role (@his_roles)
    {
        my $my_col = $self->get_col($local_roles[$his_role]);
        my $his_col = $ss->get_col($his_role);

        #
        # $his_col is the information for $his_role, indexed by
        # genome in @his_genomes.
        #
        # $my_col is hte information for my copy of $his_role,
        # indexed by genome in MY genome list.
        #

        my $my_role = $local_roles[$his_role];

#       print "merging: $self->{roles}->[$my_role] $ss->{roles}->[$his_role] his_role=$his_role my_role=$my_role\n";

        for my $his_gidx (0 .. @his_genomes - 1)
        {
            my $hisent = $his_col->[$his_gidx];

            my $my_gidx = $local_genomes[$his_gidx];


            my $myent = $my_col->[$my_gidx];

#           print "  his_gidx=$his_gidx my_gidx=$my_gidx hisent=@$hisent myent=@$myent\n";

            my %new;
            map { $new{$_}++ } @$hisent;
            map { $new{$_}++ } @$myent;

            @$myent = keys(%new);

#           print "  new entry: @$myent\n";
        }
    }

    #
    # Fix up the variant codes.
    #

    for my $his_gidx (0 .. @his_genomes - 1)
    {
        my $his_code = $ss->get_variant_code($his_gidx);
        my $my_gidx = $local_genomes[$his_gidx];

        if (!$self->get_variant_code($my_gidx))
        {
            $self->{variant_code}->[$my_gidx] = $his_code;
        }
    }

    #
    # If we are to add notes, append the other subsystem's notes text.
    #

    if ($add_notes)
    {
        my $his_notes = $ss->get_notes();

        $self->{notes} .= "\nNotes copied from $ss->{name}:\n$his_notes\n";
    }
}

sub dump
{
    my($self) = @_;

    for my $k (keys(%$self))
    {
        next if $k eq "spreadsheet" or $k eq "spreadsheet_inv";
        print "Key \"$k\": ", Dumper($self->{$k});
    }
}

#
# Increment the subsystem's version number.
#
sub incr_version {
    my($self) = @_;

    my $dir = $self->{dir};
    my $vfile = "$dir/VERSION";
    my($ver);

    if (open(my $fh,"<$vfile"))
    {
        if (defined($ver = <$fh>) && ($ver =~ /^(\S+)/))
        {
            $ver = $1;
        }
        else
        {
            $ver = 0;
        }
        close($fh);
    }
    else
    {
        $ver = 0;
    }

    $ver++;

    open(my $fh, ">$vfile") || die "could not open $vfile";
    print $fh "$ver\n";
    close($fh);

    chmod(0777, $vfile);

    $self->load_version();
}


=head3 functional_role_instances

    my @role_instances = $sub->functional_role_instances($role);

Returns the set of genes for a functional role that belong to
genomes with functional variants (> 0).

If the flag $strict is set to true,
an additional check for the correct function assignment is performed.
If the name of the functional role does not occur exaclty in the
latest function assignment of the PEG, it is not included in the
returned array. A simple index check is done.

=cut

sub functional_role_instances {

    my ($self, $role, $strict) = @_;
    my $i =0;

    my @instances;

    foreach my $cell (@{$self->get_col($self->get_role_index($role))}) {

	if ((scalar @$cell > 0) && ($self->get_variant_code($i) > 0)) {
	    foreach (@$cell) {


		unless ($strict) {
		    push @instances, $_;
		} else {
		    # check if the peg is still in sync with the role assignment
		    # will tolerate multiple role assignments but no mismatches
		    my $current_function = $self->{fig}->function_of($_);
		    if (index($current_function, $role) != -1) {
			push @instances, $_;
		    } else {
			print STDERR "[Warning] Function of $_ out of sync for role $role in subsystem ".$self->get_name()."\n";
		    }
		}
	    }
	}
	$i++;
    }


    return @instances if wantarray;
    return \@instances;

}




=head3 get_dir_from_name

    my $dirName = Subsystem::get_dir_from_name($name);

Return the name of the directory containing the SEED data for the specified
subsystem.

=over 4

=item name

Name of the subsystem whose directory is desired.

=item RETURN

Returns the fully-qualified directory name for the subsystem.

=back

=cut

sub get_dir_from_name
{
    my($name) = @_;

    my $b = $name;
    $b =~ s/ /_/g;
    my $dir = File::Spec->catfile($FIG_Config::data, 'Subsystems', $b);
    return $dir;
}

#
# Code for dealing with Bill McCune's prolog code for extending subsystems.
#
# The code here is a reconstruction of Bill's "go" script in perl with
# data pulled from the local SEED configuration.
#

sub extend_with_billogix
{
    my($self, $muser, $genomes) = @_;
    my($isMaster, $user);

    my $now = time();

    if ($muser =~ /master:(.*)/)
    {
        $isMaster = 1;
        $user = $1;
    }
    else
    {
        $isMaster = 0;
        $user = $muser;
    }

    #
    # initialize the genome list to all complete genomes, if none was passed in.
    #

    if (!$genomes)
    {
        $genomes = [$self->{fig}->genomes("complete")];
        warn "getting genome list from fig $self->{fig}";
    }

    #
    # Ensure genome list is of the right form.
    #

    if (ref($genomes) ne "ARRAY")
    {
        warn "billogix: genome list is not a list reference";
        return;
    }

    for my $g (@$genomes)
    {
        if ($g !~ /^\d+\.\d+/)
        {
            warn "billogix: genome '$g' is not of the proper form, aborting billogix run.";
            return;
        }
    }

    my $genome_list = "[" . join(", ", map { "'$_'" } @$genomes) . "]";

    warn "Genomes: $genome_list\n";
    warn Dumper($genomes);

    #
    # Find the executable.
    #

    my $exe = "$FIG_Config::bin/billogix";

    if (! -x $exe)
    {
        warn "Cannot find billogix exe at $exe\n";
        return;
    }

    my $ss_name = $self->{name};

    $ss_name =~ s/\s+/_/g;

    my $ss_dir = "$self->{dir}/";
    my $assign_dir = "$FIG_Config::data/Assignments/$user/";
    &FIG::verify_dir($assign_dir);

    my $when= strftime("%m-%d-%y:%H:%M:%S", localtime($now));
    my $job_id = "${when}:sss:$ss_name";

    my $seed = &FIG::cgi_url() . "/";
    my $export_part = "ssa.cgi?user=$muser&request=delete_or_export_ssa&export=";

    #
    # Have the prereq stuff, now start up the app.
    #

    $ENV{LOCALSZ} = "80000";
    $ENV{GLOBALSZ} = "80000";
    $ENV{TRAILSZ} = "30000";

    my $arch = &FIG::get_current_arch();

    $ENV{BILLOGIX} = "$FIG_Config::fig_disk/dist/releases/current/$arch/lib/Billogix";

    #
    # Need to ensure pl2wam is in our path
    #

    $ENV{PATH} = "${FIG_Config::ext_bin}:$ENV{PATH}";

    #
    # We're going to divide the run into $n_chunks chunks.
    #

    my $n_chunks = 10;

    my($log);
    open($log, ">$ss_dir/$job_id.log");

    for (my $this_chunk = 1; $this_chunk <= $n_chunks; $this_chunk++)
    {
        my $app_input = <<EOINP;
['\$BILLOGIX/top'].
loadup.
asserta(job_genome_list($genome_list)).
asserta(part($this_chunk, $n_chunks)).
asserta(url_default_seed('$seed')).
asserta(url_export_part('$export_part')).
asserta(ss_directory('$ss_dir')).
asserta(assign_directory('$assign_dir')).
asserta(job_id('$job_id')).
extend_test3('$ss_name').
EOINP

        print STDERR <<EOF;
Starting app

chunk $this_chunk of $n_chunks
ss_name = $ss_name
ss_dir = $ss_dir
user = $user
assign_dir = $assign_dir
exe = $exe
libdir = $ENV{BILLOGIX}
path = $ENV{PATH}

App input
$app_input
EOF
# feh, put in a block to reset perlmode indentation.
        {
            my($app_read, $app_write);

            #
            # Start the actual application with stdin and stdout redirected
            # to pipes.
            #
            # We write $app_input to the stdin pipe, and close it.
            # Then loop reading stdout, logging that output.
            #
            my $pid = open2($app_read, $app_write, $exe);

            if (!$pid)
            {
                warn "open2 $exe failed: $!\n";
                print $log "open2 $exe failed: $!\n";
                return;
            }

            print $app_write $app_input;
            close($app_write);

            #
            # Set autoflush on the logfile.
            #

            my $old = select($log);
            $| = 1;
            select(STDERR);
            $| = 1;
            select($old);

            warn "Starting $exe with pid $pid\n";
            print $log "Starting $exe with pid $pid\n";

            while (<$app_read>)
            {
                print STDERR $_;
                print $log $_;
            }

            print STDERR "App done\n";
            print $log "App done\n";

            close($app_read);

            my $ret = waitpid($pid, 0);
            my $stat = $?;
            print STDERR "Return status is $?\n";
            print $log "Return status is $?\n";

            #
            # This chunk has finished. We should see a file
            # rows.$this_chunk.$n_chunks.
            #
        }
    }
    #
    # At this point, the extension is finished (we've run the
    # $n_chunks parts of the extension job).
    #

    #
    # We read in all the individual rows files, writing the single
    # concatenation of rows.
    #

    my $ssaD = $self->{dir};

    my $rows_file = "$ssaD/rows";

    my $rowFH;
    if (!open($rowFH, ">$rows_file"))
    {
        my $err = "Cannot open rows file $ssaD/rows for writing: $!\n";
        print STDERR $err;
        print $log $err;
        return;
    }

    for (my $this_chunk = 1; $this_chunk <= $n_chunks; $this_chunk++)
    {
        my $chunkFH;
        my $cfile = "$ssaD/rows.$this_chunk.$n_chunks";
        if (!open($chunkFH, "<$cfile"))
        {
            my $err =  "Cannot open rows file $cfile for reading: $!\n";
            print STDERR $err;
            print $log $err;
            return;
        }
        while (<$chunkFH>)
        {
            print $rowFH $_;
        }
        close($chunkFH);
    }
    close($rowFH);

    #
    # Concatenate the assignments into the assignment directory.
    #

    my $assignments_file = "$assign_dir$job_id";
    my $assignFH;

    if (!open($assignFH, ">$assignments_file"))
    {
        my $err = "Cannot open assignments file $assignments_file for writing: $!\n";
        print STDERR $err;
        print $log $err;
        return;
    }

    for (my $this_chunk = 1; $this_chunk <= $n_chunks; $this_chunk++)
    {
        my $aFH;
        my $afile = "$ssaD/assignments.$this_chunk.$n_chunks";
        if (!open($aFH, "<$afile"))
        {
            my $err = "Cannot open assignments file $afile for reading: $!\n";
            print STDERR $err;
            print $log $err;
            return;
        }
        while (<$aFH>)
        {
            print $assignFH $_;
        }
        close($aFH);
    }
    close($assignFH);



    #
    # Back up the spreadsheet, and append the rows file to it.
    #

    &FIG::verify_dir("$ssaD/Backup");
    my $ts = time;
    rename("$ssaD/spreadsheet~","$ssaD/Backup/spreadsheet.$ts");
    copy("$ssaD/spreadsheet","$ssaD/spreadsheet~");
    rename("$ssaD/notes~","$ssaD/Backup/notes.$ts");

    #
    # Append the new rows to the spreadsheet.
    #

    my($ssafh, $rowsfh);
    open($ssafh, ">>$ssaD/spreadsheet") or die "Cannot open $ssaD/spreadsheet for append: $!\n";
    open($rowsfh, "<$ssaD/rows") or die "Cannot open $ssaD/rows for reading: $!\n";

    while (<$rowsfh>)
    {
        print $ssafh $_;
    }
    close($ssafh);
    close($rowsfh);

    $self->incr_version();
}

sub set_current_extend_pid
{
    my($self, $pid) = @_;

    if (open(my $fh, ">$self->{dir}/EXTEND_PID"))
    {
        print $fh "$pid\n";
    }
    else
    {
        warn "Cannot open $self->{dir}/EXTEND_PID: $!\n";
    }
}

sub get_current_extend_pid
{
    my($self) = @_;

    if (open(my $fh, "<$self->{dir}/EXTEND_PID"))
    {
        my $pid = <$fh>;
        close($fh);
        if ($pid)
        {
            chomp $pid;

            return $pid;
        }
    }
    return undef;
}

=head2 Static Utilities

These are internal static methods used by the Sprout Subsystem object
(SproutSubsys.pm). They insure that common functions are implemented with
common code.

=head3 GetDiagramIDs

    my @diagramIDs = Subsystem::GetDiagramIDs($subDir);

Return a list of the subsystem diagram IDs. The parameters are

=over 4

=item subDir

Fully-qualified directory name for the subsystem.

=item RETURN

Returns a list of the diagram IDs for this subsystem. Each diagram ID corresponds
to a diagram subdirectory in the subsystem's directory.

=back

=cut

sub GetDiagramIDs {
    # Get the parameters.
    my ($subDir) = @_;
    # Read the diagram subdirectories.
    opendir(D, "$subDir/diagrams");
    my @ids = grep { not /^\./ and -d "$subDir/diagrams/$_" } readdir(D);
    closedir(D);
    # Return the IDs.
    return @ids;

}

=head3 GetDiagramName

    my $name = Subsystem::GetDiagramName($subDir, $diagramID);

Return the name of the subsystem diagram with the specified ID.

=over 4

=item subDir

Subsystem directory name.

=item diagramID

ID of the diagram whose name is desired.

=item RETURN

Returns the name of the specified diagram, or C<undef> if the diagram does
not exist.

=back

=cut

sub GetDiagramName {
    # Get the parameters.
    my ($subDir, $diagramID) = @_;
    # Declare the return value.
    my $retVal;
    # Get the diagram's directory.
    my $ddir = "$subDir/diagrams/$diagramID";
    Trace("Looking for directory $ddir.") if T(3);
    # Only proceed if the directory exists.
    if (-d $ddir) {
        # Read the name.
        my $name = &FIG::file_head("$ddir/NAME", 1);
        # If there was no name, use the diagram ID.
        Trace("Diagram name is \"$name\".") if T(3);
        if (! $name) {
            Trace("Using default name $diagramID.") if T(3);
            $name = $diagramID;
        }
        # Lop off the line terminator.
        chomp ($name);
        # Return the result.
        $retVal = $name;
    }
    Trace("Returning diagram name \"$retVal\".") if T(3);
    return $retVal;
}

=head3 ComputeDiagramURLs

    my ($link, $imgLink) = Subsystem::ComputeDiagramURLs($self, $ssName, $diagramID, $sprout);

This is an internal static method that computes the URLs for a subsystem diagram.
It insures that both SEED and Sprout use the same rules for generating the
diagram URLs. The parameters are as follows.

=over 4

=item self

Relevant subsystem object.

=item ssName

Name of the relevant subsystem.

=item diagramID

ID of the relevant diagram.

=item sprout (optional)

If specified, indicates this should be a Sprout URL.

=item RETURN

Returns a two-element list, the first element of which is a link to the diagram
page, and the second of which is a link to the diagram image.

=back

=cut

sub ComputeDiagramURLs {
    # Get the parameters.
    my ($self, $ssName, $diagramID, $sprout) = @_;
    # Compute the CGI directory base.
    my $base = $FIG_Config::cgi_base;
    # Create the links.
    my $link;
    if ($self->is_new_diagram($diagramID)) {
        Trace("New diagram found for $diagramID in $ssName.") if T(3);
        $link = $base . "diagram.cgi?subsystem_name=$ssName&diagram=$diagramID";
    } else {
        Trace("Old diagram found for $diagramID in $ssName.") if T(3);
        $link = $base . "subsys_diagram.cgi?ssa=$ssName&diagram=$diagramID";
    }
    if ($sprout) {
        $link .= "&SPROUT=1";
    }
    my $img_link = $link . "&image=1";
    # Return them.
    return ($link, $img_link);
}

package Subsystem::Diagram;

sub new
{
    my($class, $sub, $fig, $name, $dir) = @_;

    if (!-d $dir)
    {
        return undef;
    }

    my $self = {
        fig => $fig,
        subsystem => $sub,
        name => $name,
        dir =>$ dir,
    };
    bless $self, $class;

    $self->load();

    return $self;
}

#
# Parse the diagram into internal data structure.
#

sub load
{
    my($self) = @_;

    $self->load_area();
}

sub load_area
{
    my($self) = @_;
    my $fh;

    if (!open($fh, "<$self->{dir}/area_table"))
    {
        warn "Could not load $self->{dir}/area_table: $!\n";
        return;
    }

    $self->{areas} = [];

    my $area_list = $self->{areas};

    while (<$fh>)
    {
        chomp;
        s/#.*$//;
        s/^\s+//;
        s/\s+$//;
        next if $_ eq '';
        my ($area, $tag, $value) = split(/\s+/, $_, 3);
        # print "area=$area tag=$tag value=$value\n";

        push(@$area_list, [$area, $tag, $value]);

        #
        # Do a little checking.
        #

        if ($tag eq "role")
        {
            my $idx = $self->{subsystem}->get_role_index($value);
            if (!defined($idx))
            {
                warn "Role not found for \"$value\" in subsystem $self->{subsystem}->{name}\n";
            }
        }
    }
    close($fh);
}

sub get_areas
{
    my($self) = @_;

    return @{$self->{areas}};
}


1;

=head2 Method Listing

=over 4

=item index_cell

Create the subsystem_index entries for the given cell.
(NEW).

=item delete_role(name)

Delete the given role.

=item add_role(name, abbr)

Add a new role.

=item get_subset(name)

A deprecated form of get_subsetC

=item get_subsetC(name)

Returns a given subset. A subset is an object, implemented as a blessed array
of roles.

=item add_genome(genome_id, variant_code)

=item remove_genome(genome_id)

=back

=cut

