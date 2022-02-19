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

package MetaSubsystem;

#use Carp;

#use POSIX;
use DirHandle;
#use Data::Dumper;
#use File::Copy;
#use File::Spec;
#use IPC::Open2;
#use FileHandle;
#use Tracer;
use Subsystem;
use strict;

1;

sub new {
    
    my ( $class, $name, $fig, $create, $subsystems, $genomes, $subsets, $view, $description ) = @_;

    my $msadir = get_dir_from_name( $name );
    

    if ( ! -d $msadir and not $create )    {
      return undef;
    }

    # RAE: Please do this:
    $name =~ s/^\s+//; $name =~ s/\s+$//;
    $name =~ s/ /_/g;

    my $self = {
		dir => $msadir,
		name => $name,
		fig => $fig,
	       };

    bless($self, $class);

    #
    # Check to see if the database we're running against has a variant column.
    #
#    $self->detect_db_version();

    if ($create)
    {
        my $succ = $self->create_metasubsystem( $subsystems, $genomes, $subsets, $view, $description );
	if ( $succ eq -1 ) {
	  return -1;
	}
    }
    else
    {
        $self->load();
    }

    return $self;
}


=head3 create_metasubsystem

Create a new subsystem. This creates the subsystem directory in the
correct place ($FIG_Config::data/Subsystems), and populates it with
the correct initial data.

=cut

sub create_metasubsystem {
  my( $self, $subsystems, $genomes, $subsets, $view, $description ) = @_;

  my $dir = $self->{dir};
  my $fig = $self->{fig};

  $self->{ genomes } = $genomes;
  $self->{ subsystems } = $subsystems;
  $self->{ subsets } = $subsets;
  $self->{ view } = $view;
  $self->{ description } = $description;
  
  if (-d $dir)  {
    warn "Not creating: MetaSubsystem directory $dir already exists";
    return -1;
  }

  $fig->verify_dir($dir);
  
  $self->write_metasubsystem();
}

sub write_metasubsystem {
  my( $self, $force_backup ) = @_;

  my $dir = $self->{dir};
  my $fig = $self->{fig};
  
  #
  # We first move the existing spreadsheet and notes files (if present)
  # to spreadsheet~ and notes~, and current state.
  #
  
  my $genomes_file = "$dir/genomes";
  my $genomes_bak = "$dir/genomes~";
  my $subsystems_file = "$dir/subsystems";
  my $subsystems_bak = "$dir/subsystems~";
  my $subsets_file = "$dir/subsets";
  my $subsets_bak = "$dir/subsets~";
  my $view_file = "$dir/view";
  my $view_bak = "$dir/view~";
  my $description_file = "$dir/description";
  my $description_bak = "$dir/description~";

  if ( -f $genomes_file ) {
    rename( $genomes_file, $genomes_bak);
  }
  if ( -f $subsystems_file ) {
    rename( $subsystems_file, $subsystems_bak);
  }
  if ( -f $subsets_file ) {
    rename( $subsets_file, $subsets_bak);
  }
  if ( -f $view_file ) {
    rename( $view_file, $view_bak);
  }
  if ( -f $description_file ) {
    rename( $description_file, $description_bak);
  }

  eval {
    my $fh;
    
    open( $fh, ">$genomes_file" ) or die "Cannot open $genomes_file for writing: $!\n";
    my $genomes = $self->{ genomes };
    foreach $_ ( keys %$genomes )  {
      print $fh "$_\n";
    }
    close( $fh );
    chmod( 0777,$genomes_file );
    
    open( $fh, ">$subsets_file" ) or die "Cannot open $subsets_file for writing: $!\n";
    my $subsets = $self->{ subsets };
    if ( defined( $subsets ) ) {
      foreach my $ss ( keys %$subsets )  {
	foreach my $sub ( keys %{ $subsets->{ $ss } } ) {
	  chomp $sub;
	  print $fh "$ss\t$sub\n";
	}
      }
    }
    close( $fh );
    chmod( 0777,$subsystems_file );
    
    open( $fh, ">$subsystems_file" ) or die "Cannot open $subsystems_file for writing: $!\n";
    my $subsystems = $self->{ subsystems };
    foreach my $sst ( keys %$subsystems )  {
      chomp $sst;
      print $fh "$sst\n";
    }
    close( $fh );
    chmod( 0777,$subsystems_file );
    
    open( $fh, ">$view_file" ) or die "Cannot open $view_file for writing: $!\n";
    my $view = $self->{ view };
    if ( defined( $view ) ) {
      my %subsets_view = %{ $view->{ 'Subsets' } };
      foreach my $ssname ( keys %subsets_view )  {
	my $visible = $subsets_view{ $ssname }->{ 'visible' };
	$visible = 0 if ( !defined( $visible ) );
	my $collapsed = $subsets_view{ $ssname }->{ 'collapsed' };
	$collapsed = 0 if ( !defined( $collapsed ) );
	print $fh "Subset\t$ssname\t$visible\t$collapsed\n";
      }
      if ( defined( $view->{ 'Roles' } ) ) {
	%subsets_view = %{ $view->{ 'Roles' } };
	foreach my $ssname ( keys %subsets_view )  {
	  my $visible = $subsets_view{ $ssname }->{ 'visible' };
	  $visible = 0 if ( !defined( $visible ) );
	  
	  if ( $ssname =~ /(.*)##-##(.*)/ ) {
	    my $role = $1;
	    my $subsystem = $2;
	    $subsystem = 0 if ( !defined( $subsystem ) );
	    print $fh "Role\t$ssname\t$subsystem\t$visible\n";
	  }
	}
      }
    }
    close( $fh );
    chmod( 0777,$subsystems_file );
    
    open( $fh, ">$description_file" ) or die "Cannot open $description_file for writing: $!\n";
    my $description = $self->{ description };
    chomp $description;
    print $fh "$description\n";
    close( $fh );
    chmod( 0777,$description_file );
    
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
    
    #        $self->update_backups($force_backup);
    $self->make_backup();
    
    if ($self->{version} < 100) { $self->{version} += 100 }
    open($fh, ">$dir/VERSION") or die "Cannot write $dir/VERSION: $!\n";
    print $fh "$self->{version}\n";
    close($fh);
    chmod(0777,"VERSION");
  };
  
  if ( $@ ne "" ) {
    warn "Spreadsheet write failed, reverting to backup. Error was\n$@\n";
  }
}


sub make_backup {
    my($self) = @_;

    my $dir = $self->{dir};
    my $bak = "$dir/Backup";

    $self->{fig}->verify_dir($bak);

    my $ts = time;

    rename("$dir/genomes~", "$bak/genomes.$ts");
    rename("$dir/subsystems~", "$bak/subsystems.$ts");
    rename("$dir/subsets~", "$bak/subsets.$ts");
    rename("$dir/view~", "$bak/view.$ts");
    $self->{version}++;
}



sub load {

    my($self) = @_;

    #
    # Load the subsystem.
    #

    my $ssa;
    my $genomes;
    my $subsets;
    my $view; 
    if  ( !open( $ssa,"<$self->{dir}/subsystems" ) ) {
      $self->{ empty_ss }++;
      return;
    }
    if  ( !open( $genomes,"<$self->{dir}/genomes" ) ) {
      $self->{ empty_ss }++;
      return;
    }
    if  ( !open( $subsets,"<$self->{dir}/subsets" ) ) {
      return;
    }
    if  ( !open( $view,"<$self->{dir}/view" ) ) {
      return;
    }

    $self->load_subsystems( $ssa );
    $self->load_genomes( $genomes );
    $self->load_subsets( $subsets );
    $self->load_view( $view );
    $self->load_description();
    $self->load_version();
    $self->load_curation();

    close $ssa;
    close $genomes;
    close $subsets;
    close $view;

    return 1;
}

sub load_subsystems {

  my ( $self, $ssa ) = @_;
  my $fig = $self->{ fig };

  my $sshandles;
  
  while( my $ssname = <$ssa> ) {
    chomp $ssname;
    my $sshandle = new Subsystem( $ssname, $fig, 0 );

    $self->{ subsystems }->{ $ssname } = $sshandle;
  }
  return 1;
}

sub load_genomes {
  my ( $self, $genomes ) = @_;
  
  my @genomesarr;
  while ( my $thisgenome = <$genomes> ) {
    chomp $thisgenome;
    push @genomesarr, $thisgenome;
  }
  
  my %genomeshash = map { $_ => 1 } @genomesarr;

  $self->{ genomes } = \%genomeshash;
}

sub load_subsets {
  my ( $self, $subsets ) = @_;
  
  my $subsethash;
  while ( my $thissubset = <$subsets> ) {
    chomp $thissubset;
    my ( $ssname, $ssabb ) = split( "\t", $thissubset );
    $subsethash->{ $ssname }->{ $ssabb } = 1;
  }

  $self->{ subsets } = $subsethash;
}


sub load_view {
  my ( $self, $view ) = @_;
  
  my $viewhash;
  while ( my $line = <$view> ) {
    chomp $line;
    my ( $what, $name, $third, $fourth ) = split( "\t", $line );
    if ( $what eq 'Subset' ) {
      $viewhash->{ 'Subsets' }->{ $name }->{ 'visible' } = $third;
      $viewhash->{ 'Subsets' }->{ $name }->{ 'collapsed' } = $fourth;
    }
    elsif ( $what eq 'Role' ) {
      $viewhash->{ 'Roles' }->{ $name }->{ 'visible' } = $fourth;
    }
  }
  $self->{ view } = $viewhash;
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

sub get_dir_from_name {
  my ( $name ) = @_;
  
  my $b = $name;
  $b =~ s/ /_/g;
  my $dir = File::Spec->catfile($FIG_Config::data, 'MetaSubsystems', $b);
  return $dir;
}


=head3 get_curator

    my $userName = $sub->get_curator();

Return the name of this subsystem's official curator.

=cut

sub get_curator {
  my ( $self ) = @_;
  return $self->{ 'curator' };
}


sub get_created {
  my( $self ) = @_;
  return $self->{ 'created' };
}

sub get_last_updated {
  my ( $self ) = @_;
  return $self->{ 'last_updated' };
}

sub load_curation {
  my ( $self ) = @_;
  
  if  ( open( LOG, "<$self->{dir}/curation.log" ) ) {
    my $last = 0;
    while ( defined( $_ = <LOG> ) ) {
      if ( /^(\d+)\t(\S+)\s+started/ ) {           
	my $tmpcurator = $2;
	$self->{ 'created' } = $1;
	if ( $tmpcurator =~ /master\:(.*)/ ) {
	  $self->{ 'curator' } = $1;
	}
	else {
	  $self->{ 'curator' } = $tmpcurator;
	}
      }
      if ( ( /^(\d+)/ ) && ( $1 > $last ) ) {
	$last = $1;
      }
    }
    close( LOG );
    if ($last) { $self->{last_updated} = $last; }
  }
}

sub load_version {
  my ( $self ) = @_;

  my @l = &FIG::file_head(File::Spec->catfile($self->{dir}, "VERSION"), 1);
  my $l = $l[0];
  chomp $l;
  $self->{version} = $l;
}

sub load_description {
  my ( $self ) = @_;

  if  ( open( DESC, "<$self->{dir}/description" ) ) {
    my $description;
    while ( defined( $_ = <DESC> ) ) {
      $description .= $_;
    }
    $self->{description} = $description;
  }
}

sub update_curation_log {
  my( $self ) = @_;
  
  my $fh;
  my $file = "$self->{dir}/curation.log";
  
  my $now = time;
  my $user = $self->{fig}->get_user();
  
  if ( -f $file ) {
    open( $fh, ">>$file" ) or die "Cannot open $file for writing: $!\n";
  }
  else {
    open($fh, ">$file") or die "Cannot open $file for writing: $!\n";
    print $fh "$now\t$user\tstarted\n";
  }
  print $fh "$now\t$user\tupdated\n";
  close( $fh );
}


=head3 get_description

    my $text = $sub->get_description();

Return the description for this subsystem.

=cut

sub get_description {
  my( $self ) = @_;

  return $self->{description};
}

sub set_description {
    my ( $self, $desc ) = @_;

    $self->{description} = $desc;
}

sub remove_genomes {
    my ( $self, $delgenomes ) = @_;

    my $genomes = $self->{ genomes };
    foreach $_ ( keys %$genomes )  {
      if ( defined( $delgenomes->{ $_ } ) ) {
	delete $genomes->{ $_ };
      }
    }
    $self->{ genomes } = $genomes;
}

sub add_genomes {
    my ( $self, $newgenomes ) = @_;

    my $genomes = $self->{ genomes };
    foreach $_ ( keys %$newgenomes )  {
      print STDERR $_." SHOULD BE ADDED\n";
      if ( ! defined( $genomes->{ $_ } ) ) {
	$genomes->{ $_ } = 1;
	print STDERR $_." IS BEING ADDED\n";
      }
    }
    $self->{ genomes } = $genomes;
}

#
# Increment the subsystem's version number.
#
sub incr_version {
  my ( $self ) = @_;
  
  my $dir = $self->{dir};
  my $vfile = "$dir/VERSION";
  my($ver);
  
  if ( open( my $fh,"<$vfile" ) ) {
    if ( defined( $ver = <$fh> ) && ( $ver =~ /^(\S+)/ ) ) {
      $ver = $1;
    }
    else {
      $ver = 0;
    }
    close($fh);
  }
  else {
    $ver = 0;
  }
  
  $ver++;
  
  open( my $fh, ">$vfile" ) || die "could not open $vfile";
  print $fh "$ver\n";
  close($fh);

  chmod( 0777, $vfile );

  $self->load_version();
}

sub get_curator_from_metaname {
  my ( $metaname ) = @_;
  my $msadir = get_dir_from_name( $metaname );
  
  if  ( open( LOG, "<$msadir/curation.log" ) ) {
    my $last = 0;
    while ( defined( $_ = <LOG> ) ) {
      if ( /^(\d+)\t(\S+)\s+started/ ) {           
	my $tmpcurator = $2;
	if ( $tmpcurator =~ /master\:(.*)/ ) {
	  return $1;
	}
	else {
	  return $tmpcurator;
	}
      }
    }
    close( LOG );
  }
}
