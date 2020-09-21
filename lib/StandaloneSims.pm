
#
# Code to bundle up the manipulation of standalone sims files.
#
# A StandaloneSimsSet contains up to four files: raw sims, their flips, expanded sims, and their flips.
#
# Each sims file is represente by a StandaloneSimsFile object which wraps up the low level indexing, lookup,
# and retrieval code.
#
#

package StandaloneSimsSet;

use FIG;
use strict;

my %rast_names = (flips_exp => 'expanded_similarities.flips',
		  flips_raw => 'similarities.flips',
		  fwd_raw => 'similarities',
		  fwd_exp => 'expanded_similarities');

=head3 new

    my $sims_set = StandaloneSimsSet->new($dir, %opts).

Constructor for a StandaloneSimsSet. Argument is a directory that is expected
to contain four sets of sim sets: fwd_raw, fwd_exp, flips_raw, flips_exp.

%opts may contain overrides for the filenames for the four sets listed above. If the
option rast => 1 is passed, the usual RAST names for the files are used.

Any of the four may be empty.

=cut

sub new
{
    my($class, $dir, %opts) = @_;

    my $self = {
	dir => $dir,
	opts => \%opts,
    };

    for my $set (qw(fwd_raw fwd_exp flips_raw flips_exp))
    {
	my $file = "$dir/$set";
	if ($opts{rast})
	{
	    $file = "$dir/$rast_names{$set}";
	}
	elsif ((my $alt = $opts{$set}) ne '')
	{
	    if ($alt =~ m,^/,)
	    {
		$file = $alt;
	    }
	    else
	    {
		$file = "$dir/$alt";
	    }
	}
	$self->{$set} = StandaloneSimsFile->new($file);
    }

    return bless $self, $class;
}

sub sims_fwd
{
    my ($self, $id, $maxN, $maxP, $select, $max_expand, $filters) = @_;

    return $self->get_sims('fwd', $id, $maxN, $maxP, $select, $max_expand, $filters);
}

sub sims_flip
{
    my ($self, $id, $maxN, $maxP, $select, $max_expand, $filters) = @_;

    return $self->get_sims('flips', $id, $maxN, $maxP, $select, $max_expand, $filters);
}

sub get_sims
{
    my ($self, $dir, $id_or_list, $maxN, $maxP, $select, $max_expand, $filters) = @_;

    my $cooked = ($select eq 'raw') ? 'raw' : 'exp';
    my $set_name = $dir . '_' . $cooked;

    my $set = $self->{$set_name};

    my $filter_func =  &FIG::create_sim_filter(undef, $maxP, $filters);

    my @out;
    for my $id (ref($id_or_list) ? @$id_or_list : $id_or_list)
    {
	my $start_len = @out;

	push(@out, $set->get_sims($id, $filter_func));

	if (defined($maxN) and (@out - $start_len > $maxN))
	{
	    $#out = $start_len + $maxN - 1;
	}
    }
    return @out;
}

package StandaloneSimsFile;

use FIG;
use strict;
use DB_File;

#
# A standalone sims file has the sims data in $file and
# a btree index in $file.index.
#
sub new
{
    my($class, $file, $index_file) = @_;
    
    my $hash = {};

    $index_file = "$file.index" if $index_file eq '';
    
    my $tied = tie %$hash, 'DB_File', $index_file, O_RDONLY, 0666, $DB_BTREE;

    my $fh = new FileHandle("<$file");

    my $self = {
	file => $file,
	index_file => $index_file,
	fh => $fh,
	tied => $tied,
	hash => $hash,
    };
    return bless $self, $class;
}

#
# Retrieve sims. We don't do any filtering here.
#
sub get_sims
{
    my($self, $id, $filter_func) = @_;

    my $info = $self->{hash}->{$id};

    if ($info !~ /^(\d+),(\d+)$/)
    {
	return;
    }
    
    my($seek, $len) = ($1, $2);

    my $sims_txt = &FIG::read_block($self->{fh}, $seek, $len);

    my @sims;

    if (ref($filter_func))
    {
	@sims = map { bless $_, 'Sim' } grep { &$filter_func($_) }  map { [split(/\t/)] } @$sims_txt;
    }
    else
    {
	@sims = map { bless $_, 'Sim' } map { [split(/\t/)] } @$sims_txt;
    }

    return @sims;
}

1;
