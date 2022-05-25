#
# Tools for dealing with nonredundant databases - finding sources, building NR, etc.
#

package NRTools;

use strict;

use base qw(Exporter);
use vars qw(@EXPORT);

eval {
    require Job48;
    import Job48;
};

@EXPORT = qw(scan_NR_dir scan_seed_dir scan_rast_jobs);


=head3 scan_NR_dir()

usage: scan_NR_dir(\%nr_hash, $dirname)

Scan a directory containing SEED-formatted NR directories, and fill in
%nr_hash with entries of the form $nr_hash->{name} = { name => dirname, path => full path to NR dir, size => size of fasta file}

=cut

sub scan_NR_dir
{
    my($nr_hash, $dir, $options) = @_;
    
    my $dh = new DirHandle($dir);
    $dh or die "Cannot open directory $dir: $!";
    while (defined($_ = $dh->read()))
    {
	next if /^\./;
	next if $_ eq 'SEED';	# Ignore the export of SEED data that lives here.
	
	next if $options->{skip} and /$options->{skip}/;
	my $path = "$dir/$_";
	my $fasta = "$path/fasta";
	if (-f $fasta)
	{
	    if (! -f "$path/assigned_functions")
	    {
		warn "NR directory $path missing assigned_functions\n";
	    }
	    if (! -f "$path/org.table")
	    {
		warn "NR directory $path missing org.table\n";
	    }
	    $nr_hash->{$_} = { type => "NR", name => $_, path => $path, fasta_path => $fasta, size => -s $fasta };
	}
    }
    $dh->close();
}

=head3 scan_seed_dir()

usage: @fasta = scan_seed_dir(\%nr_hash, dirname)

Scan a SEED organism directory, creating entries as in scan_NR_dir.

=cut

sub scan_seed_dir
{
    my($nr_hash, $dir, $opts) = @_;

    my $dh = new DirHandle($dir);
    my $n = 0;
    while ($_ = $dh->read())
    {
	next if /^\./;
	#next if /^9999999.\d+$/;
	
	#
	# Strip environmental sequences.
	# 
	# c.f. seed-tech mail thread of 2/8/2007 for discusson on the rationale of the following
	# logic.
	#
	# next if $fig->is_environmental($_);
	next if /^4{7}/ or /^9{7}/;

	my $path = "$dir/$_";

	next unless -d $path;
	next if (-e "$path/DELETED");

	my $fasta = "$path/Features/peg/fasta";
	if (-f $fasta)
	{
	    $nr_hash->{$_} = { type => "seed_org", name => $_, path => $path,
				   fasta_path => $fasta, size => -s _ };
	}
	last if $opts->{limit} && $n++ > $opts->{limit};
    }
    $dh->close();
}

=head3 scan_seed_dir()

usage: @jobs = scan_rast_jobs($dir)

Scan the given RAST job directory, finding all completed jobs that are marked
with import.candidate nonzero and import.action set to "import".

=cut

sub scan_rast_jobs
{
    my($jobs, $dir) = @_;

    my $dh = new DirHandle($dir);

    if (!$dh)
    {
	warn "Cannot open directory $dir: $!";
	return;
    }

    while (defined($_ = $dh->read()))
    {
	next unless /^\d+$/;

	my $job = Job48->new("$dir/$_");
	next unless $job;
	next unless $job->meta->get_metadata("status.final") eq  "complete";
	next unless $job->meta->get_metadata("import.candidate") > 0;
	next unless $job->meta->get_metadata("import.action") eq 'import';
	my $stat =  $job->meta->get_metadata('import.status');
	next if $stat eq 'computed' or $stat eq 'installed';

	push(@$jobs, $job);
    }
}

1;

