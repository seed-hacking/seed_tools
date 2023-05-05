#
# Given a set of genomes, find the subsystems that are appropriate for them and fill.
#

use strict;
use FIG;
use Data::Dumper;

use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o genome [genome, genome...]",
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die ($usage->text) if @ARGV == 0;

my @genomes = @ARGV;

my $fig = FIG->new;
my $dbh = $fig->db_handle->dbh;

my $qs = join(", ", map { "?" } @genomes);

my $res = $dbh->selectall_arrayref(qq(
SELECT DISTINCT org,subsystem 
FROM roles JOIN subsystem_index ON subsystem_index.role =roles.role 
				   WHERE org IN ($qs)), undef, @genomes);

my %ss_genomes;
push(@{$ss_genomes{$_->[1]}}, $_->[0]) foreach @$res;

for my $ss (sort keys %ss_genomes)
{
    fill_ss($ss, $ss_genomes{$ss});
    print "Filled $ss\n";

}

sub fill_ss
{
    my ($ss_name, $genomes) = @_;

    my $subsystem = $fig->get_subsystem($ss_name);
    $subsystem or die "Cannot get $ss_name\n";
    
    my $comment = "";
    
    foreach my $genome ( @$genomes ) 
    {
	my $rawgenome = $genome;
	my $loc;

	next unless$fig->is_genome( $genome );

	my $idx = $subsystem->get_genome_index($genome);
	if ( defined( $idx ) ) {
	    warn "$ss_name already has $genome\n";
	    next;
	}
      
	$subsystem->add_genome( $genome );
	my $idx = $subsystem->get_genome_index($genome);
	$subsystem->{variant_code}->[$idx] = '?';
	
	foreach my $role ( $subsystem->get_roles() ) {
	    my @inpegs = $subsystem->get_pegs_from_cell( $genome, $role );
	    next if ( @inpegs > 0 );
	    
	    my @pegs = $fig->seqs_with_role( $role, "master", $genome);
	    @pegs = grep { $subsystem->in_genome( $genome, $_ ) } @pegs;
	    my %tmppegs = map { $_ => 1 } @pegs;
	    @pegs = keys %tmppegs;
	    
	    $subsystem->set_pegs_in_cell( $genome, $role, \@pegs);
	}
    }
    $subsystem->incr_version();
    $subsystem->db_sync();
    $subsystem->write_subsystem();
}
