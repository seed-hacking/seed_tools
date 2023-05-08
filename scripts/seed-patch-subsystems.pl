#
# Given a set of genomes and subsystems, patch the variant code to ?
#

use strict;
use FIG;
use Data::Dumper;

use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o genome-file ss-file",
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die ($usage->text) if @ARGV != 2;

my $genome_file = shift;
my $ss_file = shift;

my @subsystems;
open(S, "<", $ss_file) or die "Cannot open $ss_file: $!";
while (<S>)
{
    if (/^\s*(.*?)\s*$/)
    {
	push(@subsystems, $1);
    }
}
close(S);
my @genomes;
open(G, "<", $genome_file) or die "Cannot open $genome_file: $!";
while (<G>)
{
    if (/(\d+\.\d+)/)
    {
	push(@genomes, $1);
    }
}
close(G);

my $fig = FIG->new;

for my $s (@subsystems)
{
    patch_ss($s, \@genomes);
}

sub patch_ss
{
    my ($ss_name, $genomes) = @_;

    my $subsystem = $fig->get_subsystem($ss_name);
    $subsystem or die "Cannot get $ss_name\n";
    
    foreach my $genome ( @$genomes ) 
    {
	next unless$fig->is_genome( $genome );

	my $idx = $subsystem->get_genome_index($genome);
	print "idx=$idx\n";
	if ( !defined( $idx ) ) {
	    warn "$genome not found in $ss_name, skipping\n";
	    next;
	}
      
	if ($subsystem->{variant_code}->[$idx] eq 0)
	{
	    $subsystem->{variant_code}->[$idx] = '?';
	}
    }
    $subsystem->incr_version();
    $subsystem->db_sync();
    $subsystem->write_subsystem();
}
