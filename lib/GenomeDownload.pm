package GenomeDownload;
#
# Little dancer app to enable downloads of formatted genome data
#

use JSON::XS;
use Dancer2;
use strict;
use IPC::Run qw(run);
use FIG_Config;
use FIG;

my $json = JSON::XS->new;

my %allowed = map { $_ => 1 } qw(protein_fasta genbank  patric_features);

sub get_export
{
    my($gid, $type) = @_;

    my $params = { };
    if ($type  ne "genbank")
    {
	$params->{skip_contigs} = 1;
    }

    my $fig = FIG->new;
    my $gobj = $fig->genome_id_to_genome_object($gid, $params);
    my $gobj_js = $json->encode($gobj);
    my $export;

    my $ok = run(["rast_export_genome", $type],
	"<", \$gobj_js,
	">", \$export);
    $ok or die "Error exporting";
    return $export;
}

get '/:gid/:type' => sub {
    my $gid = params->{gid};
    my $type = params->{type};

    if (!$allowed{$type}) {
	send_error("not found", 404);
    }

    my $export = get_export($gid, $type);

    send_file( \$export, content_type => 'application/x-fasta',
                   filename     => "$gid.faa" );
};

dance;
