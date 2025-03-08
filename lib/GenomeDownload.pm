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


my @allowed = qw(
	genbank gbk
	genbank_merged gbk
	spreadsheet_txt txt
	spreadsheet_xls xls
	feature_data txt
	protein_fasta faa
	contig_fasta fna
	feature_dna ffn
	seed_dir seed
	patric_features txt
	patric_specialty_genes txt
	patric_genome_metadata txt
	gff gff
	embl embl
	);
my @need_contig = qw(
	genbank
	genbank_merged
	contig_fasta
	seed_dir
	gff
	embl
	);
my %allowed;
my %need_contig = map { $_ => 1 } @need_contig;
my %suffix;
for (my $i = 0; $i < @allowed; $i += 2)
{
    $allowed{$allowed[$i]} = 1;
    $suffix{$allowed[$i]} = $allowed[$i + 1];
}

sub get_export
{
    my($gid, $type) = @_;

    my $params = { skip_contigs => !$need_contig{$type} };

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

get '/doc' => sub {
    content_type 'text/plain';
    return "Allowed types:\n" . join("\n", @allowed) . "\n";
};

get '/:gid/:type' => sub {
    my $gid = params->{gid};
    my $type = params->{type};

    if (!$allowed{$type}) {
	send_error("not found", 404);
    }

    my $export = get_export($gid, $type);

    send_file( \$export, content_type => 'application/x-fasta',
                   filename     => "$gid.$suffix{$type}" );
};

dance;
