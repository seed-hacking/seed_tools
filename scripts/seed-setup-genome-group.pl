#
# Set up a genome group from an excel spreadsheet.
# id-column is the column number that has the genome ids
# group-column is the column number that has the group name to load
#

use 5.010;
use strict;
use Spreadsheet::ParseXLSX;
use Spreadsheet::ParseExcel::Utility qw(col2int);
use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o spreadsheet.xlsx id-column group-column [group-column, ...]",
				    ["write-group" => "Write the found ids to a file named by the group."],
				    ["output-dir=s" => "Write to the given output directory."],
				    ["help|h" => "Show this help message."]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV >= 3;

my $ss_file = shift;
my $id_col = shift;
my @group_cols = @ARGV;

$id_col = col2int($id_col) if $id_col =~ /[a-z]/i;
@group_cols = map { col2int($_) if /[a-z]/i } @group_cols;

my $parser = Spreadsheet::ParseXLSX->new;
my $wb = $parser->parse($ss_file);

$wb or die "Failed to parse $ss_file";

my $ws = ($wb->worksheets)[0];

my($rstart, $rend) = $ws->row_range;

my %group;

for my $group_col (@group_cols)
{
    for my $row ($rstart .. $rend)
    {
	my $id_cell = $ws->get_cell($row, $id_col);
	my $group_cell = $ws->get_cell($row, $group_col);
	next unless $id_cell && $group_cell;
	my $id = $id_cell->value;
	my $group = $group_cell->value;
	
	next unless $group;
	next unless $id =~ /^\s*(\d+\.\d+)\s*$/;
	$id = $1;
	push(@{$group{$group}}, $id);
    }
}

if ($opt->write_group)
{
    my $dir = $opt->output_dir // ".";
    die "Output directory $dir does not exist" unless -d $dir;
    while (my($gname, $list) = each %group)
    {
	open(F, ">", "$dir/$gname") or die "Cannot write $gname: $!";
	say F $_ foreach @$list;
	close(F);
    }
}




