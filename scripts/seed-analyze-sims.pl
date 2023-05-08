#
# Look up the sims for the given peg as well as blasting with BLAST and diamond. Emit a summary chart.
#

use strict;
use 5.010;
use File::Temp;
use Data::Dumper;
use FIG;

@ARGV == 1 or die "Usage: $0 pegid\n";

my $nr = "/opt/seed/sims/004/nr";

my $fig = new FIG;

my $peg = shift;

my $trans = $fig->get_translation($peg);
my $seq = File::Temp->new;
print $seq ">$peg\n$trans\n";
close($seq);

my @sims = $fig->sims($peg, 100, 100, 'fig', 100);
print "Sims\n";
for my $s (@sims)
{
    my($q, $id, $iden, $len, $mis, $gap, $b1, $e1, $b2, $e2, $evalue) = @$s;
    print join("\t", $id, $evalue), "\n";
}

print "Diamond\n";
open(D, "-|", "diamond", "blastp", 
     "-c1",
     "--sensitive",
     "-e", "1e-4",
     "--masking", 0,
     "-k", 200,
     "--outfmt", 6,
     "-d", "$nr.dmnd",
     "-q", "$seq") or die "Diamond failed: $!";
while (<D>)
{
    chomp;
    my($id1, $id2, $iden, $ali, $mis, $gaps, $b1, $e1, $b2, $e2, $evalue) = split(/\t/);
    my($h) = $id2 =~ /md5\|(\S+)/;
    my @p = $fig->pegs_with_md5($h);
    next if grep { $_ eq $peg } @p;
    print join("\t", join(",", @p), $evalue), "\n";
}
close(D);

print "BLAST\n";
open(D, "-|", "blastall",
     "-p", "blastp",
     "-e", "1e-4",
     "-m", 8,
     "-d", $nr,
     "-i", "$seq") or die "BLAST failed: $!";
while (<D>)
{
    chomp;
    my($id1, $id2, $iden, $ali, $mis, $gaps, $b1, $e1, $b2, $e2, $evalue) = split(/\t/);
    my($h) = $id2 =~ /md5\|(\S+)/;
    my @p = $fig->pegs_with_md5($h);
    next if grep { $_ eq $peg } @p;
    print join("\t", join(",", @p), $evalue), "\n";
}
