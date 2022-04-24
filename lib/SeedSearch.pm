package SeedSearch;

use FIG_Config;
use Carp;
use strict;

our @doctypes = qw(peg
		   rna
		   atn
		   att
		   bs
		   opr
		   pbs
		   pi
		   pp
		   prm
		   pseudo
		   rsw
		   sRNA
		   trm
		   box
		   crispr_array
		   );

our %tmap;
for my $i (0..$#doctypes)
{
    $tmap{$doctypes[$i]} = $i;
    $tmap{$i} = $doctypes[$i];
}

=head1 SeedSearch

This package uses a Sphinx indexing engine to do fast lookups into
a SEED sphinx index.

=cut

sub new
{
    my($class, $params) = @_;

    if (!defined($params))
    {
	if (@FIG_Config::search_params)
	{
	    $params = \@FIG_Config::search_params;
	}
	else
	{
	    confess "SeedSearch requires a Sphinx configuration in @FIG_Config::search_params";
	}
    }

    my $sphinx = Sphinx::Search->new();

    $sphinx->SetServer($params);
    
    my $self = {
	params => $params,
	sphinx => $sphinx,
    };
    return bless $self, $class;
}

=head2 search
    
    my @results = $sphinx->search($search_terms,
			      -pagenum => i,
			      -pagesize => N)
				    

=cut    
sub start_search
{

}

sub fid_to_docid
{
    my($fid) = @_;
    
    if ($fid =~ /^fig\|(\d+)\.(\d+)\.([^.]+)\.(\d+)$/)
    {
	my ($g, $ext, $type, $num) = ($1, $2, $3, $4);
	return undef unless exists $tmap{$type};
	my $tnum = $tmap{$type};

	#
	# right to left: (cumulative)
	# 17 bits for feature number    (0)
	# 4 bits for type		(17)
	# 8 bits for ext		(21)
	# Rest for genome 		(29)
	#
	# New encoding; we ran out of bits in ext
	# 17 bits for feature number    (0)
	# 4 bits for type		(17)
	# 15 bits for ext		(21)
	# Rest for genome 		(36)

	#
	# Encoding for SEED using large values of
	# genome extension (e.g. mcSEED with 6666666.X genomes)
	#
	# 17 bits for feature number    (0)
	# 4 bits for type		(17)
	# 18 bits for ext		(21)
	# Rest for genome 		(39)
	

	my $enc;

	if ($FIG_Config::fid_to_docid_large_ext_encoding)
	{
	    $enc = $g << 39| $ext << 21 | $tnum << 17 | $num;
	}
	else
	{
	    $enc = $g << 36| $ext << 21 | $tnum << 17 | $num;
	}
	
	return $enc;
    }

    return undef;
}

sub docid_to_fid
{
    my($doc) = @_;

    my($g, $e, $t, $n);
    
    if ($FIG_Config::fid_to_docid_large_ext_encoding)
    {
	$g = $doc >> 39;
	$e = ($doc >> 21) & 0x3ffff;
	$t = ($doc >> 17) & 0xf;
	$n = $doc & 0x1ffff;
    }
    else
    {
	$g = $doc >> 36;
	$e = ($doc >> 21) & 0x7fff;
	$t = ($doc >> 17) & 0xf;
	$n = $doc & 0x1ffff;
    }
    
    my $type = $tmap{$t};
    my $genome = "$g.$e";
    my $fid = "fig|$genome.$type.$n";

    return $fid;
}


1;
