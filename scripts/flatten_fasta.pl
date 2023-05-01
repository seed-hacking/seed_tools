
$/ = "\n>";
while (defined($_ = <STDIN>))
{
    chomp;
    if ($_ =~ /^>?(\S+)[^\n]*\n(.*)/s)
    {
	$id  =  $1;
	$seq =  $2;
	$seq =~ s/\s//gs;
	print "$id\t$seq\n";
    }
}
