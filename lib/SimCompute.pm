package SimCompute;
#
# Little dancer app to allow access to sim compute stuff from a cluster.
#
# See standalone-seed-sims.pl for the interface
#

use Dancer2;
use strict;
use FIG_Config;

get '/:job/:id/query' => sub {
    my $job = params->{job};
    my $id = params->{id};
    my $path = "$FIG_Config::sims_compute_path/$job/sims.job/sims.in/seqs.added/in.$id";
    print STDERR "path=$path\n";
    return send_file($path, system_path => 1, streaming => 1);
};

get '/:job/NR' => sub {
    my $job = params->{job};
    my $path = "$FIG_Config::sims_compute_path/$job/nr";
    send_file($path, system_path => 1, streaming => 1);
};

post '/results' => sub {
    my $job = params->{job};
    my $task = params->{task};
    my $key = params->{key};
    my $path;
    if ($key eq 'output')
    {
	$path = "$FIG_Config::sims_compute_path/$job/sims.job/sims.raw/seqs.added/out.$task";
    }
    elsif ($key eq 'error')
    {
	$path = "$FIG_Config::sims_compute_path/$job/sims.job/sims.err/seqs.added/err.$task";
    }
    else
    {
	send_error("Bad path", 404);
    }
    my $upload = upload('file');
    if ($upload->copy_to($path))
    {
	return "OK\n";
    }
    else
    {
	send_error("bad write", 500);
    }
    
};
dance;
