
package P3RateLimitedUserAgent;

use Time::HiRes qw(gettimeofday usleep);
use Data::Dumper;
use strict;
use base 'LWP::UserAgent';
use URI;

sub new
{
    my($class, $limit, @rest) = @_;

    my $self = LWP::UserAgent->new(@rest);

    $self->{_limit} = $limit;
    $self->{_history} = {};

    return bless $self, $class;
}

sub request
{
    my ($self, $request, $arg, $size, $previous) = @_;
    my $u = URI->new($request->uri);
    my $host = $u->host;
    $self->{_history}->{$host} //= [];
    my $retry = 0;
    my $retry_count = 1;
    my $ret;
    do {
	my $h = $self->{_history}->{$host};
	my $now = gettimeofday;
	my $last = $self->{last_hit};
	$self->{last_hit} = $now;
	printf "ELAP %.1f\n", ($now - $last) *1000;
	$retry = 0;
	
	# print "$host $now $h <@$h> \n";
	
	@$h = grep { $_ > $now - 1 } @$h;
	# print "$self->{_limit}: >> @$h\n";
	if (@$h >= $self->{_limit})
	{
	    my $del = 1 - ($now - $h->[0]);
	    print "del=$del\n";
	    usleep(1e6 * $del) if $del > 0;
	}
	push(@$h, $now);
	$ret = $self->SUPER::request($request, $arg, $size, $previous);
	if (!$ret->is_success && $ret->code == 429)
	{
	    warn $ret->status_line;
	    #my $delay = rand(2) + $retry_count;
	    my $delay = rand(0.5);
	    print STDERR "Retry $retry_count delay=$delay\n";
	    usleep($delay * 1e6);
	    $retry_count++;
	    $retry = 1;
	}
    }
    while ($retry);
    return $ret;
}

1;
