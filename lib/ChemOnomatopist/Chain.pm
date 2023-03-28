package ChemOnomatopist::Chain;

use strict;
use warnings;

# ABSTRACT: Longest chain in a compound
# VERSION

use List::Util qw( sum0 );

sub AUTOLOAD {
	our $AUTOLOAD;
	my $call = $AUTOLOAD;
	$call =~ s/.*:://;
	return if $call eq 'DESTROY';
	if( $call =~ /^number_/ ) {
		return sum0 map { $_->can( $call )->() } @{$_[0]->{halves}};
	} else {
		return;
	}
}

sub new
{
	my( $class, @halves ) = @_;
	return bless $class, { halves => \@halves };
}

1;
