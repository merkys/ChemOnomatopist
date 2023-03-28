package ChemOnomatopist::Chain;

use strict;
use warnings;

# ABSTRACT: Longest chain in a compound
# VERSION

sub new
{
	my( $class, @halves ) = @_;
	return bless $class, { halves => \@halves };
}

1;
