package ChemOnomatopist::Group::Carbonitrile;

# ABSTRACT: Carbonitrile group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Name;

use parent ChemOnomatopist::Group::;

sub prefix { ChemOnomatopist::Name->new( 'cyano' ) }
sub suffix { ChemOnomatopist::Name->new( 'carbonitrile' ) }

1;
