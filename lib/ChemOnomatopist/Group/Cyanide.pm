package ChemOnomatopist::Group::Cyanide;

# ABSTRACT: Cyanide group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Name;

use parent ChemOnomatopist::Group::;

sub prefix() { ChemOnomatopist::Name->new( 'cyano' ) }
sub suffix() { ChemOnomatopist::Name->new( 'nitrile' ) }
sub multisuffix() { ChemOnomatopist::Name->new( 'carbonitrile' ) }

1;
