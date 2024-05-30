package ChemOnomatopist::Group::Nitramide;

# ABSTRACT: Nitramide group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub element() { 'N' }

sub prefix() { ChemOnomatopist::Name->new( 'nitramido' ) }
sub suffix() { ChemOnomatopist::Name->new( 'nitramide' ) }

1;
