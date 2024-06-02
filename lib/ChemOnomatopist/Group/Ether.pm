package ChemOnomatopist::Group::Ether;

# ABSTRACT: Ether group
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Name;

use parent ChemOnomatopist::Group::;

sub element() { 'O' }

sub is_part_of_chain() { 1 }

sub prefix() { ChemOnomatopist::Name->new( 'oxy' ) }
sub suffix() { ChemOnomatopist::Name->new }

1;
