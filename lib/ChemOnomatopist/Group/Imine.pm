package ChemOnomatopist::Group::Imine;

# ABSTRACT: Imine group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

sub element() { 'N' }

sub is_terminal() { 1 }

sub needs_multiple_bond_suffix() { '' }

sub prefix() { ChemOnomatopist::Name->new( 'imino' ) }
sub suffix() { ChemOnomatopist::Name->new( 'imine' ) }

1;
