package ChemOnomatopist::Group::Imino;

# ABSTRACT: Imino group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub element() { 'N' }

sub is_terminal() { 1 }

sub needs_multiple_bond_suffix() { '' }

sub prefix() { 'imino' }
sub suffix() { 'imine' }

1;
