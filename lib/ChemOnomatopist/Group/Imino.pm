package ChemOnomatopist::Group::Imino;

use strict;
use warnings;

# ABSTRACT: Imino group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'N' }

sub is_terminal() { return 1 }

sub needs_multiple_bond_suffix { return '' }

sub prefix { return 'imino' }
sub suffix { return 'imine' }

1;
