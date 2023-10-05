package ChemOnomatopist::Group::Ether;

use strict;
use warnings;

# ABSTRACT: Ether group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'O' }

sub is_part_of_chain() { return 1 }

sub prefix() { return '' }
sub suffix() { return '' }

1;
