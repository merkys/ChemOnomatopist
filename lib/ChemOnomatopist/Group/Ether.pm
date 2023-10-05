package ChemOnomatopist::Group::Ether;

use strict;
use warnings;

# ABSTRACT: Ether group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'O' }

sub prefix() { return 'oxy' }
sub suffix() { return 'oxy' }

1;
