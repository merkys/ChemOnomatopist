package ChemOnomatopist::Group::Carboxyl;

use strict;
use warnings;

# ABSTRACT: Carboxyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_carbon() { return 1 }

sub is_part_of_chain() { return '' }

sub prefix() { return 'carboxy' }
sub suffix() { return 'oic acid' }
sub multisuffix() { return 'carboxylic acid' }

1;
