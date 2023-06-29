package ChemOnomatopist::Group::Carboxyl;

use strict;
use warnings;

# ABSTRACT: Carboxyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'C' }

sub prefix() { return 'carboxy' }
sub suffix() { return 'oic acid' } # FIXME: Should be 'carboxylic acid' if attached to cycles
sub multisuffix() { return 'carboxylic acid' }
sub suffix_if_cycle_substituent() { return $_[0]->multisuffix }

1;
