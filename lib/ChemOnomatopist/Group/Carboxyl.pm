package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub element() { 'C' }

sub prefix() { 'carboxy' }
sub suffix() { 'oic acid' } # FIXME: Should be 'carboxylic acid' if attached to cycles
sub multisuffix() { 'carboxylic acid' }
sub suffix_if_cycle_substituent() { $_[0]->multisuffix }

1;
