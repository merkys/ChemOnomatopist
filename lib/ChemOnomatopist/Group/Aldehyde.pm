package ChemOnomatopist::Group::Aldehyde;

use strict;
use warnings;

# ABSTRACT: Aldehyde group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_carbon { return 1 }

sub is_part_of_chain() { return 1 }

sub prefix { return 'formyl' }
sub suffix { return 'al' }
sub multisuffix { return 'carbaldehyde' }
sub suffix_if_cycle_substituent { return 'carbaldehyde' }

1;
