package ChemOnomatopist::Group::Hydroxy;

use strict;
use warnings;

# ABSTRACT: Hydroxy group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_oxygen() { return 1 }

sub is_part_of_chain() { return '' }

sub prefix() { return 'hydroxy' }
sub suffix() { return 'ol' }

1;
