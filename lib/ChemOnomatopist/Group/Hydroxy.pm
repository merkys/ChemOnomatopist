package ChemOnomatopist::Group::Hydroxy;

use strict;
use warnings;

# ABSTRACT: Hydroxy group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_oxygen() { return 1 }

sub prefix() { return 'hydroxy' }
sub suffix() { return 'ol' }

1;
