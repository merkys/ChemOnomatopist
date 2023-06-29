package ChemOnomatopist::Group::Hydroperoxide;

use strict;
use warnings;

# ABSTRACT: Hydroperoxide group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return 'O' }

sub prefix { return 'hydroperoxy' }
sub suffix { return 'peroxol' }

1;
