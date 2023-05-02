package ChemOnomatopist::Group::Hydroperoxide;

use strict;
use warnings;

# ABSTRACT: Hydroperoxide group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_oxygen { return 1 }

sub prefix { return 'hydroperoxy' }
sub suffix { return 'peroxol' }

1;
