package ChemOnomatopist::Group::Carboxyl;

use strict;
use warnings;

# ABSTRACT: Carboxyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_carbon { return 1 }

sub prefix { return 'carboxy' }
sub suffix { return 'oic acid' }
sub multisuffix { return 'carboxylic acid' }

1;
