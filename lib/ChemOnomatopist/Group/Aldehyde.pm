package ChemOnomatopist::Group::Aldehyde;

use strict;
use warnings;

# ABSTRACT: Aldehyde group
# VERSION

use parent ChemOnomatopist::Group::;

sub is_carbon { return 1 }

sub prefix { return 'formyl' }
sub suffix { return 'al' }
sub multisuffix { return 'carbaldehyde' }

1;
