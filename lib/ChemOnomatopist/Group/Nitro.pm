package ChemOnomatopist::Group::Nitro;

use strict;
use warnings;

# ABSTRACT: Nitro group, or its analogue
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'nitro' }

sub is_prefix_only() { return 1 }

1;
