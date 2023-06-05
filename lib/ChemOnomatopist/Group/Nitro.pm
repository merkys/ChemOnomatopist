package ChemOnomatopist::Group::Nitro;

use strict;
use warnings;

# ABSTRACT: Nitro group
# VERSION

use parent ChemOnomatopist::Group::;

sub prefix { return 'nitro' }

sub is_part_of_chain() { return '' }
sub is_prefix_only() { return 1 }

1;
