package ChemOnomatopist::Group::Nitro;

# ABSTRACT: Nitro group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub prefix { 'nitro' }

sub is_prefix_only() { 1 }

1;
