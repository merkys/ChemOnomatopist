package ChemOnomatopist::Group::Sulfonyl;

use strict;
use warnings;

# ABSTRACT: Sulfonyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub element { return 'S' }

sub prefix { return 'sulfonyl' }

sub is_prefix_only() { return 1 }

1;
