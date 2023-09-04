package ChemOnomatopist::Group::Sulfinyl;

use strict;
use warnings;

# ABSTRACT: Sulfinyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub element { return 'S' }

sub prefix { return 'sulfinyl' }

sub is_prefix_only() { return 1 }

1;
