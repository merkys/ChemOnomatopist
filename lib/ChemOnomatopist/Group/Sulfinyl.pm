package ChemOnomatopist::Group::Sulfinyl;

use strict;
use warnings;

# ABSTRACT: Sulfinyl group
# VERSION

use parent ChemOnomatopist::Group::;

my %prefixes = (
    S  => 'sulfinyl',
    Se => 'seleninyl',
    Te => 'tellurinyl',
);

sub prefix { return $prefixes{$_[0]->element} }

sub is_prefix_only() { return 1 }

1;
