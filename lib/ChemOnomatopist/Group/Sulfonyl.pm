package ChemOnomatopist::Group::Sulfonyl;

use strict;
use warnings;

# ABSTRACT: Sulfonyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub element() { return $_[0]->C->{symbol} }

my %prefixes = (
    S  => 'sulfonyl',
    Se => 'selenonyl',
    Te => 'telluronyl',
);

sub prefix { return $prefixes{$_[0]->element} }

sub is_prefix_only() { return 1 }

1;
