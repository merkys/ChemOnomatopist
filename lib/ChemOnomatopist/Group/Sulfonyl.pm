package ChemOnomatopist::Group::Sulfonyl;

# ABSTRACT: Sulfonyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

my %prefixes = (
    S  => 'sulfonyl',
    Se => 'selenonyl',
    Te => 'telluronyl',
);

sub prefix() { ChemOnomatopist::Name->new( $prefixes{$_[0]->element} ) }

sub is_prefix_only() { 1 }

1;
