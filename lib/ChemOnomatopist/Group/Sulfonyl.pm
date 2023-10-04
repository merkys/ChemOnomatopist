package ChemOnomatopist::Group::Sulfonyl;

use strict;
use warnings;

# ABSTRACT: Sulfonyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $element ) = @_;
    return bless { element => $element }, $class;
}

sub element() { return $_[0]->{element} }

my %prefixes = (
    S  => 'sulfonyl',
    Se => 'selenonyl',
    Te => 'telluronyl',
);

sub prefix { return $prefixes{$_[0]->element} }

sub is_prefix_only() { return 1 }

1;
