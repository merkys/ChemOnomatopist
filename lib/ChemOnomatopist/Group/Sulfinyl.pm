package ChemOnomatopist::Group::Sulfinyl;

use strict;
use warnings;

# ABSTRACT: Sulfinyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub element { return 'S' }

my %prefixes = (
    S  => 'sulfinyl',
    Se => 'seleninyl',
    Te => 'tellurinyl',
);

sub prefix
{
    my( $self ) = @_;
    return $prefixes{$self->C->{symbol}};
}

sub is_prefix_only() { return 1 }

1;
