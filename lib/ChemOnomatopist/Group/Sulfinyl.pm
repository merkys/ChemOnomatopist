package ChemOnomatopist::Group::Sulfinyl;

# ABSTRACT: Sulfinyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;

my %prefixes = (
    S  => 'sulfinyl',
    Se => 'seleninyl',
    Te => 'tellurinyl',
);

sub new
{
    my( $class, $element, $ketone ) = @_;
    return bless { element => $element, ketone => $ketone }, $class;
}

sub is_prefix_only() { 1 }

sub prefix { ChemOnomatopist::Name->new( $prefixes{$_[0]->element} ) }

1;
