package ChemOnomatopist::Name;

use strict;
use warnings;

# ABSTRACT: Chemical name
# VERSION

sub new
{
    my( $class, $name ) = @_;
    $name = '' unless $name;
    return bless { name => $name }, $class;
}

sub append($)
{
    my( $self, $string ) = @_;
    $self->{name} .= $string;
}

1;
