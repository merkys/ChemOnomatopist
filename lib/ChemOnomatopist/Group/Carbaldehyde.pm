package ChemOnomatopist::Group::Carbaldehyde;

# ABSTRACT: Carbaldehyde group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $aldehyde ) = @_;
    return bless { ketone => $aldehyde->{ketone} }, $class;
}

sub element() { 'C' }

sub prefix() { 'formyl' }
sub suffix() { 'carbaldehyde' }

1;
