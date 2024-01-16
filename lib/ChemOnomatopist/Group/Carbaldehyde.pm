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

sub is_part_of_chain() { 1 }

sub prefix() { 'formyl' }
sub suffix() { 'carbaldehyde' }

1;
