package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub new()
{
    my( $class, $ketone ) = @_;
    return bless { ketone => $ketone }, $class;
}

sub element() { 'C' }

sub prefix() { 'carboxy' }
sub suffix() { 'oic acid' }
sub multisuffix() { 'carboxylic acid' }
sub suffix_if_cycle_substituent() { $_[0]->multisuffix }

1;
