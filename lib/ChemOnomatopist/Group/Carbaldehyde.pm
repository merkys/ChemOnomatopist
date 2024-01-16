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

my %suffixes = ( O => '', S => 'othi', Se => 'oselen', Te => 'otellan' );

sub prefix() { 'formyl' }
sub suffix() { 'carb' . $suffixes{$_[0]->{ketone}->element} . 'aldehyde' }

1;
