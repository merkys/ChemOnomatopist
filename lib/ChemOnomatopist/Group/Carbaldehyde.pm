package ChemOnomatopist::Group::Carbaldehyde;

# ABSTRACT: Carbaldehyde group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );

sub new
{
    my( $class, $aldehyde ) = @_;
    return bless { ketone => $aldehyde->{ketone} }, $class;
}

sub element() { 'C' }

sub prefix()
{
    my( $self ) = @_;
    return 'formyl' if $self->{ketone}->element eq 'O';

    my $name = 'methane' . $elements{$self->{ketone}->element}->{prefix};
    $name =~ s/a$//;
    $name .= 'oyl';
    return $name;
}

my %suffixes = ( O => '', S => 'othi', Se => 'oselen', Te => 'otellan' );

sub suffix() { 'carb' . $suffixes{$_[0]->{ketone}->element} . 'aldehyde' }

1;
