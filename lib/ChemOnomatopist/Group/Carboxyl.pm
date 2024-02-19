package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name;
use ChemOnomatopist::Name::Part::Multiplier;

sub new()
{
    my( $class, $hydroxy, $ketone ) = @_;
    return bless { hydroxy => $hydroxy, ketone => $ketone }, $class;
}

sub element() { 'C' }

sub prefix() { 'carboxy' }
sub suffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};
    return ChemOnomatopist::Name->new( 'oic acid' ) if $ketone->element eq 'O';

    my $element_prefix = $elements{$ketone->element}->{prefix};
    $element_prefix =~ s/a$/oic /;

    my $name = ChemOnomatopist::Name->new;
    $name .= ChemOnomatopist::Name::Part::Multiplier->new( 'di' ) if $hydroxy->element eq $ketone->element;
    $name .= $element_prefix;
    $name .= 'acid';
    return $name;
}

sub multisuffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};
    return ChemOnomatopist::Name->new( 'carboxylic acid' ) if $ketone->element eq 'O';

    my $element_prefix = $elements{$ketone->element}->{prefix};
    $element_prefix =~ s/a$/oic /;

    my $name = ChemOnomatopist::Name->new( 'carbo' );
    $name .= ChemOnomatopist::Name::Part::Multiplier->new( 'di' ) if $hydroxy->element eq $ketone->element;
    $name .= $element_prefix;
    $name .= 'acid';
    return $name;
}

sub suffix_if_cycle_substituent() { $_[0]->multisuffix }

1;
