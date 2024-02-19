package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );

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
    return 'oic acid' if $ketone->element eq 'O';

    my $name = '';
    $name .= 'di' if $hydroxy->element eq $ketone->element;
    $name .= $elements{$ketone->element}->{prefix};
    $name =~ s/a$/oic acid/;
    return $name;
}

sub multisuffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};
    return 'carboxylic acid' if $ketone->element eq 'O';

    my $name = 'carbo';
    $name .= 'di' if $hydroxy->element eq $ketone->element;
    $name .= $elements{$ketone->element}->{prefix};
    $name =~ s/a$/oic acid/;
    return $name;
}

sub suffix_if_cycle_substituent() { $_[0]->multisuffix }

1;
