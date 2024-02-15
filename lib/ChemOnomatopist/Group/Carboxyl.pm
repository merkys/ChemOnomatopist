package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );

sub new()
{
    my( $class, $ketone ) = @_;
    return bless { ketone => $ketone }, $class;
}

sub element() { 'C' }

sub prefix() { 'carboxy' }
sub suffix() { 'oic acid' }

sub multisuffix()
{
    my( $self ) = @_;
    my $ketone = $self->{ketone};
    return 'carboxylic acid' if $ketone->element eq 'O';

    my $name = 'carbo' . $elements{$self->{ketone}->element}->{prefix};
    $name =~ s/a$/oic acid/;
    return $name;
}

sub suffix_if_cycle_substituent() { $_[0]->multisuffix }

1;
