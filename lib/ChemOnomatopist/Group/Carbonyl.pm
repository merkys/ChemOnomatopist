package ChemOnomatopist::Group::Carbonyl;

use strict;
use warnings;

# ABSTRACT: Carbonyl group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $carbon, $atom ) = @_;
    return bless { C => $carbon, atom => $atom }, $class;
}

sub is_oxygen() {
    my( $self ) = @_;
    return ChemOnomatopist::is_element( $self->{atom}, 'O' );
}

my %prefixes = ( O => 'oxo', S => 'sulfanylidene' );
my %suffixes = ( O => 'one', S => 'thione' );

sub prefix
{
    my( $self ) = @_;
    return $prefixes{$self->{atom}{symbol}};
}

sub suffix
{
    my( $self ) = @_;
    return $suffixes{$self->{atom}{symbol}};
}

sub _cmp_instances
{
    my( $A, $B ) = @_;
    return $A->{atom}{symbol} cmp $B->{atom}{symbol};
}

1;
