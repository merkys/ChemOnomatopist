package ChemOnomatopist::Group::Hydroxy;

use strict;
use warnings;

# ABSTRACT: Hydroxy group
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

# From BBv2 P-63.1.5
my %prefixes = ( O => 'hydroxy', S => 'sulfanyl', Se => 'selanyl', Te => 'tellanyl' );
my %suffixes = ( O => 'ol', S => 'thiol', Se => 'selenol', Te => 'tellurol' );

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
