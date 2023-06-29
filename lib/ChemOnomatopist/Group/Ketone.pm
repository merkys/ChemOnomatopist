package ChemOnomatopist::Group::Ketone;

use strict;
use warnings;

# ABSTRACT: Ketone group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $carbon, $atom ) = @_;
    return bless { C => $carbon, atom => $atom }, $class;
}

sub element() { return ucfirst $_[0]->{atom}{symbol} }

# From BBv2 P-64.6.1
my %prefixes = ( O => 'oxo', S => 'sulfanylidene', Se => 'selanylidene', Te => 'tellanylidene' );
my %suffixes = ( O => 'one', S => 'thione', Se => 'selone', Te => 'tellone' );

sub prefix
{
    my( $self ) = @_;
    return $prefixes{$self->element};
}

sub suffix
{
    my( $self ) = @_;
    return $suffixes{$self->element};
}

sub _cmp_instances
{
    my( $A, $B ) = @_;
    return $A->element cmp $B->element;
}

1;
