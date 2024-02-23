package ChemOnomatopist::Group::Hydroxy;

# ABSTRACT: Hydroxy group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $atom ) = @_;
    return bless { atom => $atom }, $class;
}

sub element() { $_[0]->{atom}->{symbol} }

# From BBv2 P-63.1.5
my %prefixes = ( O => 'hydroxy', S => 'sulfanyl', Se => 'selanyl', Te => 'tellanyl' );
my %suffixes = ( O => 'ol', S => 'thiol', Se => 'selenol', Te => 'tellurol' );

sub prefix
{
    my( $self ) = @_;
    return $prefixes{$self->element};
}

sub suffix
{
    my( $self ) = @_;

    my $suffix = '';
    if( exists $self->{atom}->{isotope} ) {
        $suffix = '(' . $self->{atom}->{isotope} . $self->element . ')';
    }

    return $suffix . $suffixes{$self->element};
}

sub _cmp_instances
{
    my( $A, $B ) = @_;
    return $A->element cmp $B->element
}

1;
