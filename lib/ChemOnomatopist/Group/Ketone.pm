package ChemOnomatopist::Group::Ketone;

# ABSTRACT: Ketone group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

# From BBv2 P-64.6.1
my %prefixes = ( O => 'oxo', S => 'sulfanylidene', Se => 'selanylidene', Te => 'tellanylidene' );
my %suffixes = ( O => 'one', S => 'thione', Se => 'selone', Te => 'tellone' );

sub needs_multiple_bond_suffix { '' }

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
