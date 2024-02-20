package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name;
use Scalar::Util qw( blessed );

sub new()
{
    my( $class, $hydroxy, $ketone ) = @_;
    return bless { hydroxy => $hydroxy, ketone => $ketone }, $class;
}

sub element() { 'C' }

sub prefix()
{
    my( $self ) = @_;
    my $ketone = $self->{ketone};

    my $name = ChemOnomatopist::Name->new;
    if( $ketone->element ne 'O' ) {
        my $element_prefix = $elements{$ketone->element}->{prefix};
        $element_prefix =~ s/a$/o/;
        $name->append_element( $element_prefix );
    }
    return $name->append_stem( 'carboxy' );
}

sub suffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};
    return ChemOnomatopist::Name->new( 'oic acid' ) if $ketone->element eq 'O';

    my $element_prefix = $elements{$ketone->element}->{prefix};
    $element_prefix =~ s/a$/o/;

    my $name = ChemOnomatopist::Name->new;
    $name->append_multiplier( 'di' ) if ChemOnomatopist::element( $hydroxy ) eq $ketone->element;
    $name .= $element_prefix;
    $name .= 'ic acid';
    return $name;
}

sub multisuffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};
    if( $ketone->element eq 'O' ) {
        return ChemOnomatopist::Name->new( 'carboxylic acid' ) if blessed $hydroxy;
        return ChemOnomatopist::Name->new( 'carboxylate' );
    }

    my $element_prefix = $elements{$ketone->element}->{prefix};
    $element_prefix =~ s/a$/o/;

    my $name = ChemOnomatopist::Name->new( 'carbo' );
    $name->append_multiplier( 'di' ) if $hydroxy->element eq $ketone->element;
    $name .= $element_prefix;
    $name .= 'ic acid';
    return $name;
}

sub suffix_if_cycle_substituent() { $_[0]->multisuffix }

# CHECKME: Not sure if the order is correct here
sub _cmp_instances($$)
{
    my( $A, $B ) = @_;
    return $A->{ketone}->element cmp $B->{ketone}->element ||
           ChemOnomatopist::element( $A->{hydroxy} ) cmp ChemOnomatopist::element( $B->{hydroxy} );
}

1;
