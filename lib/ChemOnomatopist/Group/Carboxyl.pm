package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Name;
use List::Util qw( all any uniq );
use Scalar::Util qw( blessed );

sub new()
{
    my( $class, $hydroxy, $ketone ) = @_;
    return bless { hydroxy => $hydroxy, ketone => $ketone }, $class;
}

sub element() { 'C' }

sub element_suffix(@)
{
    my @elements = @_;

    my @prefixes = sort map  { s/a$/o/; $_ }
                        map  { $elements{$_}->{prefix} }
                        grep { $_ ne 'O' }
                             @elements;

    my $name = ChemOnomatopist::Name->new;
    if( @prefixes == 2 && uniq( @elements ) == 1 ) {
        $name->append_multiplier( 'di' );
        shift @prefixes;
    }
    for (@prefixes) {
        $name .= $_;
    }

    return $name;
}

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
    my @elements = ( ChemOnomatopist::element( $hydroxy ), $ketone->element );
    if( all { $_ eq 'O' } @elements ) {
        return ChemOnomatopist::Name->new( 'oic acid' ) if blessed $hydroxy;
        return ChemOnomatopist::Name->new( 'oate' );
    }

    return element_suffix( @elements ) . 'ate' unless blessed $hydroxy;
    return element_suffix( @elements ) . 'ic acid' if uniq( @elements ) == 1;
    return element_suffix( @elements ) . ('ic ' . $elements[0] . '-acid');
}

sub multisuffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};

    my $name = ChemOnomatopist::Name->new( 'carbo' );
    if( blessed $hydroxy &&
        $hydroxy->isa( ChemOnomatopist::Group::Hydroperoxide:: ) ) {
        $name .= 'peroxo';
        my @elements = map { ChemOnomatopist::element( $_ ) }
                           @{$hydroxy->{atoms}};
        my $element_prefix = $elements{$ketone->element}->{prefix};
        $element_prefix =~ s/a$/o/;
        $name .= $element_prefix unless $ketone->element eq 'O';

        local $" = '';
        return $name . "ic @elements-acid" if any { $_ ne 'O' } @elements;
        return $name . 'ic acid';
    } else {
        my @elements = ( ChemOnomatopist::element( $hydroxy ), $ketone->element );
        if( all { $_ eq 'O' } @elements ) {
            return ChemOnomatopist::Name->new( 'carboxylic acid' ) if blessed $hydroxy;
            return ChemOnomatopist::Name->new( 'carboxylate' );
        }

        $name .= element_suffix( @elements );
        return $name . 'ate' unless blessed $hydroxy;
        return $name . 'ic acid' if uniq( @elements ) == 1;
        return $name . ('ic ' . $elements[0] . '-acid');
    }
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
