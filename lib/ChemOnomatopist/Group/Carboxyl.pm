package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Name;
use List::Util qw( all uniq );
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
    my @elements = ( ChemOnomatopist::element( $hydroxy ), $ketone->element );
    if( all { $_ eq 'O' } @elements ) {
        return ChemOnomatopist::Name->new( 'oic acid' ) if blessed $hydroxy;
        return ChemOnomatopist::Name->new( 'oate' );
    }

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
    return $name . 'ate' unless blessed $hydroxy;
    return $name . 'ic acid' if uniq( @elements ) == 1;
    return $name . ('ic ' . $elements[0] . '-acid');
}

sub multisuffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};
    if( $ketone->element eq 'O' ) {
        # FIXME: Element of hydroxy group is not checked here
        return ChemOnomatopist::Name->new( 'carboxylic acid' ) if blessed $hydroxy;
        return ChemOnomatopist::Name->new( 'carboxylate' );
    }

    my $element_prefix = $elements{$ketone->element}->{prefix};
    $element_prefix =~ s/a$/o/;

    my $name = ChemOnomatopist::Name->new( 'carbo' );
    if( $hydroxy->isa( ChemOnomatopist::Group::Hydroperoxide:: ) ) {
        $name .= 'peroxo';
    }
    $name->append_multiplier( 'di' ) if $hydroxy->element eq $ketone->element;
    $name .= $element_prefix;
    $name .= blessed $hydroxy ? 'ic acid' : 'ate';
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
