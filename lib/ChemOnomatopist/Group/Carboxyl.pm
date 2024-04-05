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
        $name->append_stem( $_ );
    }

    return $name;
}

sub prefix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};

    my $name = element_suffix( ChemOnomatopist::element( $hydroxy ), $ketone->element );
    for (@$name) {
        $_->{value} =~ s/o$/yl/;
        $_->{value} =~ s/^selenyl$/selanyl/;
    }
    return $name->append_stem( 'carbonyl' );
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

        my @hydroxy_elements = map { ChemOnomatopist::element( $_ ) }
                                   @{$hydroxy->{atoms}};
        my $hydroxy_part = element_suffix( @hydroxy_elements );
        my $ketone_part  = element_suffix( $ketone->element );

        $hydroxy_part .= $ketone->element eq 'O' ? 'peroxoic' : 'peroxo';
        if( any { $_ ne 'O' } @hydroxy_elements ) {
            $hydroxy_part->bracket;
        }
        $ketone_part .= 'ic' if $ketone->element ne 'O';

        $name .= $hydroxy_part;
        $name .= $ketone_part;

        local $" = '';
        return $name . " @hydroxy_elements-acid" if any { $_ ne 'O' } @hydroxy_elements;
        return $name . ' acid';
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
