package ChemOnomatopist::Group::Carboxyl;

# ABSTRACT: Carboxyl group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Name;
use ChemOnomatopist::Util;
use List::Util qw( all any uniq );
use Scalar::Util qw( blessed );

sub new()
{
    my( $class, $hydroxy, $ketone ) = @_;
    die "cannot handle carboxyl group attached to nongroup atoms (most likely hydrazines)\n" unless blessed $ketone;
    die "cannot handle carboxyl group attached to anything else than ketone\n" unless $ketone->isa( ChemOnomatopist::Group::Ketone:: );
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

    my $name = ChemOnomatopist::Name->new;
    if( ChemOnomatopist::Util::element( $hydroxy ) ne 'O' ||
        $ketone->element ne 'O' ) {
        $name .= $hydroxy->prefix;
        $name->bracket unless $name->is_simple;
        $name->append_stem( 'carbono' );
        $name .= $ketone->suffix if $ketone->element ne 'O';
    } else {
        $name->append_stem( 'carboxy' );
    }
    $name->[-1]{value} =~ s/(ne|o)$/yl/;

    return $name;
}

sub suffix()
{
    my( $self ) = @_;
    my $hydroxy = $self->{hydroxy};
    my $ketone = $self->{ketone};
    my @elements = map { ChemOnomatopist::Util::element( $_ ) } ( $hydroxy, $ketone );
    if( all { $_ eq 'O' } @elements ) {
        return ChemOnomatopist::Name->new( 'oic acid' ) if blessed $hydroxy && !$hydroxy->charge;
        return ChemOnomatopist::Name->new( 'oate' );
    }

    return element_suffix( @elements ) . 'ate' unless blessed $hydroxy && !$hydroxy->charge;
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

        my @hydroxy_elements = map { ChemOnomatopist::Util::element( $_ ) } @{$hydroxy->{atoms}};
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
        return $name . " @hydroxy_elements-acid" if scalar( uniq @hydroxy_elements ) > 1;
        return $name . ' acid';
    } else {
        my @elements = ( ChemOnomatopist::Util::element( $hydroxy ), $ketone->element );
        if( all { $_ eq 'O' } @elements ) {
            return ChemOnomatopist::Name->new( 'carboxylic acid' ) if blessed $hydroxy && !$hydroxy->charge;
            return ChemOnomatopist::Name->new( 'carboxylate' );
        }

        $name .= element_suffix( @elements );
        return $name . 'ate' unless blessed $hydroxy && !$hydroxy->charge;
        return $name . 'ic acid' if uniq( @elements ) == 1;
        return $name . ('ic ' . $elements[0] . '-acid');
    }
}

sub suffix_if_cycle_substituent() { $_[0]->multisuffix }

# CHECKME: Not sure if the order is correct here
sub _cmp_instances($$)
{
    my( $A, $B ) = @_;
    return ChemOnomatopist::Util::element( $A->{ketone} )  cmp ChemOnomatopist::Util::element( $B->{ketone} ) ||
           ChemOnomatopist::Util::element( $A->{hydroxy} ) cmp ChemOnomatopist::Util::element( $B->{hydroxy} );
}

1;
