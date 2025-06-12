package ChemOnomatopist::Group::Peroxide;

# ABSTRACT: Peroxide group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Name;
use List::Util qw( any );
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $graph, @atoms ) = @_;
    return bless { graph => $graph, atoms => \@atoms }, $class;
}

sub element() { 'O' }

sub prefix()
{
    my( $self ) = @_;
    if( any  { $_->isa( ChemOnomatopist::Group::Hydroxy:: ) ||
               $_->isa( ChemOnomatopist::Group::Hydroperoxide:: ) }
        grep { blessed $_ } $self->{graph}->neighbours( $self ) ) {
        return ChemOnomatopist::Name->new( 'dioxidane' );
    } else {
        return ChemOnomatopist::Name->new( 'peroxy' );
    }
}

sub suffix() { ChemOnomatopist::Name->new( 'peroxolate' ) }

1;
