package ChemOnomatopist::Group;

use strict;
use warnings;

use ChemOnomatopist::Group::Carbonyl;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Hydroxy;
use Scalar::Util qw( blessed );

# ABSTRACT: Chemical group
# VERSION

# Order from BBv2 P-41
our @order = (
    # Radicals
    # Radical anions
    # Radical cations
    # Anions
    # Zwitterions
    # Cations
    # Acids
    ChemOnomatopist::Group::Carboxyl::,
    # Anhydrides
    # Esters
    # Acid halides and pseudohalides
    # Amides
    # Hydrazides
    # Imides
    # Nitriles
    # Aldehydes and chalcogen analogues
    ChemOnomatopist::Group::Carbonyl::,
    ChemOnomatopist::Group::Hydroxy::,
    # Hydroperoxides
    # Amines
    # Imines

    # Classes denoted by the senior atom in heterane nomenclature
);

sub new
{
    my( $class, $self ) = @_;
    $self = {} unless $self;
    return bless $self, $class;
}

# Neither of these by default
sub is_carbon { return '' }
sub is_oxygen { return '' }

sub get_name
{
    my( $class, $graph ) = @_;

    my @groups = grep { blessed $_ && $_->isa( $class ) } $graph->vertices;

    die "cannot handle multiple instances of the same group\n" if @groups > 1;
    return $class->get_mainchain_name( $graph );
}

1;
