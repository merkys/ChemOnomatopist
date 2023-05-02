package ChemOnomatopist::Group;

use strict;
use warnings;

use ChemOnomatopist::Group::Carbonyl;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;

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
    ChemOnomatopist::Group::Hydroperoxide::,
    # Amines
    # Imines

    # Classes denoted by the senior atom in heterane nomenclature
);

sub new
{
    my( $class, $carbon ) = @_;
    return bless { C => $carbon }, $class;
}

# Neither of these by default
sub is_carbon { return '' }
sub is_oxygen { return '' }

sub get_name { die "not implemented in the base class\n" }
sub prefix { return '' }
sub suffix { return '' }

# Return the attached carbon
sub C {
    my( $self ) = @_;
    return $self->is_carbon ? $self : $self->{C};
}

1;
