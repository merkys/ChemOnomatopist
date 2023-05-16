package ChemOnomatopist::Group;

use strict;
use warnings;

use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amino;
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
    ChemOnomatopist::Group::Aldehyde::,
    ChemOnomatopist::Group::Carbonyl::,
    ChemOnomatopist::Group::Hydroxy::,
    ChemOnomatopist::Group::Hydroperoxide::,
    ChemOnomatopist::Group::Amino::,
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

sub prefix { return '' }
sub suffix { return '' }
sub multisuffix { return $_[0]->suffix }

# Return the attached carbon
sub C {
    my( $self ) = @_;
    return $self->is_carbon ? $self : $self->{C};
}

1;
