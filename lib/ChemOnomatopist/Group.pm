package ChemOnomatopist::Group;

use strict;
use warnings;

use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amino;
use ChemOnomatopist::Group::Carbonyl;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Ester;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Imino;
use ChemOnomatopist::Group::Monocycle;

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
    ChemOnomatopist::Group::Ester::,
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
    ChemOnomatopist::Group::Imino::,

    # Classes denoted by the senior atom in heterane nomenclature

    ChemOnomatopist::Group::Monocycle::, # FIXME: Possibly not the right place
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
