package ChemOnomatopist::Group::Diazene;

# ABSTRACT: Diazene group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use ChemOnomatopist::Name;
use ChemOnomatopist::Name::Part::Multiplier;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;

    my @chains = ( $self, ChemOnomatopist::Group::Diazene->new( $self->graph, reverse $self->vertices ) );
    $chains[1]->{candidate_for} = $self;

    return @chains;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_suffix_locant() { $_[0]->number_of_branches != 2 }

sub prefix() { ChemOnomatopist::Name::Part::Multiplier->new( 'di' )->to_name .
               ChemOnomatopist::Name->new( 'azenyl' ) }
sub suffix() { ChemOnomatopist::Name::Part::Multiplier->new( 'di' )->to_name .
               ChemOnomatopist::Name->new( 'azene' ) }

1;
