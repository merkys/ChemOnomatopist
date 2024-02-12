package ChemOnomatopist::Group::Amidine;

# ABSTRACT: Amidine group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }
sub needs_substituent_locants { '' }

sub prefix() { 'carbamimidoyl' }
sub suffix() { 'imidamide' }

1;
