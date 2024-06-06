package ChemOnomatopist::Group::Ester;

# ABSTRACT: Ester group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use ChemOnomatopist;

sub new
{
    my( $class, $graph, $C, $ketone, $O ) = @_;
    die "cannot name esters yet\n" if $ChemOnomatopist::CAUTIOUS;
    return bless { graph => $graph, vertices => [ $C, $O ] }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }

sub needs_substituent_locants() { '' }

sub prefix()
{
    my( $self ) = @_;

    if( $self->parent && $self->graph->has_edge( $self->parent, $self->{vertices}[0] ) ) {
        die "cannot handle complicated esters\n";
    }

    return 'oxy';
}

sub suffix() { 'anoate' }

1;
