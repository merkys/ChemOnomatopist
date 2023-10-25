package ChemOnomatopist::Chain::Porphyrin;

# ABSTRACT: Porphyrin compound
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use Graph::Undirected;

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my @vertices = $graph->vertices;

    return '' unless @vertices == 24;
    return '' unless (grep { ChemOnomatopist::is_element( $_, 'C' ) } @vertices) == 20;
    return '' unless (grep { ChemOnomatopist::is_element( $_, 'N' ) } @vertices) ==  4;
}

sub ideal_graph($)
{
    my( $class ) = @_;
    my $graph = Graph::Undirected->new( refvertexed => 1 );
    my @vertices = map { { symbol => 'C' } } 1..20;
    $graph->add_cycle( @vertices );
    for (0..3) {
        $graph->add_path( $vertices[$_ * 5], { symbol => 'N' }, $vertices[$_ * 5 + 3] );
    }
    return $graph;
}

1;
