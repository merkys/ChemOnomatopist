package ChemOnomatopist::Group::Monospiro;

use strict;
use warnings;

# ABSTRACT: Monospiro compound
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Chain::VertexArray;
use Graph::Traversal::DFS;

sub new
{
    my( $class, $graph, @vertices ) = @_;

    my $subgraph = $graph->subgraph( \@vertices );
    my( $spiro_atom ) = grep { $subgraph->degree( $_ ) == 4 } @vertices;
    $subgraph->delete_vertex( $spiro_atom );

    # Graph is broken into components.
    # Each component is represented as an array of vertices in the order of traverse.
    my @components;
    for my $component (sort { @$a <=> @$b } $subgraph->connected_components) {
        my( $start ) = sort { $subgraph->degree( $a ) <=> $subgraph->degree( $b ) }
                            @$component;
        push @components,
             [ Graph::Traversal::DFS->new( $subgraph, start => $start )->dfs ];
    }

    return bless { graph => $graph, vertices => \@vertices, spiro_atom => $spiro_atom, components => \@components }, $class;
}

sub candidate_chains
{
    my( $self ) = @_;
    my( $A, $B ) = @{$self->{components}};

    # "Numbering starts in the smaller ring, if one is smaller, at a ring atom next to the spiro atom and proceeds first around that ring, then through the spiro atom and around the second ring."
    my @chains;
    push @chains,
         ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                   @$A, $self->{spiro_atom}, @$B ),
         ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                   @$A, $self->{spiro_atom}, reverse(@$B) ),
         ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                   reverse(@$A), $self->{spiro_atom}, @$B ),
         ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                   reverse(@$A), $self->{spiro_atom}, reverse(@$B) );

    if( @$A == @$B ) {
        push @chains,
             ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                       @$B, $self->{spiro_atom}, @$A ),
             ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                       @$B, $self->{spiro_atom}, reverse(@$A) ),
             ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                       reverse(@$B), $self->{spiro_atom}, @$A ),
             ChemOnomatopist::Chain::VertexArray->new( $self->{graph},
                                                       reverse(@$B), $self->{spiro_atom}, reverse(@$A) );
    }

    return @chains;
}

1;
