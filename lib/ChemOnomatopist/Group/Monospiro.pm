package ChemOnomatopist::Group::Monospiro;

use strict;
use warnings;

# ABSTRACT: Monospiro compound
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use ChemOnomatopist;
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

    # "Numbering starts in the smaller ring, if one is smaller, at a ring atom next to the spiro atom and proceeds first around that ring, then through the spiro atom and around the second ring."
    my( $A, $B ) = @components;
    my @chains;
    push @chains,
         ChemOnomatopist::Chain->new( $graph,
                                      undef,
                                      @$A, $spiro_atom, @$B ),
         ChemOnomatopist::Chain->new( $graph,
                                      undef,
                                      @$A, $spiro_atom, reverse(@$B) ),
         ChemOnomatopist::Chain->new( $graph,
                                      undef,
                                      reverse(@$A), $spiro_atom, @$B ),
         ChemOnomatopist::Chain->new( $graph,
                                      undef,
                                      reverse(@$A), $spiro_atom, reverse(@$B) );

    if( @$A == @$B ) {
        push @chains,
             ChemOnomatopist::Chain->new( $graph,
                                          undef,
                                          @$B, $spiro_atom, @$A ),
             ChemOnomatopist::Chain->new( $graph,
                                          undef,
                                          @$B, $spiro_atom, reverse(@$A) ),
             ChemOnomatopist::Chain->new( $graph,
                                          undef,
                                          reverse(@$B), $spiro_atom, @$A ),
             ChemOnomatopist::Chain->new( $graph,
                                          undef,
                                          reverse(@$B), $spiro_atom, reverse(@$A) );
    }

    my( $chain ) = ChemOnomatopist::filter_chains( @chains );

    return bless { graph => $graph, vertices => [ $chain->vertices ], spiro_atom => $spiro_atom, components => \@components }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    return ( $self );
}

sub components()
{
    my( $self ) = @_;
    return @{$self->{components}};
}

sub prefix()
{
    my( $self ) = @_;
    return 'spiro[' . join( '.', map { scalar @$_ } $self->components ) . ']' .
           ChemOnomatopist::alkane_chain_name( $self->length ) . 'ane';
}

# FIXME: This is a bit strange: class and object method with the same name
sub suffix()
{
    my( $self ) = @_;
    return '' unless ref $self;
    return 'spiro[' . join( '.', map { scalar @$_ } $self->components ) . ']' .
           ChemOnomatopist::unbranched_chain_name( $self );
}

1;
