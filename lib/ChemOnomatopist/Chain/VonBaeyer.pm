package ChemOnomatopist::Chain::VonBaeyer;

# ABSTRACT: Von Baeyer hydrocarbon
# VERSION

use strict;
use warnings;

use Graph::Traversal::DFS;
use List::Util qw( first sum );

use parent ChemOnomatopist::Chain::Circular::;

sub new
{
    my( $class, $graph, @vertices ) = @_;

    my $subgraph = $graph->subgraph( \@vertices );
    my @d3 = grep { $subgraph->degree( $_ ) == 3 } @vertices;

    $subgraph->delete_vertices( @d3 );
    my @components = sort { @$b <=> @$a } $subgraph->connected_components;

    $subgraph = $graph->subgraph( \@vertices );
    my $first_of_bridge = first { $subgraph->has_edge( $d3[0], $_ ) }
                                 @{$components[-1]};
    my $last_of_bridge  = first { $subgraph->has_edge( $d3[1], $_ ) }
                                 @{$components[-1]};
    # Disconnect the main bridge
    $subgraph->delete_edge( $d3[0], $first_of_bridge );
    $subgraph->delete_edge( $d3[1], $last_of_bridge );

    # Find the last atom of the main ring
    my $last_of_the_main = first { $subgraph->has_edge( $d3[0], $_ ) }
                                 @{$components[1]};
    $subgraph->delete_edge( $d3[0], $last_of_the_main ); # Disconnect
    # Connect in order to get the correct numbering
    $subgraph->add_edge( $last_of_the_main, $first_of_bridge );

    @vertices = reverse Graph::Traversal::DFS->new( $subgraph, start => $d3[0] )->dfs;

    return bless { graph => $graph,
                   vertices => \@vertices,
                   sizes => [ map { scalar @$_ } @components ] }, $class;
}

sub candidates()
{
    my( $self ) = @_;

    my @candidates = ( $self, $self->flipped );

    my @sizes = @{$self->{sizes}};
    if( $sizes[0] == $sizes[1] && $sizes[0] == $sizes[2] ) {
        # TODO
    } elsif( $sizes[0] == $sizes[1] ) {
        #~ push @candidates, $self->cycles_swapped( 0, 1 );
        #~ push @candidates, $self->cycles_swapped( 0, 1 )->flipped;
    } elsif( $sizes[1] == $sizes[2] ) {
        #~ push @candidates, $self->cycles_swapped( 1, 2 );
        #~ push @candidates, $self->cycles_swapped( 1, 2 )->flipped;
    }

    for (1..$#candidates) {
        $candidates[$_]->{candidate_for} = $self;
    }

    return @candidates;
}

sub flipped()
{
    my( $self ) = @_;

    my $graph = $self->graph;
    my @vertices = $self->vertices;
    my $subgraph = $graph->subgraph( \@vertices );
    my @sizes = @{$self->{sizes}};

    my @d3 = grep { $subgraph->degree( $_ ) == 3 } @vertices;
    @d3 = reverse @d3 unless $d3[0] == $vertices[0]; # CHECKME: Is this needed?

    my @vertices_now = ( $d3[1] );
    shift @vertices;
    push @vertices_now, reverse splice @vertices, 0, $sizes[0];
    shift @vertices;
    push @vertices_now, $d3[0];
    push @vertices_now, reverse splice @vertices, 0, $sizes[1];
    push @vertices_now, reverse splice @vertices, 0, $sizes[2];

    return bless { graph => $graph,
                   vertices => \@vertices_now,
                   sizes => \@sizes };
}

sub cycles_swapped($$)
{
    my( $self, $A, $B ) = @_;

    my $graph = $self->graph;
    my @vertices = $self->vertices;
    my $subgraph = $graph->subgraph( \@vertices );
    my @sizes = @{$self->{sizes}};

    my @d2 = grep { $subgraph->degree( $_ ) == 2 } @vertices;
    my @d3 = grep { $subgraph->degree( $_ ) == 3 } @vertices;
    @d3 = reverse @d3 unless $d3[0] == $vertices[0]; # CHECKME: Is this needed?

    my @A = map { $d2[$_] } (sum @sizes[0..$A-1])..(sum @sizes[0..$A])-1;
    my @B = map { $d2[$_] } (sum @sizes[0..$B-1])..(sum @sizes[0..$B])-1;

    @vertices = @d2;

    splice @vertices, sum @sizes[0..$A-1], $sizes[$A], @B;
    splice @vertices, sum @sizes[0..$B-1], $sizes[$B], @A;

    splice  @vertices, $sizes[0], $d3[1];
    unshift @vertices, $d3[0];

    return bless { graph => $graph,
                   vertices => \@vertices,
                   sizes => \@sizes };
}

sub has_form($$)
{
    my( $class, $graph ) = @_;

    my @d2 = grep { $graph->degree( $_ ) == 2 } $graph->vertices;
    my @d3 = grep { $graph->degree( $_ ) == 3 } $graph->vertices;

    return '' unless @d3 == 2;
    return '' unless @d2 + @d3 == scalar $graph->vertices;

    return '' unless $graph->is_edge_connected; # Must not have bridges
    return '' if $graph->has_edge( @d3 ); # Reject regular bicycles

    $graph = $graph->copy->delete_vertices( @d3 );
    return scalar( $graph->connected_components ) == 3;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ + 1 } @_;
}

sub prefix() { &suffix }
sub suffix()
{
    my( $self ) = @_;
    return 'bicyclo[' . join( '.', @{$self->{sizes}} ) . ']' .
           ChemOnomatopist::alkane_chain_name( $self->length ) .
           'ane';
}

1;
