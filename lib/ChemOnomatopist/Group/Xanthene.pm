package ChemOnomatopist::Group::Xanthene;

use strict;
use warnings;

# ABSTRACT: Xanthene or its close derivative
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Util::Graph qw( subgraph );

sub new
{
    my( $class, $graph, @cycles ) = @_;

    my @benzenes = grep {  $_->is_benzene } @cycles;
    my( $other ) = grep { !$_->is_benzene } @cycles;

    # Find the correct vertice order
    my $subgraph = subgraph( $graph, map { $_->vertices } @cycles );
    my @bridges = grep { $subgraph->degree( $_->[0] ) == 3 &&
                         $subgraph->degree( $_->[1] ) == 3 }
                         $subgraph->edges;
    $subgraph->delete_edges( @bridges );

    my( $heteroatom ) = $other->heteroatom_positions;
    my @other_vertices = $other->vertices;
    $subgraph->delete_vertex( $other_vertices[($heteroatom + 2) % 6] );
    my( $start ) = grep { $subgraph->has_vertex( $_ ) && $subgraph->degree( $_ ) == 1 }
                   map  { $_->vertices } @benzenes;
    my @vertices = ( Graph::Traversal::DFS->new( $subgraph, start => $start )->dfs,
                     $other_vertices[($heteroatom + 2) % 6] );

    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidates()
{
    my( $self ) = @_;
    my @candidates = ( $self );
    my @vertices = $self->vertices;

    push @candidates,
         bless { graph => $self->graph, vertices => [ reverse @vertices[11..13], @vertices[0..10] ] };
    $candidates[-1]->{candidate_for} = $self;

    if( $self->number_of_heteroatoms == 2 ) { # TODO: Add two more candidates
    }

    return @candidates;
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }

sub prefix()
{
    my( $self ) = @_;

    my( $heteroatom ) = $self->heteroatoms;
    if( $self->number_of_heteroatoms == 1 ) {
        return 'xanthene' if $heteroatom eq 'O';
        my $name = $elements{$heteroatom}->{prefix};
        $name =~ s/a$/o/;
        return $name . 'xanthene';
    } else {
        my $name = $elements{$heteroatom}->{prefix};
        $name =~ s/a$//;
        return $name . 'anthrene';
    }
}

sub suffix() { return $_[0]->prefix }

1;
