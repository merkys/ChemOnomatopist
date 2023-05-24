package ChemOnomatopist::ChainHalf;

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Util::SMILES qw( path_SMILES );
use Graph::Traversal::DFS;
use List::Util qw( sum0 );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

# ABSTRACT: Half of a longest chain
# VERSION

sub vertices();

sub new
{
    my( $class, $graph, $other_center, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph, other_center => $other_center, cache => {} };
    return bless $self, $class;
}

# Accessors

# Groups are used to check which halves of a chain can be combined together.
# If a graph contains single center, all halves will share the center.
sub group()
{
    my( $self ) = @_;
    return $self->{vertices}[1 - defined $self->{other_center}];
}

sub substituents()
{
    my( $self ) = @_;

    my $vertices = set( $self->vertices );

    my @substituents;
    for my $vertex ($self->vertices) {
        for my $neighbour ($self->{graph}->neighbours( $vertex )) {
            next if $vertices->has( $neighbour );
            push @substituents, $neighbour;
        }
    }

    return @substituents;
}

sub vertices()
{
    my( $self ) = @_;
    return @{$self->{vertices}};
}

# Properties

sub backbone_SMILES()
{
    my( $self ) = @_;
    return path_SMILES( $self->{graph}, $self->vertices );
}

sub branch_positions()
{
    my( $self ) = @_;

    return @{$self->{branch_positions}} if $self->{branch_positions};

    my $graph = $self->_disconnected_chain_graph;
    my @vertices = $self->vertices;

    my @branch_positions =
        map { ( $_ ) x $graph->degree( $vertices[$_] ) }
            grep { $graph->degree( $vertices[$_] ) }
                 0..$#vertices;

    $self->{branch_positions} = \@branch_positions;
    return @branch_positions;
}

sub group_positions
{
    my( $self, $class ) = @_;

    return @{$self->{group_positions}{$class}} if $self->{group_positions}{$class};

    my $graph = $self->_disconnected_chain_graph;
    my @vertices = $self->vertices;

    my @group_positions;
    for (0..$#vertices) {
        my $groups = grep { blessed $_ && $_->isa( $class ) }
                          $graph->neighbours( $vertices[$_] );
        next unless $groups;
        push @group_positions, ( $_ ) x $groups;
    }

    $self->{group_positions}{$class} = \@group_positions;
    return @group_positions;
}

sub heteroatom_positions()
{
    my( $self ) = @_;

    return @{$self->{heteroatom_positions}} if $self->{heteroatom_positions};

    my @vertices = $self->vertices;
    my @heteroatom_positions;
    for (0..$#vertices) {
        next if blessed $vertices[$_];
        next if ChemOnomatopist::is_element( $vertices[$_], 'C' );
        push @heteroatom_positions, $_;
    }

    $self->{heteroatom_positions} = \@heteroatom_positions;
    return @heteroatom_positions;
}

sub most_senior_group_positions()
{
    my( $self ) = @_;

    return @{$self->{most_senior_group_positions}} if $self->{most_senior_group_positions};

    my $class = ChemOnomatopist::most_senior_group( $self->vertices, $self->substituents );
    return () unless $class;

    my @vertices = $self->vertices;
    my $vertices = set( @vertices );
    my @positions;
    for (0..$#vertices) {
        my $vertex = $vertices[$_];
        push @positions, $_ if blessed $vertex && $vertex->isa( $class );
        for my $neighbour ($self->{graph}->neighbours( $vertex )) {
            next if $vertices->has( $neighbour );
            next unless blessed $neighbour;
            next unless $neighbour->isa( $class );
            push @positions, $_;
        }
    }

    $self->{most_senior_group_positions} = \@positions;
    return @positions;
}

sub heteroatoms()
{
    my( $self ) = @_;
    return map { $self->{vertices}[$_]{symbol} } $self->heteroatom_positions;
}

sub length()
{
    my( $self ) = @_;
    return scalar $self->vertices;
}

sub locant_names()
{
    my( $self ) = @_;

    return @{$self->{locant_names}} if $self->{locant_names};

    my $graph = $self->_disconnected_chain_graph->copy;

    my @locants;
    for my $vertex ($self->vertices) {
        my @current_locants;
        for my $neighbour ($graph->neighbours( $vertex )) {
            $graph->delete_edge( $vertex, $neighbour );
            if( ChemOnomatopist::is_element( $neighbour, 'C' ) ) {
                push @current_locants, ChemOnomatopist::get_sidechain_name( $graph, $neighbour );
            } elsif( blessed $neighbour ) {
                push @current_locants, $neighbour->prefix;
            }
        }
        push @locants, \@current_locants;
    }

    $self->{locant_names} = \@locants;
    return @locants;
}

sub number_of_branches_in_sidechains()
{
    my( $self ) = @_;

    return $self->{number_of_branches_in_sidechains} if exists $self->{number_of_branches_in_sidechains};

    my $graph = $self->_disconnected_chain_graph->copy;
    my @vertex_neighbours = map { $graph->neighbours( $_ ) } $self->vertices;
    $graph->delete_vertices( $self->vertices );

    my $number = sum0 map { $_ > 2 ? $_ - 2 : 0 }
                          map { $graph->degree( $_ ) }
                          map { Graph::Traversal::DFS->new( $graph, start => $_ )->dfs }
                              @vertex_neighbours;

    $self->{number_of_branches_in_sidechains} = $number;
    return $number;
}

sub number_of_carbons()
{
    my( $self ) = @_;

    return $self->{number_of_carbons} if exists $self->{number_of_carbons};

    my $graph = $self->_disconnected_chain_graph;

    my $C = grep { ChemOnomatopist::is_element( $_, 'C' ) }
            map  { Graph::Traversal::DFS->new( $graph, start => $_ )->dfs }
                 $self->vertices;

    # Since main chain carbons are included in the count, they have to be subtracted.
    $C -= $self->length;

    $self->{number_of_carbons} = $C;
    return $C;
}

sub number_of_branches()
{
    my( $self ) = @_;
    return scalar $self->branch_positions;
}

sub number_of_groups
{
    my( $self, $class ) = @_;
    return scalar $self->group_positions( $class );
}

sub number_of_heteroatoms()
{
    my( $self ) = @_;
    return scalar $self->heteroatom_positions;
}

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    return $self->{_disconnected_chain_graph} if $self->{_disconnected_chain_graph};

    my $graph = $self->{graph}->copy;
    my @vertices = $self->vertices;

    if( $self->{other_center} ) {
        # Cut the edge to the other center
        $graph->delete_edge( $vertices[0], $self->{other_center} );
    } else {
        # Cut the edges to the other candidates
        for ($graph->neighbours( $vertices[0] )) {
            $graph->delete_edge( $vertices[0], $_ );
        }
    }
    $graph->delete_path( @vertices );

    $self->{_disconnected_chain_graph} = $graph;
    return $graph;
}

1;
