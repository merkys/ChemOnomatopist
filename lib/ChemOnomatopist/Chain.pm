package ChemOnomatopist::Chain;

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::Util::SMILES qw( path_SMILES );
use Chemistry::OpenSMILES qw( %normal_valence );
use Graph::Traversal::DFS;
use List::Util qw( all sum0 );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

# ABSTRACT: Chain of atoms
# VERSION

sub vertices();

sub new
{
    my( $class, $graph, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph, cache => {} };
    return bless $self, $class;
}

# Accessors

sub graph()
{
    my( $self ) = @_;
    return $self->{graph};
}

sub substituents()
{
    my( $self ) = @_;

    my $vertices = set( $self->vertices );

    my @substituents;
    for my $vertex ($self->vertices) {
        for my $neighbour ($self->graph->neighbours( $vertex )) {
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
    return path_SMILES( $self->graph, $self->vertices );
}

sub bonds()
{
    my( $self ) = @_;

    return @{$self->{bonds}} if $self->{bonds};

    my $graph = $self->graph;
    my @vertices = $self->vertices;

    my @bonds;
    for (0..$self->length-2) {
        if( $graph->has_edge_attribute( @vertices[$_ .. $_+1], 'bond' ) ) {
            push @bonds, $graph->get_edge_attribute( @vertices[$_ .. $_+1], 'bond' );
        } else {
            push @bonds, '-';
        }
    }

    $self->{bonds} = \@bonds;
    return @bonds;
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
        next if ChemOnomatopist::is_element( $vertices[$_], 'c' ); # FIXME: is_element() should pay attention to aromaticity
        push @heteroatom_positions, $_;
    }

    $self->{heteroatom_positions} = \@heteroatom_positions;
    return @heteroatom_positions;
}

sub is_saturated()
{
    my( $self ) = @_;
    return all { $_ eq '-' } $self->bonds;
}

sub max_valence()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    my $graph = $self->graph;
    my $subgraph = $graph->subgraph( \@vertices );

    my %bond_order = (
        '-' => 1,
        ':' => 1.5, # Not in OpenSMILES
        '=' => 2,
        '#' => 3,
        '$' => 4,
    );

    my $max_valence = 0;
    for my $vertex (@vertices) {
        next if blessed $vertex; # Groups probably have full valence
        next unless exists $normal_valence{ucfirst $vertex->{symbol}};

        my $valence = 0;
        for my $neighbour ($subgraph->neighbours( $vertex )) {
            if( $graph->has_edge_attribute( $vertex, $neighbour, 'bond' ) &&
                exists $bond_order{$graph->get_edge_attribute( $vertex, $neighbour, 'bond' )} ) {
                $valence += $bond_order{$graph->get_edge_attribute( $vertex, $neighbour, 'bond' )};
            } else {
                $valence++;
            }
        }

        my( $normal_valence ) = grep { $_ >= $valence } @{$normal_valence{ucfirst $vertex->{symbol}}};
        $max_valence += $normal_valence - $valence if defined $valence;
    }

    return $max_valence;
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

sub multiple_bond_positions()
{
    my( $self ) = @_;
    my @bonds = $self->bonds;
    return grep { $bonds[$_] =~ /^[=#\$]$/ } 0..$#bonds;
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

    # has_vertex() is used to filter out neighbours within the chain
    my $number = sum0 map { $_ > 2 ? $_ - 2 : 0 }
                          map  { $graph->degree( $_ ) }
                          map  { Graph::Traversal::DFS->new( $graph, start => $_ )->dfs }
                          grep { $graph->has_vertex( $_ ) }
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

sub number_of_double_bonds()
{
    my( $self ) = @_;
    return scalar grep { $_ eq '=' } $self->bonds;
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

sub number_of_multiple_bonds()
{
    my( $self ) = @_;
    return scalar grep { $_ =~ /^[=#\$]$/ } $self->bonds;
}

sub _disconnected_chain_graph()
{
    my( $self ) = @_;

    my $graph = $self->graph->copy; # Maybe use our own graph?
    $graph->delete_path( $self->vertices );

    return $graph;
}

1;
