package ChemOnomatopist::Chain::Monospiro;

# ABSTRACT: Monospiro compound
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

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

    return bless { graph => $graph, spiro_atom => $spiro_atom, components => \@components }, $class;
}

sub candidates()
{
    my( $self ) = @_;

    # "Numbering starts in the smaller ring, if one is smaller, at a ring atom next to the spiro atom and proceeds first around that ring, then through the spiro atom and around the second ring."
    my $graph = $self->graph;
    my $spiro_atom = $self->{spiro_atom};
    my( $A, $B ) = $self->components;
    my @candidates;
    push @candidates,
         $self,
         bless( { graph => $graph, spiro_atom => $spiro_atom, components => [ $A, [ reverse @$B ] ], candidate_for => $self } ),
         bless( { graph => $graph, spiro_atom => $spiro_atom, components => [ [ reverse @$A ], $B ], candidate_for => $self } ),
         bless( { graph => $graph, spiro_atom => $spiro_atom, components => [ [ reverse @$A ], [ reverse @$B ] ], candidate_for => $self } );

    if( @$A == @$B ) {
        push @candidates,
             bless( { graph => $graph, spiro_atom => $spiro_atom, components => [ $B, $A ], candidate_for => $self } ),
             bless( { graph => $graph, spiro_atom => $spiro_atom, components => [ $B, [ reverse @$A ] ], candidate_for => $self } ),
             bless( { graph => $graph, spiro_atom => $spiro_atom, components => [ [ reverse @$B ], $A ], candidate_for => $self } ),
             bless( { graph => $graph, spiro_atom => $spiro_atom, components => [ [ reverse @$B ], [ reverse @$A ] ], candidate_for => $self } );
    }

    return @candidates;
}

sub components()
{
    my( $self ) = @_;
    return @{$self->{components}};
}

sub vertices()
{
    my( $self ) = @_;
    my( $A, $B ) = $self->components;
    my @vertices = ( @$A, $self->{spiro_atom}, @$B );
    return @vertices;
}

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my %degrees = map { $graph->degree( $_ ) => 1 } $graph->vertices;
    return '' unless join( ',', sort keys %degrees ) eq '2,4';

    my @d4 = grep { $graph->degree( $_ ) == 4 } $graph->vertices;
    return '' unless @d4 == 1;
    return 1;
}

sub needs_charge_locants() { 1 }
sub needs_isotope_locants() { 1 }
sub needs_multiple_bond_locants() { 1 }
sub needs_substituent_locants() { 1 }

sub locants() { shift; return map { $_ + 1 } @_ }

sub prefix()
{
    my( $self ) = @_;
    my $name = ChemOnomatopist::Name->new( 'spiro' );
    $name .= ChemOnomatopist::Name::Part::Fusion->new( '[' . join( '.', map { scalar @$_ } $self->components ) . ']' );
    $name .= ChemOnomatopist::alkane_chain_name( $self->length ) . 'an';

    if( $self->parent ) {
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $self->parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }

    $name .= 'yl';
    return $name;
}

sub suffix()
{
    my( $self ) = @_;
    my $name = ChemOnomatopist::Name->new( 'spiro' );
    $name .= ChemOnomatopist::Name::Part::Fusion->new( '[' . join( '.', map { scalar @$_ } $self->components ) . ']' );
    $name .= $self->SUPER::suffix;
    return $name;
}

1;
