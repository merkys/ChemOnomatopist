package ChemOnomatopist::Group::Monospiro;

use strict;
use warnings;

# ABSTRACT: Monospiro compound
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use ChemOnomatopist;
use ChemOnomatopist::Chain::Circular;
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
         ChemOnomatopist::Chain::Circular->new( $graph,
                                                @$A, $spiro_atom, @$B ),
         ChemOnomatopist::Chain::Circular->new( $graph,
                                                @$A, $spiro_atom, reverse(@$B) ),
         ChemOnomatopist::Chain::Circular->new( $graph,
                                                reverse(@$A), $spiro_atom, @$B ),
         ChemOnomatopist::Chain::Circular->new( $graph,
                                                reverse(@$A), $spiro_atom, reverse(@$B) );

    if( @$A == @$B ) {
        push @chains,
             ChemOnomatopist::Chain::Circular->new( $graph,
                                                    @$B, $spiro_atom, @$A ),
             ChemOnomatopist::Chain::Circular->new( $graph,
                                                    @$B, $spiro_atom, reverse(@$A) ),
             ChemOnomatopist::Chain::Circular->new( $graph,
                                                    reverse(@$B), $spiro_atom, @$A ),
             ChemOnomatopist::Chain::Circular->new( $graph,
                                                    reverse(@$B), $spiro_atom, reverse(@$A) );
    }

    my( $chain ) = ChemOnomatopist::filter_chains( @chains );

    return bless { graph => $graph,
                   vertices => [ $chain->vertices ],
                   spiro_atom => $spiro_atom,
                   components => \@components },
                 $class;
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

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my %degrees = map { $graph->degree( $_ ) => 1 } $graph->vertices;
    return '' unless join( ',', sort keys %degrees ) eq '2,4';

    my @d4 = grep { $graph->degree( $_ ) == 4 } $graph->vertices;
    return '' unless @d4 == 1;
    return 1;
}

sub prefix($@)
{
    my( $self, $parent ) = @_;
    my $name = ChemOnomatopist::Name->new( 'spiro' );
    $name .= ChemOnomatopist::Name::Part::Fusion->new( '[' . join( '.', map { scalar @$_ } $self->components ) . ']' );
    $name .= ChemOnomatopist::alkane_chain_name( $self->length ) . 'an';

    if( $parent ) {
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }

    $name .= 'yl';
    return $name;
}

# FIXME: This is a bit strange: class and object method with the same name
sub suffix(@)
{
    my( $self ) = @_;
    return '' unless ref $self;
    my $name = ChemOnomatopist::Name->new( 'spiro' );
    $name .= ChemOnomatopist::Name::Part::Fusion->new( '[' . join( '.', map { scalar @$_ } $self->components ) . ']' );
    $name .= ChemOnomatopist::unbranched_chain_name( $self );
    return $name;
}

1;
