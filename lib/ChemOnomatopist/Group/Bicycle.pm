package ChemOnomatopist::Group::Bicycle;

use strict;
use warnings;

# ABSTRACT: Fused bicyclic group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

use Graph::Traversal::DFS;
use Set::Object qw( set );

# From BBv2 P-25.2.1
our @names = (
    [ 'N=CN=CCC', 'NC=CN=CC=', 'pteridine' ],
    [ 'N=NC=CCC', 'CC=CC=CC=', 'cinnoline' ],
    [ 'N=CN=CCC', 'CC=CC=CC=', 'quinazoline' ],
    [ 'N=CC=NCC', 'CC=CC=CC=', 'quinoxaline' ],
    [ 'N=CC=CCC', 'NC=CC=CC=', '1,5-naphthyridine' ], # TODO: There are isomers
    [ 'C=NN=CCC', 'CC=CC=CC=', 'phthalazine' ],
    [ 'N=CC=CCC', 'CC=CC=CC=', 'quinoline' ],
    [ 'C=NC=CCC', 'CC=CC=CC=', 'isoquinoline' ],
    [ 'CC=CCNC',  'C=CC=CCN',  'quinolizine' ],

    [ 'C=NC=NC=', 'N=CNC=C', 'purine' ], # TODO: Special rules apply

    [ 'NN=CCC',  'CC=CC=CC=', 'indazole' ],
    [ 'NC=CCC',  'CC=CC=CC=', 'indole' ],
    [ 'CNC=CC=', 'C=CC=CCC',  'isoindole' ],
    [ 'CC=CNC=', 'C=CC=CCN',  'indolizine', ],
    [ 'CC=CNC',  'C=CC=CN',   '1H-pyrrolizine' ], # TODO: There are isomers
);

sub new
{
    my( $class, $graph, @vertices ) = @_;

    my $subgraph = $graph->subgraph( \@vertices );
    my @bridge = grep { $subgraph->degree( $_ ) == 3 } @vertices;
    $subgraph->delete_vertices( @bridge );

    # Graph is broken into components.
    # Each component is represented as an array of vertices in the order of traverse.
    my @components = sort { @$a <=> @$b } $subgraph->connected_components;
    for (0..1) {
        my $subgraph = $graph->subgraph( [ @{$components[$_]} ] );
        $subgraph->delete_edge( @bridge );
        $components[$_] = [ Graph::Traversal::DFS->new( $subgraph, start => $bridge[$_] )->dfs ];
    }

    # Finding a name from list
    # TODO: Unused
    my @A = @{$components[0]};
    my @B = @{$components[1]};
    my @SMILES;
    push @SMILES, cycle_SMILES( $graph, @A[1..$#A], $A[0] ) . '|' .
                  cycle_SMILES( $graph, @B[1..$#B], $B[0] );
    if( @A == @B ) {
        push @SMILES, cycle_SMILES( $graph, @B[1..$#B], $B[0] ) . '|' .
                      cycle_SMILES( $graph, @A[1..$#A], $A[0] );
    }

    my @cycles =  map { ChemOnomatopist::Group::Monocycle->new( $graph, @$_ ) }
                      @components;

    return bless { graph => $graph, vertices => \@vertices, cycles => \@cycles }, $class;
}

1;
