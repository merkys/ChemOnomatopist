package ChemOnomatopist::Group::Bicycle;

use strict;
use warnings;

# ABSTRACT: Fused bicyclic group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist;
use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Group::Monocycle::Fused;
use ChemOnomatopist::Name;
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );
use Graph::Traversal::DFS;
use List::Util qw( all any );
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

    [ 'NC=NC=CC=', 'NC=NCC', 'purine' ], # Special rules apply

    [ 'NN=CCC',  'CC=CC=CC=', '1H-indazole' ],
    [ 'NC=CC=C', 'C=CC=CC=C', '1H-indole' ],
    [ 'CNC=CC=', 'C=CC=CCC',  'isoindole' ],
    [ 'CC=CNC=', 'C=CC=CCN',  'indolizine', ],
    [ 'CC=CNC',  'C=CC=CN',   '1H-pyrrolizine' ], # TODO: There are isomers
);

# From BBv2 P-25.1.1, order of decreasing seniority
our %hydrocarbons_by_size = (
    '5,7' => 'azulene',
    '6,6' => 'naphthalene',
    '5,6' => 'indene',
);

sub new
{
    my( $class, $graph, @vertices ) = @_;

    # Has to be created early to be given for fused parts
    my $self = bless { graph => $graph, vertices => \@vertices }, $class;

    my $subgraph = $graph->subgraph( \@vertices );
    my @bridge = grep { $subgraph->degree( $_ ) == 3 } @vertices;
    $subgraph->delete_edge( @bridge );
    $self->{vertices} = [ Graph::Traversal::DFS->new( $subgraph, start => $bridge[0] )->dfs ];
    $subgraph->delete_vertices( @bridge );

    $self->_aromatise;

    # Graph is broken into components.
    # Each component is represented as an array of vertices in the order of traverse.
    my @components = sort { @$a <=> @$b } $subgraph->connected_components;
    for (0..1) {
        my $subgraph = $graph->subgraph( [ @{$components[$_]}, @bridge ] );
        $subgraph->delete_edge( @bridge );
        my @path = Graph::Traversal::DFS->new( $subgraph, start => $bridge[$_] )->dfs;
        push @path, shift @path;
        $components[$_] = \@path;
    }
    my @cycles = map { ChemOnomatopist::Group::Monocycle::Fused->new( $graph, $self, @$_ ) }
                     @components;
    $self->{cycles} = \@cycles;

    my $nbenzene = scalar grep { $_->is_benzene } @cycles;

    # The ordering should not be done if one of the cycles is benzene
    if( $nbenzene == 0 ) {
        my @flipped = map { $_->flipped } @cycles;
        # CHECKME: Additional rules from ChemOnomatopist::filter_chains() might still be needed
        my( $chain ) = sort { ChemOnomatopist::Chain::Circular::_cmp( $a, $b ) } ( @cycles, @flipped );

        if(      $chain == $cycles[1] ) {
            @cycles = reverse @cycles;
        } elsif( $chain == $flipped[0] ) {
            @cycles = @flipped;
        } elsif( $chain == $flipped[1] ) {
            @cycles = reverse @flipped;
        }
        $self->{cycles} = \@cycles;
        $self->_adjust_vertices_to_cycles;
    } elsif( $nbenzene == 1 ) {
        # Numbering has to start from cycle other than benzene
        if( $cycles[0]->is_benzene ) {
            @cycles = reverse @cycles;
            $self->{cycles} = \@cycles;
        }

        my( $chain ) = sort { ChemOnomatopist::Chain::Circular::_cmp( $a, $b ) }
                            ( $cycles[0], $cycles[0]->flipped );

        if( $chain != $cycles[0] ) {
            @cycles = map { $_->flipped } @cycles;
            $self->{cycles} = \@cycles;
        }
        $self->_adjust_vertices_to_cycles;
    }

    if( join( ',', map { $_->backbone_SMILES } @cycles ) eq 'N=CNCC,CN=CN=CC=' ) {
        @cycles = reverse map { $_->flipped } @cycles;
        $self->{cycles} = \@cycles;
    }

    return $self;
}

sub candidates()
{
    my( $self ) = @_;

    if( $self->suffix eq 'naphthalene' ) {
        # Generates all variants
        my @chains = ( $self, $self->copy, $self->copy, $self->copy );

        $chains[1]->{cycles} = [ map { $_->flipped } $chains[1]->cycles ];
        $chains[3]->{cycles} = [ map { $_->flipped } $chains[3]->cycles ];

        $chains[2]->{cycles} = [ reverse $chains[2]->cycles ];
        $chains[3]->{cycles} = [ reverse $chains[3]->cycles ];

        for (@chains) {
            $_->_adjust_vertices_to_cycles;
        }

        return @chains;
    }

    return $self;
}

sub copy()
{
    my( $self ) = @_;
    return bless { graph    => $self->graph,
                   cycles   => [ $self->cycles ],
                   vertices => [ $self->vertices ] },
                 ChemOnomatopist::Group::Bicycle::;
}

sub cycles()
{
    my( $self ) = @_;
    return @{$self->{cycles}};
}

# Tells whether the outer bonds of the bicycle qualify as aromatic
sub is_aromatic()
{
    my( $self ) = @_;
    my @outer_vertices;
    for ($self->cycles) {
        my @vertices = $_->vertices;
        pop @vertices;
        push @outer_vertices, @vertices;
    }
    return ChemOnomatopist::Chain::Circular->new( $self->graph, @outer_vertices )->is_aromatic;
}

sub is_hydrocarbon()
{
    my( $self ) = @_;
    return all { $_->is_hydrocarbon } $self->cycles;
}

# Implemented according to BBv2 P-25.3.3
sub locants(@)
{
    my $self = shift;
    my @vertices = $self->vertices;
    my @cycles = $self->cycles;
    my @locant_map = map { $_ + 1 } 0..$self->length - 1;

    if( ChemOnomatopist::is_element( $cycles[0]->{vertices}[-2], 'C' ) ) {
        splice @locant_map, $cycles[0]->length-2, 0, ($cycles[0]->length - 2) . 'a';
        pop @locant_map;
    }

    if( ChemOnomatopist::is_element( $cycles[1]->{vertices}[-2], 'C' ) ) {
        $locant_map[-1] = $locant_map[-2] . 'a';
    }

    return map { $locant_map[$_] } @_;
}

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return $self->suffix =~ /^benzo/;
}

sub needs_heteroatom_names() { return '' } # FIXME: This is not always correct

# FIXME: This is a bit strange: class and object method with the same name
sub suffix()
{
    my( $self ) = @_;
    return '' unless ref $self;
    if( $self->is_hydrocarbon ) {
        # FIXME: Check if aromatic, but with caution, as substitutions will break aromaticity
        my $cycle_sizes = join ',', map { $_->length } $self->cycles;
        return $hydrocarbons_by_size{$cycle_sizes} if exists $hydrocarbons_by_size{$cycle_sizes};

        if( $cycle_sizes =~ /^(\d+),\1$/ ) {
            return ChemOnomatopist::alkane_chain_name( $1 ) . 'alene';
        }
    }

    my @SMILES = map { $_->backbone_SMILES } $self->cycles;
    my( $retained ) = grep { ($_->[0] eq $SMILES[0] && $_->[1] eq $SMILES[1]) ||
                             ($_->[0] eq $SMILES[1] && $_->[1] eq $SMILES[0]) } @names;
    return ChemOnomatopist::Name->new( $retained->[2] ) if $retained;

    if( any { $_->is_benzene } $self->cycles ) {
        my( $other ) = grep { !$_->is_benzene } $self->cycles;
        $other = ChemOnomatopist::Group::Monocycle->new( $other->graph, $other->vertices );
        my $suffix = $other->suffix;
        $suffix =~ s/^o//;
        return 'benzo' . $suffix;
    }

    die "cannot name complex bicyclic compounds\n";
}

sub _adjust_vertices_to_cycles()
{
    my( $self ) = @_;

    my @cycles = $self->cycles;

    $self->{vertices} = [];
    push @{$self->{vertices}}, $cycles[0]->vertices;
    pop  @{$self->{vertices}};
    push @{$self->{vertices}}, $cycles[1]->vertices;
    pop  @{$self->{vertices}};

    return $self;
}

1;
