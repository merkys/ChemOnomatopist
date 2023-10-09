package ChemOnomatopist::Chain::Bicycle;

use strict;
use warnings;

# ABSTRACT: Fused bicyclic chain
# VERSION

use ChemOnomatopist;
use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Chain::Bicycle::Purine;
use ChemOnomatopist::Chain::Monocycle;
use ChemOnomatopist::Chain::Monocycle::Fused;
use ChemOnomatopist::Name;
use ChemOnomatopist::Name::Part::Fusion;
use ChemOnomatopist::Name::Part::Stem;
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );
use Chemistry::OpenSMILES qw( is_double_bond );
use Graph::Traversal::DFS;
use List::Util qw( all any min uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

use parent ChemOnomatopist::Chain::Circular::;

# From BBv2 P-25.2.1
our @names = (
    [ 'n:c:n:c:c:c:', 'n:c:c:n:c:c:', 'pteridine' ],
    [ 'c:n:c:n:c:c:', 'n:c:c:n:c:c:', 'pteridine' ],

    [ 'N=NC=CCC', 'c:c:c:c:c:c:', 'cinnoline' ],
    [ 'N=CN=CCC', 'c:c:c:c:c:c:', 'quinazoline' ],

    [ 'N=CC=NCC',     'c:c:c:c:c:c:', 'quinoxaline' ],
    [ 'n:c:c:n:c:c:', 'c:c:c:c:c:c:', 'quinoxaline' ],

    [ 'N=CC=CCC', 'NC=CC=CC=', '1,5-naphthyridine' ], # TODO: There are isomers
    [ 'C=NN=CCC', 'c:c:c:c:c:c:', 'phthalazine' ],
    [ 'n:c:c:c:c:c:', 'c:c:c:c:c:c:', 'quinoline' ],

    [ 'C=NC=CCC',     'c:c:c:c:c:c:', 'isoquinoline' ],
    [ 'c:n:c:c:c:c:', 'c:c:c:c:c:c:', 'isoquinoline' ],

    [ 'CC=CCNC',  'C=CC=CCN',  'quinolizine' ],

    [ 'c:n:c:n:c:c:', 'n:c:n:c:c:', 'purine' ], # Special rules apply

    [ 'n:n:c:c:c:', 'c:c:c:c:c:c:', '1H-indazole' ],
    [ 'n:c:c:c:c:', 'c:c:c:c:c:c:', '1H-indole' ],
    [ 'c:n:c:c:c:', 'c:c:c:c:c:c:', 'isoindole' ],
    [ 'c:c:c:c:n:', 'c:c:c:c:n:c:', 'indolizine', ],
    [ 'CC=Cn:c',    'c:c:c:c:n:',   '1H-pyrrolizine' ], # TODO: There are isomers
);

for my $name (qw( 1H-indole indolizine isoindole isoquinoline quinoline quinolizine )) {
    for (grep { $_->[2] eq $name } @names) {
        my @As_parts = @$_;
        $As_parts[0] =~ s/N/\[As\]/g;
        $As_parts[1] =~ s/N/\[As\]/g;
        $As_parts[2] =~ s/^1H-//;
        $As_parts[2] = 'ars' . $As_parts[2] unless $As_parts[2] =~ s/^iso/isoars/;
        push @names, \@As_parts;

        my @P_parts = @$_;
        $P_parts[0] =~ s/n/p/g;
        $P_parts[1] =~ s/n/p/g;
        $P_parts[2] =~ s/^1H-//;
        $P_parts[2] = 'phosph' . $P_parts[2] unless $P_parts[2] =~ s/^iso/isophosph/;
        push @names, \@P_parts;
    }
}

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
    my @cycles = map { ChemOnomatopist::Chain::Monocycle::Fused->new( $graph, $self, @$_ ) }
                     @components;
    $self->{cycles} = \@cycles;

    my $nbenzene = scalar grep { $_->is_benzene } @cycles;

    # The ordering should not be done if one of the cycles is benzene
    if( $nbenzene == 0 ) {
        # The following code is supposed to order the rings _and_ establish the traversal order
        # TODO: Maybe all traversals of both rings should be generated here?
        # FIXME: There seem to be two orders for rings: one for numbering, other for fusion naming...
        my @candidates = map { ChemOnomatopist::Chain::Monocycle->new( $_->graph, $_->vertices ) }
                             ( @cycles, map { $_->flipped } @cycles );
        for my $rule ( # P-25.3.2.4 (a): Senior heteroatom according to specific seniority order
                       \&rule_most_senior_heteroatom,
                       # TODO: P-25.3.2.4 (b): Concerns fusions of more than two rings
                       # P-25.3.2.4 (c): Second ring has to be larger
                       \&ChemOnomatopist::rule_longest_chains,
                       # P-25.3.2.4 (d): Greater number of heteroatoms of any kind
                       \&ChemOnomatopist::rule_most_heteroatoms,
                       # P-25.3.2.4 (e): Greater variety of heteroatoms
                       \&rule_greatest_variety_of_heteroatoms,
                       # P-25.3.2.4 (f): Greater number of most senior heteroatoms
                       \&ChemOnomatopist::rule_greatest_number_of_most_senior_heteroatoms,
                       # TODO: P-25.3.2.4 (g): Concerns fusions of more than two rings
                       # P-25.3.2.4 (h): Lower locants for heteroatoms
                       \&ChemOnomatopist::rule_lowest_numbered_heteroatoms,
                       # P-25.3.2.4 (i): Lower locants for senior heteroatoms
                       \&ChemOnomatopist::rule_lowest_numbered_most_senior_heteroatoms,
                       # TODO: P-25.3.2.4 (j): Concerns fusions of more than two rings
                     ) {
            my @candidates_now = $rule->( @candidates );
            if( @candidates_now == 1 ) {
                @candidates = @candidates_now;
                last;
            } elsif( @candidates ) {
                @candidates = @candidates_now;
            } else {
                last;
            }
        }
        # Here the "winning" ring is selected
        my $chain = shift @candidates;

        # Making the "winning" ring the first
        if( set( $chain->vertices ) == set( $cycles[1]->vertices ) ) {
            @cycles = reverse @cycles;
        }

        # Now we have to match the direction, I guess
        my  @bridge_in_chain = grep { set( @bridge )->has( $_ ) } $chain->vertices;
        if( $bridge_in_chain[-1] != $cycles[0]->{vertices}[-1] ) {
            @cycles = map { $_->flipped } @cycles;
        }

        $self->{cycles} = \@cycles;
        $self->_adjust_vertices_to_cycles;

        if(      join( ',', map { $_->backbone_SMILES } @cycles ) =~ /^c:n:c:n:c:c:,n:c:n:c:c:$/ ||
                 join( ',', map { $_->backbone_SMILES } @cycles ) =~ /^CNC=Nc:c,n:c:n:c:c:$/ ) {
            return ChemOnomatopist::Chain::Bicycle::Purine->new( $graph, @cycles );
        } elsif( join( ',', map { $_->backbone_SMILES } @cycles ) =~ /^N=CNCc:c,n:c:n:c:c:$/ ) {
            return ChemOnomatopist::Chain::Bicycle::Purine->new( $graph, map { $_->flipped } @cycles );
        }
    } elsif( $nbenzene == 1 ) {
        # Numbering has to start from cycle other than benzene
        if( $cycles[0]->is_benzene ) {
            @cycles = reverse @cycles;
            $self->{cycles} = \@cycles;
        }

        my( $chain ) = sort { ChemOnomatopist::Chain::Monocycle::_cmp( $a, $b ) }
                            ( $cycles[0], $cycles[0]->flipped );

        if( $chain != $cycles[0] ) {
            @cycles = map { $_->flipped } @cycles;
            $self->{cycles} = \@cycles;
        }
        $self->_adjust_vertices_to_cycles;
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
            $_->{candidate_for} = $self unless $_ == $self;
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
                 ChemOnomatopist::Chain::Bicycle::;
}

sub cycles()
{
    my( $self ) = @_;
    return @{$self->{cycles}};
}

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my %degrees = map { $graph->degree( $_ ) => 1 } $graph->vertices;
    return '' unless join( ',', sort keys %degrees ) eq '2,3';

    my @d3 = grep { $graph->degree( $_ ) == 3 } $graph->vertices;
    return '' unless @d3 == 2;
    return '' unless $graph->has_edge( @d3 );

    $graph = $graph->copy->delete_vertices( @d3 );
    return '' unless scalar( $graph->connected_components ) == 2;
    return 1;
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

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return $self->suffix =~ /^benzo/;
}

sub needs_heteroatom_names() { return '' } # FIXME: This is not always correct

sub prefix()
{
    my( $self ) = @_;

    my $name = $self->suffix;
    $name = ChemOnomatopist::Name->new( $name ) unless blessed $name;
    $name->{name}[-1] =~ s/e$//;
    if( $self->parent ) { # FIXME: Not stable for naphthalene
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $self->parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }
    $name .= 'yl';

    return $name;
}

# FIXME: This is a bit strange: class and object method with the same name
sub suffix()
{
    my( $self ) = @_;
    return '' unless ref $self;

    if( $self->is_hydrocarbon ) {
        # FIXME: Check if aromatic, but with caution, as substitutions will break aromaticity
        my $cycle_sizes = join ',', map { $_->length } $self->cycles;
        if( exists $hydrocarbons_by_size{$cycle_sizes} ) {
            return ChemOnomatopist::Name::Part::Stem->new( $hydrocarbons_by_size{$cycle_sizes} )->to_name;
        }

        if( $cycle_sizes =~ /^(\d+),\1$/ ) {
            my $name = ChemOnomatopist::alkane_chain_name( $1 ) . 'alene';
            return ChemOnomatopist::Name::Part::Stem->new( $name )->to_name;
        }
    }

    my @SMILES = map { $_->backbone_SMILES } $self->cycles;
    print STDERR "bicycle SMILES: @SMILES\n" if $ChemOnomatopist::DEBUG;
    my( $retained ) = grep { ($_->[0] eq $SMILES[0] && $_->[1] eq $SMILES[1]) ||
                             ($_->[0] eq $SMILES[1] && $_->[1] eq $SMILES[0]) } @names;
    return ChemOnomatopist::Name::Part::Stem->new( $retained->[2] )->to_name if $retained;

    if( any { $_->is_benzene } $self->cycles ) {
        my( $other ) = grep { !$_->is_benzene } $self->cycles;
        $other = ChemOnomatopist::Chain::Monocycle->new( $other->graph, $other->vertices );

        my $SMILES = $other->backbone_SMILES;
        if( $SMILES =~ /^C=C((?<el>O|S|\[Se\]|\[Te\])C|C(?<el>O|S|\[Se\]|\[Te\]))c:c$/ ) {
            # Names according to BBv2 P-25.2.1, Table 2.8, (23) and (24)
            my $element = $+{el};
            my $name = ($1 =~ /^C/ ? '2H-1-' : '1H-2-') . 'benzo';
            $element =~ s/[\[\]]//g;
            if( $element ne 'O' ) {
                $name .= $elements{$element}->{prefix};
                $name =~ s/a$/o/;
            }
            return $name . 'pyran';
        } else {
            my $name = ChemOnomatopist::Name->new( 'benzo' );
            my  $other_name = $other->suffix;
            if( $other_name->starts_with_locant ) { # Locants are moved to front
                unshift @$name, shift @$other_name;
            }
            $name .= $other_name;
            return $name;
        }
    }

    # Fusion naming according to BBv2 P-25.3.1.3

    # Find the bridge vertices
    my @bridge = $self->{cycles}[0]->vertices;
    @bridge = @bridge[-2..-1];

    # Find autosymmetric equivalents having the least locants for the bridge
    my @equiv_A = $self->{cycles}[0]->autosymmetric_equivalents;
    my @equiv_B = $self->{cycles}[1]->autosymmetric_equivalents;

    my $min_A = min map {  $_->vertex_ids( @bridge ) } @equiv_A;
    @equiv_A = grep { min( $_->vertex_ids( @bridge ) ) == $min_A } @equiv_A;

    my $min_B = min map {  $_->vertex_ids( @bridge ) } @equiv_B;
    @equiv_B = grep { min( $_->vertex_ids( @bridge ) ) == $min_B } @equiv_B;

    my $fusion = '[';
    if( @equiv_A > 1 || @equiv_B > 1 ) {
        # At least one of the rings has mirror symmetry ("flip-symmetric"), thus numeric order is ascending
        $fusion .= ($min_B+1) . ',' . ($min_B+2);
    } else {
        # Rings are rigid, thus numeric order has to be derived
        my @order_A = $equiv_A[0]->vertex_ids( @bridge );
        my @order_B = $equiv_B[0]->vertex_ids( @bridge );
        if( ($order_A[0] <=> $order_A[1]) == ($order_B[0] <=> $order_B[1]) ) {
            # Ring atoms are encountered in the same order in both of the rings
            $fusion .= ($min_B+1) . ',' . ($min_B+2);
        } else {
            # Ring atom orders differ
            $fusion .= ($min_B+2) . ',' . ($min_B+1);
        }
    }
    $fusion .= '-' . chr( 97 + $min_A ) . ']';

    my $name = ChemOnomatopist::Name->new;
    my $graph = $self->graph;

    $name->append_locants( map { $_ . 'H' } $self->locants( $self->indicated_hydrogens ) );

    my @ideal = map { ChemOnomatopist::Chain::Monocycle->new( $_->graph, $_->vertices ) }
                    $self->cycles;
    my $name_A = $ideal[1]->name;
    $name_A =~ s/^\d+H-//;

    # TODO: Complete retained prefixes from BBv2 P-25.3.2.2.3
    $name_A = 'fur'     if $name_A eq 'furan';
    $name_A = 'imidaz'  if $name_A eq 'imidazole';
    $name_A = 'pyrid'   if $name_A eq 'pyridine';
    $name_A = 'pyrimid' if $name_A eq 'pyrimidine';
    $name_A = 'thien'   if $name_A eq 'thiophene';

    $name .= $name_A;
    unless( $name->[-1] =~ s/e$/o/ ) { # BBv2 P-25.3.2.2.2
        $name->[-1] .= 'o';
    }
    $name .= ChemOnomatopist::Name::Part::Fusion->new( $fusion );

    my $name_B = $ideal[0]->name;
    $name_B->[0] =~ s/\d+H-//;
    $name .= $name_B;
    $name->bracket_numeric_locants;
    return $name;
}

sub rule_most_senior_heteroatom
{
    my( @chains ) = @_;

    # This order is taken from BBv2 P-25.3.2.4 (a) and is different from order in %elements
    my @element_order = qw( N F Cl Br I O S Se Te P As Sb Bi Si Ge Sn Pb B Al Ga In Tl );
    my %element_order = map { $element_order[$_] => $_ } 0..$#element_order;

    my( $max_value ) = sort { $element_order{$a} <=> $element_order{$b} }
                       grep { exists $element_order{$_} }
                       map  { $_->heteroatoms } @chains;
    return @chains unless $max_value;
    return grep { any { $_ eq $max_value } $_->heteroatoms } @chains;
}

sub rule_greatest_variety_of_heteroatoms
{
    my( @chains ) = @_;

    my( $max_value ) = reverse sort map { scalar uniq $_->heteroatoms } @chains;
    return @chains unless $max_value;
    return grep { scalar( uniq $_->heteroatoms ) == $max_value } @chains;
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
