package ChemOnomatopist::Chain::Bicycle;

# ABSTRACT: Fused bicyclic chain
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Chain::Bicycle::Purine;
use ChemOnomatopist::Chain::Monocycle;
use ChemOnomatopist::Chain::Monocycle::Fused;
use ChemOnomatopist::Comparable::Array::Numeric;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Name;
use ChemOnomatopist::Name::Part::Fusion;
use ChemOnomatopist::Name::Part::Stem;
use ChemOnomatopist::Util qw( all_max all_min );
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );
use Chemistry::OpenSMILES qw( is_double_bond );
use Graph::Traversal::DFS;
use List::Util qw( all any first min uniq );
use Scalar::Util qw( blessed );
use Set::Object qw( set );

use parent ChemOnomatopist::Chain::Circular::;

# From BBv2 P-25.2.1
our @names = (
    [ qw( NCNCCC NCCNCC pteridine ) ],
    [ qw( CNCNCC NCCNCC pteridine ) ],

    [ qw( NNCCCC CCCCCC cinnoline ) ],
    [ qw( NCNCCC CCCCCC quinazoline ) ],

    [ qw( NCCNCC CCCCCC quinoxaline ) ],

    [ 'c:n:n:c:c:c:', 'c:c:c:c:c:c:', 'phthalazine' ],
    [ 'n:c:c:c:c:c:', 'c:c:c:c:c:c:', 'quinoline' ],

    [ 'C=NC=CCC',     'c:c:c:c:c:c:', 'isoquinoline' ],
    [ 'c:n:c:c:c:c:', 'c:c:c:c:c:c:', 'isoquinoline' ],

    [ 'CC=CCn:c=', 'c:c:c:c:c:n:', 'quinolizine' ],
    [ 'CC=CC=c:n', 'c:c:c:c:n:c:', 'quinolizine' ],

    [ 'n:n:c:c:c:', 'c:c:c:c:c:c:', 'indazole' ],
    [ 'n:c:c:c:c:', 'c:c:c:c:c:c:', 'indole' ],
    [ 'c:n:c:c:c:', 'c:c:c:c:c:c:', 'isoindole' ],
    [ 'c:c:c:c:c:n:', 'c:c:c:n:c:', 'indolizine', ],

    [ qw( CCCCN CCCNC pyrrolizine ) ],
    [ qw( CCCNC CCCNC pyrrolizine ) ],
);

for my $name (qw( indole indolizine isoindole isoquinoline quinoline quinolizine )) {
    for (grep { $_->[2] eq $name } @names) {
        my @As_parts = @$_;
        $As_parts[0] =~ s/n/\[as\]/g;
        $As_parts[1] =~ s/n/\[as\]/g;
        $As_parts[2] =~ s/^qu//;
        $As_parts[2] = 'ars' . $As_parts[2] unless $As_parts[2] =~ s/^iso(qu)?/isoars/;
        push @names, \@As_parts;

        my @P_parts = @$_;
        $P_parts[0] =~ s/n/p/g;
        $P_parts[1] =~ s/n/p/g;
        $P_parts[2] =~ s/^qu//;
        $P_parts[2] = 'phosph' . $P_parts[2] unless $P_parts[2] =~ s/^iso(qu)?/isophosph/;
        push @names, \@P_parts;
    }
}

for (@names) {
    $_->[0] =~ s/[=:]//g;
    $_->[1] =~ s/[=:]//g;
    ( $_->[0], $_->[1] ) = map { uc } @$_[0..1];
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

    my $subgraph = $graph->subgraph( @vertices );
    my @bridge = grep { $subgraph->degree( $_ ) == 3 } @vertices;
    $subgraph->delete_edge( @bridge );
    @vertices = Graph::Traversal::DFS->new( $subgraph, start => $bridge[0] )->dfs;
    $subgraph->delete_vertices( @bridge );

    # Graph is broken into components.
    # Each component is represented as an array of vertices in the order of traverse.
    my @components = sort { @$a <=> @$b } $subgraph->connected_components;
    for (0..1) {
        if( @{$components[$_]} < 3 ) {
            die "bicycles with three and four-membered cycles are not supported yet\n"
        }

        my $subgraph = $graph->subgraph( @{$components[$_]}, @bridge );
        $subgraph->delete_edge( @bridge );
        my @path = Graph::Traversal::DFS->new( $subgraph, start => $bridge[$_] )->dfs;
        push @path, shift @path;
        $components[$_] = \@path;
    }
    my @cycles = map { ChemOnomatopist::Chain::Monocycle::Fused->new( $graph, $self, @$_ ) }
                     @components;
    $self->{cycles} = \@cycles;

    my $nbenzene = scalar grep { $_->is_benzene } @cycles;

    if( !$nbenzene && $self->is_purine ) {
        return ChemOnomatopist::Chain::Bicycle::Purine->new( $graph, @cycles );
    } elsif( !$nbenzene ) {
        # Find the senior cycle
        my @candidates = map { $_->candidates }
                         map { ChemOnomatopist::Chain::Monocycle->new( $_->graph, $_->vertices ) }
                             $self->cycles;
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
            } elsif( @candidates ) { # CHECKME: This looks strange
                @candidates = @candidates_now;
            } else {
                last;
            }
        }

        # Here the "winning" ring is selected
        my $chain = shift @candidates;

        # Making the "winning" ring the first
        if( set( $chain->vertices ) == set( $self->{cycles}[1]->vertices ) ) {
            $self->{cycles} = [ reverse $self->cycles ];
        }

        # Construct the candidates to determine the numbering order
        @candidates = ( $self,
                        $self->flipped_horizontally,
                        $self->flipped_vertically,
                        $self->flipped_horizontally->flipped_vertically );

        # Establish the order
        # FIXME: Somehow the rules from BBv3 P-14.4 interplay here
        for my $rule (
                        # P-25.3.3.1.2 (a): Lower locants for heteroatoms
                        \&ChemOnomatopist::rule_lowest_numbered_heteroatoms,
                        # P-25.3.3.1.2 (b): Lower locants for senior heteroatoms
                        \&ChemOnomatopist::rule_lowest_numbered_most_senior_heteroatoms,
                        # P-25.3.3.1.2 (c): Lower locants for fusion carbon atoms
                        \&rule_lowest_numbered_fusion_carbons,
                        # P-25.3.3.1.2 (d): Lower locants for fusion heteroatoms (rather than nonfusion)
                        \&rule_lowest_numbered_fusion_heteroatoms,
                        # TODO: P-25.3.3.1.2 (e): Lower locants for interior heteroatom
                        # P-25.3.3.1.2 (f): Lower locants for indicated hydrogen atoms
                        \&ChemOnomatopist::rule_lowest_numbered_indicated_hydrogens,

                        # P-14.4 (f): Lower locants for detachable alphabetized prefixes
                        # CHECKME: Not sure what is the relation with other rules
                        # FIXME: Does not take alphabetic order in consideration
                        \&ChemOnomatopist::rule_lowest_numbered_locants,
                     ) {
            my @candidates_now = $rule->( @candidates );
            if( @candidates_now == 1 ) {
                @candidates = @candidates_now;
                last;
            } elsif( @candidates ) { # CHECKME: This looks strange
                @candidates = @candidates_now;
            } else {
                last;
            }
        }

        $self->vertices( (shift @candidates)->vertices ); # FIXME: Simply return instead of self
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
    my @chains = ( $self );

    if( $self->is_naphthalene ) {
        # Generates all variants
        push @chains, $self->flipped_horizontally,
                      $self->flipped_vertically,
                      $self->flipped_horizontally->flipped_vertically;
        for (1..3) {
            $chains[$_]->{candidate_for} = $self;
        }
    }

    return @chains;
}

sub flipped_horizontally()
{
    my( $self ) = @_;
    my $copy = $self->copy;
    my @vertices = reverse $copy->vertices;
    push @vertices, shift @vertices;
    $copy->vertices( @vertices );
    return $copy;
}

sub flipped_vertically()
{
    my( $self ) = @_;
    my $copy = $self->copy;
    my $cycle = first { set( $_->vertices )->has( $copy->{vertices}[0] ) }
                      $copy->cycles;
    my @vertices = reverse $copy->vertices;
    for (1..$cycle->length-2) {
        unshift @vertices, pop @vertices;
    }
    $copy->vertices( @vertices );
    return $copy;
}

sub copy()
{
    my( $self ) = @_;
    return bless { graph    => $self->graph,
                   cycles   => [ $self->cycles ],
                   vertices => [ $self->vertices ],
                   parent   => $self->parent },
                 ChemOnomatopist::Chain::Bicycle::;
}

sub cycles()
{
    my( $self ) = @_;
    return @{$self->{cycles}};
}

sub fusion_positions()
{
    my( $self ) = @_;
    my $bridge = set( map { @{$_->{vertices}}[-2..-1] } $self->cycles );
    my @vertices = $self->vertices;
    return grep { $bridge->has( $vertices[$_] ) } 0..$#vertices;
}

sub fusion_carbon_positions()
{
    my( $self ) = @_;
    return grep { ChemOnomatopist::Util::element( $self->{vertices}[$_] ) eq 'C' }
                $self->fusion_positions;
}

sub fusion_heteroatom_positions()
{
    my( $self ) = @_;
    return grep { ChemOnomatopist::Util::element( $self->{vertices}[$_] ) ne 'C' }
                $self->fusion_positions;
}

sub parent(;$)
{
    my( $self, $parent ) = @_;
    my $old_parent = $self->SUPER::parent( $parent );
    return $old_parent unless $parent;
    return $old_parent if $old_parent && $parent == $old_parent;

    if( $self->is_naphthalene ) {
        my( $chain ) = ChemOnomatopist::filter_chains( $self->candidates );
        $self->vertices( $chain->vertices );
    }

    return $old_parent;
}

sub has_form($$)
{
    my( $class, $graph ) = @_;
    my %degrees = map { $graph->degree( $_ ) => 1 } $graph->vertices;
    return '' unless join( ',', sort keys %degrees ) eq '2,3';

    # Ensure subgraph has three paths between the three-degreed vertices
    my @d3 = grep { $graph->degree( $_ ) == 3 } $graph->vertices;
    return '' unless @d3 == 2;

    $graph = $graph->copy;
    my @lengths;
    for (1..3) {
        my @path = $graph->SP_Dijkstra( @d3 );
        return '' unless @path;

        $graph->delete_path( @path );
        push @lengths, scalar @path;
    }

    # Ensure a single path directly joins the three-degreed vertices
    return '' unless @lengths == 3;
    return '' unless scalar( grep { $_ > 2 } @lengths ) == 2;

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

sub is_naphthalene()
{
    my( $self ) = @_;
    return $self->is_hydrocarbon && all { $_->length == 6 } $self->cycles;
}

sub is_naphthyridine()
{
    my( $self ) = @_;
    return '' unless all { $_->length == 6 && $_->is_monoreplaced } $self->cycles;
    return join( ',', uniq $self->heteroatoms ) eq 'N';
}

sub is_purine()
{
    my( $self ) = @_;

    return '' unless $self->number_of_heteroatoms == 4;
    return '' unless join( ',', uniq $self->heteroatoms ) eq 'N';

    my @cycles = $self->cycles;
    my $pyrimidine = first { $_->length == 6 } @cycles;
    my $imidazole  = first { $_->length == 5 } @cycles;
    return '' unless $pyrimidine && $imidazole;

    return '' unless join( ',', $imidazole->heteroatom_positions ) eq '0,2';
    return '' unless join( ',', $pyrimidine->heteroatom_positions ) eq '0,2' ||
                     join( ',', $pyrimidine->heteroatom_positions ) eq '1,3';

    return 1;
}

sub needs_indicated_hydrogens() { 1 }

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return $self->suffix =~ /^benzo/ || $self->is_naphthyridine;
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;
    return $self->needs_heteroatom_locants &&
           all { !$_->is_Hantzsch_Widman } $self->cycles;
}

sub needs_substituent_locants() { 1 }

sub indicated_hydrogens_part()
{
    my( $self ) = @_;
    my $part = ChemOnomatopist::Name->new;
    return $part unless $self->number_of_indicated_hydrogens;

    return $self->SUPER::indicated_hydrogens_part if $self->is_naphthalene;

    # Use parent procedure if all the indicated hydrogens are confined to a single ring
    return $self->SUPER::indicated_hydrogens_part if any { !$_->number_of_indicated_hydrogens } $self->cycles;

    if( $self->number_of_indicated_hydrogens &&
        $self->number_of_indicated_hydrogens < $self->length ) {
        $part->append_locants( map { $_ . 'H' } $self->locants( $self->indicated_hydrogens ) );
    }

    return $part;
}

sub prefix()
{
    my( $self ) = @_;

    my $name = $self->suffix;
    $name = ChemOnomatopist::Name->new( $name ) unless blessed $name;
    $name->pop_e;
    if( $self->parent ) { # FIXME: Not stable for naphthalene
        my @vertices = $self->vertices;
        my $position = first { $self->graph->has_edge( $self->parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }
    $name .= 'yl';

    return $name;
}

sub suffix()
{
    my( $self ) = @_;

    if( $self->is_hydrocarbon ) {
        # FIXME: Check if aromatic, but with caution, as substitutions will break aromaticity
        my $cycle_sizes = join ',', sort map { $_->length } $self->cycles;
        if( exists $hydrocarbons_by_size{$cycle_sizes} ) {
            return ChemOnomatopist::Name::Part::Stem->new( $hydrocarbons_by_size{$cycle_sizes} )->to_name;
        }

        if( $cycle_sizes =~ /^(\d+),\1$/ ) {
            my $name = ChemOnomatopist::alkane_chain_name( $1 ) . 'alene';
            return ChemOnomatopist::Name::Part::Stem->new( $name )->to_name;
        }
    }

    if( $self->is_naphthyridine ) {
        return ChemOnomatopist::Name::Part::Stem->new( 'naphthyridine' )->to_name;
    }

    my @SMILES = map { $_->backbone_SMILES } $self->cycles;
    print STDERR "bicycle SMILES: @SMILES\n" if $ChemOnomatopist::DEBUG;
    my $retained = find_retained( @SMILES );
    return ChemOnomatopist::Name::Part::Stem->new( $retained->[2] )->to_name if $retained;

    if( any { $_->is_benzene } $self->cycles ) {
        my $other = first { !$_->is_benzene } $self->cycles;

        my $name = ChemOnomatopist::Name->new( 'benzo' );
        if( $other->length == 6 && $other->is_monoreplaced &&
            join( '', $other->heteroatoms ) =~ /^(O|S|Se|Te)$/ &&
            join( '', $other->heteroatom_positions ) < 4 ) {
            # Names according to BBv3 P-25.2.1, Table 2.8, (23) and (24)
            my( $element ) = $other->heteroatoms;
            if( $element ne 'O' ) {
                $name .= $elements{$element}->{prefix};
                $name->[-1] =~ s/a$/o/;
            }
            return $name . 'pyran';
        } else {
            my  $other_name = $other->suffix;
            if( $other_name->starts_with_locant ) { # Locants are moved to front
                unshift @$name, shift @$other_name;
            }
            $name->[-1] =~ s/o$// if $other_name->[0] =~ /^a/;
            return $name . $other_name;
        }
    }

    # Fusion naming according to BBv2 P-25.3.1.3

    # Find the bridge vertices
    my @bridge = $self->{cycles}[0]->vertices;
    @bridge = @bridge[-2..-1];

    # Find autosymmetric equivalents having the least locants for the bridge
    my @equiv_A = $self->{cycles}[0]->autosymmetric_equivalents;
    my @equiv_B = $self->{cycles}[1]->autosymmetric_equivalents;

    my( $A_min, $A_max ) = sort map { $_->vertex_ids( @bridge ) } @equiv_A;
    @equiv_A = grep { min( $_->vertex_ids( @bridge ) ) == $A_min } @equiv_A;

    my( $B_min, $B_max ) = sort map { $_->vertex_ids( @bridge ) } @equiv_B;
    @equiv_B = grep { min( $_->vertex_ids( @bridge ) ) == $B_min } @equiv_B;

    my $fusion;
    if( $self->{cycles}[1]->is_homogeneous &&
        $self->{cycles}[1]->number_of_branches == 2 ) {
        $fusion = '[';
    } elsif( @equiv_A > 1 || @equiv_B > 1 ) {
        # At least one of the rings has mirror symmetry ("flip-symmetric"), thus numeric order is ascending
        $fusion = '[' . ($B_min+1) . ',' . ($B_min+2) . '-';
    } else {
        # Rings are rigid, thus numeric order has to be derived
        my @order_A = $equiv_A[0]->vertex_ids( @bridge );
        my @order_B = $equiv_B[0]->vertex_ids( @bridge );
        if( ($order_A[0] <=> $order_A[1]) == ($order_B[0] <=> $order_B[1]) ) {
            # Ring atoms are encountered in the same order in both of the rings
            $fusion = '[' . ($B_min+1) . ',' . ($B_max+1);
        } else {
            # Ring atom orders differ
            $fusion = '[' . ($B_max+1) . ',' . ($B_min+1);
        }
        $fusion .= '-';
    }
    $fusion .= chr( 97 + $A_min ) . ']';

    my @ideal = map { ChemOnomatopist::Chain::Monocycle->new( $_->graph, $_->vertices ) }
                    $self->cycles;
    my $name_A = $ideal[1]->suffix;
    $name_A =~ s/^\d+H-//;

    # Retained prefixes from BBv3 P-25.3.2.2.3
    $name_A = 'anthra'    if $name_A eq 'anthracene';
    $name_A = 'naphth'    if $name_A eq 'naphthalene';
    $name_A = 'benz'      if $name_A eq 'benzene';
    $name_A = 'phenanthr' if $name_A eq 'phenanthrene';
    $name_A = 'fur'       if $name_A eq 'furan';
    $name_A = 'imidaz'    if $name_A eq 'imidazole';
    $name_A = 'pyrid'     if $name_A eq 'pyridine';
    $name_A = 'pyrimid'   if $name_A eq 'pyrimidine';
    $name_A = 'thien'     if $name_A eq 'thiophene';

    $name_A = ChemOnomatopist::Name->new( $name_A ) unless blessed $name_A;
    $name_A->pop_e;
    if( $name_A->ends_with_alkane_an_suffix ) {
        pop @$name_A;
        $name_A .= 'a';
    } else {
        $name_A .= 'o';
    }

    my $name = ChemOnomatopist::Name->new;
    $name .= $name_A;
    $name .= ChemOnomatopist::Name::Part::Fusion->new( $fusion );

    my $name_B = $ideal[0]->suffix;
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

sub rule_greatest_variety_of_heteroatoms { all_max { scalar uniq $_->heteroatoms } @_ }

sub rule_lowest_numbered_fusion_carbons { all_min { ChemOnomatopist::Comparable::Array::Numeric->new( $_->fusion_carbon_positions ) } @_ }

sub rule_lowest_numbered_fusion_heteroatoms { all_min { ChemOnomatopist::Comparable::Array::Numeric->new( $_->fusion_heteroatom_positions ) } @_ }

sub _adjust_vertices_to_cycles()
{
    my( $self ) = @_;

    my @cycles = $self->cycles;

    my @vertices;
    push @vertices, $cycles[0]->vertices;
    pop  @vertices;
    push @vertices, $cycles[1]->vertices;
    pop  @vertices;
    $self->vertices( @vertices );

    return $self;
}

sub find_retained
{
    my @SMILES = map { s/[=:]//g; uc } @_;
    return first { ($_->[0] eq $SMILES[0] && $_->[1] eq $SMILES[1]) ||
                   ($_->[0] eq $SMILES[1] && $_->[1] eq $SMILES[0]) } @names;
}

1;
