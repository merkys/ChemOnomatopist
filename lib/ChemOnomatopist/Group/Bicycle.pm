package ChemOnomatopist::Group::Bicycle;

use strict;
use warnings;

# ABSTRACT: Fused bicyclic group
# VERSION

use ChemOnomatopist;
use ChemOnomatopist::Chain::Circular;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Monocycle;
use ChemOnomatopist::Group::Monocycle::Fused;
use ChemOnomatopist::Name;
use ChemOnomatopist::Name::Part::Stem;
use ChemOnomatopist::Util::SMILES qw( cycle_SMILES );
use Graph::Traversal::DFS;
use List::Util qw( all any );
use Set::Object qw( set );

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

# From BBv2 P-25.2.1
our @names = (
    [ 'n:c:n:c:c:c:', 'n:c:c:n:c:c:', 'pteridine' ],
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

    [ 'c:n:c:n:c:c:', 'N=CNc:c', 'purine' ], # Special rules apply
    [ 'c:n:c:n:c:c:', 'NC=Nc:c', 'purine' ], # Special rules apply

    [ 'NN=Cc:c', 'c:c:c:c:c:c:', '1H-indazole' ],
    [ 'NC=Cc:c', 'c:c:c:c:c:c:', '1H-indole' ],
    [ 'CNC=CC=', 'C=CC=CCC',  'isoindole' ],
    [ 'CC=CNC=', 'C=CC=CCN',  'indolizine', ],
    [ 'CC=CNC',  'C=CC=CN',   '1H-pyrrolizine' ], # TODO: There are isomers
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
        $P_parts[0] =~ s/N/P/g;
        $P_parts[1] =~ s/N/P/g;
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
    my @cycles = map { ChemOnomatopist::Group::Monocycle::Fused->new( $graph, $self, @$_ ) }
                     @components;
    $self->{cycles} = \@cycles;

    $self->_aromatise;

    my $nbenzene = scalar grep { $_->is_benzene } @cycles;

    # The ordering should not be done if one of the cycles is benzene
    if( $nbenzene == 0 ) {
        my @flipped = map { $_->flipped } @cycles;
        # CHECKME: Additional rules from ChemOnomatopist::filter_chains() might still be needed
        my( $chain ) = sort { ChemOnomatopist::Group::Monocycle::_cmp( $a, $b ) } ( @cycles, @flipped );

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

        my( $chain ) = sort { ChemOnomatopist::Group::Monocycle::_cmp( $a, $b ) }
                            ( $cycles[0], $cycles[0]->flipped );

        if( $chain != $cycles[0] ) {
            @cycles = map { $_->flipped } @cycles;
            $self->{cycles} = \@cycles;
        }
        $self->_adjust_vertices_to_cycles;
    }

    if( join( ',', map { $_->backbone_SMILES } @cycles ) =~ /^n:c:n:c:c:c:,N(C=|=C)Nc:c$/ ) {
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

sub prefix()
{
    my( $self ) = @_;

    my $name = ChemOnomatopist::Name->new( $self->suffix );
    $name->{name}[-1] =~ s/e$//;
    if( $self->parent ) { # FIXME: Not stable for naphthalene
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $self->parent, $vertices[$_] ) } 0..$#vertices;
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
        return $hydrocarbons_by_size{$cycle_sizes} if exists $hydrocarbons_by_size{$cycle_sizes};

        if( $cycle_sizes =~ /^(\d+),\1$/ ) {
            my $name = ChemOnomatopist::alkane_chain_name( $1 ) . 'alene';
            return ChemOnomatopist::Name::Part::Stem->new( $name )->to_name;
        }
    }

    my @SMILES = map { $_->backbone_SMILES } $self->cycles;
    my( $retained ) = grep { ($_->[0] eq $SMILES[0] && $_->[1] eq $SMILES[1]) ||
                             ($_->[0] eq $SMILES[1] && $_->[1] eq $SMILES[0]) } @names;
    return ChemOnomatopist::Name::Part::Stem->new( $retained->[2] )->to_name if $retained;

    if( any { $_->is_benzene } $self->cycles ) {
        my( $other ) = grep { !$_->is_benzene } $self->cycles;
        $other = ChemOnomatopist::Group::Monocycle->new( $other->graph, $other->vertices );

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

    # TODO: Complete implementing BBv2 P-25.3.1.3 (fusion naming)
    my @cycles = $self->cycles;

    for my $rule ( # TODO: P-25.3.2.4 (a): Senior heteroatom
                   # TODO: P-25.3.2.4 (b): Concerns fusions of more than two rings
                   # P-25.3.2.4 (c): Second ring has to be larger
                   \&ChemOnomatopist::rule_longest_chains,
                   # P-25.3.2.4 (d): Greater number of heteroatoms of any kind
                   \&ChemOnomatopist::rule_most_heteroatoms,
                   # TODO: P-25.3.2.4 (e): Greater variety of heteroatoms
                   # P-25.3.2.4 (f): Greater number of most senior heteroatoms
                   \&ChemOnomatopist::rule_most_senior_heteroatoms,
                   # TODO: P-25.3.2.4 (g): Concerns fusions of more than two rings
                   # P-25.3.2.4 (h): Lower locants for heteroatoms
                   \&ChemOnomatopist::rule_lowest_numbered_heteroatoms,
                   # P-25.3.2.4 (i): Lower locants for senior heteroatoms
                   \&ChemOnomatopist::rule_lowest_numbered_most_senior_heteroatoms,
                   # TODO: P-25.3.2.4 (j): Concerns fusions of more than two rings
                 ) {
        my @cycles_now = $rule->( @cycles );
        last unless @cycles_now; # Did not succeed, quit
        if( @cycles_now == 1 ) { # Filtering completed
            @cycles = reverse @cycles if $cycles[0] == $cycles_now[0];
            last;
        }
    }

    my @ideal = map { ChemOnomatopist::Group::Monocycle->new( $_->graph, reverse $_->vertices ) } @cycles;

    # TODO: These are ad-hoc rules as for the moment generalisation is hard to make
    my $fusion = '';

    if(     $ideal[0]->{vertices}[0] == $cycles[0]->{vertices}[0] ) {
        if( $ideal[0]->{vertices}[1] == $cycles[0]->{vertices}[1] ) {
            $fusion .= '[' . join( ',', $ideal[0]->length, $ideal[0]->length - 1 );
        }
        if( $ideal[0]->{vertices}[1] == $cycles[0]->{vertices}[-1] ) {
            $fusion .= '[3,2';
        }
    } else {
        my $flipped = $cycles[0]->flipped;
        if( $ideal[0]->{vertices}[0] == $flipped->{vertices}[0] ) {
            if( $ideal[0]->{vertices}[1] == $flipped->{vertices}[1] ) {
                $fusion .= '[2,3';
            }
            if( $ideal[0]->{vertices}[1] == $cycles[0]->{vertices}[-1] ) {
                $fusion .= '[' . join( ',', $ideal[0]->length - 1, $ideal[0]->length );
            }
        }
    }

    if(     $ideal[1]->{vertices}[0] == $cycles[1]->{vertices}[0] ) {
        if( $ideal[1]->{vertices}[1] == $cycles[1]->{vertices}[1] ) {
            $fusion .= '-' . chr( 95 + $cycles[1]->length ) . ']';
        }
        if( $ideal[1]->{vertices}[1] == $cycles[1]->{vertices}[-1] ) {
            $fusion .= '-b]';
        }
    } else {
        my $flipped = $cycles[1]->flipped;
        if( $ideal[1]->{vertices}[0] == $flipped->{vertices}[0] ) {
            if( $ideal[1]->{vertices}[1] == $flipped->{vertices}[1] ) {
                $fusion .= '-' . chr( 95 + $cycles[1]->length ) . ']';
            }
        }
    }

    die "cannot name complex bicyclic compounds\n" unless $fusion =~ /^\[.+\]$/; # Full fusion is known

    my $name = $ideal[0]->name;
    # TODO: Preserve retained prefixes from BBv2 P-25.3.2.2.3
    unless( $name->[-1] =~ s/e$/o/ ) { # BBv2 P-25.3.2.2.2
        $name->[-1] .= 'o';
    }
    $name .= $fusion;
    $name .= $ideal[1]->name;
    $name->bracket_numeric_locants;
    return $name;
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

sub _aromatise()
{
    my( $self ) = @_;

    my @delocalised_cycles = grep { join( '', $_->bonds ) =~ /^((-=)+|(=-)+)$/ }
                                  ( $self, $self->cycles );
    return '' unless @delocalised_cycles;

    my $subgraph = $self->graph->subgraph( [ map { $_->vertices } @delocalised_cycles ] );
    for ($subgraph->vertices) {
        next unless $_->{symbol} =~ /^(Se|As|[BCNOPS])$/;
        $_->{symbol} = lcfirst $_->{symbol};
    }
    for ($subgraph->edges) {
        $self->graph->set_edge_attribute( @$_, 'bond', ':' );
    }
    for ($self, $self->cycles) { # Need to invalidate bond cache
        delete $_->{bonds};
    }

    return 1;
}

1;
