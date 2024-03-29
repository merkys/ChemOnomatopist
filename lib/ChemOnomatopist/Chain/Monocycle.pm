package ChemOnomatopist::Chain::Monocycle;

# ABSTRACT: Monocyclic group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Sulfinyl;
use ChemOnomatopist::Group::Sulfonyl;
use ChemOnomatopist::Name;
use ChemOnomatopist::Util qw( cmp_arrays );
use List::Util qw( all );
use Scalar::Util qw( blessed );

# From BBv2 P-22.2.1
our %names = (
    CCNCO => '1,3-oxazolidine',
    CCCNO => '1,2-oxazolidine',
    CCCCN => 'pyrrolidine',
    CCCNN => 'pyrazolidine',
    CCNCN => 'imidazolidine',

    CCCCCN => 'piperidine',
    CCNCCN => 'piperazine',
    CCNCCO => 'morpholine',

    'CCNCCS'    => 'thiomorpholine',
    'CCNCC[Se]' => 'selenomorpholine',
    'CCNCC[Te]' => 'telluromorpholine',

    # 5-membered aromatic (monoheteroatoms are handled elsewhere)

    'C=CN=CN'  => '1H-imidazole', # FIXME: Adjust for isomerism
    'C=NCCN'   => '1H-imidazole', # 4,5-dihydro-1H-imidazole

    'C=CN=CO'  => '1,3-oxazole',
    'C=CC=NO'  => '1,2-oxazole',
    'C=CC=NN'  => 'pyrazole', # FIXME: Adjust for isomerism

    'C=CC=C[Se]' => 'selenophene',
    'C=CC=C[Te]' => 'tellurophene',

    # 6-membered aromatic
    'c:c:c:c:c:c:' => 'benzene',

    'C=CC=CCO'     => '2H-pyran', # FIXME: Adjust for isomerism
    'C=CC=CCS'     => '2H-thiopyran', # FIXME: Adjust for isomerism
    'C=CC=CC[Se]'  => '2H-selenopyran', # FIXME: Adjust for isomerism
    'C=CC=CC[Te]'  => '2H-telluropyran', # FIXME: Adjust for isomerism

    'c:c:n:c:c:n:' => 'pyrazine',
    'c:c:c:c:n:n:' => 'pyridazine',
    'c:c:c:c:c:n:' => 'pyridine',
    'c:c:c:n:c:n:' => 'pyrimidine',
);

sub new
{
    my( $class, $graph, @vertices ) = @_;
    my $self = bless { graph => $graph, vertices => \@vertices }, $class;
    return $self if $self->is_homogeneous;

    ( $self ) = $self->autosymmetric_equivalents;
    return $self;
}

# FIXME: For now we generate all possible traversals of the same cycle.
#        This is not optimal, some caching could be introduced.
sub candidates()
{
    my( $self ) = @_;

    my $graph = $self->graph;
    my @vertices = $self->vertices;

    my @chains;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Monocycle->new( $graph, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Monocycle->new( $graph, @vertices );
        push @vertices, shift @vertices;
    }

    for (@chains) {
        $_->{candidate_for} = $self;
    }

    return @chains;
}

sub autosymmetric_equivalents()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    my $cycle = $self->graph->subgraph( \@vertices ); # TODO: Add attributes
    my @chains;
    for (0..$#vertices) {
        push @chains,
             ChemOnomatopist::Chain::Circular->new( $cycle, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains,
             ChemOnomatopist::Chain::Circular->new( $cycle, @vertices );
        push @vertices, shift @vertices;
    }
    # CHECKME: Additional rules from ChemOnomatopist::filter_chains() might still be needed
    @chains = sort {  ChemOnomatopist::Chain::Monocycle::_cmp( $a, $b ) } @chains;
    @chains = grep { !ChemOnomatopist::Chain::Monocycle::_cmp( $_, $chains[0] ) } @chains;
    @chains = map  { bless { graph => $self->graph, vertices => [ $_->vertices ] } } @chains;
    return ChemOnomatopist::rule_lowest_numbered_locants( @chains );
}

sub parent(;$)
{
    my( $self, $parent ) = @_;
    my $old_parent = $self->SUPER::parent( $parent );
    return $old_parent unless $parent;
    return $old_parent if $old_parent && $parent == $old_parent;

    # Addition of parent to homogeneous cycles settles the otherwise ambiguous order
    # TODO: Other autosymmetric monocycles can possibly as well be settled
    if( $self->is_homogeneous ) {
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $vertices[$_], $parent ) } 0..$#vertices;
        if( defined $position ) {
            my @chains =
                ( ChemOnomatopist::Chain::Circular->new( $self->graph,
                                                         @vertices[$position..$#vertices],
                                                         @vertices[0..$position-1] ) );
            @vertices = reverse @vertices;
            $position = $#vertices - $position;
            push @chains,
                 ChemOnomatopist::Chain::Circular->new( $self->graph,
                                                        @vertices[$position..$#vertices],
                                                        @vertices[0..$position-1] );
            for (@chains) {
                $_->{parent} = $parent;
            }

            my( $chain ) = ChemOnomatopist::filter_chains( @chains );
            $self->{vertices} = [ $chain->vertices ];
        }
    }

    return $old_parent;
}

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return '' if $self->is_hydrocarbon;
    return '' if $self->is_monoreplaced && !$self->is_substituted;
    return $self->length < 3 || $self->length > 10 || all { $_->{symbol} !~ /^[cC]$/ } $self->vertices;
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;
    return '' if $self->is_hydrocarbon;
    return $self->length < 3 || $self->length > 10 || all { $_->{symbol} !~ /^[cC]$/ } $self->vertices;
}

sub prefix()
{
    my( $self ) = @_;

    my $parent = $self->parent;
    my $name = $self->suffix;
    if( $name eq 'benzene' ) {
        if( $parent && blessed $parent &&
            ( $parent->isa( ChemOnomatopist::Group::Sulfinyl:: ) ||
              $parent->isa( ChemOnomatopist::Group::Sulfonyl:: ) ) ) {
            # Rule derived from examples in BBv2 P-63.6
            return $name;
        }
        return 'phenyl';
    }

    $name = ChemOnomatopist::Name->new( $name ) unless blessed $name;
    $name->pop_e;
    pop @$name if $name->ends_with_alkane_an_suffix;

    if( $parent && !$self->is_homogeneous ) {
        # FIXME: Order of vertices seems to be established regardless of the attachments.
        # This causes ambiguity, for example, in 1-(1,4-diazepan-1-ylsulfonyl)-8-methylisoquinoline
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }

    $name .= 'yl';

    return $name;
}

sub suffix()
{
    my( $self ) = @_;
    my $name = $self->name;
    return blessed $name ? $name : ChemOnomatopist::Name->new( $name );
}

# FIXME: Pay attention to bond orders
sub _cmp
{
    my( $A, $B ) = @_;

    my @A_heteroatoms = $A->heteroatoms;
    my @A_positions = $A->heteroatom_positions;

    @A_positions = map { $A_positions[$_] }
                       sort { $elements{$A_heteroatoms[$a]}->{seniority} <=>
                              $elements{$A_heteroatoms[$b]}->{seniority} } 0..$#A_positions;

    my @B_heteroatoms = $B->heteroatoms;
    my @B_positions = $B->heteroatom_positions;

    @B_positions = map { $B_positions[$_] }
                       sort { $elements{$B_heteroatoms[$a]}->{seniority} <=>
                              $elements{$B_heteroatoms[$b]}->{seniority} } 0..$#B_positions;

    return cmp_arrays( \@A_positions, \@B_positions )
        if cmp_arrays( \@A_positions, \@B_positions );

    return cmp_arrays( [ $A->multiple_bond_positions ], [ $B->multiple_bond_positions ] );
}

1;
