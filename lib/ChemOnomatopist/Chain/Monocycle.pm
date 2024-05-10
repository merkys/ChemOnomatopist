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
use List::Util qw( all first );
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

    'C=CN=CN'  => 'imidazole',
    'C=NCCN'   => 'imidazole', # 4,5-dihydro-1H-imidazole

    'C=CN=CO'  => '1,3-oxazole',
    'C=CC=NO'  => '1,2-oxazole',
    'C=CC=NN'  => 'pyrazole', # FIXME: Adjust for isomerism

    'C=CC=C[Se]' => 'selenophene',
    'C=CC=C[Te]' => 'tellurophene',

    # 6-membered aromatic
    'c:c:c:c:c:c:' => 'benzene',

    'C=CC=CCO'     => 'pyran',
    'C=CC=CCS'     => 'thiopyran',
    'C=CC=CC[Se]'  => 'selenopyran',
    'C=CC=CC[Te]'  => 'telluropyran',

    'c:c:n:c:c:n:' => 'pyrazine',
    'c:c:c:c:n:n:' => 'pyridazine',
    'c:c:c:c:c:n:' => 'pyridine',
    'c:c:c:n:c:n:' => 'pyrimidine',

    # Various cases of aromatisation
    'C=CCC=CN' => 'pyridine',
    'C=CC=CCN' => 'pyridine',
);

sub new
{
    my( $class, $graph, @vertices ) = @_;
    my $self = bless { graph => $graph, vertices => [ @vertices ] }, $class;
    return $self if $self->is_homogeneous;

    # TODO: This code is not optimal, but works
    my( $senior_heteroatom ) = sort { $elements{$a}->{seniority} <=>
                                      $elements{$b}->{seniority} }
                               grep { $_ ne 'C' }
                               map  { ChemOnomatopist::element( $_ ) }
                                    @vertices;
    return $self unless $senior_heteroatom;

    my @chains;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $graph, @vertices )
            if ChemOnomatopist::element( $vertices[0] ) eq $senior_heteroatom;
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $graph, @vertices )
            if ChemOnomatopist::element( $vertices[0] ) eq $senior_heteroatom;
        push @vertices, shift @vertices;
    }

    my( $first ) = sort { _cmp( $a, $b ) } @chains;
    return bless { graph => $graph, vertices => [ $first->vertices ] }, $class;
}

sub candidates()
{
    my( $self ) = @_;

    my $graph = $self->graph;
    my @vertices = $self->vertices;
    my $parent = $self->parent;
    my( $senior_heteroatom ) = $self->heteroatoms;

    my @chains;
    for (0..$#vertices) {
        push @chains, bless { graph => $graph, vertices => [ @vertices ], parent => $parent }, ChemOnomatopist::Chain::Monocycle::
            if !$senior_heteroatom || ChemOnomatopist::element( $vertices[0] ) eq $senior_heteroatom;
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, bless { graph => $graph, vertices => [ @vertices ], parent => $parent }, ChemOnomatopist::Chain::Monocycle::
            if !$senior_heteroatom || ChemOnomatopist::element( $vertices[0] ) eq $senior_heteroatom;
        push @vertices, shift @vertices;
    }

    my( $max_value ) = sort { _cmp( $a, $b ) } @chains;
    @chains = grep { !_cmp( $_, $max_value ) } @chains;

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

    # Addition/change of parent may need resetting the numbering
    my @candidates = $self->candidates;
    @candidates = rule_lowest_parent_locant( @candidates );
    @candidates = ChemOnomatopist::filter_chains( @candidates );

    $self->{vertices} = [ $candidates[0]->vertices ];
    return $old_parent;
}

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return '' if $self->is_hydrocarbon;
    return '' if $self->is_monoreplaced && !$self->is_substituted;
    return $self->length < 3 || $self->length > 10 || all { ChemOnomatopist::element( $_ ) ne 'C' } $self->vertices;
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;
    return '' if $self->is_hydrocarbon;
    return $self->length < 3 || $self->length > 10 || all { ChemOnomatopist::element( $_ ) ne 'C' } $self->vertices;
}

sub needs_indicated_hydrogens()
{
    my( $self ) = @_;
    return '' if $self->is_hydrocarbon && $self->length == 6; # BBv3 P-31.2.3.1
    return '' if $self->is_saturated; # BBv3 P-31.2.3.2
    return '' unless $self->is_Hantzsch_Widman;
    return  1;
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
        return ChemOnomatopist::Name->new( 'phenyl' );
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

sub rule_lowest_parent_locant
{
    my @chains = @_;

    my( $max_value ) = sort { $a->parent_locant <=> $b->parent_locant }
                            @chains;
    return grep { $_->parent_locant == $max_value->parent_locant } @chains;
}

# FIXME: Pay attention to bond orders
# Orders the atoms by their seniority and compares the resulting locants
sub _cmp
{
    my( $A, $B ) = @_;

    my @A_positions = $A->heteroatom_positions;
    my @B_positions = $B->heteroatom_positions;

    return cmp_arrays( \@A_positions, \@B_positions )
        if cmp_arrays( \@A_positions, \@B_positions );

    my @A_heteroatoms = $A->heteroatoms;
    my @B_heteroatoms = $B->heteroatoms;

    @A_positions = map { $A_positions[$_] }
                       sort { $elements{$A_heteroatoms[$a]}->{seniority} <=>
                              $elements{$A_heteroatoms[$b]}->{seniority} } 0..$#A_positions;

    @B_positions = map { $B_positions[$_] }
                       sort { $elements{$B_heteroatoms[$a]}->{seniority} <=>
                              $elements{$B_heteroatoms[$b]}->{seniority} } 0..$#B_positions;

    return cmp_arrays( \@A_positions, \@B_positions );
}

1;
