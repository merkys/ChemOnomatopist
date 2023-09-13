package ChemOnomatopist::Group::Monocycle;

use strict;
use warnings;

# ABSTRACT: Monocyclic group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Sulfinyl;
use ChemOnomatopist::Group::Sulfonyl;
use ChemOnomatopist::Name;
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

    # 5-membered aromatic
    'C=CC=CO' => 'furan',
    'C=COc:c' => 'furan', # For fused rings

    'C=CN=CN'  => '1H-imidazole', # FIXME: Adjust for isomerism
    'C=CN=CO'  => '1,3-oxazole',
    'C=CC=NO'  => '1,2-oxazole',
    'C=CC=NN'  => '1H-pyrazole', # FIXME: Adjust for isomerism

    'C=CC=CN'  => '1H-pyrrole', # FIXME: Adjust for isomerism
    'C=CNc:c'  => '1H-pyrrole', # For fused rings

    'C=CC=C[Se]' => 'selenophene',
    'C=CC=C[Te]' => 'tellurophene',
    'C=CC=CS' => 'thiophene',
    'C=CSc:c' => 'thiophene', # For fused rings

    # 6-membered aromatic
    'c:c:c:c:c:c:' => 'benzene',

    'C=CC=CCO'     => '2H-pyran', # FIXME: Adjust for isomerism

    'c:c:n:c:c:n:' => 'pyrazine',
    'c:c:c:c:n:n:' => 'pyridazine',
    'c:c:c:c:c:n:' => 'pyridine',
    'c:c:c:n:c:n:' => 'pyrimidine',
);

sub new
{
    my( $class, $graph, @vertices ) = @_;
    my $self = bless { graph => $graph, vertices => \@vertices }, $class;
    $self->_aromatise;
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
        push @chains, ChemOnomatopist::Group::Monocycle->new( $graph, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Group::Monocycle->new( $graph, @vertices );
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
    @chains = sort {  ChemOnomatopist::Group::Monocycle::_cmp( $a, $b ) } @chains;
    @chains = grep { !ChemOnomatopist::Group::Monocycle::_cmp( $_, $chains[0] ) } @chains;
    return map { bless { graph => $self->graph, vertices => [ $_->vertices ] } }  @chains;
}

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return $self->length < 3 || $self->length > 10 || all { $_->{symbol} !~ /^[cC]$/ } $self->vertices;
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;
    return $self->needs_heteroatom_locants;
}

sub prefix(;$)
{
    my( $self, $parent ) = @_;

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
    $name->{name}[-1] =~ s/e$//;
    pop @$name if $name->{name}[-1] eq '';
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

# FIXME: This is a bit strange: class and object method with the same name
sub suffix()
{
    my( $self ) = @_;
    return '' unless ref $self;
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

    return ChemOnomatopist::cmp_arrays( \@A_positions, \@B_positions )
        if ChemOnomatopist::cmp_arrays( \@A_positions, \@B_positions );

    return ChemOnomatopist::cmp_arrays( [ $A->multiple_bond_positions ], [ $B->multiple_bond_positions ] );
}

1;
