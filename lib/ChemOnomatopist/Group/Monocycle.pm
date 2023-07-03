package ChemOnomatopist::Group::Monocycle;

use strict;
use warnings;

# ABSTRACT: Monocyclic group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist;
use ChemOnomatopist::Elements qw( %elements );
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

    'C=CC=NC'  => '1H-imidazole', # FIXME: Adjust for isomerism
    'C=CN=CO'  => '1,3-oxazole',
    'C=CC=NO'  => '1,2-oxazole',
    'C=CC=NN'  => '1H-pyrazole', # FIXME: Adjust for isomerism
    'C=CC=CN'  => '1H-pyrrole', # FIXME: Adjust for isomerism
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

    # Select the best numbering for heteroatoms
    my @chains;
    my $cycle = $self->graph->subgraph( \@vertices ); # TODO: Add attributes
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
    my( $chain ) = sort { ChemOnomatopist::Group::Monocycle::_cmp( $a, $b ) } @chains;
    return bless { graph => $graph, vertices => [ $chain->vertices ] }, $class;
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
    return 'phenyl' if $name eq 'benzene';

    $name = ChemOnomatopist::Name->new( $name ) unless blessed $name;
    $name->{name}[-1] =~ s/(an)?e$//;

    if( $parent && !$self->is_homogeneous ) {
        my @vertices = $self->vertices;
        my( $position ) = grep { $self->graph->has_edge( $parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }

    return $name . 'yl';
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
