package ChemOnomatopist::Chain::Monocycle;

# ABSTRACT: Monocyclic group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Sulfinyl;
use ChemOnomatopist::Group::Sulfonyl;
use ChemOnomatopist::Name;
use ChemOnomatopist::Util qw( all_min circle_permutations cmp_arrays element );
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

my %five_membered_aromatic_single_heteroatom = (
    N => 'pyrrole',
    O => 'furan',
    S => 'thiophene',
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
                               map  { element( $_ ) }
                                    @vertices;
    return $self unless $senior_heteroatom;

    my @chains;
    for (circle_permutations( @vertices )) {
        next unless element( $_->[0] ) eq $senior_heteroatom;
        push @chains, ChemOnomatopist::Chain::Circular->new( $graph, @$_ );
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
            if !$senior_heteroatom || element( $vertices[0] ) eq $senior_heteroatom;
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, bless { graph => $graph, vertices => [ @vertices ], parent => $parent }, ChemOnomatopist::Chain::Monocycle::
            if !$senior_heteroatom || element( $vertices[0] ) eq $senior_heteroatom;
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

    my $cycle = $self->graph->subgraph( $self->vertices ); # TODO: Add attributes
    my @chains = map { ChemOnomatopist::Chain::Circular->new( $cycle, @$_ ) }
                     circle_permutations( $self->vertices );
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
    @candidates = rule_lowest_parent_locant( @candidates ) if defined $self->parent_locant;
    @candidates = ChemOnomatopist::filter_chains( @candidates );

    $self->{vertices} = [ $candidates[0]->vertices ];
    return $old_parent;
}

sub needs_heteroatom_locants()
{
    my( $self ) = @_;
    return '' if $self->is_hydrocarbon;
    return '' if $self->is_monoreplaced && !$self->is_substituted;
    return $self->length < 3 || $self->length > 10 || all { element( $_ ) ne 'C' } $self->vertices;
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;
    return '' if $self->is_hydrocarbon;
    return $self->length < 3 || $self->length > 10 || all { element( $_ ) ne 'C' } $self->vertices;
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
        my $position = first { $self->graph->has_edge( $parent, $vertices[$_] ) } 0..$#vertices;
        die "unknown locant in multicyclic compound\n" unless defined $position;
        $name->append_substituent_locant( $self->locants( $position ) );
    }

    $name .= 'yl';

    return $name;
}

sub suffix()
{
    my( $self ) = @_;

    my $graph = $self->graph;
    my $SMILES = $self->backbone_SMILES;

    # Check the preserved names
    if( $self->length == 5 && $self->number_of_double_bonds &&
        exists $five_membered_aromatic_single_heteroatom{join( '', $self->heteroatoms )} ) {
        return ChemOnomatopist::Name->new( $five_membered_aromatic_single_heteroatom{join( '', $self->heteroatoms )} );
    }
    return ChemOnomatopist::Name->new( $names{$SMILES} ) if exists $names{$SMILES};

    # Check for aromatic notation
    if( $SMILES =~ /:/ ) {
        $SMILES =~ s/([a-z]):/'' . uc( $1 )/ge;
        $SMILES =~ s/\[([a-z]{1,2})\]:/'[' . uc( $1 ) . ']'/ge;
        for my $SMILES_for_name (keys %names) {
            next unless $SMILES_for_name =~ /=/;
            my $name = $names{$SMILES_for_name};
            $SMILES_for_name =~ s/=//g;
            return ChemOnomatopist::Name->new( $name ) if $SMILES eq $SMILES_for_name;
        }
    }

    # Check for annulenes
    if( $self->is_hydrocarbon && $self->is_aromatic &&
        $self->length =~ /^(4|6|8|10|12|14|16)$/ ) {
        return ChemOnomatopist::Name->new(
                    'cyclo' .
                    ChemOnomatopist::IUPAC_numerical_multiplier( $self->length, 1 ) .
                    ChemOnomatopist::IUPAC_numerical_multiplier( $self->length / 2, 1 ) . 'ene' );
    }

    # Check for cycloalkanes
    if( $self->is_hydrocarbon ) {
        my $name = ChemOnomatopist::Name::Part::NondetachablePrefix->new( 'cyclo' )->to_name;
        $name .= $self->SUPER::suffix;
        return $name;
    }

    if( $self->is_Hantzsch_Widman ) {
        # Hantzsch-Widman names (BBv2 P-22.2.2.1)

        # Collect the types of heteroatoms and their attachment positions
        my %heteroatoms;
        my @vertices = $self->vertices;
        for my $i (0..$#vertices) {
            my $symbol = element( $vertices[$i] );
            next if $symbol eq 'C';
            $heteroatoms{$symbol} = [] unless $heteroatoms{$symbol};
            push @{$heteroatoms{$symbol}}, $i;
        }

        my $least_senior_element;
        my @heteroatom_locants;
        for my $element (sort { $elements{$a}->{seniority} <=> $elements{$b}->{seniority} }
                              keys %heteroatoms) {
            push @heteroatom_locants, @{$heteroatoms{$element}};
            $least_senior_element = $element;
        }

        my $name = ChemOnomatopist::Name->new;
        unless(  @heteroatom_locants == 1 ||
                (scalar keys %heteroatoms == 1 && @heteroatom_locants == $self->length - 1) ) {
            # Locants are omitted according to BBv2 P-22.2.2.1.7
            $name->append_locants( map { $_ + 1 } @heteroatom_locants );
        }

        for my $element (sort { $elements{$a}->{seniority} <=> $elements{$b}->{seniority} }
                              keys %heteroatoms) {
            if( @{$heteroatoms{$element}} > 1 ) {
                $name->append_multiplier( ChemOnomatopist::IUPAC_numerical_multiplier( scalar @{$heteroatoms{$element}} ) );
            }
            $name->append_element( exists $elements{$element}->{HantzschWidman}
                                        ? $elements{$element}->{HantzschWidman}
                                        : $elements{$element}->{prefix} );
        }
        $name->[-1] =~ s/a$//;

        if(      $self->length <= 5 ) {
            my @stems = ( 'ir', 'et', 'ol' );
            $name->append_stem( $stems[$self->length - 3] );
            if( $self->is_saturated ) {
                $name .= $heteroatoms{N} ? 'idine' : 'ane';
            } elsif( $self->length == 3 ) {
                $name .= $heteroatoms{N} ? 'ene' : 'ine';
            } else {
                $name .= 'e';
            }
            return $name;
        } elsif( $self->length == 6 ) {
            if(      ($elements{$least_senior_element}->{seniority} >= 5 &&
                      $elements{$least_senior_element}->{seniority} <= 8) || $least_senior_element eq 'Bi' ) {
                $name .= $self->is_saturated ? 'ane' : 'ine';
            } elsif( ($elements{$least_senior_element}->{seniority} >= 16 &&
                      $elements{$least_senior_element}->{seniority} <= 19) || $least_senior_element eq 'N' ) {
                $name .= $self->is_saturated ? 'inane' : 'ine';
            } else {
                $name .= $self->is_saturated ? 'inane' : 'inine';
            }
            return $name;
        } elsif( $self->length >= 7 ) {
            my @stems = ( 'ep', 'oc', 'on', 'ec' );
            $name->append_stem( $stems[$self->length - 7] );
            $name .= $self->is_saturated ? 'ane' : 'ine';
            return $name;
        }
    }

    my $name = ChemOnomatopist::Name->new( 'cyclo' );
    $name .= $self->SUPER::suffix;
    return $name;
}

sub rule_lowest_parent_locant { all_min { $_->parent_locant } @_ }

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
