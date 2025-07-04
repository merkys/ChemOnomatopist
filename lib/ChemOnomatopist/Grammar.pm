package ChemOnomatopist::Grammar;

# ABSTRACT: Grammar for chemical graphs
# VERSION

use strict;
use warnings;

use Algorithm::Combinatorics qw( combinations );
use ChemOnomatopist::Chain;
use ChemOnomatopist::Chain::ABA;
use ChemOnomatopist::Group::AcidHalide;
use ChemOnomatopist::Group::Amidine;
use ChemOnomatopist::Group::Azide;
use ChemOnomatopist::Group::Carbaldehyde;
use ChemOnomatopist::Group::Carbonitrile;
use ChemOnomatopist::Chain::Carboxamide;
use ChemOnomatopist::Chain::Chalcogen;
use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Elements qw( %elements );
use ChemOnomatopist::Group::Ether;
use ChemOnomatopist::Group::AcylHalide;
use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amide;
use ChemOnomatopist::Group::Amine;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Cyanide;
use ChemOnomatopist::Group::Diazene;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Imine;
use ChemOnomatopist::Group::Isocyanate;
use ChemOnomatopist::Group::Isocyanide;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Group::Nitramide;
use ChemOnomatopist::Group::Nitro;
use ChemOnomatopist::Group::Nitroso;
use ChemOnomatopist::Group::NoncarbonOxoacid;
use ChemOnomatopist::Group::Peroxide;
use ChemOnomatopist::Chain::Sulfimide;
use ChemOnomatopist::Group::SulfinicAcid;
use ChemOnomatopist::Group::Sulfinyl;
use ChemOnomatopist::Group::SulfonicAcid;
use ChemOnomatopist::Group::Sulfonyl;
use ChemOnomatopist::Group::Urea;
use ChemOnomatopist::Group::XO3;
use ChemOnomatopist::Util qw( element );
use Chemistry::OpenSMILES qw(
    is_double_bond
    is_triple_bond
);
use Graph::Grammar;
use List::Util qw( all any first sum );
use Scalar::Util qw( blessed );

use parent Exporter::;
our @EXPORT_OK = qw(
    parse_molecular_graph
);

sub is_nongroup      { !$_[0]->groups( $_[1] ) }
sub is_nongroup_atom { !blessed $_[1] && !$_[0]->groups( $_[1] ) && exists $_[1]->{symbol} }

sub is_C  { element( $_[1] ) && element( $_[1] ) eq 'C' }
sub is_N  { element( $_[1] ) && element( $_[1] ) eq 'N' }
sub is_O  { element( $_[1] ) && element( $_[1] ) eq 'O' }
sub is_S  { element( $_[1] ) && element( $_[1] ) eq 'S' }
sub is_Se { element( $_[1] ) && element( $_[1] ) eq 'Se' }
sub is_Te { element( $_[1] ) && element( $_[1] ) eq 'Te' }

sub is_As_N_B_P_Se_Si_Sb_S_Te { element( $_[1] ) && element( $_[1] ) =~ /^(As|Br|Cl|F|I|Sb|Se|Si|Te|N|B|P|S)$/ }
sub is_Br_Cl_F_I              { element( $_[1] ) && element( $_[1] ) =~ /^(Br|Cl|F|I)$/ }
sub is_Br_Cl_F_I_N            { element( $_[1] ) && element( $_[1] ) =~ /^(Br|Cl|F|I|N)$/ }
sub is_B_Cl_F_I               { element( $_[1] ) && element( $_[1] ) =~ /^(B|Cl|F|I)$/ }
sub is_S_Se_Te                { element( $_[1] ) && element( $_[1] ) =~ /^(S|Se|Te)$/ }
sub is_O_S_Se_Te              { element( $_[1] ) && element( $_[1] ) =~ /^(O|S|Se|Te)$/ }
sub is_C_N_O_S_Se_Te          { element( $_[1] ) && element( $_[1] ) =~ /^(C|N|O|S|Se|Te)$/ }

sub is_heteroatom { element( $_[1] ) && !&is_C }

sub charge_plus_one  {  ChemOnomatopist::charge( $_[1] ) ==  1 }
sub charge_minus_one {  ChemOnomatopist::charge( $_[1] ) == -1 }
sub no_charge        { !ChemOnomatopist::charge( $_[1] ) }

sub has_H0 { !$_[1]->{hcount} }
sub has_H1 {  exists $_[1]->{hcount} && $_[1]->{hcount} == 1 }
sub has_H2 {  exists $_[1]->{hcount} && $_[1]->{hcount} == 2 }
sub has_H3 {  exists $_[1]->{hcount} && $_[1]->{hcount} == 3 }

sub has_1_neighbour { $_[0]->degree( $_[1] ) == 1 }

sub is_aldehyde      { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Aldehyde:: ) }
sub is_amide         { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Amide:: ) }
sub is_amine         { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Amine:: ) }
sub is_cyanide       { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Cyanide:: ) }
sub is_hydroxy       { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Hydroxy:: ) }
sub is_hydroperoxide { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Hydroperoxide:: ) }
sub is_imine         { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Imine:: ) }
sub is_isocyanate    { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Isocyanate:: ) }
sub is_isocyanide    { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Isocyanide:: ) }
sub is_ketone        { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Ketone:: ) }
sub is_nitro         { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Nitro:: ) }
sub is_sulfinyl      { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Sulfinyl:: ) }
sub is_sulfonyl      { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Sulfonyl:: ) }
sub is_XO3           { blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::XO3:: ) }

sub is_benzene   { any { $_->isa( ChemOnomatopist::Chain::Monocycle:: ) && $_->is_benzene } $_[0]->groups( $_[1] ) }
sub is_chalcogen { any { $_->isa( ChemOnomatopist::Chain::Chalcogen:: ) } $_[0]->groups( $_[1] ) }
sub is_circular  { any { $_->isa( ChemOnomatopist::Chain::Circular:: ) } $_[0]->groups( $_[1] ) }
sub is_monocycle { any { $_->isa( ChemOnomatopist::Chain::Monocycle:: ) } $_[0]->groups( $_[1] ) }

sub is_hydrazine { any { $_->isa( ChemOnomatopist::Group::Hydrazine:: ) } $_[0]->groups( $_[1] ) }

sub is_ABA_chain { any { $_->isa( ChemOnomatopist::Chain::ABA:: ) } $_[0]->groups( $_[1] ) }

sub looks_like_ABA_chain
{
    my( $graph, $center ) = @_;
    # Require two neighbours
    my @neighbours = blessed $center ? $center->substituents : $graph->neighbours( $center );
    return '' unless @neighbours == 2;

    if( all { blessed $_ && $_->isa( ChemOnomatopist::Chain::ABA:: ) } @neighbours ) {
        # ABA chains on both sides
        return '' unless $neighbours[0]->outer_element eq $neighbours[1]->outer_element;
        return '' unless $neighbours[0]->inner_element eq $neighbours[1]->inner_element;
        return '' unless $neighbours[0]->inner_element eq element( $center );
    } elsif( any { blessed $_ && $_->isa( ChemOnomatopist::Chain::ABA:: ) } @neighbours ) {
        # ABA chain on one side
        @neighbours = reverse @neighbours if blessed $neighbours[1] &&
                                             $neighbours[1]->isa( ChemOnomatopist::Chain::ABA:: );
        return '' if any { blessed $_ } ( $center, $neighbours[1] );
        return '' if $neighbours[0]->inner_element eq element( $center );
        return '' if $neighbours[0]->outer_element eq element( $neighbours[1] );
    } else {
        # No ABA chain yet
        return '' if any { blessed $_ } @neighbours;
        return '' unless element( $neighbours[0] ) eq element( $neighbours[1] );
        my $outer = element( $neighbours[0] );
        my $inner = element( $center );
        return '' unless $elements{$outer}->{seniority} > $elements{$inner}->{seniority};
    }

    return 1;
}

sub anything { 1 }

my @rules = (
    # Guanidine
    [ sub { &is_nongroup_atom && &is_C && &has_H0 }, ( sub { &is_nongroup_atom && &is_N } ) x 3, NO_MORE_VERTICES,
      sub {
            my $guanidine = ChemOnomatopist::Group::Guanidine->new( $_[0], $_[1] );
            $_[0]->add_group( $guanidine );
            $_[0]->delete_vertex( $_[1] );
            for (combinations( [ @_[2..4] ], 2 )) {
                $_[0]->add_edge( @$_ );
            }
        } ],

    # Chalcogen chains
    [ ( sub { &is_nongroup_atom && &is_O  } ) x 3, NO_MORE_VERTICES,
        sub { $_[0]->add_group( ChemOnomatopist::Chain::Chalcogen->new( $_[0], undef, @_[1..3] ) ) } ],
    [ ( sub { &is_nongroup_atom && &is_S  } ) x 3, NO_MORE_VERTICES,
        sub { $_[0]->add_group( ChemOnomatopist::Chain::Chalcogen->new( $_[0], undef, @_[1..3] ) ) } ],
    [ ( sub { &is_nongroup_atom && &is_Se } ) x 3, NO_MORE_VERTICES,
        sub { $_[0]->add_group( ChemOnomatopist::Chain::Chalcogen->new( $_[0], undef, @_[1..3] ) ) } ],
    [ ( sub { &is_nongroup_atom && &is_Te } ) x 3, NO_MORE_VERTICES,
        sub { $_[0]->add_group( ChemOnomatopist::Chain::Chalcogen->new( $_[0], undef, @_[1..3] ) ) } ],

    # Carboxylic acid
    [ sub { &is_nongroup_atom && &is_C },
      sub { &is_hydroxy || &is_hydroperoxide || ( &is_nongroup_atom && &is_O && &charge_minus_one ) },
      sub { &is_ketone || &is_hydrazine || &is_imine },
      \&anything,
      NO_MORE_VERTICES,

      sub { $_[0]->replace( ChemOnomatopist::Group::Carboxyl->new( $_[2], $_[3] ), @_[1..3] ) } ],

    # Nitramide
    [ sub { &is_nongroup_atom && &is_N && &charge_plus_one }, ( sub { &is_nongroup_atom && &is_O && &has_1_neighbour } ) x 2, sub { &is_nongroup_atom && &is_N }, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Nitramide->new, @_[1..4] ) } ],

    # Azide
    [ sub { &is_nongroup_atom && &is_N && &charge_plus_one },
        EDGE { is_double_bond( @_ ) }, sub { &is_nongroup_atom && &is_N && &charge_minus_one && &has_1_neighbour },
        EDGE { is_double_bond( @_ ) }, \&is_N,
      NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Azide->new, @_[1..3] ) } ],

    # Hydrazine and diazene
    [ ( sub { &is_nongroup_atom && &is_N } ) x 2,
        sub { is_double_bond( @_ )
                ? $_[0]->add_group( ChemOnomatopist::Group::Diazene->new( @_[0..2] ) )
                : $_[0]->add_group( ChemOnomatopist::Group::Hydrazine->new( @_[0..2] ) ) } ],

    # Hydrazide
    [ sub { &is_nongroup_atom && &is_C }, \&is_hydrazine, \&is_ketone,
      sub {
            my $hydrazine = first { $_->isa( ChemOnomatopist::Group::Hydrazine:: ) }
                                    $_[0]->groups( $_[2] );
            my @vertices = $hydrazine->vertices;
            @vertices = reverse @vertices if $vertices[0] == $_[1];
            my $hydrazide = ChemOnomatopist::Group::Hydrazide->new( $_[0], $_[3], @vertices );
            $_[0]->delete_vertices( $_[3] );
            $_[0]->add_group( $hydrazide );
            $_[0]->delete_group( $hydrazine );
        } ],

    # Amide
    # CHECKME: Why ketone has to be deleted and not used in replace? It fails somewhy.
    [ sub { &is_nongroup_atom && &is_C }, \&is_amine, \&is_ketone,
      sub { $_[0]->delete_vertices( $_[3] ); $_[0]->replace( ChemOnomatopist::Group::Amide->new( $_[1], $_[3] ), $_[2] ) } ],
    [ sub { &is_sulfinyl || &is_sulfonyl }, \&is_amine,
      sub { $_[0]->replace( ChemOnomatopist::Group::Amide->new( $_[1] ), $_[2] ) } ],

    # Aldehyde
    [ sub { &is_nongroup_atom && &is_C && &has_H1 }, \&is_ketone,
      sub { $_[0]->replace( ChemOnomatopist::Group::Aldehyde->new( $_[2] ), @_[1..2] ) } ],

    # Aldehydes attached to carbon in cyclic system or a heteroatom
    [ \&is_aldehyde, sub { &is_circular || ( &is_nongroup_atom && &is_heteroatom ) }, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Carbaldehyde->new( $_[1] ), $_[1] ) } ],

    # Acid halide
    [ sub { &is_sulfinyl || &is_sulfonyl }, \&is_Br_Cl_F_I,
      sub { $_[0]->replace( ChemOnomatopist::Group::AcidHalide->new( $_[1], element( $_[2] ) ), @_[1..2] ) } ],

    # Acyl halide
    [ sub { &is_nongroup_atom && &is_C }, sub { &is_nongroup_atom && &is_Br_Cl_F_I }, sub { &is_ketone && &is_O }, \&is_C, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::AcylHalide->new( $_[2] ), @_[1..3] ) } ],
    [ sub { &is_nongroup_atom && &is_C }, sub { &is_cyanide || &is_isocyanide || &is_isocyanate }, sub { &is_ketone && &is_O }, \&is_C, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::AcylHalide->new( $_[2] ), @_[1..3] ) } ],

    # a(ba)n chain
    [ sub { &is_nongroup_atom && &is_heteroatom && &looks_like_ABA_chain }, ( sub { &is_nongroup_atom && &is_heteroatom } ) x 2, NO_MORE_VERTICES,
      sub { $_[0]->add_group( ChemOnomatopist::Chain::ABA->new( $_[0], $_[2], $_[1], $_[3] ) ) } ],
    [ sub { &is_nongroup_atom && &is_heteroatom && &looks_like_ABA_chain }, \&is_ABA_chain, sub { &is_nongroup_atom && &is_heteroatom }, NO_MORE_VERTICES,
      sub { for ($_[0]->groups( $_[2] )) { $_->add( $_[1] ); $_->add( $_[3] ) } } ],
    [ sub { &is_nongroup_atom && &is_heteroatom && &looks_like_ABA_chain }, ( \&is_ABA_chain ) x 2, NO_MORE_VERTICES,
      sub {
            my( $target ) = $_[0]->groups( $_[2] );
            my( $source ) = $_[0]->groups( $_[3] );
            $target->add( $_[1] );
            $target->add( $source );
            $_[0]->delete_group( $source );
          } ],

    # O-based groups
    [ 'hydroxy', sub { &is_nongroup_atom && &is_O && ( &has_H1 || &charge_minus_one ) }, \&anything, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Hydroxy->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_O },
        EDGE { is_double_bond( @_ ) }, \&anything,
      NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Ketone->new( element( $_[1] ) ), $_[1] ) } ],

    # Ester
    [ sub { &is_nongroup_atom && &is_C }, sub { &is_ketone && &is_O }, sub { &is_nongroup_atom && &is_O && &no_charge }, \&is_C, NO_MORE_VERTICES,
      sub { $_[0]->delete_vertex( $_[2] );
            $_[0]->add_group( ChemOnomatopist::Group::Ester->new( @_[0..3] ) ) } ],

    # Ether
    [ sub { &is_nongroup_atom && &is_O }, ( \&is_C ) x 2, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Ether->new, $_[1] ) } ],

    # Hydroxy groups and their chalcogen analogues
    [ sub { &is_nongroup_atom && &is_O_S_Se_Te && ( &has_H1 || &charge_minus_one ) }, \&is_C_N_O_S_Se_Te, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Hydroxy->new( $_[1] ), $_[1] ) } ],
    # Ketones and their chalcogen analogues
    [ sub { &is_nongroup_atom && &is_O_S_Se_Te && all { is_double_bond( @_, $_ ) } $_[0]->neighbours( $_[1] ) },
        EDGE { is_double_bond( @_ ) }, \&anything,
      NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Ketone->new( element( $_[1] ) ), $_[1] ) } ],

    # Urea
    [ sub { &is_nongroup_atom && &is_C }, \&is_ketone, ( sub { &is_nongroup_atom && &is_N } ) x 2, NO_MORE_VERTICES,
      sub { $_[0]->add_group( ChemOnomatopist::Group::Urea->new( @_ ) ) } ],

    # Isocyanide
    [ sub { &is_nongroup_atom && &is_C && &has_H0 && &charge_minus_one }, sub { &is_nongroup_atom && &is_N && &has_H0 && &charge_plus_one }, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Isocyanide->new, @_[1..2] ) } ],

    # Isocyanate
    [ sub { &is_nongroup_atom && &is_C && &has_H0 && &no_charge }, \&is_ketone, \&is_N, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Isocyanate->new( $_[2]->element ), @_[1..3] ) } ],

    # N-based groups
    [ 'cyanide without N-substitutions',
      sub { &is_nongroup_atom && &is_N },
        EDGE { is_triple_bond( @_ ) }, sub { &is_nongroup_atom && &is_C },
        NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Cyanide->new, @_[1..2] ) } ],
    [ 'cyanide with N-substitutions',
      sub { &is_nongroup_atom && &is_N },
        EDGE { is_triple_bond( @_ ) }, sub { &is_nongroup_atom && &is_C },
      sub { die "cannot name N-substituted cyano groups yet\n" } ],
    # TODO: Recognise aminium cations
    # [ sub { &is_nongroup_atom && &is_N && &has_H0 && &charge_plus_one }, ( \&anything ) x 4, NO_MORE_VERTICES,
    #   sub { $_[0]->replace( ChemOnomatopist::Group::Amine->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_N && &no_charge }, ( \&anything ) x 3, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Amine->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_N && &has_H1 }, ( \&anything ) x 2, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Amine->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_N && ( &has_H2 || ( &has_H1 && &charge_minus_one ) ) }, \&anything, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Amine->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_N && &has_H3 }, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Amine->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_N && &has_H0 },
        EDGE { is_double_bond( @_ ) && $_[0]->degree( $_[2] ) + ($_[2]->{hcount} ? $_[2]->{hcount} : 0) == 3 }, \&is_C,
        \&anything,
      NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Imine->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_N && ( &has_H1 || ( &has_H0 && &charge_minus_one ) ) },
        EDGE { is_double_bond( @_ ) && $_[0]->degree( $_[2] ) + ($_[2]->{hcount} ? $_[2]->{hcount} : 0) == 3 }, \&is_C,
      NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Imine->new( $_[1] ), $_[1] ) } ],
    [ sub { &is_nongroup_atom && &is_N && &charge_plus_one }, \&is_ketone, sub { &is_hydroxy && &is_O && &charge_minus_one }, \&anything,
      sub { $_[0]->replace( ChemOnomatopist::Group::Nitro->new, @_[1..3] ) } ],
    [ sub { &is_nongroup_atom && &is_N && &has_H0 },
        EDGE { is_double_bond( @_) && $_[0]->degree( $_[2] ) + ($_[2]->{hcount} ? $_[2]->{hcount} : 0) == 3 }, sub { &is_nongroup_atom && &is_S_Se_Te },
        \&anything,
      NO_MORE_VERTICES,
      sub { $_[0]->add_group( ChemOnomatopist::Chain::Sulfimide->new( @_[0..2] ) ) } ],

    # Esters of nitric acid and nitrous acid
    [ sub { &is_nongroup_atom && &is_N && &charge_plus_one }, \&is_ketone, sub { &is_nongroup_atom && &is_O && &charge_minus_one }, \&is_O_S_Se_Te, NO_MORE_VERTICES,
      sub { die "cannot handle nitric/nitrous acid esters yet\n" } ],
    [ sub { &is_nongroup_atom && &is_N  }, \&is_ketone, \&is_O_S_Se_Te, NO_MORE_VERTICES,
      sub { die "cannot handle nitric/nitrous acid esters yet\n" } ],

    # Carbonitrile, special case of cyanide
    [ 'cyanide attached to circular compound',
      \&is_cyanide, \&is_circular, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Carbonitrile->new, $_[1] ) } ],
    [ 'cyanide attached to heteroatom', \&is_cyanide, \&is_heteroatom, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Carbonitrile->new, $_[1] ) } ],

    # Amidines (BBv3 P-66.4.1)
    [ sub { &is_nongroup_atom && &is_C }, sub { &is_amine && &is_nongroup }, sub { &is_imine && &is_nongroup },
      sub { $_[0]->add_group( ChemOnomatopist::Group::Amidine->new( @_[0..3] ) ) } ],
    [ sub { &is_nongroup_atom && &is_S_Se_Te }, sub { &is_amine && &is_nongroup }, sub { &is_nongroup_atom && &is_N && &has_H1 }, \&anything, NO_MORE_VERTICES,
      sub { $_[0]->add_group( ChemOnomatopist::Group::Amidine->new( @_[0..3] ) ) } ],
    [ sub { &is_nongroup_atom && &is_S_Se_Te }, ( sub { &is_amine && &is_nongroup } ) x 2, sub { &is_nongroup_atom && &is_N && &has_H1 }, \&anything, NO_MORE_VERTICES,
      sub { $_[0]->add_group( ChemOnomatopist::Group::Amidine->new( @_[0..4] ) ) } ],

    # Nitroso and its analogues
    [ sub { &is_nongroup_atom && &is_Br_Cl_F_I_N }, \&is_ketone, sub { &is_C || &is_N }, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Nitroso->new( element( $_[1] ) ), @_[1..2] ) } ],

    # XO3
    [ sub { &is_nongroup_atom && &is_Br_Cl_F_I }, ( sub { &is_ketone && &is_O } ) x 3,
      sub { $_[0]->replace( ChemOnomatopist::Group::XO3->new( element( $_[1] ) ), @_[1..4] ) } ],

    # Peroxide
    [ sub { &is_nongroup_atom && &is_O }, sub { ( &is_nongroup_atom && &is_O ) || ( &is_hydroxy && &charge_minus_one ) }, \&is_C, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Peroxide->new( @_ ), @_[1..2] ) } ],

    # Hydroperoxide
    [ sub { &is_nongroup_atom && &is_O_S_Se_Te }, \&is_hydroxy, \&anything, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Hydroperoxide->new( $_[1], $_[2] ), @_[1..2] ) } ],
    [ sub { &is_nongroup_atom && &is_O_S_Se_Te }, sub { &is_nongroup_atom && &is_O && &charge_minus_one }, \&anything, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Hydroperoxide->new( $_[1], $_[2] ), @_[1..2] ) } ],

    # S-based groups
    [ sub { &is_nongroup_atom && &is_S_Se_Te }, sub { &is_ketone || ( &is_nongroup_atom && &is_N && &has_H1 ) || &is_hydrazine }, sub { &is_hydroxy || &is_hydroperoxide }, \&is_C, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::SulfinicAcid->new( element( $_[1] ), @_[2..3] ), @_[1..3] ) } ],
    [ sub { &is_nongroup_atom && &is_S_Se_Te }, ( sub { &is_ketone || ( &is_nongroup_atom && &is_N && &has_H1 ) || &is_hydrazine } ) x 2, sub { &is_hydroxy || &is_hydroperoxide }, \&is_C, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::SulfonicAcid->new( element( $_[1] ), @_[2..4] ), @_[1..4] ) } ],

    # Sulfoxide group and its analogues
    [ sub { &is_nongroup_atom && &is_S_Se_Te }, \&is_ketone, ( \&anything ) x 2, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Sulfinyl->new( element( $_[1] ), $_[2] ), @_[1..2] ) } ],
    [ sub { &is_nongroup_atom && &is_S_Se_Te }, ( \&is_ketone ) x 2, ( \&anything ) x 2, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Sulfonyl->new( element( $_[1] ), $_[2], $_[3] ), @_[1..3] ) } ],

    # Noncarbon oxoacids
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, ( sub { &is_hydroxy && &is_O } ) x 4,
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..5] ), @_[1..5] ) } ],
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, ( sub { &is_hydroxy && &is_O } ) x 2, ( sub { &is_ketone && &is_O } ) x 2,
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..5] ), @_[1..5] ) } ],
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, ( sub { &is_hydroxy && &is_O } ) x 3, sub { &is_ketone && &is_O },
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..5] ), @_[1..5] ) } ],
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, ( sub { &is_hydroxy && &is_O } ) x 3,
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..4] ), @_[1..4] ) } ],
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, ( sub { &is_hydroxy && &is_O } ) x 2, sub { &is_ketone && &is_O },
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..4] ), @_[1..4] ) } ],
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, ( sub { &is_ketone && &is_O }  ) x 2, sub { &is_hydroxy && &is_O },
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..4] ), @_[1..4] ) } ],
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, ( sub { &is_hydroxy && &is_O } ) x 2,
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..3] ), @_[1..3] ) } ],
    [ sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te },   sub { &is_hydroxy && &is_O }, sub { &is_ketone && &is_O },
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..3] ), @_[1..3] ) } ],
    [ 'noncarbon oxoacid', sub { &is_nongroup_atom && &is_As_N_B_P_Se_Si_Sb_S_Te }, sub { &is_hydroxy && &is_O },
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..2] ), @_[1..2] ) } ],
    [ \&is_XO3, sub { &is_hydroxy && &is_O },
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..2] ), @_[1..2] ) } ],
    [ sub { &is_sulfinyl || &is_sulfonyl }, ( sub { &is_hydroxy && &is_O } ) x 2,
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..3] ), @_[1..3] ) } ],
    [ \&is_nitro, sub { &is_hydroxy && &is_O },
      sub { $_[0]->replace( ChemOnomatopist::Group::NoncarbonOxoacid->new( @_[1..2] ), @_[1..2] ) } ],

    # Sulfinamides and sulfonamides
    [ \&is_amide, \&is_sulfinyl, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Sulfinamide->new( $_[2]->element, $_[2]->{ketone} ), @_[1..2] ) } ],
    [ \&is_amide, \&is_sulfonyl, NO_MORE_VERTICES,
      sub { $_[0]->replace( ChemOnomatopist::Group::Sulfonamide->new( $_[2]->element, @{$_[2]->{ketones}} ), @_[1..2] ) } ],

    # Detecting amides attached to cyclic chains
    [ sub { &is_nongroup_atom && &is_C && 1 == grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Amide:: ) && $_->{parent} == $_[1] } $_[0]->neighbours( $_[1] ) }, \&is_amide, \&is_monocycle, NO_MORE_VERTICES,
      sub {
            my( $cycle ) = $_[0]->groups( $_[3] );
            $_[0]->delete_group( $cycle );
            $_[0]->add_group( ChemOnomatopist::Chain::Carboxamide->new( $_[0], $_[2], $_[1], $cycle ) ) } ],
);

# Old unused rules
my @rules_old = (
    [ \&is_C, \&is_benzene, \&is_ketone, \&is_N, NO_MORE_VERTICES,
      sub { $_[0]->replace( { type => 'benzamide' }, @_[1..4] ) } ],

    [ \&is_benzene, \&is_hydroxy, sub { $_[0]->replace( { type => 'phenol' }, @_[1..2] ) } ],
);

sub is_mainchain() { any { $_->is_main } $_[0]->groups( $_[1] ) }

sub is_carboxamide()  { any { $_->isa( ChemOnomatopist::Chain::Carboxamide:: ) } $_[0]->groups( $_[1] ) }
sub is_purine()       { any { $_->isa( ChemOnomatopist::Chain::Bicycle::Purine:: ) } $_[0]->groups( $_[1] ) }

sub most_senior_group() { my @most_senior_groups = map { $_->most_senior_groups } $_[0]->groups( $_[1] ); return shift @most_senior_groups }
sub number_of_most_senior_groups() { scalar map { $_->most_senior_groups } $_[0]->groups( $_[1] ) }

our @mainchain_rules = (

    # Amide chains
    [ sub { &is_mainchain && !&is_carboxamide && &most_senior_group && &most_senior_group->isa( ChemOnomatopist::Group::Amide:: ) && &number_of_most_senior_groups == 1 },
      sub { &is_amide && &is_nongroup },
      sub {
            my( $chain ) = $_[0]->groups( $_[1] );
            $_[0]->delete_group( $chain );
            my $amide = ChemOnomatopist::Chain::Amide->new( $_[0], $chain, $_[2] );
            $amide->{is_main} = 1;
            $_[0]->add_group( $amide );
          } ],

    # Amine chains
    [ sub { &is_mainchain && !&is_purine && &most_senior_group && &most_senior_group->isa( ChemOnomatopist::Group::Amine:: ) && &number_of_most_senior_groups == 1 },
      sub { &is_amine && &is_nongroup },
      sub {
            my( $chain ) = $_[0]->groups( $_[1] );
            $_[0]->delete_group( $chain );
            my $amine = ChemOnomatopist::Chain::Amine->new( $_[0], $chain, $_[2] );
            $amine->{is_main} = 1;
            $_[0]->add_group( $amine );
          } ],

    # Imine chains
    [ sub { &is_mainchain && &most_senior_group && &most_senior_group->isa( ChemOnomatopist::Group::Imine:: ) && &number_of_most_senior_groups == 1 },
      sub { &is_imine && &is_nongroup },
      sub {
            my( $chain ) = $_[0]->groups( $_[1] );
            $_[0]->delete_group( $chain );
            my $imine = ChemOnomatopist::Chain::Imine->new( $_[0], $chain, $_[2] );
            $imine->{is_main} = 1;
            $_[0]->add_group( $imine );
          } ],

    # Unclaimed amine group attached to other chains
    [ sub { !&is_mainchain && &is_circular && &most_senior_group && &most_senior_group->isa( ChemOnomatopist::Group::Amine:: ) && &number_of_most_senior_groups == 1 },
      sub { &is_amine && &is_nongroup },
      sub {
            my( $chain ) = $_[0]->groups( $_[1] );
            $_[0]->delete_group( $chain );
            my $amine = ChemOnomatopist::Chain::Amine->new( $_[0], $chain, $_[2] );
            $_[0]->add_group( $amine );
          } ],

    # Unpack sidechain amides
    [ sub { !&is_mainchain && &is_amide && !&is_carboxamide },
      sub {
            $_[0]->replace( ChemOnomatopist::Group::Amine->new, $_[1] );
            $_[0]->set_edge_attribute( $_[1]->{parent}, $_[1]->{ketone}, 'bond', '=' ) if $_[1]->{ketone};
          } ],

);

sub parse_molecular_graph($)
{
    my( $graph ) = @_;
    return parse_graph( $graph, @rules );
}

1;
