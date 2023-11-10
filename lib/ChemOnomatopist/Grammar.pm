package ChemOnomatopist::Grammar;

# ABSTRACT: Grammar for chemical graphs
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Chain;
use ChemOnomatopist::Chain::Carboxamide;
use ChemOnomatopist::Chain::Circular;
use ChemOnomatopist::Chain::Ether;
use ChemOnomatopist::Group::AcylHalide;
use ChemOnomatopist::Group::Aldehyde;
use ChemOnomatopist::Group::Amide;
use ChemOnomatopist::Group::Amine;
use ChemOnomatopist::Group::Carboxyl;
use ChemOnomatopist::Group::Cyanide;
use ChemOnomatopist::Group::Hydroperoxide;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Imino;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Group::Nitro;
use ChemOnomatopist::Group::Nitroso;
use ChemOnomatopist::Group::SulfinicAcid;
use ChemOnomatopist::Group::Sulfinyl;
use ChemOnomatopist::Group::SulfonicAcid;
use ChemOnomatopist::Group::Sulfonyl;
use ChemOnomatopist::Group::XO3;
use ChemOnomatopist::Util::Graph qw(
    graph_replace
);
use Chemistry::OpenSMILES qw(
    is_double_bond
);
use Graph::Grammar;
use List::Util qw( all first sum );
use Scalar::Util qw( blessed );

use parent Exporter::;
our @EXPORT_OK = qw(
    parse_molecular_graph
);

sub is_nongroup_atom { return !blessed $_[1] && !$_[0]->groups( $_[1] ) && exists $_[1]->{symbol} }

sub is_C { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) eq 'C' }
sub is_N { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) eq 'N' }
sub is_O { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) eq 'O' }
sub is_S { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) eq 'S' }

sub is_Br_Cl_F_I     { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) =~ /^(Br|Cl|F|I)$/ }
sub is_Br_Cl_F_I_N   { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) =~ /^(Br|Cl|F|I|N)$/ }
sub is_B_Cl_F_I      { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) =~ /^(B|Cl|F|I)$/ }
sub is_S_Se_Te       { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) =~ /^(S|Se|Te)$/ }
sub is_O_S_Se_Te     { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) =~ /^(O|S|Se|Te)$/ }
sub is_C_N_O_S_Se_Te { return &is_nongroup_atom && ucfirst( $_[1]->{symbol} ) =~ /^(C|N|O|S|Se|Te)$/ }

sub is_CH1 { return &is_C &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 1 }
sub is_CH2 { return &is_C &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 2 }
sub is_CH3 { return &is_C &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 3 }
sub is_NH0 { return &is_N && !$_[1]->{hcount} }
sub is_NH1 { return &is_N &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 1 }
sub is_NH2 { return &is_N &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 2 }
sub is_NH3 { return &is_N &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 3 }
sub is_OH  { return &is_O &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 1 }
sub is_SH  { return &is_S &&  exists $_[1]->{hcount} && $_[1]->{hcount} == 1 }

sub is_Np1 { return &is_N &&  exists $_[1]->{charge} && $_[1]->{charge} ==  1 }
sub is_On1 { return &is_O &&  exists $_[1]->{charge} && $_[1]->{charge} == -1 }

sub is_any_chain { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Chain:: ) }

sub is_chain { return blessed $_[1] && blessed $_[1] eq ChemOnomatopist::Chain:: }

sub is_C_chain { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Chain:: ) }
sub is_C_chain_carboxyl { return exists $_[1]->{type} && $_[1]->{type} eq 'C_chain_carboxyl' }
sub is_carboxyl { return exists $_[1]->{type} && $_[1]->{type} eq 'carboxyl' }
sub is_headless_C_chain { return exists $_[1]->{type} && $_[1]->{type} eq 'headless_C_chain' }

sub is_amide   { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Amide:: ) }
sub is_amine   { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Amine:: ) }
sub is_hydroxy { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Hydroxy:: ) }
sub is_ketone  { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Ketone:: ) }

sub is_cyanide { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Cyanide:: ) }

sub is_circular  { return blessed $_[1] && $_->isa( ChemOnomatopist::Chain::Circular:: ) }
sub is_benzene   { return &is_monocycle && $_[1]->is_benzene }
sub is_monocycle { return &is_circular && $_->isa( ChemOnomatopist::Chain::Monocycle:: ) }

sub anything { return 1 }

my @rules = (
    # O-based groups
    [ \&is_O, \&is_chain, \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Chain::Ether->new( $_[0], undef, $_[2]->vertices, $_[1] ), @_[1..2] ) } ],

    # Rules to detect alkanes of any length
    [ \&is_CH3, \&anything, NO_MORE_VERTICES,
      sub { ChemOnomatopist::Chain->new( $_[0], undef, $_[1] ) } ],
    [ \&is_CH2, \&anything, NO_MORE_VERTICES, # Terminal alkene
      sub { ChemOnomatopist::Chain->new( $_[0], undef, $_[1] ) } ],
    [ \&is_CH2, \&is_any_chain, # Add carbon to any chain
      sub { graph_replace( $_[0], $_[2], $_[1] ); push @{$_[2]->{vertices}}, $_[1] } ],
    [ \&is_CH2, ( \&is_chain ) x 2, NO_MORE_VERTICES, # CHECKME: Why do we need this?
      sub { graph_replace( $_[0], ChemOnomatopist::Chain->new( $_[0], undef, reverse( $_[2]->vertices ), $_[1], $_[3]->vertices ), @_[1..3] ) } ],
    [ ( \&is_chain ) x 2, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Chain->new( $_[0], undef, reverse( $_[1]->vertices ), $_[2]->vertices ), @_[1..2] ) } ],
    [ \&is_chain, \&is_headless_C_chain,
      sub { graph_replace( $_[0], ChemOnomatopist::Chain->new( $_[0], undef, reverse( $_[1]->vertices ), $_[2]->vertices ), @_[1..2] ) } ],

    # Handling of headless chains
    [ \&is_CH2, ( \&anything ) x 2, NO_MORE_VERTICES, # Start a headless C chain
      sub { graph_replace( $_[0], ChemOnomatopist::Chain->new( $_[0], undef, $_[1] ), @_[1..3] ) } ],
    [ \&is_headless_C_chain, \&is_headless_C_chain, # Join two headless C chains
      sub { graph_replace( $_[0], ChemOnomatopist::Chain->new( $_[0], undef, reverse( $_[1]->vertices ), $_[2]->vertices ), @_[1..2] ) } ],

    # Carboxyl group and chains it is attached to
    [ \&is_carboxyl, \&is_headless_C_chain, NO_MORE_VERTICES,
      sub { $_[2]->{type} = 'C_chain_carboxyl'; $_[2]->{length}++; $_[0]->delete_vertex( $_[1] ) } ],
    [ \&is_carboxyl, \&is_benzene, NO_MORE_VERTICES,
      sub { $_[2]->{type} = 'benzoic acid'; $_[0]->delete_vertex( $_[1] ) } ],

    [ \&is_headless_C_chain, ( \&is_C_chain_carboxyl ) x 2, NO_MORE_VERTICES,
      sub { $_[1]->{length} = sum map { $_->{length} } @_[1..3]; $_[1]->{type} = 'C_chain_dicarboxyl'; $_[0]->delete_vertices( @_[2..3] ) } ],
    [ \&is_C_chain_carboxyl, \&is_carboxyl, NO_MORE_VERTICES,
      sub { $_[1]->{length} += 1; $_[1]->{type} = 'C_chain_dicarboxyl'; $_[0]->delete_vertices( $_[2] ) } ],

    [ \&is_C, ( \&is_N ) x 3, NO_MORE_VERTICES, # Guanidine
      sub { graph_replace( $_[0], { type => 'guanidine' }, @_[1..3] ) } ],

    [ \&is_C, \&is_benzene, \&is_ketone, \&is_N, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], { type => 'benzamide' }, @_[1..4] ) } ],

    [ \&is_NH2, \&is_NH2, NO_MORE_VERTICES, sub { graph_replace( $_[0], { type => 'hydrazine' }, @_[1..2] ) } ],

    [ \&is_SH, { type => 'sulfanyl' } ],
    [ \&is_S, \&is_ketone, ( \&anything ) x 2, NO_MORE_VERTICES, sub { graph_replace( $_[0], { type => 'sulfoxide' }, @_[1..2] ) } ],

    [ \&is_benzene, \&is_hydroxy, sub { graph_replace( $_[0], { type => 'phenol' }, @_[1..2] ) } ],

    [ \&is_C, \&is_cyanide, \&is_cycle, NO_MORE_VERTICES, sub { graph_replace( $_[0], { type => 'carbonitrile' }, @_[1..2] ) } ],
);

# Conservative rules
my @rules_conservative = (
    # Carboxylic acid
    [ \&is_C, \&is_hydroxy, \&is_ketone, \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Carboxyl->new( $_[4] ), @_[1..3] ) } ],

    # Aldehyde
    [ \&is_CH1, \&is_ketone,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Aldehyde->new( $_[2] ), @_[1..2] ) } ],

    # Amide
    [ \&is_C, \&is_amine, \&is_ketone,
      sub { $_[0]->delete_vertices( $_[3] ); graph_replace( $_[0], ChemOnomatopist::Group::Amide->new( $_[1] ), $_[2] ) } ],

    # Ester
    [ \&is_C, \&is_ketone, \&is_O, \&is_C, NO_MORE_VERTICES,
      sub {
            my $hydroxylic = first { $_ != $_[1] } $_[0]->neighbours( $_[3] );
            my $ester = ChemOnomatopist::Group::Ester->new( $hydroxylic, $_[4] );
            graph_replace( $_[0], $ester, @_[1..3] );
          } ],

    # Acyl halide
    [ \&is_C, \&is_B_Cl_F_I, \&is_ketone, \&is_C, NO_MORE_VERTICES, # FIXME: Ketone must be O
      sub { graph_replace( $_[0], ChemOnomatopist::Group::AcylHalide->new( $_[3] ), @_[1..3] ) } ],

    # O-based groups
    [ \&is_OH, \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Hydroxy->new( ChemOnomatopist::element( $_[1] ) ), $_[1] ) } ],
    [ sub { return is_O( @_ ) && all { is_double_bond( @_, $_ ) } $_[0]->neighbours( $_[1] ) }, \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Ketone->new( ChemOnomatopist::element( $_[1] ) ), $_[1] ) } ],
    [ \&is_O, ( \&is_C ) x 2, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Ether->new, $_[1] ) } ],

    # Hydroxy groups and their chalcogen analogues
    [ \&is_O_S_Se_Te, \&is_C_N_O_S_Se_Te, NO_MORE_VERTICES, # FIXME: Check for a hydrogen
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Hydroxy->new( ChemOnomatopist::element( $_[1] ) ), $_[1] ) } ],
    # Ketones and their chalcogen analogues
    [ \&is_O_S_Se_Te, \&is_C, NO_MORE_VERTICES, # FIXME: Check for double bond
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Ketone->new( ChemOnomatopist::element( $_[1] ) ), $_[1] ) } ],

    # N-based groups
    [ \&is_NH0, \&is_C, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Cyanide->new, @_[1..2] ) } ],
    [ \&is_N,   ( \&anything ) x 3, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Amine->new, $_[1] ) } ],
    [ \&is_NH1, ( \&anything ) x 2, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Amine->new, $_[1] ) } ],
    [ \&is_NH2,   \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Amine->new, $_[1] ) } ],
    [ \&is_NH3, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Amine->new, $_[1] ) } ],
    [ \&is_NH0, \&is_C, \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Imino->new, $_[1] ) } ],
    [ \&is_NH1, \&is_C, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Imino->new, $_[1] ) } ],
    [ \&is_Np1, \&is_ketone, \&is_On1, \&is_C,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Nitro->new, @_[1..3] ) } ],

    # Nitroso and its analogues
    [ \&is_Br_Cl_F_I_N, \&is_ketone, \&is_C, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Nitroso->new( ChemOnomatopist::element( $_[1] ) ), @_[1..2] ) } ],
    # XO3
    [ \&is_Br_Cl_F_I, ( \&is_ketone ) x 3, \&is_C, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::XO3->new( ChemOnomatopist::element( $_[1] ) ), @_[1..4] ) } ],

    # S-based groups
    [ \&is_S, \&is_O, \&is_hydroxy, \&is_C, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::SulfinicAcid->new( $_[4] ), @_[1..3] ) } ],
    [ \&is_S, ( \&is_O ) x 2, \&is_hydroxy, \&is_C, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::SulfonicAcid->new( $_[5] ), @_[1..4] ) } ],

    # Sulfoxide group and its analogues
    [ \&is_S_Se_Te, \&is_ketone, ( \&anything ) x 2, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Sulfinyl->new( ChemOnomatopist::element( $_[1] ) ), @_[1..2] ) } ],
    [ \&is_S_Se_Te, ( \&is_ketone ) x 2, ( \&anything ) x 2, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Sulfonyl->new( ChemOnomatopist::element( $_[1] ) ), @_[1..3] ) } ],

    # Hydroperoxide
    [ \&is_O_S_Se_Te, \&is_hydroxy, \&is_C,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Hydroperoxide->new( $_[1], $_[2] ), @_[1..2] ) } ],

    # Detecting amides attached to cyclic chains
    [ sub { return is_C( @_ ) && 1 == grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Amide:: ) && $_->{parent} == $_[1] } $_[0]->neighbours( $_[1] ) }, \&is_amide, \&is_monocycle, NO_MORE_VERTICES,
      sub { $_[0]->delete_group( $_[3] ); $_[0]->add_group( ChemOnomatopist::Chain::Carboxamide->new( $_[0], $_[2], $_[1], $_[3] ) ) } ],
);

sub parse_molecular_graph($)
{
    my( $graph ) = @_;
    return parse_graph( $graph, @rules_conservative );
}

1;
