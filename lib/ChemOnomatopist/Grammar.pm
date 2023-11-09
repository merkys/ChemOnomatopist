package ChemOnomatopist::Grammar;

# ABSTRACT: Grammar for chemical graphs
# VERSION

use strict;
use warnings;

use ChemOnomatopist::Chain;
use ChemOnomatopist::Chain::Ether;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Ketone;
use ChemOnomatopist::Util::Graph qw(
    graph_replace
);
use Chemistry::OpenSMILES::Parser;
use Graph::Grammar;
use List::Util qw( sum );
use Scalar::Util qw( blessed );

use parent Exporter::;
our @EXPORT_OK = qw(
    parse_molecular_graph
);

sub is_C { return exists $_[1]->{symbol} && ucfirst( $_[1]->{symbol} ) eq 'C' }
sub is_N { return exists $_[1]->{symbol} && ucfirst( $_[1]->{symbol} ) eq 'N' }
sub is_O { return exists $_[1]->{symbol} && ucfirst( $_[1]->{symbol} ) eq 'O' }
sub is_S { return exists $_[1]->{symbol} && ucfirst( $_[1]->{symbol} ) eq 'S' }

sub is_CH2 { return is_C( @_ ) && $_[1]->{hcount} == 2 }
sub is_CH3 { return is_C( @_ ) && $_[1]->{hcount} == 3 }
sub is_NH2 { return is_N( @_ ) && $_[1]->{hcount} == 2 }
sub is_NH3 { return is_N( @_ ) && $_[1]->{hcount} == 3 }
sub is_OH  { return is_O( @_ ) && $_[1]->{hcount} == 1 }
sub is_SH  { return is_S( @_ ) && $_[1]->{hcount} == 1 }

sub is_any_chain { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Chain:: ) }

sub is_chain { return blessed $_[1] && blessed $_[1] eq ChemOnomatopist::Chain:: }

sub is_C_chain { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Chain:: ) }
sub is_C_chain_carboxyl { return exists $_[1]->{type} && $_[1]->{type} eq 'C_chain_carboxyl' }
sub is_carboxyl { return exists $_[1]->{type} && $_[1]->{type} eq 'carboxyl' }
sub is_headless_C_chain { return exists $_[1]->{type} && $_[1]->{type} eq 'headless_C_chain' }
sub is_hydroxy { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Hydroxy:: ) }
sub is_ketone  { return blessed $_[1] && $_[1]->isa( ChemOnomatopist::Group::Ketone:: ) }

sub is_cyano { return exists $_[1]->{type} && $_[1]->{type} eq 'cyano' }

sub is_cycle { return exists $_[1]->{type} && $_[1]->{type} eq 'cycle' }
sub is_benzene { return exists $_[1]->{type} && $_[1]->{type} eq 'benzene' }

sub anything { return 1 }

my @rules = (
    # O-based groups
    [ \&is_OH, \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Hydroxy->new( ChemOnomatopist::element( $_[1] ) ), $_[1] ) } ],
    [ \&is_O,  \&anything, NO_MORE_VERTICES,
      sub { graph_replace( $_[0], ChemOnomatopist::Group::Ketone->new( ChemOnomatopist::element( $_[1] ) ), $_[1] ) } ],
    [ \&is_O,  \&is_chain, \&anything, NO_MORE_VERTICES,
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
    [ \&is_C, \&is_ketone, \&is_hydroxy, \&anything, NO_MORE_VERTICES, # Carboxyl group
      sub { graph_replace( $_[0], { type => 'carboxyl' }, @_[1..3] ) } ],
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

    [ \&is_C, \&is_NH2, \&is_ketone, \&anything, NO_MORE_VERTICES, { type => 'amide' } ],

    [ \&is_N, \&is_C, NO_MORE_VERTICES, { type => 'cyano' } ],
    [ \&is_NH2, \&is_NH2, NO_MORE_VERTICES, sub { graph_replace( $_[0], { type => 'hydrazine' }, @_[1..2] ) } ],

    [ \&is_SH, { type => 'sulfanyl' } ],
    [ \&is_S, \&is_ketone, ( \&anything ) x 2, NO_MORE_VERTICES, sub { graph_replace( $_[0], { type => 'sulfoxide' }, @_[1..2] ) } ],

    [ \&is_benzene, \&is_hydroxy, sub { graph_replace( $_[0], { type => 'phenol' }, @_[1..2] ) } ],

    [ \&is_C, \&is_cyano, \&is_cycle, NO_MORE_VERTICES, sub { graph_replace( $_[0], { type => 'carbonitrile' }, @_[1..2] ) } ],
);

sub parse_molecular_graph($)
{
    my( $graph ) = @_;
    return parse_graph( $graph, @rules );
}

1;
