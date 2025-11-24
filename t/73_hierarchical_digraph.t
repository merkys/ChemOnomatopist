#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::MolecularGraph;
use Chemistry::OpenSMILES::Parser;
use Chemistry::OpenSMILES::Writer qw( write_SMILES );
use List::Util qw( first );
use Test::More;

my @cases = (
    # From BBv3 P-92.1.4.3
    { smiles => 'O=C1CCCC12CCC=C2', atom => 5, digraph => 'CCCCC([O])(OC)C(C(C)C(C)CCC)(CCC(C)C(C)C)CCCC([O])(OC)C', AUTHOR => 1 },

    # (1R)-1-cyclopropyl-2-methylpropan-1-ol from BBv3 P-92.1.4.3
    { smiles => 'C1(CC1)[C@@H](C(C)C)O', atom => 3, digraph => 'CC(C)C(O)C(CC[C])CC[C]' },
);

eval 'use Graph::Nauty qw( are_isomorphic )';
my $has_Graph_Nauty = !$@;

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
@cases = () unless $has_Graph_Nauty;
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

my $parser = Chemistry::OpenSMILES::Parser->new;

for my $case (@cases) {
    my( $moiety ) = map { ChemOnomatopist::MolecularGraph->new( $_ ) }
                        $parser->parse( $case->{smiles} );
    my $atom = first { $_->{number} == $case->{atom} } $moiety->vertices;
    my $digraph_got = $moiety->hierarchical_digraph( $atom );
    my( $digraph_exp ) = $parser->parse( $case->{digraph} );
    ok are_isomorphic( $digraph_got, $digraph_exp, sub { $_[0]->{original} ? $_[0]->{original}{symbol} : $_[0]->{symbol} } );
}
