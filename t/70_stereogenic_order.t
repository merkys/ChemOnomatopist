#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::MolecularGraph;
use Chemistry::OpenSMILES qw( is_chiral );
use Chemistry::OpenSMILES::Parser;
use List::Util qw( first );
use Test::More;

my @cases = (
    { smiles => 'ClCC[C@H](CC(C)O)O',    order => '8,2,4,13' }, # From BBv3 P-92.1.4.1

    # From BBv3 P-92.1.4.2
    { smiles => 'C[C@H](C=C)O',          order => '4,2,0,8'  },
    { smiles => 'O[C@H](C=O)C',          order => '0,2,4,6'  },

    # From BBv3 P-92.1.4.3
    { smiles => 'C1(CCC[C@]12C=CCC2)=O', order => '0,5,8,3'  },
    { smiles => 'C1(CC1)[C@@H](C(C)C)O', order => '7,0,4,13' },

    # Derived from COD entry 1516725
    { smiles => 'C1N(C=C(CN1))Cc1ccc(cc1)[C@H]1[C@H](c2ccncc2)[C@H](c2ccc(CN3CNCC=C3)cc2)[C@H]1c1ccncc1', order => [ '14,35,10,55', '35,14,10,55' ] },
    { smiles => 'O=C1N(C=C(C(=O)N1C)C)Cc1ccc(cc1)[C@H]1[C@H](c2ccncc2)[C@H](c2ccc(CN3C(=O)N(C)C(=O)C(C)=C3)cc2)[C@H]1c1ccncc1', order => [ '18,43,14,63', '43,18,14,63' ] },

    # Synthetic test case for a center in a ring
    { smiles => '[C@H]1(Cl)OCCCCC1',     order => '1,2,7,8' },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

my $parser = Chemistry::OpenSMILES::Parser->new;

for my $case (@cases) {
    # Manual bless() as new() would perform the chirality ordering on its own
    my( $moiety ) = map { bless $_, ChemOnomatopist::MolecularGraph:: }
                        $parser->parse( $case->{smiles} );
    my $atom = first { is_chiral $_ } $moiety->vertices;
    my $digraph = $moiety->hierarchical_digraph( $atom );
    my @order = sort { ChemOnomatopist::MolecularGraph::_order_chiral_center_neighbours( $digraph, $a, $b ) }
                     $moiety->neighbours( $atom );
    my $ok;
    if( ref $case->{order} ) {
        my $regex = '^(' . join( '|', @{$case->{order}} ) . ')$';
        $ok = like join( ',', map { $_->{number} } @order ), qr/$regex/;
    } else {
        $ok = is join( ',', map { $_->{number} } @order ), $case->{order};
    }
    if( $case->{AUTHOR} && $ok ) {
        diag 'test supposed to fail with AUTHOR_TESTING' .
             ( $case->{AUTHOR} !~ /^1$/ ? ': ' . $case->{AUTHOR} : '' );
    }
}
