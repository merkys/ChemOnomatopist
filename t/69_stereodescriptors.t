#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # BBv3 P-91.3
    { smiles => 'Br[C@H](F)Cl', iupac => '(R)-bromo(chloro)(fluoro)methane' },
    { smiles => 'C\C=C/C', iupac => '(2Z)-but-2-ene', AUTHOR => 1 },
    { smiles => 'C=C[C@@H](CC)O', iupac => '(3R)-pent-1-en-3-ol' },
    { smiles => 'Cl[C@H](CC)C1=CC=CC=C1', iupac => '[(1R)-1-chloropropyl]benzene' },
    { smiles => 'Cl[C@@H](CCC(=O)O[C@@H](C)CC)CC', iupac => '(2S)-butan-2-yl (4R)-4-chlorohexanoate', AUTHOR => 1 },
    { smiles => 'Cl[C@H]([C@H](C(=O)O)O)C', iupac => '(2S,3S)-3-chloro-2-hydroxybutanoic acid' },
    { smiles => 'Cl[C@H](CC\C=C/C)C', iupac => '(2Z,6S)-6-chlorohept-2-ene', AUTHOR => 1 },
    { smiles => 'ClC1=CCCCCCCCCCC1', iupac => '(1Ξ)-1-chlorocyclododec-1-ene', AUTHOR => 1 },
    { smiles => 'C(=C\C)/C(CC=C)\C=C/C', iupac => '(5Z)-4-[(1E)-prop-1-en-1-yl]hepta-1,5-diene', AUTHOR => 1 },

    { smiles => 'C1(CC1)[C@@H](C(C)C)O', iupac => '(1R)-1-cyclopropyl-2-methylpropan-1-ol' }, # BBv3 P-92.1.4.3
    { smiles => 'C1(=CCC=C1)[C@H](C1=C[CH-]C=C1)O', iupac => '3-[(S)-(cyclopenta-1,4-dien-1-yl)(hydroxy)methyl]cyclopenta-2,4-dien-1-ide', AUTHOR => 1 }, # BBv3 P-92.1.4.4

    # BBv3 P-92.2.1.1.1
    { smiles => 'CO[C@@H](C)SC', iupac => '(1R)-1-methoxy-1-(methylsulfanyl)ethane', AUTHOR => 1 },
    { smiles => '[GeH3][C@](SC)(OC)[SiH3]', iupac => '[(R)-germyl(methoxy)(methylsulfanyl)methyl]silane', AUTHOR => 1 },
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac}, $case->{smiles} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    diag 'test supposed to fail with AUTHOR_TESTING' if $case->{AUTHOR} && $ok;
}
