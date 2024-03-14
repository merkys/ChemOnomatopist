#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-65.6.3.2.1
    { smiles => 'C(C)(=O)OCC', iupac => 'ethyl acetate', AUTHOR => 1 },
    { smiles => 'C(CCC(=O)OC)(=O)OCC', iupac => 'ethyl methyl butanedioate', AUTHOR => 1 },
    { smiles => 'C1(CCCCC1)C(=O)OC', iupac => 'methyl cyclohexanecarboxylate', AUTHOR => 1 },
    { smiles => 'C(C)C1=CC=C(C=C1)S(=O)(=O)OC', iupac => 'methyl 4-ethylbenzene-1-sulfonate', AUTHOR => 1 },

    # From BBv3 P-65.6.3.2.2
    { smiles => 'C(C)(=O)OCC(COC(C)=O)OC(C)=O', iupac => 'propane-1,2,3-triyl triacetate', AUTHOR => 1 },
    { smiles => 'C(CC(=O)OC1=CC=C(C=C1)OC(CC(=O)OC)=O)(=O)OC', iupac => 'dimethyl 1,4-phenylene dipropanedioate', AUTHOR => 1 },
    { smiles => 'C(CC(=O)OC1=CC=C(C=C1)OC(CC(=O)OC)=O)(=O)OCC', iupac => 'ethyl methyl 1,4-phenylene dipropanedioate', AUTHOR => 1 },

    # From BBv3 P-65.6.3.2.3
    { smiles => '[Br-].C(C)OC(CC[N+](C)(C)C)=O', iupac => '3-ethoxy-N,N,N-trimethyl-3-oxopropan-1-aminium bromide', AUTHOR => 1 },
    { smiles => 'C(C1=CC=CC=C1)(=O)OCCC(=O)O', iupac => '3-(benzoyloxy)propanoic acid', AUTHOR => 1 },
    { smiles => 'C(C)(=O)OCCS(=O)(=O)O', iupac => '2-(acetyloxy)ethane-1-sulfonic acid', AUTHOR => 1 },
    { smiles => 'O(C1=CC=CC=C1)S(=S)C1=CC=C(C2=CC=CC=C12)C(=O)OC', iupac => 'methyl 4-(phenoxysulfinothioyl)naphthalene-1-carboxylate', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)SS(=O)(=O)C1=CC=C(C2=CC=CC=C12)C(=O)OC', iupac => 'methyl 4-[(phenylsulfanyl)sulfonyl]naphthalene-1-carboxylate', AUTHOR => 1 },
    { smiles => 'C(C)OC(=O)OC(C(=O)OCC)C(C(C)(C)C)=O', iupac => 'ethyl 2-[(ethoxycarbonyl)oxy]-4,4-dimethyl-3-oxopentanoate', AUTHOR => 1 },
    { smiles => 'N1=CC(=CC=C1)C(=O)OCCC(=O)O', iupac => '3-[(pyridine-3-carbonyl)oxy]propanoic acid', AUTHOR => 1 },
    { smiles => 'N1=C(C=CC2=CC=CC=C12)C(=O)OCC(=O)O', iupac => '[(quinoline-2-carbonyl)oxy]acetic acid', AUTHOR => 1 },

    { smiles => 'CCCCCCCC(=O)OC(C)(C)C', iupac => 'tert-butyl octanoate', AUTHOR => 1 }, # BBv2 P-65.6.3.3.1
    { smiles => 'CCOC(=O)CC(=O)OC', iupac => 'ethyl methyl propanedioate', AUTHOR => 1 }, # BBv2 P-65.6.3.3.2.1
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
