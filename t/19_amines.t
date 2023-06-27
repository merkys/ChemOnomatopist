#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CC(CN)O', iupac => '1-aminopropan-2-ol' },
    { smiles => 'C(CC(=O)N)C(=O)C(=O)O', iupac => '5-amino-2,5-dioxopentanoic acid' },
    { smiles => 'CC(CC(C(=O)O)N)N', iupac => '2,4-diaminopentanoic acid' },
    { smiles => 'C(CN)C=O', iupac => '3-aminopropanal' },
    { smiles => 'C(CCN)CCN', iupac => 'pentane-1,5-diamine' },

    # From BBv2 P-62.2.2.1
    { smiles => 'C(C)N(CC)CC', iupac => 'N,N-diethylethanamine' },
    { smiles => 'ClCCNCCC', iupac => 'N-(2-chloroethyl)propan-1-amine' },
    { smiles => 'C(C)N(CCCC)CCC', iupac => 'N-ethyl-N-propylbutan-1-amine' },

    # From BBv2 P-62.2.2.2
    { smiles => 'CC(C#CCN(CCC)CCC)=C', iupac => '4-methyl-N,N-dipropylpent-4-en-2-yn-1-amine' },
    { smiles => 'CN(C(C)C=CC1CC=C(CC1)C)C', iupac => 'N,N-dimethyl-4-(4-methylcyclohex-3-en-1-yl)but-3-en-2-amine', AUTHOR => 1 },

    { smiles => 'CC(C)C(C)N1CCCCC(C1=O)NC', iupac => '3-(methylamino)-1-(3-methylbutan-2-yl)azepan-2-one', AUTHOR => 1 }, # PubChem 58916315 # FIXME: Misses methyl
    { smiles => 'C1CNCCC1CNCCO', iupac => '2-(piperidin-4-ylmethylamino)ethanol', AUTHOR => 1 }, # PubChem 14950460 # FIXME: Very close
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    diag 'test supposed to fail with AUTHOR_TESTING' if $case->{AUTHOR} && $ok;
}
