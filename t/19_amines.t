#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CC(CN)O', iupac => '1-aminopropan-2-ol' },
    { smiles => 'C(CC(=O)N)C(=O)C(=O)O', iupac => '5-amino-2,5-dioxopentanoic acid', AUTHOR => 1 }, # PubChem 48: seems a bit strange
    { smiles => 'CC(CC(C(=O)O)N)N', iupac => '2,4-diaminopentanoic acid' },
    { smiles => 'C(CN)C=O', iupac => '3-aminopropanal' },
    { smiles => 'C(CCN)CCN', iupac => 'pentane-1,5-diamine', AUTHOR => 1 },

    { smiles => 'CC(CC1=CC=C(NCCC(C)C)C=C1)CC', iupac => '4-(2-methylbutyl)-N-(3-methylbutyl)aniline' }, # BBv2 P-14.5.4
    { smiles => 'NC1=CC=CC=C1', iupac => 'aniline' }, # BBv2 P-34.1.1.5
    { smiles => 'COC1=CC=C(NC2=CC=CC=C2)C=C1', iupac => '4-methoxy-N-phenylaniline' }, # BBv2 P-45.2.1
    { smiles => 'NC1=C(OC2=C(NC)C=CC=C2)C=CC(=C1)C', iupac => '2-(2-amino-4-methylphenoxy)-N-methylaniline', AUTHOR => 1 }, # BBv2 P-45.2.2
    { smiles => 'CNC', iupac => 'N-methylmethanamine' }, # BBv2 P-52.1.3

    # From BBv2 P-62.2.1.2
    { smiles => 'CN', iupac => 'methanamine' },
    { smiles => 'CC(CN)C', iupac => '2-methylpropan-1-amine' },

    # From BBv2 P-62.2.2.1
    { smiles => 'C(C)N(CC)CC', iupac => 'N,N-diethylethanamine' },
    { smiles => 'ClCCNCCC', iupac => 'N-(2-chloroethyl)propan-1-amine', AUTHOR => 1 }, # FIXME: This is strange
    { smiles => 'C(C)N(CCCC)CCC', iupac => 'N-ethyl-N-propylbutan-1-amine' },
    { smiles => 'C1(=CC=CC=C1)NC=1C=NC=CC1', iupac => 'N-phenylpyridin-3-amine' },

    # From BBv2 P-62.2.2.2
    { smiles => 'CC(C#CCN(CCC)CCC)=C', iupac => '4-methyl-N,N-dipropylpent-4-en-2-yn-1-amine' },
    { smiles => 'CN(C(C)C=CC1CC=C(CC1)C)C', iupac => 'N,N-dimethyl-4-(4-methylcyclohex-3-en-1-yl)but-3-en-2-amine', AUTHOR => 1 },
    { smiles => 'CN(C(C#C)CC)C', iupac => 'N,N-dimethylpent-1-yn-3-amine', AUTHOR => 1 },
    { smiles => 'C(=C)NCCCC', iupac => 'N-ethenylbutan-1-amine' },
    { smiles => 'CC(CN(CC(=C)C)CC(=C)C)(C)C', iupac => 'N-(2,2-dimethylpropyl)-2-methyl-N-(2-methylprop-2-en-1-yl)prop-2-en-1-amine', AUTHOR => 1 },
    { smiles => 'C1(CCCCC1)NC1=CC=CC=C1', iupac => 'N-cyclohexylaniline' },
    { smiles => 'O1C(=CC=C1)NC=1NC=CC1', iupac => 'N-(furan-2-yl)-1H-pyrrol-2-amine', AUTHOR => 1 },
    { smiles => 'C(CCC)NC1CC1', iupac => 'N-butylcyclopropanamine' },
    { smiles => 'C1=C(C=CC=2CCCCC12)NC1=CC2=CC=CC=C2C=C1', iupac => 'N-(5,6,7,8-tetrahydronaphthalen-2-yl)naphthalen-2-amine', AUTHOR => 1 },

    # From BBv2 P-62.2.3
    { smiles => 'NCCC(=O)O', iupac => '3-aminopropanoic acid' },

    { smiles => 'CC(C)C(C)N1CCCCC(C1=O)NC', iupac => '3-(methylamino)-1-(3-methylbutan-2-yl)azepan-2-one', AUTHOR => 1 }, # PubChem 58916315 # FIXME: Misses methylamino
    { smiles => 'C1CNCCC1CNCCO', iupac => '2-(piperidin-4-ylmethylamino)ethanol', AUTHOR => 1 }, # PubChem 14950460 # FIXME: Very close
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
