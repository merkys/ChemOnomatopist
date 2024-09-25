#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'CC(=O)CC(=O)CCCO', iupac => '7-hydroxyheptane-2,4-dione' },
    { smiles => 'CC(C)C(CCC(C)O)CCC(=O)C', iupac => '8-hydroxy-5-propan-2-ylnonan-2-one', AUTHOR => 1 }, # PubChem 537619

    { smiles => 'C(O)OO', iupac => 'hydroperoxymethanol' },
    { smiles => 'OOCC(=O)C', iupac => '1-hydroperoxypropan-2-one' }, # Unchecked

    # From BBv3 P-61.11.1
    { smiles => 'N(=[N+]=[N-])C1=CC(=C(C=C1)F)[N+](=O)[O-]', iupac => '4-azido-1-fluoro-2-nitrobenzene' },
    { smiles => '[N+](=O)([O-])NC(N(CCC)N=O)=N', iupac => 'N\'-nitro-N-nitroso-N-propylguanidine' },

    # From BBv3 P-61.11.2
    { smiles => 'Cl[Si]([SiH2][SiH2][SiH2][Si](C)(C)C)(Cl)Cl', iupac => '1,1,1-trichloro-5,5,5-trimethylpentasilane', AUTHOR => 1 },
    { smiles => 'BrCCCl', iupac => '1-bromo-2-chloroethane' },
    { smiles => 'N(=[N+]=[N-])C1=CC=C(C=C1)N=C=O', iupac => '1-azido-4-isocyanatobenzene' },

    # From BBv2 P-65.1.2.4
    { smiles => 'CC(=O)CCCC(=O)O', iupac => '5-oxohexanoic acid' },
    { smiles => 'C(=O)(O)C(O)C(C(=O)O)C(=O)C(=O)O', iupac => '1-hydroxy-3-oxopropane-1,2,3-tricarboxylic acid' },

    { smiles => 'O=CCCC(=O)O', iupac => '4-oxobutanoic acid', AUTHOR => 1 }, # BBv2 P-66.6.1.3 - FIXME

    { smiles => 'CCCCCCCCC(C=O)C(CC)O', iupac => '2-(1-hydroxypropyl)decanal' },
    { smiles => 'S=C(CC(=O)O)C', iupac => '3-sulfanylidenebutanoic acid' },

    { smiles => '[N+](=O)([O-])C', iupac => 'nitromethane' }, # BBv2 P-61.5.1

    { smiles => 'I(=O)(=O)(=O)CCC', iupac => 'periodylpropane', AUTHOR => 1 },
    { smiles => 'I(=O)(=O)(=O)C(C)C', iupac => '2-periodylpropane' },
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
