#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-66.1.1.1.1.1
    { smiles => 'C(CCCCC)(=O)N', iupac => 'hexanamide' },
    { smiles => 'C(CCCC(=O)N)(=O)N', iupac => 'pentanediamide' },

    { smiles => 'C(C(CC(=O)N)C(=O)N)C(=O)N', iupac => 'propane-1,2,3-tricarboxamide', AUTHOR => 1 }, # BBv2 P-66.1.1.1.1.2

    # From BBv2 P-66.1.1.1.2.4
    { smiles => 'C(C=C)(=O)N', iupac => 'prop-2-enamide' },
    { smiles => 'OC(C(=O)N)C', iupac => '2-hydroxypropanamide' },

    # From BBv2 P-66.1.3
    { smiles => 'N1(CCCCC1)C(C)=O', iupac => '1-(piperidin-1-yl)ethan-1-one' },
    { smiles => 'N1(CCCC2=CC=CC=C12)C(CC)=O', iupac => '1-(3,4-dihydroquinolin-1(2H)-yl)propan-1-one', AUTHOR => 1 },

    { smiles => 'CC(CCOC)NC(=O)CCCCCCN', iupac => '7-amino-N-(4-methoxybutan-2-yl)heptanamide', AUTHOR => 1 }, # PubChem 64604850
    { smiles => 'C1CCN(CC1)CCCCC(=O)N=C(CC(=N)C2=CC=NC=C2)N', iupac => 'N-(1-amino-3-imino-3-pyridin-4-ylpropylidene)-5-piperidin-1-ylpentanamide', AUTHOR => 1 }, # PubChem 90937303

    { smiles => 'C(C)(=O)N(C(C1=CC=CC=C1)=O)C(CCCl)=O', iupac => 'N-acetyl-N-(3-chloropropanoyl)benzamide', AUTHOR => 1 }, # BBv2 P-66.1.2.1
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
