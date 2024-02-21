#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-68.3.1.3.1
    { smiles => 'N(=NC=O)C=O', iupac => 'diazenedicarbaldehyde', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)N=NC#N', iupac => 'phenyldiazenecarbonitrile' },
    { smiles => 'CN=NCC(=O)O', iupac => '(methyldiazenyl)acetic acid' },
    { smiles => 'N(=N)CCC(=O)O', iupac => '3-diazenylpropanoic acid', AUTHOR => 1 },
    { smiles => 'N(=N)C=1C=C(C(=O)O)C=CC1N=N', iupac => '3,4-bis(diazenyl)benzoic acid' },

    # From BBv3 P-68.3.1.3.2.1
    { smiles => 'CN=NC', iupac => 'dimethyldiazene' },
    { smiles => 'C1(=CC=CC=C1)N=NC1=CC=CC=C1', iupac => 'diphenyldiazene' },
    { smiles => 'ClC=1C=C(C=CC1)N=NC1=CC=C(C=C1)Cl', iupac => '(3-chlorophenyl)(4-chlorophenyl)diazene' },
    { smiles => 'C1(=CC=CC2=CC=CC=C12)N=NC1=CC2=CC=CC=C2C=C1', iupac => '(naphthalen-1-yl)(naphthalen-2-yl)diazene' },

    # From BBv3 P-68.3.1.3.2.2
    { smiles => 'C(=C)N=NC', iupac => 'ethenyl(methyl)diazene', AUTHOR => 1 },
    { smiles => 'C1=C(C=CC2=CC=CC=C12)N=NC1=CC=CC=C1', iupac => '(naphthalen-2-yl)(phenyl)diazene', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)N=NC1=CC=C(C=C1)S(=O)(=O)O', iupac => '4-(phenyldiazenyl)benzene-1-sulfonic acid' },
    { smiles => 'ClC1=CC(=C(C=C1)N=NC1=C(C=CC2=CC=CC=C12)N)C', iupac => '1-[(4-chloro-2-methylphenyl)diazenyl]naphthalen-2-amine' },
    { smiles => 'OC1=C(C2=CC=CC=C2C=C1)N=NC1=CC=C(C=C1)S(=O)(=O)O', iupac => '4-[(2-hydroxynaphthalen-1-yl)diazenyl]benzene-1-sulfonic acid' },
    { smiles => 'N(=NC1=CC=C(C(=O)O)C=C1)C1=CC=C(C(=O)O)C=C1', iupac => '4,4\'-diazenediyldibenzoic acid', AUTHOR => 1 },

    # From BBv3 P-68.3.1.3.2.3
    { smiles => 'C1=C(C=CC2=CC3=CC=CC=C3C=C12)N=NC1=CC=C2C=CC(=CC2=C1)N=NC1=CC=CC=C1', iupac => '{7-[(anthracen-2-yl)diazenyl]naphthalen-2-yl}(phenyl)diazene', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)N=NC1=C(C2=C(C(=CC=C2C=C1)N=NC1=CC=CC=C1)O)O', iupac => '2,7-bis(phenyldiazenyl)naphthalene-1,8-diol' },
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
