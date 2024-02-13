#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-66.4.1.1
    { smiles => 'C(CCCCC)(N)=N', iupac => 'hexanimidamide', AUTHOR => 1 },
    { smiles => 'C1(CCCCC1)C(N)=N', iupac => 'cyclohexanecarboximidamide' },
    { smiles => 'C(C)(N)=N', iupac => 'ethanimidamide', AUTHOR => 1 },
    { smiles => 'CS(N)=N', iupac => 'methanesulfinimidamide' },
    { smiles => 'C(N)=N', iupac => 'methanimidamide', AUTHOR => 1 },
    { smiles => 'C(CCCC(N)=N)(N)=N', iupac => 'pentanediimidamide', AUTHOR => 1 },
    { smiles => '[SiH2]([SiH2]C(N)=N)C(N)=N', iupac => 'disilane-1,2-dicarboximidamide', AUTHOR => 1 },
    { smiles => 'C(CCC(N)=N)(N)=N', iupac => 'butanediimidamide', AUTHOR => 1 },
    { smiles => 'C(C(N)=N)(N)=N', iupac => 'ethanediimidamide', AUTHOR => 1 },
    { smiles => 'C=1(C(=CC=CC1)C(N)=N)C(N)=N', iupac => 'benzene-1,2-dicarboximidamide', AUTHOR => 1 },
    { smiles => 'C1(=CC=C(C=C1)C(N)=N)C(N)=N', iupac => 'benzene-1,4-dicarboximidamide', AUTHOR => 1 },
    { smiles => 'C(C)NC(=N)C1(CCCCC1)C(N(C)C)=N', iupac => 'N\'\'1-ethyl-N1,N1-dimethylcyclohexane-1,1-dicarboximidamide', AUTHOR => 1 },

    { smiles => 'NS(=N)CCC(=O)O', iupac => '3-(S-aminosulfinimidoyl)propanoic acid' }, # From BBv3 P-66.4.1.3.4

    # From BBv3 P-66.4.1.3.5
    { smiles => 'C(C)(NC1=CC=C(C(=O)O)C=C1)=N', iupac => '4-ethanimidamidobenzoic acid', AUTHOR => 1 },
    { smiles => 'C(C)S(NC1=C(C(=O)O)C=CC=C1)(=N)=N', iupac => '2-(ethanesulfonodiimidamido)benzoic acid', AUTHOR => 1 },
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
