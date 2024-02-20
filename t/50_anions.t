#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-72.2.2.2.1.1
    { smiles => 'C(C)(=O)[O-]', iupac => 'acetate', AUTHOR => 1 },
    { smiles => 'C(CC)(=O)O[O-]', iupac => 'propaneperoxoate', AUTHOR => 1 },
    { smiles => 'C(C)(O[O-])=S', iupac => 'ethaneperoxothioate', AUTHOR => 1 },
    { smiles => 'CC(=O)O[S-]', iupac => 'ethane(OS-thioperoxoate)', AUTHOR => 1 },
    { smiles => 'C(CC)([O-])=S', iupac => 'propanethioate' },
    { smiles => 'C(C)([O-])=S', iupac => 'ethanethioate' },
    { smiles => 'C1(=CC=CC=C1)S(=O)(=O)[O-]', iupac => 'benzenesulfonate', AUTHOR => 1 },
    { smiles => 'C(C1=CC=CC=C1)P([O-])CC1=CC=CC=C1', iupac => 'dibenzylphosphinite', AUTHOR => 1 },
    { smiles => 'N1=C(C=CC=C1C(=O)[O-])C(=O)[O-]', iupac => 'pyridine-2,6-dicarboxylate' },
    { smiles => 'N1C(=CC=C1)C([O-])=N', iupac => '1H-pyrrole-2-carboximidate', AUTHOR => 1 },
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
