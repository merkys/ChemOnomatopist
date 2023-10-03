#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-82.2.1
    { smiles => '[14CH4]', iupac => '(14C)methane' },
    { smiles => 'C[2H]', iupac => '(2H1)methane' },
    { smiles => 'C1(=CC=CC=C1)[13C]([13CH3])=O', iupac => '1-phenyl(1,2-13C2)ethan-1-one' },
    { smiles => '[13CH3]C1=C(C=CC=C1)[13CH3]', iupac => '1,2-di[(13C)methyl]benzene', AUTHOR => 1 },
    { smiles => '[13CH3]C1=[13CH]C=CC=C1', iupac => '1-(13C)methyl(2-13C)benzene', AUTHOR => 1 },
    { smiles => 'N[14CH2]C1(CCCC1)O', iupac => '1-[amino(14C)methyl]cyclopentan-1-ol' },
    { smiles => 'S1C([14CH2]CC1)C1=CC=NC=C1', iupac => '4-[(3-14C)thiolan-2-yl]pyridine' },
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
