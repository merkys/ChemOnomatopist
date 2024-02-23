#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-44.4.1.11.1
    { smiles => 'C(CCCC)[2H]', iupac => '(1-2H1)pentane', AUTHOR => 1 },
    { smiles => 'C1(CCCCC1)([2H])[2H]', iupac => '(1,1-2H2)cyclohexane' },
    { smiles => '[14CH2]1CCCCC1', iupac => '(14C1)cyclohexane', AUTHOR => 1 },

    # From BBv3 P-44.4.1.11.2
    { smiles => '[14CH2]1CCCC1', iupac => '(14C1)cyclopentane', AUTHOR => 1 },
    { smiles => 'C1(CCCC1)[2H]', iupac => '(2H1)cyclopentane' },

    # From BBv3 P-44.4.1.11.3
    { smiles => 'C1(=CC=CC=C1)[3H]', iupac => '(3H1)benzene' },
    { smiles => 'C1(=CC=CC=C1)[2H]', iupac => '(2H1)benzene' },

    # From BBv3 P-44.4.1.11.4
    { smiles => 'N1=C(C=CC=C1)[2H]', iupac => '(2-2H)pyridine', AUTHOR => 1 },
    { smiles => 'N1=CC(=CC=C1)[2H]', iupac => '(3-2H)pyridine', AUTHOR => 1 },

    # From BBv3 P-82.2.1
    { smiles => '[14CH4]', iupac => '(14C)methane' },
    { smiles => 'Cl[12CH](Cl)Cl', iupac => 'trichloro(12C)methane' },
    { smiles => 'C[2H]', iupac => '(2H1)methane' },
    { smiles => 'ClC([2H])([2H])Cl', iupac => 'dichloro(2H2)methane' },
    { smiles => 'C(OC1=CC=CC=C1)([2H])([2H])[2H]', iupac => '(2H3)methoxybenzene', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)[13C]([13CH3])=O', iupac => '1-phenyl(1,2-13C2)ethan-1-one' },
    { smiles => '[13CH3]C1=C(C=CC=C1)[13CH3]', iupac => '1,2-di[(13C)methyl]benzene', AUTHOR => 1 },
    { smiles => '[13CH3]C1=[13CH]C=CC=C1', iupac => '1-(13C)methyl(2-13C)benzene', AUTHOR => 1 },
    { smiles => '[2H]CCO', iupac => '(2-H1)ethan-1-ol', AUTHOR => 1 },
    { smiles => '[12CH3]CO', iupac => '(2-C)ethan-1-ol', AUTHOR => 1 },
    { smiles => 'N[14CH2]C1(CCCC1)O', iupac => '1-[amino(14C)methyl]cyclopentan-1-ol' },
    { smiles => 'NCC1(CCCC1)[18OH]', iupac => '1-(aminomethyl)cyclopentan-1-(18O)ol', AUTHOR => 1 },
    { smiles => '[131I]C1=CC=C2C=3C=CC(=CC3CC2=C1)NC(C)=O', iupac => 'N-[7-(131I)iodo-9H-fluoren-2-yl]acetamide', AUTHOR => 1 },
    { smiles => 'C(C)OC([14CH2][14CH2]C(=O)[O-])=O.[Na+]', iupac => 'sodium 4-ethoxy-4-oxo(2,3-14C2)butanoate', AUTHOR => 1 },
    { smiles => 'S1C([14CH2]CC1)C1=CC=NC=C1', iupac => '4-[(3-14C)thiolan-2-yl]pyridine' },
    { smiles => '[35Cl]C(C[2H])C(CC)C([2H])([2H])[2H]', iupac => '2-(35Cl)chloro-3-(2H3)methyl(1-2H1)pentane', AUTHOR => 1 },

    # From BBv2 P-82.2.2.1
    { smiles => '[13CH3]C1=NC=CC=C1C', iupac => '2-(13C)methyl-3-methylpyridine' },
    { smiles => 'C(C([2H])[2H])C(CO)C(CCC)CC', iupac => '2-(2,2-2H2)ethyl-3-ethylhexan-1-ol' },

    { smiles => 'C[14CH2]C(C)[2H]', iupac => '(2-14C,3-2H1)butane' }, # From BBv3 P-82.5.2
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac}, $case->{smiles} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    if( $case->{AUTHOR} && $ok ) {
        diag 'test supposed to fail with AUTHOR_TESTING' .
             ( $case->{AUTHOR} !~ /^1$/ ? ': ' . $case->{AUTHOR} : '' );
    }
}
