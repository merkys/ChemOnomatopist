#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-44.4.1.11.1
    { smiles => 'C(CCCC)[2H]', iupac => '(1-2H1)pentane' },
    { smiles => 'C1(CCCCC1)([2H])[2H]', iupac => '(1,1-2H2)cyclohexane' },
    { smiles => '[14CH2]1CCCCC1', iupac => '(14C1)cyclohexane', AUTHOR => 1 },

    # From BBv3 P-44.4.1.11.2
    { smiles => '[14CH2]1CCCC1', iupac => '(14C1)cyclopentane', AUTHOR => 1 },
    { smiles => 'C1(CCCC1)[2H]', iupac => '(2H1)cyclopentane' },

    # From BBv3 P-44.4.1.11.3
    { smiles => 'C1(=CC=CC=C1)[3H]', iupac => '(3H1)benzene', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)[2H]', iupac => '(2H1)benzene', AUTHOR => 1 },

    # From BBv3 P-44.4.1.11.4
    { smiles => 'N1=C(C=CC=C1)[2H]', iupac => '(2-2H)pyridine' },
    { smiles => 'N1=CC(=CC=C1)[2H]', iupac => '(3-2H)pyridine' },

    # From BBv3 P-82.2.1
    { smiles => '[14CH4]', iupac => '(14C)methane' },
    { smiles => 'Cl[12CH](Cl)Cl', iupac => 'trichloro(12C)methane' },
    { smiles => 'C[2H]', iupac => '(2H1)methane' },
    { smiles => 'ClC([2H])([2H])Cl', iupac => 'dichloro(2H2)methane' },
    { smiles => 'C(OC1=CC=CC=C1)([2H])([2H])[2H]', iupac => '(2H3)methoxybenzene', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)[13C]([13CH3])=O', iupac => '1-phenyl(1,2-13C2)ethan-1-one' },
    { smiles => '[13CH3]C1=C(C=CC=C1)[13CH3]', iupac => '1,2-di[(13C)methyl]benzene', AUTHOR => 1 },
    { smiles => '[13CH3]C1=[13CH]C=CC=C1', iupac => '1-(13C)methyl(2-13C)benzene', AUTHOR => 1 },
    { smiles => '[2H]CCO', iupac => '(2-2H1)ethan-1-ol' },
    { smiles => '[12CH3]CO', iupac => '(2-12C)ethan-1-ol' },
    { smiles => 'N[14CH2]C1(CCCC1)O', iupac => '1-[amino(14C)methyl]cyclopentan-1-ol' },
    { smiles => 'NCC1(CCCC1)[18OH]', iupac => '1-(aminomethyl)cyclopentan-1-(18O)ol', AUTHOR => 1 },
    { smiles => '[131I]C1=CC=C2C=3C=CC(=CC3CC2=C1)NC(C)=O', iupac => 'N-[7-(131I)iodo-9H-fluoren-2-yl]acetamide', AUTHOR => 1 },
    { smiles => 'C(C)OC([14CH2][14CH2]C(=O)[O-])=O.[Na+]', iupac => 'sodium 4-ethoxy-4-oxo(2,3-14C2)butanoate', AUTHOR => 1 },
    { smiles => 'S1C([14CH2]CC1)C1=CC=NC=C1', iupac => '4-[(3-14C)thiolan-2-yl]pyridine' },
    { smiles => '[35Cl]C(C[2H])C(CC)C([2H])([2H])[2H]', iupac => '2-(35Cl)chloro-3-(2H3)methyl(1-2H1)pentane' },

    # From BBv2 P-82.2.2.1
    { smiles => '[13CH3]C1=NC=CC=C1C', iupac => '2-(13C)methyl-3-methylpyridine' },
    { smiles => 'C(C([2H])[2H])C(CO)C(CCC)CC', iupac => '2-(2,2-2H2)ethyl-3-ethylhexan-1-ol' },

    # From BBv3 P-82.5.1
    { smiles => 'FC(C[2H])(F)F', iupac => '1,1,1-trifluoro(2-2H1)ethane', AUTHOR => 1 },
    { smiles => 'ClC1=C(C(=CC=C1)F)[2H]', iupac => '1-chloro-3-fluoro(2-2H)benzene', AUTHOR => 1 },
    { smiles => 'COC1=C(C(=C(C(=C1[3H])[3H])[3H])[3H])O', iupac => '2-methoxy(3,4,5,6-3H4)phenol', AUTHOR => 1 },

    # From BBv3 P-82.5.2
    { smiles => 'C[14CH2]CC', iupac => '(2-14C)butane' },
    { smiles => 'CC([14CH2]C)([2H])[2H]', iupac => '(3-14C,2,2-2H2)butane' },
    { smiles => 'C[14CH2]C(C)[2H]', iupac => '(2-14C,3-2H1)butane' },
    { smiles => 'C1(=CC(=CC=C1)[3H])O', iupac => '(3-3H)phenol' },
    { smiles => 'C([C@@H](C)O)[2H]', iupac => '(2R)-(1-2H1)propan-2-ol', AUTHOR => 1 },
    { smiles => 'C[C@@H](C[C@@H](C)[2H])[3H]', iupac => '(2S,4R)-(4-2H1,2-3H1)pentane', AUTHOR => 1 },

    # From BBv3 P-82.6.3.3
    { smiles => 'C1=C(C=CC2=CC=CC=C12)[15N]=NC1=CC=CC=C1', iupac => '1-(naphthalen-2-yl)-2-phenyl(1-15N)diazene' },
    { smiles => 'C(CC)=[15N]N', iupac => '1-propylidene(1-15N)hydrazine', AUTHOR => 1 },
    { smiles => 'C(C)S[34S]SCCC(=O)O', iupac => '3-[ethyl(2-34S)trisulfanyl]propanoic acid', AUTHOR => 1 },
    { smiles => 'ClC1=C(C=CC2=CC=CC=C12)[15N]=[N+](C1=CC=CC=C1)[O-]', iupac => '1-(1-chloronaphthalen-2-yl)-2-phenyl(1-15N)diazene 2-oxide', AUTHOR => 1 },

    { smiles => '[14C](CCCC)(=O)[O][3H]', iupac => '(1-14C)pentan(3H)oic acid', AUTHOR => 1 }, # From BBv3 P-82.2.4
    { smiles => 'CCO[2H]', iupac => 'ethan(2H)ol', AUTHOR => 1 }, # From BBv3 P-82.6.1.1
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
