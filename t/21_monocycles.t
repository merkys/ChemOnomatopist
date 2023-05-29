#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C1CCNC1', iupac => 'pyrrolidine' },

    { smiles => 'C1=CC=CC=C1',  iupac => 'benzene' },
    { smiles => 'C=1C=CC=CC=1', iupac => 'benzene' },
    { smiles => 'c1ccccc1',     iupac => 'benzene' },

    # From BBv2 P-22.2.2.1.1
    { smiles => 'S1C=CC=CC=C1', iupac => 'thiepine' },
    { smiles => 'O1CCCCCCC1',   iupac => 'oxocane' },

    # From BBv2 P-22.2.2.1.5.2
    { smiles => 'O1CCC1',     iupac => 'oxetane' },
    { smiles => 'N1CCC1',     iupac => 'azetidine' },

    { smiles => 'N1NNC1', iupac => 'triazetidine' }, # Not checked

    # From BBv2 P-22.2.2.1.2
    { smiles => 'N1=CC=CN=CC=C1', iupac => '1,5-diazocine' },
    { smiles => 'O1COCC1', iupac => '1,3-dioxolane' },

    # From BBv2 P-22.2.2.1.3
    { smiles => 'S1C=NC=C1', iupac => '1,3-thiazole', AUTHOR => 1 },
    { smiles => 'O1SCCC1',   iupac => '1,2-oxathiolane', AUTHOR => 1 },
    { smiles => 'O1SCCCSC1', iupac => '1,2,6-oxadithiepane', AUTHOR => 1 },

    { smiles => 'O1CCCCCCOCCCCCCCCCC1', iupac => '1,8-dioxacyclooctadecane' }, # From BBv2 P-22.2.4

    # From BBv2 P-25.2.2.1.2
    { smiles => 'O1C=CC=CC=COC=CC=CC=CC=CC=C1', iupac => '1,8-dioxacyclooctadeca-2,4,6,9,11,13,15,17-octaene', AUTHOR => 1 },
    { smiles => 'O1CC=NC=CC=NC=CN=CC=C1', iupac => '1-oxa-4,8,11-triazacyclotetradeca-3,5,7,9,11,13-hexaene', AUTHOR => 1 },

    { smiles => 'C1=CC=CC=C1O', iupac => 'phenol' },
    { smiles => 'C(=O)(O)C1=CC=CC=C1', iupac => 'benzoic acid' },

    { smiles => 'N1C(CCCCC1)=S', iupac => 'azepane-2-thione', AUTHOR => 1 }, # From BBv2 P-64.6.1

    # From BBv2 P-31.1.3.1
    { smiles => 'C1C=CCCC1', iupac => 'cyclohexene' },
    { smiles => 'C1C=CCC=C1', iupac => 'cyclohexa-1,4-diene' },

    # From BBv2 P-31.1.3.4
    { smiles => 'C=CC1=CC=CC=C1', iupac => 'ethenylbenzene', AUTHOR => 1 },
    { smiles => 'C=C1C=CC=C1', iupac => '5-methylidenecyclopenta-1,3-diene' },

    # From BBv2 P-14.5.1
    { smiles => 'C1CCCCC1(C)CC', iupac => '1-ethyl-1-methylcyclohexane' },
    { smiles => 'CCC1CCC(C)CC1', iupac => '1-ethyl-4-methylcyclohexane' },

    { smiles => 'CC1=NC(=CC=C1)C', iupac => '2,6-dimethylpyridine' },

    { smiles => 'C1CCCCC1(C(C)(C)C)(CCCC)', iupac => '1-butyl-1-tert-butylcyclohexane' }, # Simplified version of example from BBv2 P-14.5.1

    { smiles => 'C(=O)(O)CC([Br])([Br])C1CCCCC1', iupac => '3,3-dibromo-3-cyclohexylpropanoic acid' },
    { smiles => 'ClC=1C=CC=CC=1C(F)(F)C(F)(F)F',  iupac => '1-chloro-2-(pentafluoroethyl)benzene', AUTHOR => 1 }, # From BBv2 P-14.3.4.5
    { smiles => 'CC(CCC)C1=CC=C(C=C1)C(CC)CC', iupac => '1-(pentan-2-yl)-4-(pentan-3-yl)benzene', AUTHOR => 1 }, # From BBv2 P-14.5.4

    # From BBv2 P-14.5.3
    { smiles => 'C(C)(C)(C)C=1C=CC=C(C(C)CC)C=1', iupac => '1-(butan-2-yl)-3-tert-butylbenzene', AUTHOR => 1 },

    { smiles => 'O=C1NC(=O)NC=C1C', iupac => '5-methylpyrimidine-2,4(1H,3H)-dione', AUTHOR => 1 }, # thymine
    { smiles => 'c1cc(oc1)C=O', iupac => 'furan-2-carbaldehyde', AUTHOR => 1 }, # furfural

    { smiles => 'C1CCCCC1C=O', iupac => 'cyclohexanecarbaldehyde', AUTHOR => 1 }, # From BBv2 P-66.6.1.1.3

    { smiles => 'S=C(CC1CC(CCC1)CC(CC)=O)CC', iupac => '1-[3-(2-sulfanylidenebutyl)cyclohexyl]butan-2-one', AUTHOR => 1 }, # From BBv2 P-64.7.3
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    eval { is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
}
