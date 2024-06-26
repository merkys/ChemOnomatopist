#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'C1=CC=CC=C1',  iupac => 'benzene' },
    { smiles => 'C=1C=CC=CC=1', iupac => 'benzene' },
    { smiles => 'c1ccccc1',     iupac => 'benzene' },

    { smiles => 'FC1=CC=CC=C1', iupac => 'fluorobenzene' },
    { smiles => 'FC=1C(=C(C(=C(C1)F)F)F)F', iupac => 'pentafluorobenzene' },
    { smiles => 'FC1=C(C(=C(C(=C1F)F)F)F)F', iupac => 'hexafluorobenzene' },

    { smiles => 'C1(C=CCCC1)O', iupac => 'cyclohex-2-en-1-ol' }, # From BBv2 P-14.3.2
    { smiles => 'N1=C(C=NC=C1)C(=O)O', iupac => 'pyrazinecarboxylic acid', AUTHOR => 1 }, # From BBv2 P-14.3.4.3
    { smiles => 'C1(C=CCCC1)N', iupac => 'cyclohex-2-en-1-amine' }, # From BBv2 P-14.4
    { smiles => 'C1(CCCC1)=CC(=O)O', iupac => 'cyclopentylideneacetic acid' }, # From BBv2 P-15.6.2

    # From BBv2 P-22.1.3
    { smiles => 'CC1=CC=CC=C1', iupac => 'toluene' },
    { smiles => 'C=1(C(=CC=CC1)C)C', iupac => '1,2-xylene' },
    { smiles => 'C1(=CC(=CC=C1)C)C', iupac => '1,3-xylene' },
    { smiles => 'C1(=CC=C(C=C1)C)C', iupac => '1,4-xylene' },
    { smiles => 'CC1=CC(=CC(=C1)C)C', iupac => '1,3,5-trimethylbenzene' },

    { smiles => 'N(=O)C1=CC=CC=C1',  iupac => 'nitrosobenzene' },  # From BBv2 P-59.1.9
    { smiles => 'Br(=O)C1=CC=CC=C1', iupac => 'bromosylbenzene' }, # From BBv2 P-67.1.4.5

    # From BBv3 P-22.2.1
    { smiles => 'O1C=CC=C1', iupac => 'furan' },
    { smiles => 'N1C=NC=C1', iupac => '1H-imidazole' },
    { smiles => 'O1C=NC=C1', iupac => '1,3-oxazole' },
    { smiles => 'S1C=NC=C1', iupac => '1,3-thiazole' },
    { smiles => '[Se]1C=NC=C1', iupac => '1,3-selenazole' },
    { smiles => '[Te]1C=NC=C1', iupac => '1,3-tellurazole' },
    { smiles => 'O1N=CC=C1', iupac => '1,2-oxazole' },
    { smiles => 'S1N=CC=C1', iupac => '1,2-thiazole' },
    { smiles => '[Se]1N=CC=C1', iupac => '1,2-selenazole' },
    { smiles => '[Te]1N=CC=C1', iupac => '1,2-tellurazole' },
    { smiles => 'O1CC=CC=C1', iupac => '2H-pyran' },
    { smiles => 'S1CC=CC=C1', iupac => '2H-thiopyran' },
    { smiles => '[Se]1CC=CC=C1', iupac => '2H-selenopyran' },
    { smiles => '[Te]1CC=CC=C1', iupac => '2H-telluropyran' },
    { smiles => 'N1=CC=NC=C1', iupac => 'pyrazine' },
    { smiles => 'N1N=CC=C1', iupac => '1H-pyrazole' },
    { smiles => 'N1=NC=CC=C1', iupac => 'pyridazine' },
    { smiles => 'N1=CC=CC=C1', iupac => 'pyridine' },
    { smiles => 'N1=CN=CC=C1', iupac => 'pyrimidine' },
    { smiles => 'N1C=CC=C1', iupac => '1H-pyrrole' },
    { smiles => '[Se]1C=CC=C1', iupac => 'selenophene' },
    { smiles => '[Te]1C=CC=C1', iupac => 'tellurophene' },
    { smiles => 'S1C=CC=C1', iupac => 'thiophene' },
    { smiles => 'O1CNCC1', iupac => '1,3-oxazolidine' },
    { smiles => 'O1NCCC1', iupac => '1,2-oxazolidine' },
    { smiles => 'S1CNCC1', iupac => '1,3-thiazolidine' },
    { smiles => 'S1NCCC1', iupac => '1,2-thiazolidine' },
    { smiles => '[Se]1CNCC1', iupac => '1,3-selenazolidine' },
    { smiles => '[Se]1NCCC1', iupac => '1,2-selenazolidine' },
    { smiles => '[Te]1CNCC1', iupac => '1,3-tellurazolidine' },
    { smiles => '[Te]1NCCC1', iupac => '1,2-tellurazolidine' },
    { smiles => 'N1CCCC1', iupac => 'pyrrolidine' },
    { smiles => 'N1CCOCC1', iupac => 'morpholine' },
    { smiles => 'N1CCSCC1', iupac => 'thiomorpholine' },
    { smiles => 'N1CC[Se]CC1', iupac => 'selenomorpholine' },
    { smiles => 'N1CC[Te]CC1', iupac => 'telluromorpholine' },
    { smiles => 'N1NCCC1', iupac => 'pyrazolidine' },
    { smiles => 'N1CNCC1', iupac => 'imidazolidine' },
    { smiles => 'N1CCCCC1', iupac => 'piperidine' },
    { smiles => 'N1CCNCC1', iupac => 'piperazine' },

    { smiles => 'S1CCCCCCCCCCC1', iupac => 'thiacyclododecane' }, # From BBv3 P-22.2.3.1

    # From BBv2 P-61.5.1
    { smiles => 'CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', iupac => '2-methyl-1,3,5-trinitrobenzene' },

    # From BBv2 P-22.2.2.1.1
    { smiles => 'S1C=CC=CC=C1', iupac => 'thiepine' },
    { smiles => 'O1CCCCCCC1',   iupac => 'oxocane' },

    # From BBv2 P-22.2.2.1.5.2
    { smiles => 'O1CCC1', iupac => 'oxetane' },
    { smiles => 'N1CCC1', iupac => 'azetidine' },

    { smiles => 'C1NNN1', iupac => 'triazetidine' }, # PubChem 23033474

    # From BBv2 P-22.2.2.1.2
    { smiles => 'N1=CC=CN=CC=C1', iupac => '1,5-diazocine' },
    { smiles => 'O1COCC1', iupac => '1,3-dioxolane' },

    # From BBv2 P-22.2.2.1.3
    { smiles => 'S1C=NC=C1', iupac => '1,3-thiazole' },
    { smiles => 'O1SCCC1',   iupac => '1,2-oxathiolane' },
    { smiles => 'O1SCCCSC1', iupac => '1,2,6-oxadithiepane' },
    { smiles => 'O1N=CC=P1', iupac => '1,2,5-oxazaphosphole' },

    { smiles => 'O1CCCCCCOCCCCCCCCCC1', iupac => '1,8-dioxacyclooctadecane' }, # From BBv2 P-22.2.4

    # From BBv2 P-25.2.2.1.2
    { smiles => 'O1C=CC=CC=COC=CC=CC=CC=CC=C1', iupac => '1,8-dioxacyclooctadeca-2,4,6,9,11,13,15,17-octaene' },
    { smiles => 'O1CC=NC=CC=NC=CN=CC=C1', iupac => '1-oxa-4,8,11-triazacyclotetradeca-3,5,7,9,11,13-hexaene' },

    # From BBv2 P-63.1.1.1
    { smiles => 'C1=CC=CC=C1O', iupac => 'phenol' },
    { smiles => 'BrC1=C(C=CC=C1)O', iupac => '2-bromophenol' },

    { smiles => 'C(=O)(O)C1=CC=CC=C1', iupac => 'benzoic acid' },

    { smiles => 'N1C(CCCCC1)=S', iupac => 'azepane-2-thione' }, # From BBv2 P-64.6.1

    # From BBv2 P-31.1.3.1
    { smiles => 'C1C=CCCC1', iupac => 'cyclohexene' },
    { smiles => 'C1C=CCC=C1', iupac => 'cyclohexa-1,4-diene' },

    # From BBv2 P-31.1.3.2
    { smiles => 'O1C=COCCOCCOCC1', iupac => '1,4,7,10-tetraoxacyclododec-2-ene' },
    { smiles => 'O1CC=NCCCCCCCC1', iupac => '1-oxa-4-azacyclododec-3-ene' },
    { smiles => '[SiH2]1CC#CC=CC=CCC[SiH2]CCCCCCCCC1', iupac => '1,11-disilacycloicosa-5,7-dien-3-yne' },
    { smiles => '[SiH2]1CCCCCCCC[SiH2]CC=CC=CC=CC#CC1', iupac => '1,10-disilacycloicosa-12,14,16-trien-18-yne' },

    # From BBv2 P-31.1.3.3
    { smiles => 'C=1=C=C=C=C=C=C=C=C=C=C1', iupac => 'cycloundecaundecaene' },

    # From BBv2 P-31.1.3.4
    { smiles => 'C=CC1=CC=CC=C1', iupac => 'ethenylbenzene' },
    { smiles => 'C=C1C=CC=C1', iupac => '5-methylidenecyclopenta-1,3-diene' },

    { smiles => 'S1CCNCCC1', iupac => '1,4-thiazepane' }, # From BBv2 P-31.2.3.2

    # From BBv2 P-14.5.1
    { smiles => 'C1CCCCC1(C)CC', iupac => '1-ethyl-1-methylcyclohexane' },
    { smiles => 'CCC1CCC(C)CC1', iupac => '1-ethyl-4-methylcyclohexane' },

    { smiles => 'CC1=NC(=CC=C1)C', iupac => '2,6-dimethylpyridine' },

    { smiles => 'C1CCCCC1(C(C)(C)C)(CCCC)', iupac => '1-butyl-1-tert-butylcyclohexane' }, # Simplified version of example from BBv2 P-14.5.1

    { smiles => 'C(=O)(O)CC([Br])([Br])C1CCCCC1', iupac => '3,3-dibromo-3-cyclohexylpropanoic acid' },
    { smiles => 'ClC=1C=CC=CC=1C(F)(F)C(F)(F)F',  iupac => '1-chloro-2-(pentafluoroethyl)benzene' }, # From BBv2 P-14.3.4.5
    { smiles => 'CC(CCC)C1=CC=C(C=C1)C(CC)CC', iupac => '1-(pentan-2-yl)-4-(pentan-3-yl)benzene' }, # From BBv2 P-14.5.4

    # From BBv2 P-14.5.3
    { smiles => 'C(C)(C)(C)C=1C=CC=C(C(C)CC)C=1', iupac => '1-(butan-2-yl)-3-tert-butylbenzene' },

    # From BBv3 P-44.4.1.5
    { smiles => 'O1CS[SnH2]SCOCCCC1', iupac => '1,7-dioxa-3,5-dithia-4-stannacycloundecane' },
    { smiles => 'O1CCS[SnH2]SCCOCC1', iupac => '1,9-dioxa-4,6-dithia-5-stannacycloundecane' },

    { smiles => 'O=C1NC(=O)NC=C1C', iupac => '5-methylpyrimidine-2,4(1H,3H)-dione', AUTHOR => 1 }, # thymine
    { smiles => 'c1cc(oc1)C=O', iupac => 'furan-2-carbaldehyde', AUTHOR => 1 }, # furfural

    { smiles => 'O(c1cc(cc(OC)c1OC)CCN)C', iupac => '2-(3,4,5-trimethoxyphenyl)ethanamine' }, # mescaline, SMILES from Wikipedia
    { smiles => 'COC=1C=C(C=C(C1OC)OC)CCN', iupac => '2-(3,4,5-trimethoxyphenyl)ethanamine' }, # mescaline

    { smiles => 'C1=CC=CC=CC=CC=CC#C1', iupac => 'cyclododeca-1,3,5,7,9-pentaen-11-yne' }, # From BBv3 P-54.2
    { smiles => 'C1CCCCC1C=O', iupac => 'cyclohexanecarbaldehyde' }, # From BBv2 P-66.6.1.1.3

    { smiles => 'S=C(CC1CC(CCC1)CC(CC)=O)CC', iupac => '1-[3-(2-sulfanylidenebutyl)cyclohexyl]butan-2-one' }, # From BBv2 P-64.7.3

    { smiles => 'C1=CC(=CC=C1F)Cl(=O)(=O)=O', iupac => '1-fluoro-4-perchlorylbenzene' }, # PubChem 24972904
    { smiles => 'CCC1CN1Cl(=O)(=O)=O', iupac => '2-ethyl-1-perchlorylaziridine' }, # PubChem 24973518
    { smiles => 'CC(=CC(=C(C)C)C(=C)C1=CC=CC=C1)C', iupac => '(5-methyl-3-propan-2-ylidenehexa-1,4-dien-2-yl)benzene' }, # PubChem 141889885

    { smiles => 'Cc1c(I)c(C)c(c(c1C)I)C', iupac => '1,4-diiodo-2,3,5,6-tetramethylbenzene' }, # COD 2013412

    # From BBv2 P-59.2.1.6
    { smiles => 'OC1CCC(CC1)C(CCCCCO)O', iupac => '1-(4-hydroxycyclohexyl)hexane-1,6-diol', AUTHOR => 1 },
    { smiles => 'O=C(CC1CC(C(C1)=O)=O)CC', iupac => '4-(2-oxobutyl)cyclopentane-1,2-dione' },

    { smiles => 'O=CCCCCCCC1CC(CCC1)C=O', iupac => '3-(7-oxoheptyl)cyclohexane-1-carbaldehyde', AUTHOR => 1 }, # From BBv2 P-59.2.1.7

    { smiles => 'S1C(=NC=C1)CC(=O)O', iupac => '(1,3-thiazol-2-yl)acetic acid' }, # From BBv2 P-15.6.1.2
    { smiles => 'S1C(=NC=C1)CCCCC(=O)O', iupac => '5-(1,3-thiazol-2-yl)pentanoic acid' }, # Synthetic, based on BBv2 P-15.6.1.2
    { smiles => 'CC1=C(C=C(C=C1)C(C)C)O', iupac => '2-methyl-5-(propan-2-yl)phenol' }, # From BBv2 P-63.1.1.2
    { smiles => 'S=C1OCCC1C#N', iupac => '2-sulfanylideneoxolane-3-carbonitrile' }, # From BBv2 P-65.6.3.5.1

    # From BBv2 P-64.7.1
    { smiles => 'NC1C(NCCCC1)=O', iupac => '3-aminoazepan-2-one' },
    { smiles => 'ClC=1C(C(=C(C(C1O)=O)Cl)O)=O', iupac => '2,5-dichloro-3,6-dihydroxycyclohexa-2,5-diene-1,4-dione' },
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
