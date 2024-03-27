#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-45.2.1
    { smiles => 'COC1=CC=C(NC2=CC=CC=C2)C=C1', iupac => '4-methoxy-N-phenylaniline' },
    { smiles => 'ClC1=CC(=C(C#N)C=C1)CC1=CC(=CC=C1)C#N', iupac => '4-chloro-2-[(3-cyanophenyl)methyl]benzonitrile' },
    { smiles => 'CC1=CC=C(C=C1)COC1=CC=CC=C1', iupac => '1-methyl-4-(phenoxymethyl)benzene', AUTHOR => 1 },
    { smiles => 'CN(C(C(CC1=CC(=C(C=C1)C)CC(C(=O)NC)C)C)=O)C', iupac => 'N,N,2-trimethyl-3-{4-methyl-3-[2-methyl-3-(methylamino)-3-oxopropyl]phenyl}propanamide', AUTHOR => 1 },
    { smiles => 'CS1(CC(CCC1)CSC1C[SH2]CCC1)C', iupac => '1,1-dimethyl-3-{[(1λ4-thian-3-yl)sulfanyl]methyl}-1λ4-thiane', AUTHOR => 1 },
    { smiles => 'C(C)C(C(C)C)CCC', iupac => '3-ethyl-2-methylhexane' },
    { smiles => 'CC(CC(CCC(=O)O)CCC)C', iupac => '6-methyl-4-propylheptanoic acid' },
    { smiles => 'BrC(C(C(C(=O)O)C(=C)C(CC)CBr)=C)(CC)C', iupac => '4-bromo-2-[3-(bromomethyl)pent-1-en-2-yl]-4-methyl-3-methylidenehexanoic acid', AUTHOR => 1 },
    { smiles => '[SiH2]([SiH3])SS[SiH2][Si](C)(C)C', iupac => '2-(disilanyldisulfanyl)-1,1,1-trimethyldisilane', AUTHOR => 1 },
    { smiles => 'C(=O)(O)C1=C(C=C(OC2=C(C(=C(C(=O)O)C=C2)P)[SH5])C=C1)[SH5]', iupac => '4-[4-carboxy-3-(λ6-sulfanyl)phenoxy]-2-phosphanyl-3-(λ6-sulfanyl)benzoic acid', AUTHOR => 1 },
    { smiles => '[81Br]C(C(C(CC(=O)O)C(CC)[81Br])C)C', iupac => '5-(81Br)bromo-3-[1-(81Br)bromopropyl]-4-methylhexanoic acid', AUTHOR => 1 },

    # From BBv3 P-45.2.2
    { smiles => 'NC1=C(OC2=C(NC)C=CC=C2)C=CC(=C1)C', iupac => '2-(2-amino-4-methylphenoxy)-N-methylaniline', AUTHOR => 1 },
    { smiles => 'BrC1=C(C(=CC2=CC(=CC=C12)[N+](=O)[O-])Cl)CCC1=C(C2=CC(=CC=C2C=C1F)F)F', iupac => '1-bromo-3-chloro-6-nitro-2-[2-(1,3,7-trifluoronaphthalen-2-yl)ethyl]naphthalene' },
    { smiles => 'C(=O)(O)C1=CC=C(C=C1)C(C=1C=C(C(=O)O)C=CC1)(C=1C=C(C(=O)O)C=CC1)C1=CC=C(C=C1)C(=O)O', iupac => '3,3\'-[bis(4-carboxyphenyl)methylene]dibenzoic acid', AUTHOR => 1 },
    { smiles => 'NC(C(CC=1C=CC(=C(C1)CCC(=O)NC)C)C)=O', iupac => '3-[5-(3-amino-2-methyl-3-oxopropyl)-2-methylphenyl]-N-methylpropanamide', AUTHOR => 1 },
    { smiles => 'BrC(C(CCC(=O)O)C(CC[N+](=O)[O-])Cl)C(C)Br', iupac => '5,6-dibromo-4-(1-chloro-3-nitropropyl)heptanoic acid' },
    { smiles => 'NC(CO)CCC(C(CCCO)C)CC(CCO)Cl', iupac => '2-amino-5-(2-chloro-4-hydroxybutyl)-6-methylnonane-1,9-diol' },
    { smiles => 'CC(C(CC=C)C=C(C)C)=CC', iupac => '5-methyl-4-(2-methylprop-1-en-1-yl)hepta-1,5-diene', AUTHOR => 1 },
    { smiles => 'BrC(C(C(CC(=O)O)C(CC[N+](=O)[O-])[PH4])[PH4])C', iupac => '5-bromo-3-[3-nitro-1-(λ5-phosphanyl)propyl]-4-(λ5-phosphanyl)hexanoic acid', AUTHOR => 1 },
    { smiles => 'C(=O)(O)C1=CC(=C(OC2=CC(=C(C(=O)O)C=C2)[PH4])C=C1)P', iupac => '4-(4-carboxy-2-phosphanylphenoxy)-2-(λ5-phosphanyl)benzoic acid', AUTHOR => 1 },
    { smiles => 'C(=O)(O)C1=CC(=C(OC2=CC(=C(C(=O)O)C=C2)[SH5])C=C1)[SH3]', iupac => '4-[4-carboxy-2-(λ4-sulfanyl)phenoxy]-2-(λ6-sulfanyl)benzoic acid', AUTHOR => 1 },
    { smiles => 'C(=O)(O)C1=CC(=C(OC2=CC(=C(C(=O)O)C=C2)[SH5])C=C1)P', iupac => '4-(4-carboxy-2-phosphanylphenoxy)-2-(λ6-sulfanyl)benzoic acid', AUTHOR => 1 },
    { smiles => 'N(C1C=C(C=C(C1)C)SC1=C(C(CC=C1)N([2H])[2H])C)([2H])[2H]', iupac => '3-{[3-(2H2)amino-5-methylcyclohexa-1,5-dien-1-yl]sulfanyl}-2-methylcyclohexa-2,4-dien-1-(2H2)amine', AUTHOR => 1 },
    { smiles => '[81Br]C(C(CC(=O)O)CC(C)[81Br])CC', iupac => '4-(81Br)bromo-3-[2-(81Br)bromopropyl]hexanoic acid', AUTHOR => 1 },
    { smiles => 'ClC(C(C)O)C(C(C(CC(C)O)C)CCC(C)O)C', iupac => '3-chloro-5-(3-hydroxybutyl)-4,6-dimethylnonane-2,8-diol' },

    # From BBv3 P-45.2.3
    { smiles => 'ClC=1C=NC2=CC(=CC=C2C1[N+](=O)[O-])SC1=CC=C2C(=C(C=NC2=C1)[N+](=O)[O-])Cl', iupac => '3-chloro-7-[(4-chloro-3-nitroquinolin-7-yl)sulfanyl]-4-nitroquinoline' },
    { smiles => 'BrC1=C(NC2=C(C=C(C=C2)Br)Cl)C=CC(=C1)Cl', iupac => '2-bromo-N-(4-bromo-2-chlorophenyl)-4-chloroaniline' },
    { smiles => 'C(C)C1=C(C=CC2=CC=C(C=C12)OC1=CC2=C(C(=CC=C2C=C1)CC)CCC)CCC', iupac => '1-ethyl-7-[(7-ethyl-8-propylnaphthalen-2-yl)oxy]-2-propylnaphthalene', AUTHOR => 1 },
    { smiles => 'BrC(C(CCC(=O)O)C(C(C)Br)F)C(C)F', iupac => '5-bromo-4-(2-bromo-1-fluoropropyl)-6-fluoroheptanoic acid' },
    { smiles => 'BrC(C(C(=O)O)C(CBr)O)CO', iupac => '3-bromo-2-(2-bromo-1-hydroxyethyl)-4-hydroxybutanoic acid', AUTHOR => 1 },
    { smiles => 'BrCCC(CO)CCCl', iupac => '2-(2-bromoethyl)-4-chlorobutan-1-ol' },
    { smiles => 'BrC(CCCC(CCl)CBr)Cl', iupac => '1-bromo-5-(bromomethyl)-1,6-dichlorohexane' },
    { smiles => 'C(C)C(C(CC=CC=C)C(C)C(C=C)CC)C(C=C)C', iupac => '7-ethyl-6-(3-ethylpent-4-en-2-yl)-8-methyldeca-1,3,9-triene' },
    { smiles => 'CNC(CCC1=CC(=C(C=C1)C)CCC(NCCC)=O)=O', iupac => 'N-methyl-3-{4-methyl-3-[3-oxo-3-(propylamino)propyl]phenyl}propanamide', AUTHOR => 1 },
    { smiles => 'BrC(C([PH4])C(CC(=O)O)C(C(C)Cl)[PH4])C', iupac => '3-[2-bromo-1-(λ5-phosphanyl)propyl]-5-chloro-4-(λ5-phosphanyl)hexanoic acid', AUTHOR => 1 },
    { smiles => '[81Br]C(C(CC(=O)O)C(C(C)Br)[81Br])C(C)Cl', iupac => '4-(81Br)bromo-3-[1-(81Br)bromo-2-bromopropyl]-5-chlorohexanoic acid', AUTHOR => 1 },
    { smiles => 'C(CCC)C(CC(CC)C)C(CC(CC)CC)CCC', iupac => '5-butyl-8-ethyl-3-methyl-6-propyldecane', AUTHOR => 1 },
    { smiles => 'C(C)C1=CC=C(C2=CC(=CC=C12)[Se]C1=CC2=C(C=CC(=C2C=C1)CCC)CC)CCC', iupac => '1-ethyl-6-[(8-ethyl-5-propylnaphthalen-2-yl)selanyl]-4-propylnaphthalene', AUTHOR => 1 },
    { smiles => 'BrC(CCC(C(CCl)Br)C(CBr)I)Cl', iupac => '1,5-dibromo-4-(2-bromo-1-iodoethyl)-1,6-dichlorohexane', AUTHOR => 1 },
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
