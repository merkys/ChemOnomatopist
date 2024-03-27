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
