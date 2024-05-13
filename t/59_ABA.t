#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-21.2.3.1
    { smiles => '[SnH3]O[SnH2]O[SnH2]O[SnH3]', iupac => 'tetrastannoxane' },
    { smiles => '[SiH3]S[SiH3]', iupac => 'disilathiane' },
    { smiles => 'SOS', iupac => 'dithioxane' },
    { smiles => 'P[Se]P', iupac => 'diphosphaselenane' },
    { smiles => '[SiH3]N[SiH3]', iupac => 'N-silylsilanamine', AUTHOR => 1 },

    # From BBv3 P-68.1.1.2.2
    { smiles => 'CB(OB(C)C)C', iupac => 'tetramethyldiboroxane' },
    { smiles => '[AlH2]O[AlH]O[AlH2]', iupac => 'trialuminoxane' },
    { smiles => 'CBNBC', iupac => '1-methyl-N-(methylboranyl)boranamine', AUTHOR => 1 },
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
