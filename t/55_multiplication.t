#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv3 P-45.1.1
    { smiles => 'O(C1=CC(=C(C(=O)O)C=C1)Cl)C1=CC(=C(C(=O)O)C=C1)Cl', iupac => '4,4\'-oxybis(2-chlorobenzoic acid)', AUTHOR => 1 },
    { smiles => 'C(=O)(O)C1=CC(=C(OC2=CC(=C(C(=O)O)C=C2)Cl)C=C1)Cl', iupac => '4-(4-carboxy-2-chlorophenoxy)-2-chlorobenzoic acid' },

    # From BBv3 P-45.1.2
    { smiles => 'C1(=CC=CC=C1)C(SCC1=CC=CC=C1)SCC1=CC=CC=C1', iupac => '1,1\'-[(phenylmethylene)bis(sulfanediylmethylene)]dibenzene', AUTHOR => 1 },
    { smiles => 'C(=CC1=CC=C(N)C=C1)(C1=CC=C(N)C=C1)C1=CC=C(N)C=C1', iupac => '4,4\',4\'\'-(ethene-1,1,2-triyl)trianiline', AUTHOR => 1 },
    { smiles => 'P(CP(O)(O)=O)(CP(O)(O)=O)CP(O)(O)=O', iupac => '[phosphanetriyltris(methylene)]tris(phosphonic acid)', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)C(C1=CC=CC=C1)SC(OC(C1=CC=CC=C1)(C1=CC=CC=C1)C1=CC=CC=C1)(C1=CC=CC=C1)C1=CC=CC=C1', iupac => '1,1\',1\'\'-({[(diphenylmethyl)sulfanyl]diphenylmethoxy}methanetriyl)tribenzene', AUTHOR => 1 },
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
