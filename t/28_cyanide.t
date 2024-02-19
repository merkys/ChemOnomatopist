#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-66.5.1.1.1
    { smiles => 'C(CCCCC)#N', iupac => 'hexanenitrile', AUTHOR => 1 },
    { smiles => 'C(CCCC#N)#N', iupac => 'pentanedinitrile', AUTHOR => 1 },

    { smiles => 'C(CCC)(C#N)(C#N)C#N', iupac => 'butane-1,1,1-tricarbonitrile' }, # BBv2 P-66.5.1.1.2

    # From BBv2 P-66.5.1.1.3
    { smiles => '[SiH3]C#N', iupac => 'silanecarbonitrile' },
    { smiles => 'C1(CCCCC1)C#N', iupac => 'cyclohexanecarbonitrile' },
    { smiles => 'N1(CCCCC1)C#N', iupac => 'piperidine-1-carbonitrile' },

    # From BBv2 P-66.5.1.1.4
    { smiles => 'C(#N)C1=CC=C(O1)C(=O)O', iupac => '5-cyanofuran-2-carboxylic acid' },
    { smiles => 'C(#N)CCC(=O)O', iupac => '3-cyanopropanoic acid' },
    { smiles => 'C(#N)CC(CCC#N)CCC#N', iupac => '4-(cyanomethyl)heptanedinitrile', AUTHOR => 1 },
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
