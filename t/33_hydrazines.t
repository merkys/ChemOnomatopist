#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    # From BBv2 P-68.3.1.2.1
    { smiles => 'CN(N)C', iupac => '1,1-dimethylhydrazine' },
    { smiles => 'C1(=CC=CC=C1)NN', iupac => 'phenylhydrazine' },
    { smiles => 'N(N)CN' => iupac => '1-hydrazinylmethanamine' },
    { smiles => 'N(N)C(=O)O', iupac => 'hydrazinecarboxylic acid', AUTHOR => 1 },
    { smiles => 'FN(N(F)F)F', iupac => 'tetrafluorohydrazine' },
    { smiles => 'N(N)CC#N', iupac => 'hydrazinylacetonitrile', AUTHOR => 1 },

    # From BBv3 P-68.3.1.2.2
    { smiles => 'C(CC)=NN', iupac => 'propylidenehydrazine' },
    { smiles => 'CN(N=C(C)C)C', iupac => '1,1-dimethyl-2-(propan-2-ylidene)hydrazine' },
    { smiles => 'C(C=NNC1=CC=CC=C1)=NNC1=CC=CC=C1', iupac => '1,1\'-ethanediylidenebis(2-phenylhydrazine)', AUTHOR => 1 },
    { smiles => 'CC(C)=NNC1=C(C(=O)O)C=CC=C1', iupac => '2-[(propan-2-ylidene)hydrazinyl]benzoic acid', AUTHOR => 1 },
    { smiles => 'C1(=CC=CC=C1)NN=C1CCC(CC1)C(=O)O', iupac => '4-(phenylhydrazinylidene)cyclohexane-1-carboxylic acid' },

    # From BBv3 P-68.3.1.3.5.1
    { smiles => 'C1(=CC=CC=C1)N=NC(C(C)=O)=NNC1=CC=CC=C1', iupac => '1-(phenyldiazenyl)-1-(phenylhydrazinylidene)propan-2-one' },
    { smiles => 'C1(=CC=CC=C1)C(C(=NNC1=CC=CC=C1)N=NC1=CC=CC=C1)=O', iupac => '1-phenyl-2-(phenyldiazenyl)-2-(phenylhydrazinylidene)ethan-1-one' },
    { smiles => 'C1(=CC=CC=C1)N=NC=NNC(C)=O', iupac => 'N\'-[(phenyldiazenyl)methylidene]acetohydrazide', AUTHOR => 1 },

    { smiles => 'C1(=CC=CC=C1)N=NC(CC(=O)O)=NNC1=CC=CC=C1', iupac => '3-(phenyldiazenyl)-3-(phenylhydrazinylidene)propanoic acid' }, # From BBv3 P-68.3.1.3.5.2
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
