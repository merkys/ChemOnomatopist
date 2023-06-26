#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Test::More;

my @cases = (
    { smiles => 'COCCOCCOCCOCC', iupac => '2,5,8,11-tetraoxatridecane' }, # BBv2 P-12.1
    { smiles => 'C(F)(F)(F)C(F)(F)CO', iupac => '2,2,3,3,3-pentafluoropropan-1-ol' }, # BBv2 P-14.3.4.5

    { smiles => 'ClC(C(Cl)(Cl)Cl)(Cl)Cl', iupac => 'hexachloroethane', AUTHOR => 1 },
    { smiles => 'ClC(=C(Cl)Cl)Cl', iupac => 'tetrachloroethene', AUTHOR => 1 },

    { smiles => 'C1CCCCC1CCOCC', iupac => '(3-oxapentyl)cyclohexane' },
    { smiles => 'C1CCCCC1OCCCC', iupac => 'butoxycyclohexane' }, # PubChem 13299482

    # From BBv2 P-15.4.3.1
    { smiles => 'COCSSCCOCC[Se]C', iupac => '2,8-dioxa-4,5-dithia-11-selenadodecane' },
    { smiles => '[Si]OCS[Si]', iupac => '2-oxa-4-thia-1,5-disilapentane' },

    # From BBv2 P-15.4.3.2.1
    { smiles => 'C[Si]C[Si]C[Si]CSCC', iupac => '8-thia-2,4,6-trisiladecane' },
    { smiles => 'C[Si]C[Si]C[Si]COC',  iupac => '2-oxa-4,6,8-trisilanonane' },

    # From BBv2 P-15.4.3.2.3
    { smiles => 'C[Si]C[Si]C[Si]C[Si]C(=O)O', iupac => '2,4,6,8-tetrasilanonan-1-oic acid', AUTHOR => 1 }, # FIXME
    { smiles => 'C[Si]C[Si]C[Si]C[Si]CCO', iupac => '2,4,6,8-tetrasiladecan-10-ol' },

    { smiles => 'C[SiH2]C[SiH2]C[SiH2]C[SiH2]C=C', iupac => '2,4,6,8-tetrasiladec-9-ene' }, # BBv2 P-15.4.3.2.4

    # From BBv2 P-61.5.2
    { smiles => 'ON(O)OC(C)(C)C(C(O)=O)NC(C)(C)C', iupac => '2-(tert-butylimino)-3-methyl-3-(nitrooxy)butanoic acid', AUTHOR => 1 },

    { smiles => 'FC(C(CC)F)C(CC(CCCC)CC)CCCCCC', iupac => '7-(1,2-difluorobutyl)-5-ethyltridecane' }, # BBv2 P-14.5.2
);

@cases = grep { !exists $_->{AUTHOR} } @cases unless $ENV{AUTHOR_TESTING};
plan skip_all => 'No available cases' unless @cases;
plan tests => scalar @cases;

for my $case (@cases) {
    my $ok;
    eval { $ok = is ChemOnomatopist::get_name( $case->{smiles} ), $case->{iupac} };
    $@ =~ s/\n$// if $@;
    fail $case->{smiles} . ": $@" if $@;
    diag 'test supposed to fail with AUTHOR_TESTING' if $case->{AUTHOR} && $ok;
}
