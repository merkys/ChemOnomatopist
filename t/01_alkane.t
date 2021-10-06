#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Graph::Undirected;
use Test::More;

my %SMILES_cases = (
    'CCCCC'  => 'pentane',
    'CC(C)C' => '2-methylpropane', # FIXME: 'methylpropane'
    'C1CCC1' => 'cyclobutane',
);

plan tests => 5 + scalar keys %SMILES_cases;

my $g = Graph::Undirected->new;
$g->add_path( 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H' );
$g->add_path( 'D', 'DB', 'DC', );
$g->add_edge( 'F', 'FB' );

is( ChemOnomatopist::get_name( $g ), '5-ethyl-3-methyloctane' );

$g->add_edge( 'F', 'FC' );

is( ChemOnomatopist::get_name( $g ), '3,3-dimethyl-5-ethyloctane' );

$g->add_edge( 'H', 'I' );

is( ChemOnomatopist::get_name( $g ), '4,4-dimethyl-6-ethylnonane' );

$g->add_edge( 'H', 'HB' );

is( ChemOnomatopist::get_name( $g ), '6-ethyl-2,4,4-trimethylnonane' );

$g->add_edge( 'DB', 'DBB' );

is( ChemOnomatopist::get_name( $g ), '6-1-methylethyl-2,4,4-trimethylnonane' );


is( ChemOnomatopist::get_name( 'CCCCC' ), 'pentane' );

is( ChemOnomatopist::get_name( 'CC(C)C' ), '2-methylpropane' );
# Tests with SMILES input

for my $case (sort keys %SMILES_cases) {
    is( ChemOnomatopist::get_name( $case ), $SMILES_cases{$case} );
}
