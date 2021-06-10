#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Graph;

use Test::More tests => 5;

my $g = Graph->new( undirected => 1 );
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
