#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Graph::Undirected;
use Test::More;

plan tests => 1;

my $graph;
my @paths;

$graph = Graph::Undirected->new;
$graph->add_path( 0, 11..15 );
$graph->add_path( 0, 21..25 );
$graph->add_path( 0, 31..35 );
$graph->add_path( 0, 41..45 );
$graph->add_path( 12, 16 );
$graph->add_path( 22, 26, 27 );

@paths = ChemOnomatopist::rule_most_carbon_in_side_chains_new( $graph, [ [ 0, 11..15 ] ],
                                                                       [ [ 0, 21..25 ] ],
                                                                       [ [ 0, 31..35 ] ],
                                                                       [ [ 0, 41..45 ] ] );
is scalar( @paths ), 2;
