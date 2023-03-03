#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Graph::Undirected;
use Test::More;

plan tests => 3;

my $graph;

$graph = Graph::Undirected->new;
$graph->add_path( 1, 0, 11..17 );
$graph->add_path( 1, 22..27 );
$graph->add_path( 1, 32..37 );
$graph->add_edge( 11, 18 );
$graph->add_edge( 25, 28 );
$graph->add_edge( 32, 38 );

my @paths = ChemOnomatopist::rule_lowest_numbered_locants_new( $graph, [ [ 0, 11..17 ] ],
                                                                       [ [ 0, 1, 22..27 ],
                                                                         [ 0, 1, 32..37 ] ] );
is scalar( @paths ), 2;
is join( ',', @{$paths[0]} ), '0,1,22,23,24,25,26,27';
is join( ',', @{$paths[1]} ), '0,11,12,13,14,15,16,17';
