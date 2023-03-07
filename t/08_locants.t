#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use Graph::Undirected;
use Test::More;

plan tests => 5;

my $graph;
my @paths;
my @chain;

$graph = Graph::Undirected->new;
$graph->add_path( 1, 0, 11..17 );
$graph->add_path( 1, 22..27 );
$graph->add_path( 1, 32..37 );
$graph->add_edge( 11, 18 );
$graph->add_edge( 25, 28 );
$graph->add_edge( 32, 38 );

@paths = ChemOnomatopist::rule_lowest_numbered_locants_new( $graph, [ [ 0, 11..17 ] ],
                                                                    [ [ 0, 1, 22..27 ],
                                                                      [ 0, 1, 32..37 ] ] );
is scalar( @paths ), 2;
is join( ',', @{$paths[0]} ), '0,1,22,23,24,25,26,27';
is join( ',', @{$paths[1]} ), '0,11,12,13,14,15,16,17';

@chain = ChemOnomatopist::select_main_chain_new( $graph );
is join( ',', @chain ), '27,26,25,24,23,22,1,0,11,12,13,14,15,16,17';

# Figure 7 from Urbonaitė, 2022.
# In the image, however, one branch is held as having priority over another, while in fact they are equal.
$graph = Graph::Undirected->new;
$graph->add_path( 'A'..'G' );
$graph->add_path( 'B', 'H' );
$graph->add_path( 'D', 'I'..'K' );
$graph->add_path( 'J', 'L' );
$graph->add_path( 'E', 'M' );

@paths = ChemOnomatopist::rule_lowest_numbered_locants_new( $graph, [ [ reverse 'A'..'D' ] ],
                                                                    [ [ 'D'..'G' ] ],
                                                                    [ [ 'D', 'I'..'K' ] ] );
is scalar( @paths ), 0;
