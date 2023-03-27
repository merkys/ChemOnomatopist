#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist;
use ChemOnomatopist::ChainHalf;
use Graph::Undirected;
use Test::More;

sub chain($@)
{
    my( $graph, @vertices ) = @_;
    return ChemOnomatopist::ChainHalf->new( $graph, $vertices[1], @vertices );
}

plan tests => 4;

my $graph;
my @paths;

$graph = Graph::Undirected->new;
$graph->add_path( 0, 11..15 );
$graph->add_path( 0, 21..25 );
$graph->add_path( 0, 31..35 );
$graph->add_path( 0, 41..45 );
$graph->add_path( 12, 16 );
$graph->add_path( 22, 26, 27 );

@paths = ChemOnomatopist::rule_most_carbon_in_side_chains_new( $graph, chain( $graph, 0, 11..15 ),
                                                                       chain( $graph, 0, 21..25 ),
                                                                       chain( $graph, 0, 31..35 ),
                                                                       chain( $graph, 0, 41..45 ) );
is scalar( @paths ), 2;

is join( ';', map { join ',', @$_ } @{$paths[0]} ), '0,11,12,13,14,15';
is join( ';', map { join ',', @$_ } @{$paths[1]} ), '0,21,22,23,24,25';

@paths = ChemOnomatopist::rule_most_carbon_in_side_chains_new( $graph, chain( $graph, 0, 11..15 ),
                                                                       chain( $graph, 0, 31..35 ),
                                                                       chain( $graph, 0, 41..45 ) );
is scalar( @paths ), 3;
