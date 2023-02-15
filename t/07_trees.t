#!/usr/bin/perl

use strict;
use warnings;

use ChemOnomatopist::Util::Graph qw( tree_number_of_branches );
use Graph::Undirected;
use Test::More;

plan tests => 1;

my $graph = Graph::Undirected->new;
$graph->add_path( 'A'..'E' );
is tree_number_of_branches( $graph, 'A'..'E' ), 0;
