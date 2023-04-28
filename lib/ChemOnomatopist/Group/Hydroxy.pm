package ChemOnomatopist::Group::Hydroxy;

use strict;
use warnings;

# ABSTRACT: Hydroxy group
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Chain::VertexArray;
use ChemOnomatopist::Util::Graph qw(
    graph_longest_paths_from_vertex
    graph_path_between_vertices
);
use Scalar::Util qw( blessed );
use Set::Object;

sub is_oxygen { return 1 }

sub suffix { return 'ol' }

1;
