package ChemOnomatopist::Util;

use strict;
use warnings;

# ABSTRACT: Generic utilities
# VERSION

use Exporter;
use Graph::Undirected;
use Scalar::Util qw( blessed );

use parent Exporter::;

our @EXPORT_OK = qw(
    copy
);

sub copy($)
{
    my( $object ) = @_;
    die "can only copy Graph now\n" unless blessed $object && $object->isa( Graph::Undirected:: );

    # Graphs have to be copied with the following code as Graph::copy() loses attributes.
    my $graph = $object;
    my $copy = $graph->copy;
    for my $edge ($graph->edges) {
        next unless $graph->has_edge_attributes( @$edge );
        $copy->set_edge_attributes( @$edge, $graph->get_edge_attributes( @$edge ) );
    }
    return $copy;
}

1;
