package ChemOnomatopist::Util::SMILES;

use strict;
use warnings;

# ABSTRACT: SMILES utilities
# VERSION

use parent Exporter::;

our @EXPORT_OK = qw(
    cycle_SMILES
    path_SMILES
);

sub cycle_SMILES
{
    my( $graph, @cycle ) = @_;

    my $SMILES = '';
    for my $i (0..$#cycle) {
        $SMILES .= $cycle[$i]->{symbol};
        next unless $graph->has_edge_attribute( $cycle[$i], $cycle[$i+1], 'bond' );
        $SMILES .=  $graph->get_edge_attribute( $cycle[$i], $cycle[$i+1], 'bond' );
    }
    return $SMILES;
}

sub path_SMILES
{
    my( $graph, @path ) = @_;
 
    my $SMILES = '';
    for my $i (0..$#path) {
        $SMILES .= $path[$i]->{symbol};
        next if $i == $#path;
        next unless $graph->has_edge_attribute( $path[$i], $path[$i+1], 'bond' );
        $SMILES .=  $graph->get_edge_attribute( $path[$i], $path[$i+1], 'bond' );
    }
    return $SMILES;
}

1;
