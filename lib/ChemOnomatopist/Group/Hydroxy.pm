package ChemOnomatopist::Group::Hydroxy;

use strict;
use warnings;

# ABSTRACT: Hydroxy group
# VERSION

use parent ChemOnomatopist::Group::;

use Scalar::Util qw( blessed );

sub is_oxygen { return 1 }

sub get_mainchain_name
{
    my( $class, $graph ) = @_;

    my( $hydroxy ) = grep { blessed( $_ ) && $_->isa( $class ) } $graph->vertices;
    my( $C ) = $graph->neighbours( $hydroxy );
    $graph->delete_vertex( $hydroxy );

    my $name = ChemOnomatopist::get_sidechain_name( $graph, $C );
    $name =~ s/yl$//;
    $name .= 'an' unless $name =~ /-$/;
    $name = '2-methylpropan-2-' if $name eq 'tert-butan';
    return $name .= 'ol';
}

1;
