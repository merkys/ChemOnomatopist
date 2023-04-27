package ChemOnomatopist::Group::Carbonyl;

use strict;
use warnings;

# ABSTRACT: Carbonyl group
# VERSION

use parent ChemOnomatopist::Group::;

use Scalar::Util qw( blessed );

sub is_carbon { return 1 }

sub get_name
{
    my( $class, $graph ) = @_;

    my( $ketone ) = grep { blessed( $_ ) && $_->isa( $class ) } $graph->vertices;
    my $name = ChemOnomatopist::get_sidechain_name( $graph, $ketone );
    $name =~ s/yl$//;
    $name .= 'an-1-' unless $name =~ /-$/;
    return $name . 'one';
}

sub suffix { return 'one' };

1;
