package ChemOnomatopist::Group::Ester;

# ABSTRACT: Ester group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist;
use ChemOnomatopist::Util qw( copy );
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $hydroxylic, $acid ) = @_;
    return bless { hydroxylic => $hydroxylic, acid => $acid }, $class;
}

sub element { 'C' }

sub name
{
    my( $class, $graph ) = @_;

    # TODO: Assume monoester
    my( $ester ) = grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Ester:: ) }
                        $graph->vertices;
    $graph = copy $graph;
    $graph->delete_edge( $ester, $ester->{hydroxylic} );

    my $hydroxylic_part = ChemOnomatopist::get_sidechain_name( $graph, undef, $ester->{hydroxylic} );
    my $acid_part       = ChemOnomatopist::get_sidechain_name( $graph, undef, $ester );

    $acid_part =~ s/yl$/anoate/;

    return "$hydroxylic_part $acid_part";
}

sub suffix() { 'carboxylate' }

1;
