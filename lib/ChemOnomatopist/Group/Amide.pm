package ChemOnomatopist::Group::Amide;

use strict;
use warnings;

# ABSTRACT: Amide group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $amide, $ketone, $C ) = @_;
    return bless { graph => $graph, amide => $amide, ketone => $ketone, C => $C }, $class;
}

sub vertices
{
    my( $self ) = @_;
    my @vertices = ( $self->{amide}, $self->{ketone}, $self->{C} );
    return @vertices;
}

sub prefix { return 'amido' }
sub suffix { return 'amide' }

1;
