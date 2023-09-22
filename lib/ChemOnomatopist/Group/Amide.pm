package ChemOnomatopist::Group::Amide;

use strict;
use warnings;

# ABSTRACT: Amide group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $C, $amine ) = @_;
    return bless { graph => $graph, C => $C, amine => $amine }, $class;
}

sub vertices
{
    my( $self ) = @_;
    my @vertices = ( $self->{C}, $self->{amine} );
    return @vertices;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ ? $_ : 'N' } @_;
}

sub prefix { return 'amido' } # FIXME: Not sure if really
sub suffix { return 'amide' }
sub multisuffix() { return 'carboxamide' }

1;
