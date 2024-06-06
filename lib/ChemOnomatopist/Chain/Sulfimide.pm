package ChemOnomatopist::Chain::Sulfimide;

# ABSTRACT: Sulfimide chain
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub needs_heteroatom_locants() { '' }
sub needs_heteroatom_names() { '' }

sub locants($@)
{
    my $self = shift;
    map { $_ ? $self->{vertices}[1]->{symbol} : 'N' } @_;
}

sub suffix() { 'sulfanimine' }

1;
