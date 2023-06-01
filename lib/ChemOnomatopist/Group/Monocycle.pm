package ChemOnomatopist::Group::Monocycle;

use strict;
use warnings;

# ABSTRACT: Monocyclic group
# VERSION

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::Circular::;

use ChemOnomatopist;
use List::Util qw( all );

sub new
{
    my( $class, $graph, @vertices ) = @_;

    # FIXME: For now we generate all possible traversals of the same cycle.
    #        This is not optimal, some caching could be introduced.
    my @chains;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $graph, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $graph, @vertices );
        push @vertices, shift @vertices;
    }

    my( $chain ) = ChemOnomatopist::filter_chains( @chains );

    return bless { graph => $graph, vertices => [ $chain->vertices ] }, $class;
}

sub needs_heteroatom_names()
{
    my( $self ) = @_;
    return $self->length < 3 || $self->length > 10 || all { $_->{symbol} !~ /^[cC]$/ } $self->vertices;
}

sub prefix()
{
    my( $self ) = @_;
    my $chain = ChemOnomatopist::Chain::Circular->new( $self->{graph}, @{$self->{vertices}} );
    my $name = $chain->name;
    $name =~ s/ane$/yl/;
    return $name;
}

1;
