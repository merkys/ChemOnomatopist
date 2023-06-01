package ChemOnomatopist::Group::Monocycle;

use strict;
use warnings;

# ABSTRACT: Monocyclic group
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Chain::Circular;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub candidate_chains()
{
    my( $self ) = @_;

    # FIXME: For now we generate all possible traversals of the same cycle.
    #        This is not optimal, some caching could be introduced.
    my @chains;
    my @vertices = @{$self->{vertices}};
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $self->{graph}, @vertices );
        push @vertices, shift @vertices;
    }
    @vertices = reverse @vertices;
    for (0..$#vertices) {
        push @chains, ChemOnomatopist::Chain::Circular->new( $self->{graph}, @vertices );
        push @vertices, shift @vertices;
    }

    return @chains;
}

sub prefix
{
    my( $self ) = @_;
    my $chain = ChemOnomatopist::Chain::Circular->new( $self->{graph}, @{$self->{vertices}} );
    my $name = $chain->name;
    $name =~ s/ane$/yl/;
    return $name;
}

1;
