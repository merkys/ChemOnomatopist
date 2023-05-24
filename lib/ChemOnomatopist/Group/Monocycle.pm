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

sub prefix
{
    my( $self ) = @_;
    my $chain = ChemOnomatopist::Chain::Circular->new( $self->{graph}, @{$self->{vertices}} );
    my $name = $chain->name;
    $name =~ s/ane$/yl/;
    return $name;
}

1;
