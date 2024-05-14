package ChemOnomatopist::Chain::ABA;

# ABSTRACT: a(ba)n chain as per BBv3 P-21.2.3.1
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::;

use ChemOnomatopist::Name::Part::Element;
use ChemOnomatopist::Elements qw( %elements );
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub add
{
    my( $self, @atoms ) = @_;
    die "cannot extend ABA chain into more than two sides\n" if @atoms > 2;

    my @chains = map { blessed $_ ? [ $_->vertices ] : [ $_ ] } @atoms;

    my $graph = $self->graph;
    my @vertices = $self->vertices;

    for (@chains) {
        next unless @$_;

        if(      $graph->has_edge( $self->{vertices}[ 0], $_->[ 0] ) ) {
            unshift @vertices, reverse @$_;
        } elsif( $graph->has_edge( $self->{vertices}[ 0], $_->[-1] ) ) {
            unshift @vertices, @$_;
        } elsif( $graph->has_edge( $self->{vertices}[-1], $_->[ 0] ) ) {
            push @vertices, @$_;
        } elsif( $graph->has_edge( $self->{vertices}[-1], $_->[-1] ) ) {
            push @vertices, reverse @$_;
        }
    }

    $self->{vertices} = \@vertices;
}

sub inner_element { ChemOnomatopist::element( $_[0]->{vertices}[1] ) }
sub outer_element { ChemOnomatopist::element( $_[0]->{vertices}[0] ) }

sub heteroatom_positions()
{
    my( $self ) = @_;

    return @{$self->{heteroatom_positions}} if $self->{heteroatom_positions};

    my @vertices = $self->vertices;
    my @heteroatom_positions;
    for (0..$#vertices) {
        next if ChemOnomatopist::element( $vertices[$_] ) eq 'C';
        next if ChemOnomatopist::element( $vertices[$_] ) eq $self->inner_element;
        push @heteroatom_positions, $_;
    }

    $self->{heteroatom_positions} = \@heteroatom_positions;
    return @heteroatom_positions;
}

sub needs_heteroatom_locants() { '' }

sub prefix() { ChemOnomatopist::Name::Part::Element->new( $elements{$_[0]->inner_element}->{prefix} . 'ne' )->to_name }
sub suffix() { &prefix }

1;
