package ChemOnomatopist::Chain::Ether;

# ABSTRACT: Ether chain
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Chain::;

use ChemOnomatopist;
use ChemOnomatopist::Group::Ether;
use ChemOnomatopist::Name;
use List::Util qw( first );
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $graph, $parent, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph };
    $self->{parent} = $parent if $parent;
    return bless $self, $class;
}

sub ether_position()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    return first { blessed $vertices[$_] &&
                   $vertices[$_]->isa( ChemOnomatopist::Group::Ether:: ) }
                 0..$#vertices;
}

sub heteroatom_positions()
{
    my( $self ) = @_;

    return @{$self->{heteroatom_positions}} if $self->{heteroatom_positions};

    my @vertices = $self->vertices;
    my $ether_position = $self->ether_position;
    my @heteroatom_positions;
    for (0..$#vertices) {
        next if $_ == $ether_position;
        next if ChemOnomatopist::is_element( $vertices[$_], 'C' );
        push @heteroatom_positions, $_;
    }

    $self->{heteroatom_positions} = \@heteroatom_positions;
    return @heteroatom_positions;
}

sub needs_heteroatom_locants() { 1 }
sub needs_heteroatom_names() { 1 }
sub needs_substituent_locants() { '' }

sub prefix()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    return 'oxy' if @vertices == 1;

    my $cut_position = $self->ether_position;
    if( $cut_position ) {
        my @chains = ( ChemOnomatopist::Chain->new( $self->graph,
                                                    $self->parent,
                                                    reverse @vertices[0..$cut_position-1] ),
                       ChemOnomatopist::Chain->new( $self->graph,
                                                    $self->parent,
                                                    @vertices[$cut_position+1..$#vertices] ) );
        my @prefixes = map { $_->pop_yl } map { $_->prefix } @chains;
        my $name = ChemOnomatopist::Name->new;
        $name->append_locants( $cut_position ) if $chains[0]->needs_substituent_locants;
        $name->append( $prefixes[1] );
        $name->append( 'oxy' );
        $name->append( $prefixes[0] );
        return $name;
    } else {
        my $chain = ChemOnomatopist::Chain->new( $self->graph, @vertices );
        my $name = $chain->prefix->pop_yl;
        return $name . 'oxy';
    }
}

sub suffix()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    my $cut_position = $self->ether_position;

    my @chains = ( ChemOnomatopist::Chain->new( $self->graph,
                                                $self->parent,
                                                reverse @vertices[0..$cut_position-1] ),
                   ChemOnomatopist::Chain->new( $self->graph,
                                                $self->parent,
                                                @vertices[$cut_position+1..$#vertices] ) );
    @chains = reverse @chains if $chains[0]->length > $chains[1]->length;
    my $name = $chains[0]->prefix->pop_yl;
    $name .= 'oxy';
    $name .= $chains[1]->suffix->pop_yl;
    return $name;
}

1;
