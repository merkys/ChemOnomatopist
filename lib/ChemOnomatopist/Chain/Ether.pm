package ChemOnomatopist::Chain::Ether;

use strict;
use warnings;

# ABSTRACT: Chain of atoms
# VERSION

use ChemOnomatopist;
use ChemOnomatopist::Chain;
use Scalar::Util qw( blessed );

use parent ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $parent, @vertices ) = @_;
    my $self = { vertices => \@vertices, graph => $graph };
    $self->{parent} = $parent if $parent;
    return bless $self, $class;
}

sub needs_heteroatom_locants() { return '' }
sub needs_heteroatom_names() { return '' }

sub prefix()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    my( $cut_position ) = grep { !blessed $vertices[$_] &&
                                 ChemOnomatopist::is_element( $vertices[$_], 'O' ) } 0..$#vertices;
    if( $cut_position ) {
        my @chains = ( ChemOnomatopist::Chain->new( $self->graph,
                                                    $self->parent,
                                                    reverse @vertices[0..$cut_position-1] ),
                       ChemOnomatopist::Chain->new( $self->graph,
                                                    $self->parent,
                                                    @vertices[$cut_position+1..$#vertices] ) );
        my @prefixes = map { $_->prefix } @chains;
        $prefixes[0] =~ s/ane$//;
        return $cut_position . '-' . $prefixes[0] . 'oxy' . $prefixes[1];
    } else {
        my $chain = ChemOnomatopist::Chain->new( $self->graph, @vertices );
        my $name = $chain->prefix;
        $name =~ s/ane$//;
        return $name . 'oxy';
    }
}

sub suffix()
{
    my( $self ) = @_;

    my @vertices = $self->vertices;
    my( $cut_position ) = grep { !blessed $vertices[$_] &&
                                 ChemOnomatopist::is_element( $vertices[$_], 'O' ) } 0..$#vertices;

    my @chains = ( ChemOnomatopist::Chain->new( $self->graph,
                                                $self->parent,
                                                reverse @vertices[0..$cut_position-1] ),
                   ChemOnomatopist::Chain->new( $self->graph,
                                                $self->parent,
                                                @vertices[$cut_position+1..$#vertices] ) );
    @chains = reverse @chains if $chains[0]->length > $chains[1]->length;
    return $chains[0]->prefix . 'oxy' . $chains[1]->suffix;
}

1;
