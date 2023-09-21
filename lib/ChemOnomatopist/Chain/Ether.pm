package ChemOnomatopist::Chain::Ether;

use strict;
use warnings;

# ABSTRACT: Ether chain
# VERSION

use parent ChemOnomatopist::Chain::;

use ChemOnomatopist;
use ChemOnomatopist::Name;
use Scalar::Util qw( blessed );

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
    return 'oxy' if @vertices == 1;

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
        if( $prefixes[0] =~ /ane$/ ) {
            pop @{$prefixes[0]};
            pop @{$prefixes[0]};
        }
        my $name = ChemOnomatopist::Name->new;
        $name->append_locants( $cut_position );
        $name->append( $prefixes[0] );
        $name->append( 'oxy' );
        $name->append( $prefixes[1] );
        return $name;
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
    my $name = $chains[0]->prefix;
    $name =~ s/ane$//;
    $name .= 'oxy';
    $name .= $chains[1]->suffix;
    return $name;
}

1;
