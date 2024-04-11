package ChemOnomatopist::Chain::ABA;

# ABSTRACT: a(ba)n chain as per BBv3 P-21.2.3.1
# VERSION

use parent ChemOnomatopist::Chain::;

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

    my @chain1 = ( $atoms[0] );
    my @chain2 = ( $atoms[1] );

    # Unpack the atoms in chains
    @chain1 = map { blessed $_ ? $_->vertices : $_ } @chain1;
    @chain2 = map { blessed $_ ? $_->vertices : $_ } @chain2;

    my $graph = $self->graph;
    my @vertices = $self->vertices;

    for (\@chain1, \@chain2) {
        next unless @$_;

        if( $graph->has_edge( $self->{vertices}[ 0], $_->[ 0] ) ) {
            unshift @vertices, reverse @$_;
        }
        if( $graph->has_edge( $self->{vertices}[-1], $_->[ 0] ) ) {
            push @vertices, @$_;
        }
        if( $graph->has_edge( $self->{vertices}[ 0], $_->[-1] ) ) {
            unshift @vertices, @$_;
        }
        if( $graph->has_edge( $self->{vertices}[-1], $_->[-1] ) ) {
            push @vertices, reverse @$_;
        }
    }

    $self->{vertices} = \@vertices;
}

sub needs_heteroatom_locants() { '' }

sub prefix() { $elements{ChemOnomatopist::element( $_[0]->{vertices}[1] )}->{prefix} . 'ne' }
sub suffix() { $elements{ChemOnomatopist::element( $_[0]->{vertices}[1] )}->{prefix} . 'ne' }

1;
