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

sub needs_heteroatom_locants() { '' }

sub prefix() { $elements{$_[0]->outer_element}->{prefix} . 'ne' }
sub suffix() { $elements{$_[0]->outer_element}->{prefix} . 'ne' }

1;
