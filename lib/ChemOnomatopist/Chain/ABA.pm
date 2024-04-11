package ChemOnomatopist::Chain::ABA;

# ABSTRACT: a(ba)n chain as per BBv3 P-21.2.3.1
# VERSION

use parent ChemOnomatopist::Chain::;

use ChemOnomatopist::Elements qw( %elements );

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub add
{
    my $self = shift;
    die "cannot extend ABA chain with more than two atoms\n" if @_ > 2;
    my( $first, $last ) = ( $self->{vertices}[0], $self->{vertices}[-1] );
    for (@_) {
        if( $self->graph->has_edge( $first, $_ ) ) {
            unshift @{$self->{vertices}}, $_;
        }
        if( $self->graph->has_edge( $last, $_ ) ) {
            push @{$self->{vertices}}, $_;
        }
    }
}

sub needs_heteroatom_locants() { '' }

sub prefix() { $elements{ChemOnomatopist::element( $_[0]->{vertices}[1] )}->{prefix} . 'ne' }
sub suffix() { $elements{ChemOnomatopist::element( $_[0]->{vertices}[1] )}->{prefix} . 'ne' }

1;
