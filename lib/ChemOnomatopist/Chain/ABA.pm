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
    push @{$self->{vertices}}, @_; # FIXME: Very primitive, need to see which side to push
}

sub needs_heteroatom_locants() { '' }

sub prefix() { $elements{ChemOnomatopist::element( $_[0]->{vertices}[1] )}->{prefix} . 'ne' }
sub suffix() { $elements{ChemOnomatopist::element( $_[0]->{vertices}[1] )}->{prefix} . 'ne' }

1;
