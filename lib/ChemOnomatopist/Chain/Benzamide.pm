package ChemOnomatopist::Chain::Benzamide;

use strict;
use warnings;

# ABSTRACT: Benzamide chain
# VERSION

use parent ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, $amide, $C, $benzene ) = @_;
    $benzene->parent( $C );
    return bless { graph => $graph,
                   benzene => $benzene,
                   vertices => [ $amide, $C, $benzene->vertices ] }, $class;
}

sub locants(@)
{
    my $self = shift;
    return map { $_ > 1 ? $_ - 1 : $_ ? '?' : 'N' } @_;
}

sub suffix() { return 'benz' }

1;
