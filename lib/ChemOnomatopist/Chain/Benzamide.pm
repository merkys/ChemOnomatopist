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

sub prefix() { return 'benzamido' }

sub suffix()
{
    my( $self ) = @_;
    return 'benz' if $self->{benzene}->is_benzene;

    my $suffix = $self->{benzene}->suffix;
    $suffix .= '-1-carbox';
    return $suffix;
}

1;
