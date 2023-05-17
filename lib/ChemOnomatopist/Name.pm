package ChemOnomatopist::Name;

use strict;
use warnings;

# ABSTRACT: Chemical name
# VERSION

use overload '.=' => \&append;
use overload '""' => sub { return $_[0]->{name} };

sub new
{
    my( $class, $name ) = @_;
    $name = '' unless $name;
    return bless { name => $name }, $class;
}

sub append($)
{
    my( $self, $string ) = @_;
    $self->{name} .= $string;
}

sub append_multiplier($)
{
    my( $self, $string ) = @_;
    $self->{starts_with_multiplier} = 1 unless $self->{name};
    $self->append( $string );
}

1;
