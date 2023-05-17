package ChemOnomatopist::Name;

use strict;
use warnings;

# ABSTRACT: Chemical name
# VERSION

use ChemOnomatopist;

use overload '.='  => \&append;
use overload '""'  => sub { return $_[0]->{name} };
use overload 'eq'  => sub { return  "$_[0]" eq  "$_[1]" };
use overload 'cmp' => sub { return ("$_[0]" cmp "$_[1]") * ($_[2] ? -1 : 1) };

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
    return $self;
}

sub append_multiplier($)
{
    my( $self, $string ) = @_;
    $self->{starts_with_multiplier} = 1 unless $self->{name};
    return $self->append( $string );
}

sub bracket()
{
    my( $self ) = @_;
    $self->{name} = ChemOnomatopist::bracket( $self->{name} );
}

sub starts_with_multiplier()
{
    my( $self ) = @_;
    return exists $self->{starts_with_multiplier};
}

1;
