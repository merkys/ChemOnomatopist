package ChemOnomatopist::Name;

use strict;
use warnings;

# ABSTRACT: Chemical name
# VERSION

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

sub append_locant($)
{
    my( $self, $locant ) = @_;
    $self->{has_locant} = 1;
    return $self->append( '-' . $locant . '-' );
}

sub append_multiplier($)
{
    my( $self, $string ) = @_;
    $self->{starts_with_multiplier} = 1 unless $self->{name};
    return $self->append( $string );
}

sub append_substituent_locant($)
{
    my( $self, $locant ) = @_;
    $self->{has_substituent_locant} = 1;
    return $self->append( '-' . $locant . '-' );
}

sub bracket()
{
    my( $self ) = @_;
    $self->{name} = _bracket( $self->{name} );
}

sub has_locant()
{
    my( $self ) = @_;
    return exists $self->{has_locant};
}

sub has_substituent_locant()
{
    my( $self ) = @_;
    return exists $self->{has_substituent_locant};
}

sub starts_with_multiplier()
{
    my( $self ) = @_;
    return exists $self->{starts_with_multiplier};
}

# FIXME: Implement according to BBv2 P-16.5.4: {[({[( )]})]}
sub _bracket
{
    my( $name ) = @_;
    return "($name)" if $name =~ /\{/;
    return "{$name}" if $name =~ /\[/;
    return "[$name]" if $name =~ /\(/;
    return "($name)";
}

1;
