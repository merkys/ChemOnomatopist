package ChemOnomatopist::Name;

use strict;
use warnings;

# ABSTRACT: Chemical name
# VERSION

use overload '.='  => \&append;
use overload '""'  => sub { return $_[0]->{name} };
use overload 'eq'  => sub { return  "$_[0]" eq  "$_[1]" };
use overload 'cmp' => sub { return ("$_[0]" cmp "$_[1]") * ($_[2] ? -1 : 1) };

use Scalar::Util qw( blessed );

sub new
{
    my( $class, $name ) = @_;
    $name = '' unless $name;
    return bless { name => $name }, $class;
}

# TODO: Implement vowel elision as written in BBv2 P-16.7
sub append($)
{
    my( $self, $string ) = @_;
    $self->{name} =~ s/a$// if $string =~ /^a/;

    # If names are combined and the second one starts with a number, a separator is added.
    if( $self->{name} ne '' && blessed $string && $string->isa( ChemOnomatopist::Name:: ) && $string =~ /^\d/ ) {
        $self->{name} .= '-';
    }
    $self->{name} .= $string;

    delete $self->{ends_with_multiplier};
    delete $self->{ends_with_stem};

    # Inherit locant
    if( blessed $string && $string->isa( ChemOnomatopist::Name:: ) ) {
        $self->{has_locant} = 1 if $string->has_locant;
        $self->{has_substituent_locant} = 1 if $string->has_substituent_locant;
        $self->{ends_with_multiplier} = 1 if $string->ends_with_multiplier;
        $self->{ends_with_stem} = 1 if $string->{ends_with_stem};
    }

    return $self;
}

sub append_element($)
{
    my( $self, $element ) = @_;
    return $self->append( $element );
}

sub append_locants
{
    my( $self, @locants ) = @_;
    $self->{has_locant} = 1;
    $self->append( 'a' ) if $self->{ends_with_stem} && @locants == 2;
    $self->append( '-' ) if $self->{name};
    return $self->append( join( ',', @locants ) . '-' );
}

sub append_multiplier($)
{
    my( $self, $string ) = @_;
    return $self if $string eq '';

    $self->{starts_with_multiplier} = 1 unless $self->{name};
    $self->append( $string );
    $self->{ends_with_multiplier} = 1;
    return $self;
}

sub append_stem($)
{
    my( $self, $stem ) = @_;
    $self->append( $stem );
    $self->{ends_with_stem} = 1;
    return $self;
}

sub append_substituent_locant($)
{
    my( $self, $locant ) = @_;
    $self->{has_locant} = 1;
    $self->{has_substituent_locant} = 1;
    $self->append( '-' . $locant . '-' );
    $self->{name} = 'tert-but' if $self->{name} eq '2-methylpropan-2-';
    return $self;
}

sub append_suffix($)
{
    my( $self, $suffix ) = @_;
    $self->{name} =~ s/e(-[0-9,]+-|)$/$1/ if $suffix =~ /^[aeiouy]/;              # BBv2 P-16.7.1 (a)
    $self->{name} =~ s/a$// if $self->ends_with_multiplier && $suffix =~ /^[ao]/; # BBv2 P-16.7.1 (b)
    return $self->append( $suffix );
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

sub ends_with_multiplier()
{
    my( $self ) = @_;
    return exists $self->{ends_with_multiplier};
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
