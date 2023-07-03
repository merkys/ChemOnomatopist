package ChemOnomatopist::Name;

use strict;
use warnings;

# ABSTRACT: Chemical name
# VERSION

use overload '.='  => \&append;
use overload '""'  => sub { return join '', @{$_[0]->{name}} };
use overload 'eq'  => sub { return  "$_[0]" eq  "$_[1]" };
use overload 'cmp' => sub { return ("$_[0]" cmp "$_[1]") * ($_[2] ? -1 : 1) };

use ChemOnomatopist::Name::Part::Multiplier;
use ChemOnomatopist::Name::Part::Stem;
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $name ) = @_;
    my @name_parts;
    if( defined $name && blessed $name ) {
        push @name_parts, $name;
    } elsif( defined $name && $name ne '' ) {
        push @name_parts, $name;
    }
    return bless { name => \@name_parts }, $class;
}

# TODO: Implement vowel elision as written in BBv2 P-16.7
sub append($)
{
    my( $self, $string ) = @_;
    $self->{name}[-1] =~ s/a$// if $string =~ /^a/ && @{$_[0]->{name}};

    # If names are combined and the second one starts with a number, a separator is added.
    if( @{$_[0]->{name}} && blessed $string && $string->isa( ChemOnomatopist::Name:: ) && $string =~ /^\d/ ) {
        push @{$self->{name}}, '-';
    }
    push @{$self->{name}}, blessed $string && $string->isa( ChemOnomatopist::Name:: ) ? @{$string->{name}} : $string;

    # Inherit locant
    if( blessed $string && $string->isa( ChemOnomatopist::Name:: ) ) {
        $self->{has_locant} = 1 if $string->has_locant;
        $self->{has_substituent_locant} = 1 if $string->has_substituent_locant;
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
    $self->append( 'a' ) if $self->ends_with_stem && @locants == 2;
    $self->append( '-' ) if @{$_[0]->{name}};
    return $self->append( join( ',', @locants ) . '-' );
}

sub append_multiplier($)
{
    my( $self, $string ) = @_;
    return $self if $string eq '';

    $self->append( ChemOnomatopist::Name::Part::Multiplier->new( $string ) );
    return $self;
}

sub append_stem($)
{
    my( $self, $stem ) = @_;
    $self->append( ChemOnomatopist::Name::Part::Stem->new( $stem ) );
    return $self;
}

sub append_substituent_locant($)
{
    my( $self, $locant ) = @_;
    $self->{has_locant} = 1;
    $self->{has_substituent_locant} = 1;
    $self->append( '-' . $locant . '-' );
    $self->{name} = [ 'tert-but' ] if $self eq '2-methylpropan-2-';
    return $self;
}

sub append_suffix($)
{
    my( $self, $suffix ) = @_;
    if( @{$_[0]->{name}} ) {
        if( $suffix =~ /^[aeiouy]/ ) {
            $self->{name}[-3] =~ s/e$// if $self =~ /e-[0-9,]+-$/; # BBv2 P-16.7.1 (a)
            $self->{name}[-1] =~ s/e$// if $self =~ /e$/;
        }
        $self->{name}[-1] =~ s/a$// if $self->ends_with_multiplier && $suffix =~ /^[ao]/; # BBv2 P-16.7.1 (b)
    }
    return $self->append( $suffix );
}

# FIXME: Implement according to BBv2 P-16.5.4: {[({[( )]})]}
sub bracket()
{
    my( $self ) = @_;

    if( $self =~ /\{/ ) {
        unshift @{$self->{name}}, '(';
        push    @{$self->{name}}, ')';
    } elsif( $self =~ /\[/ ) {
        unshift @{$self->{name}}, '{';
        push    @{$self->{name}}, '}';
    } elsif( $self =~ /\(/ ) {
        unshift @{$self->{name}}, '[';
        push    @{$self->{name}}, ']';
    } else {
        unshift @{$self->{name}}, '(';
        push    @{$self->{name}}, ')';
    }

    return $self;
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

sub is_enclosed()
{
    my( $self ) = @_;
    return '' unless @{$self->{name}};
    return $self->{name}[0] =~ /^[\(\[\{]/ && $self->{name}[-1] =~ /[\)\]\}]$/;
}

# FIXME: Incomplete, untested and unused
# 0 if simple, 1 if compound
sub level()
{
    my( $self ) = @_;
    return 0 + ((grep { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Stem:: ) } @{$self->{name}}) > 1);
}

sub starts_with_multiplier()
{
    my( $self ) = @_;
    return @{$self->{name}} &&
           blessed $self->{name}[0] &&
           $self->{name}[0]->isa( ChemOnomatopist::Name::Part::Multiplier:: );
}

sub ends_with_multiplier()
{
    my( $self ) = @_;
    return @{$self->{name}} &&
           blessed $self->{name}[-1] &&
           $self->{name}[-1]->isa( ChemOnomatopist::Name::Part::Multiplier:: );
}

sub ends_with_stem()
{
    my( $self ) = @_;
    return @{$self->{name}} &&
           blessed $self->{name}[-1] &&
           $self->{name}[-1]->isa( ChemOnomatopist::Name::Part::Stem:: );
}

1;
