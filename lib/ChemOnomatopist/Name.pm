package ChemOnomatopist::Name;

use strict;
use warnings;

# ABSTRACT: Chemical name
# VERSION

use overload '.='  => \&append;
use overload '""'  => sub { return join '', @{$_[0]->{name}} };
use overload 'eq'  => sub { return  "$_[0]" eq  "$_[1]" };
use overload 'cmp' => sub { return ("$_[0]" cmp "$_[1]") * ($_[2] ? -1 : 1) };
use overload '@{}' => sub { return $_[0]->{name} };

use ChemOnomatopist::Name::Part::Element;
use ChemOnomatopist::Name::Part::Locants;
use ChemOnomatopist::Name::Part::Locants::Substituent;
use ChemOnomatopist::Name::Part::Multiplier;
use ChemOnomatopist::Name::Part::Stem;
use List::Util qw( any );
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
    my( $self, $name ) = @_;

    $self->[-1] =~ s/a$// if $name =~ /^a/ && @$self;
    $self->[-1] =~ s/o$// if $name =~ /^o/ && @$self;

    # BBv2 P-16.7.1 (d)
    my $terminal_element;
    if( @$self && blessed $self->[-1] && $self->[-1]->isa( ChemOnomatopist::Name::Part::Element:: ) &&
        "$name" =~ /^[aeiouy]/ ) {
        $self->[-1] =~ s/a$//;
    }

    # If names are combined and the second one starts with a number, a separator is added.
    if( @$self && blessed $name && $name->isa( ChemOnomatopist::Name:: ) && $name =~ /^\d/ ) {
        push @$self, '-';
    }
    push @$self, blessed $name && $name->isa( ChemOnomatopist::Name:: ) ? @$name : $name;

    return $self;
}

sub append_element($)
{
    my( $self, $element ) = @_;
    return $self->append( ChemOnomatopist::Name::Part::Element->new( $element ) );
}

sub append_locants
{
    my( $self, @locants ) = @_;
    if( @$self ) {
        $self->append( ChemOnomatopist::Name::Part::Locants->new( '-' . join( ',', @locants ) . '-' ) );
    } else {
        $self->append( ChemOnomatopist::Name::Part::Locants->new( join( ',', @locants ) . '-' ) );
    }
    return $self;
}

sub append_multiplier($)
{
    my( $self, $string ) = @_;
    return $self if $string eq '';
    return $self->append( ChemOnomatopist::Name::Part::Multiplier->new( $string ) );
}

sub append_stem($)
{
    my( $self, $stem ) = @_;
    return $self->append( ChemOnomatopist::Name::Part::Stem->new( $stem ) );
}

sub append_substituent_locant($)
{
    my( $self, $locant ) = @_;
    $self->append( ChemOnomatopist::Name::Part::Locants::Substituent->new( '-' . $locant . '-' ) );
    $self->{name} = [ 'tert-but' ] if $self eq '2-methylpropan-2-';
    return $self;
}

sub append_suffix($)
{
    my( $self, $suffix ) = @_;
    if( @$self ) {
        if( $suffix =~ /^[aeiouy]/ ) {
            $self->[-2] =~ s/e$// if $self =~ /e-[0-9,]+-$/; # BBv2 P-16.7.1 (a)
            $self->[-1] =~ s/e$// if $self =~ /e$/;
        }
        $self->[-1] =~ s/a$// if $self->ends_with_multiplier && $suffix =~ /^[ao]/; # BBv2 P-16.7.1 (b)
    }
    return $self->append( $suffix );
}

# FIXME: Implement according to BBv2 P-16.5.4: {[({[( )]})]}
sub bracket()
{
    my( $self ) = @_;

    if( $self =~ /\{/ ) {
        unshift @$self, '(';
        push    @$self, ')';
    } elsif( $self =~ /\[/ ) {
        unshift @$self, '{';
        push    @$self, '}';
    } elsif( $self =~ /\(/ ) {
        unshift @$self, '[';
        push    @$self, ']';
    } else {
        unshift @$self, '(';
        push    @$self, ')';
    }

    return $self;
}

sub has_locant()
{
    my( $self ) = @_;
    return any { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Locants:: ) } @$self;
}

sub has_substituent_locant()
{
    my( $self ) = @_;
    return any { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Locants::Substituent:: ) } @$self;
}

sub is_enclosed()
{
    my( $self ) = @_;
    return '' unless @$self;
    return $self->[0] =~ /^[\(\[\{]/ && $self->[-1] =~ /[\)\]\}]$/;
}

sub is_simple()
{
    my( $self ) = @_;
    return (grep { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Stem:: ) } @$self) <= 1;
}

# FIXME: Incomplete, untested and unused
# 0 if simple, 1 if compound
sub level()
{
    my( $self ) = @_;
    return 0 + ((grep { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Stem:: ) } @$self) > 1);
}

sub starts_with_multiplier()
{
    my( $self ) = @_;
    return @$self &&
           blessed $self->[0] &&
           $self->[0]->isa( ChemOnomatopist::Name::Part::Multiplier:: );
}

sub ends_with_multiplier()
{
    my( $self ) = @_;
    return @$self &&
           blessed $self->[-1] &&
           $self->[-1]->isa( ChemOnomatopist::Name::Part::Multiplier:: );
}

sub ends_with_stem()
{
    my( $self ) = @_;
    return @$self &&
           blessed $self->[-1] &&
           $self->[-1]->isa( ChemOnomatopist::Name::Part::Stem:: );
}

1;
