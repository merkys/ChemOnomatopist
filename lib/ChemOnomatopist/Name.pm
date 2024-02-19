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

use ChemOnomatopist::Name::Part::AlkaneANSuffix;
use ChemOnomatopist::Name::Part::Element;
use ChemOnomatopist::Name::Part::Fusion;
use ChemOnomatopist::Name::Part::Locants;
use ChemOnomatopist::Name::Part::Locants::Substituent;
use ChemOnomatopist::Name::Part::Multiplier;
use ChemOnomatopist::Name::Part::Stem;
use Clone qw( clone );
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
    $self->[-1] =~ s/o$// if $name =~ /^o/ && @$self && $self->[-1] ne 'cyclo';

    # If names are combined and the second one starts with a number, a separator is added.
    if( @$self && blessed $name && $name->isa( ChemOnomatopist::Name:: ) && $name =~ /^\d/ ) {
        push @$self, '-';
    }
    # FIXME: The following needlessly converts $name into string
    $name =~ s/^-// if @$self && $self->[-1] =~ /-$/ && $name =~ /^-/;

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
    return $self unless @locants;

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

    my $name_wo_fusion = clone $self;
    # BBv2 P-16.5.4.1.2: fusion indicators are ignored
    @$name_wo_fusion =
        grep { !blessed $_ || !$_->isa( ChemOnomatopist::Name::Part::Fusion:: ) }
             @$name_wo_fusion;

    if( $name_wo_fusion =~ /\{/ ) {
        unshift @$self, '(';
        push    @$self, ')';
    } elsif( $name_wo_fusion =~ /\[/ ) {
        unshift @$self, '{';
        push    @$self, '}';
    } elsif( $name_wo_fusion =~ /\(/ ) {
        unshift @$self, '[';
        push    @$self, ']';
    } else {
        unshift @$self, '(';
        push    @$self, ')';
    }

    return $self;
}

sub bracket_numeric_locants()
{
    my( $self ) = @_;

    if( $self->starts_with_locant && $self->[0]->is_numeric ) {
        $self->[0]->{value} = '[' . $self->[0]->{value};
        $self->[0]->{value} =~ s/-$/]/;
    }

    for my $i (1..$#$self) {
        next unless blessed $self->[$i];
        next unless $self->[$i]->isa( ChemOnomatopist::Name::Part::Locants:: );
        $self->[$i-1] = '[';
        $self->[$i]->{value} =~ s/-$/]/;
    }
}

sub concatenate()
{
    my( $A, $B ) = @_;
    return clone( $A )->append( $B );
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

    my $nelements = scalar grep { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Element:: ) }
                                @$self;
    my $nstems    = scalar grep { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Stem:: ) }
                                @$self;

    return '' if $nstems >= 2;
    return '' if $nstems && $nelements;
    return  1;
}

# FIXME: Incomplete, untested and unused
# 0 if simple, 1 if compound
sub level()
{
    my( $self ) = @_;
    return 0 + ((grep { blessed $_ && $_->isa( ChemOnomatopist::Name::Part::Stem:: ) } @$self) > 1);
}

sub ends_with_alkane_an_suffix()
{
    my( $self ) = @_;
    return @$self &&
           blessed $self->[-1] &&
           $self->[-1]->isa( ChemOnomatopist::Name::Part::AlkaneANSuffix:: );
}

sub starts_with_locant()
{
    my( $self ) = @_;
    return @$self &&
           blessed $self->[0] &&
           $self->[0]->isa( ChemOnomatopist::Name::Part::Locants:: );
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
