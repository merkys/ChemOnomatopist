package ChemOnomatopist::Group::SulfonicAcid;

# ABSTRACT: Sulfonic acid group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Name;
use ChemOnomatopist::Util qw( array_frequencies );
use List::Util qw( all any first uniq );
use Scalar::Util qw( blessed );

sub new
{
    my( $class, $element, @attachments ) = @_;
    return bless { attachments => \@attachments, element => $element }, $class;
}

sub hydroxy()
{
    my( $self ) = @_;
    return first { blessed $_ && ( $_->isa( ChemOnomatopist::Group::Hydroxy:: ) ||
                                   $_->isa( ChemOnomatopist::Group::Hydroperoxide:: ) ) }
                 @{$self->{attachments}};
}

my %suffixes = ( N => 'imido', O => '', S => 'thio', Se => 'seleno', Te => 'telluro' );

sub attachments_part()
{
    my( $self ) = @_;

    my @attachments = @{$self->{attachments}};
    my $hydroxy = $self->hydroxy;
    my @non_hydroxy = grep { $_ != $hydroxy } @attachments;
    my @non_hydroxy_elements = map { ChemOnomatopist::Util::element( $_ ) } @non_hydroxy;

    my %elements = array_frequencies @non_hydroxy_elements;
    if( $hydroxy->isa( ChemOnomatopist::Group::Hydroxy:: ) ) {
        $elements{$hydroxy->element}++;
    }

    my @names;
    for (keys %elements) {
        next unless $suffixes{$_};
        my $name = ChemOnomatopist::Name->new;
        $name->append_multiplier( ChemOnomatopist::IUPAC_numerical_multiplier( $elements{$_} ) ) if $elements{$_} > 1;
        $name->append_element( $suffixes{$_} );
        push @names, $name;
    }
    if( $hydroxy->isa( ChemOnomatopist::Group::Hydroperoxide:: ) ) {
        my $suffix = $hydroxy->suffix;
        $suffix->[ 0]{value} =~ s/^-[^\-]+-//;
        $suffix->[-1]{value} =~ s/l$//;
        $suffix->bracket unless $suffix eq 'peroxo'; # non-OO needs brackets
        push @names, $suffix;
    }

    my $name = ChemOnomatopist::Name->new;
    for (sort { _cmp_names( $a, $b ) } @names) {
        $name->[-1]{value} =~ s/o$// if @$name && $_ eq 'imido';
        $name .= $_;
    }
    if( $name =~ /\)$/ ) {
        $name->[-2]{value} .= 'ic';
        $name .= ' ';
    } else {
        $name->[-1]{value} =~ s/(imid)o$/$1/;
        $name .= 'ic ';
    }
    return $name;
}

# From BBv2 P-65.3.0 and Table 4.3
# FIXME: prefix() has to enumerate elements in the attachments
sub element_prefix()
{
    my( $self ) = @_;
    my $suffix =  $suffixes{$self->element};
    $suffix = 'sulfo'    if $self->element eq 'S';
    $suffix = 'selenono' if $self->element eq 'Se';
    return ChemOnomatopist::Name->new( $suffix );
}

sub prefix() { $_[0]->element_prefix }

sub suffix()
{
    my( $self ) = @_;

    my @attachments = @{$self->{attachments}};
    my $hydroxy = first { blessed $_ && ( $_->isa( ChemOnomatopist::Group::Hydroxy:: ) ||
                                          $_->isa( ChemOnomatopist::Group::Hydroperoxide:: ) ) }
                        @attachments;
    my @non_hydroxy = grep { $_ != $hydroxy } @attachments;
    my @non_hydroxy_elements = map { ChemOnomatopist::Util::element( $_ ) } @non_hydroxy;

    if( $hydroxy->isa( ChemOnomatopist::Group::Hydroxy:: ) &&
        $hydroxy->element eq 'O' &&
        all { $_ eq 'O' } @non_hydroxy_elements ) {
        my $name = $self->element_prefix;
        $name->[-1]{value} =~ s/no$//;
        if( $hydroxy->charge ) {
            $name .= 'nate'; # BBv3 P-72.2.2.2.1.1
        } else {
            $name .= 'nic acid';
        }
        return $name;
    }

    my $name = $self->element_prefix;

    my $attachments_part = $self->attachments_part;
    $name .= $attachments_part =~ /^i/ ? 'n' : 'no';
    $name .= $self->attachments_part;

    if( $hydroxy->isa( ChemOnomatopist::Group::Hydroxy:: ) ) {
        # Needed if at least one non-hydroxy element is different (and not N)
        if( any { $_ ne 'N' && $_ ne $hydroxy->element } @non_hydroxy_elements ) {
            $name .= $hydroxy->element . '-';
        }
    } else {
        # Peroxides need explicit elements if:
        # a) elements are different
        # b) elements are OO and there is a non-N and non-O element among non-hydroxy elements
        my @elements = map { ChemOnomatopist::Util::element( $_ ) } @{$hydroxy->{atoms}};
        if( scalar( uniq @elements ) == 2 ||
            ((all { $_ eq 'O' } @elements) && any { $_ ne 'N' && $_ ne 'O' } @non_hydroxy_elements) ) {
            $name .= join( '', @elements ) . '-';
        }
    }

    return $name . 'acid';
}

sub _cmp_names
{
    my( $A, $B ) = @_;

    my @A = @$A;
    my @B = @$B;

    shift @A if $A->starts_with_multiplier;
    shift @B if $B->starts_with_multiplier;

    if( $A->is_enclosed ) {
        shift @A;
        pop @A;
    }
    if( $B->is_enclosed ) {
        shift @B;
        pop @B;
    }

    local $" = '';
    return "@A" cmp "@B";
}

1;
