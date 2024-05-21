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

my %suffixes = ( N => 'imido', O => '', S => 'thio', Se => 'seleno', Te => 'telluro' );

# From BBv2 P-65.3.0 and Table 4.3
sub prefix() { ChemOnomatopist::Name->new( $_[0]->element eq 'S' ? 'sulfo' : $suffixes{$_[0]->element} ) }
sub suffix()
{
    my( $self ) = @_;

    my @attachments = @{$self->{attachments}};
    my $hydroxy = first { blessed $_ && ( $_->isa( ChemOnomatopist::Group::Hydroxy:: ) ||
                                          $_->isa( ChemOnomatopist::Group::Hydroperoxide:: ) ) }
                        @attachments;
    my @non_hydroxy = grep { $_ != $hydroxy } @attachments;
    my @non_hydroxy_elements = map { ChemOnomatopist::element( $_ ) } @non_hydroxy;
    if( $hydroxy->isa( ChemOnomatopist::Group::Hydroxy:: ) &&
        $hydroxy->element eq 'O' &&
        all { $_ eq 'O' } @non_hydroxy_elements ) {
        my $name = $self->prefix;
        $name->[-1]{value} =~ s/no$//;
        return $name .= 'nic acid';
    }

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
    @names = sort { _cmp_names( $a, $b ) } @names;

    my $name = $self->prefix;
    $name->[-1]{value} .= 'no' unless $name =~ /no$/;
    for (sort { _cmp_names( $a, $b ) } @names) {
        $name->[-1]{value} =~ s/o$// if $_ eq 'imido';
        $name .= $_;
    }
    if( $name =~ /\)$/ ) {
        $name->[-2]{value} .= 'ic';
        $name .= ' ';
    } else {
        $name->[-1]{value} =~ s/(imid)o$/$1/;
        $name .= 'ic ';
    }

    if( $hydroxy->isa( ChemOnomatopist::Group::Hydroxy:: ) ) {
        if( any { $_ ne 'N' && $_ ne $hydroxy->element } @non_hydroxy_elements ) {
            $name .= $hydroxy->element . '-';
        }
    } else {
        my @elements = map { ChemOnomatopist::element( $_ ) } @{$hydroxy->{atoms}};
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
