package ChemOnomatopist::Group::NoncarbonOxoacid;

# ABSTRACT: Mononuclear noncarbon oxoacid as per BBv3 P-67.1.1.1
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist;
use ChemOnomatopist::Group::Hydroxy;

sub new
{
    my( $class, $atom, @attachments ) = @_;
    return bless { attachments => \@attachments, element => ChemOnomatopist::element( $atom ) }, $class;
}

my %elements = (
    As => 'ars',
    N  => 'az',
    B  => 'bor',
    Br => 'brom',
    Cl => 'chor',
    F  => 'flour',
    I  => 'iod',
    Se => 'sele',
    Si => 'silic',
    Sb => 'stib',
    S  => 'sulf',
    Te => 'tellur',
);

my %suffixes = (
    1 => {
        3 => 'oric',
        2 => 'onic',
        1 => 'inic',
        0 => 'inous',
    },
    '' => {
        3 => 'orous',
        2 => 'onous',
    },
);

sub suffix
{
    my( $self ) = @_;
    my @attachments = @{$self->{attachments}};
    my @hydroxy = grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Hydroxy:: ) } @attachments;
    my $has_ketone = @attachments != @hydroxy;

    my $name = $elements{$self->element} . $suffixes{$has_ketone}->{scalar @hydroxy} . ' acid';
    return $name;
}

1;
