package ChemOnomatopist::Group::NoncarbonOxoacid;

# ABSTRACT: Mononuclear noncarbon oxoacid as per BBv3 P-67.1.1.1
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::;

use ChemOnomatopist;
use ChemOnomatopist::Group::Hydroxy;
use ChemOnomatopist::Group::Ketone;

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
    P  => 'phosph',
    Se => 'sele',
    Si => 'silic',
    Sb => 'stib',
    S  => 'sulf',
    Te => 'tellur',
);

my %suffixes = (
    0 => { 1 => 'inous', 2 => 'onous', 3 => 'orous', 4 => 'ic' },
    1 => { 1 => 'inic',  2 => 'onic',  3 => 'oric'  },
    2 => { 1 => 'ic' },
);

sub suffix
{
    my( $self ) = @_;
    my @attachments = @{$self->{attachments}}; # print "@attachments";
    my @hydroxy = grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Hydroxy:: ) } @attachments;
    my @ketones = grep { blessed $_ && $_->isa( ChemOnomatopist::Group::Ketone:: ) }  @attachments;
    my $has_ketone = @attachments != @hydroxy;

    my $name = $elements{$self->element} . $suffixes{scalar @ketones}->{scalar @hydroxy} . ' acid';
    return $name;
}

1;
