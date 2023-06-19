package ChemOnomatopist::Group::Amine::SecondaryTertiary;

use strict;
use warnings;

# ABSTRACT: Secondary or tertiary amine
# VERSION

use parent ChemOnomatopist::Group::;

sub is_nitrogen() { return 1 }

sub C() {
    my( $self ) = @_;
    return $self;
}

sub prefix() { return 'amino' }
sub suffix() { return 'amine' }

1;
