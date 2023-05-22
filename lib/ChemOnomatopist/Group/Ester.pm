package ChemOnomatopist::Group::Ester;

use strict;
use warnings;

# ABSTRACT: Ester group
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $hydroxylic, $acid ) = @_;
    return bless { hydroxylic => $hydroxylic, acid => $acid }, $class;
}

1;
