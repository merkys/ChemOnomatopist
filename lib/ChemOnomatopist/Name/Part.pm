package ChemOnomatopist::Name::Part;

use strict;
use warnings;

# ABSTRACT: Semantic part of a chemical name
# VERSION

use ChemOnomatopist::Name;

sub new
{
    my( $class, $value ) = @_;
    return bless { value => $value }, $class;
}

sub to_name()
{
    my( $self ) = @_;
    return ChemOnomatopist::Name->new( $self );
}

1;
