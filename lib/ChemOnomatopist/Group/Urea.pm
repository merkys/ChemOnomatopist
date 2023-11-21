package ChemOnomatopist::Group::Urea;

# ABSTRACT: Urea group
# VERSION

use strict;
use warnings;

use parent ChemOnomatopist::Group::, ChemOnomatopist::Chain::;

sub new
{
    my( $class, $graph, @vertices ) = @_;
    return bless { graph => $graph, vertices => \@vertices }, $class;
}

sub needs_substituent_locants() { $_[0]->number_of_branches > 1 }

sub locants(@)
{
    my $self = shift;
    return map { $_ ? 'N' . '\'' x ($_ - 1) : 1 } @_;
}

sub prefix { return 'urea' }
sub suffix { return 'urea' }

1;
