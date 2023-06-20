package ChemOnomatopist::Group::Amine::SecondaryTertiary;

use strict;
use warnings;

# ABSTRACT: Secondary or tertiary amine
# VERSION

use parent ChemOnomatopist::Group::;

sub new
{
    my( $class, $graph ) = @_;
    return bless { graph => $graph }, $class;
}

sub C()
{
    my( $self ) = @_;
    return $self;
}

sub candidates()
{
    my( $self ) = @_;

    my @chains;
    for ($self->{graph}->neighbours( $self )) {
        my $chain = ChemOnomatopist::select_sidechain( $self->{graph}, $self, $_ );
        unshift @{$chain->{vertices}}, $self; # TODO: Invalidate cache maybe?
        push @chains, $chain;
    }

    return @chains;
}

sub is_nitrogen() { return 1 }

sub prefix() { return 'amino' }
sub suffix() { return 'amine' }

1;
