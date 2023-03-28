package ChemOnomatopist::Chain;

use strict;
use warnings;

# ABSTRACT: Longest chain in a compound
# VERSION

use List::Util qw( sum0 );

sub AUTOLOAD {
    our $AUTOLOAD;
    my $call = $AUTOLOAD;
    $call =~ s/.*:://;
    return if $call eq 'DESTROY';
    if( $call =~ /^number_/ ) {
        return sum0 map { $_->can( $call )->() } @{$_[0]->{halves}};
    } else {
        return;
    }
}

sub new
{
    my( $class, @halves ) = @_;
    return bless $class, { halves => \@halves };
}

sub locant_positions()
{
    my( $self ) = @_;
    return $self->{halves}[0]->locant_positions_backward +
           $self->{halves}[1]->locant_positions_forward;
}

1;
