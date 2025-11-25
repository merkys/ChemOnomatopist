package ChemOnomatopist::DigraphComparator;

# ABSTRACT: Hierarchical digraph comparator as per BBv3 P-92.1.4
# VERSION

use strict;
use warnings;

use Chemistry::OpenSMILES qw( is_double_bond );
use ChemOnomatopist;
use ChemOnomatopist::Util qw(
    atomic_number
    cmp_arrays
);
use Set::Object qw( set );

sub new
{
    my( $class, $graph, $root, $A, $B ) = @_;
    return bless { graph => $graph, root => $root, A => $A, B => $B,
                   A_paths => [ { vertex => $A, seen => set( $root ), prev => $root } ],
                   B_paths => [ { vertex => $B, seen => set( $root ), prev => $root } ] },
                 $class;
}

sub bump
{
    my $self = shift;

    my $graph = $self->{graph};
    my $root = $self->{root};

    my @A;
    my @B;

    my @A_paths;
    for my $path (@{$self->{A_paths}}) {
        my $vertex = $path->{vertex};
        for my $neighbour ($graph->neighbours( $vertex )) {
            if( is_double_bond( $graph, $vertex, $neighbour ) ) {
                push @A, $neighbour;
            }
            next if $neighbour == $path->{prev};
            push @A, $neighbour;
            next if $path->{seen}->has( $neighbour );
            push @A_paths, { vertex => $neighbour, seen => set( @{$path->{seen}}, $vertex ), prev => $vertex };
        }
    }
    $self->{A_paths} = \@A_paths;

    my @B_paths;
    for my $path (@{$self->{B_paths}}) {
        my $vertex = $path->{vertex};
        for my $neighbour ($graph->neighbours( $vertex )) {
            if( is_double_bond( $graph, $vertex, $neighbour ) ) {
                push @B, $neighbour;
            }
            next if $neighbour == $path->{prev};
            push @B, $neighbour;
            next if $path->{seen}->has( $neighbour );
            push @B_paths, { vertex => $neighbour, seen => set( @{$path->{seen}}, $vertex ), prev => $vertex };
        }
    }
    $self->{B_paths} = \@B_paths;

    if( @A_paths > 10000 || @B_paths > 10000 ) {
        die 'structure is too complex to rearrange its chiral centers' . "\n";
    }

    return \@A, \@B;
}

sub compare
{
    my $self = shift;

    my @A = ( $self->{A} );
    my @B = ( $self->{B} );

    while( @A || @B ) {
        my $cmp;

        $cmp = cmp_arrays( [ reverse sort map { atomic_number( $_ ) } @B ],
                           [ reverse sort map { atomic_number( $_ ) } @A ] );
        return $cmp if $cmp;

        # BBv3 P-92.3: higher atomic numbers appear first
        $cmp = cmp_arrays( [ reverse sort map { exists $_->{isotope} ? $_->{isotope} : atomic_number( $_ ) } @B ],
                           [ reverse sort map { exists $_->{isotope} ? $_->{isotope} : atomic_number( $_ ) } @A ] );
        return $cmp if $cmp;

        my( $A, $B ) = $self->bump;
        @A = @$A;
        @B = @$B;
    }

    return 0;
}

1;
