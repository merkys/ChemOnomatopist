package ChemOnomatopist::Group::Hydroxy;

use strict;
use warnings;

# ABSTRACT: Hydroxy group
# VERSION

use parent ChemOnomatopist::Group::;

use ChemOnomatopist::Util::Graph qw( graph_path_between_vertices );
use Scalar::Util qw( blessed );
use Set::Object;

sub is_oxygen { return 1 }

sub get_name
{
    my( $class, $graph ) = @_;

    my @hydroxy = grep { blessed( $_ ) && $_->isa( $class ) } $graph->vertices;
    my @C = map { $graph->neighbours( $_ ) } @hydroxy;

    if( @hydroxy == 1 ) { # The easy case
        my $hydroxy = shift @hydroxy;
        my $C = shift @C;

        $graph->delete_vertex( $hydroxy );

        my $name = ChemOnomatopist::get_sidechain_name( $graph, $C );
        $name =~ s/yl$//;
        $name .= 'an' unless $name =~ /-$/;
        $name = '2-methylpropan-2-' if $name eq 'tert-butan';
        return $name .= 'ol';
    } else {
        my @paths;
        my $max_value;
        for my $i (0..$#C) {
            for my $j (($i+1)..$#C) {
                my @path = graph_path_between_vertices( $graph, $C[$i], $C[$j] );
                my $value = (Set::Object->new( @C ) * Set::Object->new( @path ))->size;
                if(      !defined $max_value || $max_value < $value ) {
                    @paths = ( \@path );
                    $max_value = $value;
                } elsif( $max_value == $value ) {
                    push @paths, \@path;
                }
            }
        }

        die "multiple possibilities, do not know what to do\n" if @paths > 1;
        die "not implemented\n";
    }
}

1;
