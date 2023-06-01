package package ChemOnomatopist::Group::Monocycle;

use strict;
use warnings;

# ABSTRACT: Fused bicyclic group
# VERSION

# From BBv2 P-25.2.1
our @names = (
    [ 'N=CN=CCC', 'NC=CN=CC=', 'pteridine' ],
    [ 'N=NC=CCC', 'CC=CC=CC=', 'cinnoline' ],
    [ 'N=CN=CCC', 'CC=CC=CC=', 'quinazoline' ],
    [ 'N=CC=NCC', 'CC=CC=CC=', 'quinoxaline' ],
    [ 'N=CC=CCC', 'NC=CC=CC=', '1,5-naphthyridine' ], # TODO: There are isomers
    [ 'C=NN=CCC', 'CC=CC=CC=', 'phthalazine' ],
    [ 'N=CC=CCC', 'CC=CC=CC=', 'quinoline' ],
    [ 'C=NC=CCC', 'CC=CC=CC=', 'isoquinoline' ],
    [ 'CC=CCNC',  'C=CC=CCN',  'quinolizine' ],

    [ 'C=NC=NC=', 'N=CNC=C', 'purine' ], # TODO: Special rules apply

    [ 'NN=CCC',  'CC=CC=CC=', 'indazole' ],
    [ 'NC=CCC',  'CC=CC=CC=', 'indole' ],
    [ 'CNC=CC=', 'C=CC=CCC',  'isoindole' ],
    [ 'CC=CNC=', 'C=CC=CCN',  'indolizine', ],
    [ 'CC=CNC',  'C=CC=CN',   '1H-pyrrolizine' ], # TODO: There are isomers
);

1;
