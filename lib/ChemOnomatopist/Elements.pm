package ChemOnomatopist::Elements;

use strict;
use warnings;

# Replacement prefixes taken from BBv2 P-15.4.1.1

our %elements = (
    B => {
        prefix => 'bora',
        standard_bonding_number => 3,
    },
    C => {
        prefix => 'carba',
        standard_bonding_number => 4,
    },
    N => {
        prefix => 'aza',
        standard_bonding_number => 3,
    },
    O => {
        prefix => 'oxa',
        standard_bonding_number => 2,
    },
    F => {
        prefix => 'fluora',
        standard_bonding_number => 1,
    },


    Al => {
        prefix => 'alumina',
        standard_bonding_number => 3,
    },
    Si => {
        prefix => 'sila',
        standard_bonding_number => 4,
    },
    P => {
        prefix => 'phospha',
        standard_bonding_number => 3,
    },
    S => {
        prefix => 'thia',
        standard_bonding_number => 2,
    },
    Cl => {
        prefix => 'chlora',
        standard_bonding_number => 1,
    },

    Ga => {
        prefix => 'galla',
        standard_bonding_number => 3,
    },
    Ge => {
        prefix => 'germa',
        standard_bonding_number => 4,
    },
    As => {
        prefix => 'arsa',
        standard_bonding_number => 3,
    },
    Se => {
        prefix => 'selena',
        standard_bonding_number => 2,
    },
    Br => {
        prefix => 'broma',
        standard_bonding_number => 1,
    },

    In => {
        prefix => 'inda',
        standard_bonding_number => 3,
    },
    Sn => {
        prefix => 'stanna',
        standard_bonding_number => 4,
    },
    Sb => {
        prefix => 'stiba',
        standard_bonding_number => 3,
    },
    Te => {
        prefix => 'tellura',
        standard_bonding_number => 2,
    },
    I => {
        prefix => 'ioda',
        standard_bonding_number => 1,
    },

    Tl => {
        prefix => 'thalla',
        standard_bonding_number => 3,
    },
    Pb => {
        prefix => 'plumba',
        standard_bonding_number => 4,
    },
    Bi => {
        prefix => 'bisma',
        standard_bonding_number => 3,
    },
    Po => {
        prefix => 'polona',
        standard_bonding_number => 2,
    },
    At => {
        prefix => 'astata',
        standard_bonding_number => 1,
    },
);

1;
