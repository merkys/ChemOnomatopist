package ChemOnomatopist::Elements;

use strict;
use warnings;

# Replacement prefixes taken from BBv2 P-15.4.1.1

our %elements = (
    B => {
        prefix => 'bora',
        standard_bonding_number => 3,
        seniority => 20,
    },
    C => {
        prefix => 'carba',
        standard_bonding_number => 4,
        seniority => 15,
    },
    N => {
        prefix => 'aza',
        standard_bonding_number => 3,
        seniority => 10,
    },
    O => {
        prefix => 'oxa',
        standard_bonding_number => 2,
        seniority => 5,
    },
    F => {
        prefix => 'fluora',
        standard_bonding_number => 1,
        seniority => 0,
    },


    Al => {
        prefix => 'alumina',
        standard_bonding_number => 3,
        seniority => 21,
    },
    Si => {
        prefix => 'sila',
        standard_bonding_number => 4,
        seniority => 16,
    },
    P => {
        prefix => 'phospha',
        standard_bonding_number => 3,
        seniority => 11,
    },
    S => {
        prefix => 'thia',
        standard_bonding_number => 2,
        seniority => 6,
    },
    Cl => {
        prefix => 'chlora',
        standard_bonding_number => 1,
        seniority => 1,
    },

    Ga => {
        prefix => 'galla',
        standard_bonding_number => 3,
        seniority => 22,
    },
    Ge => {
        prefix => 'germa',
        standard_bonding_number => 4,
        seniority => 17,
    },
    As => {
        prefix => 'arsa',
        standard_bonding_number => 3,
        seniority => 12,
    },
    Se => {
        prefix => 'selena',
        standard_bonding_number => 2,
        seniority => 7,
    },
    Br => {
        prefix => 'broma',
        standard_bonding_number => 1,
        seniority => 2,
    },

    In => {
        prefix => 'inda',
        standard_bonding_number => 3,
        seniority => 23,
    },
    Sn => {
        prefix => 'stanna',
        standard_bonding_number => 4,
        seniority => 18,
    },
    Sb => {
        prefix => 'stiba',
        standard_bonding_number => 3,
        seniority => 13,
    },
    Te => {
        prefix => 'tellura',
        standard_bonding_number => 2,
        seniority => 8,
    },
    I => {
        prefix => 'ioda',
        standard_bonding_number => 1,
        seniority => 3,
    },

    Tl => {
        prefix => 'thalla',
        standard_bonding_number => 3,
        seniority => 24,
    },
    Pb => {
        prefix => 'plumba',
        standard_bonding_number => 4,
        seniority => 19,
    },
    Bi => {
        prefix => 'bisma',
        standard_bonding_number => 3,
        seniority => 14,
    },
    Po => {
        prefix => 'polona',
        standard_bonding_number => 2,
        seniority => 9,
    },
    At => {
        prefix => 'astata',
        standard_bonding_number => 1,
        seniority => 4,
    },
);

1;
