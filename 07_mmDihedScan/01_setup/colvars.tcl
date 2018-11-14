Colvarstrajfrequency    1000
Colvarsrestartfrequency 1000

colvar {
    name dihed182
    width 1.0

    dihedral {
        group1 { atomNumbers 1618 } 
        group2 { atomNumbers 1615 }
        group3 { atomNumbers 1613 }
        group4 { atomNumbers 1629 }
    }
}

harmonic {
    colvars        dihed182
    centers        180
    forceConstant  0.1  ;# kcal/mol/deg^2
}
