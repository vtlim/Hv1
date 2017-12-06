Colvarstrajfrequency    1000
Colvarsrestartfrequency 1000

colvar {
    name lys
    width 0.1
    lowerboundary 3.0
    upperboundary 4.0
    lowerWallConstant 100.0
    upperWallConstant 100.0

    distance {
        group1 {
            atomnumbers { 2099 }
        }
        group2 {
            atomnumbers { 613 }
        }
    }
}
