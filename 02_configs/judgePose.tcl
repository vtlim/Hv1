

# This script is meant to be a one-frame analysis of the final snapshots in 
# vanilla MD simulations for ligand pose analysis. Distinct from measuring 
# contacts over time, or RMSDs over time. This script focuses on quality > stability.

# Example usage:
# - source file.tcl
# - benzAngle 0

# After loading in list of poses in VMD session:
# - source judgePose.tcl
# - evalAll "somefile.dat" "header inside the file"

# ================= TCL PROCEDURES =================

#  [ benzAngle ] returns angles (degrees, range +-90) between two benzo planes,
#    - one on 2gbi, one on F150
#    - uses coordinates of three atoms specified on each ring
#    - extension ideas: can edit to take in arguments specifying
#      atomselect keywords (e.g. considering other Phe's, ...)
#  [ benzDist ]   returns distance (Angstroms) between two benzo centers
#  [ guanAngle ]  returns angle between NCN plane on guanidino on R211 and 2gbi
#  [ guanDist ]   returns central carbon distance of guan groups (guanidino, R211)
#  [ measDihed ]  returns 2gbi dihed angle. defined by CNCN atoms used in dihed scan
#                 NOTE / TODO: only specified for 2gbi TAUTOMER 1
#  [ rogueLipid ] returns the resname&resid of popc (if any) who's C21 or C31 is
#                 within 10 A of N214 ND2 atom.

# Further work...?
#  [ aspImidDist ]


proc planeAngle {list1 list2} {

    # Parameter: two lists, each having three selections that define a plane
    # Further work - rename vars if you want. leftover from f-->f150,g-->gbi

    set a0 [lindex $list1 0]
    set a1 [lindex $list1 1]
    set a2 [lindex $list1 2]
    set b0 [lindex $list2 0]
    set b1 [lindex $list2 1]
    set b2 [lindex $list2 2]

    # get coordinates of selected A atoms
    set fc1 [lindex [$a0 get {x y z}] 0]
    set fc2 [lindex [$a1 get {x y z}] 0]
    set fc3 [lindex [$a2 get {x y z}] 0]
    
    # compute two vectors of the Phe plane
    set fv1 [vecsub $fc1 $fc2]
    set fv2 [vecsub $fc2 $fc3]
    
    # normal vector of the Phe plane
    set fn [veccross $fv1 $fv2]

    # get coordinates of selected B atoms =====================
    set gc1 [lindex [$b0 get {x y z}] 0]
    set gc2 [lindex [$b1 get {x y z}] 0]
    set gc3 [lindex [$b2 get {x y z}] 0]
    
    # compute two vectors of the 2gbi plane
    set gv1 [vecsub $gc1 $gc2]
    set gv2 [vecsub $gc2 $gc3]
    
    # normal vector of the 2gbi plane
    set gn [veccross $gv1 $gv2]
    
    # find angle bt planes using normal vectors
    # by cos(theta)=A.B / |A||B|
    set vdot [vecdot $gn $fn]
    set v1mag [veclength $gn]
    set v2mag [veclength $fn]
    
    set radn [expr {acos($vdot/($v1mag*$v2mag))}]
    set degs [expr {$radn*180/3.14159265}]
    
    # scale to plus or minus 90
    if {$degs > 90} {set degs [expr {$degs - 180}]}
    return [format %.3f $degs]
}

proc guanAngle {molID} {
    # Parameter: VMD molID
    # Returns: angle in degrees bt the guan planes
    
    # select three atoms connected to main chain =================
    set f1 [atomselect $molID "protein and resid 211 and name NH2"]
    set f2 [atomselect $molID "protein and resid 211 and name CZ"]
    set f3 [atomselect $molID "protein and resid 211 and name NE"]
    
    
    # select three atoms on the imid H side
    set g1 [atomselect $molID "resname GBI1 and name N"]
    set g2 [atomselect $molID "resname GBI1 and name C"]
    set g3 [atomselect $molID "resname GBI1 and name N2"]

    return [format %.3f [planeAngle [list $f1 $f2 $f3] [list $g1 $g2 $g3]]]
}

proc benzAngle {molID} {
    # Parameter: VMD molID
    # Returns: angle in degrees bt the benzo ring planes
    # Reference: https://tinyurl.com/j83qxvx
    # Interpretation: less than |30| is good, i.e., mostly parallel
    
    # select three atoms connected to main chain =================
    set f1 [atomselect $molID "protein and resid 150 and name CD1"]
    set f2 [atomselect $molID "protein and resid 150 and name CG"]
    set f3 [atomselect $molID "protein and resid 150 and name CD2"]
    
    
    # select three atoms on the imid H side
    set g1 [atomselect $molID "resname GBI1 and name C5"]
    set g2 [atomselect $molID "resname GBI1 and name C6"]
    set g3 [atomselect $molID "resname GBI1 and name C7"]

    return [format %.3f [planeAngle [list $f1 $f2 $f3] [list $g1 $g2 $g3]]]

}


proc benzDist {molID} {
    set vec1 [measure center [atomselect $molID "protein and resid 150 and name CG CD2 CE2 CZ CE1 CD1"]]
    set vec2 [measure center [atomselect $molID "resname GBI1 and name C2 C3 C4 C5 C6 C7"]]
    set dist [veclength [vecsub $vec1 $vec2]]
    return [format %.3f $dist]
}


proc guanDist {molID} {
    set vec1 [measure center [atomselect $molID "protein and resid 211 and name CZ"]]
    set vec2 [measure center [atomselect $molID "resname GBI1 and name C"]]
    set dist [veclength [vecsub $vec1 $vec2]]
    return [format %.3f $dist]
}

proc measDihed {molID} {
    set a1 [[atomselect $molID "resname GBI1 and name C"] get index]  ;# guan C
    set a2 [[atomselect $molID "resname GBI1 and name N2"] get index] ;# bridge N
    set a3 [[atomselect $molID "resname GBI1 and name C1"] get index] ;# imid right C
    set a4 [[atomselect $molID "resname GBI1 and name N3"] get index] ;# unprot imid N
    set phi [measure dihed "$a1 $a2 $a3 $a4" molid $molID]
    unset a1 a2 a3 a4
    return [format %.3f $phi]
}

proc rogueLipid {molID} {
    set sel1 [atomselect $molID "name C21 C31"]
    set sel2 [atomselect $molID "protein and resid 214 and name ND2"]
    set lists [measure contacts 10 $sel1 $sel2]
    set lipindlist [lindex $lists 0]
    set asnindlist [lindex $lists 1]
    
    foreach i $lipindlist j $asnindlist {
        set dist [measure bond [list $i $j] molid $molID]
        lappend answer "[lindex [[atomselect $molID "index $i"] get {resname resid name}] 0]"
        lappend answer [format %.3f $dist]
    }
    return $answer

}

proc evalAll {outputfile header} {
    set outFile [open "$outputfile" w]
    puts $outFile "# $header"

    # measure and print benzo angle info
    puts $outFile "\n# molID\t\t\tF150-benzo angle"
    foreach i [molinfo list] {
        set x [benzAngle $i]
        puts $outFile "[molinfo $i get name]\t $x "
    }

    # measure and print benzo dist info
    puts $outFile "\n# molID\t\t\tF150-benzo dist"
    foreach i [molinfo list] {
        set x [benzDist $i]
        puts $outFile "[molinfo $i get name]\t $x "
    }

    # measure and print guan angle info
    puts $outFile "\n# molID\t\t\tR211-guan angle"
    foreach i [molinfo list] {
        set x [guanAngle $i]
        puts $outFile "[molinfo $i get name]\t $x "
    }

    # measure and print guan dist info
    puts $outFile "\n# molID\t\t\tR211-guanC dist"
    foreach i [molinfo list] {
        set x [guanDist $i]
        puts $outFile "[molinfo $i get name]\t $x "
    }

    # measure and print 2gbi dihed angle info
    puts $outFile "\n# molID\t\t\t2gbi dihed angle"
    foreach i [molinfo list] {
        set x [measDihed $i]
        puts $outFile "[molinfo $i get name]\t $x "
    }

    # measure and print potentially high lipid mols
    puts $outFile "\n# molID\t\t\tC21 C31 within 10 A of N214(ND2)"
    foreach i [molinfo list] {
        set x [rogueLipid $i]
        puts $outFile "\n[molinfo $i get name]\n  $x "
    }
    close $outFile
    
}    
