

# ====================================================================================
#
# Usage: vmdt -e dihedPercent.tcl -args /work/cluster/limvt/hv1/04_fep/1_tautomer/17041_19/c1_F150A_F182A/
#
# Calculate the percentage of each trajectory for 2GBI's C-N-C-N (see below)
# dihedral angle being in '180' (clashy with imid H and one of the guan H).
#
# angA: abs(phi) <= 70
# angB: transition region, 90+-20 and 270+-20. Not defining explicitly because neg. angles.
# angC: abs(phi)-180 <= 70 (so 110-250)
#
# Sort results: sort -t ' ' -k 1 -k 2 percentDihed.dat > percentDihed.dat1
#
# ====================================================================================

set verbose 1


# open data output file for recording dihed
if {$verbose} {
    set outFile [open "countDihed.dat" a]
    puts $outFile "# Frame\tDih. Angle\t\tRegion "
}
set summary [open "percentDihed.dat" w]
puts $summary "Direc.\tWin\t\t\t0+-70\tMid\t\t180+-70 "
puts $summary "--------------------------------------------"

set count 0
set angA 0
set angB 0
set angC 0
set home [lindex $argv 0]
cd $home
set psf [glob -dir $home/00_main *psf]


# go through all fwd and rev trajectories
foreach group {FEP_F FEP_R} {

    foreach window [glob $home/$group/lambda_*] {

        set lambda [lindex [split $window /] end]
        if {$verbose} {puts $outFile "# trajectory $group/$lambda "}

        mol new $psf
        set dcd [glob -dir $window *dcd]
        mol addfile $dcd type {dcd} first 0 last -1 step 1 waitfor -1 $count

        set a1 [[atomselect top "resname GBI1 and name C"] get index]  ;# guan C
        set a2 [[atomselect top "resname GBI1 and name N2"] get index] ;# bridge N
        set a3 [[atomselect top "resname GBI1 and name C1"] get index] ;# imid right C
        set a4 [[atomselect top "resname GBI1 and name N3"] get index] ;# unprot imid N
        set all [atomselect top all]
        
        set nframes [molinfo $count get numframes]
        for {set i 0} {$i < $nframes} {incr i} {
            $all frame $i
        
            # measure the dihedral angle
            set phi [measure dihed "$a1 $a2 $a3 $a4" frame $i]

            # classify dihedral angle
            if {$phi < 0} {set phi [expr {360+$phi}]}
            #if {$phi <= 70} {incr angA} elseif {[expr abs([expr {$phi-180}])] <= 70} {incr angC} else {incr angB}
            if {[expr abs([expr {$phi-360}])] <= 70} {
                incr angA; set code A
            } elseif {[expr abs([expr {$phi-180}])] <= 70} {
                incr angC; set code C
            } else {incr angB; set code B}
        
            # print data to output file
            if {$verbose} {puts $outFile "$i\t$phi\t$code"}

            unset phi
        } ;# done with all frames

        puts $summary "$group\t$lambda\t$angA\t\t$angB\t\t$angC"

        mol delete $count
        incr count
        set angA 0
        set angB 0
        set angC 0

    } ;# done with this window

puts $summary "--------------------------------------------"
} ;# done with this group

if {$verbose} {close $outFile}
close $summary
exit
