
# Purpose: Measure distance of selection1 to selection[2,3,4]...
# Usage: vmd -dispdev none -e calcDist.tcl -args input.psf input.dcd "selection1" "selection2" ["selection3" "selection4"...]
# Notes: If you get error: "measure center: bad weight sum, would cause divide by zero", double check your selection criteria.
#
# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set indcd [lindex $argv 1]
set pre0 [lindex $argv 2]
set pre1 [lindex $argv 3]
set pre2 [lindex $argv 4]
set pre3 [lindex $argv 5]
set pre4 [lindex $argv 6]
set pre5 [lindex $argv 7]

# processing of other arguments
set sel0 [split $pre0 {,}]
set sel1 [split $pre1 {,}]
if { $argc >= 5 } {set sel2 [split $pre2 {,}]}
if { $argc >= 6 } {set sel3 [split $pre3 {,}]}
if { $argc >= 7 } {set sel4 [split $pre4 {,}]}
if { $argc >= 8 } {set sel5 [split $pre5 {,}]} 

set skip 10
set step0 [expr $skip - 1] ;# to match start of colvars traj file 
# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD
package require pbctools

# =============================================================== #


# define output file
set outDataFile [open measDist.dat w]
puts $outDataFile "# Input PSF: $inpsf\n# Input DCD, skip $skip: $indcd\n"
puts $outDataFile "# Selection 0: $sel0\n# Selection 1: $sel1\n"

if {[info exists sel2] && [info exists sel3]} {
  puts $outDataFile "# Selection 2: $sel2\n# Selection 3: $sel3\n"
}
if {[info exists sel4] && [info exists sel5]} {
  puts $outDataFile "# Selection 4: $sel4\n# Selection 5: $sel5\n"
}
puts $outDataFile "# Frame | Distance between selections (Angstrom)"

# load files
mol new $inpsf
mol addfile $indcd type {dcd} first $step0 last -1 step $skip waitfor -1 ; # match start of colvars traj file
#mol addfile $indcd type {dcd} first 0 last -1 step 10 waitfor -1
#mol addfile $indcd type {dcd} first 0 last -1 step 1 waitfor -1
#mol addfile $indcd type {dcd} first 0 last 100 step 1 waitfor -1
pbc wrap -compound fragment -center com -centersel "protein and noh" -all 

# loop over frames
set frames [molinfo top get numframes]
set all [atomselect top all]
for {set i 0} {$i < $frames} {incr i} {
    $all frame $i
    if {[info exists sel0] && [info exists sel1]} {
      set vec0 [measure center [atomselect top "$sel0" frame $i]]
      set vec1 [measure center [atomselect top "$sel1" frame $i]]
      set distA [veclength [vecsub $vec0 $vec1]]
    }

    if {[info exists sel2] && [info exists sel3]} {
      set vec2 [measure center [atomselect top "$sel2" frame $i]]
      set vec3 [measure center [atomselect top "$sel3" frame $i]]
      set distB [veclength [vecsub $vec2 $vec3]]
    }

    if {[info exists sel4] && [info exists sel5]} {
      set vec4 [measure center [atomselect top "$sel4" frame $i]]
      set vec5 [measure center [atomselect top "$sel5" frame $i]]
      set distC [veclength [vecsub $vec4 $vec5]]
    }

    if {[info exists distC]} { puts $outDataFile "$i\t$distA\t$distB\t$distC" } \
    elseif {[info exists distB]} { puts $outDataFile "$i\t$distA\t$distB" } \
    else { puts $outDataFile "$i\t$distA\t$distB\t$distC" }
}

close $outDataFile

exit

