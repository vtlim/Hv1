
# Purpose: Write out coordinates from some NAMD trajectory.
# Usage: "vmdt -e writePDB.tcl -args input.psf input.dcd frame# output.pdb wrapBool centerBool"

# TODO (If using centerBool, use vmd and not vmdt bc can't figure out why pbc dependency is missing)
#
# wrapBoolean:   0 means no wrap, 1 means wrap with dopbc
# centerBoolean: 0 means don't alter, 1 means wrap around protein
#
# If writing out last frame, take catdcd numFrames-1

# ========================== Variables ========================= #

set inpsf  [lindex $argv 0]
set indcd  [lindex $argv 1]
set fnum   [lindex $argv 2]
set outpdb [lindex $argv 3]
set wrap   [lindex $argv 4]
set center [lindex $argv 5]

# =============================================================== #



# for dopbc
#lappend auto_path /home/victoria/Documents/tempotools/libs
package require tempoUserVMD

mol new $inpsf

if {$center} {
  dopbc -file $indcd -frames $fnum:1:$fnum -ref protein
} elseif {$wrap && !$center} {
  dopbc -file $indcd -frames $fnum:1:$fnum
} else {
  mol addfile $indcd type {dcd} first $fnum last $fnum step 1 waitfor -1
}

# if you don't have dopbc from tempotools use vmd
#if {$center} {pbc wrap -centersel "protein" -center com -compound residue -all}


animate write pdb $outpdb beg 0 end 0 skip 1
exit


