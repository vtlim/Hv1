
# Purpose: Count number of waters within x Angstrom of specified selection.

# Usage: vmdt -e ligCloseWat.tcl -args inpsf indcd selection
# Example: vmdt -e ligCloseWat.tcl -args file.psf file.dcd sidechain,and,resid,211 sidechain,and,resid,112 sidechain,and,resid,185
#
# By: Victoria Lim


# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set dcdlist [lindex $argv 1]
#set prelig [lindex $argv 2]

for {set i 2} {$i < $argc} {incr i} {
    set prelig [lindex $argv $i]
    lappend sellist [split $prelig {,}]
}
puts $sellist


set dist 4
set watlist [list]

# TEMP until this function is moved into analyzeDCD.tcl
set inskip ""
set inpdb ""

# =============================================================== #

proc average L {
    expr ([join $L +])/[llength $L].
}

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD
package require pbctools

# =============================================================== #

proc count_wat_near { outfile {dist 4} args } {
    # ============================================================
    # Count number of waters within [dist] Angstroms of the selections.
    #
    # Arguments
    #  - outfile : string
    #      Name of the output file.
    #  - dist : integer
    #      Number of Angstroms away from selection to find waters.
    #  - args : strings
    #      Variable number of arguments containing string selections for VMD.
    # Returns
    #  - (nothing)
    # Example usage
    #  - count_wat_near waters-near-selections.dat sidechain,and,resid,211 resname,GBI1
    # Notes
    #  - To specify selection, separate words with commas, not spaces. ex: protein,and,resid,112
    #
    # ========================== Variables ========================= #
    global inpsf
    global inskip
    global inpdb
    global dcdlist

    # modify the args, for example:
    # start:  {{protein and resid 211} {protein and resid 112} {}}
    # final:  {protein and resid 211} {protein and resid 112}
    set args [lindex $args 0]
    set args [lreplace $args end end]

    # define output file
    set outDataFile [open $outfile w]
    puts $outDataFile "# Number of waters (noh) within $dist Angstroms"
    puts $outDataFile "# Input PSF: $inpsf\n# Input DCD, skip $inskip: $dcdlist\n"
    set header "# Frame"
    for {set i 0} {$i < [llength $args]} {incr i} {
        puts $outDataFile "# Selection $i: [lindex $args $i]"
        append header "\t$i"
    }
    puts $outDataFile "\n$header"

    # waters calculation
    puts "Counting waters..."
    set num_steps [molinfo 0 get numframes]
    for {set frame 0} {$frame < $num_steps} {incr frame} {
        if {[expr $frame % 100 == 0]} {puts $frame}
        set curr_line "$frame\t"
        foreach x $args {
            set cw [atomselect 0 "water and oxygen within $dist of ($x and noh)" frame $frame]
            set num [$cw num]
            append curr_line "\t$num"
        }
        puts $outDataFile $curr_line
    }
    close $outDataFile
}

# load files
mol new $inpsf
mol addfile $dcdlist type {dcd} first 0 last -1 step 10 waitfor -1
#mol addfile $indcd type {dcd} first 0 last 100 step 1 waitfor -1
pbc wrap -compound fragment -center com -centersel "protein" -all

count_wat_near waters-near-selections.dat 4 $sellist

