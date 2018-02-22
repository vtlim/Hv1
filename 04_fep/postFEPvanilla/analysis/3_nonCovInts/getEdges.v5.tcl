#
# getEdges.tcl
#
# JAF and ECW
#
# jfreites@uci.edu
#
# version 1 08/14
#
# version 4 06/17  JAF
#
# no cov bond determination (as in original) but add bonds edge list
#		  an initial frame number to divide the task between different cores
#		  may use a pdb (i.e no topology or type selections) 
#
# use a psf/pdb to generate network nodes according to
# Benson and Daggett scheme
#
# Input: PSF(or PDB)/DCDs and psfile.nodes psfile.bonds (node list and bonded list)
# Output: two files
#       {outputfilename}_{framenumber}.edges one edge per line with total number of contacts in 3d column
#       {outputfilename}_{framenumber}.edgatt one edge per line with nature of interaction in 3d column
#
# Uses C-C and C-X cutoff distances as in Benson and Daggett

# System/Trajectory parameters
#------------------------------------------------------------------------
# dataDir: this is concatenated to the myPSF and file names in theFiles and myReference
# workDir: this is concatenated to the output files
# myPSF: your topology file
# trajFileType: file type for your trajectory file
# step: step length used when reading the trajectory files



set workDir /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/f1_R211S/1_S211/analysis/3_nonCovInts/
set dataDir /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/f1_R211S/1_S211/
set myPSF taut1_18629-19_s211.psf 
set trajFileType dcd
set step 1


#-----------------------------------------------------------------------------

# theFiles:
# Provide a TCL list of trajectory file names or use your TCL skills to build it

#we'll ignore the dcds with zero frames after round(numframe/1000)*1000=0
#see deamidated worksheet for details
#the current workDir is set to the "second part" first part adapted from it

set theFiles [list wrap10.dcd]


# theFileRange:
# Provide a TCL list with the first and last frame number to be analyzed in each
# trajectory file.
# Leave theFileRange empty (set to "") if you want all the frames of all the files
# to be analyzed.

#set theFileRange [list first1 last1 first2 last2 ...]
#

#set theFileRange [list 0 2500 0 2500]
set theFileRange ""

#-----------------------------------------------------
# Selection sentence for analysis
#


set mySelection "(protein or resname GBI1 or (water and same residue as name OH2 and within 6 of protein))"


#------------------------------------------------------------------------
# Cutoff Parameters 
#
#
set carbonCutoff 5.4
set Cutoff 4.6

#
#------------------------------------------------------------------------
# Output file name

set outfile taut1_S211

#
#
# Initial frame count (last frame number from previous run)

set begframe 0
#---------------------END of USER INTERFACE ---------------------------------------------------
#
#
#------------------------------------------------------------------------

set myNodeList ${workDir}${myPSF}.nodes
set myBondList ${workDir}${myPSF}.bonds

#Create a dictionary based on Eric's node list

set idf [open $myNodeList r]
while {[gets $idf line] >= 0} {
        set nindex [lindex $line 0]
        dict set network $nindex resname [lindex $line 1]
        dict set network $nindex resid [lindex $line 2]
        dict set network $nindex nloc [lindex $line 3]
        dict set network $nindex ntype [lindex $line 4]
        dict set network $nindex ndef [lindex $line 5]
	foreach thing [lrange $line 7 end] {
		dict set indices $thing $nindex
	}
}
close $idf

#Read node bonded information
set idf [open $myBondList r]
foreach thing [split [read $idf] "\n"] {
   	foreach {i j} $thing {
        	lappend bondsList [list $i $j]
        }
}
close $idf

#------------------------------------------------------------------------
#
#

set molwork [mol new ${dataDir}$myPSF]

animate delete all $molwork


set nframes $begframe

set nohs [atomselect $molwork "$mySelection and noh"]
set nocarbons [atomselect $molwork "$mySelection and not (hydrogen or carbon)"]
set carbonAtoms [atomselect $molwork "$mySelection and carbon"]



# Process beg/end trajectory files lists

if {$theFileRange == ""} {
        foreach dcdfile $theFiles {
                append theFileRange "0 -1 "
        }
} else {
        if {[llength $theFileRange] != [expr 2 * [llength $theFiles]]} {
                puts "the File Range list inconsistent with Files list"
                exit
        }
}

# Loop over files / get one frame at a time /Loop over selections
#------------------------------------------------------------------------
foreach dcdfile $theFiles {begFrame endFrame} $theFileRange {
        animate delete all $molwork
        animate read $trajFileType ${dataDir}$dcdfile beg $begFrame end $endFrame skip $step waitfor all $molwork
        set totframes [molinfo $molwork get numframes]
        for {set frame 0} {$frame < $totframes} {incr frame} {
                incr nframes
		if {$nframes < 10} {
			set framenum 00${nframes}
		} elseif {$nframes < 100} {
			set framenum 0${nframes}
		} else {
			set framenum $nframes
		}
		set idf [open ${workDir}${outfile}_${framenum}.edges w]
		set idf2 [open ${workDir}${outfile}_${framenum}.edgeatt w]
                animate goto $frame
#               prepareFrame $molwork
		$carbonAtoms update
		$nocarbons update
		$nohs update
		set carbonContacts [measure contacts $carbonCutoff $carbonAtoms]
		set Contacts [measure contacts $Cutoff $nocarbons $nohs]
		set Ai [concat [lindex $carbonContacts 0] [lindex $Contacts 0]]
		set Aj [concat [lindex $carbonContacts 1] [lindex $Contacts 1]]
#		write something
		set weight {}
		set interactType {}
		foreach thing1 $Ai thing2 $Aj {
#			puts -nonewline "**** atoms are $thing1 $thing2 "
			set Ni [dict get $indices $thing1]
			set Nj [dict get $indices $thing2]
			if {$Ni != $Nj} {
				set pair [lsort -integer [list $Ni $Nj]]
#				puts "nodes are $pair"
				if {$pair ni $bondsList} {
					dict incr weight $pair
					if {![dict exists $interactType $pair]} {
						set interaction [list [dict get $network [lindex $pair 0] ntype] [dict get $network [lindex $pair 1] ntype]]
						switch -- $interaction {
							{NOP NOP} {
								dict set interactType $pair HPHOB
							}
							{POS NEG} -
							{NEG POS} -
							{POS POS} -
							{NEG NEG} {
								dict set interactType $pair COUL
							}
							{POS DIP} -
							{DIP POS} -
							{NEG DIP} -
							{DIP NEG} -
							{DIP DIP} {
								dict set interactType $pair HBOND
							}
							default {
								dict set interactType $pair STER
							}
						}
					}
				} else {
					#dict incr weight $pair
					#dict set interactType $pair COVAL
				}
			}
		}
		dict for {key val} $weight {
			puts $idf "$key $val"
		}
		dict for {key val} $interactType {
			puts $idf2 "$key $val"
		}
		close $idf
		close $idf2
		puts "frame $nframes"
        }
}
#------------------------------------------------------------------------

exit
