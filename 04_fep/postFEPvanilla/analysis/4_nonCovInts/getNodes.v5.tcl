# 
# getNodes.tcl
#
# ECW and JAF
#
# ericw3 or jfreites @uci.edu
#
# version 2 08/14 (JAF)
# corrected some typos/omissions in the dict
# minor mods here and there to make it more efficient
#
# version 2.1 11/14
# 2.0 only works with charmm do something about it
#
# version 3 01/16
# adds node definitions for HSD and HSP
#
# use a psf to generate network nodes according to 
# Benson and Daggett scheme
#
# Input: PSF file (somename.psf)
# Output: two files
#	somename.psf.nodes is the  network node file: 
#	node index | resname | resid | node location (bb or sc)| node type | node name | total atoms | atoms psf indices
#	somename.psf.bonds contains the bonded connectivity: pairs of node indices
#	
#	


# System/Trajectory parameters
#------------------------------------------------------------------------
# dataDir: this is concatenated to the myPSF and file names in theFiles and myReference
# workDir: this is concatenated to the output files
# myPSF: your topology file

set workDir ""

set dataDir /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/f1_R211S/1_S211/
set myPSF taut1_18629-19_s211.psf


#-----------------------------------------------------------------------------

#Ligand residue name
#

set ligandres GBI1


#
#
#---------------------END of USER INTERFACE ---------------------------------------------------
#
#


set molwork [mol new ${dataDir}$myPSF]

set of [open ${workDir}${myPSF}.nodes w]
set bf [open ${workDir}${myPSF}.bonds w]

set ressel [[atomselect $molwork "protein and name CA"] get resid]

######################################################################
## Define Dictionary of Charge-Based Nodes
## NOTE: Include Hydrogens or not??? 
# ALA
lappend ndict(ALA) "CA CB"
set ntype(ALA) NOP
set ndef(ALA) "AB"
# ARG
lappend ndict(ARG) "CA CB"
lappend ndict(ARG) "CG CD"
lappend ndict(ARG) "NE CZ NH1 NH2"
set ntype(ARG) "NOP NOP POS"
set ndef(ARG) "AB GD EZ"
# ASN
lappend ndict(ASN) "CA CB"
lappend ndict(ASN) "CG OD1 ND2"
set ntype(ASN) "NOP DIP"
set ndef(ASN) "AB GD"
# ASP
lappend ndict(ASP) "CA CB"
lappend ndict(ASP) "CG"
lappend ndict(ASP) "OD1 OD2"
set ntype(ASP) "NOP NOP NEG"
set ndef(ASP) "AB G D"
# CYH (protonated Cysteine... treat as CYS)
# CYS
lappend ndict(CYS) "CA CB"
lappend ndict(CYS) "SG HG1"
set ntype(CYS) "NOP DIP"
set ndef(CYS) "AB G"
# GLY
lappend ndict(GLY) "CA"
set ntype(GLY) "NOP"
set ndef(GLY) "A"
# GLU
lappend ndict(GLU) "CA CB"
lappend ndict(GLU) "CG"
lappend ndict(GLU) "CD OE1 OE2"
set ntype(GLU) "NOP NOP NEG"
set ndef(GLU) "AB G DE"
# GLN
lappend ndict(GLN) "CA CB CG"
lappend ndict(GLN) "CD OE1 NE2" 
set ntype(GLN) "NOP DIP"
set ndef(GLN) "ABG DE"
# HIS (HSE)
lappend ndict(HSE) "CA CB"
lappend ndict(HSE) "CG CD2"
lappend ndict(HSE) "ND1"
lappend ndict(HSE) "CE1 NE2"
set ntype(HSE) "NOP NOP DIP DIP"
set ndef(HSE) "AB GD D E"
# HSD
lappend ndict(HSD) "CA CB"
lappend ndict(HSD) "CG CD2"
lappend ndict(HSD) "NE2"
lappend ndict(HSD) "CE1 ND1"
set ntype(HSD) "NOP NOP DIP DIP"
set ndef(HSD) "AB GD E D"
# HSP
lappend ndict(HSP) "CA CB"
lappend ndict(HSP) "CG CD2"
lappend ndict(HSP) "CE1 ND1 NE2"
set ntype(HSP) "NOP NOP POS"
set ndef(HSP) "AB GD DE"
# ILE
lappend ndict(ILE) "CA CB"
lappend ndict(ILE) "CG1 CG2 CD"
set ntype(ILE) "NOP NOP"
set ndef(ILE) "AB GD"
# LEU
lappend ndict(LEU) "CA CB"
lappend ndict(LEU) "CG CD1 CD2"
set ntype(LEU) "NOP NOP"
set ndef(LEU) "AB GD"
# LYS
lappend ndict(LYS) "CA CB"
lappend ndict(LYS) "CG CD CE"
lappend ndict(LYS) "NZ"
set ntype(LYS) "NOP NOP POS"
set ndef(LYS) "AB GDE Z"
# MET
lappend ndict(MET) "CA CB"
lappend ndict(MET) "CG SD CE"
set ntype(MET) "NOP NOP"
set ndef(MET) "AB GDE"
# PHE
lappend ndict(PHE) "CA CB"
lappend ndict(PHE) "CG CD1 CD2 CE1 CE2 CZ"
set ntype(PHE) "NOP NOP"
set ndef(PHE) "AB GDEZ"
# PRO
lappend ndict(PRO) "CA CB CG CD"
set ntype(PRO) NOP
set ndef(PRO) "ABGD"
# SER
lappend ndict(SER) "CA CB"
lappend ndict(SER) "OG"
set ntype(SER) "NOP DIP"
set ndef(SER) "AB G"
# THR
lappend ndict(THR) "CA CB CG2"
lappend ndict(THR) "OG1"
set ntype(THR) "NOP DIP"
set ndef(THR) "ABG G"
# TRP
lappend ndict(TRP) "CA CB CG CD1"
lappend ndict(TRP) "NE1"
lappend ndict(TRP) "CE2 CE3 CZ2 CZ3 CH2 CD2"
set ntype(TRP) "NOP DIP NOP"
set ndef(TRP) "ABGD E EZH"
# TYR
lappend ndict(TYR) "CA CB"
lappend ndict(TYR) "CG CD1 CD2 CE1 CE2 CZ"
lappend ndict(TYR) "OH"
set ntype(TYR) "NOP NOP DIP"
set ndef(TYR) "AB GDEZ H"
# VAL
lappend ndict(VAL) "CA CB CG1 CG2"
set ntype(VAL) NOP
set ndef(VAL) "ABG"
# CTR (C-terminal carboxyl)
lappend ndict(CTR) "C"
lappend ndict(CTR) "OT1 OT2"
set ntype(CTR) "NOP NEG"
set ndef(CTR) "C OT"
# TIP3
lappend ndict(TIP3) OH2
set ntype(TIP3) DIP
set ndef(TIP3) OH2

# GBI1
 lappend ndict(GBI1) "C2 C3 C4 C5 C6 C7"
 lappend ndict(GBI1) "N3"
 lappend ndict(GBI1) "C1 N4"
 lappend ndict(GBI1) "N N1 N2 C"
 set ntype(GBI1) "NOP DIP DIP POS"
 set ndef(GBI1) "BENZ NN CN GUAN"
# GBI2
 lappend ndict(GBI2) "C9 C10 C11 C12 C13 C14"
 lappend ndict(GBI2) "N8 C7 N15"
 lappend ndict(GBI2) "N3 N1 N6 C2"
 set ntype(GBI2) "NOP POS DIP"
 set ndef(GBI2) "BENZ CNN GUAN"
#
######################################################################
set nindex 1
set bbnindex 1

## Check for N-terminal residue
set res [lindex $ressel 0]
set sel [atomselect $molwork "protein and resid $res and type NH3"]
#set sel [atomselect $molwork "index 0"]
set indexsel [$sel get index]

## NOTE: Add check for no n-terminus and add regular dipolar node
if {[llength $indexsel] == 1} {
	set resname [$sel get resname]
	puts $of "${nindex}\t${resname}\t${res}\tbb\tPOS\tNT\t[llength $indexsel]\t$indexsel"
	incr nindex
	puts  -nonewline $bf "$bbnindex $nindex "
} else {
	puts "ERROR: N-terminal Amine not found"
	exit
}

set resname [$sel get resname]
$sel delete
set lista {}
foreach name $ndict($resname) type $ntype($resname) def $ndef($resname) {
	set sel [atomselect $molwork "protein and resid $res and name $name"]
	puts $of "${nindex}\t${resname}\t$res\tsc\t$type\t$def\t[llength $name]\t[$sel get index]"
	lappend lista $nindex
	incr nindex
	$sel delete
}
if {[llength $lista] > 1} {
foreach p1 [lrange $lista 0 end-1] p2 [lrange $lista 1 end] {
                        puts -nonewline $bf "$p1 $p2 "
}
puts $bf ""
}
set preCA [lindex $lista 0]


## Do the loop for middle residues
foreach res [lrang $ressel 1 end] {
	set sel [atomselect $molwork "protein and name CA and resid $res"]
	set resname [$sel get resname]
	$sel delete

	set lista {}
	# Add sidechain nodes
	foreach name $ndict($resname) type $ntype($resname) def $ndef($resname) {
		set sel [atomselect $molwork "protein and resid $res and name $name"]
		puts $of "${nindex}\t${resname}\t${res}\tsc\t${type}\t$def\t[llength $name]\t[$sel get index]"
		lappend lista $nindex
		incr nindex
		$sel delete
	}
if {[llength $lista] > 1} {
foreach p1 [lrange $lista 0 end-1] p2 [lrange $lista 1 end] {
                        puts -nonewline $bf "$p1 $p2 "
}
puts $bf ""
}

	# Add backbone dipolar node
	# residue n
	# N from residue n
	# C and O from residue n-1
	set resm1 [expr $res-1]
	set sel [atomselect $molwork "(protein and (resid $resm1 and name C O)) or (protein and (resid $res and name N))"]
	set indices [$sel get index]

	
	if {[llength $indices] == 3} {
		puts $of "$nindex\t${resname}\t${res}\tbb\tDIP\tBB\t[llength $indices]\t$indices"
		puts -nonewline $bf "$nindex [lindex $lista 0] "
		puts $bf "$nindex $preCA "
		incr nindex
	} else {
		puts "ERROR: Backbone dipolar node selection != 3 atoms"
		exit
	}
	set preCA [lindex $lista 0]
}

## Check for C-terminal residue
#set res [lindex $ressel [expr $reslength - 1]]
set res [lindex $ressel end]
set lista {}
foreach name $ndict(CTR) type $ntype(CTR) def $ndef(CTR) {
	set sel [atomselect $molwork "protein and resid $res and name $name"]
	set resname [lindex [$sel get resname] 0]
	puts $of "${nindex}\t${resname}\t${res}\tbb\t${type}\t$def\t[llength $name]\t[$sel get index]"
	lappend lista $nindex
	incr nindex
	$sel delete
}

# write to the bond list
if {[llength $lista] > 1} {
foreach p1 [lrange $lista 0 end-1] p2 [lrange $lista 1 end] {
                        puts -nonewline $bf "$p1 $p2 "
}
puts $bf ""
}
puts $bf "$preCA [lindex $lista 0]"

#add waters
set count 0
set name OH2
set def OH2
set type DIP
set resname WAT
set sel [atomselect $molwork "name $name"]
foreach ind [$sel get index] {
	puts $of "${nindex}\t${resname}\t${count}\twt\t${type}\t$def\t[llength $name]\t$ind"
	incr count
	incr nindex
}
$sel delete



set resname $ligandres
set lista {}
        # Add sidechain nodes
        foreach name $ndict($resname) type $ntype($resname) def $ndef($resname) {
                set sel [atomselect $molwork "resname $resname and name $name"]
                puts $of "${nindex}\t${resname}\t${res}\tlg\t${type}\t$def\t[llength $name]\t[$sel get index]"
                lappend lista $nindex
                incr nindex
                $sel delete
        }
if {[llength $lista] > 1} {
foreach p1 [lrange $lista 0 end-1] p2 [lrange $lista 1 end] {
                        puts -nonewline $bf "$p1 $p2 "
}
puts $bf ""
}





## Clean up
close $of
close $bf
exit
