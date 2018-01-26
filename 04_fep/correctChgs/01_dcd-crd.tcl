#!/usr/bin/tclsh


proc pdb2crd {argv} {
    set vmdargs [lindex $argv 0]
    source $vmdargs

    #Load core PSF/PDB and traj files 
    mol new "$psffile"
    #Frame 1
    mol addfile "$pdbfile" type {pdb}
    #Frame 2
    mol addfile "$dcdfile" type {dcd} first $i last $i step 1 waitfor all
    
    #Wrap cell around ligand
    package require pbctools
    pbc wrap -compound fragment -centersel "protein and resid $fep_resid" -center com -first 1 -last 1

    #Center ligand at origin
    set ligCoords [measure center [atomselect top "protein and resid $fep_resid"] weight mass]
    set toOrigin [vecinvert $ligCoords]
    [atomselect top all] moveby $toOrigin

    #Restore original resname
    [atomselect top "resid $fep_resid"] set resname $fep_resname

    #Select atoms from transformation
    #Remove trailing A/B label used for FEP
    set beta [[atomselect top "resid $fep_resid and $betasel"] get name]
    foreach {n} $beta {
        set o [string range $n 0 end-1]
        puts "$n $o"
        [atomselect top "resid $fep_resid and name $n"] set name $o
    }
    
    #Write out PDB for each frame
    set outpdb $outdir/$i.pdb
    set outcrd $outdir/$i.crd
    puts "      Converting $outpdb to $outcrd"
    [atomselect top $sel] writepdb $outpdb

    #Load the pdb file to data stream
    set idf [open $outpdb r]
    set pdbstr [read $idf]
    set pdbstr [split $pdbstr "\n"]
    close $idf

    #Create crd file
    set idf2 [open $outcrd w]
    puts $idf2 "*"
    #Total num atoms
    puts $idf2 [llength [lrange $pdbstr 1 end-2]]

    #parse pdb file and write crd
    set curresnum -99
    set currsegname "PLACEHOLDER" ;# segments w/only 1 resid fail incr
    set rserial 0
    foreach line [lrange $pdbstr 1 end-2] {
        set type [string range $line 0 5]
        set atomnum [string trim [string range $line 6 10]]
        set atomname [string range $line 12 15]
        set resname [string range $line 17 20]
        set resnum [string trim [string range $line 22 25]]
        set x [string range $line 30 37]
        set y [string range $line 38 45]
        set z [string range $line 46 53]
        set segname [string range $line 72 75]
      
        #Set residue numbering from PDB
        #Set sequential residue numbering for CHARMM
        if {$resnum != $curresnum || $segname != $currsegname} { 
            set  curresnum $resnum
            set  currsegname $segname
            incr rserial
        }

        # Patch some atom types from hybrid notation to regular notation
        set patchdict [dict create H11 "HH11" H12 "HH12" H21 "HH21" H22 "HH22"]
        dict for {key value} $patchdict {
           if { [string trim $atomname]  == $key} {
               puts "Replacing $resname $atomname atom type to $value"
               set atomname $value
           }
        }

        puts $idf2 [format "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %4s %-4s%10.5f" [string trim $atomnum] $rserial $resname [string trim $atomname] $x $y $z $segname [string trim $resnum] 0.0]
    }
    close $idf2
}

pdb2crd $argv

exit
