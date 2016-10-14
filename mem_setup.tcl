# Generate psf file for residue mutation. Usage:
# "vmdt file.psf file.pdb -e file.tcl > psfgen.out"

set outname 15183_04-F182A
package require psfgen

### separate pdb by segname for processing
mkdir "split_pdb"
cd "split_pdb"
set segs [lsort -unique [[atomselect top "all"] get segname]]
    foreach seg $segs {
        [atomselect top "segname $seg"] writepdb ${seg}.pdb
    }

### Set CHARMM protein/membrane parameters
topology /work/cluster/limvt/toppar/cgenff2b6/top_all36_cgenff.rtf
topology /work/cluster/limvt/toppar/top_all36_prot.rtf
topology /work/cluster/limvt/toppar/top_all36_lipid.rtf
topology /work/cluster/limvt/toppar/toppar_water_ions.jaf2.str
topology /work/cluster/limvt/toppar/gbi_final.str
topology /work/cluster/limvt/hv1/04_fep/top36_F2A_hybrid.inp

### Rebuild the protein segment
### Mutate the residue of interest
segment PROA {
pdb PROA.pdb
mutate 182 F2A
}

### Specify the protein coordinates
### Regnerate the protein angles/dihedrals with hybrid residue
coordpdb PROA.pdb PROA
regenerate angles dihedrals

### Rebuild the membrane/water segments
### Respecify coordinates
### auto none should be used at least for WAT
segment GBI1 {pdb GBI1.pdb}
segment MEMB {pdb MEMB.pdb}
segment WAT1 {auto none; pdb WAT1.pdb}
segment WAT2 {auto none; pdb WAT2.pdb}
segment SOD {pdb SOD.pdb}
segment CLA {pdb CLA.pdb}

coordpdb GBI1.pdb GBI1
coordpdb MEMB.pdb MEMB
coordpdb WAT1.pdb WAT1
coordpdb WAT2.pdb WAT2
coordpdb SOD.pdb SOD
coordpdb CLA.pdb CLA

### Guesses positions of newly inserted atoms from hybrid residue
### ! may raise warning of poorly guessed coordinates !
### Safe to ignore if only new atoms raise warning
guesscoord

### Write PDB and PSF of hybrid system
### Careful not to overwrite files
cd ../
writepdb ../00_main/${outname}.pdb
writepsf ../00_main/${outname}.psf

### Optional steps
# Solvate and ionize system?
# Might not be necessary with no charge change.

### Get PBC and origin
mol delete all
resetpsf
mol new ../00_main/${outname}.psf
mol addfile ../00_main/${outname}.pdb

set sel [atomselect top all]
set mm [measure minmax $sel]
set cen [measure center $sel]
puts "\n\nFYI only, use actual cell dimensions from PDB..."
puts "MinMax: $mm"
puts "Center: $cen"
exit
