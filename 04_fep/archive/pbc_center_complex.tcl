package require pbctools

set dcdfile ../alchemy01.dcd

set mol [mol new {$dcdfile} type {dcd} first 0 last -1 step 1 waitfor all]
mol addfile {1f0l_hairpin_352o349o362n.psf} type {psf} first 0 last -1 waitfor all top
molinfo top get numframes


set protmem [atomselect $mol "segid DTT or segid POPC or segid POPG"]

#Center membrane at origin, keep segments from wrapping
pbc wrap -centersel "segid POPC or segid POPG" -center origin -compound residue -all
#Center DTT protein at center of bound box
pbc wrap -centersel "segid DTT" -center bb -all

animate write pdb frames.pdb beg 0 end -1 sel $protmem

quit
