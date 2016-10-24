
# Usage: vmdt -e file.tcl
# Purpose: Get PDBs of hHv1 without ligand for FEP calculations
#   using protein poses from docking, part2, located at 
#   /tw/ageragot/hv1/mdruns-hv1+gbi.
# For the listed protein number, subtract one to get corresponding
#   frame of anton trajectory, located at
#   /pub/ageragot/hv1/anton_trajectories/depolarized
# Checked manually (15183, 19744) and are complete overlays.

set poselist {15159 15183 15183 17041 17041 17041 17120 17480 18629 19415 19460 19744}
set psf /share/pub/ageragot/hv1/anton_trajectories/depolarized/hHv1.psf
set dcd /share/pub/ageragot/hv1/anton_trajectories/depolarized/posE.hv1-every4-21us_wrapped_centered.dcd

set count 0
foreach pose $poselist {
    puts "--------- Starting on pose: $pose ----------"
    set p [expr {$pose - 1} ]

    mol new $psf type {psf} waitfor all
    mol addfile $dcd type {dcd} first $p last $p step 1 waitfor -1 $count
    animate write pdb /share/pub/limvt/hv1/02_configs/pdbs-no-gbi/${pose}.pdb beg 0 end 0 skip 1 $count
    mol delete $count
    incr count
}

exit
