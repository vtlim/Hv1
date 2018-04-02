
# View 2GBI in last frame of all FEP windows
# Arguments: psf file, last frame number of dcd
# vmd -e viewAllWins.tcl -args file.psf 250

cd /work/cluster/limvt/hv1/04_fep/17041_19/b1_F182A/
set psf [lindex $argv 0]
set count 0

color Display Background white
display projection Orthographic
display depthcue off


foreach set {FEP_F FEP_R} {
    cd $set

    foreach window [glob lambda_*] {
        cd $window
        set dcd [glob *dcd]
        mol new $psf type {psf} first 0 last -1 step 1 waitfor all
        mol addfile $dcd type {dcd} first [lindex $argv 1] last -1 step 1 waitfor -1 $count
        mol delrep 0 $count
    
        # add representation for file with GBI only
        mol addrep $count
        mol modselect 0 $count resname GBI1
        mol modstyle 0 $count Licorice 0.100000 12.000000 12.000000
        mol modcolor 0 $count Name
        
        # hide displays for ease of viewing
        mol off $count           ;# don't display anything to start
    
        incr count 
        cd ../
    }

    cd ../
}
