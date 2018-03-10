
#
# Call this script in the directory where the vanilla traj is to be run.
# sh /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/fep2Vanilla.sh > setup.out
# ________________________________________________________________________

refsys='taut2_19415-13-c1'             # noGBI, taut1_18629-19, taut2_19415-13-c1
prefix="${refsys}_a181"    # noGBI_e112, taut1_18629-19_e112, taut2_19415-13-c1_e112   
mutcode='S2A'           # D2E, V2A, etc.
endcode='ALA'           # GLU, ALA, etc.
inpsf='../00_main/19415-13-c1_S181A.psf'   # with relative path, ex: ../00_main/18629-19_V178A.psf


if [[ "$refsys" = "noGBI" ]]; then
    printf "\n\n1. copying over needed files\n"
    cp /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/mem_setup0.tcl .
    cp /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/equil01-cas.inp .
    printf "\n\n2. generating PDB from last frame of fwd40 window\n"
    /beegfs/DATA/mobley/limvt/local/vmd/vmd -dispdev none -e /work/cluster/limvt/analysis/writePDB.tcl \
        -args $inpsf ../FEP_F/lambda_40/alchemy40.dcd 249 ${refsys}_${mutcode}.pdb 1 1
    printf "\n\n3. modifying PDB at mutation residue\n"
    python /work/cluster/limvt/analysis/editFEPpdb.py ${refsys}_${mutcode}.pdb ${refsys}_${mutcode}_fromPY.pdb $mutcode $endcode
    printf "\n\n4. updating and using psfgen script for new PSF/PDB\n"
    sed -i "s/^set outname .*/set outname $prefix/" mem_setup0.tcl
    /beegfs/DATA/mobley/limvt/local/vmd/vmd -dispdev none ${refsys}_${mutcode}_fromPY.pdb -e mem_setup0.tcl > psfgen.out

elif [[ "$refsys" = "taut1_18629-19" ]]; then
    printf "\n\n1. copying over needed files\n"
    cp /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/mem_setup1.tcl .
    cp /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/equil01-cas.inp .
    printf "\n\n2. generating PDB from last frame of fwd40 window\n"
    /beegfs/DATA/mobley/limvt/local/vmd/vmd -dispdev none -e /work/cluster/limvt/analysis/writePDB.tcl -args $inpsf ../FEP_F/lambda_40/alchemy40.dcd 249 ${refsys}_${mutcode}.pdb 1 1
    printf "\n\n3. modifying PDB at mutation residue\n"
    python /work/cluster/limvt/analysis/editFEPpdb.py ${refsys}_${mutcode}.pdb ${refsys}_${mutcode}_fromPY.pdb $mutcode $endcode
    printf "\n\n4. updating and using psfgen script for new PSF/PDB\n"
    sed -i "s/^set outname .*/set outname $prefix/" mem_setup1.tcl
    /beegfs/DATA/mobley/limvt/local/vmd/vmd -dispdev none ${refsys}_${mutcode}_fromPY.pdb -e mem_setup1.tcl > psfgen.out

elif [[ "$refsys" = "taut2_19415-13-c1" ]]; then
    printf "\n\n1. copying over needed files\n"
    cp /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/mem_setup2.tcl .
    cp /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/equil01-t2-cas.inp .
    printf "\n\n2. generating PDB from last frame of fwd40 window\n"
    /beegfs/DATA/mobley/limvt/local/vmd/vmd -dispdev none -e /work/cluster/limvt/analysis/writePDB.tcl -args $inpsf ../FEP_F/lambda_40/alchemy40.dcd 249 ${refsys}_${mutcode}.pdb 1 1
    printf "\n\n3. modifying PDB at mutation residue\n"
    python /work/cluster/limvt/analysis/editFEPpdb.py ${refsys}_${mutcode}.pdb ${refsys}_${mutcode}_fromPY.pdb $mutcode $endcode
    printf "\n\n4. updating and using psfgen script for new PSF/PDB\n"
    sed -i "s/^set outname .*/set outname $prefix/" mem_setup2.tcl
    /beegfs/DATA/mobley/limvt/local/vmd/vmd -dispdev none ${refsys}_${mutcode}_fromPY.pdb -e mem_setup2.tcl > psfgen.out
fi

printf "\n\n5. Edit NAMD input file and run.\n  structure\n  coordinates\n  cellOrigin\n  runSteps\n  VIEW system and check box size.\n\n"
