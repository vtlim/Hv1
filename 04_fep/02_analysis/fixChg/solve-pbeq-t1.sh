#!/bin/bash

#SBATCH --job-name='pbeq_taut1_fwd_01'
#SBATCH --partition=mf_ilg2.3 
#SBATCH --mem=500mb
#SBATCH --time=01:00:00
#SBATCH --export=ALL 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1-250%1

#------------SIMULATION SETTINGS -----------#
i=$((SLURM_ARRAY_TASK_ID-1))
sourcedir=$SLURM_SUBMIT_DIR
scriptdir="/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/fixChg"

refsys='taut_1'   # 'taut_0, taut_1, taut_2'
window='01'       # '01', '40'
fep_resname='ARG' # 'ARG', 'SER'
fep_resid='211'
framenum='250'

parentdir="/beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/f1_R211S"
coredir="${parentdir}/00_main"
fepdir="${parentdir}/FEP_F/lambda_${window}"
dcdfile="${fepdir}/alchemy${window}.dcd"
toppardir="/work/cluster/limvt/toppar"

#-----------SELECTION SETTINGS-------------#
if [[ "$fep_resname" = "ARG" ]]; then
   betasel="beta < 0"
   invbeta="not beta > 0"
else
   betasel="beta > 0"
   invbeta="not beta < 0"
fi

if [[ "$refsys" = "taut_1" ]]; then
    #System Size (halved)
    X=60
    Y=60
    Z=60
    #Internal Box (halved)
    IX=27
    IY=27
    IZ=35
    
#    MD="0.0000"
    psffile="$coredir/18629-19_R211S.psf"
    pdbfile="$coredir/18629-19_R211S.fep"
    
    sel="(resname GBI1) or (lipid) or \
    (protein and $invbeta) or \
    (water and same residue as within 6 of \
    (protein and resid $fep_resid))"
    stage2="${scriptdir}/02_gencharmm-t1.inp"
    stage3="${scriptdir}/03_pbeq-t1.inp"

elif [[ "$refsys" = "taut_0" ]]; then   # UNVERIFIED
    #System Size (halved)
    X=60
    Y=60
    Z=60
    #Internal Box (halved)
    IX=27
    IY=28
    IZ=35

#    MD="0.0000"
    psffile="$coredir/noGBI_R211S.psf"
    pdbfile="$coredir/noGBI_R211S.fep"

    sel="(lipid) or \
    (protein and $invbeta) or \
    (water and same residue as within 6 of \
    (protein and resid $fep_resid))"
    stage2="${scriptdir}/02_gencharmm-t0.inp"
    stage3="${scriptdir}/03_pbeq-t0.inp"

elif [[ "$refsys" = "taut_2" ]]; then   # UNVERIFIED
    #System Size (halved)
    X=60
    Y=60
    Z=60
    #Internal Box (halved)
    IX=27
    IY=28
    IZ=35

#    MD="0.0000"
    psffile="$coredir/19415-13-c1_R211S.psf"
    pdbfile="$coredir/19415-13-c1_R211S.fep"

    sel="(resname GBI2) or (lipid) or \
    (protein and $invbeta) or \
    (water and same residue as within 6 of \
    (protein and resid $fep_resid))"
    stage2="${scriptdir}/02_gencharmm-t2.inp"
    stage3="${scriptdir}/03_pbeq-t2.inp"
fi

#Ns=$( awk '/OH2/ && /TIP3*/' ${pdbfile} | wc -l )
#QL=-1
#QP=0
#------------------------------------------#


prep4charmm(){
#Generalized script
mkdir -p namd && cd namd
namddir=$( pwd )
echo "###Preparing NAMD trajectory for CHARMM"
cat > $namddir/vmdargs-${i}.tcl <<EOF
set i "$i"
set fep_resid "${fep_resid}"
set fep_resname "${fep_resname}"
set betasel "${betasel}"
set invbeta "${invbeta}"
set sel "$( echo $sel )"
set psffile "${psffile}"
set pdbfile "${pdbfile}"
set dcdfile "${dcdfile}"
set outdir "$namddir"
EOF

cat $namddir/vmdargs-${i}.tcl
vmd -dispdev text -e $scriptdir/01_dcd-crd.tcl -args "$namddir/vmdargs-${i}.tcl"
}

geninput(){
#Must be tailored for system
#Check patches here
echo "###Generating CHARMM input files"

charmmdir="$workdir/charmm"
mkdir -p $charmmdir
charmm title="./charmm" \
       toppar="$toppardir" \
       indir="../namd" \
       i=$i sel="$sel" ref="$ref" \
       -i ${stage2}
}

pbeq(){
#Must be tailored for system
echo "###Starting CHARMM PBEQ"
pbeqdir="$workdir/pbeq"
mkdir -p $pbeqdir
charmm title="./pbeq" \
       toppar="$toppardir" \
       indir="./charmm" \
       i=$i sel="$sel" epsw="$epsw" \
       X="$X" Y="$Y" Z="$Z" \
       IX="$IX" IY="$IY" IZ="$IZ" \
       -i ${stage3}
}

dummysys(){
#Generalized script
echo "###Generating Dummy System"
if [[ "$JOBNAME" = "het-prot" ]]; then
    netchg=0
else
    netchg="$( head -n 1 $pbeqdir/$i.avg \
    | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' \
    | sed -e 's/^0\+//')"
fi
dumdir="$workdir/dummy"
mkdir -p $dumdir
charmm title="./dummy" \
       toppar="$toppardir" \
       i=$i X="$X" Y="$Y" Z="$Z" \
       netchg="$netchg" \
       -i $scriptdir/04_gendummy.inp
}

dummypbeq(){
#Generalized script
echo "###Starting CHARMM PBEQ for Dummy System"
chgdis="$( tail -n 1 $pbeqdir/$i.avg \
| grep -Eo '[+-]?[0-9]+([.][0-9]+)?' \
| sed -e 's/^0\+//' )"
echo $chgdis
charmm title="./pbeq" \
       toppar="$toppardir" \
       dumdir="./dummy" \
       i=$i X="$X" Y="$Y" Z="$Z" \
       IX="$IX" IY="$IY" IZ="$IZ" \
       epsw="$epsw" chgdis="$chgdis" \
       -i $scriptdir/05_dum-pbeq.inp
}

cleanup(){
rm $charmmdir/${i}.crd $charmmdir/${i}.psf
rm $dumdir/${i}_chgd.crd $dumdir/${i}_chgd.psf $dumdir/${i}.avg
find $workdir -type d -empty -exec rmdir {} \;
}

main(){
    cd $sourcedir
    mkdir -p ${JOBNAME}
    cd ${JOBNAME}
    workdir=$( pwd )
    echo "#-------Current Job Directory: ${workdir}----------------#"
    geninput 
    pbeq
    dummysys
    dummypbeq
    cd $sourcedir
}

hetlig(){
    JOBNAME="het-lig"
    epsw="97"
    ref="lig"
    main
}

hetprot(){
    JOBNAME="het-prot"
    epsw="97"
    ref="env"
    main
}

homlig(){
    JOBNAME="hom-lig"
    epsw="1"
    ref="lig"
    main
}


module load compiler/gnu/6.2.0
export PATH=$PATH:/home/limn1/bin
export PATH=$PATH:/beegfs/DATA/mobley/limvt/local/vmd


prep4charmm
hetlig
hetprot
homlig


#if [[ "$i" = "$framenum" ]]; then
#    #Command to compute correction term
#    echo "python $scriptdir/06_ana-corr.py -QL ${QL} -QP ${QP} -X ${X} -Y ${Y} -Z ${Z} -MD ${MD} -Ns ${Ns} -frames ${framenum}"
#fi

