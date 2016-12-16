### combine the Hv1+2GBI docked config pdb with that of the full Hv1 system 
### with lipids+waters. Have to align protein structures.

### Make sure to create the writepdb subdirectory (see last line) first.
### Usage: vmdt -e file.tcl

set dir "/pub/limvt/hv1/01_docking/2_tautomer"

cd $dir/pdbqt/wt_ans/x_together/
set names [lsort [glob hv1-config-*_gbi-pose-*.pdb]]

set count 0
foreach i $names {
    cd $dir/pdbqt/wt_ans
    set next [expr $count + 1]

    # open gbi+prot
    mol new ./x_together/$i
    puts "Loaded docked GBI configuration..."

    # get frame number from name
    set config [regexp -all -inline -- {[0-9]+} [lindex [split $i -] 2]]
    
    # open prot+lipids+water
    mol new $dir/pdb/full-hv1-system/hv1_depolarized_frame_${config}.pdb
    puts "Loaded corresponding protein + soln. config..."

    # align by protein CA
    set sel0 [atomselect $count "protein and name CA and resid 88 to 229"]
    set sel1 [atomselect $next "protein and name CA and resid 88 to 229"]
    set all0 [atomselect $count all]
    
    set transmat [measure fit $sel0 $sel1]

    $all0 move $transmat

    [atomselect $count "resname GBI"] writepdb "y_ligCentered/${i}"

    incr count
    incr count
    $sel0 delete
    $sel1 delete
    $all0 delete
    unset transmat
    mol delete all
}

#exit
