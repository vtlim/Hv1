
# Purpose: Draw a box of +-x, +-y, +-z centered at origin.
# Usage: while in VMD --
#    [1] "source drawBox.tcl"
#    [2] "draw_box 30 30 30" fill in 30 with desired half-x half-y half-z
#    [3] "draw delete all" to start over
# Adapted from: http://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node173.html
# Full path names for sourcing:
#    gpl: /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/fixChg/drawBox.tcl
#    cas: /home/limvt/connect/greenplanet/goto-beegfs/hv1/04_fep/analysis/fixChg/drawBox.tcl
 

proc draw_box {boxcen halfx halfy halfz} {
    # "draw_box com 30 30 30" will draw a box with center at system's center-of-mass and side lengths 60 60 60

    if {$boxcen=="origin" || $boxcen=="orig"} {
        set cenx 0.0
        set ceny 0.0
        set cenz 0.0

    } elseif {$boxcen=="com"} {
        set cen [measure center [atomselect top all]]
        set cenx [lindex $cen 0]
        set ceny [lindex $cen 1]
        set cenz [lindex $cen 2]
    }

    set minx [expr $cenx - $halfx]
    set maxx [expr $cenx + $halfx]

    set miny [expr $ceny - $halfy]
    set maxy [expr $ceny + $halfy]

    set minz [expr $cenz - $halfz]
    set maxz [expr $cenz + $halfz]


    # and draw the lines
    draw color yellow
    draw line "$minx $miny $minz" "$maxx $miny $minz"
    draw line "$minx $miny $minz" "$minx $maxy $minz"
    draw line "$minx $miny $minz" "$minx $miny $maxz"

    draw line "$maxx $miny $minz" "$maxx $maxy $minz"
    draw line "$maxx $miny $minz" "$maxx $miny $maxz"

    draw line "$minx $maxy $minz" "$maxx $maxy $minz"
    draw line "$minx $maxy $minz" "$minx $maxy $maxz"

    draw line "$minx $miny $maxz" "$maxx $miny $maxz"
    draw line "$minx $miny $maxz" "$minx $maxy $maxz"

    draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
    draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
    draw line "$maxx $maxy $maxz" "$maxx $miny $maxz"
}

proc draw_plane {halflength dimension dimvalue style} {
    # "draw_plane 40 z -45 line" will draw a square outline with side length 80 (2x40) at z=-45
    if {$dimension=="x"} {
        set miny [expr 0.0 - $halflength]
        set maxy [expr 0.0 + $halflength]

        set minz [expr 0.0 - $halflength]
        set maxz [expr 0.0 + $halflength]

        set allx $dimvalue

        if {$style=="line"} {
            puts "Drawing a plane outline at x=$dimvalue"
            draw color blue
            draw line "$allx $maxy $maxz" "$allx $miny $maxz" width 2 
            draw line "$allx $miny $maxz" "$allx $miny $minz" width 2
            draw line "$allx $miny $minz" "$allx $maxy $minz" width 2
            draw line "$allx $maxy $minz" "$allx $maxy $maxz" width 2

        } elseif {$style=="plane"} {
            puts "Drawing a plane filled at x=$dimvalue"
            draw color yellow
            draw triangle "$allx $maxy $minz" "$allx $maxy $maxz" "$allx $miny $maxz"
            draw triangle "$allx $miny $maxz" "$allx $miny $minz" "$allx $maxy $minz"
        }
    }

    if {$dimension=="y"} {
        set minx [expr 0.0 - $halflength]
        set maxx [expr 0.0 + $halflength]

        set minz [expr 0.0 - $halflength]
        set maxz [expr 0.0 + $halflength]

        set ally $dimvalue

        if {$style=="line"} {
            puts "Drawing a plane outline at y=$dimvalue"
            draw color blue
            draw line "$maxx $ally $maxz" "$minx $ally $maxz" width 2 
            draw line "$minx $ally $maxz" "$minx $ally $minz" width 2
            draw line "$minx $ally $minz" "$maxx $ally $minz" width 2
            draw line "$maxx $ally $minz" "$maxx $ally $maxz" width 2

        } elseif {$style=="plane"} {
            puts "Drawing a plane filled at y=$dimvalue"
            draw color yellow
            draw triangle "$maxx $ally $minz" "$maxx $ally $maxz" "$minx $ally $maxz"
            draw triangle "$minx $ally $maxz" "$minx $ally $minz" "$maxx $ally $minz"
        }
    }

    if {$dimension=="z"} {
        set minx [expr 0.0 - $halflength]
        set maxx [expr 0.0 + $halflength]

        set miny [expr 0.0 - $halflength]
        set maxy [expr 0.0 + $halflength]

        set allz $dimvalue

        if {$style=="line"} {
            puts "Drawing a plane outline at z=$dimvalue"
            draw color blue
            draw line "$minx $miny $allz" "$maxx $miny $allz" width 2 
            draw line "$minx $miny $allz" "$minx $maxy $allz" width 2
        
            draw line "$maxx $miny $allz" "$maxx $maxy $allz" width 2
            draw line "$minx $maxy $allz" "$maxx $maxy $allz" width 2

        } elseif {$style=="plane"} {
            puts "Drawing a plane filled at z=$dimvalue"
            draw color yellow
            draw triangle "$minx $miny $allz" "$minx $maxy $allz" "$maxx $maxy $allz"
            draw triangle "$maxx $maxy $allz" "$maxx $miny $allz" "$minx $miny $allz"
        }
    }
}

