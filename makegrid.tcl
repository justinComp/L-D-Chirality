# makegrid dKYFIL dKYFIL dKYFIL CHARMM_36_Chiral 27/27 0/27 27 3 3 3 20 6
# makegrid LnewIons Lkyfil Lkyfil CHARMM_36 64/64 0/64 64 4 4 4 20 4
# makegrid 4kyfillyfik 4kyfillyfik 4kyfillyfikout 1 0 1 1 1 1 32 5

proc makegrid {speciesL speciesD outputname ratioL ratioD totalpep gridX gridY gridZ unit_buffer solvatePad} {
    #INPUT
    # speciesL = string name of the first type of peptide being put in the grid (ex Lkyfil)
    # speciesD = string name of the second type of peptide being put in the grid (ex Dkyfil)
    # outputname = once all peptides are put in the grid, we will use this string to name the output files
    # topofile = the full topology file that will be used in making the bond info for the sim (and file type!)
    # ratioL = fration (in fraction form) of the ratio of species 1 (ex 16/64 would be 1L : 3D if species 1 was L)
    # ratioD = fration (in fraction form) of the ratio of species 2 (ex 48/64 would be 1L : 3D if species 2 was D)
    # totalpep = the total number of peptides being simulates (ex 18 or 64)
    # gridX = for the x dimensions of the grid placement, also determines how many peptides will be placed in the grid
    # gridY = for the y dimensions of the grid placement, also determines how many peptides will be placed in the grid
    # gridZ = for the z dimensions of the grid placement, also determines how many peptides will be placed in the grid
    # unit_buffer = how much space is placed between the peptides (15-32 is a good range we have found, adjust to your molarity!)
    # solvatePad = in the auto solvate line, this will dictate how far the solvent cube is from the outermost peptides! (default 6 if you cannot decide)

    puts $ratioL
    puts $ratioD

    # MOST UPDATED SCRIPT AS OF 8/13/2022 #

    # used to place 1 type of peptide into a grid, then solvate and ionize it

    # .pdb has location while .psf has charge, mass, and bond info

    # this is the file used to take structure into a psf
    package require psfgen

    # https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf
    topology top_all36_carb.rtf
    topology top_all36_lipid.rtf
    topology top_all36_prot.rtf
    topology top_all36_cgenff.rtf
    topology toppar_water_ions.str

    topology toppar_all36_prot_c36m_d_aminoacids_test.str


    # CHANGE to the file name of your pdb from avogadro
    # molname1 is for the L version
    set molname1 "$speciesL"

    # CHANGE molname2 is for the D version
    set molname2 "$speciesD"

    # CHANGE final mol name of merged peptide
    set molname "$outputname"

    # CHANGE the number of peptides you are simulating
    set numPep $totalpep

    # makes the file that will have a list of L and L peptides
    set outFile [open "makeInfoFor$molname.txt" "w"]

    # procedure to get random numbers
    proc myRand {min max} {
        expr {int(rand() * ($max + 1 - $min)) + $min}
    }

    # CHANGE divides grid into 4 sections (16 peptides per section by ratio)

    set sectionRatioL [expr ($numPep/$gridX)*($ratioL)]
    set sectionRatioD [expr ($numPep/$gridX)*($ratioD)]
    set LD []
    set LDwords []

    for {set i 1} {$i <= $gridX} {incr i} {
        set smallLcount 0
        set smallDcount 0
        for {set section 1} {$section <= [expr ($numPep/$gridX)]} {incr section} {
            set LorD [myRand 0 1]

            # for higher ratio of L (ex 3:1)
            if {$ratioL > $ratioD} {
                if {$smallLcount == $sectionRatioL && $smallDcount == $sectionRatioD} {
                    puts "Ratio of L and D has been reached. "

                    # for more L than D... $LorD > 0, selects L
                    } elseif {$LorD > 0 && $smallLcount < $sectionRatioL} {
                        lappend LD 0
                        set smallLcount [expr $smallLcount + 1]
                        puts "L"

                    # for more D than L... $LorD > 0, selects D
                    } elseif {$LorD == 0 && $smallDcount < $sectionRatioD} {
                        lappend LD 1
                        set smallDcount [expr $smallDcount + 1]
                        puts "D"

                    # if D is selected but we have maxed out how much D we can have, selects L
                    } elseif {$LorD == 0 && $smallDcount == $sectionRatioD} {
                        lappend LD 0
                        set smallLcount [expr $smallLcount + 1]
                        puts "L"

                    # if L is selected but we have maxed out how much L we can have, selects D
                    } elseif {$LorD > 0 && $smallLcount == $sectionRatioL} {
                        lappend LD 1
                        set smallDcount [expr $smallDcount + 1]
                        puts "D"
                    }

                # for higher ratio of D (ex 1:3)
                } elseif {$ratioL < $ratioD} {
                    if {$smallLcount == $sectionRatioL && $smallDcount == $sectionRatioD} {
                        puts "Ratio of L and D has been reached. "

                        # for more D than L... $LorD == 0 selects L
                        } elseif {$LorD == 0 && $smallLcount < $sectionRatioL} {
                            lappend LD 0
                            set smallLcount [expr $smallLcount + 1]
                            puts "L"

                        # for more L than D... $LorD == 0, selects D
                        } elseif {$LorD > 0 && $smallDcount < $sectionRatioD} {
                            lappend LD 1
                            set smallDcount [expr $smallDcount + 1]
                            puts "D"

                        # if D is selected but we have maxed out how much D we can have, selects L
                        } elseif {$LorD > 0 && $smallDcount == $sectionRatioD} {
                            lappend LD 0
                            set smallLcount [expr $smallLcount + 1]
                            puts "L"

                        # if L is selected but we have maxed out how much L we can have, selects D
                        } elseif {$LorD == 0 && $smallLcount == $sectionRatioL} {
                            lappend LD 1
                            set smallDcount [expr $smallDcount + 1]
                            puts "D"
                        }}}}

    puts $LD
    puts "A lil FYI: L is 0 and D is 1"

    # creates a list of IDs to merge peptides 
    set pepList []

    ##################################################
    # Loads in multiple protiens of one PDB and      #
    # rotates the peptides randomly (all at origin)  #
    ##################################################

    # adds selected number of peptides from James Tang's peptide_build.tcl

    #GET UNIT CELL.

    # CHANGE this to the ratio you need out of 64 peptides
    # ADD one so that the logic works out, ie 1 to 1 is 32 to 32 but put 33
    set Lratio [expr (($ratioL*$numPep) + 1)]
    set Dratio [expr (($ratioD*$numPep) + 1)]

    # these variables will change with iteration and are used to check when enough 
    # have been placed into the grid
    set Lcount 0
    set Dcount 0
    set i 1
    unset LorD

    puts "entering grid loop"
    # these loop through the points in the grid system
    for {set xOffset 0} {$xOffset < $gridX} {incr xOffset} {
        for {set yOffset 0} {$yOffset < $gridY} {incr yOffset} {
            for {set zOffset 0} {$zOffset < $gridZ} {incr zOffset} {

                ### makes random 0 or 1 (0 is L, 1 is D) ###
                set LorD [lindex $LD [expr $i - 1]]
                ### IF YOU ARE DOING PURE L OR D SET 0 OR 1 ACCORDINGLY ###
                # set LorD 0

                if {$Lcount == $Lratio && $Dcount == $Dratio} {
                    puts "Ratio of L and D has been reached. "
                } elseif {$LorD == 0 && $Lcount < $Lratio} {
                     
                     # loads in the single peptide file and increments 
                     mol new ${molname1}.pdb
                     mol addfile ${molname1}.psf
                     set Lcount [expr $Lcount + 1]
                     puts $Lcount
                     puts "L"

                     # selects the segment to rename
                     set sel [atomselect top "segname PROA"]
                     $sel set segname PL$i

                     segment PL$i {
                     pdb ${molname1}.pdb
                     first NTER
                     last CT2
                     }

                     lappend LDwords PL$i

                     ############################################################
                     set sel [atomselect top all]
                     $sel moveby [vecinvert [measure center $sel]]

                     set mm [measure minmax $sel]
                     #rotate unit cell (with boltzmann distribution)
                     $sel move [transaxis x [expr {360 * rand()}]]
                     $sel move [transaxis y [expr {360 * rand()}]]
                     $sel move [transaxis z [expr {360 * rand()}]]
                     #measure unit cell size, create _ angstrom buffer between peptides
                     set xCellSize [expr {$unit_buffer}]
                     set yCellSize [expr {$unit_buffer}]
                     set zCellSize [expr {$unit_buffer}]
                     # set xCellSize [expr { [lindex [lindex $mm 1] 0] - [lindex [lindex $mm 0] 0] + $unit_buffer}]
                     # set yCellSize [expr { [lindex [lindex $mm 1] 1] - [lindex [lindex $mm 0] 1] + $unit_buffer}]
                     # set zCellSize [expr { [lindex [lindex $mm 1] 2] - [lindex [lindex $mm 0] 2] + $unit_buffer}]
                     #move unit cell
                     set xTrans [expr {$xOffset * $xCellSize + rand() - 0.5 }]
                     set yTrans [expr {$yOffset * $yCellSize + rand() - 0.5 }]
                     set zTrans [expr {$zOffset * $zCellSize + rand() - 0.5 }]
                     $sel moveby [list $xTrans $yTrans $zTrans]
                     ############################################################

                     # writes a temporary .pdb that will later be merged and adds it into vmd
                     $sel writepdb "moved_$xOffset$yOffset${zOffset}_L.pdb"
                     mol delete top
                     mol new "moved_$xOffset$yOffset${zOffset}_L.pdb"
                     mol addfile ${molname1}.psf
                     
                } elseif {$LorD == 1 && $Dcount < $Dratio} {
                     mol new ${molname2}.pdb
                     mol addfile ${molname2}.psf
                     set Dcount [expr $Dcount + 1]
                     puts $Dcount
                     puts "D"

                     set sel [atomselect top "segname PROA"]
                     $sel set segname PD$i

                     segment PD$i {
                     pdb ${molname2}.pdb
                     first NTER
                     last CT2
                     }

                     lappend LDwords PD$i

                     ############################################################
                     set sel [atomselect top all]
                     $sel moveby [vecinvert [measure center $sel]]

                     set mm [measure minmax $sel]
                     #rotate unit cell (with boltzmann distribution)
                     $sel move [transaxis x [expr {360 * rand()}]]
                     $sel move [transaxis y [expr {360 * rand()}]]
                     $sel move [transaxis z [expr {360 * rand()}]]
                     #measure unit cell size, create _ angstrom buffer between peptides
                     set xCellSize [expr {$unit_buffer}]
                     set yCellSize [expr {$unit_buffer}]
                     set zCellSize [expr {$unit_buffer}]
                     # set xCellSize [expr { [lindex [lindex $mm 1] 0] - [lindex [lindex $mm 0] 0] + $unit_buffer}]
                     # set yCellSize [expr { [lindex [lindex $mm 1] 1] - [lindex [lindex $mm 0] 1] + $unit_buffer}]
                     # set zCellSize [expr { [lindex [lindex $mm 1] 2] - [lindex [lindex $mm 0] 2] + $unit_buffer}]
                     #move unit cell
                     set xTrans [expr {$xOffset * $xCellSize + rand() - 0.5 }]
                     set yTrans [expr {$yOffset * $yCellSize + rand() - 0.5 }]
                     set zTrans [expr {$zOffset * $zCellSize + rand() - 0.5 }]
                     $sel moveby [list $xTrans $yTrans $zTrans]
                     ############################################################

                     $sel writepdb "moved_$xOffset$yOffset${zOffset}_D.pdb"
                     mol delete top
                     mol new "moved_$xOffset$yOffset${zOffset}_D.pdb"
                     mol addfile ${molname2}.psf
                     
                } elseif {$LorD == 1 && $Dcount == $Dratio} {
                     mol new ${molname1}.pdb
                     mol addfile ${molname1}.psf
                     set Lcount [expr $Lcount + 1]
                     puts $Lcount
                     puts "L"

                     set sel [atomselect top "segname PROA"]
                     $sel set segname PL$i

                     segment PL$i {
                     pdb ${molname1}.pdb
                     first NTER
                     last CT2
                     }

                     lappend LDwords PL$i

                     ############################################################
                     set sel [atomselect top all]
                     $sel moveby [vecinvert [measure center $sel]]

                     set mm [measure minmax $sel]
                     #rotate unit cell (with boltzmann distribution)
                     $sel move [transaxis x [expr {360 * rand()}]]
                     $sel move [transaxis y [expr {360 * rand()}]]
                     $sel move [transaxis z [expr {360 * rand()}]]
                     #measure unit cell size, create _ angstrom buffer between peptides
                     set xCellSize [expr {$unit_buffer}]
                     set yCellSize [expr {$unit_buffer}]
                     set zCellSize [expr {$unit_buffer}]
                     # set xCellSize [expr { [lindex [lindex $mm 1] 0] - [lindex [lindex $mm 0] 0] + $unit_buffer}]
                     # set yCellSize [expr { [lindex [lindex $mm 1] 1] - [lindex [lindex $mm 0] 1] + $unit_buffer}]
                     # set zCellSize [expr { [lindex [lindex $mm 1] 2] - [lindex [lindex $mm 0] 2] + $unit_buffer}]
                     #move unit cell
                     set xTrans [expr {$xOffset * $xCellSize + rand() - 0.5 }]
                     set yTrans [expr {$yOffset * $yCellSize + rand() - 0.5 }]
                     set zTrans [expr {$zOffset * $zCellSize + rand() - 0.5 }]
                     $sel moveby [list $xTrans $yTrans $zTrans]
                     ############################################################
                     
                     $sel writepdb "moved_$xOffset$yOffset${zOffset}_L.pdb"
                     mol delete top
                     mol new "moved_$xOffset$yOffset${zOffset}_L.pdb"
                     mol addfile ${molname1}.psf

                } elseif {$LorD == 0 && $Lcount == $Lratio} {
                     mol new ${molname2}.pdb
                     mol addfile ${molname2}.psf
                     set Dcount [expr $Dcount + 1]
                     puts $Dcount
                     puts "D"

                     set sel [atomselect top "segname PROA"]
                     $sel set segname PD$i

                     segment PD$i {
                     pdb ${molname2}.pdb
                     first NTER
                     last CT2
                     }

                     lappend LDwords PD$i

                     ############################################################
                     set sel [atomselect top all]
                     $sel moveby [vecinvert [measure center $sel]]

                     set mm [measure minmax $sel]
                     #rotate unit cell (with boltzmann distribution)
                     $sel move [transaxis x [expr {360 * rand()}]]
                     $sel move [transaxis y [expr {360 * rand()}]]
                     $sel move [transaxis z [expr {360 * rand()}]]
                     #measure unit cell size, create _ angstrom buffer between peptides
                     set xCellSize [expr {$unit_buffer}]
                     set yCellSize [expr {$unit_buffer}]
                     set zCellSize [expr {$unit_buffer}]
                     # set xCellSize [expr { [lindex [lindex $mm 1] 0] - [lindex [lindex $mm 0] 0] + $unit_buffer}]
                     # set yCellSize [expr { [lindex [lindex $mm 1] 1] - [lindex [lindex $mm 0] 1] + $unit_buffer}]
                     # set zCellSize [expr { [lindex [lindex $mm 1] 2] - [lindex [lindex $mm 0] 2] + $unit_buffer}]
                     #move unit cell
                     set xTrans [expr {$xOffset * $xCellSize + rand() - 0.5 }]
                     set yTrans [expr {$yOffset * $yCellSize + rand() - 0.5 }]
                     set zTrans [expr {$zOffset * $zCellSize + rand() - 0.5 }]
                     $sel moveby [list $xTrans $yTrans $zTrans]
                     ############################################################

                     $sel writepdb "moved_$xOffset$yOffset${zOffset}_D.pdb"
                     mol delete top
                     mol new "moved_$xOffset$yOffset${zOffset}_D.pdb"
                     mol addfile ${molname2}.psf
                     
                }

                # this allows for a list of peptides to be created and merged below
                set lastpep [molinfo top]
                lappend pepList $lastpep
                puts $pepList
                set i [expr $i + 1]
            }
        }
    }

    # merges the peptides into 1 pdb and psf
    package require topotools 1.0
    namespace import TopoTools::mergemols
    set newmol [mergemols "$pepList"]
    mol rename $newmol merged_${molname}
    set all [atomselect top all]
    $all update
    $all writepdb merged_${molname}.pdb

    # delete individual 64, load in .pdb merged
    mol delete all
    mol new merged_${molname}.pdb

    # creates a .psf for the new .pdb
    writepsf merged_$molname.psf
    mol delete all

    # load in the single new pdb and psf files
    mol new merged_${molname}.psf
    mol addfile merged_${molname}.pdb


    ##################################################
    # Puts the peptides in a water box               #
    ##################################################

    #sourced from wat_sphere.tcl in box

    # Determine the center of mass of the molecule and store the coordinates
    set cen [measure center [atomselect top all] weight mass]
    set x1 [lindex $cen 0]
    set y1 [lindex $cen 1]
    set z1 [lindex $cen 2]
    set max 0

    ### Determine the distance of the farthest atom from the center of mass
    foreach atom [[atomselect top all] get index] {
      set pos [lindex [[atomselect top "index $atom"] get {x y z}] 0]
      set x2 [lindex $pos 0]
      set y2 [lindex $pos 1]
      set z2 [lindex $pos 2]
      set dist [expr pow(($x2-$x1)*($x2-$x1) + ($y2-$y1)*($y2-$y1) + ($z2-$z1)*($z2-$z1),0.5)]
      if {$dist > $max} {set max $dist}
      }

    # Solvate the molecule in a water box with enough padding (4 A).
    # One could alternatively align the molecule such that the vector 
    # from the center of mass to the farthest atom is aligned with an axis,
    # and then use no padding
    resetpsf
    package require solvate
    solvate merged_${molname}.psf merged_${molname}.pdb -t $solvatePad -o del_water

    resetpsf
    package require psfgen
    mol new del_water.psf
    mol addfile del_water.pdb
    readpsf del_water.psf
    coordpdb del_water.pdb

    # # Determine which water molecules need to be deleted and use a for loop
    # # to delete them
    # set wat [atomselect top "same residue as {water and ((x-$x1)*(x-$x1) + (y-$y1)*(y-$y1) + (z-$z1)*(z-$z1))<($max*$max)}"]
    # set del [atomselect top "water and not same residue as {water and ((x-$x1)*(x-$x1) + (y-$y1)*(y-$y1) + (z-$z1)*(z-$z1))<($max*$max)}"]
    # set seg [$del get segid]
    # set res [$del get resid]
    # set name [$del get name]
    # for {set i 0} {$i < [llength $seg]} {incr i} {
    #   delatom [lindex $seg $i] [lindex $res $i] [lindex $name $i] 
    #   }
    writepsf ${molname}_ws.psf
    writepdb ${molname}_ws.pdb

    mol delete top

    mol new ${molname}_ws.psf
    mol addfile ${molname}_ws.pdb
    puts "CENTER OF MASS OF SPHERE IS: [measure center [atomselect top all] weight mass]"
    puts "RADIUS OF SPHERE IS: $max"


    # the 106 in the below line will change for peptides with more atoms 
    set todel [atomselect top "serial 1 to [expr $totalpep * 106]"]
    set seg [$todel get segid]
    set res [$todel get resid]
    set name [$todel get name]
    for {set i 0} {$i < [llength $seg]} {incr i} {
      delatom [lindex $seg $i] [lindex $res $i] [lindex $name $i] 
      }
    # writepsf to_merge.psf
    writepdb to_merge.pdb
    mol new to_merge.pdb


    # merges the peptides into 1 pdb and psf

    # POSSIBLE CHANGE, in "set newmol [mergemols {135 136}]" the numbers might have to be changed to whatever
    # your to_merge.pdb and molname_ws.pdb index in VMD is

    # from https://www.devdungeon.com/content/how-get-user-input-tcl
    # Present a message
    puts "Enter max mol ID: "

    # Get input from standard input (keyboard) and store in $someVar
    gets stdin maxID2

    set maxID1 [expr $maxID2 - 1]

    # Output the results
    puts "Merging $maxID1 and $maxID2"

    set mergeList {}

    lappend mergeList $maxID1
    lappend mergeList $maxID2

    package require topotools 1.0
    namespace import TopoTools::mergemols
    set newmol [mergemols $mergeList]    
    #set newmol [mergemols {"$maxID1" "$maxID2"}]
    mol rename $newmol merged_water
    set all [atomselect top all]
    $all update
    $all writepdb merged_water.pdb

    # load in the single pdb and psf files
    mol delete all
    mol new merged_water.pdb
    mol addfile ${molname}_ws.psf

    # was having troubles getting autoionize to work (issues with segment names)

    # ##################################################
    # # AutoIonization and gets cell dimensions        #
    # ##################################################

    # # uses autoionize package to add NaCl (default) to solvation box
    resetpsf
    package require autoionize

    # Neutralize the system and add an additional 150mM KCl:
    # autoionize -psf input.psf -pdb input.pdb -sc 0.15 -cation POT

    autoionize -psf ${molname}_ws.psf -pdb ${molname}_ws.pdb -sc 0.150 -cation SOD -anion CLA -o ${molname}final

    # THIS IS THE WORKING OLD LINE of code
    # autoionize -psf ${molname}_ws.psf -pdb ${molname}_ws.pdb -neutralize -o ${molname}final

    # also check if auto ionization still has net charge!

    # load in the single pdb and psf files
    mol delete all
    mol new ${molname}final.pdb
    mol addfile ${molname}final.psf

    # uses findbox to present 
    # LOOK WHAT THE ATOM SELECTION NEEDS TO BE
     set all [atomselect top all]
     set minmax [measure minmax $all]
     set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]
     puts "cellBasisVector1 [lindex $vec 0] 0 0"
     puts "cellBasisVector2 0 [lindex $vec 1] 0"
     puts "cellBasisVector3 0 0 [lindex $vec 2]"
     set center [measure center $all]
     puts "cellOrigin $center"
     set volume [expr  [lindex $vec 0]*[lindex $vec 1]*[lindex $vec 2]]
     puts "volume $volume"

    # Writes this data to the output file
     puts $outFile "CENTER OF MASS OF SPHERE IS: [measure center [atomselect top all] weight mass]"
     puts $outFile "RADIUS OF SPHERE IS: $max"

     puts $outFile "cellBasisVector1 [lindex $vec 0] 0 0"
     puts $outFile "cellBasisVector2 0 [lindex $vec 1] 0"
     puts $outFile "cellBasisVector3 0 0 [lindex $vec 2]"
     puts $outFile "cellOrigin $center"
     puts $outFile "volume $volume"

     puts $outFile  $LDwords
     puts $outFile "                                    speciesL  speciesD  outputname  topofile  ratioL  ratioD  totalpep  gridX  gridY  gridZ  unit_buffer  solvatePad"
     puts $outFile "The makefile functions inputs were $speciesL $speciesD $outputname $ratioL $ratioD $totalpep $gridX $gridY $gridZ $unit_buffer $solvatePad" 
     close $outFile
     puts "DONE!"

    # should yield files called ionized.pdb and .psf, please go and rename
    # these files to something identifying (eg L or D, ratio, experiment name...)
}
