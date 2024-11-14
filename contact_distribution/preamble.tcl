mol new ../../RNA_Testset/PDB/4frg.pdb
mol addrep top

display cuemode Exp2
display cuedensity 0.18
display ambientocclusion on
display aoambient 0.75
display aodirect 0.75
display shadows on
display depthcue off

mol modcolor 0 top Chain
mol modcolor 1 top Chain
material change specular Ghost 0.050000
mol modmaterial 1 top AOChalky
mol modmaterial 0 top AOChalky
mol modselect 0 top not backbone
mol modstyle 0 top NewCartoon 0.400000 10.000000 4.100000 0
mol modstyle 1 top Tube 2.300000 12.000000

graphics top materials on
graphics top material AOChalky
axes location Off

color change rgb 0 0.01 0.24 0.42; #juelich blue
color change rgb 4 0.98 0.92 0.35; #juelich yellow
color change rgb 1 0.92 0.37 0.45; #juelich red
color change rgb 7 0.73 0.82 0.37; #juelich green
color change rgb 23 0.68 0.74 0.89; #juelich blue2
color change rgb 2 0.92 0.92 0.92; #juelich gray
color change rgb 19 0 1 0; #green2
color change rgb 29 1 0 0; #red2
color change rgb 25 0.69 0.51 0.73; #juelich violet
color Resname G blue
color Resname U red
color Resname C green
color Resname A yellow
color Chain A blue

proc geom_center {selection} {
    set gc [veczero]
    foreach coord [$selection get {x y z}] {
        set gc [vecadd $gc $coord]
    }
    return [vecscale [expr 1.0 /[$selection num]] $gc]
}
