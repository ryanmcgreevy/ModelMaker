set f [open "bad_clashes.txt" r]
set data [read $f]
close $f
set badclashes [split $data "\n"]



# start_rosetta_refine {jobname mol selections anchor cartesian mapname mapresolution score_dens bestN nstruct {cluster 0} {nPerTask 5} {scoreOnly 0} args}