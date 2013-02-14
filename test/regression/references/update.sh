#!/bin/sh

rsync -r --include='*.h' --exclude='*' ../output/r_auto/ r_auto/
rsync -r --include='*.h' --exclude='*' ../output/r_quadrature/ r_quadrature/
rsync -r --include='*.h' --exclude='*' ../output/r_quadrature_O/ r_quadrature_O/
rsync -r --include='*.out' --exclude='*' ../output/r_auto/ output/
rsync -r --include='*.json' --exclude='*' ../output/r_auto/ output/

bzr add output/*.out
bzr add output/*.json
bzr add r_auto/*.h
bzr add r_quadrature/*.h
bzr add r_quadrature_O/*.h
