#!/bin/sh

rsync -r --include='*/*.h' --include '*/*.out' --exclude='*/*' ../output/ .
bzr add */*.h */*.out
