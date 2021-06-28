#!/bin/bash
PROJECTDIR="/usr/local/AGFORTELNY/PROJECTS/Neuroblastoma"
PROJECTDIR_SIF="/media/AGFORTELNY"

singularity run \
  --bind $PROJECTDIR:$PROJECTDIR_SIF \
  scvelo.sif \
  python analyse_velocity.py
