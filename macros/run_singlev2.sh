#!/bin/bash
set -eo pipefail

PROCESSID=${1:?"usage: run_single.sh PROCESS_ID"}
NEVENTS=${NEVENTS:-10}
PARTICLE="e-"
MACRO_DIR=/sphenix/u/tomoki/work_tomoki/SiCalo/macros

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.532
set -u
export ROOT_INCLUDE_PATH=/sphenix/u/tomoki/work_tomoki/SiCalo/SiliconSeedAna/install/include:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=/sphenix/u/tomoki/work_tomoki/SiCalo/SiliconSeedAna/install/lib:$LD_LIBRARY_PATH

# Redirect stale ROOT/GenFit paths to CVMFS and supply zero-entry MVTX
# dead/hot pixel masks for the ideal MDC2 simulation.
export SICALO_MVTX_MASK_DIR=$MACRO_DIR
export LD_PRELOAD=$MACRO_DIR/sphenix_cdb_redirect.so:$MACRO_DIR/sphenix_stale_path_redirect.so${LD_PRELOAD:+:$LD_PRELOAD}

cd "$MACRO_DIR"
root -l -b -q "Fun4All_singleParticle_SiliconPositionv4.C(${PROCESSID},${NEVENTS},\"/sphenix/u/tomoki/work_tomoki/SiCalo/output\",\"${PARTICLE}\")"
