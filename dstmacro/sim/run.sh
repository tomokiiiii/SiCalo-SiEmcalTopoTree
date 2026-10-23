#!/bin/bash
set -eo pipefail

PROCESSID=${1:?process ID is required}
PROCESS6=$(printf "%06d" "${PROCESSID}")
NEVENTS=${NEVENTS:-0}
OUTPUT=/sphenix/u/tomoki/work_tomoki/SiCalo/output/ana_pythia_sicalo_position_v5/SiEmcalTopoTree_0000000029-${PROCESS6}_v5.root

if [[ -s "${OUTPUT}" ]]; then
  echo "existing v5 output (skip): ${OUTPUT}"
  exit 0
fi
WORKDIR=/sphenix/u/tomoki/work_tomoki/SiCalo/dstmacro/sim
REDIRECT_DIR=/sphenix/u/tomoki/work_tomoki/SiCalo/macros

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.532
export ROOT_INCLUDE_PATH=/sphenix/u/tomoki/work_tomoki/SiCalo/SiliconSeedAna/install/include:${ROOT_INCLUDE_PATH:-}
export LD_LIBRARY_PATH=/sphenix/u/tomoki/work_tomoki/SiCalo/SiliconSeedAna/install/lib:${LD_LIBRARY_PATH:-}

# Some ana.532 ROOT dictionaries contain absolute build-time header paths on
# retired GPFS storage. Redirect those stale paths to the active CVMFS GCC/
# GenFit headers before ROOT starts; otherwise Cling exits with code 1 and the
# Condor job is held by OnExitHold.
export LD_PRELOAD=${REDIRECT_DIR}/sphenix_stale_path_redirect.so${LD_PRELOAD:+:${LD_PRELOAD}}

cd "${WORKDIR}"
root -l -b -q "Fun4All_PYTHIA_SiliconRecoPositionv5.C(${PROCESSID},${NEVENTS})"
