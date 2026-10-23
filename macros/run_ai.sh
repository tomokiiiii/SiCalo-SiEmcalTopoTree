#!/bin/bash
set -eo pipefail

PROCESSID=$(printf "%06d" "${1:?process ID is required}")
BASE=/sphenix/u/tomoki/work_tomoki/SiCalo/output/ana_pythia_sicalo_v2
MACRODIR=/sphenix/u/tomoki/work_tomoki/SiCalo/macros
INPUT=${BASE}/ana_0000000029-${PROCESSID}.root
OUTPUT=${BASE}/SiEmcalTopoTree_0000000029-${PROCESSID}_v4.root
TEMP_OUTPUT=${OUTPUT}.tmp.$$

if [[ ! -s "${INPUT}" ]]; then
  echo "missing input: ${INPUT}" >&2
  exit 2
fi
if [[ -s "${OUTPUT}" ]]; then
  echo "existing v4 output (skip): ${OUTPUT}"
  exit 0
fi

cleanup() { rm -f "${TEMP_OUTPUT}"; }
trap cleanup EXIT

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.532
export ROOT_INCLUDE_PATH=/sphenix/u/tomoki/work_tomoki/SiCalo/SiliconSeedAna/install/include:${ROOT_INCLUDE_PATH:-}
export LD_LIBRARY_PATH=/sphenix/u/tomoki/work_tomoki/SiCalo/SiliconSeedAna/install/lib:${LD_LIBRARY_PATH:-}

cd "${MACRODIR}"

# Produce only the final v4 file. The intermediate flat tree uses a temporary
# name and is moved into place only after all v4 branches were added.
root -l -b -q \
  "make_SiEmcalTopoTree_v3.C(0,\"${INPUT}\",\"${TEMP_OUTPUT}\",true)"
test -s "${TEMP_OUTPUT}"

root -l -b -q \
  -e ".L ${MACRODIR}/add_SiEmcalTopoTreeFeatures_v4.C" \
  -e "const bool ok = SiEmcalTopoV4::AddV4Branches(\"${INPUT}\",\"${TEMP_OUTPUT}\"); gSystem->Exit(ok ? 0 : 3);"

mv "${TEMP_OUTPUT}" "${OUTPUT}"
trap - EXIT
echo "saved v4 output: ${OUTPUT}"
