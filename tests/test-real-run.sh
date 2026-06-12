#!/bin/sh

set -eu

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)

if [ "${TAFFISH_RNASEQ_REAL_OUT:-}" != "" ] && [ "${TAFFISH_RNASEQ_REAL_REFERENCE_OUT:-}" = "" ]; then
    TAFFISH_RNASEQ_REAL_REFERENCE_OUT="$TAFFISH_RNASEQ_REAL_OUT"
    export TAFFISH_RNASEQ_REAL_REFERENCE_OUT
fi

exec "$script_dir/test-real-reference.sh"
