#!/usr/bin/env bash
set -euo pipefail

URL="https://pgda.gsfc.nasa.gov/data/LOLA_20mpp/LDEM_80S_80MPP_ADJ.TIF"
OUT_FILE="LDEM_80S_80MPP_ADJ.TIF"

echo "Downloading LOLA LDEM file..."
echo "URL: ${URL}"
echo "Output: ${OUT_FILE}"

if command -v curl >/dev/null 2>&1; then
    curl -L --fail --continue-at - \
        --output "${OUT_FILE}" \
        "${URL}"
elif command -v wget >/dev/null 2>&1; then
    wget -c \
        -O "${OUT_FILE}" \
        "${URL}"
else
    echo "Error: neither curl nor wget is installed." >&2
    exit 1
fi

echo "Download completed:"
ls -lh "${OUT_FILE}"
