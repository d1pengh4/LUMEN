#!/usr/bin/env bash
# Wrapper for the Colab CLI that fixes the macOS SSL CA issue.
#
# Without this, the kernel WebSocket (wss://...prod.colab.dev) fails with
#   [SSL: CERTIFICATE_VERIFY_FAILED] unable to get local issuer certificate
# because the macOS Python.framework has no CA bundle configured for the
# websocket-client library. We point SSL_CERT_FILE at certifi's bundle.
#
# Usage:  src/generation/colabx.sh run --gpu T4 --timeout 2400 <script.py> [args]
#         src/generation/colabx.sh exec -s lumen -f <script.py> --timeout 2400
set -euo pipefail

export PATH="/Library/Frameworks/Python.framework/Versions/3.12/bin:$PATH"
CERT="$(uv tool run --from google-colab-cli python -c 'import certifi; print(certifi.where())' 2>/dev/null | tail -1)"
export SSL_CERT_FILE="$CERT"
export SSL_CERT_DIR="$(dirname "$CERT")"
exec colab "$@"
