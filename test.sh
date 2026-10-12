#!/bin/sh
set -eu
"${PYTHON:-python}" -m pytest tests "$@"
