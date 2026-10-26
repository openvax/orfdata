#!/bin/sh
set -eu
"${PYTHON:-python}" -m ruff check orfdata tests tools
