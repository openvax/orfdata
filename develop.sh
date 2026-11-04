#!/bin/sh
set -eu
"${PYTHON:-python}" -m pip install -e '.[dev,docs]'
