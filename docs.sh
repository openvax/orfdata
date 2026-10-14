#!/bin/sh
set -eu
"${PYTHON:-python}" -m mkdocs build --strict
