#!/usr/bin/env bash
# Publish only a reviewed version from clean, synchronized main.
set -euo pipefail
PYTHON="${PYTHON:-python}"
[[ "$(git branch --show-current)" == main ]] || { echo 'Release from main only' >&2; exit 1; }
[[ -z "$(git status --porcelain)" ]] || { echo 'Release requires a clean tree' >&2; exit 1; }
git fetch origin main
[[ "$(git rev-parse HEAD)" == "$(git rev-parse origin/main)" ]] || { echo 'Synchronize main first' >&2; exit 1; }
version="$($PYTHON -c 'from orfdata.version import __version__; print(__version__)')"
if git rev-parse "v$version" >/dev/null 2>&1; then echo 'Version tag already exists' >&2; exit 1; fi
./lint.sh
./test.sh
"$PYTHON" -m build
"$PYTHON" -m twine check "dist/orfdata-$version.tar.gz" "dist/orfdata-$version-py3-none-any.whl"
if [[ "${1:-}" == --dry-run ]]; then exit 0; fi
"$PYTHON" -m twine upload "dist/orfdata-$version.tar.gz" "dist/orfdata-$version-py3-none-any.whl"
git tag "v$version"
git push origin "v$version"
