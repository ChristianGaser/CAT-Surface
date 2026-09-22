#! /bin/sh
#
# Generate the Doxygen documentation into docs/doxygen/.
#
# The HTML entry point is docs/doxygen/html/index.html and the warnings are
# written to docs/doxygen/warnings.log (both set in the Doxyfile).  Run from
# anywhere; `make docs` calls this script.

set -e

cd "$(dirname "$0")/.."

if ! command -v doxygen >/dev/null 2>&1; then
    echo "doxygen not found; install it first (brew/apt-get install doxygen)." >&2
    exit 1
fi

mkdir -p docs/doxygen
doxygen Doxyfile

log=docs/doxygen/warnings.log
if [ -f "$log" ]; then
    echo "Documentation written to docs/doxygen/html/index.html"
    echo "$(grep -c 'warning:' "$log") warnings in $log"
fi
