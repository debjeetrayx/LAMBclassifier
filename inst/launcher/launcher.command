#!/bin/bash
# LAMB classifier launcher (macOS). Double-click to choose a file or folder.
# Written by LAMBclassifier::lamb_install_launcher(). Research use only.
RSCRIPT="{{RSCRIPT}}"

echo "LAMB classifier - research use only, not a clinical diagnosis."
echo

if [ ! -x "$RSCRIPT" ]; then
    echo "R was not found at:"
    echo "  $RSCRIPT"
    echo "R may have been updated or moved. Open R and run:"
    echo "  LAMBclassifier::lamb_install_launcher()"
    echo
    read -n 1 -s -r -p "Press any key to close this window."
    exit 1
fi

export LAMB_LAUNCHER=1
if [ $# -ge 1 ]; then
    LAMB_INPUT="$1" "$RSCRIPT" -e 'invisible(suppressWarnings(suppressPackageStartupMessages(library(LAMBclassifier)))); lamb_run(Sys.getenv("LAMB_INPUT"))'
else
    "$RSCRIPT" -e 'invisible(suppressWarnings(suppressPackageStartupMessages(library(LAMBclassifier)))); lamb_run()'
fi
status=$?

echo
if [ $status -ne 0 ]; then
    echo "LAMB classifier stopped with an error (see the message above)."
    echo "Nothing was changed in your data."
else
    echo "Done."
fi
read -n 1 -s -r -p "Press any key to close this window."
echo
exit $status
