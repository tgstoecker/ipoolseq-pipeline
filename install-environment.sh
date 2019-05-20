#!/bin/bash

# TODO: Deal with a non-LFS checkout (or GitHub archive download) where environment.tar.gz
# is only a LFS pointer to the actual blob

if [ -e "./environment" ] && ! [ -z "$(ls -A ./environment)"]; then
	echo "./environment alreay exists and is non-emtpty. If you want to replace it, remove the old directory first!" >&2
	exit 1
fi

echo "*** Unpacking environment.tar.gz"
test -e ./environment || mkdir ./environment || exit 1
tar -xzf environment.tar.gz -C environment || exit 1

echo "*** Cleanup prefixes in ./environment"
source ./environment/bin/activate
conda-unpack || exit 1

echo "*** DONE"
echo "The environment is now usable and can be activated with"
echo ""
echo "  source ./environment/bin/activate"
echo ""
echo "The environment must be re-activated in every terminal session"
