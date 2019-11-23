#!/bin/bash

echo

MYHOMEDEF=${HOME}/local/comer2

read -ep "Enter COMER2 install path: " -i "${MYHOMEDEF}" MYHOME

echo
echo Install path: ${MYHOME}
echo

srcdir="$(dirname $0)"

[[ -d "${srcdir}/bin" && -d "${srcdir}/var" ]] || (echo "ERROR: Source directories not found!" && exit 1)

mkdir -p "${MYHOME}" || (echo "ERROR: Failed to create destination directory!" && exit 1)

[ -d "${MYHOME}/bin" ] && (rm -fR "${MYHOME}/bin" || exit 1)
[ -d "${MYHOME}/var" ] && (rm -fR "${MYHOME}/var" || exit 1)

cp -R "${srcdir}/bin" "${MYHOME}/" || (echo "ERROR: Failed to install the package!" && exit 1)
cp -R "${srcdir}/var" "${MYHOME}/" || (echo "ERROR: Failed to install the package!" && exit 1)

echo Installation complete.

exit 0

