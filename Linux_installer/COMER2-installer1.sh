#!/bin/bash

echo

MYHOMEDEF=${HOME}/local/comer2

read -ep "Enter COMER2 install path: " -i "${MYHOMEDEF}" MYHOME

echo
echo Install path: ${MYHOME}
echo

srcdir="$(dirname $0)"

[[ -d "${srcdir}/bin" && -d "${srcdir}/var" ]] || (echo "ERROR: Source directories not found!" && exit 1)

[ -f "${srcdir}/bin/comer" ] || (echo "ERROR: Incomplete software package: executable comer missing!" && exit 1)

mkdir -p "${MYHOME}" || (echo "ERROR: Failed to create destination directory!" && exit 1)

[ -d "${MYHOME}/bin" ] && (rm -fR "${MYHOME}/bin" || exit 1)
[ -d "${MYHOME}/var" ] && (rm -fR "${MYHOME}/var" || exit 1)

cmd="${srcdir}/bin/comer --dev-list >/dev/null 2>&1"
eval ${cmd} || cat <<EOF

WARNING: The comer executable will not run on the system!
WARNING: Please make sure a GPU and appropriate NVIDIA and 
WARNING: CUDA drivers are installed. 
WARNING: If they are installed, please compile and install
WARNING: the software from the source code by typing:
WARNING: BUILD_and_INSTALL_unix.sh (gcc compiler) or
WARNING: BUILD_and_INSTALL_unix__clang.sh (clang compiler).

EOF

cp -R "${srcdir}/bin" "${MYHOME}/" || (echo "ERROR: Failed to install the package!" && exit 1)
cp -R "${srcdir}/var" "${MYHOME}/" || (echo "ERROR: Failed to install the package!" && exit 1)

echo Installation complete.

exit 0

