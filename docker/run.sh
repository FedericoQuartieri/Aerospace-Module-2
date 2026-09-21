#!/bin/sh
# Opens a shell in the container with the project mounted on /project.
#
#   ./docker/run.sh                      interactive shell
#   ./docker/run.sh ./Allwmake           compiles and exits
#
# The project is a bind mount: whatever is compiled inside ends up in
# platforms/ on the Mac, and changes made on the Mac are visible inside at once.
#
# /root instead lives in a Docker volume, not on the bind mount: that is where
# OpenFOAM installs the user libraries ($FOAM_USER_LIBBIN), and without the
# volume they would disappear at every --rm, forcing a recompilation. The volume
# is also a native filesystem of the VM, much faster than the one shared with macOS.
cd ${0%/*}/.. || exit 1

if [ -t 0 ]; then tty="-it"; else tty=""; fi

exec docker run --rm $tty \
    -v "$(pwd)":/project \
    -v aero-m2-home:/root \
    aero-m2:dev "${@:-bash}"
