#!/bin/bash
#
# compile_binaries_local.sh
#
# Build the CAT-Surface binaries for CAT on an Apple Silicon Mac with the same
# settings as .github/workflows/binaries.yml, so that the GitHub action is not
# needed for every build. The binaries are written into the build directories
# that "make cp_binaries" of CAT copies from:
#
#   -a  macOS arm64     native, runs on macOS >= 11.0       build-native-arm64
#   -x  macOS x86_64    native under Rosetta, macOS >= 10.13 build-native
#   -l  Linux x86_64    docker container manylinux_2_28 as in the CI,
#                       runs on glibc >= 2.14                build-x86_64-pc-linux
#   -w  Windows x86_64  mingw-w64 cross compiler, static    build-x86_64-w64-mingw32
#
# Without target options all four are built. Each build is checked like in the
# CI: minimum macOS version, glibc floor, and DLL dependencies on Windows.
#
# Requirements:
#   autoconf, automake, libtool    (for autoreconf)
#   Rosetta                        (-x: softwareupdate --install-rosetta)
#   docker                         (-l: e.g. OrbStack; use DOCKER_CONTEXT=orbstack
#                                   if the default context does not work)
#   mingw-w64 cross compiler       (-w: x86_64-w64-mingw32-gcc in $MINGW/bin)

SRC=$(cd "$(dirname "$0")/.." && pwd)
JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
MINGW=${MINGW:-/usr/local/x86_64-w64-mingw32}
IMAGE=${IMAGE:-quay.io/pypa/manylinux_2_28_x86_64}

MACOS_ARM64_TARGET=11.0
MACOS_X86_64_TARGET=10.13
GLIBC_FLOOR=2.14

# Compiler flags of the calling shell (e.g. -I/usr/local/include from the
# profile) must not reach the builds: they would mix headers and libraries of
# this Mac into the binaries (also into the Linux and Windows builds).
unset CFLAGS CPPFLAGS CXXFLAGS OBJCFLAGS LDFLAGS CPATH C_INCLUDE_PATH LIBRARY_PATH MACOSX_DEPLOYMENT_TARGET

targets=""
failed=""

usage ()
{
  sed -n '3,24p' "$0" | sed 's/^# \{0,1\}//'
  echo "Usage: $(basename "$0") [-a] [-x] [-l] [-w] [-j jobs]"
}

# print the tail of a log file after a failed step
fail_log ()
{
  echo "ERROR: $1 failed, see $2:"
  [ -f "$2" ] && tail -n 25 "$2"
}

# true if version $2 is not newer than version $1
not_newer ()
{
  [ "$(printf '%s\n%s\n' "$1" "$2" | sort -V | tail -n 1)" = "$1" ]
}

prepare ()
{
  echo "Running autoreconf in $SRC"
  ( cd "$SRC" && autoreconf -fi ) > /tmp/cat_autoreconf.$$.log 2>&1 || {
    fail_log autoreconf /tmp/cat_autoreconf.$$.log
    rm -f /tmp/cat_autoreconf.$$.log
    exit 1
  }
  rm -f /tmp/cat_autoreconf.$$.log
}

# build_mac <arch> <macOS target> <build dir>
build_mac ()
{
  local arch=$1 target=$2 dir=$SRC/$3 run="" bin minos n=0 status=0

  echo "--------------------------------------------------------"
  echo "macOS $arch (macOS >= $target) -> $3"
  echo "--------------------------------------------------------"
  if [ "$arch" = "x86_64" ]; then
    run="arch -x86_64"
    $run /usr/bin/true 2>/dev/null || { echo "ERROR: Rosetta is not installed"; return 1; }
  fi

  # The CI compiles with -std=gnu23, which its newer autoconf selects. The
  # autoconf on this Mac may not do that (and autoconf 2.71 would append
  # -std=gnu11), so -std=gnu23 is set here and the checks are skipped.
  rm -rf "$dir" && mkdir -p "$dir" && cd "$dir" || return 1
  $run /bin/bash -c "export MACOSX_DEPLOYMENT_TARGET=$target CC=\"clang -std=gnu23\" CFLAGS=-O2 CPPFLAGS=-O2
    ../configure --enable-static --disable-shared --disable-doc ac_cv_prog_cc_c11= ac_cv_prog_cc_c23= \
      > configure.log 2>&1" || { fail_log configure "$dir/configure.log"; return 1; }
  $run /bin/bash -c "export MACOSX_DEPLOYMENT_TARGET=$target
    make -j$JOBS > make.log 2>&1" || { fail_log make "$dir/make.log"; return 1; }

  # every binary has to run on the target macOS version
  for bin in "$dir"/Progs/CAT_*; do
    [ -x "$bin" ] && file "$bin" | grep -q 'Mach-O.*executable' || continue
    n=$((n+1))
    minos=$(otool -l "$bin" | awk '/LC_BUILD_VERSION/{b=1} b&&/minos/{print $2; exit} /LC_VERSION_MIN_MACOSX/{v=1} v&&/version/{print $2; exit}')
    if [ -z "$minos" ] || ! not_newer "$target" "$minos"; then
      echo "ERROR: $(basename "$bin") needs macOS ${minos:-?}, above $target"
      status=1
    fi
    $run "$bin" --help </dev/null >/dev/null 2>&1
    if [ $? -ge 126 ]; then
      echo "ERROR: $(basename "$bin") does not run"
      status=1
    fi
  done
  [ $n -gt 0 ] || { echo "ERROR: no binaries found in $dir/Progs"; return 1; }
  [ $status -eq 0 ] && echo "ok: $n binaries, all run on macOS >= $target"
  return $status
}

build_linux ()
{
  local dir=$SRC/build-x86_64-pc-linux

  echo "--------------------------------------------------------"
  echo "Linux x86_64 ($IMAGE, glibc >= $GLIBC_FLOOR) -> build-x86_64-pc-linux"
  echo "--------------------------------------------------------"
  command -v docker >/dev/null || { echo "ERROR: docker not found"; return 1; }
  docker version >/dev/null 2>&1 || {
    echo "ERROR: the docker engine does not answer (start it, or set DOCKER_CONTEXT, e.g. orbstack)"
    return 1
  }

  rm -rf "$dir" && mkdir -p "$dir" || return 1

  # The same steps, compiler, flags and glibc_compat.h as in the CI. autoreconf
  # has to run with the autotools of the container, because the autoconf
  # version decides on the C standard of gcc (e.g. -std=gnu23), which changes
  # the code. This is done in a copy of the sources inside the container, so
  # that the files generated for the macOS and Windows builds are not touched.
  # Only the binaries are copied back.
  cat > "$dir/docker_build.sh" <<'EOF'
#!/bin/bash
out=/src/build-x86_64-pc-linux
mkdir -p /tmp/src && cd /src &&
  tar --exclude='./build-*' --exclude='./.git' -cf - . | tar -xf - -C /tmp/src && cd /tmp/src || exit 1
command -v autoreconf >/dev/null && command -v make >/dev/null ||
  dnf install -y -q autoconf automake libtool make >/dev/null
echo "builder: $(gcc --version | sed -n 1p), $(ldd --version | sed -n 1p)"
autoreconf -fi > $out/autoreconf.log 2>&1 || { tail -n 25 $out/autoreconf.log; exit 1; }
export CFLAGS="-O2 -include $PWD/.github/glibc_compat.h" CPPFLAGS="-O2" LDFLAGS="-static-libgcc"
./configure --enable-static --disable-shared --disable-doc > $out/configure.log 2>&1 || { tail -n 25 $out/configure.log; exit 1; }
make -j$JOBS > $out/make.log 2>&1 || make -j1 > $out/make.log 2>&1 || { tail -n 25 $out/make.log; exit 1; }
mkdir -p $out/Progs && cp -p $(find Progs -maxdepth 1 -type f -executable -name 'CAT_*' ! -name '*_ui') $out/Progs/ || exit 1

# the glibc symbol versions a binary imports are a hard minimum
status=0; n=0
for f in $(find Progs -maxdepth 1 -type f -executable -name 'CAT_*' ! -name '*_ui' | sort); do
  n=$((n+1))
  highest=$(objdump -T "$f" |
            sed -nE 's/.*[ (]GLIBC_([0-9][0-9.]*)\)?[[:space:]]+([A-Za-z_][A-Za-z0-9_]*)$/\1/p' |
            sort -uV | tail -n 1)
  if [ -z "$highest" ]; then
    echo "ERROR: $(basename $f) imports no versioned glibc symbol"; status=1
  elif [ "$(printf '%s\n%s\n' "$FLOOR" "$highest" | sort -V | tail -n 1)" != "$FLOOR" ]; then
    echo "ERROR: $(basename $f) needs GLIBC_$highest, above GLIBC_$FLOOR"; status=1
  fi
done
[ $n -gt 0 ] || { echo "ERROR: no binaries found"; exit 1; }
[ $status -eq 0 ] && echo "ok: $n binaries, all run on glibc >= $FLOOR"
exit $status
EOF
  chmod +x "$dir/docker_build.sh"

  docker run --rm --platform linux/amd64 -e JOBS="$JOBS" -e FLOOR="$GLIBC_FLOOR" \
    -v "$SRC":/src "$IMAGE" /src/build-x86_64-pc-linux/docker_build.sh
}

build_windows ()
{
  local host=x86_64-w64-mingw32
  local dir=$SRC/build-$host bin deps n=0 status=0

  echo "--------------------------------------------------------"
  echo "Windows x86_64 ($MINGW, static) -> build-$host"
  echo "--------------------------------------------------------"
  [ -x "$MINGW/bin/$host-gcc" ] || { echo "ERROR: $MINGW/bin/$host-gcc not found (set MINGW)"; return 1; }

  # As in the CI (MSYS2), where configure finds no C23 support and keeps the
  # default C standard of gcc. ac_cv_prog_cc_c23=no keeps it like that with a
  # newer autoconf on this Mac.
  rm -rf "$dir" && mkdir -p "$dir" && cd "$dir" || return 1
  ( export PATH="$MINGW/bin:$PATH"
    ../configure --host=$host --enable-static --disable-shared --disable-doc CC=$host-gcc \
      CFLAGS="-O2" CPPFLAGS="-O2" LDFLAGS="-static -static-libgcc" ac_cv_prog_cc_c23=no > configure.log 2>&1
  ) || { fail_log configure "$dir/configure.log"; return 1; }
  ( export PATH="$MINGW/bin:$PATH"; make -j$JOBS > make.log 2>&1 ) || { fail_log make "$dir/make.log"; return 1; }

  # only system DLLs: the GCC, pthread and FFTW runtimes have to be linked in
  for bin in "$dir"/Progs/CAT_*.exe; do
    [ -f "$bin" ] || continue
    n=$((n+1))
    deps=$("$MINGW/bin/$host-objdump" -p "$bin" | sed -n 's/^[[:space:]]*DLL Name: //p' | tr '[:upper:]' '[:lower:]')
    if echo "$deps" | grep -qE 'fftw|gomp|winpthread|libgcc|libstdc'; then
      echo "ERROR: $(basename "$bin") depends on" $(echo "$deps" | grep -E 'fftw|gomp|winpthread|libgcc|libstdc')
      status=1
    fi
  done
  [ $n -gt 0 ] || { echo "ERROR: no binaries found in $dir/Progs"; return 1; }
  [ $status -eq 0 ] && echo "ok: $n binaries, only system DLLs"
  return $status
}

while [ $# -gt 0 ]; do
  case "$1" in
    -a) targets="$targets arm64" ;;
    -x) targets="$targets x86_64" ;;
    -l) targets="$targets linux" ;;
    -w) targets="$targets windows" ;;
    -j) JOBS=$2; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option $1"; usage; exit 1 ;;
  esac
  shift
done
[ -z "$targets" ] && targets="arm64 x86_64 linux windows"

prepare
for t in $targets; do
  case $t in
    arm64)   build_mac arm64  $MACOS_ARM64_TARGET  build-native-arm64 ;;
    x86_64)  build_mac x86_64 $MACOS_X86_64_TARGET build-native ;;
    linux)   build_linux ;;
    windows) build_windows ;;
  esac
  [ $? -ne 0 ] && failed="$failed $t"
  cd "$SRC"
done

echo "--------------------------------------------------------"
if [ -n "$failed" ]; then
  echo "FAILED:$failed"
  exit 1
fi
echo "All builds ok:$targets"
echo "Copy them to CAT with \"make cp_binaries\" in the CAT folder."
