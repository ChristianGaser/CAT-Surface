#! /bin/sh

test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.
autoreconf --force --install --verbose "$srcdir"
test -n "$NOCONFIGURE" || "$srcdir/configure" "$@"

exit
libtoolize --automake --force --copy
aclocal -I m4
autoheader
automake --add-missing --copy
autoconf
