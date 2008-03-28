#! /bin/sh

set -e
 
libtoolize --automake --copy
aclocal -I m4
autoheader
automake --add-missing --copy
autoconf
