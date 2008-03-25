#! /bin/sh

set -e
 
glibtoolize --automake --copy
aclocal -I m4
autoheader
automake --add-missing --copy
autoconf
