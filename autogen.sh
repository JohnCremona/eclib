#! /bin/sh

aclocal \
&& libtoolize \
&& automake --add-missing \
&& autoconf
