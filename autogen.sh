#! /bin/sh

aclocal

# Support for Linux and Mac OS X
# libtoolize is given the name glibtoolize
# when installing GNU libtool via 
# Homebrew or MacPorts on Mac
if hash libtoolize 2>&-
then
    libtoolize --automake
else
    glibtoolize --automake
fi

automake --add-missing
autoconf
