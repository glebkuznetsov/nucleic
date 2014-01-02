#!/bin/bash
# installprimer3unix.sh: compiles and installs the Primer3 library on UNIX
#                        systems (linux / os x).
#

install_p3(){
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    PDIR="$DIR/primer3-2.3.6"
    SDIR="$PDIR/src"
    tar -xvf "$DIR/primer3-src-2.3.6.tar.gz"
    make -C "$SDIR"
    cp $SDIR/long_seq_tm_test $SDIR/ntdpal $SDIR/ntthal $SDIR/oligotm \
       $SDIR/primer3_core /usr/local/bin/
    cp -R $SDIR/primer3_config /opt/primer3_config
    rm -rf $PDIR
}

check_p3_install(){
    ntthal -a HAIRPIN -s1 GGGGGGGGGGCCCCCCCCCC >/dev/null 2>&1
}

check_p3_install
if [ $? -eq 0 ]; then
    IDIR="$(which primer3_core)"
    IDIR=${IDIR%/*}
    echo "Primer3 is already installed in $IDIR"
else
    install_p3
    check_p3_install
    if [ $? -eq 0 ]; then
        echo "Primer3 successfully installed"
    else
        echo "Error installing primer3 (try running this script with sudo)"
    fi
fi
