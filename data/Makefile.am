# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

install:
	mkdir -p ${prefix}/data
	cp ./config ${prefix}/data
	cp -R ./i* ${prefix}/data
	chmod +r ${prefix}/data/config
	chmod +r ${prefix}/data/cs242/*
