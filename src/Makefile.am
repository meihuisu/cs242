# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${PROJ_INCLUDES}
AM_LDFLAGS = ${LDFLAGS} ${PROJ_LDFLAGS}


TARGETS = libcs242.a libcs242.so 

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcs242.so ${prefix}/lib
	cp libcs242.a ${prefix}/lib
	cp cs242.h ${prefix}/include

libcs242.a: cs242_static.o
	$(AR) rcs $@ $^

libcs242.so: cs242.o 
	$(CC) -shared $(AM_CCFLAGS) -o libcs242.so $^ $(AM_LDFLAGS)

cs242.o: cs242.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

cs242_static.o: cs242.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	

clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

