all: mbtrack-mpi

debug: mbtrack-mpi-debug

clean:
	cd src/; $(MAKE) clean
	cd doc/; $(MAKE) clean

mbtrack-mpi: FORCE
	cd src/; $(MAKE)

mbtrack-mpi-debug: FORCE
	cd src/; $(MAKE) debug

test: FORCE
	cd src/; $(MAKE) test

doc: doc/html

doc/html: FORCE
	cd doc/; $(MAKE) html

FORCE:
