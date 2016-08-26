all: mbtrack-mpi

cuda : mbtrack-cuda

debug: mbtrack-mpi-debug

clean:
	cd src/; $(MAKE) clean
	cd doc/; $(MAKE) clean

mbtrack-mpi: FORCE
	cd src/; $(MAKE)

mbtrack-cuda: FORCE
	cd src/; $(MAKE) cuda

mbtrack-mpi-debug: FORCE
	cd src/; $(MAKE) debug

test: FORCE
	cd src/; $(MAKE) test

doc: doc/html

doc/html: FORCE
	cd doc/; $(MAKE) html

remake:
	cd src/; $(MAKE) clean
	cd src/; $(MAKE) 
	cd src/; $(MAKE) clean
	cd src/; $(MAKE) cuda

FORCE:
