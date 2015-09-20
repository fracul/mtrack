########################################################
#                                                      #
#                      make_mbtrack                    #
#                                                      #
########################################################
#
#   cf)  A Makefile to create mbtrack:
#
############################################################################
#  <<Physics>>    (08/02/2002)                Last Modified: May 2011 
############################################################################
#
#===========================================================================
#

CFLAGS=-Wall -O2
CC=gcc
MPICC=mpicc -DDEBUG_INPUT
LIBS=-lm

WEAKWEAK_OBJECTS=manager_weakweak.o worker_weakweak.o bunch_weak.o transform_weak.o cyclic_array.o feedback_rf.o
WEAKSTRONG_OBJECTS=manager_weakstrong.o worker_weakstrong.o bunch_strong.o fbi_weakstrong.o rw_weakstrong.o

all : ../mbtrack-mpi

debug : MPICC += -DDEBUG -g
debug : CFLAGS = -Wall -Wextra -O0 
debug : ../mbtrack-mpi-debug

debug-log : MPICC += -DDEBUG_INPUT -DDEBUG_LOG
debug-log : debug

test : MPICC += -g
test : CFLAGS = -Wall -Wextra -O0
test : test_complexerrf test_cyclic_array test_histogram test_input test_reglin test_weak_bunch
	@for t in $^ ; do echo -n "Running $$t... " ; ./$$t ; done

#===========================================================================
# Linking
#===========================================================================

../mbtrack-mpi : mbtrack-mpi.o confmpi.o input.o test.o ${WEAKWEAK_OBJECTS} ${WEAKSTRONG_OBJECTS} bunch.o stats.o physic.o statistics.o tracking.o rw_table.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

../mbtrack-mpi-debug : mbtrack-mpi.o confmpi.o input.o test.o ${WEAKWEAK_OBJECTS} ${WEAKSTRONG_OBJECTS} bunch.o stats.o physic.o statistics.o tracking.o  rw_table.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

test_complexerrf : test.o test_complexerrf.o fbi_weakstrong.o physic.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

test_cyclic_array : test.o test_cyclic_array.o cyclic_array.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

test_histogram : test.o test_histogram.o statistics.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

test_input : test.o test_input.o input.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

test_reglin : test.o test_reglin.o physic.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

test_rw: test.o test_rw.o input.o physic.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

test_weak_bunch: test.o test_weak_bunch.o bunch_weak.o bunch.o bunch_strong.o statistics.o
	${MPICC} $(CFLAGS) $^ -o $@ $(LIBS)

# --------------------------------------------------------------------------
#  Building
# --------------------------------------------------------------------------

# objects building
%.o: %.c
	$(MPICC) -c $(CFLAGS) $< -o $@

# clean target
clean:
	rm -f *.o *.[ch]~

##############################################################################
#                                                                            #
#                        End of the Makefile: mbtrack                        #              	 
#                                                                            #
##############################################################################
