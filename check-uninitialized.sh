#!/bin/sh

#
# mbtrack binary
#
MBTRACK=./mbtrack-mpi-debug

WORK_PATH=work/tmp_`date +%Y-%m-%d`_job${PBS_JOBID}/

#
# mpirun binary
#
if [ "$MPIRUN" = "" ]
then
  MPIRUN=mpirun
fi

#
# valgrind binary
#
if [ -e $work/valgrind/bin/valgrind ]
then
  VALGRIND=$work/valgrind/bin/valgrind
else
  VALGRIND=valgrind
fi
$VALGRIND --help >/dev/null 2>&1 || { echo >&2 "Error: valgrind binary not found. For installation instruction go to http://valgrind.org/docs/manual/dist.install.html "; exit 1; }

#
# valgrind options
#
if [ -e /opt/intel/impi/3.2.2.006/etc64/valgrind.supp ]
then
  VALGRIND_SUPP="--suppressions=/opt/intel/impi/3.2.2.006/etc64/valgrind.supp"
fi
if [ -e /usr/share/openmpi/openmpi-valgrind.supp ]
then
  VALGRIND_SUPP="--suppressions=/usr/share/openmpi/openmpi-valgrind.supp"
fi
VALGRIND_MAX="--max-stackframe=208840972 --main-stacksize=64000000"
VALGRIND_OPTS="$VALGRIND_SUPP $VALGRIND_MAX --log-file=check-uninitialized_%p.log --track-origins=yes --leak-check=no --smc-check=none --vgdb=no"

#
# execution
#

mkdir "$WORK_PATH"
$MPIRUN -np 4 $MPIRUN_OPTS $VALGRIND $VALGRIND_OPTS $MBTRACK work/valgrind_initialized.conf "$WORK_PATH"
rm -rf "$WORK_PATH"
