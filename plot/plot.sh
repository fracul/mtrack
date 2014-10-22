#!/bin/bash

# Use default path if PLOT_PATH is not set
if [ "$PLOT_PATH" = "" ]
then
  export PLOT_PATH=`pwd`/plot/
fi

# Return error if directory PLOT_PATH doesn't exist
if [ ! -d "$PLOT_PATH" ]; then
  echo "Cannot find 'plot' directory. Please define environment variable PLOT_PATH"
  exit -1
fi


function gnuplot_from_file
{
  if [ -f "$1" ]
  then
    echo "Plotting $1"
    GNUPLOT_FILE=$PLOT_PATH/$(basename "$1" "$2").gnuplot
    $GNUPLOT $GNUPLOT_FILE
  else
    echo "ERROR: data file \"$1\" not found"
  fi
}

if [ -f ampinv_VER_mean.dat ]
then
  gnuplot_from_file "ampinv_VER_mean.dat" ".dat"
fi

if [ -f ampinv_HOR_mean.dat ]
then
  gnuplot_from_file "ampinv_HOR_mean.dat" ".dat"
fi
