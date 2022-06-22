#!/bin/bash

# Directories
DIRRTROVER=RTRover
DIRRTKLIB=rtklib
DIRPPPTOOLS=`pwd`

# Options
CC=g++
CFLAGS="-O3 -Wl,--stack,4194304 -DRTRover_STATIC_LIB -lm -lWs2_32 -lWinmm"
# rtklib generation
generateRtklib()
{
  echo "*****************"
  echo "RTKLIB generation"
  echo "*****************"
  cd $DIRRTKLIB ; make ; cd -
}

# RTRover generation
generateRTRover()
{
  echo "******************"  
  echo "RTRover generation"
  echo "******************"
  cd $DIRRTROVER ; make ; cd -
}

# getStream generation
generateGetStream()
{
  echo "********************"
  echo "getStream generation"
  echo "********************"
  $CC getStream.cpp -I$DIRRTKLIB $DIRRTKLIB/rtklib.a $CFLAGS -o getStream
}

# processStream generation
generateProcessStream()
{
  echo "************************"
  echo "processStream generation"
  echo "************************"
# Pour la livraison
 $CC processStream.cpp -I$DIRRTKLIB -I$DIRRTROVER $DIRRTROVER/libRTRover.a $DIRRTKLIB/rtklib.a $CFLAGS -DENAGLO -DENAGAL -DNFREQ=4 -DNEXOBS=1 -o processStream

}

# generateLowLevel generation
generateLowLevel()
{
  echo "***************************"
  echo "generateLowLevel generation"
  echo "***************************"
 $CC generateLowLevel.cpp -I$DIRRTKLIB $DIRRTKLIB/rtklib.a $CFLAGS -DENAGLO -DENAGAL -DNFREQ=4 -DNEXOBS=1 -o generateLowLevel
}

# processLowLevel generation
generateProcessLowLevel()
{
  echo "**************************"
  echo "processLowLevel generation"
  echo "**************************"	
 $CC processLowLevel.cpp -I$DIRRTROVER $DIRRTROVER/libRTRover.a $CFLAGS -o processLowLevel
}

# generation of all libraries, getStream and processStream
generateAll()
{
  generateRtklib
  generateRTRover
  generateGetStream
  generateProcessStream
  generateLowLevel
  generateProcessLowLevel
}

# Main program

# Too many options
if [ $# -gt 1 ] ; then   
  echo "Error ! Usage : compile.sh [ clean | getStream | processStream | generateLowLevel | processLowLevel | computeErrors | sendStream | all ]"
  echo "Examples :"
  echo "  ./compile.sh"
  echo "  ./compile.sh clean"
  exit
fi

# Read option
if [ $# -eq 1 ] ; then
  if [ $1 == "clean" ] ; then
    cd $DIRRTKLIB ; make clean ; cd $DIRPPPTOOLS ;
    cd $DIRRTROVER ; make clean ; cd $DIRPPPTOOLS ;
    rm -f processStream getStream processLowLevel computeErrors sendStream generateLowLevel
  elif [ $1 == "processStream" ] ; then
    generateRtklib
    generateRTRover
    generateProcessStream
  elif [ $1 == "getStream" ] ; then
    generateRtklib
    generateGetStream
  elif [ $1 == "generateLowLevel" ] ; then
    generateRtklib
    generateLowLevel
  elif [ $1 == "processLowLevel" ] ; then
    generateRTRover
    generateProcessLowLevel
  elif [ $1 == "computeErrors" ] ; then
    generateRtklib
    generateComputeErrors
    elif [ $1 == "sendStream" ] ; then
    generateRtklib
    generateSendStream
  elif [ $1 == "all" ] ; then
    generateAll
  else
    # Bad option
    echo "Error ! Usage : compile.sh [ clean | getStream | processStream | processLowLevel | generateLowLevel | computeErrors | sendStream | all ]"
    echo "Examples :"
    echo "  ./compile.sh"
    echo "  ./compile.sh clean"
    exit
  fi
fi

# No option => all
if [ $# -eq 0 ] ; then
    generateAll
fi
