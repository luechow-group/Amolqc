#! /bin/bash
echo ""
echo ""
echo "                 AMOLQC tests"
echo ""
FAILED=0
PAR=0
GEN=0
if [ -z "$SED" ] ;then
  if [[ "$OSTYPE" == "linux-gnu" ]]; then
    # Linux
    export SED=sed
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    # MacOs
    if [ -z $(which gsed) ]; then
      echo "On MacOS, gnu-sed is needed to run the testsuite: 'brew install gnu-sed'"
      exit 1
    else
      export SED=gsed
    fi
  else
    echo  -e "Amolqc/testsuite/tests.sh: operating system $OSTYPE not supported"
    exit 1
  fi
fi
export REF=ref
if [ -z "$AMOLQC_TEST" ] ;then
export AMOLQC=${PWD%testsuite/run}
else
export AMOLQC=$AMOLQC_TEST
fi
if [ -z "$AMOLQCRUN" ] ;then
export AMOLQCRUN=$AMOLQC/bin/amolqc
fi
selection=$2
if [ "$1" == "PARALLEL" ] ;then
export PAR=1
echo " Parallel test (2 CPUs)"
else
selection=$1
fi
if [ $PAR == 1 ] ;then
  export AMOLQCRUN="mpirun -n 2 $AMOLQCRUN"
  export REF=ref.par
fi

echo " Current version information: "
echo "    installation folder = "  $AMOLQC
a=`cd $AMOLQC && git describe --abbrev=6 --dirty --always --tags`
echo "    version            = " $a
if [ $GEN == 1 ] ; then
  echo " updating the references ..."
  echo "version: " $a >reference.txt
  echo "date:" `date` >>reference.txt
fi
echo ""
fcount=0
n=1
if [[ -z $selection ]]; then
  for dir in */; do fcount=$[$fcount+1] ;done
  for dir in */
  do
    cd $dir
    echo " starting test" $n "out of" $fcount
    source RUN
    n=$[$n+1]
    echo ""
    cd ..
  done
else
  for dir in */
  do
    if [[ $dir =~ $selection ]] ; then
    fcount=$[$fcount+1]
    fi
  done

  for dir in */
  do
    if [[ $dir =~ $selection ]] ; then
      cd $dir
      echo " starting test" $n "out of" $fcount
      source RUN
      n=$[$n+1]
      echo ""
      cd ..
    fi
  done
fi

if [ $FAILED == 0 ] ; then
  >&2 echo -e " All tests \033[32msuccessfully\033[39m completed."
  exit 0
else
  >&2 echo  -e " At least one test \033[31mfailed \033[39m please check the output carefully."
  exit 1
fi
