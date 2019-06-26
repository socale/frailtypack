dnl Process this file with autoconf to produce a configure script.
dnl This file is derived from Rmpi package

AC_INIT(DESCRIPTION)

## One way is to specify both include and lib paths
AC_ARG_WITH([lcmm-include],
            AC_HELP_STRING([--with-lcmm-include=INCLUDE_PATH],
                           [location of MPI header files]),
            [mpi_include_path=$withval])
if test -n "${mpi_include_path}" ; then
   MPI_INCLUDE_PATH="${mpi_include_path}"
   if test ! -f "${MPI_INCLUDE_PATH}/mpi.h"; then
     AC_MSG_ERROR([Value of --with-lcmm-include does not contain mpi.h])
   fi
   if test ! -f "${MPI_INCLUDE_PATH}/mpi.mod"; then
     AC_MSG_ERROR([Value of --with-lcmm-include does not contain mpi.mod])
   fi
else
  if test -n "${LCMM_INCLUDE}" ; then
    MPI_INCLUDE_PATH="${LCMM_INCLUDE}"
    if test ! -f "${MPI_INCLUDE_PATH}/mpi.h"; then
      AC_MSG_ERROR([Value of LCMM_INCLUDE does not contain mpi.h])
    fi
    if test ! -f "${MPI_INCLUDE_PATH}/mpi.mod"; then
      AC_MSG_ERROR([Value of LCMM_INCLUDE does not contain mpi.mod])
    fi
    echo "setting mpi include path from LCMM_INCLUDE"
  fi
fi

AC_ARG_WITH([lcmm-libpath],
            AC_HELP_STRING([--with-lcmm-libpath=LIB_PATH],
                           [location of MPI library files]),
            [mpi_lib_path=$withval])
if test -n "${mpi_lib_path}" ; then
   MPI_LIBPATH="${mpi_lib_path}"
else
  if test -n "${LCMM_LIB_PATH}" ; then
    MPI_LIBPATH="${LCMM_LIB_PATH}"
    echo "setting mpi lib path from LCMM_LIB_PATH"
  fi
fi

if test -n "${MPI_INCLUDE_PATH}" -a -z "${MPI_LIBPATH}" ; then
  AC_MSG_ERROR([Must specify both the include and lib paths])
fi
if test -z "${MPI_INCLUDE_PATH}" -a -n "${MPI_LIBPATH}" ; then
  AC_MSG_ERROR([Must specify both the include and lib paths])
fi

MPITYPE="UNKNOWN"

AC_ARG_WITH([lcmm-type],
            AC_HELP_STRING([--with-lcmm-type=MPI_TYPE],
                           [the type of MPI: OPENMPI,LAM,MPICH,MPICH2, or CRAY]),
            [mpi_type=$withval])
if test -n "${mpi_type}" ; then
   MPITYPE="${mpi_type}"
else
  if test -n "${LCMM_TYPE}" ; then
     MPITYPE="${LCMM_TYPE}"
  fi
fi

if test -z "${MPI_INCLUDE_PATH}" ; then
dnl 'paths not specified'

## another way is to specify MPI_ROOT, perhaps via --with-mpi
AC_ARG_WITH(mpi,
            AC_HELP_STRING([--with-mpi=LIB_PATH],
                           [location of top-level MPI directory]),
[   if test "${withval}" != "no"; then
        MPI_ROOT=${withval}
        MPI_INCLUDE_PATH="${MPI_ROOT}/include"
        LDFLAGS="$LDFLAGS -L${withval}/lib"
        CPPFLAGS="$CPPFLAGS -I${withval}/include"
    fi  ]
)

CC=`"${R_HOME}/bin/R" CMD config CC`
FC=`"${R_HOME}/bin/R" CMD config FC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
FCFLAGS=`"${R_HOME}/bin/R" CMD config FCFLAGS`
AC_PROG_CC
if test "X$CC" = "X" ; then
	# Our best guess...
	CC="gcc -std=gnu99"
fi

if test -z "$MPI_ROOT" ; then
	#Use pkg-config to find openmpi
	AC_PATH_PROG([MPIF], [mpif90])

	if test -n "${MPIF}"; then
		AC_MSG_CHECKING([if mpif90 has showme option])
		if "${MPIF}" --showme:version; then
			AC_MSG_RESULT([yes])
			MPITYPE="OPENMPI"
			MPI_ROOT="NONEED"
		
			OMPI_INCLUDE_PATH=`${MPIF} --showme:compile`
			OMPI_LIBS=`"${MPIF}" --showme:link`
		else
			AC_MSG_RESULT([no])
		fi
	fi
fi

if test -z "$MPI_ROOT" ; then
  for d in /opt/lib /usr/lib /usr /usr/lib/x86_64-linux-gnu/openmpi \
          /usr/local/lib /usr/local /usr/lib64/mpi/gcc/openmpi; do
    if test -f $d/include/mpi.h && test -f $d/include/mpi.mod && test -d $d/lib/openmpi; then
      echo "I am here $d and it is OpenMPI"
      MPI_ROOT=$d
      MPITYPE="OPENMPI"
      break
    elif test -f $d/include/mpi.h && test -f $d/lib/mpi.mod && test -d $d/lib/openmpi; then
      echo "I am here $d and it is OpenMPI"
      MPI_ROOT=$d
      MPITYPE="OPENMPI"
      MPI_INCLUDE_PATH=$d/lib
      break
    elif test -f $d/include/openmpi/mpi.h && test -f $d/include/openmpi/mpi.mod && test -d $d/lib/openmpi; then
      echo "I am here $d and it is OpenMPI"
      MPI_ROOT=$d
      MPITYPE="OPENMPI"
      MPI_INCLUDE_PATH=$d/include/openmpi
      break
    elif test -f $d/include/openmpi/mpi.h && test -f $d/include/openmpi/mpi.mod && test -d $d/lib64/openmpi; then
      echo "I am here $d and it is OpenMPI"
      MPI_ROOT=$d
      MPITYPE="OPENMPI"
      MPI_INCLUDE_PATH=$d/include/openmpi
      break
    elif test -f $d/openmpi/include/mpi.h && test -f $d/include/openmpi/include/mpi.mod && test -d $d/openmpi/lib; then
      echo "I am here $d/openmpi and it is OpenMPI"
      MPI_ROOT=$d/openmpi
      MPITYPE="OPENMPI"
      MPI_INCLUDE_PATH=$d/include/openmpi
      break
    elif test -f $d/include/mpi.h && test -f $d/include/mpi.mod && test -d $d/lib64; then
      echo "I am here $d and it is OpenMPI"
      MPI_ROOT=$d
      MPITYPE="OPENMPI"
      MPI_INCLUDE_PATH=$d/include
      break
   elif test -f $d/lam/include/mpi.h && test -f $d/lam/include/mpi.mod && test -f $d/lam/lib/libmpi.so; then
      echo "I am here $d/lam and it is LAM"
      MPI_ROOT=$d/lam
      MPITYPE="LAM"
      MPI_INCLUDE_PATH=$d/include/lam
      break
   elif test -f $d/include/lam/mpi.h && test -f $d/include/lam/mpi.mod && test -f $d/lib/lam/libmpi.so; then
      echo "I am here $d and it is LAM"
      MPI_ROOT=$d
      MPITYPE="LAM"
      MPI_INCLUDE_PATH=$d/include/lam
      break
    elif test -f $d/include/lam/mpi.h && test -f $d/include/lam/mpi.mod && test -f $d/lib64/lam/libmpi.so; then
      echo "I am here $d and it is LAM"
      MPI_ROOT=$d
      MPITYPE="LAM"
      MPI_INCLUDE_PATH=$d/include/lam
      break
    elif test -f $d/mpich/include/mpi.h && test -f $d/mpich/include/mpi.mod; then
      echo "I am here $d/mpich and it is MPICH"
      MPI_ROOT=$d/mpich
      MPITYPE="MPICH"
      MPI_INCLUDE_PATH=$d/mpich/include
      break
    elif test -f $d/mpi/include/mpi.h && test -f $d/mpi/include/mpi.mod; then
      echo "I am here $d/mpi"
      MPI_ROOT=$d/mpi
      MPI_INCLUDE_PATH=$d/mpi/include
      break
    elif test -f $d/include/mpi.h && test -f $d/include/mpi.mod; then
      echo "I am here $d"
      MPI_ROOT=$d
      break
    fi
  done
fi

if test "$MPI_ROOT" != "NONEED"; then
	if test -n "$MPI_ROOT"; then
		#guess MPITYPE if it is unknown and user supplies MPI_ROOT
		if test $MPITYPE = "UNKNOWN" ; then
			if test -d ${MPI_ROOT}/openmpi; then
				MPITYPE="OPENMPI"
			elif test -d ${MPI_ROOT}/lib/openmpi; then
				MPITYPE="OPENMPI"
			elif test -d ${MPI_ROOT}/lib64/openmpi; then
				MPITYPE="OPENMPI"
			elif test -f ${MPI_ROOT}/lib/liblam.so; then
				MPITYPE="LAM"
			elif test -d ${MPI_ROOT}/lib/lam; then
				MPITYPE="LAM"
			elif test -d ${MPI_ROOT}/lib64/lam; then
				MPITYPE="LAM"
			elif test -f ${MPI_ROOT}/lib/libmpich.a; then
				MPITYPE="MPICH"
			elif test -f ${MPI_ROOT}/lib64/libmpich.a; then
				MPITYPE="MPICH"
			fi
		fi
		echo "Trying to find mpi.h ..."
		if test -f ${MPI_ROOT}/include/mpi.h; then
			echo "Found in ${MPI_ROOT}/include"
			MPI_INCLUDE_PATH="${MPI_ROOT}/include"		  
		elif test -f ${MPI_ROOT}/include/openmpi/mpi.h; then
			echo "Found in ${MPI_ROOT}/include/openmpi"
			MPI_INCLUDE_PATH="${MPI_ROOT}/include/openmpi"
		elif test -f ${MPI_ROOT}/include/lam/mpi.h; then
			echo "Found in ${MPI_ROOT}/include/lam"
			MPI_INCLUDE_PATH="${MPI_ROOT}/include/lam"
		else 
			echo "Cannot find mpi.h header file"
			echo "Please use --with-mpi=/path/to/mpi or specify the environment variable MPI_ROOT"
			exit 1
		fi
		echo "Trying to find mpi.mod ..."
		if test -f ${MPI_ROOT}/include/mpi.mod; then
			echo "Found in ${MPI_ROOT}/include"
			MPI_INCLUDE_PATH="${MPI_ROOT}/include"		  
		elif test -f ${MPI_ROOT}/lib/mpi.mod; then
			echo "Found in ${MPI_ROOT}/lib"
			MPI_INCLUDE_PATH="${MPI_ROOT}/lib"
		elif test -f ${MPI_ROOT}/include/openmpi/mpi.mod; then
			echo "Found in ${MPI_ROOT}/include/openmpi"
			MPI_INCLUDE_PATH="${MPI_ROOT}/include/openmpi"
		elif test -f ${MPI_ROOT}/include/lam/mpi.mod; then
			echo "Found in ${MPI_ROOT}/include/lam"
			MPI_INCLUDE_PATH="${MPI_ROOT}/include/lam"
		else 
			echo "Cannot find mpi.mod header file"
			echo "Please use --with-mpi=/path/to/mpi or specify the environment variable MPI_ROOT"
			exit 1
		fi
	else
		AC_CHECK_HEADER(mpi.h,, AC_MSG_ERROR(["Cannot find mpi.h header file"]))
	fi 

	echo "Trying to find libmpi.so or libmpich.a ..."
	if test -f ${MPI_ROOT}/lib/libmpi.so; then
        echo "Found libmpi in ${MPI_ROOT}/lib"
        MPI_LIBPATH="${MPI_ROOT}/lib"
	elif test -f ${MPI_ROOT}/libmpi.so; then
		dnl This is what Fedora 10 had
        echo "Found libmpi in ${MPI_ROOT}"
        MPI_LIBPATH="${MPI_ROOT}"
	elif test -f ${MPI_ROOT}/lib/openmpi/libmpi.so; then
        echo "Found libmpi in ${MPI_ROOT}/lib/openmpi"
        MPI_LIBPATH="${MPI_ROOT}/lib/openmpi"
	elif test -f ${MPI_ROOT}/lib64/openmpi/libmpi.so; then
        echo "Found libmpi in ${MPI_ROOT}/lib64/openmpi"
        MPI_LIBPATH="${MPI_ROOT}/lib64/openmpi"
	elif test -f ${MPI_ROOT}/lib64/libmpi.so; then
        echo "Found libmpi in ${MPI_ROOT}/lib64"
        MPI_LIBPATH="${MPI_ROOT}/lib64"		
	elif test -f ${MPI_ROOT}/lib/lam/libmpi.so; then
        echo "Found libmpi in ${MPI_ROOT}/lib/lam"
        MPI_LIBPATH="${MPI_ROOT}/lib/lam"
		MPI_INCL2="-I${MPI_INCLUDE_PATH}/32"
	elif test -f ${MPI_ROOT}/lib64/lam/libmpi.so; then
        echo "Found libmpi in ${MPI_ROOT}/lib64/lam"
        MPI_LIBPATH="${MPI_ROOT}/lib64/lam"
		MPI_INCL2="-I${MPI_INCLUDE_PATH}/64"
	elif test -f ${MPI_ROOT}/lib/libmpich.so; then
        echo "Found libmpich in ${MPI_ROOT}/lib"
        MPI_LIBPATH="${MPI_ROOT}/lib"
		AC_CHECK_DECLS([MPICH2, MPICH2],[  MPI_DEFS="-DMPI2" ],, 
		[#include <mpi.h>])
	elif test -f ${MPI_ROOT}/lib64/libmpich.a; then
        echo "Found libmpich in ${MPI_ROOT}/lib64"
        MPI_LIBPATH="${MPI_ROOT}/lib64"
		AC_CHECK_DECLS([MPICH2, MPICH2],[  MPI_DEFS="-DMPI2" ],, 
		[#include <mpi.h>])
	else
		AC_CHECK_LIB(mpi, main, MPI_LIBPATH="${MPI_ROOT}/lib",
    	[   echo "libmpi not found. exiting..."
        	exit 1  ]
		)
	fi
fi

if test "$MPITYPE" = "OPENMPI" ; then
   AC_CHECK_PROG(ORTED, orted, yes, no)
   if test "$ORTED" = no ; then
       AC_MSG_ERROR(Cannot find orted. lcmm needs orted to run.)
   fi
fi

if test "$MPITYPE" = "LAM" ; then
	echo "Try to find liblam.so ..."
	if test -f ${MPI_ROOT}/lib/liblam.so; then
        echo "Found liblam in ${MPI_ROOT}/lib"
        MPI_LIBPATH="${MPI_ROOT}/lib"
	elif test -f ${MPI_ROOT}/lib/lam/liblam.so; then
        echo "Found liblam in ${MPI_ROOT}/lib/lam"
        MPI_LIBPATH="${MPI_ROOT}/lib/lam"
	elif test -f ${MPI_ROOT}/lib64/lam/liblam.so; then
        echo "Found liblam in ${MPI_ROOT}/lib64/lam"
        MPI_LIBPATH="${MPI_ROOT}/lib64/lam"
	else
    	AC_CHECK_LIB(lam, main, []
    	[   echo "liblam not found. Probably not LAM-MPI"]
    	)
	fi    
fi

fi dnl end of 'paths not specified'

## now we have found the include and lib paths and may know the type

AC_CHECK_LIB(util, openpty, [ MPI_LIBS="$MPI_LIBS -lutil" ])
AC_CHECK_LIB(pthread, main, [ MPI_LIBS="$MPI_LIBS -lpthread" ])

if test "$MPITYPE" = "OPENMPI" -o "$MPITYPE" = "LAM" -o "$MPITYPE" = "MPICH2" -o "$MPITYPE" = "INTELMPI" -o "$MPITYPE" = "CRAY"; then
  MPI_DEFS="-DMPI2"
fi

if test "$MPI_ROOT" = "NONEED"; then
	PKG_FCFLAGS="${OMPI_INCLUDE_PATH} -I${MPI_INCLUDE_PATH} ${MPI_DEFS} -D${MPITYPE} -Wl,-rpath ${MPI_LIBPATH}"
else
	PKG_FCFLAGS="-I${MPI_INCLUDE_PATH} ${MPI_INCL2} ${MPI_DEFS} -D${MPITYPE} -Wl,-rpath ${MPI_LIBPATH}"
fi

case "$MPITYPE" in
  OPENMPI)
    if test "$MPI_ROOT" = "NONEED"; then
		PKG_LIBS="${OMPI_LIBS} ${MPI_LIBS}"
	else
		PKG_LIBS="-L${MPI_LIBPATH} -lmpi_mpifh -lmpi ${MPI_LIBS}"
	fi
    ;;
  LAM)
    PKG_LIBS="-L${MPI_LIBPATH} -lmpi -llam ${MPI_LIBS}"
    ;;
  MPICH)
    PKG_LIBS="-L${MPI_LIBPATH} -lmpich -lmpl -lopa ${MPI_LIBS}"
    ;;
  MPICH2)
    PKG_LIBS="-L${MPI_LIBPATH} -lmpich -lmpl -lopa ${MPI_LIBS}"
    ;;
  INTELMPI)
    PKG_LIBS="-Wl,-rpath,${MPI_LIBPATH} -L${MPI_LIBPATH} -lmpich -lmpigf -lmpigi ${MPI_LIBS}" 
    ;;
  CRAY)
    PKG_LIBS="-L${MPI_LIBPATH} -lmpich_cray_default64 -lmpl -lopa ${MPI_LIBS}"
    ;;
  *)
    AC_MSG_ERROR([Unknown type of mpi: use --with-lcmm-type or LCMM_TYPE to specify it])
    ;;
esac

AC_SUBST(PKG_FCFLAGS)
AC_SUBST(PKG_LIBS)

AC_OUTPUT(src/Makevars)
