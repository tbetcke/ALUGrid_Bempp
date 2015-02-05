AC_DEFUN([DUNE_ALUGRID_CHECKS],[
  dnl check for DL malloc files
  AC_REQUIRE([ALUGRID_PATH_DLMALLOC])

  dnl check for sionlib 
  AC_REQUIRE([DUNE_PATH_SIONLIB])

  dnl check for ZOLTAN library 
  AC_REQUIRE([DUNE_PATH_ZOLTAN])

  dnl check for the METIS library (check needs to after ParMETIS check)
  AC_REQUIRE([IMMDX_LIB_METIS])

  dnl check for zlib
  AC_REQUIRE([DUNE_PATH_ZLIB])

  dnl check whether ALUGrid was found by the dune-grid module 
  dnl this conflicts with this package
  AS_IF([test "x$ALUGRID_CPPFLAGS" != "x"],
        [AC_MSG_ERROR([--with-alugrid conflicts with dune-alugrid module, remove the --with-alugrid from the configure options, use the --without-alugrid configure option, and rebuild dune-grid and dune-alugrid!])])

  dnl for some reason unknown to me the first entry appears at the bottom, the other
  dnl entries will apprear before the entries from dune-grid
  DUNE_DEFINE_ALUGRIDTYPE([ALUGRID_DUMMY_DONT_USE],[],[],[],[])

  DUNE_DEFINE_ALUGRIDTYPE([ALUGRID_CONFORM],[],[Dune::ALUGrid< dimgrid, dimworld, Dune::simplex, Dune::conforming >],[dune/alugrid/grid.hh],[dune/alugrid/dgf.hh])
  DUNE_DEFINE_ALUGRIDTYPE([ALUGRID_CUBE],[],[Dune::ALUGrid< dimgrid, dimworld, Dune::cube, Dune::nonconforming >],[dune/alugrid/grid.hh],[dune/alugrid/dgf.hh])
  DUNE_DEFINE_ALUGRIDTYPE([ALUGRID_SIMPLEX],[],[Dune::ALUGrid< dimgrid, dimworld, Dune::simplex, Dune::nonconforming >],[dune/alugrid/grid.hh],[dune/alugrid/dgf.hh])
])

AC_DEFUN([DUNE_ALUGRID_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-alugrid...])
  DUNE_CHECK_MODULES([dune-alugrid],[alugrid/grid.hh],[main])
])
