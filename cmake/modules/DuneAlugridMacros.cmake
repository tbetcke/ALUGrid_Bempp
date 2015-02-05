

include(ALUGridType)

#make libdunealugrid known locally
set(DUNE_ALUGRID_LIBRARY "${PROJECT_BINARY_DIR}/lib/libdunealugrid.a"
  CACHE FILEPATH "path to local libs in dune-alugrid" )
mark_as_advanced(DUNE_ALUGRID_LIBRARY)



#define available alugrid types
dune_define_alugridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE ALUGRID_CONFORM
    DUNETYPE "Dune::ALUGrid< dimgrid, dimworld, simplex, conforming >"
    HEADERS dune/alugrid/grid.hh dune/alugrid/dgf.hh)
dune_define_alugridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE ALUGRID_CUBE
    DUNETYPE "Dune::ALUGrid< dimgrid, dimworld, cube, nonconforming >"
    HEADERS dune/alugrid/grid.hh dune/alugrid/dgf.hh)
dune_define_alugridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE ALUGRID_SIMPLEX
    DUNETYPE "Dune::ALUGrid< dimgrid, dimworld, simplex, nonconforming >"
    HEADERS dune/alugrid/grid.hh dune/alugrid/dgf.hh)

# avoid conflicts with normal ALUGrid
if( ALUGRID_CPPFLAGS )
  message(ERROR "--with-alugrid conflicts with dune-alugrid module, 
  remove the --with-alugrid from the configure options, 
  use the --without-alugrid configure option, 
  and rebuild dune-grid and dune-alugrid!") 
  #else()
  #set(HAVE_DUNE_ALUGRID 1)
endif()

set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-DENABLE_ALUGRID")
foreach(dir ${ALUGRID_INCLUDES})
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
endforeach()

# contained in cmake system modules
find_package(ZLIB)
#set HAVE_ZLIB for config.h
set(HAVE_ZLIB ${ZLIB_FOUND})

find_package(SIONlib)
include(AddSIONlibFlags)
find_package(DLMalloc)
find_package(ZOLTAN)
include(AddZOLTANFlags)

message(AUTHOR_WARNING "TODO: Improve module test.")
