set(alucommondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/alugrid/common)
file(GLOB alucommon_HEADERS *.hh)
# installation of the header files                  
install(FILES ${alucommon_HEADERS} DESTINATION ${alucommondir})
