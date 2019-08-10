# Add executable that is built from the listed source files.
# The extensions are automatically found.
add_executable( makedb 
	makedb.cpp 
)


file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LOCALSTATEDIR}" COMERSTATEDIRNAT)
string(REGEX REPLACE "([\\])" "\\\\\\1" COMERSTATEDIRNAT "${COMERSTATEDIRNAT}")


## Specify include directories
target_include_directories( makedb PUBLIC 
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions( makedb PUBLIC
	COMER_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
	COMERSTATEDIR=${COMERSTATEDIRNAT}
)

# Link the executable to the libraries 
target_link_libraries( makedb PUBLIC pro HDP seg lib extsp )


## install targets
install(TARGETS 
	makedb 
	DESTINATION ${LOCALBINDIR}
)

