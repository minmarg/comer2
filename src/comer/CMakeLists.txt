# Add executable that is built from the listed source files.
# The extensions are automatically found.
add_executable( comer 
	comer.cpp 
)


file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LOCALSTATEDIR}" COMERSTATEDIRNAT)
string(REGEX REPLACE "([\\])" "\\\\\\1" COMERSTATEDIRNAT "${COMERSTATEDIRNAT}")


## ONE executable
## Specify include directories
target_include_directories( comer PUBLIC 
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions( comer PUBLIC
	COMER_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
	COMERSTATEDIR=${COMERSTATEDIRNAT}
)

set_target_properties( comer PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties( mycu PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# Link the executable to the libraries 
target_link_libraries( comer PUBLIC mycu pro HDP seg lib extsp Threads::Threads)


## install targets
install(TARGETS 
	comer 
	DESTINATION ${LOCALBINDIR}
)
#install(FILES 
#	comer.sh 
#	DESTINATION ${LOCALBINDIR} 
#	PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_READ WORLD_READ
#)

