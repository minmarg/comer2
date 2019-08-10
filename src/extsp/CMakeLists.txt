# Create a library with the listed source files.
add_library( extsp STATIC
	pslerror.cpp
	cheb.cpp
	exp.cpp
	zeta.cpp
	digamma.cpp
	gamma.cpp
	rng.cpp
	rv/rvnorm.cpp
	pslvector.cpp
	pslmatrix.cpp
	nsmatrix.cpp
	spdmatrix.cpp
	ivector.cpp
)

# Make sure the compiler can find include files for the library
# when other libraries or executables link to it
target_include_directories( extsp PUBLIC 
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link to the required libraries 
#target_link_libraries( extsp PUBLIC ... )

