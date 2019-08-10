# Create a library with the listed source files.
add_library( HDP STATIC
	basin.cpp
	dish.cpp
	menu.cpp
	table.cpp
	restaurant.cpp
	rntchain.cpp
	HDPscores.cpp
	HDPbase.cpp
	HDPbase_probvec.cpp
	HDPbase_probmtx.cpp
	HDPbase_rw.cpp
)

# Make sure the compiler can find include files for the library
# when other libraries or executables link to it
target_include_directories( HDP PUBLIC 
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link to the required libraries 
#target_link_libraries( HDP PUBLIC ... )

