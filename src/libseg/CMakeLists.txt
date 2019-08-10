# Create a library with the listed source files.
add_library( seg STATIC
	segdata.cpp
	Segments.cpp
	SEGAbstract.cpp
	SEGSequence.cpp
)

# Make sure the compiler can find include files for the library
# when other libraries or executables link to it
target_include_directories( seg PUBLIC 
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link to the required libraries 
#target_link_libraries( seg PUBLIC ... )

