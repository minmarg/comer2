# Create a library with the listed source files.
# Any number of sources could be listed here.
add_library( lib STATIC 
	mydirent.cpp
	msg.cpp
	myexception.cpp
	root.cpp
	alpha.cpp
	mysort.cpp
	myfiler.cpp
	mygetopt.cpp
	mystring.cpp
	logitnormal.cpp
	ConfigFile.cpp
	CtxtCoefficients.cpp
	BinarySearchStructure.cpp
)

# Make sure the compiler can find include files for the library
# when other libraries or executables link to it
target_include_directories( lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
