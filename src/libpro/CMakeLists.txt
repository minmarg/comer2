# Create a library with the listed source files.
add_library( pro STATIC
	srcpro/MOptions.cpp
	srcpro/CLOptions.cpp
	srcpro/Configuration.cpp
	srcpro/SUBSTABLE.cpp
	srcpro/TRANSPROBS.cpp
	srcpro/TCTXVECT.cpp
	srcpro/BMSequence.cpp
	srcpro/BMProtoProfile.cpp
	srcpro/PMTransModel.cpp
	srcpro/PMProfileModelBase.cpp
	srcpro/PMProfileModel.cpp
	srcpro/SEGProfile.cpp
	srcpro/VirtScores.cpp
	srcpro/SSSScores.cpp
	srcpro/CVS2Scores.cpp
#
	srcpro/DbProProbs.cpp
	srcpro/Db.cpp
	srcpro/VirtualRoDb.cpp
#
	srcsco/AbstractScoreMatrix.cpp
	srcsco/ScoresAttr.cpp
	srcsco/IntegerScoreMatrix.cpp
	srcsco/ProfileMatrix.cpp
#
	srcaln/MSA.cpp
	srcaln/MSA_rw.cpp
)

# Make sure the compiler can find include files for the library
# when other libraries or executables link to it
target_include_directories( pro PUBLIC 
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link to the required libraries 
#target_link_libraries( pro PUBLIC lib )

