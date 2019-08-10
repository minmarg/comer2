# Create a library with the listed source files.
add_library( mycu STATIC
	cupro/IOProfileModelTxt.cu
	cupro/IOProfileModelBin.cu
	cupro/CuRoDb.cpp
	cupro/CuDbReader.cpp
	cupro/PMBatchProData.cpp
	cupro/CuBatchProcessingFinalizer.cpp
	cupro/CuBatchProcessing.cu
#
	cusco/CuBatchScoreMatrix_Init.cu
	cusco/CuBatchScoreMatrix_ssss.cu
	cusco/CuBatchScoreMatrix_cvs2s.cu
	cusco/CuBatchScoreMatrix_hdp1s.cu
	cusco/CuBatchScoreMatrix_Init_hdp1s.cu
	cusco/CuBatchScoreMatrix_ssss_cvs2s.cu
	cusco/CuBatchScoreMatrix_ssss_cvs2s_hdp1s.cu
	cusco/CuBatchScoreMatrix.cu
#
	cudp/CuBatchDP_init.cu
	cudp/CuBatchDP_init_block2.cu
	cudp/CuBatchDP_init_corr.cu
	cudp/CuBatchDP_init_btck.cu
	cudp/CuBatchDP_final.cu
	cudp/CuBatchDP.cu
#
	cuss/CuBatchSS_score_probs.cu
	cuss/CuBatchSS_dynplm.cu
	cuss/CuBatchSS.cu
#
	cumapdp/CuBatchMAPDP_fwd.cu
	cumapdp/CuBatchMAPDP_bwd.cu
	cumapdp/CuBatchMAPDP_map_btck.cu
	cumapdp/CuBatchMAPDP_map_final.cu
	cumapdp/CuBatchMAPDP_prnt.cu
	cumapdp/CuBatchMAPDP.cu
#
	cualn/Devices.cu
	cualn/DevCommThread.cpp
	cualn/AlnWriter.cpp
	cualn/JobDispatcher.cpp
)

# Make sure the compiler can find include files for the library
# when other libraries or executables link to it
target_include_directories( mycu PUBLIC 
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
)

set_target_properties( mycu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Link to the required libraries 
#target_link_libraries( mycu PUBLIC pro )

