#list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/<nome file>")
#list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/<nome file>")


set(DFN_sources ${src_sources} PARENT_SCOPE)
set(DFN_headers ${src_headers} PARENT_SCOPE)
set(DFN_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
