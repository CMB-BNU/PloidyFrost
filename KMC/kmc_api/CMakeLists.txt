file(GLOB sources *.cpp )
file(GLOB headers *.h *.hpp *.hh )


add_library(kmc_static STATIC ${sources} ${headers})
add_library(kmc_dynamic SHARED ${sources} ${headers})

set_target_properties(kmc_static PROPERTIES OUTPUT_NAME "kmc")
set_target_properties(kmc_dynamic PROPERTIES OUTPUT_NAME "kmc")

target_include_directories(kmc_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(kmc_dynamic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


install(TARGETS kmc_dynamic DESTINATION lib)
install(TARGETS kmc_static DESTINATION lib)
install(FILES ${headers} DESTINATION include/kmc)


