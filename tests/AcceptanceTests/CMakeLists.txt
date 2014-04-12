file(GLOB ACCEPTANCE_SOURCES *.cpp *.h)

include_directories(../../src ../lib/googletest)
add_executable(AcceptanceTests ${ACCEPTANCE_SOURCES})
target_link_libraries(AcceptanceTests neb googletest pthread)

