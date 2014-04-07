enable_language(Fortran)

add_executable(neb_fortran neb_fortran.f90 neb.f90)
target_link_libraries(neb_fortran neb)
