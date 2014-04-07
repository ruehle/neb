module potential
contains
function get_energy_gradient(ncoords, coords, grad, userdata) result(energy) bind(c)
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: ncoords
    real(c_double), intent(in) :: coords(ncoords)
    real(c_double), intent(out) :: grad(ncoords) 
    type(c_ptr) :: userdata
    real(c_double) :: energy

    grad = 0.0
    energy = 0.0
end function
end module

program neb_example
    use neb
    use potential
    use iso_c_binding
    implicit none
    logical foo
    call neb_setup(c_funloc(get_energy_gradient), c_null_ptr)
    call neb_initialize_path(2, 2)
    call neb_start()   
    do while ( neb_step() )
        print *,"print some statistics"
    end do
    call neb_cleanup()
end program
