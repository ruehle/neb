!        typedef double potential_callback_t(int ncoords, double *coords, double
!*grad, void
!                *userdata);
!
!
!        void neb_set_image_coords(int image, int ncoords, const double *coords);
!        void neb_get_image_coords(int image, int ncoords, double *coords);
!

module neb
   use iso_c_binding
   interface
       ! void neb_setup(potential_callback_t *potential, void *userdata);
       subroutine neb_setup(potential_callback, userdata) bind(c)
           use iso_c_binding
           implicit none
           type(c_funptr), value, intent(in) :: potential_callback 
           type(c_ptr), intent(in) :: userdata
       end subroutine
 
       ! void neb_cleanup();
       subroutine neb_cleanup() bind(C)
       end subroutine 
    
       ! void neb_initialize_path(int nimages, int num_coords_per_image);
       subroutine neb_initialize_path(nimages, num_coords_per_image) bind(c)
           use iso_c_binding
           integer(c_int), value, intent(in) :: nimages
           integer(c_int), value, intent(in) :: num_coords_per_image
       end subroutine

       ! void neb_start();
       subroutine neb_start() bind(c)
       end subroutine
       
       ! bool neb_step();
       function neb_step() bind(c)
           use iso_c_binding
           logical(c_bool) :: neb_step
       end function
    end interface
end module neb

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

    print *,"potential call",coords
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
    foo = neb_step()
    call neb_cleanup()
end program
