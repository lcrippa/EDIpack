module edi2py_bindings
    use edipack
    use scifor
    use iso_c_binding
    implicit none
    
    
contains

    function i2l(var_integer) result (var_logical)
      integer    :: var_integer
      logical    :: var_logical   
    
      if (var_integer == 1) then
        var_logical = .true.
      else
        var_logical = .false.
      endif
    end function i2l
    
    function l2i(var_logical) result (var_integer)
      integer    :: var_integer
      logical    :: var_logical   
    
      if (var_logical) then
        var_integer = 1
      else
        var_integer = 0
      endif
    end function l2i

    include "edi2py_read_input.f90"
    include "edi2py_main.f90"
    include "edi2py_bath.f90"
    include "edi2py_io.f90"
    include "edi2py_bath_fit.f90"
    include "edi2py_aux_funx.f90"

end module edi2py_bindings
