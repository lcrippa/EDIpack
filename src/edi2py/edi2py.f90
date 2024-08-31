module edi2py_bindings
    use edipack
    use scifor
    use iso_c_binding
    implicit none
    
    
contains

    include "edi2py_read_input.f90"
    include "edi2py_main.f90"
    include "edi2py_bath.f90"
    include "edi2py_io.f90"
    include "edi2py_bath_fit.f90"
    include "edi2py_aux_funx.f90"

end module edi2py_bindings
