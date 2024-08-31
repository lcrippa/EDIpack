from ctypes import *
import numpy as np
import os

#################################
#link to access global variables
#################################

#dummy class, to be filled
class Link:
    pass

#function that will add a variable to the dummy class, will be called in variable definition
def add_global_variable(obj, dynamic_name, target_object, target_attribute):
    @property
    def getter(self):
        try:
            attrib = getattr(target_object, target_attribute)
            try: #this is for strings
                attrib = attrib.decode()
            except:
                pass
        except: #this is for arrays
            if(len(target_object)>1):
                return [target_object[x] for x in range(len(target_object))]
        return attrib

    @getter.setter
    def setter(self, new_value):
        try: #this is for arrays
            if(len(target_object)>1):
                minlength=min(len(target_object),len(new_value))
                target_object[0:minlength]=new_value[0:minlength]
        except:
            try:
                new_value = new_value.encode()
            except:
                pass
            setattr(target_object, target_attribute, new_value)

    # Dynamically add the property to the class
    setattr(obj.__class__, dynamic_name, getter)
    setattr(obj.__class__, dynamic_name, setter)

#create the class containing all the global variables
vars=Link()

######################################
# Load shared library with C-bindings
######################################

libpath = os.path.dirname(os.path.realpath(__file__))
libfile = os.path.join(libpath, 'libedi2py.so')
libedi2py = CDLL(libfile)

######################################
# GLOBAL VARIABLES
######################################


add_global_variable(vars, "Nbath", c_int.in_dll(libedi2py, "Nbath"), "value")
add_global_variable(vars, "Norb", c_int.in_dll(libedi2py, "Norb"), "value")
add_global_variable(vars, "Nspin", c_int.in_dll(libedi2py, "Nspin"), "value")
add_global_variable(vars, "Nloop", c_int.in_dll(libedi2py, "Nloop"), "value")
add_global_variable(vars, "Nph", c_int.in_dll(libedi2py, "Nph"), "value")
add_global_variable(vars, "Uloc", ARRAY(c_double, 5).in_dll(libedi2py, "Uloc"), "value")
add_global_variable(vars, "Ust", c_double.in_dll(libedi2py, "Ust"), "value")
add_global_variable(vars, "beta", c_double.in_dll(libedi2py, "beta"), "value")
add_global_variable(vars, "Lmats", c_int.in_dll(libedi2py, "Lmats"), "value")
add_global_variable(vars, "Lreal", c_int.in_dll(libedi2py, "Lreal"), "value")
add_global_variable(vars, "wini", c_double.in_dll(libedi2py, "wini"), "value")
add_global_variable(vars, "wfin", c_double.in_dll(libedi2py, "wfin"), "value")
add_global_variable(vars, "eps", c_double.in_dll(libedi2py, "eps"), "value")
add_global_variable(vars, "dmft_error", c_double.in_dll(libedi2py, "dmft_error"), "value")


######################################
# READ_INPUT
######################################

read_input = libedi2py.read_input
read_input.argtypes = [c_char_p]  # allows for automatic type conversion
read_input.restype = None

######################################
# ED_MAIN FUNCTIONS
######################################

# Define the function signature for the Fortran function `init_solver_site`.
init_solver_site = libedi2py.init_solver_site
init_solver_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int]  # allows for automatic type conversion
init_solver_site.restype = None

# Define the function signature for the Fortran function `init_solver_ineq`.
init_solver_ineq = libedi2py.init_solver_ineq
init_solver_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int]  # allows for automatic type conversion
init_solver_ineq.restype = None

# Define the function signature for the Fortran function `solve_site`.
solve_site = libedi2py.solve_site
solve_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,np.ctypeslib.ndpointer(ndim=4, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
solve_site.restype = None

# Define the function signature for the Fortran function `solve_ineq`.
solve_ineq = libedi2py.solve_ineq
solve_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,np.ctypeslib.ndpointer(ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
solve_ineq.restype = None

######################################
# ED_BATH FUNCTIONS
######################################

get_bath_dimension = libedi2py.get_bath_dimension
get_bath_dimension.argtypes = None  # allows for automatic type conversion
get_bath_dimension.restype = c_int

######################################
# ED_IO FUNCTIONS
######################################

get_sigma_matsubara_site = libedi2py.get_sigma_matsubara_site
get_sigma_matsubara_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_matsubara_site.restype = None

get_sigma_matsubara_ineq = libedi2py.get_sigma_matsubara_ineq
get_sigma_matsubara_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_matsubara_ineq.restype = None

get_sigma_realaxis_site = libedi2py.get_sigma_realaxis_site
get_sigma_realaxis_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_realaxis_site.restype = None

get_sigma_realaxis_ineq = libedi2py.get_sigma_matsubara_ineq
get_sigma_realaxis_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_realaxis_ineq.restype = None

######################################
# ED_BATH_FIT FUNCTIONS
######################################

chi2_fitgf_site = libedi2py.chi2_fitgf_site
chi2_fitgf_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int] 
chi2_fitgf_site.restype = None

chi2_fitgf_ineq = libedi2py.chi2_fitgf_ineq
chi2_fitgf_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int,
                            np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,
                            np.ctypeslib.ndpointer(dtype=float,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,
                            c_int] 
chi2_fitgf_ineq.restype = None

######################################
# ED_AUX_FUNX FUNCTIONS
######################################

check_convergence = libedi2py.check_convergence
check_convergence.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),c_int,c_double,c_int,c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=int,ndim=1, flags='F_CONTIGUOUS')]
check_convergence.restype = None