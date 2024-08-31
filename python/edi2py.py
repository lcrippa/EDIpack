from ctypes import *
import numpy as np
import os
import types

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


######################################
# Load shared library with C-bindings
######################################

libpath = os.path.dirname(os.path.realpath(__file__))
libfile = os.path.join(libpath, 'libedi2py.so')
libedi2py = CDLL(libfile)

######################################
# READ_INPUT
######################################

read_input_wrap = libedi2py.read_input
read_input_wrap.argtypes = [c_char_p]  # allows for automatic type conversion
read_input_wrap.restype = None

def read_input(self,input_string):
    c_string = c_char_p(input_string.encode())
    read_input_wrap(c_string)
    
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

def init_solver(self,bath):
    dim_bath=np.shape(bath)
    if len(dim_bath)<2:
        init_solver_site(bath,dim_bath[0])
    else:
        init_solver_ineq(bath,dim_bath[0],dim_bath[1])

    
def solve(self,bath,hloc,sflag=True,mpi_lanc=False):
    dim_bath=np.shape(bath)
    dim_hloc=np.shape(hloc)
    if len(dim_bath)<2:
        solve_site(bath,dim_bath[0],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],sflag)
    else:
        solve_ineq(bath,dim_bath[0],dim_bath[1],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],dim_hloc[4],mpi_lanc)
        

######################################
# ED_BATH FUNCTIONS
######################################

get_bath_dimension_wrap = libedi2py.get_bath_dimension
get_bath_dimension_wrap.argtypes = None  # allows for automatic type conversion
get_bath_dimension_wrap.restype = c_int

def get_bath_dimension(self):
    return get_bath_dimension_wrap()
    
######################################
# ED_IO FUNCTIONS
######################################

#sigma_matsubara
get_sigma_matsubara_site = libedi2py.get_sigma_matsubara_site
get_sigma_matsubara_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_matsubara_site.restype = None

get_sigma_matsubara_ineq = libedi2py.get_sigma_matsubara_ineq
get_sigma_matsubara_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_matsubara_ineq.restype = None

def get_sigma_matsubara(self,Smats):
    dim_Smats=np.shape(Smats)
    if len(dim_Smats)==5:
        get_sigma_matsubara_site(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4])
    else:
        get_sigma_matsubara_site(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4],dim_Smats[5])

#sigma_realaxis
get_sigma_realaxis_site = libedi2py.get_sigma_realaxis_site
get_sigma_realaxis_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_realaxis_site.restype = None

get_sigma_realaxis_ineq = libedi2py.get_sigma_matsubara_ineq
get_sigma_realaxis_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_realaxis_ineq.restype = None

        
def get_sigma_realaxis(self,Sreal):
    dim_Sreal=np.shape(Sreal)
    if len(dim_Sreal)==5:
        get_sigma_realaxis_site(Sreal,dim_Sreal[0],dim_Sreal[1],dim_Sreal[2],dim_Sreal[3],dim_Sreal[4])
    else:
        get_sigma_realaxis_site(Sreal,dim_Sreal[0],dim_Sreal[1],dim_Sreal[2],dim_Sreal[3],dim_Sreal[4],dim_Sreal[5])
        
#gimp_matsubara
get_gimp_matsubara_site = libedi2py.get_gimp_matsubara_site
get_gimp_matsubara_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_matsubara_site.restype = None

get_gimp_matsubara_ineq = libedi2py.get_gimp_matsubara_ineq
get_gimp_matsubara_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_matsubara_ineq.restype = None

def get_gimp_matsubara(self,Gmats):
    dim_Gmats=np.shape(Gmats)
    if len(dim_Gmats)==5:
        get_gimp_matsubara_site(Gmats,dim_Gmats[0],dim_Gmats[1],dim_Gmats[2],dim_Gmats[3],dim_Gmats[4])
    else:
        get_gimp_matsubara_site(Gmats,dim_Gmats[0],dim_Gmats[1],dim_Gmats[2],dim_Gmats[3],dim_Gmats[4],dim_Gmats[5])

#gimp_realaxis
get_gimp_realaxis_site = libedi2py.get_gimp_realaxis_site
get_gimp_realaxis_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_realaxis_site.restype = None

get_gimp_realaxis_ineq = libedi2py.get_gimp_matsubara_ineq
get_gimp_realaxis_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_realaxis_ineq.restype = None

        
def get_gimp_realaxis(self,Gmats):
    dim_Gmats=np.shape(Gmats)
    if len(dim_Gmats)==5:
        get_gimp_realaxis_site(Gmats,dim_Gmats[0],dim_Gmats[1],dim_Gmats[2],dim_Gmats[3],dim_Gmats[4])
    else:
        get_gimp_realaxis_site(Gmats,dim_Gmats[0],dim_Gmats[1],dim_Gmats[2],dim_Gmats[3],dim_Gmats[4],dim_Gmats[5])
        

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

def chi2_fitgf(self,func,bath,hloc=None,ispin=1,iorb=0):
    dim_func=np.shape(func)
    dim_bath=np.shape(bath)
    dim_hloc=np.shape(hloc)
    if(len(dim_func)==5):
        chi2_fitgf_site(func,dim_func[0],dim_func[1],dim_func[2],dim_func[3],dim_func[4],bath,dim_bath[0],ispin,iorb)
    else:
        chi2_fitgf_ineq(func,dim_func[0],dim_func[1],dim_func[2],dim_func[3],dim_func[4],dim_func[5],bath,dim_bath[0],dim_bath[1],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],dim_hloc[4],ispin)
        


######################################
# ED_AUX_FUNX FUNCTIONS
######################################

check_convergence_wrap = libedi2py.check_convergence
check_convergence_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),c_int,c_double,c_int,c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=int,ndim=1, flags='F_CONTIGUOUS')]
check_convergence_wrap.restype = None

def check_convergence(self,func,threshold,N1,N2):
    err=np.asarray([1.0])
    converged=np.asarray([0])
    dim_func=np.shape(func)
    check_convergence_wrap(func,dim_func[0],threshold,N1,N2,err,converged)
    if converged[0]==0:
        conv_bool=False
    else:
        conv_bool=True
    return err[0],conv_bool
    



####################################################################
# Create the global_env class (this is what the python module sees)
####################################################################

global_env=Link()

#variables
add_global_variable(global_env, "Nbath", c_int.in_dll(libedi2py, "Nbath"), "value")
add_global_variable(global_env, "Norb", c_int.in_dll(libedi2py, "Norb"), "value")
add_global_variable(global_env, "Nspin", c_int.in_dll(libedi2py, "Nspin"), "value")
add_global_variable(global_env, "Nloop", c_int.in_dll(libedi2py, "Nloop"), "value")
add_global_variable(global_env, "Nph", c_int.in_dll(libedi2py, "Nph"), "value")
add_global_variable(global_env, "Uloc", ARRAY(c_double, 5).in_dll(libedi2py, "Uloc"), "value")
add_global_variable(global_env, "Ust", c_double.in_dll(libedi2py, "Ust"), "value")
add_global_variable(global_env, "beta", c_double.in_dll(libedi2py, "beta"), "value")
add_global_variable(global_env, "Lmats", c_int.in_dll(libedi2py, "Lmats"), "value")
add_global_variable(global_env, "Lreal", c_int.in_dll(libedi2py, "Lreal"), "value")
add_global_variable(global_env, "wini", c_double.in_dll(libedi2py, "wini"), "value")
add_global_variable(global_env, "wfin", c_double.in_dll(libedi2py, "wfin"), "value")
add_global_variable(global_env, "eps", c_double.in_dll(libedi2py, "eps"), "value")
add_global_variable(global_env, "dmft_error", c_double.in_dll(libedi2py, "dmft_error"), "value")

#functions

global_env.read_input = types.MethodType(read_input, global_env)
global_env.init_solver = types.MethodType(init_solver, global_env)
global_env.solve = types.MethodType(solve, global_env)
global_env.get_bath_dimension = types.MethodType(get_bath_dimension, global_env)
global_env.get_sigma_matsubara = types.MethodType(get_sigma_matsubara, global_env)
global_env.get_sigma_realaxis = types.MethodType(get_sigma_realaxis, global_env)
global_env.get_gimp_matsubara = types.MethodType(get_gimp_matsubara, global_env)
global_env.get_gimp_realaxis = types.MethodType(get_gimp_realaxis, global_env)
global_env.chi2_fitgf = types.MethodType(chi2_fitgf, global_env)
global_env.check_convergence = types.MethodType(check_convergence, global_env)