from edipy import edi2py
from ctypes import *
import numpy as np

#this ed class contains all the variable trackers
vars = edi2py.vars


######################################
# READ INPUT
######################################

def read_input(input_string):
    c_string = c_char_p(input_string.encode())
    edi2py.read_input(c_string)


######################################
# ED_MAIN FUNCTIONS
######################################

def init_solver(bath):
    dim_bath=np.shape(bath)
    if len(dim_bath)<2:
        edi2py.init_solver_site(bath,dim_bath[0])
    else:
        edi2py.init_solver_ineq(bath,dim_bath[0],dim_bath[1])

    
def solve(bath,hloc,sflag=True,mpi_lanc=False):
    dim_bath=np.shape(bath)
    dim_hloc=np.shape(hloc)
    if len(dim_bath)<2:
        edi2py.solve_site(bath,dim_bath[0],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],sflag)
    else:
        edi2py.solve_ineq(bath,dim_bath[0],dim_bath[1],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],dim_hloc[4],mpi_lanc)
    
######################################
# ED_BATH FUNCTIONS
######################################

get_bath_dimension = edi2py.get_bath_dimension

######################################
# ED_IO FUNCTIONS
######################################

def get_sigma_matsubara(Smats):
    dim_Smats=np.shape(Smats)
    if len(dim_Smats)==5:
        edi2py.get_sigma_matsubara_site(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4])
    else:
        edi2py.get_sigma_matsubara_site(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4],dim_Smats[5])
        
def get_sigma_realaxis(Smats):
    dim_Smats=np.shape(Smats)
    if len(dim_Smats)==5:
        edi2py.get_sigma_realaxis_site(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4])
    else:
        edi2py.get_sigma_realaxis_site(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4],dim_Smats[5])
        
######################################
# ED_BATH_FIT FUNCTIONS
######################################

def chi2_fitgf(func,bath,hloc=None,ispin=1,iorb=0):
    dim_func=np.shape(func)
    dim_bath=np.shape(bath)
    dim_hloc=np.shape(hloc)
    if(len(dim_func)==5):
        edi2py.chi2_fitgf_site(func,dim_func[0],dim_func[1],dim_func[2],dim_func[3],dim_func[4],bath,dim_bath[0],ispin,iorb)
    else:
         edi2py.chi2_fitgf_ineq(func,dim_func[0],dim_func[1],dim_func[2],dim_func[3],dim_func[4],dim_func[5],bath,dim_bath[0],dim_bath[1],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],dim_hloc[4],ispin)
         
######################################
# ED_AUX_FUNX FUNCTIONS
######################################

def check_convergence(func,threshold,N1,N2):
    err=np.asarray([1.0])
    converged=np.asarray([0])
    dim_func=np.shape(func)
    edi2py.check_convergence(func,dim_func[0],threshold,N1,N2,err,converged)
    if converged[0]==0:
        conv_bool=False
    else:
        conv_bool=True
    return err[0],conv_bool