!ED_BATH:
integer(c_int) function get_bath_dimension_c() result(Nb) bind(c, name='get_bath_dimension')
  Nb=ed_get_bath_dimension()
end function get_bath_dimension_c

!H_REPLICA SETUP
!subroutine init_Hreplica_direct_nn(Hloc)
!  USE EDIPACK, only: ed_set_Hreplica,Norb,Nspin
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:,:,:),intent(inout) :: Hloc
!  call assert_shape(hloc,[Nspin,Nspin,Norb,Norb],"init_Hreplica_direct_nn","hloc")  
!  call ed_set_Hreplica(hloc)
!end subroutine init_Hreplica_direct_nn

!subroutine init_Hreplica_direct_so(Hloc)
!  USE EDIPACK, only: ed_set_Hreplica,Norb,Nspin
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:),intent(inout) :: Hloc
!  call assert_shape(hloc,[Nspin*Norb,Nspin*Norb],"init_Hreplica_direct_so","hloc")  
!  call ed_set_Hreplica(hloc)
!end subroutine init_Hreplica_direct_so

!subroutine init_Hreplica_symmetries_site(Hvec,lambdavec)
!  USE EDIPACK, only: ed_set_Hreplica,Norb,Nspin
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:,:,:,:),intent(inout) :: Hvec
!  real(8),dimension(:),intent(inout)         :: lambdavec
!  integer                                    :: Nsym
!  Nsym = size(lambdavec)
!  call assert_shape(Hvec,[Nspin,Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries_site","Hvec")  
!  call ed_set_Hreplica(Hvec,lambdavec)
!end subroutine init_Hreplica_symmetries_site

!subroutine init_Hreplica_symmetries_lattice(Hvec,lambdavec)
!  USE EDIPACK, only: ed_set_Hreplica,Norb,Nspin
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:,:,:,:),intent(inout) :: Hvec
!  real(8),dimension(:,:),intent(inout)       :: lambdavec
!  integer                                    :: Nsym
!  Nsym = size(lambdavec,2)
!  call assert_shape(Hvec,[Nspin,Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries_lattice","Hvec")  
!  call ed_set_Hreplica(Hvec,lambdavec)
!end subroutine init_Hreplica_symmetries_lattice




!BREAK_SYMMETRY_BATH
!subroutine break_symmetry_bath_site(bath,field,sign,save)
!  USE EDIPACK, only: ed_break_symmetry_bath
!  implicit none
!  real(8),dimension(:),intent(inout) :: bath
!  real(8),intent(in)                 :: field
!  real(8),intent(in)                 :: sign
!  logical,intent(in)                 :: save
!  call ed_break_symmetry_bath(bath,field,sign,save)
!end subroutine break_symmetry_bath_site
!
!subroutine break_symmetry_bath_ineq(bath,field,sign,save)
!  USE EDIPACK, only: ed_break_symmetry_bath
!  implicit none
!  real(8),dimension(:,:),intent(inout) :: bath
!  real(8),intent(in)                   :: field
!  real(8),intent(in)                   :: sign
!  logical,intent(in)                   :: save
!  call ed_break_symmetry_bath(bath,field,sign,save)
!end subroutine break_symmetry_bath_ineq



!SPIN_SYMMETRIZE_BATH
!subroutine spin_symmetrize_bath_site(bath,save)
!  USE EDIPACK, only: ed_spin_symmetrize_bath
!  implicit none
!  real(8),dimension(:),intent(inout) :: bath
!  logical,intent(in)                 :: save
!  call ed_spin_symmetrize_bath(bath,save)
!end subroutine spin_symmetrize_bath_site
!
!subroutine spin_symmetrize_bath_ineq(bath,save)
!  USE EDIPACK, only: ed_spin_symmetrize_bath
!  implicit none
!  real(8),dimension(:,:),intent(inout) :: bath
!  logical,intent(in)          :: save
!  call ed_spin_symmetrize_bath(bath,save)
!end subroutine spin_symmetrize_bath_ineq



!ORB_SYMMETRIZE_BATH
!subroutine orb_symmetrize_bath_site(bath,save)
!  USE EDIPACK, only: ed_orb_symmetrize_bath
!  implicit none
!  real(8),dimension(:),intent(inout) :: bath
!  logical,intent(in)        :: save
!  call ed_orb_symmetrize_bath(bath,save)
!end subroutine orb_symmetrize_bath_site
!
!subroutine orb_symmetrize_bath_ineq(bath,save)
!  USE EDIPACK, only: ed_orb_symmetrize_bath
!  implicit none
!  real(8),dimension(:,:),intent(inout) :: bath
!  logical,intent(in)          :: save
!  call ed_orb_symmetrize_bath(bath,save)
!end subroutine orb_symmetrize_bath_ineq



!ORB_EQUALITY_BATH
!subroutine orb_equality_bath_site(bath,indx,save)
!  USE EDIPACK, only: ed_orb_equality_bath
!  implicit none
!  real(8),dimension(:),intent(inout) :: bath
!  integer,intent(in)                 :: indx
!  logical,intent(in)                 :: save
!  call ed_orb_equality_bath(bath,indx,save)
!end subroutine orb_equality_bath_site
!
!subroutine orb_equality_bath_ineq(bath,indx,save)
!  USE EDIPACK, only: ed_orb_equality_bath
!  implicit none
!  real(8),dimension(:,:),intent(inout) :: bath
!  integer,intent(in)                   :: indx
!  logical,intent(in)                   :: save
!  call ed_orb_equality_bath(bath,indx,save)
!end subroutine orb_equality_bath_ineq





!PH_SYMMETRIZE_BATH
!subroutine ph_symmetrize_bath_site(bath,save)
!  USE EDIPACK, only: ed_ph_symmetrize_bath
!  implicit none
!  real(8),dimension(:),intent(inout) :: bath
!  logical,intent(in)        :: save
!  call ed_ph_symmetrize_bath(bath,save)
!end subroutine ph_symmetrize_bath_site
!
!subroutine ph_symmetrize_bath_ineq(bath,save)
!  USE EDIPACK, only: ed_ph_symmetrize_bath
!  implicit none
!  real(8),dimension(:,:),intent(inout) :: bath
!  logical,intent(in)                   :: save
!  call ed_ph_symmetrize_bath(bath,save)
!end subroutine ph_symmetrize_bath_ineq



!PH_TRANS_BATH
!subroutine ph_trans_bath_site(bath,save)
!  USE EDIPACK, only: ed_ph_trans_bath
!  implicit none
!  real(8),dimension(:),intent(inout) :: bath
!  logical,intent(in)        :: save
!  call ed_ph_trans_bath(bath,save)
!end subroutine ph_trans_bath_site
!
!subroutine ph_trans_bath_ineq(bath,save)
!  USE EDIPACK, only: ed_ph_trans_bath
!  implicit none
!  real(8),dimension(:),intent(inout) :: bath
!  logical,intent(in)        :: save
!  call ed_ph_trans_bath(bath,save)
!end subroutine ph_trans_bath_ineq






!BATH COMPONENT ROUTINES
!subroutine get_bath_component_dimension(type,ndim)
!  USE EDIPACK, only: ed_get_bath_component_dimension
!  implicit none
!  character(len=1),intent(in) :: type
!  integer,intent(out)         :: Ndim(3)
!  Ndim = ed_get_bath_component_dimension(type)
!end subroutine get_bath_component_dimension

!subroutine get_bath_component(array,bath_,type)
!  USE EDIPACK, only: ed_get_bath_component
!  implicit none
!  real(8),dimension(:,:,:),intent(inout) :: array
!  real(8),dimension(:),intent(in)        :: bath_
!  character(len=1),intent(in)            :: type
!  call ed_get_bath_component(array,bath_,type)
!end subroutine get_bath_component

!subroutine set_bath_component(array,bath_,type)
!  USE EDIPACK, only: ed_set_bath_component
!  implicit none
!  real(8),dimension(:,:,:),intent(out)   :: array
!  real(8),dimension(:),intent(in)        :: bath_
!  character(len=1),intent(in)            :: type
!  call ed_set_bath_component(array,bath_,type)
!end subroutine set_bath_component

!subroutine copy_bath_component(bathIN,bathOUT,type)
!  USE EDIPACK, only: ed_copy_bath_component
!  implicit none
!  real(8),dimension(:),intent(inout) :: bathIN
!  character(len=1),intent(in)        :: type
!  real(8),dimension(:),intent(out)   :: bathOUT
!end subroutine copy_bath_component
