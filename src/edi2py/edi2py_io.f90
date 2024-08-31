!ED_IO:
subroutine get_sigma_matsubara_site_c(sigma,dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5) bind(c, name='get_sigma_matsubara_site')
  integer(c_int),value  :: dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5
  complex(c_double_complex),dimension(dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5),intent(inout) :: sigma
  call assert_shape(sigma,[Nspin,Nspin,Norb,Norb,Lmats],"get_sigma_matsubara","Sigma")
  call ed_get_sigma_matsubara(sigma)
end subroutine get_sigma_matsubara_site_c
!
subroutine get_sigma_matsubara_ineq_c(sigma,dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5,dim_sigma_6) bind(c, name='get_sigma_matsubara_ineq')
  integer(c_int),value  :: dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5,dim_sigma_6
  integer(c_int) :: Nsites
  complex(c_double_complex),dimension(dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5,dim_sigma_6),intent(inout) :: sigma
  Nsites=size(sigma,1)
  call assert_shape(sigma,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],"get_sigma_matsubara","Sigma")
  call ed_get_sigma_matsubara(sigma,Nsites)
end subroutine get_sigma_matsubara_ineq_c



subroutine get_sigma_realaxis_site_c(sigma,dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5)  bind(c, name='get_sigma_realaxis_site')
  integer(c_int),value  :: dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5
  complex(c_double_complex),dimension(dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5),intent(inout) :: sigma
  call assert_shape(sigma,[Nspin,Nspin,Norb,Norb,Lreal],"get_sigma_realaxis","sigma")
  call ed_get_sigma_realaxis(sigma)
end subroutine get_sigma_realaxis_site_c
!
subroutine get_sigma_realaxis_ineq_c(sigma,dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5,dim_sigma_6)  bind(c, name='get_sigma_realaxis_ineq')
  integer(c_int),value  :: dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5,dim_sigma_6
  integer(c_int) :: Nsites
  complex(c_double_complex),dimension(dim_sigma_1,dim_sigma_2,dim_sigma_3,dim_sigma_4,dim_sigma_5,dim_sigma_6),intent(inout) :: sigma
  Nsites=size(sigma,1)
  call assert_shape(sigma,[Nsites,Nspin,Nspin,Norb,Norb,Lreal],"get_sigma_realaxis","sigma")
  call ed_get_sigma_realaxis(sigma,Nsites)
end subroutine get_sigma_realaxis_ineq_c



!subroutine get_gimp_matsubara_site(gimp)
!  USE EDIPACK, only: ed_get_gimp_matsubara,Nspin,Norb,Lmats
!  USE SCIFOR, only: assert_shape
!  implicit none
!  complex(8),dimension(:,:,:,:,:),intent(inout) :: gimp
!  call assert_shape(gimp,[Nspin,Nspin,Norb,Norb,Lmats],"get_gimp_matsubara","Gimp")
!  call ed_get_gimp_matsubara(gimp)
!end subroutine get_gimp_matsubara_site
!
!subroutine get_gimp_matsubara_ineq(gimp)
!  USE EDIPACK, only: ed_get_gimp_matsubara,Nspin,Norb,Lmats
!  USE SCIFOR, only: assert_shape
!  implicit none
!  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: gimp
!  integer                                         :: Nsites
!  Nsites=size(gimp,1)
!  call assert_shape(gimp,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],"get_gimp_matsubara","Gimp")
!  call ed_get_gimp_matsubara(gimp,Nsites)
!end subroutine get_gimp_matsubara_ineq



!subroutine get_gimp_realaxis_site(gimp)
!  USE EDIPACK, only: ed_get_gimp_realaxis,Nspin,Norb,Lreal
!  USE SCIFOR, only: assert_shape
!  implicit none
!  complex(8),dimension(:,:,:,:,:),intent(inout) :: gimp
!  call assert_shape(gimp,[Nspin,Nspin,Norb,Norb,Lreal],"get_gimp_realaxis","gimp")
!  call ed_get_gimp_realaxis(gimp)
!end subroutine get_gimp_realaxis_site
!
!subroutine get_gimp_realaxis_ineq(gimp)
!  USE EDIPACK, only: ed_get_gimp_realaxis,Nspin,Norb,Lreal
!  USE SCIFOR, only: assert_shape
!  implicit none
!  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: gimp
!  integer                                         :: Nsites
!  Nsites=size(gimp,1)
!  call assert_shape(gimp,[Nsites,Nspin,Nspin,Norb,Norb,Lreal],"get_gimp_realaxis","gimp")
!  call ed_get_gimp_realaxis(gimp,Nsites)
!end subroutine get_gimp_realaxis_ineq






!subroutine get_dens_site(arg)
!  USE EDIPACK, only: ed_get_dens,Norb
!  implicit none
!  real(8),dimension(:),intent(out) :: arg
!  if(size(arg)/=Norb)stop "get_dens error: size(arg)!=Norb"
!  call ed_get_dens(arg)
!end subroutine get_dens_site
!
!subroutine get_dens_ineq(arg,nlat)
!  USE EDIPACK, only: ed_get_dens,Norb
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:),intent(out) :: arg
!  integer,intent(in)                 :: nlat
!  call assert_shape(arg,[Nlat,Norb],"get_dens_ineq","arg")  
!  call ed_get_dens(arg,nlat)
!end subroutine get_dens_ineq



!subroutine get_mag_site(arg)
!  USE EDIPACK, only: ed_get_mag,Norb
!  implicit none
!  real(8),dimension(:),intent(out) :: arg
!  if(size(arg)/=Norb)stop "get_mag error: size(arg)!=Norb"
!  call ed_get_mag(arg)
!end subroutine get_mag_site
!
!subroutine get_mag_ineq(arg,nlat)
!  USE EDIPACK, only: ed_get_mag,Norb
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:),intent(out) :: arg
!  integer,intent(in)                 :: nlat
!  call assert_shape(arg,[Nlat,Norb],"get_mag_ineq","arg")  
!  call ed_get_mag(arg,nlat)
!end subroutine get_mag_ineq


!subroutine get_docc_site(arg)
!  USE EDIPACK, only: ed_get_docc,Norb
!  implicit none
!  real(8),dimension(:),intent(out) :: arg
!  if(size(arg)/=Norb)stop "get_docc error: size(arg)!=Norb"
!  call ed_get_docc(arg)
!end subroutine get_docc_site
!
!subroutine get_docc_ineq(arg,nlat)
!  USE EDIPACK, only: ed_get_docc,Norb
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:),intent(out) :: arg
!  integer,intent(in)                 :: nlat
!  call assert_shape(arg,[Nlat,Norb],"get_docc_ineq","arg")  
!  call ed_get_docc(arg,nlat)
!end subroutine get_docc_ineq



!subroutine get_eimp_site(arg)
!  USE EDIPACK, only: ed_get_eimp
!  implicit none
!  real(8),dimension(4),intent(out) :: arg
!  call ed_get_eimp(arg)
!end subroutine get_eimp_site
!
!subroutine get_eimp_ineq(arg,nlat)
!  USE EDIPACK, only: ed_get_eimp
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:),intent(out) :: arg
!  integer,intent(in)                 :: nlat
!  call assert_shape(arg,[Nlat,4],"get_eimp_ineq","arg")  
!  call ed_get_eimp(arg,nlat)
!end subroutine get_eimp_ineq



!subroutine get_doubles_site(arg)
!  USE EDIPACK, only: ed_get_doubles
!  implicit none
!  real(8),dimension(4),intent(out) :: arg
!  call ed_get_doubles(arg)
!end subroutine get_doubles_site
!
!subroutine get_doubles_ineq(arg,nlat)
!  USE EDIPACK, only: ed_get_doubles
!  USE SCIFOR, only: assert_shape
!  implicit none
!  real(8),dimension(:,:),intent(out) :: arg
!  integer,intent(in)                 :: nlat
!  call assert_shape(arg,[Nlat,4],"get_doubles_ineq","arg")  
!  call ed_get_doubles(arg,nlat)
!end subroutine get_doubles_ineq







