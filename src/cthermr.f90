!-------------------------------------------------------------------------------
! CTHERMR Module for NJOY2016
!
! Written by Hunter Belanger (hunter.belanger@gmail.com)
!
! This module produces thermal stactering data for use in continuous temperature
! simulations (such as mutliphysics applications, or when many temperatures are
! required for the same nuclide). It uses the methods outlined by Pavlou and Ji
! referenced bellow. A breif outline of the produced data, and output format is
! provided here.
!
! Coherent Elastic Scattering
!
! Incoherent Elastic Scattering
! The only data required to completely describe this reaction is the list of
! Debye-Waller integrals at their respective temperatures ( W(T) ). From this
! value, the total incoherent elastic scattering cross section may be calculated
! using Eq. 7.5 from the ENDF-6 Formats Manual. Eq. 7.4 also provides a method
! to calculated the CDF for the scattering cosine.
!
! Incoherent Inelastic Scattering
!
! Pavlou, A.T., Ji, W., 2014. On-the-fly sampling of temperature-dependent
! thermal neutron scattering data for Monte Carlo simulations. Ann Nucl Energy,
! Annals of Nuclear Energy 71, 411–426.
!-------------------------------------------------------------------------------
module cthermm
  ! Provides cthermr module
  use locale
  implicit none
  private
  public cthermr

  ! Global Variables
  integer::nendf,nin,nout,nscr
  real(kr),dimension(:),allocatable::scr
  real(kr),dimension(:),allocatable::bufo,bufn

contains

  subroutine cthermr
    use mainio
    use endf
    use util
    implicit none
    ! Local Variables
    real(kr)::time

    ! Initialize module
    call timer(time)
    write(nsyso, '(/&
      &'' cthermr...compute continuous temperature thermal scattering cross '',&
      &''section data and distributions'',f8.1,''s'')') time
    write(nsyse, '(/'' cthermr...'',59x,f8.1,''s'')') time
    
    ! End timer call
    call timer(time)
  end subroutine cthermr

end module cthermm
