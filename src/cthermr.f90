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
! -----------------------------
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

  !-------------------------------------------------------------------
  ! GLOBAL VARIABLES
  !-------------------------------------------------------------------
  integer::nendf,nin,nout,nscr
  ! material number, number of coefficients in fit, number of temperatures
  integer::mat, nfit, ntemp
  real(kr),dimension(:),allocatable::scr
  real(kr),dimension(:),allocatable::bufo,bufn
  ! temperatures, walter-debye integrals
  real(kr),dimension(:),allocatable::temps,wdi

contains

  subroutine cthermr
    use mainio
    use endf
    use util
    implicit none
    !-----------------------------------------------------------------
    ! LOCAL VARIABLES
    !-----------------------------------------------------------------
    real(kr)::time
    real(kr):: a(17)
    character(4):: ta(17)
    equivalence(a(1), ta(1))
    integer:: i, j, nb, nw, nwd, nxc, nsub
    real(kr), dimension(:), allocatable:: src
    character(66) :: string

    ! Initialize module
    call timer(time)
    write(nsyso, '(/&
      &'' cthermr...compute continuous temperature thermal scattering cross '',&
      &''section data and distributions'',f8.1,''s'')') time
    write(nsyse, '(/'' cthermr...'',59x,f8.1,''s'')') time

    ! Read in input and output tapes.
    ! Only takes input ENDF and output tapes. No pre-processing required.
    ! This is because the input ENDF file should be a thermal scattering data
    ! file which is "self contained", and does not need to be reconstructed
    ! or Doppler broadened.
    read(nsysi,*) nendf, nout, mat
    if(nout.lt.0) call error('cthermr', 'Only ASCII output allowed.',' ')

    ! Open input and output files
    call openz(nendf,0)
    call openz(nout,1)

    ! Rewind files to the beginning
    call repoz(nendf)
    call repoz(nout)

    ! Get tapeid
    call tpidio(nendf,0,0,a,nb,nw)

    ! Find the beginning of the desired mat in file (only relivant if there
    ! are multiple materials per file, which is allowed but uncommon
    call findf(mat, 1, 451, nendf)
    
    ! Read 4 CONT records (first is a HEAD)
    do i = 1, 4
      call contio(nendf,0,0,a,nb,nw)
      if (i.eq.3) then
        nsub = n1h
      end if
    end do

    ! Check nsub to make sure this is a thermal scattering neutron file
    if (nsub.ne.12) then
      string = 'Input tape is not thermal scattering neutron data'
      call error('cthermr',string,'')
    end if

    ! The number of text cards before directory is now in the n1 variable,
    ! and is refered to as nwd. The number of entries in the dictionary is
    ! now in n2h, and will be put into nxc.
    nwd = n1h
    nxc = n2h

    ! Read all text lines
    do i = 1, nwd
      call tpidio(nendf,0,0,a,nb,nw)
    end do

    ! Allocate src
    allocate(src(npage+50))
    nw = nxc 
    ! Read all nxc dictionary records.
    ! Keep notes as to what records are provided
    call dictio(nendf,0,0,src,nb,nw)

    do i = 1, nxc
      write(*,*) src((i-1)*6 + 3), src((i-1)*6 + 4), src((i-1)*6 + 5), src((i-1)*6 + 6)
    end do

    ! Read First Card
    !read(nsysi,*) mat, nfit
    !if(nfit.lt.10) call error('cthermr', 'Must have a fit of at least 10.',' ')
    
    ! End timer call
    call timer(time)
  end subroutine cthermr

end module cthermm
