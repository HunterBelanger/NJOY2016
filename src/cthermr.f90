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
! ---------------------------
!
! Incoherent Elastic Scattering
! -----------------------------
! The only data required to completely describe this reaction is the list of
! Debye-Waller integrals at their respective temperatures ( W(T) ). From this
! value, the total incoherent elastic scattering cross section may be calculated
! using Eq. 7.5 from the ENDF-6 Manual. Eq. 7.4 also provides a method to
! calculated the CDF for the scattering cosine.
!
! Incoherent Inelastic Scattering
! -------------------------------
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
  ! PRIVATE MODULE VARIABLES
  !-------------------------------------------------------------------
  integer:: mat, nin, nout
  integer:: lat, lasym, lln, ns, ni, ntemps, nbetas, nalphas
  integer:: lthr ! Determines coherent elastic (1) or incoherent elastic (2)
  real(kr):: b(24) ! Max len of B(N) is 24 from ENDF manual p. 151
  real(kr), dimension(:), allocatable:: alpha, beta, temps, energies
  real(kr), dimension(:,:,:), allocatable:: s_a_b_t

contains
  !============================================================================
  ! Main CTHERMR subroutine
  subroutine cthermr
    use mainio
    use endf
    use util
    !-----------------------------------------------------------------
    ! LOCAL VARIABLES
    !-----------------------------------------------------------------
    real(kr):: time
    real(kr):: a(17)
    character(4):: ta(17)
    equivalence(a(1), ta(1))
    integer:: nwd, nxc, nsub ! Numbers from ENDF format
    integer:: i, nb, nw ! nb and nw from NJOY dev info in manual
    real(kr), dimension(:), allocatable:: dict
    logical :: elastic = .false., inco_inelastic = .false.
    

    ! Initialize module
    call timer(time)
    write(nsyso, '(/&
      &'' cthermr...computes continuous temperature thermal scattering '',&
      &''data'',f8.1,''s'')') time
    write(nsyse, '(/'' cthermr...'',59x,f8.1,''s'')') time

    ! Read in input and output tapes.
    ! Only takes input ENDF and output tapes. No pre-processing required.
    ! This is because the input ENDF file should be a thermal scattering data
    ! file which is "self contained", and does not need to be reconstructed
    ! or Doppler broadened.
    read(nsysi,*) nin, nout, mat
    if(nout.lt.0) call error('cthermr', 'Only ASCII output allowed.',' ')
    ! TODO write input, output, and mat to nsyso

    ! Open input and output files
    call openz(nin,0)
    call openz(nout,1)

    ! Rewind files to the beginning
    call repoz(nin)
    call repoz(nout)

    ! Get tapeid
    call tpidio(nin,0,0,a,nb,nw)

    ! Find the beginning of the desired mat in file (only relivant if there
    ! are multiple materials per file, which is allowed but uncommon
    call findf(mat, 1, 451, nin)
    
    ! Read 4 CONT records (first is a HEAD)
    do i = 1, 4
      call contio(nin,0,0,a,nb,nw)
      if (i.eq.3) then
        nsub = n1h
      end if
    end do

    ! Check nsub to make sure this is a thermal scattering neutron file
    if (nsub.ne.12) then
      call error('cthermr','Input tape does not contain thermal neutron scattering data','')
    end if

    ! The number of text cards before directory is now in the n1 variable,
    ! and is refered to as nwd. The number of entries in the dictionary is
    ! now in n2h, and will be put into nxc.
    nwd = n1h
    nxc = n2h

    ! Read all text lines
    do i = 1, nwd
      call tpidio(nin,0,0,a,nb,nw)
    end do

    ! Allocate dict
    allocate(dict(nxc * 6))
    nw = nxc 
    ! Read all nxc dictionary records (use 6*nxc for zeros).
    ! Keep notes as to what records are provided
    call dictio(nin,0,0,dict,nb,nw)

    ! check for reaction types in dictionary
    do i = 1, nxc
      if ((dict((i-1)*6 + 3) .eq. 7) .and. (dict((i-1)*6 + 4) .eq. 4)) then
        ! Check for incoherent inelastic scattering MF=7, MT=4. This should
        ! ALWAYS be present, but I make sure just in case.
        inco_inelastic = .true.
      else if ((dict((i-1)*6 + 3) .eq. 7) .and. (dict((i-1)*6 + 4) .eq. 2)) then
        ! Coherent and Incoherent elastic are both in MT=2, because there never
        ! are both present. It is only one or the other. It must be checked latter
        ! if it is Coherent or Incoherent
        elastic = .true.
      end if
    end do
    ! TODO write which reactions are contained in file to nsyso

    if (elastic) then
      call process_elastic_data()
    end if

    if (inco_inelastic) then
      call read_s_a_b_t_table()
    end if

    ! Read First Card
    !read(nsysi,*) mat, nfit
    !if(nfit.lt.10) call error('cthermr', 'Must have a fit of at least 10.',' ')

    ! Deallocate arrays
    if (allocated(dict)) deallocate(dict)
    if (allocated(alpha)) deallocate(alpha)
    if (allocated(beta)) deallocate(beta)
    if (allocated(temps)) deallocate(temps)
    if (allocated(energies)) deallocate(energies)
    if (allocated(s_a_b_t)) deallocate(s_a_b_t)
    
    ! End timer call
    call timer(time)
  end subroutine cthermr
  
  !============================================================================
  ! Subroutine to determine Coherent or Incoherent Elastic
  subroutine process_elastic_data()
    use mainio
    use endf
    use util
    integer:: i, nb, nw
    real(kr):: a(17)
    
    ! Find MF 7 MT 4 in input file
    call findf(mat, 7, 2, nin)

    ! Read in head
    call contio(nin, 0, 0, a, nb, nw)
    lthr = l1h

    ! Process data depending on representation

  end subroutine process_elastic_data

  !============================================================================
  ! Subourinte to read in entire S(a,b,T) table, with a, b, and T grids
  subroutine read_s_a_b_t_table()
    use mainio
    use endf
    use util
    integer:: i, b_i, a_i, t_i, nb, nw, loc, low
    real(kr):: tmp(5*npage)

    ! Find MF 7 MT 4 in input file
    call findf(mat, 7, 4, nin)

    ! Read head
    call contio(nin, 0, 0, tmp, nb, nw)
    lat = l2h
    lasym = n1h

    ! Read list of B(N)
    call listio(nin, 0, 0, tmp, nb, nw)
    if (nb .ne. 0) then
      call error('cthermr','B(N) list longer than 24. Should not be possible','')
    else
      ni = n1h
      ns = n2h
      do i = 1, ni
        b(i) = tmp(6+i)
      end do
    end if
    

    ! Read tab2 with info on beta
    call tab2io(nin, 0, 0, tmp, nb, nw)
    nbetas = n2h
    allocate(beta(nbetas))

    ! Read in data for all betas
    do b_i = 1, nbetas
      loc = 1
      call tab1io(nin, 0, 0, tmp(loc), nb, nw)
      if(b_i .eq. 1) then
        nalphas = n2h
        ntemps = l1h + 1
        low = 6 + 2*n1h
        allocate(alpha(nalphas))
        allocate(temps(ntemps))
        allocate(s_a_b_t(nbetas,nalphas,ntemps))
      end if
      ! TODO should not use such an adhoc method. Should fix one day
      if (nalphas > 2*5*npage) then
        call error('cthermr','Cant fit all alphas in tmp','')
      end if
      loc = loc + nw
      do while (nb .ne. 0)
        call moreio(nin, 0, 0, tmp(loc), nb, nw)
        loc = loc + nw
      end do
      beta(b_i) = c2h

      ! Assign all alphas for T0
      a_i = 1
      do i = 1, 2*nalphas, 2
        ! Only save alphas for first run
        if (b_i .eq. 1) then
          alpha(a_i) = tmp(low+i)
          temps(1) = c1h
        else
          if(alpha(a_i) .ne. tmp(low+i)) &
            & call error('cthermr','Apha dissagrement','')

          if(temps(1) .ne. c1h) &
            & call error('cthermr','Temperature dissagrement','')
        end if
        s_a_b_t(b_i, a_i, 1) = tmp(low+i+1)
        a_i = a_i + 1
      end do

      !read in other lists for all temperatures
      do t_i = 2, ntemps
        loc = 1
        call listio(nin, 0, 0, tmp(loc), nb, nw)
        loc = loc + nw
        do while (nb .ne. 0)
          call moreio(nin, 0, 0, tmp(loc), nb, nw)
          loc = loc + nw
        end do

        if (b_i .eq. 1) then
          temps(t_i) = c1h
        else
          if(temps(t_i) .ne. c1h) &
            & call error('cthermr','Temperature dissagrement','')
        endif

        do a_i = 1, nalphas
          s_a_b_t(b_i, a_i, t_i) = tmp(6+a_i)
        end do
      end do
    end do
    
    write(nsyso, '(/''Number of Beta values  : '',i4,'''')') nbetas
    write(nsyso, '(''Number of Alpha values : '',i4,'''')') nalphas
    write(nsyso, '(''Number of Temperatures : '',i4,'''')') ntemps
    write(nsyso, '(/''Betas:''/6(E12.6,X))') (beta(i), i = 1, nbetas)
    write(nsyso, '(/''Alphas:''/6(E12.6,X))') (alpha(i), i = 1, nalphas)
    write(nsyso, '(/''Temperatures:''/6(E12.6,X))') (temps(i), i = 1, ntemps)
  end subroutine read_s_a_b_t_table

end module cthermm
