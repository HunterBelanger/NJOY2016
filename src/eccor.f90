module eccom
  ! provides subroutine eccor for NJOY2016
  use locale
  use mainio ! provides nsysi,contio,nsyso,nsyse
  use endf   ! provides endf routines and variables
  use util   ! provides timer,openz,repoz,error
  implicit none
  private
  public eccor

  ! file units
  integer::ning,ninr,nine,nour,nout

  ! error message string
  character(len=131) hsmg

  ! .. scalars in parameter ..
  !
  ! ngrmax - maximum number of neutron groups
  ! nggmax - maximum number of gamma groups
  ! tmpmax - maximum number of temperatures
  ! maxord - maximum order of sub-group data
  ! maxmtx - maximum number of matrix types
  ! maxleg - maximum number of legendre orders
  ! maxxs  - maximum number of cross-section types
  ! maxnfr - maximum number of response function types
  ! maxel  - maximum number of elements per library
  ! maxthr - maximum number of thermal scattering groups
  ! maxmth - maximum number of reaction number of thermal scattering data
  ! ngblk  - maximum number of group blocks
  ! maxnbf - size of buffer array of integer numbers
  ! maxcbf - size of buffer array of characters
  ! maxbuf - size of buffer array of real numbers
  !
  integer    ngrmax,nggmax,     tmpmax,            maxord,maxmtx, &
             maxleg,            maxxs,maxnfr,      maxel, &
             maxthr,ngblk,      msize,             psize,maxcbf, &
             maxnbf,            maxbuf,maxmth
  parameter (ngrmax=1970,nggmax=200,tmpmax=15,maxord=12,maxmtx=3,maxleg=5, &
  & maxxs=15,maxnfr=40,maxel=50,maxthr=200,ngblk=150,msize=10000000, &
  & psize=msize*maxleg,maxcbf=maxel*2+maxxs+maxmtx+maxnfr,maxmth=10, &
  & maxnbf=3000000,maxbuf=3000000)

  ! common values
  real(kr) xsec(maxxs,tmpmax,ngrmax),mubar(ngrmax),spec(ngrmax), &
  & weight(maxord,tmpmax,ngrmax),prosig(maxxs,maxord,tmpmax,ngrmax), &
  & rf(maxnfr,ngrmax),chid(ngrmax),nud(ngrmax),nu5(ngrmax),egrid(ngrmax+1)
  integer lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  character listfr(maxnfr)*16,listmi(maxxs)*16,listmx(maxmtx)*16,sycor(maxel)*16
  integer iprint ! print level
  logical delay,prompt
  integer nel,nelp,nfg,ngi,nfre,nsmc,npmc,matno,nsmic,nsmtx,nfr
  integer ntemp(maxel),iel(maxel),ismc(maxxs),ifre(maxnfr),ipmc(maxxs), &
  & ismtx(maxmtx),npne(maxmtx),inmin(maxmtx,ngrmax),inmax(maxmtx,ngrmax), &
  & ishld(tmpmax,ngrmax),nsmx(maxel),itemp(maxel)
  integer ntty
  integer ntotal,nelas,ninel,nnufis,ntran,nfissn,ncapt,melas,mnxn,minel, &
  & selas,sinel,mtherm,nnxn,rn2nd,rn2n,rn3n,rn4n,rn2np,rn3np,rn2na,rn3na, &
  & rn2n2a,rnany,ntherm(maxmth)
  logical pmt16,pmt103,pmt104,pmt105,pmt106,pmt107,found

contains

  subroutine eccor
  !-----------------------------------------------------------------
  !
  !     Produce an eccolib interface file from Njoy intermediate
  !   cross-section library.
  !
  !     The eccolib format provide an efficient way to store multi-
  !   group isotopic nuclear data to be used in the ecco code,
  !
  !     Working from thermr and groupr output tape,this module
  !   produce the standard eccolib interface file.
  !
  ! Copyright:
  !  Copyright (C) 2025 Newcleo LFR
  !  This library is free software; you can redistribute it and/or
  !  modify it under the terms of the BSD license.
  !
  !---input specifications (free format)----
  !
  ! card 1 file units
  !   ning      input gendf unit (>0 formatted / <0 unformatted)
  !   ninr      input unit of reference table (=0 for create option)
  !   nine      input unit of eccolib (=0 for create option)
  !   nour      output unit of reference table
  !   nout      output unit of eccolib
  !   iprint    long print option (0/1=minimum/maximum) (default=0)
  !
  ! card 2 eccolib structure
  !   lrecl     record length (bytes) for library (usually 4096)
  !   lwordi    number of bytes in an integer word (usually 4)
  !   lwordc    cumber of characters in a string (usually 16)
  !   lwordr    number of bytes in a real(kr) word (usually 4)
  !
  ! card 3 group block structure (2 cards)
  !   nreci     number of group blocks
  !   ing       array of integer, number of groups stored together in each
  !             group block. The sum of nreci values must equal the total
  !             number of energy groups in the eccolib.
  !
  ! card 4 title (if nine = 0)
  !   title     char(len=80) description of this library
  !
  ! card 5 data source for element
  !   sysour    char(len=16), data source for this element
  !
  ! card 6 material number
  !   matno     material number to be processed from GENDF tape
  !
  ! card 7 material number
  !   type      type of nuclide (either 'FERTILE', 'FISSILE' or
  !             'NON-FISSILE')
  !
  ! card 8 reaction number of thermal scattering matrices (2 cards)
  !   mtherm    number of MT indices
  !   ntherm    mtherm values of MT incices (221 is the value for free
  !             gas scattering)
  !
  ! card 9 isotope card
  !   iel       identifier number for this isotope
  !   awt       atomic weight
  !   sycor     char(len=16) isotope name
  !   sybel     char(len=16) element name
  !
  ! card 10 isotope properties
  !   ecng      non-gamma capture energy (mev)
  !   ecg       gamma capture energy (mev)
  !   efng      non-gamma fission energy (mev)
  !   efg       gamma fission energy (Mev)
  !   dis       disintegration constant (s-1)
  !
  ! card11 vector cross section names (2 cards) (if nine = 0)
  !   nsmic2    number of x-section type names. If nsmic2=0, use 'TOTAL',
  !             'TRANSPORT','CAPTURE','ELASTIC','INELASTIC','FISSION',
  !             'NU*FISSION','N,XN','AXIAL_TRANSPORT','RADIAL_TRANSPORT'
  !   listmi    user-defined char(len=16) names (if nsmic2 > 0)
  !
  ! card12 matrix cross section names (2 cards) (if nine = 0)
  !   nsmtx     number of x-section type names. If nsmtx=0, use 'ELASTIC',
  !             'INELASTIC','N,XN'
  !   listmx    user-defined char(len=16) names (if nsmtx > 0)
  !
  ! card13 response function names (2 cards) (if nine = 0)
  !   nfr       number of response function names. If nsmtx=0, use 'N,P',
  !             'N,ALPHA','N,2N','N,3N','N,4N','N,GAMMA','N,D','N,T',
  !             'N,HE3','N,2ALPHA','N,3ALPHA','N,2P',N,P+ALPHA',
  !             'N,T+2ALPHA','N,D+2ALPHA','N,F','N,N+F','N,2N+F','N,3N+F',
  !             'N,N+ALPHA','N,N+3ALPHA','N,2N+ALPHA','N,3N+ALPHA',
  !             'N,N+P','N,N+2ALPHA','N,2N+2ALPHA','N,N+D','N,N+T',
  !             'N,N+HE3','N,N+D+2ALPHA','N,N+T+2ALPHA','N,2NP','N,3NP',
  !             'N,PD','N,PT,’N,ANYTHING’
  !   listfr    user-defined char(len=16) names (if nfr > 0)
  !
  ! card14 response function names (2 cards) (if nine = 0)
  !   nsgp      number of gamma production x-section names
  !   listgp    user-defined char(len=16) names (if nsgp > 0)
  !
  ! card15 response function names (if nine = 0)
  !   nmpn      maximum order of pn
  !
  ! card16      gamma energy mesh (2 cards) (if nine = 0)
  !   nfgg      number of gamma energy groups
  !   eg        energy limits (nfgg+1 values)
  !  ...
  !
  !            cards 9 to 16 are set for the first isotope in the library.
  !            Repeat cards 9 and 10 for all other isotopes present in the
  !            library.
  !
  !-----------------------------------------------------------------
  ! .. local scalars ..
  integer nbuff(1024),i,ifis,irect,lrecl,lsblk,lwordc,lwordi,lwordr, &
  & maxblk,maxm,maxp,mbuf,mreci,nfis,nreci,nrect,ntfis,thdim1,thrmin
  logical lcreate
  character fertil*4,fissil*4,nonfis*4,chdate*8,sysour*16,type*16,title*80, &
  & text6*6
  ! .. local variables and arrays ..
  real(kr) time
  integer lenrec(3),alloc_ok
  ! .. intrinsic functions ..
  intrinsic max
  ! .. data statements ..
  data fertil/'FERT'/
  data fissil/'FISS'/
  data nonfis/'NON-'/
  character word*10
  integer idate(8)
  ! .. allocatable arrays ..
  integer,allocatable,dimension(:) :: ing,lblock,mblock,pblock,sblock, &
  & mstart,pstart,istart,lunit
  integer,allocatable,dimension(:,:) :: tnmin,tnmax
  real(kr),allocatable,dimension(:) :: eng,matrix,matrx2,porder,pordr2, &
  & pordr3
  real(kr),allocatable,dimension(:,:) :: totinm
  real(kr),allocatable,dimension(:,:,:) :: therma
  real(kr),allocatable,dimension(:,:,:,:) :: torder
  !
  call timer(time)
  write(nsyso,'(/'' eccor...produce an eccolib format output file'',23x,f8.1, &
  & ''s'')') time
  write(nsyse,'(/'' eccor...'',60x,f8.1,''s'')') time
  !----
  !  scratch storage allocation
  !----
  allocate(ing(ngblk),lunit(maxmtx),tnmin(tmpmax,maxthr), &
  & tnmax(tmpmax,maxthr),stat=alloc_ok)
  if(alloc_ok /= 0) call error('eccor','unable to allocate int arrays',' ')
  allocate(eng(ngrmax+1),therma(maxthr,maxthr,tmpmax),totinm(ngrmax,maxmtx), &
  & torder(maxthr,maxthr,maxleg,tmpmax),stat=alloc_ok)
  if(alloc_ok /= 0) call error('eccor','unable to allocate real arrays',' ')
  ing(:ngblk) = 0
  tnmin(:tmpmax,:maxthr) = 0
  tnmax(:tmpmax,:maxthr) = 0
  totinm(:ngrmax,:maxmtx) = 0.0d0
  therma(:maxthr,:maxthr,:tmpmax) = 0.0d0
  torder(:maxthr,:maxthr,:maxleg,:tmpmax) = 0.0d0
  !----
  !  initialise variables
  !----
  matno = 0
  melas = 0
  minel = 0
  mnxn = 0
  mtherm = 0
  nel = 0
  nelas = 0
  nfg = 0
  nfissn = 0
  ncapt = 0
  nfr = 0
  nfre = 0
  ngi = 0
  ninel = 0
  nnxn = 0
  npmc = 0
  nsc = 0
  nsh = 0
  nsmc = 0
  nsmic = 0
  nsmtx = 0
  nsp = 0
  ntty = 0
  nnufis = 0
  ntotal = 0
  ntran = 0
  rn2nd = 0
  rn2n = 0
  rn3n = 0
  rn4n = 0
  rn2np = 0
  rn3np = 0
  rn2na = 0
  rn3na = 0
  rn2n2a = 0
  rnany = 0
  selas = 0
  sinel = 0
  nelp = 0
  prompt = .false.
  delay = .false.
  !----
  !  initialise arrays
  !----
  chid(:ngrmax)=0.0d0
  egrid(:ngrmax+1)=0.0d0
  mubar(:ngrmax)=0.0d0
  nud(:ngrmax)=0.0d0
  nu5(:ngrmax)=0.0d0
  prosig(:maxxs,:maxord,:tmpmax,:ngrmax)=0.0d0
  rf(:maxnfr,:ngrmax)=0.0d0
  spec(:ngrmax)=0.0d0
  weight(:maxord,:tmpmax,:ngrmax)=0.0d0
  xsec(:maxxs,:tmpmax,:ngrmax)=0.0d0
  ntherm(:maxmth) = 0
  ntemp(:maxel)=0
  iel(:maxel)=0
  ifre(:maxnfr)=0
  inmax(:maxmtx,:ngrmax)=0
  inmin(:maxmtx,:ngrmax)=0
  ipmc(:maxxs)=0
  ishld(:tmpmax,:ngrmax)=0
  ismc(:maxxs)=0
  ismtx(:maxmtx)=0
  npne(:maxmtx)=0
  lnrec1(:ngblk*5,:2,:maxmtx)=0
  lnrec2(:ngblk*5,:2,:maxmtx)=0
  nsmx(:maxel)=0
  itemp(:maxel)=0
  !----
  !  initialise local variables
  !----
  ifis = 0
  irect = 0
  lrecl = 0
  lsblk = 0
  lwordc = 0
  lwordi = 0
  lwordr = 0
  maxblk = 0
  maxm = 0
  maxp = 0
  mbuf = 0
  mreci = 0
  nfis = 0
  nine = 0
  ning = 0
  ninr = 0
  nour = 0
  nout = 0
  nreci = 0
  nrect = 0
  ntfis = 0
  thdim1 = 0
  iprint = 0
  ! set defaults
  iverf = 5
  thrmin = ngrmax
  ! start run
  !
  ! read input data
  read(nsysi,*) ning,ninr,nine,nour,nout,iprint
  if(iprint>0) then
    write(nsyso,9020)
    write(nsyso,9040) ngrmax,maxbuf,tmpmax,maxord,maxmtx,maxleg,maxxs, &
    maxnfr,maxel,maxthr,ngblk,msize,psize,maxnbf,maxcbf
  endif
  call openz(ning,0)
  if(nine==0) then
    lcreate = .true.
    if(iprint>0) write(nsyso,9170)
  else
    lcreate = .false.
    if(iprint>0) write(nsyso,9180)
  endif
  ! read in record length (bytes) for library and word size
  ! for int, real and char
  read(nsysi,*) lrecl,lwordi,lwordc,lwordr
  if((iprint>0).and.(nine/=0)) then
    write(nsyso,9160) ninr
    write(nsyso,9150) nine
  endif
  ! record length (in words) for each data type
  lenrec(1) = lrecl/lwordi
  lenrec(2) = lrecl/lwordc
  lenrec(3) = lrecl/lwordr
  if(iprint>0) then
    write(nsyso,9060) nour
    write(nsyso,9050) nout,lrecl,lenrec
  endif
  ! read in blocking structure for library
  read(nsysi,*) nreci
  if(nreci>ngblk) then
    ! dimensions exceeded
    write(nsyso,9010) ngblk,nreci
    call error('eccor','maximum number of group blocks exceeded',' ')
  endif
  maxblk = 0
  mreci = 0
  read(nsysi,*) (ing(i),i=1,nreci)
  do i = 1,nreci
    maxblk = max(maxblk,ing(i))
    mreci = mreci + ing(i)
  enddo
  if(lcreate) then
    ! read in title for library
    read(nsysi,*) title
    if(iprint>0) write(nsyso,9130) title
  endif
  ! read in data source for element
  read(nsysi,*) sysour
  if(iprint>0) write(nsyso,9140) sysour
  ! read in material number
  read(nsysi,*) matno
  if(iprint>0) write(nsyso,9110) matno
  ! nfis: number of fissile isotopes
  ! read in type of input nuclide
  read(nsysi,*) type
  if(iprint>0) write(nsyso,9090) type
  if(type(1:4)==fertil) then
    ifis = 2
    nfis = 1
    ntfis = 1
  else if(type(1:4)==fissil) then
    ifis = 1
    nfis = 1
    ntfis = 1
  else if(type(1:4)==nonfis) then
    ifis = 0
  else
    write(nsyso,9100) type
    call error('eccor','keyword not recognised',' ')
  endif
  ! read in reaction number of thermal scatter
  read(nsysi,*) mtherm
  if(mtherm/=0) then
    if(mtherm>maxmth) call error('eccor','maxmth overflow',' ')
    read(nsysi,*) (ntherm(i), i=1,mtherm)
    if(iprint>0) write(nsyso,9120) (ntherm(i), i=1,mtherm)
  endif
  ! check if partials of parasitic absortpion reactions are present
  call ectpmt
  !
  !if(iprint>0) write(nsyso,9000) msize,psize
  if(lcreate) then
    ! no input ecco library so create new library
    write(text6,'(4htape,i2.2)') nout
    open(nout,file=text6,access='direct',status='unknown',recl=lrecl, &
    & form='unformatted')
    write(text6,'(4htape,i2.2)') nour
    open(nour,file=text6,access='direct',status='unknown',recl=lrecl, &
    & form='unformatted')
    ! ecco disk record reference array
    ! lrec(1,n): data block type (1-integer,2-character,3-real,4-double p.)
    ! lrec(2,n): start record for the nth data block
    ! lrec(3,n): number of words in the nth data block
    ! nrect is the current ecco direct access data record
    ! create the reference data package
    call ecfdat(nrect,nour,lenrec)
    !
    ! scratch storage allocation
    allocate(matrix(msize),matrx2(msize),porder(psize),pordr2(psize))
    !
    ! create the ecco general contents package
    maxm = 0
    maxp = 0
    call ecncon(nrect,lenrec,eng,nfis,ifis,ntfis,title,sysour,nreci,ing, &
    & matrix,matrx2,porder,pordr2,maxm,maxp,lunit,mreci,lwordr,thrmin, &
    & therma,torder,tnmin,tnmax,thdim1)
    !
    ! scratch storage deallocation
    deallocate(pordr2,porder,matrx2,matrix)
    !
    ! scratch storage allocation
    allocate(lblock(maxblk),mstart(maxblk+1),istart(maxblk+1),mblock(maxblk), &
    & pblock(maxblk),sblock(maxblk),pstart(maxblk))
    allocate(matrix(maxm),porder(maxp),pordr2(maxp),pordr3(maxp))
    !
    ! create the ecco x-section package
    call ecosec(nrect,lenrec,nreci,ing,matrix,maxm,porder,maxp,pordr2,pordr3, &
    & lunit,maxblk,lblock,mstart,istart,mblock,pblock,sblock,pstart,totinm, &
    & thrmin,therma,torder,tnmin,tnmax,thdim1)
    !
    ! scratch storage deallocation
    deallocate(pordr3,pordr2,porder,matrix)
    deallocate(pstart,sblock,pblock,mblock,istart,mstart,lblock)
  else
    ! input ecco library exists so add data to it
    ! the input reference file on unit ninr
    ! the output reference file on unit nour
    ! the input cross-section file on unit nine
    ! the output cross-section file on unit nout
    write(text6,'(4htape,i2.2)') ninr
    open(ninr,file=text6,access='direct',status='old',recl=lrecl, &
    form='unformatted')
    write(text6,'(4htape,i2.2)') nour
    open(nour,file=text6,access='direct',status='unknown',recl=lrecl, &
    form='unformatted')
    write(text6,'(4htape,i2.2)') nine
    open(nine,file=text6,access='direct',status='old',recl=lrecl, &
    form='unformatted')
    write(text6,'(4htape,i2.2)') nout
    open(nout,file=text6,access='direct',status='unknown',recl=lrecl, &
    form='unformatted')
    ! read the reference data package
    nrect = 0
    call refin(nrect,lenrec)
    !
    ! scratch storage allocation
    allocate(matrix(msize),matrx2(msize),porder(psize),pordr2(psize))
    maxm = 0
    maxp = 0
    call ecnadd(irect,nrect,lenrec,eng,nfis,ifis,ntfis,sysour,nreci,ing, &
    & matrix,matrx2,porder,pordr2,maxm,maxp,lunit,lwordr,thrmin,therma, &
    & torder,tnmin,tnmax,thdim1)
    !
    ! scratch storage deallocation
    deallocate(pordr2,porder,matrx2,matrix)
    !
    ! scratch storage allocation
    allocate(lblock(maxblk),mstart(maxblk+1),istart(maxblk+1),mblock(maxblk), &
    & pblock(maxblk),sblock(maxblk),pstart(maxblk))
    allocate(matrix(maxm),porder(maxp),pordr2(maxp),pordr3(maxp))
    !
    ! extend the ecco general contents package
    call ecosad(irect,nrect,lenrec,nreci,ing,matrix,maxm,porder,maxp,pordr2, &
    & pordr3,lunit,maxblk,lblock,mstart,istart,mblock,pblock,sblock,pstart, &
    & totinm,thrmin,therma,torder,tnmin,tnmax,thdim1)
    !
    ! scratch storage deallocation
    deallocate(pordr3,pordr2,porder,matrix)
    deallocate(pstart,sblock,pblock,mblock,istart,mstart,lblock)
  endif
  !
  ! scratch storage deallocation
  deallocate(torder,totinm,therma,eng)
  deallocate(tnmax,tnmin,lunit,ing)
  !
  ! write a zero record at end of library
  nbuff(:1024)=0
  call ecwint(nout,lenrec(1),nbuff,lenrec(1),nrect,1)
  call timer(time)
  if(iprint>0) write(nsyso,9070) nrect
  !
  ! eccor is complete
  call repoz(ning)
  call closz(ning)
  call closz(nour)
  call closz(nout)
  if(ninr/=0) call closz(ninr)
  if(nine/=0) call closz(nine)
  write(nsyso,'(69x,f8.1,''s''/1x,77(''*''))') time
  return
  !----
  !  formats
  !----
  9000 format (' msize=',i10,' psize =',i10)
  9010 format (' error - maximum number of group blocks exceeded',/, &
  ' ngblk = ',i6,' nreci = ',i6)
  9020 format (/,10x,'gendf* to ecco library conversion ',/)
  9040 format (/,10x,'values of parameters used in this version',/,7x, &
  'ngrmax     maxbuf         tmpmax     maxord     maxmtx',/,5i12,/,7x, &
  'maxleg     maxxs          maxnfr     maxel      maxthr',/,5i12,/,7x, &
  'ngblk      msize          psize      maxnbf     maxcbf',/,5i12)
  9050 format (' output library written to                         : tape', &
  i2.2,/, ' record length                                     :',i6,/, &
  ' record length in words for each data type (i,c,r) :',3i6)
  9060 format (' output reference file written to                  : tape',i2.2)
  9070 format (1x,'maximum eccolib record reached:',i8)
  9090 format (' isotope type                                      :',a)
  9100 format (' *** keyword not recognised ***',/,' keyword read was ',a, &
  ' but expected either fertile,',' fissile or non-fissile')
  9110 format (' material number                                   :',i5)
  9120 format (' thermal scattering data on reaction number        :',i5)
  9130 format (' title                                             :',a)
  9140 format (' element origin                                    :',a)
  9150 format (' input library read from                           : tape',i2.2)
  9160 format (' input reference file read from                    : tape',i2.2)
  9170 format (/,' an ecco library is to be created',/,1x,32 ('='),/)
  9180 format (/,' a new isotope is to be added to an existing ecco', &
  ' library',/,1x,56('='),/)
  end subroutine eccor
  !
  subroutine ectpmt
  !
  !  program checks tape ning to see if it contains components of
  !  parasitic absorption or (n,2n) thus:
  !    name      combined         partials
  !    (n,2n)       16             876-891
  !    (n,p)       103             600-649
  !    (n,d)       104             650-699
  !    (n,t)       105             700-749
  !    (n,he3)     106             750-799
  !    (n,alpha)   107             800-849
  !
  !    assumes if levels are present the lower partial reaction is present
  !    but it is possible to just include continuum - the higher partial.
  !    under both circumstances the variable associated with the "combined"
  !    should be set to true
  !
  call repoz(ning)
  !
  !  components of (n,p)
  call find(matno,3,600,ning,found)
  pmt103=found
  if(.not.found)then
    call find(matno,3,649,ning,found)
    pmt103=found
  endif
  !
  !  components of (n,d)
  call find(matno,3,650,ning,found)
  pmt104=found
  if(.not.found)then
    call find(matno,3,699,ning,found)
    pmt104=found
  endif
  !
  !  components of (n,t)
  call find(matno,3,700,ning,found)
  pmt105=found
  if(.not.found)then
    call find(matno,3,749,ning,found)
    pmt105=found
  endif
  !
  !  components of (n,he3)
  call find(matno,3,750,ning,found)
  pmt106=found
  if(.not.found)then
    call find(matno,3,799,ning,found)
    pmt106=found
  endif
  !
  !  components of (n,alpha)
  call find(matno,3,800,ning,found)
  pmt107=found
  if(.not.found)then
    call find(matno,3,849,ning,found)
    pmt107=found
  endif
  !
  !  components of (n,2n)
  call find(matno,3,875,ning,found)
  pmt16=found
  if(.not.found)then
    call find(matno,3,891,ning,found)
    pmt16=found
  endif
  call repoz(ning)
  end subroutine ectpmt
  !
  subroutine ecfdat(nrect,nour,lenrec)
  ! code to create the reference data package
  ! read the ecco reference file from card images
  ! .. scalar arguments ..
  integer nour,nrect
  ! .. array arguments ..
  integer lenrec(3)
  ! .. local scalars ..
  real(kr) avog,er,etop,fcwe,rn,testeg,width
  integer i,ig,ipos,j,lbuff,nbl,ndata,next,nfgg,nfrecs,ngopt,ngr,nmpn, &
  nopt,npri,nrange,ns,nsgp,ntopg,loc,nb,nw,nz,nfr2,nmisc2,nsmtx2,nsmic2
  logical found
  ! .. local arrays ..
  integer nbuff(maxnbf)
  real(kr) buff(maxnbf),awt(maxel),dis(maxel),ecg(maxel),ecng(maxel), &
  & efg(maxel),efng(maxel),eg(nggmax)
  integer lfr(maxnfr),lrec(3,4)
  character listgp(5)*16,sybel(maxel)*16
  character(len=16) cbuff(maxcbf)
  ! .. intrinsic functions ..
  intrinsic abs,exp,float,log
  character*16 midat(10),mxdat(3),frdat(maxnfr)
  integer  lfrdat(maxnfr)
  data midat /'TOTAL','TRANSPORT','CAPTURE','ELASTIC','INELASTIC','FISSION', &
  'NU*FISSION','N,XN','AXIAL_TRANSPORT','RADIAL_TRANSPORT'/
  data mxdat /'ELASTIC','INELASTIC','N,XN'/
  data frdat /'N,P','N,ALPHA','N,2N','N,3N','N,4N','N,GAMMA','N,D','N,T', &
  & 'N,HE3','N,2ALPHA','N,3ALPHA','N,2P','N,P+ALPHA','N,T+2ALPHA', &
  & 'N,D+2ALPHA','N,F','N,N+F','N,2N+F','N,3N+F','N,N+ALPHA','N,N+3ALPHA', &
  & 'N,2N+ALPHA','N,3N+ALPHA','N,N+P','N,N+2ALPHA','N,2N+2ALPHA','N,N+D', &
  & 'N,N+T','N,N+HE3','N,N+D+2ALPHA','N,N+T+2ALPHA','N,2NP','N,3NP','N,PD', &
  & 'N,PT','N,ANYTHING','N,2N+D','N,N+2P','N,N+P+ALPHA','N,D+ALPHA'/
  data lfrdat /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
  0,0,0,0,0,0,0,0,0/
  !----
  !  read input data
  !----
  nel = 1
  !  iel   - isotope identifier
  !  awt   - atomic weight
  ! sycor  - isotope name
  ! sybel  - element name
  read(nsysi,*) iel(nel),awt(nel),sycor(nel),sybel(nel)
  if(iprint>0) then
    write(nsyso,9040) nel
    write(nsyso,9050) iel(nel),awt(nel),sycor(nel),sybel(nel)
  endif
  ! ecng - non-gamma capture energy (Mev)
  ! ecg  - gamma capture energy (Mev)
  ! efng - non-gamma fission energy (Mev)
  ! efg  - gamma fission energy (Mev)
  ! dis  - disintegration constant (s-1)
  read(nsysi,*) ecng(nel),ecg(nel),efng(nel),efg(nel),dis(nel)
  if(iprint>0) then
    write(nsyso,9070)
    write(nsyso,9080) ecng(nel),ecg(nel),efng(nel),efg(nel),dis(nel)
  endif
  ! avog  - avogadro's number
  ! fcwe  - mev to watts conversion factor
  ! er    - slowing down energy (mev)
  avog = 0.602295
  fcwe = 1.602E-13
  er = 0.0
  if(iprint>0) write(nsyso,9090) avog,fcwe,er
  ! nsmic  - number of x-section type names
  ! listmi - x-section type name
  read(nsysi,*) nsmic2
  if(nsmic2==0) then
    nsmic = 10
    if(iprint>0) write(nsyso,9100) nsmic
    do ns = 1,nsmic
      listmi(ns) = midat(ns)
      if(iprint>0) write(nsyso,9110) ns,listmi(ns)
    enddo
  else
    nsmic = nsmic2
    if(iprint>0) write(nsyso,9100) nsmic
    do ns = 1,nsmic
      read(nsysi,*) listmi(ns)
      if(iprint>0) write(nsyso,9110) ns,listmi(ns)
    enddo
  endif
  ! nsmtx  - number of x-section matrix type names
  ! listmx - x-section type matrix name
  read(nsysi,*) nsmtx2
  if(nsmtx2==0) then
    nsmtx = 3
    if(iprint>0) write(nsyso,9120) nsmtx
    do ns = 1,nsmtx
      listmx(ns) = mxdat(ns)
      if(iprint>0) write(nsyso,9110) ns,listmx(ns)
    enddo
  else
    nsmtx = nsmtx2
    if(iprint>0) write(nsyso,9120) nsmtx
    do ns = 1,nsmtx
      read(nsysi,*) listmx(ns)
      if(iprint>0) write(nsyso,9110) ns,listmx(ns)
    enddo
  endif
  ! nfr    - number of response function names
  ! listfr - response function name
  ! lfr    - r.f. is a combination of the next lfr types
  read(nsysi,*) nfr2
  if(nfr2==0) then
    nfr = maxnfr
    if(iprint>0) write(nsyso,9130) nfr
    do ns = 1,nfr
      listfr(ns) = frdat(ns)
      lfr(ns) = lfrdat(ns)
      if(iprint>0) write(nsyso,9140) ns,listfr(ns),lfr(ns)
    enddo
  else
    nfr = nfr2
    if(nfr>maxnfr)then
      write(nsyso,9131) nfr,maxnfr
      call error('ecfdat','dimensions exceeded',' ')
    else
      read(nsysi,*) (listfr(ns),ns=1,nfr)
      read(nsysi,*) (lfr(ns),ns=1,nfr)
      if(iprint>0) then
        write(nsyso,9130) nfr
        do ns = 1,nfr
          write(nsyso,9140) ns,listfr(ns),lfr(ns)
        enddo
      endif
    endif
  endif
  ! nsgp   - number of gamma production x-section names
  read(nsysi,*) nsgp
  if(iprint>0) write(nsyso,9150) nsgp
  if(nsgp>0) then
    ! listgp - gamma production x-section name
    if(iprint>0) write(nsyso,9155)
    do ns = 1,nsgp
      read(nsysi,*) listgp(ns)
      if(iprint>0) write(nsyso,9110) ns,listgp(ns)
   enddo
  endif
  ! nmpn  - maximum order of pn
  npri=1
  read(nsysi,*) nmpn
  if(iprint>0) write(nsyso,9160) nmpn
  !----
  !  recover neutron energy mesh from gendf file
  !----
  call skiprz(ning,1)
  call find(iel(nel),1,451,ning,found)
  if(.not.found) call error('ecfdat','fail to recover the energ'//'y mesh',' ')
  call contio(ning,0,0,buff,nb,nw)
  nz=nint(buff(4))
  call listio(ning,0,0,buff(1),nb,nw)
  loc=1+nw
  do while (nb/=0)
    if(loc+302>maxnbf) call error('ecfdat','maxnbf overflow(1)',' ')
    call moreio(ning,0,0,buff(loc),nb,nw)
    loc=loc+nw
  enddo
  nfg=nint(buff(3))
  do ig=1,nfg+1
    egrid(nfg-ig+2)=buff(7+nz+ig)
  enddo
  ndata = nfg + 1
  if(iprint>0) then
    if(npri/=0) then
      write(nsyso,9000)
      write(nsyso,9010) (i,egrid(i),log(egrid(i)/egrid(i+1)),i=1,nfg)
    endif
  endif
  egrid(:ndata) = egrid(:ndata)/1.e6
  if(iprint>0) write(nsyso,9180) nfg
  !----
  !  recover gamma energy mesh from input
  !----
  ! nfgg   - number of fine gamma groups
  read(nsysi,*) nfgg
  if(nfgg>nggmax) call error('ecfdat','nggmax overflow',' ')
  if(iprint>0) write(nsyso,9190) nfgg
  if(nfgg>0) then
    ! eg(i): upper energy (mev) in group i (gammas)
    ! eg(nfgg+1) is the lower bound of the group scheme
    read(nsysi,9170) (eg(i),i=1,nfgg + 1)
  endif
  ! record reference array entries:
  ! record block 1: identifier and record reference array
  ! record block 2: integer data
  ! record block 3: character data
  ! record block 4: real data
  ! record block 1: identifier and record reference array
  ! data block type: (integer)
  lrec(1,1) = 1
  ! start record:
  lrec(2,1) = 1
  ! number of words stored:
  lrec(3,1) = 9 + 3*4
  ! record block 2: integer data
  nbuff(1) = nel
  lbuff = 1
  nbuff(lbuff+1:lbuff+nel) = iel(:nel)
  lbuff = lbuff + nel
  nbuff(lbuff+1) = nsmic
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nsmtx
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nfr
  lbuff = lbuff + 1
  nbuff(lbuff+1:lbuff+nfr) = lfr(:nfr)
  lbuff = lbuff + nfr
  nbuff(lbuff+1) = nsgp
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nmpn
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nfg
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nfgg
  lbuff = lbuff + 1
  if(lbuff>maxnbf) call error('ecfdat','maxnbf overflow(2)',' ')
  ! block 2: (integer)
  lrec(1,2) = 1
  ! block 2: start record
  lrec(2,2) = 2
  ! block 2: number of words
  lrec(3,2) = lbuff
  nrect     = 1
  call ecwint(nour,lenrec(1),nbuff,lbuff,nrect,1)
  ! record block 3: character data
  lbuff=0
  cbuff(lbuff+1:lbuff+nel) = sycor(:nel)
  lbuff = nel
  cbuff(lbuff+1:lbuff+nel) = sybel(:nel)
  lbuff = lbuff + nel
  cbuff(lbuff+1:lbuff+nsmic) = listmi(:nsmic)
  lbuff = lbuff + nsmic
  cbuff(lbuff+1:lbuff+nsmtx) = listmx(:nsmtx)
  lbuff = lbuff + nsmtx
  cbuff(lbuff+1:lbuff+nfr) = listfr(:nfr)
  lbuff = lbuff + nfr
  if(nfgg>0) then
    cbuff(lbuff+1:lbuff+nfgg) = listgp(:nfgg)
    lbuff = lbuff + nfgg
  endif
  ! block 3: (character)
  lrec(1,3) = 2
  ! block 3: start record
  lrec(2,3) = nrect + 1
  ! block 3: number of data words
  lrec(3,3) = lbuff
  call ecwc16(nour,lenrec(2),cbuff,lbuff,nrect,1)
  ! record block 4: real data
  buff(:nel) = awt(:nel)
  lbuff = nel
  buff(lbuff+1:lbuff+nel) = ecng(:nel)
  lbuff = lbuff + nel
  buff(lbuff+1:lbuff+nel) = ecg(:nel)
  lbuff = lbuff + nel
  buff(lbuff+1:lbuff+nel) = efng(:nel)
  lbuff = lbuff + nel
  buff(lbuff+1:lbuff+nel) = efg(:nel)
  lbuff = lbuff + nel
  buff(lbuff+1:lbuff+nel) = dis(:nel)
  lbuff = lbuff + nel
  buff(lbuff+1) = avog
  buff(lbuff+2) = er
  buff(lbuff+3) = fcwe
  lbuff = lbuff + 3
  buff(lbuff+1:lbuff+nfg+1) = egrid(:nfg+1)
  lbuff = lbuff + nfg + 1
  if(nfgg>0) then
    buff(lbuff+1:lbuff+nfgg+1) = eg(:nfgg+1)
    lbuff = lbuff + nfgg + 1
  endif
  ! block 4: (reals)
  lrec(1,4) = 3
  ! block 4: start record
  lrec(2,4) = nrect + 1
  ! block 4: number of words
  lrec(3,4) = lbuff
  call ecwr4(nour,lenrec(3),buff,lbuff,nrect,1)
  ! finally, create identifier record block
  ! type number of this package
  nbuff(1) = 1
  ! name number of this package
  nbuff(2) = 1
  ! father number of this package
  nbuff(3) = 1
  ! structure number
  nbuff(4) = 1
  ! number of physical records in this package
  nbuff(5) = nrect
  ! number of fortran records in this package
  nfrecs = 4
  nbuff(6) = nfrecs
  ! library origin - 2 = jef
  nbuff(7) = 2
  ! zero unused variables
  nbuff(8) = 0
  nbuff(9) = 0
  lbuff = 9
  do j = 1,nfrecs
    nbuff(lbuff+1) = lrec(1,j)
    nbuff(lbuff+2) = lrec(2,j)
    nbuff(lbuff+3) = lrec(3,j)
    lbuff = lbuff + 3
  enddo
  nrect = 0
  call ecwint(nour,lenrec(1),nbuff,lbuff,nrect,1)
  ! reset nrect to maximum value
  nrect = nbuff(5)
  if(iprint>0) then
    write(nsyso,9200) nrect
    ! print out data block structure
    write(nsyso,9210)
    write(nsyso,9220)
    do nbl = 1,nfrecs
      write(nsyso,9230) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    enddo
  endif
  return
  !----
  !  formats
  !----
  9000 format (/,30x,'reference group structure',/,30x,25 ('='),/, &
  ' group      energy         width  ',' group      energy         width  ', &
  ' group      energy         width',/)
  9010 format (3 (1x,i5,f16.5,f12.5))
  9020 format (1x,'error calculating fine group energy bounds ',1p,2e14.5)
  9040 format (/,5x,'ecco reference file - number of elements:',i4,/,1x, &
  'element',7x,'atomic',9x,'isotope',9x,'element',/,1x,' number',7x,'weight', &
  9x,'  name ',9x,'  name ',/)
  9050 format (2x,i6,f12.3,12x,2a16)
  9070 format (/,3x,'non-gamma',5x,' gamma ',3x,'non-gamma',5x,' gamma ',5x, &
  'disint-',/,3x,' capture ',5x,'capture',3x,' fission ',5x,'fission',5x, &
  'egration',/,3x,'  energy ',5x,' energy',3x,'  energy ',5x,' energy',5x, &
  'constant',/,3x,'  (mev)  ',5x,' (mev) ',3x,'  (mev)  ',5x,' (mev) ',5x, &
  ' (s-1)  ',/)
  9080 format (1x,1p,5e12.5)
  9090 format (/,1x,32havogadro's number              :,1p,e14.6,/,1x, &
  'mev to watts conversion factor :',e14.6,/,1x, &
  'slowing down energy (mev)      :',e14.6)
  9100 format (/,1x,'number of x-section type names:',i6,/,1x,' type ',10x, &
  ' name',/,1x,'number',/)
  9110 format (1x,i4,12x,a16)
  9120 format (/,1x,'number of x-section matrix type names:',i6,/,1x,' type ', &
  10x,' name',/,1x,'number',/)
  9130 format (/,1x,'number of response function names:',i6,/,1x,' type ',10x, &
  ' name',15x,'lfr code',/,1x,'number',/)
  9131 format (/,1x,'dimensions exceeded',/, &
  ' number of response function names:',i6,/,1x,' dimensions allow: ',i6)
  9140 format (1x,i4,12x,a16,2x,i8)
  9150 format (/,1x,'number of gamma production x-section names:',i6)
  9155 format (/,1x,' type ',10x,' name',/,1x,'number',/)
  9160 format (/,1x,'maximum order of pn:',i6)
  9170 format (/,1x,' calculation of fine group structure ',/, &
  ' 1/lethargy interval:',i6,/,1x,'groups above 10 mev:',i6,/,1x, &
  'groups below 10 mev:',i6)
  9180 format (/,1x,'number of fine neutron groups:',i6)
  9190 format (/,1x,'number of fine gamma groups:  ',i6)
  9200 format (/,1x,'last record of reference data:',i6)
  9210 format (/,1x,'reference data record block structure',/)
  9220 format (2x,'f rec',4x,'data',4x,'start',4x,' words ',/,2x,' no. ',4x, &
  'type',3x,'record',3x,' stored ',/)
  9230 format (2x,i5,3i9)
  9240 format (1x,1p,e14.6,3x,i4,7x,i1,4x,e14.6)
  9250 format (' energy ranges',/,'   top energy    no. groups  option')
  end subroutine ecfdat
  !
  subroutine ecncon(nrect,lenrec,eng,nfis,ifis,ntfis,title,sysour,nreci,ing, &
  & matrix,matrx2,porder,pordr2,maxm,maxp,lunit,mreci,lwordr,thrmin,therma, &
  & torder,tnmin,tnmax,thdim1)
  ! code to create the general contents package
  !
  ! arguments
  !  nrect - record number on general constants file
  !  eng - broad group energy boundaries
  ! general contents data block structure:
  ! record block 1: integers; (identifier)
  ! record block 2: integers;
  ! record block 3: integers;
  ! record block 4: characters
  ! record block 5: reals
  ! record block 5+ifis: reals
  ! last record block 5+nfis: reals
  ! .. scalar arguments ..
  integer ifis,lwordr,maxm,maxp,mreci,nfis,nreci,nrect,ntfis,thdim1,thrmin
  character sysour*16,title*80
  ! .. array arguments ..
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize), &
  & eng(ngrmax+1),therma(maxthr,maxthr,tmpmax), &
  & torder(maxthr,maxthr,maxleg,tmpmax)
  integer ing(nreci),lenrec(3),lunit(maxmtx),tnmax(tmpmax,maxthr), &
  tnmin(tmpmax,maxthr)
  ! .. local scalars ..
  real(kr) buff(maxnbf),atemp,comp,emax,emin,fe,totspc,za
  integer nbuff(maxnbf),i,ielp,ielr,if,itempj,j,l,lbuff,msmc,msmtx,nb,nbl, &
  & ncards,ndil,nfrecs,ngg,npnfp,nreg,nsmxj,nspge,ntempj,ntw,nw,nword2,mg, &
  & mbuf,mat
  logical fision,found
  data fision/.false./,found/.false./
  ! .. local arrays ..
  real(kr) b(ngrmax),btemp(10),flux(ngrmax),sdil(10)
  integer lrec(3,maxel),nmax(ngrmax+1)
  character(len=16) cbuff(maxcbf)
  ! .. intrinsic functions ..
  intrinsic abs
  !----
  !  initialise local scalars and arrays
  !----
  mbuf = maxnbf
  comp = 0.0
  atemp = 0.0
  emax = 0.0
  emin = 0.0
  fe = 0.0
  totspc = 0.0
  za = 0.0
  i = 0
  ielp = 0
  ielr = 0
  if = 0
  itempj = 0
  j = 0
  l = 0
  lbuff = 0
  msmc = 0
  msmtx = 0
  nb = 0
  nbl = 0
  ncards = 0
  ndil = 0
  nfrecs = 0
  ngg = 0
  npnfp = 0
  nreg = 0
  nsmxj = 0
  nspge = 0
  ntempj = 0
  ntw = 0
  nw = 0
  nword2 = 0
  b(:ngrmax)=0.0d0
  btemp(:10)=0.0d0
  flux(:ngrmax)=0.0d0
  sdil(:10)=0.0d0
  lrec(:3,:maxel)=0
  nmax(:ngrmax+1)=0
  ! zero array spec
  spec(:ngi)=0.0d0
  nelp = 1
  ntemp(nelp) = 0
  itemp(nelp) = 0
  !----
  !  read input gendf
  !  read mf=1 on gendf from unit ning
  !----
  call repoz(ning)
  call tpidio(ning,0,0,buff,nb,nw)
  10 call contio(ning,0,0,buff,nb,nw)
  if(mfh/=0) then
    ntemp(nelp) = ntemp(nelp) + 1
    ntempj = ntemp(nelp)
    ! file 1 only read for first temperature
    za = c1h
    if(n1h/=-1) then
      write(nsyso,9000)
      call error('ecncon','** not a gendf tape **',' ')
    endif
    ntw = n2h
    call ecfil1(sdil,ndil,ntw,eng,btemp(ntempj),ngrmax)
    if(ngi/=mreci) then
      write(nsyso,9010) mreci,ngi
      call error('ecncon','wrong number of groups on input gendf*',' ')
    endif
    if(iprint>0) write(nsyso,9020) btemp(ntempj)
    if(ntempj==1) then
      emin = eng(1)
      emax = eng(ngi+1)
      if(iprint>0) then
        write(nsyso,9030) ning,za,math,ndil,ngi,emin,emax
        write(nsyso,9040) (sdil(l),l=1,ndil)
      endif
      !  calculate array nmax for the general contents file
      !  nmax(l) - fine group number with the same upper neutron
      !            energy boundary as the lth group of the present
      !            group structure
      ! eccolib group structure is stored in mev so convert gendf structure
      eng(:ngi+1) = eng(:ngi+1)/1.e6
      i  = ngi + 1
      j  = 1
      !  modified not to check the bottom boundary
      ! it is 1.1e-4 in uk and 1.00001e-5 in france!
      do l = 1,nfg
        if(abs((egrid(l)-eng(i))/egrid(l))<0.0005) then
          nmax(j) = l
          i  = i - 1
          j  = j + 1
        else if(egrid(l)<eng(i)) then
          write(nsyso,9450) l,egrid(l),eng(i)
          call error('ecncon','input group boundary is not a fine group boun' &
          & //'dary',' ')
        endif
      enddo
    endif
    ! read the gendf to set constants
    ! read cross-section data
    thdim1 = ngi - maxthr + 1
    if(thdim1<=0) thdim1 = 1
    mat = matno
    call find(mat,3,0,ning,found)
    if(found) then
      call ecfil3(flux,eng)
    else
      write(nsyso,9480)
      call error('ecncon','no file 3 (cross section) data found on gendf* ' &
      & //'tape', ' ')
    endif
    itempj = itemp(nelp)
    if(ntempj==1) then
      ! read delayed chi from file 5 if required
      mat = matno
      call find(mat,5,0,ning,found)
      if(found) call ecfil5
      ! read matrices from file 6
      mat = matno
      call find(mat,6,0,ning,found)
      if(found) then
        call ecfil6(nreci,ing,flux,lenrec,lwordr,matrix,matrx2,porder, &
        pordr2,maxm,maxp,lunit,itempj,thrmin,fision)
        if(iprint>0) write(nsyso,9510) itempj,thrmin
        if(ngi-thrmin>maxthr) then
          write(nsyso,9500) maxthr,ngi - thrmin
          call error('ecncon','ngi-thrmin>maxthr',' ')
        endif
      else if(iprint>0) then
        write(nsyso,9460)
      endif
      ! read temperature dependent thermal scattering matrix if present
    else if(itempj>0 .and. ntempj>1) then
      mat = matno
      call find(mat,6,0,ning,found)
      if(found) then
        thdim1 = ngi - maxthr + 1
        if(thdim1<=0) thdim1 = 1
        call ecfil6t(therma,thrmin,tnmin,tnmax,torder,thdim1)
      else
        ! thermal scattering data not present
        write(nsyso,9470)
        call error('ecncon','** temperature dependant thermal scattering matr' &
        & //'ices have been requested but are not available on gendf * **',' ')
      endif
    endif
    ! read sub-group data from ecfil50
    mat = matno
    call find(mat,50,0,ning,found)
    if(found) call ecfil50
    call tomend(ning,0,0,buff)
    call skiprz(ning,-1)
    call contio(ning,0,0,buff,nb,nw)
  endif
  if(math==-1) then
    go to 40
  else
    go to 10
  endif
  ! end of gendf
  !----
  !  assign the general contents parameters (fine group structure)
  !----
  40 nrect = 0
  ! nreg: number of cell regions
  ! nelp: number of elements in this group structure
  ! ncards: number of cards read as input
  ! ngi: number of groups in this structure
  nreg = 1
  ncards = 0
  ! store start record for identifier record block
  lrec(2,1) = nrect + 1
  ! allow 1 record for identifier record block storage
  nrect = nrect + 1
  ! nmax(ngi): fine group number with same upper boundary
  ! ielr: identifier reference array by region
  ! ielp: identifier reference array by position
  ielp = nelp
  ielr = ielp
  if(iprint>0) then
    write(nsyso,9070) ielr
    write(nsyso,9080) ielp
  endif
  ! enter data for record block number 2:
  ! nreg: number of cell regions
  ! nelp: number of elements in this group structure
  ! ncards: number of cards read as input
  ! ngi: number of groups in this structure
  ! put into buffer array ready for disk write
  nbuff(1) = nreg
  nbuff(2) = nelp
  nbuff(3) = ncards
  nbuff(4) = ngi
  if(iprint>0) write(nsyso,9060) nreg,nelp,ncards,ngi
  ! nmax(ngi): fine group number with same upper boundary
  ! put data into disk buffer
  nbuff(5:4+ngi) = nmax(:ngi)
  ! ntemp: number of temperatures for each element
  if(iprint>0) write(nsyso,9090) ntempj
  ! if thermal scattering data has been requested but there is only
  ! one temperature on the input gendf then itemp(nelp) should be
  ! set to 0 on the output library
  if(ntempj==1) then
    itempj = 0
    itemp(nelp) = 0
  endif
  nsmxj = nsmx(nelp)
  if(iprint>0) then
    write(nsyso,9100) itempj
    write(nsyso,9110) nfis
    ! ifis: table of fissile isotope indicators
    write(nsyso,9120) ifis
    ! ntfis: number of fission spectra for each element
    write(nsyso,9130) ntfis
    ! nsmc: number of x-section types per element
    write(nsyso,9140) nsmc
    ! npmc: number of x-section types per element with subgroups
    write(nsyso,9150) npmc
    ! nsmx: number of x-section matrix types per element
    write(nsyso,9160) nsmxj
    ! msmc: number of macroscopic x-section types
    ! msmtx: number of macroscopic x-section matrices
  endif
  msmc = 0
  msmtx = 0
  ngg = 0
  nspge = 0
  npnfp = 0
  if(iprint>0) then
    write(nsyso,9170) msmc,msmtx
    write(nsyso,9180) nfre ! number of response function types per element
    write(nsyso,9190) ngg ! number of gamma groups
    write(nsyso,9200) nspge
    write(nsyso,9210) npnfp ! pn order (+1) to which fluxes are stored
  endif
  ! update buffer array
  nbuff(ngi+5) = ielr
  nbuff(ngi+6) = ielp
  nbuff(ngi+7) = ntempj
  nbuff(ngi+8) = itempj
  nbuff(ngi+9) = nfis
  nbuff(ngi+10) = ifis
  nbuff(ngi+11) = ntfis
  nbuff(ngi+12) = nsmc
  nbuff(ngi+13) = npmc
  nbuff(ngi+14) = nsmxj
  nbuff(ngi+15) = msmc
  nbuff(ngi+16) = msmtx
  nbuff(ngi+17) = nfre
  nbuff(ngi+18) = npnfp
  nbuff(ngi+19) = ngg
  nbuff(ngi+20) = nspge
  ! block:  2 (integers)
  ! number of data words in block
  nword2 = ngi + 20
  ! data block 2: (integer)
  lrec(1,2) = 1
  ! data block 2: start record
  lrec(2,2) = nrect + 1
  ! data block 2: number of words
  lrec(3,2) = nword2
  call ecwint(nout,lenrec(1),nbuff,nword2,nrect,1)
  ! block:  3 (integers)
  ! start address in nbuff
  lbuff = 0
  ! ismc(nsmc): reference list number of x-section type per element
  if(nsmc>0) then
    if(iprint>0) write(nsyso,9220) (ismc(i),i=1,nsmc)
    ! update library buffer
    nbuff(lbuff+1:lbuff+nsmc) = ismc(:nsmc)
    lbuff = lbuff + nsmc
  endif
  ! ipmc(npmc): reference list number of x-section types with sub-
  !             groups per element
  if(npmc>0) then
    if(iprint>0) write(nsyso,9230) (ipmc(i),i=1,npmc)
    ! update library buffer
    nbuff(lbuff+1:lbuff+npmc) = ipmc(:npmc)
    lbuff = lbuff + npmc
  endif
  ! ismtx(nsmxj): reference list number of x-section type per elem.
  if(nsmxj>0) then
    if(iprint>0) write(nsyso,9240) (ismtx(i),i=1,nsmxj)
    ! update disk buffer
    nbuff(lbuff+1:lbuff+nsmxj) = ismtx(:nsmxj)
    lbuff = lbuff + nsmxj
  endif
  ! ifre(nfre): reference list number of x-section type per elem.
  if(nfre>0) then
    if(iprint>0) write(nsyso,9250) (ifre(i),i=1,nfre)
    ! update disk buffer
    nbuff(lbuff+1:lbuff+nfre) = ifre(:nfre)
    lbuff = lbuff + nfre
  endif
  ! jfr(nfre): value type indicator for response functions
  ! set to value of ifre.
  if(nfre>0) then
    if(iprint>0) write(nsyso,9260) (ifre(i),i=1,nfre)
    nbuff(lbuff+1:lbuff+nfre) = ifre(:nfre)
    lbuff = lbuff + nfre
  endif
  ! npne(nsmxj): pn orders; per matrix type, per element
  if(nsmxj>0) then
    if(iprint>0) write(nsyso,9270) (npne(i),i=1,nsmxj)
    ! update disk buffer
    nbuff(lbuff+1:lbuff+nsmxj) = npne(:nsmxj)
    lbuff = lbuff + nsmxj
  endif
  ! data block 3: (integer)
  lrec(1,3) = 1
  ! start record:
  lrec(2,3) = nrect + 1
  ! number of words:
  lrec(3,3) = lbuff
  call ecwint(nout,lenrec(1),nbuff,lbuff,nrect,1)
  ! block:  4 (characters)
  if(iprint>0) write(nsyso,9280) title
  ! transfer data to ecco buffer
  cbuff(1) = title(1:16)
  cbuff(2) = title(17:32)
  cbuff(3) = title(33:48)
  cbuff(4) = title(49:64)
  cbuff(5) = title(65:80)
  lbuff = 5
  ! read x-section data origin for each element
  if(iprint>0) write(nsyso,9290)
  if(iprint>0) write(nsyso,9300) sycor(ielp),sysour
  ! transfer to ecco buffer
  cbuff(lbuff+1) = sysour
  lbuff = lbuff + 1
  ! data block 4: (character)
  lrec(1,4) = 2
  ! start record:
  lrec(2,4) = nrect + 1
  ! number of words:
  lrec(3,4) = lbuff
  call ecwc16(nout,lenrec(2),cbuff,lbuff,nrect,1)
  ! block:  5 (reals)
  ! start address in nbuff
  lbuff = 0
  ! comp: compositions of nuclides in region 1
  comp = 1.0
  if(iprint>0) write(nsyso,9310) comp
  ! update disk buffer
  buff(lbuff+1) = comp
  lbuff = lbuff + 1
  ! atemp: temperature of this region
  atemp = 0.0
  if(iprint>0) write(nsyso,9320) atemp
  ! update disk buffer
  buff(lbuff+1) = atemp
  lbuff = lbuff + 1
  ! b(ngi): =1.0 for reference group structure
  b(:ngi) = 1.0
  ! update disk buffer
  buff(lbuff+1:lbuff+ngi) = b(:ngi)
  lbuff = lbuff + ngi
  ! btemp(ntemp): tabulated temperatures for each element
  if(iprint>0) write(nsyso,9340)
  ! number of temperatures for this nuclide
  if(iprint>0) write(nsyso,9350) ielp, (btemp(i),i=1,ntempj)
  ! update disk buffer
  buff(lbuff+1:lbuff+ntempj) = btemp(:ntempj)
  lbuff = lbuff + ntempj
  ! data block 5:
  lrec(1,5) = 3
  ! start record:
  lrec(2,5) = nrect + 1
  ! number of words:
  lrec(3,5) = lbuff
  call ecwr4(nout,lenrec(3),buff,lbuff,nrect,1)
  ! data for ecco blocks: 6 -> 5 + nfis
  ! transfer each fission spectrum separately
  ! currently only one spectrum stored per element
  if(nfis>0) then
    if(nnufis/=0) then
      if = 1
      ! find 4ev and set the fission spectrum below it to 1e-20
      do mg = ngi,1,-1
        if(egrid(mg)>4.0001e-6) go to 147
        spec(mg)=1e-20
      enddo
      147 totspc = 0.0
      do i = 1,ngi
        buff(i) = spec(i)
        totspc = totspc + buff(i)
      enddo
      if(iprint>0) write(nsyso,9440) totspc
      if(abs(totspc-1.0)>0.0001) then
        if(iprint>0) write(nsyso,9430) totspc
        ! normalise fission spectra to sum to 1.0
        buff(:ngi) = buff(:ngi)/totspc
      endif
      if(iprint>0) write(nsyso,9360) if
      if(iprint>0) write(nsyso,9370) (buff(i),i=1,ngi)
      ! ef: fission energy
      fe = 0.0
      if(.not.fision) fe = 1.e+6
      if(iprint>0) write(nsyso,9380) fe
      ! update disk buffer
      buff(ngi+1) = fe
      lbuff = ngi + 1
      ! each fission spectrum is in a new data block
      ! block 5+if: (real)
      lrec(1,5+if) = 3
      ! start record:
      lrec(2,5+if) = nrect + 1
      ! number of words:
      lrec(3,5+if) = lbuff
      call ecwr4(nout,lenrec(3),buff,lbuff,nrect,1)
    else
      write(nsyso,9490)
      call error('ecncon','nu bar not on input gendf* tape so nu*' &
      & //'fission can not be formed',' ')
    endif
  endif
  !----
  !  finally, create identifier record block
  !----
  ! data type is integer:
  lrec(1,1) = 1
  ! number of words in first record block
  lrec(3,1) = 9 + 3* (5+nfis)
  ! type number of this package
  nbuff(1) = 2
  ! name number of this package
  nbuff(2) = 2
  ! father number of this package
  nbuff(3) = 1
  ! structure number
  nbuff(4) = 1
  ! number of physical records in this package
  nbuff(5) = nrect - lrec(2,1) + 1
  ! number of fortran records in this package
  nfrecs = 5 + nfis
  nbuff(6) = nfrecs
  ! library origin - 2 = jef
  nbuff(7) = 2
  ! zero unused variables
  nbuff(8) = 0
  nbuff(9) = 0
  lbuff = 9
  do j = 1,nfrecs
    nbuff(lbuff+1:lbuff+3) = lrec(:3,j)
    lbuff = lbuff + 3
  enddo
  ! reset record counter to beginning of package
  nrect = lrec(2,1) - 1
  call ecwint(nout,lenrec(1),nbuff,lbuff,nrect,1)
  ! reset nrect to maximum value
  nrect = lrec(2,1) + nbuff(5) - 1
  if(iprint>0) then
    write(nsyso,9390) nrect
    ! print out data block structure
    write(nsyso,9400)
    write(nsyso,9410)
    do nbl = 1,nfrecs
      write(nsyso,9420) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    enddo
  endif
  return
  !----
  !  formats
  !----
  9000 format (/,'***** error - not a gendf tape *****')
  9010 format (/,' sum of groups in group blocks - ',i5,' differs', &
  ' from number of groups on input gendf* - ',i5)
  9020 format (/,25x,'temperature ',f8.2,/,25x,20 ('='))
  9030 format (/,' gendf read from unit              :',i7,/, &
  ' material identifier (za)          :',f7.0,/, &
  ' material identifier (mat)         :',i7,/, &
  ' number of dilutions               :',i7,/, &
  ' number of groups                  :',i7,4x,1p,' from ',e10.4,' to ',e10.4, &
  & ' ev')
  9040 format (/,' dilutions                            :',1p,10e12.4)
  9050 format (/,' group structure                      :',1p,/,(8e14.6))
  9060 format (/,' nreg,nelp,ncards,ngi                 :',1x,4i6,/)
  9070 format (/,' ielr                                 :', (1x,10i6))
  9080 format (/,' ielp                                 :', (1x,10i6))
  9090 format (/,' number of temperatures - ntemp       :', (1x,10i6))
  9100 format (/,' group number indicating that from this group',' onwards',/, &
  & '  the elastic scattering matrix is temperature',' dependent - itemp:',i6)
  9110 format (/,' number of fissile elements - nfis   :',1x,i6)
  9120 format (/,' type number of element (0=non-fissile, 1=fissile, ', &
  & '2=fertile) - ifis :', (1x,10i6))
  9130 format (/,' number of fission spectra - ntfis   :', (1x,10i6))
  9140 format (/,' number of cross-section types - nsmc:', (1x,10i6))
  9150 format (/,' number of sub-group types - npmc    :', (1x,10i6))
  9160 format (/,' number of matrix - nsmx             :', (1x,10i6))
  9170 format (/,' msmc,msmtx:',2i6)
  9180 format (/,' number of response functions - nfre :', (1x,10i6))
  9190 format (/,' number of gamma groups - ngg        :',1x,i6)
  9200 format (/,1x,'number of gamma production matrices - nspge:',i6)
  9210 format (/,1x,'pn order +1 to which the fluxes are stored',' - npnfp:',i6)
  9220 format (/,' ismc                                :', (10i6))
  9230 format (/,' ipmc                                :', (10i6))
  9240 format (/,' ismtx                               :', (10i6))
  9250 format (/,' ifre                                :', (10i6))
  9260 format (/,' jfr                                 :', (10i6))
  9270 format (/,' npne                                :', (10i6))
  9280 format (/,' title                               :',5a16)
  9290 format (/,' element',10x,'origin')
  9300 format (1x,2a16)
  9310 format (/,1x,'comp                              : ',/,(1p,10e12.5))
  9320 format (/,1x,'atemp                             : ',1p,e12.5)
  9330 format (/,1x,'b                                 :',/,(5x,1p,10e12.5))
  9340 format (/,' element',10x,'temperatures')
  9350 format (1x,i4,10x,10f12.3)
  9360 format (/,' spectrum                            :',i6,/)
  9370 format (1x,1p,10e12.5)
  9380 format (/,1x,'fe                                :',1p,e12.5)
  9390 format (/,1x,'last record of general contents   :',i6)
  9400 format (/,1x,'general contents record block structure')
  9410 format (2x,'f rec',4x,'data',4x,'start',4x,' words ',/,2x,' no. ',4x, &
  & 'type',3x,'record',3x,' stored ')
  9420 format (2x,i5,3i9)
  9430 format (/,' ****** warning - ',' fission spectra does not sum to 1.0', &
  '  value is ',f10.4,' ******',/, &
  '  fission spectra will be normalised to sum to 1.0')
  9440 format (/,' sum of fission spectrum = ',f10.5)
  9450 format (' input group boundary is not a fine group boundary',/, &
  & ' group ',i5,' energies ',1p,2e12.5)
  9460 format (' ****** warning - no file 6 (matrix) data found on gendf* tape')
  9470 format (/,' ***** error - temperature dependant thermal', &
  ' scattering matrices have been requested but are', &
  ' not available on the gendf* *****')
  9480 format (' ****** error - no file 3 (cross section) data found', &
  & //' on gendf* tape')
  9490 format (/,' ***** error -',' nu bar not on input gendf* tape so nu*' &
  & //'fission ','can not be formed *****')
  9500 format (/,' ***** error - array dimension exceeded. the ', &
  'parameter maxthr should be increased from ',i6,' to ',i6)
  9510 format (/,' thermal upscatter starts in group          ',i6,/, &
  & ' lowest group to which scattering occurs is ',i6)
  end subroutine ecncon
  !
  subroutine ecnadd(irect,nrect,lenrec,eng,nfis,ifis,ntfis,sysour,nreci, &
  & ing,matrix,matrx2,porder,pordr2,maxm,maxp,lunit,lwordr,thrmin,therma, &
  & torder,tnmin,tnmax,thdim1)
  !
  ! add element to the general contents package
  !
  ! arguments
  !  irect - record number on input general constants file
  !  nrect - record number on output general constants file
  !  lenrec - length of physical records for each data type
  !  eng - broad group energy boundaries
  !
  ! general contents data block structure:
  ! record block 1: integers; (identifier)
  ! record block 2: integers;
  ! record block 3: integers;
  ! record block 4: characters
  ! record block 5: reals
  ! record block 5+ifis: reals
  ! last record block 5+nfis: reals
  ! .. scalar arguments ..
  integer ifis,irect,lwordr,maxm,maxp,nfis,nreci,nrect,ntfis,thdim1,thrmin
  character sysour*16
  ! .. array arguments ..
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  real(kr) eng(ngrmax+1),therma(maxthr,maxthr,tmpmax),torder(maxthr,maxthr, &
  maxleg,tmpmax)
  integer ing(nreci),lenrec(3),lunit(maxmtx),tnmax(tmpmax,maxthr), &
  tnmin(tmpmax,maxthr)
  ! .. local scalars ..
  real(kr) comp,emax,emin,fe,totspc,za
  integer i,ielp,ielr,if,itempj,j,l,lbuff,lielp,lielr,lifis,litemp,lnfre, &
  lngg,lnpmc,lnsmc,lnsmx,lnspge,lntemp,lntfis,loc,mbuff,msmc,msmtx,n,nb, &
  nbl,ncards,ndil,nfis1,nfrecs,ngg,npnfp,nreg,nsmxj,nspge,ntempj,ntw,nw, &
  nword2,mg,mbuf,mat
  logical fision,found
  data    fision/.false./,found/.false./
  ! .. local arrays ..
  real(kr) buff(maxnbf),btemp(10),flux(ngrmax),sdil(10)
  integer nbuff(maxnbf),ielp1(maxel),ielr1(maxel),ifis1(maxel),lrec(3,maxel), &
  nfre1(maxel),ngmax1(ngrmax),nmax(ngrmax+1),npmc1(maxel),nsmc1(maxel), &
  nspge1(maxel),ntfis1(maxel)
  character(len=16) cbuff(maxcbf)
  ! .. intrinsic functions ..
  intrinsic abs
  !
  ! initialise local arrays and variables
  mbuf = maxnbf
  comp = 0.0
  emax = 0.0
  emin = 0.0
  fe = 0.0
  totspc = 0.0
  za = 0.0
  i = 0
  ielp = 0
  ielr = 0
  if = 0
  itempj = 0
  j = 0
  l = 0
  lbuff = 0
  lielp = 0
  lielr = 0
  lifis = 0
  litemp = 0
  lnfre = 0
  lngg = 0
  lnpmc = 0
  lnsmc = 0
  lnsmx = 0
  lnspge = 0
  lntemp = 0
  lntfis = 0
  loc = 0
  mbuff = 0
  msmc = 0
  msmtx = 0
  n = 0
  nb = 0
  nbl = 0
  ncards = 0
  ndil = 0
  nfis1 = 0
  nfrecs = 0
  ngg = 0
  npnfp = 0
  nreg = 0
  nsmxj = 0
  nspge = 0
  ntempj = 0
  ntw = 0
  nw = 0
  nword2 = 0
  btemp(:10)=0.0d0
  flux(:ngrmax)=0.0d0
  sdil(:10)=0.0d0
  ielp1(:maxel)=0
  ielr1(:maxel)=0
  ifis1(:maxel)=0
  lrec(:3,:maxel)=0
  nfre1(:maxel)=0
  ngmax1(:ngrmax)=0
  nmax(:ngrmax+1)=0
  npmc1(:maxel)=0
  nsmc1(:maxel)=0
  nspge1(:maxel)=0
  ntfis1(:maxel)=0
  ! zero array spec
  spec(:ngi)=0.0d0
  !----
  !  read input general contents package
  !----
  nrect = 0
  irect = 0
  10 call ecrint(nine,lenrec(1),nbuff,lenrec(1),irect,1)
  if(nbuff(1)==0) then
    ! end of input ecco library
    ! general contents not found
    write(nsyso,9440) nine
    call error('ecnadd','general contents data not found',' ')
  else if(nbuff(1)/=2) then
    ! read next package
    go to 10
  endif
  ! number of fortran records in this package
  nfrecs = nbuff(6)
  ! library origin - 2 = jef
  if(nbuff(7)/=2) then
    write(nsyso,9450) nbuff(7)
    call error('ecnadd','wrong library origin',' ')
  endif
  ! print out data block structure
  if(iprint>0) write(nsyso,9380)
  if(iprint>0) write(nsyso,9390)
  lbuff = 9
  do j = 1,nfrecs
    lrec(:3,j) = nbuff(lbuff+1:lbuff+3)
    lbuff = lbuff + 3
    if(iprint>0) write(nsyso,9400) j,lrec(1,j),lrec(2,j),lrec(3,j)
  enddo
  loc = lrec(3,1)
  ! read record 1 - block 2
  irect = lrec(2,2) - 1
  if(iprint>0) write(nsyso,'(//,6(a,i8),//)') ' nine: ', nine,' lenrec(1): ', &
  & lenrec(1),' nbuff(loc): ',nbuff(loc),' lrec(3,2): ', lrec(3,2),' irect: ', &
  & irect
  call ecrint(nine,lenrec(1),nbuff(loc),lrec(3,2),irect,loc)
  ! nreg: number of cell regions
  ! nelp: number of elements in this group structure
  ! ncards: number of cards read as input
  ! ngi: number of groups in this structure
  nreg = nbuff(loc)
  ! if nreg = 1 then gecco library being added to
  nelp = nbuff(loc+1)
  ncards = nbuff(loc+2)
  ngi = nbuff(loc+3)
  ! nmax(ngi): fine group number with same upper boundary
  nmax(:ngi) = nbuff(loc+4:loc+3+ngi)
  ! ielr1: identifier reference array by region
  ! ielp1: identifier reference array by position
  lielr = loc + 3 + ngi
  lielp = lielr + nelp
  lntemp = lielp + nelp
  litemp = lntemp + nelp
  lifis = litemp + nelp + 1
  lntfis = lifis + nelp
  lnsmc = lntfis + nelp
  lnpmc = lnsmc + nelp
  lnsmx = lnpmc + nelp
  lnfre = lnsmx + nelp + 2
  do i = 1,nelp
    ielr1(i) = nbuff(lielr+i)
    ielp1(i) = nbuff(lielp+i)
    ntemp(i) = nbuff(lntemp+i)
    itemp(i) = nbuff(litemp+i)
    ifis1(i) = nbuff(lifis+i)
    ntfis1(i) = nbuff(lntfis+i)
    nsmc1(i) = nbuff(lnsmc+i)
    npmc1(i) = nbuff(lnpmc+i)
    nsmx(i) = nbuff(lnsmx+i)
    nfre1(i) = nbuff(lnfre+i)
  enddo
  nfis1 = nbuff(litemp+nelp+1)
  msmc = nbuff(lnsmx+nelp+1)
  msmtx = nbuff(lnsmx+nelp+2)
  if(iprint>0) write(nsyso,'(//,a,/,(2(a,i8),/))') ' ecnadd avant',' npnfp: ', &
  & npnfp,' lnfre: ', lnfre,' nelp: ',  nelp,' nbuff(lnfre+nelp+1): ', &
  & nbuff(lnfre+nelp+1)
  npnfp = nbuff(lnfre+nelp+1)
  if(iprint>0) write(nsyso,'(//,a,/,(2(a,i8),/))') ' ecnadd apres',' npnfp: ', &
  & npnfp, ' lnfre: ', lnfre,' nelp: ',  nelp,' nbuff(lnfre+nelp+1): ', &
  & nbuff(lnfre+nelp+1)
  ngg = nbuff(lnfre+nelp+2)
  if(ngg/=0) then
    lngg = lnfre + nelp + 2
    ngmax1(:ngg) = nbuff(lngg+1:lngg+ngg)
  endif
  lnspge = lnfre + nelp + 2 + ngg
  nspge1(:nelp) = nbuff(lnspge+1:lnspge+nelp)
  ! increase number of elements in package by 1
  nelp = nelp + 1
  itemp(nelp) = 0
  ! read input gendf
  ! read mf=1 on gendf from unit ning
  call repoz(ning)
  call tpidio(ning,0,0,buff,nb,nw)
  20 call contio(ning,0,0,buff,nb,nw)
  if(mfh/=0) then
    ntemp(nelp) = ntemp(nelp) + 1
    ntempj = ntemp(nelp)
    ! file 1 only read for first temperature
    za = c1h
    if(n1h/=-1) then
      write(nsyso,9000)
      call error('ecnadd',' ******  not a gendf tape ******',' ')
    endif
    ntw = n2h
    call ecfil1(sdil,ndil,ntw,eng,btemp(ntempj),ngrmax)
    if(iprint>0) write(nsyso,9010) btemp(ntempj)
    if(ntempj==1) then
      emin = eng(1)
      emax = eng(ngi+1)
      if(iprint>0) write(nsyso,9020) ning,za,math,ndil,ngi,emin,emax
      if(iprint>0) write(nsyso,9030) (sdil(l),l=1,ndil)
      !  calculate array nmax for the general contents file
      !  nmax(l) - fine group number with the same upper neutron
      !            energy boundary as the lth group of the present
      !            group structure
      ! ecco group structure is stored in mev so convert gendf structure
      eng(:ngi+1) = eng(:ngi+1)/1.d+6
      i  = ngi + 1
      j  = 1
      do l = 1,nfg
        if(abs((egrid(l)-eng(i))/egrid(l))<0.0005)then
          if(l/=nmax(j) .and. j<=ngi) then
            ! different group structure
            write(nsyso,9430) l,nmax(j)
            call error('ecnadd','group structures differ',' ')
          endif
          i  = i - 1
          j  = j + 1
        endif
      enddo
    endif
    ! read the gendf to set constants
    ! read cross-section data
    thdim1 = ngi - maxthr + 1
    if(thdim1<=0) thdim1 = 1
    mat = matno
    call find(mat,3,0,ning,found)
    if(found) then
      call ecfil3(flux,eng)
    else
      write(nsyso,9480)
      call error('ecnadd','** &
      & warning - no file 3 (cross section)data found on gendf* tape',' ')
    endif
    itempj = itemp(nelp)
    if(ntempj==1) then
      ! read delayed chi from file 5 if required
      mat = matno
      call find(mat,5,0,ning,found)
      if(found) call ecfil5
      ! read matrices from file 6
      mat = matno
      call find(mat,6,0,ning,found)
      if(found) then
        call ecfil6(nreci,ing,flux,lenrec,lwordr,matrix,matrx2,porder, &
        & pordr2,maxm,maxp,lunit,itempj,thrmin,fision)
        if((itempj/=0).and.(iprint>0))  write(nsyso,9510) itempj,thrmin
        if(ngi-thrmin>maxthr) then
          write(nsyso,9500) maxthr,ngi - thrmin
          call error('ecnadd','ngi-thrmin>maxthr','')
        endif
      else if(iprint>0) then
        write(nsyso,9460)
      endif
      ! read temperature dependent thermal scattering matrix if present
    else if(itempj>0 .and. ntempj>1) then
      mat = matno
      call find(mat,6,0,ning,found)
      if(found) then
        call ecfil6t(therma,thrmin,tnmin,tnmax,torder,thdim1)
      else
        ! temperature dependant data requested but not present
        write(nsyso,9470)
        call error('ecnadd',' ** temperature dependant thermal scattering mat' &
        & //'rices have been requested but are not available on gendf* **',' ')
      endif
    endif
    ! read sub-group data from ecfil50
    mat = matno
    call find(mat,50,0,ning,found)
    if(found) call ecfil50
    call tomend(ning,0,0,buff)
    call skiprz(ning,-1)
    call contio(ning,0,0,buff,nb,nw)
  endif
  if(math==-1) then
    go to 100
  else
    go to 20
  endif
  ! end of gendf
  ! ntemp: number of temperatures for each element
  !----
  !  assign the general contents parameters (fine group structure)
  !----
  100 nrect = 0
  ! store start record for identifier record block
  lrec(2,1) = nrect + 1
  ! allow 1 record for identifier record block storage
  nrect = nrect + 1
  ! enter data for record block number 2:
  ! nreg: number of cell regions
  ! nelp: number of elements in this group structure
  ! ncards: number of cards read as input
  ! ngi: number of groups in this structure
  if(iprint>0) write(nsyso,9050) nreg,nelp,ncards,ngi
  ! put into buffer array ready for disk write
  nbuff(1) = 1
  nbuff(2) = nelp
  nbuff(3) = ncards
  nbuff(4) = ngi
  ! nmax(ngi): fine group number with same upper boundary
  ! put data into disk buffer
  nbuff(5:4+ngi) = nmax(:ngi)
  ! ielp: identifier reference array by position
  ielp = nelp
  ielr = ielp
  if(iprint>0) then
    write(nsyso,9060) ielr
    write(nsyso,9070) ielp
    write(nsyso,9080) ntempj
  endif
  ! if thermal scattering data has been requested but there is only
  ! one temperature on the input gendf then itemp(nelp) should be
  ! set to 0 on the output library
  if(ntempj==1) then
    itempj = 0
    itemp(nelp) = 0
  endif
  nsmxj = nsmx(nelp)
  ngg = 0
  if(iprint>0) then
    write(nsyso,9090) itempj
    write(nsyso,9100) nfis1 + nfis
    write(nsyso,9110) ifis ! table of fissile isotope indicators
    write(nsyso,9120) ntfis ! number of fission spectra for each element
    write(nsyso,9130) nsmc ! number of x-section types per element
    write(nsyso,9140) npmc ! number of x-section types per element with subgroup
    write(nsyso,9150) nsmxj ! number of x-section matrix types per element
    ! msmc: number of macroscopic x-section types
    ! msmtx: number of macroscopic x-section matrices
    write(nsyso,9160) msmc,msmtx
    write(nsyso,9170) nfre ! number of response function types
    write(nsyso,9180) ngg ! number of gamma groups
    write(nsyso,9190) npnfp ! pn order (+1) to which fluxes are stored
  endif
  ! update buffer array
  lielr = ngi + 4
  lielp = lielr + nelp
  lntemp = lielp + nelp
  litemp = lntemp + nelp
  lifis = litemp + nelp + 1
  lntfis = lifis + nelp
  lnsmc = lntfis + nelp
  lnpmc = lnsmc + nelp
  lnsmx = lnpmc + nelp
  lnfre = lnsmx + nelp + 2
  do i = 1,nelp - 1
    nbuff(lielr+i) = ielr1(i)
    nbuff(lielp+i) = ielp1(i)
    nbuff(lntemp+i) = ntemp(i)
    nbuff(litemp+i) = itemp(i)
    nbuff(lifis+i) = ifis1(i)
    nbuff(lntfis+i) = ntfis1(i)
    nbuff(lnsmc+i) = nsmc1(i)
    nbuff(lnpmc+i) = npmc1(i)
    nbuff(lnsmx+i) = nsmx(i)
    nbuff(lnfre+i) = nfre1(i)
  enddo
  nbuff(lielr+nelp) = ielr
  nbuff(lielp+nelp) = ielp
  nbuff(lntemp+nelp) = ntemp(nelp)
  nbuff(litemp+nelp) = itemp(nelp)
  nbuff(lifis+nelp) = ifis
  nbuff(lntfis+nelp) = ntfis
  nbuff(lnsmc+nelp) = nsmc
  nbuff(lnpmc+nelp) = npmc
  nbuff(lnsmx+nelp) = nsmxj
  nbuff(lnfre+nelp) = nfre
  nbuff(litemp+nelp+1) = nfis1 + nfis
  nbuff(lnsmx+nelp+1) = msmc
  nbuff(lnsmx+nelp+2) = msmtx
  nbuff(lnfre+nelp+1) = npnfp
  nbuff(lnfre+nelp+2) = ngg
  if(ngg/=0) then
    lngg = lnfre + nelp + 2
    nbuff(lngg+1:lngg+ngg) = ngmax1(:ngg)
  endif
  lnspge = lnfre + nelp + 2 + ngg
  nbuff(lnspge+1:lnspge+nelp-1) = nspge1(:nelp-1)
  nspge = 0
  nbuff(lnspge+nelp) = nspge
  ! block:  2 (integers)
  ! number of data words in block
  nword2 = 9 + ngi + ngg + nreg*nelp + 10*nelp
  ! data block 2: (integer)
  lrec(1,2) = 1
  ! data block 2: start record
  lrec(2,2) = nrect + 1
  ! data block 2: number of words
  lrec(3,2) = nword2
  ! write output record
  call ecwint(nout,lenrec(1),nbuff,nword2,nrect,1)
  ! block:  3 (integers)
  ! read input record
  irect = lrec(2,3) - 1
  call ecrint(nine,lenrec(1),nbuff,lrec(3,3),irect,1)
  ! start address in nbuff
  lbuff = lrec(3,3)
  mbuff = 0
  ! ismc(nsmc): reference list number of x-section type per element
  ! update library buffer
  do j = 1,nelp - 1
    if(nsmc1(j)>0) then
      nbuff(lbuff+1:lbuff+nsmc1(j)) = nbuff(mbuff+1:mbuff+nsmc1(j))
      lbuff = lbuff + nsmc1(j)
      mbuff = mbuff + nsmc1(j)
    endif
  enddo
  if(nsmc>0) then
    if(iprint>0) write(nsyso,9200) (ismc(i),i=1,nsmc)
    nbuff(lbuff+1:lbuff+nsmc) = ismc(:nsmc)
    lbuff = lbuff + nsmc
  endif
  ! ipmc(npmc): reference list number of x-section types with sub-
  !             groups per element
  ! update library buffer
  do j = 1,nelp - 1
    if(npmc1(j)>0) then
      nbuff(lbuff+1:lbuff+npmc1(j)) = nbuff(mbuff+1:mbuff+npmc1(j))
      lbuff = lbuff + npmc1(j)
      mbuff = mbuff + npmc1(j)
    endif
  enddo
  if(npmc>0) then
    if(iprint>0) write(nsyso,9210) (ipmc(i),i=1,npmc)
    nbuff(lbuff+1:lbuff+npmc) = ipmc(:npmc)
    lbuff = lbuff + npmc
  endif
  ! ismtx(nsmxj): reference list number of x-section type per elem.
  ! update disk buffer
  do j = 1,nelp - 1
    if(nsmx(j)>0) then
      nbuff(lbuff+1:lbuff+nsmx(j)) = nbuff(mbuff+1:mbuff+nsmx(j))
      lbuff = lbuff + nsmx(j)
      mbuff = mbuff + nsmx(j)
    endif
  enddo
  if(nsmxj>0) then
    if(iprint>0) write(nsyso,9220) (ismtx(i),i=1,nsmxj)
    nbuff(lbuff+1:lbuff+nsmxj) = ismtx(:nsmxj)
    lbuff = lbuff + nsmxj
  endif
  ! ifre(nfre): reference list number of x-section type per elem.
  ! update disk buffer
  do j = 1,nelp - 1
    if(nfre1(j)>0) then
      nbuff(lbuff+1:lbuff+nfre1(j)) = nbuff(mbuff+1:mbuff+nfre1(j))
      lbuff = lbuff + nfre1(j)
      mbuff = mbuff + nfre1(j)
    endif
  enddo
  if(nfre>0) then
    if(iprint>0) write(nsyso,9230) (ifre(i),i=1,nfre)
    nbuff(lbuff+1:lbuff+nfre) = ifre(:nfre)
    lbuff = lbuff + nfre
  endif
  ! jfr(nfre): value type indicator for response functions
  !  set equal to ifre
  ! update disk buffer
  do j = 1,nelp - 1
    if(nfre1(j)>0) then
      nbuff(lbuff+1:lbuff+nfre1(j)) = nbuff(mbuff+1:mbuff+nfre1(j))
      lbuff = lbuff + nfre1(j)
      mbuff = mbuff + nfre1(j)
    endif
  enddo
  if(nfre>0) then
    if(iprint>0) write(nsyso,9240) (ifre(i),i=1,nfre)
    nbuff(lbuff+1:lbuff+nfre) = ifre(:nfre)
    lbuff = lbuff + nfre
  endif
  ! npne(nsmx): pn orders; per matrix type, per element
  ! update disk buffer
  do j = 1,nelp - 1
    if(nsmx(j)>0) then
      nbuff(lbuff+1:lbuff+nsmx(j)) = nbuff(mbuff+1:mbuff+nsmx(j))
      lbuff = lbuff + nsmx(j)
      mbuff = mbuff + nsmx(j)
    endif
  enddo
  if(nsmxj>0) then
    if(iprint>0) write(nsyso,9250) (npne(i),i=1,nsmxj)
    nbuff(lbuff+1:lbuff+nsmx(nelp)) = npne(:nsmx(nelp))
    lbuff = lbuff + nsmxj
  endif
  ! data block 3: (integer)
  lrec(1,3) = 1
  ! start record:
  lrec(2,3) = nrect + 1
  ! number of words:
  loc = lrec(3,3) + 1
  lbuff = lbuff - loc + 1
  lrec(3,3) = lbuff
  ! write output record
  call ecwint(nout,lenrec(1),nbuff(loc),lbuff,nrect,1)
  ! block:  4 (characters)
  ! read input record
  irect = lrec(2,4) - 1
  call ecrc16(nine,lenrec(2),cbuff,lrec(3,4),irect,1)
  ! read x-section data origin for this element
  if(iprint>0) then
    write(nsyso,9260) (cbuff(i),i=1,5)
    write(nsyso,9270)
    write(nsyso,9280) sycor(ielp),sysour
  endif
  ! transfer to ecco buffer
  lbuff = lrec(3,4)
  cbuff(lbuff+1) = sysour
  lbuff = lbuff + 1
  ! data block 4: (character)
  lrec(1,4) = 2
  ! start record:
  lrec(2,4) = nrect + 1
  ! number of words:
  lrec(3,4) = lbuff
  ! write output record
  call ecwc16(nout,lenrec(2),cbuff,lbuff,nrect,1)
  ! block:  5 (reals)
  ! read input record
  irect = lrec(2,5) - 1
  call ecrr4(nine,lenrec(3),buff,lrec(3,5),irect,1)
  ! start address in nbuff
  lbuff = lrec(3,5)
  ! update disk buffer
  ! note . nreg should be =0 for an input library
  ! however gecco incorrectly set nreg=1
  ! comp: compositions of nuclides in region 1
  comp = 1.0
  if(iprint>0) write(nsyso,9290) comp
  ! update disk buffer
  buff(lbuff+1:lbuff+nelp-1) = buff(:nelp-1)
  mbuff = nelp
  lbuff = lbuff + nelp
  buff(lbuff) = comp
  lbuff = lbuff + 1
  ! atemp: temperature of this region
  ! update disk buffer
  buff(lbuff) = buff(mbuff)
  mbuff = mbuff + 1
  ! update disk buffer
  do i = 1,ngi
    buff(lbuff+i) = buff(mbuff)
    mbuff = mbuff + 1
  enddo
  lbuff = lbuff + ngi
  ! btemp(ntemp): tabulated temperatures for each element
  if(iprint>0) write(nsyso,9310)
  ! number of temperatures for this nuclide
  if(iprint>0) write(nsyso,9320) ielp, (btemp(i),i=1,ntempj)
  ! update disk buffer
  do i = 1,nelp - 1
    do n = 1,ntemp(i)
      buff(lbuff+n) = buff(mbuff)
      mbuff = mbuff + 1
    enddo
    lbuff = lbuff + ntemp(i)
  enddo
  buff(lbuff+1:lbuff+ntempj) = btemp(:ntempj)
  lbuff = lbuff + ntempj
  ! data block 5:
  lrec(1,5) = 3
  ! start record:
  lrec(2,5) = nrect + 1
  ! number of words:
  loc = lrec(3,5) + 1
  lbuff = lbuff - lrec(3,5)
  lrec(3,5) = lbuff
  call ecwr4(nout,lenrec(3),buff(loc),lbuff,nrect,1)
  ! data for ecco blocks: 6 -> 5 + nfis
  ! calculate fission spectrum
  ! transfer each fission spectrum separately
  if(nfis1/=0) then
    ! input library contains fissile elements
    do if = 1,nfis1
      irect = lrec(2,5+if) - 1
      lrec(2,5+if) = nrect + 1
      call ecrr4(nine,lenrec(3),buff,lrec(3,5+if),irect,1)
      call ecwr4(nout,lenrec(3),buff,lrec(3,5+if),nrect,1)
    enddo
  endif
  if(nfis/=0) then
    if(nnufis/=0) then
      if = nfis1 + nfis
      nfis = nfis1 + nfis
      ! find 4ev and set the fission spectrum below it to 1e-20
      do mg = ngi,1,-1
        if(egrid(mg)>4.0001e-6) go to 387
        spec(mg)=1e-20
      enddo
      ! fissile so read in spectra
      387 totspc = 0.0
      do i = 1,ngi
        buff(i) = spec(i)
        totspc = totspc + buff(i)
      enddo
      if(iprint>0) write(nsyso,9420) totspc
      if(abs(totspc-1.0)>0.0001) then
        if(iprint>0) write(nsyso,9410) totspc
        ! normalise fission spectra to sum to 1.0
        buff(:ngi) = buff(:ngi)/totspc
      endif
      if(iprint>0) write(nsyso,9330) if
      if(iprint>0) write(nsyso,9340) (buff(i),i=1,ngi)
      ! ef: fission energy
      fe = 0.0
      if(.not.fision) fe = 1.e+6
      if(iprint>0) write(nsyso,9350) fe
      ! update disk buffer
      buff(ngi+1) = fe
      lbuff = ngi + 1
      ! each fission spectrum is in a new data block
      ! block 5+if: (real)
      lrec(1,5+if) = 3
      ! start record:
      lrec(2,5+if) = nrect + 1
      ! number of words:
      lrec(3,5+if) = lbuff
      ! write output record
      call ecwr4(nout,lenrec(3),buff,lbuff,nrect,1)
    else
      ! nu bar not on input tape
      write(nsyso,9490)
      call error('ecnadd',' nu bar not on input gendf* tape so nu* &
      fission cannot be formed *****',' ')
    endif
  else
    nfis = nfis1
  endif
  !----
  !  finally, create identifier record block
  !----
  ! data type is integer:
  lrec(1,1) = 1
  ! number of words in first record block
  lrec(3,1) = 9 + 3* (5+nfis)
  ! type number of this package
  nbuff(1) = 2
  ! name number of this package
  nbuff(2) = 2
  ! father number of this package
  nbuff(3) = 1
  ! structure number
  nbuff(4) = 1
  ! number of physical records in this package
  nbuff(5) = nrect - lrec(2,1) + 1
  ! number of fortran records in this package
  nfrecs = nfrecs + nfis - nfis1
  nbuff(6) = nfrecs
  ! library origin - 2 = jef
  nbuff(7) = 2
  ! zero unused variables
  nbuff(8) = 0
  nbuff(9) = 0
  lbuff = 9
  do j = 1,nfrecs
    nbuff(lbuff+1:lbuff+3) = lrec(:3,j)
    lbuff = lbuff + 3
  enddo
  ! reset record counter to beginning of package
  nrect = lrec(2,1) - 1
  call ecwint(nout,lenrec(1),nbuff,lbuff,nrect,1)
  ! reset nrect to maximum value
  nrect = lrec(2,1) + nbuff(5) - 1
  ! print out data block structure
  if(iprint>0) then
    write(nsyso,9360) nrect
    write(nsyso,9370)
    write(nsyso,9390)
    do nbl = 1,nfrecs
      write(nsyso,9400) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    enddo
  endif
  return
  !----
  !  formats
  !----
  9000 format (' ****** error - not a gendf tape ******')
  9010 format (/,25x,'temperature ',f8.2,/,25x,20 ('='))
  9020 format (/,' gendf read from unit               :',i7,/, &
  ' material identifier (za)          :',f7.0,/, &
  ' material identifier (mat)         :',i7,/, &
  ' number of dilutions               :',i7,/, &
  ' number of groups                  :',i7,4x,1p,' from ',e10.4,' to ',e10.4, &
  ' ev')
  9030 format (/,' dilutions                         :',1p,10e12.4)
  9040 format (/,' group structure                      :',1p,/,(8e14.6))
  9050 format (/,' nreg,nelp,ncards,ngi                 :',4i6,/)
  9060 format (/,1x,'ielr:', (1x,10i6))
  9070 format (/,' ielp                                 :', (1x,10i6))
  9080 format (/,' number of temperatures - ntemp       :', (1x,10i6))
  9090 format (/,' group number indicating that from this group',' onwards',/ &
  ,'  the elastic scattering matrix is temperature',' dependent - itemp:',i6)
  9100 format (/,' number of fissile elements - nfis   :',i6)
  9110 format (/,' type number of element (0=non-fissile',',1=fissile, &
  2=fertile) - ifis            :', (1x,10i6))
  9120 format (/,' number of fission spectra - ntfis   :', (1x,10i6))
  9130 format (/,' number of cross-section types - nsmc:', (1x,10i6))
  9140 format (/,' number of sub-group types - npmc    :', (1x,10i6))
  9150 format (/,' number of matrix - nsmx             :', (1x,10i6))
  9160 format (/,' msmc,msmtx:',2i6)
  9170 format (/,' number of response functions - nfre :', (1x,10i6))
  9180 format (/,' number of gamma groups - ngg        :',i6)
  9190 format (/,' pn order +1 to which the fluxes are stored',' - npnfp:',i6)
  9200 format (/,' ismc                                :', (10i6))
  9210 format (/,' ipmc                                :', (10i6))
  9220 format (/,' ismtx                               :', (10i6))
  9230 format (/,' ifre                                :', (10i6))
  9240 format (/,' jfr                                 :', (10i6))
  9250 format (/,' npne                                :', (10i6))
  9260 format (/,' title                               :',5a16)
  9270 format (/,' element',10x,'origin')
  9280 format (1x,2a16)
  9290 format (/,1x,'comp                              :',(1x,1p,10e12.5))
  9300 format (/,1x,'b                                 :',/,(5x,1p,10e12.5))
  9310 format (/,' element',10x,'temperatures')
  9320 format (1x,i4,(10x,10f12.3))
  9330 format (/,' normalised spectrum         :',i6,/)
  9340 format (1x,1p,10e12.5)
  9350 format (/,1x,'fe                                :',1p,e12.5)
  9360 format (/,1x,'last record of general contents   :',i6)
  9370 format (/,1x,'output general contents record block structure')
  9380 format (/,1x,'input general contents record block structure')
  9390 format (2x,'f rec',4x,'data',4x,'start',4x,' words ',/,2x,' no. ',4x, &
  'type',3x,'record',3x,' stored ')
  9400 format (2x,i5,3i9)
  9410 format (/,' ****** warning - fission spectra does not sum', &
  ' to 1.0 value is ',f10.4,' *****',/, &
  '  fission spectra will be normalised to sum to 1.0')
  9420 format (/,' sum of fission spectrum = ',f10.5)
  9430 format (/,' ****** error - group structures differ ',2i5,' ******',/)
  9440 format (/,' ****** error - in ecnadd: general contents data', &
  ' not found on unit:',i6,' ******',/)
  9450 format (/,' ****** error - in ecnadd: library origin is not',' jef',i6 &
  ,' ******',/)
  9460 format (' ****** warning - no file 6 (matrix) data found on',' gendf* &
   tape')
  9470 format (/,' ***** error - temperature dependant thermal ', &
  ' scattering matrices have been requested but are not', &
  ' available on the gendf* *****')
  9480 format (' ****** warning - no file 3 (cross section) data found', &
  ' on gendf* tape')
  9490 format (/,' ***** error ',' nu bar not on input gendf* tape so nu* &
  fission','cannot be formed *****')
  9500 format (/,' ***** error - array dimension exceeded. the ', &
  'parameter maxthr should be increased from ',i6,' to ',i6)
  9510 format (/,' thermal upscatter starts in group          ',i6,/, &
  ' lowest group to which scattering occurs is ',i6)
  end subroutine ecnadd
  !
  subroutine ecosad(irect,nrect,lenrec,nreci,ing,matrix,mdim,porder,pdim, &
  & pordr2,pordr3,lunit,maxblk,lblock,mstart,istart,mblock,pblock,sblock, &
  & pstart,totinm,thrmin,therma,torder,tnmin,tnmax,thdim1)
  ! add element to the ecco broad group cross section package
  ! arguments
  ! irect - record number on input general constants file
  ! nrect - record number on output general constants file
  ! lenrec - length of physical records for each data type
  ! eng - broad group energy boundaries
  ! maxblk - maximum number of groups in group block
  ! .. scalar arguments ..
  integer irect,maxblk,mdim,nreci,nrect,pdim,thdim1,thrmin
  ! .. array arguments ..
  real(kr) matrix(mdim),porder(pdim),pordr2(pdim),pordr3(pdim)
  real(kr) therma(thdim1:ngi,thdim1:ngi,tmpmax),torder(thdim1:ngi,thdim1:ngi, &
  maxleg,tmpmax),totinm(ngrmax,maxmtx)
  integer ing(nreci),istart(maxblk+1),lblock(maxblk),lenrec(3),lunit(maxmtx), &
  & mblock(maxblk),mstart(maxblk+1),pblock(maxblk),pstart(maxblk), &
  & sblock(maxblk),tnmax(tmpmax,thdim1:ngi),tnmin(tmpmax,thdim1:ngi)
  ! .. local scalars ..
  real(kr) buff(maxnbf),delta,diff,elterm,interm,sumelg,suming,sumnxn,totnxn, &
  & xnterm
  integer i,ig,ighi,iglo,im,in,ingmax,ingmin,ingp,inrec,ip,ir,iseg,istp, &
  istrt,it,itempj,j,jpos,k,lbuff,loc,m,mbuff,mmg,mpos,mstop,mxsize,n,nbl, &
  nfrecs,ng,nm,nor,npne1,npos,nrecs,nsmxj,nt,ntempj,numrec,nwords,outrec,pos, &
  pp,ppos,ppos2,ppos3,tot1,tot2,totmat,totord,totrf,totsub,totxs,xpos
  logical ifail
  ! .. local arrays ..
  integer nbuff(maxnbf),lrec(3,9000),morder(6),ntot(ngblk),orec(3,9000)
  ! .. intrinsic functions ..
  intrinsic abs
  !
  ! zero local arrays and set defaults
  totinm(:ngrmax,:maxmtx)=0.0d0
  lrec(:3,9000)=0
  orec(3,:9000)=0
  morder(:6)=0
  ntot(:ngblk)=0
  delta = 0.0
  diff = 0.0
  elterm = 0.0
  interm = 0.0
  sumelg = 0.0
  suming = 0.0
  sumnxn = 0.0
  totnxn = 0.0
  xnterm = 0.0
  ig = 0
  ighi = 0
  iglo = 0
  im = 0
  in = 0
  ingmax = 0
  ingmin = 0
  ingp = 0
  inrec = 0
  ip = 0
  ir = 0
  iseg = 0
  istp = 0
  istrt = 0
  it = 0
  itempj = 0
  j = 0
  jpos = 0
  k = 0
  lbuff = 0
  loc = 0
  m = 0
  mbuff = 0
  mmg = 0
  mpos = 0
  mstop = 0
  mxsize = 0
  n = 0
  nbl = 0
  nfrecs = 0
  ng = 0
  nm = 0
  nor = 0
  npne1 = 0
  npos = 0
  nrecs = 0
  nsmxj = 0
  nt = 0
  ntempj = 0
  numrec = 0
  nwords = 0
  outrec = 0
  pos = 0
  pp = 0
  ppos = 0
  ppos2 = 0
  ppos3 = 0
  tot1 = 0
  tot2 = 0
  totmat = 0
  totord = 0
  totrf = 0
  totsub = 0
  totxs = 0
  xpos = 0
  do m = 1,nsmx(nelp)
    tot1 = 0
    tot2 = 0
    if(iprint>0) write(nsyso,9000) listmx(ismtx(m))
    do n = 1,nreci
      if(iprint>0) write(nsyso,9010) n,lnrec1(n,1,m),lnrec1(n,2,m)
      tot1 = tot1 + lnrec1(n,2,m)
      if(iprint>0) write(nsyso,9020) lnrec2(n,1,m),lnrec2(n,2,m)
      tot2 = tot2 + lnrec2(n,2,m)
    enddo
    if(iprint>0) write(nsyso,9030) tot1,tot2
  enddo
  ! read input cross-section record
  irect = 0
  30 call ecrint(nine,lenrec(1),nbuff,lenrec(1),irect,1)
  ! type number of this package
  if(nbuff(1)==0 .and. nbuff(2)==0) then
    ! end of library
    write(nsyso,9190)
    call error('ecosad',' ****** error - end of input ecco library reached c' &
    & //'ross-section record expected but not found ******',' ')
  else if(nbuff(1)/=3) then
    ! not cross-section record
    ! move to next package
    irect = irect + nbuff(5) - 1
    go to 30
  endif
  if(nbuff(12)>lenrec(1)) then
    ! first fortran record is greater 1 physical record
    ! reread last record
    irect = irect - 1
    nwords = nbuff(12)
    call ecrint(nine,lenrec(1),nbuff,nwords,irect,1)
  endif
  ! number of fortran records in this package
  nfrecs = nbuff(6)
  ! library origin - 2 = jef
  if(nbuff(7)/=2) then
    ! not jef library origin
    call error('ecosad',' ****** error - library origin not jef ******',' ')
  endif
  ! print out data block structure
  if(iprint>0) write(nsyso,9130)
  if(iprint>0) write(nsyso,9150)
  lbuff = 9
  do j = 1,nfrecs
    lrec(1,j) = nbuff(lbuff+1)
    lrec(2,j) = nbuff(lbuff+2)
    lrec(3,j) = nbuff(lbuff+3)
    if(iprint>0) write(nsyso,9160) j,lrec(1,j),lrec(2,j),lrec(3,j)
    lbuff = lbuff + 3
  enddo
  !  identifier block: (number 1)
  !  block type: (integer)
  !  store start record for identifier record block (block 1)
  orec(2,1) = nrect + 1
  orec(1,1) = 1
  !  number of words  = number of words on input library +6 ( lrec values for
  !                     records b and c for new element )
  orec(3,1) = lrec(3,1) + 6*nreci
  ! number of records which will be occupied by the identifier record
  nrecs = (orec(3,1)+ (lenrec(1)-1))/lenrec(1)
  nrect = nrect + nrecs
  !  transfer data to ecco library
  !  data block 2:
  !  data type: (integer)
  orec(1,2) = 1
  orec(2,2) = nrect + 1
  !  copy record 1 from input library
  call ecrint(nine,lenrec(1),nbuff,lrec(3,2),irect,1)
  orec(3,2) = lrec(3,2)
  call ecwint(nout,lenrec(1),nbuff,orec(3,2),nrect,1)
  if(nreci/=nbuff(1)) then
    ! different blocking structure on input library than requested in input data
    if(iprint>0) write(nsyso,9120) nbuff(1),nreci
    call error('ecosad','different blocking structure',' ')
  endif
  outrec = 3
  inrec = 3
  ! loop over group blocks
  ingmax = 0
  do n = 1,nreci
    ifail = .false.
    ingmin = ingmax + 1
    ingmax = ingmin + ing(n) - 1
    ! read next block record a from input ecco library
    ! integer data
    call ecrint(nine,lenrec(1),nbuff,lrec(3,inrec),irect,1)
    lbuff = lrec(3,inrec) + 1
    npos = 1
    ! put ishld/inmin/inmax data in buffer and write to disk
    ! store the current value of nrect
    ! block: 3+(n-1)*3
    ! data type: integer
    ! start record
    orec(2,outrec) = nrect + 1
    nsmxj = nsmx(nelp)
    ntempj = ntemp(nelp)
    itempj = itemp(nelp)
    ! update buffer array with  ishld values
    ! copy ishld data from input library
    do ig = ingmin,ingmax
      do i = 1,nelp - 1
        do nt = 1,ntemp(i)
          nbuff(lbuff) = nbuff(npos)
          lbuff = lbuff + 1
          npos = npos + 1
        enddo
      enddo
      ! add data for new element
      do nt = 1,ntempj
        nbuff(lbuff) = ishld(nt,ig)
        lbuff = lbuff + 1
      enddo
    enddo
    ! update buffer array with inmin values
    ! copy inmin data from input library
    do ig = ingmin,ingmax
      do i = 1,nelp - 1
        do nm = 1,nsmx(i)
          nbuff(lbuff) = nbuff(npos)
          lbuff = lbuff + 1
          npos = npos + 1
        enddo
      enddo
      ! add inmin data for new element
      do k = 1,nsmxj
        nbuff(lbuff) = inmin(k,ig)
        lbuff = lbuff + 1
      enddo
    enddo
    ! update buffer array with inmax values
    ! copy inmax data from input library
    do ig = ingmin,ingmax
      do i = 1,nelp - 1
        do nm = 1,nsmx(i)
          nbuff(lbuff) = nbuff(npos)
          lbuff = lbuff + 1
          npos = npos + 1
        enddo
      enddo
      ! add inmax data for new element
      do k = 1,nsmxj
        nbuff(lbuff) = inmax(k,ig)
        lbuff = lbuff + 1
      enddo
    enddo
    ! temperature dependent thermal elastic matrix
    ! copy tnmin data from input library
    do ig = ingmin,ingmax
      do i = 1,nelp - 1
        if((itemp(i)>0) .and. (ig>=itemp(i))) then
          do nt = 2,ntemp(i)
            nbuff(lbuff) = nbuff(npos)
            lbuff = lbuff + 1
            npos = npos + 1
          enddo
        endif
      enddo
      ! add tnmin data for new element*
      if((itempj>0) .and. (ig>=itempj)) then
        do nt = 2,ntempj
          nbuff(lbuff) = tnmin(nt,ig)
          lbuff = lbuff + 1
        enddo
      endif
    enddo
    ! copy tnmax data from input library
    do ig = ingmin,ingmax
      do i = 1,nelp - 1
        if((itemp(i)>0) .and. (ig>=itemp(i))) then
          do nt = 2,ntemp(i)
            nbuff(lbuff) = nbuff(npos)
            lbuff = lbuff + 1
            npos = npos + 1
          enddo
        endif
      enddo
      ! add tnmax data for new element
      if((itempj>0) .and. (ig>=itempj)) then
        do nt = 2,ntempj
          nbuff(lbuff) = tnmax(nt,ig)
          lbuff = lbuff + 1
        enddo
      endif
    enddo
    ! number of data words in this block
    loc = lrec(3,inrec) + 1
    lbuff = lbuff - loc
    orec(3,outrec) = lbuff
    orec(1,outrec) = lrec(1,inrec)
    nrecs = (orec(3,outrec)+ (lenrec(1)-1))/lenrec(1)
    call ecwint(nout,lenrec(1),nbuff(loc),lbuff,nrect,loc)
    ! copy cross-sections records b and c from input ecco library
    inrec = inrec + 1
    outrec = outrec + 1
    numrec = (nelp-1)*2
    do j = 1,numrec
      ! copy cross-section records from input library to output library
      ! real data
      irect = lrec(2,inrec) - 1
      call ecrr4(nine,lenrec(3),buff,lrec(3,inrec),irect,1)
      orec(1,outrec) = lrec(1,inrec)
      orec(2,outrec) = nrect + 1
      orec(3,outrec) = lrec(3,inrec)
      call ecwr4(nout,lenrec(3),buff,orec(3,outrec),nrect,1)
      inrec = inrec + 1
      outrec = outrec + 1
    enddo
    ! calculate size of each group block
    ! mblock(ig) = storage required for matrices
    ! sblock(ig) = storage required for sub-group data
    ! lbuff = length of cross=section record
    ! istart(ig) = starting position of each group in cross-section record
    lbuff = 0
    morder(:nsmxj) = 0
    ! sub group data
    ng = 0
    do ig = ingmin,ingmax
      ng = ng + 1
      sblock(ng) = 0
      do it = 1,ntempj
        sblock(ng) = sblock(ng) + (npmc+1)*ishld(it,ig)
      enddo
      ! matrices
      mblock(ng) = 0
      pblock(ng) = 0
      do i = 1,nsmxj
        iglo = inmin(i,ig)
        ighi = inmax(i,ig)
        if((iglo/=0) .and. (ighi/=0)) then
          ingp = ighi - iglo + 1
          mblock(ng) = mblock(ng) + ingp
          morder(i) = morder(i) + ingp* (npne(i)-1)
          pblock(ng) = pblock(ng) + ingp* (npne(i)-1)
        endif
      enddo
      if((itempj>0) .and. (ig>=itempj)) then
        ! temperature dependent thermal matrices
        do nt = 2,ntempj
          iglo = tnmin(nt,ig)
          ighi = tnmax(nt,ig)
          if((iglo/=0) .and. (ighi/=0)) then
            ingp = ighi - iglo + 1
            mblock(ng) = mblock(ng) + ingp
            morder(melas) = morder(melas) +ingp* (npne(melas)-1)
            pblock(ng) = pblock(ng) +ingp* (npne(melas)-1)
          endif
        enddo
      endif
      totsub = totsub + sblock(ng)
      totmat = totmat + mblock(ng)
      totord = totord + pblock(ng)
      lblock(ng) = mblock(ng) + nsmc*ntempj + sblock(ng) + nfre
      lbuff = lbuff + lblock(ng)
    enddo
    mbuff = lbuff + 1
    istart(1) = 1
    pstart(1) = mbuff
    do ng = 1,ing(n) - 1
      istart(ng+1) = istart(ng) + lblock(ng)
      pstart(ng+1) = pstart(ng) + pblock(ng)
    enddo
    ! set up space for higher order scattering matrices
    mxsize = 0
    do im = 1,nsmxj
      mxsize = mxsize + morder(im)
    enddo
    mstop = mbuff + mxsize - 1
    if(iprint>0) write(nsyso,9040) maxnbf,mstop
    if(mstop>maxnbf) then
      write(nsyso,9050) maxnbf,mstop
      call error('ecosad','buffer array exceeded',' ')
    endif
    ! zeroise buffer
    buff(:mstop)=0.0
    ! check matrix terms
    call ectsum(totinm,nreci,ing,ismtx,listmx,nsmxj,ntempj,itempj,tnmin,tnmax, &
    therma,lenrec,matrix,porder,lnrec1,lnrec2,npne,lunit,inmin,inmax,mdim, &
    pdim,thdim1,ngi,ingmin,ingmax,n)
    ! read in matrix terms
    mpos = mbuff
    ng = 0
    if(minel>0) then
      call buffin(minel,n,lenrec(3),xpos,ppos,matrix,mdim,porder,pdim,lnrec1, &
      lnrec2,npne(minel),lunit(minel))
      ppos = ppos + 1
    endif
    if(melas>0) then
      call buffin(melas,n,lenrec(3),xpos,ppos2,matrix,mdim,pordr2,pdim,lnrec1 &
      ,lnrec2,npne(melas),lunit(melas))
      ppos2 = ppos2 + 1
    endif
    if(mnxn>0) then
      call buffin(mnxn,n,lenrec(3),xpos,ppos3,matrix,mdim,pordr3,pdim,lnrec1, &
      lnrec2,npne(mnxn),lunit(mnxn))
      ppos3 = ppos3 + 1
    endif
    ! loop over groups
    do ig = ingmin,ingmax
      ng = ng + 1
      jpos = istart(ng)
      ! sum inelastic order 1 over g'
      suming = 0.0
      if(minel>0) then
        if(npne(minel)>1) then
          iglo = inmin(minel,ig)
          ighi = inmax(minel,ig)
          if(iglo/=0 .and. ighi/=0) then
            do j = iglo,ighi
              ppos = ppos - (npne(minel)-1)
              suming = suming + porder(ppos)
            enddo
          endif
        endif
      endif
      ! sum nxn order 1 over g'
      sumnxn = 0.0
      if(mnxn>0) then
        if(npne(mnxn)>1) then
          iglo = inmin(mnxn,ig)
          ighi = inmax(mnxn,ig)
          if(iglo/=0 .and. ighi/=0) then
            do j = iglo,ighi
              ppos3 = ppos3 - (npne(mnxn)-1)
              sumnxn = sumnxn + pordr3(ppos3)
            enddo
          endif
        endif
      endif
      ! loop over temperature
      sumelg = 0.0
      do nt = 1,ntempj
        ! calculate (sum elastic order 1 over g')/elastic  = mubar
        if(melas>0) then
          if(npne(melas)>1) then
            if((itempj>0) .and. (ig>=itempj) .and.nt>1) then
              ! temperature dependent matrices
              iglo = tnmin(nt,ig)
              ighi = tnmax(nt,ig)
              if((iglo/=0) .and. (ighi/=0)) then
                sumelg = 0.0
                do j = iglo,ighi
                  sumelg = sumelg +torder(ig,j,1,nt)
                enddo
                sumelg = sumelg/xsec(nelas,nt,ig)
              endif
            else if(nt==1) then
              iglo = inmin(melas,ig)
              ighi = inmax(melas,ig)
              if(iglo/=0 .and. ighi/=0) then
                do j = iglo,ighi
                  ppos2 = ppos2 - (npne(melas)-1)
                  sumelg = sumelg + pordr2(ppos2)
                enddo
                sumelg = sumelg/xsec(nelas,1,ig)
              endif
            endif
            ! check (sum of elastic order 1 over g')/elastic  = mubar
            if(mubar(ig)/=0 .and. nt==1 .and.sumelg/=0.0) then
              diff = (sumelg-mubar(ig))/sumelg
              if(diff>0.0001) then
                if(ifail) then
                  if(iprint>0) write(nsyso,9220) ig,sumelg,mubar(ig),diff
                else
                  if(iprint>0) write(nsyso,9210)
                  if(iprint>0) write(nsyso,9220) ig,sumelg,mubar(ig),diff
                  ifail = .true.
                endif
              endif
            endif
          endif
        endif
        ! calculate transport cross-section
        ! elastic term
        if(sumelg/=0) then
          if(nelas>0) then
            elterm = sumelg*xsec(nelas,nt,ig)
          else
            elterm = 0.0
          endif
          ! inelastic term
          interm = 0.0
          if(ninel>0) then
            if(xsec(ninel,1,ig)/=0.0) then
              interm = suming*xsec(ninel,nt,ig)/xsec(ninel,1,ig)
            endif
          endif
          ! nxn term
          xnterm = 0.0
          if(nnxn>0) then
            if(xsec(nnxn,1,ig)/=0.0) then
              xnterm = sumnxn*xsec(nnxn,nt,ig)/xsec(nnxn,1,ig)
            endif
          endif
          ! transport cross-section=total-elastic term - inelastic term
          xsec(ntran,nt,ig) = xsec(ntotal,nt,ig) - elterm -interm - xnterm
        else
          ! no elastic scattering matrix - calculate transport as 1/3 total xsec
          if(iprint>0) write(nsyso,9230)
          xsec(ntran,nt,ig) = xsec(ntotal,nt,ig)/3.0
        endif
        ! cross-sections
        do in = 1,nsmc
          buff(jpos) = xsec(in,nt,ig)
          jpos = jpos + 1
        enddo
        ! sub group data
        ! weights
        nor = ishld(nt,ig)
        do im = 1,nor
          ! weights
          buff(jpos) = weight(im,nt,ig)
          jpos = jpos + 1
          ! total
          buff(jpos) = prosig(1,im,nt,ig)
          jpos = jpos + 1
          ! calculate transport
          elterm = 0.0
          if(selas>0) then
            elterm = sumelg*prosig(selas,im,nt,ig)
          endif
          interm = 0.0
          if(ninel>0) then
            if(xsec(ninel,1,ig)/=0.0) then
              interm = suming/xsec(ninel,1,ig)*prosig(sinel,im,nt,ig)
            endif
          endif
          prosig(npmc,im,nt,ig) = prosig(1,im,nt,ig) -elterm - interm
          ! other reactions
          do in = 2,npmc
            buff(jpos) = prosig(in,im,nt,ig)*weight(im,nt,ig)
            jpos = jpos + 1
          enddo
        enddo
        mstart(ng) = jpos
      enddo
      ! end of group loop
    enddo
    ! matrices
    do i = 1,nsmxj
      ! read in matrix and higher orders from buffer from buffer
      call buffin(i,n,lenrec(3),xpos,ppos,matrix,mdim,porder,pdim,lnrec1, &
      lnrec2,npne(i),lunit(i))
      ng = 0
      do ig = ingmin,ingmax
        ng = ng + 1
        jpos = mstart(ng)
        mpos = pstart(ng)
        iglo = inmin(i,ig)
        ighi = inmax(i,ig)
        totinm(ig,i) = 0.0
        if(iglo/=0 .and. ighi/=0) then
          do j = iglo,ighi
            ! sum matrices
            totinm(ig,i) = totinm(ig,i) + matrix(xpos)
            buff(jpos) = matrix(xpos)
            xpos = xpos - 1
            jpos = jpos + 1
          enddo
          ! legendre orders
          if(npne(i)>1) then
            iseg = ighi - iglo + 1
            npne1 = npne(i) - 1
            pos = ppos - npne1
            do ip = 1,npne1
              pp = pos + ip
              do j = iglo,ighi
                buff(mpos) = porder(pp)
                pp = pp - npne1
                mpos = mpos + 1
              enddo
            enddo
            ppos = ppos - npne1*iseg
          endif
        endif
        mstart(ng) = jpos
        pstart(ng) = mpos
      enddo
    enddo
    ng = 0
    do ig = ingmin,ingmax
      ng = ng + 1
      jpos = mstart(ng)
      mpos = pstart(ng)
      ! temperature dependent matrices
      if((itempj>0) .and. (ig>=itempj)) then
        do nt = 2,ntempj
          istrt = tnmin(nt,ig)
          istp = tnmax(nt,ig)
          if((istrt/=0) .and. (istp/=0)) then
            do j = istrt,istp
              buff(jpos) = therma(ig,j,nt)
              jpos = jpos + 1
            enddo
            ! legendre orders
            if(npne(melas)>1) then
              do ip = 1,npne(melas) - 1
                do j = istrt,istp
                  buff(mpos) = torder(ig,j,ip,nt)
                  mpos = mpos + 1
                enddo
              enddo
            endif
          endif
        enddo
      endif
      ! response functions
      do ir = 1,nfre
        buff(jpos) = rf(ir,ig)
        jpos = jpos + 1
      enddo
      ! end of group loop
    enddo
    ! write cross-sections away to library
    nfrecs = nfrecs + 1
    orec(1,outrec) = 3
    orec(2,outrec) = nrect + 1
    orec(3,outrec) = lbuff
    call ecwr4(nout,lenrec(3),buff,lbuff,nrect,1)
    nfrecs = nfrecs + 1
    outrec = outrec + 1
    orec(1,outrec) = 3
    orec(2,outrec) = nrect + 1
    orec(3,outrec) = mstop - mbuff + 1
    nrecs = (orec(3,outrec)+ (lenrec(3)-1))/lenrec(3)
    call ecwr4(nout,lenrec(3),buff(mbuff),orec(3,outrec),nrect,mbuff)
    outrec = outrec + 1
    ! end of group block loop
  enddo
  ! check elastic matrix sums to cross-section
  if(melas/=0) then
    ifail = .false.
    if(iprint>0) write(nsyso,9070)
    do ig = 1,ngi
      if(xsec(nelas,1,ig)/=0.0) then
        delta = abs(xsec(nelas,1,ig)-totinm(ig,melas))
        diff = (delta/xsec(nelas,1,ig))*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9080) ig,totinm(ig,melas),xsec(nelas,1,ig),diff
        ifail = .true.
      endif
      mmg = mmg - 1
    enddo
    if(.not.ifail) then
      if(iprint>0) write(nsyso,9100)
    endif
  endif
  ! check n,xn matrix sums to cross-section
  if(mnxn/=0) then
    ifail = .false.
    if(iprint>0) write(nsyso,9090)
    do ig = 1,ngi
      totnxn = 0.0
      if(rn2nd/=0) totnxn = totnxn + rf(rn2nd,ig)*2.0
      if(rn2n/=0) totnxn = totnxn + rf(rn2n,ig)*2.0
      if(rn3n/=0) totnxn = totnxn + rf(rn3n,ig)*3.0
      if(rn4n/=0) totnxn = totnxn + rf(rn4n,ig)*4.0
      if(rn2np/=0) totnxn = totnxn + rf(rn2np,ig)*2.0
      if(rn3np/=0) totnxn = totnxn + rf(rn3np,ig)*3.0
      if(rn2na/=0) totnxn = totnxn + rf(rn2na,ig)*2.0
      if(rn3na/=0) totnxn = totnxn + rf(rn3na,ig)*3.0
      if(rn2n2a/=0) totnxn = totnxn + rf(rn2n2a,ig)*2.0
      if(rnany/=0) totnxn = totnxn + rf(rnany,ig)*nu5(ig)
      if(totnxn/=0.0) then
        delta = abs(totnxn-totinm(ig,mnxn))
        diff = (delta/totnxn)*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9080) ig,totinm(ig,mnxn),totnxn,diff
        ifail = .true.
      endif
      mmg = mmg - 1
    enddo
    if(.not.ifail) then
      if(iprint>0) write(nsyso,9100)
    endif
  endif
  ! check inelastic matrix sums to cross-section
  if(minel/=0 .and. ninel>0) then
    ifail = .false.
    if(iprint>0) write(nsyso,9110)
    do ig = 1,ngi
      if(xsec(ninel,1,ig)/=0.0) then
        delta = abs(xsec(ninel,1,ig)-totinm(ig,minel))
        diff = (delta/xsec(ninel,1,ig))*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9080) ig,totinm(ig,minel),xsec(ninel,1,ig),diff
        ifail = .true.
      endif
      mmg = mmg - 1
    enddo
    if(.not.ifail) then
      if(iprint>0) write(nsyso,9100)
    endif
  endif
  totxs = nsmc*ngi*ntempj
  totrf = nfre*ngi
  if(iprint>0) write(nsyso,9060) totxs,totrf,totsub,totmat,totord
  ! create first data block
  ! type number of this package
  nbuff(1) = 3
  ! name number of this package
  nbuff(2) = 3
  ! father number of this package
  nbuff(3) = 2
  ! structure number
  nbuff(4) = 1
  ! number of physical records in this package
  nbuff(5) = nrect - orec(2,1) + 1
  ! number of fortran records in this package
  nbuff(6) = nfrecs
  ! library origin - 2 = jef
  nbuff(7) = 2
  ! zero unused variables
  nbuff(8) = 0
  nbuff(9) = 0
  lbuff = 9
  do j = 1,nfrecs
    nbuff(lbuff+1) = orec(1,j)
    nbuff(lbuff+2) = orec(2,j)
    nbuff(lbuff+3) = orec(3,j)
    lbuff = lbuff + 3
  enddo
  nrect = orec(2,1) - 1
  ! records occupied by data block record reference array
  nwords = orec(3,1)
  ! write header block (i.e. number 1)
  call ecwint(nout,lenrec(1),nbuff,nwords,nrect,1)
  ! reset record counter to maximum record reached
  nrect = orec(2,1) + nbuff(5) - 1
  ! print out data block structure
  if(iprint>0) then
    write(nsyso,9140)
    write(nsyso,9150)
    nbl = 1
    write(nsyso,9160) nbl,orec(1,nbl),orec(2,nbl),orec(3,nbl)
    nbl = 2
    write(nsyso,9160) nbl,orec(1,nbl),orec(2,nbl),orec(3,nbl)
  endif
  do i = 1,nreci
    ! zero total word per block scaler (n.b. summed over elements)
    ntot(i) = 0
    ! record reference for ishld/inmin/inmax data
    nbl = nbl + 1
    if(iprint>0) write(nsyso,9160) nbl,orec(1,nbl),orec(2,nbl),orec(3,nbl)
    ! record reference for x-section data
    do ir = 1,nelp
      nbl = nbl + 1
      if(iprint>0) write(nsyso,9170) nbl,orec(1,nbl),orec(2,nbl),orec(3,nbl)
      ntot(i) = ntot(i) + orec(3,nbl)
      ! record reference for p > 0 terms
      nbl = nbl + 1
      if(iprint>0) write(nsyso,9170) nbl,orec(1,nbl),orec(2,nbl),orec(3,nbl)
      ntot(i) = ntot(i) + orec(3,nbl)
    enddo
  enddo
  ! output total words per block
  if(iprint>0) write(nsyso,9180) (ntot(i),i=1,nreci)
  return
  !----
  !  formats
  !----
  9000 format (/,10x,a,/,' block  record   size')
  9010 format (i5,i6,4x,i6)
  9020 format (5x,i6,4x,i6)
  9030 format (' total size of p0 matrix = ',i8,/, &
  ' total size of p1 - pn matrix = ',i8)
  9040 format (/,' buffer array usage - space available is ',i8, &
  ' space required is ',i8,/)
  9050 format (' buffer array exceeded. space available is ',i8, &
  ' space required is ',i8)
  9060 format (/,' space required for cross-section data     : ',i10,/, &
  ' space required for response function data : ',i10,/, &
  ' space required for sub-group data         : ',i10,/, &
  ' space required for matrix p0 data         : ',i10,/, &
  ' space required for matrix p1 - p5 data    : ',i10,/)
  9070 format (/,' total elastic has been formed by summing the', &
  ' matrix terms',/,1x,59 ('='),/,/,' group ','  calculated    gendf* &
       difference(%)',/,' -----   ----------    ------     -------------')
  9080 format (1x,i4,1x,1p,2e13.5,0p,f10.3)
  9090 format (/,' total n,xn has been formed by summing the',' matrix terms' &
  ,/,1x,56 ('='),/,/,' group ','  calculated    gendf*     difference(%)',/, &
  ' -----   ----------    ------     -------------')
  9100 format (' no differences greater than 0.01 % found ')
  9110 format (/,' total inelastic has been formed by summing the', &
  ' matrix terms',/,1x,61 ('='),/,/,' group ','  calculated    gendf* &
       difference(%)',/,' -----   ----------    ------     -------------')
  9120 format (/,' ****** error - different blocking structure on', &
  ' input library (',i5,') than requested in input data (',i5,')',' ******')
  9130 format (/,1x,'input cross-section record block structure')
  9140 format (/,1x,'output cross-section record block structure')
  9150 format (2x,'f rec',4x,'data',4x,'start ',4x,'  words ',/,2x,' no. ',4x &
  ,'type',3x,'record',6x,' stored ')
  9160 format (2x,i5,3i9)
  9170 format (2x,i5,3i9)
  9180 format (/,1x,'ntot:', (1x,10i10))
  9190 format (' ****** error - end of input ecco library reached', &
  ' cross-section record expected but not found ******')
  9210 format (' ****** warning - ', &
  ' difference > 0.0001 found between mubar on the gendf* tape ',/, &
  ' and mubar calculated from the elastic cross-section',' ******',/,/, &
  ' group    calculated    gendf*    difference ',/, &
  ' -----    ----------    ------    ----------')
  9220 format (i5,4x,1p,e12.5,e12.5,0p,f10.5)
  9230 format (/,'***note no elastic scattering matrix available for ', &
  'this nuclide. transport will be calculated as 1/3 total')
  end subroutine ecosad
  !
  subroutine ecosec(nrect,lenrec,nreci,ing,matrix,mdim,porder,pdim,pordr2, &
  & pordr3,lunit,maxblk,lblock,mstart,istart,mblock,pblock,sblock,pstart, &
  & totinm,thrmin,therma,torder,tnmin,tnmax,thdim1)
  !
  ! routine to create the ecco broad group cross section package
  !
  !  arguments
  !  nrect - record number on output general constants file
  !  lenrec - length of physical records for each data type
  !  maxblk - maximum number of groups in group block
  ! .. scalar arguments ..
  integer maxblk,mdim,nreci,nrect,pdim,thdim1,thrmin
  ! .. array arguments ..
  real(kr) matrix(mdim),porder(pdim),pordr2(pdim),pordr3(pdim)
  real(kr) therma(thdim1:ngi,thdim1:ngi,tmpmax),torder(thdim1:ngi,thdim1:ngi, &
  maxleg,tmpmax),totinm(ngrmax,maxmtx)
  integer ing(nreci),istart(maxblk+1),lblock(maxblk),lenrec(3),lunit(maxmtx), &
  & mblock(maxblk),mstart(maxblk+1),pblock(maxblk),pstart(maxblk), &
  & sblock(maxblk),tnmax(tmpmax,thdim1:ngi),tnmin(tmpmax,thdim1:ngi)
  ! .. local scalars ..
  real(kr) delta,diff,elterm,interm,sumelg,suming,sumnxn,totnxn,xnterm
  integer i,ig,ighi,iglo,im,in,ingmax,ingmin,ingp,ip,ir,iseg,istp,istrt,it, &
  itempj,j,jpos,k,lbuff,m,mbuff,mmg,mpos,mstop,mxsize,n,nbl,nfrecs,ng,nor, &
  npne1,nrecs,nsmxj,nt,ntempj,nwords,pos,pp,ppos,ppos2,ppos3,tot1,tot2,totmat &
  ,totord,totrf,totsub,totxs,xpos
  logical ifail
  ! .. local arrays ..
  integer nbuff(maxnbf),lrec(3,9000),morder(6),ntot(ngblk)
  real(kr) buff(maxnbf)
  ! .. intrinsic functions ..
  intrinsic abs
  !----
  !  initialise local scalars
  !----
  delta = 0.0
  diff = 0.0
  elterm = 0.0
  interm = 0.0
  sumelg = 0.0
  suming = 0.0
  sumnxn = 0.0
  totnxn = 0.0
  xnterm = 0.0
  i = 0
  ig = 0
  ighi = 0
  iglo = 0
  im = 0
  in = 0
  ingmax = 0
  ingmin = 0
  ingp = 0
  ip = 0
  ir = 0
  iseg = 0
  istp = 0
  istrt = 0
  it = 0
  itempj = 0
  j = 0
  jpos = 0
  k = 0
  lbuff = 0
  m = 0
  mbuff = 0
  mmg = 0
  mpos = 0
  mstop = 0
  mxsize = 0
  n = 0
  nbl = 0
  nfrecs = 0
  ng = 0
  nor = 0
  npne1 = 0
  nrecs = 0
  nsmxj = 0
  nt = 0
  ntempj = 0
  nwords = 0
  pos = 0
  pp = 0
  ppos = 0
  ppos2 = 0
  ppos3 = 0
  tot1 = 0
  tot2 = 0
  totmat = 0
  totord = 0
  totrf = 0
  totsub = 0
  totxs = 0
  xpos = 0
  !----
  !  zero arrays
  !----
  totinm(:ngrmax,:maxmtx)=0.0d0
  lrec(:3,:9000)=0
  morder(:6)=0
  ntot(:ngblk)=0
  ifail = .false.
  do m = 1,nsmx(nelp)
    tot1 = 0
    tot2 = 0
    if(iprint>0) write(nsyso,9010) listmx(ismtx(m))
    do n = 1,nreci
      if(iprint>0) write(nsyso,9020) n,lnrec1(n,1,m),lnrec1(n,2,m)
      tot1 = tot1 + lnrec1(n,2,m)
      if(iprint>0) write(nsyso,9030) lnrec2(n,1,m),lnrec2(n,2,m), &
      & lnrec1(n,2,m) + lnrec2(n,2,m),tot1 + tot2
      tot2 = tot2 + lnrec2(n,2,m)
    enddo
    if(iprint>0) write(nsyso,9040) tot1,tot2
  enddo
  !----
  !  store start record for identifier record block (block 1)
  !----
  lrec(2,1) = nrect + 1
  if(iprint>0) write(nsyso,9140) nreci
  nbuff(1) = nreci
  if(iprint>0) write(nsyso,9150) (ing(i),i=1,nreci)
  ! update buffer
  nbuff(2:nreci+1) = ing(:nreci)
  ! number of data words in this record block
  lbuff = 1 + nreci
  ! identifier block: (number 1)
  ! block type: (integer)
  lrec(1,1) = 1
  ! number of words
  lrec(3,1) = 9 + 3* (2+ (2+1)*nreci)
  ! number of records which will be occupied by the identifier record
  nrecs = (lrec(3,1)+ (lenrec(1)-1))/lenrec(1)
  ! record reached after identifier record
  nrect = nrect + nrecs
  !----
  !  transfer data to ecco library
  !  data block 2:
  !  data type: (integer)
  !----
  lrec(1,2) = 1
  ! start record
  lrec(2,2) = nrect + 1
  ! number of words
  lrec(3,2) = lbuff
  call ecwint(nout,lenrec(1),nbuff,lbuff,nrect,1)
  !----
  !  loop over group blocks
  !----
  totsub = 0
  totmat = 0
  totord = 0
  ingmax = 0
  do n = 1,nreci
    ingmin = ingmax + 1
    ingmax = ingmin + ing(n) - 1
    if(iprint>0) write(nsyso,9050) n,ingmin,ingmax
    ! block: 3
    ! data type: (integer)
    lrec(1,3+ (n-1)*3) = 1
    lrec(2,3+ (n-1)*3) = nrect + 1
    !----
    !  update buffers
    !----
    !  put ishld/inmin/inmax data in buffer and write to disk
    !  store the current value of nrect
    !  block: 3
    !  start record
    !----
    ntempj = ntemp(nelp)
    itempj = itemp(nelp)
    nsmxj = nsmx(nelp)
    lbuff = 0
    !  update buffer array with ing(ibl) ishld values
    do ig = ingmin,ingmax
      ! number of temperatures for this nuclide
      do nt = 1,ntempj
        nbuff(lbuff+nt) = ishld(nt,ig)
      enddo
      lbuff = lbuff + ntempj
    enddo
    !  update buffer array with ing(ibl) inmin values
    do ig = ingmin,ingmax
      do k = 1,nsmxj
        lbuff = lbuff + 1
        nbuff(lbuff) = inmin(k,ig)
      enddo
    enddo
    !  update buffer array with ing(ibl) inmax values
    do ig = ingmin,ingmax
      do k = 1,nsmxj
        lbuff = lbuff + 1
        nbuff(lbuff) = inmax(k,ig)
      enddo
    enddo
    !----
    !  temperature dependent thermal elastic matrix
    !----
    do ig = ingmin,ingmax
      if((itempj>0) .and. (ig>=itempj)) then
        do nt = 2,ntempj
          lbuff = lbuff + 1
          nbuff(lbuff) = tnmin(nt,ig)
        enddo
      endif
    enddo
    do ig = ingmin,ingmax
      if((itempj>0) .and. (ig>=itempj)) then
        do nt = 2,ntempj
          lbuff = lbuff + 1
          nbuff(lbuff) = tnmax(nt,ig)
        enddo
      endif
    enddo
    ! number of data words in this block
    lrec(3,3+ (n-1)*3) = lbuff
    call ecwint(nout,lenrec(1),nbuff,lbuff,nrect,1)
    !----
    !  calculate size of each group block
    !  mblock(ig) = storage required for matrices
    !  sblock(ig) = storage required for sub-group data
    !  lbuff = length of cross=section record
    !  istart(ig) = starting position of each group in cross-section record
    !----
    lbuff = 0
    morder(:nsmxj) = 0
    !----
    ! sub group data
    !----
    ng = 0
    do ig = ingmin,ingmax
      ng = ng + 1
      sblock(ng) = 0
      do it = 1,ntempj
        sblock(ng) = sblock(ng) + (npmc+1)*ishld(it,ig)
      enddo
      ! matrices
      mblock(ng) = 0
      pblock(ng) = 0
      do i = 1,nsmxj
        iglo = inmin(i,ig)
        ighi = inmax(i,ig)
        if((iglo/=0) .and. (ighi/=0)) then
          ingp = ighi - iglo + 1
          mblock(ng) = mblock(ng) + ingp
          morder(i) = morder(i) + ingp* (npne(i)-1)
          pblock(ng) = pblock(ng) + ingp* (npne(i)-1)
        endif
      enddo
      if((itempj>0) .and. (ig>=itempj)) then
        ! temperature dependent thermal matrices
        do nt = 2,ntempj
          iglo = tnmin(nt,ig)
          ighi = tnmax(nt,ig)
          if((iglo/=0) .and. (ighi/=0)) then
            ingp = ighi - iglo + 1
            mblock(ng) = mblock(ng) + ingp
            morder(melas) = morder(melas) +ingp* (npne(melas)-1)
            pblock(ng) = pblock(ng) +ingp* (npne(melas)-1)
          endif
        enddo
      endif
      totsub = totsub + sblock(ng)
      totmat = totmat + mblock(ng)
      totord = totord + pblock(ng)
      lblock(ng) = mblock(ng) + nsmc*ntempj + sblock(ng) + nfre
      lbuff = lbuff + lblock(ng)
    enddo
    mbuff = lbuff + 1
    istart(1) = 1
    pstart(1) = mbuff
    do ng = 1,ing(n) - 1
      istart(ng+1) = istart(ng) + lblock(ng)
      pstart(ng+1) = pstart(ng) + pblock(ng)
    enddo
    !----
    !  set up space for higher order scattering matrices
    !----
    mxsize = 0
    do im = 1,nsmxj
      mxsize = mxsize + morder(im)
    enddo
    mstop = mbuff + mxsize - 1
    if(iprint>0) write(nsyso,9060) maxnbf,mstop
    if(mstop>maxnbf) then
      write(nsyso,9070) maxnbf,mstop
      call error('ecosec','buffer array exceeded',' ')
    endif
    !  zeroise buffer
    buff(:mstop)=0.0
    !----
    !  check matrix terms
    !----
    call ectsum(totinm,nreci,ing,ismtx,listmx,nsmxj,ntempj,itempj,tnmin,tnmax, &
    therma,lenrec,matrix,porder,lnrec1,lnrec2,npne,lunit,inmin,inmax,mdim, &
    pdim,thdim1,ngi,ingmin,ingmax,n)
    !----
    !  read in matrix terms
    !----
    mpos = mbuff
    ng = 0
    if(minel>0) then
      call buffin(minel,n,lenrec(3),xpos,ppos,matrix,mdim,porder,pdim,lnrec1, &
      lnrec2,npne(minel),lunit(minel))
      ppos = ppos + 1
    endif
    if(melas>0) then
      call buffin(melas,n,lenrec(3),xpos,ppos2,matrix,mdim,pordr2,pdim,lnrec1 &
      ,lnrec2,npne(melas),lunit(melas))
      ppos2 = ppos2 + 1
    endif
    if(mnxn>0) then
      call buffin(mnxn,n,lenrec(3),xpos,ppos3,matrix,mdim,pordr3,pdim,lnrec1, &
      lnrec2,npne(mnxn),lunit(mnxn))
      ppos3 = ppos3 + 1
    endif
    !----
    ! loop over groups
    !----
    do ig = ingmin,ingmax
      ng = ng + 1
      jpos = istart(ng)
      ! sum inelastic order 1 over g'
      suming = 0.0
      if(minel>0) then
        if(npne(minel)>1) then
          iglo = inmin(minel,ig)
          ighi = inmax(minel,ig)
          if(iglo/=0 .and. ighi/=0) then
            do j = iglo,ighi
              ppos = ppos - (npne(minel)-1)
              suming = suming + porder(ppos)
            enddo
          endif
        endif
      endif
      !----
      ! sum nxn order 1 over g'
      !----
      sumnxn = 0.0
      if(mnxn>0) then
        if(npne(mnxn)>1) then
          iglo = inmin(mnxn,ig)
          ighi = inmax(mnxn,ig)
          if(iglo/=0 .and. ighi/=0) then
            do j = iglo,ighi
              ppos3 = ppos3 - (npne(mnxn)-1)
              sumnxn = sumnxn + pordr3(ppos3)
            enddo
          endif
        endif
      endif
      !----
      !  loop over temperature
      !----
      sumelg = 0.0
      do nt = 1,ntempj
        ! calculate (sum elastic order 1 over g')/elastic  = mubar
        if(melas>0) then
          if(npne(melas)>1) then
            if((itempj>0) .and. (ig>=itempj) .and.nt>1) then
              ! temperature dependent matrices
              iglo = tnmin(nt,ig)
              ighi = tnmax(nt,ig)
              if((iglo/=0) .and. (ighi/=0)) then
                sumelg = 0.0
                do j = iglo,ighi
                  sumelg = sumelg +torder(ig,j,1,nt)
                enddo
                sumelg = sumelg/xsec(nelas,nt,ig)
              endif
            else if(nt==1) then
              iglo = inmin(melas,ig)
              ighi = inmax(melas,ig)
              if(iglo/=0 .and. ighi/=0) then
                do j = iglo,ighi
                  ppos2 = ppos2 - (npne(melas)-1)
                  sumelg = sumelg + pordr2(ppos2)
                enddo
                sumelg = sumelg/xsec(nelas,1,ig)
              endif
            endif
            ! check (sum of elastic order 1 over g')/elastic  = mubar
            if(mubar(ig)/=0 .and. nt==1 .and.sumelg/=0.0) then
              diff = (sumelg-mubar(ig))/sumelg
              if(diff>0.0001) then
                if(ifail) then
                  if(iprint>0) write(nsyso,9230) ig,sumelg,mubar(ig),diff
                else
                  if(iprint>0) write(nsyso,9220)
                  if(iprint>0) write(nsyso,9230) ig,sumelg,mubar(ig),diff
                  ifail = .true.
                endif
              endif
            endif
          endif
        endif
        !----
        !  calculate transport cross-section
        !  elastic term
        !----
        if(sumelg/=0) then
          if(nelas>0) then
            elterm = sumelg*xsec(nelas,nt,ig)
          else
            elterm = 0.0
          endif
          ! inelastic term
          interm = 0.0
          if(ninel>0) then
            if(xsec(ninel,1,ig)/=0.0) then
              interm = suming*xsec(ninel,nt,ig)/xsec(ninel,1,ig)
            endif
          endif
          ! nxn term
          xnterm = 0.0
          if(nnxn>0) then
            if(xsec(nnxn,1,ig)/=0.0) then
              xnterm = sumnxn*xsec(nnxn,nt,ig)/xsec(nnxn,1,ig)
            endif
          endif
          ! transport cross-section=total-elastic term - inelastic term
          xsec(ntran,nt,ig) = xsec(ntotal,nt,ig) - elterm -interm - xnterm
        else
          ! no elastic scattering matrix - calculate transport as 1/3 total xsec
          if(iprint>0) write(nsyso,9240)
          xsec(ntran,nt,ig) = xsec(ntotal,nt,ig)/3.0
        endif
        ! cross-sections
        do in = 1,nsmc
          buff(jpos) = xsec(in,nt,ig)
          jpos = jpos + 1
        enddo
        !----
        !  sub group data
        !  weights
        !----
        nor = ishld(nt,ig)
        do im = 1,nor
          ! weights
          buff(jpos) = weight(im,nt,ig)
          jpos = jpos + 1
          ! total
          buff(jpos) = prosig(1,im,nt,ig)
          jpos = jpos + 1
          ! calculate transport
          elterm = 0.0
          if(selas>0) then
            elterm = sumelg*prosig(selas,im,nt,ig)
          endif
          interm = 0.0
          if(ninel>0) then
            if(xsec(ninel,1,ig)/=0.0) then
              interm = suming/xsec(ninel,1,ig)*prosig(sinel,im,nt,ig)
            endif
          endif
          prosig(npmc,im,nt,ig) = prosig(1,im,nt,ig) -elterm - interm
          ! other reactions
          do in = 2,npmc
            buff(jpos) = prosig(in,im,nt,ig)*weight(im,nt,ig)
            jpos = jpos + 1
          enddo
        enddo
        mstart(ng) = jpos
      enddo
      ! end of group loop
    enddo
    !----
    !  matrices
    !----
    do i = 1,nsmxj
      ! read in matrix and higher orders from buffer
      call buffin(i,n,lenrec(3),xpos,ppos,matrix,mdim,porder,pdim,lnrec1, &
      lnrec2,npne(i),lunit(i))
      ng = 0
      do ig = ingmin,ingmax
        ng = ng + 1
        jpos = mstart(ng)
        mpos = pstart(ng)
        iglo = inmin(i,ig)
        ighi = inmax(i,ig)
        totinm(ig,i) = 0.0
        if(iglo/=0 .and. ighi/=0) then
          do j = iglo,ighi
            ! sum matrices
            totinm(ig,i) = totinm(ig,i) + matrix(xpos)
            buff(jpos) = matrix(xpos)
            xpos = xpos - 1
            jpos = jpos + 1
          enddo
          ! legendre orders
          if(npne(i)>1) then
            iseg = ighi - iglo + 1
            npne1 = npne(i) - 1
            pos = ppos - npne1
            do ip = 1,npne1
              pp = pos + ip
              do j = iglo,ighi
                buff(mpos) = porder(pp)
                pp = pp - npne1
                mpos = mpos + 1
              enddo
            enddo
            ppos = ppos - npne1*iseg
          endif
        endif
        mstart(ng) = jpos
        pstart(ng) = mpos
      enddo
    enddo
    ng = 0
    do ig = ingmin,ingmax
      ng = ng + 1
      jpos = mstart(ng)
      mpos = pstart(ng)
      ! temperature dependent matrices
      if((itempj>0) .and. (ig>=itempj)) then
        do nt = 2,ntempj
          istrt = tnmin(nt,ig)
          istp = tnmax(nt,ig)
          if((istrt/=0) .and. (istp/=0)) then
            do j = istrt,istp
              buff(jpos) = therma(ig,j,nt)
              jpos = jpos + 1
            enddo
            ! legendre orders
            if(npne(melas)>1) then
              do ip = 1,npne(melas) - 1
                do j = istrt,istp
                  buff(mpos) = torder(ig,j,ip,nt)
                  mpos = mpos + 1
                enddo
              enddo
            endif
          endif
        enddo
      endif
      ! response functions
      do ir = 1,nfre
        buff(jpos) = rf(ir,ig)
        jpos = jpos + 1
      enddo
      ! end of group loop
    enddo
    !----
    !  write cross-sections away to library
    !----
    lrec(1,4+ (n-1)*3) = 3
    lrec(2,4+ (n-1)*3) = nrect + 1
    lrec(3,4+ (n-1)*3) = lbuff
    call ecwr4(nout,lenrec(3),buff,lbuff,nrect,1)
    lrec(1,5+ (n-1)*3) = 3
    lrec(2,5+ (n-1)*3) = nrect + 1
    lrec(3,5+ (n-1)*3) = mstop - mbuff + 1
    call ecwr4(nout,lenrec(3),buff(mbuff),lrec(3,5+ (n-1)*3),nrect,mbuff)
    ! end of group block loop
  enddo
  !----
  !  check elastic matrix sums to cross-section
  !----
  if(melas/=0) then
    ifail = .false.
    if(iprint>0) write(nsyso,9090)
    do ig = 1,ngi
      if(xsec(nelas,1,ig)/=0.0) then
        delta = abs(xsec(nelas,1,ig)-totinm(ig,melas))
        diff = (delta/xsec(nelas,1,ig))*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9100) ig,totinm(ig,melas),xsec(nelas,1,ig),diff
        ifail = .true.
      endif
      mmg = mmg - 1
    enddo
    if(.not.ifail) then
      if(iprint>0) write(nsyso,9120)
    endif
  endif
  !----
  !  check n,xn matrix sums to cross-section
  !----
  if(mnxn/=0) then
    ifail = .false.
    if(iprint>0) write(nsyso,9110)
    do ig = 1,ngi
      totnxn = 0.0
      if(rn2nd/=0) totnxn = totnxn + rf(rn2nd,ig)*2.0
      if(rn2n/=0) totnxn = totnxn + rf(rn2n,ig)*2.0
      if(rn3n/=0) totnxn = totnxn + rf(rn3n,ig)*3.0
      if(rn4n/=0) totnxn = totnxn + rf(rn4n,ig)*4.0
      if(rn2np/=0) totnxn = totnxn + rf(rn2np,ig)*2.0
      if(rn3np/=0) totnxn = totnxn + rf(rn3np,ig)*3.0
      if(rn2na/=0) totnxn = totnxn + rf(rn2na,ig)*2.0
      if(rn3na/=0) totnxn = totnxn + rf(rn3na,ig)*3.0
      if(rn2n2a/=0) totnxn = totnxn + rf(rn2n2a,ig)*2.0
      if(rnany/=0) totnxn = totnxn + rf(rnany,ig)*nu5(ig)
      if(totnxn/=0.0) then
        delta = abs(totnxn-totinm(ig,mnxn))
        diff = (delta/totnxn)*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9100) ig,totinm(ig,mnxn),totnxn,diff
        ifail = .true.
      endif
      mmg = mmg - 1
    enddo
    if((.not.ifail).and.(iprint>0)) write(nsyso,9120)
  endif
  !----
  !  check inelastic matrix sums to cross-section
  !----
  if(minel/=0 .and. ninel>0) then
    ifail = .false.
    if(iprint>0) write(nsyso,9130)
    do ig = 1,ngi
      if(xsec(ninel,1,ig)/=0.0) then
        delta = abs(xsec(ninel,1,ig)-totinm(ig,minel))
        diff = (delta/xsec(ninel,1,ig))*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9100) ig,totinm(ig,minel),xsec(ninel,1,ig),diff
        ifail = .true.
      endif
      mmg = mmg - 1
    enddo
    if((.not.ifail).and.(iprint>0)) write(nsyso,9120)
  endif
  totxs = nsmc*ngi*ntempj
  totrf = nfre*ngi
  if(iprint>0) write(nsyso,9080) totxs,totrf,totsub,totmat,totord
  !----
  !  create first data block
  !  type number of this package
  !----
  nbuff(1) = 3
  ! name number of this package
  nbuff(2) = 3
  ! father number of this package
  nbuff(3) = 2
  ! structure number
  nbuff(4) = 1
  ! number of physical records in this package
  nbuff(5) = nrect - lrec(2,1) + 1
  ! number of fortran records in this package
  nfrecs = 2 + (2+1)*nreci
  nbuff(6) = nfrecs
  ! library origin - 2 = jef
  nbuff(7) = 2
  ! zero unused variables
  nbuff(8) = 0
  nbuff(9) = 0
  lbuff = 9
  do j = 1,nfrecs
    nbuff(lbuff+1) = lrec(1,j)
    nbuff(lbuff+2) = lrec(2,j)
    nbuff(lbuff+3) = lrec(3,j)
    lbuff = lbuff + 3
  enddo
  nrect = lrec(2,1) - 1
  ! records occupied by data block record reference array
  nwords = lrec(3,1)
  !----
  !   write header block (i.e. number 1)
  !----
  call ecwint(nout,lenrec(1),nbuff,nwords,nrect,1)
  !  reset record counter to maximum record reached
  nrect = lrec(2,1) + nbuff(5) - 1
  ! print out data block structure
  if(iprint>0) then
    write(nsyso,9160)
    write(nsyso,9170)
    nbl = 1
    write(nsyso,9180) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    nbl = 2
    write(nsyso,9180) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
  endif
  ntotal = 0
  do i = 1,nreci
    ! zero total word per block scaler (n.b. summed over elements)
    ntot(i) = 0
    ! record reference for ishld/inmin/inmax data
    nbl = nbl + 1
    if(iprint>0) write(nsyso,9180) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    ! record reference for x-section data
    nbl = nbl + 1
    if(iprint>0) write(nsyso,9190) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    ntot(i) = ntot(i) + lrec(3,nbl)
    ! record reference for p > 0 terms
    nbl = nbl + 1
    if(iprint>0) write(nsyso,9190) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    ntot(i) = ntot(i) + lrec(3,nbl)
    ntotal = ntotal + ntot(i)
  enddo
  !----
  !  output total words per block
  !----
  if(iprint>0) write(nsyso,9200) (ntot(i),i=1,nreci)
  if(iprint>0) write(nsyso,9210) ntotal
  return
  !----
  !  formats
  !----
  9000 format (1x,' temp  group  tnmin tnmax ',/, (4i6))
  9010 format (/,10x,a,/,' block  record   size')
  9020 format (i5,i8,4x,i8)
  9030 format (5x,i8,4x,i8)
  9040 format (' total size of p0 matrix = ',i8,/, &
  ' total size of p1 - pn matrix = ',i8)
  9050 format (' group block ',i5,' from group ',i5,' to ',i5)
  9060 format (/,' buffer array usage - space available is ',i8, &
  ' space required is ',i8,/)
  9070 format (' buffer array exceeded. space available is ',i8, &
  ' space required is ',i8)
  9080 format (/,' space required for cross-section data     : ',i10,/, &
  ' space required for response function data : ',i10,/, &
  ' space required for sub-group data         : ',i10,/, &
  ' space required for matrix p0 data         : ',i10,/, &
  ' space required for matrix p1 - p5 data    : ',i10,/)
  9090 format (/,' total elastic has been formed by summing the', &
  ' matrix terms',/,1x,59 ('='),/,/,' group ','  calculated    gendf*', &
  ' -----   ----------    ------     -------------')
  9100 format (1x,i4,1x,1p,2e13.5,0p,f10.3)
  9110 format (/,' total n,xn has been formed by summing the',' matrix terms' &
  ,/,1x,56 ('='),/,/,' group ','  calculated    gendf*     difference(%)',/, &
  ' -----   ----------    ------     -------------')
  9120 format (' no differences greater than 0.01 % found ')
  9130 format (/,' total inelastic has been formed by summing the', &
  ' matrix terms',/,1x,61 ('='),/,/,' group ','  calculated    gendf* &
       difference(%)',/,' -----   ----------    ------     -------------')
  9140 format (/,1x,'nreci:',i6)
  9150 format (/,1x,'ing:', (1x,21i6))
  9160 format (/,1x,'cross-section record block structure')
  9170 format (2x,'f rec',4x,'data',4x,'start ',4x,'  words ',/,2x,' no','. ' &
  ,4x,'type',3x,'record',6x,' stored ')
  9180 format (2x,i5,3i9)
  9190 format (2x,i5,3i9)
  9200 format (/,1x,'total words per block:'/ (1x,10i10))
  9210 format (/,1x,'total words:',1x,i12)
  9220 format (/,' ****** warning -', &
  ' difference > 0.0001 found between mubar on the gendf* tape ',/, &
  ' and mubar calculated from the elastic cross-section',' ******',/,/, &
  ' group    calculated    gendf*    difference ',/, &
  ' -----    ----------    ------    ----------')
  9230 format (i5,4x,1p,e12.5,e12.5,0p,f10.5)
  9240 format (/,'***note no elastic scattering matrix available for ', &
  'this nuclide. transport will be calculated as 1/3 total')
  end subroutine ecosec
  !
  subroutine ectsum(totinm,nreci,ing,ismtx,listmx,nsmxj,ntempj,itempj,tnmin, &
  tnmax,therma,lenrec,matrix,porder,lnrec1,lnrec2,npne,lunit,inmin,inmax,mdim, &
  pdim,thdim1,ngi,ingmin,ingmax,n)
  integer mdim,pdim,nreci,thdim1,ngi,ingmin,ingmax,n,nsmxj,ntempj,itempj
  real(kr) fact
  !  .. array arguments ..
  real(kr) matrix(mdim),porder(pdim)
  real(kr) totinm(ngrmax,maxmtx),therma(thdim1:ngi,thdim1:ngi,tmpmax)
  integer inmin(maxmtx,ngrmax),inmax(maxmtx,ngrmax),npne(maxmtx), &
  lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx),ing(nreci),lenrec(3), &
  lunit(maxmtx),tnmax(tmpmax,thdim1:ngi),tnmin(tmpmax,thdim1:ngi),ismtx(maxmtx)
  character listmx(maxmtx)*16
  !  .. local scalars ..
  integer mpos,ppos,ypos,i,j,ig,ighi,iglo,istp,istrt,ix,ng,nt
  real(kr) totnxn
  !----
  !  matrices
  !----
  do i = 1,nsmxj
    if(i==melas)then
      ix=nelas
    else if(i==mnxn)then
      ix=nnxn
    else if(i==minel)then
      ix=ninel
    else
      write(nsyso,'(//,5(a,i8),//)') ' i: ',i,' nsmxj: ', nsmxj, ' melas: ', &
       melas,' mnxn: ',  mnxn,  ' minel: ', minel
      call error('ectsum','i/=melas,mnxn,minel',' ')
    endif
    !----
    !  read in matrix and higher orders from buffer
    !----
    call buffin(i,n,lenrec(3),mpos,ppos,matrix,mdim,porder,pdim,lnrec1,lnrec2, &
    npne(i),lunit(i))
    ypos = mpos
    ng = 0
    do ig = ingmin,ingmax
      ng = ng + 1
      iglo = inmin(i,ig)
      ighi = inmax(i,ig)
      totinm(ig,i)=0.0
      fact = 1.0
      if(iglo/=0 .and. ighi/=0) then
        do j = iglo,ighi
          !----
          !  sum matrices
          !----
          totinm(ig,i) = totinm(ig,i)+matrix(mpos)
          mpos = mpos - 1
        enddo
        if(totinm(ig,i)/=0)fact = xsec(ix,1,ig)/totinm(ig,i)
        totnxn = 0.0
        if(i==mnxn)then
          !----
          !  calculate n,xn cross-section from partials
          !----
          if(rn2nd/=0)totnxn=totnxn+rf(rn2nd,ig)*2.0
          if(rn2n/=0)totnxn=totnxn+rf(rn2n,ig)*2.0
          if(rn3n/=0)totnxn=totnxn+rf(rn3n,ig)*3.0
          if(rn4n/=0)totnxn=totnxn+rf(rn4n,ig)*4.0
          if(rn2np/=0)totnxn=totnxn+rf(rn2np,ig)*2.0
          if(rn3np/=0)totnxn=totnxn+rf(rn3np,ig)*3.0
          if(rn2na/=0)totnxn=totnxn+rf(rn2na,ig)*2.0
          if(rn3na/=0)totnxn=totnxn+rf(rn3na,ig)*3.0
          if(rn2n2a/=0)totnxn=totnxn+rf(rn2n2a,ig)*2.0
          if(rnany/=0)totnxn=totnxn+rf(rnany,ig)*nu5(ig)
          if(totnxn/=0.0)fact = totnxn/totinm(ig,i)
        endif
        do j = iglo,ighi
          !----
          ! normalise matrices
          !----
          matrix(ypos) = matrix(ypos)*fact
          ypos = ypos - 1
        enddo
      endif
      !----
      !  temperature dependent matrices
      !----
      if((itempj>0).and.(ig>=itempj).and.i==melas) then
        do nt = 2,ntempj
          istrt = tnmin(nt,ig)
          istp = tnmax(nt,ig)
          totinm(ig,i)=0.0
          if((istrt/=0) .and. (istp/=0)) then
            do j = istrt,istp
              totinm(ig,i) = totinm(ig,i)+therma(ig,j,nt)
            enddo
            if(totinm(ig,i)>0.0) then
              fact = xsec(ix,nt,ig)/totinm(ig,i)
            else if(iprint>0) then
              write(nsyso,1000)
              write(nsyso,1100)i,ig,nt,xsec(ix,nt,ig),totinm(ig,i),fact
            endif
            therma(ig,istrt:istp,nt) = therma(ig,istrt:istp,nt) * fact
          endif
        enddo
      endif
      totinm(ig,i)=fact
    enddo ! end of group loop
    call ecfwrt(i,n,lenrec(3),mpos,ppos,matrix,mdim,porder,pdim,lnrec1,lnrec2 &
    ,npne(i),lunit(i))
  enddo ! end of matrix loop
  if(iprint>0) then
    if(nsmxj==1)write(nsyso,1220)(listmx(ismtx(i)),i=1,nsmxj),(j,(totinm(j,i), &
    i=1,nsmxj),j=ingmin,ingmax)
    if(nsmxj==2)write(nsyso,1210)(listmx(ismtx(i)),i=1,nsmxj),(j,(totinm(j,i), &
    i=1,nsmxj),j=ingmin,ingmax)
    if(nsmxj==3)write(nsyso,1200)(listmx(ismtx(i)),i=1,nsmxj),(j,(totinm(j,i), &
    i=1,nsmxj),j=ingmin,ingmax)
  endif
  return
  !----
  !  formats
  !----
  1000 format(' matrix group temp  cross-section  matrix terms  factor')
  1100 format(1x,3i5,3x,1p,3e14.5)
  1200 format(/,' factors used to scale matrices',/, &
  ' ------------------------------',/,' group        ',3a16,/,(1x,i4,3x,1p, &
  3e16.5))
  1210 format(/,' factors used to scale matrices',/, &
  ' ------------------------------',/,' group        ',2a16,/,(1x,i4,3x,1p, &
  2e16.5))
  1220 format(/,' factors used to scale matrices',/, &
  ' ------------------------------',/,' group        ',a16,/,(1x,i4,3x,1p, &
  e16.5))
  end subroutine ectsum
  !
  subroutine ecfil1(sig0,nsig0,ntw,en,btemp,ngrmax)
  !***************************************************************************
  ! this subroutine reads mf=1 and saves the sigma-0 values and the neutron
  ! group structure
  ! sig0 -
  ! nsig0
  ! en
  ! btemp
  ! ngrmax - maximum number of groups
  !***************************************************************************
  !.. scalar arguments ..
  real(kr) btemp
  integer nsig0,ntw,ngrmax
  !.. array arguments ..
  real(kr) en(ngrmax+1),sig0(10)
  real(kr) buff(maxnbf)
  !.. local scalars ..
  integer igr,isig0,l,ll,nb,ng2,nw
  !
  !***input.
  if(iprint>0) write(nsyso,1030)
  nsig0 = l2h
  if(nsig0>10) call error('ecfil1',' *** stop  ===> increase sig0',' ')
  l = 1
  call listio(ning,0,0,buff(1),nb,nw)
  if(math/=matno)then
    write(hsmg,'(15hmaterial number,i5,14h requested but,i5,6h found)')
    call error('ecfil1',hsmg,' ')
  endif
  ngi = l1h
  btemp = c1h
  10 if(nb==0) go to 20
  l = l + nw
  if(l>maxnbf) call error('ecfil1',' *** buffer size exceeded',' ')
  call moreio(ning,0,0,buff(l),nb,nw)
  go to 10
  20 continue
  ll = 6 + ntw
  do isig0 = 1,nsig0
    ll = ll + 1
    sig0(isig0) = buff(ll)
  enddo
  ng2 = ngi + 1
  do igr = 1,ng2
    ll = ll + 1
    en(igr) = buff(ll)
  enddo
  call tosend(ning,0,0,buff(1))
  return
  !----
  !  formats
  !----
  1030 format (/,30('-'),' group structure ',30('-'),/)
  end subroutine ecfil1
  !
  subroutine ecfil3(flux,eng)
  ! this subroutine reads the averaged cross section data on unit ning
  ! .. array arguments ..
  real(kr) eng(ngrmax+1),flux(ngrmax)
  ! .. local scalars ..
  real(kr) delta,diff
  integer i,ifail,ig,in,in2n,ipos,j,jpos,ll,mg,mmg,nb,ng2,nleg,nprint,nsig, &
  nsmci,ntempj,nthres,nw,sumfis,sumine
  character blank*16,reac*16
  ! .. local arrays ..
  real(kr) buff(maxnbf),total(ngrmax),totfis(ngrmax),totine(ngrmax)
  !     temp for thermal scatter
  real(kr) tempscat(ngrmax)
  integer ismci(maxel)
  ! .. intrinsic functions ..
  intrinsic abs
  integer trouve
  !
  if(iprint>0) write(nsyso,9120)
  ! initialise local scalars
  nb = 0
  ifail = 0
  in2n = 0
  ipos = 0
  jpos = 0
  nleg = 0
  nsig = 0
  nsmci = 0
  nthres = 0
  nw = 0
  sumfis = 0
  sumine = 0
  ! set defaults
  delay = .false.
  blank = '   -            '
  ntempj = ntemp(nelp)
  nfissn = 0
  ntotal = 0
  ntran = 0
  nelas = 0
  ninel = 0
  nnxn = 0
  ncapt = 0
  nprint = 0
  if(iprint>0) write(nsyso,9000)
  ! zeroise arrays for summing
  ismci(:maxel)=0
  total(:ngi)=0.0d0
  totfis(:ngi)=0.0d0
  totine(:ngi)=0.0d0
  tempscat(:ngrmax)=0.0d0
  if(ntempj==1) then
    nfre = 0
    nud(:ngi)=0.0d0
    rf(:nfr,:ngi) = 0.0d0
  endif
  xsec(:nsmic,ntempj,:ngi) = 0.0d0
  !----
  !  start loop over reactions on tape ning
  !----
  50 call contio(ning,0,0,buff,nb,nw)
  nleg = l1h
  nsig = l2h
  if(mfh==3) then
    ! cross-sections
    if(mth>19 .and. nprint==0 .and. nfissn/=0) then
      nprint = 1
      if(iprint>0) write(nsyso,9010) mfh,18,listmi(ismci(nfissn)), &
      listmi(ismci(ntotal)),blank
    endif
    if(mth==1) then
      ! total
      nsmci = nsmci + 1
      reac = 'TOTAL'
      ! search for reaction on reference table
      call search(reac,ismci(nsmci),listmi,nsmic)
      ntotal = nsmci
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank,blank, &
      & listmi(ismci(ntotal))
    else if(mth==2) then
      ! elastic
      nsmci = nsmci + 1
      reac = 'ELASTIC'
      call search(reac,ismci(nsmci),listmi,nsmic)
      nelas = nsmci
      if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(nsmci)), &
      listmi(ismci(ntotal)),blank
    else if(mth==4) then
      ! inelastic
      sumine = 1
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank,blank, &
      & 'INELASTIC       '
    else if(mth==18) then
      ! fission
      nsmci = nsmci + 1
      reac = 'FISSION'
      call search(reac,ismci(nsmci),listmi,nsmic)
      nfissn = nsmci
    else if(mth==19 .or. mth==20 .or. mth==21 .or.mth==38) then
      sumfis = 1
      if(nfissn==0) then
        ! sum partials to form total fission
        nsmci = nsmci + 1
        reac = 'FISSION'
        call search(reac,ismci(nsmci),listmi,nsmic)
        nfissn = nsmci
      endif
      if(ntempj==1) then
        nfre = nfre + 1
        if(mth==19) then
          reac = 'N,F'
        else if(mth==20) then
          reac = 'N,N+F'
        else if(mth==21) then
          reac = 'N,2N+F'
        else if(mth==38) then
          reac = 'N,3N+F'
        endif
        call search(reac,ifre(nfre),listfr,nfr)
        if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(nfissn)), &
        listmi(ismci(ntotal)),listfr(ifre(nfre))
      else
        if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(nfissn)), &
        listmi(ismci(ntotal)),blank
      endif
    else if(ntrouve(mth,mtherm,ntherm)/=0) then
      ! thermal scatter
      if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(nelas)), &
      & listmi(ismci(ntotal)),blank,'thermal scatter '
    else if(mth==452) then
      ! nubar
      nsmci = nsmci + 1
      reac = 'NU*FISSION'
      nnufis = nsmci
      call search(reac,ismci(nsmci),listmi,nsmic)
      if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(nsmci)),blank,blank
    else if(mth==455) then
      ! nubar - delayed
      delay = .true.
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank,blank, &
      'delayed-nubar us','ed to form the f','ission spectrum '
    else if(mth==456) then
      ! nubar - prompt
      ! prompt = .true.
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank,blank, &
      'prompt-nubar use','d to form the fi','ssion spectrum  '
    else if(mth==251) then
      ! mubar
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank,blank, &
      & 'transport        '
      ! calculate n,xn
    else if(mth==11 .or. mth==16 .or. mth==17 .or. mth==24 .or. mth==25 .or. &
           & mth==30 .or. mth==37 .or. mth==41 .or. mth==42 .or. mth==5  .or. &
           & (mth>= 875 .and. mth<=891)) then
      if(nnxn==0) then
        nsmci = nsmci + 1
        nnxn = nsmci
        reac = 'N,XN'
        call search(reac,ismci(nsmci),listmi,nsmic)
        in2n = nsmci
      endif
      if(ntempj==1) then
        if(mth<876) then
          nfre = nfre + 1
        else if(reac/='n,2n') then
          nfre = nfre + 1
        endif
        if(mth==11) then
          rn2nd = nfre
          reac = 'N,2N+D'
        else if(mth==16) then
          rn2n = nfre
          reac = 'N,2N'
          if(pmt16)then
            if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
            & 'replaced by partials'
            call tosend(ning,0,0,buff)
            nfre = nfre - 1
            go to 50
          endif
        else if(mth==17) then
          rn3n = nfre
          reac = 'N,3N'
        else if(mth==24) then
          rn2na = nfre
          reac = 'N,2N+ALPHA'
        else if(mth==25) then
          rn3na = nfre
          reac = 'N,3N+ALPHA'
        else if(mth==30) then
          rn2n2a = nfre
          reac = 'N,2N+2ALPHA'
        else if(mth==37) then
          rn4n = nfre
          reac = 'N,4N'
        else if(mth==41) then
          rn2np = nfre
          reac = 'N,2NP'
        else if(mth==42) then
          rn3np = nfre
          reac = 'N,3NP'
        else if(mth==5) then
          rnany = nfre
          reac = 'N,ANYTHING'
        else if(mth>= 875 .and. mth<=891) then
          rn2n = nfre
          reac = 'N,2N'
        endif
        call search(reac,ifre(nfre),listfr,nfr)
        if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(nnxn)), &
        listmi(ismci(ntotal)),listfr(ifre(nfre))
      else
        if(pmt16)then
          if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
          & 'replaced by partials'
          call tosend(ning,0,0,buff)
          go to 50
        else
          if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(nnxn)), &
          listmi(ismci(ntotal)),blank
        endif
      endif
    else if((mth>=102 .and. mth<=109) .or.(mth>=111 .and. mth<=117) .or. &
           & (mth>=600 .and. mth<=849)) then
      ! capture
      if(ncapt==0) then
        nsmci = nsmci + 1
        reac = 'CAPTURE'
        call search(reac,ismci(nsmci),listmi,nsmic)
        ncapt = nsmci
      endif
      if(ntempj==1) then
        nfre = nfre + 1
        if(mth==102) then
          reac = 'N,GAMMA'
        else if(mth==103) then
          reac = 'N,P'
          if(pmt103) then
            if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
            & 'replaced by partials'
            call tosend(ning,0,0,buff)
            nfre = nfre - 1
            go to 50
          endif
        else if(mth==104) then
          reac = 'N,D'
          if(pmt104) then
            if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
            & 'replaced by partials'
            call tosend(ning,0,0,buff)
            nfre = nfre - 1
            go to 50
          endif
        else if(mth==105) then
          reac = 'N,T'
          if(pmt105) then
            if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
            & 'replaced by partials'
            call tosend(ning,0,0,buff)
            nfre = nfre - 1
            go to 50
          endif
        else if(mth==106) then
          reac = 'N,HE3'
          if(pmt106) then
            if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
            & 'replaced by partials'
            call tosend(ning,0,0,buff)
            nfre = nfre - 1
            go to 50
          endif
        else if(mth==107) then
          reac = 'N,ALPHA'
          if(pmt107) then
            if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
            & 'replaced by partials'
            call tosend(ning,0,0,buff)
            nfre = nfre - 1
            go to 50
          endif
        else if(mth==108) then
          reac = 'N,2ALPHA'
        else if(mth==109) then
          reac = 'N,3ALPHA'
        else if(mth==111) then
          reac = 'N,2P'
        else if(mth==112) then
          reac = 'N,P+ALPHA'
        else if(mth==113) then
          reac = 'N,T+2ALPHA'
        else if(mth==114) then
          reac = 'N,D+2ALPHA'
        else if(mth==115) then
          reac = 'N,PD'
        else if(mth==116) then
          reac = 'N,PT'
        else if(mth==117) then
          reac = 'N,D+ALPHA'
        else if(mth>=600 .and. mth<=649) then
          if(reac=='n,p') nfre = nfre - 1
          reac = 'N,P'
        else if(mth>=650 .and. mth<=699) then
          if(reac=='n,d') nfre = nfre - 1
          reac = 'N,D'
        else if(mth>=700 .and. mth<=749) then
          if(reac=='n,t') nfre = nfre - 1
          reac = 'N,T'
        else if(mth>=750 .and. mth<=799) then
          if(reac=='n,he3') nfre = nfre - 1
          reac = 'N,HE3'
        else if(mth>=800 .and. mth<=849) then
          if(reac=='n,alpha') nfre = nfre - 1
          reac = 'N,ALPHA'
        endif
        call search(reac,ifre(nfre),listfr,nfr)
        if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(ncapt)), &
        listmi(ismci(ntotal)),listfr(ifre(nfre))
      else
        if((mth==103.and.pmt103).or.(mth==104.and.pmt104).or. &
          & (mth==105.and.pmt105).or.(mth==106.and.pmt106).or. &
          & (mth==107.and.pmt107)) then
          if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank, &
          & 'replaced by partials'
          call tosend(ning,0,0,buff)
          go to 50
        else
          if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(ncapt)), &
          listmi(ismci(ntotal)),blank
        endif
      endif
      ! calculate total inelastic
    else if((mth>=22.and.mth<=23) .or.(mth>=28.and.mth<=29) .or. &
           & (mth>=32.and.mth<=36) .or.(mth>=44.and.mth<=45) .or. &
           & (mth>=51.and.mth<=91)) then
      if(ninel==0) then
        nsmci = nsmci + 1
        reac = 'INELASTIC'
        call search(reac,ismci(nsmci),listmi,nsmic)
        ninel = nsmci
      endif
      if((ntempj==1) .and. (mth<51)) then
        nfre = nfre + 1
        if(mth==22) then
          reac = 'N,N+ALPHA'
        else if(mth==23) then
          reac = 'N,N+3ALPHA'
        else if(mth==28) then
          reac = 'N,N+P'
        else if(mth==29) then
          reac = 'N,N+2ALPHA'
        else if(mth==32) then
          reac = 'N,N+D'
        else if(mth==33) then
          reac = 'N,N+T'
        else if(mth==34) then
          reac = 'N,N+HE3'
        else if(mth==35) then
          reac = 'N,N+D+2ALPHA'
        else if(mth==36) then
          reac = 'N,N+T+2ALPHA'
        else if(mth==44) then
          reac = 'N,N+2P'
        else if(mth==45) then
          reac = 'N,N+P+ALPHA'
        endif
        call search(reac,ifre(nfre),listfr,nfr)
        if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(ninel)), &
        listmi(ismci(ntotal)),listfr(ifre(nfre))
      else if(mth>=51 .and. mth<=91) then
        if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(ninel)), &
        listmi(ismci(ntotal)),blank
      else if(ntempj>1) then
        if(iprint>0) write(nsyso,9010) mfh,mth,listmi(ismci(ninel)), &
        listmi(ismci(ntotal)),blank
      endif
    else if(mth==26) then
      ! (n,abs)
      write(nsyso,9180)
      stop
    else if(mth==0) then
      go to 50
    else
      !  reaction not recognised
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,blank,'not recognized  '
      call tosend(ning,0,0,buff)
      go to 50
    endif
    !----
    !  start loop over groups on tape ning
    !----
    60 ll = 1
    call listio(ning,0,0,buff(1),nb,nw)
    70 if(nb/=0) then
      if((ll+nw)>maxbuf) then
        write(nsyso,9080)
        call error('ecfil3',' workspace exceeded reading gendf* **** ',' ')
      else
        ll = ll + nw
        call moreio(ning,0,0,buff(ll),nb,nw)
        go to 70
      endif
    endif
    ng2 = l1h
    nw = n1h
    ig = n2h
    ipos = (nleg*nsig* (ng2-1)) + 7
    mg = ngi - ig + 1
    if(mth==1) then
      ! total
      total(mg) = total(mg) + buff(ipos)
    else if(mth==4) then
      ! inelastic
      totine(mg) = totine(mg) + buff(ipos)
    else if(mth==18) then
      ! fission
      totfis(mg) = totfis(mg) + buff(ipos)
    else if(mth==2) then
      ! elastic
      xsec(nsmci,ntempj,mg) = buff(ipos)
      xsec(ntotal,ntempj,mg) = xsec(ntotal,ntempj,mg) +buff(ipos)
    else if(ntrouve(mth,mtherm,ntherm)/=0) then
      ! thermal scatter
      ! cut out test on small thermal scatter believe wrong
      ! caused a problem in jeff3 be9 where coherent elastic is 1e-14 in
      ! bottom 2 groups of 1968 structure.
      ! test needed becaue the terminator record contains a cross section of
      ! zero and number of thermal groups is set to 1968.
      if(buff(ipos)>0.0) then
        ! correct total
        xsec(ntotal,ntempj,mg) = xsec(ntotal,ntempj,mg) + buff(ipos) - &
        & xsec(nelas,ntempj,mg)
        ! replace elastic
        tempscat(mg)=tempscat(mg)+buff(ipos)
        xsec(nelas,ntempj,mg) = 0.0
        itemp(nelp) = mg
      endif
    else if(mth==452) then
      ! nu
      jpos = nleg*nsig + 7
      xsec(nsmci,ntempj,mg) = buff(jpos)
      if(ntempj==1) flux(mg) = buff(7)
    else if(ntempj==1) then
      if(mth==251) then
        ! mubar
        jpos = nleg*nsig + 7
        mubar(mg) = buff(jpos)
      else if(mth==455) then
        ! nubar -delayed
        jpos = nleg*nsig + 7
        nud(mg) = buff(jpos)
      endif
    endif
    ! calculate n,xn
    if(mth==16 .or. mth==17 .or. mth==24 .or. mth==25 .or.mth==30 .or. mth==37 &
      & .or. mth==41 .or. mth==42 .or. mth==5 .or. (mth>=875 .and. mth<=891)) &
      & then
      xsec(in2n,ntempj,mg) = xsec(in2n,ntempj,mg) + buff(ipos)
      xsec(ntotal,ntempj,mg) = xsec(ntotal,ntempj,mg) +buff(ipos)
      if(ntempj==1) rf(nfre,mg) = rf(nfre,mg)+buff(ipos)
    endif
    ! calculate total fission
    if(mth==19 .or. mth==20 .or. mth==21 .or. mth==38) then
      ! sum partials
      xsec(nfissn,ntempj,mg) = xsec(nfissn,ntempj,mg) +buff(ipos)
      xsec(ntotal,ntempj,mg) = xsec(ntotal,ntempj,mg) +buff(ipos)
      if(ntempj==1) rf(nfre,mg) = buff(ipos)
    endif
    ! calculate capture
    if((mth>=102 .and. mth<=109) .or.(mth>=111 .and. mth<=116) .or. &
      & (mth>=600 .and. mth<=849)) then
      xsec(ncapt,ntempj,mg) = xsec(ncapt,ntempj,mg) + buff(ipos)
      xsec(ntotal,ntempj,mg) = xsec(ntotal,ntempj,mg) +buff(ipos)
      if(ntempj==1) rf(nfre,mg) = rf(nfre,mg)+buff(ipos)
    endif
    ! calculate total inelastic
    if((mth>=22.and.mth<=23) .or. (mth>=28.and.mth<=29) .or. &
      & (mth>=32.and.mth<=36) .or. (mth==44) .or.(mth==45) .or. &
      & (mth>=51.and.mth<=91)) then
      xsec(ninel,ntempj,mg) = xsec(ninel,ntempj,mg) + buff(ipos)
      xsec(ntotal,ntempj,mg) = xsec(ntotal,ntempj,mg) +buff(ipos)
      if(mth<51 .or. mth>91) then
        totine(mg) = totine(mg) + buff(ipos)
        if(ntempj==1) rf(nfre,mg) = buff(ipos)
      endif
    endif
    if(ig<ngi) then
      ll = 1
      call listio(ning,0,0,buff(1),nb,nw)
      go to 70
    endif
    call tosend(ning,0,0,buff)
    go to 50
  endif
  if(mfh/=0) call skiprz(ning,-1)
  ! maj elastic
  do ig=1,ngrmax
    if(tempscat(ig)>0.0) xsec(nelas,ntempj,ig) = tempscat(ig)
  enddo
  ! add transport to list of primary cross-sections
  nsmci = nsmci + 1
  reac = 'TRANSPORT'
  call search(reac,ismci(nsmci),listmi,nsmic)
  ntran = nsmci
  if(ntempj==1) then
    nsmc = nsmci
    ismc(:nsmc) = ismci(:nsmc)
  else if(nsmc/=nsmci) then
    ! error. different number of reactions per temperature
    write(nsyso,9070) nsmc,nsmci
    write(nsyso,*) (listmi(ismc(i)),i=1,nsmc)
    write(nsyso,*) (listmi(ismci(i)),i=1,nsmci)
    call error('ecfil3','different number of primary cross sections',' ')
  endif
  ! write out total fission if calculated
  if(sumfis==1) then
    ! calculate difference
    if(iprint>0) then
      write(nsyso,9020)
      write(nsyso,9045)
    endif
    ifail = 0
    mmg = ngi + 1
    do ig = 1,ngi
      delta = abs(xsec(nfissn,ntempj,ig)-totfis(ig))
      if(xsec(nfissn,ntempj,ig)/=0.0) then
        diff = (delta/xsec(nfissn,ntempj,ig))*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9050) ig,eng(mmg-1),eng(mmg),totfis(ig), &
        xsec(nfissn,ntempj,ig),diff
        ifail = 1
      endif
      mmg = mmg - 1
    enddo
    if((ifail==0).and.(iprint>0)) write(nsyso,9150)
  else if(nfissn/=0) then
    ! no partial fission on tape so total fission should be used
    do ig = 1,ngi
      xsec(nfissn,ntempj,ig) = totfis(ig)
      xsec(ntotal,ntempj,ig) = xsec(ntotal,ntempj,ig) +xsec(nfissn,ntempj,ig)
    enddo
  endif
  !  calculate nu*fission
  if(nnufis/=0) then
    do ig = 1,ngi
      xsec(nnufis,ntempj,ig) = xsec(nnufis,ntempj,ig)*xsec(nfissn,ntempj,ig)
    enddo
  endif
  ! write out total inelastic if calculated
  if(sumine==1) then
    ! calculate difference
    if(iprint>0) then
      write(nsyso,9030)
      write(nsyso,9045)
    endif
    ifail = 0
    mmg = ngi + 1
    do ig = 1,ngi
      delta = abs(xsec(ninel,ntempj,ig)-totine(ig))
      if(xsec(ninel,ntempj,ig)/=0.0) then
        diff = (delta/xsec(ninel,ntempj,ig))*100.0
      else
        diff = 0.0
      endif
      ! if difference greater than .01 %
      if(diff>0.01) then
        if(iprint>0) write(nsyso,9050) ig,eng(mmg-1),eng(mmg),totine(ig), &
        xsec(ninel,ntempj,ig),diff
        ifail = 1
      endif
      mmg = mmg - 1
    enddo
    if((ifail==0).and.(iprint>0)) write(nsyso,9150)
  endif
  ! write out total
  if(iprint>0) then
    write(nsyso,9040)
    write(nsyso,9045)
  endif
  mmg = ngi + 1
  ifail = 0
  do ig = 1,ngi
    delta = abs(xsec(ntotal,ntempj,ig)-total(ig))
    diff = (delta/xsec(ntotal,ntempj,ig))*100.0
    ! if difference greater than 1 %
    if(diff>1.0) then
      if(iprint>0) write(nsyso,9060) ig,eng(mmg-1),eng(mmg),total(ig), &
      xsec(ntotal,ntempj,ig),diff,delta
      ifail = 1
    else if(diff>0.01) then
      if(iprint>0) write(nsyso,9130) ig,eng(mmg-1),eng(mmg),total(ig), &
      xsec(ntotal,ntempj,ig),diff,delta
      if(ifail/=1)ifail = 2
    endif
    mmg = mmg - 1
  enddo
  if(iprint>0) then
    if(ifail==1) then
      write(nsyso,9140)
    else if(ifail==0) then
      write(nsyso,9150)
    endif
    ! list reactions to be processed
    write(nsyso,9090)
    write(nsyso,9100) (listmi(ismc(i)),i=1,nsmc)
    write(nsyso,9110)
    write(nsyso,9100) (listfr(ifre(i)),i=1,nfre)
  endif
  ! check array bounds not exceeded for tnmin and tnmax
  if(itemp(nelp)/=0) then
    nthres = ngi - itemp(nelp) + 1
    if(iprint>0) write(nsyso,9170) nthres
    if(nthres>maxthr) then
      write(nsyso,9160) nthres,maxthr
      call error('ecfil3','array bound exceeded',' ')
    endif
  endif
  return
  !----
  !  formats
  !----
  9000 format ('  mf  mt',16x,'primary',11x,'response',8x, &
  'cross-check/comment',/,'  ==  ==',16x,'=======',11x,'========',8x, &
  '===================',/)
  9010 format (1x,i3,1x,i3,2x,8a16,/)
  9020 format (/,30x,' fission cross-check',/,31x,22 ('='),/)
  9030 format (/,30x,' inelastic cross-check',/,31x,22 ('='),/)
  9040 format (/,30x,' total cross-check',/,31x,17 ('='),/, &
  ' note: the calculated total will be written ',' to the ecco library',/)
  9045 format (' group  lower energy  upper ','energy     gendf*      ', &
  & 'calculated    diff(%)    diff',/, &
  & ' -----  ------------  ------------', &
  & '   ----------    ----------    -------------',/)
  9050 format (1x,i4,1x,1p,4e14.6,0p,f10.3)
  9060 format (1x,i4,1x,1p,4e14.6,0p,2f10.4,'*****')
  9070 format (/,' ***** error - ,',' different number of primary cross-', &
  'sections. found ',i5,' for this temp and ',i5,' for the first temp ****')
  9080 format (/,' ***** error - ',' workspace exceeded reading gendf* **** ')
  9090 format (/,15x,' primary cross-sections to be processed',/,16x,37 ('-'))
  9100 format (1x,5a16,/)
  9110 format (/,15x,' response functions to be processed',/,16x,33 ('-'))
  9120 format (/,30 ('-'),' cross-section data ',30 ('-'),/)
  9130 format (1x,i4,1x,1p,4e14.6,0p,2f10.4)
  9140 format (/,1x,6 ('**** error ***'),/, &
  ' differences greater than 1% found between the total ','cross-section on', &
  /,' the gendf* and the sum of the ','partials. ',/,1x,6 ('**** error ***'))
  9150 format (' no differences greater than 0.01% found ')
  9160 format (' array bound exceeded. number of groups with thermal', &
  ' scattering data is ',i5,' gecco maximum is ',i5,' (parameter maxthr)')
  9170 format (/,' number of groups with thermal',' scattering data is ',i5)
  9180 format (' absorption cross section (mt 27) present on gendf*', &
  ' but gecco not programmed to process it correctly ')
  contains
    !
    integer function ntrouve(mth,mtherm,ntherm)
    integer mth,mtherm,ii,ntherm(mtherm)
    ntrouve = 0
    if(mtherm/=0) then
      do ii=1,mtherm
        if(mth==ntherm(ii)) then
          ntrouve = ii
          return
        endif
      enddo
    endif
  end function ntrouve
  end subroutine ecfil3
  !
  subroutine ecfil5
  ! .. local scalars ..
  integer ig,ipos,it,ll,mg,nb,ng2,ntime,nw
  ! .. local arrays ..
  real(kr) buff(maxnbf)
  ! ..
  if(iprint>0) write(nsyso,9030)
  if(iprint>0) write(nsyso,9000)
  ! zeroise chid
  chid(:ngi)=0.0d0
  ! start loop over sections on tape ning
  10 call contio(ning,0,0,buff,nb,nw)
  if(mfh==5) then
    ntime = l1h
    if(mth==455) then
      ! delayed chi
      if(iprint>0) write(nsyso,9010) mfh,mth,'delayed-chi used', &
      & ' to form the fis','sion spectrum   '
      ll = 1
      call listio(ning,0,0,buff(1),nb,nw)
      20 if(nb/=0) then
        if((ll+nw)>maxbuf) then
          call error('ecfil5','workspace exceeded reading gendf* (1)',' ')
        endif
        ll = ll + nw
        call moreio(ning,0,0,buff(ll),nb,nw)
        go to 20
      endif
      ng2 = l1h
      nw = n1h
      ipos = 7 + ntime
      do ig = 1,ng2 - 1
        mg = ngi - ig + 1
        do it = 1,ntime
          chid(mg) = chid(mg) + buff(ipos)
          ipos = ipos + 1
        enddo
      enddo
    else if(mth==452) then
      ! 1 mev fission spectrum
      if(iprint>0) write(nsyso,9010) mfh,mth,'fission spectrum', &
      ' at 1 mev (only',' used when fissi','on matrix not pr','esent)          '
      ll = 1
      call listio(ning,0,0,buff(1),nb,nw)
      50 if(nb/=0) then
        if((ll+nw)>maxbuf) then
          call error('ecfil5','workspace exceeded reading gendf* (2)',' ')
        endif
        ll = ll + nw
        call moreio(ning,0,0,buff(ll),nb,nw)
        go to 50
      endif
      ng2 = l1h
      nw = n1h
      ipos = 8
      do ig = 1,ng2 - 1
        mg = ngi - ig + 1
        spec(mg) = spec(mg) + buff(ipos)
        ipos = ipos + 1
      enddo
    else
      !  reaction not recognised
      if(iprint>0) write(nsyso,9010) mfh,mth,'not recognized  '
    endif
    call tosend(ning,0,0,buff)
    go to 10
  else
    call skiprz(ning,-1)
  endif
  ! end loop over sections on tape ning
  return
  !----
  !  formats
  !----
  9000 format ('  mf  mt',13x,'comment',/,'  ==  ==',13x,'=======',/)
  9010 format (1x,i3,1x,i3,2x,8a16,/)
  9020 format (/,' ****** error - workspace exceeded reading gendf* ',' ******')
  9030 format (/,30 ('-'),' energy distribution data ',30 ('-'),/)
  end subroutine ecfil5
  !
  subroutine ecfil6(nreci,ing,flux,lenrec,lwordr,matrix,matrx2,porder,pordr2, &
  maxm,maxp,lunit,itempj,thrmin,fision)
  ! this subroutine reads the neutron scattering matrices on unit ning
  ! and sets the arrays inmin and inmax
  ! .. scalar arguments ..
  integer itempj,lwordr,maxm,maxp,nreci,thrmin
  logical fision
  logical pmt16,pmt103,pmt104,pmt105,pmt106,pmt107
  ! .. array arguments ..
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  real(kr) flux(ngrmax)
  integer ing(nreci),lenrec(3),lunit(maxmtx)
  ! .. local scalars ..
  real(kr) bot,diff,sum
  integer i,ig,lrecl,mg,nb,nleg,npne1,nsig,nsmxj,nw,totf
  logical ifail
  character blank,reac*16
  ! .. local arrays ..
  real(kr) buff(maxnbf),fispec(ngrmax),profis(ngrmax),xdelay(ngrmax)
  real(kr) tbuff(maxnbf),tmatrix(msize),tmatrx2(msize),tporder(psize), &
  & tpordr2(psize)
  ! .. intrinsic functions ..
  intrinsic abs
  !
  blank = '   -            '
  if(iprint>0) write(nsyso,9040)
  if(iprint>0) write(nsyso,9000)
  ! set defaults
  fision = .false.
  melas = 0
  minel = 0
  mnxn = 0
  nleg = 0
  bot = 0.0
  diff = 0.0
  sum = 0.0
  lrecl = 0
  nsmxj = 0
  totf = 0
  npne1 = 0
  nsig = 0
  ! zeroise arrays
  fispec(:ngi)=0.0d0
  profis(:ngi)=0.0d0
  tbuff(:maxnbf)=0.0d0
  tmatrix(:msize)=0.0d0
  tmatrx2(:msize)=0.0d0
  tporder(:psize)=0.0d0
  tpordr2(:psize)=0.0d0
  xdelay(:ngi)=0.0d0
  inmax(:maxmtx,:ngi)=0
  inmin(:maxmtx,:ngi)=0
  ! open scratch files
  lrecl = lenrec(1)*lwordr
  do i = 1,nsmtx
    open(80+i,access='direct',status='scratch',recl=lrecl,form='unformatted')
    open(85+i,access='direct',status='scratch',recl=lrecl,form='unformatted')
    open(90+i,access='direct',status='scratch',recl=lrecl,form='unformatted')
    open(95+i,access='direct',status='scratch',recl=lrecl,form='unformatted')
  enddo
  ! read gendf
  40 continue
  call contio(ning,0,0,buff,nb,nw)
  if(mfh==6) then
    nleg = l1h
    nsig = l2h
    if(mth==2) then
      reac = 'ELASTIC'
      call ecfout(nsmxj,ing,nreci,nleg,nsig,lenrec(3),reac,maxm,maxp, &
      matrix,porder,ngi,buff,maxnbf,inmin,inmax,lnrec1,lnrec2,npne,ismtx, &
      listmx,nsmtx,thrmin)
      lunit(nsmxj) = 80 + nsmxj
      melas = nsmxj
    else if(mth==5 .or. mth==11 .or. mth==16 .or. mth==17 .or. mth==24 .or. &
          &  mth==25 .or.mth==30 .or. mth==37 .or. mth==41 .or. mth==42 .or. &
          & (mth>=875 .and. mth<= 891)) then
      if(mnxn==0) then
        reac = 'N,XN'
        call ecfout(nsmxj,ing,nreci,nleg,nsig,lenrec(3),reac,maxm,maxp, &
        matrix,porder,ngi,buff,maxnbf,inmin,inmax,lnrec1,lnrec2,npne,ismtx, &
        listmx,nsmtx,thrmin)
        lunit(nsmxj) = 80 + nsmxj
        mnxn = nsmxj
        if(iprint>0) write(nsyso,9010) mfh,mth,listmx(ismtx(mnxn))
      else
        if(iprint>0) write(nsyso,9010) mfh,mth,listmx(ismtx(mnxn))
        npne1 = npne(mnxn)
        if(npne1==1 .and. nleg==1) then
          call ecfad1(mnxn,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix,buff, &
          maxnbf,lunit(mnxn),lnrec1,inmin,inmax,maxm)
        else if(npne1>1 .and. npne1==nleg) then
          call ecfad2(mnxn,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(mnxn),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,maxm,maxp)
        else if(npne1>1 .and. nleg==1) then
          call ecfad3(mnxn,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(mnxn),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,npne1,maxm,maxp)
        else if(npne1==1 .and. nleg>1) then
          npne(mnxn) = nleg
          call ecfad4(mnxn,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(mnxn),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,npne1,maxm,maxp)
        else
          npne(mnxn) = nleg
          call ecfad5(mnxn,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(mnxn),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,npne1,maxm,maxp)
        endif
      endif
    else if((mth>=22.and.mth<=23) .or.(mth>=28.and.mth<=29) .or. &
           & (mth>=32.and.mth<=36) .or.(mth>=44.and.mth<=45) .or. &
           & (mth>=51.and.mth<=91)) then
      if(minel==0) then
        reac = 'INELASTIC'
        call ecfout(nsmxj,ing,nreci,nleg,nsig,lenrec(3),reac,maxm,maxp, &
        matrix,porder,ngi,buff,maxnbf,inmin,inmax,lnrec1,lnrec2,npne,ismtx, &
        listmx,nsmtx,thrmin)
        lunit(nsmxj) = 80 + nsmxj
        minel = nsmxj
        if(iprint>0) write(nsyso,9010) mfh,mth,listmx(ismtx(minel))
      else
        if(iprint>0) write(nsyso,9010) mfh,mth,listmx(ismtx(minel))
        npne1 = npne(minel)
        if(npne1==1 .and. nleg==1) then
          call ecfad1(minel,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix,buff, &
          maxnbf,lunit(minel),lnrec1,inmin,inmax,maxm)
        else if(npne1>1 .and. npne1==nleg) then
          call ecfad2(minel,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(minel),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,maxm,maxp)
        else if(npne1>1 .and. nleg==1) then
          call ecfad3(minel,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(minel),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,npne1,maxm,maxp)
        else if(npne1==1 .and. nleg>1) then
          npne(minel) = nleg
          call ecfad4(minel,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(minel),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,npne1,maxm,maxp)
        else
          npne(minel) = nleg
          call ecfad5(minel,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(minel),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,npne1,maxm,maxp)
        endif
      endif
    else if(mth==18) then
      fision = .true.
      call ecfissi(nleg,nsig,flux,totf,profis,fispec,xdelay)
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,'chi used to form', &
      ' the fission spe','ctrum           '
    else if(mth==19 .or. mth==20 .or. mth==21 .or.mth==38) then
      if(mth==19) spec(:ngrmax)=0.0d0
      fision = .true.
      call ecrfis(nleg,nsig,flux,totf,profis,fispec,xdelay)
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,'chi used to form', &
      ' the fission spe','ctrum           '
    else if(ntrouve(mth,mtherm,ntherm)/=0) then
      if(melas==0) then
        reac = 'ELASTIC'
        call ecfout(nsmxj,ing,nreci,nleg,nsig,lenrec(3),reac,maxm,maxp, &
        matrix,porder,ngi,buff,maxnbf,inmin,inmax,lnrec1,lnrec2,npne,ismtx, &
        listmx,nsmtx,thrmin)
        lunit(nsmxj) = 80 + nsmxj
        melas = nsmxj
      else
        if(ntrouve(mth,mtherm,ntherm)==1) then
          call ecfthr(melas,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix, &
          buff,maxnbf,lunit(melas),lnrec1,lnrec2,inmin,inmax,porder,pordr2, &
          nleg,maxm,maxp,itempj,thrmin)
        else
          call ecfadth(melas,ing,nreci,nsig,lenrec(3),matrx2,ngi,matrix,buff, &
          maxnbf,lunit(melas),lnrec1,lnrec2,inmin,inmax,porder,pordr2,nleg, &
          maxm,maxp,itempj,thrmin)
        endif
      endif
      if(iprint>0) write(nsyso,9010) mfh,mth,listmx(ismtx(melas)), &
      & 'thermal scatter '
      !----
      !  print elastic matrix coding from ectsum
      !  read in matrix and higher orders from buffer
      !----
    else
      if(iprint>0) write(nsyso,9010) mfh,mth,blank,' not recognized'
    endif
    if(nleg-1>maxleg) then
      write(nsyso,9080) nleg - 1,maxleg
      call error('ecfil6','maximum number of legendre orders exceeded',' ')
    endif
    if(mth/=0) call tosend(ning,0,0,buff)
    go to 40
  else
    call skiprz(ning,-1)
  endif
  nsmx(nelp) = nsmxj
  if(iprint>0) write(nsyso,9020) maxm,maxp
  ! list reactions to be processed
  if(iprint>0) write(nsyso,9030) (listmx(ismtx(i)),i=1,nsmxj)
  !
  if(fision) then
    ! normalise fission spectrum
    bot = 0.0
    do mg = 1,ngi
      bot = bot + xsec(nnufis,1,mg)*flux(mg)
    enddo
    if(iprint>0) write(nsyso,9060)
    sum = 0.0
    ifail = .false.
    do mg = 1,ngi
      diff = 0.0
      if(xsec(nnufis,1,mg)>0.0) then
        fispec(mg) = fispec(mg) +(profis(mg)+xdelay(mg))/xsec(nnufis,1,mg)
      endif
      sum = sum + fispec(mg)
      spec(mg) = spec(mg)/bot
      ! check if chi gg' summed over g' = 1.0
      if(fispec(mg)/=0.0) then
        diff = abs(1.0-fispec(mg))
      endif
      if(diff>0.0001) then
        ifail = .true.
        if(iprint>0) write(nsyso,9050) mg,fispec(mg),profis(mg),xdelay(mg), &
        xsec(nnufis,1,mg),xsec(nfissn,1,mg),diff
      endif
    enddo
    if((.not.ifail).and.(iprint>0)) write(nsyso,9070)
  endif
  return
  !----
  !  formats
  !----
  9000 format ('  mf  mt',5x,'reaction',5x,'cross-check/comment',/,'  ==  ==' &
  ,5x,'========',5x,'===================',/)
  9010 format (1x,i3,1x,i3,2x,8a16,/)
  9020 format (/,' maximum array size required for matrices is',i6,/, &
  ' maximum array size required for orders is',i6)
  9030 format (/,15x,' matrices to be processed',/,16x,24 ('-'),/,1x,5a16,/,/)
  9040 format (/,30 ('-'),' matrix data ',30 ('-'),/)
  9050 format (i5,1p,e15.5,e14.5,e13.5,e13.5,e12.5,0p,f10.4)
  9060 format (/,28x,'fission spectrum chi gg dash summed over g',/,28x, &
  38 ('='),/,10x,' differences greater than 0.0001',/,10x, &
  ' ===============================',/,/, &
  '  group   chi gg dash     prompt      delayed    ','nu* &
  fission    fission   difference',/, &
  '                          matrix      matrix',/)
  9070 format (' all differences less than 0.0001',/)
  9080 format (' ***error**** maximum number of legendre orders',' exceeded', &
  /,i6,' found on tape, maximum allowed in',' this version of gecco ',i6)
  contains
    !
    integer function ntrouve(mth,mtherm,ntherm)
    integer mth,mtherm,ii,ntherm(mtherm)
    ntrouve = 0
    if(mtherm/=0) then
      do ii=1,mtherm
        if(mth==ntherm(ii)) then
          ntrouve = ii
          return
        endif
      enddo
    endif
  end function ntrouve
  end subroutine ecfil6
  !
  subroutine ecfil6t(therma,thrmin,tnmin,tnmax,torder,thdim1)
  !-----------------------------------------------------------------------
  !-----reads the thermal neutron scattering matrices on unit ning
  !----- and sets the arrays tnmin and tnmax
  !-----------------------------------------------------------------------
  !     .. scalar arguments ..
  integer thrmin,thdim1,tnmin(tmpmax,thdim1:ngi),tnmax(tmpmax,thdim1:ngi)
  ! .. array arguments ..
  real(kr) therma(thdim1:ngi,thdim1:ngi,tmpmax),torder(thdim1:ngi,thdim1:ngi, &
  maxleg,tmpmax)
  ! ..
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ll,nb,nleg,nsig,nw,ipfin,ipos1,i,j,k,ip,ipos,istart, &
  & istop,itempj,mg,mg2hi,mg2lo,ntempj
  ! .. local arrays ..
  real(kr) buff(maxnbf)
  !
  if(iprint>0) write(nsyso,900)
  ntempj = ntemp(nelp)
  itempj = itemp(nelp)
  !----
  !  read gendf
  !----
  10 call contio(ning,0,0,buff,nb,nw)
  if(mfh==6) then
    if(ntrouve(mth,mtherm,ntherm)/=0) then
      nleg = l1h
      nsig = l2h
      !
      ! thermal upscatter matrix
      if((nleg)/=npne(melas))then
        ! different number of legendre orders on thermal scattering
        ! matrix and elastic matrix
        write(hsmg,'(52hdifferent number of Legendre orders for thermal matr, &
        & 22hix and elastic matrix:,2i5)') nleg,npne(melas)
        call error('ecfil6t',hsmg,' ')
      endif
      if(iprint>0) write(nsyso,910) mfh,mth,'elastic         ', &
      & 'thermal scatteri','ng data          '
      20 ll = 1
      call listio(ning,0,0,buff,nb,nw)
      30 if(nb/=0) then
        ll = ll + nw
        call moreio(ning,0,0,buff(ll),nb,nw)
        go to 30
      endif
      ig2lo = l2h
      nw = n1h
      ig = n2h
      ig2hi = nw/(nleg*nsig)-1+ig2lo-1
      mg=ngi-ig+1
      if(mg>=itempj)then
        mg2lo=ngi-ig2hi+1
        mg2hi=ngi-ig2lo+1
        if(tnmin(ntempj,mg)==0.0) then
          tnmin(ntempj,mg) = mg2lo
        else if(tnmin(ntempj,mg) >mg2lo) then
          tnmin(ntempj,mg) = mg2lo
        endif
        if(tnmax(ntempj,mg)==0.0) then
          tnmax(ntempj,mg) = mg2hi
        else if(tnmax(ntempj,mg)<mg2hi) then
          tnmax(ntempj,mg) = mg2hi
        endif
        ipos = nleg*nsig + 7
        istart=mg2hi
        istop=mg2lo
        do j=istart,istop,-1
          therma(mg,j,ntempj)=therma(mg,j,ntempj)+buff(ipos)
          ipfin=npne(melas)-1
          do ip=1,npne(melas)-1
            torder(mg,j,ip,ntempj)=torder(mg,j,ip,ntempj)+buff(ipos+ip)
            if(ip==1) ipos1=ipos+ip
          enddo
          ipos=ipos+nleg*nsig
        enddo
      else if(mg/=1)then
        ! error - number of thermal groups inconsistent
        write(hsmg,'(37hinconsistent umber of thermal groups:,2i5)') mg,itempj
        call error('ecfil6t',hsmg,' ')
      endif
      if(n2h<ngi) go to 20
    endif
    call tosend(ning,0,0,buff)
    go to 10
  else
    call skiprz(ning,-1)
  endif
  return
  !----
  !  formats
  !----
  900 format (/,1x,30('-'),' matrix data ',30('-')//'  mf  mt',5x, &
  & 'reaction',5x,'cross-check/comment'/,'  ==  ==',5x,'========',5x, &
  & '==================='/)
  910 format (1x,i3,1x,i3,2x,8a16,/)
  contains
    !
    integer function ntrouve(mth,mtherm,ntherm)
    integer mth,mtherm,ii,ntherm(mtherm)
    ntrouve = 0
    if(mtherm/=0) then
      do ii=1,mtherm
        if(mth==ntherm(ii)) then
          ntrouve = ii
          return
        endif
      enddo
    endif
  end function ntrouve
  end subroutine ecfil6t
  !
  subroutine ecfil50
  ! .. local scalars ..
  integer nb,nw,scapt,i,ig,im,ll,mg,ngipt,nor,ntempj,jc
  character*16 reac,blank
  real(kr) totcap,diff,totela,totsum
  logical ifail,mt101
  ! .. local arrays ..
  real(kr) buff(maxnbf)
  ! ..
  if(iprint>0) write(nsyso,1070)
  !----
  !  set defaults
  !----
  blank = '   -            '
  if(iprint>0) write(nsyso,900)
  mt101 = .false.
  ifail = .false.
  ntempj = ntemp(nelp)
  npmc=0
  scapt=0
  !----
  !  zeroise arrays
  !----
  if(ntempj==1) then
    ishld(:tmpmax,:ngrmax)=0
    prosig(:maxxs,:maxord,:tmpmax,:ngrmax)=0.d0
    weight(:maxord,:tmpmax,:ngrmax)=0.d0
  endif
  !----
  !  start loop over sections on tape ning
  !----
  10 call contio(ning,0,0,buff,nb,nw)
  if(mfh==50) then
    npmc=npmc+1
    if(l2h>maxord) then
      write(hsmg,'(35hmaximum order of probability table=,i5, &
      & 18h; maximum allowed=,i5)') l2h,maxord
      call error('ecfil50',hsmg,' ')
    endif
    140 ll = 1
    call listio(ning,0,0,buff(1),nb,nw)
    150 if(nb/=0) then
      if((ll+nw)>maxnbf) then
        write(nsyso,1020)
        call error('ecfil50','work space exceeded reading subgroup data',' ')
      endif
      ll = ll + nw
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 150
    endif
    if(mth==101) mt101=.true.
    nor = l2h
    ig = n2h
    mg=ngi-ig+1
    !----
    !  write this group away to prosig
    !  nor=1 for final group
    !----
    if(nor>1) then
      ishld(ntempj,mg)=nor
      if( (mth==1).or.(mth==2).or.(mth==4).or.(mth==16).or.(mth==17).or. &
      & (mth==18).or.(mth==101).or.(.not.mt101.and.mth==102) ) then
        jc=6
        !----
        !  read weight
        !----
        if(mth==1) weight(:nor,ntempj,mg) = buff(jc+1:jc+nor)
        !----
        !  read cross-section
        !----
        jc=jc+nor
        if((mth==16).or.(mth==17)) then
          ! sumup (n,2n) and (n,3n)
          prosig(npmc,:nor,ntempj,mg) = prosig(npmc,:nor,ntempj,mg) + &
          & buff(jc+1:jc+nor)
        else
          prosig(npmc,:nor,ntempj,mg) = buff(jc+1:jc+nor)
        endif
      endif
    endif
    if(ig<ngi) go to 140
    if(mth==1) then
      reac='TOTAL'
      call search(reac,ipmc(npmc),listmi,nsmic)
      if(iprint>0) write(nsyso,910) mth,blank,blank,listmi(ipmc(npmc))
    else if(mth==2) then
      reac='ELASTIC'
      call search(reac,ipmc(npmc),listmi,nsmic)
      selas=npmc
      if(iprint>0) write(nsyso,910) mth,listmi(ipmc(npmc)),'total           '
    else if(mth==4) then
      reac='INELASTIC'
      call search(reac,ipmc(npmc),listmi,nsmic)
      sinel=npmc
      if(iprint>0) write(nsyso,910) mth,listmi(ipmc(npmc)),'total           '
    else if((mth==16).or.(mth==17)) then
      reac='N,XN'
      call search(reac,ipmc(npmc),listmi,nsmic)
      if(iprint>0) write(nsyso,910) mth,listmi(ipmc(npmc)),'total           '
    else if(mth==18) then
      reac='FISSION'
      call search(reac,ipmc(npmc),listmi,nsmic)
      if(iprint>0) write(nsyso,910) mth,listmi(ipmc(npmc)),'total           '
    else if(mth==101.or.(.not.mt101.and.mth==102)) then
      reac='CAPTURE'
      scapt = npmc
      call search(reac,ipmc(npmc),listmi,nsmic)
      if(iprint>0) write(nsyso,910) mth,listmi(ipmc(npmc)),'total           '
    else
      !----
      !  reaction not recognised
      !----
      if(iprint>0) write(nsyso,fmt=1010) mfh,mth,ntempj
      npmc = npmc - 1
    endif
    call tosend(ning,0,0,buff)
    go to 10
    !
  else
    !----
    !  check if end of mfh has been reached
    !----
    if(mfh/=0) call skiprz(ning,-1)
    go to 20
  endif
  !-----end loop over sections on tape ning
  20 continue
  !******************************************************************
  !*  compare capture calculated from the sub group data with total *
  !*  capture from summing the capture partials.                    *
  !*  replace primary capture cross sections with values from       *
  !* sub group data when present.                                   *
  !******************************************************************
  if(scapt/=0)then
    if(iprint>0) write(nsyso,1140)
    do mg = 1,ngi
     totcap = 0.0
      nor = ishld(ntempj,mg)
      if( nor>=1 ) then
        do im = 1,nor
          totcap = totcap + weight(im,ntempj,mg)*prosig(scapt,im,ntempj,mg)
        enddo
        diff = abs((xsec(ncapt,ntempj,mg)-totcap))/totcap*100.
        if((diff>1.0).and.(iprint>0)) then
          write(nsyso,1090)mg,xsec(ncapt,ntempj,mg),totcap,diff
        endif
        !----
        !  correct total
        !----
        xsec(ntotal,ntempj,mg) = xsec(ntotal,ntempj,mg)- xsec(ncapt,ntempj,mg) &
        & + totcap
        !----
        !  replace capture
        !----
        xsec(ncapt,ntempj,mg)=totcap
      endif
    enddo
  endif
  !----
  !  check subgroup elastic if thermal elastic present
  !----
  if(selas/=0.and.mtherm/=0)then
    if(iprint>0) write(nsyso,1150)
    do mg = itemp(nelp),ngi
      totela = 0.0
      nor = ishld(ntempj,mg)
      if( nor>=1 ) then
        do im = 1,nor
          totela = totela + weight(im,ntempj,mg)*prosig(selas,im,ntempj,mg)
        enddo
        diff = abs((xsec(nelas,ntempj,mg)-totela))/totela*100.
        write(nsyso,1090)mg,xsec(nelas,ntempj,mg),totela,diff
        !----
        !  correct  subgroup data
        !----
        do im=1,nor
          ! correct total
          prosig(1,im,ntempj,mg) = prosig(1,im,ntempj,mg) - &
          & prosig(selas,im,ntempj,mg) +prosig(selas,im,ntempj,mg)* &
          & xsec(nelas,ntempj,mg)/totela
          ! correct elastic
          prosig(selas,im,ntempj,mg) = prosig(selas,im,ntempj,mg)* &
          xsec(nelas,ntempj,mg)/totela
        enddo
      endif
    enddo
  endif
  !----
  !  check if totals are consistent
  !----
  if(iprint>0) write(nsyso,1100)
  if(iprint>0) write(nsyso,1120)
  do mg = 1,ngi
    totsum = 0.0
    nor = ishld(ntempj,mg)
    if( nor>=1 ) then
      do im = 1,nor
        totsum = totsum + weight(im,ntempj,mg)*prosig(1,im,ntempj,mg)
      enddo
      diff = abs((xsec(ntotal,ntempj,mg)-totsum))/totsum*100.
      if(diff>1.0)then
        ifail = .true.
        if(iprint>0) write(nsyso,1090) mg,xsec(ntotal,ntempj,mg),totsum,diff
      endif
      prosig(1,im,ntempj,mg) = prosig(1,im,ntempj,mg)*xsec(ntotal,ntempj,mg)/ &
      & totsum
    endif
  enddo
  if(iprint>0) then
    if(ifail)then
      write(nsyso,1130)
    else
      write(nsyso,1110)
    endif
  endif
  !----
  !  set up ipmc to include transport
  !----
  npmc=npmc+1
  reac='TRANSPORT'
  call search(reac,ipmc(npmc),listmi,nsmic)
  !----
  !  print reactions processed
  !----
  if(iprint>0) write(nsyso,1050)
  if(iprint>0) write(nsyso,1060)(listmi(ipmc(i)),i=1,npmc)
  return
  !
  900 format ('  mt',3x,'reaction',22x,'cross-check/comment'/'  ==',3x, &
  & '========',22x,'==================='/)
  910 format (1x,i3,2x,8a16,/)
  1010 format (/' ****** warning - mf ',i3,' mt ',i3,' temperature ',i3, &
  ' has been ignored ******')
  1020 format (/,' ***** error - work space exceeded reading ', &
  & ' sub-group data *****')
  1050 format (/,15x,'Sub-group data to be processed',/,15x,30('-'),/)
  1060 format (1x,5a16)
  1070 format (/,30('-'),' sub-group data ',30('-'),/)
  1090 format (1x,i4,2x,1p,e13.5,3x,e13.5,2x,0p,f10.4)
  1100 format (/,' TOTAL has been calculated from the sub-group data', &
   &    ' and compared ','with the TOTAL calculated from file 3',/)
  1110 format (' No differences greater than 1% were found')
  1120 format (/,' NOTE: The sub-group total data will be scaled to be', &
  &        ' consistent with the primary total cross-section.',/,/, &
  &        ' Group File 3  TOTAL   Sub-group TOTAL    diff(%)',/, &
  &        ' ----- -------------   ---------------    -------')
  1130 format (/,1x,6('*** ERROR ***'),/, &
  &        ' Differences greater than 1% found  *****', &
  &        /,1x,6('*** ERROR ***'))
  1140 format (/,' CAPTURE has been calculated from the sub-group data', &
  &       ' and compared with the CAPTURE calculated from ',/, &
  &       ' file 3. The sub group value is written to the ECCO', &
  &       ' library when present',/,/, &
  &       ' Group File 3  CAPTURE  Sub-group CAPTURE  diff(%)',/, &
  &       ' ----- -------------   ---------------    -------')
  1150 format (/,' ELASTIC has been calculated from the sub-group data', &
  &       ' and compared with the ELASTIC calculated from ',/, &
  &       ' file 3. The sub group value is adjusted when', &
  &       ' thermal scattering data is present',/,/, &
  &       ' Group Library ELASTIC  Sub-group ELASTIC  diff(%)',/, &
  &       ' ----- -------------   ---------------    -------')
  end subroutine ecfil50
  !
  subroutine ecfad1(nsmxj,ing,nreci,nsig,lenrec,matrx2,ngi,matrix,buff, &
  & bufsiz,lunit,lnrec1,inmin,inmax,maxm)
  !
  ! add matrices - nleg = 1 (no higher order matrices present)
  !
  ! no higher order data in reactions added so far
  ! .. scalar arguments ..
  integer bufsiz,lenrec,lunit,ngi,nreci,nsig,nsmxj,maxm
  ! .. array arguments ..
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),matrx2(msize)
  integer ing(nreci),inmax(maxmtx,ngrmax),inmin(maxmtx,ngrmax), &
  lnrec1(ngblk*5,2,maxmtx)
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ilower,imgmax,imgmin,ipos,j,ll,mg,mg2hi, &
  mg2lo,mpos,mpos2,mrec1,n,m,nb,nrec1,nw,minlo,maxhi,ijk,is,istop, &
  lunit1,lunit3,mend,mglast,mtop
  !
  lunit1 = lunit
  if(lunit<90)then
    lunit3 = lunit1 + 10
  else
    lunit3 = lunit1 - 10
  endif
  nrec1 = 0
  mpos = 0
  mpos2 = 0
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  mglast = ngi + 1
  !----
  !  read in first block
  !----
  mrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(mpos>0)then
    call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
  else
    matrix(:msize)=0.0
  endif
  mpos = 0
  10 continue
  ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 continue
  if(nb/=0) then
    ll = ll + nw
    if(ll>bufsiz) then
      call error('ecfad1','** error reading gendf, workspace exceeded **',' ')
    else
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
  endif
  nw = n1h
  ig = n2h
  ig2lo = l2h
  ig2hi = nw/nsig - 1 + ig2lo - 1
  mg = ngi - ig + 1
  mg2lo = ngi - ig2hi + 1
  mg2hi = ngi - ig2lo + 1
  ipos = nsig + 7
  if(inmin(nsmxj,mg)==0)then
    maxhi = 0
    minlo = ngi + 1
  else
    minlo = ngi - inmax(nsmxj,mg) + 1
    maxhi = ngi - inmin(nsmxj,mg) + 1
  endif
  !----
  !  find correct group block for current group
  !----
  30 continue
  if(mg<imgmin) then
    !----
    !  group is not in current group block
    !  check if old buffer has been copied to new buffer
    !----
    if(mpos<lnrec1(n,2,nsmxj))then
      !----
      !  copy rest of old buffer
      !----
      mend = lnrec1(n,2,nsmxj)
      do m = mpos+1,mend
        mpos2 = mpos2 + 1
        matrx2(mpos2) = matrix(m)
      enddo
    endif
    !----
    ! write new block to buffer
    !----
    if(mpos2>0)then
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos2
      maxm = max(mpos2,maxm)
      if(maxm>msize) then
        write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxm,msize
        call error('ecfad1',hsmg,' ')
      endif
      call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
      mpos2 = 0
    endif
    !----
    !  read in next block
    !----
    n = n - 1
    imgmax = imgmin - 1
    imgmin = imgmax - ing(n) + 1
    mrec1 = lnrec1(n,1,nsmxj) - 1
    mpos = lnrec1(n,2,nsmxj)
    if(mpos>0) then
      call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
      mpos = 0
    else
      matrix(:msize)=0.0
    endif
    go to 30
  endif
  !----
  !  check for dummy last group
  !----
  if((ig/=ngi).or.(ig2lo/=ig2hi).or.(buff(ipos)/=0.0))then
    !----
    !  if threshold reaction and mg is not the first group in the group block
    !  or gaps in incident energy groups
    !  check whether the input buffer contains data for previous groups
    !----
    if(mg/=imgmax) then
      if(mg<mglast-1) then
        mtop = min(mglast-1,imgmax)
        do is = mg+1,mtop
          if(inmin(nsmxj,is)/=0) then
            do j = inmin(nsmxj,is),inmax(nsmxj,is)
              mpos = mpos + 1
              mpos2 = mpos2 + 1
              matrx2(mpos2) = matrix(mpos)
            enddo
          endif
        enddo
      endif
    endif
    if(inmin(nsmxj,mg)==0) then
      !----
      !  no data on buffer , add new data
      !----
      inmin(nsmxj,mg) = mg2lo
      inmax(nsmxj,mg) = mg2hi
      do j = mg2lo,mg2hi
        mpos2 = mpos2 + 1
        matrx2(mpos2) = buff(ipos)
        ipos = ipos + nsig
      enddo
    else
      if(minlo<ig2lo) then
        !----
        !  more data on buffer so copy first
        !----
        ilower = minlo
        !----
        !  check for gaps
        !----
        if(ig2lo>maxhi)then
          !----
          !  copy data from buffer and and zeroise gap
          !----
          do j=ilower,maxhi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
          enddo
          do j=maxhi+1,ig2lo-1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = 0.0
          enddo
          maxhi = ig2lo-1
        else
          !----
          !  copy from buffer
          !----
          do j = ilower,ig2lo - 1
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
          enddo
        endif
        ilower = ig2lo
      else if(minlo>ig2lo) then
        !----
        !  more data on gendf so copy first
        !----
        istop = minlo - 1
        if(ig2hi<minlo) istop = ig2hi
        ilower = ig2lo
        do j = ilower,istop
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          ipos = ipos + nsig
        enddo
       ilower = istop + 1
      else
        !----
        !  lower groups equal
        !----
        ilower = ig2lo
      endif
      !----
      !  check for gaps in data
      !----
      if(ilower<minlo) then
        do j = ilower,minlo-1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = 0.0
        enddo
        ilower = minlo
      endif
      if(maxhi<=ig2hi) then
        !----
        !  add data
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos) + buff(ipos)
          ipos = ipos + nsig
        enddo
        ilower = maxhi + 1
        if(ig2hi>maxhi) then
          do j = ilower,ig2hi
            mpos2 = mpos2 + 1
            matrx2(mpos2) = buff(ipos)
            ipos = ipos + nsig
          enddo
          ilower = ig2hi + 1
        endif
      else
        if(ig2hi>=ilower) then
          do j = ilower,ig2hi
           mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos) + buff(ipos)
            ipos = ipos + nsig
          enddo
          ilower = ig2hi+ 1
        endif
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
        enddo
      endif
    endif
    inmin(nsmxj,mg) = min(inmin(nsmxj,mg),mg2lo)
    inmax(nsmxj,mg) = max(inmax(nsmxj,mg),mg2hi)
    mglast = mg
    if(ig<ngi) go to 10
  endif
  !
  ! copy any remaining data
  !
  if(mpos+1<lnrec1(n,2,nsmxj))then
    istop = lnrec1(n,2,nsmxj)-1
    do ijk=mpos,istop
      mpos=mpos+1
      mpos2=mpos2+1
      matrx2(mpos2) = matrix(mpos)
    enddo
  endif
  maxm = max(mpos2,maxm)
  if(maxm>msize) then
    write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
    & i9,13h actual size=,i8)') mth,maxm,msize
    call error('ecfad1',hsmg,' ')
  endif
  !----
  !  write last block away to buffer
  !----
  if(mpos2>0)then
    maxm = max(mpos2,maxm)
    if(maxm>msize) then
      write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxm,msize
      call error('ecfad1',hsmg,' ')
    endif
    lnrec1(n,1,nsmxj) = nrec1 + 1
    lnrec1(n,2,nsmxj) = mpos2
    call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
  endif
  lunit = lunit3
  end subroutine ecfad1
  !
  subroutine ecfad2(nsmxj,ing,nreci,nsig,lenrec,matrx2,ngi,matrix,buff,bufsiz, &
  & lunit,lnrec1,lnrec2,inmin,inmax,porder,pordr2,nleg,maxm,maxp)
  !
  ! add matrices - nleg = npne(nmsxj) > 1  (higher order matrices present)
  !
  ! buffer contains data for higher order data but this reaction does not
  !
  ! .. scalar arguments ..
  integer bufsiz,lenrec,lunit,ngi,nleg,nreci,nsig,nsmxj,maxm,maxp
  ! .. array arguments ..
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  integer ing(nreci),inmax(maxmtx,ngrmax),inmin(maxmtx,ngrmax), &
  lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ilower,imgmax,imgmin,ipos,j,ll,lunit2,lunit4,mg, &
  mg2hi,mg2lo,mpos,mpos2,mrec1,mrec2,n,nb,nrec1,nrec2,nw,ppos,ppos2,minlo, &
  maxhi,p,pend,ijk,ip,is,istop,lunit1,lunit3,m,mend,mglast,mtop
  !
  nrec1 = 0
  nrec2 = 0
  mpos = 0
  ppos = 0
  mpos2 = 0
  ppos2 = 0
  lunit1 = lunit
  lunit2 = lunit1 + 5
  if(lunit<90)then
    lunit3 = lunit1 + 10
    lunit4 = lunit2 + 10
  else
    lunit3 = lunit1 - 10
    lunit4 = lunit2 - 10
  endif
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  mglast = ngi + 1
  !----
  !  read in first block
  !----
  mrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(mpos>0)then
    call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
    mrec2 = lnrec2(n,1,nsmxj) - 1
    ppos = lnrec2(n,2,nsmxj)
    if(ppos>0)then
      call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
      ppos = 0
    else
      porder(:psize)=0.0
    endif
  else
    matrix(:msize)=0.0
  endif
  mpos = 0
  10 continue
  ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 continue
  if(nb/=0) then
    ll = ll + nw
    if(ll>bufsiz) then
      call error('ecfad2',' ** error reading gendf, workspace exceeded',' ')
    else
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
  endif
  nw = n1h
  ig = n2h
  ig2lo = l2h
  ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
  mg = ngi - ig + 1
  mg2lo = ngi - ig2hi + 1
  mg2hi = ngi - ig2lo + 1
  ipos = nleg*nsig + 7
  if(inmin(nsmxj,mg)==0)then
    maxhi = 0
    minlo = ngi + 1
  else
    minlo = ngi - inmax(nsmxj,mg) + 1
    maxhi = ngi - inmin(nsmxj,mg) + 1
  endif
  !----
  !  find correct group block for current group
  !----
  30 if(mg<imgmin) then
    !
    ! group is not in current group block
    !
    ! check if old buffer has been copied to new buffer
    !
    if(mpos<lnrec1(n,2,nsmxj))then
      !
      ! copy rest of old buffer
      !
      mend = lnrec1(n,2,nsmxj)
      do m = mpos+1,mend
        mpos2 = mpos2 + 1
        matrx2(mpos2) = matrix(m)
      enddo
      if(ppos<lnrec2(n,2,nsmxj))then
        pend = lnrec2(n,2,nsmxj)
        do p = ppos+1,pend
          ppos2 = ppos2 + 1
          pordr2(ppos2) = porder(p)
        enddo
      endif
    endif
    !----
    !  write new block to buffer
    !----
    if(mpos2>0)then
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos2
      maxm = max(mpos2,maxm)
      if(maxm>msize) then
        write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxm,msize
        call error('ecfad2',hsmg,' ')
      endif
      call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
      lnrec2(n,1,nsmxj) = nrec2 + 1
      lnrec2(n,2,nsmxj) = ppos2
      maxp = max(ppos2,maxp)
      if(maxp>psize) then
        write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxp,psize
        call error('ecfad2',hsmg,' ')
      endif
      call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
      mpos2 = 0
      ppos2 = 0
    endif
    !----
    !  read in next block
    !----
    n = n - 1
    imgmax = imgmin - 1
    imgmin = imgmax - ing(n) + 1
    mrec1 = lnrec1(n,1,nsmxj) - 1
    mpos = lnrec1(n,2,nsmxj)
    if(mpos>0) then
      call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
      mpos = 0
      mrec2 = lnrec2(n,1,nsmxj) - 1
      ppos = lnrec2(n,2,nsmxj)
      if(ppos>0)then
        call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
        ppos = 0
      else
        porder(:psize)=0.0
      endif
    else
      matrix(:msize)=0.0
    endif
    go to 30
  endif
  !----
  !  check for dummy last group
  !----
  if((ig/=ngi).or.(ig2lo/=ig2hi).or.(buff(ipos)/=0.0))then
    !----
    !  if threshold reaction and mg is not the first group in the group block
    !  or gaps in incident energy groups
    !  check whether the input buffer contains data for previous groups
    !----
    if(mg/=imgmax) then
      if(mg<mglast-1) then
        mtop = min(mglast-1,imgmax)
        do is = mg+1,mtop
          if(inmin(nsmxj,is)/=0) then
            do j = inmin(nsmxj,is),inmax(nsmxj,is)
              mpos = mpos + 1
              mpos2 = mpos2 + 1
              matrx2(mpos2) = matrix(mpos)
              do ip = 1,nleg - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            enddo
          endif
        enddo
      endif
    endif
    if(inmin(nsmxj,mg)==0) then
      !----
      !  no data on buffer , add new data
      !----
      inmin(nsmxj,mg) = mg2lo
      inmax(nsmxj,mg) = mg2hi
      do j = mg2lo,mg2hi
        mpos2 = mpos2 + 1
        matrx2(mpos2) = buff(ipos)
        do ip = 1,nleg - 1
          ppos2 = ppos2 + 1
          pordr2(ppos2) = buff(ipos+ip)
        enddo
        ipos = ipos + nleg*nsig
      enddo
    else
      if(minlo<ig2lo) then
        !----
        !  more data on buffer so copy first
        !----
        ilower = minlo
        !----
        !  check for gaps
        !----
        if(ig2lo>maxhi)then
          !----
          !  copy data from buffer and and zeroise gap
          !----
          do j=ilower,maxhi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            do ip = 1,nleg - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos)
            enddo
          enddo
          do j=maxhi+1,ig2lo-1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = 0.0
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = 0.0
            enddo
          enddo
          maxhi = ig2lo-1
        else
          do j = ilower,ig2lo - 1
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            do ip = 1,nleg - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos)
            enddo
          enddo
        endif
        ilower = ig2lo
      else if(minlo>ig2lo) then
        !----
        !  more data on gendf so copy first
        !----
        ilower = ig2lo
        istop = minlo -1
        if(ig2hi<minlo) istop = ig2hi
        do j = ilower,istop
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = istop + 1
      else
        !----
        !  lower groups equal
        !----
        ilower = ig2lo
      endif
      !----
      !  check for gaps in data
      !----
      if(ilower<minlo) then
        do j = ilower,minlo-1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = 0.0
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
          enddo
        enddo
        ilower = minlo
      endif
      if(maxhi<=ig2hi) then
        !----
        !  add data
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos) + buff(ipos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos) + buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = maxhi + 1
        if(ig2hi>maxhi) then
          do j = ilower,ig2hi
            mpos2 = mpos2 + 1
            matrx2(mpos2) = buff(ipos)
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = buff(ipos+ip)
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
      else
        if(ig2hi>=ilower)then
          do j = ilower,ig2hi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos) + buff(ipos)
            do ip = 1,nleg - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos) + buff(ipos+ip)
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos)
          enddo
        enddo
        ilower = maxhi + 1
      endif
    endif
    inmin(nsmxj,mg) = min(inmin(nsmxj,mg),mg2lo)
    inmax(nsmxj,mg) = max(inmax(nsmxj,mg),mg2hi)
    mglast = mg
    if(ig<ngi) go to 10
  endif
  !----
  !  copy any remaining data
  !----
  if(mpos+1<lnrec1(n,2,nsmxj))then
    istop = lnrec1(n,2,nsmxj) -1
    do ijk=mpos,istop
      mpos=mpos+1
      mpos2=mpos2+1
      matrx2(mpos2) = matrix(mpos)
    enddo
  endif
  if(ppos+1<lnrec2(n,2,nsmxj))then
    istop = lnrec2(n,2,nsmxj) -1
    do ijk=ppos,istop
      ppos=ppos+1
      ppos2=ppos2+1
      pordr2(ppos2) = porder(ppos)
    enddo
  endif
  if(mpos2>0)then
    maxm = max(mpos2,maxm)
    if(maxm>msize) then
      write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxm,msize
      call error('ecfad2',hsmg,' ')
    endif
    lnrec1(n,1,nsmxj) = nrec1 + 1
    lnrec1(n,2,nsmxj) = mpos2
    call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
    maxp = max(ppos2,maxp)
    if(maxp>psize) then
      write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxp,psize
      call error('ecfad2',hsmg,' ')
    endif
    lnrec2(n,1,nsmxj) = nrec2 + 1
    lnrec2(n,2,nsmxj) = ppos2
    call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
  endif
  lunit = lunit3
  end subroutine ecfad2
  !
  subroutine ecfad3(nsmxj,ing,nreci,nsig,lenrec,matrx2,ngi,matrix,buff, &
  & bufsiz,lunit,lnrec1,lnrec2,inmin,inmax,porder,pordr2,nleg,npne1,maxm,maxp)
  !
  ! add matrices - nleg = 1 and  npne1  > 1. (higher orders)
  !
  !  higher order data on buffer but not present for this reaction
  ! .. scalar arguments ..
  integer bufsiz,lenrec,lunit,ngi,nleg,npne1,nreci,nsig,nsmxj,maxm,maxp

  ! .. array arguments ..
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  integer ing(nreci),inmax(maxmtx,ngrmax),inmin(maxmtx,ngrmax), &
  lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ilower,imgmax,imgmin,ipos,j,ll,lunit2,lunit4,mg, &
  mg2hi,mg2lo,mpos,mpos2,mrec1,mrec2,n,nb,nrec1,nrec2,nw,ppos,ppos2, &
  minlo,maxhi,totzer,p,pend,ijk,ip,is,istop,lunit1,lunit3,m,mend,mglast,mtop
  !
  totzer = 0
  nrec1 = 0
  nrec2 = 0
  mpos = 0
  ppos = 0
  mpos2 = 0
  ppos2 = 0
  lunit1 = lunit
  lunit2 = lunit1 + 5
  if(lunit<90)then
    lunit3 = lunit1 + 10
    lunit4 = lunit2 + 10
  else
    lunit3 = lunit1 - 10
    lunit4 = lunit2 - 10
  endif
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  mglast = ngi + 1
  !
  !  read in first block
  !
  mrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(mpos>0)then
    call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
    mrec2 = lnrec2(n,1,nsmxj) - 1
    ppos = lnrec2(n,2,nsmxj)
    if(ppos>0)then
      call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
      ppos = 0
    else
      porder(:psize)=0.0
    endif
  else
    matrix(:msize)=0.0
  endif
  mpos = 0
  10 continue
  ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 continue
  if(nb/=0) then
    ll = ll + nw
    if(ll>bufsiz) then
      call error('ecfad3',' ** error reading gendf, workspace exceeded',' ')
    else
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
  endif
  nw = n1h
  ig = n2h
  ig2lo = l2h
  ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
  mg = ngi - ig + 1
  mg2lo = ngi - ig2hi + 1
  mg2hi = ngi - ig2lo + 1
  ipos = nleg*nsig + 7
  if(inmin(nsmxj,mg)==0)then
    maxhi = 0
    minlo = ngi + 1
  else
    minlo = ngi - inmax(nsmxj,mg) + 1
    maxhi = ngi - inmin(nsmxj,mg) + 1
  endif
  30 continue
  if(mg<imgmin) then
    !----
    !  group is not in current group block
    !  check if old buffer has been copied to new buffer
    !----
    if(mpos<lnrec1(n,2,nsmxj))then
      !----
      !  copy rest of old buffer
      !----
      mend = lnrec1(n,2,nsmxj)
      do m = mpos+1,mend
        mpos2 = mpos2 + 1
        matrx2(mpos2) = matrix(m)
      enddo
      if(ppos<lnrec2(n,2,nsmxj))then
        pend = lnrec2(n,2,nsmxj)
        do p = ppos+1,pend
          ppos2 = ppos2 + 1
          pordr2(ppos2) = porder(p)
        enddo
      endif
    endif
    !----
    !  write new block to buffer
    !----
    if(mpos2>0)then
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos2
      maxm = max(mpos2,maxm)
      if(maxm>msize) then
        write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=,i9, &
        13h actual size=,i8)') mth,maxm,msize
        call error('ecfad3',hsmg,' ')
      endif
      call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
      lnrec2(n,1,nsmxj) = nrec2 + 1
      lnrec2(n,2,nsmxj) = ppos2
      maxp = max(ppos2,maxp)
      if(maxp>psize) then
        write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=,i9, &
        13h actual size=,i8)') mth,maxp,psize
        call error('ecfad3',hsmg,' ')
      endif
      call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
      mpos2 = 0
      ppos2 = 0
    endif
    !----
    !  read in next block
    !----
    n = n - 1
    imgmax = imgmin - 1
    imgmin = imgmax - ing(n) + 1
    mrec1 = lnrec1(n,1,nsmxj) - 1
    mpos = lnrec1(n,2,nsmxj)
    if(mpos>0) then
      call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
      mpos = 0
      mrec2 = lnrec2(n,1,nsmxj) - 1
      ppos = lnrec2(n,2,nsmxj)
      if(ppos>0)then
        call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
        ppos = 0
      else
        porder(:psize)=0.0
      endif
    else
      matrix(:msize)=0.0
    endif
    go to 30
  endif
  !----
  !  check for dummy last group
  !----
  if((ig/=ngi).or.(ig2lo/=ig2hi).or.(buff(ipos)/=0.0))then
    !----
    !  if threshold reaction and mg is not the first group in the group block
    !  or gaps in incident energy groups
    !  check whether the input buffer contains data for previous groups
    !----
    if(mg/=imgmax) then
      if(mg<mglast-1) then
        mtop = min(mglast-1,imgmax)
        do is = mg+1,mtop
          if(inmin(nsmxj,is)/=0) then
            do j = inmin(nsmxj,is),inmax(nsmxj,is)
              mpos = mpos + 1
              mpos2 = mpos2 + 1
              matrx2(mpos2) = matrix(mpos)
              do ip = 1,npne1 - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            enddo
          endif
        enddo
      endif
    endif
    if(inmin(nsmxj,mg)==0) then
      !----
      !  no data on buffer , add new data
      !----
      inmin(nsmxj,mg) = mg2lo
      inmax(nsmxj,mg) = mg2hi
      do j = mg2lo,mg2hi
        mpos2 = mpos2 + 1
        matrx2(mpos2) = buff(ipos)
        do ip = 1,npne1 - 1
          ppos2 = ppos2 + 1
          pordr2(ppos2) = 0.0
          totzer = totzer + 1
        enddo
        ipos = ipos + nleg*nsig
      enddo
    else
      if(minlo<ig2lo) then
        !----
        !  more data on buffer so copy first
        !----
        ilower = minlo
        !----
        !  check for gaps
        !----
        if(ig2lo>maxhi)then
          !----
          !  copy data from buffer and and zeroise gap
          !----
          do j=ilower,maxhi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            if(npne1>1) then
              do ip = 1,npne1 - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            endif
          enddo
          do j=maxhi+1,ig2lo-1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = 0.0
            if(npne1>1) then
              do ip = 1,npne1 - 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = 0.0
              enddo
            endif
          enddo
          maxhi = ig2lo-1
        else
          !----
          !  copy data from buffer
          !----
          do j = ilower,ig2lo - 1
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            if(npne1>1) then
              do ip = 1,npne1 - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            endif
          enddo
        endif
        ilower = ig2lo
      else if(minlo>ig2lo) then
        !----
        !  more data on gendf so copy first
        !----
        ilower = ig2lo
        istop = minlo - 1
        if( ig2hi<minlo) then
          istop = ig2hi
        endif
        do j = ilower,istop
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          do ip = 1,npne1 - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
            totzer = totzer + 1
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = istop + 1
      else
        !----
        !  lower groups equal
        !----
        ilower = ig2lo
      endif
      !----
      !  check for gaps in data
      !----
      if(ilower<minlo) then
        do j = ilower,minlo-1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = 0.0
          do ip = 1,npne1-1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
          enddo
        enddo
        ilower = minlo
      endif
      if(maxhi<=ig2hi) then
        !----
        !  add data
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos) + buff(ipos)
          do ip = 1,npne1 - 1
            ppos2 = ppos2 + 1
            ppos = ppos + 1
            pordr2(ppos2) = porder(ppos)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = maxhi + 1
        if(ig2hi>maxhi) then
          !----
          !  more data on gendf so copy
          !----
          do j = ilower,ig2hi
            mpos2 = mpos2 + 1
            matrx2(mpos2) = buff(ipos)
            do ip = 1,npne1 - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = 0.0
              totzer = totzer + 1
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
      else
        if(ig2hi>=ilower)then
          !----
          !*    add data *
          !----
          do j = ilower,ig2hi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos) + buff(ipos)
            do ip = 1,npne1 - 1
              ppos2 = ppos2 + 1
              ppos = ppos + 1
              pordr2(ppos2) = porder(ppos)
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
        !----
        !  more data on buffer so copy
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
          do ip = 1,npne1 - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos)
          enddo
        enddo
        ilower = maxhi + 1
      endif
    endif
    inmin(nsmxj,mg) = min(inmin(nsmxj,mg),mg2lo)
    inmax(nsmxj,mg) = max(inmax(nsmxj,mg),mg2hi)
    mglast = mg
    if(ig<ngi) go to 10
  else
    !----
    !  dummy last group
    !----
  endif
  !----
  !  copy any remaining data
  !----
  if(mpos+1<lnrec1(n,2,nsmxj))then
    istop = lnrec1(n,2,nsmxj) -1
    do ijk=mpos,istop
      mpos=mpos+1
      mpos2=mpos2+1
      matrx2(mpos2) = matrix(mpos)
    enddo
  endif
  if(ppos+1<lnrec2(n,2,nsmxj))then
    istop = lnrec2(n,2,nsmxj) -1
    do ijk=ppos,istop
      ppos=ppos+1
      ppos2=ppos2+1
      pordr2(ppos2) = porder(ppos)
    enddo
  endif
  if(mpos2>0)then
    maxm = max(maxm,mpos2)
    if(maxm>msize) then
      write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=,i9, &
      13h actual size=,i8)') mth,maxm,msize
      call error('ecfad3',hsmg,' ')
    endif
    lnrec1(n,1,nsmxj) = nrec1 + 1
    lnrec1(n,2,nsmxj) = mpos2
    call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
    maxp = max(maxp,ppos2)
    if(maxp>psize) then
      write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=,i9, &
      13h actual size=,i8)') mth,maxp,psize
      call error('ecfad2',hsmg,' ')
    endif
    lnrec2(n,1,nsmxj) = nrec2 + 1
    lnrec2(n,2,nsmxj) = ppos2
    call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
  endif
  lunit = lunit3
  if(iprint>0) write(nsyso,"(' number of zero in inelastic matrix is ',i8)") &
  & totzer
  end subroutine ecfad3
  !
  subroutine ecfad4(nsmxj,ing,nreci,nsig,lenrec,matrx2,ngi,matrix,buff, &
  & bufsiz,lunit,lnrec1,lnrec2,inmin,inmax,porder,pordr2,nleg,npne1,maxm,maxp)
  !
  ! add matrices - nleg > 1 and npne1 = 1. (higher orders)
  !
  !  no higher order data on buffer but present for this reaction
  ! .. scalar arguments ..
  integer bufsiz,lenrec,lunit,ngi,nleg,npne1,nreci,nsig,nsmxj,maxm,maxp
  ! .. array arguments ..
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  integer ing(nreci),inmax(maxmtx,ngrmax),inmin(maxmtx,ngrmax), &
  lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  !  .. local scalars ..
  integer ig,ig2hi,ig2lo,ilower,imgmax,imgmin,ipos,j,ll,lunit2,lunit4, &
  mg,mg2hi,mg2lo,mpos,mpos2,mrec1,mrec2,n,nb,nrec1,nrec2,nw,ppos2,minlo, &
  maxhi,totzer,ijk,ip,is,istop,lunit1,lunit3,m,mend,mglast,mtop
  !
  mglast = 0
  totzer = 0
  nrec1 = 0
  nrec2 = 0
  mpos = 0
  mpos2 = 0
  ppos2 = 0
  lunit1 = lunit
  lunit2 = lunit1 + 5
  if(lunit<90)then
    lunit3 = lunit1 + 10
    lunit4 = lunit2 + 10
  else
    lunit3 = lunit1 - 10
    lunit4 = lunit2 - 10
  endif
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  !----
  !  read in first block
  !----
  mrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(mpos>0)then
    call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
    mrec2 = lnrec2(n,1,nsmxj) - 1
    pordr2(:psize)=0.0
  else
    matrix(:msize)=0.0
    pordr2(:psize)=0.0
  endif
  mpos = 0
  10 continue
  ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 continue
  if(nb/=0) then
    ll = ll + nw
    if(ll>bufsiz) then
      call error('ecfad4',' ** error reading gendf, workspace exceeded',' ')
    else
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
  endif
  nw = n1h
  ig = n2h
  ig2lo = l2h
  ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
  ipos = nleg*nsig + 7
  mg = ngi - ig + 1
  mg2lo = ngi - ig2hi + 1
  mg2hi = ngi - ig2lo + 1
  if(inmin(nsmxj,mg)==0)then
    maxhi = 0
    minlo = ngi + 1
  else
    minlo = ngi - inmax(nsmxj,mg) + 1
    maxhi = ngi - inmin(nsmxj,mg) + 1
  endif
  30 continue
  if(mg<imgmin) then
    !----
    !  group is not in current group block
    !  check if old buffer has been copied to new buffer
    !----
    if(mpos<lnrec1(n,2,nsmxj))then
      !----
      !  copy rest of old buffer
      !----
      mend = lnrec1(n,2,nsmxj)
      do m = mpos+1,mend
        mpos2 = mpos2 + 1
        matrx2(mpos2) = matrix(m)
      enddo
      !----
      !  no higher order data on buffer
      !----
      ppos2 = mpos2*(nleg-1)
    endif
    !----
    !  write new block to buffer
    !----
    if(mpos2>0)then
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos2
      maxm = max(mpos2,maxm)
      if(maxm>msize) then
        write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxm,msize
        call error('ecfad4',hsmg,' ')
      endif
      call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
      lnrec2(n,1,nsmxj) = nrec2 + 1
      lnrec2(n,2,nsmxj) = ppos2
      maxp = max(ppos2,maxp)
      if(maxp>psize) then
        write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxp,psize
        call error('ecfad4',hsmg,' ')
      endif
      call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
      mpos2 = 0
      ppos2 = 0
    endif
    !----
    !  read in next block
    !----
    n = n - 1
    imgmax = imgmin - 1
    imgmin = imgmax - ing(n) + 1
    mrec1 = lnrec1(n,1,nsmxj) - 1
    mpos = lnrec1(n,2,nsmxj)
    if(mpos>0) then
      call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
      mpos = 0
      mrec2 = lnrec2(n,1,nsmxj) - 1
      pordr2(:psize)=0.0
    else
      matrix(:msize)=0.0
    endif
    go to 30
  endif
  !----
  !  check for dummy last group
  !----
  if((ig/=ngi).or.(ig2lo/=ig2hi).or.(buff(ipos)/=0.0))then
    !----
    !  if threshold reaction and mg is not the first group in the group block
    !  or gaps in incident energy groups
    !  check whether the input buffer contains data for previous groups
    !----
    if(mg/=imgmax) then
      if(mg<mglast-1) then
        mtop = min(mglast-1,imgmax)
        do is = mg+1,mtop
          if(inmin(nsmxj,is)/=0) then
            do j = inmin(nsmxj,is),inmax(nsmxj,is)
              mpos = mpos + 1
              mpos2 = mpos2 + 1
              matrx2(mpos2) = matrix(mpos)
            enddo
          endif
        enddo
        ppos2 = mpos2 * (nleg-1)
      endif
    endif
    if(inmin(nsmxj,mg)==0) then
      !----
      !  no data on buffer , add new data
      !----
      inmin(nsmxj,mg) = mg2lo
      inmax(nsmxj,mg) = mg2hi
      do j = mg2lo,mg2hi
        mpos2 = mpos2 + 1
        matrx2(mpos2) = buff(ipos)
        do ip = 1,nleg - 1
          ppos2 = ppos2 + 1
          pordr2(ppos2) = buff(ipos+ip)
        enddo
        ipos = ipos + nleg*nsig
      enddo
    else
      if(minlo<ig2lo) then
        !----
        !  more data on buffer so copy first
        !----
        ilower = minlo
        !
        ! check for gaps
        !
        if(ig2lo>maxhi)then
          !----
          !  copy data from buffer and and zeroise gap
          !----
          do j=ilower,maxhi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
          enddo
          do j=maxhi+1,ig2lo-1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = 0.0
          enddo
          maxhi = ig2lo-1
          ppos2 = mpos2 * (nleg-1)
        else
          do j = ilower,ig2lo - 1
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
          enddo
          ppos2 = mpos2 * (nleg-1)
        endif
        ilower = ig2lo
      else if(minlo>ig2lo) then
        !----
        !  more data on gendf so copy first
        !----
        ilower = ig2lo
        istop = minlo - 1
        if( ig2hi<minlo) then
          istop = ig2hi
        endif
        do j = ilower,istop
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = istop + 1
      else
        !----
        !  lower groups equal
        !----
        ilower = ig2lo
      endif
      !----
      !  check for gaps in data
      !----
      if(ilower<minlo) then
        do j = ilower,minlo-1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = 0.0
          do ip = 1,nleg-1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
          enddo
        enddo
        ilower = minlo
      endif
      if(maxhi<=ig2hi) then
        !----
        !  add data
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos) + buff(ipos)
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = maxhi + 1
        if(ig2hi>maxhi) then
          !----
          !  more data on gendf*
          !----
          do j = ilower,ig2hi
            mpos2 = mpos2 + 1
            matrx2(mpos2) = buff(ipos)
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = buff(ipos+ip)
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
      else
        if(ig2hi>=ilower)then
          do j = ilower,ig2hi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos) + buff(ipos)
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = buff(ipos+ip)
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
        enddo
        ppos2 = mpos2 * (nleg -1)
        ilower = maxhi + 1
      endif
    endif
    inmin(nsmxj,mg) = min(inmin(nsmxj,mg),mg2lo)
    inmax(nsmxj,mg) = max(inmax(nsmxj,mg),mg2hi)
    mglast = mg
    if(ig<ngi) go to 10
  else
    !----
    !  dummy last group
    !----
  endif
  !----
  !  copy any remaining data
  !----
  if(mpos+1<lnrec1(n,2,nsmxj))then
    istop = lnrec1(n,2,nsmxj) -1
    do ijk=mpos,istop
      mpos=mpos+1
      mpos2=mpos2+1
      matrx2(mpos2) = matrix(mpos)
    enddo
  endif
  ppos2 = mpos2 * (nleg -1)
  if(mpos2>0)then
    maxm = max(maxm,mpos2)
    if(maxm>msize) then
      write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxm,msize
      call error('ecfad4',hsmg,' ')
    endif
    lnrec1(n,1,nsmxj) = nrec1 + 1
    lnrec1(n,2,nsmxj) = mpos2
    call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
    maxp = max(maxp,ppos2)
    if(maxp>psize) then
      write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxp,psize
      call error('ecfad4',hsmg,' ')
    endif
    lnrec2(n,1,nsmxj) = nrec2 + 1
    lnrec2(n,2,nsmxj) = ppos2
    call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
  endif
  lunit = lunit3
  if(iprint>0) write(nsyso,"(' number of zero in inelastic matrix is ',i8)") &
  & totzer
  end subroutine ecfad4
  !
  subroutine ecfad5(nsmxj,ing,nreci,nsig,lenrec,matrx2,ngi,matrix,buff,bufsiz, &
  lunit,lnrec1,lnrec2,inmin,inmax,porder,pordr2,nleg,npne1,maxm,maxp)
  !
  ! add matrices - nleg /= npne1 .and > 1. (higher orders)
  !
  ! higher order data present on buffer and reaction has different number
  ! of orders
  !  npne1 = number of orders on current buffer
  !  nleg = number of orders in reaction being processed
  ! .. scalar arguments ..
  integer bufsiz,lenrec,lunit,ngi,nleg,npne1,nreci,nsig,nsmxj,maxm,maxp
  ! .. array arguments
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  integer ing(nreci),inmax(maxmtx,ngrmax),inmin(maxmtx,ngrmax), &
  & lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ilower,imgmax,imgmin,ipos,j,ll,lunit2,lunit4,mg, &
  mg2hi,mg2lo,mpos,mpos2,mrec1,mrec2,n,nb,nrec1,nrec2,nw,ppos,ppos2,minlo, &
  maxhi,totzer,p,pend,ijk,ip,is,istop,lunit1,lunit3,m,mend,mglast,mtop,nlmax
  !
  nlmax = max0(nleg,npne1)-1
  totzer = 0
  nrec1 = 0
  nrec2 = 0
  mpos = 0
  ppos = 0
  mpos2 = 0
  ppos2 = 0
  lunit1 = lunit
  lunit2 = lunit1 + 5
  if(lunit<90)then
    lunit3 = lunit1 + 10
    lunit4 = lunit2 + 10
  else
    lunit3 = lunit1 - 10
    lunit4 = lunit2 - 10
  endif
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  mglast = ngi + 1
  !----
  !  read in first block
  !----
  mrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(mpos>0)then
    call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
    mrec2 = lnrec2(n,1,nsmxj) - 1
    ppos = lnrec2(n,2,nsmxj)
    if(ppos>0)then
      call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
      ppos = 0
    else
      porder(:psize)=0.0
    endif
  else
    matrix(:msize)=0.0
  endif
  mpos = 0
  10 continue
  ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 continue
  if(nb/=0) then
    ll = ll + nw
    if(ll>bufsiz) then
      call error('ecfad5',' ** error reading gendf, workspace exceeded',' ')
    else
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
  endif
  nw = n1h
  ig = n2h
  ig2lo = l2h
  ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
  ipos = nleg*nsig + 7
  mg = ngi - ig + 1
  mg2lo = ngi - ig2hi + 1
  mg2hi = ngi - ig2lo + 1
  if(inmin(nsmxj,mg)==0)then
    maxhi = 0
    minlo = ngi + 1
  else
    minlo = ngi - inmax(nsmxj,mg) + 1
    maxhi = ngi - inmin(nsmxj,mg) + 1
  endif
  30 continue
  if(mg<imgmin) then
    !----
    !  group is not in current group block
    !  check if old buffer has been copied to new buffer
    !----
    if(mpos<lnrec1(n,2,nsmxj))then
      !----
      !  copy rest of old buffer
      !----
      mend = lnrec1(n,2,nsmxj)
      do m = mpos+1,mend
        mpos2 = mpos2 + 1
        matrx2(mpos2) = matrix(m)
      enddo
      if(ppos<lnrec2(n,2,nsmxj))then
        pend = lnrec2(n,2,nsmxj)
        do p = ppos+1,pend
          ppos2 = ppos2 + 1
          pordr2(ppos2) = porder(p)
        enddo
      endif
    endif
    !----
    !  write new block to buffer
    !----
    if(mpos2>0)then
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos2
      maxm = max(mpos2,maxm)
      if(maxm>msize) then
        write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxm,msize
        call error('ecfad5',hsmg,' ')
      endif
      call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
      lnrec2(n,1,nsmxj) = nrec2 + 1
      lnrec2(n,2,nsmxj) = ppos2
      maxp = max(ppos2,maxp)
      if(maxp>psize) then
        write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxp,psize
        call error('ecfad4',hsmg,' ')
      endif
      call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
      mpos2 = 0
      ppos2 = 0
    endif
    !----
    !  read in next block
    !----
    n = n - 1
    imgmax = imgmin - 1
    imgmin = imgmax - ing(n) + 1
    mrec1 = lnrec1(n,1,nsmxj) - 1
    mpos = lnrec1(n,2,nsmxj)
    if(mpos>0) then
      call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
      mpos = 0
      mrec2 = lnrec2(n,1,nsmxj) - 1
      ppos = lnrec2(n,2,nsmxj)
      if(ppos>0)then
        call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
        ppos = 0
      else
        porder(:psize)=0.0
      endif
    else
      matrix(:msize)=0.0
    endif
    go to 30
  endif
  !----
  !  check for dummy last group
  !----
  if((ig/=ngi).or.(ig2lo/=ig2hi).or.(buff(ipos)/=0.0))then
    !----
    !  if threshold reaction and mg is not the first group in the group block
    !  or gaps in incident energy groups
    !  check whether the input buffer contains data for previous groups
    !----
    if(mg/=imgmax) then
      if(mg<mglast-1) then
        mtop = min(mglast-1,imgmax)
        do is = mg+1,mtop
          if(inmin(nsmxj,is)/=0) then
            do j = inmin(nsmxj,is),inmax(nsmxj,is)
              mpos = mpos + 1
              mpos2 = mpos2 + 1
              matrx2(mpos2) = matrix(mpos)
              do ip = 1,npne1 - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            enddo
          endif
        enddo
      endif
    endif
    if(inmin(nsmxj,mg)==0) then
      !----
      !  no data on buffer , add new data
      !----
      inmin(nsmxj,mg) = mg2lo
      inmax(nsmxj,mg) = mg2hi
      do j = mg2lo,mg2hi
        mpos2 = mpos2 + 1
        matrx2(mpos2) = buff(ipos)
        if(nleg>1) then
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
        endif
        if(npne1>nleg) then
          do ip = nleg,npne1 - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
            totzer = totzer + 1
          enddo
        endif
        ipos = ipos + nleg*nsig
      enddo
    else
      if(minlo<ig2lo) then
        !----
        !  more data on buffer so copy first
        !----
        ilower = minlo
        !----
        ! check for gaps
        !----
        if(ig2lo>maxhi)then
          !----
          !  copy data from buffer and and zeroise gap
          !----
          do j=ilower,maxhi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            if(npne1>1) then
              do ip = 1,npne1 - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            endif
          enddo
          do j=maxhi+1,ig2lo-1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = 0.0
            if(npne1>1) then
              do ip = 1,npne1 - 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = 0.0
                totzer = totzer + 1
              enddo
            endif
          enddo
          maxhi = ig2lo-1
        else
          do j = ilower,ig2lo - 1
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            if(npne1>1) then
              do ip = 1,npne1 - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            endif
          enddo
        endif
        ilower = ig2lo
      else if(minlo>ig2lo) then
        !----
        !  more data on gendf so copy first
        !----
        ilower = ig2lo
        istop = minlo - 1
        if( ig2hi<minlo) then
          istop = ig2hi
        endif
        do j = ilower,istop
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          if(nleg>1) then
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = buff(ipos+ip)
            enddo
          endif
          if(npne1>nleg) then
            do ip = 1,npne1 - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = 0.0
              totzer = totzer + 1
            enddo
          endif
          ipos = ipos + nleg*nsig
        enddo
        ilower = istop + 1
      else
        !----
        !  lower groups equal
        !----
        ilower = ig2lo
      endif
      !----
      !  check for gaps in data
      !----
      if(ilower<minlo) then
        do j = ilower,minlo-1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = 0.0
          do ip = 1,nlmax
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
            totzer = totzer + 1
          enddo
        enddo
        ilower = minlo
      endif
      if(maxhi<=ig2hi) then
        !----
        !  add data
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos) + buff(ipos)
          if(nleg>1) then
            do ip = 1,nleg - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos) + buff(ipos+ip)
            enddo
          endif
          if(npne1>nleg) then
            do ip = nleg,npne1 - 1
              ppos2 = ppos2 + 1
              ppos = ppos + 1
              pordr2(ppos2) = porder(ppos)
            enddo
          endif
          ipos = ipos + nleg*nsig
        enddo
        ilower = maxhi + 1
        if(ig2hi>maxhi) then
          do j = ilower,ig2hi
            mpos2 = mpos2 + 1
            matrx2(mpos2) = buff(ipos)
            if(nleg>1) then
              do ip = 1,nleg - 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = buff(ipos+ip)
              enddo
            endif
            if(npne1>nleg) then
              do ip = nleg,npne1 - 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = 0.0
                totzer = totzer + 1
              enddo
            endif
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
      else
        if(ig2hi>=ilower)then
          do j = ilower,ig2hi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos) + buff(ipos)
            if(nleg>1) then
              do ip = 1,nleg - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos) + buff(ipos+ip)
              enddo
            endif
            if(npne1>nleg) then
              do ip = nleg,npne1 - 1
                ppos2 = ppos2 + 1
                ppos = ppos + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            endif
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
          if(npne1>1) then
            do ip = 1,npne1 - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos)
            enddo
          endif
        enddo
        ilower = maxhi + 1
      endif
    endif
    inmin(nsmxj,mg) = min(inmin(nsmxj,mg),mg2lo)
    inmax(nsmxj,mg) = max(inmax(nsmxj,mg),mg2hi)
    mglast = mg
    if(ig<ngi) go to 10
  else
    !----
    !  dummy last group
    !----
  endif
  !----
  !  copy any remaining data
  !----
  if(mpos+1<lnrec1(n,2,nsmxj))then
    istop = lnrec1(n,2,nsmxj) -1
    do ijk=mpos,istop
      mpos=mpos+1
      mpos2=mpos2+1
      matrx2(mpos2) = matrix(mpos)
    enddo
  endif
  if(ppos+1<lnrec2(n,2,nsmxj))then
    istop = lnrec2(n,2,nsmxj) -1
    do ijk=ppos,istop
      ppos=ppos+1
      ppos2=ppos2+1
      pordr2(ppos2) = porder(ppos)
    enddo
  endif
  if(mpos2>0)then
    maxm = max(maxm,mpos2)
    if(maxm>msize) then
      write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxm,msize
      call error('ecfad5',hsmg,' ')
    endif
    lnrec1(n,1,nsmxj) = nrec1 + 1
    lnrec1(n,2,nsmxj) = mpos2
    call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
    maxp = max(maxp,ppos2)
    if(maxp>psize) then
      write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxp,psize
      call error('ecfad5',hsmg,' ')
    endif
    lnrec2(n,1,nsmxj) = nrec2 + 1
    lnrec2(n,2,nsmxj) = ppos2
    call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
  endif
  lunit = lunit3
  if(iprint>0) write(nsyso,"(' number of zero in inelastic matrix is ',i8)") &
  & totzer
  end subroutine ecfad5
  !
  subroutine ecfadth(nsmxj,ing,nreci,nsig,lenrec,matrx2,ngi,matrix,buff, &
  bufsiz,lunit,lnrec1,lnrec2,inmin,inmax,porder,pordr2,nleg,maxm,maxp,itempj, &
  thrmin)
  !
  ! add matrices - nleg = npne(nmsxj) > 1  (higher order matrices present)
  !
  ! .. scalar arguments ..
  integer bufsiz,lenrec,lunit,ngi,nleg,nreci,nsig,nsmxj,thrmin,maxm,maxp, &
  & itempj
  ! .. array arguments ..
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  integer ing(nreci),inmax(maxmtx,ngrmax),inmin(maxmtx,ngrmax), &
  lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ilower,imgmax,imgmin,ipos,j,ll,lunit2,lunit4,mg, &
  mg2hi,mg2lo,mpos,mpos2,mrec1,mrec2,n,nb,nrec1,nrec2,nw,ppos,ppos2,minlo, &
  maxhi,p,pend,ijk,ip,is,istop,lunit1,lunit3,m,mend,mtop,mglast
  !
  nrec1 = 0
  nrec2 = 0
  mpos = 0
  ppos = 0
  mpos2 = 0
  ppos2 = 0
  lunit1 = lunit
  lunit2 = lunit1 + 5
  if(lunit<90)then
    lunit3 = lunit1 + 10
    lunit4 = lunit2 + 10
  else
    lunit3 = lunit1 - 10
    lunit4 = lunit2 - 10
  endif
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  mglast = ngi + 1
  !----
  !  read in first block
  !----
  mrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(mpos>0)then
    call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
    mrec2 = lnrec2(n,1,nsmxj) - 1
    ppos = lnrec2(n,2,nsmxj)
    if(ppos>0)then
      call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
      ppos = 0
    else
      porder(:psize)=0.0
    endif
  else
    matrix(:msize)=0.0
  endif
  mpos = 0
  10 continue
  ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 continue
  if(nb/=0) then
    ll = ll + nw
    if(ll>bufsiz) then
      call error('ecfadth',' ** error reading gendf, workspace exceeded',' ')
    else
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
  endif
  nw = n1h
  ig = n2h
  ig2lo = l2h
  ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
  mg = ngi - ig + 1
  mg2lo = ngi - ig2hi + 1
  mg2hi = ngi - ig2lo + 1
  ipos = nleg*nsig + 7
  if(inmin(nsmxj,mg)==0)then
    maxhi = 0
    minlo = ngi + 1
  else
    minlo = ngi - inmax(nsmxj,mg) + 1
    maxhi = ngi - inmin(nsmxj,mg) + 1
  endif
  !----
  !  find correct group block for current group
  !----
  30 continue
  if(mg<imgmin) then
    !----
    !  group is not in current group block
    !  check if old buffer has been copied to new buffer
    !----
    if(mpos<lnrec1(n,2,nsmxj))then
      !----
      !  copy rest of old buffer
      !----
      mend = lnrec1(n,2,nsmxj)
      do m = mpos+1,mend
        mpos2 = mpos2 + 1
        matrx2(mpos2) = matrix(m)
      enddo
      if(ppos<lnrec2(n,2,nsmxj))then
        pend = lnrec2(n,2,nsmxj)
        do p = ppos+1,pend
          ppos2 = ppos2 + 1
          pordr2(ppos2) = porder(p)
        enddo
      endif
    endif
    !----
    !  write new block to buffer
    !----
    if(mpos2>0)then
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos2
      maxm = max(mpos2,maxm)
      if(maxm>msize) then
        write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxm,msize
        call error('ecfadth',hsmg,' ')
      endif
      call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
      lnrec2(n,1,nsmxj) = nrec2 + 1
      lnrec2(n,2,nsmxj) = ppos2
      maxp = max(ppos2,maxp)
      if(maxp>psize) then
        write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxp,psize
        call error('ecfadth',hsmg,' ')
      endif
      call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
      mpos2 = 0
      ppos2 = 0
    endif
    !----
    !  read in next block
    !----
    n = n - 1
    imgmax = imgmin - 1
    imgmin = imgmax - ing(n) + 1
    mrec1 = lnrec1(n,1,nsmxj) - 1
    mpos = lnrec1(n,2,nsmxj)
    if(mpos>0) then
      call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
      mpos = 0
      mrec2 = lnrec2(n,1,nsmxj) - 1
      ppos = lnrec2(n,2,nsmxj)
      if(ppos>0)then
        call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
        ppos = 0
      else
        porder(:psize)=0.0
      endif
    else
      matrix(:msize)=0.0
    endif
    go to 30
  endif
  !----
  !  check for dummy last group
  !----
  if((ig/=ngi).or.(ig2lo/=ig2hi).or.(buff(ipos)/=0.0))then
    !----
    !  if threshold reaction and mg is not the first group in the group block
    !  or gaps in incident energy groups
    !  check whether the input buffer contains data for previous groups
    !----
    if(mg/=imgmax) then
      if(mg<mglast-1) then
        mtop = min(mglast-1,imgmax)
        do is = mg+1,mtop
          if(inmin(nsmxj,is)/=0) then
            do j = inmin(nsmxj,is),inmax(nsmxj,is)
              mpos = mpos + 1
              mpos2 = mpos2 + 1
              matrx2(mpos2) = matrix(mpos)
              do ip = 1,nleg - 1
                ppos = ppos + 1
                ppos2 = ppos2 + 1
                pordr2(ppos2) = porder(ppos)
              enddo
            enddo
          endif
        enddo
      endif
    endif
    if(inmin(nsmxj,mg)==0) then
      !----
      !  no data on buffer , add new data
      !----
      inmin(nsmxj,mg) = mg2lo
      inmax(nsmxj,mg) = mg2hi
      thrmin = min(mg2lo,thrmin)
      do j = mg2lo,mg2hi
        mpos2 = mpos2 + 1
        matrx2(mpos2) = buff(ipos)
        do ip = 1,nleg - 1
          ppos2 = ppos2 + 1
          pordr2(ppos2) = buff(ipos+ip)
        enddo
        ipos = ipos + nleg*nsig
      enddo
    else
      if(minlo<ig2lo) then
        !----
        !  more data on buffer so copy first
        !----
        ilower = minlo
        !----
        !  check for gaps
        !----
        if(ig2lo>maxhi)then
          !----
          !  copy data from buffer and and zeroise gap
          !----
          do j=ilower,maxhi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            do ip = 1,nleg - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos)
            enddo
          enddo
          do j=maxhi+1,ig2lo-1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = 0.0
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = 0.0
            enddo
          enddo
          maxhi = ig2lo-1
        else
          do j = ilower,ig2lo - 1
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos)
            do ip = 1,nleg - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos)
            enddo
          enddo
        endif
        ilower = ig2lo
      else if(minlo>ig2lo) then
        !----
        !  more data on gendf so copy first
        !----
        ilower = ig2lo
        istop = minlo -1
        if(ig2hi<minlo)then
          istop = ig2hi
        endif
        do j = ilower,istop
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = istop + 1
      else
        !----
        !  lower groups equal
        !----
        ilower = ig2lo
      endif
      !----
      !  check for gaps in data
      !----
      if(ilower<minlo) then
        do j = ilower,minlo-1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = 0.0
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
          enddo
        enddo
        ilower = minlo
      endif
      if(maxhi<=ig2hi) then
        !----
        !  add data
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos) + buff(ipos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos) + buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        ilower = maxhi + 1
        if(ig2hi>maxhi) then
          do j = ilower,ig2hi
            mpos2 = mpos2 + 1
            matrx2(mpos2) = buff(ipos)
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = buff(ipos+ip)
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
      else
        if(ig2hi>=ilower)then
          do j = ilower,ig2hi
            mpos = mpos + 1
            mpos2 = mpos2 + 1
            matrx2(mpos2) = matrix(mpos) + buff(ipos)
            do ip = 1,nleg - 1
              ppos = ppos + 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = porder(ppos) + buff(ipos+ip)
            enddo
            ipos = ipos + nleg*nsig
          enddo
          ilower = ig2hi + 1
        endif
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos)
          enddo
        enddo
        ilower = maxhi + 1
      endif
    endif
    inmin(nsmxj,mg) = min(inmin(nsmxj,mg),mg2lo)
    inmax(nsmxj,mg) = max(inmax(nsmxj,mg),mg2hi)
    thrmin = min(inmin(nsmxj,mg),thrmin)
    mglast=mg
    if(ig<ngi) go to 10
  else
    !----
    !  dummy last group
    !----
  endif
  !----
  !  last group *
  !  copy any remaining data
  !----
  if(mpos+1<lnrec1(n,2,nsmxj))then
    istop = lnrec1(n,2,nsmxj) -1
    do ijk=mpos,istop
      mpos=mpos+1
      mpos2=mpos2+1
      matrx2(mpos2) = matrix(mpos)
    enddo
  endif
  if(ppos+1<lnrec2(n,2,nsmxj))then
    istop = lnrec2(n,2,nsmxj) -1
    do ijk=ppos,istop
      ppos=ppos+1
      ppos2=ppos2+1
      pordr2(ppos2) = porder(ppos)
    enddo
  endif
  if(mpos2>0)then
    maxm = max(mpos2,maxm)
    if(maxm>msize) then
      write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxm,msize
      call error('ecfadth',hsmg,' ')
    endif
    lnrec1(n,1,nsmxj) = nrec1 + 1
    lnrec1(n,2,nsmxj) = mpos2
    call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
    maxp = max(ppos2,maxp)
    if(maxp>psize) then
      write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxp,psize
      call error('ecfadth',hsmg,' ')
    endif
    lnrec2(n,1,nsmxj) = nrec2 + 1
    lnrec2(n,2,nsmxj) = ppos2
    call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
  endif
  lunit = lunit3
  end subroutine ecfadth
  !
  subroutine ecfout(nsmxj,ing,nreci,nleg,nsig,lenrec,reac,maxm,maxp,matrix, &
  porder,ngi,buff,bufsiz,inmin,inmax,lnrec1,lnrec2,npne,ismtx,listmx, &
  nsmtx,thrmin)
  ! .. scalar arguments ..
  integer nsmxj,nreci,lenrec,maxm,maxp,bufsiz,ngi,nsmtx,thrmin
  character*16 reac
  ! .. array arguments
  integer ing(nreci),ismtx(maxmtx),inmin(maxmtx,ngi),inmax(maxmtx,ngi), &
  lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx),npne(maxmtx)
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),porder(psize)
  character listmx(maxmtx)*16
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ll,nb,nleg,nsig,nw,mpos,ppos,i,j,imgmax,imgmin, &
  ip,ipos,ipos5,lunit1,lunit2,mg,mg2hi,mg2lo,n,nlegnsig,nrec1,nrec2
  real(kr) nu5sum
  logical first
  !
  first = .true.
  nsmxj = nsmxj + 1
  call search(reac,ismtx(nsmxj),listmx,nsmtx)
  npne(nsmxj)=nleg
  nrec1 = 0
  nrec2 = 0
  mpos = 0
  ppos = 0
  lunit1 = 80 + nsmxj
  lunit2 = 85 + nsmxj
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  20 ll = 1
  call listio(ning,0,0,buff,nb,nw)
  30 if(nb/=0) then
    ll = ll + nw
    call moreio(ning,0,0,buff(ll),nb,nw)
    go to 30
  endif
  nw = n1h
  ig = n2h
  ig2lo = l2h
  ig2hi = nw/(nleg*nsig)-1+ig2lo-1
  mg=ngi-ig+1
  mg2lo=ngi-ig2hi+1
  mg2hi=ngi-ig2lo+1
  ipos = nleg*nsig + 7
  !----
  !  check for dummy last group
  !----
  if((ig/=ngi).or.(ig2lo/=ig2hi).or.(buff(ipos)/=0.0))then
    !----
    !  threshold reaction ?
    !----
    if(first) then
      first = .false.
      25 if(mg<imgmin) then
        n = n - 1
        imgmax = imgmin -1
        imgmin = imgmax - ing(n) + 1
        go to 25
      endif
    else if(mg<imgmin) then
      !----
      !  write block away to buffer
      !----
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos
      maxm = max(mpos,maxm)
      if(maxm>msize) then
        write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxm,msize
        call error('ecfout',hsmg,' ')
      endif
      call ecwr4(lunit1,lenrec,matrix,mpos,nrec1,1)
      if( nleg>1) then
        lnrec2(n,1,nsmxj) = nrec2 + 1
        lnrec2(n,2,nsmxj) = ppos
        maxp = max(maxp,ppos)
        if(maxp>psize) then
          write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
          & i9,13h actual size=,i8)') mth,maxp,psize
          call error('ecfout',hsmg,' ')
        endif
        call ecwr4(lunit2,lenrec,porder,ppos,nrec2,1)
        ppos = 0
      endif
      mpos = 0
      !----
      !  find which group block current group is in
      !----
      35 n = n - 1
      imgmax = imgmin - 1
      imgmin = imgmax - ing(n) + 1
      if(mg<imgmin) go to 35
    endif
    !----
    !  position of group in current group block (reverse order)
    !----
    inmin(nsmxj,mg) = mg2lo
    inmax(nsmxj,mg) = mg2hi
    do i=1,mtherm
      if(mth==ntherm(i)) thrmin = min(mg2lo,thrmin)
    enddo
    !----
    !  store the number of neutrons per (n,anything) reaction
    !----
    if(mth==5) then
      nu5sum=0.0
      ipos5=ipos
      nlegnsig=nleg*nsig
      do j=mg2lo,mg2hi
        nu5sum=nu5sum+buff(ipos5)
        ipos5=ipos5+nlegnsig
      enddo
      nu5(mg)=nu5sum/rf(rnany,mg)
      ! maximum neutron production for reaction other than fission should be 4
      if((nu5(mg)>4.0).and.(iprint>0)) write(nsyso,1020) nu5(mg),mg
    endif
    nlegnsig=nleg*nsig
    do j=mg2lo,mg2hi
      mpos = mpos + 1
      matrix(mpos)=buff(ipos)
      if( nleg>1) then
        do ip=1,nleg-1
          ppos = ppos + 1
          porder(ppos)=buff(ipos+ip)
        enddo
      endif
      ipos=ipos+nlegnsig
    enddo
    if(ig<ngi) go to 20
  endif
  !----
  !  write final block to buffer
  !----
  if(mpos>0) then
    lnrec1(n,1,nsmxj) = nrec1 + 1
    lnrec1(n,2,nsmxj) = mpos
    maxm = max(mpos,maxm)
    if(maxm>msize) then
      write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=, &
      & i9,13h actual size=,i8)') mth,maxm,msize
      call error('ecfout',hsmg,' ')
    endif
    call ecwr4(lunit1,lenrec,matrix,mpos,nrec1,1)
    if( nleg>1) then
      lnrec2(n,1,nsmxj) = nrec2 + 1
      lnrec2(n,2,nsmxj) = ppos
      maxp = max(ppos,maxp)
      if(maxp>psize) then
        write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=, &
        & i9,13h actual size=,i8)') mth,maxp,psize
        call error('ecfout',hsmg,' ')
      endif
      call ecwr4(lunit2,lenrec,porder,ppos,nrec2,1)
    endif
  endif
  return
  !
  1020  format(' **warning*** number of neutrons form (n,anything)>4',/ &
  & ' there are', 1pe12.5, 'neutrons in group',i6)
  end subroutine ecfout
  !
  subroutine buffin(nsmxj,n,lenrec,mpos,ppos,matrix,maxm,porder,maxp,lnrec1, &
  lnrec2,npne,lunit)
  ! .. scalar arguments ..
  integer lenrec,nsmxj,n,mpos,ppos,maxm,maxp,npne,lunit
  ! .. array arguments ..
  integer lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  real(kr) matrix(maxm),porder(maxp)
  ! .. local scalars ..
  integer lunit1,lunit2,nrec1,nrec2
  !
  lunit1 = lunit
  lunit2 = lunit1 + 5
  nrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(nrec1>=0) call ecrr4(lunit1,lenrec,matrix,mpos,nrec1,1)
  if(npne>1) then
    nrec2 = lnrec2(n,1,nsmxj) - 1
    ppos = lnrec2(n,2,nsmxj)
    if(nrec2>=0) call ecrr4(lunit2,lenrec,porder,ppos,nrec2,1)
  endif
  end subroutine buffin
  !
  subroutine ecrfis(nleg,nsig,flux,totf,profis,fispec,xdelay)
  ! .. scalar arguments ..
  integer nleg,nsig,totf
  ! .. array arguments ..
  real(kr) flux(ngrmax),fispec(ngrmax),profis(ngrmax),xdelay(ngrmax)
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ll,nb,nw,ipos,j,mg,mg2hi,mg2lo
  real(kr) nudig,fiss,fp,fp1,fp2,fprod
  ! .. local arrays ..
  integer, parameter :: nbuf=10000
  real(kr) buff(nbuf)
  !----
  !  partial fission
  !  modified to read compressed fission matrix - see njoy-groupr-gendf
  !----
  !  zero array spec on first entry only
  !----
  if(totf==0) spec(:ngrmax)=0.0d0
  if(totf<=0) then
    if (iprint>0) write(nsyso,1000) mth
    !----
    !  read a list record and consider its type
    !
    !  a record with the group index, ig=0, contains a spectrum.
    !
    !  a record with the index to the lowest nonzero group, ig2lo=0
    !  contains flux,production cross section for group ig.
    !
    !  a record with ig non zero and iglo non zero is a matrix record.
    !----
    !  set default values
    !----
    fprod=0
    mg2lo = 1
    mg2hi = 1
    !----
    !  read record
    !----
    10 ll = 1
    call listio(ning,0,0,buff,nb,nw)
    20 if(nb/=0) then
      ll = ll + nw
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
    ig2lo = l2h
    nw = n1h
    ig = n2h
    ig2hi = nw/ (nleg*nsig)  + ig2lo - 1
    ipos = nleg*nsig + 6
    mg = ngi - ig + 1
    !----
    !  spectrum record
    !----
    if(ig==0) then
      mg2lo = ngi - ig2hi + 1
      mg2hi = ngi - ig2lo + 1
      do j = mg2hi,mg2lo,-1
        spec(j) = spec(j) + buff(ipos)
        ipos = ipos + nleg*nsig
      enddo
      go to 10
    endif
    if(ig2lo==0) then
      !----
      !  flux,production cross section record
      !----
      fp = buff(ipos+1)
      profis(mg) = profis(mg) + fp
      fprod = fprod +fp* flux(mg)
      if(n2h<ngi) then
        go to 10
      else
        do j = mg2hi,mg2lo,-1
          spec(j) = spec(j) * fprod
        enddo
        go to 90
      endif
    endif
    !----
    !  matrix record
    !----
    do j = mg2hi,mg2lo,-1
      spec(j) = spec(j) * fprod
    enddo
    mg2lo = ngi - ig2hi + 2
    mg2hi = ngi - ig2lo + 1
    ipos=ipos+1
    fp1=flux(mg)
    fp2=profis(mg)
    do j = mg2hi,mg2lo,-1
      fp=buff(ipos)
      spec(j) = spec(j) + fp * fp1
      fp2 = fp2 + fp
      ipos = ipos + nleg*nsig
    enddo
    profis(mg)=fp2
    !----
    !  2nd or later matrix record
    !  once the first matrix record processed all remaining records
    !  must be matrix - coded to avoid tests on record type.
    !----
    60 if(n2h<ngi) then
      ll = 1
      call listio(ning,0,0,buff,nb,nw)
      70 if(nb/=0) then
        ll = ll + nw
        call moreio(ning,0,0,buff(ll),nb,nw)
        go to 70
      endif
      ig2lo = l2h
      nw = n1h
      ig = n2h
      ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
      mg = ngi - ig + 1
      mg2lo = ngi - ig2hi + 1
      mg2hi = ngi - ig2lo + 1
      ipos = nleg*nsig + 7
      fp1=flux(mg)
      fp2=profis(mg)
      do j = mg2hi,mg2lo,-1
        fp=buff(ipos)
        spec(j) = spec(j) + fp * fp1
        fp2 = fp2 + fp
        ipos = ipos + nleg*nsig
      enddo
      profis(mg)=fp2
      go to 60
    endif
  endif
  !----
  ! add prompt and delayed (if present) fission matrices
  !
  !  nudig : nu bar delayed
  !  fiss  :  fission cross-section
  !  buff(ipos) : prompt fission matrix
  !  chid : fission spectrum for delayed neutrons (summed over time group i)
  !----
  90 if(totf==0)then
    if(delay)then
      do mg = 1,ngi
        nudig = nud(mg)
        fiss = xsec(nfissn,1,mg)
        do j = 1,ngi
          xdelay(mg) = xdelay(mg) + nudig*fiss*chid(j)
          spec(j) = spec(j) + nudig*fiss*chid(j)*flux(mg)
        enddo
      enddo
    endif
  endif
  totf = -1
  return
  1000 format (/' reaction ',i3,' has been summed to form the fission spectrum')
  end subroutine ecrfis
  !
  subroutine ecfthr(nsmxj,ing,nreci,nsig,lenrec,matrx2,ngi,matrix,buff, &
  bufsiz,lunit,lnrec1,lnrec2,inmin,inmax,porder,pordr2,nleg,maxm,maxp,itempj, &
  thrmin)
  !
  ! add matrices - nleg = npne(nmsxj) > 1  (higher order matrices present)
  !
  ! .. scalar arguments ..
  integer bufsiz,lenrec,lunit,ngi,nleg,nreci,nsig,nsmxj,maxm,maxp,itempj, &
  thrmin
  ! .. array arguments
  real(kr) buff(bufsiz)
  real(kr) matrix(msize),matrx2(msize),porder(psize),pordr2(psize)
  !
  integer ing(nreci),inmax(maxmtx,ngrmax),inmin(maxmtx,ngrmax), &
  lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  ! .. local scalars ..
  integer ig,ig2hi,ig2lo,ilower,imgmax,imgmin,ipos,j,ll,lunit2,lunit4,mg, &
  mg2hi,mg2lo,mpos,mpos2,mrec1,mrec2,n,nb,nrec1,nrec2,nw,ppos,ppos2,minlo, &
  maxhi,npne,i,ip,lunit1,lunit3,m,nn
  !
  nrec1 = 0
  nrec2 = 0
  mpos = 0
  ppos = 0
  mpos2 = 0
  ppos2 = 0
  lunit1 = lunit
  lunit2 = lunit1 + 5
  if(lunit<90)then
    lunit3 = lunit1 + 10
    lunit4 = lunit2 + 10
  else
    lunit3 = lunit1 - 10
    lunit4 = lunit2 - 10
  endif
  n = nreci
  imgmax = ngi
  imgmin = imgmax - ing(n) + 1
  !----
  !  read in buffer for first group in this reaction
  !----
  mrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(mpos>0)then
    call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
    mpos = 0
    mrec2 = lnrec2(n,1,nsmxj) - 1
    ppos = lnrec2(n,2,nsmxj)
    if(ppos>0)then
      call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
      ppos = 0
    else
      porder(:psize)=0.0
    endif
  else
    matrix(:msize)=0.0
  endif
  10 continue
  ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 continue
  if(nb/=0) then
    ll = ll + nw
    if(ll<=bufsiz) then
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 20
    endif
  else
    go to 30
  endif
  write(nsyso,9000)
  call error('ecfad2',' ** error reading gendf, workspace exceeded',' ')
  30 continue
  nw = n1h
  ig = n2h
  if(ig<=itempj)then
    ig2lo = l2h
    ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
    mg = ngi - ig + 1
    mg2lo = ngi - ig2hi + 1
    mg2hi = ngi - ig2lo + 1
    if(inmin(nsmxj,mg)==0)then
      maxhi = 0
      minlo = ngi + 1
    else
      minlo = ngi - inmax(nsmxj,mg) + 1
      maxhi = ngi - inmin(nsmxj,mg) + 1
    endif
    if(mg<imgmin) then
      !----
      !  write block away to buffer
      !----
      lnrec1(n,1,nsmxj) = nrec1 + 1
      lnrec1(n,2,nsmxj) = mpos2
      maxm = max(mpos2,maxm)
      call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
      lnrec2(n,1,nsmxj) = nrec2 + 1
      lnrec2(n,2,nsmxj) = ppos2
      maxp = max(ppos2,maxp)
      call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
      mpos2 = 0
      ppos2 = 0
      !----
      !  read in next block
      !----
      n = n - 1
      imgmax = imgmin - 1
      imgmin = imgmax - ing(n) + 1
      mrec1 = lnrec1(n,1,nsmxj) - 1
      mpos = lnrec1(n,2,nsmxj)
      if(mpos>0)then
        call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
        mpos = 0
        mrec2 = lnrec2(n,1,nsmxj) - 1
        ppos = lnrec2(n,2,nsmxj)
        if(ppos>0)then
          call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
          ppos = 0
        else
          porder(:psize)=0.0
        endif
      else
        matrix(:msize)=0.0
      endif
    endif
    ipos = nleg*nsig + 7
    if(inmin(nsmxj,mg)==0) then
      !----
      !  no data on buffer , add new data
      !----
      inmin(nsmxj,mg) = mg2lo
      inmax(nsmxj,mg) = mg2hi
      thrmin = min(mg2lo,thrmin)
      do j = mg2lo,mg2hi
        mpos2 = mpos2 + 1
        matrx2(mpos2) = buff(ipos)
        do ip = 1,nleg - 1
          ppos2 = ppos2 + 1
          pordr2(ppos2) = buff(ipos+ip)
        enddo
        ipos = ipos + nleg*nsig
      enddo
    else
      if(minlo<ig2lo) then
        !----
        !  more data on buffer so copy first
        !----
        ilower = ig2lo
        do j = minlo,ig2lo - 1
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos)
          enddo
        enddo
      else if(minlo>ig2lo) then
        !----
        !  more data on gendf so copy first
        !----
        ilower = minlo
        do j = ig2lo,minlo - 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
      else
        !----
        !  lower groups equal
        !----
        ilower = ig2lo
      endif
      !----
      !  check for gaps in data
      !----
      if(ilower<minlo) then
        do j = ilower,minlo-1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = 0.0
          do ip = 1,nleg - 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = 0.0
          enddo
        enddo
        ilower = minlo
      endif
      if(maxhi<=ig2hi) then
        !----
        !  add data
        !----
        do j = ilower,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        if(ig2hi>maxhi) then
          do j = maxhi + 1,ig2hi
            mpos2 = mpos2 + 1
            matrx2(mpos2) = buff(ipos)
            do ip = 1,nleg - 1
              ppos2 = ppos2 + 1
              pordr2(ppos2) = buff(ipos+ip)
            enddo
            ipos = ipos + nleg*nsig
          enddo
        endif
      else
        do j = ilower,ig2hi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = buff(ipos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = buff(ipos+ip)
          enddo
          ipos = ipos + nleg*nsig
        enddo
        do j = ig2hi + 1,maxhi
          mpos = mpos + 1
          mpos2 = mpos2 + 1
          matrx2(mpos2) = matrix(mpos)
          do ip = 1,nleg - 1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos)
          enddo
        enddo
      endif
    endif
    inmin(nsmxj,mg) = min(inmin(nsmxj,mg),mg2lo)
    inmax(nsmxj,mg) = max(inmax(nsmxj,mg),mg2hi)
    thrmin = min(inmin(nsmxj,mg),thrmin)
    if(ig<ngi) go to 10
  endif
  !----
  !  last group
  !----
  if( itempj/=imgmin) then
    !----
    !  remainder of group block to be copied
    !----
    do i=imgmin,itempj-1
      mg2lo = inmin(nsmxj,i)
      mg2hi = inmax(nsmxj,i)
      do m = mg2lo,mg2hi
        mpos = mpos + 1
        mpos2 = mpos2 + 1
        matrx2(mpos2) = matrix(mpos)
        if(nleg>1)then
          do ip = 1,nleg-1
            ppos = ppos + 1
            ppos2 = ppos2 + 1
            pordr2(ppos2) = porder(ppos)
          enddo
        endif
      enddo
    enddo
  endif
  !----
  !  write block away to buffer
  !----
  lnrec1(n,1,nsmxj) = nrec1 + 1
  lnrec1(n,2,nsmxj) = mpos2
  maxm = max(mpos2,maxm)
  if(maxm>msize) then
    write(hsmg,'(34harray matrix exceeded for reaction,i5,9h: needed=,i9, &
    13h actual size=,i8)') mth,maxm,msize
    call error('ecfthr',hsmg,' ')
  endif
  call ecwr4(lunit3,lenrec,matrx2,mpos2,nrec1,1)
  lnrec2(n,1,nsmxj) = nrec2 + 1
  lnrec2(n,2,nsmxj) = ppos2
  maxp = max(ppos2,maxp)
  if(maxp>psize) then
    write(hsmg,'(34harray porder exceeded for reaction,i5,9h: needed=,i9, &
    13h actual size=,i8)') mth,maxp,psize
    call error('ecfthr',hsmg,' ')
  endif
  call ecwr4(lunit4,lenrec,pordr2,ppos2,nrec2,1)
  mpos2 = 0
  ppos2 = 0
  n = n - 1
  !----
  !  copy from input buffer to output buffer
  !----
  if(n/=0) then
    do nn = n,1,-1
      mrec1 = lnrec1(nn,1,nsmxj) - 1
      mpos = lnrec1(nn,2,nsmxj)
      if(mpos>0)then
        call ecrr4(lunit1,lenrec,matrix,mpos,mrec1,1)
        lnrec1(nn,1,nsmxj) = nrec1 + 1
        call ecwr4(lunit3,lenrec,matrix,mpos,nrec1,1)
        mrec2 = lnrec2(nn,1,nsmxj) - 1
        ppos = lnrec2(nn,2,nsmxj)
        if(ppos>0) then
          call ecrr4(lunit2,lenrec,porder,ppos,mrec2,1)
          lnrec2(nn,1,nsmxj) = nrec2 + 1
          call ecwr4(lunit4,lenrec,porder,ppos,nrec2,1)
        endif
      endif
    enddo
  endif
  !
  lunit = lunit3
  return
  !
  9000 format (/,' ****** error reading gendf, workspace exceeded')
  end subroutine ecfthr
  !
  subroutine ecfwrt(nsmxj,n,lenrec,mpos,ppos,matrix,maxm,porder,maxp,lnrec1, &
  lnrec2,npne,lunit)
  ! .. scalar arguments ..
  integer lenrec,nsmxj,n,mpos,ppos,maxm,maxp,npne,lunit
  ! .. array arguments ..
  integer lnrec1(ngblk*5,2,maxmtx),lnrec2(ngblk*5,2,maxmtx)
  real(kr) matrix(maxm),porder(maxp)
  ! .. local scalars ..
  integer lunit1,lunit2,nrec1,nrec2
  !
  lunit1 = lunit
  lunit2 = lunit1 + 5
  nrec1 = lnrec1(n,1,nsmxj) - 1
  mpos = lnrec1(n,2,nsmxj)
  if(nrec1>=0) call ecwr4(lunit1,lenrec,matrix,mpos,nrec1,1)
  if(npne>1) then
    nrec2 = lnrec2(n,1,nsmxj) - 1
    ppos = lnrec2(n,2,nsmxj)
    if(nrec2>=0) call ecwr4(lunit2,lenrec,porder,ppos,nrec2,1)
  endif
  end subroutine ecfwrt
  !
  subroutine refin(nrect,lenrec)
  !  code to read the reference data package and add data for a new substance
  ! .. scalar arguments ..
  integer nrect
  ! .. array arguments ..
  integer lenrec(3)
  ! .. local scalars ..
  integer, parameter :: nbuf=10000
  real(kr) avog,er,fcwe
  integer i,j,la,lbuff,ldis,le,lecg,lecng,lefg,lefng,liel,llfr,llstfr,llstgp, &
  & llstmi,llstmx,lsybel,lsycor,mrect,nbl,nfgg,nfrecs,nmpn,ns,nsgp
  ! .. local arrays ..
  real(kr) buff(nbuf),awt(maxel),dis(maxel),ecg(maxel),ecng(maxel),efg(maxel), &
  & efng(maxel),eg(ngrmax)
  integer nbuff(nbuf),lfr(maxnfr),lrec(3,4)
  character listgp(5)*16,sybel(maxel)*16
  character(len=16) cbuff(maxcbf)
  !
  mrect = 1
  !----
  !  read identifier record for reference package
  !----
  call ecrint(ninr,lenrec(1),nbuff,lenrec(1),nrect,1)
  ! number of fortran records in this package
  nfrecs = nbuff(6)
  ! print out data block structure
  if (iprint>0) write(nsyso,9130)
  if (iprint>0) write(nsyso,9140)
  lbuff = 9
  do j = 1,nfrecs
    lrec(1,j) = nbuff(lbuff+1)
    lrec(2,j) = nbuff(lbuff+2)
    lrec(3,j) = nbuff(lbuff+3)
    if (iprint>0) write(nsyso,9150) j,lrec(1,j),lrec(2,j),lrec(3,j)
    lbuff = lbuff + 3
  enddo
  if (iprint>0) write(nsyso,9170) (nbuff(i),i=1,12)
  lrec(1,1) = 1
  lrec(2,1) = mrect + 1
  call ecrint(ninr,lenrec(1),nbuff,lrec(3,2),nrect,1)
  ! nel is the number of reference elements
  nel = nbuff(1)
  ! set up extra addresses for block 2
  ! address of iel
  liel = 2
  do i = 1,nel
    iel(i) = nbuff(liel)
    liel = liel + 1
  enddo
  ! number of x-sections
  nsmic = nbuff(liel)
  ! number of matrices
  nsmtx = nbuff(liel+1)
  ! number of response functions
  nfr = nbuff(liel+2)
  ! address of lfr
  llfr = liel + 3
  do i = 1,nfr
    lfr(i) = nbuff(llfr)
    llfr = llfr + 1
  enddo
  ! nsgp   - number of gamma production x-section names
  nsgp = nbuff(llfr)
  ! nmpn - maximum number of pn
  nmpn = nbuff(llfr+1)
  ! number of fine groups
  nfg = nbuff(llfr+2)
  ! number of gamma groups
  nfgg = nbuff(llfr+3)
  ! add to reference file
  nel = nel + 1
  !  iel   - isotope identifier
  !  awt   - atomic weight
  ! sycor  - isotope name
  ! sybel  - element name
  read(nsysi,*) iel(nel),awt(nel),sycor(nel),sybel(nel)
  read(nsysi,*) ecng(nel),ecg(nel),efng(nel),efg(nel),dis(nel)
  ! record block 2: integer data
  nbuff(1) = nel
  lbuff = 1
  nbuff(lbuff+1:lbuff+nel) = iel(:nel)
  lbuff = lbuff + nel
  nbuff(lbuff+1) = nsmic
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nsmtx
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nfr
  lbuff = lbuff + 1
  nbuff(lbuff+1:lbuff+nfr) = lfr(:nfr)
  lbuff = lbuff + nfr
  nbuff(lbuff+1) = nsgp
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nmpn
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nfg
  lbuff = lbuff + 1
  nbuff(lbuff+1) = nfgg
  lbuff = lbuff + 1
  ! block 2: start record
  lrec(2,2) = mrect + 1
  ! block 2: number of words
  lrec(3,2) = lbuff
  call ecwint(nour,lenrec(1),nbuff,lbuff,mrect,1)
  call ecrc16(ninr,lenrec(2),cbuff,lrec(3,3),nrect,1)
  ! set up extra addresses for block 3
  ! start address of sycor array
  lsycor = 1
  ! start address of sybel array
  lsybel = lsycor + nel - 1
  do i = 1,nel - 1
    sycor(i) = cbuff(lsycor)
    lsycor = lsycor + 1
    sybel(i) = cbuff(lsybel)
    lsybel = lsybel + 1
  enddo
  ! start address of listmi array
  llstmi = lsybel
  do i = 1,nsmic
    listmi(i) = cbuff(llstmi)
    llstmi = llstmi + 1
  enddo
  ! start address of listmx array
  llstmx = llstmi
  do i = 1,nsmtx
    listmx(i) = cbuff(llstmx)
    llstmx = llstmx + 1
  enddo
  ! start address of listfr array
  llstfr = llstmx
  do i = 1,nfr
    listfr(i) = cbuff(llstfr)
    llstfr = llstfr + 1
  enddo
  ! start address of listgp array
  llstgp = llstfr
  do i = 1,nsgp
    listgp(i) = cbuff(llstgp)
    llstgp = llstgp + 1
  enddo
  ! record block 3: character data
  do i = 1,nel
    cbuff(i) = sycor(i)
  enddo
  lbuff = nel
  do i = 1,nel
    cbuff(lbuff+i) = sybel(i)
  enddo
  lbuff = lbuff + nel
  do j = 1,nsmic
    cbuff(lbuff+j) = listmi(j)
  enddo
  lbuff = lbuff + nsmic
  do j = 1,nsmtx
    cbuff(lbuff+j) = listmx(j)
  enddo
  lbuff = lbuff + nsmtx
  do j = 1,nfr
    cbuff(lbuff+j) = listfr(j)
  enddo
  lbuff = lbuff + nfr
  if(nfgg>0) then
    do j = 1,nfgg
      cbuff(lbuff+j) = listgp(j)
    enddo
    lbuff = lbuff + nfgg
  endif
  ! block 3: start record
  lrec(2,3) = mrect + 1
  ! block 3: number of data words
  lrec(3,3) = lbuff
  call ecwc16(nour,lenrec(2),cbuff,lrec(3,3),mrect,1)
  ! read block 4
  call ecrr4(ninr,lenrec(3),buff,lrec(3,4),nrect,1)
  ! set up extra addresses for block 4
  ! start address of a array
  la = 1
  ! start address of ecng array
  lecng = la + nel - 1
  ! start address of ecg array
  lecg = lecng + nel - 1
  ! start address of efng array
  lefng = lecg + nel - 1
  ! start address of efg array
  lefg = lefng + nel - 1
  ! start address of dis array
  ldis = lefg + nel - 1
  ! start address of neutron energy fine group array
  le = ldis + nel - 1 + 3
  do i = 1,nel - 1
    !   a    - atomic weight
    awt(i) = buff(la)
    la = la + 1
    ! ecng - non-gamma capture energy (mev)
    ecng(i) = buff(lecng)
    lecng = lecng + 1
    ! ecg  - gamma capture energy (mev)
    ecg(i) = buff(lecg)
    lecg = lecg + 1
    ! efng - non-gamma fission energy (mev)
    efng(i) = buff(lefng)
    lefng = lefng + 1
    ! efg  - gamma fission energy (mev)
    efg(i) = buff(lefg)
    lefg = lefg + 1
    ! dis  - disintegration constant (s-1)
    dis(i) = buff(ldis)
    ldis = ldis + 1
  enddo
  lbuff = ldis
  avog = buff(lbuff)
  er = buff(lbuff+1)
  fcwe = buff(lbuff+2)
  do i = 1,nfg + 1
    egrid(i) = buff(le)
    le = le + 1
  enddo
  lbuff = le
  if(nfgg>0) then
    do j = 1,nfgg + 1
      eg(j) = buff(lbuff)
      lbuff = lbuff + 1
    enddo
  endif
  ! write block 4
  ! record block 4: real data
  do i = 1,nel
    buff(i) = awt(i)
  enddo
  lbuff = nel
  do i = 1,nel
    buff(lbuff+i) = ecng(i)
  enddo
  lbuff = lbuff + nel
  do i = 1,nel
    buff(lbuff+i) = ecg(i)
  enddo
  lbuff = lbuff + nel
  do i = 1,nel
    buff(lbuff+i) = efng(i)
  enddo
  lbuff = lbuff + nel
  do i = 1,nel
    buff(lbuff+i) = efg(i)
  enddo
  lbuff = lbuff + nel
  do i = 1,nel
    buff(lbuff+i) = dis(i)
  enddo
  lbuff = lbuff + nel
  buff(lbuff+1) = avog
  buff(lbuff+2) = er
  buff(lbuff+3) = fcwe
  lbuff = lbuff + 3
  do j = 1,nfg + 1
    buff(lbuff+j) = egrid(j)
  enddo
  lbuff = lbuff + nfg + 1
  if(nfgg>0) then
    do j = 1,nfgg + 1
      buff(lbuff+j) = eg(j)
    enddo
    lbuff = lbuff + nfgg + 1
  endif
  lrec(2,4) = mrect + 1
  lrec(3,4) = lbuff
  call ecwr4(nour,lenrec(3),buff,lrec(3,4),mrect,1)
  if (iprint>0) write(nsyso,9000) nel
  !  iel   - isotope identifier
  !   a    - atomic weight
  ! sycor  - isotope name
  ! sybel  - element name
  do i = 1,nel
    write(nsyso,9010) iel(i),awt(i),sycor(i),sybel(i)
  enddo
  ! ecng - non-gamma capture energy (mev)
  ! ecg  - gamma capture energy (mev)
  ! efng - non-gamma fission energy (mev)
  ! efg  - gamma fission energy (mev)
  ! dis  - disintegration constant (s-1)
  if (iprint>0) then
    write(nsyso,9020)
    do i = 1,nel
      write(nsyso,9030) ecng(i),ecg(i),efng(i),efg(i),dis(i)
    enddo
    ! nsmic  - number of x-section type names
    write(nsyso,9040) nsmic
    ! listmi - x-section type name
    write(nsyso,9050) (ns,listmi(ns),ns=1,nsmic)
    ! nsmtx  - number of x-section matrix type names
    write(nsyso,9060) nsmtx
    ! listmx - x-section type matrix name
    write(nsyso,9050) (ns,listmx(ns),ns=1,nsmtx)
    ! nfr    - number of response function names
    write(nsyso,9070) nfr
    ! listfr - response function name
    ! lfr    - r.f. is a combination of the next lfr types
    if(nfr>0) write(nsyso,9080) (ns,listfr(ns),lfr(ns),ns=1,nfr)
    ! nsgp   - number of gamma production x-section names
    write(nsyso,9090) nsgp
    ! listgp - gamma production x-section name
    if(nsgp>0) write(nsyso,9050) (ns,listgp(ns),ns=1,nsgp)
    ! nmpn  - maximum order of pn
    write(nsyso,9100) nmpn
    ! fine group scheme:-
    ! egrid(nfg+1) is the lower bound of the group scheme
    write(nsyso,9110) nfg
    ! nfgg   - number of fine gamma groups
    write(nsyso,9120) nfgg
  endif
  ! finally, create identifier record block
  ! type number of this package
  nbuff(1) = 1
  ! name number of this package
  nbuff(2) = 1
  ! father number of this package
  nbuff(3) = 1
  ! structure number
  nbuff(4) = 1
  ! number of physical records in this package
  nbuff(5) = mrect
  ! number of fortran records in this package
  nfrecs = 4
  nbuff(6) = nfrecs
  ! library origin - 2 = jef
  nbuff(7) = 2
  !----
  !  zero unused variables
  !----
  nbuff(8) = 0
  nbuff(9) = 0
  lrec(3,1) = 9 + nfrecs*3
  lrec(2,1) = 1
  lbuff = 9
  do j = 1,nfrecs
    nbuff(lbuff+1) = lrec(1,j)
    nbuff(lbuff+2) = lrec(2,j)
    nbuff(lbuff+3) = lrec(3,j)
    lbuff = lbuff + 3
  enddo
  mrect = 0
  call ecwint(nour,lenrec(1),nbuff,lbuff,mrect,1)
  !----
  !  print out data block structure
  !----
  if (iprint>0) then
    write(nsyso,9180) mrect
    write(nsyso,9160)
    write(nsyso,9140)
    do nbl = 1,nfrecs
      write(nsyso,9150) nbl,lrec(1,nbl),lrec(2,nbl),lrec(3,nbl)
    enddo
  endif
  return
  !----
  !  formats
  !----
  9000 format (/,5x,'ecco reference file - number of elements:',i4,/,1x, &
  'element',7x,'atomic',9x,'isotope',9x,'element',/,1x,' number',7x,'weight', &
  9x,'  name ',9x,'  name ',/)
  9010 format (2x,i6,f12.3,12x,2a16)
  9020 format (/,3x,'non-gamma',5x,' gamma ',3x,'non-gamma',5x,' gamma ',5x, &
  'disint-',/,3x,' capture ',5x,'capture',3x,' fission ',5x,'fission',5x, &
  'egration',/,3x,'  energy ',5x,' energy',3x,'  energy ',5x,' energy',5x, &
  'constant',/,3x,'  (mev)  ',5x,' (mev) ',3x,'  (mev)  ',5x,' (mev) ',5x, &
  ' (s-1)  ',/)
  9030 format (1x,1p,5e12.5)
  9040 format (/,1x,'number of x-section type names:',i6,/,1x,' type ',10x, &
  ' name',/,1x,'number',/)
  9050 format (1x,i4,12x,a16)
  9060 format (/,1x,'number of x-section matrix type names:',i6,/,1x,' type ' &
  ,10x,' name',/,1x,'number',/)
  9070 format (/,1x,'number of response function names:',i6,/,1x,' type ',10x &
  ,' name',15x,'lfr code',/,1x,'number',/)
  9080 format (1x,i4,12x,a16,2x,i8)
  9090 format (/,1x,'number of gamma production x-section names:',i6,/,1x, &
  ' type ',10x,' name',/,1x,'number',/)
  9100 format (/,1x,'maximum order of pn:',i6)
  9110 format (/,1x,'number of fine groups:',i6)
  9120 format (/,1x,'number of fine gamma groups:',i6)
  9130 format (/,1x,'input reference data record block structure')
  9140 format (2x,'f rec',4x,'data',4x,'start ',4x,'  words ',/,2x,' no. ',4x &
  ,'type',3x,'record',6x,' stored ')
  9150 format (2x,i5,3i9)
  9160 format (/,1x,'output reference data record block structure')
  9170 format (1x,'identifier record:',12i6)
  9180 format (/,1x,'last record of reference data:',i6)
  end subroutine refin
  !
  subroutine ecfissi(nleg,nsig,flux,totf,profis,fispec,xdelay)
  ! .. scalar arguments ..
  integer totf
  ! .. array arguments
  real(kr) flux(ngrmax),fispec(ngrmax),profis(ngrmax),xdelay(ngrmax)
  ! .. local scalars ..
  integer, parameter :: nbuf=10000
  integer ig,ig2hi,ig2lo,ll,nb,nleg,nsig,nw,ipos,j,mg,mg2hi,mg2lo
  real(kr) buff(nbuf),nudig,fiss,fp,fp1,fp2,fprod
  !----
  !  fission
  !  modified to read compressed fission matrix - see njoy-groupr-gendf
  !----
  !  zero array spec
  !----
  spec(:ngrmax)=0.0d0
  totf = 1
  !----
  !  read a list record and consider its type
  !
  !  a record with the group index, ig=0, contains a spectrum.
  !
  !  a record with the index to the lowest nonzero group, ig2lo=0
  !  contains flux,production cross section for group ig.
  !
  !  a record with ig non zero and iglo non zero is a matrix record.
  !
  !  set default values
  !----
  fprod=0
  mg2lo = 1
  mg2hi = 1
  ! read record
  10 ll = 1
  call listio(ning,0,0,buff,nb,nw)
  20 if(nb/=0) then
    ll = ll + nw
    call moreio(ning,0,0,buff(ll),nb,nw)
    go to 20
  endif
  ig2lo = l2h
  nw = n1h
  ig = n2h
  ig2hi = nw/ (nleg*nsig)  + ig2lo - 1
  ipos = nleg*nsig + 6
  mg = ngi - ig + 1
  !----
  ! spectrum record
  !----
  if(ig==0) then
    mg2lo = ngi - ig2hi + 1
    mg2hi = ngi - ig2lo + 1
    do j = mg2hi,mg2lo,-1
      spec(j) = spec(j) + buff(ipos)
      ipos = ipos + nleg*nsig
    enddo
    go to 10
  endif
  if(ig2lo==0) then
    !----
    !  flux,production cross section record
    !----
    fp = buff(ipos+1)
    profis(mg) = profis(mg) + fp
    fprod = fprod +fp* flux(mg)
    if(n2h<ngi) then
      go to 10
    else
      do j = mg2hi,mg2lo,-1
        spec(j) = spec(j) * fprod
      enddo
      go to 90
    endif
  endif
  !----
  !  matrix record
  !----
  do j = mg2hi,mg2lo,-1
    spec(j) = spec(j) * fprod
  enddo
  mg2lo = ngi - ig2hi + 2
  mg2hi = ngi - ig2lo + 1
  ipos=ipos+1
  fp1=flux(mg)
  fp2=profis(mg)
  do j = mg2hi,mg2lo,-1
    fp=buff(ipos)
    spec(j) = spec(j) + fp * fp1
    fp2 = fp2 + fp
    ipos = ipos + nleg*nsig
  enddo
  profis(mg)=fp2
  !----
  !  2nd or later matrix record
  !  once the first matrix record processed all remaining records
  !  must be matrix - coded to avoid tests on record type.
  !----
  60 if(n2h<ngi) then
    ll = 1
    call listio(ning,0,0,buff,nb,nw)
    70 if(nb/=0) then
      ll = ll + nw
      call moreio(ning,0,0,buff(ll),nb,nw)
      go to 70
    endif
    ig2lo = l2h
    nw = n1h
    ig = n2h
    ig2hi = nw/ (nleg*nsig) - 1 + ig2lo - 1
    mg = ngi - ig + 1
    mg2lo = ngi - ig2hi + 1
    mg2hi = ngi - ig2lo + 1
    ipos = nleg*nsig + 7
    fp1=flux(mg)
    fp2=profis(mg)
    do j = mg2hi,mg2lo,-1
      fp=buff(ipos)
      spec(j) = spec(j) + fp * fp1
      fp2 = fp2 + fp
      ipos = ipos + nleg*nsig
    enddo
    profis(mg)=fp2
    go to 60
  endif
  !----
  ! add prompt and delayed (if present) fission matrices
  !
  !  nudig : nu bar delayed
  !  fiss  :  fission cross-section
  !  buff(ipos) : prompt fission matrix
  !  chid : fission spectrum for delayed neutrons (summed over time group i)
  !----
  90 if(delay)then
    do mg = 1,ngi
      nudig = nud(mg)
      fiss = xsec(nfissn,1,mg)
      do j = 1,ngi
        xdelay(mg) = xdelay(mg) + nudig*fiss*chid(j)
        spec(j) = spec(j) + nudig*fiss*chid(j)*flux(mg)
      enddo
    enddo
  endif
  end subroutine ecfissi
  !
  subroutine ecrc16(lunit,lenrec,ch,nwords,irec,nextc)
  ! code to read a direct access dataset in character*16 form
  ! .. scalar arguments ..
  integer irec,lenrec,lunit,nwords,nextc,maxcbf,dimabs
  ! .. array arguments ..
  character(len=16) ch(nwords)
  ! .. local scalars ..
  integer j,jl,ju,nr,nrecs
  !
  ju = 0
  ! number of records to be read
  nrecs = (nwords+lenrec-1)/lenrec
  do nr = 1,nrecs
    jl = ju + 1
    ju = ju + lenrec
    if(ju>nwords) ju = nwords
    irec = irec + 1
    read(lunit,rec=irec) (ch(j),j=jl,ju)
  enddo
  return
  end subroutine ecrc16
  !
  subroutine ecrint(lunit,lenrec,nc,nwords,irec,nextn)
  ! code to read a direct access dataset in integer form
  ! .. scalar arguments ..
  integer irec,lenrec,lunit,nwords,nextn,dimabs
  ! ..
  ! .. array arguments ..
  integer nc(nwords)
  ! ..
  ! .. local scalars ..
  integer j,jl,ju,nr,nrecs
  character*4096 ctemp
  !
  ju = 0
  ! number of records to be read
  nrecs = (nwords+lenrec-1)/lenrec
  do nr = 1,nrecs
    jl = ju + 1
    ju = ju + lenrec
    if(ju>nwords) ju = nwords
    irec = irec + 1
    read(lunit,rec=irec) (nc(j),j=jl,ju)
  enddo
  return
  end subroutine ecrint
  !
  subroutine ecrr4(lunit,lenrec,bc,nwords,irec,nextr)
  ! code to read a direct access block from unit 'lunit'
  ! with a record length of lenrec real words
  ! .. scalar arguments ..
  integer irec,lenrec,lunit,nwords,nextr,dimabs
  ! .. array arguments ..
  real(kr) bc(nwords)
  ! .. local scalars ..
  integer j,jl,ju,nr,nrecs
  real, dimension(:), allocatable :: dummy
  !
  nrecs = (nwords+lenrec-1)/lenrec
  ju = 0
  allocate(dummy(nwords))
  do nr = 1,nrecs
    jl = ju + 1
    ju = ju + lenrec
    if(ju>nwords) ju = nwords
    irec = irec + 1
    read(lunit,rec=irec) (dummy(j),j=jl,ju)
  enddo
  bc(:nwords)=dummy(:nwords) ! single to double precision
  deallocate(dummy)
  return
  end subroutine ecrr4
  !
  subroutine ecwc16(lunit,lenrec,nc,nwords,irec,nextc)
  ! code to write a direct access block to unit 'lunit'
  ! with record length of lenrec 16 byte characters
  ! total number of records to be written in this block
  ! .. scalar arguments ..
  integer irec,lenrec,lunit,nwords,nextc,dimabs
  ! .. array arguments ..
  character(len=16) nc(nwords)
  ! .. local scalars ..
  integer j,jl,ju,nr,nrecs
  !
  nrecs = (nwords+lenrec-1)/lenrec
  ju = 0
  do nr = 1,nrecs
    jl = ju + 1
    ju = ju + lenrec
    if(ju>nwords) ju = nwords
    irec = irec + 1
    write(lunit,rec=irec) (nc(j),j=jl,ju)
  enddo
  return
  end subroutine ecwc16
  !
  subroutine ecwint(lunit,lenrec,nc,nwords,irec,nextn)
  ! code to write a direct access block to unit 'lunit'
  ! with a record length of lenrec integer words
  ! total number of records to be written in this block
  ! .. scalar arguments ..
  integer irec,lenrec,lunit,nwords,nextn,dimabs
  ! .. array arguments ..
  integer nc(nwords)
  ! .. local scalars ..
  integer j,jl,ju,nr,nrecs
  !
  nrecs = (nwords+lenrec-1)/lenrec
  ju = 0
  do nr = 1,nrecs
    jl = ju + 1
    ju = ju + lenrec
    if(ju>nwords) ju = nwords
    irec = irec + 1
    write(lunit,rec=irec) (nc(j),j=jl,ju)
  enddo
  return
  end subroutine ecwint
  !
  subroutine ecwr4(lunit,lenrec,bc,nwords,irec,nextr)
  ! code to write a direct access block to unit 'lunit'
  ! with a record length of lenrec real words
  ! total number of records to be written in this block
  ! .. scalar arguments ..
  integer irec,lenrec,lunit,nwords,nextr,dimabs
  ! .. array arguments ..
  real(kr) bc(nwords)
  ! .. local scalars ..
  integer j,jl,ju,nr,nrecs
  real, dimension(:), allocatable :: dummy
  !
  allocate(dummy(nwords))
  dummy(:nwords)=real(bc(:nwords)) ! double to single precision
  nrecs = (nwords+lenrec-1)/lenrec
  ju = 0
  do nr = 1,nrecs
    jl = ju + 1
    ju = ju + lenrec
    if(ju>nwords) ju = nwords
    irec = irec + 1
    write(lunit,rec=irec) (dummy(j),j=jl,ju)
  enddo
  deallocate(dummy)
  return
  end subroutine ecwr4
  !
  subroutine find (mat,mf,mt,ntape,found)
  !     based on njoy routine findf but search forward only
  !     find specified section on an endf/b format tape.
  !     if ntape lt 0, it is assumed to be in binary mode.
  !     if mt=0, find the first section in mat,mf.
  !     if mf=0, find the first file in mat.
  !     routine searches down and leaves tape positioned
  !     to read the *head* card for the section requested.
  integer mat,mf,mt,ntape,itape
  logical found,there
  !
  found=.false.
  ! ***read first card.
  itape=iabs(ntape)
  inquire(itape,exist=there)
  100 continue
  if(ntape>0) then
    read(itape,'(66x,i4,i2,i3)') math,mfh,mth
  else
    read(itape) math,mfh,mth
  endif
  ! ***test for mat
  if((math<=0) .or. (math>mat)) then
    call skiprz(ntape,-1)
    return
    if(math<mat) go to 100
  else
    ! ***test for file
    if(mf==0) go to 300
    if(mfh==mf) go to 130
    if(mfh<mf) go to 100
    call skiprz(ntape,-1)
    return
    ! ***test for section
    130 if(mt==0) go to 300
    if(mth==mt) go to 300
    if(mth<mt) go to 100
    call skiprz(ntape,-1)
    return
    ! ***desired section has been found
    ! ***backspace to before head card
    300 call skiprz(ntape,-1)
    found= .true.
    return
  endif
  end subroutine find
  !
  subroutine search(reac,number,list,isize)
  integer number,isize,in
  character*16 reac,list(isize)
  !----
  ! search for reaction on reference table
  !----
  do in=1,isize
    if(reac==list(in))then
      number=in
      return
    endif
  enddo
  !----
  ! reaction not found
  !----
  write(nsyso,'(/21h available reactions:/(5x,8a17))') (list(in),in=1,isize)
  write(hsmg,'(9hreaction ,a,10h not found)') trim(reac)
  call error('search',hsmg,' ')
  end subroutine search
end module eccom
