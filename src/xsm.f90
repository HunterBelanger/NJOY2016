module kdim
   !
   !-----------------------------------------------------------------------
   !
   ! kdi file support module.
   !
   ! this module provide the possibility to use a direct access (da) file
   ! to simulate memory like support. maxbuf static buffers are used.
   !
   !  contains:
   !  . kdiget    extract data from the kdi file.
   !  . kdiput    save data into the kdi file.
   !  . kdicl     terminate operations on the kdi file.
   !
   !  remark : file opening must be done using the open statement
   !           provided in the Fortran-90 standard:
   !           open(ifile,access='direct',recl=lrecl*iunit) where iunit
   !           is the number of memory units in a single precision word.
   !
   ! saved variables:
   !  nbuf           : number of buffers used. nbuf<=maxbuf
   !  ibufs          : index of the last buffer used.
   !  nrecs(ibuf)    : relative record number corresponding to the buffer.
   !  ifiles(ibuf)   : unit number corresponding to the buffer.
   !  iage(ibuf)     : age of the buffer.
   !  modif(ibuf)    : modification flag for the buffer.
   !  record(1,ibuf) : static buffer. dimension record(lrecl,maxbuf)
   !
   ! Author: Alain Hebert (original version -- 1994)
   !
   !-----------------------------------------------------------------------
   !
   implicit integer(a-z)
   private
   public :: kdiget,kdiput,kdicl
   integer,parameter :: maxbuf=10,lrecl=1024
   integer,save :: nbuf=0
   integer,save :: ibufs,nrecs(maxbuf),ifiles(maxbuf),iage(maxbuf)
   logical,save :: modif(maxbuf)
   integer,save :: record(lrecl,maxbuf)
   !
   interface kdiget
      !
      !----------------------------------------------------------------------
      !
      ! extract data from the kdi file.
      !
      ! input parameters:
      !  ifile  : unit number of the kdi file.
      !  iofset : offset, from start of kdi file where data is extracted.
      !  length : number of addressable units to read from file.
      !  irc    : return code (zero if no error is detected).
      !
      ! output parameter:
      !  data   : data to be copied from the kdi file.
      !
      !----------------------------------------------------------------------
      !
      module procedure kdiget_i,kdiget_r,kdiget_h,kdiget_d,kdiget_l,kdiget_c
   end interface
   !
   interface kdiput
      !
      !----------------------------------------------------------------------
      !
      ! save data into the kdi file.
      !
      ! input parameters:
      !  ifile  : unit number of the kdi file.
      !  data   : data to be copied from the kdi file.
      !  iofset : offset, from start of kdi file where data is written.
      !  length : number of addressable units to write into file.
      !  irc    : return code (zero if no error is detected).
      !
      !----------------------------------------------------------------------
      !
      module procedure kdiput_i,kdiput_r,kdiput_h,kdiput_d,kdiput_l,kdiput_c
   end interface
contains
   subroutine kdiget_i(ifile,data,iofset,length,irc)
   integer,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-90
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=record(n,ibuf)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_i
   !
   subroutine kdiget_r(ifile,data,iofset,length,irc)
   real,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-91
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=transfer(record(n,ibuf),0.0)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_r
   !
   subroutine kdiget_h(ifile,data,iofset,length,irc)
   character(len=*),dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-92
      return
   else if(len(data(1))/=4) then
      irc=-93
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=transfer(record(n,ibuf),'aaaa')
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_h
   !
   subroutine kdiget_d(ifile,data,iofset,length,irc)
   double precision,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-94
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            data(l)=transfer(record(n-1:n,ibuf),0.0d0)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_d
   !
   subroutine kdiget_l(ifile,data,iofset,length,irc)
   logical,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-95
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            data(l)=transfer(record(n,ibuf),.true.)
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_l
   !
   subroutine kdiget_c(ifile,data,iofset,length,irc)
   complex,dimension(:) :: data
   !
   if(lbound(data,1)/=1) then
      irc=-94
      return
   endif
   irc=0
   kbuf=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            data(l)=transfer(record(n-1:n,ibuf),(0.0,0.0))
         enddo
         do i=1,nbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=20,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   ! ----------------------------------------------------------------
   read(ifile,rec=nrec,err=20,iostat=irc) (record(i,ibufs),i=1,lrecl)
   ! ----------------------------------------------------------------
   go to 10
   20 return
   end subroutine kdiget_c
   !
   subroutine kdiput_i(ifile,data,iofset,length,irc)
   integer,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-97
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=data(l)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_i
   !
   subroutine kdiput_r(ifile,data,iofset,length,irc)
   real,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-98
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=transfer(data(l),0)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_r
   !
   subroutine kdiput_h(ifile,data,iofset,length,irc)
   character(len=*),dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-99
      return
   else if(len(data(1))/=4) then
      irc=-100
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=transfer(data(l),0)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_h
   !
   subroutine kdiput_d(ifile,data,iofset,length,irc)
   double precision,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-101
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            record(n-1:n,ibuf)=transfer(data(l),(/0/))
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_d
   !
   subroutine kdiput_l(ifile,data,iofset,length,irc)
   logical,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-102
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=lmin,ngro
            n=n+1
            record(n,ibuf)=transfer(data(l),0)
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_l
   !
   subroutine kdiput_c(ifile,data,iofset,length,irc)
   complex,dimension(:) :: data
   character(len=72) :: text72
   !
   if(lbound(data,1)/=1) then
      irc=-103
      return
   endif
   kbuf=0
   irc=0
   nrec=1+iofset/lrecl
   n=mod(iofset,lrecl)
   lmin=1
   10 ngro=min(lrecl+lmin-n-1,length)
   kdm=-1
   do iiii=1,nbuf
      ibuf=1+mod(iiii+ibufs-2,nbuf)
      if((ifile==ifiles(ibuf)).and.(nrec==nrecs(ibuf))) then
         do l=(lmin+1)/2,(ngro+1)/2
            n=n+2
            record(n-1:n,ibuf)=transfer(data(l),(/0/))
         enddo
         do i=1,maxbuf
            if(i/=ibuf) iage(i)=iage(i)+1
         enddo
         modif(ibuf)=.true.
         ibufs=ibuf
         if(ngro==length) return
         nrec=nrec+1
         n=0
         lmin=ngro+1
         ibufs=ibuf+1
         go to 10
      endif
      if(iage(ibuf)>kdm) then
         kdm=iage(ibuf)
         kbuf=ibuf
      endif
   enddo
   if(nbuf.lt.maxbuf) then
      nbuf=nbuf+1
      ibufs=nbuf
   else
      ibufs=kbuf
      if(modif(ibufs)) then
         ! -------------------------------------------------------
         write(ifiles(ibufs),rec=nrecs(ibufs),err=30,iostat=irc) &
         & (record(i,ibufs),i=1,lrecl)
         ! -------------------------------------------------------
      endif
   endif
   nrecs(ibufs)=nrec
   ifiles(ibufs)=ifile
   iage(ibufs)=0
   modif(ibufs)=.false.
   if((n/=0).or.(lrecl+lmin-n-1>length)) then
      ! ---------------------------------------------------------------
      read(ifile,rec=nrec,err=20,iostat=ii) (record(i,ibufs),i=1,lrecl)
      ! ---------------------------------------------------------------
   endif
   go to 10
   20 inquire(ifile,name=text72,recl=nn)
   close(ifile,status='keep')
   open(ifile,file=text72,access='direct',recl=nn,status='old')
   go to 10
   30 return
   end subroutine kdiput_c
   !
   subroutine kdicl(ifile,istatu,jrc)
   kbuf=0
   if(istatu==1) then
      do ibuf=1,nbuf
         if(modif(ibuf).and.(ifile==ifiles(ibuf))) then
            ! -----------------------------------------------
            write(ifile,rec=nrecs(ibuf),err=130,iostat=jrc) &
            & (record(i,ibuf),i=1,lrecl)
            ! -----------------------------------------------
            modif(ibuf)=.false.
         endif
         if(ifile==ifiles(ibuf)) ifiles(ibuf)=0
      enddo
      close(ifile,err=130,status='keep',iostat=jrc)
   else if(istatu==2) then
      do ibuf=1,nbuf
         if(ifile==ifiles(ibuf)) then
            iage(ibuf)=999999
            modif(ibuf)=.false.
            ifiles(ibuf)=0
         endif
      enddo
      close(ifile,err=130,status='delete',iostat=jrc)
   else
      jrc=-999
   endif
   130 return
   end subroutine kdicl
end module kdim

module filem
   !
   !-----------------------------------------------------------------------
   !
   ! file allocation/release support module.
   !
   ! allocate and release file units associated to a given file name.
   ! direct access (kdi), sequential (formatted or not) and direct
   ! access (da) files are permitted
   !
   !  contains:
   !  . kdropn    open file and allocate unit number.
   !  . kdrcls    close file and release unit number.
   !  . kdrfst    return information about a file.
   !
   ! saved variables:
   !  nbtape   : number of units available                     i
   !  nreser   : number of reserved names                      i
   !  ndummy   : number of dummy units                         i
   !  creser   : reserved names                     c(nreser)*12
   !  cdummy   : names of dummy files               c(ndummy)*12
   !  ctapen   : names allocated to units           c(nbtape)*16
   !  itapet   : type of files for all units           i(nbtape)
   !             =0  no type associated to file
   !             =iutype  file of type iutype opened
   !             =-iutype  file of type iutype closed
   !  ilenda   : length of da file                     i(nbtape)
   !             =0       not a da file
   !             =lrda    length of da file
   !
   ! Author: Alain Hebert (original version -- 1994)
   !
   !-----------------------------------------------------------------------
   !
   use kdim
   private
   public :: kdropn,kdrcls,kdrfst
   integer,parameter :: nbtape=63,nreser=2,ndummy=4,nb=nbtape-6,lrecld=1024
   character(len=16),save,dimension(nbtape) :: ctapen=               &
    (/ ('        ',i=1,4),'ft05f001','ft06f001',                     &
    ('        ',i=1,nb) /)
   character(len=12),save,dimension(ndummy) :: cdummy=               &
    (/ 'dummyda','dummysq','dummyca','dummyin' /)
   character(len=12),save,dimension(nreser) :: creser=               &
    (/ 'ft05f001','ft06f001' /)
   character(len=22),save,dimension(ndummy) :: ctype=                &
    (/ 'direct access kdi     ','sequential unformatted',            &
    'sequential character  ','direct access da      ' /)
   integer,save,dimension(nbtape) :: itapet=0,ilenda=0
   !
   interface kdrfst
      !
      !-----------------------------------------------------------------------
      !
      ! return the unit allocated to file name or the file name allocated
      ! to unit. also return the iutype of that file.
      !
      ! input parameters:
      !  cuname/ifile : character file name or file unit number.
      !
      ! output parameters (optional):
      !  name  : character file name.
      !  unit  : file unit number.
      !  ldra  : number of words in da file. returned only if iutype=4.
      !  type  : error flag or iutype of the file.
      !          > 0  file iutype
      !          =-2  file not oppened
      !          =-3  this file has been opened by routines other than kdropn
      !          =-4  file unit is reserved (5,6)
      !          =-5  illegal unit number
      !
      !-----------------------------------------------------------------------
      !
      module procedure kdrfst_h,kdrfst_i
   end interface
   !
contains
   !
   function kdropn(osname,iactio,iutype,lrda)
   !
   !-----------------------------------------------------------------------
   !
   ! open file and allocate unit number. allocate unit number to file
   ! name if unit is already opened, returns unit number.
   !
   ! input parameters:
   !     osname   : filename                                     c
   !                if osname=' ', use a default name
   !     iactio   : action on file                               i
   !                =0 to create a new file
   !                =1 to access and modify an existing file
   !                =2 to access an existing file in read-only mode
   !                =3 unknown
   !     iutype   : file type                                    i
   !                =1  kdi file type
   !                =2  sequential unformatted
   !                =3  sequential formatted
   !                =4  direct access (da) unformatted file
   !     ldra     : number of words in da file                   i
   !                required for  iutype=4 only
   !
   ! output parameter:
   !     kdropn   : unit number/error status                     i
   !                >0  unit number allocated
   !               =-1  no more unit available
   !               =-2  file type requested inconsistent with file type
   !                    of this file
   !               =-3  this file has been already opened
   !               =-4  file name is reserved
   !                      sysin,ft06f00,ft07f00
   !               =-5  illegal file type 0 < iutype < 4
   !               =-6  error on open of kdi file
   !               =-7  error on open of unformatted sequential file
   !               =-8  error on open of formatted sequential file
   !               =-9  error on open of direct access file
   !               =-10 invalid number of word in direct access file
   !                    lrda must be > 0 and given iif iutype=4
   !               =-11 da file length inconsistent with previous length
   !                    for this file
   !
   !-----------------------------------------------------------------------
   !
   character(len=*) :: osname
   integer,optional :: lrda
   character(len=16) :: crdnam
   character(len=12) :: cform,cstatu
   logical :: lfilop
   integer,dimension(:),allocatable :: data
   !----
   !  check if iutype is valid
   !----
   kdropn=-9999
   if((iutype>4).or.(iutype<1)) then
     kdropn=-5
     return
   endif
   !----
   !  check if lrda is valid
   !----
   if(present(lrda)) then
     if((iutype/=4).or.(lrda<1)) then
       kdropn=-10
       return
     endif
   else
     if(iutype==4) then
       kdropn=-10
       return
     endif
   endif
!----
!  check if file name not forbidden
!----
   do ireser=1,nreser
     if(osname==creser(ireser)) then
       kdropn=-4
       return
     endif
   enddo
   !----
   !  check for dummy file name/allocate dummy file name if requested
   !----
   do idummy=1,ndummy
     if(osname==cdummy(idummy)) then
       if(idummy/=iutype) then
         kdropn=-2
         return
       endif
     endif
   enddo
   if(osname==' ') then
      crdnam=cdummy(iutype)
   else
      crdnam=osname
   endif
   !----
   !  check if file opened/permitted
   !----
   inquire(file=crdnam,opened=lfilop,number=iboucl)
   if(lfilop) then
     kdropn=-3
     return
   endif
   !----
   !  look for a new unit number for file crdnam
   !  look if file name is in table of previously allocated units
   !----
   do jboucl=nbtape,1,-1
     if(ctapen(jboucl)==crdnam) then
       iboucl=jboucl
       if(iutype/=-itapet(iboucl)) then
         kdropn=-2
         return
       else if(iutype==4) then
         kdropn=-11
         return
         if(lrda/=ilenda(iboucl)) then
           kdropn=-11
           return
         endif
       endif
       go to 121
     endif
   enddo
   !----
   !  look for never allocated unit location
   !----
   do jboucl=nbtape,1,-1
     if(ctapen(jboucl)==' ') then
       iboucl=jboucl
       inquire(unit=iboucl,opened=lfilop)
       if(.not.lfilop) go to 121
     endif
   enddo
   !----
   !  look for a unit previously allocated but released.
   !----
   do jboucl=nbtape,1,-1
     if(itapet(jboucl)<0) then
       iboucl=jboucl
       go to 121
     endif
   enddo
   !----
   !  error - no unit number available
   !----
   kdropn=-1
   return
   !
   121 lrda2=0
   if(iutype==1) then
      !----
      !  open kdi file
      !----
     kdropn=-6
     lrda2=lrecld
     cform='unformatted'
   else if(iutype==2) then
     !----
     !  open sequential unformatted file
     !----
     kdropn=-7
     cform='unformatted'
   else if(iutype==3) then
     !----
     !  open sequential formatted file
     !----
      kdropn=-8
      cform='formatted'
   else if(iutype==4) then
     !----
     !  open da file
     !----
     kdropn=-9
     lrda2=lrda
     cform='unformatted'
   endif
   if(iactio==0) then
     cstatu='new'
   else if(iactio==1) then
     cstatu='old'
   else if(iactio==2) then
     cstatu='old'
   else
     cstatu='unknown'
   endif
   if(((iutype==1).or.(iutype==4)).and.(iactio==2)) then
     allocate(data(lrda2))
     data=0
     inquire(iolength=lrecl) (data(i),i=1,lrda2)
     deallocate(data)
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='direct',recl=lrecl,status=cstatu,action='read')
   else if((iutype==1).or.(iutype==4)) then
     allocate(data(lrda2))
     data=0
     inquire(iolength=lrecl) (data(i),i=1,lrda2)
     deallocate(data)
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='direct',recl=lrecl,status=cstatu)
   else if(((iutype==2).or.(iutype==3)).and.(iactio==2)) then
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='sequential',status=cstatu,action='read')
     rewind(iboucl)
   else if((iutype==2).or.(iutype==3)) then
     open(unit=iboucl,file=crdnam,err=9000,iostat=iercod,form=cform, &
     access='sequential',status=cstatu)
     rewind(iboucl)
   endif
   kdropn=iboucl
   ctapen(iboucl)=crdnam
   itapet(iboucl)=iutype
   if(iutype==4) ilenda(iboucl)=lrda
   return
   9000 write(6,8000) crdnam, ctype(iutype),iercod
   !----
   !  error format
   !----
   8000 format('1',5x,'error in opening of file in kdropn'/ &
   6x,'file name=',a7/6x,'file type=',a22/                  &
   6x,'   iercod=',i7)
   return
   end function kdropn
   !
   function kdrcls(itapno,iactio)
   !
   !-----------------------------------------------------------------------
   !
   ! close file and release unit number.
   !
   ! input parameters:
   !  itapno   : unit number                                  i
   !            =0    close all units
   !             > 0    unit to close
   !  iactio   : action on file                               i
   !            =1    to keep the file;
   !            =2    to delete the file.
   !
   ! output parameter:
   !  kdrcls   : error status                                 i
   !            = 0  unit closed
   !            =-2  file not oppened
   !            =-3  this file has been opened by routines
   !                   other than kdropn
   !            =-4  file unit is reserved (5,6)
   !            =-5  illegal unit number
   !            =-6  error on close of kdi file
   !            =-7  error on close of unformatted
   !                   sequential file
   !            =-8  error on close of formatted sequential file
   !            =-9  error on close of da file
   !            =-10  invalid close action iactio=1,2 permitted only
   !
   !-----------------------------------------------------------------------
   !
   logical :: lfilop
   !
   if(itapno==0) then
   !----
   !  unit=0 close all opened units
   !----
     idebtp=1
     ifintp=nbtape
     go to 200
   else if((itapno>nbtape).or.(itapno<0)) then
   !----
   !  invalid unit number
   !----
     kdrcls=-5
     return
   endif
   inquire(unit=itapno,opened=lfilop)
   if(lfilop) then
     if(itapet(itapno)==0) then
   !----
   !  reserved unit
   !----
       kdrcls=-4
       return
     endif
   !----
   !  unit permitted
   !----
     idebtp=itapno
     ifintp=itapno
     go to 200
   else
   !----
   !  unit not opened
   !----
     kdrcls=-2
     return
   endif
   200 do iboucl=idebtp,ifintp
     if(itapet(iboucl)==1) then
   !----
   !  close kdi type file
   !----
       kdrcls=-6
       call kdicl(iboucl,iactio,iercod)
       if(iercod/=0) go to 9001
       if(iactio==1) then
         itapet(iboucl)=-itapet(iboucl)
       else if(iactio==2) then
         itapet(iboucl)=0
         ctapen(iboucl)=' '
       else
         kdrcls=-10
       endif
     else if(itapet(iboucl)>0) then
   !----
   !  close other types of file
   !----
       kdrcls=-(5+itapet(iboucl))
       if(iactio==1) then
         close(iboucl,iostat=iercod,status='keep',err=9001)
         itapet(iboucl)=-itapet(iboucl)
       else if(iactio==2) then
         close(iboucl,iostat=iercod,status='delete',err=9001)
         itapet(iboucl)=0
         ctapen(iboucl)=' '
         ilenda(iboucl)=0
       else
         kdrcls=-10
       endif
       if(iercod/=0) go to 9001
     endif
   enddo
   kdrcls=0
   return
   9001 write(6,8001) iboucl,ctapen(iboucl),ctype(itapet(iboucl)),iercod
   !----
   !  error format
   !----
   8001 format('1',5x,'error in close of file in kdropn'/            &
   6x,'unit nb. =',i10/6x,'file name=',a7/6x,'file type=',a22/       &
   6x,'   iercod=',i10 )
   return
   end function kdrcls
   !
   subroutine kdrfst_h(cuname,name,unit,type,lrda)
   character(len=*) :: cuname
   character(len=*),optional :: name
   integer,optional :: unit,type,lrda
   character(len=16) :: crdnam
   logical :: lfilop
   !----
   !  return file number and type
   !----
   if(present(name)) name=cuname
   if(present(lrda)) lrda=0
   crdnam=cuname
   inquire(file=crdnam,opened=lfilop,number=iunit)
   if(present(unit)) unit=iunit
   if(lfilop) then
      if((iunit<1).or.(iunit>nbtape)) then
        if(present(type)) type=-5
      else if((ctapen(iunit)==' ').or.(itapet(iunit)<=0)) then
        if(present(type)) type=-3
      else
        if(present(type)) type=itapet(iunit)
        if((itapet(iunit)==4).and.present(lrda)) lrda=ilenda(iunit)
      endif
   else
      if(present(type)) type=-2
   endif
   return
   end subroutine kdrfst_h
   !
   subroutine kdrfst_i(ifile,name,unit,type,lrda)
   character(len=*),optional :: name
   integer,optional :: unit,type,lrda
   logical :: lfilop
   !
   if(present(name)) name=' '
   if(present(unit)) unit=ifile
   if(present(lrda)) lrda=0
   !----
   !  illegal unit
   !----
   if((ifile<1).or.(ifile>nbtape)) then
       if(present(type)) type=-5
       return
   endif
   !----
   !  return file name and type
   !----
   inquire(ifile,opened=lfilop)
   if(lfilop) then
      if((ctapen(ifile)==' ').or.(itapet(ifile)<=0)) then
        if(present(type)) type=-3
      else
        if(present(name)) name=ctapen(ifile)
        if(present(type)) type=itapet(ifile)
        if((itapet(ifile)==4).and.present(lrda)) lrda=ilenda(ifile)
      endif
   else
      if(present(type)) type=-2
   endif
   return
   end subroutine kdrfst_i
end module filem

module xsmm
   !
   !-----------------------------------------------------------------------
   !
   ! xsm file support module.
   !
   ! the xsm database results from the juxtaposition of a hierarchical
   ! logical structure into a direct access file with fixe length records.
   ! the direct access file is fortran-90 compatible and is managed by
   ! kdiget/put/cl. xsmop/put/get/inf/nod/cl entries provide a set of
   ! methods to access an xsm file.
   !
   ! the logical structure of an xsm file is made of a root directory fol-
   ! lowed by variable-length nodes containing the usefull information.
   ! each time a directory is full, an extent is automatically created at
   ! the end of the file, so that the total number of nodes in a direc-
   ! tory is only limited by the maximum size of the direct access file.
   ! any node can contain a sub-directory in order to create a hierar-
   ! chical structure.
   !input parameters:
   !  imp    : type of access.   0: new file mode;
   !                             1: modification mode;
   !                             2: read only mode.
   !  impx   : if impx=0, we suppress printing on xsmop.
   !  namp   : character*12 name of the current block.
   !  nammy  : charecter*12 name of the active directory.
   !  ilong  : number of information elements stored in the current block.
   !  itype  : type of information elements stored in the current block.
   !           0: directory                1: integer
   !           2: single precision         3: character*4
   !           4: double precision         5: logical
   !           6: complex                  7: undefined
   !  data1 :  information elements. dimension data1(ilong)
   !  iact  :  type of movement in the hierarchical structure.
   !           0: return to the root directory;
   !           1: move to a son directory;
   !           2: move to the parent directory.
   !  istatu : =1 to keep the xsm file at close ; =2 to destroy it.
   !
   ! entry points:
   !  entry xsmop  : open a xsm file.
   !  entry xsmput : copy a block from memory into the xsm file.
   !  entry xsmget : copy a block from the xsm file into memory.
   !  entry xsmlen : return the length and type of a block. return 0 if
   !                 the block doed not exists.
   !  entry xsmsix : move in the hierarchical structure.
   !  entry xsmnxt : find the name of the next block stored in the active
   !                 directory. if namp=' ' at input, find any name for
   !                 any block stored in this directory.
   !  entry xsmcl  : close the xsm file.
   !
   ! type(xsm_file):
   !  cmt    : character*12 names of the blocks stored in the active
   !           directory extent.
   !  idir   : offset of the active directory extent.
   !  iofmax : dimension of vectors iofset,jlong jtype, and nmt
   !  nmt    : number of blocks stored in the active directory extent.
   !  link   : offset of the next directory extent. the last extent is
   !           always linked to the first one.
   !  iroot  : offset of an extent that belong to the parent directory.
   !           equal to -1 for the root directory.
   !  myname : name of the active directory. myname='/' for the root level.
   !  iofset : offset table for the active directory extent (position of
   !           the first element of each block in the xsm file).
   !  jlong  : number of information elements contained in each block of
   !           the active directory extent.
   !  jtype  : type of information elements contained in each block of the
   !           active directory extent.
   !
   ! Author: Alain Hebert (original version -- 1995)
   !
   !-----------------------------------------------------------------------
   !
   use kdim
   use filem
   implicit none
   private
   public :: xsm_file,xsmop,xsmget,xsmput,xsmcl,xsmsix,xsmlib,xsmexp,xsmlen
   integer,parameter :: iofmax=30,iprim=3,iwrd=3,klong=5+iwrd+(3+iwrd)*iofmax
   ! iprim=address of the root directory.
   ! klong=length of a directory extent in memory.
   type xsm_file
      private
      character(len=12) :: myname,cmt(iofmax)
      logical modif
      integer :: ifile,nmt,link,iroot,iofset(iofmax),jlong(iofmax), &
      & jtype(iofmax),idir,ioftop,impfil
   end type xsm_file
   character(len=131) :: hsmg,tex131
   character(len=8) :: circ
   integer,dimension(5) :: itemp
   !
   interface xsmget
     !
     !----------------------------------------------------------------------
     !
     ! copy a node of data or a directory from a xsm file into memory.
     !
     ! input parameters:
     !   pfxsm : xsm_file object handle.
     !    namp : character*lnod name of the current node.
     !
     ! output parameter:
     !   data1 : information elements.
     !
     !-----------------------------------------------------------------------
     !
     module procedure xsmget_i1,xsmget_r1,xsmget_h1,xsmget_d1
   end interface
   !
   interface xsmput
     !
     !-----------------------------------------------------------------------
     !
     ! copy a node of data or a directory from memory into a xsm file.
     !
     ! input parameters:
     !   pfxsm : xsm_file object handle.
     !    namp : character*lnod name of the current node.
     !   data1 : information elements.
     !
     !-----------------------------------------------------------------------
     !
     module procedure xsmput_i1,xsmput_r1,xsmput_h1,xsmput_d1
   end interface
   !
contains
   !
   subroutine xsmop(pfxsm,fname,mode,impx)
   !-----------------------------------------------------------------------
   !
   ! open an existing or create a new xsm file.
   !
   ! input parameters:
   !    fname: character*lnod name of the xsm_file object. an xsm file may
   !           contains one or many xsm_file objects.
   !    mode : type of access.  =0: new file mode;
   !                            =1: modification mode;
   !                            =2: read only mode (default value).
   !    impx : if impx=0, we suppress printing on xsmop (optional).
   !
   ! output parameter:
   !  pfxsm : xsm_file object handle.
   !
   !-----------------------------------------------------------------------
   use util   ! provides error
   use mainio ! provides nsysi,contio,nsyso,nsyse
   !
   type(xsm_file),pointer :: pfxsm
   character(len=*) :: fname
   integer,optional :: mode,impx
   !
   integer :: imp,impy,irc
   character(len=4) :: header
   !
   ! open the xsm file.
   imp=2
   if(present(mode)) imp=mode
   impy=0
   if(present(impx)) impy=impx
   pfxsm%impfil=imp
   pfxsm%ifile=kdropn(fname,imp,1)
   if(pfxsm%ifile<=0) call error('xsmop','unable to open xsm file',' ')
   if(pfxsm%impfil.ge.1) then
      ! recover the root directory if the xsm file already exists.
      if(pfxsm%ifile<=0) go to 160
      call kdiget(pfxsm%ifile,itemp,0,3,irc)
      if(irc.ne.0) go to 140
      write (header,'(a4)') itemp(1)
      if(header.ne.'$xsm') then
         call error('xsmop','unable to import xsm file',' ')
      endif
      pfxsm%ioftop=itemp(2)
      pfxsm%idir=itemp(3)
      call xsmdir(1,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link,pfxsm%iroot, &
      & pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir)
      pfxsm%modif=.false.
      if(impy.gt.0) write (nsyso,185) fname,pfxsm%ioftop,pfxsm%myname
   else
      ! the xsm file is new.
      pfxsm%ioftop=iprim+klong
      pfxsm%idir=iprim
      pfxsm%iroot=-1
      pfxsm%nmt=0
      pfxsm%link=pfxsm%idir
      pfxsm%myname='/'
      pfxsm%modif=.true.
      header='$xsm'
      read (header,'(a4)') itemp(1)
      itemp(2)=pfxsm%ioftop
      itemp(3)=pfxsm%idir
      call kdiput(pfxsm%ifile,itemp,0,3,irc)
      if(irc.ne.0) go to 150
   endif
   return
   140 write(circ,'(i8)') irc
   call error('xsmop','(kdiget) iostat='//circ,' ')
   150 write(circ,'(i8)') irc
   call error('xsmop','(kdiput) iostat='//circ,' ')
   160 write(hsmg,'("unable to open xsm file ''",a,"'' (ier =",i3,&
   &").")') pfxsm%myname,pfxsm%ifile
   call error('xsmop',hsmg,' ')
   185 format (/35h xsmop: xsm file recovery. file = ',a,1h'/27x, &
   & 28hhighest attainable address =,i10/27x,20hactive directory = ',&
   & a12,1h'/)
   end subroutine xsmop
   !
   subroutine xsmput_i1(pfxsm,namp,data1)
   use util   ! provides error
   type(xsm_file),pointer :: pfxsm
   integer,dimension(:) :: data1
   character(len=*) :: namp
   character(len=12) :: namt
   integer :: ilong,itype,iii,irc,ipos
   !
   ! copy a node from memory into the xsm file.
   if(pfxsm%impfil.eq.2) then
      call error('xsmput',' forbidden operation in read-only mode.',' ')
   endif
   namt=namp
   ilong=size(data1)
   itype=1
   call xsmrep(namt,2,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block already exists. replace it.
      pfxsm%modif=.true.
      if(ilong.gt.pfxsm%jlong(iii)) then
         pfxsm%iofset(iii)=pfxsm%ioftop
         pfxsm%ioftop=pfxsm%ioftop+ilong
      endif
      pfxsm%jlong(iii)=ilong
      pfxsm%jtype(iii)=itype
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(iii),ilong,irc)
      if(irc.ne.0) go to 150
      return
   else
      ! the block does not already exists. create it.
      if(pfxsm%nmt.eq.iofmax) then
         ! the active directory is full. create an extent.
         ipos=pfxsm%link
         pfxsm%link=pfxsm%ioftop
         if(pfxsm%modif) then
            call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
            & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype, &
            & pfxsm%idir)
         else
            itemp(1)=pfxsm%link
            call kdiput(pfxsm%ifile,itemp,pfxsm%idir+3,1,irc)
            if(irc.ne.0) go to 150
         endif
         pfxsm%idir=pfxsm%link
         pfxsm%link=ipos
         pfxsm%ioftop=pfxsm%ioftop+klong
         pfxsm%nmt=0
      endif
      pfxsm%nmt=pfxsm%nmt+1
      pfxsm%cmt(pfxsm%nmt)=namt
      pfxsm%jlong(pfxsm%nmt)=ilong
      pfxsm%jtype(pfxsm%nmt)=itype
      pfxsm%iofset(pfxsm%nmt)=pfxsm%ioftop
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(pfxsm%nmt),ilong,irc)
      if(irc.ne.0) go to 150
      pfxsm%modif=.true.
      pfxsm%ioftop=pfxsm%ioftop+ilong
   endif
   return
   150 write(circ,'(i8)') irc
   call error('xsmput_i1','(kdiput) iostat='//circ,' ')
   end subroutine xsmput_i1
   !
   subroutine xsmput_r1(pfxsm,namp,data1)
   use util   ! provides error
   type(xsm_file),pointer :: pfxsm
   real,dimension(:) :: data1
   character(len=*) :: namp
   character(len=12) :: namt
   integer :: ilong,itype,iii,irc,ipos
   !
   ! copy a node from memory into the xsm file.
   if(pfxsm%impfil.eq.2) then
      call error('xsmput',' forbidden operation in read-only mode.',' ')
   endif
   namt=namp
   ilong=size(data1)
   itype=2
   call xsmrep(namt,2,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block already exists. replace it.
      pfxsm%modif=.true.
      if(ilong.gt.pfxsm%jlong(iii)) then
         pfxsm%iofset(iii)=pfxsm%ioftop
         pfxsm%ioftop=pfxsm%ioftop+ilong
      endif
      pfxsm%jlong(iii)=ilong
      pfxsm%jtype(iii)=itype
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(iii),ilong,irc)
      if(irc.ne.0) go to 150
      return
   else
      ! the block does not already exists. create it.
      if(pfxsm%nmt.eq.iofmax) then
         ! the active directory is full. create an extent.
         ipos=pfxsm%link
         pfxsm%link=pfxsm%ioftop
         if(pfxsm%modif) then
            call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link,& 
            & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype, &
            & pfxsm%idir)
         else
            itemp(1)=pfxsm%link
            call kdiput(pfxsm%ifile,itemp,pfxsm%idir+3,1,irc)
            if(irc.ne.0) go to 150
         endif
         pfxsm%idir=pfxsm%link
         pfxsm%link=ipos
         pfxsm%ioftop=pfxsm%ioftop+klong
         pfxsm%nmt=0
      endif
      pfxsm%nmt=pfxsm%nmt+1
      pfxsm%cmt(pfxsm%nmt)=namt
      pfxsm%jlong(pfxsm%nmt)=ilong
      pfxsm%jtype(pfxsm%nmt)=itype
      pfxsm%iofset(pfxsm%nmt)=pfxsm%ioftop
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(pfxsm%nmt),ilong,irc)
      if(irc.ne.0) go to 150
      pfxsm%modif=.true.
      pfxsm%ioftop=pfxsm%ioftop+ilong
   endif
   return
   150 write(circ,'(i8)') irc
   call error('xsmput_i1','(kdiput) iostat='//circ,' ')
   end subroutine xsmput_r1
   !
   subroutine xsmput_h1(pfxsm,namp,data1)
   use util   ! provides error
   type(xsm_file),pointer :: pfxsm
   character(len=4),dimension(:) :: data1
   character(len=*) :: namp
   character(len=12) :: namt
   integer :: ilong,itype,iii,irc,ipos
   !
   ! copy a node from memory into the xsm file.
   if(pfxsm%impfil.eq.2) then
      call error('xsmput',' forbidden operation in read-only mode.',' ')
   endif
   namt=namp
   ilong=size(data1)
   itype=3
   call xsmrep(namt,2,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block already exists. replace it.
      pfxsm%modif=.true.
      if(ilong.gt.pfxsm%jlong(iii)) then
         pfxsm%iofset(iii)=pfxsm%ioftop
         pfxsm%ioftop=pfxsm%ioftop+ilong
      endif
      pfxsm%jlong(iii)=ilong
      pfxsm%jtype(iii)=itype
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(iii),ilong,irc)
      if(irc.ne.0) go to 150
      return
   else
      ! the block does not already exists. create it.
      if(pfxsm%nmt.eq.iofmax) then
         ! the active directory is full. create an extent.
         ipos=pfxsm%link
         pfxsm%link=pfxsm%ioftop
         if(pfxsm%modif) then
            call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
            & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype, &
            & pfxsm%idir)
         else
            itemp(1)=pfxsm%link
            call kdiput(pfxsm%ifile,itemp,pfxsm%idir+3,1,irc)
            if(irc.ne.0) go to 150
         endif
         pfxsm%idir=pfxsm%link
         pfxsm%link=ipos
         pfxsm%ioftop=pfxsm%ioftop+klong
         pfxsm%nmt=0
      endif
      pfxsm%nmt=pfxsm%nmt+1
      pfxsm%cmt(pfxsm%nmt)=namt
      pfxsm%jlong(pfxsm%nmt)=ilong
      pfxsm%jtype(pfxsm%nmt)=itype
      pfxsm%iofset(pfxsm%nmt)=pfxsm%ioftop
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(pfxsm%nmt),ilong,irc)
      if(irc.ne.0) go to 150
      pfxsm%modif=.true.
      pfxsm%ioftop=pfxsm%ioftop+ilong
   endif
   return
   150 write(circ,'(i8)') irc
   call error('xsmput_i1','(kdiput) iostat='//circ,' ')
   end subroutine xsmput_h1
   !
   subroutine xsmput_d1(pfxsm,namp,data1)
   use util   ! provides error
   type(xsm_file),pointer :: pfxsm
   double precision,dimension(:) :: data1
   character(len=*) :: namp
   character(len=12) :: namt
   integer :: ilong,itype,iii,irc,ipos
   !
   ! copy a node from memory into the xsm file.
   if(pfxsm%impfil.eq.2) then
      call error('xsmput',' forbidden operation in read-only mode.',' ')
   endif
   namt=namp
   ilong=2*size(data1)
   itype=4
   call xsmrep(namt,2,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block already exists. replace it.
      pfxsm%modif=.true.
      if(ilong.gt.pfxsm%jlong(iii)) then
         pfxsm%iofset(iii)=pfxsm%ioftop
         pfxsm%ioftop=pfxsm%ioftop+ilong
      endif
      pfxsm%jlong(iii)=ilong
      pfxsm%jtype(iii)=itype
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(iii),ilong,irc)
      if(irc.ne.0) go to 150
      return
   else
      ! the block does not already exists. create it.
      if(pfxsm%nmt.eq.iofmax) then
         ! the active directory is full. create an extent.
         ipos=pfxsm%link
         pfxsm%link=pfxsm%ioftop
         if(pfxsm%modif) then
            call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
            & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype, &
            & pfxsm%idir)
         else
            itemp(1)=pfxsm%link
            call kdiput(pfxsm%ifile,itemp,pfxsm%idir+3,1,irc)
            if(irc.ne.0) go to 150
         endif
         pfxsm%idir=pfxsm%link
         pfxsm%link=ipos
         pfxsm%ioftop=pfxsm%ioftop+klong
         pfxsm%nmt=0
      endif
      pfxsm%nmt=pfxsm%nmt+1
      pfxsm%cmt(pfxsm%nmt)=namt
      pfxsm%jlong(pfxsm%nmt)=ilong
      pfxsm%jtype(pfxsm%nmt)=itype
      pfxsm%iofset(pfxsm%nmt)=pfxsm%ioftop
      call kdiput(pfxsm%ifile,data1,pfxsm%iofset(pfxsm%nmt),ilong,irc)
      if(irc.ne.0) go to 150
      pfxsm%modif=.true.
      pfxsm%ioftop=pfxsm%ioftop+ilong
   endif
   return
   150 write(circ,'(i8)') irc
   call error('xsmput_i1','(kdiput) iostat='//circ,' ')
   end subroutine xsmput_d1
   !
   subroutine xsmget_i1(pfxsm,namp,data2)
   use util   ! provides error
   !
   type(xsm_file),pointer :: pfxsm
   integer,dimension(:) :: data2
   character(len=*) :: namp
   !
   character(len=12) :: namt
   integer :: iii,irc
   !
   ! copy a block from xsm file into memory.
   namt=namp
   call xsmrep(namt,1,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block exists. recover it.
      call kdiget(pfxsm%ifile,data2,pfxsm%iofset(iii),pfxsm%jlong(iii),irc)
      if(irc.ne.0) go to 140
   else
   ! the block does not exists.
      go to 110
   endif
   return
   110 write(hsmg,190) namt,pfxsm%myname
   call error('xsmget_i1',hsmg,' ')
   140 write(circ,'(i8)') irc
   call error('xsmget_i1','(kdiget) iostat='//circ,' ')
   190 format (7hblock ',a12,19h' is not on level ',a12,18h' of the xsm file.)
   end subroutine xsmget_i1
   !
   subroutine xsmget_r1(pfxsm,namp,data2)
   use util   ! provides error
   !
   type(xsm_file),pointer :: pfxsm
   real,dimension(:) :: data2
   character(len=*) :: namp
   !
   character(len=12) :: namt
   integer :: iii,irc
   !
   ! copy a block from xsm file into memory.
   namt=namp
   call xsmrep(namt,1,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block exists. recover it.
      call kdiget(pfxsm%ifile,data2,pfxsm%iofset(iii),pfxsm%jlong(iii),irc)
      if(irc.ne.0) go to 140
   else
   ! the block does not exists.
      go to 110
   endif
   return
   110 write(hsmg,190) namt,pfxsm%myname
   call error('xsmget_i1',hsmg,' ')
   140 write(circ,'(i8)') irc
   call error('xsmget_i1','(kdiget) iostat='//circ,' ')
   190 format (7hblock ',a12,19h' is not on level ',a12,18h' of the xsm file.)
   end subroutine xsmget_r1
   !
   subroutine xsmget_h1(pfxsm,namp,data2)
   use util   ! provides error
   !
   type(xsm_file),pointer :: pfxsm
   character(len=4),dimension(:) :: data2
   character(len=*) :: namp
   !
   character(len=12) :: namt
   integer :: iii,irc
   !
   ! copy a block from xsm file into memory.
   namt=namp
   call xsmrep(namt,1,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block exists. recover it.
      call kdiget(pfxsm%ifile,data2,pfxsm%iofset(iii),pfxsm%jlong(iii),irc)
      if(irc.ne.0) go to 140
   else
   ! the block does not exists.
      go to 110
   endif
   return
   110 write(hsmg,190) namt,pfxsm%myname
   call error('xsmget_i1',hsmg,' ')
   140 write(circ,'(i8)') irc
   call error('xsmget_i1','(kdiget) iostat='//circ,' ')
   190 format (7hblock ',a12,19h' is not on level ',a12,18h' of the xsm file.)
   end subroutine xsmget_h1
   !
   subroutine xsmget_d1(pfxsm,namp,data2)
   use util   ! provides error
   !
   type(xsm_file),pointer :: pfxsm
   double precision,dimension(:) :: data2
   character(len=*) :: namp
   !
   character(len=12) :: namt
   integer :: iii,irc
   !
   ! copy a block from xsm file into memory.
   namt=namp
   call xsmrep(namt,1,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block exists. recover it.
      call kdiget(pfxsm%ifile,data2,pfxsm%iofset(iii),pfxsm%jlong(iii),irc)
      if(irc.ne.0) go to 140
   else
   ! the block does not exists.
      go to 110
   endif
   return
   110 write(hsmg,190) namt,pfxsm%myname
   call error('xsmget_i1',hsmg,' ')
   140 write(circ,'(i8)') irc
   call error('xsmget_i1','(kdiget) iostat='//circ,' ')
   190 format (7hblock ',a12,19h' is not on level ',a12,18h' of the xsm file.)
   end subroutine xsmget_d1
   !
   subroutine xsmlen(pfxsm,namp,ilong,itype)
   !----------------------------------------------------------------------
   !
   ! obtain length and type of a record
   !
   ! input parameters:
   !   pfxsm : xsm_file object handle.
   !    namp : character*lnod name of the current node.
   !
   ! output parameter:
   !   ilong : length of node
   !  ityxsm : type of information elements contained in node namp.
   !
   !----------------------------------------------------------------------
   type(xsm_file),pointer :: pfxsm
   character(len=*) :: namp
   integer :: ilong,itype
   character(len=12) :: namt
   integer :: iii
   !
   namt=namp
   ilong=0
   itype=7
   call xsmrep(namt,1,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.gt.0) then
      ! the block exists.
      ilong=pfxsm%jlong(iii)
      itype=pfxsm%jtype(iii)
   endif
   return
   end subroutine xsmlen
   !
   subroutine xsmsix(pfxsm,namp,iact)
   !----------------------------------------------------------------------
   !
   ! move in the hierarchical structure.
   !
   ! input parameters:
   !   pfxsm : xsm_file object handle.
   !    namp : character*lnod name of the daughter directory.
   !     ind : =1: move up to daughter directory namp; =2 move down.
   !
   !----------------------------------------------------------------------
   use util   ! provides error
   use mainio ! provides nsysi,contio,nsyso,nsyse
   type(xsm_file),pointer :: pfxsm
   character(len=*) :: namp
   integer :: iact
   character(len=12) :: namt
   integer :: iii,irc,ipos
   !
      if(iact.eq.0) then
         if(pfxsm%modif) then
            call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
            & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,&
            & pfxsm%idir)
            pfxsm%modif=.false.
         endif
         pfxsm%idir=iprim
      else if(iact.eq.2) then
         if(pfxsm%modif) then
            call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
            & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype, &
            & pfxsm%idir)
            pfxsm%modif=.false.
         endif
         pfxsm%idir=pfxsm%iroot
      else if(iact.eq.1) then
         namt=namp
         call xsmrep(namt,2,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt, &
         & pfxsm%link,pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong, &
         & pfxsm%jtype,pfxsm%idir,iii)
         if(iii.gt.0) then
            ! the son level already exists.
            if(pfxsm%jtype(iii).ne.0) go to 120
            if(pfxsm%modif) then
               call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
               & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong, &
               & pfxsm%jtype,pfxsm%idir)
               pfxsm%modif=.false.
            endif
            pfxsm%idir=pfxsm%iofset(iii)
         else
            ! the son level does not exists and is created.
            if(pfxsm%impfil.eq.2) then
               write(nsyso,'(/25h try to create directory ,a12, &
               & 16h from directory ,a12)') namt,pfxsm%myname
               call error('xsmsix',' forbidden operation in read-only mode.', &
               & ' ')
            endif
            if(pfxsm%nmt.eq.iofmax) then
               ! the active directory is full. create an extent.
               ipos=pfxsm%link
               pfxsm%link=pfxsm%ioftop
               if(pfxsm%modif) then
                  call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
                  & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong, &
                  & pfxsm%jtype,pfxsm%idir)
               else
                  itemp(1)=pfxsm%link
                  call kdiput(pfxsm%ifile,itemp,pfxsm%idir+3,1,irc)
                  if(irc.ne.0) go to 150
               endif
               pfxsm%idir=pfxsm%link
               pfxsm%link=ipos
               pfxsm%ioftop=pfxsm%ioftop+klong
               pfxsm%nmt=0
            endif
            pfxsm%nmt=pfxsm%nmt+1
            pfxsm%cmt(pfxsm%nmt)=namt
            pfxsm%jlong(pfxsm%nmt)=-1
            pfxsm%jtype(pfxsm%nmt)=0
            pfxsm%iofset(pfxsm%nmt)=pfxsm%ioftop
            pfxsm%ioftop=pfxsm%ioftop+klong
            call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
            & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype, &
            & pfxsm%idir)
            ! create a new directory.
            pfxsm%iroot=pfxsm%idir
            pfxsm%myname=namt
            pfxsm%idir=pfxsm%iofset(pfxsm%nmt)
            pfxsm%nmt=0
            pfxsm%link=pfxsm%idir
            pfxsm%modif=.true.
            return
         endif
      endif
      call xsmdir(1,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link,pfxsm%iroot, &
      & pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir)
   return
   120 write(tex131,200) namt
   call error('xsmsix',tex131,' ')
   150 write(circ,'(i8)') irc
   call error('xsmsix','(kdiput) iostat='//circ,' ')
   200 format (7hblock ',a12,28h' is not a directory of the ,9hxsm file.)
   end subroutine xsmsix
   !
   subroutine xsmnxt(pfxsm,namp,nammy)
   !----------------------------------------------------------------------
   !
   ! find the name of the next block stored into the active directory.
   !
   ! input parameters:
   !   pfxsm : xsm_file object handle.
   !    namp : character*lnod name of the current node.
   !
   ! output parameters:
   !    namp : character*lnod name of the next node.
   !   nammy : name of xsm object
   !
   !----------------------------------------------------------------------
   use util   ! provides error
   character(len=12) :: namt
   integer :: iii
   type(xsm_file),pointer :: pfxsm
   character(len=*) :: namp,nammy
   !
   nammy=pfxsm%myname
   namt=namp
   if(namt.eq.' ') then
      if(pfxsm%nmt.eq.0) return
      namp=pfxsm%cmt(1)
      return
   else if(pfxsm%nmt.eq.0) then
      go to 110
   endif
   call xsmrep(namt,1,pfxsm%ifile,pfxsm%modif,pfxsm%cmt,pfxsm%nmt,pfxsm%link, &
   & pfxsm%iroot,pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir, &
   & iii)
   if(iii.eq.0) then
      go to 110
   else if(iii+1.le.pfxsm%nmt) then
      namp=pfxsm%cmt(iii+1)
      return
   endif
   !
   ! switch to the next directory.
   if(pfxsm%idir.ne.pfxsm%link) then
      if(pfxsm%modif) then
         call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link,pfxsm%iroot, &
         & pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir)
         pfxsm%modif=.false.
      endif
      pfxsm%idir=pfxsm%link
      ! recover the next directory.
      call xsmdir(1,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link,pfxsm%iroot, &
      & pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir)
   endif
   namp=pfxsm%cmt(1)
   return
   110 write(tex131,190) namt,pfxsm%myname
   call error('xsmnxt',tex131,' ')
   190 format (7hblock ',a12,19h' is not on level ',a12,18h' of the xsm file.)
   end subroutine xsmnxt
   !
   subroutine xsmcl(pfxsm,istatu)
   !-----------------------------------------------------------------------
   !
   ! close the xsm file.
   !
   ! input parameters:
   !  pfxsm  : xsm_file object handle.
   !  istatu : =1 to keep the file at close ; =2 to destroy it.
   !
   !-----------------------------------------------------------------------
   use util   ! provides error
   type(xsm_file),pointer :: pfxsm
   integer :: irc,istatu,iii
   !
   if(pfxsm%modif) then
      if(pfxsm%impfil.eq.2) then
         call error('xsmcl',' forbidden operation in read-only mode.',' ')
      endif
      call xsmdir(2,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link,pfxsm%iroot, &
      & pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,pfxsm%idir)
   endif
   if(pfxsm%impfil.eq.2) then
      call xsmdir(1,pfxsm%ifile,pfxsm%cmt,pfxsm%nmt,pfxsm%link,pfxsm%iroot, &
      & pfxsm%myname,pfxsm%iofset,pfxsm%jlong,pfxsm%jtype,iprim)
   endif
   call kdiget(pfxsm%ifile,itemp,1,1,irc)
   if(irc.ne.0) go to 140
   iii=itemp(1)
   if(pfxsm%ioftop.gt.iii) then
      if(pfxsm%impfil.eq.2) then
         call error('xsmcl',' forbidden operation in read-only mode.',' ')
      endif
      itemp(1)=pfxsm%ioftop
      call kdiput(pfxsm%ifile,itemp,1,1,irc)
      if(irc.ne.0) go to 150
   endif
   call kdiget(pfxsm%ifile,itemp,2,1,irc)
   if(irc.ne.0) go to 140
   iii=itemp(1)
   if((pfxsm%idir.ne.iii).and.(pfxsm%impfil.le.1)) then
      itemp(1)=pfxsm%idir
      call kdiput(pfxsm%ifile,itemp,2,1,irc)
      if(irc.ne.0) go to 150
   endif
   irc=kdrcls(pfxsm%ifile,istatu)
   if(irc.ne.0) go to 160
   pfxsm%ifile=0
   return
   140 write(circ,'(i8)') irc
   call error(' xsmcl','(kdiget) iostat='//circ,' ')
   150 write(circ,'(i8)') irc
   call error('xsmcl','(kdiput) iostat='//circ,' ')
   160 write(circ,'(i8)') irc
   call error('xsmcl','(kdrcls) iostat='//circ,' ')
   end subroutine xsmcl
   !
   subroutine xsmdir(ind,ifile,cmt,nmt,link,iroot,myname,iofset,jlong,jtype, &
   & idir)
   !-----------------------------------------------------------------------
   !
   ! import or export a directory using the kdi utility.
   !
   ! input parameters:
   !  ind   : =1 for import ; =2 for export.
   !  ifile  : unit number of the xsm file.
   !  cmt    : list of character*12 names of each block (record or
   !           directory) that belong to the active directory extent.
   !  nmt    : number of blocks stored on the active directory extent.
   !  link   : offset of the next directory extent.
   !  iroot  : offset of any parent directory extent.
   !  myname : character*12 name of the active directory. myname='/' for
   !           the root level.
   !  iofset : offset list (position of the first element of each block
   !           that belong to the active directory extent).
   !  jlong  : length of each record (jlong=0 for a directory) that belong
   !           to the active directory extent.
   !  jtype  : type of each block that belong to the active directory
   !           extent.
   !  idir   : offset of the active directory extent.
   !-----------------------------------------------------------------------
   use util   ! provides error
   integer :: ind,ifile,nmt,link,iroot,idir
   character(len=*),dimension(:) :: cmt
   character(len=*) :: myname
   character header*4,circ*8,hgar*12
   integer,parameter :: iofmax=30,iwrd=3
   integer,dimension(:) :: iofset,jlong,jtype
   integer :: i,ibase,iofma2,irc
   !
   if(ind.eq.1) then
      call kdiget(ifile,itemp,idir,5,irc)
      if(irc.ne.0) go to 40
      write (header,'(a4)') itemp(1)
      if(header.ne.'$$$$') go to 30
      iofma2=itemp(2)
      nmt=itemp(3)
      if(nmt.gt.iofmax) go to 30
      link=itemp(4)
      iroot=itemp(5)
      call kdiget(ifile,itemp,idir+5,iwrd,irc)
      write (hgar(:4),'(a4)') itemp(1)
      write (hgar(5:8),'(a4)') itemp(2)
      write (hgar(9:),'(a4)') itemp(3)
      myname=hgar
      if(nmt.eq.0) return
      ibase=idir+5+iwrd
      call kdiget(ifile,iofset,ibase,nmt,irc)
      call kdiget(ifile,jlong,ibase+iofma2,nmt,irc)
      call kdiget(ifile,jtype,ibase+2*iofma2,nmt,irc)
      do i=1,nmt
         call kdiget(ifile,itemp,ibase+3*iofma2+(i-1)*3,3,irc)
         if(irc.ne.0) go to 40
         write (hgar(:4),'(a4)') itemp(1)
         write (hgar(5:8),'(a4)') itemp(2)
         write (hgar(9:),'(a4)') itemp(3)
         cmt(i)=hgar
      enddo
   else if(ind.eq.2) then
      header='$$$$'
      read (header,'(a4)') itemp(1)
      itemp(2)=iofmax
      itemp(3)=nmt
      itemp(4)=link
      itemp(5)=iroot
      call kdiput(ifile,itemp,idir,5,irc)
      if(irc.ne.0) go to 50
      hgar=myname
      read (hgar(:4),'(a4)') itemp(1)
      read (hgar(5:8),'(a4)') itemp(2)
      read (hgar(9:),'(a4)') itemp(3)
      call kdiput(ifile,itemp,idir+5,iwrd,irc)
      if(nmt.eq.0) return
      ibase=idir+5+iwrd
      call kdiput(ifile,iofset,ibase,nmt,irc)
      call kdiput(ifile,jlong,ibase+iofmax,nmt,irc)
      call kdiput(ifile,jtype,ibase+2*iofmax,nmt,irc)
      do i=1,nmt
         hgar=cmt(i)
         read (hgar(:4),'(a4)') itemp(1)
         read (hgar(5:8),'(a4)') itemp(2)
         read (hgar(9:),'(a4)') itemp(3)
         call kdiput(ifile,itemp,ibase+3*iofmax+(i-1)*3,3,irc)
         if(irc.ne.0) go to 50
      enddo
   endif
   return
   !
   !     abort on fatal errors
   !
   30 call error('xsmdir','unable to recover directory.',' ')
   40 write(circ,'(i8)') irc
   call error('xsmdir','(kdiget) error no. '//circ,' ')
   50 write(circ,'(i8)') irc
   call error('xsmdir','(kdiput) error no. '//circ,' ')
   end subroutine xsmdir
   !
   subroutine xsmrep(namt,ind,ifile,modif,cmt,nmt,link,iroot,myname, &
   & iofset,jlong,jtype,idir,iii)
   !-----------------------------------------------------------------------
   !
   ! find a block (record or directory) position in the active directory
   ! and related extents.
   !
   ! input parameters:
   !  namt   : character*12 name of the required block.
   !  ind    : =1 search namt ; =2 search and positionning in an empty
   !           slot of the active directory if namt does not exists.
   !  ifile  : unit number of the xsm file.
   !  modif  : =.true. if the active directory extent have been modified.
   !  cmt    : list of character*12 names of each block (record or
   !           directory) that belong to the active directory extent.
   !  nmt    : number of blocks stored on the active directory extent.
   !  link   : offset of the next directory extent.
   !  iroot  : offset of any parent directory extent.
   !  myname : character*12 name of the active directory. myname='/' for
   !           the root level.
   !  iofset : offset list (position of the first element of each block
   !           that belong to the active directory extent).
   !  jlong  : length of each record (jlong=0 for a directory) that belong
   !           to the active directory extent.
   !  jtype  : type of each block that belong to the active directory
   !           extent.
   !  idir   : offset of the active directory extent.
   !
   ! output parameter:
   !  iii    : return code. =0 if the block named namt does not exists;
   !           =position in the active directory extent if namt extsts.
   !-----------------------------------------------------------------------
   integer, parameter :: iofmax=30
   integer :: ind,ifile,nmt,link,iroot,idir,iii
   character(len=*),dimension(:) :: cmt
   character(len=*) :: namt,myname
   integer :: i,ipos,istart
   logical :: modif
   integer,dimension(:) :: iofset,jlong,jtype
   !
   ipos=-1
   if(nmt.lt.iofmax) ipos=idir
   if(nmt.eq.0) go to 50
   do 10 i=1,nmt
   if(namt.eq.cmt(i)) then
   ! the block already exists.
      iii=i
      return
   endif
   10 continue
   !
   ! the block namt does not exists in the active directory extent. we
   ! search in other extents that belong to the active directory.
   if(idir.ne.link) then
      ! recover a new directory extent.
      istart=link
      if(modif) then
         call xsmdir(2,ifile,cmt,nmt,link,iroot,myname,iofset,jlong,jtype,idir)
         modif=.false.
      endif
      idir=istart
      30 call xsmdir(1,ifile,cmt,nmt,link,iroot,myname,iofset,jlong,jtype,idir)
      if(nmt.lt.iofmax) ipos=idir
      do 40 i=1,nmt
      if(namt.eq.cmt(i)) then
         ! the block namt was found in the active directory extent.
         iii=i
         return
      endif
      40 continue
      if(link.eq.istart) go to 50
      idir=link
      go to 30
   endif
   50 iii=0
   if(ind.eq.1) return
   !
   ! positionning in an extent with empty slots.
   if((ipos.ge.0).and.(ipos.ne.idir)) then
      ! an extent with an empty slot was found.
      idir=ipos
      call xsmdir(1,ifile,cmt,nmt,link,iroot,myname,iofset,jlong,jtype,idir)
   endif
   return
   end subroutine xsmrep
   !
   subroutine xsmlib(pfxsm)
   !-----------------------------------------------------------------------
   !
   ! list the content of every extents that belong to the active directory.
   !
   !-----------------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   type(xsm_file),pointer :: pfxsm
   integer, parameter :: ntype=8
   integer :: inmt,ilong,ityxsm,itot,ity
   character namt*12,myname*12,first*12,ctype(ntype)*16
   save ctype
   data (ctype(ity),ity=1,ntype)/'directory','integer','real', &
   & 'character','double precision','logical','complex','undefined'/
   !
   namt=' '
   call xsmnxt(pfxsm,namt,myname)
   if(namt.eq.' ') then
      write (nsyso,90) myname
      return
   endif
   first=namt
   write(nsyso,100) myname
   itot=0
   inmt=0
   !
   10 inmt=inmt+1
   call xsmlen(pfxsm,namt,ilong,ityxsm)
   if(ityxsm.eq.0) then
      write (nsyso,110) inmt,namt,ctype(1)
   else
      write (nsyso,120) inmt,namt,ilong,ctype(ityxsm+1)
      itot=itot+ilong
   endif
   call xsmnxt(pfxsm,namt,myname)
   if(namt.eq.first) go to 20
   go to 10
   !
   20 write(nsyso,130) myname,itot
   return
   !
   90 format (/10h xsmlib: ',a12,24h' is an empty directory.)
   100 format (//38h xsmlib: content of active directory ',a12,2h':// &
   & 11h block name,10(1h-),4x,6hlength,4x,4htype/)
   110 format (1x,i4,3h  ',a12,1h',14x,a16)
   120 format (1x,i4,3h  ',a12,1h',i10,4x,a16)
   130 format (//37h total number of words on directory ',a12,3h' =,i10/)
   end subroutine xsmlib
   !
   subroutine xsmexp(pfxsm,nunit,imode,idir,impx)
   !-----------------------------------------------------------------------
   !
   ! export/import  the content of the xsm file using the contour method.
   ! export start from the active directory. esope version.
   !
   !  nunit  : unit number of the file where the export is performed.
   !  imode  : type of export/import file:
   !           =1 sequential unformatted; =2 sequential formatted (ascii).
   !  idir   : =1 to export ; =2 to import.
   !  impx   : equal to zero for no print.
   !
   !-----------------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use util   ! provides error
   type(xsm_file),pointer :: pfxsm
   integer :: nunit,imode,idir,impx
   integer, parameter :: nlevel=50,nblk=24,lendat=4
   integer :: i,j,k,ilevel,jlevel,ilong,jlong,ityxsm,iofset,itot,jmin,lennam
   character namt*12,myname*12,path(nlevel)*12,first(nlevel)*12,hsmg*131
   integer,allocatable,dimension(:) :: ibase
   real,allocatable,dimension(:) :: rbase
   character(len=4),allocatable,dimension(:) :: hbase
   double precision,allocatable,dimension(:) :: dbase
   !
   if((imode.lt.1).or.(imode.gt.2)) then
   write(hsmg,'(44hxsmexp: invalid file type on the xsm file lo, &
      & 13hcated on unit,i3,1h.)') nunit
      call error('xsmexp',hsmg,' ')
   else if((idir.lt.1).or.(idir.gt.2)) then
      write(hsmg,'(44hxsmexp: invalid action on the xsm file locat, &
      & 10hed on unit,i3,1h.)') idir
      call error('xsmexp',hsmg,' ')
   endif
   itot=0
   ilevel=1
   namt=' '
   call xsmnxt(pfxsm,namt,myname)
   if(idir.eq.2) go to 170
   !
   ! xsm file export.
   if(impx.gt.0) write(nsyso,300) 'export on',nunit,myname
   10 namt=' '
   lennam=12
   call xsmnxt(pfxsm,namt,myname)
   if(namt.eq.' ') then
      if(ilevel.eq.1) return
      namt=path(ilevel)
      ilevel=ilevel-1
      call xsmsix(pfxsm,' ',2)
      if((nunit.ne.0).and.(imode.eq.1)) then
         write(nunit) 0,0,0,0
      else if((nunit.ne.0).and.(imode.eq.2)) then
         write(nunit,310) 0,0,0,0
      endif
      if(impx.gt.0) write(nsyso,350) ilevel
      go to 130
   endif
   first(ilevel)=namt
   !
   20 call xsmlen(pfxsm,namt,ilong,ityxsm)
   jlong=ilong
   if((ityxsm.eq.4).or.(ityxsm.eq.6)) jlong=ilong/2
   if(impx.gt.0) write(nsyso,320) ilevel,namt,ityxsm,abs(jlong)
   if((nunit.ne.0).and.(imode.eq.1)) then
      write(nunit) ilevel,lennam,ityxsm,abs(jlong)
      if(lennam.gt.0) write(nunit) namt
   else if((nunit.ne.0).and.(imode.eq.2)) then
      write(nunit,310) ilevel,lennam,ityxsm,abs(jlong)
      if(lennam.gt.0) write(nunit,'(a12,68(1h ))') namt
   endif
   if(ityxsm.eq.0) then
      ! directory data.
      ilevel=ilevel+1
      if(ilevel.gt.nlevel) call error('xsmexp','too many directory levels.',' ')
      call xsmsix(pfxsm,namt,1)
      path(ilevel)=namt
      go to 10
   else if((ilong.ne.0).and.(ityxsm.le.6)) then
      itot=itot+ilong
      ! export a node.
      if(ityxsm.eq.1) then
         ! integer data.
         allocate(ibase(jlong))
         call xsmget(pfxsm,namt,ibase)
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(jlong-1)/nblk
               jmin=min(nblk,jlong-(i-1)*nblk)
               iofset=(i-1)*nblk
               write(nunit) (ibase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            write(nunit,'(8i10)') (ibase(i),i=1,jlong)
         endif
         deallocate(ibase)
      else if(ityxsm.eq.2) then
         ! single precision data.
         allocate(rbase(jlong))
         call xsmget(pfxsm,namt,rbase)
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(jlong-1)/nblk
               jmin=min(nblk,jlong-(i-1)*nblk)
               iofset=(i-1)*nblk
               write(nunit) (rbase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            write(nunit,'(1p,5e16.8)') (rbase(i),i=1,jlong)
         endif
         deallocate(rbase)
      else if(ityxsm.eq.3) then
         ! character*4 data.
         allocate(hbase(jlong))
         call xsmget(pfxsm,namt,hbase)
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(jlong-1)/nblk
               jmin=min(nblk,jlong-(i-1)*nblk)
               write(nunit) (lendat,k=1,jmin)
            enddo
            do i=1,1+(jlong-1)/nblk
               jmin=min(nblk,jlong-(i-1)*nblk)
               iofset=(i-1)*nblk
               write(nunit) (hbase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            write(nunit,'(8i10)') (lendat,i=1,jlong)
            write(nunit,'(20a4)') (hbase(i),i=1,jlong)
         endif
         deallocate(hbase)
      else if(ityxsm.eq.4) then
         ! double precision data.
         allocate(dbase(jlong))
         call xsmget(pfxsm,namt,dbase)
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(jlong-1)/nblk
               jmin=min(nblk,jlong-(i-1)*nblk)
               iofset=(i-1)*nblk
               write(nunit) (dbase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            write(nunit,'(1p,4d20.12)') (dbase(i),i=1,jlong)
         endif
         deallocate(dbase)
      endif
   endif
   130 call xsmnxt(pfxsm,namt,myname)
   if(namt.eq.first(ilevel)) then
      if((nunit.ne.0).and.(imode.eq.1)) then
         write(nunit) -ilevel,0,0,0
      else if((nunit.ne.0).and.(imode.eq.2)) then
         write(nunit,310) -ilevel,0,0,0
      endif
      if(impx.gt.0) write(nsyso,350) ilevel
      if(ilevel.eq.1) go to 140
      namt=path(ilevel)
      ilevel=ilevel-1
      call xsmsix(pfxsm,' ',2)
      go to 130
   endif
   go to 20
   140 if(impx.gt.0) write(nsyso,330) 'exported',itot
   return
   !
   ! file import.
   170 if(impx.gt.0) write(nsyso,300) 'import from',nunit,myname
   180 if((nunit.ne.0).and.(imode.eq.1)) then
      read(nunit,end=270) jlevel,lennam,ityxsm,ilong
   else if((nunit.ne.0).and.(imode.eq.2)) then
      read(nunit,340,end=270) jlevel,lennam,ityxsm,ilong
   endif
   if(jlevel.eq.-1) go to 270
   if(jlevel.eq.-ilevel) go to 260
   if(lennam.gt.12) call error('xsmexp','a record name is greater th" &
   & //"an 12 characters.',' ')
   if(lennam.eq.0) then
      namt=' '
   else if((nunit.ne.0).and.(imode.eq.1)) then
      read(nunit) namt
   else if((nunit.ne.0).and.(imode.eq.2)) then
      read(nunit,'(a12)') namt
   endif
   !
   if(impx.gt.0) write(nsyso,320) jlevel,namt,ityxsm,ilong
   if(jlevel.ne.ilevel) then
      write(hsmg,'(44hxsmexp: unable to import the xsm file locate, &
      & 9hd on unit,i3,1h.)') nunit
      call error('xsmexp',hsmg,' ')
   endif
   if(ityxsm.eq.0) then
      ! import associative table data.
      ilevel=ilevel+1
      call xsmsix(pfxsm,namt,1)
      go to 180
   else if(ityxsm.le.6) then
      ! import a node
      if(ityxsm.eq.1) then
         ! integer data.
         allocate(ibase(ilong))
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(ilong-1)/nblk
               jmin=min(nblk,ilong-(i-1)*nblk)
               iofset=(i-1)*nblk
               read(nunit) (ibase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            read(nunit,'(8i10)') (ibase(i),i=1,ilong)
         endif
         call xsmput(pfxsm,namt,ibase)
         deallocate(ibase)
      else if(ityxsm.eq.2) then
         ! single precision data.
         allocate(rbase(ilong))
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(ilong-1)/nblk
               jmin=min(nblk,ilong-(i-1)*nblk)
               iofset=(i-1)*nblk
               read(nunit) (rbase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            read(nunit,'(5e16.0)') (rbase(i),i=1,ilong)
         endif
         call xsmput(pfxsm,namt,rbase)
         deallocate(rbase)
      else if(ityxsm.eq.3) then
         ! character*4 data.
         allocate(ibase(ilong))
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(ilong-1)/nblk
               jmin=min(nblk,ilong-(i-1)*nblk)
               iofset=(i-1)*nblk
               read(nunit) (ibase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            read(nunit,'(8i10)') (ibase(i),i=1,ilong)
         endif
         ilong=(sum(ibase)+3)/4
         deallocate(ibase)
         allocate(hbase(ilong))
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(ilong-1)/nblk
               jmin=min(nblk,ilong-(i-1)*nblk)
               iofset=(i-1)*nblk
               read(nunit) (hbase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            read(nunit,'(20a4)') (hbase(i),i=1,ilong)
         endif
         call xsmput(pfxsm,namt,hbase)
         deallocate(hbase)
      else if(ityxsm.eq.4) then
         ! double precision data.
         allocate(dbase(ilong))
         if((nunit.ne.0).and.(imode.eq.1)) then
            do i=1,1+(ilong-1)/nblk
               jmin=min(nblk,ilong-(i-1)*nblk)
               iofset=(i-1)*nblk
               read(nunit) (dbase(iofset+j),j=1,jmin)
            enddo
         else if((nunit.ne.0).and.(imode.eq.2)) then
            read(nunit,'(4d20.0)') (dbase(i),i=1,ilong)
         endif
         call xsmput(pfxsm,namt,dbase)
         deallocate(dbase)
      endif
      jlong=ilong
      if((ityxsm.eq.4).or.(ityxsm.eq.6)) jlong=2*ilong
      itot=itot+jlong
   else
      write(hsmg,'(2a,i3,a)') 'try to import unknown type record f', &
      & 'rom the xsm file located on unit ''',nunit,'''(1).'
      call error('xsmexp',hsmg,' ')
   endif
   go to 180
   !
   260 if(impx.gt.0) write(nsyso,350) ilevel
   ilevel=ilevel-1
   call xsmsix(pfxsm,' ',2)
   go to 180
   !
   270 if(ilevel.ne.1) call error('xsmexp','truncated import file.',' ')
   if(impx.gt.0) write(nsyso,330) 'imported',itot
   return
   !
   300 format (//18h xsmexp: xsm file ,a,5h unit,i3, &
   & 17h from directory ',a12,2h'.//18h level  block name, &
   & 4(1h-),4x,12htype  length/)
   310 format ('->',4i8,32(1h ),' <-')
   320 format (1x,i5,3h  ',a12,1h',2i8)
   330 format (//31h xsmexp: total number of words ,a,2h =,i10/)
   340 format (2x,4i8)
   350 format (1x,i5,2x,14('-'))
   end subroutine xsmexp
end module xsmm

