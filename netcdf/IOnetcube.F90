 subroutine create_nc_cube(ncid, name, titlerun,  n1, n2, n3, dx, dy, dz, time, dtime, va2max,  ier)
 implicit none
 include 'netcdf.inc'

 character(*) name
 character(*) titlerun
 
 integer :: ier, ncid
 integer :: n1, n2, n3

 integer :: status

 integer :: dimx, dimy, dimz
 integer :: dims(3), varid 
 integer :: numberdims, length 

 ! --- attributes
 real(kind=4) time, dtime, va2max
 real(kind=4) dx, dy, dz



 
 ier = 0 
 status = nf_create( name, OR(NF_CLOBBER,NF_64BIT_OFFSET), ncid )


!   --- include, run-title, time, time-step, resolution: 

 length = len(titlerun)
 status = nf_put_att_text( ncid, NF_GLOBAL, 'run-title', length, titlerun)
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_put_att_real( ncid, NF_GLOBAL, 'time', NF_REAL, 1, time) 
 if ( status .ne. NF_NOERR ) call  handle_error(status)


 status = nf_put_att_real( ncid, NF_GLOBAL, 'tstep', NF_REAL, 1, dtime)
 if ( status .ne. NF_NOERR ) call  handle_error(status)


 status = nf_put_att_real( ncid, NF_GLOBAL, 'dx', NF_REAL, 1, dx)
 if ( status .ne. NF_NOERR ) call  handle_error(status)


 status = nf_put_att_real( ncid, NF_GLOBAL, 'dy', NF_REAL, 1, dy)
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_put_att_real( ncid, NF_GLOBAL, 'dz', NF_REAL, 1, dz)
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_put_att_real( ncid, NF_GLOBAL, 'va2max', NF_REAL, 1, va2max)
 if ( status .ne. NF_NOERR ) call  handle_error(status)


! --------------------
! ---- start defining the dimensions of the variables: 

 status = nf_def_dim( ncid, 'nx', n1, dimx )
 if ( status .ne. NF_NOERR ) call  handle_error(status) 

 status = nf_def_dim( ncid, 'ny', n2, dimy )
 if ( status .ne. NF_NOERR ) call  handle_error(status)     

 if (n3 .eq. 1 ) then 
   numberdims = 2
 else 
   status = nf_def_dim( ncid, 'nz', n3, dimz )
   if ( status .ne. NF_NOERR ) call  handle_error(status)     
   numberdims = 3
 endif 

 dims(1) = dimx
 dims(2) = dimy
 if (numberdims .eq. 3)  dims(3) = dimz

 ! -------- create variables

 status = nf_def_var( ncid, 'R',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'P',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'T',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)
 
 status = nf_def_var( ncid, 'E',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'U',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'V',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'W',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'Bx',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'By',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)

 status = nf_def_var( ncid, 'Bz',  NF_REAL, numberdims ,dims, varid )
 if ( status .ne. NF_NOERR ) call  handle_error(status)
 
! nf_enddef puts named file out of "define" mode
  status = nf_enddef(ncid)
  if ( status .ne. NF_NOERR ) call  handle_error(status)

! nf_sync synchronises disk writes with memory buffers (immediately 
! available after writing)
  status = nf_sync(ncid)
  if ( status .ne. NF_NOERR ) call  handle_error(status)


 end subroutine 



!  subroutine close_netcdf(ncid, ier)
! implicit none

! include 'netcdf.inc'
! integer ier, ncid, status

! ier = 0

!      status = nf_sync(ncid)
!      if ( status .ne. NF_NOERR ) call  handle_error(status)

! Closes off named NetCDF file
!      status = nf_close(ncid)
!      if ( status .ne. NF_NOERR ) call handle_error(status)


!      return
 

! end subroutine




 subroutine write_nc_cube(ncid, myrank, nproc, rho, T, P, eps, u, v, w, bx, by, bz, nl,  n1, n2, n3, comm,  ier)

  implicit none

  include 'netcdf.inc'
#ifdef MPI
  include 'mpif.h'
  real(kind=4) dwk(nl,n2, n3)
  integer ncount, vartag

#endif
   integer myrank, comm, nproc
   integer nl
   integer ier, ncid, n1, n2, n3
   integer status
   integer varid
   integer icount3(3), istart3(3)

   integer i, j
   real(kind = 4) rho(nl, n2,  n3), T(nl, n2,  n3), P(nl, n2, n3), eps(nl, n2, n3)
   real(kind= 4) u(nl, n2, n3), v(nl, n2, n3), w(nl, n2, n3), bx(nl, n2, n3), by(nl, n2, n3), bz(nl, n2, n3)
 
  
   ier = 0

#ifdef MPI
!   Rho 

   vartag = 11
   call MPI_Barrier(comm, ier) 

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, rho, dwk, 'R' , nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  Density '
   else
     ncount = n2*nl*n3
     call MPI_Send(rho, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif

! -- Temp 
   vartag = 12
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, T, dwk, 'T', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  Temp '
   else
     ncount = n2*nl*n3
     call MPI_Send(T, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif

! -- Pressure  
   vartag = 13
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, P, dwk, 'P', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finished  gather write  Pres '
   else
     ncount = n2*nl*n3
     call MPI_Send(P, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif

! ---- energy 

   vartag = 14
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, eps, dwk, 'E', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  eps  '
   else
     ncount = n2*nl*n3
     call MPI_Send(eps, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif


! --- vx 
   vartag = 15
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, u, dwk, 'U', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  u-velocity '
   else
     ncount = n2*nl*n3
     call MPI_Send(u, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif

! ---- vy  
   vartag = 16
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, v, dwk, 'V', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  v-velocity '
   else
     ncount = n2*nl*n3
     call MPI_Send(v, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif

! --- vz
   vartag = 17
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, w, dwk, 'W', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  w-velocity '
   else
     ncount = n2*nl*n3
     call MPI_Send(w, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif

! ---- Bx
   vartag = 18
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, bx, dwk, 'Bx', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  Bx field '
   else
     ncount = n2*nl*n3
     call MPI_Send(bx, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif

! --- By 
   vartag = 19
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, by, dwk, 'By', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write  By fields '
   else
     ncount = n2*nl*n3
     call MPI_Send(By, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif
! --- Bz 
   vartag = 20
   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call gatherwrite(ncid, comm, nproc, bz, dwk, 'Bz', nl,  n1, n2, n3, vartag, ier )
     print*, 'Finish gather write Bz fields  '
   else
     ncount = n2*nl*n3
     call MPI_Send(bz, ncount, MPI_REAL, 0, vartag, comm, ier)
   endif




#else

!  Rho 

   status = nf_inq_varid( ncid,  'R', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)

   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, rho(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do


!  Temp  
   
   status = nf_inq_varid( ncid,  'T', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)

   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, T(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do


!  pres   

   status = nf_inq_varid( ncid,  'P', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)
   
   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, P(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do

!  EPS   

   status = nf_inq_varid( ncid,  'E', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)
   
   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, eps(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do

! velocities: 



!  vx - u velocity   

   status = nf_inq_varid( ncid,  'U', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)
   
   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, u(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do


!  vy - v velocity   

   status = nf_inq_varid( ncid,  'V', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)

   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, v(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do


!  vx - u velocity   

   status = nf_inq_varid( ncid,  'W', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)

   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, w(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do

! ---- magnetic fields: 


!  Bx -    

   status = nf_inq_varid( ncid,  'Bx', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)

   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, bx(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do



   status = nf_inq_varid( ncid,  'By', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)
   
   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, by(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do



   status = nf_inq_varid( ncid,  'Bz', varid )
   if ( status .ne. NF_NOERR) call handle_error(status)
   
   do i = 1, n3
      do j = 1, n2

         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = n1
         icount3(2) = 1
         icount3(3) = 1

          status = nf_put_vara_real(ncid,varid,istart3, icount3, bz(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

         end do
      end do




   return
#endif


 end subroutine


 subroutine open_nc_cube(ncid, name, n1, n2, n3, dx, dy, dz, time, ier ) 

!
      implicit none
      include 'netcdf.inc'
      integer          ier
      character(len=*)   name
      integer::        n1, n2, n3 
      integer::        n1id, n2id, n3id
      integer::        status, varid, ncid


      integer :: numberdims, length

 ! --- attributes
      real(kind=4) rtime, rdx, rdy, rdz
      real(kind=8) time, dx, dy, dz
      integer:: dxid, dyid, dzid, timeid
      logical:: oldcubes

      oldcubes = .false. 



      ier = -1
      ncid = -1

! NF_OPEN reads named file, NF_WRITE implies file is writeable
! ncid is returned NetCDF file ID.

      status =  NF_OPEN( trim(name), NF_WRITE, ncid )
      if ( status .ne. NF_NOERR ) then
         print*,'ERROR: Unable to open file ', name
         call  handle_error(status)
         return
      end if



!.......Read the dimension id's
! nf_inq_dimid returns dimension ID (integer) of "string" in ncid

      status = nf_inq_dimid( ncid, 'nx', n1id )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'ny', n2id )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      numberdims = 3
      status = nf_inq_dimid( ncid, 'nz', n3id )
      if ( status .ne. NF_NOERR ) call  handle_error(status)



!.......Read the actual dimensions
! nf_inq_dimlen returns length of dimension
!
      status = nf_inq_dimlen( ncid, n1id, n1 )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, n2id, n2 )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_inq_dimlen( ncid, n3id, n3 )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      ier = 0





      if ( .not. oldcubes) then 

!.......Read the actual attributes 

      status = nf_get_att_real( ncid, NF_GLOBAL, 'dx', rdx )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_get_att_real( ncid, NF_GLOBAL, 'dy', rdy )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_get_att_real( ncid, NF_GLOBAL, 'dz', rdz )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_get_att_real( ncid, NF_GLOBAL, 'time', rtime  )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
    
      dx = dble(rdx)
      dy = dble(rdy)
      dz = dble(rdz)
      time = dble(rtime) 

      ier = 0 

      else 

        time = 0.0d0 
        dx = 1757812.500000d0 
        dy = dx
        dz = 1000000.0d0 
 
      endif 

      return


 end subroutine 
!---------------------------------------------------------------------------------------------------

 subroutine read_nc_cube_one(ncid, myrank, nproc, varname, var, nl,  n1, n2, n3, comm,  ier)

  implicit none

  include 'netcdf.inc'
#ifdef MPI
  include 'mpif.h'
  integer ncount, vartag
  integer ireq
  integer mstat(MPI_STATUS_SIZE,1)
#endif
   real(kind=4) dwk(nl, n2, n3) 
   integer myrank, comm, nproc
   integer nl
   integer ier, ncid, n1, n2, n3
   integer ln1, ln2, ln3
   integer nxid, nyid, nzid
   integer status
   integer varid
   integer icount(3), istart(3)

   integer i, j, k
   real(kind = 4) var(nl, n2, n3)
   character(*) varname
!-------------------------------------------
   ier = 0
!-----------------------------------------------------------------
      if ( myrank .eq. 0 ) then
!-----------------------------------------------------------------

!
!.......Read the dimensions id's
!
        status = nf_inq_dimid( ncid, 'nx', nxid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimid( ncid, 'ny', nyid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimid( ncid, 'nz', nzid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

!.......Read the dimensions 
! 
        status = nf_inq_dimlen( ncid, nxid, ln1 )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimlen( ncid, nyid, ln2 )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimlen( ncid, nzid, ln3 )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

!
!.......Paranoid check
!

        if (( ln1 .ne. n1 ).or.( ln2.ne.n2).or.(ln3.ne.n3)) then
          write(*,*)'ERROR: Dimensions disagree ', ln1, n1,ln2,n2,ln3,n3
          ier = -1
          return
        end if

!-----------------------------------------------------------------
      end if


#ifdef MPI

   vartag = 11

   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call ReadBcast(ncid, comm, nproc, var, dwk, trim(varname) , nl,  n1, n2, n3, vartag, ier )

   else
     ncount = n2*nl*n3
     call MPI_Irecv( var, ncount, MPI_REAL, 0, vartag, comm, ireq, ier )
     call MPI_Waitall( 1, ireq, mstat, ier )
   endif 

#else 

        status = nf_inq_varid( ncid, trim(varname),  varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, var(1,j,k) )
           end do
         end do


#endif 



   return

 end subroutine read_nc_cube_one 


! --------------------------------------------------------------------------------------------------
 subroutine read_nc_cube(ncid, myrank, nproc, rho, T, P, eps, u, v, w, bx, by, bz, nl,  n1, n2, n3, comm,  ier)

  implicit none

  include 'netcdf.inc'
#ifdef MPI
  include 'mpif.h'
  integer ncount, vartag(10) 
  integer ireq(10)
  integer mstat(MPI_STATUS_SIZE,10)
#endif
   integer myrank, comm, nproc
   integer nl
   integer ier, ncid, n1, n2, n3
   integer ln1, ln2, ln3
   integer nxid, nyid, nzid
   integer status
   integer varid
   integer icount(3), istart(3)

   integer i, j, k 
   real(kind=4) dwk(nl,n2, n3)
   real(kind = 4) rho(nl, n2,  n3), T(nl, n2,  n3), P(nl, n2, n3), eps(nl, n2, n3)
   real(kind= 4) u(nl, n2, n3), v(nl, n2, n3), w(nl, n2, n3), bx(nl, n2, n3), by(nl, n2, n3), bz(nl, n2, n3)


   ier = 0


      ier = 0
!-----------------------------------------------------------------
      if ( myrank .eq. 0 ) then
!-----------------------------------------------------------------

!
!.......Read the dimensions id's
!
        status = nf_inq_dimid( ncid, 'nx', nxid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimid( ncid, 'ny', nyid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimid( ncid, 'nz', nzid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

!.......Read the dimensions 
! 
        status = nf_inq_dimlen( ncid, nxid, ln1 )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimlen( ncid, nyid, ln2 )
        if ( status .ne. NF_NOERR ) call  handle_error(status)
        status = nf_inq_dimlen( ncid, nzid, ln3 )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

!
!.......Paranoid check
!

        if (( ln1 .ne. n1 ).or.( ln2.ne.n2).or.(ln3.ne.n3)) then
          write(*,*)'ERROR: Dimensions disagree ', ln1, n1,ln2,n2,ln3,n3
          ier = -1
          return
        end if

!-----------------------------------------------------------------
      end if
!-----------------------------------------------------------------



#ifdef MPI

   do i = 1, 10 
     vartag(i) = 10+i
   end do 

   call MPI_Barrier(comm, ier)

   if (myrank .eq. 0) then
     call ReadBcast(ncid, comm, nproc, rho, dwk, 'R' , nl,  n1, n2, n3, vartag(1), ier )

     call ReadBcast(ncid, comm, nproc, T, dwk, 'T' , nl,  n1, n2, n3, vartag(2), ier )
  
     call ReadBcast(ncid, comm, nproc, P, dwk, 'P' , nl,  n1, n2, n3, vartag(3), ier )

     call ReadBcast(ncid, comm, nproc, eps , dwk, 'E' , nl,  n1, n2, n3, vartag(4), ier )

     call ReadBcast(ncid, comm, nproc, u, dwk, 'U' , nl,  n1, n2, n3, vartag(5), ier )

     call ReadBcast(ncid, comm, nproc, v, dwk, 'V' , nl,  n1, n2, n3, vartag(6), ier )

     call ReadBcast(ncid, comm, nproc, w, dwk, 'W' , nl,  n1, n2, n3, vartag(7), ier )

     call ReadBcast(ncid, comm, nproc, bx, dwk, 'Bx' , nl,  n1, n2, n3, vartag(8), ier )

     call ReadBcast(ncid, comm, nproc, by, dwk, 'By' , nl,  n1, n2, n3, vartag(9), ier )

     call ReadBcast(ncid, comm, nproc, bz, dwk, 'Bz' , nl,  n1, n2, n3, vartag(10), ier )

     print*, ' all quantites read and bcast' 

   else
     ncount = n2*nl*n3
     call MPI_Irecv( rho, ncount, MPI_REAL, 0, vartag(1), comm, ireq(1), ier )

     call MPI_Irecv ( T, ncount, MPI_REAL, 0, vartag(2), comm, ireq(2), ier)

     call MPI_Irecv ( P, ncount, MPI_REAL, 0, vartag(3), comm, ireq(3), ier)

     call MPI_Irecv ( eps, ncount, MPI_REAL, 0, vartag(4), comm, ireq(4), ier)

     call MPI_Irecv ( u, ncount, MPI_REAL, 0, vartag(5), comm, ireq(5), ier)

     call MPI_Irecv ( v, ncount, MPI_REAL, 0, vartag(6), comm, ireq(6), ier)

     call MPI_Irecv ( w, ncount, MPI_REAL, 0, vartag(7), comm, ireq(7), ier)

     call MPI_Irecv ( bx, ncount, MPI_REAL, 0, vartag(8), comm, ireq(8), ier)

     call MPI_Irecv ( by, ncount, MPI_REAL, 0, vartag(9), comm, ireq(9), ier)

     call MPI_Irecv ( bz, ncount, MPI_REAL, 0, vartag(10), comm, ireq(10), ier)

     call MPI_Waitall( 10, ireq, mstat, ier )

     
   endif






#else


         !-------------------
         ! Read data as REAL
         !-------------------


         status = nf_inq_varid( ncid, 'R', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1 
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, rho(1, j, k ))

           end do
         end do

         status = nf_inq_varid( ncid, 'T', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3 
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1 
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, t(1, j, k) )
           end do
         end do


         status = nf_inq_varid( ncid, 'P', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, p(1,j,k) )

           end do
         end do


         status = nf_inq_varid( ncid, 'E', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, eps(1,j,k) )

           end do
         end do



        status = nf_inq_varid( ncid, 'U', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3 
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1 
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, u(1,j,k) )
           end do
         end do


         status = nf_inq_varid( ncid, 'V', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3 
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1  
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, v(1,j,k) )
           end do
         end do

         status = nf_inq_varid( ncid, 'W', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3 
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, w(1,j,k) )
           end do
         end do


         status = nf_inq_varid( ncid, 'Bx', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, bx(1,j,k) )
           end do
         end do

         status = nf_inq_varid( ncid, 'By', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3 
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, by(1,j,k) )
           end do
         end do

         status = nf_inq_varid( ncid, 'Bz', varid )
         if ( status .ne. NF_NOERR ) call  handle_error(status)
         do k=1,n3 
           do j=1,n2
             istart(1) = 1
             istart(2) = j
             istart(3) = k

             icount(1) = n1
             icount(2) = 1
             icount(3) = 1

             status = nf_get_vara_real(ncid, varid, istart, icount, bz(1,j,k) )
           end do
         end do

#endif 

    return 

 end subroutine read_nc_cube 
 

!---------------------------------------------------------------------------------------------------


 
! subroutine handle_error(status)
!      implicit none

!      integer  status
!      include 'netcdf.inc'

!      if ( status .ne. NF_NOERR ) then
!        print*, NF_STRERROR(STATUS)
!      end if
!end  subroutine handle_error

#ifdef MPI
 subroutine gatherwrite(ncid, comm, nproc, var, dwk, varname, nlayer,  n1, n2, n3, vartag, ier)
  implicit none
  include 'mpif.h'
  include 'netcdf.inc'

  integer ncid, comm, nrpoc, n1, n2, nlayer, n3, ier
  integer vartag, varid, nproc
  integer status(MPI_STATUS_SIZE)

  real(kind=4) var( nlayer,n2, n3)
  real(kind=4) dwk( nlayer,n2, n3)

  character(*) varname

  integer i,j,k,l
  integer istart3(3), icount3(3), iproc, ncount
  integer istt

  ier = nf_inq_varid(ncid, trim(varname), varid)
  if ( ier .ne. NF_NOERR ) call  handle_error(ier)

! loop over all procs

  do iproc = 0, nproc-1

    istt = (iproc * n1)/nproc
    ncount = n2*nlayer*n3

   if(iproc .eq. 0) then
     do i = 1, nlayer
      do j = 1, n2
       do k = 1, n3
        dwk(i,j,k) = var(i,j,k)
       end do
      end do
     end do
   else
     call MPI_Recv( dwk, ncount, MPI_REAL, iproc, vartag, comm, status, ier )
   end if

    do i = 1, n3
      do j = 1, n2

         istart3(1) = istt+1
         istart3(2) = j
         istart3(3) = i
         icount3(1) = nlayer
         icount3(2) = 1
         icount3(3) = 1

          ier  = nf_put_vara_real(ncid,varid,istart3, icount3, dwk(1,j,i) )
          if ( ier .ne. NF_NOERR) call handle_error(ier)


         end do
      end do
! close loop over all procs
  end do

 end subroutine gatherwrite

 subroutine ReadBcast (ncid, comm, nproc, var, dwk, varname, nlayer,  n1, n2, n3, vartag, ier)

  implicit none
  include 'mpif.h'
  include 'netcdf.inc'

  integer ncid, comm, nrpoc, n1, n2, nlayer, n3, ier
  integer vartag, varid, nproc
  integer status(MPI_STATUS_SIZE)

  real(kind=4) var( nlayer,n2, n3)
  real(kind=4) dwk( nlayer,n2, n3)

  character(*) varname

  integer i,j,k,l
  integer istart3(3), icount3(3), iproc, ncount
  integer istt



     ier = nf_inq_varid( ncid, trim(varname), varid )
     if ( ier .ne. NF_NOERR ) call  handle_error(ier )

     do iproc=0,nproc-1
       istt = ((iproc)*n1)/nproc +1

       ncount = nlayer*n2*n3

       do k= 1, n3 
         do j=1,n2
           istart3(1) = istt
           istart3(2) = j
           istart3(3) = k
           icount3(1) = nlayer
           icount3(2) = 1
           icount3(3) = 1

             ier  = nf_get_vara_real(ncid, varid, istart3, icount3, dwk(1,j,k))
         end do
       end do



       if ( iproc .eq. 0 ) then
         do i= 1, nlayer
           do j=1,n2
             do k=1,n3
               var(i,j,k) = dwk(i,j,k)
             end do
           end do
         end do
       else
         call MPI_Send( dwk, ncount, MPI_REAL, iproc, vartag, comm, ier )
       end if
     end do
 
   ier = 0 
   return



 end subroutine ReadBcast



#endif

 
