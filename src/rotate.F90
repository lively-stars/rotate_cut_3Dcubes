  module rotate

  use arrays


  CONTAINS

  subroutine rotate_cube(mu, pivot, nx, ny, nz, dx, dy, dz)

   implicit none 
   integer, intent(in) :: nx, ny, nz

   integer  Nzcut 
   integer  i,j,jj,k, m,  xn, xnp, x_l, z_l 
   real(kind=8), intent(in) :: mu, pivot, dx, dy, dz 
   real(kind=8)  theta, pivotdx,  newdx, newx, newz 
   real(kind=8)   z_f,  x_f 
   real(kind=8)   zold, znew, Tnewd 
   real(kind=8)   time1, time2, timeall 
   zold = 0.0d0
   znew = 0.0d0 
 
 
   theta = acos(mu)
   pivotdx = tan(theta) * pivot
   newdx = dz * sin(theta)
   Nzcut = Nz
   if (mu .le. 0.9) Nzcut = int(Nz*(0.9d0/mu))
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   zgrid(1)=0.0d0
   do k=2,Nzcut
      zgrid(k)=zgrid(k-1)+dz
   end do

! if mu not equal 1 we need new z indices! 

 

!---- start to do the whole cube 
!--- over the x-y plan

      call cpu_time(time1) 

      do i = 1, Nx
        do k = 1, Ny
           do j = 1, Nzcut 
             newx = (k-1)*dx - pivotdx +(j-1)*newdx 
             newz = (j-1)*(dz)*mu;

             if (newx .lt. 0.0) then 
                  newx = newx - floor(newx/(dble(Ny)*dx)) * (Ny*dx)  
                  do while( newx .lt. 0.0)  
                    newx = newx+(Ny*dx)
                  enddo 
                 
             elseif (newx .ge. (Ny*dx)) then 
                  newx = newx - (floor(newx/(dble(Ny)*dx))*(Ny*dx))  
                  do while(newx .ge. (Ny*dx))  
                    newx = newx - (floor(newx/(dble(Ny)*dx))*(Ny*dx))   
                  enddo 
             endif 

             x_f = newx/dx
             x_l = int(x_f)
             x_f = x_f-dble(x_l)
             z_f = newz/dz
             z_l = int(z_f)
             z_f = z_f-dble(z_l)  
             z_l = z_l +1
             xn = x_l + 1
             xnp = xn+1 
             if (xnp .gt. Ny) xnp = mod(xnp, Ny)  

             newT(i,k,j)  = ((T(i, xn, z_l)) * (1.0d0-x_f) + x_f*(T(i, xnp,  z_l))) *(1.0d0- z_f)
             newT(i,k,j)  = newT(i,k,j) +  z_f*((T(i, xn, z_l+1))*(1.0d0-x_f)+ x_f*(T(i, xnp,z_l+1)))

             newP(i,k,j)  = ((P(i, xn, z_l)) * (1.0d0-x_f) + x_f*(P(i, xnp,  z_l))) *(1.0d0- z_f)
             newP(i,k,j)  = newP(i,k,j) +  z_f*((P(i, xn, z_l+1))*(1.0d0-x_f)+ x_f*(P(i, xnp,z_l+1)))


             newrho(i,k,j)  = ((rho(i, xn, z_l)) * (1.0d0-x_f) + x_f*(rho(i, xnp,  z_l))) *(1.0d0- z_f)
             newrho(i,k,j)  = newrho(i,k,j) +  z_f*((rho(i, xn, z_l+1))*(1.0d0-x_f)+ x_f*(rho(i, xnp,z_l+1)))

#ifdef VELO
             newVtot(i,k,j)  = ((Vtot(i, xn, z_l)) * (1.0d0-x_f) + x_f*(Vtot(i, xnp,  z_l))) *(1.0d0- z_f)
             newVtot(i,k,j)  = newVtot(i,k,j) +  z_f*((Vtot(i, xn, z_l+1))*(1.0d0-x_f)+ x_f*(Vtot(i, xnp,z_l+1)))
#endif 

            znew = zold + dz*(j-1)*1.0d-5 

           end do
         end do 
       end do 

       call cpu_time(time2) 
       timeall = time2-time1
       print*, ' time spend in loop =', timeall 

  end subroutine 



  subroutine rotate_cube_azimuth(mu, phi, pivot, nx, ny, nz, dx, dy, dz)

   implicit none
   integer, intent(in) :: nx, ny, nz

   integer  Nzcut
   integer  i,j,jj,k, m,  yn, ynp, y_l, z_l
   integer  xn, xnp, x_l

   real(kind=8), intent(in) :: mu, pivot, dx, dy, dz
   real(kind=8)  theta, pivotdx,  newdy, newy, newz
   real(kind=8) phi , pi_n
   real(kind=8)  newdx, newx

   real(kind=8)   z_f,  y_f, x_f
   real(kind=8)   zold, znew, Tnewd
   real(kind=8)   time1, time2, timeall
   zold = 0.0d0
   znew = 0.0d0

   pi_n = 4.0*atan(1.d0)
   phi = pi_n/180.0d0*phi

   theta = acos(mu)
   pivotdx = tan(theta) * pivot
   newdy = dz * sin(theta) *cos(phi)
   newdx = dz *sin(theta)*sin(phi)

   Nzcut = Nz
   if (mu .le. 0.9) Nzcut = int(Nz*(0.9d0/mu))

!!!!!!!!!!!!!!!!!!!!!!!!!!!
   zgrid(1)=0.0d0
   do k=2,Nzcut
      zgrid(k)=zgrid(k-1)+dz
   end do

! if mu not equal 1 we need new z indices!



!---- start to do the whole cube
!--- over the x-y plan

      call cpu_time(time1)

      do i = 1, Nx
        do k = 1, Ny
           do j = 1, Nzcut
             newy = (k-1)*dy - pivotdx +(j-1)*newdy
             newx = (i-1)*dx - pivotdx + (j-1)*newdx

             newz = (j-1)*(dz)*mu;

             if (newy .lt. 0.0) then
                  newy = newy - floor(newy/(dble(Ny)*dy)) * (Ny*dy)
                  do while( newy .lt. 0.0)
                    newy = newy+(Ny*dy)
                  enddo

             elseif (newy .ge. (Ny*dy)) then
                  newy = newy - (floor(newy/(dble(Ny)*dy))*(Ny*dy))
                  do while(newy .ge. (Ny*dy))
                    newy = newy - (floor(newy/(dble(Ny)*dy))*(Ny*dy))
                  enddo
             endif

             if (newx .lt. 0.0) then
                  newx = newx - floor(newx/(dble(Nx)*dx)) * (Nx*dx)
                  do while( newx .lt. 0.0)
                    newx = newx+(Nx*dx)
                  enddo

             elseif (newx .ge. (Ny*dx)) then
                  newx = newx - (floor(newx/(dble(Nx)*dx))*(Nx*dx))
                  do while(newx .ge. (Nx*dx))
                    newx = newx - (floor(newx/(dble(Nx)*dx))*(Nx*dx))
                  enddo
             endif
! ----- iterpolate between 8 points!
! ----
             y_f = newy/dy
             y_l = int(y_f)
             y_f = y_f-dble(y_l)

             z_f = newz/dz
             z_l = int(z_f)
             z_f = z_f-dble(z_l)
             z_l = z_l +1

             yn = y_l + 1
             ynp = yn+1

             x_f = newx/dx
             x_l = int(x_f)
             x_f = x_f-dble(x_l)

             xn = x_l +1
             xnp = xn+1


             if (ynp .gt. Ny) ynp = mod(ynp, Ny)
             if (xnp .gt. Nx) xnp = mod(xnp, Nx)

! temperature
             newT(i,k,j)  = ((T(xn, yn, z_l)) * (1.0d0-y_f) + y_f*(T(xn, ynp,  z_l))) *(1.0d0- z_f)
             newT(i,k,j)  = newT(i,k,j) +  z_f*((T(xn, yn, z_l+1))*(1.0d0-y_f)+ y_f*(T(nx, ynp,z_l+1)))
             newT(i,k,j)  = newT(i,k,j) * (x_f) +(1.0d0-x_f)* ((T(xnp, yn, z_l)) * (1.0d0-y_f) + y_f*(T(xnp, ynp,  z_l))) *(1.0d0- z_f) &
     &                     + (1.0d0-x_f) *  z_f*((T(xnp, yn, z_l+1))*(1.0d0-y_f)+ y_f*(T(xnp, ynp,z_l+1)))


! pressure
            newP(i,k,j)  = ((P(xn, yn, z_l)) * (1.0d0-y_f) + y_f*(P(xn, ynp,  z_l))) *(1.0d0- z_f)
             newP(i,k,j)  = newP(i,k,j) +  z_f*((P(xn, yn, z_l+1))*(1.0d0-y_f)+ y_f*(P(nx, ynp,z_l+1)))
             newP(i,k,j)  = newP(i,k,j) * (x_f) +(1.0d0-x_f)* ((P(xnp, yn, z_l)) * (1.0d0-y_f) + y_f*(P(xnp, ynp,  z_l))) *(1.0d0- z_f) &
     &                     + (1.0d0-x_f) *  z_f*((P(xnp, yn, z_l+1))*(1.0d0-y_f)+ y_f*(P(xnp, ynp,z_l+1)))


! density

            newrho(i,k,j)  = ((rho(xn, yn, z_l)) * (1.0d0-y_f) + y_f*(rho(xn, ynp,  z_l))) *(1.0d0- z_f)
             newrho(i,k,j)  = newrho(i,k,j) +  z_f*((rho(xn, yn, z_l+1))*(1.0d0-y_f)+ y_f*(rho(nx, ynp,z_l+1)))
             newrho(i,k,j)  = newrho(i,k,j) * (x_f) +(1.0d0-x_f)* ((rho(xnp, yn, z_l)) * (1.0d0-y_f) + y_f*(rho(xnp, ynp,  z_l))) *(1.0d0- z_f) &
     &                     + (1.0d0-x_f) *  z_f*((rho(xnp, yn, z_l+1))*(1.0d0-y_f)+ y_f*(rho(xnp, ynp,z_l+1)))



#ifdef VELO

            newVtot(i,k,j)  = ((Vtot(xn, yn, z_l)) * (1.0d0-y_f) + y_f*(Vtot(xn, ynp,  z_l))) *(1.0d0- z_f)
             newVtot(i,k,j)  = newVtot(i,k,j) +  z_f*((Vtot(xn, yn, z_l+1))*(1.0d0-y_f)+ y_f*(Vtot(nx, ynp,z_l+1)))
             newVtot(i,k,j)  = newVtot(i,k,j) * (x_f) +(1.0d0-x_f)* ((Vtot(xnp, yn, z_l)) * (1.0d0-y_f) + y_f*(Vtot(xnp, ynp,  z_l))) *(1.0d0- z_f) &
     &                     + (1.0d0-x_f) *  z_f*((Vtot(xnp, yn, z_l+1))*(1.0d0-y_f)+ y_f*(Vtot(xnp, ynp,z_l+1)))


#endif

            znew = zold + dz*(j-1)*1.0d-5

           end do
         end do
       end do

       call cpu_time(time2)
       timeall = time2-time1
       print*, ' time spend in loop =', timeall

  end subroutine  rotate_cube_azimuth

  end module 
