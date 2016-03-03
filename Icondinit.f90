!================================================================
!================================================================
!================================================================
!================================================================
!THIS IS LIKE CONDINIT, BUT NOW IT KEEPS ONLY THE VELOCITY FIELD FROM THE INPUT GRAFIC FILE AND THE DENSITY, PRESSURE, AND POSSIBLY OTHER ICs ARE HERE GENERATED.
subroutine condinit_cylinder(x,u,dx,nn,rho_center)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

INTERFACE
	subroutine rho_profile(d, rho, rho_center, k)
        	double precision, INTENT(IN) :: d,rho_center
        	double precision, INTENT(OUT) :: rho
		integer, intent(in) :: k
       	END SUBROUTINE rho_profile
END INTERFACE

INTERFACE
	subroutine distance_new(x01,x02,x03,x1x,x1y,x1z,x2x,x2y,x2z,d,n)
		double precision, intent(in)::x1x,x1y,x1z 		! 1st point on fil given
		double precision, intent(in)::x2x,x2y,x2z 		! 2nd point on fil given
		double precision, intent(in)::x01,x02,x03 		! coords of point in space given		! 
		double precision, intent(out)::d,n		!distances
	END SUBROUTINE distance_new
END INTERFACE

  real(dp),dimension(1:MAXREGION)::rho_center, edge !Inhereted values from init_flow_fine
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  integer::i,ivar,k, check			!i goes from 1 to nn. check is either 0 or 1.
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::d1, d2                  	!distance from point to fil1 and fil2
  real(dp),dimension(3)::x11 		!one point defining the first fil
  real(dp),dimension(3)::x12 		!another point defining the first fil
  real(dp),dimension(3)::x21, x22		!two points defining the second filament
  real(dp)::rho_tot,n			
  real(dp),dimension(1:MAXREGION)::rho=0.
  integer:: p, odyn, dva, try  
  real(dp),allocatable::array(:,:) !Good name eh? First coord should be 2**(lvl + 1)**NDIM 
  integer::m		!Number of filaments	
  real(dp)::boltz, G, pi,m_h, big_K	!thermal Boltzmann, Gravity, pi, mass of hydrogen, 
                                !scratch for calc time saving
  integer::first_coord, allocatestatus, deallocatestatus, read_index
!  real(dp)::turb_pwr_multiplicative_factor    !function to multiply power of turbulence wrt to density
!  real(dp):: factor_exponent    !exponent for the turb_pwr_multiplicative_factor
  real(dp)::turb_pwr_multiplicative_factor !
!parameter (boltz = 1.380648813e-56)	!L**2M/T*t**2
parameter (G = 6.67384d0)		
parameter (pi = 3d0*atan(sqrt(3d0)))
!parameter (m_h=1.673534e-57)		!M
!factor_exponent = 1d0!this should be a namelist parameter at some point
first_coord= 2**(levelmin)
!write(*,*) "FIRST first cooord", first_coord
first_coord = first_coord**NDIM !first coordinate for ARRAY. Second should be NDIM
!write(*,*) "SECND fircoord=", first_coord, "NDIM=", NDIM, "lvlmin", levelmin
  allocate (array(first_coord,NDIM), stat=allocatestatus) !stat will tell you if
!  write(*,*)  "sizeof(array)", sizeof(array)                                                    !there's enough memory
IF (AllocateStatus /= 0) then
        write(*,*) "Not enough memory"
        STOP 
endif
  ! Set some (tiny) default values in case n_region=0
  q(1:nn,1)=d_background
  q(1:nn,2)=0.0d0
#if NDIM>1
  q(1:nn,3)=0.0d0
#endif
#if NDIM>2
  q(1:nn,4)=0.0d0
#endif
  q(1:nn,ndim+2)=smallr*smallc**2/gamma
#if NVAR > NDIM + 2
  do ivar=ndim+3,nvar
     q(1:nn,ivar)=0.0d0
  end do
#endif
k=0
check = 0	!check==0 means you are in background, check==1 means in at least one filament
IF (nn .eq. 0) RETURN

rho_tot = 0.0d0
#if NDIM<3
write(*,*) "ICondinit can only be used for NDIM==3" 
Stop
#endif
 
  open(unit=101,file='vxfile',access='stream',status='old',form='unformatted')
  read(101) array(:,1) 
  open(unit=102,file='vyfile',access='stream',form='unformatted',status='old')
  read(102) array(:,2)
  open(unit=103,file='vzfile',access='stream',form='unformatted',status='old')
  read(103) array(:,3)
close(101)
close(102)
close(103)

if(d_background .le. 1e-10) then

  iloop: DO i=1,nn !Should put something here that stops the distance_ etc and
                        !rho_profile from being called if NDIM <3
    odyn = int(x(i,1)/dx)
    dva = int(x(i,2)/dx)
    try = int(x(i,3)/dx)
    if (odyn<0 .or. odyn>(2**levelmin)) write(*,*) "odyn OOB ==", odyn
    if (dva<0 .or. dva>(2**levelmin)) write(*,*) "dva OOB ==", dva
    if (try<0 .or. try>(2**levelmin))  write(*,*)"try OOB ==", try
    p = odyn*(2**levelmin)**2 + dva*(2**levelmin)+ try + 1
    if(p<1 .or. p>first_coord) write(*,*)"p OOB ==", p
    kloop:	do k=1,nregion
      call distance_new(x(i,1),x(i,2),x(i,3),x1_center(k),y1_center(k),z1_center(k),x2_center(k),y2_center(k),z2_center(k),d1,n)
      call rho_profile(d1,rho(k), rho_center(k), k)  !calc density given d, store val in rho.
      turb_pwr_multiplicative_factor = (rho(k)/2.0d0)**(-factor_exponent) !denom is ave density
      q(i,ndim+2)=p_region(k)
#if NVAR>NDIM+2
      do ivar=ndim+3,nvar
        q(i,ivar)=var_region(k,ivar-ndim-2)
      end do
#endif
      rho_tot = rho(k) + rho_tot
    ENDDO kloop
    q(i,1) = rho_tot !Set the primitive variables from namelist values, DENSITY
    rho_tot = 0.0d0
    q(i,2) =  turb_pwr_multiplicative_factor*array(p,1) !x velocity
#if NDIM>1
    q(i,3) = turb_pwr_multiplicative_factor*array(p,2) !y velocity
#endif
#if NDIM>2
    q(i,4) =  turb_pwr_multiplicative_factor*array(p,3) !z velocity
#endif
  enddo iloop

else !there is a background

  ilooptwo: DO i=1,nn
    odyn = int(x(i,1)/dx)
    dva = int(x(i,2)/dx)
    try = int(x(i,3)/dx)
    if (odyn<0 .or. odyn>(2**levelmin)) write(*,*) "odyn OOB ==", odyn
    if (dva<0 .or. dva>(2**levelmin)) write(*,*) "dva OOB ==", dva
    if (try<0 .or. try>(2**levelmin))  write(*,*)"try OOB ==", try
    p = odyn*(2**levelmin)**2 + dva*(2**levelmin)+ try + 1
    if(p<1 .or. p>first_coord) write(*,*)"p OOB ==", p

   if(filamentary .eq. .true.) then
         klooptwo:	do k=1,nregion

      call distance_new(x(i,1),x(i,2),x(i,3),x1_center(k),y1_center(k),z1_center(k),x2_center(k),y2_center(k),z2_center(k),d1,n)
      call rho_profile(d1,rho(k), rho_center(k), k)  !calc density given d, store val in rho.
      turb_pwr_multiplicative_factor = (rho(k)/2.0d0)**(-factor_exponent) !denom is ave density
      if(rho(k) .ge. d_background) then

      rho_tot = rho(k) + rho_tot
      endif !rho(k) >= d_back if

      q(i,ndim+2)=p_region(k)
#if NVAR>NDIM+2
      do ivar=ndim+3,nvar
        q(i,ivar)=var_region(k,ivar-ndim-2)
      end do
#endif
    ENDDO klooptwo
endif !filamentary .eq. .true. if
if(rho_tot .gt. 0d0 )then
        q(i,1) = rho_tot
endif

    !Set the primitive variables from namelist values
    !    print *, "q(i,1) = ", q(i,1), "rho_tot", rho_tot, "rho(1)", rho(1)
    rho_tot = 0.0d0
    q(i,2) =  array(p,1) !x velocity
#if NDIM>1
    q(i,3) = array(p,2) !y velocity
#endif
#if NDIM>2
    q(i,4) =  array(p,3) !z velocity
#endif

  enddo ilooptwo
endif !background if




     
  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  	u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  	u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  	u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  	u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,ndim+2+ivar)=q(1:nn,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+ivar)
  enddo
#endif
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do


   DEALLOCATE (array, STAT = DeAllocateStatus)

if(deallocatestatus /= 0) write(*,*) "Deallocate of ARRAY unsuccesful"
end subroutine condinit_cylinder
