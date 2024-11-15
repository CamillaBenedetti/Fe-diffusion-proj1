program diffusion
!constants
real*8,parameter::pi=3.14159265359
real*8,parameter::mSun = 1.989d33 !!Sun mass in g
real*8,parameter::mEarth = 5.972d27 !!Earth mass in g
real*8,parameter::kpc_to_cm = 3.084d21      !! 1 kpc !!
real*8,parameter::au_to_cm=1.495978707d13   !! 1 AU !!
real*8,parameter::years_to_seconds=3.156d7
real*8,parameter::mu=0.61 !H mass fraction in the gas
real*8,parameter::boltz=1.38066e-16
real*8,parameter::guniv=6.6720e-8 !cgs: cm^3/(g*s^2)
real*8,parameter::mp=1.67265e-24 !proton mass in g

!problem parameters for density and mass calculations
real*8,parameter::rho0nfw = 7.35d-26 !central density for Navarro Frenk White density distrib, g/cm^3
real*8,parameter::rs=435.7*kpc_to_cm  !scale length, kpc to cm
real*8,parameter::rho0= 4.d-26   !central gas density
!! 5000 points, with BCG, isothermal gas !!old: 2.882d-26
!real*8,parameter:: rho0=9.d-26   !! 5000 points, with BCG !!
!real*8,parameter::rho0=1.6d-25   !! 5000 points, with BCG and dT/dr !! second part, modified density?
real*8,parameter::ticm=8.9e7 !temperature of isothermal gas, first considered constant

real*8,parameter::rvir=2797.*kpc_to_cm !virial radius
real*8,parameter::r500=rvir/2. !normalization constant for radius, effective radius deviding by half mass radius
real*8,parameter::fc=1.138799
real*8,parameter::mvir=1.3e15*mSun !virial mass
real*8,parameter::mbcg=1.d12*mSun   !mass of biggest central galaxy 
real*8,parameter::ahern=12.*kpc_to_cm/(1.+2.**(0.5)) !scale lenght in Hernquist profile for BCG mass

!declaration of variables 
parameter(jmax=5000) !total # of points in the grid

real*8 :: rmin,rmax
 

!declaration of arrays of dimension jmax for first part, mass and density calculations
real*8, dimension(jmax) :: r(jmax),rr(jmax), vol(jmax), mnfw(jmax),mdark(jmax), mhern(jmax), &
rhonfw(jmax), rhost(jmax), grvnfw(jmax), grvnfw_noBCG(jmax), lnd(jmax), lnd_iso(jmax), rho(jmax), &
rho_iso(jmax), mgas(jmax), fbarr(jmax), fgasr(jmax), tem(jmax)


!declarations for second part, Fe diffusion
real*8, dimension(jmax) :: zfeobs(jmax), zfe(jmax), zfest(jmax), rhofe(jmax), rhofeobs(jmax), &
amfeiniz(jmax), amfeobs(jmax), gradzfe(jmax), amfe(jmax), rhofedot(jmax)

!problem parameters for Fe diffusion calculations
real*8,parameter::aml=7.5    !! this is the mass to light ratio 
real*8,parameter::zfesol=1.8e-3 !Solar metallicity
real*8,parameter::zfesn=0.744/1.4
real*8,parameter::snu=0.5
real*8,parameter::tnow=13.7*1.e9*years_to_seconds
real*8::lturb, kappa, vturb, slope


!step 1: build the numerical grid
!we build 2 grids, shifted by dr/2
rmin = 0*kpc_to_cm
rmax = 3000.*kpc_to_cm

do j=1, jmax
    r(j) = rmin+(j-1)*rmax/(jmax-1) !build first grid array r
    rr(j)= r(j)+0.5*rmax/(jmax-1) !build second grid array rr, half a step ahead
enddo

open(20,file = 'grid.dat')
do j=1, jmax
   write(10, *) real(r(j)/kpc_to_cm), real(rr(j)/kpc_to_cm) !check the grid
enddo
close(20)
!alternative: non uniform grid



!step 2: define density with DM NFW density profile, analytical formula
do j=1,jmax
   x=rr(j)/rs !rho centered with r at r(j+1/2)
   rhonfw(j)=rho0nfw/(x*(1.+x)**2) !DM density, Navarro Frenk White profile
   !rhost(j)=mbcg/(2.*pi)*(ahern/rr(j))/(ahern+rr(j))**3 !for second part, stellar density
enddo
!change rho0 central gas density to get Mbaryonic/Mtot = 0.16-0.17 (try to automatize)
!correct hydrostatic eq equation with Temperature gradient


!step 3: integrate density to get mass of DM, variables are centered and calculated with rr
open(20, file='mass.dat')
vol(1) = 4/3*pi*r(j)**3 !calculate volume of the shell at radius rr(j-1), array count starts at 1
mnfw(1) = 0.
do j=2, jmax
   vol(j) = 4/3*pi*(r(j)**3-r(j-1)**3) !calculate volume of the shell centered at radius rr(j-1)
   x= r(j)/rs
   mnfw(j) = mnfw(j-1)+ rhonfw(j-1)*vol(j) !numerical solution, Navarro Frenk White profile for DM
   mdark(j) = mvir*(log(1.+x)-x/(1.+x))/fc !analytican solution for DM to compare
   !make it more realistic adding the BCG mass
   mhern(j) = mbcg*r(j)**2/(r(j)+ahern)**2 !Mass of stars in the BCG, here the mass in stars of other cluster galaxies is neglected

   write(20, 1001)r(j)/kpc_to_cm, mnfw(j)/mSun, mdark(j)/mSun, mhern(j)/mSun
enddo
1001 format(4(1pe12.4))
close(20)

!step 4: calculate grav potential to calculate gas density
open(20,file='grv.dat')
grvnfw(1)=0.          !! ok per alone NFW, isotermo o beta-model
grvnfw_noBCG(1)=0.
do j=2,jmax
   grvnfw_noBCG(j)=guniv*(mnfw(j))/r(j)**2 !first only with dark matter potential 
   grvnfw(j)=guniv*(mnfw(j)+mhern(j))/r(j)**2 !then also with bcg, stars neglected
   write(20,1002)r(j)/kpc_to_cm,grvnfw(j)/mSun, grvnfw_noBCG(j)/mSun
enddo
1002 format(2(1pe12.4))
close(20)

! Temperature profile
open(20,file='temperature.dat',status='unknown')
do j=1,jmax
  x=rr(j)/r500
  xx=x/0.045
  tem(j)=ticm*1.35*(xx**1.9+0.45)/(xx**1.9+1.)* &   !! this is for Perseus !!
      1./(1.+(x/0.6)**2)**0.45
  write(20,1003)rr(j)/kpc_to_cm, tem(j), ticm
enddo
close(20)
1003 format(4(1pe12.4))


! calculate the gas density, assuming ticm (temperature intercluster medium)
lnd(1)=log(rho0)          !! equation for gas in hydro eq. with the grav potential, lnd = ln of density
lnd_iso(1)=lnd(1)
rho(1)= rho0
rho_iso(1)=rho(1)
do j=2,jmax
   gg=grvnfw(j) !grav acceleration from DM halo
   gg_noBCG=grvnfw_noBCG(j)
   temmed=0.5*(tem(j)+tem(j-1))
   lnd_iso(j)=lnd_iso(j-1)-gg_noBCG*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm) !ln of density !! isoth and without BCG !!   
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*temmed) &
          - (log(tem(j)) - log(tem(j-1))) ! with BCG and Temperature gradient
   rho(j)=exp(lnd(j)) !density
   rho_iso(j)=exp(lnd_iso(j))
enddo

!write the different densities in a file
open(20,file='density.dat',status='unknown')
do j=1,jmax
   write(20,1000)rr(j)/kpc_to_cm,rho(j),rhonfw(j),rho_iso(j)
enddo
close(20)
1000 format(4(1pe12.4))

!calculate the mass of the gas integrating the density
open(20,file='mgas.dat',status='unknown')
mgas(1)=rho(1)*vol(1)
do j=2,jmax
   mgas(j)=mgas(j-1)+rho(j-1)*vol(j)
   write(20,1100)r(j)/kpc_to_cm, mgas(j)/mSun, mnfw(j)/mSun
enddo
close(20)

!step 5: calculate the barion fraction fb(r) in jmax-1 (not on the border) to check value for rho0
fbar=(mhern(jmax-1)+mgas(jmax-1))/(mnfw(jmax-1)+mgas(jmax-1)+mhern(jmax-1)) !baryonic mass = stars + gas
fgas=(mgas(jmax-1))/(mnfw(jmax-1)+mgas(jmax-1)+mhern(jmax-1))
print*,'fbar= ',real(fbar), 'fgas= ', real(fgas)

open(20,file='barfrac.dat')
do j=2,jmax-1
   fbarr(j)=(mhern(j)+mgas(j))/(mnfw(j)+mgas(j)+mhern(j))
   fgasr(j)=(mgas(j))/(mnfw(j)+mgas(j)+mhern(j))
   write(20,1100)r(j)/ckpc_to_cm, fbarr(j), fgasr(j)
enddo
close(20)
1100 format(4(1pe12.4))

!addendum: thermal pressure ok, add turbulent pressure



!***********************************************************************
!! Second part: we have the gas density profile and we can proceed
!! with the integration of the diffusion equation for rhofe
!***********************************************************************

!! Set the initial abundance profile
zfeout=0.4*zfesol   !! this is the background abundance !!

do j=1,jmax
   x=rr(j)/(80.*kpc_to_cm)
   zfeobs(j)=zfesol*0.3*1.4*1.15*(2.2+x**3)/(1+x**3)/1.15  !Perseus! First approach: observed Fe profile from XMM Newton Data, reduced to fit Chandra obs slide 39
   zfeobs(j)=zfeobs(j) - zfeout   !! subtract background abundance to get the excess 
   zfeobs(j)=max(zfeobs(j),0.)
   zfest(j)=1.*zfesol    !! set the stellar abundance !! for source 
   zfe(j)=zfeobs(j)   !!0. !!zfeout !!zfeobs(j)  !! which initial zfe? !! used for adding diffusion and source
   rhofe(j)=rho(j)*zfe(j)/1.4
   rhofeobs(j)=rho(j)*zfeobs(j)/1.4
enddo

!! Calculate the initial excess of iron mass

amfeiniz(1)=rhofe(1)*vol(1)
amfeobs(1)=rhofeobs(1)*vol(1)
do j=2,jmax
   amfeiniz(j)=amfeiniz(j-1)+rhofe(j-1)*vol(j)
   amfeobs(j)=amfeobs(j-1)+rhofeobs(j-1)*vol(j)
enddo

open(10,file='zfe_initial.dat')
open(20,file='initial.dat',status='unknown')
do j=1,jmax
   write(10,1500)rr(j)/kpc_to_cm,zfe(j)/zfesol,zfeobs(j)/zfesol, &
                 r(j)/kpc_to_cm,amfeiniz(j)/mSun,amfeobs(j)/mSun
   write(20,3001)rr(j)/kpc_to_cm+0.001,zfe(j)/zfesol,zfeobs(j)/zfesol
enddo
close(10)
close(20)
1500 format(6(1pe12.4))
3001  format(4(1pe12.4))

!have to overplot quantities at pag 42 and compare

!first check diffusin, then source alone, then combine
! for now don't add the source term

!! boundary conditions (outflows) at the grid border, before time integration
zfe(1)=zfe(2)
zfe(jmax)=zfe(jmax-1)
rhofe(1)=rho(1)*zfe(1)/1.4 !but it has already been calculated, what changes??
rhofe(jmax)=rho(jmax)*zfe(jmax)/1.4


!!  set the diffusion coefficient kappa = C*v*l (for now constant)
!open(20,file='kappa.dat',status='unknown')
vturb=260.e5   !! like Perseus !! turbulent velocity
lturb=15.*kpc_to_cm  !! this is quite uncertain !! turbulence lengthscale
!rscala=30.*kpc_to_cm  !needed for variable kappa
kappa=0.11*vturb*lturb

!! do j=1,jmax    !! for making kappa non constant
!!!    kappa(j)=0.11*vturb*lturb   !! constant !!
!!    kappa(j)=rhost(j)  !!0.333*vturb*lturb   !! constant !!
!!    kappa(j)=kappa(j)-0.6*kappa(j)*exp(-(r(j)/rscala)**2)
!!    write(20,*)real(r(j)/kpc_to_cm),kappa(j)
!! enddo
!close(20)


!! Here start the time integration (use FTCS method)
print*,'start time integration cycle'
tend=tnow
time0=tnow-5.*1.e9*years_to_seconds
time=time0
!! calculate the timestep (to be modified if the grid is non-uniform)
dt=0.4*(r(5)-r(4))**2/(2.*kappa)  !! ok for Delta_r costant !!
!print*,real(dt/years_to_seconds) !!check time step
!print*,real(time/years_to_seconds)


!! write the source terms (SNIa and stellar winds)
slope=1.1
alphast=4.7e-20*(time/tnow)**(-1.26)
alphasn=4.436e-20*(snu/aml)*(time/tnow)**(-slope)

!!  print*,'alphast,sn = ',alphast,alphasn

do j=2,jmax-2
   rhofedot(j)=(alphast*zfest(j)/1.4+alphasn*zfesn)*rhost(j)
enddo

do while (time>0 .and. time<=tend)    !!start the main time cycle
   gradzfe(1) = 0. !!boundary conditions, reset at each iteration?
   gradzfe(jmax) = 0.
   !! the equation to be solved is d(n*zfe)/dt = div(kappa*n*grad(zfe)) + S
   !! (according to Rebusco et al. 2006)
   !! Use the FTCS scheme.

!goto 776 !to skip source and check diffusion
!! source step
do j=2,jmax-1
   !!!    write(70,*)rhofe(j),dt*rhofedot(j)
   !!!    if(j.eq.5)print*,'azz ',dt,rhofe(j),dt*rhofedot(j),rhofedot(j)
       rhofe(j)=rhofe(j) + dt*rhofedot(j)
       zfe(j)=rhofe(j)/rho(j) * 1.4
    enddo
   
   !! set the boundary conditions (outflows)
   
         zfe(1)=zfe(2)
         zfe(jmax)=zfe(jmax-1)
         rhofe(1)=rhofe(2)
         rhofe(jmax)=rhofe(jmax-1)

!776   continue

  
!  diffusive step   !  check the Fe conservation !
   do j=2,jmax-1 !don't touch boundary
      zfe(j)=1.4*rhofe(j)/rho(j) !source part

      !goto777 !to check source and skip diffusion
      gradzfe(j)=(zfe(j)-zfe(j-1))/(rr(j)-rr(j-1))  !! dZ/dr centered at "j" !!
      rhojp1=0.5*(rho(j+1)+rho(j))  !! rho centered at "j+1" !!
      rhoj=0.5*(rho(j-1)+rho(j))    !! rho centered at "j" !!
      rhofe(j)=rhofe(j) &
              + (dt/1.4)*(r(j+1)**2*kappa*rhojp1*gradzfe(j+1) &
              -r(j)**2*kappa*rhoj*gradzfe(j)) / (0.33333333*(r(j+1)**3-r(j)**3))
      zfe(j)=1.4*rhofe(j)/rho(j)  !! update Z_Fe with the new rho_Fe !!
    enddo
    
    zfe(1)=zfe(2) !reset boundary conditions for the diffusion test (outflows)
    zfe(jmax)=zfe(jmax-1)
    rhofe(1)=rhofe(2)
    rhofe(jmax)=rhofe(jmax-1)
    !777   continue

    if(time>=(time0+1e9*years_to_seconds - 5.1780151e12) .and. time<=(time0+ 1e9*years_to_seconds + 5.1780151e12)) then
      open(20,file='Fe_1Gyr.dat',status='unknown')
      !need to calculate total mass of Fe and check conservation
      amfe(1)=rhofe(1)*vol(1)
      do j=2,jmax
         amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)
      enddo
      do j=1,jmax
         write(20,1700)rr(j)/kpc_to_cm,zfe(j)/zfesol, rhofe(j), amfe(j)/mSun
      enddo
      close(20)
      1700 format(4(1pe12.4))
      !print*,'masses at 2 Gyr'
      !write(6,3002)amfe(jmax)/mSun,amfeiniz(jmax)/mSun,amfeobs(jmax)/mSun
      !write(6,3003)amfe(180)/mSun,amfeiniz(180)/mSun,amfeobs(180)/mSun


    else if(time>=(time0+2e9*years_to_seconds - 5.1780151e12) .and. time<=(time0+ 2e9*years_to_seconds + 5.1780151e12)) then
         open(20,file='Fe_2Gyr.dat',status='unknown')
      !need to calculate total mass of Fe and check conservation
         amfe(1)=rhofe(1)*vol(1)
         do j=2,jmax
            amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)
         enddo
         do j=1,jmax
            write(20,1700)rr(j)/kpc_to_cm,zfe(j)/zfesol, rhofe(j), amfe(j)/mSun
         enddo
         close(20)
         !print*,'masses at 2 Gyr'
         !write(6,3002)amfe(jmax)/mSun,amfeiniz(jmax)/mSun,amfeobs(jmax)/mSun
         !write(6,3003)amfe(180)/mSun,amfeiniz(180)/mSun,amfeobs(180)/mSun

   else if(time>=(time0+5e9*years_to_seconds - 5.1780151e12) .and. time<=(time0+ 5e9*years_to_seconds + 5.1780151e12)) then
      open(20,file='Fe_5Gyr.dat',status='unknown')
      !need to calculate total mass of Fe and check conservation
      amfe(1)=rhofe(1)*vol(1)
      do j=2,jmax
         amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)
      enddo
      do j=1,jmax
         write(20,1700)rr(j)/kpc_to_cm,zfe(j)/zfesol, rhofe(j), amfe(j)/mSun
      enddo
      close(20)
      !print*,'masses at 5 Gyr'
      !write(6,3002)amfe(jmax)/mSun,amfeiniz(jmax)/mSun,amfeobs(jmax)/mSun
      !write(6,3003)amfe(180)/mSun,amfeiniz(180)/mSun,amfeobs(180)/mSun
    end if

    time=time+dt
    !!print*,'dt (yr),time (Gyr) = ', real(dt/years_to_seconds), &
            !!real(time/1.e9/years_to_seconds)
end do
print*,'end of cycle'

amfe(1)=rhofe(1)*vol(1)
do j=2,jmax
   amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)
enddo

!need to calculate total mass of Fe and check conservation
write(6,3002)amfe(jmax)/mSun,amfeiniz(jmax)/mSun,amfeobs(jmax)/mSun
!mass of Fe at 100 kpc
write(6,3003)amfe(180)/mSun,amfeiniz(180)/mSun,amfeobs(180)/mSun
3002  format('M_Fe(tot), M_Fein(tot) (Msol) = ',3(1pe12.4))
3003  format('M_Fe(100kpc), M_Fein(100kpc) (Msol) = ',3(1pe12.4))


stop
end 
