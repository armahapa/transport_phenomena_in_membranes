!-----------------------------------------------------------------------------------------------------------------
! Numerical simulation of a square lipid bilayer with curvature induced-proteins
!
! Simulation code by: Arijit Mahapatra, David Saintillan, Padmini Ranganami (UCSD)
!
! For details of the model, refer to:
! "Transport phenomena in fluid films with curvature elasticity" by Mahapatra et al.
! (https://arxiv.org/abs/2001.07539)
!-----------------------------------------------------------------------------------------------------------------
program sumingup
implicit none
integer:: N,M,P,i,j,k,ptc,itn,lamit,pvti,pvtj,velit,ii,jj,ilg,jlg,lgp
double precision:: sig0, lam0, bsig_s, ell, kp, beta,c,err,Tst,Ted,Tstp,stth,flg,dt1
double precision:: pr, L, D1, D2, D3, D4, D5, D6, D7, D8, D9,start,finish,dlr,fc,kBT,D
double precision:: D10,lr,tT,lX,lY,dx,dy,dt,Dmf,xa,xb,ya,yb,dialn,lgsigdx,lgsigdy
double precision:: cth,lct,tol,sigtol,ur,zr,zrr,sigrr,err_zeta,err_lam,err_z
double precision:: zer,zerr,lamrr,lamr,rZe,zenu,znu,vrr,diln2,test
double precision:: signu,dzetae,dzetaw,dzetan,dzetas,zcon,zecon,velerr,utol
double precision:: zetadx,zetady,azE,azW,azN,azS,azP,uer,ver,xd,yd
double precision:: sigdx,sigdy,Lsig,Lsigsq,sigo,uW,uE,uN,uS,urr,uP,sigsqdy,D11
double precision:: zeo,zo,azeE,azeW,azeN,azeS,sigsqdx,dvy,dvx,duy,dux,d2zy,d2zxy,d2zx
double precision:: dsign,dsige,dsigw,dsigs,umean,vmean,usum,vsum,xt,yt,r,PI,ddx,ddy,wbc 
double precision:: laW,laS,laP,laN,lamo,laE,lamnu,zera,uup,vup,nu,uo,vo,delta,divm,xc,yc
double precision:: asiP,zx,zy,wold,wdx,wdy,hxx,hxy,hyy,wzedx,wzedy,Lzeta,Lzew,wxx,wyy,wxy
double precision, dimension(0:500):: x,y,siggs,sigge,siggw,siggn
double precision, dimension(0:5000):: t,tsig,Je,Jw,Jn,Js
double precision, dimension(0:20000000)::xLg,yLg
double precision, dimension(0:500,0:500):: sigr,azeP,Sze,d2h,Sz,Ssi,Sla,shst,uold,vold
double precision, dimension(0:500,0:500):: zeta,z,sig,lam,sigp,u,v,Sv,Su,zold,zetaold
double precision, dimension(0:500,0:500)::asiE,asiW,asiN,asiS,g1,g2,fx,fy,div
double precision, dimension(0:500,0:500):: oldu,oldv,oldlam,oldsig,w,zp,g0 

!! Parameters
sig0=1d-3
lam0=1d-4
bsig_s=10.d0
ell=5d-1
kp=84.d0
kBT=4.14d-2
D=1e5
nu=1d-5
pr=0.d0
L=1000.d0


!!======================= value of pi
PI=2.D0*DASIN(1.D0)

!! basic non dimensional number
D6=L**2.d0*kBT*sig0/kp        !!! C
D10=kp*ell/(kBT*L)       !!! B
D8=2.d0*L**2.d0*lam0/kp            !!! T
D11=(lam0*L**2.d0)/(nu*D)     !! Pe
D9=2.d0*L**3.d0*pr/kp              !! P


!!! Derived nondimensional number 
D1= kp*ell*sig0/(L*lam0)     !!! 2 B C /T
D2=2.d0*(ell)**2.d0*kp*sig0**2.d0/lam0   !!! 4 C^2 B^2/ T
D3=kBT*sig0/lam0         !!   2 C/ T
D5=2.d0*L*ell*sig0      !! 2 B C
D7=2.d0*kp*ell**2.d0*sig0/kBT     !!! 2 B^2 C
!!! When you specify the dimensionless numbers independently
! D1= 2*D10*D6/D8
! D2= (4.d0*D6**2.d0*D10**2.d0)/D8
! D3= 2*D6/D8
! D5=2*D6*D10
! D7=2*D10**2.d0*D6

!! Grid size (N along the x, M along the y, P along time)
N=101
M=101
P=1001
lr=1.d0             !! aspect ratio
tT=1d-1             !! totat non dimensional time

!!=========================================================================
lX=1.d0
lY=lX/lr
dx=lX/(N-1)
dy=lY/(M-1)
delta=(dx**2.d0+dy**2.d0)**0.5d0
!delta=dx
y(1)=-lY/2.d0
do j=2,M
y(j)=y(j-1)+dy
enddo
x(1)=-lX/2.d0
do i=2,N
x(i)=x(i-1)+dx
enddo
dt=tT/(P-1)
t(1)=0.d0
do k=2,P
t(k)=t(k-1)+dt
enddo

!! initiations of the variable
do i=0,N+1
do j=0,M+1
zeta(i,j)=0.d0
sig(i,j)=0.d0
z(i,j)=0.d0
lam(i,j)=1.d0
d2h(i,j)=0.d0
shst(i,j)=0
u(i,j)=0.d0
v(i,j)=0.d0
fx(i,j)=0.d0
fy(i,j)=0.d0
w(i,j)=0.d0
enddo
enddo


!! Boundary conditions
do j=1,M
z(1,j)=0.d0
z(N,j)=0.d0
lam(N,j)=1.d0
lam(1,j)=1.d0
enddo
do i=1,N
z(i,1)=0.d0
z(i,M)=0.d0
lam(i,M)=1.d0
lam(i,1)=1.d0
enddo

Tst=tT/10.d0;
Ted=tT/5.d0;
Tstp=tT/8.d0;
stth=10.d0;


!! prescribing end flux
do k=2,M
    !Je(k)=-0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))
    !Jw(k)=0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))  
    !Jn(k)=-0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))
    !Js(k)=0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))  
    Je(k)=0.d0
    Jw(k)=0.d0 
    Jn(k)=0.d0
    Js(k)=0.d0  
enddo


!! calculation of log (for the Greens function)
ilg=1

          do i=1,N   ! open 53
          if (i==1 .OR. i==N) then !open 53a
          do j=1,M    ! open 54
        
          do ii=2,N-1    !open 55
          do jj=2,M-1    ! open 56
          !if (ii/=i .OR. jj/=j)then
          xt=(x(i)-x(ii))
          yt=(y(j)-y(jj))
          r=(xt**2.d0+yt**2.d0)**(0.5d0)
          xlg(ilg)=DLOG(r)
          ilg=ilg+1
          !endif
          enddo     ! close 56
          enddo    ! close 55
          enddo    ! close 54
          endif ! close 53a
          enddo    ! close 53
 jlg=1
 
        do j=1,M   ! open 53
          if (j==1 .OR. j==M) then !open 53a
          do i=2,N-1     ! open 54
          do ii=2,N-1    !open 55
          do jj=2,M-1    ! open 56
          !if (ii/=i .OR. jj/=j)then
          xt=(x(i)-x(ii))
          yt=(y(j)-y(jj))
          r=(xt**2.d0+yt**2.d0)**(0.5d0)
          ylg(jlg)=DLOG(r)
          jlg=jlg+1
          !endif
          enddo     ! close 56
          enddo    ! close 55
          enddo    ! close 54
          endif ! close 53a
          enddo    ! close 53       



!!! Initial condition for sigma;
!! for 1 patch use lct=1.78d-1 activate xa=0, ya=0 and first line of sig exp
!!!for 2 patches lct=1.26d-1 activate xa=ya=-0.25, xb=yb=0.25 and first 2 lines of sig exp
!!! for 4 patches lct=0.89d-1 activate xa,xb,xc,xd and ya,yb,yc,yd and first 4 lines of 
!!!sig exp
cth=20.d0
lct=1.26d-1

!xa=0.d0
!ya=0.d0

xa=-0.25d0
ya=-0.25d0

xb=0.25d0
yb=0.25d0

!xc=-0.25d0
!yc=0.25d0

!xd=0.25d0
!yd=-0.25d0
do i=1,N
do j=1,M
sig(i,j)=0.5d0*(1.d0-tanh(cth*(((x(i)-xa)**2.d0+(y(j)-ya)**2.d0)**0.5d0-lct))) + &
0.5d0*(1.d0-tanh(cth*(((x(i)-xb)**2.d0+(y(j)-yb)**2.d0)**0.5d0-lct)))  !+ &
!0.5d0*(1.d0-tanh(cth*(((x(i)-xc)**2.d0+(y(j)-yc)**2.d0)**0.5d0-lct))) + &
!0.5d0*(1.d0-tanh(cth*(((x(i)-xd)**2.d0+(y(j)-yd)**2.d0)**0.5d0-lct)))
!sig(i,j)=1.d0-tanh(9.d0*((x(i)**2.d0+y(j)**2.d0)**0.5d0-2d-1))
!sig(i,j)=(tanh(cth*(x(i)+lct))-tanh(cth*(x(i)-lct)))*(tanh(cth*(y(j)+lct))- &
!tanh(cth*(y(j)-lct)))
!sig(i,j)=0.d0
enddo
enddo

tsig(1)=0.d0
          do i=1,N
          do j=1,M
          tsig(1)=tsig(1)+dx*dy*sig(i,j)
          enddo
          enddo

!!! tolerence value for convergenece
tol=5d-7
sigtol=1d-7
!utol=1d-10
ur=7d-1
!! coefficients of Poissons solver
 		 azeE=1.d0/dx**2.d0
         azeW=1.d0/dx**2.d0
         azeN=1.d0/dy**2.d0
         azeS=1.d0/dy**2.d0
         
          azE=1.d0/dx**2.d0
          azW=1.d0/dx**2.d0
          azN=1.d0/dy**2.d0
          azS=1.d0/dy**2.d0
          azP=-(2.d0/dx**2.d0+2.d0/dy**2.d0)
		asiP=-(2.d0/dx**2.d0+2.d0/dy**2.d0+1.d0/dt)
        
        laE=1.d0/dx**2.d0
          laW=1.d0/dx**2.d0
          laN=1.d0/dy**2.d0
          laS=1.d0/dy**2.d0
          laP=-2.d0*(1.d0/dx**2.d0+1.d0/dy**2.d0)
          
          uE=1.d0/(dx**2.d0)
   		 uW=1.d0/(dx**2.d0)
   		 uN=1.d0/(dy**2.d0)
   		 uS=1.d0/(dy**2.d0)
   		 uP=-2.d0*(1.d0/(dx**2.d0)+1.d0/(dy**2.d0))




  	fc=0.d0

do k=1,P
!if (k>1) then   ! open 007          ! just to write file

  	if(k>1)then
  	do i=1,N
  	do j=1,M
	sigp(i,j)=sig(i,j)
	zp(i,j)=z(i,j)
	fc=1.d0
	 enddo
	 enddo
	endif
	itn=1
	err=5.d0
	!if (k<3)then
	!dt=0.05d0
	!else
	!dt=dt1
	!endif
	
	!asiP=-(2.d0/dx**2.d0+2.d0/dy**2.d0+1.d0/dt)
        do while (err>tol .OR. itn<5)
         !do itn=1,20
        itn=itn+1
  !=====================zeta and z iteration=========================================
  !============================================================================
         zerr=0.d0
        
          
          
         do i=1,N
         do j=1,M
         zetaold(i,j)=zeta(i,j)
         enddo
         enddo
         
          do i=1,N
         do j=1,M
         zold(i,j)=z(i,j)
         enddo
         enddo
         
          
         
         do i=1,N
         do j=1,M
         azeP(i,j)=-(2.d0/dx**2.d0+2.d0/dy**2.d0)
         rZe=(sig(i+1,j)-2.d0*sig(i,j)+sig(i-1,j))/dx**2.d0+   &
        (sig(i,j+1)-2.d0*sig(i,j)+sig(i,j-1))/dy**2.d0
         Sze(i,j)=D5*rZe+D8*shst(i,j) +  &
    D6*zeta(i,j)*(2.d0*fc*sig(i,j)*(DLOG(abs(sig(i,j)/bsig_s))-1.d0)+  & 
    D7*(sig(i,j))**2.d0)+  &
    D8*zeta(i,j)*lam(i,j)
         enddo
         enddo
         
         zcon=10.d0
         do while (zcon>sigtol)  ! open 552
         zcon=0.d0
         zecon=0.d0
         !zeta bc
         
         do j=1,M
          zeta(0,j)=4.d0*z(2,j)/(dx**2.d0)-zeta(2,j)
          zeta(N+1,j)=4.d0*z(N-1,j)/(dx**2.d0)-zeta(N-1,j)
          !zeta(1,j)=0.d0
          !zeta(N,j)=0.d0
          enddo
          do i=1,N
          zeta(i,0)=4.d0*z(i,2)/(dy**2.d0)-zeta(i,2)
          zeta(i,M+1)=4.d0*z(i,M-1)/(dy**2.d0)-zeta(i,M-1)
          !zeta(i,1)=0.d0
          !zeta(i,M)=0.d0
          enddo
         
         
         
         
         do i=1,N
         do j=1,M
         zeo=zeta(i,j)
         zenu=Sze(i,j)-azeE*zeta(i+1,j)-   &
         azeW*zeta(i-1,j)-azeN*zeta(i,j+1)-azeS*zeta(i,j-1)
         zeta(i,j)=(1.d0-ur)*zeta(i,j)+ur*zenu/azeP(i,j)
         !zera=zeo-zeta(i,j)
         !zer=ABS(zera)
         if (abs(zeo-zeta(i,j))>zcon)then
         zcon=zeo-zeta(i,j)
         end if


         enddo
         enddo
  

!z iteration

          do i=2,N-1
          do j=2,M-1

          Sz(i,j)=zeta(i,j)
          enddo
          enddo
		
          do i=2,N-1
          do j=2,M-1
          zo=z(i,j)
         znu=Sz(i,j)-azE*z(i+1,j)-azW*z(i-1,j)-        &
        azN*z(i,j+1)-azS*z(i,j-1)
          z(i,j)=(1.d0-ur)*z(i,j)+ur*znu/azP
          !zr=abs(z(i,j)-zo)
          if (abs(z(i,j)-zo)>zcon)then
          zcon=abs(z(i,j)-zo)
          endif
          enddo
          enddo
          test=-1.2d0
   !print*, "time step",k, "z iteration error is",zcon,"fc",fc,"D3",D3

     enddo        ! close 552
     
     ! now check convergence
     err=0.d0
     err_zeta=0.d0
        do i=1,N
        do j=1,M
        if(abs(zeta(i,j)-zetaold(i,j))>err_zeta) then
        err_zeta=abs(zeta(i,j)-zetaold(i,j))
        endif
        
        enddo
        enddo
        err=err_zeta
        err_z=0.d0
          do i=2,N-1
        do j=2,M-1
        if(abs(z(i,j)-zold(i,j))>err_z) then
        err_z=abs(z(i,j)-zold(i,j))
        endif
        
        enddo
        enddo
        
        if (err_z>err)then
        err=err_z
        endif
        
       if (k>1)then
        do i=1,N
        do j=1,M
        wold=w(i,j)
        w(i,j)=(ur)*w(i,j)+(1.d0-ur)*((z(i,j)-zp(i,j))/dt)/D11
        err=max(abs(wold-w(i,j)),err)
        enddo
        enddo
       endif
        
         
     
 
!sig iteration starts==========================================================
!==============================================================================
if (k>1)then

         do i=1,N
  		 do j=1,M
  		 oldsig(i,j)=sig(i,j)
  		 enddo
  		 enddo
          
          !!! Coefficient for Poisson's solver
          do i=1,N
          do j=1,M
          asiE(i,j)=1.d0/dx**2.d0-u(i,j)*D11/(2.d0*dx)
          asiW(i,j)=1.d0/dx**2.d0+u(i,j)*D11/(2.d0*dx)
          asiN(i,j)=1.d0/dy**2.d0-v(i,j)*D11/(2.d0*dy)
          asiS(i,j)=1.d0/dy**2.d0+v(i,j)*D11/(2.d0*dy)   
          enddo
          enddo
         
         do i=1,N
          do j=1,M
Ssi(i,j)=D10*sig(i,j)*(zeta(i,j)*(D6*(2.d0*sig(i,j)*(DLOG(abs(sig(i,j)/bsig_s))-1.d0)+ & 
           D7*(sig(i,j))**2.d0)+  &
           D8*lam(i,j))+D9+D8*shst(i,j))-sigp(i,j)/dt+D11*w(i,j)*zeta(i,j)*sig(i,j) + &
           -D7*((sig(i+1,j)-sig(i-1,j))**2.d0/(4.d0*dx**2.d0)+  & 
           (sig(i,j+1)-sig(i,j-1))**2.d0/(4.d0*dy**2.d0)) + &
           D10*((sig(i+1,j)-sig(i-1,j))*(zeta(i+1,j)-zeta(i-1,j))/(4.d0*dx**2.d0)+  & 
           (sig(i,j+1)-sig(i,j-1))*(zeta(i,j+1)-zeta(i,j-1))/(4.d0*dy**2.d0))  
          enddo
          enddo
         
          
           

 		
          
          sigrr=10.d0
		do while (sigrr>sigtol)
          sigrr=0.d0
          
         
         
          
         !! sig bc 
         do j=1,M
          dzetaw=(zeta(2,j)-zeta(1,j))/dx
          dsigw=(D10*dzetaw*sig(1,j)+D11*u(1,j)*sig(1,j)-Jw(k))/(1.d0+sig(1,j)*D7)
          siggw(j)=sig(2,j)-2.d0*dx*dsigw
          sig(0,j)=siggw(j)
          enddo



          do i=1,N
          dzetas=(zeta(i,2)-zeta(i,1))/dy
          dsigs=(D10*dzetas*sig(i,1)+D11*v(i,1)*sig(i,1)-Js(k))/(1.d0+sig(i,1)*D7)
          siggs(i)=sig(i,2)-2.d0*dy*dsigs
          sig(i,0)=siggs(i)
          enddo





          do j=1,M
          dzetae=(zeta(N,j)-zeta(N-1,j))/dx
          dsige=(D10*dzetae*sig(N,j)+D11*u(N,j)*sig(N,j)-Je(k))/(1.d0+sig(N,j)*D7)
          sigge(j)=sig(N-1,j)+2.d0*dx*dsige
          sig(N+1,j)=sigge(j)
          enddo


          do i=1,N
          dzetan=(zeta(i,M)-zeta(i,M-1))/dy
          dsign=(D10*dzetan*sig(i,M)+D11*sig(i,M)*v(i,M)-Jn(k))/(1.d0+sig(i,M)*D7)
          siggn(i)=sig(i,M-1)+2.d0*dy*dsign
          sig(i,M+1)=siggn(i)
          enddo

         
  		  do i=1,N
          do j=1,M
          sigo=sig(i,j)
          sig(i,j)=(1.d0-ur)*sig(i,j)+ur*(Ssi(i,j)-asiE(i,j)*sig(i+1,j)-      &
          asiW(i,j)*sig(i-1,j)-asiN(i,j)*sig(i,j+1)-asiS(i,j)*sig(i,j-1))/asiP
          sigr(i,j)=abs(sig(i,j)-sigo)
          if (sigr(i,j)>sigrr)then
          sigrr=sigr(i,j)
          endif
          enddo
          enddo

!print*, "time step",k, "sig iteration error is",sigrr,"sig",sig(1,10)
       enddo ! end of sigtol 
       
       do i=1,N
    	 do j=1,M
    	 if(abs(oldsig(i,j)-sig(i,j))>err)then
    	 err=abs(oldsig(i,j)-sig(i,j))
    	 endif
    	 enddo
    	 enddo  
endif ! end of  checking for solving for diffusion equation
       
		 
       

 

  !!! Extrapolated variable used velocity iteration       
          do i=1,N
          zeta(i,0)=2.d0*zeta(i,1)-zeta(i,2)
          zeta(i,M+1)=2.d0*zeta(i,M)-zeta(i,M-1)
           sig(i,0)=2.d0*sig(i,1)-sig(i,2)
          sig(i,M+1)=2.d0*sig(i,M)-sig(i,M-1)
          w(i,0)=2.d0*w(i,1)-w(i,2)
          w(i,M+1)=2.d0*w(i,M)-w(i,M-1)
          enddo
          
          do j=1,M
          zeta(0,j)=2.d0*zeta(1,j)-zeta(2,j)
          zeta(N+1,j)=2.d0*zeta(N,j)-zeta(N-1,j)
         sig(0,j)=2.d0*sig(1,j)-sig(2,j)
          sig(N+1,j)=2.d0*sig(N,j)-sig(N-1,j)
           w(0,j)=2.d0*w(1,j)-w(2,j)
          w(N+1,j)=2.d0*w(N,j)-w(N-1,j)
          enddo
          
   ! lam iteration starts========================================================
 !=============================================================================
          do i=1,N
  		 do j=1,M
  		 oldlam(i,j)=lam(i,j)
  		 enddo
  		 enddo       
          
          
          do i=2,N-1
          do j=2,M-1
 Lzew=(w(i+1,j)*zeta(i+1,j)-2.d0*w(i,j)*zeta(i,j)+w(i-1,j)*zeta(i-1,j))/(dx**2.d0)+  &
 (w(i,j+1)*zeta(i,j+1)-2.d0*w(i,j)*zeta(i,j)+w(i,j-1)*zeta(i,j-1))/(dy**2.d0)
 Lzeta=(zeta(i+1,j)-2.d0*zeta(i,j)+zeta(i-1,j))/(dx**2.d0)+  &
         (zeta(i,j+1)-2.d0*zeta(i,j)+zeta(i,j-1))/(dy**2.d0)
         hxx=(z(i+1,j)-2.d0*z(i,j)+z(i-1,j))/(dx**2.d0)
          hyy=(z(i,j+1)-2.d0*z(i,j)+z(i,j-1))/(dy**2.d0)
          hxy=(z(i+1,j+1)-z(i+1,j-1)-z(i-1,j+1)+z(i-1,j-1))/(4.d0*dx*dy)
          wxx=(w(i+1,j)-2.d0*w(i,j)+w(i-1,j))/(dx**2.d0)
          wyy=(w(i,j+1)-2.d0*w(i,j)+w(i,j-1))/(dy**2.d0)
          wxy=(w(i+1,j+1)-w(i+1,j-1)-w(i-1,j+1)+w(i-1,j-1))/(4.d0*dx*dy)
          wbc=hxx*wxx+2.d0*hxy*wxy+hyy*wyy
          zetadx=(zeta(i+1,j)-zeta(i-1,j))/(2.d0*dx)
          zetady=(zeta(i,j+1)-zeta(i,j-1))/(2.d0*dy)
          wdx=(w(i+1,j)-w(i-1,j))/(2.d0*dx)
          wdy=(w(i,j+1)-w(i,j-1))/(2.d0*dy)
          sigdx=(sig(i+1,j)-sig(i-1,j))/(2.d0*dx)
          sigdy=(sig(i,j+1)-sig(i,j-1))/(2.d0*dy)
          lgsigdx=(DLOG(abs(sig(i+1,j)))-DLOG(abs(sig(i-1,j))))/(2.d0*dx)
          lgsigdy=(DLOG(abs(sig(i,j+1)))-DLOG(abs(sig(i,j-1))))/(2.d0*dy)
          
          Lsig=(sig(i+1,j)-2.d0*sig(i,j)+sig(i-1,j))/dx**2.d0+  &
         (sig(i,j+1)-2.d0*sig(i,j)+sig(i,j-1))/dy**2.d0
        Lsigsq=((sig(i+1,j))**2.d0-2.d0*(sig(i,j))**2.d0+(sig(i-1,j))**2.d0)/dx**2.d0+ &
         ((sig(i,j+1))**2.d0-2.d0*(sig(i,j))**2.d0+(sig(i,j-1))**2.d0)/dy**2.d0
          Sla(i,j)= (D1*zeta(i,j)-fc*D3*DLOG(abs(sig(i,j)/bsig_s)))*Lsig +   & 
          D1*(zetadx*sigdx+zetady*sigdy)- fc*D3*(lgsigdx*sigdx+lgsigdy*sigdy) -  &
    (5d-1)*(D2)*Lsigsq-2.d0*Lzew+4.d0*(zetadx*wdx+zetady*wdy)+2.d0*w(i,j)*Lzeta+ &
        2.d0*wbc
          end do
          end do
          
          
          
          do i=1,N !open 51
          do j=1,M ! open 52
          if (i==1 .OR. i==N)then
          ddx=dx/2.d0
          else
          ddx=dx
          endif
          if (j==1 .OR. j==M)then
          ddy=dy/2.d0
          else
          ddy=dy
          endif
          g0(i,j)=Sla(i,j)*ddx*ddy
          enddo    ! close 52
          enddo    ! close 51
          
          lamrr=10.d0
          lamit=0
          ptc=0
          
          call cpu_time(start)
          ilg=1
          jlg=1
          do i=1,N   ! open 53
          if (i==1 .OR. i==N) then !open 53a
          do j=1,M    ! open 54
          lam(i,j)=1.d0
          do ii=2,N-1    !open 55
          do jj=2,M-1    ! open 56
          !if (ii/=i .OR. jj/=j)then
          xt=(x(i)-x(ii))
          yt=(y(j)-y(jj))
          !r=(xt**2.d0+yt**2.d0)**(0.5d0)
          lam(i,j)=lam(i,j)+(1.d0/(2.d0*PI))*g0(ii,jj)*xlg(ilg)
          ilg=ilg+1
          !endif
          
          enddo     ! close 56
          enddo    ! close 55
          
          
          enddo    ! close 54
          endif ! close 53a
          enddo    ! close 53
          !!---------- rest two boundaries
          do j=1,M   ! open 53
          if (j==1 .OR. j==M) then !open 53a
          do i=2,N-1     ! open 54
          lam(i,j)=1.d0
          do ii=2,N-1    !open 55
          do jj=2,M-1    ! open 56
          !if (ii/=i .OR. jj/=j)then
          xt=(x(i)-x(ii))
          yt=(y(j)-y(jj))
          !r=(xt**2.d0+yt**2.d0)**(0.5d0)
          lam(i,j)=lam(i,j)+(1.d0/(2.d0*PI))*g0(ii,jj)*ylg(jlg)
          jlg=jlg+1
          !end if   !! if you change anything dont forget to change in log calc
          
          enddo     ! close 56
          enddo    ! close 55
          
          
          enddo    ! close 54
          endif ! close 53a
          enddo    ! close 53
    call cpu_time(finish)
          
          lamrr=10.d0
          lamit=0
          ptc=0
          
  		 do while (lamrr>sigtol)

          lamrr=0.d0


          do i=2,N-1
          do j=2,M-1
          lamo=lam(i,j)
          lamnu=Sla(i,j)-laE*lam(i+1,j)-       &
          laW*lam(i-1,j)-laN*lam(i,j+1)-laS*lam(i,j-1)
          lam(i,j)=(1.d0-ur)*lam(i,j)+ur*lamnu/laP
          lamr=abs(lam(i,j)-lamo)
          if (lamr>lamrr)then
          lamrr=lamr
          endif
          enddo
          enddo
          lamit=lamit+1
          ptc=ptc+1
        ! if (ptc>10)then
         !print*, "time step",k, "lam iteration error is",lamrr,"lam",lam(1,1)
         ! ptc=0
        !endif
          enddo ! lam iteration stops here        
          err_lam=0.d0
         do i=1,N
    	 do j=1,M
    	 if(abs(oldlam(i,j)-lam(i,j))>err_lam)then
    	 err_lam=abs(oldlam(i,j)-lam(i,j))
    	 endif
    	 enddo
    	 enddo  
          
          if(err_lam>err)then
          err=err_lam
          endif
          
         
!!==============================================================================
!!============================velocity iteration===============================
if(k>1)then  ! start v iteration

!!! calculating SU and SV




velerr=10.d0
  do i=1,N
  lam(i,M+1)=2.d0*lam(i,M)-lam(i,M-1)
  lam(i,0)=2.d0*lam(i,1)-lam(i,2)
  z(i,0)=z(i,2)
  z(i,M+1)=z(i,M-1)
  enddo


  do j=1,M
  lam(N+1,j)=2.d0*lam(N,j)-lam(N-1,j)
  lam(0,j)=2.d0*lam(1,j)-lam(2,j)
  z(0,j)=z(2,j)
  z(N+1,j)=z(N-1,j)
  enddo

do i=1,N
   do j=1,M
   		  zetadx=(zeta(i+1,j)-zeta(i-1,j))/(2.d0*dx)
          zetady=(zeta(i,j+1)-zeta(i,j-1))/(2.d0*dy)
          wdx=(w(i+1,j)-w(i-1,j))/(2.d0*dx)
          wdy=(w(i,j+1)-w(i,j-1))/(2.d0*dy)
          wzedx=(zeta(i+1,j)*w(i+1,j)-zeta(i-1,j)*w(i-1,j))/(2.d0*dx)
          wzedy=(zeta(i,j+1)*w(i,j+1)-zeta(i,j-1)*w(i,j-1))/(2.d0*dy)
          hxx=(z(i+1,j)-2.d0*z(i,j)+z(i-1,j))/(dx**2.d0)
          hyy=(z(i,j+1)-2.d0*z(i,j)+z(i,j-1))/(dy**2.d0)
          hxy=(z(i+1,j+1)-z(i+1,j-1)-z(i-1,j+1)+z(i-1,j-1))/(4.d0*dx*dy)
   sigdx=(sig(i+1,j)-sig(i-1,j))/(2.d0*dx)
   sigdy=(sig(i,j+1)-sig(i,j-1))/(2.d0*dy)
   sigsqdx=((sig(i+1,j))**2.d0-(sig(i-1,j))**2.d0)/(2.d0*dx)
   sigsqdy=((sig(i,j+1))**2.d0-(sig(i,j-1))**2.d0)/(2.d0*dy)
   Su(i,j)= (D1*zeta(i,j)-D3*DLOG(abs(sig(i,j)/bsig_s)))*sigdx-(5d-1)*(D2)*sigsqdx -  &
   (lam(i+1,j)-lam(i-1,j))/(2.d0*dx)+2.d0*w(i,j)*zetadx+ &
  2.d0*(wdx*hxx+wdy*hxy) -wzedx
   Sv(i,j)= (D1*zeta(i,j)-D3*DLOG(abs(sig(i,j)/bsig_s)))*sigdy-(5d-1)*(D2)*sigsqdy - &
    (lam(i,j+1)-lam(i,j-1))/(2.d0*dy)+2.d0*w(i,j)*zetady+ &
  2.d0*(wdx*hxy+wdy*hyy) -wzedy
   enddo
   enddo



  do i=1,N
  do j=1,M
  oldu(i,j)=u(i,j)
  oldv(i,j)=v(i,j)
  enddo
  enddo





do while (velerr>sigtol) ! open 59
		 do i=1,N
          do j=1,M
          uold(i,j)=u(i,j)
          vold(i,j)=v(i,j)
          enddo
          enddo


!!! the vel bc was here
do i=1,N
  u(i,0)=u(i,1)
  v(i,0)=v(i,1)
  u(i,M+1)=u(i,M)
  v(i,M+1)=v(i,M)
  enddo


  do j=1,M
  u(0,j)=u(1,j)
  v(0,j)=v(1,j)
  u(N+1,j)=u(N,j)
  v(N+1,j)=v(N,j)
  enddo
  


   
   
velit=0
urr=10.d0
do while (urr>sigtol)




  

urr=0.d0
  do i=1,N
  do j=1,M
   uo=u(i,j)
   u(i,j)=(1.d0-ur)*u(i,j)+ur*(Su(i,j)-uE*u(i+1,j)-      &
   uW*u(i-1,j)-uN*u(i,j+1)-uS*u(i,j-1))/uP
   uer=abs(u(i,j)-uo)
   if (uer>urr)then
   urr=uer
   endif
 enddo
enddo
  !if (urr>err) then
   !err=urr
  !endif

  vrr=0.d0
    do i=1,N
    do j=1,M
     vo=v(i,j)
     v(i,j)=(1.d0-ur)*v(i,j)+ur*(Sv(i,j)-uE*v(i+1,j)-      &
     uW*v(i-1,j)-uN*v(i,j+1)-uS*v(i,j-1))/uP
     ver=abs(v(i,j)-vo)
     if (ver>vrr)then
     vrr=ver
     endif
   enddo
   enddo
   if (vrr>urr) then
     urr=vrr
    endif
    
    
     
    
    velit=velit+1
    
    !print*, "time step",k, "vel it error is",urr
  enddo  !end of do while for u and v
   
    dialn=0.d0
usum=0.d0
vsum=0.d0
  do i=1,N
  do j=1,M
if (i>1 .AND. i<N) then
ddx=dx
else
ddx=dx/2.d0
endif
if (j>1 .AND. i<M) then
ddy=dy
else
ddy=dy/2.d0
endif

  usum=usum+u(i,j)*ddx*ddy
  vsum=vsum+v(i,j)*ddx*ddy
  dialn=dialn+w(i,j)*zeta(i,j)*ddx*ddy
  enddo
  enddo
  umean=usum/(Lx*Ly)
  vmean=vsum/(Lx*Ly)
  
  do i=1,N
  do j=1,M
  u(i,j)=u(i,j)-umean
  v(i,j)=v(i,j)-vmean
  enddo
  enddo
  
  velerr=0.d0
  do i=1,N
  do j=1,M
  if(abs(u(i,j)-uold(i,j))>velerr)then
  velerr=abs(u(i,j)-uold(i,j))
  endif
  
  enddo
  enddo
  do i=1,N
  do j=1,M
  if(abs(v(i,j)-vold(i,j))>velerr)then
  velerr=abs(v(i,j)-vold(i,j))
  endif
  
  enddo
  enddo
  
  !print*, "time step",k, "velerr",velerr
  enddo     ! close 59
    
    do i=1,N
    do j=1,M
    if(abs(oldu(i,j)-u(i,j))>err)then
    err=abs(oldu(i,j)-u(i,j))
    endif
    enddo
    enddo
    
    do i=1,N
    do j=1,M
    if(abs(oldv(i,j)-v(i,j))>err)then
    err=abs(oldv(i,j)-v(i,j))
    endif
    enddo
    enddo
    
    
    
z(0,0)=z(2,2)
z(N+1,M+1)=z(N-1,M-1)
z(0,M+1)=z(2,M-1)
z(N+1,0)=z(N-1,2)

do i=1,N
do j=1,M
d2zx=(z(i+1,j)-2.d0*z(i,j)+z(i-1,j))/(dx**2.d0)
d2zy=(z(i,j+1)-2.d0*z(i,j)+z(i,j-1))/(dy**2.d0)
d2zxy=(z(i+1,j+1)-z(i+1,j-1)-z(i-1,j+1)+z(i-1,j-1))/(4.d0*dx*dy)
dux=(u(i+1,j)-u(i-1,j))/(2.d0*dx)
dvy=(v(i,j+1)-v(i,j-1))/(2.d0*dy)
duy=(u(i,j+1)-u(i,j-1))/(2.d0*dy)
dvx=(v(i+1,j)-v(i-1,j))/(2.d0*dx)
shst(i,j)=2.d0*(d2zx*dux+d2zy*dvy+d2zxy*(duy+dvx))
enddo
enddo




endif   ! end z iteration
         
          
       
!-======== div cal
divm=0.d0
   do i=1,N
   do j=1,M
   div(i,j)=(u(i+1,j)-u(i-1,j))/(2.d0*dx)+(v(i,j+1)-v(i,j-1))/(2.d0*dy)-zeta(i,j)*w(i,j)
   if(abs(div(i,j))>=divm)then
   divm=abs(div(i,j))
   endif
   enddo
   enddo
 !-======== div cal
    
diln2=(v(1,M)+v(N,M)-v(1,1)-v(N,1))*dx/2.d0+ (u(N,1)+u(N,M)-u(1,1)-u(1,M))*dy/2.d0
do i=2,N-1
diln2=diln2+v(i,M)*dx-v(i,1)*dx
enddo     
do j=2,M-1
diln2=diln2+u(N,j)*dy-u(1,j)*dy
enddo

         pvti=1+(N-1)/2
         pvtj=1+(M-1)/2
print*,"time step",k, "itn no",itn,"error is",err,"z",z(pvti,pvtj),"v",u(1,1),"div",divm,&
"w",w(pvti,pvtj),"diln2 ", diln2,"dialn",dialn
          !itn=itn+1
 

enddo ! for do while ending

 d2h(1,1)=((z(2,1)-z(1,1))/dx)**2.d0+  &
         ((z(1,2)-z(1,1))/dy)**2.d0
          d2h(N,1)=((z(N,1)-z(N-1,1))/dx)**2.d0+   &
        ((z(N,2)-z(N,1))/dy)**2.d0
          d2h(N,M)=((z(N,M)-z(N-1,M))/dx)**2.d0+   &
         ((z(N,M)-z(N,M-1))/dy)**2.d0
          d2h(1,M)=((z(2,M)-z(1,M))/dx)**2.d0+  &
         ((z(1,M)-z(1,M-1))/dy)**2.d0
          do i=1,N-1
          d2h(i,1)=((z(i+1,1)-z(i-1,1))/(2.d0*dx))**2.d0+ &
         ((z(i,2)-z(i,1))/dy)**2.d0
          d2h(i,M)=((z(i+1,M)-z(i-1,M))/(2*dx))**2.d0+ &
         ((z(i,M)-z(i,M-1))/dy)**2.d0
          enddo
          do j=1,M-1
          d2h(1,j)=((z(2,j)-z(1,j))/(dx))**2.d0+ &
        ((z(1,j+1)-z(1,j-1))/(2.d0*dy))**2.d0
          d2h(N,j)=((z(N,j)-z(N-1,j))/(dx))**2.d0+ &
         ((z(N,j+1)-z(N,j-1))/(2.d0*dy))**2.d0
          enddo
          do i=2,N-1
          do j=2,M-1
          d2h(i,j)=((z(i+1,j)-z(i-1,j))/(2.d0*dx))**2.d0+ &
          ((z(i,j+1)-z(i,j-1))/(2.d0*dy))**2.d0
          enddo
          enddo

 tsig(k)=0.d0
          do i=1,N
  do j=1,M
if (i>1 .AND. i<N) then
ddx=dx
else
ddx=dx/2.d0
endif
if (j>1 .AND. i<M) then
ddy=dy
else
ddy=dy/2.d0
endif
          tsig(k)=tsig(k)+dx*dy*sig(i,j)
          enddo
          enddo




!endif  ! close 007


   
open(unit=27,file="zeta.txt")
         do i=1,N
         write(27,32)(zeta(i,j),j=1,M)
         enddo
  32     format(1x,5000f16.8)

      open(unit=28,file="sig.txt")
         do i=1,N
         write(28,33)(sig(i,j),j=1,M)
         enddo
  33     format(1x,5000f16.8)

  open(unit=29,file="z.txt")
         do i=1,N
         write(29,34)(z(i,j),j=1,M)
         enddo
  34     format(1x,5000f16.8)

  open(unit=30,file="lam.txt")
         do i=1,N
         write(30,35)(lam(i,j),j=1,M)
         enddo
  35     format(1x,5000f16.8)
  open(unit=31,file="u.txt")
         do i=1,N
         write(31,36)(u(i,j),j=1,M)
         enddo
  36     format(1x,5000f16.8)
  open(unit=32,file="v.txt")
         do i=1,N
         write(32,37)(v(i,j),j=1,M)
         enddo
  37     format(1x,5000f16.8)
  
  
  open(unit=36,file="div.txt")
         do i=1,N
         write(36,41)(div(i,j),j=1,M)
         enddo
  41     format(1x,5000f16.8)
  open(unit=56,file="w.txt")
         do i=1,N
         write(56,61)(w(i,j),j=1,M)
         enddo
  61     format(1x,5000f16.8)
  open(unit=69,file="zmax.txt")
         write(69,74)(z(pvti,pvtj))
  74     format(1x,5000f16.8)


enddo  ! time loop ends

open(unit=33,file="totalsig.txt")
         do i=1,P
         write(33,38)(tsig(i))
         enddo
  38     format(1x,5000f16.8)
  open(unit=37,file="t.txt")
         do i=1,P
         write(37,42)(t(i))
         enddo
  42     format(1x,5000f16.8)
open(unit=38,file="x.txt")
         do i=1,N
         write(38,43)(x(i))
         enddo
  43     format(1x,5000f16.8)
   open(unit=39,file="y.txt")
         do j=1,M
         write(39,44)(y(j))
         enddo
  44     format(1x,5000f16.8)


    

!print*,"lamP=",laP,"z=",z(10,10),"divm",divm

end program sumingup
