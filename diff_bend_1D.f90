!-----------------------------------------------------------------------------------------------------------------
! Numerical simulation of in one dimension with diffusion of curvature induced-proteins
!
! Simulation code by: Arijit Mahapatra, David Saintillan, Padmini Ranganami (UCSD)
!
! For details of the model, refer to:
! "Transport phenomena in fluid films with curvature elasticity" by Mahapatra et al.
! (https://arxiv.org/abs/2001.07539)
!-----------------------------------------------------------------------------------------------------------------
program sumingup
implicit none
integer:: N,M,P,i,j,k,ptc,itn,lamit,pvti,pvtj,velit,ii,jj,kN
double precision:: sig0, lam0, alpha, ell, kp, beta,c,err,fc
double precision:: pr, L, D1, D2, D3, D4, D5, D6, D7, D8, D9
double precision:: D10,lr,tT,lX,dx,dt, bsig_s, kBT, D
double precision:: cth,lct,tol,sigtol,ur,zr,zrr,sigrr
double precision:: zer,zerr,lamrr,lamr,rZe,zenu,znu,vrr
double precision:: signu,dzetae,dzetaw,dzetan,dzetas,zcon,zecon
double precision:: zetadx,zetady,azE,azW,azN,azS,azP
double precision:: sigdx,sigdy,Lsig,Lsigsq,sigo
double precision:: zeo,zo,azeE,azeW,azeN,azeS
double precision:: dsign,dsige,dsigw,dsigs,PI,RHSlam
double precision:: laW,laS,laP,laN,lamo,laE,lamnu,zera,uup,vup,nu,uo,vo,delta,divm
double precision:: asiP,asiE,asiW,zetam,sigm
double precision, dimension(0:500):: x,y,siggs,sigge,siggw,siggn,azeP,Sze,Sz,zetaold,zold 
double precision, dimension(0:500)::Ssi,d2h,sigold
double precision, dimension(0:5000):: t,tsig
double precision, dimension(0:500,0:2000):: sigr,Sla,shst,uold,vold
double precision, dimension(0:500,0:2000):: zeta,z,sig,lam,s
sig0=5d-3
lam0=1d-4
bsig_s=10.d0
ell=5d-1
kp=84.d0
!beta=0.1d0*alpha*sig0
kBT=4.14d-2
D=1e5
!nu=1d-4
pr=0.d0
L=1000.d0
PI=2.D0*DASIN(1.D0)


D1= kp*ell*sig0/(L*lam0)
D2=2.d0*(ell)**2.d0*kp*sig0**2.d0/lam0
D3=kBT*sig0/lam0
!D4=beta/(alpha*sig0)
D5=2.d0*L*ell*sig0
D6=L**2.d0*kBT*sig0/kp
D7=2.d0*kp*ell**2.d0*sig0/kBT
D8=2.d0*L**2.d0*lam0/kp
D9=2.d0*L**3.d0*pr/kp
D10=kp*ell/(kBT*L)

N=201
P=1001
tT=1d0
lX=1.d0
dx=lX/(N-1)
dt=tT/(P-1)
!delta=dx
x(1)=-lX/2.d0
do i=2,N
x(i)=x(i-1)+dx
enddo
t(1)=0.d0
do k=2,P
t(k)=t(k-1)+dt
enddo

cth=20.d0
lct=1d-1

do i=1,N
sig(i,1)=0.5*(tanh(cth*(x(i)+lct))-tanh(cth*(x(i)-lct)))
enddo


do i=0,N+1
do k=1,P
zeta(i,k)=0.d0
sig(i,k)=sig(i,1)
z(i,k)=0.d0
lam(i,k)=1.d0
enddo
enddo


do k=1,P
z(1,k)=0.d0
z(N,k)=0.d0
lam(N,k)=1.d0
!lam(1,k)=1.d0
enddo

! coefficients 

 zerr=0.d0
         azeE=1.d0/dx**2.d0
         azeW=1.d0/dx**2.d0  
          azE=1.d0/dx**2.d0
          azW=1.d0/dx**2.d0
         
          azP=-(2.d0/dx**2.d0)
          
          
          asiP=-(2.d0/dx**2.d0+1.d0/dt)
          asiE=1.d0/dx**2.d0
          asiW=1.d0/dx**2.d0  

fc=0.d0
tol=5d-7
sigtol=1d-7
!utol=1d-10
ur=9d-1
do k=1,P
  	
  	if(k>1)then
  	fc=1.d0
  	do i=0,N+1
	z(i,k)=z(i,k-1)
	zeta(i,k)=zeta(i,k-1)
	sig(i,k)=sig(i,k-1)
	lam(i,k)=lam(i,k-1)
	enddo
  	endif
  	
	itn=1
	err=5.d0
        do while (err>tol .OR. itn<20)
         !do itn=1,20
        itn=itn+1
  !=====================zeta and z iteration=========================================
  !============================================================================
           
         do i=1,N
         zetaold(i)=zeta(i,k)
         enddo
         
          do i=1,N
         zold(i)=z(i,k)
         enddo
         
         do i=2,N-1
         azeP(i)=-(2.d0/dx**2.d0)-  &
        D6*(fc*2.d0*fc*sig(i,k)*(log(abs(sig(i,k)/bsig_s))-1.d0)+  & 
        D7*(sig(i,k))**2.d0)-D8*lam(i,k)
         rZe=(sig(i+1,k)-2.d0*sig(i,k)+sig(i-1,k))/(dx**2.d0)
         Sze(i)=D5*rZe
         enddo
         
         
         zcon=10.d0
         do while (zcon>sigtol)  ! open 552
         zcon=0.d0
         zecon=0.d0
         !zeta bc
        
          zeta(1,k)=2.d0*z(2,k)/(dx**2.d0)
          zeta(N,k)=2.d0*z(N-1,k)/(dx**2.d0)
         
         
         

         
         do i=2,N-1
         zeo=zeta(i,k)
         zenu=Sze(i)-azeE*zeta(i+1,k)-azeW*zeta(i-1,k)
         zeta(i,k)=(1.d0-ur)*zeta(i,k)+ur*zenu/azeP(i)
         zera=zeo-zeta(i,k)
         zer=ABS(zera)
         if (zer>zecon)then
         zecon=zer
         end if
         enddo
  !z iteration

          do i=2,N-1
           Sz(i)=zeta(i,k)
          enddo
         do i=2,N-1
          zo=z(i,k)
         znu=Sz(i)-azE*z(i+1,k)-azW*z(i-1,k)
          z(i,k)=(1.d0-ur)*z(i,k)+ur*znu/azP
          zr=abs(z(i,k)-zo)
          if (zr>zcon)then
          zcon=zr
          endif
          enddo
   


     enddo        ! close 552
     
     ! now check convergence
     zerr=0.d0
        do i=1,N
        if(abs(zeta(i,k)-zetaold(i))>zerr) then
        zerr=abs(zeta(i,k)-zetaold(i))
        endif
        
        enddo
        
        err=zerr
          do i=1,N
      
        if(abs(z(i,k)-zold(i))>err) then
        err=abs(z(i,k)-zold(i))
        endif
        
        enddo
      
        
    
 
!sig iteration starts==========================================================
!==============================================================================
if (k>1)then
   do i=1,N
   sigold(i)=sig(i,k)
   enddo    
           
          sigrr=10.d0
		do while (sigrr>sigtol)
          sigrr=0.d0
          
          do i=1,N
 Ssi(i)=D10*sig(i,k)*(zeta(i,k)*(D6*(2.d0*fc*sig(i,k)*(log(abs(sig(i,k)/bsig_s))-1.d0)+ &
          D7*(sig(i,k))**2.d0)+  &
           D8*lam(i,k))+D9)-sig(i,k-1)/dt + & 
           D10*(zeta(i+1,k)-zeta(i-1,k))*(sig(i+1,k)-sig(i-1,k))/(4.d0*dx**2.d0) - &
           D7*(sig(i+1,k)-sig(i-1,k))*(sig(i+1,k)-sig(i-1,k))/(4.d0*dx**2.d0)
          enddo
      
          dzetaw=(zeta(2,k)-zeta(1,k))/dx
          dsigw=(D1*dzetaw*sig(1,k))/(D2*sig(1,k)+D3)
          sig(0,k)=sig(2,k)-2.d0*dx*dsigw
        
   
          dzetae=(zeta(N,k)-zeta(N-1,k))/dx
          dsige=(D1*dzetae*sig(N,k))/(D2*sig(N,k)+D3)
          sig(N+1,k)=sig(N-1,k)+2.d0*dx*dsige
          
   
 		
  		  do i=1,N
          sigo=sig(i,k)
          sig(i,k)=(1.d0-ur)*sig(i,k)+ur*(Ssi(i)-asiE*sig(i+1,k)-      &
          asiW*sig(i-1,k))/asiP
          if (abs(sig(i,k)-sigo)>sigrr)then
          sigrr=abs(sig(i,k)-sigo)
          endif
          enddo

       enddo ! end of sigtol 
       
       do i=1,N
       if(abs(sig(i,k)-sigold(i))>err)then
       err=abs(sig(i,k)-sigold(i))
       endif
       enddo
endif ! end of  checking for solving for diffusion equation

 ! lam iteration starts========================================================
 !=============================================================================
         
         lamrr=0
    do i=N-1,1,-1
        lamo=lam(i,k)
        zetam=(zeta(i+1,k)+zeta(i,k))/2
        sigm=(sig(i+1,k)+sig(i,k))/2
       RHSlam=((sig(i+1,k)-sig(i,k))/(dx))*(D1*zetam-D2*sigm-  & 
       fc*D3*log(abs(sig(i,k)/bsig_s)))
       lam(i,k)=lam(i+1,k)-RHSlam*dx
       lamr=abs(lam(i,k)-lamo)
       if (lamr>lamrr)then
       lamrr=lamr
       endif
    enddo 
  
    if (lamrr>err)then
      err=lamrr
    endif
     
   
         z(N+1,k)=z(N-1,k)
         z(0,k)=z(2,k)
         

           kN=(N-1)/2+1
          s(kN,k)=0
          
          do i=kN+1,N
             s(i,k)=s(i-1,k)+(dx**2.d0+(z(i,k)-z(i-1,k))**2.d0)**0.5d0
          enddo
          
          do i=kN-1,1,-1
          s(i,k)=s(i+1,k)-(dx**2.d0+(z(i,k)-z(i+1,k))**2)**0.5d0
          enddo
          
        if (k==1)then

         
          
          do i=1,N
            sig(i,1)=0.5d0*(tanh(cth*(s(i,1)+lct))-tanh(cth*(s(i,1)-lct))) 
          enddo
       
	sig(0,k)=2.d0*sig(1,k)-sig(2,k)
   sig(N+1,k)=2.d0*sig(N,k)-sig(N-1,k)


         endif

          tsig(k)=0.d0
          do i=1,N
          tsig(k)=tsig(k)+abs(s(i,k)-s(i-1,k))*sig(i,k)
          enddo


         
          !itn=itn+1
 pvti=1+(N-1)/2
          print*,"time step",k, "itn no",itn,"error is",err,"z",z(pvti,k)

enddo ! for do while ending

enddo  ! time loop ends

   
    

print*,"lamP=",laP,"z=",z(10,10),"divm",divm
open(unit=27,file="zeta.txt")
         do i=1,N
         write(27,32)(zeta(i,j),j=1,P)
         enddo
  32     format(1x,5000f16.8)

      open(unit=28,file="sig.txt")
         do i=1,N
         write(28,33)(sig(i,j),j=1,P)
         enddo
  33     format(1x,5000f16.8)

  open(unit=29,file="z.txt")
         do i=1,N
         write(29,34)(z(i,j),j=1,P)
         enddo
  34     format(1x,5000f16.8)

  open(unit=30,file="lam.txt")
         do i=1,N
         write(30,35)(lam(i,j),j=1,P)
         enddo
  35     format(1x,5000f16.8)
  open(unit=33,file="totalsig.txt")
         do i=1,P
         write(33,38)(tsig(i))
         enddo
  38     format(1x,5000f16.8)
 open(unit=34,file="t.txt")
         do i=1,P
         write(34,39)(t(i))
         enddo
  39     format(1x,5000f16.8)
  open(unit=35,file="x.txt")
         do i=1,N
         write(35,40)(x(i))
         enddo
  40     format(1x,5000f16.8)
end program sumingup
