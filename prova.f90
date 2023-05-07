PROGRAM MAIN                                                     
IMPLICIT REAL*8(A-H,O-Z)
dimension xr(1000),frev(1000), alphavec(4)
a0=0.00433  
n1=700 
step=0.015 
aa=1000000. 
time=0.00005 
alphavec(1)=0.3
alphavec(2)=0.4
alphavec(3)=0.5
alphavec(4)=0.8
iter=70000
open(1,file='den.dat')
do j=1,4
    alpha=alphavec(j)
    pi=4.d0*DATAN(1.d0)
    piin=1.d0/(4.d0*pi)
    pi2in=dsqrt(piin)
    alpha2=alpha*alpha
    cvar=2.d0*dsqrt(alpha)**3/dsqrt(dsqrt(pi))

    do i=1,n1
    xr(i)=step*dfloat(i-1)
    xr2=xr(i)*xr(i)
    frev(i)=cvar*xr(i)*dexp(-0.5d0*alpha2*xr2)
    write(1,*) xr(i),frev(i)
    enddo
    write(1,*) ' '
    write(1,*) ' '
enddo
close(1)
END program MAIN