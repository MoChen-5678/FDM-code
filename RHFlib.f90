!*************************************************************************************************!
!                                                                                                 !
Module RHFlib                                                                                     !
!                                                                                                 !
!*************************************************************************************************!
!
!--- Contain: rtbrent, zriddr, ration, simps, integ2, DET, sdiag
!
!--- Subroutine and function archieve 
!
Use Define
implicit none

    integer, parameter, private :: itmax=100
    double precision, parameter, private :: epss=1.0d-16

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine NUCLEI(IS, NPRO, TE)                                                                   !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!     is = 1 determines the symbol for a given proton number npro
!          2 determines the proton number for a given symbol te
!
!-----------------------------------------------------------------------
!
    INTEGER :: IS, NPRO, NP
    CHARACTER, DIMENSION(2) :: TE
    INTEGER, PARAMETER :: MAXZ = 210 

    CHARACTER, DIMENSION(2*MAXZ+2) :: T

    DATA T(  1: 40)/ ' ',' ',' ','H','H','e','L','i','B','e',' ','B',' ','C',' ','N',' ','O',' ','F', & 
                   & 'N','e','N','a','M','g','A','l','S','i',' ','P',' ','S','C','l','A','r',' ','K' /
    DATA T( 41: 80)/ 'C','a','S','c','T','i',' ','V','C','r','M','n','F','e','C','o','N','i','C','u', & 
                   & 'Z','n','G','a','G','e','A','s','S','e','B','r','K','r','R','b','S','r',' ','Y' /
    DATA T( 81:120)/ 'Z','r','N','b','M','o','T','c','R','u','R','h','P','d','A','g','C','d','I','n', & 
                   & 'S','n','S','b','T','e',' ','I','X','e','C','s','B','a','L','a','C','e','P','r' /
    DATA T(121:160)/ 'N','d','P','m','S','m','E','u','G','d','T','b','D','y','H','o','E','r','T','m', & 
                   & 'Y','b','L','u','H','f','T','a',' ','W','R','e','O','s','I','r','P','t','A','u' /
    DATA T(161:200)/ 'H','g','T','l','P','b','B','i','P','o','A','t','R','n','F','r','R','a','A','c', & 
                   & 'T','h','P','a',' ','U','N','p','P','u','A','m','C','m','B','k','C','f','E','s' /
    DATA T(201:240)/ 'F','m','M','d','N','o','L','r','R','f','D','b','S','g','B','h','H','s','M','t', & 
                   & 'A','0','A','1','A','2','A','3','A','4','A','5','A','6','A','7','A','8','A','9' /
    DATA T(241:280)/ 'B','0','B','1','B','2','B','3','B','4','B','5','B','6','B','7','B','8','B','9', & 
                   & 'C','0','C','1','C','2','C','3','C','4','C','5','C','6','C','7','C','8','C','9' /
    DATA T(281:320)/ 'D','0','D','1','D','2','D','3','D','4','D','5','D','6','D','7','D','8','D','9', & 
                   & 'E','0','E','1','E','2','E','3','E','4','E','5','E','6','E','7','E','8','E','9' /
    DATA T(321:360)/ 'F','0','F','1','F','2','F','3','F','4','F','5','F','6','F','7','F','8','F','9', & 
                   & 'G','0','G','1','G','2','G','3','G','4','G','5','G','6','G','7','G','8','G','9' /
                                 !
    IF (IS.EQ.1) THEN
        IF (NPRO.LT.0.OR.NPRO.GT.MAXZ) STOP 'IN NUCLEUS: NPRO WRONG' 
        TE = T(2*NPRO+1:2*NPRO+2)
        RETURN
    ELSE
        DO NP = 0, MAXZ
            IF (TE(1) .EQ. T(2*NP+1) .AND. TE(2).EQ.T(2*NP+2)) THEN
                NPRO = NP
                RETURN
            ENDIF
         ENDDO
         WRITE(6,100) TE
  100    FORMAT(//,' NUCLEUS ',A2,'  UNKNOWN')
    end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
END Subroutine NUCLEI                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
SUBROUTINE Gauleg(x1,x2,x,w,n)                                                                    !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Given the lower and upper limits of integration x1 and x2, and given n
!--- this routine returns array x(1:n) and w(1:n) of length n, containing the 
!--- abscissas and weight of the Gauss-Legendr n-point quadrature formula. 
!
    INTEGER, intent(in) :: n
    double precision, intent(in) :: x1,x2
    double precision, dimension(n), intent(out) :: x,w
    DOUBLE PRECISION, parameter :: EPS=3.d-14

    INTEGER :: i,j,m
    DOUBLE PRECISION :: p1,p2,p3,pp,xl,xm,z,z1

    m   = (n+1)/2
    xm  = 0.5d0*(x2+x1)
    xl  = 0.5d0*(x2-x1)
    do i= 1, m
        z           = cos(pi*(i -.25d0)/(n+.5d0))
1       continue
        p1          = 1.d0
        p2          = 0.d0
        do j= 1,n
            p3  = p2
            p2  = p1
            p1  = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do
        pp          = n*(z*p1-p2)/(z*z-1.d0)
        z1          = z
        z           = z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)        = xm-xl*z
        x(n+1-i)    = xm+xl*z
        w(i)        = 2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)    = w(i)
    end do
    return
END Subroutine Gauleg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Function DOT(x1, x2, cosx)                                                                        !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the dot product of vectors x1 and x2, and cosin of separation angle
!
    double precision :: DOT, DOTP
    double precision, dimension(3), intent(in) :: x1, x2
    double precision, intent(out) :: cosx
    
    integer :: i
    double precision :: r1, r2
    
    DOTP = zero;     r1  = zero;     r2  = zero
    do i = 1, 3
        DOTP    = DOTP + x1(i)*x2(i)
        r1      = x1(i)**2 + r1
        r2      = x2(i)**2 + r2
    end do
    r1  = dsqrt(r1);     r2  = dsqrt(r2)

    cosx    = DOTP/(r1*r2)
    DOT     = DOTP
    return
    
    DOTP    = (DOT - x1(3)*x2(3))/(r1*r2)
    
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Function DOT                                                                                  !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine VECTP(x1, x2, x)                                                                       !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the vector product of the vectors x1, x2
!
    double precision, dimension(3), intent(in) :: x1, x2
    double precision, dimension(3), intent(out) :: x
    
    x(1)    = x1(2)*x2(3) - x1(3)*x2(2)
    x(2)    = x1(3)*x2(1) - x1(1)*x2(3)
    x(3)    = x1(1)*x2(2) - x1(2)*x2(1)
    
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine VECTP                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Function seconds_to_string( seconds ) result( t )                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        character(len=40) :: t
        double precision  :: seconds
        integer           :: s, m, h, d
        character(len=64) :: f
        double precision  :: sec_per_day = 86400d0
        double precision  :: sec_per_hou = 3600d0
        double precision  :: sec_per_min = 60.0d0
        double precision  :: r
        ! input: time in seconds
        ! output: string with the time formatted in days, hours minutes and seconds
        r = seconds
        d = floor( r / sec_per_day )
        r = r - dble(d)*sec_per_day

        h = floor( r / sec_per_hou )
        r = r - dble(h)*sec_per_hou

        m = floor( r / sec_per_min )
        r = r - dble(m)*sec_per_min
        
        s = floor( r )
        
        f = '(1x,i2," day(s), ",i2.2,":",i2.2,":",i2.2,5x,"(",f12.6,")")'
        write(t,f) d, h, m, s, seconds
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end function seconds_to_string                                                                    !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine deriv(f,df,num,step)                                                                   !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to calcuate the first derivative for
!--- one dimension function f, df is the derivative and num is
!--- the total point number for the arry, step is distance between
!--- two points.

    integer, intent(in) :: num
    double precision, dimension(num) :: f, df
    double precision, intent(in) :: step

    double precision, dimension(5,5) :: A
    double precision :: emfact,sum
    integer :: i,j,k,jj,nmx

    data A(1,1:5)/-50.,96.,-72.,32.,-6./
    data A(2,1:5)/-6.,-20.,36.,-12.,2./
    data A(3,1:5)/2.,-16.,0.,16.,-2./
    data A(4,1:5)/-2.,12.,-36.,20.,6./
    data A(5,1:5)/6.,-32.,72.,-96.,50./
    data emfact/24./

    nmx = num - 2
    do j = 1, num
        k = 3
        if(j.lt.3)      k = j
        if(j.gt.nmx)    k = j - num + 5
        sum = 0.
        do i = 1,5
            jj  = j + i - k
            sum = sum + A(k,i)*f(jj)
        end do
        df(j)  = sum/(step*emfact)
    end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine deriv                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine deriv2(f,df,num,sstep)                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to calcuate the second derivative for
!--- one dimension function f, df is the derivative and num is
!--- the total point number for the arry, step is distance between
!--- two points.

    integer, intent(in) :: num
    double precision, dimension(num) :: f, df
    double precision, intent(in) :: sstep

    double precision, dimension(5,5) :: A
    double precision :: emfact,sum,step
    integer :: i,j,k,jj,nmx

    data A(1,1:5)/35.,-104.,114.,-56.,11./
    data A(2,1:5)/11.,-20.,6.,4.,-1./
    data A(3,1:5)/-1.,16.,-30.,16.,-1./
    data A(4,1:5)/-1.,4.,6.,-20.,11./
    data A(5,1:5)/11.,-56.,114.,-104.,35./
    data emfact/12./

    nmx  = num - 2
    step = sstep**2
    do j = 1, num
        k = 3
        if(j.lt.3)      k = j
        if(j.gt.nmx)    k = j - num + 5
        sum = 0.
        do i = 1,5
            jj  = j + i - k
            sum = sum + A(k,i)*f(jj)
        end do
        df(j)  = sum/(step*emfact)
    end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine deriv2                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
subroutine tderiv(f,df,num,step)                                                                  !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to calcuate the first derivative for
!--- one dimension function f, df is the derivative and num is
!--- the total point number for the arry, step is distance between
!--- two points.

   integer, intent(in) :: num
   double precision, dimension(num) :: f, df
   double precision, intent(in) :: step

   double precision, dimension(6,6) :: A
   double precision :: emfact,sum
   integer :: i,j,k,jj,nmx

   data A(1,1:6)/-274., 600.,-600., 400.,-150.,  24./
   data A(2,1:6)/ -24.,-130., 240.,-120.,  40.,  -6./
   data A(3,1:6)/   6., -60., -40., 120., -30.,   4./
   data A(4,1:6)/  -4.,  30.,-120.,  40.,  60.,  -6./
   data A(5,1:6)/   6., -40., 120.,-240., 130.,  24./
   data A(6,1:6)/ -24., 150.,-400., 600.,-600., 274./
   data emfact/120./

   nmx = num - 3
   do j = 1, num
      k = 4
      if(j.lt.4)     k = j
      if(j.gt.nmx)   k = j - num + 6
      sum = 0.
      do i = 1,6
         jj  = j + i - k
         sum = sum + A(k,i)*f(jj)
      end do
      df(j)  = sum/(step*emfact)
   end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine tderiv                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
subroutine tderiv2(f,df,num,sstep)                                                                !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to calcuate the second derivative for
!--- one dimension function f, df is the derivative and num is
!--- the total point number for the arry, step is distance between
!--- two points.

   integer, intent(in) :: num
   double precision, dimension(num) :: f, df
   double precision, intent(in) :: sstep

   double precision, dimension(6,6) :: A
   double precision :: emfact,sum,step
   integer :: i,j,k,jj,nmx

   data A(1,1:6)/ 225.,-770.,1070.,-780., 305., -50./
   data A(2,1:6)/  50., -75., -20.,  70., -30.,   5./
   data A(3,1:6)/  -5.,  80.,-150.,  80.,  -5.,   0./
   data A(4,1:6)/   0.,  -5.,  80.,-150.,  80.,  -5./
   data A(5,1:6)/   5., -30.,  70., -20., -75.,  50./
   data A(6,1:6)/ -50., 305.,-780.,1070.,-770., 225./
   data emfact/60./

   nmx  = num - 3
   step = sstep**2
   do j = 1, num
      k = 4
      if(j.lt.4)     k = j
      if(j.gt.nmx)   k = j - num + 6
      sum = 0.
      do i = 1,6
         jj  = j + i - k
         sum = sum + A(k,i)*f(jj)
      end do
      df(j)  = sum/(step*emfact)
   end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine tderiv2                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
function rtbrent(func,x1,x2,y1,y2,tol)                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!     Using Brent's method find the root of the function FUNC(X) 
!     known to lie between X1 and X2.
!     Y1 = FUNC(X1), Y2 = FUNC(X2) 
!     The root will be returned as RTBRENT with accuracy TOL
!
!     from: NUMERICAL RECIPIES, 9.3
!
!----------------------------------------------------------------------c 
    double precision :: rtbrent,func
    double precision, intent(in) :: x1,x2,y1,y2
    double precision, intent(in) :: tol
    external func

    double precision :: a,b,c,d,e,p,q,r,s,xm,fa,fb,fc,tol1
    integer :: iter

    a  = x1
    b  = x2
    fa = y1
    fb = y2
    if (fa*fb.gt.zero) stop ' in RTBRENT: root must be bracketed'
    fc = fb
    do  iter = 1,itmax
       if (fb*fc.gt.zero) then
!--- rename a,b,c and adjust bounding interval
         c  = a  
         fc = fa
         d  = b - a
         e  = d
       endif
       if (dabs(fc).lt.dabs(fb)) then
         a  = b
         b  = c
         c  = a
         fa = fb
         fb = fc
         fc = fa
       endif

!---  convergence check
       tol1 = 2*epss*dabs(b)+half*tol
       xm = half*(c-b)
       if (dabs(xm).le.tol1 .or. fb.eq.zero) then
          rtbrent = b
          return
       endif

       if (dabs(e).ge.tol1 .and. dabs(fa).gt.dabs(fb)) then
!---  attempt inverse quadratic interpolation
         s = fb/fa
         if (a.eq.c) then
            p = 2*xm*s
            q = one - s
         else
            q = fa/fc
            r = fb/fc 
            p = s*(2*xm*q*(q-r) - (b-a)*(r-one))
            q = (q-one)*(r-one)*(s-one)
         endif
         if (p.gt.zero) q = -q
!--- check whether in bounds
         p = abs(p)
         if (2*p.lt.dmin1(3*xm*q-dabs(tol1*q),dabs(e*q))) then
!--- accept interpolation
            e = d
            d = p/q
         else
!--- interpolation failed, use besection
            d = xm
            e = d
         endif
       else
!--- bounds decreasing too slowly, use bisection
         d = xm
         e = d
       endif
!--- move last best guess to a
       a  = b
       fa = fb
       if (abs(d).gt.tol1) then
!---  evaluate new trial root
          b = b + d
       else
          b = b + sign(tol1,xm)
       endif
       fb = func(b)
    end do
    stop ' in RTBRENT: exceeding maximum number of iterations'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end function rtbrent                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
subroutine ration(x,eps,ka,func,l)                                                                !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to solve the nonlinear equation 
!--- by the rational numerical methods.
!
    double precision :: func
    double precision, intent(inout) :: x
    double precision, intent(in) :: eps
    integer, intent(out) :: l
    integer, intent(in) :: ka

    external func
    
    double precision, dimension(10) :: A,Y
    double precision :: q,h2,z
    integer :: i, j, k, m,ll

    l   = 10
    q   = 1.0d+35
    ll  = 0
    do while(l.ne.0.and.ll.eq.0) 
        j   = 0
        do while(j.le.7)
            if(j.le.2) then
                z  = x + 0.01*j
            else
                z  = h2
            end if
            y(j+1) = func(ka,z)
            h2     = z
            if(j.eq.0) then
                A(1)  = z
            else
                m     = 0
                do i  = 1, j 
                    if(m.eq.0) then
                        if(dabs(h2-A(i))+1.0.eq.1.0) then
                            m  = m + 1
                        else
                            h2 = (Y(j+1) - y(i))/(h2 - A(i))
                        end if
                    end if
                end do
                A(j+1) = h2
                if(m.ne.0) then
                    A(j+1) = q
                    Stop 'A(j+1)'
                end if
                h2     = 0.0
                do i = j,1,-1
                    if(dabs(A(i+1)+h2)+1.0.eq.1.0) then
                        h2  = q
                        Stop 'h2'
                    else
                        h2  = -Y(i)/(A(i+1)+h2)
                    end if
                end do
                h2     = h2 + A(1)
            end if
            if(dabs(Y(j+1)).gt.eps) then
                j = j+1
                k = j
            else
                k = j+1
                j = 8
            end if
        end do
        if(dabs(Y(k)).gt.eps) then
            x  = h2
            l  = l - 1
            ll = 0
        else
            ll = 1
        end if
    end do
    x  = h2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine ration                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
subroutine simps(f,n,h,results)                                                                   !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!     This routine performs Simpson's rule integration of a
!     function defined by a table of equispaced values.
!     Parameters are:
! F  ---   Array of values of the function
! N  ---   Number of points
! H  ---   The uniform spacing between X values
! RESULTS  ---   Estimate of the integral that is returned to caller.
!
!----------------------------------------------------------------------c
!
    integer, intent(in) :: n
    double precision, dimension(n) ,intent(in) ::  f
    double precision,intent(out) ::  results
    double precision, intent(in) :: h


    double precision, parameter ::c3d8=0.375d0
    integer :: npanel,nhalf, nbegin,i,nend
    double precision :: x

!--- Check to see if number of panels is even.  Number of panels is N-1.
    npanel  = n - 1
    nhalf   = npanel/2
    nbegin  = 1
    results = zero

!--- Number of panels is odd.  Use 3/8 rule on first three panels, 1/3 rule on rest of them.
!
    if ((npanel-2*nhalf).ne.0) then
        results = h*c3d8*( f(1) + 3*(f(2)+f(3)) + f(4) )
        if (n.eq.4) return
        nbegin=4
    endif

!--- Apply 1/3 rule - add in first, second, last values
!
    results = results + h*third*( f(nbegin) + 4*f(nbegin+1) + f(n) )
    nbegin  = nbegin+2
    if (nbegin.eq.n) then
        return
    else
        x = zero
        nend = n-2
        do i = nbegin,nend,2
            x = x + f(i) + 2*f(i+1) 
        enddo
        results = results + h*two*third*x

        return
    endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine simps                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine integ2(fun,npt,xstep,match,results,flo,fhi)                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!     This subroutine integrates the function FUN which has a
!     discontinuity at MATCH.
!
!----------------------------------------------------------------------c
    integer, intent(in) :: match,npt
    double precision, dimension(match), intent(out) :: fun
    double precision, intent(out) :: results
    double precision, intent(in) :: xstep,flo,fhi

    double precision :: fmatch,res1,res2
    fmatch     = fun(match)

!---- Integrate up until match point
    fun(match) = flo
    call simps(fun,match,xstep,res1)

!---- Integrate from the match point to the end
    fun(match) = fhi
    call simps(fun(match),npt-match+1,xstep,res2)

!---- Add the resulting integrals
    results    = res1+res2

    fun(match) = fmatch

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine integ2                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
function DET(ma, n, aa)                                                                           !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the determinant
!
    double precision :: DET
    integer, intent(in) :: ma,n
    double precision, dimension(ma,n), intent(in) :: aa
    double precision, dimension(ma,n) :: a

    double precision, parameter :: tol = 1.0d-20
    integer :: k, i, j,l
    double precision :: s, p, q,cq,cp

    a   = aa
!--- Normalization of the coulomn of the matrix a
    do k = 1,n
        s = zero
        do i = 1, n
            s      = s + a(i,k)**2
        end do
        s = one/dsqrt(s)
        do i = 1, n
            a(i,k) = s*a(i,k)
        end do
    end do

!    
    DET    = one
    do k = 1, n
        p  = zero
        do j = k, n
            q  = dabs(a(j,k))
            if(q.gt.p) then
                p  = q
                i  = j
            end if
        end do
        if(p.lt.tol) then
            DET  = 0.0d0
            return
        end if
        if(i.ne.k) then
            DET  = -DET
            do l = k,n
                cq     = a(i,l)
                a(i,l) = a(k,l)
                a(k,l) = cq
            end do
        end if
        DET  = DET*a(k,k)
        cp   = one/a(k,k)
        do i = k+1, n
            cq  = a(i,k)*cp
            do l = k+1, n
                a(i,l) = a(i,l) - cq*a(k,l)
            end do
        end do
    end do
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end function DET                                                                                  !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
subroutine sdiag(nmax,n,a,d,x,e,is)                                                               !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!=======================================================================
!
!     A   matrix to be diagonalized
!     D   eigenvalues    
!     X   eigenvectors
!     E   auxiliary field
!     IS = 1  eigenvalues are ordered and major component of X is positiv
!          0  eigenvalues are not ordered            
!-----------------------------------------------------------------------
    integer, intent(in) :: nmax, n,is
    double precision, dimension(nmax,nmax) :: a,x
    double precision, dimension(n) :: e,d

    double precision, parameter ::  tol = 1.d-32,eps=1.d-10    
    integer :: i, j,l,j1,k, im
    double precision :: f,g,h,hi,s,p,b,pr,r,c
                                  
!
!    write(*,*) tol, eps
      if (n.eq.1) then
         d(1)   = a(1,1)  
         x(1,1) = one
         return
      endif
!
      do i = 1,n 
         do j = 1,i 
            x(i,j)=a(i,j)
         enddo
      enddo
!
!cc   Householder-reduction
      do i = n,2,-1
         l = i - 2
         f = x(i,i-1)
         g = f            
         h = zero
         do k = 1,l
            h = h + x(i,k)*x(i,k)
         enddo
         s = h + f*f
         if (s.lt.tol) h = zero
         if (h.gt.zero) then         
            l = l+1                        
            g = dsqrt(s)
            if (f.ge.zero) g = -g        
            h = s - f*g                 
            hi = one/h                
            x(i,i-1) = f - g          
            f = zero                
            do j = 1,l        
               x(j,i) = x(i,j)*hi  
               s = zero             
               do k = 1,j     
                  s = s + x(j,k)*x(i,k)                      
               enddo
               j1 = j+1                                
               do k = j1,l                        
                  s = s + x(k,j)*x(i,k)                  
               enddo
               e(j) = s*hi                         
               f = f + s*x(j,i)                     
            enddo
            f = f*hi*half
            do j = 1,l                  
               s    = x(i,j)                    
               e(j) = e(j) - f*s              
               p    = e(j)                    
               do k = 1,j              
                  x(j,k) = x(j,k) - s*e(k) - x(i,k)*p        
               enddo
            enddo
         endif
         d(i)   = h                            
         e(i-1) = g                         
      enddo
!            
!cc   transformation matrix 
      d(1) = zero                               
      e(n) = zero                             
      b    = zero                               
      f    = zero                              
      do i = 1,n                      
         l = i-1                            
         if (d(i).ne.zero) then        
            do j = 1,l                  
               s = zero                         
               do k = 1,l                
                  s = s + x(i,k)*x(k,j)          
               enddo
               do k = 1,l              
                  x(k,j) = x(k,j) - s*x(k,i)   
               enddo
            enddo
         endif
         d(i)   = x(i,i)            
         x(i,i) = one              
         do j = 1,l        
            x(i,j) = zero         
            x(j,i) = zero         
         enddo
      enddo
!
!cc   diagonalizition of tri-diagonal-matrix
      do l = 1,n                     
         h = eps*( abs(d(l))+ abs(e(l)))
         if (h.gt.b) b = h             
!
!cc      test for splitting        
         do j = l,n              
            if ( abs(e(j)).le.b) goto 10
         enddo
!
!cc      test for convergence    
   10    if (j.gt.l) then   
   20       p  = (d(l+1)-d(l))/(2*e(l))          
            r  = dsqrt(p*p+1.d0)
            pr = p + r                           
            if (p.lt.zero) pr = p - r             
            h = d(l) - e(l)/pr                 
            do i = l,n                  
               d(i) = d(i) - h
            enddo                  
            f = f + h                       
!
!cc         QR-transformation          
            p = d(j)                    
            c = one                     
            s = zero                    
            do i = j-1,l,-1 
               g = c*e(i)            
               h = c*p              
               if ( abs(p).lt.abs(e(i))) then
                  c = p/e(i)                     
                  r = dsqrt(c*c+one)
                  e(i+1) = s*e(i)*r             
                  s = one/r                      
                  c = c/r                     
               else
                  c = e(i)/p                          
                  r = dsqrt(c*c+one)
                  e(i+1) = s*p*r                     
                  s = c/r                           
                  c = one/r                         
               endif
               p = c*d(i) - s*g             
               d(i+1) = h + s*(c*g+s*d(i)) 
               do k = 1,n           
                  h        = x(k,i+1)            
                  x(k,i+1) = x(k,i)*s + h*c
                  x(k,i)   = x(k,i)*c - h*s 
               enddo
            enddo          
            e(l) = s*p          
            d(l) = c*p         
            if ( abs(e(l)).gt.b) goto 20
         endif
!
!cc      convergence      
         d(l) = d(l) + f    
      enddo
!
      if (is.eq.0) return
!cc   ordering of eigenvalues    
      do i = 1,n            
         k  = i                    
         p  = d(i)                
         do j = i+1,n          
            if (d(j).lt.p) then 
                k = j                    
                p = d(j)
            endif              
         enddo             
         if (k.ne.i) then
            d(k) = d(i)          
            d(i) = p            
            do j = 1,n     
               p      = x(j,i)        
               x(j,i) = x(j,k)  
               x(j,k) = p
            enddo
         endif
      enddo      
!                 
!     signum
      do  k = 1,n
         s = zero
         do i = 1,n
            h = abs(x(i,k))
            if (h.gt.s) then
               s  = h
               im = i
            endif
         enddo
         if (x(im,k).lt.zero) then
            do i = 1,n
               x(i,k) = -x(i,k)
            enddo
         endif
      enddo
! 
      return
!-end-SDIAG
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine sdiag                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
FUNCTION zriddr(func,x1,x2, f1, f2, xacc, lg, lf, gf)                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
    double precision :: zriddr, func
    double precision, intent(in) :: x1, x2, xacc, f1, f2
    integer, intent(in) :: lg, lf
    double precision, dimension(2,2) :: gf
    EXTERNAL func

    INTEGER, parameter :: MAXIT = 60
    double precision, parameter :: UNUSED = -1.11d30
    
    INTEGER :: j
    double precision :: fh,fl,fm,fnew,s,xh,xl,xm,xnew
    fl   =  f1
    fh   =  f2
    if((fl.gt.zero.and.fh.lt.zero).or.(fl.lt.zero.and.fh.gt.zero))then
        xl     =  x1
        xh     =  x2
        zriddr =  UNUSED
        do j=  1,MAXIT
            xm     =  0.5*(xl+xh)
            fm     =  func(lg, lf, xm, gf)
            s      =  dsqrt(fm**2-fl*fh)
            if(s.eq.0.)return
            xnew   =  xm+(xm-xl)*(sign(1.d0,fl-fh)*fm/s)
            if (dabs(xnew-zriddr).le.xacc) return
            zriddr =  xnew
            fnew   =  func(lg, lf, zriddr, gf)
            if (fnew.eq.0.) return
            if(dsign(fm,fnew).ne.fm) then
                xl   =  xm
                fl   =  fm
                xh   =  zriddr
                fh   =  fnew
            else if(dsign(fl,fnew).ne.fl) then
                xh   =  zriddr
                fh   =  fnew
            else if(dsign(fh,fnew).ne.fh) then
                xl   =  zriddr
                fl   =  fnew
            else
                Stop 'never get here in zriddr'
            endif
            if(dabs(xh-xl).le.xacc) return
        end do
        Stop 'zriddr exceed maximum iterations'
    else if (fl.eq.zero) then
        zriddr =  x1
    else if (fh.eq.zero) then
        zriddr =  x2
    else
        Stop 'root must be bracketed in zriddr'
    endif
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Function zriddr                                                                               !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine intpol6(v, vh, n)                                                                      !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!---- Six point Lagrange interpolation
!
    integer, intent(in) :: n
    double precision, dimension(n), intent(in)  :: v
    double precision, dimension(n), intent(out) :: vh

    double precision, dimension(5, 6) :: A
    integer :: i, j

    data A(1, 1:6) / 0.24609375d0,   1.23046875d0,  -0.8203125d0,   0.4921875d0,  -0.17578125d0,   0.02734375d0/
    data A(2, 1:6) /-0.02734375d0,   0.41015625d0,   0.8203125d0,  -0.2734375d0,   0.08203125d0,  -0.01171875d0/
    data A(3, 1:6) / 0.01171875d0,  -0.09765625d0,   0.5859375d0,   0.5859375d0,  -0.09765625d0,   0.01171875d0/
    data A(4, 1:6) /-0.01171875d0,   0.08203125d0,  -0.2734375d0,   0.8203125d0,   0.41015625d0,  -0.02734375d0/
    data A(5, 1:6) / 0.02734375d0,  -0.17578125d0,   0.4921875d0,  -0.8203125d0,   1.23046875d0,   0.24609375d0/

    do i = 1, 2
        vh(i)   = v(1)*A(i, 1) + v(2)*A(i, 2) + v(3)*A(i, 3) + v(4)*A(i, 4) + v(5)*A(i, 5) + v(6)*A(i, 6)
    end do

    do i = 3, n-3
        vh(i)   = v(i-2)*A(3,1) + v(i-1)*A(3,2) + v(i)*A(3,3) + v(i+1)*A(3,4) + v(i+2)*A(3,5) + v(i+3)*A(3,6)
    end do

    do i = n-2, n-1
        j       = i + 6 - n
        vh(i)   = v(n-5)*A(j,1) + v(n-4)*A(j,2) + v(n-3)*A(j,3) + v(n-2)*A(j,4) + v(n-1)*A(j,5) + v(n)*A(j,6)
    end do

    vh(n)       = zero
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine intpol6                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine intpolm(v,vh,n)                                                                        !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!     Interpolates potentials between mesh points  VH(i) = V(i+1/2)
!     Lagrange 4-point interpolation
!
!----------------------------------------------------------------------c
    integer, intent(in) :: n
    integer :: ie,i
    double precision , dimension(n) :: v,vh
    double precision, parameter :: c1=0.0625d0,c5=0.3125d0,c9=0.5625d0,c15=0.9375d0

    ie = n-2
    do i = 2,ie
       vh(i)  = (v(i)+v(i+1))*c9 - (v(i-1)+v(i+2))*c1
    enddo
    vh(1)   = (v(1)-v(3)  )*c5 + v(2)  *c15 + v(4)  *c1
    vh(n-1) = (v(n)-v(n-2))*c5 + v(n-1)*c15 + v(n-3)*c1
    vh(n)   = 0.d0
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine intpolm                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine intpolp(v,vh,n)                                                                        !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!     Interpolates potentials between mesh points  VH(i) = V(i-1/2)
!     Lagrange 4-point interpolation
!
!----------------------------------------------------------------------c
    integer, intent(in) :: n
    integer :: ie,i
    double precision , dimension(n) :: v,vh
    double precision, parameter :: c1=0.0625d0,c5=0.3125d0,c9=0.5625d0,c15=0.9375d0

    ie = n-1
    do i = 3,ie
       vh(i)  = (v(i)+v(i-1))*c9 - (v(i+1)+v(i-2))*c1
    enddo
    vh(2)   = (v(1)-v(3)  )*c5 + v(2)  *c15 + v(4)  *c1
    vh(n)   = (v(n)-v(n-2))*c5 + v(n-1)*c15 + v(n-3)*c1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine intpolp                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
SUBROUTINE ratint(xa,ya,n,x,y,dy)                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Given arrays xa, ya, each of length n, and given a value of x, this routine
!--- returns a value of y and an accuracy dy. The value returned is that of the 
!--- diagonal rational function, evaluated at x, which passes through the n points
!--- (xa_i, ya_i), i = 1, ..., n
!
    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: xa, ya
    double precision, intent(in) :: x
    double precision, intent(out) :: dy, y
    integer, parameter :: NMAX=10
    double precision, parameter :: TINY=1.e-25

    integer :: i,m,ns
    double precision :: dd,h,hh,t,w,c(NMAX),d(NMAX)
    
    ns  =  1
    hh  =  dabs(x-xa(1))
    do i =  1,n
        h       =  dabs(x-xa(i))
        if (h.eq.0.)then
            y   =  ya(i)
            dy  =  0.0
            return
        else if (h.lt.hh) then
            ns  =  i
            hh  =  h
        endif
        c(i)    =  ya(i)
        d(i)    =  ya(i)+TINY
    end do
    y   =  ya(ns)
    ns  =  ns-1
    do m =  1,n-1
        do i =  1,n-m
            w   =  c(i+1)-d(i)
            h   =  xa(i+m)-x
            t   =  (xa(i)-x)*d(i)/h
            dd  =  t-c(i+1)
            if(dd.eq.0.)    Stop 'failure in ratint'
            dd  =  w/dd
            d(i)=  c(i+1)*dd
            c(i)=  t*dd
        end do
        if (2*ns.lt.n-m)then
            dy  =  c(ns+1)
        else
            dy  =  d(ns)
            ns  =  ns-1
        endif
        y   =  y+dy
    end do
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
END subroutine ratint                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
SUBROUTINE polint(xa,ya,n,x,y,dy)                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: xa, ya
    double precision, intent(in) :: x
    double precision, intent(out) :: dy, y
    integer, parameter :: NMAX=10

    integer :: i,m,ns
    double precision :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns  =  1
    dif =  abs(x-xa(1))
    do i =  1,n
        dift    =  abs(x-xa(i))
        if (dift.lt.dif) then
            ns  =  i
            dif =  dift
        endif
        c(i)    =  ya(i)
        d(i)    =  ya(i)
    end do
    y   =  ya(ns)
    ns  =  ns-1
    do  m =  1,n-1
        do  i=  1,n-m
            ho  =  xa(i)-x
            hp  =  xa(i+m)-x
            w   =  c(i+1)-d(i)
            den =  ho-hp
            if(den.eq.0.)Stop 'failure in polint'
            den     =  w/den
            d(i)    =  hp*den
            c(i)    =  ho*den
        end do
        if (2*ns.lt.n-m)then
            dy  =  c(ns+1)
        else
            dy  =  d(ns)
            ns  =  ns-1
        endif
        y   =  y+dy
    end do
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
END subroutine polint                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!                                                                                                 !
end module RHFlib                                                                                 !
!                                                                                                 !
!*************************************************************************************************!
