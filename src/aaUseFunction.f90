
! ******************** BGOS ********************************

    subroutine bgos(sx,id,x1,x2,ro)

! ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)
    implicit none

    integer::id
    double precision::ro,sx,f,v1,v2,s,dls,ro2,x1,x2,uniran

5     continue

    x1=uniran()
    x2=uniran()

    if (id.ne.1) goto 10
    f=2.*sqrt(3.)
    x1=(x1-0.5)*f
    x2=(x2-0.5)*f
    goto 20
10    continue
    v1=2.*x1-1
    v2=2.*x2-1
    s=v1*v1+v2*v2
    if (s.ge.1.) goto 5
    dls=sqrt(-2.*log(s)/s)
    x1=v1*dls
    x2=v2*dls
20    continue
    ro2=ro*ro
    if (abs(ro).gt.1.E-10) x2=(x1+x2*sqrt(1./ro2-1.))*ro
    x1=x1*sx
    x2=x2*sx

    return

    end subroutine bgos


!------------------- FIN SUBROUTINE BGOS -----------------

! ------------------------------------------------------

    double precision function uniran()
!
!     Random number generator(RCARRY), adapted from F. James
!     "A Review of Random Number Generators"
!      Comp. Phys. Comm. 60(1990), pp. 329-344.
!
    double precision,save::carry
    double precision,dimension(24),save::seeds
    double precision,parameter::one=1
    double precision,parameter::twom24 = ONE/16777216
    integer,save::i,j
    data i, j, carry / 24, 10, 0.0 /
    data seeds / &
    0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635, &
    0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
    0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
    0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/

    uniran = seeds(i) - seeds(j) - carry

    if (uniran .lt. 0) then
        uniran = uniran + 1
        carry = twom24
    else
        carry = 0
    end if

    seeds(I) = uniran
    I = 24 - MOD( 25-I, 24 )
    J = 24 - MOD( 25-J, 24 )

    end function uniran



    subroutine percentile3(t,n,p,out)

    implicit none

    integer::n,ib,i,indd
    double precision::a,b,c,temp,p
    double precision,dimension(n)::t
    double precision,intent(out)::out

    n=size(t)

    ! tri des temps
    indd=1
    do while (indd.eq.1)
        indd=0
        do i=1,(n-1)
            if (t(i).gt.t(i+1)) then
                temp=t(i)
                t(i)=t(i+1)
                t(i+1)=temp
                indd=1
            end if
        end do
    end do

    ! quantile d'ordre p%
    a=(n-1)*p
    ! b=mod(a,1.0d0) ! scl: 23/11/2018
    b = a - int(a) ! scl: 23/11/2018
    c=a-b
    ib=int(c)
    
    !out= (1-b)*t(ib+1)+b*t(ib+2) ! scl: 23/11/2018
    if(ib <= n-2)then ! scl: 23/11/2018
    out = (1-b)*t(ib+1)+b*t(ib+2)
    else
    out = t(n) ! l'on suppose ici qu'on cherche le 100th percentile de la serie, ce qui est normale car dans ce cas, q est tres proche de 1
    endif

    end subroutine percentile3


    subroutine percentile2(t,n,t25,t975)

    implicit none
    ! Based on discussion in this paper http://www.haiweb.org/medicineprices/manual/quartiles_iTSS.pdf

    integer::n,ib,i,indd
    double precision::a,b,c,temp
    double precision,dimension(n)::t
    double precision,intent(out)::t25,t975

    n=size(t)

    ! tri des temps
    indd=1
    do while (indd.eq.1)
        indd=0
        do i=1,(n-1)
            if (t(i).gt.t(i+1)) then
                temp=t(i)
                t(i)=t(i+1)
                t(i+1)=temp
                indd=1
            end if
        end do
    end do

    ! quantile d'ordre 2.5%
    a=(n-1)*0.025d0
    ! b=mod(a,1.0d0) ! scl: 23/11/2018
    b = a - int(a) ! scl: 23/11/2018
    c=a-b
    ib=int(c)
    t25= (1-b)*t(ib+1)+b*t(ib+2)

    ! quantile d'ordre 97.5%
    a=(n-1)*0.975d0
    ! b=mod(a,1.0d0) ! scl: 23/11/2018
    b = a - int(a) ! scl: 23/11/2018
    c=a-b
    ib=int(c)
    t975= (1-b)*t(ib+1)+b*t(ib+2)

    end subroutine percentile2


    subroutine percentile(t,t25,t975) ! pour le MCMC : version qui ne marche qu'avec 1000 echantillons

    integer::indd,i
    double precision::t25,t975,temp
    double precision,dimension(1000)::t

    indd=1
    do while (indd.eq.1)
        indd=0
        do i=1,999
            if (t(i).gt.t(i+1)) then
                temp=t(i)
                t(i)=t(i+1)
                t(i+1)=temp
                indd=1
            end if
        end do
    end do
    
        !t25=t(25) ! scl: 23/11/2018
        t25 = 0.25d0*t(250)+0.75d0*t(251) ! scl: 23/11/2018
        !t975=t(975) ! scl: 23/11/2018
        t975 = 0.975d0*t(975)+0.025d0*t(976) ! scl: 23/11/2018

    end subroutine percentile



    subroutine bb(i1, i2, i3, y, newknots, result, dumsub)
    ! This subroutine calculates i2 th basis of spline of
    ! degree (i3-1).
    IMPLICIT NONE
    integer(kind=4) i1, i2, i3
    double precision y, newknots(i1), temp1, temp2, result, result1, result2
    external dumsub

    if(i3.eq.1) then ! ordre des B-splines egal a 1
        if((y.ge.newknots(i2)).and.(y.lt.newknots(i2+1))) then
            result=1.d0
        else
            result=0.d0
        endif
    else ! ordre de B-splines superieur a 1
        call dumsub(i1, i2, (i3-1), y, newknots, result1, dumsub)
        temp1=(y-newknots(i2))*result1/(newknots(i2+i3-1)-newknots(i2))

        if(temp1.ne.temp1) temp1=0.d0

        call dumsub(i1, (i2+1), (i3-1), y, newknots, result2, dumsub)
        temp2=(newknots(i2+i3)-y)*result2/(newknots(i2+i3)-newknots(i2+1))

        if(temp2.ne.temp2) temp2=0.d0

        result=temp1+temp2
    endif

    return

    end
    !If one wants to call this subroutine from R following is an example code.
    !We assume that the above Fortran subroutines were saved in a file named �spline.f�.




    subroutine splinebasisIndiv(d, m, m1, k, x, innerknots,boundaryknots, basis)
! This subroutine generates Bspline basis functions.
! x(n) is a n by 1 input vector for which B-spline basis
! function will be evaluated.
! innerknots(m1) set of m1 innerknot points.
! newknots is the entire set of knots, of length m=m1+2(d+1)
! where d is the degree of the splines.
! k=number of spline basis=m1+d+1
    IMPLICIT NONE
    integer(kind=4) d, k, m, m1
    double precision x, innerknots(m1), boundaryknots(2)
    double precision newknots(m), basis(k), result
    external bb
    integer(kind=4) i1, j

    do i1=1, (d+1)
        newknots(i1)=boundaryknots(1)
    end do

    do i1=(d+2), (m1+d+1)
        newknots(i1)=innerknots(i1-d-1)
    end do

    do i1=(m1+d+2), m
        newknots(i1)=boundaryknots(2)
    end do

!    do i=1, n
        if(x.eq.boundaryknots(2)) then
            basis(k)=1.d0
            do j=1,(k-1)
                basis(j)=0.d0
            end do
        else
            do j=1,k
                call bb(m, j, (d+1), x, newknots, result, bb)
                basis(j)=result
            end do
        endif
!    end do

    return
    end

!====================================================================
!====================================================================
!     double precision function Basissplines(j,q,t,zii,m)
!
!     integer::j,q,m
!     double precision::t,bbb
!     double precision,dimension(-q+1:m+q)::zii
!     bbb=0.d0
!
!     if (q.eq.3) then
!         bbb=0.d0
!         if (t.ge.zii(j).and.t.lt.zii(j+1))then
!             bbb =((t-zii(j))**2.d0)/((zii(j+2)-zii(j))*(zii(j+1)-zii(j)))!ok
!         endif
!         if (t.ge.zii(j+1).and.t.lt.zii(j+2))then
!             bbb =((t-zii(j))*(zii(j+2)-t))/((zii(j+2)-zii(j))*(zii(j+2)-zii(j+1))) &
!             + ((zii(j+3)-t)*(t-zii(j+1)))/((zii(j+3)-zii(j+1))*(zii(j+2)-zii(j+1)))
!         endif
!         if (t.ge.zii(j+2).and.t.lt.zii(j+3))then
!             bbb=((zii(j+3)-t)**2.d0)/((zii(j+3)-zii(j+1))*(zii(j+3)-zii(j+2)))
!         endif
!          if (t.eq.zii(m+1).and.j.eq.m)then
! ! !            bbb=((zii(m+q)-t)**2.d0)/((zii(m+q)-zii(m+q-2))*(zii(m+q)-zii(m+q-1)))
! ! !            bbb=((t-zii(m))**2.d0)/((zii(m+2)-zii(m))*(zii(m+1)-zii(m)))
!              bbb=1.d0
!         else
!              bbb=0.d0
!          endif
!     endif
!
!     if (q.eq.2) then
!         bbb=0.d0
!         if (t.ge.zii(j).and.t.lt.zii(j+1))then
!             bbb =(t-zii(j))/(zii(j+1)-zii(j))!ok
!         endif
!         if (t.ge.zii(j+1).and.t.lt.zii(j+2))then
!             bbb =(zii(j+2)-t)/(zii(j+2)-zii(j+1))
!         endif
!          if (t.eq.zii(m+1))then
! !              t=t-(t/1.d5)
! !              bbb =(t-zii(m))/(zii(m+1)-zii(m))!ok
!              bbb=1.d0
!          endif
!     endif
!
!     if (q.eq.1) then
!         bbb=0.d0
!          if (t.ge.zii(j).and.t.lt.zii(j+1))then
!                  bbb =1.d0
!          endif
!           if (t.eq.zii(m+1))then
!
!              bbb=1.d0
!              bbb=1.d0
!           endif
!     endif
!
!     Basissplines=bbb
!     return
!
!     end function Basissplines

!===============================================================================================================
!===============================================================================================================
! quadratic B-splines

!     innerknots(1:nbinnerknots)=knotsTPS(1:nbinnerknots)
!     boundaryknots(1)=knotsTPS(0)
!     boundaryknots(2)=knotsTPS(nbinnerknots+1)
! !d=degree , n=nsujet, m=m1+2*(d+1), m1=nb noeuds internes, k=m1+d+1 (nb de bases de splines), x: vecteur tps, innerknots: noeuds interne, boundaryknots: bornes de l'intervalle de tps de suivi - 0 et cens, basis: base pour chaque tps
!     call splinebasis(qorder-1, nsujet,nbinnerknots+2*qorder, nbinnerknots, nbinnerknots+qorder, t1, innerknots,boundaryknots, basis)



    !Following is the Fortran code to generate B-spline basis function.
    subroutine splinebasis(d, n, m, m1, k, x, innerknots,boundaryknots, basis)
    ! This subroutine generates Bspline basis functions.
    ! x(n) is a n by 1 input vector for which B-spline basis
    ! function will be evaluated.
    ! innerknots(m1) set of m1 innerknot points.
    ! newknots is the entire set of knots, of length m=m1+2(d+1)
    ! where d is the degree of the splines.
    ! k=number of spline basis=m1+d+1
    IMPLICIT NONE
    integer(kind=4) d, k, m, m1, n
    double precision x(n), innerknots(m1), boundaryknots(2)
    double precision newknots(m), basis(n, k), result
    external bb
    integer(kind=4) i1, i, j

    do i1=1, (d+1)
        newknots(i1)=boundaryknots(1)
    end do

    do i1=(d+2), (m1+d+1)
        newknots(i1)=innerknots(i1-d-1)
    end do

    do i1=(m1+d+2), m
        newknots(i1)=boundaryknots(2)
    end do

    do i=1,n
        if(x(i).eq.boundaryknots(2)) then
            basis(i, k)=1.d0
            do j=1, (k-1)
                basis(i, j)=0.d0
            end do
        else
            do j=1, k
                call bb(m, j, (d+1), x(i), newknots, result, bb)
                basis(i, j)=result
            end do
        endif
    end do

    return
    end


!====================================================================!
!=============  SUBROUTINE D'INTEGRATION SUR A ET B =================!
!====================================================================!


       subroutine integration(f,a,b,result,abserr,resabs,resasc,i,bh,np)
!***begin prologue  dqk15
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           double precision version
!
!           parameters
!            on entry
!              f      - double precision
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - double precision
!                       lower limit of integration
!
!              b      - double precision
!                       upper limit of integration
!
!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the7-point gauss rule(resg).
!
!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - double precision
!                       approximation to the integral j
!
!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  dqk15
!
       double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,d1mach(5),epmach,f,fc,fsum,fval1,fval2,&
       fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk
       integer :: j,jtw,jtwm1
       external :: f
       integer ::i,np
       double precision :: bh(np)

       dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!          wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    wg(1)=0.129484966168869693270611432679082d0
    wg(2)=0.279705391489276667901467771423780d0
    wg(3)=0.381830050505118944950369775488975d0
    wg(4)=0.417959183673469387755102040816327d0

    xgk(1)=0.991455371120812639206854697526329d0
    xgk(2)=0.949107912342758524526189684047851d0
    xgk(3)=0.864864423359769072789712788640926d0
    xgk(4)=0.741531185599394439863864773280788d0
    xgk(5)=0.586087235467691130294144838258730d0
    xgk(6)=0.405845151377397166906606412076961d0
    xgk(7)=0.207784955007898467600689403773245d0
    xgk(8)=0.000000000000000000000000000000000d0

    wgk(1)=0.022935322010529224963732008058970d0
    wgk(2)=0.063092092629978553290700663189204d0
    wgk(3)=0.104790010322250183839876322541518d0
    wgk(4)=0.140653259715525918745189590510238d0
    wgk(5)=0.169004726639267902826583426598550d0
    wgk(6)=0.190350578064785409913256402421014d0
    wgk(7)=0.204432940075298892414161999234649d0
    wgk(8)=0.209482141084727828012999174891714d0


!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!
!***first executable statement  dqk15

    D1MACH(1)=2.23D-308
    D1MACH(2)=1.79D+308
    D1MACH(3)=1.11D-16
    D1MACH(4)=2.22D-16
    D1MACH(5)=0.301029995663981195D0

    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.


    fc = f(centr,i,bh,np)

    resg = fc*wg(4)
    resk = fc*wgk(8)
    resabs = dabs(resk)

    do j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,i,bh,np)
        fval2 = f(centr+absc,i,bh,np)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    enddo

    do j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,i,bh,np)
        fval2 = f(centr+absc,i,bh,np)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    enddo

    reskh = resk*0.5d+00
    resasc = wgk(8)*dabs(fc-reskh)
    do j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    enddo

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

    return

    end


!====================================================================!
!=============  SUBROUTINE D'INTEGRATION SUR A ET B =================!
!====================================================================!

       subroutine integrationdc(f,a,b,result,abserr,resabs,resasc,i,bh,np,frail)
!***begin prologue  dqk15
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           double precision version
!
!           parameters
!            on entry
!              f      - double precision
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - double precision
!                       lower limit of integration
!
!              b      - double precision
!                       upper limit of integration
!
!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the7-point gauss rule(resg).
!
!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - double precision
!                       approximation to the integral j
!
!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  dqk15
    use comon,only:nea

       double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,d1mach(5),epmach,f,fc,fsum,fval1,fval2,&
       fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,uflow,wg,wgk,xgk
       double precision,dimension(nea)::frail
       integer :: j,jtw,jtwm1
       external :: f
       integer ::i,np
       double precision :: bh(np)

       dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!          wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    wg(1)=0.129484966168869693270611432679082d0
    wg(2)=0.279705391489276667901467771423780d0
    wg(3)=0.381830050505118944950369775488975d0
    wg(4)=0.417959183673469387755102040816327d0

    xgk(1)=0.991455371120812639206854697526329d0
    xgk(2)=0.949107912342758524526189684047851d0
    xgk(3)=0.864864423359769072789712788640926d0
    xgk(4)=0.741531185599394439863864773280788d0
    xgk(5)=0.586087235467691130294144838258730d0
    xgk(6)=0.405845151377397166906606412076961d0
    xgk(7)=0.207784955007898467600689403773245d0
    xgk(8)=0.000000000000000000000000000000000d0

    wgk(1)=0.022935322010529224963732008058970d0
    wgk(2)=0.063092092629978553290700663189204d0
    wgk(3)=0.104790010322250183839876322541518d0
    wgk(4)=0.140653259715525918745189590510238d0
    wgk(5)=0.169004726639267902826583426598550d0
    wgk(6)=0.190350578064785409913256402421014d0
    wgk(7)=0.204432940075298892414161999234649d0
    wgk(8)=0.209482141084727828012999174891714d0


!           list of major variables
!           -----------------------
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
!           machine dependent constants
!           ---------------------------
!
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!
!***first executable statement  dqk15

    D1MACH(1)=2.23D-308
    D1MACH(2)=1.79D+308
    D1MACH(3)=1.11D-16
    D1MACH(4)=2.22D-16
    D1MACH(5)=0.301029995663981195D0

    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.


    fc = f(centr,i,bh,np,frail)

    resg = fc*wg(4)
    resk = fc*wgk(8)
    resabs = dabs(resk)

    do j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)

        fval1 = f(centr-absc,i,bh,np,frail)

        fval2 = f(centr+absc,i,bh,np,frail)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    enddo

    do j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,i,bh,np,frail)
        fval2 = f(centr+absc,i,bh,np,frail)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    enddo

    reskh = resk*0.5d+00
    resasc = wgk(8)*dabs(fc-reskh)
    do j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    enddo
!write(*,*)resk,hlgth,a,b
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

    return

    end


!****************************************************************
    subroutine searchknotstps(tps0,knots0,nbinnerknots0,qorder0,nsujetmax0,equidistantTPS0,c0,begin)

    use comon

    integer::ii,jj,nsujetmax0,nbrecu0,indd,entTPS,nbinnerknots0,qorder0,equidistantTPS0
    integer,dimension(nsujetmax0)::c0
    double precision::temp
    double precision,dimension(-qorder0+1:nbinnerknots0+qorder0)::knots0
    double precision,dimension(nsujetmax0)::tps0
    double precision::begin

    equidistantTPS0=0
    jj=0
    do ii=1,nsujetmax0
        if(tps0(ii).ne.(0.d0).and.c0(ii).eq.1)then
            jj=jj+1
        endif
    end do
    nbrecu0=jj

    allocate(t2(nbrecu0))

! remplissage du vecteur de temps
    jj=0

    do ii=1,nsujetmax0
        if(tps0(ii).ne.(0.d0).and.c0(ii).eq.1)then
            jj=jj+1
            t2(jj)=tps0(ii)
        endif
    end do

    indd=1
    do while (indd.eq.1)
        indd=0
        do i=1,nbrecu0-1
            if (t2(i).gt.t2(i+1)) then
                temp=t2(i)
                t2(i)=t2(i+1)
                t2(i+1)=temp
                indd=1
            end if
        end do
    end do

    entTPS=int(nbrecu0/(nbinnerknots0+1))

    j=0

! ajout beta tps 24/02/2012
    do j=1,nbinnerknots0
        !if (equidistantTPS0.eq.0) then
            knots0(j)=(t2(entTPS*j)+t2(entTPS*j+1))/(2.d0)
        !else
        !    knots0(j)=(cens/(nbinnerknots0+1))*j
        !endif
    end do
    knots0(-qorder0+1:0)=begin !0.d0
    knots0(nbinnerknots0+1:nbinnerknots0+qorder0)=cens

    deallocate(t2)

    end subroutine searchknotstps

! ==============================================

    subroutine drawTimeCoef(np,b,nvar,filtre,Out)

    use comon
    use betatttps

    implicit none

    integer::p,i,iii,ptps,jj,jjj,np
    double precision::tps,iiii
    integer::nvar
    double precision,dimension(np)::b
    integer,dimension(nvar)::filtre
    double precision,dimension(0:100,0:4*sum(filtre)),intent(out)::Out

    allocate(betatpsX(0:100),betatpsminX(0:100),betatpsmaxX(0:100),varBetatps(0:100))
    allocate(BasisSinhaTPS(nbinnerknots+qorder))

    p=1
    ptps=0
    do i=1,nvar
        if (filtre(i).eq.1) then ! beta dependant du temps
            betatpsX=0.d0
            betatpsminX=0.d0
            betatpsmaxX=0.d0
            varBetatps=0.d0
            do iii=0,100
                iiii=dble(iii)
                tps=knotsTPS(0)+(knotsTPS(nbinnerknots+1)-knotsTPS(0))*(iiii/100.d0)

                Out(iii,0) = tps

                call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,nbinnerknots+qorder, &
                tps,innerknots,boundaryknots,basisSinhaTPS)

                ! on estime les beta dependant du temps sur 100 points (arbitraire)
                do jj=-qorder+1,nbinnerknots
                    betatpsX(iii)=betatpsX(iii)+b(np-(nva+npbetatps)+(p-1)+jj+qorder)*BasisSinhaTPS(jj+qorder)
                end do

                ! on estime les variances des beta dependant du temps
                do jj=-qorder+1,nbinnerknots
                    do jjj=-qorder+1,nbinnerknots
                        varBetatps(iii)=varBetatps(iii)+BasisSinhaTPS(jj+qorder)*BasisSinhaTPS(jjj+qorder) &
                        *(H_hess(np-(nva+npbetatps)+(p-1)+jj+qorder,np-(nva+npbetatps)+(p-1)+jjj+qorder))
                    end do
                end do

                ! bornes de l'intervalle de confiance
                betatpsminX(iii)=betatpsX(iii)-1.96*dsqrt(varBetatps(iii))
                betatpsmaxX(iii)=betatpsX(iii)+1.96*dsqrt(varBetatps(iii))

                ! matrice avec estimations/ICmin/ICmax/ecarts-type
                Out(iii,1+ptps)=betatpsX(iii)
                Out(iii,2+ptps)=betatpsminX(iii)
                Out(iii,3+ptps)=betatpsmaxX(iii)
                Out(iii,4+ptps)=dsqrt(varBetatps(iii))
            end do
        ptps = ptps + 4 ! decalage de 4 colonnes
        endif
        p=p+filtre(i)*(nbinnerknots+qorder-1)+1
    end do

    deallocate(betatpsX,betatpsminX,betatpsmaxX,varBetatps)
    deallocate(BasisSinhaTPS)

    end subroutine drawTimeCoef


    subroutine drawTimeCoefdc(np,b,nvar,filtre,Out)

    use comon
    use betatttps

    implicit none

    integer::p,i,iii,ptps,jj,jjj,np,p1
    double precision::tps,iiii
    integer::nvar
    double precision,dimension(np)::b
    integer,dimension(nvar)::filtre
    double precision,dimension(0:100,0:4*sum(filtre)),intent(out)::Out

    allocate(betatpsX(0:100),betatpsminX(0:100),betatpsmaxX(0:100),varBetatps(0:100))
    allocate(BasisSinhaTPS(nbinnerknots+qorder))

    p1=nva1+npbetatps1
    p=1
    ptps=0
    do i=1,nvar
        if (filtre(i).eq.1) then ! beta dependant du temps
            betatpsX=0.d0
            betatpsminX=0.d0
            betatpsmaxX=0.d0
            varBetatps=0.d0
            do iii=0,100
                iiii=dble(iii)
                tps=knotsTPS(0)+(knotsTPS(nbinnerknots+1)-knotsTPS(0))*(iiii/100.d0)

                Out(iii,0) = tps

                call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,nbinnerknots+qorder, &
                tps,innerknots,boundaryknots,basisSinhaTPS)

                do jj=-qorder+1,nbinnerknots
                    betatpsX(iii)=betatpsX(iii)+b(np-(nva+npbetatps)+(p-1)+p1+jj+qorder)*BasisSinhaTPS(jj+qorder)
                end do

                do jj=-qorder+1,nbinnerknots
                    do jjj=-qorder+1,nbinnerknots
                        varBetatps(iii)=varBetatps(iii)+BasisSinhaTPS(jj+qorder)*BasisSinhaTPS(jjj+qorder) &
                        *(H_hess(np-(nva+npbetatps)+(p-1)+p1+jj+qorder,np-(nva+npbetatps)+(p-1)+p1+jjj+qorder))
                    end do
                end do

                betatpsminX(iii)=betatpsX(iii)-1.96*dsqrt(varBetatps(iii))
                betatpsmaxX(iii)=betatpsX(iii)+1.96*dsqrt(varBetatps(iii))

                Out(iii,1+ptps)=betatpsX(iii)
                Out(iii,2+ptps)=betatpsminX(iii)
                Out(iii,3+ptps)=betatpsmaxX(iii)
                Out(iii,4+ptps)=dsqrt(varBetatps(iii))
            end do
        ptps = ptps + 4
        endif
        p=p+filtre(i)*(nbinnerknots+qorder-1)+1
    end do

    deallocate(betatpsX,betatpsminX,betatpsmaxX,varBetatps)
    deallocate(BasisSinhaTPS)

    end subroutine drawTimeCoefdc

    !**************************************
!*********** Determinant of a Matrix *********
!**************************************

    !Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
    double precision FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    double precision :: m, temp
    INTEGER :: i, j, k, l
    double precision, DIMENSION(n,n) :: matrix
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO

    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO

    END FUNCTION FindDet



!C==========================================================================
!C
!C
!C
!C     Fonctions pour calculer une fonction de repartition issue d'une
!C     loi normale
!C                                   ajout le 28/04/2014
!C
!C===========================================================================




      double precision function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    David Hill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
      implicit none

      real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
      real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
      real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
      real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
      real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
      real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
      real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
      real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
      real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
      real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
      real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
      real ( kind = 8 ), parameter :: con = 1.28D+00
      real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
      real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
      real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
      real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
      real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
      real ( kind = 8 ), parameter :: ltone = 7.0D+00
      real ( kind = 8 ), parameter :: p = 0.398942280444D+00
      real ( kind = 8 ), parameter :: q = 0.39990348504D+00
      real ( kind = 8 ), parameter :: r = 0.398942280385D+00
      logical up
      logical upper
      real ( kind = 8 ), parameter :: utzero = 18.66D+00
      real ( kind = 8 ) x
      real ( kind = 8 ) y
      real ( kind = 8 ) z

      up = upper
      z = x

      if ( z < 0.0D+00 ) then
         up = .not. up
         z = - z
      end if

      if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

         if ( up ) then
            alnorm = 0.0D+00
         else
            alnorm = 1.0D+00
         end if

         return

      end if

      y = 0.5D+00 * z * z

      if ( z <= con ) then

         alnorm = 0.5D+00 - z * ( p - q * y     &
   / ( y + a1 + b1/ ( y + a2 + b2 / ( y + a3 ))))

      else

         alnorm = r * exp ( - y )     &
      / ( z + c1 + d1 / ( z + c2 + d2 / ( z + c3 + d3     &
      / ( z + c4 + d4 / ( z + c5 + d5 / ( z + c6 ))))))

      end if

      if ( .not. up ) then
         alnorm = 1.0D+00 - alnorm
      end if

      return
      end function alnorm


! ------------------------------------------------
!     Cholesky decomposition.
     subroutine cholesky_sub(A,n)

    implicit none

  ! formal vars
      integer :: n      ! number of rows/cols in matrix
      double precision    :: A(n,n) ! matrix to be decomposed

  ! local vars
      integer :: j,i      ! iteration counter

  ! begin loop
   do j = 1,n

    ! perform diagonal component
    A(j,j) = sqrt(A(j,j) - dot_product(A(j,1:j-1),A(j,1:j-1)))

    ! perform off-diagonal component
    if (j < n) A(j+1:n,j) = (A(j+1:n,j) - matmul(A(j+1:n,1:j-1),A(j,1:j-1))) / &
              A(j,j)

     end do
     do i = 1,n
        do j = i+1 , n
          A(i,j) = 0.d0
    end do
  end do
end subroutine cholesky_sub
