
! Calcul de la fonction de survie a un point t donne


!AM:modification 12/07/2012-baseline hazard

!AM:add argument lam
    subroutine survival_frailty(t,the_s,the1_s,nz,zi_s,su,lam,nst)

    implicit none

    double precision,intent(in)::t
    integer,intent(in)::nz,nst
    double precision,dimension(nz),intent(in)::the_s,the1_s
    double precision,dimension(nz+4),intent(in)::zi_s
!AM:add lam
    double precision,dimension(2),intent(out)::su,lam
    integer::j,k,i
    double precision,dimension(-2:nz-3)::the,the1
    double precision,dimension(-2:nz+1)::zi
    double precision::ht,ht2,h2,som,htm,h2t,h3,h2n,hn,im,im1,&
    im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,hh,gl,som1,gl1


    the=the_s
    zi=zi_s
    the1=the1_s

    som = 0.d0
    som1 = 0.d0
    gl = 0.d0 
    gl1 = 0.d0 
    su = 0.d0
!AM: add lam
    lam = 0.d0
    
    do k = 2,nz ! avant y avait nz-1...
        if ((t.ge.zi(k-1)).and.(t.lt.zi(k)))then
            j = k-1
            if (j.gt.1)then
                do i=2,j
                    som = som+the(i-4)
                end do  
                do i=2,j
                    som1 = som1+the1(i-4)
                end do  
            endif

            ht = t-zi(j)
            htm= t-zi(j-1)
            h2t= t-zi(j+2)
            ht2 = zi(j+1)-t
            ht3 = zi(j+3)-t
            hht = t-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)

            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(t-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
!AM: add lam
            lam(1) = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
            if (nst == 2) then
                gl1 = som1 +(the1(j-3)*im3)+(the1(j-2)*im2)+(the1(j-1)*im1)+(the1(j)*im)
!AM: add lam
                lam(2) = (the1(j-3)*mm3)+(the1(j-2)*mm2)+(the1(j-1)*mm1)+(the1(j)*mm)
            end if
        endif
    end do
   
    if(t.ge.zi(nz))then
        som = 0.d0
        do i=1,nz
            som = som+the(i-3)
            som1 = som1+the1(i-3)
        end do
        gl = som
        gl1 = som1
    endif

    su(1) = dexp(-gl)
    if(nst == 2) su(2) = dexp(-gl1)

    end subroutine survival_frailty


   subroutine survival2(t,the0,nz,zi0,su,nst)

    implicit none

    double precision,intent(in)::t
    integer,intent(in)::nst,nz
    double precision,dimension(nz,nst),intent(in)::the0
    double precision,dimension(nz+4),intent(in)::zi0
    double precision,dimension(nst),intent(out)::su
    integer::j,k,i,l
    double precision,dimension(-2:nz-3,nst)::the
    double precision,dimension(nst)::somT,gl
    double precision,dimension(-2:nz+1)::zi
    double precision::ht,ht2,h2,htm,h2t,h3,h2n,hn,im,im1,&
    im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,hh


    the = the0
    zi = zi0

    somT = 0.d0
    gl = 0.d0
    su = 0.d0

    do k = 2,nz ! avant y avait nz-1...
        if ((t.ge.zi(k-1)).and.(t.lt.zi(k)))then
            j = k-1
            if (j.gt.1)then
                do i=2,j
                    do l=1,nst
                        somT(l) = somT(l)+the(i-4,l)
                    end do
                end do
            endif

            ht = t-zi(j)
            htm= t-zi(j-1)
            h2t= t-zi(j+2)
            ht2 = zi(j+1)-t
            ht3 = zi(j+3)-t
            hht = t-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)

            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(t-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0

            do l=1,nst
                gl(l) = somT(l) +(the(j-3,l)*im3)+(the(j-2,l)*im2)+(the(j-1,l)*im1)+(the(j,l)*im)
            end do

        endif
    end do

    if(t.ge.zi(nz))then
        somT = 0.d0
        do i=1,nz
            do l=1,nst
                somT(l) = somT(l)+the(i-3,l)
            end do
        end do
        do l=1,nst
            gl(l) = somT(l)
        end do
    endif

    do l=1,nst
        su(l) = dexp(-gl(l))
    end do

    end subroutine survival2


    subroutine survival_cpm(t,b,nst,nbintervR,time,surv)
    
    use optim
    
    implicit none

    integer::nbintervR,i,j,nst
    double precision::t,som1,som2,su
    double precision,dimension(nbintervR+1)::time
    double precision,dimension(0:nbintervR)::ttt
    double precision,dimension(2)::surv
    double precision,dimension(nbintervR*nst)::b

    ttt = time

    som1=0.d0
    som2=0.d0
    su=0.d0
    surv=0.d0

    do i=1,nbintervR
        if ((t.ge.(ttt(i-1))).and.(t.lt.(ttt(i)))) then
            som1=(b(i)**2)*(t-ttt(i-1))
            if (i.ge.2)then
                do j=1,i-1
                    som2=som2+(b(j)**2)*(ttt(j)-ttt(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
        if ((t.eq.(ttt(nbintervR)))) then
            som1=(b(nbintervR)**2)*(t-ttt(nbintervR-1))
            if (nbintervR.ge.2)then
                do j=1,nbintervR-1
                    som2=som2+(b(j)**2)*(ttt(j)-ttt(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
    end do

    surv(1) = su

    if(surv(1).lt.0.d0)then
        surv(1) = 0.d0
    endif

    if(surv(1).gt.1.d0)then
        surv(1) = 1.d0 
    endif

    if(nst == 2) then
!!! Ajout
        som1=0.d0
        som2=0.d0
        su=0.d0
!!!
        do i=1,nbintervR
            if ((t.ge.(ttt(i-1))).and.(t.lt.(ttt(i)))) then
                som1=(b(i+nbintervR)**2)*(t-ttt(i-1))
                if (i.ge.2)then
                    do j=1,i-1
                        som2=som2+(b(j+nbintervR)**2)*(ttt(j)-ttt(j-1))
                    end do
                endif
                su=dexp(-(som1+som2))
            endif
            if ((t.eq.(ttt(nbintervR)))) then
                som1=(b(2*nbintervR)**2)*(t-ttt(nbintervR-1))
                if (nbintervR.ge.2)then
                    do j=1,nbintervR-1
                        som2=som2+(b(j+nbintervR)**2)*(ttt(j)-ttt(j-1))
                    end do
                endif
                su=dexp(-(som1+som2))
            endif
        end do    

        surv(2) = su
    
        if(surv(2).lt.0.d0)then
            surv(2) = 0.d0
        endif
    
        if(surv(2).gt.1.d0)then
            surv(2) = 1.d0 
        endif
    end if

    end subroutine survival_cpm


    subroutine survivalj_cpm(t,b,nbintervR,nbintervDC,time,timedc,surv)
    
    use optim
    
    implicit none

    integer::nbintervR,nbintervDC,i,j
    double precision::t,som1,som2,su
    double precision,dimension(nbintervR+1)::time
    double precision,dimension(nbintervDC+1)::timedc
    double precision,dimension(0:nbintervR)::ttt
    double precision,dimension(0:nbintervDC)::tttdc
    double precision,dimension(2)::surv
    double precision,dimension(nbintervR+nbintervDC)::b

    ttt = time
    tttdc = timedc
    som1=0.d0
    som2=0.d0
    su=0.d0
    surv=0.d0

!recurrent
    do i=1,nbintervR
        if ((t.ge.(ttt(i-1))).and.(t.lt.(ttt(i)))) then
            som1=(b(i)**2)*(t-ttt(i-1))
            if (i.ge.2)then
                do j=1,i-1
                    som2=som2+(b(j)**2)*(ttt(j)-ttt(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
        if ((t.eq.(ttt(nbintervR)))) then
            som1=(b(nbintervR)**2)*(t-ttt(nbintervR-1))
            if (nbintervR.ge.2)then
                do j=1,nbintervR-1
                    som2=som2+(b(j)**2)*(ttt(j)-ttt(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
    end do    

    surv(1) = su

    if(surv(1).lt.0.d0)then
        surv(1) = 0.d0
    endif

    if(surv(1).gt.1.d0)then
        surv(1) = 1.d0 
    endif

!death
!!! Ajout
    som1=0.d0
    som2=0.d0
    su=0.d0
!!!
    do i=1,nbintervDC
        if ((t.ge.(tttdc(i-1))).and.(t.lt.(tttdc(i)))) then
            som1=(b(i+nbintervR)**2)*(t-tttdc(i-1))
            if (i.ge.2)then
                do j=1,i-1
                    som2=som2+(b(j+nbintervR)**2)*(tttdc(j)-tttdc(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
        if ((t.eq.(tttdc(nbintervDC)))) then
            som1=(b(nbintervDC+nbintervR)**2)*(t-tttdc(nbintervDC-1))
            if (nbintervDC.ge.2)then
                do j=1,nbintervDC-1
                    som2=som2+(b(j+nbintervR)**2)*(tttdc(j)-tttdc(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
    end do    

    surv(2) = su

    if(surv(2).lt.0.d0)then
        surv(2) = 0.d0
    endif

    if(surv(2).gt.1.d0)then
        surv(2) = 1.d0 
    endif


    end subroutine survivalj_cpm

! nouvelles fonctions survie pour plus de strates

    subroutine survival_cpm2(t,b,nst,nbintervR,time,surv)

    use optim

    implicit none

    integer::nbintervR,i,j,k,nst
    double precision::t,som1,som2,su
    double precision,dimension(nbintervR+1)::time
    double precision,dimension(0:nbintervR)::ttt
    double precision,dimension(nst)::surv
    double precision,dimension(nbintervR*nst)::b

    ttt = time

    surv=0.d0

    do k=1,nst
        som1=0.d0
        som2=0.d0
        su=0.d0
        do i=1,nbintervR
            if ((t.ge.(ttt(i-1))).and.(t.lt.(ttt(i)))) then
                som1=(b(i+(k-1)*nbintervR)**2)*(t-ttt(i-1))
                if (i.ge.2)then
                    do j=1,i-1
                        som2=som2+(b(j+(k-1)*nbintervR)**2)*(ttt(j)-ttt(j-1))
                    end do
                endif
                su=dexp(-(som1+som2))
            endif
            if ((t.eq.(ttt(nbintervR)))) then
                som1=(b(k*nbintervR)**2)*(t-ttt(nbintervR-1))
                if (nbintervR.ge.2)then
                    do j=1,nbintervR-1
                        som2=som2+(b(j+(k-1)*nbintervR)**2)*(ttt(j)-ttt(j-1))
                    end do
                endif
                su=dexp(-(som1+som2))
            endif
        end do

        surv(k) = su

        if(surv(k).lt.0.d0)then
            surv(k) = 0.d0
        endif

        if(surv(k).gt.1.d0)then
            surv(k) = 1.d0
        endif

    end do

    end subroutine survival_cpm2


    subroutine survivalj_cpm2(t,b,nst,nbintervR,nbintervDC,time,timedc,surv)

    use optim

    implicit none

    integer::nbintervR,nbintervDC,i,j,k,nst
    double precision::t,som1,som2,su
    double precision,dimension(nbintervR+1)::time
    double precision,dimension(nbintervDC+1)::timedc
    double precision,dimension(0:nbintervR)::ttt
    double precision,dimension(0:nbintervDC)::tttdc
    double precision,dimension(nst)::surv
    double precision,dimension(nst*nbintervR+nbintervDC)::b

    ttt = time
    tttdc = timedc

    surv=0.d0

!recurrent
    do k=1,nst-1
        som1=0.d0
        som2=0.d0
        su=0.d0
        do i=1,nbintervR
            if ((t.ge.(ttt(i-1))).and.(t.lt.(ttt(i)))) then
                som1=(b(i+(k-1)*nbintervR)**2)*(t-ttt(i-1))
                if (i.ge.2)then
                    do j=1,i-1
                        som2=som2+(b(j+(k-1)*nbintervR)**2)*(ttt(j)-ttt(j-1))
                    end do
                endif
                su=dexp(-(som1+som2))
            endif
            if ((t.eq.(ttt(nbintervR)))) then
                som1=(b(k*nbintervR)**2)*(t-ttt(nbintervR-1))
                if (nbintervR.ge.2)then
                    do j=1,nbintervR-1
                        som2=som2+(b(j+(k-1)*nbintervR)**2)*(ttt(j)-ttt(j-1))
                    end do
                endif
                su=dexp(-(som1+som2))
            endif
        end do

        surv(k) = su

        if(surv(k).lt.0.d0)then
            surv(k) = 0.d0
        endif

        if(surv(k).gt.1.d0)then
            surv(k) = 1.d0
        endif

    end do

!death
    som1=0.d0
    som2=0.d0
    su=0.d0

    do i=1,nbintervDC
        if ((t.ge.(tttdc(i-1))).and.(t.lt.(tttdc(i)))) then
            som1=(b(i+(nst-1)*nbintervR)**2)*(t-tttdc(i-1))
            if (i.ge.2)then
                do j=1,i-1
                    som2=som2+(b(j+(nst-1)*nbintervR)**2)*(tttdc(j)-tttdc(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
        if ((t.eq.(tttdc(nbintervDC)))) then
            som1=(b(nbintervDC+(nst-1)*nbintervR)**2)*(t-tttdc(nbintervDC-1))
            if (nbintervDC.ge.2)then
                do j=1,nbintervDC-1
                    som2=som2+(b(j+(nst-1)*nbintervR)**2)*(tttdc(j)-tttdc(j-1))
                end do
            endif
            su=dexp(-(som1+som2))
        endif
    end do

    surv(nst) = su

    if(surv(nst).lt.0.d0)then
        surv(nst) = 0.d0
    endif

    if(surv(nst).gt.1.d0)then
        surv(nst) = 1.d0
    endif


    end subroutine survivalj_cpm2

