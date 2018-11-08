
! Calcul te de la fonction de risque a un point t donne

!========================================================
!
!            Splines
!
!========================================================

    subroutine risque(t,the_s,the1_s,nz,zi_s,lam,nst)

    implicit none

    double precision,intent(in)::t
    integer,intent(in)::nz,nst
    double precision,dimension(nz),intent(in)::the_s,the1_s
    double precision,dimension(nz+4),intent(in)::zi_s
    double precision,dimension(2),intent(out)::lam
    integer::j,k
    double precision,dimension(-2:nz-3)::the,the1
    double precision,dimension(-2:nz+1)::zi
    double precision::ht,ht2,h2,htm,h2t,h3,h2n,hn,im,im1,im2,&
    mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,hh

    the=the_s
    zi=zi_s
    the1=the1_s
    lam=0.d0
    do k = 2,nz-1
        if ((t.ge.zi(k-1)).and.(t.lt.zi(k)))then
            j = k-1

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
            lam(1) = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)

            if (nst==2) then
                lam(2) = (the1(j-3)*mm3)+(the1(j-2)*mm2)+(the1(j-1)*mm1)+(the1(j)*mm)
            end if
        endif
    end do

    end subroutine risque

! nouvelle fonction risque pour plus de strates

    subroutine risque2(t,the0,nz,zi0,lam,nst)

    implicit none

    double precision,intent(in)::t
    integer,intent(in)::nz,nst
    double precision,dimension(nz,nst),intent(in)::the0
    double precision,dimension(nz+4),intent(in)::zi0
    double precision,dimension(nst),intent(out)::lam
    integer::j,k,l
    double precision,dimension(-2:nz-3,nst)::the
    double precision,dimension(-2:nz+1)::zi
    double precision::ht,ht2,h2,htm,h2t,h3,h2n,hn,im,im1,im2,&
    mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,hh

    the = the0
    zi = zi0
    lam = 0.d0
    do k = 2,nz-1
        if ((t.ge.zi(k-1)).and.(t.lt.zi(k)))then
            j = k-1

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
                lam(l) = (the(j-3,l)*mm3)+(the(j-2,l)*mm2)+(the(j-1,l)*mm1)+(the(j,l)*mm)
            end do
        endif
    end do

    end subroutine risque2
