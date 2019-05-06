!==========================  RISQUE   ====================================
!== permet avoir fonction de risque (lam) en tous points
!== aussi survie su01 ou a modifier pour avoir risq cumul� gl

! calcul les points pour les fonctions 
! et leur bandes de confiance

    subroutine risqueG(x,the,n,zi,su,lam)
    use tailles

    implicit none

    integer::j,k,n,i
    double precision::x,ht,ht2,h2,som,lam,su,h1,htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3,&
    ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,gl,hh
    double precision,dimension(-2:npmax)::the,zi
!      do i=1,n
!         the(i-3)=(b(i))*(b(i))
!      end do

    som = 0.d0
    ht3 = 0.d0
    hht = 0.d0
    h4 = 0.d0
    h3m = 0.d0
    hh3 = 0.d0
    hh2 = 0.d0
    mm = 0.d0
    im3 = 0.d0
    mm2 = 0.d0
    h = 0.d0
    gl = 0.d0
    hh = 0.d0
    do k = 2,n-1
        if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            j = k-1
            if (j.gt.1)then
                do i=2,j
                    som = som+the(i-4)
                end do  
            endif   
            ht = x-zi(j)
            htm= x-zi(j-1)
            h2t= x-zi(j+2)
            ht2 = zi(j+1)-x
            ht3 = zi(j+3)-x
            hht = x-zi(j-2)
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
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        endif
    end do

!--------- ATTENTION 
!=== si x est plus grand que le dernier noeud zi(), onsupposera
! lam,cumulhaz,su constant

    if(x.ge.zi(n))then
        som = 0.d0
        do i=1,n
            som = som+the(i-3)
        end do
        gl = som
!------------- for hazard
        i = n-2
        h1 = (zi(i)-zi(i-1))
        lam = (4.d0*the(i-1)/h1)
    endif

    su  = dexp(-gl)

    return

    end subroutine risqueG
!=============fin RISQUE




    double precision function integrant3(frail,varx,vary) 
! calcul de l integrant, pour un effet aleatoire donn� frail 
! sachant          un temps de progression varx et temps dc vary

    use tailles
    !use comon,only:AG,auxig,c,cdc,effet,ndate,ndatedc,nst,nsujet,nt0,nt1,&
    !nt0dc,nt1dc,nva,nva1,nva2,t0,t1,t0dc,t1dc,pe,aux1,aux2,g,nig,res1,res3
    use comon,only:theta,alpha,zi,nz1,nz2
    
    !use comongroup,only:gsuj,nigdc
    use comongroup,only:expb1,expb2,the1,the2

    implicit none

    double precision::hazrec,cumhazrec,hazdc,cumhazdc,frail,varx,vary,&
    su,lam,logGammaJ

! frail=u
! varx=valeur du temps x
! vary=valeur du temps y
! risque_recurrences(x)=hazrec 
! risque_deces(x)=hazdc
! risqcumul�e_recurrences(x)=cumhazrec 
! risqcumul�e_deces(x)=cumhazdc
! ajurec= exp(beta Z) for recurrences<<<<<<
! ajudc= exp(beta Z) for deaths
      
! for recurrences, calcul des fonctions de risques
    call risqueG(varx,the1,nz1+2,zi,su,lam)
    hazrec = lam 
    cumhazrec = -dlog(su)
!             open(4,file='outjointclust')
!             write(4,*)'dans RISQUE integran>t3',varx,su,lam,frail
! for deaths, calcul des fonctions de risques
    call risqueG(vary,the2,nz2+2,zi,su,lam)
    hazdc = lam 
    cumhazdc = -dlog(su)
!             open(4,file='outjointclust')
!         write(4,*)'integrant3: hazdc ',hazdc,cumhazdc
!         write(4,*)'integrant3: vary',vary
!         write(4,*)''

    integrant3 = frail*hazrec*expb1*dexp(-frail*cumhazrec*expb1) &
    *(frail**alpha)*hazdc*expb2*dexp((-frail**alpha)*cumhazdc*expb2) &
    *(frail**(1.d0/theta-1.d0))*(dexp(-frail/theta)) &
    /(dexp(logGammaJ(1.d0/theta))*theta**(1.d0/theta))


    if(integrant3.lt.0.d0)then
!         write(*,*)'dans integrant3',integrant3,frail,theta,alpha ,expb1,expb2,&
!         cumhazrec,cumhazdc,hazrec,hazdc,logGammaJ(1.d0/theta)
    endif

    return
    
    end function integrant3


    double precision function integrant4(frail,varx,vary) 
! calcul de l integrant, pour un effet aleatoire donn� frail et un groupe donne auxig (cf funcpa)      
    use tailles
    !use comon,only:AG,g,nig,auxig,c,cdc,nt0,nt1,nt0dc,&
    !nt1dc,nsujet,nva,nva1,nva2,ndate,ndatedc,nst,pe,effet,&
    !t0,t1,t0dc,t1dc,aux1,aux2,date,datedc,res1,res3
    use comon,only:theta,alpha,zi,nz1,nz2
    
    !use comongroup,only:vet,vet2,gsuj,nigdc
    use comongroup,only:expb1,expb2,the1,the2

    implicit none

    double precision::hazrec,cumhazrec,hazdc,cumhazdc,frail,varx,vary,su,lam
    double precision::logGammaJ

! frail=u
! varx=valeur du temps x
! vary=valeur du temps y
! risque_recurrences(x)=hazrec 
! risque_deces(x)=hazdc
! risqcumulee_recurrences(x)=cumhazrec 
! risqcumulee_deces(x)=cumhazdc
! ajurec= exp(beta Z) for recurrences
! ajudc= exp(beta Z) for deaths

! for recurrences, calcul des fonctions de risques
    call risqueG(varx,the1,nz1+2,zi,su,lam)
    hazrec = lam 
    cumhazrec = -dlog(su)
! for deaths, calcul des fonctions de risques
    call risqueG(vary,the2,nz2+2,zi,su,lam)
    hazdc = lam 
    cumhazdc = -dlog(su)
    integrant4 = -frail*cumhazrec *expb1-(frail**alpha)*cumhazdc*expb2

    integrant4 = exp(integrant4)
    integrant4 =  integrant4 * (frail**(1.d0/theta-1.d0))*(dexp(-frail/theta)) &
    /(dexp(logGammaJ(1.d0/theta))*theta**(1.d0/theta))

    return

    end function integrant4

! 20 points / sur (0,+infty) pour tau kendall
!== integrale sur xx (temps progression)
!== integrale sur yy (temps deces)
!== ss : donne le resultat de l'int�grale
!== gaulagkend34 fait appelle � deux int�grales

      subroutine gaulagKend34(ss,xx,yy,choix) 

    use tailles
    !use comon,only:auxig
    use donnees,only:w,x

    implicit none

    double precision::ss,xx,yy,auxfunca
    integer::j,choix,nan
    double precision,external::integrant3, integrant4
! gauss laguerre
! integrant3 ou 4 est l int�grant, ss le r�sultat de l integrale sur 0 ,  +infty
    ss=0.d0 
! Will be twice the average value of the function,since the ten
! wei hts (five numbers above each used twice) sum to 2.
    do j=1,20
        if (choix.eq.3) then   !integrale 3
            auxfunca=integrant3(x(j),xx,yy)
            ss = ss+w(j)*(auxfunca)
        endif
        if (choix.eq.4) then   !integrale 4
            auxfunca=integrant4(x(j),xx,yy)
            ss = ss+w(j)*(auxfunca)
            if(ss.eq.nan)then 
!     write(*,*)'----int1',ss ,w(j),auxfunca,j
!     stop
            endif
        endif
    end do

    return
    
    end subroutine gaulagKend34


! 20 points / sur (0,+infty) pour tau kendall
! integrale sur x (temps progression)
! pour un temps de deces yy donn�

    subroutine gaulagKend2(ss,yy) 

    use tailles
    !use comon,only:auxig
    use donnees,only:w,x

    implicit none

    integer::j
    double precision::ss,yy,ss3,ss4,auxfunca

    ss=0.d0 
! Will be twice the average value of the function,since the ten
! wei hts (five numbers above each used twice) sum to 2.
    do j=1,20
        call gaulagKend34(ss3,x(j),yy,3)
        call gaulagKend34(ss4,x(j),yy,4)
        auxfunca=ss3*ss4
        ss = ss+w(j)*(auxfunca)
    end do

    return
    
    end subroutine gaulagKend2


!==================================================================
!==================================================================
! 20 points / sur (0,+infty) pour tau kendall
!== integrale sur y (temps deces)

    subroutine gaulagKend1(ss) 

    use tailles
    !use comon,only:auxig
    use donnees,only:w,x

    implicit none

    integer::j
    double precision::ss,ss2,auxfunca
! gauss laguerre
! func1 est l int�grant, ss le r�sultat de l integrale sur 0 ,  +infty


    ss=0.d0 
! Will be twice the average value of the function,since the ten
! wei hts (five numbers above each used twice) sum to 2.
    do j=1,20
        call gaulagKend2(ss2,x(j)) 
!pour une valeur des temps de deces fix�e (xj)
        auxfunca=ss2 
        ss = ss+w(j)*(auxfunca)
    end do

    return
    
    end subroutine gaulagKend1


! 20 points / sur (0,+infty) pour tau kendall
!== integrale sur u_j (une fragilite)
!== pour une autre fragilit� donn�e u_i

    subroutine gaulagKend2bis(ss,ui) 

    use tailles
    !use comon,only:AG,auxig,nva,nva1,nva2,ndate,ndatedc,nst,pe,effet,t0,t1,t0dc,&
    !t1dc,c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet,date,datedc,
    use comon,only:zi,nz1,nz2,alpha,theta
    !use comongroup,only:expb1,expb2,vet,vet2
    use comongroup,only:the1,the2
    use donnees,only:w,x

    implicit none

    integer::j
    double precision::hazrec,cumhazrec,hazdc,cumhazdc,ss,ui,su,lam,auxfunca,logGammaJ
! gauss laguerre

    ss=0.d0 
! Will be twice the average value of the function,since the ten
! wei hts (five numbers above each used twice) sum to 2.
    do j=1,20
! for deaths, hazard functions
        call risqueG(x(j),the2,nz2+2,zi,su,lam)
        hazdc = lam 
        cumhazdc = -dlog(su)
! for TTP, hazard functions
        call risqueG(x(j),the1,nz1+2,zi,su,lam)
        hazrec = lam 
        cumhazrec = -dlog(su)
        auxfunca=(1.d0/(ui+x(j)))*(1.d0/(ui**alpha+x(j)**alpha))* &
        (ui**(alpha+1)+x(j)**(alpha+1))*(x(j)**(1.d0/theta-1.d0))*(dexp(-x(j)/theta)) &
        *(ui**(1.d0/theta-1.d0))*(dexp(-ui/theta))/(dexp(logGammaJ(1.d0/theta))*theta**(1.d0/theta))**2

        ss = ss+w(j)*(auxfunca)
    end do

    return

    end subroutine gaulagKend2bis




!============================================================
!==================================================================
! 20 points / sur (0,+infty) pour tau kendall
!== integrale sur ui (un effet aleatoire)

    subroutine gaulagKend1bis(ss) 

    use tailles
    !use comon,only:auxig
    use donnees,only:w,x

    implicit none

    integer::j
    double precision::ss,ss2,auxfunca
! gauss laguerre
! func1 est l int�grant, ss le r�sultat de l integrale sur 0 ,  +infty
    ss=0.d0 
! Will be twice the average value of the function,since the ten
! wei hts (five numbers above each used twice) sum to 2.
    do j=1,20
        call gaulagKend2bis(ss2,x(j)) 
!pour une valeur de la fragilit� fix�e (u_i)
        auxfunca=ss2 
        ss = ss+w(j)*(auxfunca)
    end do

    return

    end subroutine gaulagKend1bis
