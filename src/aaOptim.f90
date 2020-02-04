
!- Version fortran 90
!
!
      module type

      interface verif1
      subroutine marq98j(k0,b,m,ni,v,rl,ier,istop,effet,ca,cb,dd,fctnames)
         integer,intent(in) :: m,effet
         integer,intent(inout)::ni,ier,istop
         double precision,dimension(m*(m+3)/2),intent(out)::v
     double precision,dimension(2)::k0
         double precision,intent(out)::rl
         double precision,dimension(m),intent(inout)::b
     double precision,intent(out)::ca,cb,dd
     external::fctnames
     double precision::fctnames

      end subroutine marq98j
      
      ! IJ 2018 NCC design
      subroutine marq98jAlt(k0,b,m,ni,v,rl,ier,istop,effet,ca,cb,dd,fctnames,fctnamesindiv)
         integer,intent(in) :: m,effet
         integer,intent(inout)::ni,ier,istop
         double precision,dimension(m*(m+3)/2),intent(out)::v
     double precision,dimension(2)::k0
         double precision,intent(out)::rl
         double precision,dimension(m),intent(inout)::b
     double precision,intent(out)::ca,cb,dd
     external::fctnames, fctnamesindiv
     double precision::fctnames, fctnamesindiv
      end subroutine marq98jAlt

      subroutine derivaJ(b,m,v,rl,k0,fctnames)
        integer,intent(in)::m
    double precision,dimension(2)::k0
        double precision,intent(inout)::rl
        double precision,dimension(m),intent(in)::b
        double precision,dimension((m*(m+3)/2)),intent(out)::v
    double precision,external::fctnames
      end subroutine derivaJ

      !IJ 2018 NCC design
          subroutine derivaJindiv(b,m,vindiv,rl,k0,index,fctnamesindiv)
        integer,intent(in)::m,index
    double precision,dimension(2)::k0
        double precision,intent(inout)::rl
        double precision,dimension(m),intent(in)::b
        double precision,dimension(m),intent(out)::vindiv
    double precision,external::fctnamesindiv
      end subroutine derivaJindiv

      subroutine searpasj(vw,step,b,bh,m,delta,fim,k0,fctnames)
        integer,intent(in)::m
    double precision,dimension(2)::k0
        double precision,dimension(m),intent(in)::b
        double precision,dimension(m),intent(inout)::bh,delta
        double precision,intent(inout)::vw,fim,step
    double precision,external::fctnames
      end subroutine searpasj

      subroutine dmfsdj(a,n,eps,ier)
        integer,intent(in)::n
        integer,intent(inout)::ier
        double precision,intent(inout)::eps
        double precision,dimension(n*(n+1)/2),intent(inout)::A
      end subroutine dmfsdj

      subroutine valfpaj(vw,fi,b,bk,m,delta,k0,fctnames)
        integer,intent(in)::m
    double precision,dimension(2)::k0
    double precision,intent(in)::vw
        double precision,dimension(m),intent(in)::b,delta
        double precision,dimension(m),intent(out)::bk
        double precision,intent(out)::fi
    double precision,external::fctnames
      end subroutine valfpaj

      subroutine dmaxt(maxt,delta,m)
        integer,intent(in)::m
        double precision,dimension(m),intent(in)::delta
        double precision,intent(out)::maxt
      end subroutine dmaxt
      end interface verif1

      interface verif2
      subroutine dsinvj(A,N,EPS,IER)
        integer,intent(in)::n
        integer,intent(inout)::ier
        double precision,intent(inout)::eps
        double precision,dimension(n*(n+1)/2),intent(inout)::A
      end subroutine dsinvj

      subroutine dcholej(a,k,nq,idpos)
      integer,intent(in)::k,nq
      integer,intent(inout)::idpos
      double precision,dimension(k*(k+3)/2),intent(inout)::a
      end subroutine dcholej
      end interface verif2

      end module type
!-----------------------------------------------------------
! Derniere mis a jour : 09/02/2011
!-----------------------------------------------------------

!-------------------------------------------------------------
!                   MARQ98
!-------------------------------------------------------------
      module optim

      implicit none
! -Interface permettant la verification des type des arguments
      interface verif1
           module procedure marq98j,derivaJ,marq98jAlt,derivaJindiv,searpasj,dmfsdj,valfpaj
      end interface verif1

      interface verif2
        module procedure dsinvj,dcholej,dmaxt
      end interface verif2

      CONTAINS
!------------------------------------------------------------
!                   MARQ98
!-------------------------------------------------------------


    subroutine marq98j(k0,b,m,ni,v,rl,ier,istop,effet,ca,cb,dd,fctnames)
!
!  fu = matrice des derivees secondes et premieres
!
!  istop: raison de l'arret
!  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
!  2: nb max d'iterations atteints
!  4: Erreur
    !use residusM,only:indg
    use parameters
    !use comon,only: c,cdc,Hspl_hess,indic_eta,ndate,ndatedc,nst,nsujet,&
    !nt0,nt0dc,nt1,nt1dc,nva1,nva2,PEN_deri,t0,t0dc,t1,t1dc
    use comon,only:nva,model,I_hess,H_hess,hess,indic_ALPHA,&
    typeof,vvv!,indic_tronc

!add additive
    !use additiv,only:correl

    IMPLICIT NONE
!   variables globales
    integer,intent(in) :: m,effet
    integer,intent(inout)::ni,ier,istop
    double precision,dimension(m*(m+3)/2),intent(out)::v
    double precision,intent(out)::rl
    double precision,dimension(m),intent(inout)::b
    double precision,intent(out)::ca,cb,dd
    double precision,dimension(2)::k0
    double precision,dimension(2)::zero
!   variables locales
    integer::nql,ii,nfmax,idpos,ncount,id,jd,m1,j,i,ij,k 
    character(len=100)::bar ! Pour la progressBar
    double precision,dimension(m*(m+3)/2)::fu,v1,vnonpen
    double precision,dimension(m)::delta,b1,bh
    double precision::da,dm,ga,tr
    double precision::GHG,step,eps,vw,fi,maxt, &
    z,rl1,th,ep
    external::fctnames
    double precision::fctnames
!---------- ajout
    integer::kkk
    
    zero=0.d0
    id=0
    jd=0
    z=0.d0
    th=1.d-5
    eps=1.d-7!1.d-6
    nfmax=m*(m+1)/2
    ca=epsa+1.d0
    cb=epsb+1.d0
    rl1=-1.d+10
    ni=0
    istop=0
    da=0.01d0
    dm=5.d0
    nql=1
    m1=m*(m+1)/2
    ep=1.d-20
    
    
    Main:Do       
        
        call derivaJ(b,m,v,rl,k0,fctnames)         
        
    rl1=rl
    
        !print *,"rl1 = ", rl1
        if(rl.eq.-1.D9) then
            istop=4
            goto 110
        endif
        !write(*,*)'iteration***',ni,'vrais',rl
        
        !Myriam ! Construction de la progressBar
        if (((effet == 2).and.(model == 1)).or.(model == 7)) then 
            if(maxiter < 100) then
                bar(1:3) = "0%|"    
                do k=1, maxiter
                    bar(3+k:3+k)="-"
                enddo
                bar(maxiter+5:maxiter+5) = "|"
                do k=maxiter+6, 100
                    bar(k:k+1)=" "
                enddo
                
                do k=1, ni
                    bar(3+k:3+k)="*"
                enddo
                call intpr(bar, -1, ni, 0)
                call intpr('Iteration:', -1, ni, 1)
            else
                bar(1:3) = "0%|"    
                do k=1, 70
                    bar(3+k:3+k)="-"
                enddo
                bar(73:100) = "|                          "
                
                if (modulo(ni,3) == 0) then 
                    do k=1, ni/3
                        bar(3+k:3+k)="*"
                    enddo
                    if (ni < 300) then 
                        call intpr(bar, -1, ni, 0)
                        call intpr('Iteration:', -1, ni, 1) 
                    else 
                        call intpr('Iteration:', -1, ni, 1)
                    endif
                endif
            endif            
        endif
    
    ! Fin de la progressBar
    
    !write(*,*) 'iteration ***', ni, 'vecteur b ', b
    !write(*,*) 'm', m, 'iun', iun, 'z', z, 'k0', k0

    dd = 0.d0
    fu = 0.d0

    do i=1,m
        do j=i,m
            ij=(j-1)*j/2+i
            fu(ij)=v(ij)
        end do
    end do

        call dsinvj(fu,m,ep,ier)

    if (ier.eq.-1) then ! hessienne non inversible
        !print*,"here"
        dd=epsd+1.d0
    else
        GHG = 0.d0
        do i=1,m
            do j=1,m
                if(j.ge.i) then
                    ij=(j-1)*j/2+i
                else
                    ij=(i-1)*i/2+j
                end if
                GHG = GHG + v(m1+i)*fu(ij)*V(m1+j)
            end do
        end do
        dd=GHG/dble(m)
    end if

!write(*,*)'epsa',ca,'epsb',cb,'epsd',dd
!   write(*,*)'aaOptim: b(knots)',b(1:12)
! write(*,*) 'aaOptim: m',m, nva, nva1, nva2
  !! write(*,*)'aaOptim: b(random):', b((m-nva-3):(m-nva))
   !!write(*,*)'aaOptim: b(covar)',b((m-nva+1):m)

!     print*,ca,cb,dd
!     print*,b
!    print*,"-------------------------"
!     if (ni.eq.50) then
!         print*,"====================> FAUSSE CONVERGENCE"
!         ca=1.d-5
!         cb=1.d-5
!         dd=1.d-5
!     endif

    if(ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) exit main

        tr=0.d0
        do i=1,m
            ii=i*(i+1)/2
            tr=tr+dabs(v(ii))
        end do
        tr=tr/dble(m)

        ncount=0
        ga=0.01d0

 400    do i=1,nfmax+m
           fu(i)=v(i)
        end do
        
        do i=1,m
            ii=i*(i+1)/2
            if (v(ii).ne.0) then
                fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
            else
                fu(ii)=da*ga*tr
            endif
        end do
        
        call dcholej(fu,m,nql,idpos)

        if (idpos.ne.0) then
            ncount=ncount+1
            if (ncount.le.3.or.ga.ge.1.d0) then
                da=da*dm
            else
                ga=ga*dm
                if (ga.gt.1.d0) ga=1.d0
            endif

            goto 400

        else
        
            do i=1,m
            
                delta(i)=fu(nfmax+i)
                b1(i)=b(i)+delta(i)
            end do
            
            rl=fctnames(b1,m,id,z,jd,z,k0)
            
            if(rl.eq.-1.D9) then
                istop=4
                goto 110
            end if
            if (rl1.lt.rl) then
                if(da.lt.eps) then
                    da=eps
                else
                    da=da/(dm+2.d0)
                endif
                goto 800
            endif
        endif
!      write(6,*) 'loglikelihood not improved '
        call dmaxt(maxt,delta,m)
                
        if(maxt.eq.0.D0) then
            vw=th
        else
            !call dmaxt(maxt,delta,m)
            vw=th/maxt

        endif
        step=dlog(1.5d0)
        call searpasj(vw,step,b,bh,m,delta,fi,k0,fctnames)
        
        rl=-fi
        if(rl.eq.-1.D9) then
            istop=4
            goto 110
        end if
        do i=1,m
            delta(i)=vw*delta(i)
        end do
        da=(dm-3.d0)*da

 800     cb=dabs(rl1-rl)

        ca=0.d0
        do i=1,m
            ca=ca+delta(i)*delta(i)
        end do
        !write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd

        !print*,"rl1=",rl1
        !print*,"ca =",ca," cb =",cb," dd =",dd
        do i=1,m
            b(i)=b(i)+delta(i)
        end do
        
        ni=ni+1
        if (ni.ge.maxiter) then
            istop=2
!            write(6,*) 'maximum number of iteration reached'
            goto 110
        end if
    End do Main
    
    v=0.D0

    v(1:m*(m+1)/2)=fu(1:m*(m+1)/2)

    istop=1

!================ pour les bandes de confiance
!==== on ne retient que les para des splines
        
    call derivaJ(b,m,v,rl,k0,fctnames)
     !print*,"dans aaOptim, appel de derivaJ  with GAP vec funcpa... "
    if(rl.eq.-1.D9) then
        istop=4
        goto 110
    end if

    do i=1,(m*(m+3)/2)
        v1(i)=0.d0
    end do
!---- Choix du model

    select case(model)
        case(1,7)
            m1=m-nva-effet-indic_alpha !joint
        case(6)  !Uni 
            m1 = m - nva - effet*2 - indic_alpha*2
        !cas supplementaire rajout? pour joint general
        case(5)
            m1=m-nva-effet-indic_alpha !joint general (du moment que indic_eta=indic_alpha=1)
        case(2)
            m1=m-nva-effet*2 !additive
        case(3)
            m1=m-nva-effet !nested
        case(4)
            m1=m-nva-effet !shared
    end select

    kkk=m1*(m1+1)/2

    do i=1,m1
        kkk=kkk+1
        do j=i,m1
            k = (((j-1)*j)/2) +i
            v1(k)=v(k)/(4.d0*b(i)*b(j))
        end do
        v1(kkk)=v1(kkk)+(v(kkk)/(4.d0*b(i)*b(i)*b(i)))
    end do
    
    ep=10.d-10
    call dsinvJ(v1,m1,ep,ier)

    if (ier.eq.-1) then
!         write(*,*)   'echec inversion matrice information'
        istop=3
    endif



    do i=1,m1
        do j=i,m1
            hess(i,j)=v1((j-1)*j/2+i)
        end do
    end do

    do i=2,m1
        do j=1,i-1
            hess(i,j)=hess(j,i)
        end do
    end do

    ep=10.d-10
    call dsinvJ(v,m,ep,ier)

    if (ier.eq.-1) then
        istop=3

!AD:
!        call dsinvj(v1,m1,ep,ier)
!        if (ier.eq.-1) then
!             write(*,*)'echec inversion matrice information
!     & prms fixes'
!            istop=31
!        else
!            DO k=1,m1*(m1+1)/2
!                v(k)=v1(k)
!            END DO
!        end if
! fin ajout amadou
    endif
    
    ep=10.d-10
    call derivaJ(b,m,vnonpen,rl,zero,fctnames)

    do i=1,m
        do j=i,m
            I_hess(i,j)=vnonpen((j-1)*j/2+i)
        end do
    end do

    do i=2,m
        do j=1,i-1
            I_hess(i,j)=I_hess(j,i)
        end do
    end do

!========================================================

!   H_hess est moins la hessienne inverse sur la vraisemblance penalisee
    do i=1,m
        do j=i,m
            H_hess(i,j)=v((j-1)*j/2+i)
        end do
    end do
!       write(*,*) 'H_hess(16,16) fin marq',H_hess(16,16),m,((j-1)*j/2+i)

    do i=2,m
        do j=1,i-1
            H_hess(i,j)=H_hess(j,i)
        end do
    end do

 !AD:
    if (typeof .ne. 0) then
        do i=1,m*(m+1)/2
            vvv(i)=v(i)
        end do
    end if

 110   continue
 
       return
       end subroutine marq98j

!------------------------------------------------------------
!                          DERIVA
!------------------------------------------------------------

    subroutine derivaJ(b,m,v,rl,k0,fctnames)
    !use comon,only:nst,k0T
    use comon,only:model!,nstRec, indic_tronc
    implicit none

    integer,intent(in)::m
    double precision,intent(inout)::rl
    double precision,dimension(2)::k0
    double precision,dimension(m),intent(in)::b
    double precision,dimension((m*(m+3)/2)),intent(inout)::v
    double precision,dimension(m)::fcith
    integer ::i0,m1,ll,i,k,j,iun
    double precision::fctnames,thn,th,z,vl,th2,vaux
    external::fctnames
    logical::endDeriva

    endDeriva=.false.
    fcith=0.d0
    select case(model)
    case(1,7)
        th=1.d-3 !joint
    case(5)
        th=1.d-5 !joint general !cas supplementaire rajout?
    case(2)
        th=5.d-3 !additive
    case(3)
        th=1.d-5 !nested
    case(4)
        th=1.d-5 !shared
    end select

    thn=-th
    th2=th*th
    z=0.d0
    i0=0
    iun =1

        !print *,"k0 = ",k0
        !print *,"b = ", b
        !size(b) =32 , m =32
        !print *,"z = ", z            ! =0.000000000
        !print *,"iun", iun           ! =1
        !print *,"rl = ",rl           ! =0.000000000
    
    rl=fctnames(b,m,iun,z,iun,z,k0)
    
    if(rl.ne.-1.d9) then
        do i=1,m
            if((fcith(i).ne.-1.d9).and.(endDeriva.eqv..false.)) then
                fcith(i)=fctnames(b,m,i,th,i0,z,k0)
            else
                rl=-1.d9
                endDeriva=.true.
            end if
        end do
                
        if (endDeriva.eqv..false.) then
            k=0
            m1= (m*(m+1)/2)
            ll=m1

            do i=1,m
                !print*, "deuxieme boucle DERIVAj, tour",i," :"
                ll=ll+1
                !write(*,*) 'b', b
                vaux=fctnames(b,m,i,thn,i0,z,k0)
                !write(*,*) 'vaux', vaux
                if(vaux.eq.-1.d9) then
                    rl=-1.d9
                    endDeriva=.true.
                end if
                if (endDeriva.eqv..false.) then
                    vl=(fcith(i)-vaux)/(2.d0*th)
                    v(ll)=vl
                    do j=1,i
                        k=k+1
                        v(k)=-(fctnames(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2
                !        write(*,*)k,fctnames(b,m,i,th,j,th,k0),b,m,i,th,j,k0
                    end do
                end if
            end do
        end if
    end if
        
    !stop
    return

    end subroutine derivaJ
!------------------------------------------------------------
!                        SEARPAS
!------------------------------------------------------------

      subroutine searpasj(vw,step,b,bh,m,delta,fim,k0,fctnames)
!
!  MINIMISATION UNIDIMENSIONNELLE
!
      implicit none

      integer,intent(in)::m
      double precision,dimension(m),intent(in)::b
      double precision,intent(inout)::vw
      double precision,dimension(m),intent(inout)::bh,delta
      double precision,intent(inout)::fim,step
      double precision::vlw,vlw1,vlw2,vlw3,vm,fi1,fi2,fi3
      double precision,dimension(2)::k0
      double precision::fctnames
      external::fctnames
      integer::i

       vlw1=dlog(vw)
       vlw2=vlw1+step
       call valfpaj(vlw1,fi1,b,bh,m,delta,k0,fctnames)
       call valfpaj(vlw2,fi2,b,bh,m,delta,k0,fctnames)

       if(fi2.ge.fi1) then
      vlw3=vlw2
      vlw2=vlw1
      fi3=fi2
      fi2=fi1
      step=-step

          vlw1=vlw2+step
          call valfpaj(vlw1,fi1,b,bh,m,delta,k0,fctnames)
          if(fi1.gt.fi2) goto 50
       else
          vlw=vlw1
          vlw1=vlw2
          vlw2=vlw
          fim=fi1
          fi1=fi2
          fi2=fim
       end if

       do i=1,40
          vlw3=vlw2
          vlw2=vlw1
          fi3=fi2
          fi2=fi1

          vlw1=vlw2+step
          call valfpaj(vlw1,fi1,b,bh,m,delta,k0,fctnames)
          if(fi1.gt.fi2) goto 50
          if(fi1.eq.fi2) then
             fim=fi2
             vm=vlw2
             goto 100
          end if
       end do
!
!  PHASE 2 APPROXIMATION PAR QUADRIQUE
!
50     continue
!
!  CALCUL MINIMUM QUADRIQUE
!
      vm=vlw2-step*(fi1-fi3)/(2.d0*(fi1-2.d0*fi2+fi3))
      call valfpaj(vm,fim,b,bh,m,delta,k0,fctnames)
      if(fim.le.fi2) goto 100
      vm=vlw2
      fim=fi2
100   continue
      vw=dexp(vm)

      return

      end subroutine searpasj

!------------------------------------------------------------
!                         DCHOLE
!------------------------------------------------------------

      subroutine dcholej(a,k,nq,idpos)

      implicit none

      integer,intent(in)::k,nq
      integer,intent(inout)::idpos
      double precision,dimension(k*(k+3)/2),intent(inout)::a

      integer::i,ii,i1,i2,i3,m,j,k2,jmk
      integer::ijm,irm,jji,jjj,l,jj,iil,jjl,il
      integer,dimension(k)::is
      double precision ::term,xn,diag,p
      equivalence (term,xn)


!      ss programme de resolution d'un systeme lineaire symetrique
!
!       k ordre du systeme /
!       nq nombre de seconds membres
!
!       en sortie les seconds membres sont remplaces par les solutions
!       correspondantes
!

      i2=0
      ii=0
      idpos=0
      k2=k+nq
!     calcul des elements de la matrice
      do i=1,k
         ii=i*(i+1)/2
!       elements diagonaux
         diag=a(ii)
         i1=ii-i
         if(i-1.ne.0) goto 1
         if(i-1.eq.0) goto 4
1        i2=i-1
         do l=1,i2
             m=i1+l
             p=a(m)
             p=p*p
             if(is(l).lt.0) goto 2
             if(is(l).ge.0) goto 3
2            p=-p
3            diag=diag-p
         end do

4        if(diag.lt.0) goto 5
         if(diag.eq.0) goto 50
         if(diag.gt.0) goto 6
5        is(i)=-1
         idpos=idpos+1
         diag=-dsqrt(-diag)
         a(ii)=-diag
         goto 7
6        is(i)=1
         diag=dsqrt(diag)
         a(ii)=diag
!       elements non diagonaux
7        i3=i+1
         do j=i3,k2
            jj=j*(j-1)/2+i
            jmk=j-k-1
            if(jmk.le.0) goto 9
            if(jmk.gt.0) goto 8
8           jj=jj-jmk*(jmk+1)/2
9           term=a(jj)
            if(i-1.ne.0) goto 10
            if(i-1.eq.0) goto 13
10          do l=1,i2
               iil=ii-l
               jjl=jj-l
               p=a(iil)*a(jjl)
               il=i-l
               if(is(il).lt.0) goto 11
               if(is(il).ge.0) goto 12
11             p=-p
12             term=term-p
            end do
13            a(jj)=term/diag
       end do
      end do

!       calcul des solutions
      jj=ii-k+1
      do l=1,nq
         jj=jj+k
         i=k-1
14       jji=jj+i
         xn=a(jji)
         if(i-k+1.lt.0) goto 20
         if(i-k+1.ge.0) goto 22
20       j=k-1
21       jjj=jj+j
         ijm=i+1+j*(j+1)/2
         xn=xn-a(jjj)*a(ijm)
         if(j-i-1.le.0) goto 22
         if(j-i-1.gt.0) goto 30
30       j=j-1
         goto 21
22       irm=(i+1)*(i+2)/2
         a(jji)=xn/a(irm)
         if(i.le.0) cycle
         if(i.gt.0) goto 40
40       i=i-1
         go to 14
      end do
50    continue
      return
      end subroutine dcholej


      subroutine dmfsdj(a,n,eps,ier)
!
!   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
!   MATRICE = TRANSPOSEE(T)*T
!   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
!            PAR COLONNE DE LA METRICE A FACTORISER
!   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
!
!   SUBROUTINE APPELE PAR DSINV
!
!   N : DIM. MATRICE
!   EPS : SEUIL DE TOLERANCE
!   IER = 0 PAS D'ERREUR
!   IER = -1 ERREUR
!   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
!
      implicit none

      integer,intent(in)::n
      integer,intent(inout)::ier
      double precision,intent(inout)::eps
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision :: dpiv,dsum,tol
      integer::i,k,l,kpiv,ind,lend,lanf,lind

!
!   TEST ON WRONG INPUT PARAMETER N
!
      dpiv=0.d0
      if (n-1.lt.0) goto 12
      if (n-1.ge.0) ier=0
!
!   INITIALIZE DIAGONAL-LOOP
!
      kpiv=0
      do k=1,n
          kpiv=kpiv+k
          ind=kpiv
          lend=k-1
!
!   CALCULATE TOLERANCE
!
          tol=dabs(eps*sngl(A(kpiv)))
!
!   START FACTORIZATION-LOOP OVER K-TH ROW
!
         do i=k,n
        dsum=0.d0
            if (lend.lt.0) goto 2
            if (lend.eq.0) goto 4
            if (lend.gt.0) goto 2
!
!   START INNER LOOP
!
2           do l=1,lend
               lanf=kpiv-l
               lind=ind-l
           dsum=dsum+A(lanf)*A(lind)
            end do

!
!   END OF INNEF LOOP
!
!   TRANSFORM ELEMENT A(IND)
!
4           dsum=A(ind)-dsum
            if (i-k.ne.0) goto 10
            if (i-k.eq.0) goto 5
!
!   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
!
5           if (sngl(dsum)-tol.le.0) goto 6
            if (sngl(dsum)-tol.gt.0) goto 9
6           if (dsum.le.0) goto 12
            if (dsum.gt.0) goto 7
7           if (ier.le.0) goto 8
            if (ier.gt.0) goto 9
8           ier=k-1
!
!   COMPUTE PIVOT ELEMENT
!
9           dpiv=dsqrt(dsum)
            A(kpiv)=dpiv
            dpiv=1.D0/dpiv
            goto 11
!
!   CALCULATE TERMS IN ROW
!
10          A(ind)=dsum*dpiv
11          ind=ind+i
         end do
      end do

!
!   END OF DIAGONAL-LOOP
!
      return
12    ier=-1
      return

      end subroutine dmfsdj


!------------------------------------------------------------
!                            DSINV
!------------------------------------------------------------


    subroutine dsinvj(A,N,EPS,IER)

!
!     INVERSION D'UNE MATRICE SYMETRIQUE DEFINIE POSITIVE :
!
!     MATRICE = TRANSPOSEE(T)*T
!     INERSE(MATRICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
!
!     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA MATRICE A INVERSER
!         STOCKEE COLONNE PAR COLONNE
!     DIM. MATRICE A INVERSER = N
!     DIM. TABLEAU A = N*(N+1)/2
!
!     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
!           COMME NUL
!
!     IER : CODE D'ERREUR
!         IER=0 PAS D'ERREUR
!         IER=-1 ERREUR SUR LA DIM.N OU MATRICE PAS DEFINIE POSITIVE
!         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
!
    implicit none

    integer,intent(in)::n
    integer,intent(inout)::ier
    double precision,intent(inout)::eps
    double precision,dimension(n*(n+1)/2),intent(inout)::A
    double precision::din,work
    integer::ind,ipiv,i,j,k,l,min,kend,lhor,lver,lanf

!
!     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
!     A=TRANSPOSE(T) * T
!

    call dmfsdj(A,n,eps,ier)



    if (ier.lt.0) goto 9
!      if (ier.ge.0) det=0.d0
!
!     INVERT UPPER TRIANGULAR MATRIX T
!     PREPARE INVERSION-LOOP
!
!
! calcul du log du determinant

!      do i=1,n
!         det=det+dlog(A(i*(i+1)/2))
!      end do
!      det=2*det
    ipiv=n*(n+1)/2
    ind=ipiv
!
!     INITIALIZE INVERSION-LOOP
!
    do i=1,n
        din=1.d0/A(ipiv)
        A(ipiv)=din
        min=n
        kend=i-1
        lanf=n-kend
        if (kend.le.0) goto 5
        if (kend.gt.0) j=ind
    !
!     INITIALIZE ROW-LOOP
!
        do k=1,kend
            work=0.d0
            min=min-1
            lhor=ipiv
            lver=j
!
!     START INNER LOOP
!
            do l=lanf,min
                lver=lver+1
                lhor=lhor+l
                work=work+A(lver)*A(lhor)
            end do
!
!     END OF INNER LOOP
!
            A(j)=-work*din
            j=j-min
        end do

!
!     END OF ROW-LOOP
!
5        ipiv=ipiv-min
        ind=ind-1
    end do

!
!     END OF INVERSION-LOOP
!
!     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
!     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
!     INITIALIZE MULTIPLICATION-LOOP
!
    do i=1,n
        ipiv=ipiv+i
        j=ipiv
!
!     INITIALIZE ROW-LOOP
!
        do k=i,n
            work=0.d0
            lhor=j
!
!     START INNER LOOP
!
            do l=k,n
                lver=lhor+k-i
                work=work+A(lhor)*A(lver)
                lhor=lhor+l
            end do
!
!     END OF INNER LOOP
!
            A(j)=work
            j=j+k
        end do
    end do

!
!     END OF ROW-AND MULTIPLICATION-LOOP
!
9     return
      end subroutine dsinvj

!------------------------------------------------------------
!                          VALFPA
!------------------------------------------------------------

    subroutine valfpaj(vw,fi,b,bk,m,delta,k0,fctnames)

    implicit none

    integer,intent(in)::m
    double precision,dimension(m),intent(in)::b,delta
    double precision,dimension(m),intent(out)::bk
    double precision,intent(out)::fi
    double precision::vw,fctnames,z
    double precision,dimension(2)::k0
    integer::i0,i
    external::fctnames

    z=0.d0
    i0=1
    do i=1,m
    bk(i)=b(i)+dexp(vw)*delta(i)
    end do
    fi=-fctnames(bk,m,i0,z,i0,z,k0)
    if(fi.eq.-1.D9) then
        goto 1
    end if
1       continue
    return

    end subroutine valfpaj

!------------------------------------------------------------
!                            MAXT
!------------------------------------------------------------
    subroutine dmaxt(maxt,delta,m)

    implicit none

    integer,intent(in)::m
    double precision,dimension(m),intent(in)::delta
    double precision,intent(out)::maxt
    integer::i

    maxt=Dabs(delta(1))

    do i=2,m
        if(Dabs(delta(i)).gt.maxt)then
            maxt=Dabs(delta(i))
        end if

    end do

    return
    end subroutine dmaxt

!----------- add for multiv Yassin

    subroutine marq98(k0,b,m,ni,v,rl,ier,istop,effet,ca,cb,dd,fctnames)

!
!  fu = matrice des derivees secondes et premieres
!
!  istop: raison de l'arret
!  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
!  2: nb max d'iterations atteints
!  4: Erreur
    !use residusMmultiv,only:indg
    use parametersmultiv
    !use comonmultiv,only:c,cdc,Hspl_hess,indic_eta,ndate,&
    !ndatedc,nst,nsujet,nt0,nt0dc,nt1,nt1dc,nva1,nva2,&
    !PEN_deri,t0,t0dc,t1,t1dc
    use comonmultiv,only: nva,model,I_hess,H_hess,hess,indic_ALPHA,&
    typeof,vvv

!add additive
    !use additiv,only:correl

    IMPLICIT NONE
!   variables globales
    integer,intent(in) :: m,effet
    integer,intent(inout)::ni,ier,istop
    double precision,dimension(m*(m+3)/2),intent(out)::v
    double precision,intent(out)::rl
    double precision,dimension(m),intent(inout)::b
    double precision,intent(out)::ca,cb,dd
    double precision,dimension(3)::k0
        double precision,dimension(3)::zero
!   variables locales
    integer::nql,ii,nfmax,idpos,ncount,id,jd,m1,j,i,ij,k
    double precision,dimension(m*(m+3)/2)::fu,v1,vnonpen
    double precision,dimension(m)::delta,b1,bh
    double precision::da,dm,ga,tr
    double precision::GHG,step,eps,vw,fi,maxt, &
    z,rl1,th,ep
    external::fctnames
    double precision::fctnames
!---------- ajout
    integer::kkk

    zero=0.d0
    id=0
    jd=0
    z=0.d0
    th=1.d-5
    eps=1.d-7!1.d-6
    nfmax=m*(m+1)/2
    ca=epsa+1.d0
    cb=epsb+1.d0
    rl1=-1.d+10
    ni=0
    istop=0
    da=0.01d0
    dm=5.d0
    nql=1
    m1=m*(m+1)/2
    ep=1.d-20
    Main:Do

        call deriva(b,m,v,rl,k0,fctnames)
        rl1=rl
    if(rl.eq.-1.D9) then
        istop=4
        goto 110
    end if
!    write(*,*)'parameters***',b
!     write(*,*)'Optim iteration***',ni,'vrais',rl

        dd = 0.d0

        fu=0.D0

    do i=1,m
        do j=i,m
            ij=(j-1)*j/2+i
            fu(ij)=v(ij)
        end do
    end do

        call dsinvj(fu,m,ep,ier)

    if (ier.eq.-1) then
        dd=epsd+1.d0
    else
        GHG = 0.d0
        do i=1,m
            do j=1,m
                if(j.ge.i) then
                    ij=(j-1)*j/2+i
                else
                    ij=(i-1)*i/2+j
                end if
                GHG = GHG + v(m1+i)*fu(ij)*V(m1+j)
            end do
        end do
        dd=GHG/dble(m)
    end if

!    print*,ca,cb,dd
    if(ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) exit main

        tr=0.d0
        do i=1,m
            ii=i*(i+1)/2
            tr=tr+dabs(v(ii))
        end do
        tr=tr/dble(m)

        ncount=0
        ga=0.01d0

 400    do i=1,nfmax+m
           fu(i)=v(i)
        end do

        do i=1,m
            ii=i*(i+1)/2
            if (v(ii).ne.0) then
                fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
            else
                fu(ii)=da*ga*tr
            endif
        end do

        call dcholej(fu,m,nql,idpos)
        if (idpos.ne.0) then
            ncount=ncount+1
            if (ncount.le.3.or.ga.ge.1.d0) then
                da=da*dm
            else
                ga=ga*dm
                if (ga.gt.1.d0) ga=1.d0
            endif

            goto 400

        else
            do i=1,m
                delta(i)=fu(nfmax+i)
                b1(i)=b(i)+delta(i)
            end do
            rl=fctnames(b1,m,id,z,jd,z,k0)
            if(rl.eq.-1.D9) then
                istop=4
                goto 110
            end if
            if (rl1.lt.rl) then
                if(da.lt.eps) then
                    da=eps
                else
                    da=da/(dm+2.d0)
                endif
                goto 800
            endif
        endif
!      write(6,*) 'loglikelihood not improved '
        call dmaxt(maxt,delta,m)

        if(maxt.eq.0.D0) then
            vw=th
        else
            !call dmaxt(maxt,delta,m)
            vw=th/maxt

        endif
        step=dlog(1.5d0)
!      write(*,*) 'searpas'
        call searpas(vw,step,b,bh,m,delta,fi,k0,fctnames)
        rl=-fi
        if(rl.eq.-1.D9) then
            istop=4
            goto 110
        end if

        do i=1,m
            delta(i)=vw*delta(i)
        end do
        da=(dm-3.d0)*da

 800     cb=dabs(rl1-rl)
        ca=0.d0
        do i=1,m
            ca=ca+delta(i)*delta(i)
        end do
!         write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd
        do i=1,m
            b(i)=b(i)+delta(i)
        end do

        ni=ni+1

        if (ni.ge.maxiter) then
            istop=2
 !           write(6,*) 'maximum number of iteration reached'
            goto 110
        end if
    End do Main

    v=0.D0

    v(1:m*(m+1)/2)=fu(1:m*(m+1)/2)

    istop=1

!================ pour les bandes de confiance
!==== on ne retient que les para des splines

    call deriva(b,m,v,rl,k0,fctnames)
    if(rl.eq.-1.D9) then
        istop=4
        goto 110
    end if

    do i=1,(m*(m+3)/2)
        v1(i)=0.d0
    end do
!---- Choix du model

    select case(model)
        case(1,7)
            m1=m-nva-effet-indic_alpha !joint

        !cas supplementaire rajout? pour joint general
        case(5)
            m1=m-nva-effet-indic_alpha !joint general
        case(2)
            m1=m-nva-effet*2 !additive
        case(3)
            m1=m-nva-effet !nested
        case(4)
            m1=m-nva-effet !shared
    end select

    kkk=m1*(m1+1)/2

    do i=1,m1
        kkk=kkk+1
        do j=i,m1
            k = (((j-1)*j)/2) +i
            v1(k)=v(k)/(4.d0*b(i)*b(j))
        end do
        v1(kkk)=v1(kkk)+(v(kkk)/(4.d0*b(i)*b(i)*b(i)))
    end do

    ep=10.d-10
    call dsinvJ(v1,m1,ep,ier)

    if (ier.eq.-1) then
!         write(*,*)   'echec inversion matrice information'
        istop=3
    endif



    do i=1,m1
        do j=i,m1
            hess(i,j)=v1((j-1)*j/2+i)
        end do
    end do

    do i=2,m1
        do j=1,i-1
            hess(i,j)=hess(j,i)
        end do
    end do


    ep=10.d-10
    call dsinvJ(v,m,ep,ier)

    if (ier.eq.-1) then
        istop=3

!AD:
!        call dsinvj(v1,m1,ep,ier)
!        if (ier.eq.-1) then
!             write(*,*)'echec inversion matrice information
!     & prms fixes'
!            istop=31
!        else
!            DO k=1,m1*(m1+1)/2
!                v(k)=v1(k)
!            END DO
!        end if
! fin ajout amadou
    endif


    ep=10.d-10
    call deriva(b,m,vnonpen,rl,zero,fctnames)

    do i=1,m
        do j=i,m
            I_hess(i,j)=vnonpen((j-1)*j/2+i)
        end do
    end do

    do i=2,m
        do j=1,i-1
            I_hess(i,j)=I_hess(j,i)
        end do
    end do

!========================================================

!   H_hess est moins la hessienne inverse sur la vraisemblance penalisee
    do i=1,m
        do j=i,m
            H_hess(i,j)=v((j-1)*j/2+i)
        end do
    end do
!       write(*,*) 'H_hess(16,16) fin marq',H_hess(16,16),m,((j-1)*j/2+i)

    do i=2,m
        do j=1,i-1
            H_hess(i,j)=H_hess(j,i)
        end do
    end do

 !AD:
    if (typeof .ne. 0) then
        do i=1,m*(m+1)/2
            vvv(i)=v(i)
        end do
    end if

 110   continue

       return

       end subroutine marq98


    subroutine deriva(b,m,v,rl,k0,fctnames)
    use comonmultiv,only:model
    implicit none

    integer,intent(in)::m
    double precision,intent(inout)::rl
    double precision,dimension(3)::k0
    double precision,dimension(m),intent(in)::b
    double precision,dimension((m*(m+3)/2)),intent(out)::v
    double precision,dimension(m)::fcith
    integer ::i0,m1,ll,i,k,j,iun
    double precision::fctnames,thn,th,z,vl,th2,vaux
    external::fctnames

    select case(model)
    case(1,7)
        th=1.d-3 !joint
    case(5)
        th=1.d-5 !joint general
    case(2)
        th=5.d-3 !additive
    case(3)
        th=1.d-5 !nested
    case(4)
        th=1.d-5 !shared
    end select

    thn=-th
    th2=th*th
    z=0.d0
    i0=0
    iun =1

    rl=fctnames(b,m,iun,z,iun,z,k0)
        if(rl.eq.-1.d9) then
            rl=-1.d9
            goto 123
        end if
    do i=1,m
        fcith(i)=fctnames(b,m,i,th,i0,z,k0)
                if(fcith(i).eq.-1.d9) then
                    rl=-1.d9
                    goto 123
                end if
    end do

    k=0
    m1=m*(m+1)/2
    ll=m1

    do i=1,m
        ll=ll+1
        vaux=fctnames(b,m,i,thn,i0,z,k0)
                if(vaux.eq.-1.d9) then
                    rl=-1.d9
                    goto 123
                end if
        vl=(fcith(i)-vaux)/(2.d0*th)
        v(ll)=vl
        do j=1,i
            k=k+1
            v(k)=-(fctnames(b,m,i,th,j,th,k0)-fcith(j)-fcith(i)+rl)/th2
        end do
    end do

123   continue
    return

    end subroutine deriva


      subroutine searpas(vw,step,b,bh,m,delta,fim,k0,fctnames)
!
!  MINIMISATION UNIDIMENSIONNELLE
!
      implicit none

      integer,intent(in)::m
      double precision,dimension(m),intent(in)::b
      double precision,intent(inout)::vw
      double precision,dimension(m),intent(inout)::bh,delta
      double precision,intent(inout)::fim,step
      double precision::vlw,vlw1,vlw2,vlw3,vm,fi1,fi2,fi3
      double precision,dimension(3)::k0
      double precision::fctnames
      external::fctnames
      integer::i

       vlw1=dlog(vw)
       vlw2=vlw1+step
       call valfpa(vlw1,fi1,b,bh,m,delta,k0,fctnames)
       call valfpa(vlw2,fi2,b,bh,m,delta,k0,fctnames)

       if(fi2.ge.fi1) then
      vlw3=vlw2
      vlw2=vlw1
      fi3=fi2
      fi2=fi1
      step=-step

          vlw1=vlw2+step
          call valfpa(vlw1,fi1,b,bh,m,delta,k0,fctnames)
          if(fi1.gt.fi2) goto 50
       else
          vlw=vlw1
          vlw1=vlw2
          vlw2=vlw
          fim=fi1
          fi1=fi2
          fi2=fim
       end if

       do i=1,40
          vlw3=vlw2
          vlw2=vlw1
          fi3=fi2
          fi2=fi1

          vlw1=vlw2+step
          call valfpa(vlw1,fi1,b,bh,m,delta,k0,fctnames)
          if(fi1.gt.fi2) goto 50
          if(fi1.eq.fi2) then
             fim=fi2
             vm=vlw2
             goto 100
          end if
       end do
!
!  PHASE 2 APPROXIMATION PAR QUADRIQUE
!
50     continue
!
!  CALCUL MINIMUM QUADRIQUE
!
      vm=vlw2-step*(fi1-fi3)/(2.d0*(fi1-2.d0*fi2+fi3))
      call valfpa(vm,fim,b,bh,m,delta,k0,fctnames)
      if(fim.le.fi2) goto 100
      vm=vlw2
      fim=fi2
100   continue
      vw=dexp(vm)

      return

      end subroutine searpas


    subroutine valfpa(vw,fi,b,bk,m,delta,k0,fctnames)

    implicit none

    integer,intent(in)::m
    double precision,dimension(m),intent(in)::b,delta
    double precision,dimension(m),intent(out)::bk
    double precision,intent(out)::fi
    double precision::vw,fctnames,z
    double precision,dimension(3)::k0
    integer::i0,i
    external::fctnames

    z=0.d0
    i0=1
    do i=1,m
    bk(i)=b(i)+dexp(vw)*delta(i)
    end do
    fi=-fctnames(bk,m,i0,z,i0,z,k0)
    if(fi.eq.-1.D9) then
        goto 1
    end if
1       continue
    return

    end subroutine valfpa


!                   MARQ98JAlt - IJ
!-------------------------------------------------------------


    subroutine marq98jAlt(k0,b,m,ni,v,rl,ier,istop,effet,ca,cb,dd,fctnames,fctnamesindiv)
    
!
!  fu = matrice des derivees secondes et premieres
!
!  istop: raison de l'arret
!  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
!  2: nb max d'iterations atteints
!  4: Erreur
    !use residusM,only:indg
    use parameters
    !use comon,only: c,cdc,Hspl_hess,indic_eta,ndate,ndatedc,nst,nsujet,&
    !nt0,nt0dc,nt1,nt1dc,nva1,nva2,PEN_deri,t0,t0dc,t1,t1dc
    use comon,only:nva,model,I_hess,H_hess,hess,indic_ALPHA,&
    typeof,vvv,ng !indic_tronc

!add additive
    !use additiv,only:correl

    IMPLICIT NONE
!   variables globales
    integer,intent(in) :: m,effet
    integer,intent(inout)::ni,ier,istop
    double precision,dimension(m*(m+3)/2),intent(out)::v
    double precision,intent(out)::rl
    double precision,dimension(m),intent(inout)::b
    double precision,intent(out)::ca,cb,dd
    double precision,dimension(2)::k0
    double precision,dimension(2)::zero
!   variables locales
    integer::nql,ii,nfmax,idpos,ncount,id,jd,m1,j,i,ij,k,sub !iun
    character(len=100)::bar ! Pour la progressBar
    double precision,dimension(m*(m+3)/2)::fu,v1 !vnonpen
    double precision,dimension(m)::delta,b1,bh,vindiv ! last one IJ
    double precision::da,dm,ga,tr
    double precision::GHG,step,eps,vw,fi,maxt, &
    z,rl1,th,ep
    external::fctnames,fctnamesindiv
    double precision::fctnames,fctnamesindiv
!---------- ajout
    integer::kkk
    
    zero=0.d0
    id=0
    jd=0
    z=0.d0
    th=1.d-5
    eps=1.d-7!1.d-6
    nfmax=m*(m+1)/2
    ca=epsa+1.d0
    cb=epsb+1.d0
    rl1=-1.d+10
    ni=0
    istop=0
    da=0.01d0
    dm=5.d0
    nql=1
    m1=m*(m+1)/2
    ep=1.d-20
    
    !write(*,*) "started marq98jAlt"
    
    
    Main:Do       
        
        call derivaJ(b,m,v,rl,k0,fctnames)     
        
        !write(*,*) "took derivative of full"    
        
    rl1=rl
    
        !print *,"rl1 = ", rl1
        if(rl.eq.-1.D9) then
            istop=4
            goto 110
        endif
        !write(*,*)'iteration***',ni,'vrais',rl
        
        !Myriam ! Construction de la progressBar
        if (((effet == 2).and.(model == 1)).or.(model == 7)) then 
            if(maxiter < 100) then
                bar(1:3) = "0%|"    
                do k=1, maxiter
                    bar(3+k:3+k)="-"
                enddo
                bar(maxiter+5:maxiter+5) = "|"
                do k=maxiter+6, 100
                    bar(k:k+1)=" "
                enddo
                
                do k=1, ni
                    bar(3+k:3+k)="*"
                enddo
                call intpr(bar, -1, ni, 0)
                call intpr('Iteration:', -1, ni, 1)
            else
                bar(1:3) = "0%|"    
                do k=1, 70
                    bar(3+k:3+k)="-"
                enddo
                bar(73:100) = "|                          "
                
                if (modulo(ni,3) == 0) then 
                    do k=1, ni/3
                        bar(3+k:3+k)="*"
                    enddo
                    if (ni < 300) then 
                        call intpr(bar, -1, ni, 0)
                        call intpr('Iteration:', -1, ni, 1) 
                    else 
                        call intpr('Iteration:', -1, ni, 1)
                    endif
                endif
            endif            
        endif
    
    ! Fin de la progressBar
    
    !write(*,*) 'iteration ***', ni, 'vecteur b ', b
    !write(*,*) 'm', m, 'iun', iun, 'z', z, 'k0', k0

    dd = 0.d0
    fu = 0.d0

    do i=1,m
        do j=i,m
            ij=(j-1)*j/2+i
            fu(ij)=v(ij)
        end do
    end do

        call dsinvj(fu,m,ep,ier)
        
        !write(*,*) "inversion line"   

    if (ier.eq.-1) then ! hessienne non inversible
        !print*,"here"
        dd=epsd+1.d0
    else
        GHG = 0.d0
        do i=1,m
            do j=1,m
                if(j.ge.i) then
                    ij=(j-1)*j/2+i
                else
                    ij=(i-1)*i/2+j
                end if
                GHG = GHG + v(m1+i)*fu(ij)*V(m1+j)
            end do
        end do
        dd=GHG/dble(m)
    end if

!write(*,*)'epsa',ca,'epsb',cb,'epsd',dd
!   write(*,*)'aaOptim: b(knots)',b(1:12)
! write(*,*) 'aaOptim: m',m, nva, nva1, nva2
  !! write(*,*)'aaOptim: b(random):', b((m-nva-3):(m-nva))
   !!write(*,*)'aaOptim: b(covar)',b((m-nva+1):m)

!     print*,ca,cb,dd
!     print*,b
!    print*,"-------------------------"
!     if (ni.eq.50) then
!         print*,"====================> FAUSSE CONVERGENCE"
!         ca=1.d-5
!         cb=1.d-5
!         dd=1.d-5
!     endif

 !write(*,*) ca
  !write(*,*) cb
   !write(*,*) dd

    if(ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) exit main

        tr=0.d0
        do i=1,m
            ii=i*(i+1)/2
            tr=tr+dabs(v(ii))
        end do
        tr=tr/dble(m)

        ncount=0
        ga=0.01d0

 400    do i=1,nfmax+m
           fu(i)=v(i)
        end do
        
        do i=1,m
            ii=i*(i+1)/2
            if (v(ii).ne.0) then
                fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
            else
                fu(ii)=da*ga*tr
            endif
        end do
        
        call dcholej(fu,m,nql,idpos)

        if (idpos.ne.0) then
            ncount=ncount+1
            if (ncount.le.3.or.ga.ge.1.d0) then
                da=da*dm
            else
                ga=ga*dm
                if (ga.gt.1.d0) ga=1.d0
            endif

            goto 400

        else
        
            do i=1,m
            
                delta(i)=fu(nfmax+i)
                b1(i)=b(i)+delta(i)
            end do
            
            rl=fctnames(b1,m,id,z,jd,z,k0)
            
            if(rl.eq.-1.D9) then
                istop=4
                goto 110
            end if
            if (rl1.lt.rl) then
                if(da.lt.eps) then
                    da=eps
                else
                    da=da/(dm+2.d0)
                endif
                goto 800
            endif
        endif
!      write(6,*) 'loglikelihood not improved '
        call dmaxt(maxt,delta,m)
                
        if(maxt.eq.0.D0) then
            vw=th
        else
            !call dmaxt(maxt,delta,m)
            vw=th/maxt

        endif
        step=dlog(1.5d0)
        call searpasj(vw,step,b,bh,m,delta,fi,k0,fctnames)
        
        rl=-fi
        if(rl.eq.-1.D9) then
            istop=4
            goto 110
        end if
        do i=1,m
            delta(i)=vw*delta(i)
        end do
        da=(dm-3.d0)*da

 800     cb=dabs(rl1-rl)

        ca=0.d0
        do i=1,m
            ca=ca+delta(i)*delta(i)
        end do
        !write(6,*) 'ca =',ca,' cb =',cb,' dd =',dd

        !print*,"rl1=",rl1
        !print*,"ca =",ca," cb =",cb," dd =",dd
        do i=1,m
            b(i)=b(i)+delta(i)
        end do
        
        !write(*,*) ni 
        
        ni=ni+1
        if (ni.ge.maxiter) then
            istop=2
!            write(6,*) 'maximum number of iteration reached'
            goto 110
        end if
    End do Main
    
    v=0.D0

    v(1:m*(m+1)/2)=fu(1:m*(m+1)/2)

    istop=1

!================ pour les bandes de confiance
!==== on ne retient que les para des splines
        
    call derivaJ(b,m,v,rl,k0,fctnames)
     !print*,"dans aaOptim, appel de derivaJ  with GAP vec funcpa... "
    if(rl.eq.-1.D9) then
        istop=4
        goto 110
    end if

    do i=1,(m*(m+3)/2)
        v1(i)=0.d0
    end do
!---- Choix du model

    select case(model)
        case(1,7)
            m1=m-nva-effet-indic_alpha !joint
        case(6)  !Uni 
            m1 = m - nva - effet*2 - indic_alpha*2
        !cas supplementaire rajout? pour joint general
        case(5)
            m1=m-nva-effet-indic_alpha !joint general (du moment que indic_eta=indic_alpha=1)
        case(2)
            m1=m-nva-effet*2 !additive
        case(3)
            m1=m-nva-effet !nested
        case(4)
            m1=m-nva-effet !shared
    end select

    kkk=m1*(m1+1)/2

    do i=1,m1
        kkk=kkk+1
        do j=i,m1
            k = (((j-1)*j)/2) +i
            v1(k)=v(k)/(4.d0*b(i)*b(j))
        end do
        v1(kkk)=v1(kkk)+(v(kkk)/(4.d0*b(i)*b(i)*b(i)))
    end do
    
    ep=10.d-10
    call dsinvJ(v1,m1,ep,ier)

    if (ier.eq.-1) then
!         write(*,*)   'echec inversion matrice information'
        istop=3
    endif



    do i=1,m1
        do j=i,m1
            hess(i,j)=v1((j-1)*j/2+i)
        end do
    end do

    do i=2,m1
        do j=1,i-1
            hess(i,j)=hess(j,i)
        end do
    end do

    ep=10.d-10
    call dsinvJ(v,m,ep,ier)

    if (ier.eq.-1) then
        istop=3

!AD:
!        call dsinvj(v1,m1,ep,ier)
!        if (ier.eq.-1) then
!             write(*,*)'echec inversion matrice information
!     & prms fixes'
!            istop=31
!        else
!            DO k=1,m1*(m1+1)/2
!                v(k)=v1(k)
!            END DO
!        end if
! fin ajout amadou
    endif
    
    ep=10.d-10
    
    !!!! IJ edit: changing I_hess
    !call derivaJ(b,m,vnonpen,rl,zero,fctnames)

    !do i=1,m
    !    do j=i,m
    !        I_hess(i,j)=vnonpen((j-1)*j/2+i)
    !    end do
    !end do
    
    !write(*,*) "starting new portion"
    
    ! IJ: initialize I_hess (for summing each element)
    do i=1,m
        do j=1,m
            I_hess(i,j)=0.d0
        end do
    end do
    
    !write(*,*) "initialized I_hess"
    
    ! IJ: compute matrix summing the outer products of the score vectors
    ! for each individual (group)
    do sub=1,ng 
        !write(*,*) sub
        ! compute ONLY first derivatives of the individual contribution to the likelihood for subject sub 
        call derivaJindiv(b,m,vindiv,rl,zero,sub,fctnamesindiv)
        !write(*,*) "indiv deriv computed"
        ! add each term of the outer product to the appropriate slot (upper triangular matrix)
        do i=1,m
            do j=i,m
                I_hess(i,j)=I_hess(i,j)+(vindiv(i)*vindiv(j))
            end do
        end do    
        !write(*,*) I_hess(1,1)
    end do        
    
    !write(*,*) "finished upper triangle of matrix"

    ! fill in lower triangle of I_hess
    do i=2,m
        do j=1,i-1
            I_hess(i,j)=I_hess(j,i)
        end do
    end do
    
    !write(*,*) "reflected matrix"

!========================================================

!   H_hess est moins la hessienne inverse sur la vraisemblance penalisee
    do i=1,m
        do j=i,m
            H_hess(i,j)=v((j-1)*j/2+i)
        end do
    end do
!       write(*,*) 'H_hess(16,16) fin marq',H_hess(16,16),m,((j-1)*j/2+i)

    do i=2,m
        do j=1,i-1
            H_hess(i,j)=H_hess(j,i)
        end do
    end do

 !AD:
    if (typeof .ne. 0) then
        do i=1,m*(m+1)/2
            vvv(i)=v(i)
        end do
    end if

 110   continue
 
       return
       end subroutine marq98jAlt
 
 
     
!------------------------------------------------------------
!                          DERIVA INDIV - IJ
!------------------------------------------------------------

        ! first derivatives only
    subroutine derivaJindiv(b,m,vindiv,rl,k0,index,fctnamesindiv)
    !use comon,only:nst,k0T
    use comon,only:model !nstRec, indic_tronc
    implicit none

    integer,intent(in)::m,index
    double precision,intent(inout)::rl
    double precision,dimension(2)::k0
    double precision,dimension(m),intent(in)::b
    double precision,dimension(m),intent(inout)::vindiv
    double precision,dimension(m)::fcith
    integer ::i0,i,iun !ll,m1,k,j
    double precision::fctnamesindiv,thn,th,z,vl,th2,vaux
    external::fctnamesindiv
    logical::endDeriva

    endDeriva=.false.
    
    select case(model)
    case(1,7)
        th=1.d-3 !joint
    case(5)
        th=1.d-5 !joint general !cas supplementaire rajout?
    case(2)
        th=5.d-3 !additive
    case(3)
        th=1.d-5 !nested
    case(4)
        th=1.d-5 !shared
    end select

    thn=-th
    th2=th*th
    z=0.d0
    i0=0
    iun =1

        !print *,"k0 = ",k0
        !print *,"b = ", b
        !size(b) =32 , m =32
        !print *,"z = ", z            ! =0.000000000
        !print *,"iun", iun           ! =1
        !print *,"rl = ",rl           ! =0.000000000
    
    rl=fctnamesindiv(b,m,iun,z,iun,z,k0,index)
    
    if(rl.ne.-1.d9) then
        do i=1,m
            if((fcith(i).ne.-1.d9).and.(endDeriva.eqv..false.)) then
                fcith(i)=fctnamesindiv(b,m,i,th,i0,z,k0,index)
            else
                rl=-1.d9
                endDeriva=.true.
            end if
        end do
                
        if (endDeriva.eqv..false.) then
            !k=0
            !m1=m*(m+1)/2
            !ll=m1
            do i=1,m
                !print*, "deuxieme boucle DERIVAj, tour",i," :"
                !ll=ll+1
                !write(*,*) 'b', b
                vaux=fctnamesindiv(b,m,i,thn,i0,z,k0,index)
                !write(*,*) 'vaux', vaux
                if(vaux.eq.-1.d9) then
                    rl=-1.d9
                    endDeriva=.true.
                end if
                if (endDeriva.eqv..false.) then
                    vl=(fcith(i)-vaux)/(2.d0*th)
                    vindiv(i)=vl
                    !do j=1,i
                    !   k=k+1
                    !   vindiv(k)=-(fctnamesindiv(b,m,i,th,j,th,k0,index)-fcith(j)-fcith(i)+rl)/th2
                !   !     write(*,*)k,fctnames(b,m,i,th,j,th,k0),b,m,i,th,j,k0
                    !end do
                end if
            end do
        end if
    end if
        
    !stop
    return

    end subroutine derivaJindiv

    

    
    
    end module optim


