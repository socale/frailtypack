


!========================          funcpaMultivSplines         ====================
    double precision function funcpaMultivSplines(b,np,id,thi,jd,thj,k0)

    use taillesmultiv
    use comonmultiv
    use residusMmultiv

    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :

    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(3),intent(in)::k0
    double precision,intent(in)::thi,thj
    integer::n,i,j,k,vj,ig,choix
    integer,dimension(ngmax)::cpt,cptmeta
    double precision::pe1,pe2,pe3,sum,inv,som1,som2,som3,res,vet,vet2,vet3,h1,pi
    double precision,dimension(-2:(nzloco+2)):: the1
    double precision,dimension(-2:(nzdc+2)):: the2
    double precision,dimension(-2:(nzmeta+2)):: the3
    double precision,dimension(np)::bh
    !double precision,dimension(0:ndatemax)::imdc
    double precision,dimension(ngmax)::res2,res1dc,res2dc &
    ,res3dc,integrale1,integrale2,integrale3,res2meta
!AD: for death,change dimension
    double precision,dimension(ndatemax)::dut1
    double precision,dimension(ndatemaxdc)::dut2
    double precision,dimension(ndatemeta)::dut3
!AD:end
    double precision,dimension(0:ndatemax)::ut1
    double precision,dimension(0:ndatemaxdc)::ut2
    double precision,dimension(0:ndatemeta)::ut3
    double precision::int

    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    pi=3.141592653589793d0

!    n = (np-nva-effet-indic_ALPHA)/3
    the1=0.d0
    the2=0.d0
    the3=0.d0

    do i=1,nzloco+2
        the1(i-3)=(bh(i))*(bh(i))
    end do

    do i=1,nzdc+2
        the2(i-3)=(bh(nzloco+2+i))*(bh(nzloco+2+i))
    end do

    do i=1,nzmeta+2
        the3(i-3)=(bh(nzloco+nzdc+4+i))*(bh(nzloco+nzdc+4+i))
    end do

    if(effet.eq.1) then
        theta = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
!        alpha = bh(np-nva)
        alpha =bh(np-nva-indic_rho) !rho
        alpha1=bh(np-nva-indic_a1)
        alpha2=bh(np-nva-indic_a2)
        eta = bh(np-nva-indic_eta)*bh(np-nva-indic_eta)
    endif


!----------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!AD:modify
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
    dut2(1) = (the2(-2)*4.d0/(zidc(2)-zidc(1)))
    dut3(1) = (the3(-2)*4.d0/(zimeta(2)-zimeta(1)))

    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
    ut2(1) = the2(-2)*dut2(1)*0.25d0*(zidc(1)-zidc(-2))
    ut3(1) = the3(-2)*dut3(1)*0.25d0*(zimeta(1)-zimeta(-2))

    ut1(0) = 0.d0
    ut2(0) = 0.d0
    ut3(0) = 0.d0
!//// NEW AMADOU vvv :

!--- strate1
    som1 = 0.d0
    vj = 0
    do i=2,ndate-1
        do k = 2,nzloco!n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som1 = som1 + the1(j-4)
                vj  = j
                endif
            endif
        end do

        ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
        +(the1(j-1)*im1(i))+(the1(j)*im(i))

        dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
        +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

    end do

!--- strate2
    vj = 0
    som2 = 0.d0
    do i=2,ndatedc-1
        do k = 2,nzdc!n-2
            if (((datedc(i)).ge.(zidc(k-1))).and.(datedc(i).lt.zidc(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som2 = som2 + the2(j-4)
                vj  = j
                endif
            endif
        end do

        ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
        +(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))

        dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
        +(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))

    end do

!--- strate3
    som3 = 0.d0
    vj = 0
    do i=2,ndatemeta-1
        do k = 2,nzmeta
            if (((datemeta(i)).ge.(zimeta(k-1))).and.(datemeta(i).lt.zimeta(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som3 = som3 + the3(j-4)
                vj  = j
                endif
            endif
        end do

        ut3(i) = som3 +(the3(j-3)*im3meta(i))+(the3(j-2)*im2meta(i)) &
        +(the3(j-1)*im1meta(i))+(the3(j)*immeta(i))

        dut3(i) = (the3(j-3)*mm3meta(i))+(the3(j-2)*mm2meta(i)) &
        +(the3(j-1)*mm1meta(i))+(the3(j)*mmmeta(i))

    end do
!-------------fin strate3

    i = nzloco
    h1 = (zi(i)-zi(i-1))
    ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)

    i = nzdc
    h1 = (zidc(i)-zidc(i-1))
    ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
    dut2(ndatedc) = (4.d0*the2(i-1)/h1)

    i = nzmeta
     h1 = (zimeta(i)-zimeta(i-1))
    ut3(ndatemeta)=som3+the3(i-4)+the3(i-3)+the3(i-2)+the3(i-1)
    dut3(ndatemeta) = (4.d0*the3(i-1)/h1)
!//// fin NEW AMADOU-vvv
!AD:end
!-------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------

    do k=1,ng
        res1(k) = 0.d0
        res2(k) = 0.d0
        res3(k) = 0.d0
        res1dc(k) = 0.d0
        res2dc(k) = 0.d0
        res3dc(k) = 0.d0
        res1meta(k) = 0.d0
        res2meta(k) = 0.d0
        res3meta(k) = 0.d0
        cpt(k) = 0
        cptmeta(k) = 0
        integrale1(k) = 0.d0
        integrale2(k) = 0.d0
        integrale3(k) = 0.d0
        aux1(k)=0.d0
        aux2(k)=0.d0
    end do

!*********************************************
!-----avec un effet aleatoire dans le modele
!*********************************************

    inv = 1.d0/theta

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc


    do i=1,nsujet
        cpt(g(i))=cpt(g(i))+1
        if(nva1.gt.0)then
            vet = 0.d0
            do j=1,nva1
                vet =vet + bh(np-nva+j)*dble(ve(i,j))
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif

        if((c(i).eq.1))then
            res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
        endif

        if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
            funcpaMultivSplines=-1.d9
            goto 123
        end if

!     nouvelle version
        res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            funcpaMultivSplines=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature:
        res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            funcpaMultivSplines=-1.d9
            goto 123
        end if
    end do

!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces
!ccccccccccccccccccccccccccccccccccccccccc

    do k=1,ng
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 =vet2 + bh(np-nva3-nva2+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2=1.d0
        endif

        if(cdc(k).eq.1)then
            res2dc(k) = dlog(dut2(nt1dc(k))*vet2)

            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                funcpaMultivSplines=-1.d9
                goto 123
            end if
        endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
        aux1(k)=ut2(nt1dc(k))*vet2
!        aux2(k)=aux2(k)+ut2(nt0(k))*vet2 !vraie troncature
        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            funcpaMultivSplines=-1.d9
            goto 123
        end if
!        if ((aux2(k).ne.aux2(k)).or.(abs(aux2(k)).ge. 1.d30)) then
!            funcpaMultivSplines=-1.d9
!            goto 123
!        end if
    end do


    do i=1,nsujetmeta
        cptmeta(gmeta(i))=cptmeta(gmeta(i))+1
        if(nva3.gt.0)then
            vet3 = 0.d0
            do j=1,nva3
                vet3 =vet3 + bh(np-nva3+j)*dble(vemeta(i,j))
            end do
            vet3 = dexp(vet3)
        else
            vet3=1.d0
        endif

        if((cmeta(i).eq.1))then
            res2meta(gmeta(i)) = res2meta(gmeta(i))+dlog(dut3(nt1meta(i))*vet3)
        endif

! Fonction de risque cumulée de recidive au tepms T_ij
        res1meta(gmeta(i)) = res1meta(gmeta(i))+ut3(nt1meta(i))*vet3
! Fonction de risque de recidive au tepms T_i(j-1)
        res3meta(gmeta(i)) = res3meta(gmeta(i))+ut3(nt0meta(i))*vet3

    end do

!**************INTEGRALES ****************************
    do ig=1,ng
        auxig=ig
        call gausshermiteBIS2011(int,30)
        integrale3(ig) = int !moins bon
    end do

!************* FIN INTEGRALES **************************

    res = 0.d0
    do k=1,ng
        sum=0.d0
!        if(cpt(k).gt.0.or. cptmeta(k).gt.0)then

        res= res + res2(k)+res2meta(k)+res2dc(k)+dlog(integrale3(k)) &
        - dlog((2.d0)*pi*sqrt(theta*eta* &
        (1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2)))

        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpaMultivSplines=-1.d9
            goto 123
        end if

    end do

!---------- calcul de la penalisation -------------------

    pe1 = 0.d0
    pe2 = 0.d0
    pe3 = 0.d0

    do i=1,nzloco-1
        pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
        *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
        the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
        m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
        the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
        m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
        *the1(i)*m1m(i))
    end do
    do i=1,nzdc-1
        pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3b(i))+(the2(i-2) &
        *the2(i-2)*m2m2b(i))+(the2(i-1)*the2(i-1)*m1m1b(i))+( &
        the2(i)*the2(i)*mmmb(i))+(2.d0*the2(i-3)*the2(i-2)* &
        m3m2b(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1b(i))+(2.d0* &
        the2(i-3)*the2(i)*m3mb(i))+(2.d0*the2(i-2)*the2(i-1)* &
        m2m1b(i))+(2.d0*the2(i-2)*the2(i)*m2mb(i))+(2.d0*the2(i-1) &
        *the2(i)*m1mb(i))
    end do
    do i=1,nzmeta-1
        pe3 = pe3+(the3(i-3)*the3(i-3)*m3m3c(i))+(the3(i-2) &
        *the3(i-2)*m2m2c(i))+(the3(i-1)*the3(i-1)*m1m1c(i))+(  &
        the3(i)*the3(i)*mmmc(i))+(2.d0*the3(i-3)*the3(i-2)*  &
        m3m2c(i))+(2.d0*the3(i-3)*the3(i-1)*m3m1c(i))+(2.d0*  &
        the3(i-3)*the3(i)*m3mc(i))+(2.d0*the3(i-2)*the3(i-1)*  &
        m2m1c(i))+(2.d0*the3(i-2)*the3(i)*m2mc(i))+(2.d0*the3(i-1)  &
        *the3(i)*m1mc(i))
    end do

    pe = k0(1)*pe1 + k0(2)*pe2 + k0(3)*pe3

    resnonpen = res
    res = res - pe

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaMultivSplines=-1.d9
        do k=1,ng
            Rrec(k)=0.d0
            Nrec(k)=0

            Rdc(k)=0.d0
            Ndc(k)=0

            Rrec2(k)=0.d0
            Nrec2(k)=0
        end do
        goto 123
    else
        funcpaMultivSplines = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)

            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)

            Rrec2(k)=res1meta(k)
            Nrec2(k)=nigmeta(k)
        end do
    end if

!Ad:
123     continue

    return

    end function funcpaMultivSplines



!==========================  DISTANCE   =================================

    subroutine distanceJ_splines(nzloco,nzdc,nzmeta,b,mt1,mt2,mt3,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out,x3Out,lam3Out,su3Out)

    use taillesmultiv
    !use comonmultiv,only:c,cdc,date,datedc,hess,Hspl_hess,I_hess,ndate,ndatedc,&
    !nt0,nt1,nt0dc,nt1dc,nsujet,nva,nva1,nva2,t0,t1,t0dc,t1dc,PEN_deri
    use comonmultiv,only:zi,zidc,zimeta,H_hess !,nst &    

    Implicit none

    integer,intent(in):: nzloco,nzdc,nzmeta,mt1,mt2,mt3
    double precision ,dimension(npmax),intent(in):: b
    integer::i,j,n,k,l      
    double precision::x1,x2,h,hdc,hmeta,su,bsup,binf,lam,lbinf,lbsup
    double precision ,dimension(npmax,npmax)::hes1,hes2,hes3
    double precision ,dimension(-2:npmax):: the1,the2,the3     
    double precision,dimension(mt1)::x1Out
    double precision,dimension(mt2)::x2Out
    double precision,dimension(mt3)::x3Out
    double precision,dimension(mt1,3):: lamOut,suOut
    double precision,dimension(mt2,3):: lam2Out,su2Out
    double precision,dimension(mt3,3):: lam3Out,su3Out    

    n  = nzloco+2
    !loco
    do i=1,nzloco+2
        do j=1,nzloco+2
            hes1(i,j)=h_Hess(i,j)
        end do
    end do

    !dc
    k = 0
    do i=nzloco+3,nzloco+2+nzdc+2
        k = k + 1
        l = 0
        do j=nzloco+3,nzloco+2+nzdc+2
            l = l + 1
            hes2(k,l)=H_hess(i,j)
        end do
    end do

    !meta
    k = 0
    do i=nzloco+nzdc+5,nzloco+2+nzdc+2+nzmeta+2
        k = k + 1
        l = 0
        do j=nzloco+nzdc+5,nzloco+2+nzdc+2+nzmeta+2
            l = l + 1
            hes2(k,l)=H_hess(i,j)
        end do
    end do


    do i=1,nzloco+2
        the1(i-3)=(b(i))*(b(i))
    end do

    do i=1,nzdc+2
        j = nzloco+2+i
        the2(i-3)=(b(j))*(b(j))
    end do

    do i=1,nzmeta+2
        j = nzloco+2+nzdc+2+i
        the3(i-3)=(b(j))*(b(j))
    end do

    h = (zi(n)-zi(1))*0.01d0
    hdc = (zidc(n)-zidc(1))*0.01d0
    hmeta = (zimeta(n)-zimeta(1))*0.01d0

! Recurrent
    x1 = zi(1)
    x2 = zi(1)
    do i=1,mt1
        if(i .ne.1)then
            x1 = x1 + h
        end if
        call cosp(x1,the1,nzloco+2,hes1,zi,binf,su,bsup,lbinf,lam,lbsup)

        if(bsup.lt.0.d0)then
            bsup = 0.d0
        endif
        if(binf.gt.1.d0)then
            binf = 1.d0
        endif
        if(lbinf.lt.0.d0)then
            lbinf = 0.d0
        endif 
!!!   Replaced by next sentences and add new ones JRG January 05

        x1Out(i)=x1
        lamOut(i,1)=lam
        lamOut(i,2)=lbinf
        lamOut(i,3)=lbsup
        suOut(i,1)=su
        suOut(i,2)=bsup
        suOut(i,3)=binf
    end do

! Death
    x1 = zidc(1)
    x2 = zidc(1)
    do i=1,mt2
!!!   Replaced by next sentences and add new ones JRG January 05
        if(i.ne.1)then
            x2 = x2 + hdc
        endif
        call cosp(x2,the2,nzdc+2,hes2,zidc,binf,su,bsup,lbinf,lam,lbsup)
        if(bsup.lt.0.d0)then
            bsup = 0.d0
        endif
        if(binf.gt.1.d0)then
            binf = 1.d0 
        endif
        if(lbinf.lt.0.d0)then
            lbinf = 0.d0
        endif 

        x2Out(i)=x2
        lam2Out(i,1)=lam
        lam2Out(i,2)=lbinf
        lam2Out(i,3)=lbsup
        su2Out(i,1)=su
        su2Out(i,2)=binf
        su2Out(i,3)=bsup 
    end do

! Meta
    x1 = zimeta(1)
    x2 = zimeta(1)
    do i=1,mt3
        if(i .ne.1)then
            x1 = x1 + hmeta 
        end if
        call cosp(x1,the3,nzmeta+2,hes3,zimeta,binf,su,bsup,lbinf,lam,lbsup)

        if(bsup.lt.0.d0)then
            bsup = 0.d0
        endif
        if(binf.gt.1.d0)then
            binf = 1.d0
        endif
        if(lbinf.lt.0.d0)then
            lbinf = 0.d0
        endif
!!!   Replaced by next sentences and add new ones JRG January 05

        x3Out(i)=x1
        lam3Out(i,1)=lam
        lam3Out(i,2)=lbinf
        lam3Out(i,3)=lbsup
        su3Out(i,1)=su
        su3Out(i,2)=bsup
        su3Out(i,3)=binf
    end do
    return

    end subroutine distanceJ_splines

