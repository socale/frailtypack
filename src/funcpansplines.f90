
!========================          FUNCPA  NESTED_SPLINES       ====================


    double precision function funcpansplines(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:nz2,t0,t1,
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m &
    ,m2m1,m2m,m1m,mm3,mm2,mm1,mm,im3,im2,im1,im &
    ,date,zi,c,nt0,nt1,nsujet,nva,ndate,nst &
    ,stra,pe,effet,nz1,ve &
    ,g,nig,indictronq,AG,auxig,alpha,eta,resnonpen,nb_gl
    !use commun,only:nssgexact,
    use commun,only:ngexact,mij,mid,ssg,aux1,aux2,mij_ind
    use residusM

    Implicit none

    integer::nb,n,np,id,jd,i,j,k,vj,cptg,l,ig,ip,issg,choix
    integer,dimension(ngexact)::cpt
    double precision::thi,thj,pe1,pe2,dnb,sum,theta,inv &
    ,som1,som2,res,vet,h1,int,gammlnN
    double precision,dimension(-2:np)::the1,the2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngexact)::res1,res2,res3 &
    ,integrale1,integrale2,integrale3,sum1
    double precision,dimension(2)::k0
    double precision,dimension(ndatemax)::dut1,dut2
    double precision,dimension(0:ndatemax)::ut1,ut2

    dut1=0.d0
    dut2=0.d0
    ut1=0.d0
    ut2=0.d0

    bh=b
    j=0
    theta=0.d0
    res=0.d0

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    n = nz1+2

    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i 
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do

    if(effet.eq.1) then
        theta = (bh(np-nva)*bh(np-nva)) ! variance effet groupe
    endif

    if(effet.eq.2) then
        alpha = (bh(np-nva-1)*bh(np-nva-1)) ! variance effet groupe
        eta = (bh(np-nva)*bh(np-nva))  ! variance effet sous groupe
    endif

    vj = 0
    som1 = 0.d0
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
    ut1(0) = 0.d0
    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))

    if (nst.eq.2) then
        som2 = 0.d0
        dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
        ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
        ut2(0) = 0.d0
    endif

    do i=2,ndate-1
        do k = 2,n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                    som1 = som1+the1(j-4)
                    som2 = som2+the2(j-4)
                    vj  = j
                endif
            endif
        end do

        ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
        +(the1(j-1)*im1(i))+(the1(j)*im(i))
        dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
        +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

        if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3(i))+(the2(j-2)*im2(i)) &
            +(the2(j-1)*im1(i))+(the2(j)*im(i))
            dut2(i) = (the2(j-3)*mm3(i))+(the2(j-2)*mm2(i)) &
            +(the2(j-1)*mm1(i))+(the2(j)*mm(i))
        endif
    end do

    i = n-2
    h1 = (zi(i)-zi(i-1))
    ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)
    if(nst.eq.2)then
        ut2(ndate)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
        dut2(ndate) = (4.d0*the2(i-1)/h1)
    endif

!--------------------------------------------------------
!----------calcul de la vraisemblance ------------------
!---------------------------------------------------------

    res1 = 0.d0
    res2 = 0.d0
    cpt = 0

!*******************************************
!----- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then
        do i=1,nsujetmax
            cpt(g(i))=cpt(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*ve(i,j)
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if((c(i).eq.1).and.(stra(i).eq.1))then
                res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
            endif

            if((c(i).eq.1).and.(stra(i).eq.2))then
                res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
                res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet-ut1(nt0(i))*vet
            endif

            if(stra(i).eq.2)then
                res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet-ut2(nt0(i))*vet
            endif

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if
        end do

        res = 0.d0
        cptg = 0
! k indice les groupes
        do k=1,ngexact
            if(cpt(k).gt.0)then !nb de sujets dans un gpe=nig()
                res = res-res1(k)+ res2(k)
                cptg = cptg + 1
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpansplines=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle effet=0

!**************************************************
!-----avec un seul  effet aleatoire dans le modele
!**************************************************

    if (effet.eq.1) then
!    write(*,*)'AVEC 1 EFFET ALEATOIRE'

        inv = 1.d0/theta
        cpt=0
        res1=0.d0
        res2=0.d0
        res3=0.d0

        do i=1,nsujetmax

            cpt(g(i))=cpt(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if((c(i).eq.1).and.(stra(i).eq.1))then
                res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
            endif

            if((c(i).eq.1).and.(stra(i).eq.2))then
                res2(g(i)) = res2(g(i))+dlog(dut2(nt1(i))*vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
! nouvelle version
                res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
            endif

            if(stra(i).eq.2)then
! nouvelle version
                res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
            endif

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
            endif

            if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ut2(nt0(i))*vet
            endif

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if
        end do 

        res = 0.d0
        cptg = 0
        mid = 0
!     gam2 = gamma(inv)
! k indice les groupes

        do k=1,nsujet
            if(c(k).eq.1)then
                mid(g(k))=mid(g(k))+1
            endif
        end do

        do k=1,ngexact!exact
            sum=0.d0
            if(cpt(k).gt.0)then
                nb = mid(k)!nb de deces par groupe
                dnb = dble(nb)

                if (dnb.gt.1.d0) then
                    do l=1,nb
                        sum=sum+dlog(1.d0+theta*dble(nb-l))
                    end do
                endif
                if(theta.gt.(1.d-5)) then
!cccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        + res2(k) + sum  
!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res= res-(inv+dnb)*dlog(theta*(res1(k))+1.d0) &
                        +(inv)*dlog(theta*res3(k)+1.d0) &
                        + res2(k) + sum
                    endif

                    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        funcpansplines=-1.d9
                        goto 123
                    end if
                else
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                        +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0) &
                        +res2(k)+sum
!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-dnb*dlog(theta*res1(k)+1.d0) &
                        -res1(k)*(1.d0-theta*res1(k)/2.d0 &
                        +theta*theta*res1(k)*res1(k)/3.d0) &
                        +res2(k)+sum &
                        +res3(k)*(1.d0-theta*res3(k)/2.d0 &
                        +theta*theta*res3(k)*res3(k)/3.d0)
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpansplines=-1.d9
                    goto 123
                end if
            endif

        end do

        if (indic_cumul==1) then
            do i= 1,ngmax
                do j=1,n_ssgbygrp(i)
                    cumulhaz1(i,j) = res1((i-1)*n_ssgbygrp(i)+j)
                    cumulhaz0(i,j) = res3((i-1)*n_ssgbygrp(i)+j)
                end do
            end do
        end if
    endif !fin boucle effet=1

!************************************************
!-----avec deux effets aleatoires dans le modele
!************************************************

    if (effet.eq.2) then

        mid=0
        mij=0
        mij_ind = 0
        res1=0.d0
        res2=0.d0
        aux1=0.d0
        aux2=0.d0
        integrale1=0.d0
        integrale2=0.d0
        integrale3=0.d0

!===== MODIFICATION DE LA VRAISEMBLANCE POUR LE NESTED FRAILTY MODEL

        do k=1,nsujet
		
            if(c(k).eq.1)then
                mid(g(k))=mid(g(k))+1
                mij(g(k),ssg(k,g(k)))=mij(g(k),ssg(k,g(k)))+1
                mij_ind(g(k)) = mij_ind(g(k))+1
!nb de dc ds ss gpe ssg(k)
            endif
        end do
        do k=1,nsujet
            if(nva.gt.0)then
                vet = 0.d0
                do ip=1,nva
                    vet =vet + bh(np-nva+ip)*dble(ve(k,ip))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if((c(k).eq.1).and.(stra(k).eq.1))then
                res2(g(k)) = res2(g(k))+dlog(dut1(nt1(k))*vet)
            endif
            if((c(k).eq.1).and.(stra(k).eq.2))then
                res2(g(k)) = res2(g(k))+dlog(dut2(nt1(k))*vet)
            endif

            if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if

            if(stra(k).eq.1)then
                aux1(g(k),ssg(k,g(k)))=aux1(g(k),ssg(k,g(k)))+ut1(nt1(k))*vet
                aux2(g(k),ssg(k,g(k)))=aux2(g(k),ssg(k,g(k)))+ut1(nt0(k))*vet
            endif
            if(stra(k).eq.2)then
                aux1(g(k),ssg(k,g(k)))=aux1(g(k),ssg(k,g(k)))+ut2(nt1(k))*vet
                aux2(g(k),ssg(k,g(k)))=aux2(g(k),ssg(k,g(k)))+ut2(nt0(k))*vet
            endif

            if ((aux1(g(k),ssg(k,g(k))).ne.aux1(g(k),ssg(k,g(k)))).or.(abs(aux1(g(k),ssg(k,g(k)))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if

            if ((aux2(g(k),ssg(k,g(k))).ne.aux2(g(k),ssg(k,g(k)))).or.(abs(aux2(g(k),ssg(k,g(k)))).ge. 1.d30)) then
                funcpansplines=-1.d9
                goto 123
            end if
        end do

!================== calcul des integrales par Gauss LAGUERRE
!     memes points et poids dans chq groupe
        do ig=1,ngexact
            auxig=ig
            choix=1
            call gaulagN(int,choix,nb_gl)
            integrale1(auxig)=int
!     integrale sur la troncature:
            if(indictronq.eq.1)then
                if(AG.eq.1)then !andersen gill
                    choix=3
                    call gaulagN(int,choix,nb_gl)
                    integrale3(auxig)=int
                else !troncature classique
                    choix=2
                    call gaulagN(int,choix,nb_gl)
                    integrale2(auxig)=int
                endif
            endif
        end do

!======================================================================

        do ig=1,ngexact
            sum1(ig)=0.d0
            do issg=1,n_ssgbygrp(ig) !nssgbyg !!! NON ICI NSSGBYG
                if(mij(ig,issg).gt.1) then
                    do l=1,mij(ig,issg)
                sum1(ig)=sum1(ig)+dlog(1.d0+eta*dble(mij(ig,issg)-l))
                    end do
                endif
            end do
        end do

        res = 0.d0

        do k=1,ngexact  
            if(nig(k).gt.0)then
                if(indictronq.eq.0)then
                    res = res+res2(k)+sum1(k) &
                    -dlog(alpha)/(alpha)-gammlnN(1.d0/alpha) &
                !    -dlog(alpha)/(alpha)-dble(gammlnN(real(1./alpha))) &
                    +dlog(integrale1(k))

                endif
                if(indictronq.eq.1)then
                    if(AG.eq.1)then
!cccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                        res = res+res2(k)+sum1(k) &
                        -dlog(alpha)/(alpha)-gammlnN(1.d0/alpha) &
                    !    -dlog(alpha)/(alpha)-dble(gammlnN(real(1./alpha))) &
                        +dlog(integrale3(k))
                    else
! vraisemblance pr donnees censurees dte et tronquees a gauche
                        res = res+res2(k)+sum1(k)+dlog(integrale1(k)) &
                        -dlog(integrale2(k))
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpansplines=-1.d9
                    goto 123
                end if
            endif
        end do

    endif !fin boucle effet=2

!----------calcul de la penalisation -------------------

    pe1 = 0.d0
    pe2 = 0.d0
    do i=1,n-3
        pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
        *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
        the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
        m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
        the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
        m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
        *the1(i)*m1m(i))

        if (nst.eq.2) then
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
            *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
            the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
            m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
            the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
            m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
            *the2(i)*m1m(i))
        endif
    end do

    if (nst.eq.2) then
        pe = k0(1)*pe1 + k0(2)*pe2
    else
        pe = k0(1)*pe1
    endif

    resnonpen = res

    res = res - pe

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpansplines=-1.d9
        goto 123
    end if

    funcpansplines = res

123     continue

    return

    end function funcpansplines


!==========================  DISTANCE   =================================

    subroutine distancensplines(nz1,nz2,b,effet,mt,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)

    use tailles
    !use comon,only:I_hess,Hspl_hess,t0,t1,nsujet,nt0,nt1,date,c,ndate
    use comon,only:zi,nva,nst,H_hess

    implicit none

    integer::nz1,nz2,i,j,n,np,k,l,effet,mt
    double precision::x1,x2,h,su,bsup,binf,lam,lbinf,lbsup
    double precision,dimension(npmax,npmax)::hes1,hes2
    double precision,dimension(-2:npmax)::the1,the2
    double precision,dimension(npmax)::b
    double precision,dimension(mt),intent(out)::x1Out,x2Out
    double precision,dimension(mt,3),intent(out)::lamOut,suOut,lam2Out,su2Out

    n  = nz1+2

    if(nst.eq.2)then
        np  = nz1+2+nz2+2+effet+nva
    else
        np  = nz1+2+effet+nva
    endif

    do i=1,nz1+2
        do j=1,nz1+2
            hes1(i,j)=h_Hess(i,j)
        end do
    end do

    if(nst.eq.2)then
        k = 0
        do i=nz1+3,nz1+2+nz2+2
            k = k + 1 
            l = 0
            do j=nz1+3,nz1+2+nz2+2
                l = l + 1
                hes2(k,l)=h_Hess(i,j)
            end do
        end do
    endif

    do i=1,nz1+2
        the1(i-3)=(b(i))*(b(i))
    end do

    if(nst.eq.2)then
        do i=1,nz2+2
            j = nz1+2+i
            the2(i-3)=(b(j))*(b(j))
        end do
    endif

    h = (zi(n)-zi(1))*0.01d0
    x1 = zi(1)
    x2 = zi(1)

    do i=1,mt
        if(i .ne.1)then
            x1 = x1 + h
        end if
        call cospN(x1,the1,nz1+2,hes1,zi,bsup,su,binf,lbinf,lam,lbsup)
        if(bsup.lt.0.d0)then
            bsup = 0.d0
        endif
        if(binf.gt.1.d0)then
            binf = 1.d0
        endif
        if(lbinf.lt.0.d0)then
            lbinf = 0.d0
        endif

        x1Out(i)=x1
        lamOut(i,1)=lam
        lamOut(i,2)=lbinf
        lamOut(i,3)=lbsup 
        suOut(i,1)=su
        suOut(i,2)=binf! en fait tracé de la borne sup
        suOut(i,3)=bsup ! en fait tracé de la borne inf 

        if(nst.eq.2)then
            if(i.ne.1)then
                x2 = x2 + h
            endif 
            call cospN(x2,the2,nz2+2,hes2,zi,bsup,su,binf,lbinf,lam,lbsup)

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
        endif
    end do

    return

    end subroutine distancensplines
