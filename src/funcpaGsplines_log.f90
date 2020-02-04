


!========================          FUNCPAG_SPLINES         ====================
    double precision function funcpaGsplines_log(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:kkapa,t0,t1,t0dc,t1dc
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m,mm3,mm2,mm1,mm,&
    im3,im2,im1,im,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc,date,datedc,zi,&
    c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet,nva,nva1,nva2,ndate,ndatedc,nst, &
    effet,stra,ve,vedc,pe,ng,g,nig,AG,indic_ALPHA,ALPHA,sig2, &
    auxig,aux1,aux2,res1,res3,res5,resnonpen,indictronq,nb_gh
    use residusM
    use comongroup

    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :
    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    integer::n,i,j,k,vj,ig,choix
    integer,dimension(ngmax)::cpt
    double precision::pe1,pe2,som1,som2,res,h1
    double precision,dimension(np)::bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc &
    ,res3dc,integrale1,integrale2,integrale3
!AD: for death,change dimension
    double precision,dimension(ndatemax)::dut1
    double precision,dimension(ndatemaxdc)::dut2
!AD:end
    double precision,dimension(0:ndatemax)::ut1
    double precision,dimension(0:ndatemaxdc)::ut2
    double precision::int
    double precision,parameter::pi=3.141592653589793d0

    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0
    ut1=0.d0
    ut2=0.d0
    dut2=0.d0
    dut1=0.d0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    n = (np-nva-effet-indic_ALPHA)/nst

    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i 
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do

    if(effet.eq.1) then
        sig2 = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
        if (indic_alpha.eq.1) then ! new : joint more flexible alpha = 1 
            alpha = bh(np-nva)
        else
            alpha = 1.d0
        endif
    endif

!----------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!AD:modify
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
    dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))

    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
    ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))

    ut1(0) = 0.d0
    ut2(0) = 0.d0

!//// NEW AMADOU vvv :
!--- strate1
    som1 = 0.d0
    vj = 0
    do i=2,ndate-1
        do k = 2,n-2
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
        do k = 2,n-2
            if (((datedc(i)).ge.(zi(k-1))).and.(datedc(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som2 = som2 + the2(j-4)
                vj  = j
                endif
            endif
        end do

        if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
            +(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))
            dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
            +(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))
        endif

    end do

!-------------fin strate2
    i = n-2
    h1 = (zi(i)-zi(i-1))

    ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)

    ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
    dut2(ndatedc) = (4.d0*the2(i-1)/h1)

!//// fin NEW AMADOU-vvv
!AD:end
!-------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    res1dc = 0.d0
    res2dc = 0.d0
    res3dc = 0.d0
    cpt = 0
    integrale1 = 0.d0
    integrale2 = 0.d0
    integrale3 = 0.d0
    aux1 = 0.d0
    aux2 = 0.d0
    res5 = 0.d0

!*******************************************
!----- A SIMPLE SHARED FRAILTY  MODEL
!      write(*,*)'SIMPLE SHARED FRAILTY MODEL'
!*******************************************

    if(indic_joint.eq.0)then

        do i=1,nsujet

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

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
                funcpaGsplines_log=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
                res5(i) = ut1(nt1(i))*vet
                res1(g(i)) = res1(g(i)) + res5(i) ! pour les residus
            endif

            if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ut2(nt0(i))*vet
                res5(i) = ut2(nt1(i))*vet
                res1(g(i)) = res1(g(i)) + res5(i) ! pour les residus
            endif

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
                funcpaGsplines_log=-1.d9
                goto 123
            end if

            if ((res5(i).ne.res5(i)).or.(abs(res5(i)).ge.1.d30)) then
                funcpaGsplines_log=-1.d9
                goto 123
            end if
        end do

!**************INTEGRALES ****************************
        do ig=1,ng
            auxig = ig
            choix = 1
            call gauherS(int,choix,nb_gh)
            integrale1(ig) = int
            if (AG.eq.1) then
                choix = 3
                call gauherS(int,choix,nb_gh)
                integrale3(ig) = int
            endif
            if (indictronq.eq.1) then
                choix = 2
                call gauherS(int,choix,nb_gh)
                integrale2(ig) = int
            endif
        end do
!************* FIN INTEGRALES ************************

        res = 0.d0
        do k=1,ng
            if(AG.EQ.1)then
                res = res+res2(k)-dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0+ &
                dlog(integrale3(k))
            else
                res = res+res2(k)+dlog(integrale1(k)) &
                -dlog(integrale2(k))
            endif
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                     funcpaGsplines_log=-1.d9
                     goto 123
            end if
        enddo
    else ! passage au modele conjoint

!********************************************
!*******************************************
!----- JOINT FRAILTY MODEL
!*********************************************

!     pour les donnees recurrentes

        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

            if(nva1.gt.0)then
                vet = 0.d0
                do j=1,nva1
                    vet = vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet = 1.d0
            endif

            if (c(i).eq.1) then
                res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
            endif
            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                !print*,"here7"
                funcpaGsplines_log=-1.d9
                goto 123
            end if

!     nouvelle version
            res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                !print*,"here6"
                funcpaGsplines_log=-1.d9
                goto 123
            end if

!     modification pour nouvelle vraisemblance / troncature:
            res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                !print*,"here5"
                funcpaGsplines_log=-1.d9
                goto 123
            end if
        end do

!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces
!ccccccccccccccccccccccccccccccccccccccccc

        do k=1,lignedc !ng
            if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
                    vet2 = vet2 + bh(np-nva2+j)*dble(vedc(k,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2 = 1.d0
            endif

            if(cdc(k).eq.1)then
                res2dc(gsuj(k)) = res2dc(gsuj(k)) + dlog(dut2(nt1dc(k))*vet2)
                if ((res2dc(gsuj(k)).ne.res2dc(gsuj(k))).or.(abs(res2dc(gsuj(k))).ge. 1.d30)) then
                    !print*,"here3",res2dc(gsuj(k))
                    funcpaGsplines_log=-1.d9
                    goto 123
                end if
            endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
            aux1(gsuj(k)) = aux1(gsuj(k)) + ut2(nt1dc(k))*vet2
            aux2(gsuj(k)) = aux2(gsuj(k)) + ut2(nt0dc(k))*vet2 !vraie troncature

            if ((aux1(gsuj(k)).ne.aux1(gsuj(k))).or.(abs(aux1(gsuj(k))).ge. 1.d30)) then
                funcpaGsplines_log=-1.d9
                goto 123
            end if
            if ((aux2(gsuj(k)).ne.aux2(gsuj(k))).or.(abs(aux2(gsuj(k))).ge. 1.d30)) then
                funcpaGsplines_log=-1.d9
                goto 123
            end if
        end do

!**************INTEGRALES ****************************
        do ig=1,ng
            auxig = ig
            choix = 3
            call gauherJ(int,choix,nb_gh)
            integrale3(ig) = int
        end do
!************* FIN INTEGRALES **************************

        res = 0.d0
        do k=1,ng
            if(cpt(k).gt.0)then
                !if(theta.gt.(1.d-8)) then
!cccc ancienne vraisemblance : pour calendar sans vrai troncature cccccccc
                    res= res + res2(k) &
!--      pour le deces:
                    + res2dc(k)  &
                    - dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0  &
                    + dlog(integrale3(k))
                !else
!*************************************************************************
!     developpement de taylor d ordre 3
!*************************************************************************
!                   write(*,*)'************** TAYLOR *************'
                !    res= res + res2(k) &
                !    + res2dc(k)  &
                !    - logGammaJ(1./theta)-dlog(theta)/theta  &
                !    + dlog(integrale3(k))
                !endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    !print*,"here"
                    funcpaGsplines_log=-1.d9
                    goto 123
                end if
            endif
        end do
    endif ! fin boucle indic_joint=0 ou 1

!---------- calcul de la penalisation -------------------

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
        if(nst.eq.1)then
            pe2=0.d0
        else
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
            *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
            the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
            m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
            the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
            m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
            *the2(i)*m1m(i))
        endif
    end do

    pe = k0(1)*pe1 + k0(2)*pe2

    resnonpen = res

    res = res - pe

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaGsplines_log=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpaGsplines_log = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
        end do
    end if

123     continue

    return

    end function funcpaGsplines_log

