

!========================          FUNCPA_SPLINES_INTCENS          ====================
    double precision function funcpassplines_intcens(b,np,id,thi,jd,thj,k0)
    use tailles
    !use comon,only:AG,nig,nz1,nz2,t0,t1,tU
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
    mm3,mm2,mm1,mm,im3,im2,im1,im,date,zi,c,nt0,nt1,nsujet,nva,ndate, &
    nst,stra,ve,pe,effet,ng,g,resnonpen,theta,ntU,d,dmax
        use residusM
    
        implicit none

! *** NOUVELLLE DECLARATION F90 :
    integer::n,np,id,jd,i,j,k,vj,cptg,l,r
    integer,dimension(ngmax)::cpt,cptC
    double precision::thi,thj,pe1,pe2,inv,som1,som2,res,vet,h1,num,den
    double precision,dimension(-2:npmax)::the1,the2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res1,res2,res3
    double precision,dimension(2)::k0
    double precision,dimension(ndatemax)::dut1,dut2
    double precision,dimension(0:ndatemax)::ut1,ut2
    double precision,dimension(1:2**dmax,1:ng)::p,p2 ! matrice
    integer,dimension(1:2**dmax,1:ng)::nbU,nbU2

    j=0  
    
    theta=0.d0
    do i=1,np
        bh(i)=b(i)
    end do 
    dut1=0.d0
    dut2=0.d0
    ut1=0.d0
    ut2=0.d0
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    n = (np-nva-effet)/nst ! nombre de eta (coef des splines)

    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do


    if(effet.eq.1) then
        theta = bh(np-nva)*bh(np-nva)
    endif

!---------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

    vj = 0
    som1 = 0.d0
    som2 = 0.d0
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))

    dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
    ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
    ut1(0) = 0.d0
    ut2(0) = 0.d0
    do i=2,ndate-1 ! pour toutes les dates du jeu de données
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
    ut2(ndate)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)
    dut2(ndate) = (4.d0*the2(i-1)/h1)

!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!--------------------------------------------------------

!--- avec ou sans variable explicative  ------cc
    res1= 0.d0
    res2= 0.d0
    res3= 0.d0
    cpt = 0
    cptC = 0
    p = 1.d0
    p2 = p
    nbU = 0
    nbU2 = 0

!*******************************************
!---- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then
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

            if((c(i).eq.0).and.(stra(i).eq.1))then
                res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
            endif
            if((c(i).eq.0).and.(stra(i).eq.2))then
                res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge.1.d30)) then
                funcpassplines_intcens=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
                RisqCumul(i) = ut1(nt1(i))*vet
            endif
            if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ut2(nt0(i))*vet 
                RisqCumul(i) = ut2(nt1(i))*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge.1.d30)) then
                funcpassplines_intcens=-1.d9
                goto 123
            end if

            if (c(i).eq.1) then ! c = censure par intervalle
                r = 1
                do k=1,2**cptC(g(i))
                    if (stra(i).eq.1) then
                        p2(r,g(i)) = p(k,g(i))*dexp(ut1(nt1(i))*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(ut1(ntU(i))*vet)
                    endif
                    if (stra(i).eq.2) then
                        p2(r,g(i)) = p(k,g(i))*dexp(ut2(nt1(i))*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(ut2(ntU(i))*vet)
                    endif
                    r = r + 2
                enddo
                cptC(g(i)) = cptC(g(i)) + 1
            endif
            p = p2

        end do

        res = 0.d0
        cptg = 0

! k indice les groupes
        do k=1,ng
            if(cpt(k).gt.0)then
                ! matrice du nombre de U
                if (d(k).gt.0) then
                    r = 1
                    do l=1,2**(d(k)-1)
                        nbU2(r,k) = nbU(l,k)
                        nbU2(r+1,k) = nbU(l,k) + 1
                        r = r + 2
                        nbU = nbU2
                    enddo
                endif
                ! produit de Kronecker
                do r=1,2**d(k)
                    res2(k) = res2(k) + ((-1)**nbU(r,k))*dexp(-res1(k)-dlog(p(r,k)))
                enddo
                res = res + dlog(res2(k)) + res3(k)
                cptg = cptg + 1
                if ((res.ne.res).or.(abs(res).ge.1.d30)) then
                    funcpassplines_intcens=-1.d9
                    goto 123
                endif
            endif
        end do

!*********************************************
!----avec un effet aleatoire dans le modele
!*********************************************

    else
! write(*,*)'debut funcpa interval censoring'
!      write(*,*)'AVEC EFFET ALEATOIRE'
        inv = 1.d0/theta

        do i=1,nsujet ! i indice les sujets

            cpt(g(i))=cpt(g(i))+1

! creation de l'exp(bX) de l'individu
            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet = vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet = 1.d0
            endif

! t0 = troncature, t1 = L, tU = U
! pour les i qui sont censurés à droite : tU le meme que t1
            if ((c(i).eq.0).and.(stra(i).eq.1)) then ! censure a droite
                res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
            endif
            if ((c(i).eq.0).and.(stra(i).eq.2)) then
                res1(g(i)) = res1(g(i)) + ut2(nt1(i))*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge.1.d30)) then
                funcpassplines_intcens=-1.d9
                goto 123
            end if

! produit de Kronecker
! création de la matrice p : 2**dmax lignes et ng colones
            if (c(i).eq.1) then ! c = censure par intervalle
                r = 1
                do k=1,2**cptC(g(i))
                    if (stra(i).eq.1) then
                        p2(r,g(i)) = p(k,g(i))*dexp(ut1(nt1(i))*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(ut1(ntU(i))*vet)
                    endif
                    if (stra(i).eq.2) then
                        p2(r,g(i)) = p(k,g(i))*dexp(ut2(nt1(i))*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(ut2(ntU(i))*vet)
                    endif
                    r = r + 2
                enddo
                cptC(g(i)) = cptC(g(i)) + 1
            endif
            p = p2

! modification pour nouvelle vraisemblance / troncature:
            if (stra(i).eq.1) then
                res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet
            endif
            if (stra(i).eq.2) then
                res3(g(i)) = res3(g(i)) + ut2(nt0(i))*vet
            endif
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
              funcpassplines_intcens=-1.d9
              goto 123
            end if

        enddo

        res = 0.d0
        cptg = 0

! k indice les groupes
        do k=1,ng
            if(cpt(k).gt.0)then
! création de nbU : nombre de U dans chaque élément de la matrice p
                if (d(k).gt.0) then
                    r = 1
                    do l=1,2**(d(k)-1)
                        nbU2(r,k) = nbU(l,k)
                        nbU2(r+1,k) = nbU(l,k) + 1
                        r = r + 2
                        nbU = nbU2
                    enddo
                endif
! r indice les éléments dans le produit de Kronecker
! res3 null if no left truncation
                do r=1,2**d(k)
                    num = ((-1)**nbU(r,k))*((1.d0+theta*res3(k))**inv)
                    den = (1.d0+theta*res1(k)+theta*dlog(p(r,k)))**inv
                    res2(k) = res2(k) + num/den
                enddo
                if ((res2(k).ne.res2(k)).or.(abs(res2(k)).ge.1.d30)) then
                  funcpassplines_intcens=-1.d9
                  goto 123
                end if

                res = res + dlog(res2(k))
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                  funcpassplines_intcens=-1.d9
                  goto 123
                end if
            endif
        end do

    endif !fin boucle effet=0

!--------- calcul de la penalisation -------------------

    pe1 = 0.d0
    pe2 = 0.d0
    pe=0.d0
    do i=1,n-3

        pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
        *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
        the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
        m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
        the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
        m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
        *the1(i)*m1m(i))
        pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
        *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
        the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
        m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
        the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
        m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
        *the2(i)*m1m(i))

    end do


!    Changed JRG 25 May 05
    if (nst.eq.1) then
        pe2=0.d0
    end if

    pe = k0(1)*pe1 + k0(2)*pe2 

    resnonpen = res

    res = res - pe

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpassplines_intcens=-1.d9
        goto 123
    end if

    funcpassplines_intcens = res

    do k=1,ng
        cumulhaz(k)=res1(k)
    end do

123     continue

    return

    end function funcpassplines_intcens


!==========================  DISTANCE   =================================
! fonction supprimée car déjà définie dans funcpassplines.f90


