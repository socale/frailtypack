
!========================          FUNCPAG_CPM         ====================
    double precision function funcpaG_tps(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:im,im2,im3,im1,mm,mm1,mm2,mm3,nva1,nva2,res4,stra,ve,vedc
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
    t0,t1,t0dc,t1dc,c,cdc,nsujet,nva, &
    nst,effet,ng,g,nig,AG,kkapa,indic_alpha,alpha,theta, &
    auxig,aux1,aux2,res1,res3,typeof,pe,resnonpen,nb_gl
    use residusM
    use comongroup
    use betatttps

    implicit none

    integer::np,id,jd,i,k,ig,choix,n,j,l,nb,cptg
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,sum,res,pe1,pe2,inv,dnb
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3,integrale3gap
    double precision,dimension(ngmax)::integrale4
    double precision::logGammaJ,int
    double precision,dimension(2)::k0
    double precision::result,result3,abserr,resabs,resasc,resultdc
    double precision,external::risqindivrec,risqindivdc,risqindiv

    kkapa=k0

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if(effet.eq.1) then
        theta = bh(np-(nva+npbetatps)-indic_alpha)**(2.d0)
        if (indic_alpha.eq.1) then ! new : joint more flexible alpha = 1
            alpha = bh(np-(nva+npbetatps))
        else
            alpha = 1.d0
        endif
    endif

!-------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------cc

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0

    res1dc = 0.d0
    res2dc = 0.d0

    cpt = 0
    integrale1 = 0.d0
    integrale2 = 0.d0
    integrale3 = 0.d0

    integrale4 = 0.d0
    integrale3gap = 0.d0
    aux1 = 0.d0
    aux2 = 0.d0

!*********************************************
!----- A SIMPLE SHARED FRAILTY  MODEL
!      write(*,*)'SIMPLE SHARED FRAILTY MODEL'
!*********************************************
    if(indic_joint.eq.0)then

        inv = 1.d0/theta
        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

!cccccccccccccccccccccc
! Fonction de risque de base recidive au temps T_ij
!cccccccccccccccccccccc
            if(c(i).eq.1)then
                res2(g(i)) = res2(g(i)) + dlog(risqindiv(t1(i),i,bh,np))
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaG_tps=-1.d9
                goto 123
            end if

!!cccccccccccccccccccccc
!! Fonction de risque cumulée de recidive au temps T_ij
!!cccccccccccccccccccccc
            call integration(risqindiv,0.d0,t1(i),result,abserr,resabs,resasc,i,bh,np)
            res1(g(i)) = res1(g(i)) + result

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpaG_tps=-1.d9
                goto 123
            end if

!!cccccccccccccccccccccc
!! Fonction de risque cumulée de recidive au temps T_i(j-1)
!!cccccccccccccccccccccc
            call integration(risqindiv,0.d0,t0(i),result3,abserr,resabs,resasc,i,bh,np)
            res3(g(i)) = res3(g(i)) + result3

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                funcpaG_tps=-1.d9
                goto 123
            end if

            if (res1(g(i)).lt.res3(g(i))) then
                funcpaG_tps=-1.d9
                goto 123
            end if

        end do

        res = 0.d0
        cptg = 0

!     gam2 = gamma(inv)
! k indice les groupes
!!cccccccccccccccccccccc
!! CALCUL DE LOG-VRAISEMBLANCE
!!cccccccccccccccccccccc
        do k=1,ng
            sum=0.d0
            if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))
                if (dnb.gt.1.d0) then
                    do l=1,nb
                        sum=sum+dlog(1.d0+theta*dble(nb-l))
                    end do
                endif
                if(theta.gt.(1.d-5)) then
!ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res= res-(inv+dnb)*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        + res2(k) + sum
!ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-(inv+dnb)*dlog(theta*res1(k)+1.d0)  &
                        +(inv)*dlog(theta*res3(k)+1.d0)+ res2(k) + sum
                    endif
                else
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
                        -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                        +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0)+res2(k)+sum
!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res-dnb*dlog(theta*res1(k)+1.d0)-res1(k)*(1.d0-theta*res1(k) &
                        /2.d0+theta*theta*res1(k)*res1(k)/3.d0) &
                        +res2(k)+sum &
                        +res3(k)*(1.d0-theta*res3(k)/2.d0 &
                        +theta*theta*res3(k)*res3(k)/3.d0)
                    endif
                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    !print*,"here4",res
                    funcpaG_tps=-1.d9
                    goto 123
                end if
            endif
        end do
    else !passage au modele conjoint

!*******************************************
!----- JOINT FRAILTY MODEL
!*******************************************

! pour les donnees recurrentes

        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

!cccccccccccccccccccccc
! Fonction de risque de base recidive au temps T_ij
!cccccccccccccccccccccc

            if (c(i).eq.1) then
                res2(g(i)) = res2(g(i)) + dlog(risqindivrec(t1(i),i,bh,np))
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaG_tps=-1.d9
                goto 123
            end if

!cccccccccccccccccccccc
! Fonction de risque cumulée de recidive au temps T_ij
!cccccccccccccccccccccc
            call integration(risqindivrec,0.d0,t1(i),result,abserr,resabs,resasc,i,bh,np)
            res1(g(i)) = res1(g(i)) + result

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                !print*,"here2"
                funcpaG_tps=-1.d9
                goto 123
            end if

!cccccccccccccccccc
! Fonction de risque de recidive au temps T_i(j-1)
!cccccccccccccccccc
            call integration(risqindivrec,0.d0,t0(i),result3,abserr,resabs,resasc,i,bh,np)
            res3(g(i)) = res3(g(i)) + result3

            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
                !print*,"here3"
                funcpaG_tps=-1.d9
                goto 123
            end if

        end do

! pour le deces

        do k=1,lignedc!ng
!cccccccccccccccccc
! Fonction de risque de deces au temps T_i*
!cccccccccccccccccc
            if (cdc(k).eq.1) then
                res2dc(gsuj(k)) = dlog(risqindivdc(t1dc(k),k,bh,np))
            endif

            if ((res2dc(gsuj(k)).ne.res2dc(gsuj(k))).or.(abs(res2dc(gsuj(k))).ge. 1.d30)) then
                !print*,"here4"
                funcpaG_tps=-1.d9
                goto 123
            end if

!cccccccccccccccccc
! Fonction de risque cumulée de dcd au temps T_i*
!cccccccccccccccccc
            call integration(risqindivdc,t0dc(k),t1dc(k),resultdc,abserr,resabs,resasc,k,bh,np)
            aux1(gsuj(k)) = aux1(gsuj(k)) + resultdc

            if ((aux1(gsuj(k)).ne.aux1(gsuj(k))).or.(abs(aux1(gsuj(k))).ge. 1.d30)) then
                !print*,"here5"
                funcpaG_tps=-1.d9
                goto 123
            end if
        end do

!***************INTEGRALES ****************************
        do ig=1,ng
            auxig = ig
            choix = 3
            call gaulagj(int,choix,nb_gl)
            integrale3(ig) = int !moins bon
            if ((integrale3(ig).eq.0.d0).and.(typeof.ne.2)) then
                integrale3(ig) = 1.d-300
            endif
        end do
!************** FIN INTEGRALES ************************

        res = 0.d0
        do k=1,ng
            sum = 0.d0
            if(cpt(k).gt.0)then
                if(theta.gt.(1.d-8)) then
!cccc ancienne vraisemblance : pour calendar sans vrai troncature cccccccc
                    res= res + res2(k) &
!--      pour le deces:
                    + res2dc(k)  &
                    - logGammaJ(1./theta)-dlog(theta)/theta  &
                    + dlog(integrale3(k))
                else
!*************************************************************************
!     developpement de taylor d ordre 3
!*************************************************************************
!                   write(*,*)'************** TAYLOR *************'
                    res= res + res2(k) &
                    + res2dc(k) &
                    - logGammaJ(1./theta)-dlog(theta)/theta  &
                    + dlog(integrale3(k))

                endif
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    !print*,"here6",k,res,res2(k),res2dc(k),integrale3(k)
                    funcpaG_tps=-1.d9
                    goto 123
                end if
            endif
        end do
    endif !fin boucle indic_joint=0 or 1

    if (typeof.eq.0) then ! penalisation pour les splines
        n = (np-(nva+npbetatps)-effet)/nst
        do k=1,n
            the1(k-3)=(bh(k))*(bh(k))
            j = n+k
            if (nst.eq.2) then
                the2(k-3)=(bh(j))*(bh(j))
            endif
        end do
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

        if (nst.eq.1) then
            pe2=0.d0
        end if

        pe = k0(1)*pe1 + k0(2)*pe2

        resnonpen = res

        res = res - pe
    endif

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        !print*,"here7"
        funcpaG_tps =-1.d9
        goto 123
    endif

    funcpaG_tps = res

    ! pour les martingales
    Rrec = res1
    Nrec = nig
    Rdc = aux1
    Ndc = cdc

123     continue

    return

    end function funcpaG_tps
