

!========================          FUNCPA_WEIB_INTCENS          ====================
    double precision function funcpasweib_intcens(b,np,id,thi,jd,thj,k0)
    use tailles
    !use comon,only:AG,nig
    use comon,only:t0,t1,c,nsujet,nva, &
    nst,stra,ve,effet,ng,g,etaR,etaD,betaR,betaD,kkapa,theta,tU,d,dmax ! rajouts
        use residusM
    
    
    implicit none
    
! *** NOUVELLLE DECLARATION F90 :
    
    integer::np,id,jd,i,j,k,cptg,l,r
    integer,dimension(ngmax)::cpt,cptC
    double precision::thi,thj,inv,res,vet,num,den
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res1,res2,res3
    double precision,dimension(2)::k0
    double precision,dimension(1:2**dmax,1:ng)::p,p2 ! matrice
    integer,dimension(1:2**dmax,1:ng)::nbU,nbU2

    kkapa=k0
    j=0
    theta=0.d0

    bh = b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if (nst == 1) then
        betaR= bh(1)**2
        etaR= bh(2)**2
        etaD= 0.d0
        betaD= 0.d0
    else
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= bh(3)**2
        etaD= bh(4)**2
    end if

    if(effet.eq.1) then
        theta = bh(np-nva)*bh(np-nva)
    endif

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
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet
            endif
            if((c(i).eq.0).and.(stra(i).eq.2))then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge.1.d30)) then
                          funcpasweib_intcens=-1.d9
                          goto 123
            end if

            if(stra(i).eq.1)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
            endif
            if(stra(i).eq.2)then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaD)**betaD)*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge.1.d30)) then
                          funcpasweib_intcens=-1.d9
                          goto 123
            end if

            if (c(i).eq.1) then ! c = censure par intervalle
                r = 1
                do k=1,2**cptC(g(i))
                    if (stra(i).eq.1) then
                        p2(r,g(i)) = p(k,g(i))*dexp(((t1(i)/etaR)**betaR)*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(((tU(i)/etaR)**betaR)*vet)
                    endif
                    if (stra(i).eq.2) then
                        p2(r,g(i)) = p(k,g(i))*dexp(((t1(i)/etaD)**betaD)*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(((tU(i)/etaD)**betaD)*vet)
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
                          funcpasweib_intcens=-1.d9
                          goto 123
                endif
            endif
        end do

!*********************************************
!----avec un effet aleatoire dans le modele
!*********************************************

    else
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
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet
            endif
            if ((c(i).eq.0).and.(stra(i).eq.2)) then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet
            endif
            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge.1.d30)) then
                          funcpasweib_intcens=-1.d9
                          goto 123
                        end if

! produit de Kronecker
! création de la matrice p : 2**dmax lignes et ng colones
            if (c(i).eq.1) then ! c = censure par intervalle
                r = 1
                do k=1,2**cptC(g(i))
                    if (stra(i).eq.1) then
                        p2(r,g(i)) = p(k,g(i))*dexp(((t1(i)/etaR)**betaR)*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(((tU(i)/etaR)**betaR)*vet)
                    endif
                    if (stra(i).eq.2) then
                        p2(r,g(i)) = p(k,g(i))*dexp(((t1(i)/etaD)**betaD)*vet)
                        p2(r+1,g(i)) = p(k,g(i))*dexp(((tU(i)/etaD)**betaD)*vet)
                    endif
                    r = r + 2
                enddo
                cptC(g(i)) = cptC(g(i)) + 1
            endif
            p = p2

! modification pour nouvelle vraisemblance / troncature:
            if (stra(i).eq.1) then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaR)**betaR)*vet
            endif
            if (stra(i).eq.2) then
                res3(g(i)) = res3(g(i)) + ((t0(i)/etaD)**betaD)*vet
            endif
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
              funcpasweib_intcens=-1.d9
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
                do r=1,2**d(k)
                    num = ((-1)**nbU(r,k))*((1.d0+theta*res3(k))**inv)
                    den = (1.d0+theta*res1(k)+theta*dlog(p(r,k)))**inv
                    res2(k) = res2(k) + num/den
                enddo
                if (res2(k).eq.0.d0) then
                    res2(k) = 1.d-300
                endif

                if ((res2(k).ne.res2(k)).or.(abs(res2(k)).ge.1.d30)) then
                  funcpasweib_intcens=-1.d9
                  goto 123
                end if

                res = res + dlog(res2(k))
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    !print*,"here",res2(k),dlog(res2(k)),inv*dlog(theta)
                  funcpasweib_intcens=-1.d9
                  goto 123
                end if
            endif
        end do

    endif !fin boucle effet=0

!--------------------------------------------------

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpasweib_intcens=-1.d9
        goto 123
    end if

    funcpasweib_intcens = res

    do k=1,ng
        cumulhaz(k)=res1(k)
    end do

123     continue

    return

    end function funcpasweib_intcens
