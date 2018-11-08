
!========================          funcpaMultivCpm         ====================
    double precision function funcpaMultivCpm(b,np,id,thi,jd,thj,k0)

    use taillesmultiv
    ! use comonmultiv,only:AG,cens,nst,stra,t0dc
    use comonmultiv,only:nbintervR,nbintervDC,t0,t1,t1dc,c,cdc,nsujet,nva,nva1,nva2, &
    effet,ve,vedc,ng,g,nig,indic_ALPHA,ALPHA,theta, &
    auxig,aux1,aux2,res1,res3,res4,ttt,tttdc,betacoef,kkapa, &
!add meta
    nbintervM,t0meta,t1meta,cmeta,nva3,vemeta,tttmeta,eta,alpha1,alpha2,nsujetmeta,gmeta,&
    indic_a1,indic_a2,indic_eta,res1meta,res3meta,nigmeta,indic_rho
    use residusMmultiv

    implicit none

    integer::np,id,jd,i,j,k,ig
    integer,dimension(ngmax)::cpt,cptmeta
    double precision::som11,som21,somf
    integer::jj,gg,gg2
    double precision::thi,thj,res,vet,vet2,vet3
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc,res2meta
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3,integrale3gap
    double precision,dimension(ngmax)::integrale4
    double precision::int
    double precision,dimension(3)::k0
    double precision,parameter::pi=3.141592653589793d0

    kkapa=k0
    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

! Allocation des parametre des fcts de risque cte par morceaux
!cpm
    betacoef = 0.d0
    do i = 1,(nbintervR+nbintervDC+nbintervM)
        betacoef(i)=bh(i)**2
    end do
!cpm

    if(effet.eq.1) then
        theta = bh(np-nva-indic_ALPHA)*bh(np-nva-indic_ALPHA)
!        alpha = bh(np-nva)
!add
        alpha =bh(np-nva-indic_rho) !rho
        alpha1=bh(np-nva-indic_a1)
        alpha2=bh(np-nva-indic_a2)
        eta = bh(np-nva-indic_eta)*bh(np-nva-indic_eta)
    endif


!-------------------------------------------------------
!---------- calcul de la vraisemblance ------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------cc

    do k=1,ng
        res1(k) = 0.d0
        res2(k) = 0.d0
        res3(k) = 0.d0
        res4(k) = 0.d0
        res1dc(k) = 0.d0
        res2dc(k) = 0.d0
        cpt(k) = 0
        integrale1(k) = 0.d0
        integrale2(k) = 0.d0
        integrale3(k) = 0.d0
        integrale4(k) = 0.d0
        integrale3gap(k) = 0.d0
        aux1(k)=0.d0
        aux2(k)=0.d0
        res1meta(k) = 0.d0
        res2meta(k) = 0.d0
        res3meta(k) = 0.d0
        cptmeta(k) = 0
    end do

!**********************************************
!-----avec un effet aleatoire dans le modele
!**********************************************


!==========================================================================
!     pour les donnees loco
!==========================================================================
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
            do gg=1,nbintervR
                if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                    res2(g(i)) =  res2(g(i))+ dlog(betacoef(gg)*vet)
                endif
            end do
        endif  
        if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
            funcpaMultivCpm=-1.d9
            goto 123
        end if

!cccccccccccccccccccccc
! Fonction de risque cumulée de recidive au tepms T_ij
!cccccccccccccccccccccc

        som11=0.d0
        som21=0.d0
        gg2=0
        do gg=1,nbintervR
            if((t1(i).gt.(ttt(gg-1))).and.(t1(i).le.(ttt(gg))))then
                som11=betacoef(gg)*(t1(i)-ttt(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                    end do
                endif
                res1(g(i)) = res1(g(i))+vet*(som11+som21)
            endif
        end do
        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            funcpaMultivCpm=-1.d9
            goto 123
        end if

!cccccccccccccccccc
! Fonction de risque de recidive au tepms T_i(j-1)
!cccccccccccccccccc

        som11=0.d0
        som21=0.d0
        gg2=0
        do gg=1,nbintervR
            if ((t0(i).ge.(ttt(gg-1))).and.(t0(i).lt.(ttt(gg)))) then
                som11=betacoef(gg)*(t0(i)-ttt(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                    end do
                endif
                res3(g(i)) = res3(g(i))+vet*(som11+som21)
            endif
        end do
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            funcpaMultivCpm=-1.d9
            goto 123
        end if
    end do

!==========================================================================
!     pour les donnees dc
!==========================================================================

    do k=1,ng
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2-nva3+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2=1.d0
        endif

!cccccccccccccccccc
! Fonction de risque de deces au tepms T_i*
!cccccccccccccccccc

        if(cdc(k).eq.1)then

            do gg=1,nbintervDC
                if ((t1dc(k).gt.(tttdc(gg-1))).and.(t1dc(k).le.(tttdc(gg)))) then
                    res2dc(k) = dlog(betacoef(nbintervR+gg)*vet2)
                endif
            end do
            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                funcpaMultivCpm=-1.d9
                goto 123
            end if
        endif

!cccccccccccccccccc
! Fonction de risque cumulée de dcd au tepms T_i*
!cccccccccccccccccc

        som11=0.d0
        som21=0.d0
        somf=0.d0
        gg2=0
        do gg=1,nbintervDC
            if((t1dc(k).gt.(tttdc(gg-1))).and.(t1dc(k).le.(tttdc(gg))))then
                som11=betacoef(gg+nbintervR)*(t1dc(k)-tttdc(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(jj+nbintervR)*(tttdc(jj)-tttdc(jj-1))
                    end do
                endif
            endif
            aux1(k)=(som11+som21)*vet2
        end do

        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            funcpaMultivCpm=-1.d9
            goto 123
        end if
    end do

!==========================================================================
!     pour les donnees Meta
!==========================================================================

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
            do gg=1,nbintervM
                if((t1meta(i).ge.(tttmeta(gg-1))).and.(t1meta(i).lt.(tttmeta(gg))))then
                    res2meta(gmeta(i)) =  res2meta(gmeta(i))+ dlog(betacoef(gg+nbintervR+nbintervDC)*vet3)
                endif
            end do
        endif
        if ((res2meta(gmeta(i)).ne.res2meta(gmeta(i))).or.(abs(res2meta(gmeta(i))).ge. 1.d30)) then
            funcpaMultivCpm=-1.d9
            goto 123
        end if

!cccccccccccccccccccccc
! Fonction de risque cumulée de recidive au tepms T_ij
!cccccccccccccccccccccc

        som11=0.d0
        som21=0.d0
        gg2=0
        do gg=1,nbintervM
            if((t1meta(i).gt.(tttmeta(gg-1))).and.(t1meta(i).le.(tttmeta(gg))))then
                som11=betacoef(gg+nbintervR+nbintervDC)*(t1meta(i)-tttmeta(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(jj+nbintervR+nbintervDC)*(tttmeta(jj)-tttmeta(jj-1))
                    end do
                endif
                res1meta(gmeta(i)) = res1meta(gmeta(i))+vet3*(som11+som21)
            endif
        end do
        if ((res1meta(gmeta(i)).ne.res1meta(gmeta(i))).or.(abs(res1meta(gmeta(i))).ge. 1.d30)) then
            funcpaMultivCpm=-1.d9
            goto 123
        end if

!cccccccccccccccccc
! Fonction de risque de recidive au tepms T_i(j-1)
!cccccccccccccccccc

        som11=0.d0
        som21=0.d0
        gg2=0
        do gg=1,nbintervM
            if ((t0meta(i).ge.(tttmeta(gg-1))).and.(t0meta(i).lt.(tttmeta(gg)))) then
                som11=betacoef(gg+nbintervR+nbintervDC)*(t0meta(i)-tttmeta(gg-1))
                gg2=gg
                if (gg2.ge.2)then
                    do jj=1,gg2-1
                        som21=som21+betacoef(jj+nbintervR+nbintervDC)*(tttmeta(jj)-tttmeta(jj-1))
                    end do
                endif
                res3meta(gmeta(i)) = res3meta(gmeta(i))+vet3*(som11+som21)
            endif
        end do
        if ((res3meta(gmeta(i)).ne.res3meta(gmeta(i))).or.(abs(res3meta(gmeta(i))).ge. 1.d30)) then
            funcpaMultivCpm=-1.d9
            goto 123
        end if
    end do

!==========================================================================
!     fin pour les donnees Meta
!==========================================================================

!***************INTEGRALES ****************************
    do ig=1,ng
        auxig=ig
!        choix = 3
!        call gaulagj(int,choix)
        call gausshermiteBIS2011(int,30)!!!!! garder 30!!!!
        integrale3(ig) = int !moins bon
    end do

!************** FIN INTEGRALES **************************

    res = 0.d0
    do k=1,ng
        res= res + res2(k)+res2meta(k)+res2dc(k)+dlog(integrale3(k))-dlog((2.d0)*pi*sqrt(theta*eta* &
             (1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2)))
    end do

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaMultivCpm =-1.d9
        Rrec=0.d0
        Nrec=0
        Rdc=0.d0
        Ndc=0
        Rrec2=0.d0
        Nrec2=0
        goto 123
    else
        funcpaMultivCpm = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
            Rrec2(k)=res1meta(k)
            Nrec2(k)=nigmeta(k)
        end do
    end if

123     continue

    return

    end function funcpaMultivCpm


