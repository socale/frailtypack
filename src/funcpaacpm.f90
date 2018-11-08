!===================================================================
!========================          FUNCPA       ====================
    double precision function funcpaacpm(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:alpha,im,
    use comon,only:effet,stra,t0,t1,c,nsujet,nva,nst,auxig,ng,ve,g,nig,indictronq,&
    nbintervR,ttt,betacoef,kkapa,indic_tronc
    !use additiv,only:invD,rho
    use additiv,only:correl,ngexact,ve2,sigma2,tau2,cov,mid,betaaux
    use residusM,only:som_Xbeta,indic_cumul,cumulhaz

    implicit none

    integer::np,id,jd,i,j,k,cptg,ig,ip
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,som1,som2,res,vet,somm1,somm2
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::integrale1,funcaux
    double precision,dimension(2)::k0

    real,parameter::pi = 3.1415926535
    integer::restar,nf
!****** pour la maximisation avec marq98aux
    integer :: npaux,niaux,ieraux,istopaux
    DOUBLE PRECISION::resaux
    double precision,dimension((npmax*(npmax+3)/2))::vaux
    double precision,dimension(npmax)::baux
!****** derivanal
    double precision , dimension(ngmax)::res1,res2,res3,res4
    double precision , dimension(ngmax)::res5,res6
!****** u_tilde
    double precision  :: u_tilde,v_tilde
    integer::gg,jj
    double precision,external::funcpaocpm
!******************

    kkapa=k0
    somm1=0.d0
    somm2=0.d0
    vet=0.d0
    j=0
    res=0.d0
    som2=0.d0
    restar = 0
    nf = 1

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    betacoef = 0.d0

    do i=1,nst*nbintervR
        betacoef(i)=bh(i)**2
    end do

    if(effet.eq.1) then
!terme de correlation entre 2 frailties/avec contraintes
        sigma2 = (bh(np-nva-1)*bh(np-nva-1)) ! variance intercept
        tau2 = (bh(np-nva)*bh(np-nva))  ! variance traitement * groupe
        cov=dsqrt(sigma2*tau2)*(2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)
    endif

!---------------------------------------------------------
!----------calcul de la vraisemblance --------------------
!---------------------------------------------------------

!---- avec ou sans variable explicative  ------cc

    do ig=1,ngexact!ng!
        res1(ig) = 0.d0
        res2(ig) = 0.d0
    end do
    cpt=0

!*******************************************
!----- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then

        do i=1,nsujet

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
                do gg=1,nbintervR !!
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                         res2(g(i)) = res2(g(i))+dlog(betacoef(gg)*vet)
                    end if!!
                end do !!
            endif

            if((c(i).eq.1).and.(stra(i).eq.2))then
                do gg=1,nbintervR !!
                    if((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg))))then
                          res2(g(i)) = res2(g(i))+dlog(betacoef(nbintervR+gg)*vet)
                    end if!!
                end do !!
            endif

           if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
               funcpaacpm=-1.d9
               goto 123
           end if

            if(stra(i).eq.1)then
                som1=0.d0
                som2=0.d0
                somm1=0.d0
                somm2=0.d0
                do gg=1,nbintervR
                    if ((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg)))) then
                        som1=betacoef(gg)*(t1(i)-ttt(gg-1))
                        if (gg.ge.2)then
                            do jj=1,gg-1
                                som2=som2+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                            end do!!
                        endif!!

                        res1(g(i)) = res1(g(i)) + (som1+som2)*vet

                        if (indic_cumul==1) then
                            cumulhaz(g(i)) = (som1+som2)
                        end if
                    end if!!

                    if ((t0(i).ge.(ttt(gg-1))).and.(t0(i).lt.(ttt(gg)))) then
                        somm1=betacoef(gg)*(t0(i)-ttt(gg-1))
                        if (gg.ge.2)then
                            do jj=1,gg-1
                                somm2=somm2+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                            end do!!
                        endif!!
                        res1(g(i)) = res1(g(i)) - (somm1+somm2)*vet
                    end if
                end do!!
                if (indic_cumul==1) then
                    cumulhaz(g(i)) = (som1+som2)
                end if
            endif

            if(stra(i).eq.2)then
                som1=0.d0
                som2=0.d0
                somm1=0.d0
                somm2=0.d0

                do gg=1,nbintervR
                    if ((t1(i).ge.(ttt(gg-1))).and.(t1(i).lt.(ttt(gg)))) then
                        som1=betacoef(nbintervR+gg)*(t1(i)-ttt(gg-1))
                        if (gg.ge.2)then
                            do jj=1,gg-1
                                som2=som2+betacoef(nbintervR+jj)*(ttt(jj)-ttt(jj-1))
                            end do!!
                        endif!!
                        res1(g(i)) = res1(g(i)) + (som1+som2)*vet 
                    end if!!

                    if ((t0(i).ge.(ttt(gg-1))).and.(t0(i).lt.(ttt(gg)))) then
                        somm1=betacoef(nbintervR+gg)*(t0(i)-ttt(gg-1))
                        if (gg.ge.2)then
                            do jj=1,gg-1
                                somm2=somm2+betacoef(nbintervR+jj)*(ttt(jj)-ttt(jj-1))
                            end do!!
                        endif!!
                        res1(g(i)) = res1(g(i)) - (somm1+somm2)*vet
                    end if
                end do
                if (indic_cumul==1) then
                    cumulhaz(g(i)) = (som1+som2)
                end if
            endif

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpaacpm=-1.d9
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
                    funcpaacpm=-1.d9
                    goto 123
                end if
            endif
        end do

    else
!*************************************************************
!-----avec deux effets aleatoires dans le modele et correl=0
!************************************************************
        if(effet.eq.1.and.correl.eq.0)then

            mid=0
            integrale1=0.d0

            do k=1,nsujet
                if(c(k).eq.1)then
                    mid(g(k))=mid(g(k))+1
                endif
            end do 

            do k=1,nsujet

                if(nva.gt.0)then
                    vet = 0.d0
                    do ip=1,nva
                        vet =vet + bh(np-nva+ip)*ve(k,ip)
                    end do
                    vet = dexp(vet)
                else
                    vet=1.d0
                endif

                if((c(k).eq.1).and.(stra(k).eq.1))then
                    do gg=1,nbintervR
                        if((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg))))then
                            res2(g(k)) = res2(g(k))+dlog(betacoef(gg)*vet)
                        end if
                    end do
                endif

                if((c(k).eq.1).and.(stra(k).eq.2))then
                    do gg=1,nbintervR
                        if((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg))))then
                            res2(g(k)) = res2(g(k))+dlog(betacoef(nbintervR+gg)*vet)
                        end if
                    end do
                endif

                if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                    funcpaacpm=-1.d9
                    goto 123
                end if

            end do

!=========================================================================
!==== calcul des integrales par transformation de LAPLACE pour chq gpe ===
!=========================================================================

            do ig=1,ng
                res3(ig)=0.d0
                res4(ig)=0.d0
                res5(ig)=0.d0
                res6(ig)=0.d0
                auxig=ig
                baux(1)=0.05d0      !initialisation de u_tilde
                baux(2)=0.05d0      !initialisation de v_tilde
!======================================================================
! maximisation pour deux parametres, les 2 effets aleatoires : u et v
!======================================================================
                npaux=2
                niaux=0
                ieraux=0
                istopaux=0
                vaux=0.d0!vecteur derivees 2nd et 1ERES
                funcaux=0.d0!vraisemblance
                resaux=0.d0!vraisemblance
                do ip=1,nva
                    betaaux(ip)= bh(np-nva+ip)
                end do

    call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux,funcpaocpm)
       if (ieraux .eq.-1) then
          funcpaacpm=-1.d9
          goto 123
       end if
                u_tilde = baux(1)!u_tilde
                v_tilde = baux(2)!v_tilde

                do k=1,nsujet

                    if(nva.gt.0.and.g(k).eq.ig)then
                        vet = 0.d0 
                        do ip=1,nva
                            vet =vet + bh(np-nva+ip)*ve(k,ip)
                        end do
                        vet = dexp(vet)
                    else
                        vet=1.d0
                    endif

                    if(g(k).eq.ig)then

                        if(c(k).eq.1)then
                            res3(ig) = res3(ig)+u_tilde+v_tilde*ve2(k,1)
                        endif

                        if(stra(k).eq.1)then
                            somm1=0.d0
                            som1=0.d0
                            som2=0.d0
                            do gg=1,nbintervR
                                if ((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg)))) then
                                    som1=betacoef(gg)*(t1(k)-ttt(gg-1))
                                    if (gg.ge.2)then
                                        do jj=1,gg-1
                                            som2=som2+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                                        end do!!
                                    endif!!
                                    res4(ig) = res4(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1)) 
                                    res5(ig) = res5(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                                    res6(ig) = res6(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                                end if!!
                            end do

                        endif

                        if(stra(k).eq.2)then
                            som1=0.d0
                            som2=0.d0
                            somm2=0.d0
                            do gg=1,nbintervR
                                if ((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg)))) then
                                    som1=betacoef(nbintervR+gg)*(t1(k)-ttt(gg-1))
                                    if (gg.ge.2)then
                                        do jj=1,gg-1
                                            som2=som2+betacoef(nbintervR+jj)*(ttt(jj)-ttt(jj-1))
                                        end do!!
                                    endif!!
                                    res4(ig) = res4(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1)) 
                                    res5(ig) = res5(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                                    res6(ig) = res6(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                                end if!!
                            end do

                        endif

                    endif
                end do

                som_Xbeta(ig) = vet

!=====fin maximisation aux

                integrale1(ig)= -0.5d0*dlog(sigma2*tau2) &
                -0.5d0* &
                dlog(res4(ig)*res6(ig)+res4(ig)/tau2+res6(ig)/sigma2 &
                +1.d0/(sigma2*tau2)-res5(ig)*res5(ig)) & !det(K"(b))
                +res2(ig)+res3(ig)-res4(ig) &
                -0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2)!-k(b)
            end do

!======= fin calcul de log integrale
!======================================================================

            res = 0.d0
            do k=1,ng!ngexact
                if(nig(k).gt.0)then

                    if(indictronq.eq.0)then
                        res = res+integrale1(k)!integrale1 donne le log de I directement
                        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                            funcpaacpm=-1.d9
                            goto 123
                        end if
                    else !troncature
                        indic_tronc = 1
                        funcpaacpm=-1.d9
                        goto 123
                    endif
                endif
            end do
        endif                !fin boucle effet= 1 and correl = 0

!=======================================================================

!************************************************************
!-----avec deux effets aleatoires dans le modele et correl=1
!************************************************************

        if(effet.eq.1.and.correl.eq.1)then

            mid=0
            integrale1=0.d0

            do k=1,nsujet
                if(c(k).eq.1)then
                    mid(g(k))=mid(g(k))+1
                endif
            end do 

            do k=1,nsujet
                if(nva.gt.0)then
                    vet = 0.d0 
                    do ip=1,nva
                        vet =vet + bh(np-nva+ip)*ve(k,ip)
                    end do
                    vet = dexp(vet)
                else
                    vet=1.d0
                endif

                if((c(k).eq.1).and.(stra(k).eq.1))then
                    do gg=1,nbintervR !!
                        if((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg))))then
                            res2(g(k)) = res2(g(k))+dlog(betacoef(gg)*vet)
                        end if!!
                    end do !!
                endif
                if((c(k).eq.1).and.(stra(k).eq.2))then
                    do gg=1,nbintervR !!
                        if((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg))))then
                            res2(g(k)) = res2(g(k))+dlog(betacoef(nbintervR+gg)*vet)
                        end if!!
                    end do !!
                endif 
                if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                    funcpaacpm=-1.d9
                    goto 123
                end if
            end do

!=========================================================================
!==== calcul des integrales par transformation de LAPLACE pour chq gpe ===
!=========================================================================

            do ig=1,ngexact!ng!
                res3(ig)=0.d0
                res4(ig)=0.d0
                res5(ig)=0.d0
                res6(ig)=0.d0
                auxig=ig
                baux(1)=0.05d0   !initialisation de u_tilde
                baux(2)=0.05d0   !initialisation de v_tilde

!======================================================================
! maximisation pour deux parametres, les 2 effets aleatoires : u et v
!======================================================================
                npaux=2
                niaux=0
                vaux=0.d0!vecteur derivees 2nd et 1ERES
                funcaux=0.d0!vraisemblance
                resaux=0.d0!vraisemblance
                do ip=1,nva
                    betaaux(ip)= bh(np-nva+ip)
                end do
    call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux,funcpaocpm)
       if (ieraux .eq.-1) then
          funcpaacpm=-1.d9
          goto 123
       end if
                u_tilde = baux(1)!u_tilde(ig)
                v_tilde = baux(2)!v_tilde

                if(effet.eq.1) then
                end if
                do k=1,nsujet
                    if(nva.gt.0.and.g(k).eq.ig)then
                        vet = 0.d0 
                        do ip=1,nva
                            vet =vet + bh(np-nva+ip)*ve(k,ip)
                        end do
                        vet = dexp(vet)
                    else
                        vet=1.d0
                    endif

                    if(g(k).eq.ig)then
                        if(c(k).eq.1)then
                            res3(ig) = res3(ig)+u_tilde+v_tilde*ve2(k,1)
                        endif

                        if(stra(k).eq.1)then
                            som1=0.d0
                            som2=0.d0

                            do gg=1,nbintervR
                                if ((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg)))) then
                                    som1=betacoef(gg)*(t1(k)-ttt(gg-1))
                                    if (gg.ge.2)then
                                        do jj=1,gg-1
                                            som2=som2+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                                        end do!!
                                    endif!!
                                    res4(ig) = res4(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1)) 
                                    res5(ig) = res5(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                                    res6(ig) = res6(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                                end if!!
                            end do

                        endif
                        if(stra(k).eq.2)then
                            som1=0.d0
                            som2=0.d0
                            do gg=1,nbintervR
                                if ((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg)))) then
                                    som1=betacoef(nbintervR+gg)*(t1(k)-ttt(gg-1))
                                    if (gg.ge.2)then
                                        do jj=1,gg-1
                                            som2=som2+betacoef(nbintervR+jj)*(ttt(jj)-ttt(jj-1))
                                        end do!!
                                    endif!!
                                    res4(ig) = res4(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1)) 
                                    res5(ig) = res5(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                                    res6(ig) = res6(ig) + (som1+som2)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
                                end if!!
                            end do
                        endif
                    endif
                end do

!=====fin maximisation aux
                som_Xbeta(ig) = vet

                cov=dsqrt(sigma2*tau2)* & !avec contrainte
                (2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)

                integrale1(ig)= -0.5d0*dlog(sigma2*tau2-cov*cov) &!-0.5logdet(D)
                -0.5d0* &
                dlog((res4(ig)+tau2/(tau2*sigma2-cov**2)) &
                *(res6(ig)+sigma2/(tau2*sigma2-cov**2)) &
                -(res5(ig)-cov/(tau2*sigma2-cov**2))**2) &!-0.5*log(det(ka2)
                +res2(ig)+res3(ig)-res4(ig) &
                -0.5d0*((u_tilde**2)/sigma2+(v_tilde**2)/tau2 &
                -2.d0*u_tilde*v_tilde*cov/(sigma2*tau2)) &
                /(1.d0-(cov**2)/(sigma2*tau2))         !-ka

            end do
!======= fin calcul de log integrale
!======================================================================

!======= fin calcul de integrale
!======================================================================

            res = 0.d0
            do k=1,ngexact  !ng!
                if(nig(k).gt.0)then
                    if(indictronq.eq.0)then
                        res = res + integrale1(k)!integrale1 donne le log de I directement
                        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                            funcpaacpm=-1.d9
                            goto 123
                        end if
                    else !troncature
                        indic_tronc = 1
                        funcpaacpm=-1.d9
                        goto 123
                    endif
                endif
            end do
        endif                     !fin boucle effet= 1 and correl = 1
    endif                     !fin boucle globale effet=0

    funcpaacpm = res
    if ((funcpaacpm.ne.funcpaacpm).or.(abs(funcpaacpm).ge. 1.d30)) then
        funcpaacpm=-1.d9
        goto 123
    end if
123     continue
    return

    end function funcpaacpm


!========================    DEBUT FUNCPAO_cpm       ====================
    double precision function funcpaocpm(b,np,id,thi,jd,thj)

    use tailles
    !use comon,only:alpha,nig,nst,t0,
    use comon,only:g,c,t1,nsujet,nva, &
    stra,ve,auxig,nbintervR,ttt,betacoef
    !use additiv,only:rho,ve2,
    use additiv,only:betaaux,sigma2,tau2,cov

    implicit none

    integer::np,id,jd,k,i,gg,jj
    double precision,dimension(np)::bhaux,b
    double precision::thi,thj,res,vet
!****** u_tilde
    double precision::u_tilde,v_tilde,som1,som2
!****** derivanal
    double precision,dimension(ngmax)::res3,res4
    double precision,dimension(ngmax)::res5,res6,res8

!==============================================
!================POUR UN GROUPE AUXIG donne !
!==============================================

    res3=0.d0
    res4=0.d0
    res5=0.d0
    res6=0.d0
    res8=0.d0

    bhaux=b

    if (id.ne.0) bhaux(id)=bhaux(id)+thi 
    if (jd.ne.0) bhaux(jd)=bhaux(jd)+thj

    u_tilde=bhaux(1)
    v_tilde=bhaux(2)

!--------------------------------------------------------
!----------calcul de la vraisemblance ------------------
!---------------------------------------------------------

    do k=1,nsujet
        if(nva.gt.0.and.g(k).eq.auxig)then
            vet = 0.d0 
            do i=1,nva
                vet =vet + betaaux(i)*ve(k,i)
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif

        if(g(k).eq.auxig)then
            if(c(k).eq.1)then
                res3(auxig) = res3(auxig)+u_tilde+v_tilde*ve(k,1)
                res8(auxig) = res8(auxig)+ve(k,1)
            endif
            if(stra(k).eq.1)then
                som1=0.d0
                som2=0.d0
                do gg=1,nbintervR
                    if ((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg)))) then
                        som1=betacoef(gg)*(t1(k)-ttt(gg-1))
                        if (gg.ge.2)then
                            do jj=1,gg-1
                                som2=som2+betacoef(jj)*(ttt(jj)-ttt(jj-1))
                            end do
                        endif
                    end if
                    res4(auxig) =res4(auxig)+(som1+som2)*vet*dexp(u_tilde+v_tilde*ve(k,1))
                    res5(auxig) = res5(auxig)+(som1+som2)*vet*dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)
                    res6(auxig) = res6(auxig)+(som1+som2)*vet*dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
                end do

            endif

            if(stra(k).eq.2)then
                som1=0.d0
                som2=0.d0
                do gg=1,nbintervR
                    if ((t1(k).ge.(ttt(gg-1))).and.(t1(k).lt.(ttt(gg)))) then
                        som1=betacoef(nbintervR+gg)*(t1(k)-ttt(gg-1))
                        if (gg.ge.2)then
                            do jj=1,gg-1
                                som2=som2+betacoef(nbintervR+jj)*(ttt(jj)-ttt(jj-1))
                            end do
                        endif
                    end if
                    res4(auxig) =res4(auxig)+(som1+som2)*vet*dexp(u_tilde+v_tilde*ve(k,1))
                    res5(auxig) = res5(auxig)+(som1+som2)*vet*dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1) 
                    res6(auxig) = res6(auxig)+(som1+som2)*vet*dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
                end do
            endif
        endif
    end do

    res = - res3(auxig)+res4(auxig)+0.5d0*(((u_tilde)**2)/sigma2+((v_tilde)**2)/tau2 &
    -2.d0*u_tilde*v_tilde*cov/(sigma2*tau2))/(1.d0-(cov**2)/(sigma2*tau2))

    funcpaocpm=-res

    return

    end function funcpaocpm




