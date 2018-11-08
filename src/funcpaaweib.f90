!===================================================================
!========================          FUNCPA       ====================
    double precision function funcpaaweib(b,np,id,thi,jd,thj,k0)

    use tailles
    use comon,only:effet,stra,t0,t1,c,nsujet,nva,nst,auxig,ng,ve,g,nig,indictronq, &
    etaR,etaD,betaR,betaD,kkapa,indic_tronc
    !use additiv,only:invD,rho,
    use additiv,only:correl,ngexact,ve2,sigma2,tau2,cov,mid,betaaux
    use residusM,only:som_Xbeta,indic_cumul,cumulhaz

    implicit none

    integer::np,id,jd,i,j,k,cptg,ig,ip
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,som2,res,vet
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
!      common /derivanal/res1,res2,res3,res4,res5,res6,res8

!****** u_tilde
    double precision  :: u_tilde,v_tilde
    double precision,external::funcpaoweib
!      common /utilde/u_tilde,v_tilde
!******************

    kkapa=k0

    j=0
    res=0.d0
    vet=0.d0
    som2=0.d0
    restar = 0
    nf = 1

    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    if (nst == 1) then
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= 0.d0
        etaD= 0.d0
    else
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= bh(3)**2
        etaD= bh(4)**2
    end if

    if(effet.eq.1) then
!terme de correlation entre 2 frailties/avec contraintes
        sigma2 = (bh(np-nva-1)*bh(np-nva-1)) ! variance intercept
        tau2 = (bh(np-nva)*bh(np-nva))  ! variance traitement * groupe
        cov=dsqrt(sigma2*tau2)*(2.d0 * dexp(bh(np-nva-2))/(1.d0+dexp(bh(np-nva-2)))-1.d0)
    endif

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
                res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
            endif  

            if((c(i).eq.1).and.(stra(i).eq.2))then
                res2(g(i)) = res2(g(i))+(betaD-1.d0)*dlog(t1(i))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpaaweib=-1.d9
                goto 123
            end if

            if(stra(i).eq.1)then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaR)**betaR)*vet - ((t0(i)/etaR)**betaR)*vet
                if (indic_cumul==1) then
                    cumulhaz(g(i)) = ((t1(i)/etaR)**betaR)
                end if
            endif

            if(stra(i).eq.2)then
                res1(g(i)) = res1(g(i)) + ((t1(i)/etaD)**betaD)*vet - ((t0(i)/etaD)**betaD)*vet
                if (indic_cumul==1) then
                    cumulhaz(g(i)) = ((t1(i)/etaD)**betaD)
                end if
            endif

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpaaweib=-1.d9
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
                    funcpaaweib=-1.d9
                    goto 123
                end if
            endif
        end do

    else

!*************************************************************
!-----avec deux effets aleatoires dans le modele et correl=0
!*************************************************************
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
                    res2(g(k)) = res2(g(k))+(betaR-1.d0)*dlog(t1(k))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
                endif

                if((c(k).eq.1).and.(stra(k).eq.2))then
                    res2(g(k)) = res2(g(k))+(betaD-1.d0)*dlog(t1(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
                endif

                if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                    funcpaaweib=-1.d9
                    goto 123
                end if

            end do

!=========================================================================
!==== calcul des integrales par transformation de LAPLACE pour chq gpe ===
!=========================================================================

            do ig=1,ng!ngexact
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

    call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux,funcpaoweib)
       if (ieraux .eq.-1) then
          funcpaaweib=-1.d9
          goto 123
       end if
                u_tilde = baux(1)!u_tilde
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
                            res4(ig) = res4(ig)+((t1(k)/etaR)**betaR)*vet* &
                            dexp(u_tilde+v_tilde*ve2(k,1))

                            res5(ig) = res5(ig)+((t1(k)/etaR)**betaR)*vet* &
                            dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)

                            res6(ig) = res6(ig)+((t1(k)/etaR)**betaR)*vet* &
                            dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2

                        endif

                        if(stra(k).eq.2)then
                            res4(ig) = res4(ig)+((t1(k)/etaD)**betaD)*vet* &
                            dexp(u_tilde+v_tilde*ve2(k,1))
                            res5(ig) = res5(ig)+((t1(k)/etaD)**betaD)*vet* &
                            dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                            res6(ig) = res6(ig)+((t1(k)/etaD)**betaD)*vet* &
                            dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2
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
                            funcpaaweib=-1.d9
                            goto 123
                        end if
                    else !troncature
                        indic_tronc = 1
                        funcpaaweib=-1.d9
                        goto 123
                !        write(*,*)'***TRAITER TRONCATURE**'
                    endif
                endif
            end do
        endif                !fin boucle effet= 1 and correl = 0

!=======================================================================

!*************************************************************
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
                    res2(g(k)) = res2(g(k))+(betaR-1.d0)*dlog(t1(k))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
                endif
                if((c(k).eq.1).and.(stra(k).eq.2))then
                    res2(g(k)) = res2(g(k))+(betaD-1.d0)*dlog(t1(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet)
                endif
                if ((res2(g(k)).ne.res2(g(k))).or.(abs(res2(g(k))).ge. 1.d30)) then
                    funcpaaweib=-1.d9
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
    call marq98o(baux,npaux,niaux,vaux,resaux,ieraux,istopaux,funcpaoweib)
       if (ieraux .eq.-1) then
          funcpaaweib=-1.d9
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
                            res3(ig) = res3(ig) &
                            +u_tilde+v_tilde*ve2(k,1)
                        endif

                        if(stra(k).eq.1)then
                            res4(ig) = res4(ig) &
                            +((t1(k)/etaR)**betaR)*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                            res5(ig) = res5(ig) &
                            +((t1(k)/etaR)**betaR)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                            res6(ig) = res6(ig) &
                            +((t1(k)/etaR)**betaR)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2

                        endif
                        if(stra(k).eq.2)then
                            res4(ig) = res4(ig) &
                            +((t1(k)/etaD)**betaD)*vet*dexp(u_tilde+v_tilde*ve2(k,1))
                            res5(ig) = res5(ig) &
                            +((t1(k)/etaD)**betaD)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*ve2(k,1)
                            res6(ig) = res6(ig) &
                            +((t1(k)/etaD)**betaD)*vet*dexp(u_tilde+v_tilde*ve2(k,1))*(ve2(k,1))**2

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

            res = 0.d0
            do k=1,ngexact  !ng!
                if(nig(k).gt.0)then
                    if(indictronq.eq.0)then
                        res = res+integrale1(k)!integrale1 donne le log de I directement
                        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                            funcpaaweib=-1.d9
                            goto 123
                        end if
                    else !troncature
                        indic_tronc = 1
                        funcpaaweib=-1.d9
                        goto 123
                    endif
                endif
            end do

        endif                     !fin boucle effet= 1 and correl = 1
    endif                     !fin boucle globale effet=0 

    funcpaaweib = res
    if ((funcpaaweib.ne.funcpaaweib).or.(abs(funcpaaweib).ge. 1.d30)) then
        funcpaaweib=-1.d9
        goto 123
    end if
123     continue
    return

    end function funcpaaweib


!========================    DEBUT FUNCPAO_cpm       ====================
    double precision function funcpaoweib(b,np,id,thi,jd,thj)

    use tailles
    !use comon,only:ndate,nig,nst,t0
    use comon,only:g,t1,c,nsujet,nva, &
    stra,ve,auxig,etaR,etaD,betaR,betaD
    !use additiv,only:rho,ve2
    use additiv,only:betaaux,sigma2,tau2,cov

    implicit none

    integer::np,id,jd,i,k,ip
    double precision,dimension(np)::bhaux,b
    double precision::thi,thj,res,vet
!****** u_tilde
    double precision  :: u_tilde,v_tilde
!****** derivanal
    double precision , dimension(ngmax)::res3,res4
    double precision , dimension(ngmax)::res5,res6,res8

!==============================================
!================POUR UN GROUPE AUXIG donne !
!==============================================

    res3=0.d0
    res4=0.d0
    res5=0.d0
    res6=0.d0
    res8=0.d0

    do i=1,np
        bhaux(i)=b(i)
    end do

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
            do ip=1,nva
                vet =vet + betaaux(ip)*ve(k,ip)
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
                res4(auxig) =res4(auxig)+((t1(k)/etaR)**betaR)*vet*dexp(u_tilde+v_tilde*ve(k,1))

                res5(auxig) = res5(auxig)+((t1(k)/etaR)**betaR)*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)

                res6(auxig) = res6(auxig)+((t1(k)/etaR)**betaR)*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
            endif
            if(stra(k).eq.2)then
                res4(auxig) = res4(auxig) &
                +((t1(k)/etaD)**betaD)*vet*dexp(u_tilde+v_tilde*ve(k,1))

                res5(auxig) = res5(auxig)+((t1(k)/etaD)**betaD)*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)

                res6(auxig) = res6(auxig)+((t1(k)/etaD)**betaD)*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
            endif

        endif
    end do

    res=-res3(auxig) + res4(auxig)+0.5d0*(((u_tilde)**2)/sigma2+((v_tilde)**2)/tau2 &
    -2.d0*u_tilde*v_tilde*cov/(sigma2*tau2))/(1.d0-(cov**2)/(sigma2*tau2))         !-ka

    funcpaoweib= -res

    return

    end function funcpaoweib
