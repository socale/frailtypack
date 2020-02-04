


!========================          FUNCPAJ_WEIBULL FOR JOINT NESTED         ====================
    double precision function funcpajweib_fam(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:AG,auxif,auxig,betaR,etaR,fsize,fam,res4,nst,t0dc
    use comon,only: etaD,betaD, etaT, betaT, nstRec, &
    t0,t1,t1dc,c,cdc,nsujet,nva,nva1,nva2,& 
    effet,stra,ve,vedc,ng,g,nig,indic_ALPHA,ALPHA,theta,& 
    nfam,fam,fsize,indic_xi,xi,eta, & !for family 
    aux1,aux2,res1,res3,kkapa,nb_gl
    use residusM
    !use comongroup,only:the1,the2
    use comongroup,only:vet,vet2

    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :

    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    integer::n,i,j,k,vj,ig,choix,jj
    integer,dimension(ngmax)::cpt
    double precision :: sum, res
    
    double precision,dimension(np)::bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc &
    ,res3dc,integrale1,integrale2,integrale3
    double precision,dimension(nfam)::integrale3f
    double precision::int

    kkapa=k0
    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0
    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj
    
    do jj=1,nstRec !en plus strates A.Lafourcade 07/2014
        betaT(jj)=bh((jj-1)*2+1)**2     
    etaT(jj)= bh((jj-1)*2+2)**2
    end do
    betaD= bh(2*nstRec+1)**2
    etaD= bh(2*nstRec+2)**2

    if(effet.ge.1) then
        theta = bh(np-nva-indic_ALPHA-indic_xi-1)*bh(np-nva-indic_ALPHA-indic_xi-1)
        eta = bh(np-nva-indic_ALPHA-indic_xi)*bh(np-nva-indic_ALPHA-indic_xi)
        if (indic_alpha.eq.1) then ! new : joint more flexible alpha = 1 
            alpha = bh(np-nva-indic_xi)
        else
            alpha = 1.d0
        endif
        if (indic_xi.eq.1) then ! new : joint more flexible xi = 1 
            xi = bh(np-nva)
        else
            xi = 1.d0
        endif
    endif

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
        cpt(k) = 0
        integrale1(k) = 0.d0
        integrale2(k) = 0.d0
        integrale3(k) = 0.d0
        aux1(k)=0.d0
        aux2(k)=0.d0
    end do

    do i=1,nfam
        integrale3f(i) = 0.d0
    end do

!*********************************************
!-----avec un effet aleatoire dans le modele
!*********************************************

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

    do i=1,nsujet
        cpt(g(i))=cpt(g(i))+1
        if(nva1.gt.0)then
            vet = 0.d0
            do j=1,nva1
                vet = vet + bh(np-nva+j)*dble(ve(i,j))
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif

        if((c(i).eq.1))then
        res2(g(i)) = res2(g(i))+(betaT(stra(i))-1.d0)*dlog(t1(i))+ &
            dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpajweib_fam=-1.d9
                goto 123
            end if        
        endif
        
        res1(g(i)) = res1(g(i))+((t1(i)/etaT(stra(i)))**betaT(stra(i)))*vet
        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
            funcpajweib_fam=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature:
        res3(g(i)) = res3(g(i))+((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
            funcpajweib_fam=-1.d9
            goto 123
        end if    
    end do

!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces
!ccccccccccccccccccccccccccccccccccccccccc

    do k=1,ng ! dans Joint ng=nb individus
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2=1.d0
        endif
        if(cdc(k).eq.1)then
            res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2)
            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                funcpajweib_fam=-1.d9
                goto 123
            end if
        endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
        aux1(k)=((t1dc(k)/etaD)**betaD)*vet2

        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            funcpajweib_fam=-1.d9
            goto 123
        end if
    end do

    choix = 3
    call gaulagJf(int,nb_gl)
    res=int
    
!************* FIN INTEGRALES **************************
    do k=1,ng
        sum=0.d0
        if(cpt(k).gt.0)then        
            res= res + res2(k) + res2dc(k)
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajweib_fam=-1.d9
                goto 123
            end if
        endif
    end do

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpajweib_fam=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123
    else
        funcpajweib_fam = res
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
        end do
        
        k=0
        do i= 1,ng
            do j=1,fsize(fam(i))
                cumulhaz1(i,j) = res1(k+j)
                cumulhaz0(i,j) = res3(k+j)
                cumulhazdc(i,j) = aux1(k+j)
            end do
            if(i.lt.ng.and.fam(i).ne.fam(i+1)) then 
                k = k +fsize(fam(i))
            end if
        end do
    end if
!Ad:
123     continue
    return

    end function funcpajweib_fam