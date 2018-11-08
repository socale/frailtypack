    

    
    module propred         
        integer,dimension(:),allocatable,save::stracross!pour crossvalidation   
        double precision,dimension(:),allocatable,save::aux
        double precision,dimension(:,:),allocatable,save::y
        double precision,dimension(:),allocatable,save::v        
        double precision,dimension(:,:),allocatable,save::I1_hess,H1_hess
        double precision,dimension(:,:),allocatable,save::I2_hess,H2_hess
        double precision,dimension(:,:),allocatable,save::HI1,HI2
        double precision,dimension(:,:),allocatable,save::HIH,IH,HI    
    end module propred
    
!
!
!#########################################################################################################
!
!
    subroutine additive(ns0,ng0,nst0,nz0,xmin10,xmin20,tt00,tt10,ic0,groupe0,nva0, &
    str0,vax0,interaction,ag0,noVar,maxiter0,irep10,correl0,np,b,coef,varcoef,varcoef2, &
    rhoEnd,covEnd,varcovEnd,varSigma2,varTau2,ni,res,LCV,k0,x1Out,lamOut,xSu1,suOut,x2Out, &
    lam2Out,xSu2,su2Out,typeof0,equidistant,nbintervR0,mt,ier,ddl,istop,shapeweib,scaleweib,mt1,trunc,ziOut,time,&
    Resmartingale,frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtycov,linearpred,EPS)

    use parameters        
    use tailles
    use propred
    use comon
    use additiv
    use optim
    use residusM
    
    implicit none    
     
    integer::noVar,interaction,mt,mt1,trunc
    integer::groupe,j,k,nz,n,np,cpt,ii,iii,iiii,ver, &
    cptstr1,cptstr2,i,ic,ni,ier,istop,l,str, &
        nb_echec,nb_echecor,auxng, &
    irep1,nvacross,nstcross,effetcross
!----------------- ajout      
    integer::ns0,ng0,nst0,correl0,maxiter0,nva0,ag0,nz0,irep10
    double precision,dimension(nz0+6),intent(out)::ziOut
    double precision,dimension(ns0)::tt00,tt10
    integer,dimension(ns0)::groupe0,ic0,str0
    double precision, dimension(ns0,nva0)::vax0
    double precision::xmin10,xmin20
    double precision,dimension(np)::b
!------------------
    double precision::tt0,tt1
    integer,dimension(nva0)::filtre,filtre2
    double precision::xmin1,xmin2,res,min,max,maxt,pord, &
    lrs,trace,trace1,trace2,auxi,ax,bx,cx,tol,ddl, &
    fa,fb,fc,goldenadd,estimvadd,f1,f2,f3,varcov    
    double precision,dimension(2):: auxkappa,k0 
    double precision,dimension(:),allocatable::vax
    double precision::rhoEnd,covEnd,varcovEnd
    double precision,dimension(nva0)::coef,varcoef,varcoef2    
    double precision,dimension(2)::varSigma2,varTau2
    double precision,dimension(mt)::x1Out,x2Out
    double precision,dimension(mt,3)::lamOut,lam2Out
    double precision,dimension(mt1,3)::suOut,su2Out
    double precision,dimension(mt1)::xSu1,xSu2
!AD: add traceLCV    
    double precision,dimension(2),intent(out)::LCV,shapeweib,scaleweib
!AD: add for new marq
    double precision::ca,cb,dd,funcpaasplines,funcpaacpm,funcpaaweib
    external::funcpaasplines,funcpaacpm,funcpaaweib

!Cpm
    integer::typeof0,nbintervR0,equidistant,ent,indd
    double precision::temp
    double precision,dimension(nbintervR0+1)::time
    integer,dimension(3)::istopp
    double precision,dimension(ng0)::frailtysd,frailtysd2
!predictor
    double precision,dimension(ng0),intent(out)::Resmartingale,frailtypred,frailtypred2,frailtyvar,&
    frailtyvar2,frailtycov
    double precision,external::funcpaares
    double precision,dimension(ns0),intent(out)::linearpred
    double precision,dimension(1,ns0)::XBeta1
    double precision,dimension(1,nva0)::Xcoef
    double precision,dimension(2,2)::sigma
    double precision,dimension(3),intent(inout)::EPS ! seuils de convergence
        
!cpm
    indic_cumul=0
    istopp=0
    ca=0.d0
    cb=0.d0
    dd=0.d0
    filtre=0
    filtre2=0
!AD:end    
    epsa=EPS(1) !1.d-4
    epsb=EPS(2) !1.d-4
    epsd=EPS(3) !1.d-3
    typeof = typeof0
    model=2
!------------------------
    if (typeof .ne. 0) then
        nbintervR = nbintervR0
    end if    
    shapeweib = 0.d0
    scaleweib = 0.d0
    auxng=0
    maxiter=maxiter0
!----------------
    str=0
!---------------- 
    allocate(invd(2,2))
     
!=== initialisations 
        lrs=0.d0    
        nb_echec=0 
        nb_echecor=0

!=== fin initialisations 
    
    xmin2=xmin20
    xmin1=xmin10
!-----------------------------------------------------------------------------------------
    
    nsujet=ns0
    nsujetmax=nsujet
    allocate(t0(nsujetmax),t1(nsujetmax),c(nsujetmax),nt0(nsujetmax),nt1(nsujetmax), &
    stra(nsujetmax),g(nsujetmax),stracross(nsujetmax),aux(2*nsujetmax))

    ndatemax=2*nsujet
    allocate(date(ndatemax))
    if (typeof == 0) then
        allocate(mm3(ndatemax),mm2(ndatemax),mm1(ndatemax),mm(ndatemax), &
        im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax),dut1(ndatemax),dut2(ndatemax), &
        ut1(0:ndatemax),ut2(0:ndatemax))
    end if
!-----------------------------------------------------------------------------------------
    ng=ng0
    ngmax=ng
    allocate(nig(ngmax),mid(ngmax))    


    allocate(cumulhaz1(ngmax,2))
!-----------------------------------------------------------------------------------------
    nst=nst0
    correl=correl0  
    correlini=correl
!-----------------------------------------------------------------------------------------
    ver=nva0
    nvarmax=ver
    allocate(ve(nsujetmax,nvarmax),ve2(nsujetmax,nvarmax),betaaux(nvarmax))
    allocate(som_Xbeta(ngmax))
    allocate(vax(nvarmax))

    if (noVar.eq.1) then 
        filtre=0
        nva=0  
    else
        filtre=1
        filtre2(interaction)=1
        nva=nva0  
    end if        
    ag=ag0
    

    nz=nz0
    
    effet = 1
    nig = 0
    g=0  
!------------  lecture fichier -----------------------
    maxt = 0.d0
    mint = 0.d0
    cpt = 0
    k = 0
    cptstr1 = 0
    cptstr2 = 0
    
    do i = 1,nsujet

        if (i.eq.1) then
            mint = tt00(i) ! affectation du min juste une fois
        endif

        if(nst.eq.2)then
            tt0=tt00(i)
            tt1=tt10(i)
            groupe=groupe0(i)
            ic=ic0(i)
            str=str0(i)
            do j=1,nva
                vax(j)=vax0(i,j)
            end do
        else
            tt0=tt00(i)
            tt1=tt10(i)
            groupe=groupe0(i)
            ic=ic0(i)
            do j=1,nva
                vax(j)=vax0(i,j)
            end do            
        endif
        k = k +1
        if(k.eq.1)then
            auxng=groupe
            ngexact=1
            g(k)=1
        else
            g(k)=groupe
        endif
!------------------   observation c=1
        if(ic.eq.1)then
            cpt = cpt + 1
            c(k)=1
            if(str.eq.1.and.nst.eq.2)then
                stra(k) = 1
                cptstr1 = cptstr1 + 1
            endif
            if(str.eq.2.and.nst.eq.2)then
                stra(k) = 2
                cptstr2 = cptstr2 + 1
            endif

            if(nst.eq.1)then
                stra(k) = 1
                cptstr1 = cptstr1 + 1
            endif
            t0(k) = tt0
            t1(k) = tt1

            if(auxng.ne.groupe)then !chgt de groupes
                ngexact=ngexact+1
                auxng=groupe

                if(k.ne.1)then
                    g(k)=g(k-1)+1
                endif
                goto 100             
            endif
  
100  continue

            nig(g(k)) = nig(g(k))+1
            iii = 0
            iiii = 0
            do ii = 1,ver
                if(filtre(ii).eq.1)then
                    iii = iii + 1
                    ve(i,iii) = dble(vax(ii))
                    ve2(i,iii) = ve(i,iii) 
!====================================================================
                endif
                if(filtre2(ii).eq.1)then
                    iiii = iiii + 1
                    ve(k,iiii) = dble(vax(ii))
                    ve2(k,iiii) = ve(k,iiii)  
                endif
!====================================================================
            end do   
        else 
!------------------   censure a droite  c=0
            if(ic.eq.0)then
                c(k) = 0 
                if(str.eq.1.and.nst.eq.2)then
                    stra(k) = 1
                    cptstr1 = cptstr1 + 1
                endif
                if(str.eq.2.and.nst.eq.2)then
                    stra(k) = 2
                    cptstr2 = cptstr2 + 1
                endif
                if(nst.eq.1)then
                    stra(k) = 1
                    cptstr1 = cptstr1 + 1
                endif
                iii = 0
                iiii = 0
                do ii = 1,ver
                    if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(k,iii) = dble(vax(ii))
                        ve2(k,iii) =  ve(k,iii) 
!==================recodage en +-1/2 ===============================
                    endif
                    if(filtre2(ii).eq.1)then
                        iiii = iiii + 1
                        ve(k,iiii) = dble(vax(ii))
                        ve2(k,iiii) =  ve(k,iiii) 
                    endif
                end do 

                t0(k) =  tt0
                t1(k) = tt1
                if(auxng.ne.groupe)then !chgt de groupes
                    ngexact=ngexact+1
                    auxng=groupe
                    if(k.ne.1)then
                        g(k)=g(k-1)+1
                    endif
                    goto 101
                endif   
101   continue
                nig(g(k)) = nig(g(k))+1
            endif
        endif
        if (maxt.lt.t1(k))then
            maxt = t1(k)
        endif
        if (mint.gt.t0(k)) then
            mint = t0(k)
        endif
    end do 
    
!AD:
    if (typeof .ne. 0) then 
        cens = maxt
    end if
!Ad    
        
    if (typeof == 0) then    
        nz1=nz
        nz2=nz
        if(nz.gt.20)then
            nz = 20
        endif
        if(nz.lt.4)then
            nz = 4
        endif
    end if

!***************************************************
!--------------- zi- ----------------------------------

!      construire vecteur zi (des noeuds)

    
    min = 1.d-10
    max = maxt
    
    do i = 1,2*nsujet
        do k = 1,nsujet
            if((t0(k).ge.min))then
                if(t0(k).lt.max)then
                    max = t0(k)
                endif
            endif
            if((t1(k).ge.min))then
                if(t1(k).lt.max)then
                    max = t1(k)
                endif
            endif
        end do   
        aux(i) = max
        min = max + 1.d-12
        max = maxt
    end do
    
    date(1) = aux(1)
    k = 1
    do i=2,2*nsujet
        if(aux(i).gt.aux(i-1))then
            k = k+1
            date(k) = aux(i)
        endif 
    end do 
    
    if(typeof == 0) then
    

! Al:10/03/2014 emplacement des noeuds splines en percentile (sans censure par intervalle)
        i=0
        j=0
!----------> taille - nb de recu
        do i=1,nsujet
            if(t1(i).ne.(0.d0).and.c(i).eq.1) then
                j=j+1
            endif
        end do
        nbrecu=j

!----------> allocation des vecteur temps
        allocate(t2(nbrecu))
        
!----------> remplissage du vecteur de temps
        j=0
        do i=1,nsujet
            if (t1(i).ne.(0.d0).and.c(i).eq.1) then
                j=j+1
                t2(j)=t1(i)
            endif
        end do
            
        nzmax=nz+3
        allocate(zi(-2:nzmax))
        ndate = k

        zi(-2) = mint
        zi(-1) = mint !date(1)
        zi(0) = mint !date(1)
        zi(1) = mint !date(1) 
        j=0
        do j=1,nz-2
            pord = dble(j)/(dble(nz)-1.d0)
            call percentile3(t2,nbrecu,pord,zi(j+1))
        end do
        zi(nz) = maxt !date(ndate)
        zi(nz+1) = maxt !zi(nz)
        zi(nz+2) = maxt !zi(nz)
        zi(nz+3) = maxt !zi(nz)
        ziOut = zi
        deallocate(t2)
!         nzmax=nz+3
!         allocate(zi(-2:nzmax))
!     
!         ndate = k
! 
!         zi(-2) = date(1)
!         zi(-1) = date(1)
!         zi(0) = date(1)
!         zi(1) = date(1)
!         h = (date(ndate)-date(1))/dble(nz-1)
!         do i=2,nz-1
!             zi(i) =zi(i-1) + h   
!         end do
! 
!         zi(nz) = date(ndate)
!         zi(nz+1)=zi(nz)
!         zi(nz+2)=zi(nz)
!         zi(nz+3)=zi(nz)
!         ziOut = zi


    end if
!---------- affectation nt0,nt1----------------------------

    indictronq=0
    do i=1,nsujet 
        if (typeof == 0) then
            if(t0(i).eq.0.d0)then
                nt0(i) = 0
            endif
        end if

        if(t0(i).ne.0.d0)then
            indictronq=1
        endif

        if (typeof == 0) then
            do j=1,ndate
                if(date(j).eq.t0(i))then
                    nt0(i)=j
                endif
                if(date(j).eq.t1(i))then
                    nt1(i)=j
                endif
            end do
        end if

    end do 

!---------- affectation des vecteurs de splines -----------------
    if (typeof == 0) then
        n  = nz+2
    
        call vecspliadd(n,ndate)
    
        
        allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax), &
        m3m1(nzmax),m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax))
        
        call vecpenadd(n)  
            
        np = nst*n + nva + correl + 2*effet 
    end if 
!-----------------------
!    if(typeof == 1) then 
!        np = nst*nbintervR + nva + correl + 2*effet 
!    end if
!    
!    if(typeof == 2) then
!        np = nst*2 + nva + correl + 2*effet 
!    end if
    
    npmax=np
    
    allocate(b_temp(np))
    allocate(I_hess(npmax,npmax),H_hess(npmax,npmax),Hspl_hess(npmax,npmax) &
    ,hess(npmax,npmax),PEN_deri(npmax,1))
    allocate(y(npmax,npmax),v((npmax*(npmax+3)/2)),&
    I1_hess(npmax,npmax),H1_hess(npmax,npmax),I2_hess(npmax,npmax),H2_hess(npmax,npmax), &
    HI1(npmax,npmax),HI2(npmax,npmax),HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax))
!---------------------    
! term de correlation entre frailty + 2 termes de var de la meme frailty   
    nbpara =np
            
!------- initialisation des parametres                  
    do i=1,np
        b(i)= 1.d-1
    end do
! traitement, intercept et pente initialis�s plus loin
    if (typeof == 1) then
!------- RECHERCHE DES NOEUDS
!----------> Enlever les zeros dans le vecteur de temps
        i=0
        j=0

!----------> taille - nb de recu
        do i=1,nsujet
            if(t1(i).ne.(0.d0).and.c(i).eq.1) then
                j=j+1
            endif
        end do
        nbrecu=j

    !    n = nbintervR
!----------> allocation des vecteur temps
        allocate(t2(nbrecu))
        
!----------> remplissage du vecteur de temps
        j=0
        do i=1,nsujet
            if (t1(i).ne.(0.d0).and.c(i).eq.1) then
                j=j+1
                t2(j)=t1(i)
            endif
        end do
        
!----------> tri du vecteur de temps
        indd=1
        do while (indd.eq.1)
            indd=0
            do i=1,nbrecu-1
                if (t2(i).gt.t2(i+1)) then
                    temp=t2(i)
                    t2(i)=t2(i+1)
                    t2(i+1)=temp
                    indd=1        
                end if
            end do
        end do    
        
        ent=int(nbrecu/(nbintervR))
        
        allocate(ttt(0:nbintervR))
        
        ttt(0)=0.d0
        
        ttt(nbintervR)=cens
        
        j=0
        do j=1,nbintervR-1
            if (equidistant.eq.0) then
                ttt(j)=(t2(ent*j)+t2(ent*j+1))/(2.d0)
            else
                ttt(j)=(cens/(nbintervR))*j
            endif
        end do
        time = ttt
        deallocate(t2)
!------- FIN RECHERCHE DES NOEUDS    
    end if    

!***********************************************************
!************** NEW : cross validation  ***********************
!!!!! sur une seule strate, sans var expli , sans frailties ****
!*****************************************************************
    if (typeof == 0) then
        nvacross=nva !pour la recherche du parametre de lissage sans var expli
        nva=0
        effetcross=effet
        effet=0
        nstcross=nst
        nst=1
        correl=0

        do l=1,nsujet  
            stracross(l)=stra(l)
        end do
        do l=1,nsujet  
            stra(l)=1
        end do

    
        irep1=irep10

!        if(irep1.eq.0.and.nst.ge.2)then
! ne se produit jamais maintenant (1 seule strate pour CV): changement de juin 2009
!            stop
!        endif
    
        xmin1=xmin10
        
        if(xmin1.le.0.d0)then
            xmin1 = 0.d0
        endif  

    
!*************************************************
!    nvat=nva
        nva=0
    end if


!on  travaille d'abord sur une seule strate
    if(typeof == 0) then
        if(irep1.eq.1)then   !pas recherche du parametre de lissage
            xmin1 = dsqrt(xmin1)

            auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

            if (ni.ge.250) then

                do i=1,nz+2
                    b(i)=1.d-1
                end do     
                xmin1 = sqrt(10.d0)*xmin1
                auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

                if (ni.lt.250) then

                else
                    do i=1,nz+2
                        b(i)=1.d-1
                    end do     
                    xmin1 = sqrt(10.d0)*xmin1
                    auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)
     
                endif
            else

            endif
    !----------------------------------------------------
        else                   !recherche du parametre de lissage

            if(xmin1.le.0.d0)then
                xmin1 = 1.d0
            endif  

            xmin1 = dsqrt(xmin1)
            auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

            if(ddl.gt.-2.5d0)then
                xmin1 = dsqrt(xmin1)
                auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)
    
                if(ddl.gt.-2.5d0)then
                    xmin1 = dsqrt(xmin1)
                    auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)


                    if(ddl.gt.-2.5d0)then
                        xmin1 = dsqrt(xmin1)
                        auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

                        if(ddl.gt.-2.5d0)then
                            xmin1 = dsqrt(xmin1)
                            auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

                            if(ddl.gt.-2.5d0)then
                                xmin1 = dsqrt(xmin1)
                            endif   
                        endif   
                    endif   
                endif
            endif 

            if (ni.ge.250) then
                do i=1,nz+2
                    b(i)=1.d-1
                end do     
                xmin1 = sqrt(10.d0)*xmin1
                auxi = estimvadd(xmin1,n,b,y,ddl,ni,res)

                if (ni.ge.250) then
                    do i=1,nz+2
                        b(i)=1.d-1
                    end do     
                    xmin1 = sqrt(10.d0)*xmin1
                endif
            endif 
            ax = xmin1
            bx = xmin1*dsqrt(1.5d0)  

            call mnbrakadd(ax,bx,cx,fa,fb,fc,b,n)            
            tol = 0.001d0
            res = goldenadd(ax,bx,cx,tol,xmin1,n,b,y,ddl)

            effet=0
            correl=0

            auxkappa(1)=xmin1*xmin1
            auxkappa(2)=0.d0
            call marq98J(auxkappa,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaasplines)
            if (istop .ne. 1) then
                istopp(1)=1
                goto 1000
            end if    
        end if
    end if
    
    
    if(typeof == 0) then
        nva=nvacross ! pour la recherche des parametres de regression
        nst=nstcross ! avec stratification si n�cessaire
        effet=effetcross ! avec effet initial
        do l=1,nsujet  
            stra(l)=stracross(l) !r�tablissement stratification
        end do
    
!********************************************************************
        if(nst.eq.2)then
            xmin2=xmin20
        endif
    
        k0(1) = xmin1*xmin1
        k0(2) = 0.d0
        
        if(nst.eq.2)then
            k0(2) = xmin2
        endif
    end if

    if (typeof .ne. 0) then
        allocate(vvv((npmax*(npmax+1)/2))) !,kkapa(2))
    end if    
    
!=============================- fin cross validation
    
!===== initialisation des parametres de regression/pas effets aleatopires
!    write(*,*)'====================================='
!    write(*,*)'== avec var explicatives============='
!    write(*,*)'==================================='
    
    effet=0

!    write(*,*)'===avant marq98==========',n,np,nva,nst,effet
    indic_cumul=1
    allocate(cumulhaz(ngmax))

    select case(typeof)
        case(0)
            np = nst*n + nva
            b(nst*n+1)=-0.15d0 !initialisation traitementc
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaasplines)
        case(1)
            np = nst*nbintervR + nva
            b(nst*nbintervR+1)=-0.15d0
            allocate(betacoef(nst*nbintervR))
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaacpm)
            
        case(2)
            np = nst*2 + nva
            b(nst*2+1)=-0.15d0
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaaweib)
    end select
    if (istop .ne. 1) then
        istopp(2)=1
        goto 1000
    end if            
!===== recherche de l'ensemble des parametres
    
!    write(*,*)'====================================='
!    write(*,*)'== ensemble des parametres ======='
!    write(*,*)'====================================='
    
    effet=1
    correl=correlini
    do i=1,nva
        b(np-i+2+correl+1)=b(np-i+1)
    end do
    np=nbpara 
    
    if(correl.eq.1)then
        b(np-nva-2)=0.1d0!-0.69d0!-1.61d0!-0.64d0!initialisation cov avec contrainte 
    endif
    

    b(np-nva-1)=0.5d0!0.5477d0      !      initialisation intercept
    b(np-nva)=0.15d0!0.5477d0  !0.13d0!  !   0.15d0!   initialisation pente
    
    
!    write(*,*)' ============> marq2 typeof :',typeof
!    write(*,*)'effet ',effet,' correl ',correl
    indic_cumul=0
    select case(typeof)
        case(0)
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaasplines)
        case(1)            
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaacpm)            
        case(2)
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaaweib)
    end select
    !    write(*,*)'fin marquardt'
!Al:
    EPS(1) = ca
    EPS(2) = cb
    EPS(3) = dd
!Al:
    if (istop .ne. 1) then
        istopp(3)=1
        goto 1000
    end if
    
    j=(np-nva)*(np-nva+1)/2
    
    trace=0
    trace1=0
    trace2=0
    
    select case(typeof)
        case(0)    
! que lorsque la matrice totale est inversible, ie sans echec inversion
! strate1 : 
            do i=1,nz1+2
                do j=1,nz1+2
                    H1_hess(i,j)=H_hess(i,j)
                    I1_hess(i,j)=I_hess(i,j)
                end do
            end do 
            call multiJ(H1_hess,I1_hess,nz1+2,nz1+2,nz1+2,HI1)
            do i =1,nz1+2
                trace1=trace1+HI1(i,i)
            end do

! strate2 :
            if(nst.eq.2)then
                do i=1,nz2+2
                    k=nz1+2+i
                    do j=1,nz2+2
                        l=nz1+2+j                   
                        H2_hess(i,j)=H_hess(k,l)
                        I2_hess(i,j)=I_hess(k,l)
                    end do
                end do
                call multiJ(H2_hess,I2_hess,nz2+2,nz2+2,nz2+2,HI2)
                do i =1,nz2+2
                    trace2=trace2+HI2(i,i)
                end do
            endif ! pour la deuxieme strate
            
        case(1)
            
            do i=1,nbintervR
                do j=1,nbintervR
                    H1_hess(i,j)=H_hess(i,j)
                    I1_hess(i,j)=I_hess(i,j)
                end do
            end do 
            call multiJ(H1_hess,I1_hess,nbintervR,nbintervR,nbintervR,HI1)
            do i =1,nbintervR
                trace1=trace1+HI1(i,i)
            end do

! strate2 :
            if(nst.eq.2)then
                do i=1,nbintervR
                    k=nbintervR+i
                    do j=1,nbintervR
                        l=nbintervR+j                   
                        H2_hess(i,j)=H_hess(k,l)
                        I2_hess(i,j)=I_hess(k,l)
                    end do
                end do
                call multiJ(H2_hess,I2_hess,nbintervR,nbintervR,nbintervR,HI2)
                do i =1,nbintervR
                    trace2=trace2+HI2(i,i)
                end do
            endif ! pour la deuxieme strate
            
        case(2)
        
            do i=1,2
                do j=1,2
                    H1_hess(i,j)=H_hess(i,j)
                    I1_hess(i,j)=I_hess(i,j)
                end do
            end do 
            call multiJ(H1_hess,I1_hess,2,2,2,HI1)
            do i =1,2
                trace1=trace1+HI1(i,i)
            end do

! strate2 :
            if(nst.eq.2)then
                do i=1,2
                    k=2+i
                    do j=1,2
                        l=2+j                   
                        H2_hess(i,j)=H_hess(k,l)
                        I2_hess(i,j)=I_hess(k,l)
                    end do
                end do
                call multiJ(H2_hess,I2_hess,2,2,2,HI2)
                do i =1,2
                    trace2=trace2+HI2(i,i)
                end do
            endif ! pour la deuxieme strate
    end select
    
    call multiJ(I_hess,H_hess,np,np,np,IH)
    call multiJ(H_hess,IH,np,np,np,HIH)
! Salida Juan Aug'07
    
    f1 = (b(np-nva))*(2.d0* dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2)))-1.d0)
    f2 = (b(np-nva-1))*(2.d0* dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2)))-1.d0)
    f3= 2.d0*(b(np-nva-1))*(b(np-nva))*(dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2)))- &
    (dexp(b(np-nva-2))/(1.d0+dexp(b(np-nva-2))))**2)

    varcov=f1*f1*H_hess(np-nva-1,np-nva-1)+f2*f2*H_hess(np-nva,np-nva)+ &
    f3*f3*H_hess(np-nva-2,np-nva-2)+2.d0*f1*f3*H_hess(np-nva-1,np-nva-2)+ &
    2.d0*f2*f3*H_hess(np-nva,np-nva-2)+2.d0*f1*f2*H_hess(np-nva-1,np-nva)
     
    
    if (correl.eq.1)then
        rho=cov/(dsqrt(b(np-nva-1)*b(np-nva-1)*b(np-nva)*b(np-nva)))
    else
        rho=-1
    end if

    
    rhoEnd=rho
    covEnd=cov
    varcovEnd=varcov



    varSigma2(1)=((2.d0*b(np-nva-1))**2)*H_hess(np-nva-1,np-nva-1)
    varSigma2(2)=((2.d0*b(np-nva-1))**2)*HIH(np-nva-1,np-nva-1)
    
    varTau2(1)=((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva)
    varTau2(2)=((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva)


! Covariates

    do i=1,nva
        coef(i)=b(np-nva+i)
        varcoef(i)=H_hess(np-nva+i,np-nva+i) 
        varcoef2(i)=HIH(np-nva+i,np-nva+i)
    end do

    
    if(effet.eq.1.and.ier.eq.-1)then
        v((np-nva)*(np-nva+1)/2)=10.d10
    endif
        
    j=(np-nva)*(np-nva+1)/2
    b_temp = b
    
! --------------  Lambda and survival and cumulative hazard estimates 
    
!    write(*,*)'===========> avant distance <==========='
    select case(typeof)
        case(0)
            call distanceasplines(nz1,nz2,b,effet,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
        case(1)
            Call distancecpm(b,nst*nbintervR,mt,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
        case(2)
            if (nst == 1) then
                typeof2 = 1
            else
                typeof2 = 2
            end if
            Call distanceweib(b,np,mt,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
    end select
        
!    write(*,*)'===========> apres distance <==========='
    if (nst == 1) then
        scaleweib(1) = etaR !betaR
        shapeweib(1) = betaR !etaR
        scaleweib(2) = 0.d0
        shapeweib(2) = 0.d0
    else
        scaleweib(1) = etaR !betaR
        shapeweib(1) = betaR !etaR
        scaleweib(2) = etaD !betaD
        shapeweib(2) = betaD !etaD
    end if    
!AD:add LCV
!     calcul de la trace, pour le LCV (likelihood cross validation)
    LCV=0.d0
    if (typeof == 0) then    
!        write(*,*)'The approximate like cross-validation Criterion in the non parametric case'
        call multiJ(H_hess,I_hess,np,np,np,HI)
        do i =1,np
            LCV(1) = LCV(1) + HI(i,i)
        end do 
        LCV(1) = (LCV(1) - resnonpen) / nsujet
    else        
!        write(*,*)'=========> Akaike information Criterion <========='
!        LCV(2) = 2.d0 * np - 2.d0 * res
        LCV(2) = (1.d0 / nsujet) * (np - res)
!        write(*,*)'======== AIC :',LCV(2)
    end if
    
!AD:end     
1000    continue
    trunc = indic_tronc

    Resmartingale=0.d0
    frailtypred=0.d0
    frailtyvar=0.d0
    frailtysd=0.d0
    frailtypred2=0.d0
    frailtyvar2=0.d0
    frailtysd2=0.d0
    frailtycov=0.d0    

    if(nva .gt. 0) then
!        write(*,*)'========= Call martingale =========='
        allocate(Xbeta(1,ns0),invsigma(2,2))
        sigma(1,1)=b(np-nva-1)**2
        sigma(2,2)=b(np-nva)**2    

        if (correl == 1) then
            sigma(1,2)=dsqrt(sigma(1,1))*dsqrt(sigma(2,2))*rhoEnd
        else
            sigma(1,2)=0.d0
        end if
        sigma(2,1)=sigma(1,2)

        detsigma=sigma(1,1)*sigma(2,2)-sigma(1,2)**2
        invsigma(1,1)=sigma(2,2)
        invsigma(2,2)=sigma(1,1)
        invsigma(1,2)=-sigma(1,2)
        invsigma(2,1)=-sigma(2,1)

        do i= 1,2
            do j=1,2
                invsigma(i,j)= 1.d0/detsigma*invsigma(i,j)
            end do
        end do
        
        Xcoef(1,:) = coef
        Xbeta1 = matmul(Xcoef,transpose(ve))        

        Call ResidusMartingalea(b,np,funcpaares,Resmartingale,frailtypred,frailtyvar,frailtysd,&
        frailtypred2,frailtyvar2,frailtysd2,frailtycov)

        do i=1,nsujet
            linearpred(i)=Xbeta1(1,i) + frailtypred(g(i)) + frailtypred2(g(i)) * ve2(i,1)
        end do
    end if

!    write(*,*)'<======== Fin ========> 1'
    if(typeof==0) then
        deallocate(mm3,mm2,mm1,mm,im3,im2,im1,im,dut1,dut2,ut1,ut2)
        deallocate(zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
    else
        deallocate(vvv) !,kkapa)
        if (typeof == 1) then
            deallocate(ttt,betacoef)
        end if
    end if

    deallocate(cumulhaz1,som_Xbeta)

    deallocate(XBeta,invsigma,b_temp)

    deallocate(t0,t1,c,nt0,nt1,stra,g,stracross,aux,invd,nig,mid,date,ve,ve2,betaaux,vax,cumulhaz)    
    deallocate(I_hess,H_hess,Hspl_hess,hess,PEN_deri,y,v,I1_hess,H1_hess,I2_hess,H2_hess, &
    HI1,HI2,HIH,IH,HI)
    
    end subroutine additive

    
    
    
    !========================== VECSPLI =====================
    subroutine vecspliadd(n,ndate) 
    
    use tailles
    use comon,only:zi,date,mm3,mm2,mm1,mm,im3,im2,im1,im
    
    implicit none
    
    integer::n,ndate,i,j,k
    double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2
    double precision::h3,h4,h3m,h2n,hn,hh3,hh2
    
    
!----------  calcul de u(ti) ---------------------------
    j=0
    do i=1,ndate-1
        do k = 2,n-2
            if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
                j = k-1
            endif
        end do 
        ht = date(i)-zi(j)
        htm= date(i)-zi(j-1)
        h2t= date(i)-zi(j+2)
        ht2 = zi(j+1)-date(i)
        ht3 = zi(j+3)-date(i)
        hht = date(i)-zi(j-2)
        h = zi(j+1)-zi(j)
        hh= zi(j+1)-zi(j-1)
        h2= zi(j+2)-zi(j)
        h3= zi(j+3)-zi(j)
        h4= zi(j+4)-zi(j)
        h3m= zi(j+3)-zi(j-1)
        h2n=zi(j+2)-zi(j-1)
        hn= zi(j+1)-zi(j-2)
        hh3 = zi(j+1)-zi(j-3)
        hh2 = zi(j+2)-zi(j-2)
        mm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
        *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
        h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        im3(i) = (0.25d0*(date(i)-zi(j-3))*mm3(i))+(0.25d0*hh2 &
        *mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
        im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
        im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
        im(i)  = ht*mm(i)*0.25d0
    end do
    
    end subroutine vecspliadd
    
!========================== VECPEN ==============================
    subroutine vecpenadd(n) 
    
    use tailles
    use comon,only:zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m!date
    
    implicit none
    
    integer::n,i
    double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2, &
    a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x
    
    do i=1,n-3
        h = zi(i+1)-zi(i)
        hh= zi(i+1)-zi(i-1)
        h2= zi(i+2)-zi(i)
        h3= zi(i+3)-zi(i)
        h4= zi(i+4)-zi(i)
        h3m= zi(i+3)-zi(i-1)
        h2n=zi(i+2)-zi(i-1)
        hn= zi(i+1)-zi(i-2)
        hh3 = zi(i+1)-zi(i-3)
        hh2 = zi(i+2)-zi(i-2)
        a3 = h*hh*hn*hh3
        a2 = hh2*hh*h*hn
        b2 = hh2*h2n*hh*h
        c2 = hh2*h2*h*h2n
        a1 = h3m*h2n*hh*h
        b1 = h3m*h2*h*h2n
        c1 = h3m*h3*h2*h
        a0 = h4*h3*h2*h
        x3 = zi(i+1)*zi(i+1)*zi(i+1)-zi(i)*zi(i)*zi(i)
        x2 = zi(i+1)*zi(i+1)-zi(i)*zi(i)
        x  = zi(i+1)-zi(i)
        m3m3(i) = (192.d0*h/(hh*hn*hh3*hh*hn*hh3))
        m2m2(i) = 64.d0*(((3.d0*x3-(3.d0*x2*(2.d0*zi(i+1)+zi(i-2) &
        ))+x*(4.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i-2)+4.d0*zi(i+1) &
        *zi(i-2)))/(a2*a2)))
        m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2) &
        +zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1) &
        +zi(i+1)*zi(i+1)+2.d0*zi(i+2)*zi(i-1)+2.d0*zi(i+2) &
        *zi(i+1)+2.d0*zi(i-1)*zi(i+1)))/(b2*b2)))
        m2m2(i) = m2m2(i) +64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i+2) &
        +zi(i)))+x*(4.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i)+4.d0*zi(i+2) &
        *zi(i)))/(c2*c2))
        m2m2(i) = m2m2(i) +128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+2) &
        +zi(i-1)+3.d0*zi(i+1)+zi(i-2)))+x*(2.d0*zi(i+1)*zi(i+2) &
        +2.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i+2) &
        +zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*b2))
        m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0* &
        x2*(2.d0*zi(i+2)+zi(i)+2.d0*zi(i+1)+zi(i-2)))+x* &
        (4.d0*zi(i+1)*zi(i+2)+2.d0*zi(i+1)*zi(i)+2.d0*zi(i-2) &
        *zi(i+2)+zi(i-2)*zi(i)))/(a2*c2))
        m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0*x2 &
        *(3.d0*zi(i+2)+zi(i)+zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i)+ &
        2.d0*zi(i-1)*zi(i+2)+zi(i)*zi(i-1)+2.d0*zi(i+1)*zi(i+2) &
        +zi(i+1)*zi(i)+2.d0*zi(i+2)*zi(i+2)))/(b2*c2))
        m1m1(i) = 64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i-1)+zi(i+1))) &
        +x*(4.d0*zi(i-1)*zi(i-1)+zi(i+1)*zi(i+1)+4.d0*zi(i-1) &
        *zi(i+1)))/(a1*a1))
        m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i)   &  
        +zi(i+2)))+x*(zi(i-1)*zi(i-1)+zi(i)*zi(i)+zi(i+2)* &
        zi(i+2)+2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+2.d0* &
        zi(i)*zi(i+2)))/(b1*b1))
        m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i+3) &
        +2.d0*zi(i)))+x*(zi(i+3)*zi(i+3)+4.d0*zi(i)*zi(i) &
        +4.d0*zi(i+3)*zi(i)))/(c1*c1)) 
        m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(3.d0 &
        *zi(i-1)+zi(i)+zi(i+2)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i-1) &
        +2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1) &
        +zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(a1*b1))
        m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+ &
        2.d0*zi(i)+2.d0*zi(i-1)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i+3) &
        +4.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i))) &
        /(a1*c1))    
        m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+3.d0 &
        *zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1)  &   
        *zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3) &
        +2.d0*zi(i+2)*zi(i)))/(b1*c1))
        mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
        m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2) &
        ))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2)) &
        +((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x* &
        (zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2)) &
        +((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x* &
        (2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)))
        m3m1(i) = 192.d0*(((x3-(0.5d0*x2*(4.d0*zi(i+1)+2.d0*zi(i-1) &
        ))+x*(2.d0*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*a1)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+2)+zi(i-1)+zi(i))) &
        +x*(zi(i+1)*zi(i-1)+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(b1*a3)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+3)+2.d0*zi(i)))+x*(zi(i+1) &
        *zi(i+3)+2.d0*zi(i+1)*zi(i)))/(c1*a3)) )
        m3m(i) = 576.d0*((-(x3/3.d0)+(0.5d0*x2*(zi(i+1)+zi(i))) &
        -x*zi(i+1)*zi(i))/(a3*a0))
        m2m1(i) = 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1)+3.d0* &
        zi(i+1)+zi(i-2)))-x*(4.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1) &
        *zi(i+1)+2.d0*zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*a1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1)+ &
        zi(i)+zi(i+2)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i-1) &
        +2.d0*zi(i+1)*zi(i)+2.d0*zi(i+1)*zi(i+2)+zi(i-2)*zi(i-1)+ &
        zi(i-2)*zi(i)+zi(i-2)*zi(i+2)))/(a2*b1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)+2.d0 &
        *zi(i)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i+3)+4.d0 &
        *zi(i+1)*zi(i)+zi(i-2)*zi(i+3)+2.d0*zi(i-2)*zi(i)))/(a2*c1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2* &
        (3.d0*zi(i-1)+2.d0*zi(i+1)+zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1) &
        +zi(i+2)*zi(i+1)+2.d0*zi(i-1)*zi(i-1)+3.d0 &
        *zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(b2*a1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0 &
        *zi(i-1)+zi(i)+2.d0*zi(i+2)+zi(i+1)))-x*(zi(i+2)*zi(i-1) &
        +zi(i+2)*zi(i)+zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1)+zi(i-1) &
        *zi(i)+zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i) &
        +zi(i+1)*zi(i+2)))/(b2*b1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
        +2.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1)))-x*(zi(i+2)*zi(i+3) &
        +2.d0*zi(i+2)*zi(i)+zi(i-1)*zi(i+3)+2.d0*zi(i-1)*zi(i) &
        +zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i)))/(b2*c1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1) &
        +zi(i+1)+2.d0*zi(i+2)+zi(i)))-x*(4.d0*zi(i+2)*zi(i-1)+2.d0* &
        zi(i+2)*zi(i+1)+2.d0*zi(i)*zi(i-1)+zi(i)*zi(i+1)))/(c2*a1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1) &
        +2.d0*zi(i)+3.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1)+2.d0 &
        *zi(i+2)*zi(i)+2.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i-1)+zi(i) &
        *zi(i)+zi(i)*zi(i+2)))/(c2*b1))
        m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
        +3.d0*zi(i)+2.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i+3)+4.d0 &
        *zi(i+2)*zi(i)+zi(i)*zi(i+3)+2.d0*zi(i)*zi(i)))/(c2*c1))
        m2m(i) = 192.d0*(((x3-(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i+1) &
        +zi(i-2)))+x*(2.d0*zi(i+1)*zi(i)+zi(i-2)*zi(i)))/(a2*a0)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1))) &
        +x*(zi(i+2)*zi(i)+zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(b2*a0)) &
        +((x3-(0.5d0*x2*(4.d0*zi(i)+2.d0*zi(i+2)))+x*(2.d0*zi(i+2) &
        *zi(i)+zi(i)*zi(i)))/(c2*a0)) )
        m1m(i) = 192.d0*(((-x3+(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i-1) &
        +zi(i+1)))-x*(2.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(a1*a0)) &
        +((-x3+(0.5d0*x2*(4.d0*zi(i)+zi(i-1)+zi(i+2))) &
        -x*(zi(i-1)*zi(i)+zi(i)*zi(i)+zi(i+2)*zi(i)))/(b1*a0)) &
        +((-x3+(0.5d0*x2*(5.d0*zi(i)+zi(i+3)))-x*(zi(i+3)*zi(i) &
        +2.d0*zi(i)*zi(i)))/(c1*a0)) )
    end do
    
    end subroutine vecpenadd

!==========================  SUSP  ====================================
    subroutine suspadd(x,the,n,su,lam,zi)
    
    use tailles
    
    implicit none
    integer::j,k,n,i
    double precision::x,ht,ht2,h2,som,lam,su,htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3
    double precision::ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,gl,hh 
    double precision,dimension(-2:npmax)::zi,the
    
    som = 0.d0
    gl=0.d0
    do k = 2,n+1
        if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            j = k-1
            if (j.gt.1)then
                do i=2,j
                    som = som+the(i-4)
                end do  
            endif   
            ht = x-zi(j)
            htm= x-zi(j-1)
            h2t= x-zi(j+2)
            ht2 = zi(j+1)-x
            ht3 = zi(j+3)-x
            hht = x-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        endif
    end do
    
    if(x.ge.zi(n))then
        som = 0.d0
        do i=1,n+1
            som = som+the(i-3)
        end do
        gl = som
    endif
    su  = dexp(-gl)
    
    return
    
    end subroutine suspadd
    
!==========================  COSP  ====================================
    
    subroutine cospadd(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)
    
    use tailles
    
    implicit none
    
    integer::j,k,n,i
    double precision::x,ht,ht2,h2,som,lam,su,binf,bsup,lbinf,lbsup,pm, &
    htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3,ht3,hht,h4, &
    h3m,hh3,hh2,mm,im3,mm2,h,gl,hh
    double precision,dimension(-2:npmax)::the,zi
    double precision,dimension(npmax,npmax)::y

    j=0
    gl=0.d0    
    som = 0.d0
    do k = 2,n-1
        if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            j = k-1
            if (j.gt.1)then
                do i=2,j
                    som = som+the(i-4)
                end do  
            endif   
            ht = x-zi(j)
            htm= x-zi(j-1)
            h2t= x-zi(j+2)
            ht2 = zi(j+1)-x
            ht3 = zi(j+3)-x
            hht = x-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        endif
    end do
    
    if(x.ge.zi(n))then
        som = 0.d0
        do i=1,n
            som = som+the(i-3)
        end do
        gl = som
    endif
    
    call confadd(x,j,n,y,pm,zi)
    
    binf = dexp(-gl - 1.96d0*pm)
    su  = dexp(-gl)
    bsup = dexp(-gl + 1.96d0*pm)
    
    call conf1add(x,j,n,y,pm,zi)
    lbinf = lam - 1.96d0*pm
    lbsup = lam + 1.96d0*pm
    
    return
    
    end subroutine cospadd
    
!=====================  CONF1  =============================
    subroutine conf1add(x,ni,n,y,pm,zi)
    
    use tailles
    
    implicit none
    
    integer::ni,i,n,j
    double precision::mmspadd,x,pm,res
    double precision,dimension(-2:npmax) :: zi
    double precision,dimension(npmax) :: vecti,aux
    double precision,dimension(npmax,npmax) :: y
    
    do i=1,n
        vecti(i) = mmspadd(x,ni,i,zi)
    end do   
    
    do i=1,n
        aux(i) = 0.d0
        do j=1,n
            aux(i) = aux(i) - y(i,j)*vecti(j)
        end do
    end do   
    
    res = 0.d0
    do i=1,n
        res = res + aux(i)*vecti(i)
    end do
    
    if (res.lt.0)then 
        res=-res
    endif 
    
    pm = dsqrt(res) 
    
    end subroutine conf1add
    

!=====================  CONF  =============================
    subroutine confadd(x,ni,n,y,pm,zi)
    
    use tailles
    
    implicit none
    
    integer::ni,i,n,j
    double precision::ispadd,x,pm,res
    double precision,dimension(-2:npmax) :: zi
    double precision,dimension(52) :: vecti,aux
    double precision,dimension(npmax,npmax) :: y
    
    
    do i=1,n
        vecti(i) = ispadd(x,ni,i,zi)
    end do   
    
    do i=1,n
        aux(i) = 0.d0
        do j=1,n
            aux(i) = aux(i) - y(i,j)*vecti(j)
        end do
    end do   
    
    res = 0.d0
    do i=1,n
        res = res + aux(i)*vecti(i)
    end do
    
    if (res.lt.0)then 
        res=-res
    endif 
    
    pm = dsqrt(res)
    
    end subroutine confadd
    

!==========================   ISP   ==================================
    double precision function ispadd(x,ni,ns,zi)
    
    use tailles
    
    implicit none
    
    integer::ni,ns
    double precision::val,mmspadd,x
    double precision,dimension(-2:npmax)::zi
    
    if(x.eq.zi(ni))then
        if(ni.le.ns-3)then
            val = 0.d0
        else
            if(ni.le.ns-2)then
                val = ((zi(ni)-zi(ni-1))*mmspadd(x,ni,ns,zi))*0.25d0
            else
                if (ni.eq.ns-1)then
                    val = ((zi(ni)-zi(ni-2))*mmspadd(x,ni,ns,zi)+ &
                    (zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+1,zi))*0.25d0
                else
                    if(ni.eq.ns)then
                        val = ((zi(ni)-zi(ni-3))*mmspadd(x,ni,ns,zi)+ &
                        (zi(ni+2)-zi(ni-2))*mmspadd(x,ni,ns+1,zi) &
                        +(zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+2,zi))*0.25d0
                    else
                        val = 1.d0
                    endif
                endif
            endif   
        endif
    else   
        if(ni.lt.ns-3)then
            val = 0.d0
        else
            if(ni.eq.ns-3)then
                val = (x-zi(ni))*mmspadd(x,ni,ns,zi)*0.25d0
            else  
                if(ni.eq.ns-2)then
                    val = ((x-zi(ni-1))*mmspadd(x,ni,ns,zi)+ &
                    (zi(ni+4)-zi(ni))*mmspadd(x,ni,ns+1,zi))*0.25d0
                else   
                    if (ni.eq.ns-1)then
                        val =((x-zi(ni-2))*mmspadd(x,ni,ns,zi)+ &
                        (zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+1,zi) &
                        +(zi(ni+4)-zi(ni))*mmspadd(x,ni,ns+2,zi))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val =((x-zi(ni-3))*mmspadd(x,ni,ns,zi)+ &
                            (zi(ni+2)-zi(ni-2))*mmspadd(x,ni,ns+1,zi) &
                            +(zi(ni+3)-zi(ni-1))*mmspadd(x,ni,ns+2,zi) &
                            +(zi(ni+4)-zi(ni))*mmspadd(x,ni,ns+3,zi))*0.25d0
                        else
                            val = 1.d0
                        endif
                    endif
                endif
            endif
        endif 
    endif
    
    ispadd = val
    
    return
    
    end function ispadd
    
!==========================  MMSP   ==================================
    
    double precision function mmspadd(x,ni,ns,zi)
    
    use tailles
    
    implicit none
    
    integer::ni,ns
    double precision::val,x
    double precision,dimension(-2:npmax)::zi

    
    if(ni.lt.ns-3)then
        val = 0.d0
    else
        if(ns-3.eq.ni)then
            if(x.eq.zi(ni))then
                val = 0.d0
            else  
                val = (4.d0*(x-zi(ni))*(x-zi(ni)) &
                *(x-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3) &
                -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
            endif
        else 
            if(ns-2.eq.ni)then
                if(x.eq.zi(ni))then
                    val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1))) &
                    /((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1)) &
                    *(zi(ni+1)-zi(ni-1)))
                else  
                    val = (4.d0*(x-zi(ni-1))*(x-zi(ni-1)) &
                    *(zi(ni+1)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    -zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni))) &
                    +   (4.d0*(x-zi(ni-1))*(x-zi(ni)) &
                    *(zi(ni+2)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1)))  &
                    +   (4.d0*(x-zi(ni))*(x-zi(ni)) &
                    *(zi(ni+3)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+3) &
                    -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                endif
            else   
                if (ns-1.eq.ni)then
                    if(x.eq.zi(ni))then
                        val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1) &
                        -zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1) &
                        -zi(ni-1))*(zi(ni+1)-zi(ni-2)))) &
                        +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni))) &
                        /((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1)) &
                        *(zi(ni+1)-zi(ni-1)))))
                    else
                        val = (4.d0*((x-zi(ni-2))*(zi(ni+1) &
                        -x)*(zi(ni+1)-x))/((zi(ni+2) &
                        -zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)- &
                        zi(ni))*(zi(ni+1)-zi(ni-2)))) &
                        +((4.d0*((x-zi(ni-1))*(zi(ni+2)-x)  &
                        *(zi(ni+1)-x))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni))))) &
                        +((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x)  &
                        *(x-zi(ni)))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni)))))
                    endif 
                else
                    if(ni.eq.ns)then
                        if(x.eq.zi(ni))then
                            val =(4.d0*(x-zi(ni+1))*(x &
                            -zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1) &
                            -zi(ni-2))*(zi(ni+1)-zi(ni-3))))
                        else   
                            val =(4.d0*(x-zi(ni+1))*(x &
                            -zi(ni+1))*(zi(ni+1)-x)/((zi(ni+1) &
                            -zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1) &
                            -zi(ni))*(zi(ni+1)-zi(ni-3))))
                        endif
                    else
                        val = 0.d0
                    endif
                endif
            endif
        endif
    endif
    
    mmspadd = val
    
    return
    
    end function mmspadd



    




!=====================cross validation

!========================          MNBRAK         ===================
    subroutine mnbrakadd(ax,bx,cx,fa,fb,fc,b,n)
    use tailles
    implicit none
    
    double precision::ax,bx,cx,fa,fb,fc,aux,res,dum,fu,q,r,u,ulim
    double precision,dimension(npmax)::b
    double precision,dimension(npmax,npmax)::y
    double precision::estimvadd,gold,glimit,tiny
    parameter (gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
    integer::n,ni
    
    
    fa = estimvadd(ax,n,b,y,aux,ni,res)
    
    fb = estimvadd(bx,n,b,y,aux,ni,res)
    
    
    if(fb.gt.fa)then
        dum = ax
        ax = bx
        bx = dum
        dum = fb
        fb = fa
        fa = dum
    endif
    cx = bx + gold*(bx-ax)
    
    fc = estimvadd(cx,n,b,y,aux,ni,res)
    
1       if(fb.ge.fc)then
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx-((bx-cx)*q-(bx-ax)*r)/ &
        (2.d0*sign(max(abs(q-r),tiny),q-r))
        ulim = bx + glimit*(cx-bx)
        if((bx-u)*(u-cx).gt.0.d0)then
            fu = estimvadd(u,n,b,y,aux,ni,res)
            if(fu.lt.fc)then
                ax = bx
                fa = fb
                bx = u
                fb = fu
                return
            else
                if(fu.gt.fb)then
                    cx = u
                    fc = fu
                    return
                endif   
            endif
            u = cx + gold*(cx-bx)
            fu = estimvadd(u,n,b,y,aux,ni,res)
        else
            if((cx-u)*(u-ulim).gt.0.d0)then
                fu = estimvadd(u,n,b,y,aux,ni,res)
                if(fu.lt.fc)then
                    bx = cx
                    cx = u
                    u = cx + gold*(cx-bx)
                    fb = fc
                    fc = fu
                    fu = estimvadd(u,n,b,y,aux,ni,res)
                endif  
            else
                if((u-ulim)*(ulim-cx).ge.0.d0)then
                    u = ulim
                    fu = estimvadd(u,n,b,y,aux,ni,res)
                else
                    u = cx + gold*(cx-bx)
                    fu = estimvadd(u,n,b,y,aux,ni,res)
                endif
            endif   
        endif
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        goto 1
    endif
    
    return 
    
    end subroutine mnbrakadd

!========================      GOLDEN   =========================
    double precision function goldenadd(ax,bx,cx,tol,xmin,n,b,y,aux)
    
    use tailles
    
    implicit none
    
    double precision,dimension(npmax,npmax)::y
    double precision,dimension(npmax)::b
    double precision ax,bx,cx,tol,xmin,r,c,aux,res
    parameter (r=0.61803399d0,c=1.d0-r)
    double precision::f1,f2,x0,x1,x2,x3,estimvadd
    integer::n,ni
    
    x0 = ax
    x3 = cx
    if(abs(cx-bx).gt.abs(bx-ax))then
        x1 = bx
        x2 = bx + c*(cx-bx)
    else
        x2 = bx
        x1 = bx - c*(bx-ax)
    endif
    f1 = estimvadd(x1,n,b,y,aux,ni,res)
    f2 = estimvadd(x2,n,b,y,aux,ni,res)
    
1       if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
            if(f2.lt.f1)then
                x0 = x1
                x1 = x2
                x2 = r*x1 + c*x3
                f1 = f2
                f2 = estimvadd(x2,n,b,y,aux,ni,res)
    
            else
                x3 = x2
                x2 = x1
                x1 = r*x2+c*x0
                f2 = f1
                f1 = estimvadd(x1,n,b,y,aux,ni,res)
    
            endif
            go to 1
        endif
        if(f1.lt.f2)then
            goldenadd = f1
            xmin = x1
        else
            goldenadd = f2
            xmin = x2
        endif
        return
    
    end function goldenadd
    

!========================          ESTIMV         ===================

    double precision function estimvadd(k00,n,b,y,aux,ni,res)
    
    use tailles
    !use comon,only:c,nst,nsujet,nt0,nt1,nva,nz1,nz2,t0,t1,
    use comon,only:ndate,date,zi,pe,effet,mm3,mm2,mm1,mm,im3,im2,im1,im
    !use additiv,only:correl
    use optim
    
    implicit none
    integer n,ij,i,k,j,vj,ier,istop,ni
    double precision,dimension((n*(n+3)/2))::v
    double precision,dimension(n,n)::y
    double precision,dimension(ndatemax)::ut,dut
    double precision,dimension(n)::bh,b
    double precision,dimension(-2:npmax)::the    
    double precision::res,k00,som,h1,aux
    double precision,dimension(2)::k0    
    double precision::ca,cb,dd,funcpaasplines
    external::funcpaasplines
    
    ca=0.d0
    cb=0.d0
    dd=0.d0
    j=0
    estimvadd=0.d0
    k0(1) = k00*k00
    k0(2)=0.d0
    
    call marq98J(k0,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaasplines)    
!AD:    
    if (istop.eq.4) goto 50
!AD:        
    if(k0(1).gt.0.d0)then
        do ij=1,n
            the(ij-3)=(b(ij))*(b(ij))
            bh(ij) = (b(ij))*(b(ij))
        end do
    
        vj = 0
        som = 0.d0
        dut(1) = (the(-2)*4.d0/(zi(2)-zi(1)))
        ut(1) = the(-2)*dut(1)*0.25d0*(zi(1)-zi(-2))
        do i=2,ndate-1
            do k = 2,n-2
                if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
                    j = k-1
                    if ((j.gt.1).and.(j.gt.vj))then
                        som = som+the(j-4)
                        vj  = j
                    endif   
                endif
            end do 
            ut(i) = som +(the(j-3)*im3(i))+(the(j-2)*im2(i)) &
            +(the(j-1)*im1(i))+(the(j)*im(i))
            dut(i) = (the(j-3)*mm3(i))+(the(j-2)*mm2(i)) &
            +(the(j-1)*mm1(i))+(the(j)*mm(i))
        end do
        i = n-2
        h1 = (zi(i)-zi(i-1))
        ut(ndate) = som+ the(i-4) + the(i-3)+the(i-2)+the(i-1)
        dut(ndate) = (4.d0*the(i-1)/h1)
    
        call testadd(dut,k0,n,aux,y)
        estimvadd = - ((res-pe)) - aux
    
    else
        aux = -n
    endif
!AD:
50    continue      
!AD:
    return
    end function estimvadd
      
!=================calcul de la hessienne  et de omega  ==============
    subroutine testadd(dut,k0,n,res,y)
    
    use tailles
    !use comon,only:c,nst,nsujet,date,nt0,nt1,nva,t0,t1,ndate,zi
    implicit none
    
    integer::n,i,j,np
    double precision::res,tra,d
    double precision,dimension(npmax,npmax)::hessh,hess,omeg,y
    double precision,dimension(ndatemax)::dut
    integer,dimension(npmax)::indx     
    double precision,dimension(2)::k0 
    
    do i = 1,n
        do j = 1,n
        hess(i,j) = 0.d0 
        end do
    end do    
    
    do i = 1,n
        do j = i,n
            call matadd(hess(i,j),dut,i,j,n)
        end do
    end do
    do i = 2,n
        do j = 1,i-1
            hess(i,j)=hess(j,i)
        end do
    end do    
    
    call calcomegadd(n,omeg)
    
    do i = 1,n
        do j = 1,n
            hessh(i,j)=-hess(i,j)
            hess(i,j) = hess(i,j) - (2.d0*k0(1)*omeg(i,j)) 
        end do   
    end do
    
    np = n
    do i=1,n
        do j=1,n
            y(i,j)=0.d0
        end do
        y(i,i)=1.d0
    end do
    
    call ludcmpadd(hess,n,indx,d)
    
    do j=1,n
        call lubksbadd(hess,n,indx,y(1,j))
    end do
    
    tra = 0.d0
    do i=1,n
        do j=1,n
            tra = tra + y(i,j)*hessh(j,i)
        end do
    end do
    
    res = (tra)
    
    end subroutine testadd

!======================  LUBKSB  ======================================
    subroutine lubksbadd(a,n,indx,b)
    
    use tailles
    
    implicit none
    
    integer::n,i,ii,j,ll
    double precision::sum
    integer,dimension(npmax)::indx   
    double precision,dimension(npmax)::b
    double precision,dimension(npmax,npmax)::a
    
    ii = 0
    do i=1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if(ii.ne.0)then
            do j=ii,i-1
                sum = sum -a(i,j)*b(j)
            end do
        else
            if(sum.ne.0.d0)then
                ii=i
            endif
        endif
        b(i)=sum
    end do
    do i=n,1,-1
        sum = b(i)
        do j = i+1,n
            sum = sum-a(i,j)*b(j)
        end do
        b(i)=sum/a(i,i)
    end do
    
    return
    
    end subroutine lubksbadd    

!======================  LUDCMP  ======================================
    subroutine ludcmpadd(a,n,indx,d)
    
    use tailles
    
    implicit none
    
    integer::n,nmax,i,imax,j,k
    integer,dimension(npmax)::indx 
    double precision,dimension(npmax,npmax)::a         
    double precision::d,tiny,aamax,dum,sum
    parameter (nmax=500,tiny=1.d-20)
    double precision,dimension(nmax)::vv
    
    d = 1.d0
    imax=0
    do i=1,n
        aamax=0.d0
        do j=1,n
            if (dabs(a(i,j)).gt.aamax)then
                aamax=dabs(a(i,j))
            endif
        end do
        if (aamax.eq.0.d0) then
        end if
        vv(i) = 1.d0/aamax
    end do
    do j = 1,n
        do i=1,j-1
            sum = a(i,j)
            do k=1,i-1
                sum = sum - a(i,k)*a(k,j)
            end do
            a(i,j) = sum
        end do
        aamax = 0.d0
        do i = j,n
            sum = a(i,j)
            do k=1,j-1
                sum = sum -a(i,k)*a(k,j)
            end do
            a(i,j) = sum
            dum = vv(i)*dabs(sum)
            if (dum.ge.aamax) then
                imax = i
                aamax = dum
            endif
        end do
        if(j.ne.imax)then
            do k=1,n
                dum = a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k) = dum
            end do
            d = -d
            vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)then
            a(j,j)=tiny
        endif
        if(j.ne.n)then
            dum = 1.d0/a(j,j)
            do i = j+1,n
                a(i,j) = a(i,j)*dum
            end do
        endif
    end do
    
    return
        
    end subroutine ludcmpadd
!=======================  CALOMEG  ===========================
    subroutine calcomegadd(n,omeg)
    
    use tailles
    !use comon,only:date,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m2m1,m2m,m1m
    use comon,only:m3m
    
    implicit none
    
    integer::n,i,j
    double precision::calc00add,calc01add,calc02add
    double precision,dimension(npmax,npmax)::omeg
    
    
    do i=1,n
        do j=1,n
        omeg(i,j)=0.d0
        end do
    end do
    
    omeg(1,1)=calc00add(1,n)
    omeg(1,2)=calc01add(1,n)
    omeg(1,3)=calc02add(1,n)
    omeg(1,4)=m3m(1)
    omeg(2,1)=omeg(1,2)
    omeg(2,2)=calc00add(2,n)
    omeg(2,3)=calc01add(2,n)
    omeg(2,4)=calc02add(2,n)
    omeg(2,5)=m3m(2)
    omeg(3,1)=omeg(1,3)
    omeg(3,2)=omeg(2,3)
    omeg(3,3)=calc00add(3,n)
    omeg(3,4)=calc01add(3,n)
    omeg(3,5)=calc02add(3,n)
    omeg(3,6)=m3m(3)
    do i=4,n-3
        omeg(i,i-3)=omeg(i-3,i)
        omeg(i,i-2)=omeg(i-2,i)
        omeg(i,i-1)=omeg(i-1,i)
        omeg(i,i)=calc00add(i,n)
        omeg(i,i+1)=calc01add(i,n)
        omeg(i,i+2)=calc02add(i,n)
        omeg(i,i+3)=m3m(i)
    end do   
    omeg(n-2,n-5)=omeg(n-5,n-2)
    omeg(n-2,n-4)=omeg(n-4,n-2)
    omeg(n-2,n-3)=omeg(n-3,n-2)
    omeg(n-2,n-2)=calc00add(n-2,n)
    omeg(n-2,n-1)=calc01add(n-2,n)
    omeg(n-2,n)=calc02add(n-2,n)
    omeg(n-1,n-4)=omeg(n-4,n-1)
    omeg(n-1,n-3)=omeg(n-3,n-1)
    omeg(n-1,n-2)=omeg(n-2,n-1)
    omeg(n-1,n-1)=calc00add(n-1,n)
    omeg(n-1,n)=calc01add(n-1,n)
    omeg(n,n-3)=omeg(n-3,n)
    omeg(n,n-2)=omeg(n-2,n)
    omeg(n,n-1)=omeg(n-1,n)
    omeg(n,n)=calc00add(n,n)
    
    end subroutine calcomegadd


!====================  MAT  ==================================
    subroutine matadd(res,dut,k,l,n)
    
    use tailles
    !use comon,only:nst,t0,t1,nt0,nva,ndate
    use comon,only:date,zi,c,nt1,nsujet

    implicit none
    
    integer::k,l,j,ni,n,i
    double precision,dimension(ndatemax)::dut
    double precision::res,res1,mspadd,aux2,u2
    
!---------- calcul de la hessienne ij ------------------
    res = 0.d0
    res1 = 0.d0
    do i=1,nsujet
        if(c(i).eq.1)then  !event
            u2 = dut(nt1(i)) 
            do j = 2,n-2
            if((date(nt1(i)).ge.zi(j-1)).and.(date(nt1(i)).lt.zi(j)))then
                ni = j-1
            endif
        end do 
            if(date(nt1(i)).eq.zi(n-2))then
                ni = n-2
            endif   
!-------attention numero spline 
            aux2 = mspadd(nt1(i),ni,k)*mspadd(nt1(i),ni,l)
            if (u2.le.0.d0)then
                res1 = 0.d0
            else   
                res1 = - aux2/(u2*u2)
            endif  
        else !censure  
            res1 = 0.d0
        endif 
        res = res + res1
    end do   
    
    end subroutine matadd

!==========================  MSP   ==================================
    double precision function mspadd(i,ni,ns)
    
    use tailles
    use comon,only:date,zi
    
    implicit none
    
    integer::ni,ns,i
    double precision::val
    
    
    if(ni.lt.ns-3)then
        val = 0.d0
    else
        if(ns-3.eq.ni)then
            if(date(i).eq.zi(ni))then
                val = 0.d0
            else  
                val = (4.d0*(date(i)-zi(ni))*(date(i)-zi(ni)) &
                *(date(i)-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3) &
                -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
            endif
        else 
            if(ns-2.eq.ni)then
                if(date(i).eq.zi(ni))then
                    val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1))) &
                    /((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1)) &
                    *(zi(ni+1)-zi(ni-1)))
                else  
                    val = (4.d0*(date(i)-zi(ni-1))*(date(i)-zi(ni-1)) &
                    *(zi(ni+1)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    -zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni))) &
                    +   (4.d0*(date(i)-zi(ni-1))*(date(i)-zi(ni)) &
                    *(zi(ni+2)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1)))  &
                    +   (4.d0*(date(i)-zi(ni))*(date(i)-zi(ni)) &
                    *(zi(ni+3)-date(i)))/((zi(ni+3)-zi(ni-1))*(zi(ni+3) &
                    -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                endif
            else   
                if (ns-1.eq.ni)then
                    if(date(i).eq.zi(ni))then
                        val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1) &
                        -zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1) &
                        -zi(ni-1))*(zi(ni+1)-zi(ni-2)))) &
                        +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni)))  &
                        /((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1)) &
                        *(zi(ni+1)-zi(ni-1)))))
                    else
                        val = (4.d0*((date(i)-zi(ni-2))*(zi(ni+1) &
                        -date(i))*(zi(ni+1)-date(i)))/((zi(ni+2) &
                        -zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)- &
                        zi(ni))*(zi(ni+1)-zi(ni-2)))) &
                        +((4.d0*((date(i)-zi(ni-1))*(zi(ni+2)-date(i))  &
                        *(zi(ni+1)-date(i)))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni))))) &
                        +((4.d0*((zi(ni+2)-date(i))*(zi(ni+2)-date(i))  &
                        *(date(i)-zi(ni)))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni)))))
                    endif 
                else
                    if(ni.eq.ns)then
                        if(date(i).eq.zi(ni))then
                            val =(4.d0*(date(i)-zi(ni+1))*(date(i) &
                            -zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1) &
                            -zi(ni-2))*(zi(ni+1)-zi(ni-3))))
                        else   
                            val =(4.d0*(date(i)-zi(ni+1))*(date(i) &
                            -zi(ni+1))*(zi(ni+1)-date(i))/((zi(ni+1) &
                            -zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1) &
                            -zi(ni))*(zi(ni+1)-zi(ni-3))))
                        endif
                    else
                        val = 0.d0
                    endif
                endif
            endif
        endif
    endif
    
    mspadd = val
        
    return
        
    end function mspadd

!==========================   SP   ==================================
    double precision function spadd(i,ni,ns)
    
    use tailles
    use comon,only:date,zi
    
    implicit none
    
    integer::ni,ns,i
    double precision::val,mspadd
    
    if(date(i).eq.zi(ni))then
        if(ni.le.ns-3)then
            val = 0.d0
        else
            if(ni.le.ns-2)then
                val = ((zi(ni)-zi(ni-1))*mspadd(i,ni,ns))*0.25d0
            else
                if (ni.eq.ns-1)then
                    val = ((zi(ni)-zi(ni-2))*mspadd(i,ni,ns)+ &
                    (zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+1))*0.25d0
                else
                    if(ni.eq.ns)then
                        val = ((zi(ni)-zi(ni-3))*mspadd(i,ni,ns)+ &
                        (zi(ni+2)-zi(ni-2))*mspadd(i,ni,ns+1) &
                        +(zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+2))*0.25d0
                    else
                        val = 1.d0
                    endif
                endif
            endif   
        endif
    else   
        if(ni.lt.ns-3)then
            val = 0.d0
        else
            if(ni.eq.ns-3)then
                val = (date(i)-zi(ni))*mspadd(i,ni,ns)*0.25d0
            else  
                if(ni.eq.ns-2)then
                    val = ((date(i)-zi(ni-1))*mspadd(i,ni,ns)+ &
                    (zi(ni+4)-zi(ni))*mspadd(i,ni,ns+1))*0.25d0
                else   
                    if (ni.eq.ns-1)then
                        val =((date(i)-zi(ni-2))*mspadd(i,ni,ns)+ &
                        (zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+1) &
                        +(zi(ni+4)-zi(ni))*mspadd(i,ni,ns+2))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val =((date(i)-zi(ni-3))*mspadd(i,ni,ns)+ &
                            (zi(ni+2)-zi(ni-2))*mspadd(i,ni,ns+1) &
                            +(zi(ni+3)-zi(ni-1))*mspadd(i,ni,ns+2) &
                            +(zi(ni+4)-zi(ni))*mspadd(i,ni,ns+3))*0.25d0
                        else
                            val = 1.d0
                        endif
                    endif
                endif
            endif
        endif 
    endif
    spadd = val
        
    return
        
    end function spadd
!================

!=========================  CALC00  =========================
    double precision function calc00add(j,n) 
    
    use tailles
    !use comon,only:m3m2,m3m1,m3m,m2m1,m2m,m1m
    use comon,only:m3m3,m2m2,m1m1,mmm
        
    implicit none
    
    double precision::part
    integer::j,n
    
    
    if(j.eq.1)then
        part = m3m3(j)
    else
        if(j.eq.2)then
            part = m3m3(j) + m2m2(j-1)
        else
            if(j.eq.3)then
                part = m3m3(j) + m2m2(j-1) + m1m1(j-2)
            else
                if(j.eq.n-2)then
                    part = m2m2(j-1) + m1m1(j-2) + mmm(j-3)
                else   
                    if(j.eq.n-1)then
                        part = mmm(j-3) + m1m1(j-2)
                    else
                        if(j.eq.n)then
                            part = mmm(j-3)
                        else   
                            part=mmm(j-3)+m1m1(j-2)+m2m2(j-1)+m3m3(j)
                        endif
                    endif
                endif   
            endif   
        endif   
    endif 
    
    calc00add = part
    
    return
    
    end function calc00add
    
!=========================  CALC01  =========================
    
    double precision function calc01add(j,n)
    
    use tailles
    !use comon,only:m3m3,m2m2,m1m1,mmm,m3m1,m3m,m2m
    use comon,only:m3m2,m2m1,m1m
        
    implicit none
    
    double precision::part
    integer::j,n
    
    
    
    if(j.eq.1)then
        part = m3m2(j)
    else   
        if(j.eq.2)then
            part = m3m2(j) + m2m1(j-1) 
        else
            if(j.eq.n-2)then
                part = m1m(j-2) + m2m1(j-1) 
            else
                if(j.ne.n-1)then
                    part = m3m2(j) + m2m1(j-1) + m1m(j-2)
                else
                    part = m1m(j-2)
                endif
            endif   
        endif
    endif   
    
    calc01add = part
    
    return
    
    end function calc01add

!=========================  CALC02  =========================
    double precision function calc02add(j,n)
    
    use tailles
    !use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m,m2m1,m1m
    use comon,only:m3m1,m2m
        
    implicit none
    
    
    
    double precision::part
    integer::j,n
    
    
    
    if(j.eq.1)then
        part = m3m1(j)
    else   
        if(j.ne.n-2)then
            part = m3m1(j) + m2m(j-1) 
        else
            part = m2m(j-1)
        endif
    endif   
    
    calc02add = part
    
    return
    
    end function calc02add


!===============================    MARQ98AUX =========================