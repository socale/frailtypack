
    subroutine frailpenal(nsujetAux,ngAux,icenAux,nstAux,effetAux, &
    nzAux,axT,tt0Aux,tt1Aux,icAux,groupeAux,nvaAux,strAux,vaxAux, &
    AGAux,noVar,maxitAux,irep1,np,b,H_hessOut,HIHOut,resOut,LCV, &
    xTOut,lamTOut,xSuT,suTOut,typeof0,equidistant,nbintervR0,mt, &
    ni,cpt,ier,k0,ddl,istop,shapeweib,scaleweib,mt1,ziOut,Resmartingale,martingaleCox,&
    frailtypred,frailtyvar,frailtysd,linearpred,time,intcensAux,ttUAux,logNormal0, &
    timedep0,nbinnerknots0,qorder0,filtretps0,BetaTpsMat,EPS) ! rajout


!
! Obs: noVar=0 indicates no variables but I need to pass all 0's in vaxAux
    use parameters
    use tailles
    use comon
    use optim
    use residusM
    use betatttps

    implicit none

    integer::groupe,ij,kk,j,k,nz,n,np,cpt,ii,iii,ver,cptstr1,cptstr2, &
    i,ic,ni,ier,istop,cptni,cptni1,cptni2,ibou,nbou2,id,cptbiais,l, &
    icen,irep1,nvacross,nstcross,effetcross,mt,mt1,p

    integer::ss,sss,noVar,AGAux,maxitAux,nsujetAux,ngAux,icenAux,nstAux,effetAux, &
    nzAux,nvaAux,intcensAux,str

    integer,dimension(nvaAux)::filtre
    integer,dimension(nsujetAux)::stracross
    double precision,dimension(nvaAux)::vax
    double precision::tt0,tt1,ttU
    double precision::h,res,min,max,maxt,pord !,scoreT,scoretest
    double precision,dimension(2)::auxkappa,k0,res01
    double precision,dimension(3*nsujetAux)::aux
    double precision,dimension(np*(np+3)/2)::v
    double precision,dimension(np)::b
    double precision,dimension(np,np)::HIH,IH,HI
    !******************************************   Add JRG January 05
    
    double precision,dimension(nzAux+6),intent(out)::ziOut
    double precision::resOut
    double precision,dimension(nsujetAux)::tt0Aux,tt1Aux,ttUAux
    integer,dimension(nsujetAux)::icAux,groupeAux
    double precision,dimension(nsujetAux)::strAux
    double precision,dimension(nsujetAux,nvaAux)::vaxAux
    double precision,dimension(np,np)::H_hessOut,HIHOut
    double precision,dimension(nstAux),intent(out)::shapeweib ,scaleweib
    double precision,dimension(nstAux),intent(in):: axT
    double precision,dimension(nstAux):: xminT
    double precision,dimension(mt,nstAux)::xTOut
    double precision,dimension(mt,3,nstAux)::lamTOut
    double precision,dimension(mt1,3,nstAux)::suTOut
    double precision,dimension(100,nstAux)::xSuT
    !*******************************************  Add JRG May 05 (Cross-validation)    
    double precision::auxi,ax,bx,cx,tol,ddl,fa,fb,fc,goldens,estimvs, &
    xmin1,xmin2
    double precision,dimension(np,np)::y
!AD:add
    double precision,dimension(2),intent(out)::LCV
    double precision::ca,cb,dd
    double precision,external::funcpassplines,funcpascpm,funcpasweib
    double precision,external::funcpascpm_intcens,funcpassplines_intcens,funcpasweib_intcens
    double precision,external::funcpassplines_log,funcpasweib_log,funcpascpm_log
    double precision,external::funcpas_tps
!AD:
!Cpm
    integer::typeof0,nbintervR0,equidistant,ent,indd
    double precision::temp
!predictor
    double precision,dimension(ngAux),intent(out)::Resmartingale,frailtypred,frailtysd,frailtyvar
    double precision,external::funcpasres
    double precision,dimension(nsujetAux),intent(out)::linearpred,martingaleCox
    double precision,dimension(1,nvaAux)::coefBeta
    double precision::coefBeta2
    double precision,dimension(1,nsujetAux)::XBeta
    double precision,dimension(nbintervR0+1)::time
    integer,dimension(2)::istopp
    integer,intent(in)::logNormal0
    integer,intent(in)::timedep0,nbinnerknots0,qorder0
    integer,dimension(nvaAux),intent(in)::filtretps0
    double precision,dimension(0:100,0:4*sum(filtretps0))::BetaTpsMat ! matrice des effets depedants du temps (4 colonnes par effet dep du tps)
    double precision,dimension(nbinnerknots0+qorder0)::basis
    double precision,dimension(3),intent(inout)::EPS ! seuils de convergence

!verification
    !print*,nsujetAux,ngAux,icenAux,nstAux,effetAux,nzAux,ax1,ax2,nvaAux
    !do i =1,5
    !    print*,i,tt0Aux(i),tt1Aux(i),ttUAux(i),icAux(i),groupeAux(i)
    !end do
    !do i=1,5
    !    print*,(vaxAux(i,j),j=1,nvaAux)
    !    print*,'--'
    !end do
    !STOP
    
    
!cpm
    istopp=0
    time = 0.d0

!AD:add for new marq
    epsa=EPS(1) !1.d-4
    epsb=EPS(2) !1.d-4
    epsd=EPS(3) !1.d-4
!AD:end

    timedep = timedep0
    npbetatps = (nbinnerknots0+qorder0-1)*sum(filtretps0)
    nbinnerknots = nbinnerknots0
    qorder = qorder0
    logNormal = logNormal0
    intcens = intcensAux
    typeof = typeof0
    model=4
    npmax=np
    if (typeof .ne. 0) then
        nbintervR = nbintervR0
    end if
!     shapeweib = 0.d0
!     scaleweib = 0.d0

    NSUJETMAX=nsujetAux
    allocate(t0(nsujetmax),t1(nsujetmax),tU(nsujetmax),c(nsujetmax),stra(nsujetmax),g(nsujetmax))

    allocate(RisqCumul(nsujetmax))
    if (typeof == 0) then
        allocate(nt0(nsujetmax),nt1(nsujetmax),ntU(nsujetmax))
        nt0=0
        nt1=0
        ntU=0
    end if
    c=0
    g=0

    ndatemax=2*nsujetAux+sum(icAux) ! on ajoute le nombre de temps d'entree en plus : les tU
                                    ! c'est a dire le nombre de censures par intervalle

    allocate(date(0:ndatemax))
    date=0.d0

    allocate(mm3(ndatemax),mm2(ndatemax),mm1(ndatemax),mm(ndatemax),im3(ndatemax), &
    im2(ndatemax),im1(ndatemax),im(ndatemax))    

    mm=0.d0
    mm1=0.d0
    mm2=0.d0
    mm3=0.d0

    ngmax=ngAux

!Al: utile pour le calcul d'integrale avec distribution log normale
    allocate(res3(ngmax),res5(nsujetmax),res1(ngmax))
!Al

    allocate(nig(ngmax))
    nig=0
    allocate(Residus(ngmax),cumulhaz(ngmax),vuu(1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    j=0
    nbou2=0
    id=1    
    cptni=0
    cptni1=0
    cptni2=0
    cptbiais=0
!
!  --------- Added January 05 JRG  
!

    pe=0.d0
    nsujet=0

    ndate=0
    nst=0
    ibou=1
    ij=0
    kk=0
    ni=0
    cpt=0

!    b=0.d0 on veut pouvoir laisser l'utilisateur initialiser
    resOut=0.d0

    nsujet=nsujetAux
    ng=ngAux
    icen=icenAux
    nst=nstAux
    effet=effetAux 
    nz=nzAux

    nvarmax=nvaAux
    allocate(etaT(nstAux),betaT(nstAux))
    nstRec=0 ! dans marq98j de distinguer le modele shared du modele joint
    shapeweib = 0.d0
    scaleweib = 0.d0
    lamTOut=0.d0
    suTOut=0.d0

    allocate(filtretps(nvarmax),betatps(nvarmax),betatps3(nvarmax))
    filtretps = filtretps0

    if (noVar.eq.1) then 
        do i=1,nvaAux
            filtre(i)=0
        enddo  
        nva=0  
    else
        do i=1,nvaAux
            filtre(i)=1
        enddo  
        nva=nvaAux  
    end if  

    ver=nvaAux
    nvarmax=ver
    allocate(ve(nsujetmax,nvarmax))

    ve=0.d0
    res=0.d0
    v=0.d0
 
    AG=AGAux
    maxiter=maxitAux 

!**************************************************
!**************************************************
!********************* prog spline****************

    res01=0.d0
     
!------------  lecture fichier -----------------------

    maxt = 0.d0
    mint = 0.d0
    
    cpt = 0
    k = 0
    cptstr1 = 0
    cptstr2 = 0

    do i = 1,nsujet
        str=1

        if (i.eq.1) then
            mint = tt0Aux(i) ! affectation du min juste une fois
        endif

        if(nst.ge.2)then
            tt0=tt0Aux(i)
            tt1=tt1Aux(i)
            ttU=ttUAux(i)
            ic=icAux(i)
            groupe=groupeAux(i)
            str=INT(strAux(i))
            
            do j=1,nva
                vax(j)=vaxAux(i,j)
            enddo
        else
            
            tt0=tt0Aux(i)
            tt1=tt1Aux(i)
            ttU=ttUAux(i)
            ic=icAux(i)
            groupe=groupeAux(i)
            
            do j=1,nva
                vax(j)=vaxAux(i,j)
            enddo
        endif
        k = k +1

!     essai sans troncature
!------------------   observation c=1
        if(ic.eq.1)then
            cpt = cpt + 1
            c(k)=1

!             if(str.eq.1)then
!                 stra(k) = 1
!                 cptstr1 = cptstr1 + 1
!             else
!                 if(str.eq.2)then
!                     stra(k) = 2
!                     cptstr2 = cptstr2 + 1
!                 endif
!             endif
            stra(k)=str ! pas besoin de if et des compteurs

            t0(k) = tt0
            t1(k) = tt1
            tU(k) = ttU
            g(k) = groupe

! nb de dc dans un groupe
! attention on a un nombre de groupes, mais le numero de groupe peut depasser la liste
            nig(groupe) = nig(groupe)+1

            iii = 0
            do ii = 1,ver
                if(filtre(ii).eq.1)then
                iii = iii + 1
                ve(i,iii) = vax(ii)
                endif
            end do

        else 
!------------------   censure a droite  c=0
            if(ic.eq.0)then
                c(k) = 0 
!                 if(str.eq.1)then
!                     stra(k) = 1
!                     cptstr1 = cptstr1 + 1
!                 else
!                     if(str.eq.2)then
!                         stra(k) = 2
!                         cptstr2 = cptstr2 + 1
!                     endif
!                 endif
                stra(k)=str !en plus strates A.Lafourcade 05/2014
                iii = 0

                do ii = 1,ver
                    if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) =vax(ii)
                    endif
                end do 

                t0(k) = tt0
                t1(k) = tt1
                tU(k) = ttU
                g(k) = groupe
            endif
        endif

        if (maxt.lt.t1(k))then
            maxt = t1(k)
        endif

        if ((maxt.lt.tU(k)).and.(tU(k).ne.t1(k))) then
            maxt = tU(k)
        endif

        if (mint.gt.t0(k)) then
            mint = t0(k)
        endif
        
    end do

!AD:
!    if (typeof .ne. 0) then
        cens = maxt
!    end if
!Ad

! rajout Alexandre 18/05/2012
! creation de d : nombre de censures par intervalle dans chaque groupe
    allocate(d(ngmax))
    d = 0
    do i=1,nsujet
        if (c(i).eq.1) then
            d(g(i)) = d(g(i)) + 1
        endif
    enddo
    
    dmax = 0
    do i=1,ng
        if (dmax.lt.d(i)) then
            dmax = d(i) ! max number of subject avec event=1( interval censored) by group
        endif
    enddo
!    write(*,*)'dmax dans frailtypack.f90',dmax
    

! %%%%%%%%%%%%% SANS EFFET ALEATOIRE %%%%%%%%%%%%%%%%%%%%%%%%%

    nsujet = k
    
    if (typeof == 0) then    
        nz1=nz
        nz2=nz
        
        if(nz.gt.20)then
            nz = 20
        endif
        
        if(nz.lt.4)then
            nz = 4
        endif
!***************************************************
!-------------- zi- ----------------------------------
!      construire vecteur zi (des noeuds)
    end if
    
    min = 1.d-10
    max = maxt

    do i = 1,(2*nsujet+sum(icAux)) ! changement comme indique plus haut
        do k = 1,nsujet
            if (t0(k).ge.min) then
                if(t0(k).lt.max)then
                    max = t0(k)
                endif
            endif
            
            if (t1(k).ge.min) then
                if(t1(k).lt.max)then
                    max = t1(k)
                endif
            endif
            
            if (tU(k).ne.t1(k)) then
                if((tU(k).ge.min))then
                    if (tU(k).lt.max) then
                        max = tU(k)
                    endif
                endif
            endif
        end do
        
        aux(i) = max
        min = max + 1.d-12 ! pour virer les doublons
        max = maxt
    end do
    
    date(1) = aux(1)

    k = 1
    do i=2,(2*nsujet+sum(icAux))
        if(aux(i).gt.aux(i-1))then
            k = k+1
            date(k) = aux(i)
        endif
    end do
    
    if(typeof == 0) then

! Al:10/03/2014 emplacement des noeuds splines en percentile (sans censure par intervalle)
        if(equidistant.eq.0) then ! percentile
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
        else ! equidistant
            nzmax=nz+3
            allocate(zi(-2:nzmax))
            ndate = k

            zi(-2) = mint! VR 20-fev-15 date(1)
            zi(-1) =  mint!date(1)
            zi(0) =  mint!date(1)
            zi(1) =  mint!date(1)
            h = (date(ndate)- mint)/dble(nz-1)
            do i=2,nz-1
                zi(i) =zi(i-1) + h
            end do  


            zi(nz) = maxt !date(ndate)
            zi(nz+1) = maxt !zi(nz)
            zi(nz+2) = maxt !zi(nz)
            zi(nz+3) = maxt !zi(nz)
            ziOut = zi
        endif
!--------- affectation nt0,nt1,ntU----------------------------
        indictronq=0
        do i=1,nsujet

            if(t0(i).eq.0.d0)then
                nt0(i) = 0
            endif
            if(t0(i).ne.0.d0)then
                indictronq=1
            endif
!Al:
            !if ((intcens.eq.1).and.(t1(i).eq.0.d0)) then
            !    nt1(i) = 0
            !endif
!Al:
            do j=1,ndate

                if(date(j).eq.t0(i))then
                    nt0(i)=j
                endif

                if(date(j).eq.t1(i))then
                    nt1(i)=j
                endif

                if(date(j).eq.tU(i))then
                    ntU(i)=j
                endif
            end do
        end do

!--------- affectation des vecteurs de splines -----------------
        n = nz+2

        call vecsplis(n,ndate)

        allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax),m3m1(nzmax),m3m(nzmax), &
        m2m1(nzmax),m2m(nzmax),m1m(nzmax))

        call vecpens(n)
    end if

    allocate(k0T(nst)) !k0T en plus strates A.Lafourcade 05/2014
    allocate(H_hess(npmax,npmax),Hspl_hess(npmax,npmax), &
    hess(npmax,npmax),I_hess(npmax,npmax))

    H_hess=0.d0
    I_hess=0.d0
    Hspl_hess=0.d0
    hess=0.d0
    HI=0.d0

!------ initialisation des parametres

    ! savoir si l'utilisateur entre des parametres initiaux
    if (sum(b).eq.0.d0) then
        b = 1.d-1
    else if (sum(b(np-nva+1:np)).eq.0.d0) then
        b(1:np-nva-1) = 1.d-1
        b(np-nva+1:np) = 1.d-1
    else if (b(np-nva).eq.0.d0) then
        b(1:np-nva-1) = 1.d-1
        b(np-nva) = 1.d-1
    else
        b(1:np-nva-1) = 1.d-1
    endif

    if (typeof == 1) then
!         if (nst == 2) then
!             b(nbintervR+1:2*nbintervR)=0.2d0
!         end if
        b(nbintervR+1:nst*nbintervR)=0.2d0
    end if
!    Esto se cambia ya que ahora xmin1 es kappa1

    xmin1=axT(1)
    if (nst.ge.2) then
        xmin2=axT(2)
    else
        xmin2=0.d0
    endif
    if(nst.ge.2)then
        do l=2,nst
            xminT(l)=axT(l)
        end do
    end if !fin strates

    ent = 0
    entTPS = 0
    indd = 1

    allocate(knotsTPS(-qorder+1:nbinnerknots+qorder))
    allocate(innerknots(nbinnerknots))

    if (typeof == 1) then

        if (intcens.eq.0) then !! enlever le piecewise-per pour la censure par intervalle
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

            n = nbintervR
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

            ent=int(nbrecu/nbintervR)
            entTPS=int(nbrecu/(nbinnerknots+1))! ajout 24/02/2012
        endif

        allocate(ttt(0:nbintervR))

        ttt(0) = mint !0.d0 pour prendre en compte une eventuelle troncature

        ttt(nbintervR)=cens
        j=0
        do j=1,nbintervR-1
            if (equidistant.eq.0) then ! ici se fait la difference entre piecewise-per et equi
                ttt(j) = (t2(ent*j)+t2(ent*j+1))/(2.d0)
            else
                ttt(j) = mint + ((cens-mint)/nbintervR)*j
            endif
        end do

        if (intcens.eq.0) then
            ! ajout beta tps 24/02/2012
            do j=1,nbinnerknots
!                if (equidistantTPS.eq.0) then
                    knotsTPS(j)=(t2(entTPS*j)+t2(entTPS*j+1))/(2.d0)
!                else
!                    knotsTPS(j)=(cens/(nbinnerknots+1))*j
!                endif
            end do

            knotsTPS(-qorder+1:0)=ttt(0)
            knotsTPS(nbinnerknots+1:nbinnerknots+qorder)=cens
!            write(*,*),'knotsTPS',knotsTPS
            time = ttt

            deallocate(t2)
        endif
!------- FIN RECHERCHE DES NOEUDS

    end if

    if (typeof.ne.1) then
        
        if (intcens.eq.0) then !! enlever le piecewise-per pour la censure par intervalle
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

            !n = nbintervR
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
        
            entTPS=int(nbrecu/(nbinnerknots+1))! ajout 24/02/2012
!        endif

        ! ajout beta tps 24/02/2012
        do j=1,nbinnerknots
!            if (equidistant.eq.0) then
                !print*,knotsTPS(j)
                knotsTPS(j)=(t2(entTPS*j)+t2(entTPS*j+1))/(2.d0)
!            else
!                knotsTPS(j)=(cens/(nbinnerknots+1))*j
!            endif
        end do
!        censTPS=cens
        knotsTPS(-qorder+1:0)= mint !t2(1)
        knotsTPS(nbinnerknots+1:nbinnerknots+qorder)=cens
!        write(*,*),'knotsTPS',knotsTPS

!        if (intcens.eq.0) then
            deallocate(t2)
        endif
!------- FIN RECHERCHE DES NOEUDS
    end if

!***********************************************************
!************** NEW : cross validation  ***********************
!       sur une seule strate, sans var expli , sans frailties ****
!***********************************************************
   

    nvacross=nva !pour la recherche du parametre de lissage sans var expli
    nva=0
    effetcross=effet
    effet=0
    nstcross=nst
    nst=1

    do l=1,nsujet  
        stracross(l)=stra(l)
    end do
    
    do l=1,nsujet
        stra(l)=1
    end do

    if(typeof == 0) then
        if(irep1.eq.1)then   !pas recherche du parametre de lissage
            xmin1 = dsqrt(xmin1)
            auxi = estimvs(xmin1,n,b,y,ddl,ni,res)

            if (ni.ge.maxiter) then
    !            write(*,*) ' '
    !            write(*,*) 'no convergence with the chosen smoothing parameter'
                do i=1,nz+2
                    b(i)=1.d-1
                end do
                xmin1 = sqrt(10.d0)*xmin1
                auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
                if (ni.lt.maxiter) then
    !                write(*,*)' '
    !                write(*,*)'Value of the smoothing parameter :' &
    !                ,real(xmin1*xmin1), '  DoF :',-ddl
    !                write(*,*)' '
    !                write(*,*)'Log-vraisemblance :',res
                else
                    do i=1,nz+2
                        b(i)=1.d-1
                    end do
                    
                    xmin1 = sqrt(10.d0)*xmin1
                    auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
                    if (ni.lt.maxiter) then
    !                    write(*,*)' '
    !                    write(*,*)'Value of the smoothing parameter :' &
    !                    ,real(xmin1*xmin1), '  DoF :',-ddl
    !                    write(*,*)' '
    !                    write(*,*)'Log-vraisemblance :',res    
                    endif
                endif
            else
    !            write(*,*)' '
    !            write(4,*)'Value of the smoothing parameter :' &
    !            ,real(xmin1*xmin1), '  DoF :',-ddl
    !            write(4,*)' '
    !            write(*,*)' '
    !            write(*,*)'                Searching smoothing parameter'
    !            write(*,*)' Value of the smoothing parameter :' &
    !            ,real(xmin1*xmin1), '  DoF :',-ddl
    !            write(*,*)' '
    !            write(*,*)'Log-vraisemblance :',res
            endif
!---------------------------------------------------
        else                   !recherche du parametre de lissage
        
            if(xmin1.le.0.d0)then
                xmin1 = 1.d0
            endif 
    !        write(*,*)' '
    !        write(*,*)'                Searching smoothing parameter'
    !        write(*,*)' '
            xmin1 = dsqrt(xmin1)
            auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
    
            if(ddl.gt.-2.5d0)then
    
                xmin1 = dsqrt(xmin1)
            
                auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
                if(ddl.gt.-2.5d0)then
                    xmin1 = dsqrt(xmin1)
                    auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
    
                    if(ddl.gt.-2.5d0)then
                        xmin1 = dsqrt(xmin1)
                        auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
    
                        if(ddl.gt.-2.5d0)then
                            xmin1 = dsqrt(xmin1)
                            auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
                
                            if(ddl.gt.-2.5d0)then
                                xmin1 = dsqrt(xmin1)
                            endif
                        endif
                    endif
                endif
            endif
        
            if (ni.ge.maxiter) then
                do i=1,nz+2
                    b(i)=1.d-1
                end do
                xmin1 = sqrt(10.d0)*xmin1
                auxi = estimvs(xmin1,n,b,y,ddl,ni,res)
    
                if (ni.ge.maxiter) then
                    do i=1,nz+2
                    b(i)=1.d-1
                end do
                    xmin1 = sqrt(10.d0)*xmin1
                endif
            endif
            
            ax = xmin1
            bx = xmin1*dsqrt(1.5d0)
            call mnbraks(ax,bx,cx,fa,fb,fc,b,n)
        
            tol = 0.001d0
        
            res = goldens(ax,bx,cx,tol,xmin1,n,b,y,ddl)
    !        write(4,*)'******************************************* '
    !        write(4,*)'Best smoothing parameter',real(xmin1*xmin1), '  DoF :',-ddl
    !        write(4,*)'******************************************* '
    !        write(*,*)'Best smoothing parameter',real(xmin1*xmin1), '  DoF :',-ddl
            
            auxkappa(1)=xmin1*xmin1
            auxkappa(2)=0.d0
            if (logNormal.eq.0) then
                if (intcens.eq.1) then
                    call marq98j(auxkappa,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines_intcens)
                else
                    call marq98j(auxkappa,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines)
                endif
            else
                call marq98j(auxkappa,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines_log)
            endif

            if (istop.ne.1) then
                istopp(1)=1
                goto 1000
            end if
            if (ni.ge.maxiter) then
    !            write(*,*)' '
    !            write(*,*) 'non convergence'
            else
    !            write(*,*)' '
    !            write(*,*)'Log-vraisemblance :',res
            endif

        endif
    !    write(4,*)'nombre de noeuds:',nz
    end if

    nva=nvacross ! pour la recherche des parametres de regression
    nst=nstcross ! avec stratification si necessaire
    effet=effetcross ! avec effet initial

    do l=1,nsujet
        stra(l)=stracross(l) !retablissement stratification
    end do

    if (typeof .ne. 0) then
        allocate(vvv((npmax*(npmax+1)/2)))
    end if
    
    if (typeof == 0) then
        k0(1) = xmin1*xmin1
        if(nst.eq.2)then
            k0(2) = xmin2
        endif
        k0T(1) = xmin1*xmin1 !en plus strates A.Lafourcade 05/2014
        if(nst.ge.2)then
            do l=2,nst
                k0T(l) = xminT(l)
            end do
        endif
    end if
!------------  indicateur d'effet aleatoire ou non dans le modele
    !write(*,*)'np',np
    ca=0.d0
    cb=0.d0
    dd=0.d0
!    if (typeof .ne. 0) then
!        allocate(kkapa(2))
!    end if
!    print*,"appel de funcpa"
    select case(typeof)
        case(0)
            if (timedep.eq.0) then
                if (logNormal.eq.0) then
                    if (intcens.eq.1) then
                        call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines_intcens)
                    else
                        call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines)
                    endif
                else
                    call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines_log)
                endif
            else
                call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpas_tps)
            endif
        case(1)
            allocate(betacoef(nst*nbintervR))

            if (timedep.eq.0) then
                if (logNormal.eq.0) then
                    if (intcens.eq.1) then
                        call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpascpm_intcens)
                    else
                        call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpascpm)
                    endif
                else
                    call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpascpm_log)
                endif
            else
                call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpas_tps)
            endif
        case(2)
            if (timedep.eq.0) then
                if (logNormal.eq.0) then
                    if (intcens.eq.1) then
                        call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpasweib_intcens)
                    else
                        call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpasweib)
                    endif
                else
                    call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpasweib_log)
                endif
            else
                call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpas_tps)
            endif
    end select

!    if (typeof .ne. 0) then
!        deallocate(kkapa)
!    end if
!Al:

    EPS(1) = ca
    EPS(2) = cb
    EPS(3) = dd
!Al:

    if (istop.ne.1) then
        istopp(2)=1
        goto 1000
    end if


!    j=(np-nva)*(np-nva+1)/2

    call multis(I_hess,H_hess,np,np,np,IH)
    call multis(H_hess,IH,np,np,np,HIH)

!    if(effet.eq.1)then
!        j=(np-nva)*(np-nva+1)/2
    !    write(*,*)''
    !    write(*,*)'variance non corrigee H pour theta:'
    !    write(*,*)((2.d0*b(np-nva))**2)*H_hess(np-nva,np-nva)
    !    write(*,*)'variance corrigee HIH pour theta:'
    !    write(*,*)((2.d0*b(np-nva))**2)*HIH(np-nva,np-nva)
!    endif

!    if(effet.eq.1)then
!        j=(np-nva)*(np-nva+1)/2
!    endif

    if(effet.eq.1.and.ier.eq.-1)then
        v((np-nva)*(np-nva+1)/2)=10.d10
    endif

    res01(effet+1)=res

    !write(4,*)'valeur de ni',ni 
    !write(*,*)'valeur de ni',ni

!    j=(np-nva)*(np-nva+1)/2

! --------------  Lambda and survival estimates JRG January 05

    select case(typeof)
        case(0)
!             call distancessplines(nz1,nz2,b,effet,mt,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
            call distancessplines(nz1,b,effet,mt,xTOut,lamTOut,suTOut)!en plus strates A.Laf 05/2014
        case(1)
!             Call distancecpm(b,nbintervR*nst,mt,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
            Call distanceScpm(b,nbintervR*nst,mt,xTOut,lamTOut,xSuT,suTOut)!en plus strates A.Laf 05/2014
        case(2)
            typeof2 = 1            
!             Call distanceweib(b,np,mt,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
            Call distanceSweib(b,np,mt,xTOut,lamTOut,xSuT,suTOut)!en plus strates A.Laf 05/2014            
    end select

    resOut=res

!     if (nst == 1) then
!         scaleweib(1) = etaR !betaR
!         shapeweib(1) = betaR !etaR
!         scaleweib(2) = 0.d0
!         shapeweib(2) = 0.d0
!     else
!         scaleweib(1) = etaR !betaR
!         shapeweib(1) = betaR !etaR
!         scaleweib(2) = etaD !betaD
!         shapeweib(2) = betaD !etaD
!     end if
    do l=1,nst !en plus strates A.Lafourcade 05/2014
        scaleweib(l) = etaT(l)
        shapeweib(l) = betaT(l)
    end do
    

!AD:add LCV
!LCV(1) The approximate like cross-validation Criterion
!LCV(2) Akaike information Criterion
!     calcul de la trace, pour le LCV (likelihood cross validation)
    LCV=0.d0
    if (typeof == 0) then
!        write(*,*)'The approximate like cross-validation Criterion in the non parametric case'
        call multis(H_hess,I_hess,np,np,np,HI)
        do i =1,np
            LCV(1) = LCV(1) + HI(i,i)
        end do
        LCV(1) = (LCV(1)-resnonpen) / nsujet
    else
!        write(*,*)'=========> Akaike information Criterion <========='
!        LCV(2) = 2.d0 * np - 2.d0 * resOut
        LCV(2) = (1.d0 / nsujet) *(np - resOut)
!        write(*,*)'======== AIC :',LCV(2)
    end if


!AD:end

    do ss=1,npmax
        do sss=1,npmax
            HIHOut(ss,sss) = HIH(ss,sss)
            H_hessOut(ss,sss)= H_hess(ss,sss)
        end do
    end do


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!    Bias and Var eliminated
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!AD:
1000    continue

    if (timedep.eq.1) then

        BetaTpsMat = 0.d0

! CALCUL DES ESTIMATIONS DES beta DEPENDANTS DU TEMPS

        if(nva.gt.0)then
            call drawTimeCoef(np,b,nva,filtretps,BetaTpsMat)
        endif

    endif

!---------------------Calcul residus de martingales

    deallocate(I_hess,H_hess)

    coefBeta = 0.d0
    Xbeta = 0.d0

    if (timedep.eq.0) then
        coefBeta(1,:) = b((np-(nva+npbetatps)+1):np)
        Xbeta = matmul(coefBeta,transpose(ve))
    else
        do j=1,nsujet
            p=1
            call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,nbinnerknots+qorder,t1(j), &
            innerknots,boundaryknots,basis)
            do i=1,nva
                coefBeta2 = 0.d0
                if (filtretps(i).eq.1) then
                    do k=-qorder+1,nbinnerknots
                        coefBeta2 = coefBeta2 + b(np-(nva+npbetatps)+p-1+k+qorder)*basis(k+qorder)
                    end do
                else
                    coefBeta2 = b(np-(nva+npbetatps)+p)
                endif
                Xbeta(1,j) = Xbeta(1,j) + coefBeta2*dble(ve(j,i))
                p=p+filtretps(i)*(nbinnerknots+qorder-1)+1
            end do
        enddo
    endif

    if(typeof==0) then

        if((istopp(1) == 0).and.(istopp(2) == 0)) then
            !deallocate(I_hess,H_hess)

            allocate(vecuiRes(ng),vres((1*(1+3)/2)),I_hess(1,1),H_hess(1,1),post_esp(ng),post_SD(ng))
            effetres = effet

            if (effet == 1) then
                Call ResidusMartingale(b,np,funcpasres,Resmartingale,frailtypred,frailtyvar,frailtysd)
            else
                Resmartingale=0.d0
                frailtypred=0.d0
                frailtyvar=0.d0
                frailtysd=0.d0
            end if

            do i=1,nsujet
                if (effet == 1) then
                    if (logNormal.eq.0) then
                        linearpred(i)=Xbeta(1,i)+dlog(frailtypred(g(i)))
                    else
                        linearpred(i)=Xbeta(1,i)+frailtypred(g(i))
                    endif
                else
                    linearpred(i)=Xbeta(1,i)
                    martingaleCox(i) = c(i) - RisqCumul(i)
                end if
            end do

            deallocate(I_hess,H_hess,vres,vecuiRes,post_esp,post_SD)
        else
            !deallocate(I_hess,H_hess)
            Resmartingale=0.d0
            frailtypred=0.d0
            linearpred=0.d0
            martingaleCox=0.d0
        end if

    else
        if((istopp(2) == 0)) then

            !deallocate(I_hess,H_hess)

            allocate(vecuiRes(ng),vres((1*(1+3)/2)),I_hess(1,1),H_hess(1,1),post_esp(ng),post_SD(ng))    
            effetres = effet

            if (effet == 1) then
                Call ResidusMartingale(b,np,funcpasres,Resmartingale,frailtypred,frailtyvar,frailtysd)
            else
                Resmartingale=0.d0
                frailtypred=0.d0
                frailtyvar=0.d0
                frailtysd=0.d0
            end if

            do i=1,nsujet
                if (effet == 1) then
                    linearpred(i)=Xbeta(1,i)+dlog(frailtypred(g(i)))
                else
                    linearpred(i)=Xbeta(1,i)
                    martingaleCox(i) = c(i) - RisqCumul(i)
                end if
            end do

            deallocate(I_hess,H_hess,vres,vecuiRes,post_esp,post_SD)
        else
            !deallocate(I_hess,H_hess)
            Resmartingale=0.d0
            frailtypred=0.d0
            linearpred=0.d0
            martingaleCox=0.d0
        end if
    end if

! ! Al : 25/03/2014 Score Test for frailties (Commenges/Andersen)
!     if (effet.eq.1) then    
!         scoreT = scoretest(b,frailtypred)
!         print*,scoreT
!     endif

    deallocate(t0,t1,tU,c,d,stra,g,nig,ve,Hspl_hess,hess, &
    mm3,mm2,mm1,mm,im3,im2,im1,im,date,cumulhaz,Residus,vuu, &
    RisqCumul,res5,res3,res1,filtretps,knotsTPS,innerknots,betatps,betatps3)

    if (typeof == 0) then
        deallocate(nt0,nt1,ntU,zi)
        deallocate(m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
    end if

    if (typeof == 1) then
        deallocate(ttt,vvv,betacoef)
    end if

    if (typeof == 2) then
        deallocate(vvv)
    end if
    deallocate(etaT,betaT,k0T)!en plus strates A.Lafourcade 05/2014

    return

    end subroutine frailpenal


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!CCCCCCCCCCC!**********SUBROUTINES******  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



!========================== VECSPLI ==============================
    subroutine vecsplis(n,ndate)

    !use tailles,only:ndatemax,NSUJETMAX,npmax,nvarmax
    !use comon,only:ve
    use comon,only:date,zi,mm2,mm3,mm1,mm,im3,im2,im1,im

    implicit none
    
    integer::n,ndate,i,j,k
    double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2

!---------  calcul de u(ti) ---------------------------

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
    
    end subroutine vecsplis

!========================== VECPEN ==============================
    subroutine vecpens(n)
    !use tailles,only:ndatemax,npmax
    !use comon,only:date
    use comon,only:zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
    
    implicit none
        
    integer::n,i
    double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
    double precision::a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x
     
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
        ))+x*(4.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i-2)+4.d0*zi(i+1)*zi(i-2)))/(a2*a2)))
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
        m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i) &   
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
        *zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1) &   
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
        +2.d0*zi(i)*zi(i)))/(c1*a0)))
    
    end do
    
    end subroutine vecpens



!==========================  SUSP  ====================================
    subroutine susps(x,the,n,su,lam,zi)

    !use tailles,only:ndatemax
    use tailles,only:npmax

    implicit none

    integer::j,k,n,i
    double precision::x,ht,ht2,h2,som,lam,su,htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3, &
    ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,gl,hh
    double precision,dimension(-2:npmax)::zi,the


    som = 0.d0
    gl = 0.d0
    
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

    end subroutine susps

!==========================  COSP  ====================================
! calcul les points pour les fonctions
! et leur bandes de confiance

    subroutine cosps(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)

    !use tailles,only:ndatemax
    use tailles,only:npmax

    implicit none

    integer::j,k,n,i
    double precision::x,ht,ht2,h2,som,lam,su,binf,bsup,lbinf,lbsup,pm, &
    htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2, &
    h,gl,hh
    double precision,dimension(-2:npmax)::the,zi
    double precision,dimension(npmax,npmax)::y


    j=0
    som = 0.d0
    gl = 0.d0 

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

    call confs(x,j,n,y,pm,zi)

    binf = dexp(-gl - 1.96d0*pm)
    su  = dexp(-gl)
    bsup = dexp(-gl + 1.96d0*pm)

    call conf1s(x,j,n,y,pm,zi)
    lbinf = lam - 1.96d0*pm
    lbsup = lam + 1.96d0*pm
    
    return
    
    end subroutine cosps


!=====================  CONF1  =============================
    subroutine conf1s(x,ni,n,y,pm,zi)

    !use tailles,only:ndatemax
    use tailles,only:npmax

    implicit none

    integer::ni,i,n,j
    double precision::mmsps,x,pm,res
    double precision,dimension(-2:npmax)::zi
    double precision,dimension(npmax)::vecti,aux
    double precision,dimension(npmax,npmax)::y

    do i=1,n
        vecti(i) = mmsps(x,ni,i,zi)
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

    res=-res
    pm = dsqrt(res)

    end subroutine conf1s


!=====================  CONF  =============================
    subroutine confs(x,ni,n,y,pm,zi)

    !use tailles,only:ndatemax
    use tailles,only:npmax

    implicit none

    integer::ni,i,n,j
    double precision::isps,x,pm,res
    double precision,dimension(-2:npmax)::zi
    double precision,dimension(npmax)::vecti,aux
    double precision,dimension(npmax,npmax)::y

    do i=1,n
        vecti(i) = isps(x,ni,i,zi)
    end do

    do i=1,n
        aux(i) = 0.d0
        do j=1,n
            aux(i) = aux(i) - y(i,j)*vecti(j) ! multiplication matrice H*BaseSpl=aux
        end do
    end do

    res = 0.d0
    do i=1,n
        res = res + aux(i)*vecti(i) ! BaseSpl*aux
    end do

    res=-res
    pm = dsqrt(res)    

    end subroutine confs


!==========================   ISP   ==================================
    double precision function isps(x,ni,ns,zi)

    !use tailles,only:ndatemax
    use tailles,only:npmax

    implicit none

    integer::ni,ns
    double precision::val,mmsps,x
    double precision,dimension(-2:npmax)::zi



    if(x.eq.zi(ni))then
        if(ni.le.ns-3)then
            val = 0.d0
        else
            if(ni.le.ns-2)then
                val = ((zi(ni)-zi(ni-1))*mmsps(x,ni,ns,zi))*0.25d0
            else
                if (ni.eq.ns-1)then
                    val = ((zi(ni)-zi(ni-2))*mmsps(x,ni,ns,zi)+ &
                      (zi(ni+3)-zi(ni-1))*mmsps(x,ni,ns+1,zi))*0.25d0
                else
                    if(ni.eq.ns)then
                        val = ((zi(ni)-zi(ni-3))*mmsps(x,ni,ns,zi)+ &
                        (zi(ni+2)-zi(ni-2))*mmsps(x,ni,ns+1,zi) &
                        +(zi(ni+3)-zi(ni-1))*mmsps(x,ni,ns+2,zi))*0.25d0
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
                val = (x-zi(ni))*mmsps(x,ni,ns,zi)*0.25d0
            else
                if(ni.eq.ns-2)then
                    val = ((x-zi(ni-1))*mmsps(x,ni,ns,zi)+ &
                    (zi(ni+4)-zi(ni))*mmsps(x,ni,ns+1,zi))*0.25d0
                else
                    if (ni.eq.ns-1)then
                        val =((x-zi(ni-2))*mmsps(x,ni,ns,zi)+ &
                        (zi(ni+3)-zi(ni-1))*mmsps(x,ni,ns+1,zi) &
                        +(zi(ni+4)-zi(ni))*mmsps(x,ni,ns+2,zi))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val =((x-zi(ni-3))*mmsps(x,ni,ns,zi)+ &
                            (zi(ni+2)-zi(ni-2))*mmsps(x,ni,ns+1,zi) &
                            +(zi(ni+3)-zi(ni-1))*mmsps(x,ni,ns+2,zi) &
                            +(zi(ni+4)-zi(ni))*mmsps(x,ni,ns+3,zi))*0.25d0
                        else
                            val = 1.d0
                        endif
                    endif
                endif
            endif
        endif
    endif

    isps = val

    return

    end function isps


!==========================  MMSP   ==================================
    double precision function mmsps(x,ni,ns,zi)

    !use tailles,only:ndatemax
    use tailles,only:npmax

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
                        +((4.d0*((x-zi(ni-1))*(zi(ni+2)-x) & 
                        *(zi(ni+1)-x))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni))))) &
                        +((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x) & 
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

    mmsps = val

    return

    end function mmsps



!========================          MNBRAK         ===================
    subroutine mnbraks(ax,bx,cx,fa,fb,fc,b,n)

    use tailles,only:npmax

    implicit none

    double precision::ax,bx,cx,fa,fb,fc,aux,res
    double precision,dimension(npmax)::b
    double precision,dimension(npmax,npmax)::y
    double precision::estimvs,gold,glimit,tiny
    parameter (gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
    double precision::dum,fu,q,r,u,ulim
    integer::n,ni

    fa = estimvs(ax,n,b,y,aux,ni,res)
    fb = estimvs(bx,n,b,y,aux,ni,res)

    if(fb.gt.fa)then
        dum = ax
        ax = bx
        bx = dum
        dum = fb
        fb = fa
        fa = dum
    endif

    cx = bx + gold*(bx-ax)
    fc = estimvs(cx,n,b,y,aux,ni,res)

 1       if(fb.ge.fc)then
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),tiny),q-r))
        ulim = bx + glimit*(cx-bx)

        if((bx-u)*(u-cx).gt.0.d0)then
            fu = estimvs(u,n,b,y,aux,ni,res)
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
            fu = estimvs(u,n,b,y,aux,ni,res)
        else
            if((cx-u)*(u-ulim).gt.0.d0)then
                fu = estimvs(u,n,b,y,aux,ni,res)
                if(fu.lt.fc)then
                    bx = cx
                    cx = u
                    u = cx + gold*(cx-bx)
                    fb = fc
                    fc = fu
                    fu = estimvs(u,n,b,y,aux,ni,res)
                endif
            else
                if((u-ulim)*(ulim-cx).ge.0.d0)then
                    u = ulim
                    fu = estimvs(u,n,b,y,aux,ni,res)
                else
                    u = cx + gold*(cx-bx)
                    fu = estimvs(u,n,b,y,aux,ni,res)
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

    end subroutine mnbraks

!========================      GOLDEN   =========================
    double precision function goldens(ax,bx,cx,tol,xmin,n,b,y,aux)

    !use tailles,only:ndatemax
    use tailles,only:npmax

    implicit none
    
    double precision,dimension(npmax,npmax)::y
    double precision,dimension(npmax)::b
    double precision::ax,bx,cx,tol,xmin,r,c,aux,res
    parameter (r=0.61803399d0,c=1.d0-r)
    double precision::f1,f2,x0,x1,x2,x3,estimvs
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

         f1 = estimvs(x1,n,b,y,aux,ni,res)
         f2 = estimvs(x2,n,b,y,aux,ni,res)
 
 1       if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
        if(f2.lt.f1)then
            x0 = x1
            x1 = x2
            x2 = r*x1 + c*x3
            f1 = f2
            f2 = estimvs(x2,n,b,y,aux,ni,res)
        else
            x3 = x2
            x2 = x1
            x1 = r*x2+c*x0
            f2 = f1
            f1 = estimvs(x1,n,b,y,aux,ni,res)
        endif
        go to 1
    endif
    if(f1.lt.f2)then
        goldens = f1
        xmin = x1
    else
        goldens = f2
        xmin = x2
    endif

    return

    end function goldens


!========================          ESTIMV         ===================

    double precision function estimvs(k00,n,b,y,aux,ni,res)

    !use tailles,only:ndatemax,NSUJETMAX
    use tailles,only:npmax
    !use comon,only:t0,t1,c,nt0,nt1,nsujet,nva,nst,nz1,nz2,typeof
    use comon,only:ndate,k0T,date,zi,pe,effet,mm3,mm2,mm1,mm,im3,&
    im2,im1,im,intcens,logNormal

    use optim
    implicit none

    double precision,dimension(npmax,npmax)::y
    double precision,dimension((npmax*(npmax+3)/2))::v
    double precision,dimension(-2:npmax)::the
    double precision,dimension(2)::k0
    double precision,dimension(ndate)::ut,dut
    double precision,dimension(npmax)::bh,b
    double precision::res,k00,som,h1
    double precision::aux
    integer::n,ij,i,k,j,vj,ier,istop,ni
    double precision::ca,cb,dd
    double precision,external::funcpassplines,funcpassplines_intcens,funcpassplines_log

    j=0
    estimvs=0.d0
    k0T(1)= k00*k00 !en plus strates A.Lafourcade 05/2014
    k0(1) = k00*k00
    k0(2) = 0.d0

    if (logNormal.eq.0) then
        if (intcens.eq.1) then
                      
            call marq98j(k0,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines_intcens)
        else
            call marq98j(k0,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines)
        endif
    else
        call marq98j(k0,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpassplines_log)
    endif

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

        call tests(dut,k0,n,aux,y)
        estimvs = - ((res-pe)) - aux
    else
        aux = -n
    endif
!AD:
50    continue
!AD:

    return

    end function estimvs

!=================calcul de la hessienne  et de omega  ==============
    subroutine tests(dut,k0,n,res,y)

    !use tailles,only:ndatemax,NSUJETMAX
    use tailles,only:npmax
    !use comon,only:date,zi,t0,t1,c,nt0,nt1,nsujet,nva,nst
    use comon,only:ndate

    implicit none

    double precision,dimension(npmax,npmax)::hessh,hess,omeg,y
    integer,dimension(npmax)::indx
    integer::n,i,j,np
    double precision,dimension(2)::k0
    double precision,dimension(ndate)::dut
    double precision::d,res,tra


    do i = 1,n
        do j = 1,n
            hess(i,j) = 0.d0
        end do
    end do

    do i = 1,n
        do j = i,n
            call mats(hess(i,j),dut,i,j,n)
        end do
    end do

    do i = 2,n
        do j = 1,i-1
            hess(i,j)=hess(j,i)
        end do
    end do

    call calcomegs(n,omeg)

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
    
    call ludcmps(hess,n,indx,d) ! decomposition LU

    do j=1,n
        call lubksbs(hess,n,indx,y(1,j)) ! resolution AX=B
    end do

    tra = 0.d0
    do i=1,n
        do j=1,n
            tra = tra + y(i,j)*hessh(j,i)
        end do
    end do

    res = (tra)

    end subroutine tests


!=======================  CALOMEG  ===========================
    subroutine calcomegs(n,omeg)

    !use tailles,only:ndatemax
    use tailles,only:npmax
    !use comon,only:date,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m2m1,m2m,m1m
    use comon,only:m3m
      
    implicit none

    double precision,dimension(npmax,npmax)::omeg
    integer::n,i,j
    double precision::calc00s,calc01s,calc02s

    do i=1,n
        do j=1,n
            omeg(i,j)=0.d0
        end do
    end do

    omeg(1,1)=calc00s(1,n)
    omeg(1,2)=calc01s(1,n)
    omeg(1,3)=calc02s(1,n)
    omeg(1,4)=m3m(1)
    omeg(2,1)=omeg(1,2)
    omeg(2,2)=calc00s(2,n)
    omeg(2,3)=calc01s(2,n)
    omeg(2,4)=calc02s(2,n)
    omeg(2,5)=m3m(2)
    omeg(3,1)=omeg(1,3)
    omeg(3,2)=omeg(2,3)
    omeg(3,3)=calc00s(3,n)
    omeg(3,4)=calc01s(3,n)
    omeg(3,5)=calc02s(3,n)
    omeg(3,6)=m3m(3)

    do i=4,n-3
        omeg(i,i-3)=omeg(i-3,i)
        omeg(i,i-2)=omeg(i-2,i)
        omeg(i,i-1)=omeg(i-1,i)
        omeg(i,i)=calc00s(i,n)
        omeg(i,i+1)=calc01s(i,n)
        omeg(i,i+2)=calc02s(i,n)
        omeg(i,i+3)=m3m(i)
    end do

    omeg(n-2,n-5)=omeg(n-5,n-2)
    omeg(n-2,n-4)=omeg(n-4,n-2)
    omeg(n-2,n-3)=omeg(n-3,n-2)
    omeg(n-2,n-2)=calc00s(n-2,n)
    omeg(n-2,n-1)=calc01s(n-2,n)
    omeg(n-2,n)=calc02s(n-2,n)
    omeg(n-1,n-4)=omeg(n-4,n-1)
    omeg(n-1,n-3)=omeg(n-3,n-1)
    omeg(n-1,n-2)=omeg(n-2,n-1)
    omeg(n-1,n-1)=calc00s(n-1,n)
    omeg(n-1,n)=calc01s(n-1,n)
    omeg(n,n-3)=omeg(n-3,n)
    omeg(n,n-2)=omeg(n-2,n)
    omeg(n,n-1)=omeg(n-1,n)
    omeg(n,n)=calc00s(n,n)

    end subroutine calcomegs


!====================  MAT  ==================================
    subroutine mats(res,dut,k,l,n)

    !use tailles,only:npmax,ndatemax,NSUJETMAX
    !use comon,only:t0,t1,nt0,nva,nst
    use comon,only:date,zi,c,nt1,nsujet,ndate    

    implicit none

    double precision::res,res1,msps,aux2,u2
    double precision,dimension(ndate)::dut
    integer::k,l,j,ni,n,i

!--------- calcul de la hessienne ij ------------------

    res = 0.d0
    res1 = 0.d0
    do i=1,nsujet
        if(c(i).eq.1)then  !event
            u2 = dut(nt1(i))
            do j = 2,n-2
                if((date(nt1(i)).ge.zi(j-1)).and. &
                    (date(nt1(i)).lt.zi(j)))then
                    ni = j-1
                endif
            end do
            if(date(nt1(i)).eq.zi(n-2))then
                ni = n-2
            endif
!------attention numero spline
            aux2 = msps(nt1(i),ni,k)*msps(nt1(i),ni,l)
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
    end subroutine mats

!==========================  MSP   ==================================
    double precision function msps(i,ni,ns)

    !use tailles,only:npmax,ndatemax
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
                    -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) &
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
                        +((4.d0*((date(i)-zi(ni-1))*(zi(ni+2)-date(i)) & 
                        *(zi(ni+1)-date(i)))/((zi(ni+2)-zi(ni-2)) &
                        *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
                        (zi(ni+1)-zi(ni))))) &
                        +((4.d0*((zi(ni+2)-date(i))*(zi(ni+2)-date(i)) & 
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

    msps = val

    return

    end function msps


!=========================  CALC00  =========================
    double precision function calc00s(j,n)

    !use tailles,only:npmax
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

    calc00s = part

    return

    end function calc00s


!=========================  CALC01  =========================
    double precision function calc01s(j,n)

    !use tailles,only:npmax
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

    calc01s = part

    return

    end function calc01s

!=========================  CALC02  =========================
    double precision function calc02s(j,n)

    !use tailles,only:npmax
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

    calc02s = part

    return

    end function calc02s


!================== multiplication de matrice  ==================

! multiplie A par B avec le resultat dans C

    subroutine multis(A,B,IrowA,JcolA,JcolB,C)

    use tailles,only:npmax

    implicit none

    integer::IrowA,JcolA,JcolB,i,j,k
    double precision::sum
    double precision,dimension(npmax,npmax) ::A,B,C

    do I=1,IrowA
        do J=1,JcolB
            sum=0
            do K=1,JcolA
                sum=sum+A(I,K)*B(K,J)
            end do
            C(I,J)=sum
        end do
    end do

    return

    end subroutine multis

!======================  LUBKSB  ======================================
    subroutine lubksbs(a,n,indx,b)

    use tailles,only:npmax

    implicit none

    integer::n,i,ii,j,ll
    integer,dimension(npmax)::indx
    double precision,dimension(npmax,npmax)::a
    double precision,dimension(npmax)::b
    double precision::sum

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

    end subroutine lubksbs

!======================  LUDCMP  ======================================
       subroutine ludcmps(a,n,indx,d)

    use tailles,only:npmax

    implicit none

    integer::n,i,imax,j,k
    integer,dimension(n)::indx
    double precision::d
    double precision,dimension(npmax,npmax)::a
    integer,parameter::nmax=500
    double precision,parameter::tiny=1.d-20
    double precision aamax,dum,sum
    double precision,dimension(nmax)::vv

    imax=0
    d = 1.d0
    do i=1,n
        aamax=0.d0
        do j=1,n
            if (dabs(a(i,j)).gt.aamax)then
                aamax=dabs(a(i,j))
            endif
        end do
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
            do i=j+1,n
                a(i,j) = a(i,j)*dum
            end do
        endif
    end do

    return

    end subroutine ludcmps

!=============================================================
! gauss hermite
! func est l integrant, ss le resultat de l integrale sur -infty , +infty

    SUBROUTINE gauherS(ss,choix)

    use tailles
    use donnees,only:x2,w2,x3,w3
    !use comon,only:auxig
    use comon,only:typeof

    Implicit none

    double precision,intent(out)::ss
    integer,intent(in)::choix

    double precision::auxfunca,func1S,func2S,func3S
    external::func1S,func2S,func3S
    integer::j

    ss=0.d0
    if (typeof.eq.0) then
        do j=1,20
            if (choix.eq.1) then
                auxfunca=func1S(x2(j))
                ss = ss+w2(j)*(auxfunca)
            else                   !choix=2, troncature
                if (choix.eq.2) then
                    auxfunca=func2S(x2(j))
                    ss = ss+w2(j)*(auxfunca)
                else
                    if (choix.eq.3) then
                        auxfunca=func3S(x2(j))
                        ss = ss+w2(j)*(auxfunca)
                    endif
                endif
            endif
        end do
    else
        do j=1,32
            if (choix.eq.1) then
                auxfunca=func1S(x3(j))
                ss = ss+w3(j)*(auxfunca)
            else                   !choix=2, troncature
                if (choix.eq.2) then
                    auxfunca=func2S(x3(j))
                    ss = ss+w3(j)*(auxfunca)
                else
                    if (choix.eq.3) then
                        auxfunca=func3S(x3(j))
                        ss = ss+w3(j)*(auxfunca)
                    endif
                endif
            endif
        end do
    endif

    return

    END SUBROUTINE gauherS

!=====================================================================

    double precision function func1S(frail)

    use tailles
    !use comon,only:stra
    use comon,only:auxig,sig2,g,res5,c

    IMPLICIT NONE

    double precision,intent(in)::frail

    integer::i
    double precision::prod
    double precision,parameter::pi=3.141592653589793d0

    prod = 1.d0
    do i=1,nsujetmax
        if (g(i).eq.auxig) then
            prod = prod * (dexp(frail)**c(i)) * &
            dexp(-dexp(frail)*res5(i))
        endif
    enddo

    func1S = prod*(1.d0/dsqrt(2.d0*pi*sig2))* &
    dexp(-(frail**2.d0)/(2.d0*sig2))

    return

    end function func1S

!=====================================================================

    double precision function func2S(frail)

    use comon,only:auxig,sig2,res3

    IMPLICIT NONE

    double precision,intent(in)::frail
    double precision,parameter::pi=3.141592653589793d0

    func2S = dexp(-dexp(frail)*res3(auxig))* &
    (1.d0/dsqrt(2.d0*pi*sig2))* &
    dexp(-(frail**2.d0)/(2.d0*sig2))

    return

    end function func2S

!=====================================================================

    double precision function func3S(frail)

    use comon,only:auxig,sig2,res1,res3,nig
    
    IMPLICIT NONE

    double precision,intent(in)::frail

!    func3S = dexp(frail*nig(auxig))* &
!    dexp(-dexp(frail)*(res1(auxig)-res3(auxig))-(frail**2.d0)/(2.d0*sig2))

    func3S = dexp(frail*nig(auxig) &
    -dexp(frail)*(res1(auxig)-res3(auxig))-(frail**2.d0)/(2.d0*sig2))

    return

    end function func3S
