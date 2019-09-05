    module perso
    implicit none
    double precision,dimension(:),allocatable,save::aux
    integer,dimension(:),allocatable,save::gaux,gnew
    integer,dimension(:),allocatable,save::stracross,filtre
    double precision,dimension(:),allocatable,save::vax
    double precision,dimension(:,:),allocatable,save:: I1_hess &
    ,H1_hess,I2_hess,H2_hess,HI1,HI2,HIH,IH,HI
    double precision,dimension(:),allocatable,save::v
    double precision,dimension(:,:),allocatable,save:: y
    end module perso

    subroutine nested(ns0,ng0,nssgbyg0,nst0,nz0,axT,tt00,tt10,ic0,groupe0, &
    ssgroupe0,nva0,str0,vax0,AG0,noVar,maxiter0,irep1,np,maxngg,b,H_hessOut,HIHOut,resOut, &
    LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out,typeof0,equidistant,nbintervR0,mt,ni, &
    cpt,ier,k0,ddl,istop,shapeweib,scaleweib,mt1,ziOut,time, &
    Resmartingale,frailtypred,frailtypredg,frailtyvar,frailtyvarg,frailtysd,frailtysdg,linearpred,EPS,nbgl)

    use tailles
    use parameters
    use optim
    use commun
    use comon,only:date,zi,t0,t1,c,nt0,nt1,ve,stra,effet,nz1,nz2,I_hess,H_hess,&
    Hspl_hess,g,nig,indictronq,ag,mm3,mm2,mm1,mm,im3,im2,im1,im,m3m3,m2m2,m1m1,mmm,m3m2, &
    m3m1,m3m,m2m1,m2m,m1m,resnonpen,nsujet,nva,ndate,nst,model,hess,typeof, &
    ttt,betacoef,typeof2,t2,vvv,nbintervR,cens,nbrecu,etaR,etaD,betaR,betaD,nb_gl
    !alpha,auxig,pe,eta,kkapa
    use perso
    use residusM

    Implicit none

    integer::ns0,ng0,nssgbyg0,nst0,np,nz0,nva0,AG0,noVar,maxiter0,mt,mt1,maxngg
    double precision,dimension(nz0+6),intent(out)::ziOut
!     double precision,intent(in)::ax1,ax2
    double precision,dimension(2),intent(in)::axT
    double precision,dimension(ns0)::tt00,tt10
    integer,dimension(ns0)::ic0,groupe0,ssgroupe0
    double precision,dimension(ns0)::str0
    double precision,dimension(ns0,nva0)::vax0
    double precision,intent(out)::resOut
    double precision,dimension(np,np),intent(out)::HIHOut,H_hessOut
    integer::ni,ier,nz,k,j,ii,iii,l
    double precision,dimension(2)::k0,res01
    double precision::h,xmin1,xmin2,maxtt,min,max,pord,mint
    double precision,dimension(np)::b
    !declaration pour R
    double precision,dimension(mt),intent(out)::x1Out,x2Out
    double precision,dimension(100)::xSu1,xSu2

    double precision,dimension(mt,3),intent(out)::lamOut,lam2Out
    double precision,dimension(mt1,3),intent(out)::suOut,su2Out
    integer ::cptni,cptni1,cptni2,ic,n,cpt,nvacross,effetcross,nstcross,irep1 &
    ,ss,sss,ngaux,istop
    integer::i,cptstr1,cptstr2,groupe,ssgroupe,ver
    double precision::str,ddl
!    real::tt0,tt1
    double precision::tt0,tt1
    integer::auxng,auxssng
    double precision::auxi,ax,bx,cx,tol,fa,fb,fc,res,goldenN,estimvN
    double precision:: bgpe,bssgpe
    double precision ::trace,trace1,trace2
    double precision,dimension(2)::auxkappa
!AD:add
    double precision,dimension(2),intent(out)::LCV,shapeweib,scaleweib
    double precision::ca,cb,dd,funcpansplines,funcpancpm,funcpanweib
    external::funcpansplines,funcpancpm,funcpanweib
!Cpm
    integer::typeof0,nbintervR0,equidistant,ent,indd
    integer,intent(in)::nbgl
    double precision::temp
    double precision,dimension(nbintervR0+1)::time
!cpm
!predictor
    double precision,dimension(nssgbyg0),intent(out)::Resmartingale
    double precision,dimension(ng0),intent(out)::frailtypred,frailtysd,frailtyvar
    double precision,external::funcpanres
    double precision,dimension(ns0),intent(out)::linearpred
    double precision,dimension(1,nva0)::coefBeta
    double precision,dimension(1,ns0)::XBeta
    integer,dimension(5)::istopp
    double precision,dimension(ng0,maxngg),intent(out)::frailtypredg,frailtysdg,frailtyvarg
    double precision,dimension(3),intent(inout)::EPS ! seuils de convergence

    istopp = 0
    indic_cumul=0
    Resmartingale=0.d0
    frailtypred=0.d0
    frailtyvar=0.d0
    frailtysd=0.d0
    linearpred=0.d0
!AD:end

    time = 0.d0
    lamOut=0.d0
    lam2Out=0.d0
    suOut=0.d0
    su2Out=0.d0
    auxng=0
    auxssng=0

    maxiter=maxiter0

    epsa=EPS(1) !1.d-4
    epsb=EPS(2) !1.d-4
    epsd=EPS(3) !1.d-4

    nb_gl = nbgl
    
    ca=0.d0
    cb=0.d0
    dd=0.d0
    typeof = typeof0
!----- Type of model joint==1
    model=3

    if (typeof == 1) then
        nbintervR = nbintervR0
    end if
    shapeweib = 0.d0
    scaleweib = 0.d0

    nsujetmax=ns0
    nsujet=ns0

    ngmax=ng0

    nssgbyg=nssgbyg0
    ngexact=nssgbyg0

    xmin1 = axT(1) !ax1
    xmin2 = axT(2) !ax2

    ver=nva0
    nvarmax=nva0

    ndatemax=2*ns0

    allocate(filtre(ver),date(ndatemax),aux(2*ns0))

    if (noVar.eq.1) then

        do i=1,nva0
            filtre(i)=0
        end do
        nva=0
    else
        do i=1,nva0
            filtre(i)=1
        end do
        nva=nva0
    end if

    if (typeof == 0) then
        allocate(mm3(ndatemax),mm2(ndatemax) &
        ,mm1(ndatemax),mm(ndatemax),im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax))
        allocate(nt0(ns0),nt1(ns0))
    end if


    allocate(t0(ns0),t1(ns0),c(ns0),stra(ns0),stracross(ns0),g(ns0),gaux(ns0),gnew(ns0))   
    allocate(ve(ns0,nvarmax))

    allocate(ssg(ns0,ngexact),nig(ngexact),mid(ngexact),n_ssgbygrp(ngexact))
    AG=AG0

    nz=nz0
    cptni=0
    cptni1=0
    cptni2=0

    res01(1)=0.d0
    res01(2)=0.d0

    effet = 1
    nig = 0
    ssg=0
    g=0

    mint = 0.d0
    maxtt = 0.d0
    cpt = 0
    k = 0
    cptstr1 = 0
    cptstr2 = 0

    allocate(vax(nva))
    nst=nst0
    do i = 1,ns0

        if (i.eq.1) then
            mint = tt00(i) ! affectation du min juste une fois
        endif

        if(nst.eq.2)then
            tt0=tt00(i)
            tt1=tt10(i)
            ic=ic0(i)
            ssgroupe=ssgroupe0(i)
            groupe=groupe0(i)
            str=str0(i)
            do j=1,nva
                vax(j)=vax0(i,j)  
            enddo
        else
            tt0=tt00(i)
            tt1=tt10(i)
            ic=ic0(i)
            ssgroupe=ssgroupe0(i)
            groupe=groupe0(i)
            str=str0(i)
            do j=1,nva
                vax(j)=vax0(i,j)  
            enddo
        endif

        k = k +1

        if(k.eq.1)then
            auxng=groupe
            auxssng=ssgroupe
            ngexact=1
            nssgexact=1
            g(k)=1
            ssg(k,g(k))=1
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
            t0(k) = tt0!dble(tt0)
            t1(k) = tt1!dble(tt1)
            if(auxng.ne.groupe)then !chgt de groupes
                ngexact=ngexact+1
                auxng=groupe
                nssgexact=nssgexact+1
                auxssng=ssgroupe
                if(k.ne.1)then
                    g(k)=g(k-1)+1
                    ssg(k,g(k))=1  !ssg(k-1,g(k-1))+1 
                endif

                goto 100
            endif

            if(auxssng.ne.ssgroupe.and.auxng.eq.groupe)then 
!     chgt de ssgroupe mais pas de groupe
                nssgexact=nssgexact+1
                auxssng=ssgroupe
                if(k.ne.1)then
                    g(k)=g(k-1)
                    ssg(k,g(k))=ssg(k-1,g(k-1))+1
                endif
                goto 100
            endif

            if(k.ne.1)then
                if(auxssng.eq.ssgroupe.and.auxng.eq.groupe)then 
!     aucun chgt de ssgroupe ni de gpe
                    if(k.ne.1)then
                        g(k)=g(k-1)
                        ssg(k,g(k))=ssg(k-1,g(k-1))
                    endif
                    goto 100
                endif
            endif
 100       continue

            nig(g(k)) = nig(g(k))+1
            iii = 0
            do ii = 1,nva
                if(filtre(ii).eq.1)then
                    iii = iii + 1
                    ve(i,iii) = vax(ii)!dble(vax(ii))
                endif
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
                do ii = 1,nva
                    if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) = vax(ii)!dble(vax(ii))

                    endif
                end do 
                t0(k) =  tt0!dble(tt0)
                t1(k) = tt1!dble(tt1)

                if(auxng.ne.groupe)then !chgt de groupes
                    ngexact=ngexact+1
                    auxng=groupe
                    nssgexact=nssgexact+1
                    auxssng=ssgroupe
                    if(k.ne.1)then
                        g(k)=g(k-1)+1
                        ssg(k,g(k))=1   !ssg(k-1,g(k-1))+1 
                    endif
!     si on suppose vraiment une structure NESTED
                    goto 101
                endif

                if(auxssng.ne.ssgroupe.and.auxng.eq.groupe)then 
!     chgt de ssgroupe mais pas de groupe
                    nssgexact=nssgexact+1
                    auxssng=ssgroupe
                    if(k.ne.1)then
                        g(k)=g(k-1)
                        ssg(k,g(k))=ssg(k-1,g(k-1))+1
                    endif
!      write(*,*)'** chgt ssgp',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
                    goto 101
                endif

                if(k.ne.1)then
                    if(auxssng.eq.ssgroupe.and.auxng.eq.groupe)then 
!     aucun chgt de ssgroupe ni de gpe
                        if(k.ne.1)then
                            g(k)=g(k-1)
                            ssg(k,g(k))=ssg(k-1,g(k-1))
                        endif
!      write(*,*)'** pas chgt',k,ssg(k,g(k)),nssgexact,g(k),ngexact,ic
                        goto 101
                    endif
                endif
101     continue
                 nig(g(k)) = nig(g(k))+1
            endif
        endif

        if (maxtt.lt.t1(k))then
            maxtt = t1(k)
        endif

        if (mint.gt.t0(k)) then
            mint = t0(k)
        endif

    end do 
!    print*,ngexact
!AD:
    if (typeof .ne. 0) then 
        cens = maxtt
    end if
!Ad
    nsujet = k
    nssgbyg=nssgexact

!--------------------------- fin lecture du fichier


    if (typeof == 0) then

        nz1=nz
        if(nst.eq.2)then
            nz2=nz
        endif
        if(nz.gt.20)then
            nz = 20
        endif 
        if(nz.lt.4)then
            nz = 4
        endif
!***************************************************
!--------------- zi- ----------------------------------

!      construire vecteur zi (des noeuds)
    end if

    min = 1.d-10
    max = maxtt

    do i=1,2*nsujet
        do k=1,nsujet
            if((t0(k) .ge. min))then
                if(t0(k) .lt. max)then
                max = t0(k)
                endif
            endif
            if((t1(k) .ge. min))then
                if(t1(k) .lt. max)then
                max = t1(k)
                endif
            endif
        end do
        aux(i) = max
        min = max + 1.d-12
        max = maxtt
    end do

    date(1) = aux(1)

    k = 1
    do i=2,2*ns0
        if(aux(i).gt.aux(i-1))then
            k = k+1
            date(k) = aux(i)
        endif
    end do

    if(typeof == 0) then

! Al:10/03/2014 emplacement des noeuds splines en percentile
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
            zi(-1) = mint
            zi(0) = mint
            zi(1) = mint
            j=0
            do j=1,nz-2
                pord = dble(j)/(dble(nz)-1.d0)
                call percentile3(t2,nbrecu,pord,zi(j+1))
            end do
            zi(nz) = maxtt
            zi(nz+1) = maxtt
            zi(nz+2) = maxtt
            zi(nz+3) = maxtt
            ziOut = zi
            deallocate(t2)
        else ! equidistant
            nzmax=nz+3
            allocate(zi(-2:nzmax))
            ndate = k

            zi(-2) = date(1)
            zi(-1) = date(1)
            zi(0) = date(1)
            zi(1) = date(1)
            h = (date(ndate)-date(1))/(nz-1)
            do i=2,nz-1
                zi(i) =zi(i-1) + h
            end do
            zi(nz) = date(ndate)
            zi(nz+1)=zi(nz)
            zi(nz+2)=zi(nz)
            zi(nz+3)=zi(nz)
            ziOut = zi
        endif

    endif

!---------- affectation nt0,nt1----------------------------

    indictronq=0
    do i=1,ns0
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

    if (typeof == 0) then

!---------- affectation des vecteurs de splines -----------------
        n = nz+2

        call vecspliN(n,ndate)

        allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax), &
        m3m2(nzmax),m3m1(nzmax),m3m(nzmax),m2m1(nzmax), &
        m2m(nzmax),m1m(nzmax))

        call vecpenN(n)
    end if

    nbpara =np
    npmax=np

!===== Allocation
    allocate(I1_hess(np,np),H1_hess(np,np),I2_hess(np,np) &
    ,H2_hess(np,np),HI1(np,np),HI2(np,np),HIH(np,np),IH(np,np),HI(np,np) &
    ,y(np,np),Hspl_hess(np,np)) 


    b=3.d-1

!***********************************************************
!************** NEW : cross validation  ***********************
!***********************************************************

    nvacross=nva
    nva=0
    effetcross=effet
    effet=0
    nstcross=nst
    nst=1
    do l=1,ns0
        stracross(l)=stra(l)
    end do
    do l=1,nsujet
        stra(l)=1
    end do
!*************************************************
    if(xmin1.le.0.d0)then
        xmin1 = 0.d0
    endif


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

        allocate(ttt(0:nbintervR))

        ttt(0)=0.d0

        ttt(nbintervR)=cens

        j=0
        do j=1,nbintervR-1
            if (equidistant.eq.0) then
                ttt(j)=(t2(ent*j)+t2(ent*j+1))/(2.d0)
            else
                ttt(j)=(cens/nbintervR)*j
            endif
        end do
        time = ttt
        deallocate(t2)
!------- FIN RECHERCHE DES NOEUDS
    end if


    if(typeof == 0) then
        allocate(hess(n,n),I_hess(n,n),H_hess(n,n),v(n*(n+3)/2))
        if(irep1.eq.1)then   !pas recherche du parametre de lissage
            xmin1 = dsqrt(xmin1)
            auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
            if (ni.ge.250) then
                do i=1,nz+2
                    b(i)=1.d-1 !3.d-1
                end do
                xmin1 = sqrt(10.d0)*xmin1
                auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
                if (ni.lt.250) then

                else
                    do i=1,nz+2
                        b(i)=1.d-1
                    end do
                    xmin1 = sqrt(10.d0)*xmin1
                    auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
                    if (ni.lt.250) then

                    endif
                endif
            endif
        else                   !recherche du parametre de lissage

            if(xmin1.le.0.d0)then
                xmin1 = 1.d0
            endif
!          write(*,*)' '
!          write(*,*)'                Searching smoothing parameter'
!          write(*,*)' '

            xmin1 = dsqrt(xmin1)
            auxi = estimvN(xmin1,n,b,y,ddl,ni,res)

            if(ddl.gt.-2.5d0)then
                xmin1 = dsqrt(xmin1)
                auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
                if(ddl.gt.-2.5d0)then
                    xmin1 = dsqrt(xmin1)
                    auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
                    if(ddl.gt.-2.5d0)then
                        xmin1 = dsqrt(xmin1)
                        auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
                        if(ddl.gt.-2.5d0)then
                            xmin1 = dsqrt(xmin1)
                            auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
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
                auxi = estimvN(xmin1,n,b,y,ddl,ni,res)
                if (ni.ge.250) then
                    do i=1,nz+2
                        b(i)=1.d-1
                    end do
                    xmin1 = sqrt(10.d0)*xmin1
                endif
            endif

            ax = xmin1
            bx = xmin1*dsqrt(1.5d0)
            call mnbrakN(ax,bx,cx,fa,fb,fc,b,n)
            tol = 0.001d0
            res = goldenN(ax,bx,cx,tol,xmin1,n,b,y,ddl)
            effet=0
            auxkappa(1)=xmin1*xmin1
            auxkappa(2)=0.d0
            call marq98j(auxkappa,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpansplines)

            if (istop .ne. 1) then
                istopp(1)=1
                goto 1000
            end if
        endif

    end if

    nva=nvacross ! pour la recherche des parametres de regression
    nst=nstcross ! avec stratification si necessaire
    effet=effetcross ! avec effet initial

    do l=1,ns0
        stra(l)=stracross(l) !retablissement stratification
    end do

!     if (typeof .ne. 0) then
!         allocate(vvv((npmax*(npmax+1)/2)))
!     end if

    if (typeof == 0) then 
        k0(1) = xmin1*xmin1
        k0(2) = 0.d0
        if(nst.eq.2)then
            k0(2) = xmin2
        endif
    end if

!=============================- fin cross validation
!===== initialisation des parametres de regression/pas effets aleatoires

!    write(*,*),'====================================='
!    write(*,*),'== avec var explicatives !t ============='
!    write(*,*),'====================================='

    ca=0.d0
    cb=0.d0
    dd=0.d0

    select case(typeof)
        case(1)
            b(1:nst*nbintervR) = 0.08d0
        case(2)
            b(1:nst*2) = 0.8d0
    end select 

    if (typeof .ne. 0) then
        !allocate(kkapa(2))
        allocate(betacoef(nst*nbintervR))
    end if

    effet=0

    select case(typeof)
        case(0)
            deallocate(hess,I_hess,H_hess,v)
            np = nst*n + nva
            allocate(hess(np,np),I_hess(np,np),H_hess(np,np),v(np*(np+3)/2))
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpansplines) 
        case(1)
            np = nst*nbintervR + nva
            allocate(hess(np,np),I_hess(np,np),H_hess(np,np),v(np*(np+3)/2),vvv(np*(np+1)/2))
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpancpm) 
        case(2)
            np = nst*2 + nva
            allocate(hess(np,np),I_hess(np,np),H_hess(np,np),v(np*(np+3)/2),vvv(np*(np+1)/2))
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpanweib) 

    end select

!    write(*,*)'istop ==========>',istop,' np ',np
    if (istop .ne. 1) then
        istopp(2)=1
        goto 1000
    end if

!=================================   debut effet = 1 groupe =======================================
!    write(*,*),'================================================'
!    write(*,*),'== avec var explicatives + effet groupe ========='
!    write(*,*),'================================================'
    do i=1,nva
        b(np-i+2)=b(np-i+1)
    end do
!    b(np-nva)=0.1d0
    b(np-nva)=0.5d0
    effet=1
    deallocate(hess,I_hess,H_hess,v)
    select case(typeof)
        case(0)
            np = nst*n + nva + effet
            allocate(hess(np,np),I_hess(np,np),H_hess(np,np),v(np*(np+3)/2))
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpansplines) 
        case(1)
            np = nst*nbintervR + nva + effet
            deallocate(vvv)
            allocate(hess(np,np),I_hess(np,np),H_hess(np,np),v(np*(np+3)/2),vvv(np*(np+1)/2))
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpancpm) 
        case(2)
            np = nst*2 + nva + effet
            deallocate(vvv)
            allocate(hess(np,np),I_hess(np,np),H_hess(np,np),v(np*(np+3)/2),vvv(np*(np+1)/2))
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpanweib) 
    end select

!    write(*,*)'istop ==========>',istop,' np ',np
    if (istop .ne. 1) then
        istopp(3)=1
        goto 1000
    end if  
    bgpe=b(np-nva)

!=================================   Fin effet = 1 groupe =======================================

!=================================   debut effet = 1 sub groupe =======================================
!    write(*,*),'================================================'
!    write(*,*),'== avec var explicatives + effet sub group ====='
!    write(*,*),'================================================'

    gaux=g !numero de groupe stokes
    gnew=0
    do i=1,ngexact
        do j=1,nsujet
            if(g(j).eq.i)then
                gnew(j)=ssg(j,i)
            endif
        end do
    end do
    g=gnew

    g=ssgroupe0

    n_ssgbygrp=0
    nssgmax=0
    do i=1,ng0
        do k=1,nsujet
            if(gaux(k) == i) then
                if(k==1)then
                    n_ssgbygrp(i) =  n_ssgbygrp(i)+1
                else
                    if(g(k) .ne. g(k-1))n_ssgbygrp(i) =  n_ssgbygrp(i)+1
                end if
            end if
        end do
        if(nssgmax <= n_ssgbygrp(i)) then
            nssgmax = n_ssgbygrp(i)
        end if
    end do

    indic_cumul=1
    allocate(mij(ng0,nssgmax),mij_ind(nssgbyg0),aux1(ngexact,nssgmax),aux2(ngexact,nssgmax))
    allocate(cumulhaz1(ngmax,nssgmax),cumulhaz0(ngmax,nssgmax))

    effet=1
    ngaux=ngexact
    ngexact=nssgexact
    b(np-nva)=0.5d0
    select case(typeof)
        case(0)
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpansplines) 
        case(1)
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpancpm) 
        case(2)
            call marq98j(k0,b(1:np),np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpanweib) 
    end select

!    write(*,*)'istop ==========>',istop,' np ',np 
    if (istop .ne. 1) then
        istopp(4)=1
        goto 1000
    end if

    bssgpe=b(np-nva)

!=================================   Fin effet = 1 sub groupe =======================================

        indic_cumul=0 
!=================================   debut  ensemble de parametres  =================================
!    write(*,*),'=========================================='
!    write(*,*),'== ensemble des parametres ============='!,k0(1)
!    write(*,*),'====================================='

    effet=2
    ngexact = ngaux
    g=gaux
    do i=1,nva
        b(np-i+2)=b(np-i+1)
    end do

    np=nbpara

    b(np-nva-1)=1.d-1 !bgpe! !initialisation alpha(groupe)
    b(np-nva)=1.d-1 !bssgpe!0.15d0   !initialisation eta(sous groupe)

!    write(*,*)'np ',np
!    write(*,*)'b ',b
    deallocate(hess,I_hess,H_hess,v)
    allocate(hess(np,np),I_hess(np,np),H_hess(np,np),v(np*(np+3)/2))
    select case(typeof)
        case(0)            
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpansplines) 
        case(1)
            deallocate(vvv)
            allocate(vvv(np*(np+1)/2))
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpancpm) 
        case(2)
            deallocate(vvv)
            allocate(vvv(np*(np+1)/2))
            call marq98j(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpanweib)
    end select

!    write(*,*)'istop ==========>',istop,' np ',np
!Al:
    EPS(1) = ca
    EPS(2) = cb
    EPS(3) = dd
!Al:
    if (istop .ne. 1) then
        istopp(5)=1
        goto 1000
    end if

!=================================   fin  ensemble de parametres  =================================

    j=(np-nva)*(np-nva+1)/2

    trace=0
    trace1=0
    trace2=0

    select case(typeof)

        case(0)
! strate1 :
            do i=1,nz1+2
                do j=1,nz1+2
                    H1_hess(i,j)=H_hess(i,j)
                    I1_hess(i,j)=I_hess(i,j)
                end do
            end do

            call multiN(H1_hess,I1_hess,nz1+2,nz1+2,nz1+2,HI1)

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
                call multiN(H2_hess,I2_hess,nz2+2,nz2+2,nz2+2,HI2)
                do i =1,nz2+2
                    trace2=trace2+HI2(i,i)
                end do
            endif ! pour la deuxieme strate

        case(1)
! strate1 :
            do i=1,nbintervR
                do j=1,nbintervR
                    H1_hess(i,j)=H_hess(i,j)
                    I1_hess(i,j)=I_hess(i,j)
                end do
            end do

            call multiN(H1_hess,I1_hess,nbintervR,nbintervR,nbintervR,HI1)

            do i =1,nbintervR
                trace1=trace1+HI1(i,i)
            end do
! strate2 :
            if (nst.eq.2) then
                do i=1,nbintervR
                    k=nbintervR+i
                    do j=1,nbintervR
                        l=nbintervR+j
                        H2_hess(i,j)=H_hess(k,l)
                        I2_hess(i,j)=I_hess(k,l)
                    end do
                end do
                call multiN(H2_hess,I2_hess,nbintervR,nbintervR,nbintervR,HI2)
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

            call multiN(H1_hess,I1_hess,2,2,2,HI1)

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
                call multiN(H2_hess,I2_hess,2,2,2,HI2)
                do i =1,2
                    trace2=trace2+HI2(i,i)
                end do
            endif ! pour la deuxieme strate
    end select

    call multiN(I_hess,H_hess,np,np,np,IH)
    call multiN(H_hess,IH,np,np,np,HIH)

    if(effet.eq.2.and.ier.eq.-1)then
        v((np-nva)*(np-nva+1)/2)=10.d10
    endif

    resOut=res

    do ss=1,npmax
        do sss=1,npmax
            HIHOut(ss,sss) = HIH(ss,sss)
            H_hessOut(ss,sss)= H_hess(ss,sss)
        end do
    end do
! --------------  Lambda and survival and cumulative hazard estimates

    select case(typeof)
        case(0)
            call distancensplines(nz1,nz2,b,effet,mt,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
        case(1)
            call distancecpm(b,nbintervR*nst,mt,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
        case(2)
            if (nst == 1) then
                typeof2 = 1
            else
                typeof2 = 2
            end if
            call distanceweib(b,np,mt,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
    end select

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
        call multiN(H_hess,I_hess,np,np,np,HI)
        do i =1,np
            LCV(1) = LCV(1) + HI(i,i)
        end do 
        LCV(1) = (LCV(1) - resnonpen) / nsujet
    else
!        write(*,*)'=========> Akaike information Criterion <========='
!        LCV(2) = 2.d0 * np - 2.d0 * resOut
        LCV(2) = (1.d0 / nsujet) *(np - resOut)
!        write(*,*)'======== AIC :',LCV(2)
    end if
!AD:end



1000    continue


    linearpred=0.d0
    Resmartingale=0.d0
    frailtypred=0.d0
    frailtypredg=0.d0
    frailtyvar=0.d0
    frailtyvarg=0.d0
    frailtysd=0.d0
    frailtysdg=0.d0

!---------------------Calcul residus de martingal
!    write(*,*),'=========================================='
!    write(*,*),'=============== Martingale =============',nva
!    write(*,*),'====================================='

    if(nva .gt. 0) then
        do i=1,nva !2
            coefBeta(1,i)=b(np-nva-effet+i)
        end do

        Xbeta = matmul(coefBeta,transpose(ve))

        Call ResidusMartingalen(funcpanres,Resmartingale,frailtypred,maxngg,frailtypredg,&
        frailtyvar,frailtyvarg,frailtysd,frailtysdg)

        do i=1,nsujet
            linearpred(i)=Xbeta(1,i)+dlog(frailtypred(g(i))*frailtypredg(g(i),ssg(i,g(i))))
        end do
    end if

    if((istopp(2)== 0).and.(istopp(3)== 0)) then
        deallocate(mij,mij_ind,aux1,aux2,cumulhaz1,cumulhaz0)
    end if

    deallocate(H1_hess,I2_hess,H2_hess,HI1,HI2,HIH,IH,HI,y,Hspl_hess,hess,gaux,gnew,aux,H_hess,I_hess)
    deallocate(date,t0,t1,c,stra,stracross,g,ssg,nig,mid,ve,vax,filtre,v,I1_hess,n_ssgbygrp)

    if (typeof == 0) then
        deallocate(zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
        mm3,mm2,mm1,mm,im3,im2,im1,im,nt0,nt1)
    end if

    if (typeof == 1) then
        deallocate(ttt,vvv,betacoef)
    end if

    if (typeof == 2) then
        deallocate(vvv,betacoef)
    end if

!    if (typeof .ne. 0) then
!        deallocate(kkapa)
!    end if
!    write(*,*) '====Fin nested  === '

    end subroutine nested!FIN prog principal



!========================== VECSPLI =====================
    subroutine vecspliN(n,ndate)

    use tailles
    use comon,only:date,zi,mm3,mm2,mm1,mm,im3,im2,im1,im 
    Implicit none

    integer::n,ndate,i,j,k
    double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2

!----------  calcul de u(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)
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

    end subroutine vecspliN

!========================== VECPEN ==============================
    subroutine vecpenN(n) 

    use tailles
    use comon,only:zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m!date

    Implicit none

    integer::n,i
    double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2,a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x

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
        *zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1)   &  
        *zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3) &
        +2.d0*zi(i+2)*zi(i)))/(b1*c1))
        mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
        m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2) &
        ))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2)) &
        +((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x* &
        (zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2)) &
        +((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x* &
        (2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)) )
        m3m1(i) = 192.d0*(((x3-(0.5d0*x2*(4.d0*zi(i+1)+2.d0*zi(i-1) &
        ))+x*(2.d0*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*a1)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+2)+zi(i-1)+zi(i))) &
        +x*(zi(i+1)*zi(i-1)+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(b1*a3)) &
        +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+3)+2.d0*zi(i)))+x*(zi(i+1) &
        *zi(i+3)+2.d0*zi(i+1)*zi(i)))/(c1*a3)))
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
        *zi(i)+zi(i)*zi(i)))/(c2*a0)))
        m1m(i) = 192.d0*(((-x3+(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i-1) &
        +zi(i+1)))-x*(2.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(a1*a0)) &
        +((-x3+(0.5d0*x2*(4.d0*zi(i)+zi(i-1)+zi(i+2))) &
        -x*(zi(i-1)*zi(i)+zi(i)*zi(i)+zi(i+2)*zi(i)))/(b1*a0)) &
        +((-x3+(0.5d0*x2*(5.d0*zi(i)+zi(i+3)))-x*(zi(i+3)*zi(i) &
        +2.d0*zi(i)*zi(i)))/(c1*a0))) 
    end do

    end subroutine vecpenN



!==========================  COSP  ====================================

    subroutine cospN(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)
    use tailles 
    Implicit none

    integer :: j,k,n,i
    double precision::x,ht,ht2,h2,som,lam,su,binf,bsup,lbinf,lbsup &
    ,pm,htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3 &
    ,hh2,mm,im3,mm2,h,gl,hh
    double precision,dimension(-2:npmax):: the,zi
    double precision,dimension(npmax,npmax):: y

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

    call confN(x,j,n,y,pm,zi)

    binf = dexp(-gl - 1.96d0*pm)
    su  = dexp(-gl)
    bsup = dexp(-gl + 1.96d0*pm)

    call conf1N(x,j,n,y,pm,zi)
    lbinf = lam - 1.96d0*pm
    lbsup = lam + 1.96d0*pm

    return

    end subroutine cospN
!=====================  CONF1  =============================
    subroutine conf1N(x,ni,n,y,pm,zi)
    use tailles
    Implicit none

    integer :: ni,i,n,j
    double precision :: mmspN,x,pm,res
    double precision,dimension(-2:npmax)::zi
    double precision,dimension(npmax)::vecti,aux
    double precision,dimension(npmax,npmax)::y

    do i=1,n
        vecti(i) = mmspN(x,ni,i,zi)
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

    end subroutine conf1N


!=====================  CONF  =============================
    subroutine confN(x,ni,n,y,pm,zi)

    use tailles
    Implicit none

    integer::ni,i,n,j
    double precision::ispN,x,pm,res
    double precision,dimension(-2:npmax)::zi
    double precision,dimension(52)::vecti,aux
    double precision,dimension(npmax,npmax)::y

    do i=1,n
        vecti(i) = ispN(x,ni,i,zi)
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

    end subroutine confN

!==========================   ISP   ==================================

    double precision function ispN(x,ni,ns,zi)

    use tailles
    Implicit none

    integer::ni,ns
    double precision::val,mmspN,x
    double precision,dimension(-2:npmax)::zi

    if(x.eq.zi(ni))then
        if(ni.le.ns-3)then
            val = 0.d0
        else
            if(ni.le.ns-2)then
                val = ((zi(ni)-zi(ni-1))*mmspN(x,ni,ns,zi))*0.25d0
            else
                if (ni.eq.ns-1)then
                    val = ((zi(ni)-zi(ni-2))*mmspN(x,ni,ns,zi)+ &
                    (zi(ni+3)-zi(ni-1))*mmspN(x,ni,ns+1,zi))*0.25d0
                else
                    if(ni.eq.ns)then
                        val = ((zi(ni)-zi(ni-3))*mmspN(x,ni,ns,zi)+ &
                        (zi(ni+2)-zi(ni-2))*mmspN(x,ni,ns+1,zi) &
                        +(zi(ni+3)-zi(ni-1))*mmspN(x,ni,ns+2,zi))*0.25d0
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
                val = (x-zi(ni))*mmspN(x,ni,ns,zi)*0.25d0
            else  
                if(ni.eq.ns-2)then
                    val = ((x-zi(ni-1))*mmspN(x,ni,ns,zi)+ &
                    (zi(ni+4)-zi(ni))*mmspN(x,ni,ns+1,zi))*0.25d0
                else   
                    if (ni.eq.ns-1)then
                                   val =((x-zi(ni-2))*mmspN(x,ni,ns,zi)+ &
                        (zi(ni+3)-zi(ni-1))*mmspN(x,ni,ns+1,zi) &
                        +(zi(ni+4)-zi(ni))*mmspN(x,ni,ns+2,zi))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val =((x-zi(ni-3))*mmspN(x,ni,ns,zi)+ &
                            (zi(ni+2)-zi(ni-2))*mmspN(x,ni,ns+1,zi) &
                            +(zi(ni+3)-zi(ni-1))*mmspN(x,ni,ns+2,zi) &
                            +(zi(ni+4)-zi(ni))*mmspN(x,ni,ns+3,zi))*0.25d0
                        else
                            val = 1.d0
                        endif
                    endif
                endif
            endif
        endif
    endif
        ispN = val

    return

    end function ispN
!==========================  MMSP   ==================================
    double precision function mmspN(x,ni,ns,zi)
    use tailles
    Implicit none

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
                *(x-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3)&
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
                        +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni)))  &
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

    mmspN = val
    return

    end function mmspN

!================== multiplication de matrice  ==================


    subroutine multiN(A,B,IrowA,JcolA,JcolB,C)

    use tailles
    Implicit none

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

    end subroutine multiN

!=======cross validation

!========================          MNBRAK         ===================
    subroutine mnbrakN(ax,bx,cx,fa,fb,fc,b,n)

    use tailles
    Implicit none

    double precision::ax,bx,cx,fa,fb,fc,res,estimvN,gold &
    ,glimit,tiny,dum,fu,q,r,u,ulim,aux
    double precision,dimension(npmax)::b
    double precision,dimension(npmax,npmax)::y
    parameter (gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
    integer::n,ni


    fa = estimvN(ax,n,b,y,aux,ni,res)
    fb = estimvN(bx,n,b,y,aux,ni,res)

    if(fb.gt.fa)then
        dum = ax
        ax = bx
        bx = dum
        dum = fb
        fb = fa
        fa = dum
    endif

    cx = bx + gold*(bx-ax)
    fc = estimvN(cx,n,b,y,aux,ni,res)
 1       if(fb.ge.fc)then
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),tiny),q-r))
        ulim = bx + glimit*(cx-bx)
        if((bx-u)*(u-cx).gt.0.d0)then
            fu = estimvN(u,n,b,y,aux,ni,res)
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
            fu = estimvN(u,n,b,y,aux,ni,res)
        else
            if((cx-u)*(u-ulim).gt.0.d0)then
                fu = estimvN(u,n,b,y,aux,ni,res)
                if(fu.lt.fc)then
                    bx = cx
                    cx = u
                    u = cx + gold*(cx-bx)
                    fb = fc
                    fc = fu
                    fu = estimvN(u,n,b,y,aux,ni,res)
                endif
            else
                if((u-ulim)*(ulim-cx).ge.0.d0)then
                u = ulim
                fu = estimvN(u,n,b,y,aux,ni,res)
                else
                u = cx + gold*(cx-bx)
                fu = estimvN(u,n,b,y,aux,ni,res)
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

    end subroutine mnbrakN

!========================      GOLDEN   =========================
    double precision function goldenN(ax,bx,cx,tol,xmin,n,b,y,aux)

    use tailles
    Implicit none

    double precision,dimension(npmax,npmax)::y
    double precision,dimension(npmax)::b
    double precision::ax,bx,cx,tol,xmin,r,c,res,f1,f2,x0 &
    ,x1,x2,x3,estimvN,aux
    parameter (r=0.61803399d0,c=1.d0-r)
    integer n,ni

    x0 = ax
    x3 = cx
    if(abs(cx-bx).gt.abs(bx-ax))then
        x1 = bx
        x2 = bx + c*(cx-bx)
    else
        x2 = bx
        x1 = bx - c*(bx-ax)
    endif
    f1 = estimvN(x1,n,b,y,aux,ni,res)
    f2 = estimvN(x2,n,b,y,aux,ni,res)

 1    if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
        if(f2.lt.f1)then
            x0 = x1
            x1 = x2
            x2 = r*x1 + c*x3
            f1 = f2
            f2 = estimvN(x2,n,b,y,aux,ni,res)
        else
            x3 = x2
            x2 = x1
            x1 = r*x2+c*x0
            f2 = f1
            f1 = estimvN(x1,n,b,y,aux,ni,res)
        endif
        go to 1
    endif
    if(f1.lt.f2)then
        goldenN = f1
        xmin = x1
    else
        goldenN = f2
        xmin = x2
    endif

    return

    end function goldenN


!========================          ESTIMV         ===================

    double precision function estimvN(k00,n,b,y,aux,ni,res)

    use tailles
    use optim
    use comon,only:ndate,date,zi,pe,effet,mm3,mm2,mm1,mm,im3,im2,im1,im
    !nt0,nt1,nva,t0,t1,c,nsujet,nst,nz1,nz2

    Implicit none

    double precision::res,som,h1,k00,aux
    double precision,dimension(npmax*(npmax+3)/2)::v
    double precision,dimension(npmax,npmax)::y
    double precision,dimension(npmax)::bh,b
    double precision,dimension(ndatemax)::ut,dut
    double precision,dimension(2)::k0
    double precision,dimension(-2:npmax)::the
    integer::n,ij,i,k,j,vj,ier,istop,ni
    double precision::ca,cb,dd,funcpansplines
    external::funcpansplines
    ca=0.d0
    cb=0.d0
    dd=0.d0
    estimvN=0.d0
    j=0
    k0(1) = k00*k00
    k0(2)=0.d0

    call marq98j(k0,b,n,ni,v,res,ier,istop,effet,ca,cb,dd,funcpansplines)
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

        call testN(dut,k0,n,aux,y)
        estimvN = - ((res-pe)) - aux
    else
        aux = -n
    endif
!AD:
50    continue
!AD:
    return

    end function estimvN

!=================calcul de la hessienne  et de omega  ==============
    subroutine testN(dut,k0,n,res,y)

    use tailles
    !use comon,only:date,zi,t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst

    Implicit none

    double precision,dimension(npmax,npmax)::hessh,hess,omeg,y
    integer::n,i,j,np,indx(npmax)
    double precision,dimension(2)::k0
    double precision,dimension(ndatemax)::dut
    double precision::d,tra,res

    res=0.d0

    do i = 1,n
        do j = 1,n
            hess(i,j) = 0.d0 
        end do
    end do

    do i = 1,n
        do j = i,n
            call matN(hess(i,j),dut,i,j,n)
        end do
    end do
    do i = 2,n
        do j = 1,i-1
            hess(i,j)=hess(j,i)
        end do
    end do

    call calcomegN(n,omeg)

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
    call ludcmpN(hess,n,indx,d)

    do j=1,n
        call lubksbN(hess,n,indx,y(1,j))
    end do

    tra = 0.d0
    do i=1,n
        do j=1,n
            tra = tra + y(i,j)*hessh(j,i)
        end do
    end do

    res = (tra)

    end subroutine testN

!======================  LUBKSB  ======================================
    subroutine lubksbN(a,n,indx,b)

    use tailles

    Implicit none

    integer::n,i,ii,j,ll
    integer,dimension(npmax)::indx
    double precision,dimension(npmax)::b
    double precision,dimension(npmax,npmax)::a
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

    end subroutine lubksbN

!======================  LUDCMP  ======================================
    subroutine ludcmpN(a,n,indx,d)

    use tailles

    Implicit none

    integer::n,i,imax,j,k
    integer,dimension(n)::indx
    double precision,dimension(npmax,npmax)::a
    integer,parameter::nmax=500
    double precision,parameter::tiny=1.d-20
    double precision::aamax,dum,sum,d
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

    end subroutine ludcmpN

!=======================  CALOMEG  ===========================

    subroutine calcomegN(n,omeg)

    use tailles
    use comon,only:m3m!,m3m1,m3m3,mmm,m3m2,m2m1,m2m,m2m2,m1m,m1m1,date,zi

    Implicit none

    double precision,dimension(npmax,npmax)::omeg
    integer::n,i,j
    double precision:: calc00N,calc01N,calc02N


    do i=1,n
        do j=1,n
            omeg(i,j)=0.d0
        end do
    end do

    omeg(1,1)=calc00N(1,n)
    omeg(1,2)=calc01N(1,n)
    omeg(1,3)=calc02N(1,n)
    omeg(1,4)=m3m(1)
    omeg(2,1)=omeg(1,2)
    omeg(2,2)=calc00N(2,n)
    omeg(2,3)=calc01N(2,n)
    omeg(2,4)=calc02N(2,n)
    omeg(2,5)=m3m(2)
    omeg(3,1)=omeg(1,3)
    omeg(3,2)=omeg(2,3)
    omeg(3,3)=calc00N(3,n)
    omeg(3,4)=calc01N(3,n)
    omeg(3,5)=calc02N(3,n)
    omeg(3,6)=m3m(3)

    do i=4,n-3
        omeg(i,i-3)=omeg(i-3,i)
        omeg(i,i-2)=omeg(i-2,i)
        omeg(i,i-1)=omeg(i-1,i)
        omeg(i,i)=calc00N(i,n)
        omeg(i,i+1)=calc01N(i,n)
        omeg(i,i+2)=calc02N(i,n)
        omeg(i,i+3)=m3m(i)
    end do

    omeg(n-2,n-5)=omeg(n-5,n-2)
    omeg(n-2,n-4)=omeg(n-4,n-2)
    omeg(n-2,n-3)=omeg(n-3,n-2)
    omeg(n-2,n-2)=calc00N(n-2,n)
    omeg(n-2,n-1)=calc01N(n-2,n)
    omeg(n-2,n)=calc02N(n-2,n)
    omeg(n-1,n-4)=omeg(n-4,n-1)
    omeg(n-1,n-3)=omeg(n-3,n-1)
    omeg(n-1,n-2)=omeg(n-2,n-1)
    omeg(n-1,n-1)=calc00N(n-1,n)
    omeg(n-1,n)=calc01N(n-1,n)
    omeg(n,n-3)=omeg(n-3,n)
    omeg(n,n-2)=omeg(n-2,n)
    omeg(n,n-1)=omeg(n-1,n)
    omeg(n,n)=calc00N(n,n)

    end subroutine calcomegN


!====================  MAT  ==================================
    subroutine matN(res,dut,k,l,n)

    use tailles
    use comon,only:date,zi,c,nt1,nsujet!,nva,ndate,nstnt0,t0,t1

    Implicit none

    double precision::res,res1,mspN,aux2,u2
    double precision,dimension(ndatemax)::dut
    integer::k,l,j,ni,n,i

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
            aux2 = mspN(nt1(i),ni,k)*mspN(nt1(i),ni,l)
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

    end subroutine matN

!==========================  MSP   ==================================
    double precision function mspN(i,ni,ns)

    use tailles
    use comon,only:date,zi
    Implicit none

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
                *(date(i)-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3)&
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

    mspN = val

    return

    end function mspN


!=========================  CALC00  =========================
    double precision function calc00N(j,n)

    use tailles
    use comon,only:m3m3,m2m2,m1m1,mmm!,m3m2,m3m1,m2m1,m3m,m2m,m1m

    Implicit none

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

    calc00N = part

    return

    end function calc00N

!=========================  CALC01  =========================

    double precision function calc01N(j,n)

    use tailles
    use comon,only:m1m,m2m1,m3m2!,m2m,m3m,m1m1,m2m2,m3m3,m3m1,mmm

    Implicit none

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

    calc01N = part

    return

    end function calc01N

!=========================  CALC02  =========================

    double precision function calc02N(j,n)

    use tailles
    use comon,only:m3m1,m2m!,m1m,m3m,m2m1,m3m3,m2m2,m1m1,mmm,m3m2

    Implicit none

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

    calc02N = part

    return

    end function calc02N

!===================================================================

    double precision function gammlnN(xx)

    double precision::xx
!     Returns the value ln[gamma(xx)] for xx > 0.
    INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
!     Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
        24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
        -.5395239384953d-5,2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*dlog(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
    end do

    gammlnN = tmp + dlog(stp*ser/x)

    return

    END FUNCTION gammlnN

!==================================================================

!==================================================================

    SUBROUTINE qgauss1N(a,b,ss) ! sans troncature

    use tailles
    !use comon,only:auxig
    !use commun,only:mij,mid,ngexact,nssgexact
    Implicit none

    double precision :: a,b,ss
    double precision ::auxfunc1a,auxfunc1b
!      external :: func1N
!      Returns as ss the integral of the function func1 between a and b, by ten-point Gauss-Legendre integration: the function is evaluated exactly ten times at interior points in the range of integration.
! func1 est l integrant, ss le resultat de l integrale

    integer::j
    double precision::func1N
    double precision,dimension(5)::w,x
    double precision::dx,xm,xr  !The abscissas and weights.
    SAVE  w,x
    DATA w/.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
    DATA x/.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/

    xm=5.d-1*(b+a)
    xr=5.d-1*(b-a)
    ss=0.d0 ! Will be twice the average value of the function,since the ten
    do j=1,5  !weights (five numbers above each used twice) sum to 2.
        dx=xr*x(j)
        auxfunc1a=func1N(xm+dx)
        auxfunc1b=func1N(xm-dx)
        ss = ss+w(j)*( auxfunc1a+ auxfunc1b)
    end do
    ss=xr*ss            !  Scale the answer to the range of integration.
    return

    END SUBROUTINE qgauss1N

!================================================
!==================================================================

    SUBROUTINE gaulagN(ss,choix,nnodes) 

    use tailles
    use donnees
    !use comon,only:auxig
    Implicit none

    double precision::ss,auxfunca,func0N,func1N,func2N,func3N,func4N &
    ,func5N,func6N
    external::func0N,func1N,func2N,func3N,func4N,func5N,func6N
! gauss laguerre
! func1 est l integrant, ss le resultat de l integrale sur 0 ,  +infty
    integer::j,choix,nnodes
    double precision,dimension(nnodes):: xx,ww

    if(nnodes.eq.20) then
      xx(1:nnodes) = x(1:nnodes)
      ww(1:nnodes) = w(1:nnodes)
    else if (nnodes.eq.32) then
      xx(1:nnodes) = x1(1:nnodes)
      ww(1:nnodes) = w1(1:nnodes)
    end if

    ss=0.d0
! Will be twice the average value of the function,since the ten
! weights (five numbers above each used twice) sum to 2.
    do j=1,nnodes
        if (choix.eq.1) then 
            auxfunca=func1N(xx(j))
            ss = ss+ww(j)*(auxfunca)
        else                   !choix=2, troncature
            if (choix.eq.2) then 
                auxfunca=func2N(xx(j))
                ss = ss+ww(j)*(auxfunca)
            else                   !choix=3,essai, res = -1
                if (choix.eq.3) then
                    auxfunca=func3N(xx(j)) !dexp(-x(j))
                    ss = ss+ww(j)*(auxfunca)
                endif
            endif
        endif
    end do
    
    return

    END SUBROUTINE gaulagN

!================================================

    double precision function func1N(frail)
! calcul de l integrant, pour un effet aleatoire donne frail et un groupe donne auxig (cf funcpa)

    use tailles
    use comon,only:auxig,g,nsujet,alpha,eta
    !t0,t1,c,nt0,nt1,indictronq,stra,nva,nig,ndate,nst
    
    use commun,only:ngexact,aux1,ssg,mij,mid!nssgexact,aux2
    use residusM,only:n_ssgbygrp
    Implicit none


    double precision::frail
    integer::ig,issg,k
    double precision,dimension(ngexact)::prod1 !

    ig=auxig
!     initialisation de prod1 et prod2 par le numerateur de l integrant

    ! prod1(ig)=(dexp(-frail/alpha))*(frail**(1.d0/alpha-1.d0+mid(ig)))
    ! reecriture du numerateur pour eviter les bugs quand alpha trop petit
    prod1(ig) = dexp((1.d0/alpha-1.d0+mid(ig))*dlog(frail)-(frail/alpha))

!    write(*,*)'groupe',ig,'mid',mid(ig),'n_ssgbygrp(ig)',n_ssgbygrp(ig)
!    write(*,*)'eta',eta,'alpha',alpha,'frail',frail

    do issg=1,n_ssgbygrp(ig) !!! attention sous gpe pour un gpe donne
        do k=1,nsujet
            if((g(k).eq.ig).and.(ssg(k,g(k)).eq.issg))then
                prod1(ig)=prod1(ig) &
                *(1.d0+eta*frail*aux1(g(k),ssg(k,g(k)))) &
                **(-(1.d0/eta)-mij(g(k),ssg(k,g(k))))
            !    write(*,*)'groupe',ig,'eta',eta,'alpha',alpha,'frail',frail
!                write(*,*)'mij(g(k),ssg(k,g(k)))',mij(g(k),ssg(k,g(k)))
!                write(*,*)'(1.d0+eta*frail*aux1(g(k),ssg(k,g(k))))',(1.d0+eta*frail*aux1(g(k),ssg(k,g(k))))
!                write(*,*)'(-(1.d0/eta)-mij(g(k),ssg(k,g(k))))',(-(1.d0/eta)-mij(g(k),ssg(k,g(k))))
            !    write(*,*)'** prod1 **',prod1(ig),'**value',(1.d0+eta*frail*aux1(g(k),ssg(k,g(k)))) &
            !    **(-(1.d0/eta)-mij(g(k),ssg(k,g(k))))
                exit
            endif
        end do
    end do
!    write(*,*)'** prod1 **',prod1(ig),frail,eta,alpha,ig
    func1N = prod1(ig)

    return

    end function func1N
!==================================================================

    SUBROUTINE qgauss2N(a,b,ss) ! avec troncature
    !use comon,only:auxig
    implicit none

    double precision::a,b,ss,func2N
    double precision::auxfunc2a,auxfunc2b
!      Returns as ss the integral of the function func between a and b, by ten-point Gauss-Legendre integration: the function is evaluated exactly ten times at interior points in the range of integration.
! func2 est l integrant
    integer::j
    double precision,dimension(5)::w,x
    double precision::dx,xm,xr  !     The abscissas and weights.
    SAVE w,x
    DATA w/.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
    DATA x/.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
    xm=5.d-1*(b+a)
    xr=5.d-1*(b-a)
    ss=0.d0 ! Will be twice the average value of the function,since the ten
    do j=1,5  !weights (five numbers above each used twice) sum to 2.
        dx=xr*x(j)
        auxfunc2a=func2N(xm+dx)
        auxfunc2b=func2N(xm-dx)
        ss = ss+w(j)*( auxfunc2a+ auxfunc2b) 
    end do
    ss=xr*ss            !  Scale the answer to the range of integration.
    return

    END SUBROUTINE qgauss2N

!================================================

    double precision function func2N(frail) ! calcul de l integrant, pour le calcul d intregrale avec troncature
! calcul de l integrant, pour un effet aleatoire donne frail et un groupe donne auxig (cf funcpa)

    use tailles
    use comon,only:auxig,g &
    ,alpha,eta,indictronq, &
    nsujet!c,ndate,nig,nst,nt0,nt1,nva,t0,t1,stra
    use commun,only:ssg,ngexact,aux2!aux1,mid,mij,nssgexact
    use residusM,only:n_ssgbygrp

    Implicit none

    integer::ig,issg,k
    double precision::frail
    double precision,dimension(ngexact)::prod2 !ngmax

    ig=auxig

!    PROD2(IG)=(DEXP(DBLE(-frail/ALPHA))) &
!    *(DBLE(frail)**(1.d0/alpha-1.d0))

!    PROD2(IG)=(DEXP(-frail/ALPHA)) &
!    *(frail**(1.d0/alpha-1.d0))

    ! reecriture du numerateur pour eviter les bugs quand alpha trop petit
    prod2(ig) = dexp((1.d0/alpha-1.d0)*dlog(frail)-(frail/alpha))

    do issg=1,n_ssgbygrp(ig)
        do k=1,nsujet
        if((g(k).eq.ig).and.(ssg(k,g(k)).eq.issg))then
            if(indictronq.eq.1)then
            prod2(ig)=prod2(ig)*((1.d0 &
        +(eta*frail*aux2(g(k),ssg(k,g(k)))))**(-1.d0/eta)) 
    !    +(eta*dble(frail)*aux2(g(k),ssg(k,g(k)))))**(-1.d0/eta)) 
            endif
            exit
        endif
        end do
    end do

    func2N= prod2(ig)
    return

    end function func2N

!================================================
!================================================

    double precision function func3N(frail)

    use tailles
    use comon,only:auxig,g,alpha,eta,indictronq,nsujet
    !t0,t1,c,nt0,nt1,nva,ndate,nst,stra,nig
    
    use commun,only:ssg,ngexact,aux1,aux2,mij,mid!nssgexact
    use residusM,only:n_ssgbygrp

    Implicit none

    integer::ig,issg,k
    double precision::frail
    double precision,dimension(ngexact)::prod3 !ngmax

    ig = auxig
    ! prod3(ig) = (frail**(1.d0/alpha-1.d0+mid(ig)))*dexp(-frail/alpha)
    ! reecriture du numerateur pour eviter les bugs quand alpha trop petit
    prod3(ig) = dexp((1.d0/alpha-1.d0+mid(ig))*dlog(frail)-(frail/alpha))

    do issg=1,n_ssgbygrp(ig)
        do k=1,nsujet
        if ((g(k).eq.ig).and.(ssg(k,g(k)).eq.issg)) then
            if (indictronq.eq.1) then
                prod3(ig) = prod3(ig)*(eta*frail*(aux1(g(k),ssg(k,g(k)))-aux2(g(k),ssg(k,g(k))))+1.d0) &
                **(-(1.d0/eta)-mij(g(k),ssg(k,g(k))))
            endif
            exit
        endif
        end do
    end do

    func3N = prod3(ig)

    return

    end function func3N

!===============================================================================
!===================================== END =====================================
!===============================================================================


