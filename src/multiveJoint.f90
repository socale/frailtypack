    
!-------------------------------------
!             Joint multive
!-------------------------------------

!nobsEvent=c(nsujet0,nsujetmeta0,ng0)
!nbvar=c(nva10,nva20,nva30)
!noVarEvent=c(noVar1,noVar2,noVar3)
!maxIteration=c(maxit0,maxit00)
!nbIntervEvent=c(nbintervR0,nbintervDC0,nbintervM0)
!mtEvent=c(mt1,mt2,mt3)
!cptEvent=c(cpt,cpt_dc,cptmeta)
!critCV=c(ier,istop,istopshared(1),istopshared(2),istopshared(3))
!mt1Event=c(mt11,mt12,mt13)
!ResMartingaleEvent=c(Res_martingale,Res_martingaledc,Res_martingale2)
!frailtyEstimates=c(frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr)

    subroutine joint_multiv(nobsEvent,nz0,k0,tt00,tt10,tt0meta0,tt1meta0,ic0,icmeta0, &
    groupe0,groupe0meta,groupe0dc,tt0dc0,tt1dc0,icdc0,nbvar,vax0,vaxmeta0,vaxdc0,noVarEvent, &
    maxIteration,initialize,np,b,H_hessOut,HIHOut,resOut,LCV,critCV,x1Out,lamOut,xSu1,suOut,x2Out, & 
    lam2Out,xSu2,su2Out,x3Out,lam3Out,xSu3,su3Out,typeof0,equidistant0,nbIntervEvent, &
    mtEvent,ni,cptEvent,shape_weib,scale_weib,mt1Event,irep,ag0,ResMartingaleEvent,frailtyEstimates, &
    linearpred,linearpreddc,linearpredM,ziOut1,ziOutdc,ziOutmeta,time,timedc,timeM)

    use parametersmultiv
    use residusMmultiv
    use comonmultiv
    use taillesmultiv
    use splines
    use optim
    use sortiemultive

    IMPLICIT NONE  
    
    integer,dimension(3),intent(in)::nobsEvent,nbvar,noVarEvent,nbIntervEvent,mtEvent,mt1Event
    integer,dimension(2),intent(in)::maxIteration
    integer,dimension(3),intent(out)::cptEvent
    integer,dimension(5),intent(out)::critCV
    integer::maxit0,maxit00,mt11,mt12,mt13,nsujetmeta0,nva30,ag0,initialize
    integer,dimension(3),intent(in)::nz0
    integer::nsujet0,ng0,nva10,nva20,mt1,mt2,mt3,irep
    integer::np,equidistant,equidistant0
    integer,dimension(nobsEvent(1))::groupe0,ic0
    integer,dimension(nobsEvent(2))::groupe0meta,icmeta0 
    integer,dimension(nobsEvent(3)),intent(in)::icdc0,groupe0dc
    double precision,dimension(nobsEvent(3))::tt0dc0,tt1dc0
    double precision,dimension(nobsEvent(1))::tt00,tt10
    double precision,dimension(nobsEvent(2))::tt0meta0,tt1meta0
    double precision,dimension(3)::k0,kappaCV
    double precision,dimension(nobsEvent(1),nbvar(1)),intent(in):: vax0
    double precision,dimension(nobsEvent(3),nbvar(2)),intent(in):: vaxdc0
    double precision,dimension(nobsEvent(2),nbvar(3)),intent(in):: vaxmeta0
    double precision,dimension(np,np)::H_hessOut,HIHOut
    double precision::resOut
    double precision,dimension(mtEvent(1))::x1Out
    double precision,dimension(mtEvent(2))::x2Out
    double precision,dimension(mtEvent(3))::x3Out
    double precision,dimension(mtEvent(1),3)::lamOut
    double precision,dimension(mt1Event(1),3)::suOut
    double precision,dimension(mtEvent(2),3)::lam2Out
    double precision,dimension(mt1Event(2),3)::su2Out
    double precision,dimension(mtEvent(3),3)::lam3Out
    double precision,dimension(mt1Event(3),3)::su3Out
    integer::ss,sss
    double precision,dimension(np):: b
    
    double precision,dimension(:),allocatable::b01,b02,b03
    integer::np1,np2,np3        
    double precision,dimension(2)::LCV
    double precision,dimension(3)::shape_weib,scale_weib    
    integer::noVar1,noVar2,noVar3
    integer::cpt,cptmeta,cpt_dc,ier,ni
    integer::groupe,groupemeta,ij,kk,j,k,n,ii,iii,iii2,cptstr1,cptstr2   &
    ,i,ic,icmeta,icdc,istop,cptni,cptni1,cptni2,nb_echec,nb_echecor,id,cptbiais &
    ,cptauxdc   
    double precision::tt0,tt0meta,tt0dc,tt1,tt1meta,tt1dc,h,res,min,mindc,max,pord, &
    maxdc,maxt,maxtmeta,maxtdc,moy_peh0,moy_peh1,lrs,BIAIS_moy,mint,mintdc,mintmeta
    double precision,dimension(2)::res01
!AD: add for new marq
    double precision::ca,cb,dd
    double precision,external::funcpaMultivSplines,funcpaMultivCpm,funcpaMultivWeib
    double precision,dimension(mt1Event(1))::xSu1
    double precision,dimension(mt1Event(2))::xSu2
    double precision,dimension(mt1Event(3))::xSu3
!cpm
    integer::indd,ent,entdc,entmeta,typeof0,nzsha
    double precision::temp,kappa1    
    integer::nbintervR0,nbintervDC0,nbintervM0    

!predictor
    double precision,dimension(nobsEvent(3))::Res_martingale,Res_martingaledc,Res_martingale2,&
    frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr
    double precision,dimension(nobsEvent(3),3)::ResMartingaleEvent
    double precision,dimension(nobsEvent(3),5)::frailtyEstimates
    double precision,external::funcpamultires
    double precision,dimension(nobsEvent(1))::linearpred
    double precision,dimension(nobsEvent(3))::linearpreddc
    double precision,dimension(nobsEvent(2))::linearpredM    
    double precision,dimension(1,nbvar(1))::coefBeta
    double precision,dimension(1,nbvar(2))::coefBetadc    
    double precision,dimension(1,nbvar(3))::coefBetaM        
    double precision,dimension(1,nobsEvent(1))::XBeta
    double precision,dimension(1,nobsEvent(3))::XBetadc
    double precision,dimension(1,nobsEvent(2))::XBetaM
    double precision,dimension(:,:),allocatable::ve1,ve2,ve3
    double precision,dimension(nbIntervEvent(1)+1)::time
    double precision,dimension(nbIntervEvent(2)+1)::timedc
    double precision,dimension(nbIntervEvent(3)+1)::timeM
    double precision,dimension(nz0(1)+6)::ziOut1
    double precision,dimension(nz0(2)+6)::ziOutdc
    double precision,dimension(nz0(3)+6)::ziOutmeta
!-------- Parametres shared
    double precision::ddls
    double precision,dimension(:),allocatable::str00
    integer::cpts
    double precision,dimension(:,:),allocatable::H_hess0,HIH0
    double precision,dimension(2)::LCVs,shapeweibs,scaleweibs,k0s
!     double precision,dimension(:),allocatable::x1Outs,x2Outs
    double precision,dimension(:,:),allocatable::xTOuts !en plus
!     double precision,dimension(:,:),allocatable::lamOuts,lam2Outs
    double precision,dimension(:,:,:),allocatable::lamTOuts !en plus
!     double precision,dimension(100,3)::suOuts,su2Outs
    double precision,dimension(100,3,1)::suTOuts !en plus
!     double precision,dimension(100)::xSu1s,xSu2s
    double precision,dimension(100,1)::xSuTs !en plus
    double precision,dimension(:),allocatable::zis
    double precision,dimension(nobsEvent(3))::Resmartingales,frailtypreds,frailtysds,frailtyvars
    double precision,dimension(:),allocatable::linearpreds,martingaleCoxs,times
    integer::timedepMul,nbinnerknots0Mul,qorder0Mul
    double precision,dimension(:),allocatable::filtretps0Mul
    double precision,dimension(0:100,1)::BetaTpsMatMul
    double precision,dimension(3)::EPS


!     do i=1,5
!         write(*,*)groupe0(i),groupe0meta(i),groupe0dc(i)
!     enddo
    
!     do i=1,5
!         write(*,*)(vax0(i,j),j=1,nbvar(1))
!     enddo
!     write(*,*)'-----------------------------'    
!     do i=1,5
!         write(*,*)(vaxmeta0(i,j),j=1,nbvar(3))
!     enddo
!     write(*,*)'-----------------------------'
!     do i=1,5
!         write(*,*)(vaxdc0(i,j),j=1,nbvar(2))
!     enddo        

!     do i=1,5
!         write(*,*)ic0(i),groupe0(i),icmeta0(i),groupe0meta(i)
!     enddo        
    
!      write(*,11)'dataR',(tt00(i),tt10(i),ic0(i),groupe0(i),(vax0(i,j),j=1,nbvar(1)),i=1,20)
!     format(11)
! !     write(*,*)'dataDC',(tt00(i),tt10(i),ic0(i),groupe0(i),(vax0(i,j),j=1,nbvar(1)),i=1,20)    
!      write(*,*)'dataM',(tt0meta0(i),tt1meta0(i),icmeta0(i),groupe0meta(i),(vaxmeta0(i,j),j=1,nbvar(3)),i=1,20)    
!      stop
        
    nsujet0=nobsEvent(1)
    nsujetmeta0=nobsEvent(2)
    ng0=nobsEvent(3)
    nva10=nbvar(1)
    nva20=nbvar(2)
    nva30=nbvar(3)
    noVar1=noVarEvent(1)
    noVar2=noVarEvent(2)
    noVar3=noVarEvent(3)
    maxit0=maxIteration(1)
    maxit00=maxIteration(2)

    nbintervR0=nbIntervEvent(1)
    nbintervDC0=nbIntervEvent(2)
    nbintervM0=nbIntervEvent(3)
    mt1=mtEvent(1)
    mt2=mtEvent(2)
    mt3=mtEvent(3)
    mt11=mt1Event(1)
    mt12=mt1Event(2)
    mt13=mt1Event(3)
    kappaCV = 0.d0


!-------- end Parametres shared

    equidistant=equidistant0
    maxiter = maxit0 !maxiter joint
    typeof = typeof0


!!!!!   INITIALISATION DES B A PARTIR DES SHARED 

    select case(typeof)
        case(0)
            np1 = nz0(1) + 2 + nva10 + 1 
            np2 = nz0(2) + 2 + nva20 + 1 ! 0 Cox pour le dc
            np3 = nz0(3) + 2 + nva30 + 1
        case(1)
            k0=0.d0
            np1 = nbintervR0 + nva10 + 1 
            np2 = nbintervDC0 + nva20 + 1 ! 0 Cox
            np3 = nbintervM0 + nva30 + 1
        case(2)
            k0=0.d0
            np1 = 2 + nva10 +1 
            np2 = 2 + nva20 +1 ! 0 Cox
            np3 = 2 + nva30 +1
    end select
    
    allocate(b01(np1),b02(np2),b03(np3))
    
!================>
!================> initialisation parametres frailtypack
!================>
    if (initialize == 1) then
        b01=0.d0
        b02=0.d0
        b03=0.d0
        timedepMul=0
        nbinnerknots0Mul=0
        qorder0Mul=0
        BetaTpsMatMul=0.d0
!==================================> LOCO
!         write(*,*)''
!         write(*,*)'=============== estimation loco ================='
!         write(*,*)'========= typeof = (0:Splines, 1:Piecewise, 2:Weibull) ',typeof0    
!         write(*,*)'========= nombre de parametres loco ',np1
!         write(*,*)'========= nombre de sujet ',nsujet0
!         write(*,*)'========= nombre de groupe ',ng0,nva10,size(vax0),noVar1,ag0

        b01 = 0.25d0
        nzsha = nz0(1)
        kappa1 = k0(1)
        EPS(1) = 1.d-3
        EPS(2) = 1.d-3
        EPS(3) = 1.d-3
        allocate(str00(nsujet0))
        str00=1
        allocate(H_hess0(np1,np1),HIH0(np1,np1),zis(nzsha+6))
!         allocate(x1Outs(mt1),x2Outs(mt1),lamOuts(mt1,3),lam2Outs(mt1,3))
        allocate(xTOuts(mt1,1),lamTOuts(mt1,3,1))!en plus
        allocate(linearpreds(nsujet0),martingaleCoxs(nsujet0),times(nbintervR0+1))

        allocate(filtretps0Mul(nva10))
        filtretps0Mul = 0
        
        Call frailpenal(nsujet0,ng0,1,1,1, &
        nzsha,kappa1,tt00,tt10,ic0,groupe0,nva10,str00,vax0, &
        ag0,noVar1,maxit00,irep,np1,b01,H_hess0,HIH0,resOut,LCVs, &
        xTOuts,lamTOuts,xSuTs,suTOuts,typeof0,equidistant0,nbintervR0,mt1, &    
        ni,cpts,ier,k0s,ddls,istop,shapeweibs,scaleweibs,100,zis,Resmartingales,martingaleCoxs,&
        frailtypreds,frailtyvars,frailtysds,linearpreds,times,0,tt10,0, &
        timedepMul,nbinnerknots0Mul,qorder0Mul,filtretps0Mul,BetaTpsMatMul,EPS)


        deallocate(filtretps0Mul)
        
        deallocate(H_hess0,HIH0,zis,str00)
        
!         deallocate(x1Outs,x2Outs,lamOuts,lam2Outs)
        deallocate(xTOuts,lamTOuts)!en plus
        
        deallocate(linearpreds,martingaleCoxs,times)

        critCV(3) = istop
!        write(*,*)''
!        write(*,*)' Critere de convergence istop loco',istop
!        write(*,*)''
        if(typeof == 0) then
            if(irep.eq.0) then
                kappaCV(1)=k0s(1) ! kappa branche loco
            else
                kappaCV(1)=kappa1
            endif
 !            write(*,*)'Smoothing parameter loco ',kappaCV(1)
 !            write(*,*)'K0s ',k0s
 !            write(*,*)'irep',irep    
            
        endif

!==================================> META
!         write(*,*)''
!          write(*,*)'=============== estimation meta ================='
!         write(*,*)'========= typeof = (0:Splines, 1:Piecewise, 2:Weibull) ',typeof0    
!         write(*,*)'========= nombre de parametres ',np3
!         write(*,*)'========= nombre de sujet ',nsujetmeta0
!         write(*,*)'========= nombre de groupe ',ng0

        b03 = 0.25d0
        nzsha = nz0(3)
        kappa1 = k0(3)
        EPS(1) = 1.d-3
        EPS(2) = 1.d-3
        EPS(3) = 1.d-3
      allocate(str00(nsujetmeta0))
        str00=1
       allocate(H_hess0(np3,np3),HIH0(np3,np3),zis(nzsha+6))
    !         allocate(x1Outs(mt3),x2Outs(mt3),lamOuts(mt3,3),lam2Outs(mt3,3))
        allocate(xTOuts(mt1,1),lamTOuts(mt1,3,1))!en plus
       allocate(linearpreds(nsujetmeta0),martingaleCoxs(nsujetmeta0),times(nbintervM0+1))
        allocate(filtretps0Mul(nva30))
        filtretps0Mul = 0

        Call frailpenal(nsujetmeta0,ng0,1,1,1, &
        nzsha,kappa1,tt0meta0,tt1meta0,icmeta0,groupe0meta,nva30,str00,vaxmeta0, &
        ag0,noVar3,maxit00,irep,np3,b03,H_hess0,HIH0,resOut,LCVs, &
        xTOuts,lamTOuts,xSuTs,suTOuts,typeof0,equidistant0,nbintervM0,mt3, &    
        ni,cpts,ier,k0s,ddls,istop,shapeweibs,scaleweibs,100,zis,Resmartingales,martingaleCoxs,&
        frailtypreds,frailtyvars,frailtysds,linearpreds,times,0,tt10,0, &
        timedepMul,nbinnerknots0Mul,qorder0Mul,filtretps0Mul,BetaTpsMatMul,EPS)

        deallocate(filtretps0Mul)
        deallocate(H_hess0,HIH0,zis,str00)
!         deallocate(x1Outs,x2Outs,lamOuts,lam2Outs)
        deallocate(xTOuts,lamTOuts)!en plus
        deallocate(linearpreds,martingaleCoxs,times)

         critCV(5) = istop

!        write(*,*)''
!        write(*,*)' Critere de convergence istop meta',istop
!        write(*,*)''

        if(typeof == 0) then
            if(irep.eq.0) then
                kappaCV(3)=k0s(1) ! kappa branche meta
            else
                kappaCV(3)=kappa1
            end if
 !            write(*,*)'Smoothing parameter Meta ',kappaCV(3)
 !            write(*,*)'K0s ',k0s
 !            write(*,*)'irep',irep
        end if

! ==================================> DC
!        write(*,*)''
!        write(*,*)'=============== estimation dc ================='
!        write(*,*)'========= typeof = (0:Splines, 1:Piecewise, 2:Weibull) ',typeof0    
!        write(*,*)'========= nombre de parametres ',np2
!        write(*,*)'========= nombre de sujet ',ng0
!        write(*,*)'========= nombre de groupe ',ng0

        b02 = 0.25d0
        nzsha = nz0(2)
        kappa1 = k0(2)
        EPS(1) = 1.d-3
        EPS(2) = 1.d-3
        EPS(3) = 1.d-3
        allocate(str00(ng0))
        str00=1
        allocate(H_hess0(np2,np2),HIH0(np2,np2),zis(nzsha+6))
!         allocate(x1Outs(mt2),x2Outs(mt2),lamOuts(mt2,3),lam2Outs(mt2,3))
        allocate(xTOuts(mt1,1),lamTOuts(mt1,3,1))!en plus
        allocate(linearpreds(ng0),martingaleCoxs(ng0),times(nbintervDC0+1))
        allocate(filtretps0Mul(nva20))
        filtretps0Mul = 0


        Call frailpenal(ng0,ng0,1,1,1, & ! 0 Cox proportional hazards model
        nzsha,kappa1,tt0dc0,tt1dc0,icdc0,groupe0dc,nva20,str00,vaxdc0, &
        ag0,noVar2,maxit00,irep,np2,b02,H_hess0,HIH0,resOut,LCVs, &
        xTOuts,lamTOuts,xSuTs,suTOuts,typeof0,equidistant0,nbintervDC0,mt2, &    
        ni,cpts,ier,k0s,ddls,istop,shapeweibs,scaleweibs,mt12,zis,Resmartingales,martingaleCoxs,&
        frailtypreds,frailtyvars,frailtysds,linearpreds,times,0,tt10,0, &
        timedepMul,nbinnerknots0Mul,qorder0Mul,filtretps0Mul,BetaTpsMatMul,EPS)

        deallocate(filtretps0Mul)
        deallocate(H_hess0,HIH0,zis,str00)
!         deallocate(x1Outs,x2Outs,lamOuts,lam2Outs)
        deallocate(xTOuts,lamTOuts)!en plus
        deallocate(linearpreds,martingaleCoxs,times)

         critCV(4) = istop

!        write(*,*)''
!        write(*,*)' Critere de convergence istop dc',istop
!        write(*,*)''
        if(typeof == 0) then
            if(irep.eq.0) then
                kappaCV(2)=k0s(1) ! kappa branche dc
            else
                kappaCV(2)=kappa1
            endif
 !            write(*,*)'Smoothing parameter DC ',kappaCV(2)
 !            write(*,*)'K0s ',k0s
 !            write(*,*)'irep',irep
 !            write(*,*)'kappaCV',kappaCV
        end if
    else
        b01 = 0.25d0 ! recurrent
        b03 = 0.25d0 ! meta
        b02 = 0.25d0 ! loco
    end if

    if(typeof.ne.0) kappaCV=1.d0 

    !k0 = kappaCV

!!!!!  FIN  INITIALISATION DES B A PARTIR DES SHARED             

    Res_martingale=0.d0
    Res_martingaledc=0.d0
    Res_martingale2=0.d0
    frailtypred=0.d0
    frailtypred2=0.d0
    frailtyvar=0.d0
    frailtyvar2=0.d0
    frailtyCorr=0.d0
    linearpred=0.d0
    linearpreddc=0.d0
    linearpredM=0.d0

    if (typeof .eq. 0) then
        ziOut1=0.d0
        ziOutdc=0.d0
        ziOutmeta=0.d0
    endif

    model = 1
    
    indic_alpha=4
    indic_eta=3
    indic_rho=0
    indic_a1=1
    indic_a2=2
    
    if (typeof .ne. 0) then
        nbintervR = nbintervR0
        nbintervDC = nbintervDC0
        nbintervM = nbintervM0        
    end if

    allocate(vectn(3))
    
    epsa = 1.d-3
    epsb = 1.d-3
    epsd = 1.d-3

                
    lrs = 0.d0
    moy_peh0 = 0.d0
    moy_peh1 = 0.d0
    
    nb_echec = 0
    nb_echecor = 0
    nb0recu = 0
    moyrecu =0.d0             
        
    ngmax=ng0
    ng=ng0

    allocate(ResidusRec(ngmax),Residusdc(ngmax),ResidusRec2(ngmax),Rrec(ngmax),Nrec(ngmax),Rdc(ngmax),&
    Ndc(ngmax),Rrec2(ngmax),Nrec2(ngmax),vuu(2),nig(ngmax),nigmeta(ngmax),cdc(ngmax),t0dc(ngmax),t1dc(ngmax),&
    aux1(ngmax),aux2(ngmax),res1(ngmax),res4(ngmax),res3(ngmax),mi(ngmax),res1meta(ngmax),res3meta(ngmax))

    shape_weib = 0.d0
    scale_weib = 0.d0
    
    nsujetmax=nsujet0
    nsujet=nsujet0

    allocate(t0(nsujetmax),t1(nsujetmax),c(nsujetmax),stra(nsujetmax),g(nsujetmax),aux(2*nsujetmax))

    nsujetmeta=nsujetmeta0
    nsujetmetamax=nsujetmeta0
    allocate(t0meta(nsujetmeta),t1meta(nsujetmeta),cmeta(nsujetmeta),gmeta(nsujetmeta),auxmeta(2*nsujetmeta))
    ndatemeta=2*nsujetmeta 
    ndatemaxdc=2*ng0     

    if (typeof == 0) then
        allocate(nt0dc(ngmax),nt1dc(ngmax),nt0(nsujetmax),nt1(nsujetmax),mm3dc(ndatemaxdc),mm2dc(ndatemaxdc),&
        mm1dc(ndatemaxdc),mmdc(ndatemaxdc),im3dc(ndatemaxdc),im2dc(ndatemaxdc),im1dc(ndatemaxdc),imdc(ndatemaxdc),&
        nt0meta(nsujetmeta),nt1meta(nsujetmeta),mm3meta(ndatemeta),immeta(ndatemeta),mm2meta(ndatemeta),&
        mm1meta(ndatemeta),mmmeta(ndatemeta),im3meta(ndatemeta),im2meta(ndatemeta),im1meta(ndatemeta))        
    end if    

    nst=3    
    ni=0      
!---debut des iterations de simulations       
    id=1
    cptni=0
    cptni1=0
    cptni2=0
    biais_moy=0.d0
    cptbiais=0
    cptaux=0
    cptauxdc=0
    ij=0
    kk=0
    groupe=0
    groupemeta=0
    n=0
    nzloco=0
    nzdc=0
    nzmeta=0

    effet=1
    res01(1)=0.d0
    res01(2)=0.d0
!------------  entre non fichier et nombre sujet -----        
    nvarmax=ver

    allocate(vax(nva10),vaxdc(nva20),vaxmeta(nva30))
            
    nva1=nva10
    nva2=nva20
    nva3=nva30
    nva = nva1+nva2+nva3
    nvarmax=nva

    allocate(ve(nsujetmax,nvarmax),vedc(ngmax,nvarmax),vemeta(nsujetmeta,nvarmax),ve1(nsujetmax,nva1),&
    ve2(ngmax,nva2),ve3(nsujetmeta,nva3),filtre(nva10),filtre2(nva20),filtre3(nva30))
    
    nig=0
    nigmeta=0 
    
! AD: recurrent
    if (noVar1.eq.1) then 
!        write(*,*)'filtre recurrent desactive'
        filtre=0
        nva1=0
    else
!        write(*,*)'filtre recurrent active'
        filtre=1
    end if    
!AD:death
    if (noVar2.eq.1) then 
!        write(*,*)'filtre death desactive'
        filtre2=0
        nva2=0
    else
!        write(*,*)'filtre death active'
        filtre2=1
    end if    
    
    if ((noVar1.eq.1).or.(noVar2.eq.1)) then
        nva = nva1+nva2 
    end if
! AD: recurrent meta
    if (noVar3.eq.1) then 
!        write(*,*)'filtre meta desactive'
        filtre3=0
        nva3=0
    else
!        write(*,*)'filtre meta active'
        filtre3=1
    end if

!------------  lecture fichier -----------------------
    maxt = 0.d0
    mint = 0.d0

    maxtmeta = 0.d0
    mintmeta = 0.d0

    maxtdc = 0.d0
    mintdc = 0.d0

    cpt = 0
    cptmeta =0 
    cptcens = 0
    cptcensmeta = 0
    cpt_dc = 0
    k = 0
    cptstr1 = 0
    cptstr2 = 0

!ccccccccccccccccccccc
! pour le deces
!cccccccccccccccccccc
    do k = 1,ng 

        if (k.eq.1) then
            mintdc = tt0dc0(k) ! affectation du min juste une fois
        endif

        tt0dc=tt0dc0(k)
        tt1dc=tt1dc0(k)
        icdc=icdc0(k)
        groupe=groupe0dc(k)
        do j=1,nva20        
            vaxdc(j)=vaxdc0(k,j)
        enddo                
        if(tt0dc.gt.0.d0)then
            cptauxdc=cptauxdc+1
        endif                  
!------------------   deces c=1 pour donnees de survie
        if(icdc.eq.1)then
            cpt_dc = cpt_dc + 1
            cdc(k)=1
            t0dc(k) = tt0dc      !/100.d0
            t1dc(k) = tt1dc      !+0.000000001
            iii = 0
            iii2 = 0           
            do ii = 1,nva20
                if(filtre2(ii).eq.1)then
                    iii2 = iii2 + 1
                    vedc(k,iii2) = dble(vaxdc(ii))
                endif
            end do   
        else 
!------------------   censure a droite ou event recurr  c=0 
            if(icdc.eq.0)then
                cdc(k) = 0 
                iii = 0
                iii2 = 0
                do ii = 1,nva20
                    if(filtre2(ii).eq.1)then
                    iii2 = iii2 + 1
                    vedc(k,iii2) = dble(vaxdc(ii))
                    endif
                end do 
                t0dc(k) =  tt0dc 
                t1dc(k) = tt1dc    
            endif
        endif
        if (maxtdc.lt.t1dc(k))then
            maxtdc = t1dc(k)
        endif
        if (mintdc.gt.t0dc(k)) then
            mintdc = t0dc(k)
        endif
    end do

!AD:
    if (typeof .ne. 0) then 
        cens = maxtdc
    end if
!Ad    
    k = 0
    cptstr1 = 0
    cptstr2 = 0

!cccccccccccccccccccccccccccccccccc
! pour les donnees recurrentes  loco
!cccccccccccccccccccccccccccccccccc
    do i = 1,nsujet     !sur les observations

        if (i.eq.1) then
            mint = tt00(i) ! affectation du min juste une fois
        endif

        tt0=tt00(i)
        tt1=tt10(i)
        ic=ic0(i)
        groupe=groupe0(i)

        do j=1,nva10
            vax(j)=vax0(i,j)
        enddo

        if(tt0.gt.0.d0)then
            cptaux=cptaux+1
        endif
!-----------------------------------------------------
!    essai sans troncature
!    tt0=0.
!------------------   observation c=1 pour donnï¿½es recurrentes
        if(ic.eq.1)then
            cpt = cpt + 1
            c(i)=1
            t0(i) = tt0
            t1(i) = tt1
            t1(i) = t1(i)
            g(i) = groupe
            nig(groupe) = nig(groupe)+1 ! nb d event recurr dans un groupe
            iii = 0
            iii2 = 0
  
            do ii = 1,nva10
                if(filtre(ii).eq.1)then
                    iii = iii + 1
                    ve(i,iii) = dble(vax(ii)) !ici sur les observations
                endif
            end do
        else 
!------------------   censure a droite  c=0 pour donnees recurrentes
            if(ic.eq.0)then
                cptcens=cptcens+1
                c(i) = 0 
                iii = 0
                iii2 = 0

                do ii = 1,nva10
                    if(filtre(ii).eq.1)then
                    iii = iii + 1
                    ve(i,iii) = dble(vax(ii))
                    endif
                end do 
                t0(i) =  tt0
                t1(i) = tt1
                t1(i) = t1(i)
                g(i) = groupe
            endif
        endif
        if (maxt.lt.t1(i))then
            maxt = t1(i)
        endif
        if (mint.gt.t0(i)) then
            mint = t0(i)
        endif
    end do 

    nsujet=i-1
!cccccccccccccccccccccccccccccccccc
! pour les donnees recurrentes  meta
!cccccccccccccccccccccccccccccccccc
    do i = 1,nsujetmeta     !sur les observations

        if (i.eq.1) then
            mintmeta = tt0meta0(i) ! affectation du min juste une fois
        endif

        tt0meta=tt0meta0(i)
        tt1meta=tt1meta0(i)
        icmeta=icmeta0(i)
        groupemeta=groupe0meta(i)

        do j=1,nva30
            vaxmeta(j)=vaxmeta0(i,j)
        enddo

        if(tt0meta.gt.0.d0)then
            cptauxmeta=cptauxmeta+1
        endif
!-----------------------------------------------------
!     essai sans troncature
!     tt0=0.
!------------------   observation c=1 pour donnees recurrentes
        if(icmeta.eq.1)then
            cptmeta = cptmeta + 1
            cmeta(i)=1
            t0meta(i) = tt0meta 
            t1meta(i) = tt1meta  
            t1meta(i) = t1meta(i)
            gmeta(i) = groupemeta
            nigmeta(groupemeta) = nigmeta(groupemeta)+1 ! nb d event recurr dans un groupe
            iii = 0
            iii2 = 0

            do ii = 1,nva30
                if(filtre3(ii).eq.1)then
                    iii = iii + 1
                    vemeta(i,iii) = dble(vaxmeta(ii)) !ici sur les observations
                endif
            end do
        else 
!------------------   censure a droite  c=0 pour donnees recurrentes
            if(icmeta.eq.0)then
                cptcens=cptcens+1
                cmeta(i) = 0 
                iii = 0
                iii2 = 0
                do ii = 1,nva30
                    if(filtre3(ii).eq.1)then
                    iii = iii + 1
                    vemeta(i,iii) = dble(vaxmeta(ii))
                    endif
                end do 
                t0meta(i) =  tt0meta
                t1meta(i) = tt1meta
                t1meta(i) = t1meta(i)
                gmeta(i) = groupemeta
            endif
        endif
        if (maxtmeta.lt.t1meta(i))then
            maxtmeta = t1meta(i)
        endif
        if (mintmeta.gt.t0meta(i)) then
            mintmeta = t0meta(i)
        endif
    end do 

    if (typeof == 0) then    
        nzloco=nz0(1)
        nzdc=nz0(2)
        nzmeta=nz0(3)
        vectn(1)=nz0(1)+2
        vectn(2)=nz0(2)+2
        vectn(3)=nz0(3)+2
        
        nz1=nzloco
        nz2=nzdc
        nz3=nzmeta
        
        if(nzloco.gt.20)then
            nzloco = 20
        endif
        if(nzloco.lt.4)then
            nzloco = 4
        endif
        if(nzdc.gt.20)then
            nzdc = 20
        endif
        if(nzdc.lt.4)then
            nzdc = 4
        endif
        if(nzmeta.gt.20)then
            nzmeta = 20
        endif
        if(nzmeta.lt.4)then
            nzmeta = 4
        endif
    end if

    ndatemax=2*nsujet

    allocate(date(ndatemax),datedc(ndatemax),datemeta(ndatemeta))
    
    if(typeof == 0) then
        allocate(mm3(ndatemax),mm2(ndatemax) &
        ,mm1(ndatemax),mm(ndatemax),im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax))
    end if

!!!  DONNEES DECES

    mindc = 0.d0
    maxdc = maxtdc
    do i = 1,2*ng     
        do k = 1,ng
            if((t0dc(k).ge.mindc))then
                if(t0dc(k).lt.maxdc)then
                    maxdc = t0dc(k)
                endif
            endif
            if((t1dc(k).ge.mindc))then
                if(t1dc(k).lt.maxdc)then
                    maxdc = t1dc(k)
                endif
            endif
        end do   
        aux(i) = maxdc
        mindc = maxdc + 1.d-12
        maxdc = maxtdc
    end do
! 
    datedc(1) = aux(1)
    k = 1
    do i=2,2*ng
        if(aux(i).gt.aux(i-1))then
            k = k+1
            datedc(k) = aux(i)
        endif 
    end do 
    
    if(typeof == 0) then
        ndatedc = k   
    end if   

!--------------- zi- ----------------------------------
!      construire vecteur zi (des noeuds)

!    DONNEES RECURRENTES loco
    min = 0.d0
    aux =0.d0
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

    if(typeof==0) then
        ndate = k
    end if

!   DONNEES RECURRENTES meta

    min = 0.d0
    auxmeta =0.d0
    max = maxtmeta


    do i = 1,2*nsujetmeta
        do k = 1,nsujetmeta
            if((t0meta(k).ge.min))then
                if(t0meta(k).lt.max)then
                    max = t0meta(k)
                endif
            endif
            if((t1meta(k).ge.min))then
                if(t1meta(k).lt.max)then
                    max = t1meta(k)
                endif
            endif
        end do   
        auxmeta(i) = max
        min = max + 1.d-12
        max = maxtmeta
    end do

    datemeta(1) = auxmeta(1)
    k = 1
    do i=2,2*nsujetmeta
        if(auxmeta(i).gt.auxmeta(i-1))then
            k = k+1
            datemeta(k) = auxmeta(i)
        endif 
    end do 

    if(typeof==0) then
        ndatemeta = k
    end if

!==========================>
!==== Construction des Zi pour les splines
!==========================>
    if(typeof == 0) then
        if (equidistant.eq.0) then ! percentile

            ! recurrent
            i=0
            j=0
            do i=1,nsujet
                if(t1(i).ne.(0.d0).and.c(i).eq.1) then
                    j=j+1
                endif
            end do
            nbrecu=j

            allocate(t2(nbrecu))
            j=0
            do i=1,nsujet
                if (t1(i).ne.(0.d0).and.c(i).eq.1) then
                    j=j+1
                    t2(j)=t1(i)
                endif
            end do

            allocate(zi(-2:(nzloco+3)))
            ndate = k
            zi(-2) = mint
            zi(-1) = mint
            zi(0) = mint
            zi(1) = mint
            j=0
            do j=1,nzloco-2
                pord = dble(j)/(dble(nzloco)-1.d0)
                call percentile3(t2,nbrecu,pord,zi(j+1))
            end do
            zi(nzloco) = maxt
            zi(nzloco+1) = maxt
            zi(nzloco+2) = maxt
            zi(nzloco+3) = maxt
            ziOut1 = zi
            deallocate(t2)

            ! death
            i=0
            j=0
            do i=1,ng
                if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1) then
                    j=j+1
                endif
            end do
            nbdeces=j

            allocate(t3(nbdeces))
            j=0
            do i=1,ng
                if (t1dc(i).ne.(0.d0).and.cdc(i).eq.1) then
                    j=j+1
                    t3(j)=t1dc(i)
                endif
            end do

            allocate(zidc(-2:(nzdc+3)))
            zidc(-2) = mintdc
            zidc(-1) = mintdc
            zidc(0) = mintdc
            zidc(1) = mintdc
            j=0
            do j=1,nzdc-2
                pord = dble(j)/(dble(nzdc)-1.d0)
                call percentile3(t3,nbdeces,pord,zidc(j+1))
            end do
            zi(nzdc) = maxtdc
            zi(nzdc+1) = maxtdc
            zi(nzdc+2) = maxtdc
            zi(nzdc+3) = maxtdc
            ziOutdc = zidc
            deallocate(t3)

            ! meta
            i=0
            j=0
            do i=1,nsujetmeta
                if(t1meta(i).ne.(0.d0).and.cmeta(i).eq.1) then
                    j=j+1
                endif
            end do
            nbrecumeta=j

            allocate(t4(nbrecumeta))
            j=0
            do i=1,nsujetmeta
                if (t1meta(i).ne.(0.d0).and.cmeta(i).eq.1) then
                    j=j+1
                    t4(j)=t1meta(i)
                endif
            end do

            allocate(zimeta(-2:(nzmeta+3)))
            zimeta(-2) = mintmeta
            zimeta(-1) = mintmeta
            zimeta(0) = mintmeta
            zimeta(1) = mintmeta
            j=0
            do j=1,nzmeta-2
                pord = dble(j)/(dble(nzmeta)-1.d0)
                call percentile3(t4,nbrecumeta,pord,zimeta(j+1))
            end do
            zimeta(nzmeta) = maxtmeta
            zimeta(nzmeta+1) = maxtmeta
            zimeta(nzmeta+2) = maxtmeta
            zimeta(nzmeta+3) = maxtmeta
            ziOutmeta = zimeta
            deallocate(t4)

        else ! equidistant

            ! recurrent
            nzmax=nzloco+3
            allocate(zi(-2:nzmax))
            zi(-2) = date(1)
            zi(-1) = date(1)
            zi(0) = date(1)
            zi(1) = date(1)
            h = (date(ndate)-date(1))/dble(nzloco-1)
            do i=2,nzloco-1
                zi(i) =zi(i-1) + h
            end do
            zi(nzloco) = date(ndate)
            zi(nzloco+1)=zi(nzloco)
            zi(nzloco+2)=zi(nzloco)
            zi(nzloco+3)=zi(nzloco)
            ziOut1 = zi

            ! death
            allocate(zidc(-2:(nzdc+3)))
            zidc(-2) = datedc(1)
            zidc(-1)= datedc(1)
            zidc(0)= datedc(1)
            zidc(1)= datedc(1)
            h =(datedc(ndatedc)-datedc(1))/dble(nzdc-1)
            do i=2,nzdc-1
                zidc(i) =zidc(i-1)+h
            end do
            zidc(nzdc) = datedc(ndatedc)
            zidc(nzdc+1)=zidc(nzdc)
            zidc(nzdc+2)=zidc(nzdc)
            zidc(nzdc+3)=zidc(nzdc)
            ziOutdc = zidc

            ! meta
            allocate(zimeta(-2:(nzmeta+3)))
            zimeta(-2) = datemeta(1)
            zimeta(-1)= datemeta(1)
            zimeta(0)= datemeta(1)
            zimeta(1)= datemeta(1)
            h =(datemeta(ndatemeta)-datemeta(1))/dble(nzmeta-1)
            do i=2,nzmeta-1
                zimeta(i) =zimeta(i-1)+h
            end do
            zimeta(nzmeta) = datemeta(ndatemeta)
            zimeta(nzmeta+1)=zimeta(nzmeta)
            zimeta(nzmeta+2)=zimeta(nzmeta)
            zimeta(nzmeta+3)=zimeta(nzmeta)
            ziOutmeta = zimeta
         end if
    endif

!---------- affectation nt0dc,nt1dc DECES ----------------------------

    indictronqdc=0
    do k=1,ng
        if (typeof == 0) then
            if(nig(k).eq.0.d0)then
                nb0recu = nb0recu + 1 !donne nb sujet sans event recu
            endif
            moyrecu =  moyrecu + dble(nig(k))

            if(t0dc(k).eq.0.d0)then
                nt0dc(k) = 0
            endif
        end if

        if(t0dc(k).ne.0.d0)then
            indictronqdc=1
        endif

        if (typeof == 0) then
            do j=1,ndatedc
                if(datedc(j).eq.t0dc(k))then
                    nt0dc(k)=j
                endif
                if(datedc(j).eq.t1dc(k))then
                    nt1dc(k)=j
                endif
            end do
        end if
    end do

!---------- affectation nt0,nt1 RECURRENTS----------------------------

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
    indictronqmeta=0
    do i=1,nsujetmeta

        if (typeof == 0) then
            if(t0meta(i).eq.0.d0)then
                nt0meta(i) = 0
            endif
        end if
        if(t0meta(i).ne.0.d0)then
            indictronqmeta=1
        endif
        if (typeof == 0) then
            do j=1,ndatemeta
                if(datemeta(j).eq.t0meta(i))then
                    nt0meta(i)=j
                endif
                if(datemeta(j).eq.t1meta(i))then
                    nt1meta(i)=j
                endif
            end do
        end if
    end do

     if (typeof == 0) then
!---------- affectation des vecteurs de splines -----------------
        call vecspli(ndate,ndatedc,ndatemeta)
        allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax),m3m1(nzmax),&
        m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax),m3m3b(nzdc),m2m2b(nzdc),m1m1b(nzdc), &
        mmmb(nzdc),m3m2b(nzdc),m3m1b(nzdc),m3mb(nzdc),m2m1b(nzdc),m2mb(nzdc),m1mb(nzdc), &
        m3m3c(nzmeta),m2m2c(nzmeta),m1m1c(nzmeta),mmmc(nzmeta),m3m2c(nzmeta),m3m1c(nzmeta), &
        m3mc(nzmeta),m2m1c(nzmeta),m2mc(nzmeta),m1mc(nzmeta))    

        call vecpenP(nzloco+2,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
        call vecpenP(nzdc+2,zidc,m3m3b,m2m2b,m1m1b,mmmb,m3m2b,m3m1b,m3mb,m2m1b,m2mb,m1mb)
        call vecpenP(nzmeta+2,zimeta,m3m3c,m2m2c,m1m1c,mmmc,m3m2c,m3m1c,m3mc,m2m1c,m2mc,m1mc)
    end if

    npmax=np

    allocate(I_hess(npmax,npmax),H_hess(npmax,npmax),Hspl_hess(npmax,npmax) &
    ,PEN_deri(npmax,1),hess(npmax,npmax),v((npmax*(npmax+3)/2)),I1_hess(npmax,npmax) &
    ,H1_hess(npmax,npmax),I2_hess(npmax,npmax),H2_hess(npmax,npmax),HI2(npmax,npmax) & 
    ,HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax),BIAIS(npmax,1))

    if (typeof .ne. 0) then
        allocate(vvv((npmax*(npmax+1)/2)))
    end if

     if (typeof == 1) then
!============================================================================================>
!=======================================>   LOCO  <===========================================
!============================================================================================>
        j=0
        do i=1,nsujet
            if(t1(i).ne.(0.d0).and.c(i).eq.1) then
                j=j+1
            endif
        end do
        nbrecu=j
        allocate(t2(nbrecu))

        j=0
        do i=1,nsujet
            if (t1(i).ne.(0.d0).and.c(i).eq.1) then
                j=j+1
                t2(j)=t1(i)
            endif
        end do
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
!============================================================================================>
!=======================================>   DECES  <==========================================
!============================================================================================>
!
        j=0
        do i=1,ngmax
            if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
                j=j+1
            endif
        end do
        nbdeces=j

        allocate(t3(nbdeces))
        j=0
        do i=1,ngmax
            if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
                j=j+1
                t3(j)=t1dc(i)
            endif
        end do

        indd=1
        do while (indd.eq.1)
            indd=0
            do i=1,nbdeces-1
                if (t3(i).gt.t3(i+1)) then
                    temp=t3(i)
                    t3(i)=t3(i+1)
                    t3(i+1)=temp
                    indd=1
                end if
            end do
        end do

        entdc=int(nbdeces/nbintervDC)
        allocate(tttdc(0:nbintervDC))
        tttdc(0)=0.d0
        tttdc(nbintervDC)=cens

        j=0
        do j=1,nbintervDC-1
            if (equidistant.eq.0) then
                tttdc(j)=(t3(entdc*j)+t3(entdc*j+1))/(2.d0)
            else
                tttdc(j)=(cens/nbintervDC)*j
            endif
        end do     
        timedc = tttdc

!============================================================================================>
!=======================================>   META  <===========================================
!============================================================================================>


        j=0

        do i=1,nsujetmeta
            if(t1meta(i).ne.(0.d0).and.cmeta(i).eq.1) then
                j=j+1
            endif
        end do
        nbrecumeta=j
        allocate(t4(nbrecumeta))

        j=0
        do i=1,nsujetmeta
            if (t1meta(i).ne.(0.d0).and.cmeta(i).eq.1) then
                j=j+1
                t4(j)=t1meta(i)
            endif
        end do

        indd=1
        do while (indd.eq.1)
            indd=0
            do i=1,nbrecumeta-1
                if (t4(i).gt.t4(i+1)) then
                    temp=t4(i)
                    t4(i)=t4(i+1)
                    t4(i+1)=temp
                    indd=1
                end if
            end do
        end do

        entmeta=int(nbrecumeta/nbintervM)

        allocate(tttmeta(0:nbintervM))

        tttmeta(0)=0.d0
        tttmeta(nbintervM)=cens

        j=0
        do j=1,nbintervM-1
            if (equidistant.eq.0) then
                tttmeta(j)=(t4(entmeta*j)+t4(entmeta*j+1))/(2.d0)
            else
                tttmeta(j)=(cens/nbintervM)*j
            endif
        end do
        timeM = tttmeta
        deallocate(t2,t3,t4)
    end if

    ca=0.d0
    cb=0.d0
    dd=0.d0
    if (typeof .ne. 0) then
        allocate(kkapa(3))
    end if

!     write(*,*)''
!     write(*,*)'=============== Modele final ================='
!     write(*,*)'========= typeof = (0:Splines, 1:Piecewise, 2:Weibull) ',typeof    
!     write(*,*)'========= nombre de parametres ',np
!     write(*,*)'========= nombre de sujet, d obs meta, d obs loco ',ng,nsujetmeta,nsujet

    b=1.d-1

!---> initialisation des parametres par shared pour initialize == 1
    if(initialize == 1) then
        b(1:(np1-nva1-1))=b01(1:(np1-nva1-1))
        b((1+(np1-nva1-1)):(np1-nva1+np2-nva2-2))=b02(1:(np2-nva2-1)) !b02(1:(np2-nva2))
        b((1+(np1-nva1+np2-nva2-2)):(np1-nva1+np2-nva2+np3-nva3-3))=b03(1:(np3-nva3-1))
    end if

!---> initialisation des parametres associes aux var expli
    do i=1,nva1
        b(i+5+(np1-nva1+np2-nva2+np3-nva3-3))=b01(np1-nva1+i)
    end do
    do i=1,nva2
        b(i+5+(np1-nva1+np2-nva2+np3-nva3-3)+nva1)=b02(np2-nva2+i)
    end do
    do i=1,nva3
        b(i+5+(np1-nva1+np2-nva2+np3-nva3-3)+nva2+nva1)=b03(np3-nva3+i)
    end do        

    if((typeof == 0).and.(initialize==1)) then
        k0 = kappaCV    
    end if

! parametres de dependance
    b(np-nva-indic_alpha)=1.d0
    b(np-nva-indic_eta)=1.d0
    b(np-nva-indic_a1)=0.5d0
    b(np-nva-indic_a2)=0.5d0
    b(np-nva-indic_rho)=0.5d0
! fin parametres de dependance

! NEW initialisation to compare with frailtypack
!    b=0.5d0
! NEW initialisation to compare with frailtypack
!    write(*,*)'param init',b

    res=0.d0
    select case(typeof)
        case(0)
            call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivSplines)
        case(1)
            allocate(betacoef(nbintervR + nbintervDC + nbintervM))
            call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivCpm)
        case(2)
            call marq98(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaMultivWeib)
    end select

    if (typeof .ne. 0) then
        deallocate(kkapa)
    end if

    resOut=res
    critCV(1)=ier
    critCV(2)=istop

    if (istop .ne. 1) then
        goto 1000
    end if

    call multi(I_hess,H_hess,np,np,np,IH)
    call multi(H_hess,IH,np,np,np,HIH)

    if(effet.eq.1.and.ier.eq.-1)then
        v((np-nva-indic_alpha)*(np-nva-indic_alpha+1)/2)=10.d10
    endif

    res01(effet+1)=res

! --------------  Lambda and survival estimates JRG January 05
    select case(typeof)
        case(0)
            call distanceJ_splines(nzloco,nzdc,nzmeta,b,mt1,mt2,mt3,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out,&
            x3Out,lam3Out,su3Out)
        case(1)
            Call distanceJ_cpm(b,nbintervR+nbintervDC+nbintervM,mt1,mt2,mt3,x1Out,lamOut,xSu1,suOut,x2Out, &
            lam2Out,xSu2,su2Out,x3Out,lam3Out,xSu3,su3Out)
        case(2)
            Call distanceJ_weib(b,np,mt1,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out,x3Out,lam3Out,xSu3,su3Out)
            scale_weib(1) = etaR
            shape_weib(1) = betaR
            scale_weib(2) = etaD
            shape_weib(2) = betaD
            scale_weib(3) = etaM
            shape_weib(3) = betaM
    end select

    do ss=1,npmax
        do sss=1,npmax
            HIHOut(ss,sss) = HIH(ss,sss)
            H_hessOut(ss,sss)= H_hess(ss,sss)
        end do  
    end do

!AD:add LCV
!LCV(1) The approximate like cross-validation Criterion
!LCV(2) Akaike information Criterion 
!     calcul de la trace, pour le LCV (likelihood cross validation)
    LCV=0.d0
    if (typeof == 0) then
!        write(*,*)'The approximate like cross-validation Criterion in the non parametric case'
        call multi(H_hess,I_hess,np,np,np,HI)    
        do i =1,np
            LCV(1) = LCV(1) + HI(i,i)
        end do
        LCV(1) = (LCV(1) - resnonpen) / nsujet
    else
!        write(*,*)'=========> Akaike information Criterion <========='
        LCV(2) = (1.d0 / nsujet) *(np - resOut)
!        write(*,*)'======== AIC :',LCV(2)
    end if


1000 continue

!    write(*,*)'=========== coefBeta loco =========='
    coefBeta(1,:) = b((np-nva+1):(np-nva+nva1))
!    print*,coefBeta

!    write(*,*)'=========== coefBeta dc =========='
    coefBetadc(1,:) = b((np-nva+nva1+1):(np-nva+nva1+nva2))
!    print*,coefBetadc

!    write(*,*)'=========== coefBeta meta =========='
    coefBetaM(1,:) = b((np-nva+nva1+nva2+1):np)
!    print*,coefBetaM(1,:)

    do i=1,nsujet
        do j=1,nva1
            ve1(i,j)=ve(i,j)
        end do
    end do

    do i=1,ng
        do j=1,nva2
            ve2(i,j)=vedc(i,j)
        end do
    end do

    do i=1,nsujetmeta
        do j=1,nva3
            ve3(i,j)=vemeta(i,j)
        end do
    end do

    Xbeta = matmul(coefBeta,transpose(ve1))
    Xbetadc = matmul(coefBetadc,transpose(ve2))
    XbetaM = matmul(coefBetaM,transpose(ve3))

    if((istop == 1) .and. (effet == 1)) then
!        print*,'======== Call Residus Martingale ==========='
        deallocate(I_hess,H_hess)

        allocate(vres((2*(2+3)/2)),I_hess(2,2),H_hess(2,2))

        effetres = effet

        Call Residus_Martingale_multive(b,np,funcpamultires,Res_martingale,Res_martingaledc,Res_martingale2,&
        frailtypred,frailtypred2,frailtyvar,frailtyvar2,frailtyCorr)

        do i=1,nsujet
            linearpred(i)=Xbeta(1,i)+frailtypred(g(i))
        end do

        do i=1,ng
            linearpreddc(i)=Xbetadc(1,i)+alpha1*frailtypred(g(i))+alpha2*frailtypred2(gmeta(i))
        end do

        do i=1,nsujetmeta
            linearpredM(i)=XbetaM(1,i)+frailtypred2(gmeta(i))
        end do

        deallocate(I_hess,H_hess,vres)
    else
        deallocate(I_hess,H_hess)
    end if
!AD:end
    deallocate(ResidusRec,Residusdc,ResidusRec2,Rrec,Nrec,Rdc,Ndc,vuu,t0meta,t1meta,&
    cmeta,gmeta,auxmeta,vemeta,vectn,res1meta,res3meta,ve1,ve2,ve3,nig,nigmeta,cdc,t0dc,&
    t1dc,aux1,aux2,res1,res4,res3,mi,t0,t1,c,stra,g,aux,vax,vaxdc,vaxmeta,ve,vedc,filtre,&
    filtre2,filtre3,Hspl_hess,PEN_deri,hess,v,I1_hess,H1_hess,I2_hess,H2_hess,HI2,HIH,IH,&
    HI,BIAIS,date,datedc,datemeta,b01,b02,b03,Rrec2,Nrec2)

    if (typeof == 0) then
        deallocate(mm3meta,immeta,mm2meta,mm1meta,mmmeta,im3meta,im2meta,im1meta,&
        nt0dc,nt1dc,nt0,nt1,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc,mm3,mm2,&
        mm1,mm,im3,im2,im1,im,zi,zidc,zimeta,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,&
        m2m,m1m,nt0meta,nt1meta)
        deallocate(m3m3b,m2m2b,m1m1b,mmmb,m3m2b,m3m1b,m3mb,m2m1b,m2mb,m1mb,m3m3c,&
        m2m2c,m1m1c,mmmc,m3m2c,m3m1c,m3mc,m2m1c,m2mc,m1mc)
    end if

    if (typeof .ne. 0) then
        deallocate(vvv)
    end if

    if (typeof == 1) then
        deallocate(ttt,tttdc,tttmeta,betacoef)
    end if

    cptEvent(1)=cpt
    cptEvent(2)=cpt_dc
    cptEvent(3)=cptmeta

    ResMartingaleEvent(1:nobsEvent(3),1)=Res_martingale(1:nobsEvent(3))
    ResMartingaleEvent(1:nobsEvent(3),2)=Res_martingaledc(1:nobsEvent(3))
    ResMartingaleEvent(1:nobsEvent(3),3)=Res_martingale2(1:nobsEvent(3))

    frailtyEstimates(1:nobsEvent(3),1)=frailtypred(1:nobsEvent(3))
    frailtyEstimates(1:nobsEvent(3),2)=frailtypred2(1:nobsEvent(3))
    frailtyEstimates(1:nobsEvent(3),3)=frailtyvar(1:nobsEvent(3))
    frailtyEstimates(1:nobsEvent(3),4)=frailtyvar2(1:nobsEvent(3))
    frailtyEstimates(1:nobsEvent(3),5)=frailtyCorr(1:nobsEvent(3))

    return

    end subroutine joint_multiv


!========================== VECSPLI ==============================
!AD:add argument:ndatedc 
    subroutine vecspli(ndate,ndatedc,ndatemeta)
!AD:end
    use taillesmultiv
!AD:
    use comonmultiv,only:date,datedc,datemeta,zi,zimeta,mm3,mm2,mm1,mm,im3,im2,im1,im &
    ,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc,mm3meta,mm2meta,mm1meta,&
    mmmeta,im3meta,im2meta,im1meta,immeta,vectn!zidc,

!AD:end
    IMPLICIT NONE

    integer,intent(in)::ndate,ndatedc,ndatemeta
    integer::i,j,k
    double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2

!----------  calcul de u(ti) :  STRATE1 ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

    j=0
    do i=1,ndate-1
        do k = 2,vectn(1)-2
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
        im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0) &
            +(h4*mm(i)*0.25d0)
        im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
        im(i)  = ht*mm(i)*0.25d0

    end do
!AD: add for death 
!----------  calcul de u(ti) :  STRATE2 ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)
    j=0
    do i=1,ndatedc-1
        do k = 2,vectn(2)-2
            if ((datedc(i).ge.zi(k-1)).and.(datedc(i).lt.zi(k)))then
                j = k-1
            endif
        end do 
        ht = datedc(i)-zi(j)
        htm= datedc(i)-zi(j-1)
        h2t= datedc(i)-zi(j+2)
        ht2 = zi(j+1)-datedc(i)
        ht3 = zi(j+3)-datedc(i)
        hht = datedc(i)-zi(j-2)
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
        mm3dc(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        mm2dc(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
        *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        mm1dc(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
        h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        mmdc(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        im3dc(i) = (0.25d0*(datedc(i)-zi(j-3))*mm3dc(i))+(0.25d0*hh2 &
        *mm2dc(i))+(0.25d0*h3m*mm1dc(i))+(0.25d0*h4*mmdc(i))
        im2dc(i) = (0.25d0*hht*mm2dc(i))+(h3m*mm1dc(i)*0.25d0) &
            +(h4*mmdc(i)*0.25d0)
        im1dc(i) = (htm*mm1dc(i)*0.25d0)+(h4*mmdc(i)*0.25d0)
        imdc(i)  = ht*mmdc(i)*0.25d0

    end do

!AD:end
    j=0
    do i=1,ndatemeta-1
        do k = 2,vectn(3)-2
            if ((datemeta(i).ge.zimeta(k-1)).and.(datemeta(i).lt.zimeta(k)))then
                j = k-1
            endif
        end do 
        ht = datemeta(i)-zimeta(j)
        htm= datemeta(i)-zimeta(j-1)
        h2t= datemeta(i)-zimeta(j+2)
        ht2 = zimeta(j+1)-datemeta(i)
        ht3 = zimeta(j+3)-datemeta(i)
        hht = datemeta(i)-zimeta(j-2)
        h = zimeta(j+1)-zimeta(j)
        hh= zimeta(j+1)-zimeta(j-1)
        h2= zimeta(j+2)-zimeta(j)
        h3= zimeta(j+3)-zimeta(j)
        h4= zimeta(j+4)-zimeta(j)
        h3m= zimeta(j+3)-zimeta(j-1)
        h2n=zimeta(j+2)-zimeta(j-1)
        hn= zimeta(j+1)-zimeta(j-2)
        hh3 = zimeta(j+1)-zimeta(j-3)
        hh2 = zimeta(j+2)-zimeta(j-2)
        mm3meta(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        mm2meta(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
        *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n)) 
        mm1meta(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
        h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        mmmeta(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        im3meta(i) = (0.25d0*(datemeta(i)-zimeta(j-3))*mm3meta(i))+(0.25d0*hh2  &
        *mm2meta(i))+(0.25d0*h3m*mm1meta(i))+(0.25d0*h4*mmmeta(i))
        im2meta(i) = (0.25d0*hht*mm2meta(i))+(h3m*mm1meta(i)*0.25d0)+(h4*mmmeta(i)*0.25d0)
        im1meta(i) = (htm*mm1meta(i)*0.25d0)+(h4*mmmeta(i)*0.25d0)
        immeta(i)  = ht*mmmeta(i)*0.25d0
         end do

    end subroutine vecspli


    subroutine vecpenP(n,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m) 

    IMPLICIT NONE

    integer,intent(in)::n
    integer::i
    double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2,a3,a2,b2 &
    ,c2,a1,b1,c1,a0,x3,x2,x
    double precision,dimension(-2:(n+1))::zi
    double precision,dimension(n-2),intent(inout)::m3m3,m2m2,m1m1,mmm,m3m2, &
    m3m1,m3m,m2m1,m2m,m1m

!*********************************************************************
        m3m3=0.d0
    m2m2=0.d0
    m1m1=0.d0
    mmm=0.d0
    m3m2=0.d0
    m3m1=0.d0
    m3m=0.d0
    m2m1=0.d0
    m2m=0.d0
    m1m=0.d0
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
        m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2)  &
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
        (2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)) )
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

    end subroutine vecpenP




!==========================  SUSP  ====================================
    subroutine susp(x,the,n,su,lam,zi)

    use taillesmultiv

    IMPLICIT NONE 

    integer,intent(in)::n
    double precision,intent(out)::lam,su
    double precision,dimension(-2:npmax),intent(in)::zi,the
    double precision,intent(in)::x
    integer::j,k,i
    double precision::ht,ht2,h2,som,htm,h2t,h3,h2n,hn, &
    im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2 &
    ,h,gl,hh

    gl=0.d0
    som = 0.d0
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

    end subroutine susp

!==========================  COSP  ====================================
! calcul les points pour les fonctions
! et leur bandes de confiance

    subroutine cosp(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::n
    double precision,intent(in)::x
    double precision,intent(out)::lam,su
    double precision,intent(out)::binf,bsup,lbinf,lbsup
    double precision,dimension(npmax,npmax),intent(in)::y
    double precision,dimension(-2:npmax),intent(in)::the,zi
    integer::j,k,i
    double precision::ht,ht2,h2,som,pm,htm,h2t,h3,h2n,hn, &
    im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2, &
    h,gl,hh

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

    call conf(x,j,n,y,pm,zi)

    binf = dexp(-gl + 1.96d0*pm)
    su  = dexp(-gl)
    bsup = dexp(-gl - 1.96d0*pm)

    call conf1(x,j,n,y,pm,zi)
    lbinf = lam - 1.96d0*pm
    lbsup = lam + 1.96d0*pm
!         write(*,*)'lbinf apres conf1',lbinf,lam,pm

    return

    end subroutine cosp


!=====================  CONF1  =============================


    subroutine  conf1(x,ni,n,y,pm,zi)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::ni,n
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision,dimension(npmax,npmax),intent(in)::y
    double precision,intent(out)::pm
    integer::i,j
    double precision::res,mmsp
    double precision,dimension(npmax)::vecti,aux

    do i=1,n
        vecti(i) = mmsp(x,ni,i,zi)
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

    res = -res
    pm = dsqrt(res)

    end subroutine  conf1

!=====================  CONF  =============================

    subroutine  conf(x,ni,n,y,pm,zi)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::ni,n
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision,dimension(npmax,npmax),intent(in)::y
    double precision,intent(out)::pm
    integer::i,j
    double precision::res,isp
    double precision,dimension(52)::vecti,aux

    do i=1,n
        vecti(i) = isp(x,ni,i,zi)
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

    end subroutine  conf


!==========================   ISP   ==================================

    double precision function isp(x,ni,ns,zi)

    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::ni,ns
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision::val,mmsp



    if(x.eq.zi(ni))then
        if(ni.le.ns-3)then
            val = 0.d0
            else
                if(ni.le.ns-2)then
                    val = ((zi(ni)-zi(ni-1))*mmsp(x,ni,ns,zi))*0.25d0
                else
                    if (ni.eq.ns-1)then
                        val = ((zi(ni)-zi(ni-2))*mmsp(x,ni,ns,zi)+ &
                        (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val = ((zi(ni)-zi(ni-3))*mmsp(x,ni,ns,zi)+ &
                            (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi) &
                            +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi))*0.25d0
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
                val = (x-zi(ni))*mmsp(x,ni,ns,zi)*0.25d0
            else
                if(ni.eq.ns-2)then
                    val = ((x-zi(ni-1))*mmsp(x,ni,ns,zi)+ &
                    (zi(ni+4)-zi(ni))*mmsp(x,ni,ns+1,zi))*0.25d0
                else
                    if (ni.eq.ns-1)then
                        val =((x-zi(ni-2))*mmsp(x,ni,ns,zi)+ &
                        (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi) &
                        +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+2,zi))*0.25d0
                    else
                        if(ni.eq.ns)then
                            val =((x-zi(ni-3))*mmsp(x,ni,ns,zi)+ &
                            (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi) &
                            +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi) &
                            +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+3,zi))*0.25d0
                        else
                            val = 1.d0
                        endif
                    endif
                endif
            endif
        endif 
    endif

    isp = val

    return

    end function isp

!==========================  MMSP   ==================================

    double precision function mmsp(x,ni,ns,zi)

    use taillesmultiv

    IMPLICIT NONE 

    integer,intent(in)::ni,ns
    double precision,intent(in)::x
    double precision,dimension(-2:npmax),intent(in)::zi
    double precision::val

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
                    -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) &
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

    mmsp = val

    return

    end function mmsp


!================== multiplication de matrice  ==================

! multiplie A par B avec le resultat dans C

    subroutine multi(A,B,IrowA,JcolA,JcolB,C)
!     remarque :  jcolA=IrowB
    use taillesmultiv

    IMPLICIT NONE

    integer,intent(in)::IrowA,JcolA,JcolB
    double precision,dimension(npmax,npmax),intent(in):: A,B
    double precision,dimension(npmax,npmax),intent(out)::C
    integer::i,j,k
    double precision::sum

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

    end subroutine multi


    double precision  function func30(frail1,frail2)
! calcul de l integrant, pour un effet aleatoire donne frail et un groupe donne auxig (cf funcpa)

    use taillesmultiv
    use comonmultiv,only:nigmeta,alpha1,alpha2,res1meta,res3meta,&
    nig,auxig,alpha,theta,eta,aux1,res1,res3,cdc 

    implicit none

    double precision::frail1,frail2

    func30=0.d0
    func30 = frail1*(cdc(auxig)*alpha1+nig(auxig)) &
    +frail2*(cdc(auxig)*alpha2+nigmeta(auxig)) &
    -dexp(frail1)*(res1(auxig)-res3(auxig))-dexp(frail2)*(res1meta(auxig)-res3meta(auxig)) &
    -dexp(frail1*alpha1+frail2*alpha2)*aux1(auxig) &
    +(2.d0*((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0) &
    *frail1*frail2/sqrt(theta*eta)  &
    -(frail1**2.d0)/theta -(frail2**2.d0)/eta) &
    /(2.d0*(1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2.d0))
    func30 = dexp(func30)

    return

    end function func30

