!mtaille =c(mt1,mt2,mt11,mt12)
!paraweib =c(shapeweib(1),shapeweib(2),scaleweib(1),scaleweib(2))
!tempdc=matrice(tt0dc0,tt1dc0) pour les deces
!kendall: 1ere colonne ss0, 2eme colonne tau
!paratps = c(timedep0,nbinnerknots,qorder0)

!--entête pour fortran
    subroutine joint_surrogate(nsujet0,ng0,ntrials0,lignedc0,nz0,nst0,k0,tt00,tt10,ic0,groupe0,trials,pourtrial0,nigs0,&
                               cdcs0,groupe00,tt0dc0,tt1dc0,icdc0,tempdc,icdc00,nva10,vax0,nva20,vaxdc0,vaxdc00,noVar1,&
                               noVar2,ag0,maxit0,np,b,H_hessOut,HIHOut,resOut,LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,&
                               xSu2,su2Out,typeof0,equidistant,nbintervR0,nbintervDC0,mtaille,ni,cpt,cpt_dc,ier,istop,&
                               paraweib,MartinGales,linearpred,linearpreddc,ziOut,time,timedc,linearpredG,typeJoint0,&
                               intcens0,ttU0,logNormal0,paratps,filtretps0,BetaTpsMat,BetaTpsMatDc,EPS,nsim_nodes,indice_esti,&
                               indice_covST0,paGH,param_weibull0)

!AD: add for new marq
    use parameters
    use splines
    use comon
    use tailles

    use optim_scl_0
!AD:pour fortran
    use sortie
    use residusM
    use comongroup
    use betatttps
    use var_surrogate ! pour le surrogate
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    !use mpi ! module pour l'environnement MPI
!AD:
    !use mod_Adaptative !pour l'estimation des fragilites
    implicit none
    
    integer, dimension(4), intent(in)::indice_esti
    integer::maxit0,mt1,mt2,mt11,mt12,nnodes,aaa,control !nn,npinit,nvatmp 
    integer,dimension(4),intent(in)::mtaille
    integer,intent(in)::nsujet0,ng0,nz0,nva10,nva20,lignedc0,ag0,ntrials0,nst0
    integer,intent(in),dimension(ntrials0)::trials
    double precision,dimension(nz0+6),intent(out)::ziOut
    integer::np,equidistant
    integer,dimension(nsujet0),intent(in)::groupe0,ic0
    integer,dimension(ng0),intent(in)::icdc0
    integer,dimension(lignedc0),intent(in)::groupe00,icdc00

    double precision,dimension(ng0)::tt0dc0,tt1dc0
    double precision,dimension(lignedc0,2)::tempdc
    double precision,dimension(nsujet0)::tt00,tt10,ttU0 !! rajout
    double precision,dimension(2)::k0
    double precision,dimension(nsujet0,nva10),intent(in):: vax0
    double precision,dimension(ng0,nva20),intent(in):: vaxdc0
    double precision,dimension(lignedc0,nva20),intent(in):: vaxdc00
    double precision,dimension(np,np)::H_hessOut,HIHOut ! H_hessOut=matrice des variances-covariances
    double precision::resOut,tp2,tp1
    double precision,dimension(mtaille(1))::x1Out
    double precision,dimension(mtaille(2))::x2Out
    double precision,dimension(mtaille(1),3)::lamOut
    double precision,dimension(mtaille(3),3)::suOut
    double precision,dimension(mtaille(2),3)::lam2Out
    double precision,dimension(mtaille(4),3)::su2Out
    integer::ss,sss
    double precision,dimension(np), intent(inout):: b!,b_save
    double precision,dimension(2),intent(out)::LCV
    double precision,dimension(2)::shapeweib,scaleweib
    double precision,dimension(4),intent(out)::paraweib
    
    integer,intent(in)::noVar1,noVar2,intcens0,param_weibull0 !param_weibull! parametrisation de la weibull utilisee: 0= parametrisation par defaut dans le programme de Virginie, 1= parametrisation 
    !                                                         a l'aide de la fonction de weibull donnee dans le cours de Pierre
    integer,intent(out)::cpt,cpt_dc,ier,ni
    integer::groupe,ij,kk,j,k,nz,n,ii,iii,iii2,cptstr1,cptstr2   & !code
    ,i,ic,icdc,istop,cptni,cptni1,cptni2,nb_echec,nb_echecor,id,cptbiais &
    ,cptauxdc,p !rang,erreur
    double precision::tt0,tt0dc,tt1,tt1dc,h,hdc,res,min,mindc,max_, &
    maxdc,maxt,maxtdc,moy_peh0,moy_peh1,lrs,BIAIS_moy,ttU,mintdc !! rajout
    double precision,dimension(2)::res01
!AD: add for new marq
    double precision::ca,cb,dd !result_
    double precision,external::funcpajsplines_surrogate,funcpajsplines_surrogate_1,funcpajcpm,funcpajweib
    double precision,external::funcpajsplines_intcens,funcpajweib_intcens
    double precision,external::funcpajsplines_log,funcpajcpm_log,funcpajweib_log
    double precision,external::funcpaGsplines,funcpaGcpm,funcpaGweib
    double precision,external::funcpaGsplines_intcens,funcpaGcpm_intcens,funcpaGweib_intcens
    double precision,external::funcpaGsplines_log,funcpaGcpm_log,funcpaGweib_log
    double precision,external::funcpaj_tps,funcpaG_tps,funcpajsplines_copule_surrogate
    double precision,dimension(100)::xSu1,xSu2
!cpm
    integer::indd,ent,entdc,typeof0,nbintervR0,nbintervDC0, np_2
    double precision::temp    
!predictor
    double precision,dimension(ng0)::Resmartingale,Resmartingaledc,frailtypred,frailtyvar
    double precision,dimension(ng0,4),intent(out)::MartinGales

    double precision,external::funcpajres,funcpajres_log
    double precision,dimension(nsujet0),intent(out)::linearpred
    double precision,dimension(ng0),intent(out)::linearpreddc
    double precision,dimension(lignedc0),intent(out)::linearpredG
    double precision,dimension(1,nva10)::coefBeta
    double precision,dimension(1,nva20)::coefBetadc
    double precision::coefBeta2
    double precision,dimension(1,nsujet0)::XBeta

    double precision,dimension(1,ng0)::XBetadc
    double precision,dimension(1,lignedc0)::XBetaG

    double precision,dimension(nbintervR0+1)::time
    double precision,dimension(nbintervDC0+1)::timedc


!    double precision,dimension(4,2),intent(out)::kendall!ss0,tau
!    double precision::sstmp
    integer,intent(in)::typeJoint0 !initialisation
!    double precision,dimension(nn)::Bshared
    integer::ngtemp
    integer,intent(in)::logNormal0,indice_covST0

    integer,dimension(3),intent(in)::paratps
    integer,dimension(nva10+nva20),intent(in)::filtretps0
    double precision,dimension(0:100,0:4*sum(filtretps0(1:nva10)))::BetaTpsMat !!! a refaire
    double precision,dimension(0:100,0:4*sum(filtretps0(nva10+1:nva10+nva20)))::BetaTpsMatDc
    double precision,dimension(paratps(2)+paratps(3))::basis
    double precision,dimension(3),intent(inout)::EPS ! seuils de convergence : on recupere les valeurs obtenues lors de l'algorithme a la fin
    integer, dimension(13),intent(in)::nsim_nodes !scl nsim_nodes: vecteur contenant le nbre de simulation(1) pour le MC et de noeud(2) pour la quadrature,le troisieme element indique si on fait l'adaptative(1) ou la non adaptative(0), le quatrieme indique la methode d'integration
                                                 !0=Monte carlo,1= MC+quadrature, 2=quadrature, le cinquieme le nombre de parametres associes a la fragilite
                                                 ! le septieme indique le nombre d'effet aleatoire dans le cas de la quadrature adaptative
                                                 ! et le huitieme le type de modele a estimer (0=joint surrogate classique,1=joint surrogate complet)
    double precision,dimension(ng0,nsim_nodes(7)+1+nsim_nodes(7) + (nsim_nodes(7)*(nsim_nodes(7)-1))/2),intent(in):: paGH ! parametre pour l'adaptative: en ligne les individus, en colone on a respectivement: les ui_cham,racine carree du determinant de l'inverse de la cholesky,variance des ui_chap,les covariances estimees des fragilites pour chaque individu, sachant que la matrice de variances covariance est bien la cholesky                                                        

    integer, dimension(ntrials0,2),intent(in)::nigs0,cdcs0
    integer, dimension(nsujet0),intent(in)::pourtrial0
    double precision::bi,bs,wres
    integer::rangparam,nb_ree,nb0,indice_eta0,n_par_pro,suplement ! nb0= nobre effet aleatoire pour adaptative
    
    !recherche du rang du processus courant
    !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
    
    param_weibull=param_weibull0
    estim_wij_chap=0
    control_affichage = 0
    control_adaptative_laplace = 0
    
   ! 100 continue
    !rang=0
    control=0  ! donnee de controle de l'estimation des effets aleatoires
    
    cpteu=0 ! variable globale utilisee juste pour des tests
    mt1=mtaille(1)
    mt2=mtaille(2)
    mt11=mtaille(3)
    mt12=mtaille(4)

    allocate(vaxdc(nva20),vax(nva10))

    timedep = paratps(1)
    nbinnerknots = paratps(2)
    qorder = paratps(3)

    allocate(filtretps(nva10),filtre2tps(nva20))
    allocate(betatps(nva10),betatps2(nva20))
    filtretps = filtretps0(1:nva10)
    filtre2tps = filtretps0(nva10+1:nva10+nva20)

    npbetatps1 = (nbinnerknots+qorder-1)*sum(filtretps)
    npbetatps2 = (nbinnerknots+qorder-1)*sum(filtre2tps)
    npbetatps = npbetatps1 + npbetatps2

    !!print*,"logNormal0=",logNormal0
    logNormal = logNormal0
    type_mod_surr=logNormal0 !scl
    intcens = intcens0
    typeJoint = typeJoint0
    ag = ag0
    typeof = typeof0
    !model = 1  !indique le type de modele utilise
    model = 8 !scl pour le model surrogate
    indic_alpha = 0
    type_joint=nsim_nodes(8)
    !!print*,"type_mod_surr=",type_mod_surr
    !!print*,"type_mod_surr=",type_mod_surr
    !stop
!    if (type_mod_surr.eq.1) then ! scl ici on est dans le cas du modele avec tous les effets aleatoires gaussiens
      aaa=0
      indice_eta0=indice_esti(1)
      frailt_base=indice_esti(2)
        if(type_joint==2)then
            indice_eta = 0
        else
            if(indice_eta0==0)then
                indice_eta = 0 ! on fixe eta à 1
            else
                indice_eta = 1 
            endif
        endif      
      
      !!print*,"indice_eta=",indice_eta
      if(type_joint .ne.3) indice_theta=1
      indice_alpha=1
      indice_sigma=1
      indice_varS=1
      indice_varT=1
      indice_theta_t=1
      indice_theta_st=1
    if(frailt_base==1) then
        indice_gamma=1
        indice_gamma_t=1
        indice_gamma_st=indice_esti(4)
        indice_alpha_ui=indice_esti(3)
        if(type_joint==1 .or. type_joint==3) then ! modele a fragilites partages
            allocate(chol(3,3))
        else ! modele complet
            allocate(chol(6,6))
        endif
    else
        allocate(chol(2,2))
    endif
      if(indice_covST0==1)then
        indice_covST=1
      else
        indice_covST=0
      endif
      Nmax=nsujet0
      nsim=nsim_nodes(1)
      methodInt=nsim_nodes(4) ! on suppose estimer les integrales multidimensionnelles par monte carlo, bien evidemment il faudrait que le nombre de simulation soit positif ie nsim0>0
      ntrials=ntrials0
      allocate(nsujeti(ntrials0),delta(nsujet0),deltastar(ng0),const_res4(nsujet0),const_res5(ng0),const_res1(ntrials)&
      ,const_aux1(ntrials),nigs(ntrials),cdcs(ntrials),nigts(ntrials),cdcts(ntrials),pourtrial(nsujet0),res2s(ntrials)&
      ,res2_dcs(ntrials),res2s_sujet(nsujet0),res2_dcs_sujet(nsujet0))      
      nsujeti=trials ! nombre de sujet par cluster
      sujet_essai_max=MAXVAL(nsujeti(1:ntrials))
      nparamfrail=nsim_nodes(5)
      !if(rang==0) !print*,"Nombre de sujet dans le plus grand cluster:",sujet_essai_max
       if(methodInt .eq.0 .or. methodInt .eq.2 .or. methodInt .eq.4) then ! on alloue le vecteur des elements a simuler pour MC seulement si on fait du MC
        !allocate(Vect_sim_MC(nsim,sujet_essai_max+1))! +1 car on integre au plus sur le nombre de sujet +1 pour l'effet aleatoire niveau essai
        allocate(Vect_sim_MC(nsim,nparamfrail))
        a_deja_simul=0
       endif
       
    if(methodInt==3) then !integration par laplace
        allocate(wij_chap(nsujet0,1),control_wij_chap(nsujet0))
        control_wij_chap=0
        if(frailt_base == 1) then
            np_2 = 3
            
        else
            np_2 = 2
        endif
        allocate(IhessLaplace(np_2,np_2),H_hess_laplace(np_2,np_2),&
                b_i_laplace(np_2),v_i_laplace(np_2*(np_2+3)/2),hess_laplace(np_2,np_2),vvv_laplace(np_2*(np_2+1)/2))
        ! !print*,"suis dans joint",size(wij_chap),size(wij_chap,1),size(wij_chap,2)
    endif
    
      !!print*,"suis dans le joint"
      nigs=nigs0(:,1) ! nombre de sujets avec une progression par essai
      cdcs=cdcs0(:,1) ! nombre de sujets decedes par essai
      nigts=nigs0(:,2) ! nombre de sujets avec une progression par essai en interaction avec le traitement
      cdcts=cdcs0(:,2) ! nombre de sujets decedes par essai en interaction avec le traitement
      !!print*,"nigs=",nigs
      !!print*,"nigts=",nigts
      !!print*,"cdcts=",cdcts
      !!print*,"cdcs=",cdcs
      
      pourtrial=pourtrial0 ! indice des essais par sujet
      delta=ic0
      deltastar=icdc0
      !!print*,"size(deltastar)=",size(deltastar),"size(icdc0)=",size(icdc0)
      !!print*,"size(delta)=",size(delta),"size(ic0)=",size(ic0)
      npoint=nsim_nodes(2) 
      vectorisation=nsim_nodes(6) !permet de reduire le temps de calcul
      const_res4=0.d0
      const_res5=0.d0
      allocate(xx1(npoint))
      allocate(ww1(npoint))
      nnodes=npoint
        if(logNormal==1)then
                if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
                else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
                else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
                else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
                else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
                else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
                else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
                end if
            else ! guass laguerre
                if (nnodes.eq.20) then
                    xx1(1:nnodes) = x(1:nnodes)
                    ww1(1:nnodes) = w(1:nnodes)
                else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x1(1:nnodes)
                    ww1(1:nnodes) = w1(1:nnodes)
                endif
            endif
            
      !!print*,"nsim_nodes(3)",nsim_nodes(3)
      ! stop
 
        if(nsim_nodes(3).ne.0) then
            adaptative=.true.
        else
            adaptative=.false.
        end if
        ! controle si l'on fait de l'adaptative ou de la pseudo adaptative
        if(nsim_nodes(3).ne.0)then
            control_adaptative=1
            switch_adaptative=1 ! en cas de non convergence de la procedure d'estimation des effets aleatoires a posteriorie, on met cette variable a 0 pour dire de ne pas considerer la pseudo adaptative
        else
            control_adaptative=0
            switch_adaptative=0
        endif
    !endif
    
    if (typeof .ne. 0) then
        nbintervR = nbintervR0
        nbintervDC = nbintervDC0
    end if

    lignedc = lignedc0

    maxiter = maxit0
!AD:add for new marq
    epsa = EPS(1) !1.d-4
    epsb = EPS(2) !1.d-4
    epsd = EPS(3) !1.d-4
!     if (intcens0.eq.1) then
!         epsa = 1.d-2
!         epsb = 1.d-2
!         epsd = 1.d-2
!     endif
!AD:end

    lrs = 0.d0
    moy_peh0 = 0.d0
    moy_peh1 = 0.d0

    nb_echec = 0
    nb_echecor = 0
    nb0recu = 0
    moyrecu =0.d0

    ngmax=ng0
    ng=ng0
    if(typeJoint==1) then
        ngtemp=ng
    else
        ngtemp=lignedc
    end if
!Prise en compte des temps dc
    allocate(temps0dc(ngtemp),temps1dc(ngtemp),ictemp(ngtemp),variable(ngtemp,nva20))

    if(typeJoint==1)then
        temps0dc=tt0dc0
        temps1dc=tt1dc0
        ictemp=icdc0
        variable=vaxdc0
    else
        temps0dc=tempdc(:,1)
        temps1dc=tempdc(:,2)
        ictemp=icdc00
        variable=vaxdc00
    endif
!end
    allocate(ResidusRec(ngtemp),Residusdc(ngtemp),Rrec(ngtemp),Nrec(ngtemp),Rdc(ngtemp),Ndc(ngtemp),vuu(1))
    allocate(cdc(ngtemp),t0dc(ngtemp),t1dc(ngtemp),aux1(ngtemp),aux2(ngtemp) &
    ,res1(ngtemp),res4(ngtemp),res3(ngtemp),mi(ngtemp))

    allocate(nig(ngtemp),nigdc(ng))
    shapeweib = 0.d0
    scaleweib = 0.d0
    nsujetmax=nsujet0
    nsujet=nsujet0

!Al: utile pour le calcul d'integrale avec distribution log normale
    allocate(res5(nsujetmax))
!Al

    allocate(t0(nsujetmax),t1(nsujetmax),tU(nsujetmax),c(nsujetmax),stra(nsujetmax),g(nsujetmax),aux(3*nsujetmax),gsuj(nsujetmax)) !! chgt dimension aux

! censure par intervalle
    allocate(resL(nsujetmax),resU(nsujetmax))
    
    if(typeJoint==1) then
        ndatemaxdc=2*ng0
    else
        ndatemaxdc=2*lignedc
    end if

    if (typeof == 0) then
        allocate(nt0dc(ngtemp),nt1dc(ngtemp),nt0(nsujetmax),nt1(nsujetmax),ntU(nsujetmax))!! rajout
        allocate(mm3dc(ndatemaxdc),mm2dc(ndatemaxdc),mm1dc(ndatemaxdc),mmdc(ndatemaxdc) &
        ,im3dc(ndatemaxdc),im2dc(ndatemaxdc),im1dc(ndatemaxdc),imdc(ndatemaxdc))
    end if
    nst=2
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
    n=0
    nz=0
!**************************************************
!********************* prog spline****************
    effet=1
    res01(1)=0.d0
    res01(2)=0.d0
!------------  entre non fichier et nombre sujet -----
    nvarmax=ver
    nva1=nva10
    nva2=nva20
    nva = nva1+nva2
    !===========Fin scl: 12/04:2019================
    !nvarmax=nva
    !allocate(ve(nsujetmax,nvarmax),vedc(ngtemp,nvarmax))
    allocate(ve(nsujetmax,nva1),vedc(ngtemp,nva2))
    !===========Fin scl: 12/04:2019================
    allocate(ve1(nsujetmax,nva1),ve2(ngtemp,nva2))
    allocate(filtre(nva10),filtre2(nva20))
    nig=0
! AD: recurrent
    if (noVar1.eq.1) then 
!        !write(*,*)'filtre 1 desactive'
        filtre=0
        nva1=0
    else
        filtre=1
    end if
!AD:death
    if (noVar2.eq.1) then 
!        !write(*,*)'filtre 2 desactive'
        filtre2=0
        nva2=0
    else
        filtre2=1
    end if
    
    if ((noVar1.eq.1).or.(noVar2.eq.1)) then
        nva = nva1+nva2
    end if

!AD:end

!------------  lecture fichier -----------------------
    maxt = 0.d0
    mint = 0.d0
    
    maxtdc = 0.d0
    mintdc = 0.d0
    
    cpt = 0
    cptcens = 0
    cpt_dc = 0
    k = 0
    cptstr1 = 0
    cptstr2 = 0
    nigdc = 0

!ccccccccccccccccccccc
! pour le deces
!cccccccccccccccccccc
    do k = 1,ngtemp
    
        if (k.eq.1) then
            mintdc = temps0dc(k) ! affectation du min juste une fois
        endif
        
        tt0dc=temps0dc(k)
        tt1dc=temps1dc(k)
        icdc=ictemp(k)
        if(typeJoint==1)then
            !groupe=groupe0(k) !!!! inutile et surtout faux car groupe0 contient les groupes
            !!!! pour les evenements et non les deces (les lignes ne sont pas les memes)
        else
            groupe=groupe00(k)
        end if
        do j=1,nva20
            vaxdc(j)=variable(k,j)
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
            if(typeJoint.ne.1) then
                gsuj(k) = groupe
                nigdc(groupe) = nigdc(groupe)+1 ! nb de deces par groupe
            endif
            iii = 0
            iii2 = 0
            !===========scl: 12/04:2019================
            ! do ii = 1,nva20
                ! if(filtre2(ii).eq.1)then
                    ! iii2 = iii2 + 1
                    ! vedc(k,iii2) = dble(vaxdc(ii))
                ! endif
            ! end do
            vedc(k,:) = dble(vaxdc)
            !===========Fin scl: 12/04:2019================
        else
!------------------   censure a droite ou event recurr  c=0
            if(icdc.eq.0)then
                cdc(k) = 0
                iii = 0
                iii2 = 0
                !===========scl: 12/04:2019================
                ! do ii = 1,nva20
                    ! if(filtre2(ii).eq.1)then
                    ! iii2 = iii2 + 1
                    ! vedc(k,iii2) = dble(vaxdc(ii))
                    ! endif
                ! end do 
                vedc(k,:) = dble(vaxdc) 
                !===========Fin scl: 12/04:2019================
                t0dc(k) = tt0dc
                t1dc(k) = tt1dc
                if(typeJoint.ne.1) gsuj(k) = groupe
            endif
        endif
        
        if (maxtdc.lt.t1dc(k))then
            maxtdc = t1dc(k)
        endif
        
        if (mintdc.gt.t0dc(k)) then
            mintdc = t0dc(k)
        endif
    end do

    deallocate(temps0dc,temps1dc,ictemp,variable)
!AD:
    if (typeof .ne. 0) then
        cens = maxtdc
    end if
!Ad
    k = 0
    cptstr1 = 0
    cptstr2 = 0

!cccccccccccccccccccccccccccccccccc
! pour les donnees recurrentes
!cccccccccccccccccccccccccccccccccc
    do i = 1,nsujet     !sur les observations
    
        if (i.eq.1) then
            mint = tt00(i) ! affectation du min juste une fois
        endif
        
        tt0=tt00(i)
        tt1=tt10(i)
        ttU=ttU0(i) !! rajout
        ic=ic0(i)
        groupe=groupe0(i)
        !call dblepr("suis danc funcpan vax0=", -1, dble(vax0(i,:)), size(vax0,2))
!------------------
        do j=1,nva10
            vax(j)=vax0(i,j) ! ensemble des observation du sujet i associees au surrrogate
        enddo
!--------------
        if(tt0.gt.0.d0)then
            cptaux=cptaux+1
        endif
!-----------------------------------------------------

!     essai sans troncature
!     tt0=0
!------------------   observation c=1 pour donnees recurrentes

        if(ic.eq.1)then
            cpt = cpt + 1
            c(i)=1
            t0(i) = tt0
            t1(i) = tt1
            tU(i) = ttU !! rajout
            t1(i) = t1(i)
            g(i) = groupe
            nig(groupe) = nig(groupe)+1 ! nb d event recurr dans un groupe
            iii = 0
            iii2 = 0
            !===========scl: 12/04:2019================
            ! do ii = 1,nva10
                ! if(filtre(ii).eq.1)then
                    ! iii = iii + 1
                    ! ve(i,iii) = dble(vax(ii)) !ici sur les observations
                ! endif
            ! end do
            ve(i,:) = dble(vax)
            !===========Fin scl: 12/04:2019================
        else
!------------------   censure a droite  c=0 pour donnees recurrentes
            if(ic.eq.0)then
                cptcens=cptcens+1
                c(i) = 0
                iii = 0
                iii2 = 0
                !===========scl: 12/04:2019================
                ! do ii = 1,nva10
                    ! if(filtre(ii).eq.1)then
                    ! iii = iii + 1
                    ! ve(i,iii) = dble(vax(ii))
                    ! endif
                ! end do
                ve(i,:) = dble(vax)
                !===========Fin scl: 12/04:2019================
                t0(i) =  tt0
                t1(i) = tt1
                tU(i) = ttU !! rajout
                t1(i) = t1(i)
                g(i) = groupe
            endif
        endif
        if (maxt.lt.t1(i))then
            maxt = t1(i)
        endif

        if ((maxt.lt.tU(i)).and.(tU(i).ne.t1(i))) then
            maxt = tU(i)
        endif
        
        if (mint.gt.t0(i)) then
            mint = t0(i)
        endif
    end do


    deallocate(filtre,filtre2)
    nsujet=i-1

    if (typeof == 0) then
        nz=nz0
        nz1=nz
        nz2=nz

        nzloco=nz1
        nzdc=nz2

        if(nz.gt.20)then
            nz = 20
        endif
        if(nz.lt.4)then
            nz = 4
        endif
    end if

    ndatemax=2*nsujet+sum(ic0) !! rajout
    allocate(date(ndatemax),datedc(ndatemax))
    !!print*,"ndatemax=",ndatemax,"sum(ic0)=",sum(ic0)
    !stop
    if(typeof == 0) then
        allocate(mm3(ndatemax),mm2(ndatemax) &
        ,mm1(ndatemax),mm(ndatemax),im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax))
    end if

!!!  DONNEES DECES

    mindc = 0.d0
    maxdc = maxtdc
    do i = 1,2*ngtemp
        do k = 1,ngtemp
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

    datedc(1) = aux(1)
    k = 1
    do i=2,2*ngtemp
        if(aux(i).gt.aux(i-1))then
            k = k+1
            datedc(k) = aux(i)
        endif
    end do
    
    if(typeof == 0) then
        ndatedc = k
    end if
!!         !write(*,*)'** ndatemax,ndatemaxdc',ndatemax,ndatemaxdc
!--------------- zi- ----------------------------------

!      construire vecteur zi (des noeuds)
!!! DONNEES RECURRENTES

    min = 0.d0
    aux = 0.d0
    max_ = maxt

    do i = 1,(2*nsujet+sum(ic0)) !! rajout
        do k = 1,nsujet
            if (t0(k).ge.min) then
                if (t0(k).lt.max_) then
                    max_ = t0(k)
                endif
            endif
            if (t1(k).ge.min) then
                if (t1(k).lt.max_) then
                    max_ = t1(k)
                endif
            endif
            if (tU(k).ne.t1(k)) then
                if (tU(k).ge.min) then
                    if (tU(k).lt.max_) then !! rajout
                        max_ = tU(k)
                    endif
                endif
            endif
        end do
        aux(i) = max_
        min = max_ + 1.d-12
        max_ = maxt
    end do

    date(1) = aux(1)
    k = 1
    do i=2,(2*nsujet+sum(ic0))
        if(aux(i).gt.aux(i-1))then
            k = k+1
            date(k) = aux(i)
        endif
    end do

    if(typeof == 0) then
        ndate = k
        nzmax=nz+3
            
        allocate(zi(-2:nzmax))
        zi(-2) =0
        zi(-2) = date(1)
        zi(-1) = date(1)
        zi(0) = date(1)
        zi(1) = date(1)

        h = (date(ndate)-date(1))/dble(nz-1)

        do i=2,nz-1
            zi(i) =zi(i-1) + h
        end do
        zi(nz) = date(ndate)
        zi(nz+1)=zi(nz)
        zi(nz+2)=zi(nz)
        zi(nz+3)=zi(nz)
        !!print*,"size(ziout)",size(ziout),"size(zi)",size(zi)
        !!print*,"ziout(1:3)",ziout(1:3)
        !!print*,"zi(1:3)",zi(1:3)
        ziOut = zi
! ajout TPS
    allocate(zidc(-2:(nzdc+3)))

!-------------------------------------------------------------------
    zidc(-2) = datedc(1) 
    zidc(-1)= datedc(1) 
    zidc(0)= datedc(1) 
    zidc(1)= datedc(1) 
!
    hdc =(datedc(ndatedc)-datedc(1))/dble(nzdc-1)
    do i=2,nzdc-1
        zidc(i) =zidc(i-1)+hdc
    end do 
!      hh=0
!     do i=2,nzdc-1
!      hh=int(i*ndatedc/nzdc)
!          zidc(i)=datedc(hh)
!     end do

    zidc(nzdc) = datedc(ndatedc)
    zidc(nzdc+1)=zidc(nzdc)
    zidc(nzdc+2)=zidc(nzdc)
    zidc(nzdc+3)=zidc(nzdc)
!
! fin ajout TPS

    end if

!---------- affectation nt0dc,nt1dc DECES ----------------------------

    indictronqdc=0
    do k=1,ngtemp
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
    !!print*,"datedc(",datedc
    !!print*,"nt1dc",nt1dc
    !stop

!---------- affectation nt0,nt1,ntU RECURRENTS----------------------------

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
                if(date(j).eq.tU(i))then !! rajout
                    ntU(i)=j
                endif
            end do
        end if
    end do
    
!!print*,"groupe=",groupe0
!stop
!     !print*,"event"
!     do i=20,40
!         !print*,i,t1(i),tU(i),c(i),ve(i,1)
!     enddo
!     !print*,"death"
!     do i=1,20
!         !print*,i,t1dc(i),cdc(i),ve(i,1)
!     enddo

    if (typeof == 0) then
!---------- affectation des vecteurs de splines -----------------
        n  = nz+2
!AD:add argument:ndatedc
        call vecspliJ(n,ndate,ndatedc)
!AD:end
        allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax) &
        ,m3m1(nzmax),m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax))

        call vecpenJ(n)
    end if

    npmax=np
    allocate(the1(-2:npmax),the2(-2:npmax))
    allocate(hess(npmax,npmax),I1_hess(npmax,npmax) &
    ,H1_hess(npmax,npmax),I2_hess(npmax,npmax),H2_hess(npmax,npmax),HI2(npmax,npmax) & 
    ,HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax),BIAIS(npmax,1))
    !Hspl_hess(npmax,npmax),PEN_deri(npmax,1),

!------- initialisation des parametres

    ! savoir si l'utilisateur a entre des parametres initiaux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (sum(b).eq.0.d0) then
        aaa=1
        !!print*,"npmax=",npmax
        do i=1,(npmax-nva)
            b(i)=5.d-1
        end do
        !b(npmax-nva+1)=-0.5d0
        !b(npmax-nva+2)=-0.25d0
        b(npmax-nva+1)=0.5d0
        b(npmax-nva+2)=0.5d0
    endif

    if(typeof ==1) then
         b(1:nbintervR) = 0.8d0!1.d-2!
!         b((nbintervR+1):(nbintervR+nbintervDC)) = 0.8d0!1.d-2
!         b(np-nva-indic_alpha)=5.d-1 ! pour theta
    end if

!    !write(*,*)'typeof',typeof

    if (typeof == 2) then
        b(1:4)=1.d-1!0.8d0
!        b(np-nva-indic_alpha)=5.d-1 ! pour theta
!        b(np-nva-indic_alpha)=1.d0 ! pour theta
    end if

    if ((typeof == 0) .and. (aaa.eq.1)) then
        !b(np-nva-indic_alpha)=1.d0 ! pour sigma2
        !b(np-nva)=1.d0 ! pour alpha
        b(np-nva-nparamfrail+1:np-nva)=1.d0
        !!print*,b
        !stop
    end if
    !!print*,b
    !!print*,"size(b)",size(b)
    !stop




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

        ttt(0)=mint !0.d0
        ttt(nbintervR)=cens

        j=0
        do j=1,nbintervR-1
            if (equidistant.eq.0) then
                ttt(j)=(t2(ent*j)+t2(ent*j+1))/(2.d0)
            else
                ttt(j)=(cens/nbintervR)*j !! ici c'est surement faux, car la troncature n'est pas prise en compte
            endif
        end do
        time = ttt
!----------> taille - nb de deces
        j=0
        do i=1,ngtemp
            if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
                j=j+1
            endif
        end do
        nbdeces=j

!----------> allocation des vecteur temps
        allocate(t3(nbdeces))
!----------> remplissage du vecteur de temps
        j=0
        do i=1,ngtemp
            if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1)then
                j=j+1
                t3(j)=t1dc(i)
            endif
        end do
!----------> tri du vecteur de temps
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
        tttdc(0)=mintdc !0.d0
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
        deallocate(t2,t3)
!------- FIN RECHERCHE DES NOEUDS
    end if

    ca=0.d0
    cb=0.d0
    dd=0.d0

!     allocate(knotsTPS(-qorder+1:nbinnerknots+qorder))
!     allocate(knotsdcTPS(-qorder+1:nbinnerknots+qorder))
!     allocate(innerknots(nbinnerknots))
!     allocate(innerknotsdc(nbinnerknots))
! 
! ! ajout TPS
!     call searchknotstps(t1,knotsTPS,nbinnerknots,qorder,nsujet,equidistantTPS,c,mint)
!     call searchknotstps(t1dc,knotsdcTPS,nbinnerknots,qorder,ng,equidistantTPS,cdc,mintdc) !! pas ng ?
! 
!     innerknots(1:nbinnerknots)=knotsTPS(1:nbinnerknots)
!     innerknotsdc(1:nbinnerknots)=knotsdcTPS(1:nbinnerknots)
!     boundaryknots(1)=knotsdcTPS(0)
!     boundaryknots(2)=knotsdcTPS(nbinnerknots+1)
! fin ajout TPS

    !if (typeof.ne.0)allocate(kkapa(2))
!!print*,"youhou"
!     if(typeJoint .ne. 1) then
! !====== INITIALISATION SUR SHARED Frailty MODEL
!         indic_joint=0
!         indic_alpha=0
!         nvatmp=nva
!         nva=nva1
!         select case(typeof)
!             case(0)
!                 npinit = nz+2+nva1+effet
!             case(1)
!                 npinit = nbintervR+nva1+effet
!             case(2)
!                 npinit = 2+nva1+effet
!         end select
! 
!         allocate(Binit(npinit))
!         if (typeof .ne. 0)allocate(vvv((npinit*(npinit+1)/2)))
!         if(typeof==1)allocate(betacoef(nbintervR))
! 
!         nst=1
!         stra=1
! !        select case(initialisation)
! !            case(1)
!                 !=======> initialisation par shared
!                 Binit=1.d-1 !5.d-1
!                 if(typeof==0) then
!                     Binit((nz+2+1):(nz+2+nva1))=1.d-1
!                     Binit(nz+2+nva1+effet)=1.d0
!                 end if
! !    !write(*,*),'===================================',effet,npinit
! !    !write(*,*),'== sur un SHARED FRAILTY model ====',(nz+2),nva1,nva2
! !    !write(*,*),'=== donnees recurrentes uniquement ='
! !    !write(*,*),'==================================='
! !     !write(*,*),(Binit(i),i=1,npinit)
! 
!                 allocate(I_hess(npinit,npinit),H_hess(npinit,npinit),v((npinit*(npinit+3)/2)))
! 
!                 select case(typeof)
!                     case(0)
!                         if (timedep.eq.0) then
!                             if (logNormal.eq.0) then
!                                 call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaGsplines)
!                             else
!                                 call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaGsplines_log)
!                             endif
!                         else
!                             call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaG_tps)
!                         endif
!                     case(1)
!                         if (timedep.eq.0) then
!                             call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaGcpm)
!                         else
!                             call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaG_tps)
!                         endif
!                     case(2)
!                         if (timedep.eq.0) then
!                             call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaGweib)
!                         else
!                             call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaG_tps)
!                         endif
!                 end select
!             deallocate(I_hess,H_hess,v)
! 
! !    !write(*,*)'======= Fin Shared ======='
! !    !write(*,*),(Binit(i),i=1,npinit)
! 
!             b(1:npinit) = Binit
! 
! !            case(2)
! !            !=======> initialisation par utilisateur
! !                do i=1,npinit
! !                    b(i) = Bshared(i)
! !                end do
! !            case(3)
! !            !=======> initialisation par defaut
! !            b(1:npinit)=1.d-5
! !        end select
! 
!         deallocate(Binit)
! 
!         do j=1,nva1
!             b(np-nva1+j)=b(nz+2+1+j)
!         end do
!         select case(typeof)
!             case(0)
!                 b((np-nva1+1):np)=b((nz+4):(nz+3+nva1))
!             case(1)
!                 b((np-nva1+1):np)=b((nz+4):(nz+3+nva1))
!             case(2)
!                 b((np-nva1+1):np)=b((nz+4):(nz+3+nva1))
!         end select
! 
!         do j=1,nva1
!             b(np-nva1-nva2+j)=b(np-nva2+j)
!         end do
! 
!         nva=nvatmp
! 
!         select case(typeof)
!             case(0)
!                 b(np-nva-indic_alpha)=b(nz+2+effet)
!             case(1)
!                 b(np-nva-indic_alpha)=b(nbintervR+effet)
!             case(2)
!             b(np-nva-indic_alpha)=b(2+effet)
!         end select
! 
!         b(np-nva)=1.d0         !2.2d0 ! pour alpha
! !=====initialisation des splines
! 
!         select case(typeof)
!             case(0)
!                 b((nz+3):(2*(nz+2)))=b(1:(nz+2))
!             case(1)
!                 b((nbintervR+1):2*(nbintervR))=b(1:nbintervR)
!             case(2)
!                 b(3:4)=b(1:2)
!         end select
! 
!         nst=2
!         indic_joint = 1
!         indic_alpha = 1
! 
!         if (typeof.ne.0) deallocate(vvv)
!         if(typeof==1)deallocate(betacoef)
! 
!     end if

    allocate(I_hess(np,np),H_hess(np,np),v((np*(np+3)/2)))

    if (typeof==1)allocate(betacoef(nbintervR+nbintervDC))
    if (typeof .ne. 0)allocate(vvv((np*(np+1)/2)))

    !if (istop .ne. 1)goto 1000 ! si l'initialisation ne marche pas, ne pas faire le modele

!     if(typeJoint == 1) then
                ! !print*,"suis dans joint:paramètres initiaux pour appel Marquard"
                ! !print*,"k0=",k0
                ! !print*,"b=",b
                ! !print*,"np=",np
                ! !print*,"ni=",ni
                ! !print*,"v=",v
                ! !print*,"res=",res
                ! !print*,"ier=",ier
                ! !print*,"istop=",istop
                ! !print*,"effet=",effet
                ! !print*,"ca,cb,dd =",ca,cb,dd
    if(adaptative .and. estim_wij_chap.eq.0)then
        nb0=nsim_nodes(7)
        nb_ree = nb0 + (nb0*(nb0-1))/2.d0 ! les variances + covariances estimes de la cholesky
        ! Parametres pour GH pseudo-adaptative
        allocate(ui_chap(ng0,1),invBi_cholDet(ng0),invBi_chol(ng0,nb_ree),invBi_cholDet_Essai(ntrials))
        allocate(invBi_chol_Essai(ntrials*9),invBi_chol_Individuel(ng0),ui_chap_Essai(ntrials,3))
                        
    endif
    
    !call intpr("type_joint =", -1, type_joint, 1)    
    call cpu_time(tp1)
    select case(typeof)
        case(0) ! fonction de risque de base approximee par des splines
            select case(type_joint)
                case(0)
                    call marq98j_SCL_0(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajsplines_surrogate_1)
                case(1)                
                    if(methodInt==3) then !integration par laplace
                        !!call MPI_ABORT(MPI_COMM_WORLD,erreur,code)! on stop tous les programmes appartenant au communicateur code, equivalent de l'instruction stop en sequantiel
                        !========= gestion du nombre d'essai a manipuler par processus dans le cas de laplace=========== 
                        n_par_pro=INT(ntrials/nb_procs)
                        suplement=ntrials-n_par_pro*nb_procs ! donne le nombre de simulation a partager entre les premiers processus seulement
                        
                        ! remplissage du table du nombre de simulation a effectuer par processus
                        !!print*,nb_procs,n_sim,suplement
                        allocate(table_par_pro(nb_procs))
                        table_par_pro(1:nb_procs)=n_par_pro
                        !!print*,tableNsim
                        table_par_pro(1:suplement)=n_par_pro+1 ! tous les essais jusqu'au rang supplement-1 recoivent une tâche supplementaire a realiser
                    endif
                    ! if(rang==0) then
                        call marq98j_SCL_0(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajsplines_surrogate)
                    ! endif
                    
                    if(methodInt==3) then !integration par laplace
                        deallocate(table_par_pro)
                    endif
                    
                    ! !call MPI_FINALIZE(code) 
                
                case(2)
                    call marq98j_SCL_0(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajsplines_surrogate)    
                case(3) ! the joint frailty-copula model
                    call marq98j_SCL_0(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajsplines_copule_surrogate)    
            endselect
        case(1) ! fonctions de risque de base supposees constantes par morceau
!                 if (timedep.eq.0) then
!                     if (logNormal.eq.0) then
!                        call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajcpm)
!                     else
!                         call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajcpm_log)
!                     endif
!                 else
!                     call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaj_tps)
!                 endif
        case(2) ! fonctions de risque de base approchees par Weibull
!                 if (timedep.eq.0) then
!                     if (logNormal.eq.0) then
!                         if (intcens.eq.1) then
!                             call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajweib_intcens)
!                         else
!                            call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajweib)
!                         endif
!                     else
!                         call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajweib_log)
!                     endif
!                 else
!                     call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaj_tps)
!                 endif
    end select
        call cpu_time(tp2)
!     end if

!if(rang==0) !print*,"============fin d'optimisation à l'aide de Marq98J==============="

    resOut=res
!Al:
    EPS(1) = ca
    EPS(2) = cb
    EPS(3) = dd
!Al:
    if (istop.ne.1) goto 1000

    call multiJ(I_hess,H_hess,np,np,np,IH) ! IH pour le calcul du LCV, I_hess= matrice hessienne , H_hess= matrice hessienne (variances-covariances des parametres estimes)
    call multiJ(H_hess,IH,np,np,np,HIH)
!!print*,"covar========================1"
    if(effet.eq.1.and.ier.eq.-1)then
        v((np-nva-indic_alpha)*(np-nva-indic_alpha+1)/2)=10.d10
    endif
    
        res01(effet+1)=res

! --------------  Lambda and survival estimates JRG January 05
    nstRec = 1 ! scl 30/10/2018 
    select case(typeof)
        case(0)
            call distanceJsplines(nz1,nz2,b,mt1,mt2,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
        case(1)
            Call distanceJcpm(b,nbintervR+nbintervDC,mt1,mt2,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
        case(2)
            Call distanceJweib(b,np,mt1,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
    end select
!!print*,"covar========================2"
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
    paraweib(1) = shapeweib(1)
    paraweib(2) = shapeweib(2)
    paraweib(3) = scaleweib(1)
    paraweib(4) = scaleweib(2)

    do ss=1,npmax
        do sss=1,npmax
            HIHOut(ss,sss) = HIH(ss,sss)
            H_hessOut(ss,sss)= H_hess(ss,sss)
        end do
    end do
    ! call dblepr("H_hessOut(17,17):", -1, H_hessOut(17,17), 1)

!!print*,"covar=======================3="    
!Ouverture du fichiers resultat

! open(5,file='OutJoint_Result_surrogate.txt')
if (istop.ne.1) then
    !write(*,*)"ERREUR : LE MODELE N'A PAS CONVERGE. Istop=",istop
    !write(5,*)"ERREUR : LE MODELE N'A PAS CONVERGE. Istop=",istop
else
    ! !write(*,*)"LE MODELE A BIEN CONVERGE"
    ! !write(*,*)" "
    ! !write(*,*)"***************************"
    ! !write(*,*)"       Resultats bruts"
    ! !write(*,*)"***************************"
    ! !write(*,*)" "
    ! !write(*,*)'nombre total de parametres',np
    ! !write(*,*)'kappa1',k0(1)
    ! !write(*,*)'kappa2',k0(2)
    ! !write(*,*)'nombre de noeuds:',nz         
    
    ! !write(*,*)'***Estimates des effets aléatoires:***'
    ! !write(*,*)" "
    ! !write(*,*)"=======================Niveau individuel=================="
    ! rangparam=np-nva-nparamfrail+indice_eta+indice_theta
    ! !write(*,*)'THETA^2 = variance de w_ij(w_ij suivent une log-normale)=',&
    ! b(rangparam)**2.d0
! !!print*,"covar========================"
    ! ! pour tenir compte du changt de var dans l'estimation de la variance (car on estime theta et pas theta^2), on utilise la delta methode pour lecart type
    ! !write(*,*)'SE theta (=H)',dsqrt(((2.d0*b(rangparam))**2)* &
    ! H_hessOut(rangparam,rangparam))
    ! !write(*,*)'SE theta (=HIH)',dsqrt(((2.d0*b(rangparam))**2)* &
    ! HIHOut(rangparam,rangparam))
    
    ! !write(*,*)" "
    ! !!write(*,*)'Parametres de fragilité, eta dans eta*w_ij, pour les deces=',&
    ! !b(np-nva-nparamfrail+indice_eta)
    ! !!write(*,*)'SE eta (=H)',dsqrt(H_hessOut(np-nva-nparamfrail+indice_eta,np-nva-nparamfrail+indice_eta))
    ! !!write(*,*)'SE eta (=HIH)',dsqrt(HIHOut(np-nva-nparamfrail+indice_eta,np-nva-nparamfrail+indice_eta))
    ! !write(*,*)" "
    ! !write(*,*)'*** Effets estimes du traitement ***'
    ! !write(*,*)" "
    ! ! Effets des covariables avec les intervalles de confiances
    ! do i=1,nva
        ! rangparam=np-nva+i
        ! !Intervalle de confiance
        ! bi = b(rangparam) - 1.96*dsqrt(H_hessOut(rangparam,rangparam))
        ! bs = b(rangparam) + 1.96*dsqrt(H_hessOut(rangparam,rangparam))
                
        ! if(i.eq.1)then
            ! !write(*,*)'*** For surrogate ***'
        ! endif
                
        ! if(i.eq.nva1+1)then
            ! !write(*,*)'*** For true endpoint ***',nva,i
        ! endif
                
        ! !write(*,*)'**************** '
        ! !write(*,*)'Variable traitement : '
                
        ! !write(*,*)i,')','beta=',b(rangparam)
        ! !write(*,*)' '
        ! !write(*,*)i,')',' SE (=H)',dsqrt(H_hessOut(rangparam,rangparam))
        ! wres=(b(rangparam))/dsqrt(H_hessOut(rangparam,rangparam))
        ! !write(*,*)'---> WALD',wres
        ! !write(*,*)' '
        ! !write(*,*)i,')',' SE (=HIH)',dsqrt(HIHOut(rangparam,rangparam))
        ! !write(*,*)'---> WALD',(b(rangparam))/dsqrt(HIHOut(rangparam,rangparam))
        ! !write(*,*)' '
        ! !write(*,*)'RR : ',dexp(b(rangparam)),'  IC',dexp(bi),dexp(bs) 
        ! !write(*,*)'**************** ' 
        ! !write(*,*)" "        
    ! end do 
            
    ! !write(*,*)'*** Parametres de convergence ***'    
    ! !write(*,*)" "
    ! !write(*,*)'valeur de ni(Nombre ditterations): ni=',ni
    ! !write(*,*)"log-vraisemblance:res=",res     
    ! !write(*,*)"critere sur les coefficients: ca=",ca
    ! !write(*,*)"critere sur la vraisemblance: cb=",cb
    ! !write(*,*)"critere sur le gradient: dd=",dd
    ! !write(*,*)"joint ligne 1344 pour exploiter les resultats du joint enlever le stop"
    ! !write(*,*)'************joint ligne 1221 temps_CPU mis pour lestimation des parametres (en secondes)******* =',tp2-tp1
    
    ! ecriture dans le fichier
    !write(5,*)"LE MODELE A BIEN CONVERGE"
    !write(5,*)" "
    !write(5,*)"***************************"
    !write(5,*)"       Resultats bruts"
    !write(5,*)"***************************"
    !write(5,*)" "
    !write(5,*)'nombre total de parametres',np
    !write(5,*)'kappa1',k0(1)
    !write(5,*)'kappa2',k0(2)
    !write(5,*)'nombre de noeuds:',nz         
    
    !write(5,*)'***Estimates des effets aleatoires:***'
    !write(5,*)" "
    !write(5,*)"=======================Niveau individuel=================="
    rangparam=np-nva-nparamfrail+indice_eta+indice_theta
    !write(5,*)'THETA^2 = variance de w_ij(w_ij suivent une log-normale)=',&
    !b(rangparam)**2.d0
    
    ! pour tenir compte du changt de var dans l'estimation de la variance (car on estime theta et pas theta^2), on utilise la delta methode pour lecart type
    !write(5,*)'SE theta (=H)',dsqrt(((2.d0*b(rangparam))**2)* &
    !H_hessOut(rangparam,rangparam))
    !write(5,*)'SE theta (=HIH)',dsqrt(((2.d0*b(rangparam))**2)* &
   ! HIHOut(rangparam,rangparam))
    
    !!write(5,*)'Paramètres de fragilité, eta dans eta*w_ij, pour les deces=',&
    !b(np-nva-nparamfrail+indice_eta)
    !!write(5,*)'SE eta (=H)',dsqrt(H_hessOut(np-nva-nparamfrail+indice_eta,np-nva-nparamfrail+indice_eta))
    !!write(5,*)'SE eta (=HIH)',dsqrt(HIHOut(np-nva-nparamfrail+indice_eta,np-nva-nparamfrail+indice_eta))
    !write(5,*)" "
    !write(5,*)'*** Effets estimes du traitement ***'
    !write(5,*)" "
    ! Effets des covariables avec les intervalles de confiances
    do i=1,nva
        rangparam=np-nva+i
        !Intervalle de confiance
        bi = b(rangparam) - 1.96*dsqrt(H_hessOut(rangparam,rangparam))
        bs = b(rangparam) + 1.96*dsqrt(H_hessOut(rangparam,rangparam))
                
        if(i.eq.1)then
            !write(5,*)'*** For surrogate ***'
        endif
                
        if(i.eq.nva1+1)then
            !write(5,*)'*** For true endpoint ***',nva,i
        endif
                
        !write(5,*)'**************** '
        !write(5,*)'Variable traitement : '
                
        !write(5,*)i,')','beta=',b(rangparam)
        !write(5,*)' '
        !write(5,*)i,')',' SE (=H)',dsqrt(H_hessOut(rangparam,rangparam))
        wres=(b(np-nva+i))/dsqrt(H_hessOut(rangparam,rangparam))
        !write(5,*)'---> WALD',wres
        !write(5,*)' '
        !write(5,*)i,')',' SE (=HIH)',dsqrt(HIHOut(rangparam,rangparam))
        !write(5,*)'---> WALD',(b(rangparam))/dsqrt(HIHOut(rangparam,rangparam))
        !write(5,*)' '
        !write(5,*)'RR : ',dexp(b(rangparam)),'  IC',dexp(bi),dexp(bs) 
        !write(5,*)'**************** '   
        !write(5,*)" "        
    end do 
    !write(5,*)" "
    !write(5,*)'*** Parametres de convergence ***'    
    !write(5,*)" "
    !write(5,*)'valeur de ni(Nombre dittérations): ni=',ni
    !write(5,*)"log-vraisemblance:res=",res     
    !write(5,*)"critere sur les coefficients: ca=",ca
    !write(5,*)"critere sur la vraisemblance: cb=",cb
    !write(5,*)"critere sur le gradient: dd=",dd
    !write(5,*)"joint ligne 1418 pour exploiter les resultats du joint enlever le stop"
    !write(5,*)'************joint ligne 1226 temps_CPU mis pour lestimation des parametres (en secondes)******* =',tp2-tp1
endif

!AD:add LCV
!LCV(1) The approximate like cross-validation Criterion
!LCV(2) Akaike information Criterion
!     calcul de la trace, pour le LCV (likelihood cross validation)
    LCV=0.d0
    if (typeof == 0) then
        !write(5,*)'The approximate like cross-validation Criterion in the non parametric case'
        call multiJ(H_hess,I_hess,np,np,np,HI)
        do i =1,np
            LCV(1) = LCV(1) + HI(i,i)
        end do     
        LCV(1) = (LCV(1) - resnonpen) / nsujet
        !write(*,*)'======== LCV :',LCV(1)
    else
        !write(5,*)'=========> Akaike information Criterion <========='
        LCV(2) = (1.d0 / nsujet) *(np - resOut)
        !write(5,*)'======== AIC :',LCV(2)
    end if

! close(5)
goto 12345

1000 continue
!AD:end
!!print*,"covar========================252"
    if (timedep.eq.1) then

        BetaTpsMat = 0.d0
        BetaTpsMatDc = 0.d0
    
! CALCUL DES ESTIMATIONS DES beta DEPENDANTS DU TEMPS

        if (nva1.gt.0) then
            call drawTimeCoef(np,b,nva1,filtretps,BetaTpsMat)
        end if

        if (nva2.gt.0) then
            call drawTimeCoefdc(np,b,nva2,filtre2tps,BetaTpsMatDc)
        end if

    endif

!    !write(*,*)'======== Calcul des residus de martingale ========'
!!print*,"covar========================12"
    deallocate(I_hess,H_hess)
!!print*,"covar========================13"
    coefBeta = 0.d0
    coefBetadc = 0.d0
    Xbeta = 0.d0
    Xbetadc = 0.d0
    XbetaG = 0.d0

    do i=1,nsujet
        do j=1,nva1
            ve1(i,j)=ve(i,j)
        end do
    end do
    do i=1,ngtemp
        do j=1,nva2
            ve2(i,j)=vedc(i,j)
        end do
    end do

    if (timedep.eq.0) then
        coefBeta(1,:) = b((np-nva+effet):(np-nva+nva1)) !b((np-nva+effet):(np-nva+effet+nva1))
        coefBetadc(1,:) = b((np-nva+effet+nva1):np) !b((np-nva+effet+nva1+1):np)
        Xbeta = matmul(coefBeta,transpose(ve1))
        if(typeJoint==1)then
            Xbetadc = matmul(coefBetadc,transpose(ve2))
        else
            XbetaG = matmul(coefBetadc,transpose(ve2))
        endif
    else
        do j=1,nsujet
            p=1
            call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,nbinnerknots+qorder,t1(j), &
            innerknots,boundaryknots,basis)
            do i=1,nva1
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
        do j=1,ngtemp
            p=1
            call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,nbinnerknots+qorder,t1dc(j), &
            innerknots,boundaryknots,basis)
            do i=1,nva2
                coefBeta2 = 0.d0
                if (filtre2tps(i).eq.1) then
                    do k=-qorder+1,nbinnerknots
                        coefBeta2 = coefBeta2 + b(np-(nva+npbetatps)+nva1+p-1+k+qorder)*basis(k+qorder)
                    end do
                else
                    coefBeta2 = b(np-(nva+npbetatps)+nva1+p)
                endif
                if (typeJoint==1) Xbetadc(1,j) = Xbetadc(1,j) + coefBeta2*dble(vedc(j,i))
                if (typeJoint==0) XbetaG(1,j) = XbetaG(1,j) + coefBeta2*dble(vedc(j,i))
                p=p+filtre2tps(i)*(nbinnerknots+qorder-1)+1
            end do
        enddo
    endif
!!print*,"covar========================"
    !deallocate(I_hess,H_hess)

    if((istop == 1) .and. (effet == 1)) then

!        if(typeJoint == 1) then
            allocate(vecuiRes(ng),vres((1*(1+3)/2)),I_hess(1,1),H_hess(1,1))
            effetres = effet

            if (logNormal.eq.0) then
                            !print*,"suis là 1 nsujet=",nsujet
                Call ResidusMartingalej(b,np,funcpajres,Resmartingale,Resmartingaledc,frailtypred,frailtyvar)
            else    
                !!print*,"suis là 2 nsujet=",b,np,size(b)            
                Call ResidusMartingalej(b,np,funcpajres_log,Resmartingale,Resmartingaledc,frailtypred,frailtyvar)
            endif

            if (istopres.eq.1) then
                do i=1,nsujet
                    if (logNormal.eq.0) then
                        linearpred(i)=Xbeta(1,i)+dlog(frailtypred(g(i)))
                    else
                        linearpred(i)=Xbeta(1,i)+frailtypred(g(i))
                    endif
                end do

                if (typeJoint.eq.1) then
                    do i=1,ng
                        if (logNormal.eq.0) then
                            linearpreddc(i)=Xbetadc(1,i)+alpha*dlog(frailtypred(i))
                        else
                            linearpreddc(i)=Xbetadc(1,i)+alpha*frailtypred(i)
                        endif
                    end do
                else
                    do i=1,lignedc
                        if (logNormal.eq.0) then
                            linearpredG(i)=XbetaG(1,i)+alpha*dlog(frailtypred(gsuj(i)))
                        else
                            linearpredG(i)=XbetaG(1,i)+alpha*frailtypred(gsuj(i))
                        endif
                    end do
                endif
            endif
            

!        else
!            allocate(vecuiRes(lignedc),vres((1*(1+3)/2)),I_hess(1,1),H_hess(1,1))
!            effetres = effet
!
!            Call ResidusMartingalej(b,np,funcpajres,Resmartingale,Resmartingaledc,frailtypred,frailtyvar)
!
!            if (istopres.eq.1) then
!                do i=1,nsujet
!                    linearpred(i)=Xbeta(1,i)+dlog(frailtypred(g(i)))
!                end do
!    
!                do i=1,ngtemp
!                    linearpreddc(i)=Xbetadc(1,i)+alpha*dlog(frailtypred(i))
!                end do
!            endif
!        end if

        MartinGales(:,1)=Resmartingale
        MartinGales(:,2)=Resmartingaledc
        MartinGales(:,3)=frailtypred
        MartinGales(:,4)=frailtyvar
!!print*,"covar========================251"
        deallocate(I_hess,H_hess,vres,vecuiRes)
!!print*,"covar========================252"
    else
        !deallocate(I_hess,H_hess)

! les 4 variables suivantes sont maintenant dans MartinGales
!         Resmartingale=0.d0
!         Resmartingaledc=0.d0
!         frailtypred=0.d0
!         frailtyvar=0.d0

        MartinGales=0.d0
        linearpred=0.d0
        linearpreddc=0.d0
        linearpredG=0.d0
    end if
!!print*,"covar========================3333"
!============================= (Audrey)
! Prediction : fait dans prediction.f90 (fonction a part entiere dans le package)
!=============================


!     if(typeJoint.ne.1)then
! !add for joint group
! ! ATTENTION le tau de Kendall est calcule
! ! pour une valeur des variables explicatives
! !'************ TAU DE KENDALL 1 *************'
!         if (typeof==0) then
!             if(nva.gt.0)then 
!                 expb1=0.d0
!                 expb2=0.d0
!                 do i=1,nva1
!                     expb1=expb1+b(np-nva+i)
!                 end do
!                 expb1=dexp(expb1)
! 
!                 do i=1,nva2
!                     expb2=expb2+b(np-nva2+i)
!                 end do
!                 expb2=dexp(expb2)
! 
!                 call gaulagKend1(sstmp)
! 
!                 kendall(1,1)=sstmp
!                 kendall(1,2)=4.d0*sstmp-1.d0
!             endif
! 
!             expb1=1.d0
!             expb2=1.d0
!             call gaulagKend1(sstmp)
!             kendall(2,1)=sstmp
!             kendall(2,2)=4.d0*sstmp-1.d0
! 
! !'************ TAU DE KENDALL 2 *************'
!             if(nva.gt.0)then 
!                 expb1=0.d0
!                 expb2=0.d0
!                 do i=1,nva1
!                     expb1=expb1+b(np-nva+i)
!                 end do 
!                 expb1=dexp(expb1)
! 
!                 do i=1,nva2
!                     expb2=expb2+b(np-nva2+i)
!                 end do
!                 expb2=dexp(expb2)
! 
!                 call gaulagKend1bis(sstmp)
!                 kendall(3,1)=sstmp
!                 kendall(3,2)=2.d0*sstmp-1.d0
!             endif
!             expb1=1.d0
!             expb2=1.d0
!             call gaulagKend1bis(sstmp)
!             kendall(4,1)=sstmp
!             kendall(4,2)=2.d0*sstmp-1.d0
!         else
!             kendall=0.d0
!         end if
!     endif

12345 continue
!!print*,"covar========================1"
    if (istop.eq.1) deallocate(I_hess,H_hess)
!!print*,"covar========================2"
deallocate(res2s_sujet,res2_dcs_sujet)
    deallocate(nig,cdc,t0dc,t1dc,aux1,aux2,res1,res4,res3,mi,t0,t1,tU,c,stra,g,resL,resU,res5)
    deallocate(aux,vax,vaxdc,ve,vedc)
    deallocate(hess,v,I1_hess,H1_hess,I2_hess,H2_hess)
    deallocate(HI2,HIH,IH,HI)
    deallocate(BIAIS)
    deallocate(datedc)
    !!print*,"date=",date
    deallocate(date)
    !deallocate(hess,v,I1_hess,H1_hess,I2_hess,H2_hess,HI2,HIH,IH,HI,BIAIS,date,datedc)
    deallocate(ResidusRec,Residusdc,Rrec,Nrec,Rdc,Ndc,vuu,ve1,ve2)
    deallocate(the1,the2,nigdc,gsuj)
    deallocate(filtretps,filtre2tps)
    deallocate(betatps,betatps2)
    deallocate(const_res1,const_aux1,nigs,cdcs,pourtrial,res2s,res2_dcs,nigts,cdcts,nsujeti)
    deallocate(delta,deltastar)
    deallocate(const_res4,const_res5)
    deallocate(xx1)
    deallocate(ww1,chol)
    if(nsim_nodes(3).ne.0) then
        deallocate(ui_chap,invBi_cholDet,invBi_chol,invBi_cholDet_Essai,invBi_chol_Essai,invBi_chol_Individuel,&
                    ui_chap_Essai)
    endif
    
    if(methodInt==3) then !integration par laplace
        deallocate(wij_chap,control_wij_chap)
    endif
    
!     deallocate(knotsTPS,knotsdcTPS,innerknots,innerknotsdc)
    if(methodInt.eq.0 .or. methodInt.eq.2 .or. methodInt.eq.4) then ! scl: on desalloue le vecteur des elements a simuler pour MC seulement si on fait du MC
       deallocate(Vect_sim_MC)
    endif
! !print*,"covar========================21"
    if (typeof == 0) then
        deallocate(nt0dc,nt1dc,nt0,nt1,ntU,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc, &
        mm3,mm2,mm1,mm,im3,im2,im1,im,zi,zidc,m3m3,m2m2,m1m1,mmm,&
        m3m2,m3m1,m3m,m2m1,m2m,m1m)
    end if
    
    if(methodInt==3) deallocate(IhessLaplace,H_hess_laplace,hess_laplace,vvv_laplace,b_i_laplace,v_i_laplace)

    if (typeof .ne. 0)deallocate(vvv) !,kkapa)
    
    if (typeof == 1) then
        deallocate(ttt,tttdc,betacoef)
    end if
    !if(control.eq.1) goto 100
    return
    end subroutine joint_surrogate


!========================== VECSPLI ==============================
!AD:add argument:ndatedc
    ! subroutine vecspliJ(n,ndate,ndatedc)
    ! use tailles

    ! use comon,only:date,datedc,zi,mm3,mm2,mm1,mm,im3,im2,im1,im &
    ! ,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc

    ! IMPLICIT NONE
    
    ! integer,intent(in)::n,ndate,ndatedc
    ! integer::i,j,k
    ! double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2

! ! ----------  calcul de u(ti) :  STRATE1 ---------------------------
   ! ! attention the(1)  sont en nz=1
       ! ! donc en ti on a the(i)
    ! j=0
    ! do i=1,ndate-1
        ! do k = 2,n-2
            ! if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
                ! j = k-1
            ! endif
        ! end do 
        ! ht = date(i)-zi(j)
        ! htm= date(i)-zi(j-1)
        ! h2t= date(i)-zi(j+2)
        ! ht2 = zi(j+1)-date(i)
        ! ht3 = zi(j+3)-date(i)
        ! hht = date(i)-zi(j-2)
        ! h = zi(j+1)-zi(j)
        ! hh= zi(j+1)-zi(j-1)
        ! h2= zi(j+2)-zi(j)
        ! h3= zi(j+3)-zi(j)
        ! h4= zi(j+4)-zi(j)
        ! h3m= zi(j+3)-zi(j-1)
        ! h2n=zi(j+2)-zi(j-1)
        ! hn= zi(j+1)-zi(j-2)
        ! hh3 = zi(j+1)-zi(j-3)
        ! hh2 = zi(j+2)-zi(j-2)
        ! mm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        ! mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
        ! *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        ! mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
        ! h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        ! mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        ! im3(i) = (0.25d0*(date(i)-zi(j-3))*mm3(i))+(0.25d0*hh2 &
        ! *mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
        ! im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0) &
            ! +(h4*mm(i)*0.25d0)
        ! im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
        ! im(i)  = ht*mm(i)*0.25d0

    ! end do
! ! AD: add for death
! ! ----------  calcul de u(ti) :  STRATE2 ---------------------------
   ! ! attention the(1)  sont en nz=1
       ! ! donc en ti on a the(i)

    ! do i=1,ndatedc-1
        ! do k = 2,n-2
            ! if ((datedc(i).ge.zi(k-1)).and.(datedc(i).lt.zi(k)))then
                ! j = k-1
            ! endif
        ! end do 
        ! ht = datedc(i)-zi(j)
        ! htm= datedc(i)-zi(j-1)
        ! h2t= datedc(i)-zi(j+2)
        ! ht2 = zi(j+1)-datedc(i)
        ! ht3 = zi(j+3)-datedc(i)
        ! hht = datedc(i)-zi(j-2)
        ! h = zi(j+1)-zi(j)
        ! hh= zi(j+1)-zi(j-1)
        ! h2= zi(j+2)-zi(j)
        ! h3= zi(j+3)-zi(j)
        ! h4= zi(j+4)-zi(j)
        ! h3m= zi(j+3)-zi(j-1)
        ! h2n=zi(j+2)-zi(j-1)
        ! hn= zi(j+1)-zi(j-2)
        ! hh3 = zi(j+1)-zi(j-3)
        ! hh2 = zi(j+2)-zi(j-2)
        ! mm3dc(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        ! mm2dc(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
        ! *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        ! mm1dc(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
        ! h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        ! mmdc(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        ! im3dc(i) = (0.25d0*(datedc(i)-zi(j-3))*mm3dc(i))+(0.25d0*hh2 &
        ! *mm2dc(i))+(0.25d0*h3m*mm1dc(i))+(0.25d0*h4*mmdc(i))
        ! im2dc(i) = (0.25d0*hht*mm2dc(i))+(h3m*mm1dc(i)*0.25d0) &
            ! +(h4*mmdc(i)*0.25d0)
        ! im1dc(i) = (htm*mm1dc(i)*0.25d0)+(h4*mmdc(i)*0.25d0)
        ! imdc(i)  = ht*mmdc(i)*0.25d0

    ! end do
! ! AD:end        
    ! end subroutine vecspliJ

! !========================== VECPEN ==============================
    ! subroutine vecpenJ(n) 
    
    ! use tailles
    
    ! use comon,only:date,datedc,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
    
    ! IMPLICIT NONE
    
    ! integer,intent(in)::n
    ! integer::i
    ! double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2,a3,a2,b2 &
    ! ,c2,a1,b1,c1,a0,x3,x2,x


! !*********************************************************************
         
    ! do i=1,n-3
        ! h = zi(i+1)-zi(i)
        
        ! hh= zi(i+1)-zi(i-1)
        ! h2= zi(i+2)-zi(i)
        ! h3= zi(i+3)-zi(i)
        ! h4= zi(i+4)-zi(i)
        ! h3m= zi(i+3)-zi(i-1)
        ! h2n=zi(i+2)-zi(i-1)
        ! hn= zi(i+1)-zi(i-2)
        ! hh3 = zi(i+1)-zi(i-3)
        ! hh2 = zi(i+2)-zi(i-2)
        ! a3 = h*hh*hn*hh3
        ! a2 = hh2*hh*h*hn
        ! b2 = hh2*h2n*hh*h
        ! c2 = hh2*h2*h*h2n
        ! a1 = h3m*h2n*hh*h
        ! b1 = h3m*h2*h*h2n
        ! c1 = h3m*h3*h2*h
        ! a0 = h4*h3*h2*h
        ! x3 = zi(i+1)*zi(i+1)*zi(i+1)-zi(i)*zi(i)*zi(i)
        ! x2 = zi(i+1)*zi(i+1)-zi(i)*zi(i)
        ! x  = zi(i+1)-zi(i)
        
        ! m3m3(i) = (192.d0*h/(hh*hn*hh3*hh*hn*hh3))
        ! m2m2(i) = 64.d0*(((3.d0*x3-(3.d0*x2*(2.d0*zi(i+1)+zi(i-2) &
        ! ))+x*(4.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i-2)+4.d0*zi(i+1) &
        ! *zi(i-2)))/(a2*a2)))
        ! m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2)  &
        ! +zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1) &
        ! +zi(i+1)*zi(i+1)+2.d0*zi(i+2)*zi(i-1)+2.d0*zi(i+2) &
        ! *zi(i+1)+2.d0*zi(i-1)*zi(i+1)))/(b2*b2)))
        ! m2m2(i) = m2m2(i) +64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i+2) &
        ! +zi(i)))+x*(4.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i)+4.d0*zi(i+2) &
        ! *zi(i)))/(c2*c2))
        
        ! m2m2(i) = m2m2(i) +128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+2) &
        ! +zi(i-1)+3.d0*zi(i+1)+zi(i-2)))+x*(2.d0*zi(i+1)*zi(i+2) &
        ! +2.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i+2) &
        ! +zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*b2))
        ! m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0* & 
        ! x2*(2.d0*zi(i+2)+zi(i)+2.d0*zi(i+1)+zi(i-2)))+x* &
        ! (4.d0*zi(i+1)*zi(i+2)+2.d0*zi(i+1)*zi(i)+2.d0*zi(i-2) &
        ! *zi(i+2)+zi(i-2)*zi(i)))/(a2*c2))
        ! m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0*x2 &
        ! *(3.d0*zi(i+2)+zi(i)+zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i)+ &
        ! 2.d0*zi(i-1)*zi(i+2)+zi(i)*zi(i-1)+2.d0*zi(i+1)*zi(i+2) &
        ! +zi(i+1)*zi(i)+2.d0*zi(i+2)*zi(i+2)))/(b2*c2))
        ! m1m1(i) = 64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i-1)+zi(i+1))) &
        ! +x*(4.d0*zi(i-1)*zi(i-1)+zi(i+1)*zi(i+1)+4.d0*zi(i-1) &
        ! *zi(i+1)))/(a1*a1))
        ! m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i) &     
        ! +zi(i+2)))+x*(zi(i-1)*zi(i-1)+zi(i)*zi(i)+zi(i+2)* &
        ! zi(i+2)+2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+2.d0* &
        ! zi(i)*zi(i+2)))/(b1*b1))
        ! m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i+3) &
        ! +2.d0*zi(i)))+x*(zi(i+3)*zi(i+3)+4.d0*zi(i)*zi(i) &
        ! +4.d0*zi(i+3)*zi(i)))/(c1*c1)) 
        ! m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(3.d0 &
        ! *zi(i-1)+zi(i)+zi(i+2)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i-1) &
        ! +2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1) &
        ! +zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(a1*b1))
        ! m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+ &
        ! 2.d0*zi(i)+2.d0*zi(i-1)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i+3) &
        ! +4.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i))) &
        ! /(a1*c1))    
        ! m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+3.d0 &
        ! *zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1) &    
        ! *zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3) &
        ! +2.d0*zi(i+2)*zi(i)))/(b1*c1))
        ! mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
        ! m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2) &
        ! ))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2)) &
        ! +((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x* &
        ! (zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2)) &
        ! +((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x* &
        ! (2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)) )
        ! m3m1(i) = 192.d0*(((x3-(0.5d0*x2*(4.d0*zi(i+1)+2.d0*zi(i-1) &
        ! ))+x*(2.d0*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*a1)) &
        ! +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+2)+zi(i-1)+zi(i))) &
        ! +x*(zi(i+1)*zi(i-1)+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(b1*a3)) &
        ! +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+3)+2.d0*zi(i)))+x*(zi(i+1) &
        ! *zi(i+3)+2.d0*zi(i+1)*zi(i)))/(c1*a3)) )
        ! m3m(i) = 576.d0*((-(x3/3.d0)+(0.5d0*x2*(zi(i+1)+zi(i))) &
        ! -x*zi(i+1)*zi(i))/(a3*a0))
        ! m2m1(i) = 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1)+3.d0* &
        ! zi(i+1)+zi(i-2)))-x*(4.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1) &
        ! *zi(i+1)+2.d0*zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*a1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1)+ &
        ! zi(i)+zi(i+2)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i-1) &
        ! +2.d0*zi(i+1)*zi(i)+2.d0*zi(i+1)*zi(i+2)+zi(i-2)*zi(i-1)+ &
        ! zi(i-2)*zi(i)+zi(i-2)*zi(i+2)))/(a2*b1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)+2.d0 &
        ! *zi(i)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i+3)+4.d0 &
        ! *zi(i+1)*zi(i)+zi(i-2)*zi(i+3)+2.d0*zi(i-2)*zi(i)))/(a2*c1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2* &
        ! (3.d0*zi(i-1)+2.d0*zi(i+1)+zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1) &
        ! +zi(i+2)*zi(i+1)+2.d0*zi(i-1)*zi(i-1)+3.d0 &
        ! *zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(b2*a1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0 &
        ! *zi(i-1)+zi(i)+2.d0*zi(i+2)+zi(i+1)))-x*(zi(i+2)*zi(i-1) &
        ! +zi(i+2)*zi(i)+zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1)+zi(i-1) &
        ! *zi(i)+zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i) &
        ! +zi(i+1)*zi(i+2)))/(b2*b1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
        ! +2.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1)))-x*(zi(i+2)*zi(i+3) &
        ! +2.d0*zi(i+2)*zi(i)+zi(i-1)*zi(i+3)+2.d0*zi(i-1)*zi(i) &
        ! +zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i)))/(b2*c1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1) &
        ! +zi(i+1)+2.d0*zi(i+2)+zi(i)))-x*(4.d0*zi(i+2)*zi(i-1)+2.d0* &
        ! zi(i+2)*zi(i+1)+2.d0*zi(i)*zi(i-1)+zi(i)*zi(i+1)))/(c2*a1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1) &
        ! +2.d0*zi(i)+3.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1)+2.d0 &
        ! *zi(i+2)*zi(i)+2.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i-1)+zi(i) &
        ! *zi(i)+zi(i)*zi(i+2)))/(c2*b1))
        ! m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
        ! +3.d0*zi(i)+2.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i+3)+4.d0 &
        ! *zi(i+2)*zi(i)+zi(i)*zi(i+3)+2.d0*zi(i)*zi(i)))/(c2*c1))
        ! m2m(i) = 192.d0*(((x3-(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i+1) &
        ! +zi(i-2)))+x*(2.d0*zi(i+1)*zi(i)+zi(i-2)*zi(i)))/(a2*a0)) &
        ! +((x3-(0.5d0*x2*(3.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1))) &
        ! +x*(zi(i+2)*zi(i)+zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(b2*a0)) &
        ! +((x3-(0.5d0*x2*(4.d0*zi(i)+2.d0*zi(i+2)))+x*(2.d0*zi(i+2) &
        ! *zi(i)+zi(i)*zi(i)))/(c2*a0)) )
        ! m1m(i) = 192.d0*(((-x3+(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i-1) &
        ! +zi(i+1)))-x*(2.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(a1*a0)) &
        ! +((-x3+(0.5d0*x2*(4.d0*zi(i)+zi(i-1)+zi(i+2))) &
        ! -x*(zi(i-1)*zi(i)+zi(i)*zi(i)+zi(i+2)*zi(i)))/(b1*a0)) &
        ! +((-x3+(0.5d0*x2*(5.d0*zi(i)+zi(i+3)))-x*(zi(i+3)*zi(i) &
        ! +2.d0*zi(i)*zi(i)))/(c1*a0)) )

    ! end do

    ! end subroutine vecpenJ






! !==========================  SUSP  ====================================
    ! subroutine suspJ(x,the,n,su,lam,zi)
    
    ! use tailles
    
    ! IMPLICIT NONE 
    
    ! integer,intent(in)::n
    ! double precision,intent(out)::lam,su
    ! double precision,dimension(-2:npmax),intent(in)::zi,the
    ! double precision,intent(in)::x
    ! integer::j,k,i
    ! double precision::ht,ht2,h2,som,htm,h2t,h3,h2n,hn, &
    ! im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2 &
    ! ,h,gl,hh

    ! gl=0.d0
    ! som = 0.d0
    ! do k = 2,n+1
        ! if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            ! j = k-1
            ! if (j.gt.1)then
                ! do i=2,j
                    ! som = som+the(i-4)
                ! end do  
            ! endif   
            ! ht = x-zi(j)
            ! htm= x-zi(j-1)
            ! h2t= x-zi(j+2)
            ! ht2 = zi(j+1)-x
            ! ht3 = zi(j+3)-x
            ! hht = x-zi(j-2)
            ! h = zi(j+1)-zi(j)
            ! hh= zi(j+1)-zi(j-1)
            ! h2= zi(j+2)-zi(j)
            ! h3= zi(j+3)-zi(j)
            ! h4= zi(j+4)-zi(j)
            ! h3m= zi(j+3)-zi(j-1)
            ! h2n=zi(j+2)-zi(j-1)
            ! hn= zi(j+1)-zi(j-2)
            ! hh3 = zi(j+1)-zi(j-3)
            ! hh2 = zi(j+2)-zi(j-2)
            ! mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            ! mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            ! *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            ! mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            ! h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            ! mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            ! im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            ! +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            ! im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            ! im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            ! im  = ht*mm*0.25d0
            ! gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            ! lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        ! endif
    ! end do
   
    ! if(x.ge.zi(n))then
        ! som = 0.d0
        ! do i=1,n+1
            ! som = som+the(i-3)
        ! end do
        ! gl = som
    ! endif
    
    ! su  = dexp(-gl)
    
    ! return
     
    ! end subroutine suspJ

! !==========================  COSP  ====================================
! ! calcul les points pour les fonctions 
! ! et leur bandes de confiance

    ! subroutine cospJ(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)
    
    ! use tailles
    
    ! IMPLICIT NONE
    
    ! integer,intent(in)::n
    ! double precision,intent(in)::x
    ! double precision,intent(out)::lam,su
    ! double precision,intent(out)::binf,bsup,lbinf,lbsup
    ! double precision,dimension(npmax,npmax),intent(in)::y
    ! double precision,dimension(-2:npmax),intent(in)::the,zi
    ! integer::j,k,i
    ! double precision::ht,ht2,h2,som,pm,htm,h2t,h3,h2n,hn, &
    ! im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2, &
    ! h,gl,hh
    
    ! j=0
    ! gl=0.d0
    ! som = 0.d0
    ! do k = 2,n-1
        ! if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
            ! j = k-1
            ! if (j.gt.1)then
                ! do i=2,j
                ! som = som+the(i-4)
                ! end do  
            ! endif   
            ! ht = x-zi(j)
            ! htm= x-zi(j-1)
            ! h2t= x-zi(j+2)
            ! ht2 = zi(j+1)-x
            ! ht3 = zi(j+3)-x
            ! hht = x-zi(j-2)
            ! h = zi(j+1)-zi(j)
            ! hh= zi(j+1)-zi(j-1)
            ! h2= zi(j+2)-zi(j)
            ! h3= zi(j+3)-zi(j)
            ! h4= zi(j+4)-zi(j)
            ! h3m= zi(j+3)-zi(j-1)
            ! h2n=zi(j+2)-zi(j-1)
            ! hn= zi(j+1)-zi(j-2)
            ! hh3 = zi(j+1)-zi(j-3)
            ! hh2 = zi(j+2)-zi(j-2)
            ! mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            ! mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            ! *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            ! mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            ! h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            ! mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            ! im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            ! +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            ! im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            ! im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            ! im  = ht*mm*0.25d0
            ! gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
            ! lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
        ! endif
    ! end do
   
    ! if(x.ge.zi(n))then
        ! som = 0.d0
        ! do i=1,n
            ! som = som+the(i-3)
        ! end do
        ! gl = som
    ! endif

    ! call confJ(x,j,n,y,pm,zi)

    ! binf = dexp(-gl + 1.96d0*pm)
    ! su  = dexp(-gl)
    ! bsup = dexp(-gl - 1.96d0*pm)

    ! call conf1J(x,j,n,y,pm,zi)
    ! lbinf = lam - 1.96d0*pm
    ! lbsup = lam + 1.96d0*pm
! !         !write(*,*)'lbinf apres conf1',lbinf,lam,pm

    ! return

    ! end subroutine cospJ
     
     
! !=====================  CONF1  =============================


    ! subroutine  conf1J(x,ni,n,y,pm,zi)
    
    ! use tailles
    
    ! IMPLICIT NONE  
    
    ! integer,intent(in)::ni,n
    ! double precision,intent(in)::x
    ! double precision,dimension(-2:npmax),intent(in)::zi
    ! double precision,dimension(npmax,npmax),intent(in)::y
    ! double precision,intent(out)::pm
    ! integer::i,j
    ! double precision::res,mmspJ
    ! double precision,dimension(npmax)::vecti,aux

      
           
    ! do i=1,n
        ! vecti(i) = mmspJ(x,ni,i,zi)
    ! end do
    
    ! do i=1,n
        ! aux(i) = 0.d0
        ! do j=1,n
            ! aux(i) = aux(i) - y(i,j)*vecti(j)
        ! end do
    ! end do 


    ! res = 0.d0
    ! do i=1,n
        ! res = res + aux(i)*vecti(i)
    ! end do
    
    ! res=-res 
    ! pm = dsqrt(res)
    
    ! end subroutine  conf1J
     
! !=====================  CONF  =============================

    ! subroutine  confJ(x,ni,n,y,pm,zi)
    
    ! use tailles
    
    ! IMPLICIT NONE  
    
    ! integer,intent(in)::ni,n
    ! double precision,intent(in)::x
    ! double precision,dimension(-2:npmax),intent(in)::zi
    ! double precision,dimension(npmax,npmax),intent(in)::y
    ! double precision,intent(out)::pm
    ! integer::i,j
    ! double precision::res,ispJ
    ! double precision,dimension(52)::vecti,aux
    
    ! do i=1,n
    ! vecti(i) = ispJ(x,ni,i,zi)
    ! end do   

    ! do i=1,n
    ! aux(i) = 0.d0
    ! do j=1,n
        ! aux(i) = aux(i) - y(i,j)*vecti(j)
    ! end do
    ! end do   

    ! res = 0.d0
    ! do i=1,n
    ! res = res + aux(i)*vecti(i)
    ! end do
    ! res=-res
    ! pm = dsqrt(res)
               
    ! end subroutine  confJ


! !==========================   ISP   ==================================

    ! double precision function ispJ(x,ni,ns,zi)
    
    ! use tailles
    
    ! IMPLICIT NONE  
    
    ! integer,intent(in)::ni,ns
    ! double precision,intent(in)::x
    ! double precision,dimension(-2:npmax),intent(in)::zi
    ! double precision::val,mmspJ



    ! if(x.eq.zi(ni))then
        ! if(ni.le.ns-3)then
            ! val = 0.d0
            ! else
                ! if(ni.le.ns-2)then
                    ! val = ((zi(ni)-zi(ni-1))*mmspJ(x,ni,ns,zi))*0.25d0
                ! else
                    ! if (ni.eq.ns-1)then
                        ! val = ((zi(ni)-zi(ni-2))*mmspJ(x,ni,ns,zi)+ &
                        ! (zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+1,zi))*0.25d0
                    ! else
                        ! if(ni.eq.ns)then
                            ! val = ((zi(ni)-zi(ni-3))*mmspJ(x,ni,ns,zi)+ &
                            ! (zi(ni+2)-zi(ni-2))*mmspJ(x,ni,ns+1,zi) &
                            ! +(zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+2,zi))*0.25d0
                            ! else
                                ! val = 1.d0
                            ! endif
                        ! endif
                ! endif   
        ! endif
    ! else   
        ! if(ni.lt.ns-3)then
            ! val = 0.d0
        ! else
            ! if(ni.eq.ns-3)then
                ! val = (x-zi(ni))*mmspJ(x,ni,ns,zi)*0.25d0
            ! else  
                ! if(ni.eq.ns-2)then
                    ! val = ((x-zi(ni-1))*mmspJ(x,ni,ns,zi)+ &
                    ! (zi(ni+4)-zi(ni))*mmspJ(x,ni,ns+1,zi))*0.25d0
                ! else   
                    ! if (ni.eq.ns-1)then
                        ! val =((x-zi(ni-2))*mmspJ(x,ni,ns,zi)+ &
                        ! (zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+1,zi) &
                        ! +(zi(ni+4)-zi(ni))*mmspJ(x,ni,ns+2,zi))*0.25d0
                    ! else
                        ! if(ni.eq.ns)then
                            ! val =((x-zi(ni-3))*mmspJ(x,ni,ns,zi)+ &
                            ! (zi(ni+2)-zi(ni-2))*mmspJ(x,ni,ns+1,zi) &
                            ! +(zi(ni+3)-zi(ni-1))*mmspJ(x,ni,ns+2,zi) &
                            ! +(zi(ni+4)-zi(ni))*mmspJ(x,ni,ns+3,zi))*0.25d0
                        ! else
                            ! val = 1.d0
                        ! endif
                    ! endif
                ! endif
            ! endif
        ! endif 
    ! endif
      
    ! ispJ = val
    
    ! return
    
    ! end function ispJ
    
! !==========================  MMSP   ==================================

    ! double precision function mmspJ(x,ni,ns,zi)
    
    ! use tailles
    
    ! IMPLICIT NONE 
    
    ! integer,intent(in)::ni,ns
    ! double precision,intent(in)::x
    ! double precision,dimension(-2:npmax),intent(in)::zi
    ! double precision::val

    ! if(ni.lt.ns-3)then
        ! val = 0.d0
    ! else
        ! if(ns-3.eq.ni)then
            ! if(x.eq.zi(ni))then
                ! val = 0.d0
            ! else  
                ! val = (4.d0*(x-zi(ni))*(x-zi(ni)) &
                ! *(x-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3) &
                ! -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
            ! endif
        ! else 
            ! if(ns-2.eq.ni)then
                ! if(x.eq.zi(ni))then
                    ! val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1))) &
                    ! /((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1)) &
                    ! *(zi(ni+1)-zi(ni-1)))
                ! else  
                    ! val = (4.d0*(x-zi(ni-1))*(x-zi(ni-1)) &
                    ! *(zi(ni+1)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    ! -zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni))) &
                    ! +   (4.d0*(x-zi(ni-1))*(x-zi(ni)) &
                    ! *(zi(ni+2)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                    ! -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) &
                    ! +   (4.d0*(x-zi(ni))*(x-zi(ni)) &
                    ! *(zi(ni+3)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+3) &
                    ! -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                ! endif
            ! else   
                ! if (ns-1.eq.ni)then
                    ! if(x.eq.zi(ni))then
                        ! val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1) &
                        ! -zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1) &
                        ! -zi(ni-1))*(zi(ni+1)-zi(ni-2)))) &
                        ! +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni))) &
                        ! /((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1)) &
                        ! *(zi(ni+1)-zi(ni-1)))))
                    ! else
                        ! val = (4.d0*((x-zi(ni-2))*(zi(ni+1) &
                        ! -x)*(zi(ni+1)-x))/((zi(ni+2) &
                        ! -zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)- &
                        ! zi(ni))*(zi(ni+1)-zi(ni-2)))) &
                        ! +((4.d0*((x-zi(ni-1))*(zi(ni+2)-x)  &
                        ! *(zi(ni+1)-x))/((zi(ni+2)-zi(ni-2)) &
                        ! *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
                        ! (zi(ni+1)-zi(ni))))) &
                        ! +((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x) &
                        ! *(x-zi(ni)))/((zi(ni+2)-zi(ni-2)) &
                        ! *(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))* &
                        ! (zi(ni+1)-zi(ni)))))
                    ! endif 
                ! else
                    ! if(ni.eq.ns)then
                            ! if(x.eq.zi(ni))then
                            ! val =(4.d0*(x-zi(ni+1))*(x &
                            ! -zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1) &
                            ! -zi(ni-2))*(zi(ni+1)-zi(ni-3))))
                        ! else   
                            ! val =(4.d0*(x-zi(ni+1))*(x &
                            ! -zi(ni+1))*(zi(ni+1)-x)/((zi(ni+1) &
                            ! -zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1) &
                            ! -zi(ni))*(zi(ni+1)-zi(ni-3))))
                        ! endif
                    ! else
                        ! val = 0.d0
                    ! endif
                ! endif
            ! endif
        ! endif
    ! endif

    ! mmspJ = val
    
    ! return
    
    ! end function mmspJ

! !================== multiplication de matrice  ==================

! ! multiplie A par B avec le resultat dans C

    ! subroutine multiJ(A,B,IrowA,JcolA,JcolB,C)
! !     remarque :  jcolA=IrowB
    ! use tailles
    
    ! IMPLICIT NONE
    
    ! integer,intent(in)::IrowA,JcolA,JcolB
    ! double precision,dimension(IrowA,JcolA),intent(in):: A
    ! double precision,dimension(JcolA,JcolB),intent(in):: B
    ! double precision,dimension(IrowA,JcolB),intent(out)::C       
    ! integer::i,j,k
    ! double precision::sum
    
    ! !!print*,"A=",A
    ! !!print*,"B=",B
    ! do I=1,IrowA
        ! do J=1,JcolB
            ! sum=0
            ! do K=1,JcolA
                ! sum=sum+A(I,K)*B(k,J)
            ! end do
            ! C(I,J)=sum
            ! !!print*,"sum=",sum
        ! end do
    ! end do
    
    ! return
    
    ! end subroutine multiJ
! !====================================================================
! !============================    GAMMA      ==============================

! !       function qui calcule le log de  Gamma
    ! double precision function logGammaJ(xx)
    
    ! use donnees,only:cof,stp,half,one,fpf 
    
    ! implicit none
    
    ! integer::j
    ! double precision,intent(in)::xx
    
    
    ! double precision::x,tmp,ser
    
    ! x = xx - one
    ! tmp = x + fpf
    ! tmp = (x+half)*dlog(tmp) - tmp
    ! ser = one
    ! do j = 1,6
        ! x = x + one
        ! ser = ser + cof(j)/x
    ! end do
    ! logGammaJ = tmp + dlog(stp*ser)
    
    ! return
    
    ! end function logGammaJ


! !==================================================================
! ! calcul de l integrant, pour un effet aleatoire donne frail et un groupe donne auxig (cf funcpa)
! !==================================================================

    ! double precision function func1J(frail)

    ! use tailles
    ! use comon,only:auxig,alpha,theta,aux1,aux1,mi

    ! IMPLICIT NONE

    ! double precision,intent(in)::frail

    ! func1J = (frail**(alpha*mi(auxig)+1./theta-1.))* &
    ! dexp(-(frail**alpha) *aux1(auxig))*dexp(-frail/theta)
    
    ! return
    
    ! end function func1J

! !==================================================================

    ! double precision function func2J(frail)

    ! use tailles
    ! use comon,only:auxig,ALPHA,theta,aux2

    ! IMPLICIT NONE

    ! double precision,intent(in)::frail
    ! double precision::logGammaJ

    ! func2J = dexp(-(frail**alpha)*aux2(auxig))*dexp(-frail/theta)*(frail) &
    ! /(exp(logGammaJ(1.d0/theta))*(theta**(1./theta)))

    ! return

    ! end function func2J

! !==================================================================

    ! double precision function func3J(frail)

    ! use tailles
    ! use comon,only:nig,auxig,alpha,theta,res1,res3,aux1,cdc

    ! IMPLICIT NONE

    ! double precision,intent(in)::frail

    ! func3J = (nig(auxig)+ alpha*cdc(auxig)+ 1./theta-1.)*dlog(frail) &
    ! - frail*(res1(auxig)-res3(auxig)) & !res3=0 si AG=0
    ! - (frail**alpha)*(aux1(auxig))- frail/theta

    ! func3J = exp(func3J)

    ! return

    ! end function func3J

! !==================================================================

    ! double precision function func3bis(frail)

    ! use tailles
    ! use comon,only:g,nig,auxig,aux1,aux2,res1,res3,t0,t1,t0dc,&
    ! t1dc,c,cdc,nt0,nt1,nt0dc,nt1dc,nsujet,nva,nva1,nva2,ndate,ndatedc,nst,AG,&
    ! alpha,theta
    ! use comongroup,only:nigdc,gsuj
    
    ! implicit none

    ! double precision::frail

    ! func3bis = (nig(auxig)+ alpha*nigdc(auxig)+ 1.d0/theta-1.d0)*dlog(frail) &
    ! - frail*(res1(auxig)-res3(auxig))-(frail**alpha)*(aux1(auxig)) &
    ! - frail/theta

    ! func3bis = exp(func3bis)

    ! return

    ! end function func3bis

! !==================================================================
! ! gauss laguerre
! ! func1 est l integrant, ss le resultat de l integrale sur 0 ,  +infty

    ! subroutine gaulagJ(ss,choix)

    ! use tailles
    ! use comon,only:auxig,typeof,typeJoint
    ! use donnees,only:w,x,w1,x1
    
    ! implicit none

    ! integer,intent(in)::choix
    ! double precision,intent(out):: ss
    ! double precision ::auxfunca,func1J,func2J,func3J,func3bis
    ! external :: func1J,func2J,func3J,func3bis
    ! integer :: j

    ! ss=0.d0
! !    if (typeof == 0)then
! ! Will be twice the average value of the function,since the ten
! ! wei hts (five numbers above each used twice) sum to 2.
        ! do j=1,20
            ! if (choix.eq.1) then !integrale 1
                ! auxfunca=func1J(x(j))
                ! ss = ss+w(j)*(auxfunca)
            ! else                   !choix=2, survie marginale, vraie troncature
                ! if (choix.eq.2) then 
                    ! auxfunca=func2J(x(j))
                    ! ss = ss+w(j)*(auxfunca)
                ! else                   !choix=3, AG model
                    ! if (choix.eq.3) then
                        ! if(typeJoint==1)then
                            ! auxfunca=func3J(x(j))
                        ! else
                            ! auxfunca=func3bis(x(j))
                        ! endif
                        ! ss = ss+w(j)*(auxfunca)
                    ! endif
                ! endif
            ! endif
        ! end do
! !     else
! !         do j=1,32
! !             if (choix.eq.1) then !integrale 1
! !                 auxfunca=func1j(x1(j))
! !                 ss = ss+w1(j)*(auxfunca)
! !             else                   !choix=2, survie marginale, vraie troncature
! !                 if (choix.eq.2) then 
! !                     auxfunca=func2j(x1(j))
! !                     ss = ss+w1(j)*(auxfunca)
! !                 else                   !choix=3, AG model
! !                     if (choix.eq.3) then
! !                         if(typeJoint==1)then
! !                             auxfunca=func3j(x1(j))
! !                         else
! !                             auxfunca=func3bis(x1(j))
! !                         endif
! !                         ss = ss+w1(j)*(auxfunca)
! !                         !print*,ss
! !                     endif
! !                 endif
! !             endif
! !         end do
! !     end if

    ! return

    ! end subroutine gaulagJ

! !==================================================================

    ! subroutine gaulagJ_intcens(ss,choix)
    
    ! use tailles
    ! use comon,only:auxig,typeof
    ! use donnees,only:w,x,w1,x1
    
    ! implicit none
    
    ! integer,intent(in)::choix
    ! double precision,intent(out)::ss
    ! double precision::auxfunca,func4J,func5J
    ! external::func4J,func5J

    ! integer::j

    ! ss = 0.d0
    ! !if (typeof.eq.0) then
    ! !    do j=1,20
    ! !        if (choix.eq.1) then
    ! !            auxfunca=func4J(x(j))
    ! !            ss = ss+w(j)*(auxfunca)
    ! !        else
    ! !            if (choix.eq.2) then
    ! !                auxfunca=func5J(x(j))
    ! !                ss = ss+w(j)*(auxfunca)
    ! !            endif
    ! !        endif
    ! !    end do
    ! !else
        ! do j=1,32
            ! if (choix.eq.1) then
                ! auxfunca=func4J(x1(j))
                ! ss = ss+w1(j)*(auxfunca)
            ! else
                ! if (choix.eq.2) then
                    ! auxfunca=func5J(x1(j))
                    ! ss = ss+w1(j)*(auxfunca)
                ! endif
            ! endif
        ! end do
    ! !endif

    ! return

    ! end subroutine gaulagJ_intcens

! !==================================================================

    ! double precision function func4J(frail)

    ! use tailles
    ! use comon,only:g,c,auxig,alpha,theta,resL,resU, &
    ! res1,aux1,typeJoint,cdc
    ! use comongroup,only:nigdc

    ! IMPLICIT NONE

    ! double precision,intent(in)::frail
    ! double precision::prod
    ! integer::i,mi

    ! prod = 1.d0
    ! do i=1,nsujetmax
        ! if ((g(i).eq.auxig).and.(c(i).eq.1)) then
            ! prod = prod*(dexp(-frail*resL(i))-dexp(-frail*resU(i)))
        ! endif
    ! enddo

    ! if (typeJoint.eq.1) then
        ! mi = cdc(auxig)
    ! else
        ! mi = nigdc(auxig)
    ! endif

    ! func4J = dexp((alpha*mi+1.d0/theta-1.d0)*dlog(frail) &
    ! -frail/theta-frail*res1(auxig)-(frail**alpha)*aux1(auxig)) &
    ! *prod

! !    func4J = frail**(alpha*nigdc(auxig)+1.d0/theta-1.d0)* &
! !    dexp(-frail/theta-frail*res1(auxig)-(frail**alpha)*aux1(auxig)) &
! !    *prod

    ! if (func4J.lt.0.d0) then
! !        !print*,"func4J",func4J,prod,auxig
    ! endif

    ! return

    ! end function func4J

! !==================================================================

    ! double precision function func5J(frail)

    ! use tailles
    ! use comon,only:auxig,alpha,theta,res3,aux2

    ! IMPLICIT NONE

    ! double precision,intent(in)::frail

    ! func5J = dexp((1.d0/theta-1.d0)*dlog(frail) &
    ! -frail/theta-frail*res3(auxig)-(frail**alpha)*aux2(auxig))

! !    func5J = frail**(1.d0/theta-1.d0)* &
! !    dexp(-frail/theta-frail*res3(auxig)-(frail**alpha)*aux2(auxig))

    ! return

    ! end function func5J

! !=============================================================
! ! gauss hermite
! ! func est l integrant, ss le resultat de l integrale sur -infty , +infty

    ! SUBROUTINE gauherJ(ss,choix)

    ! use tailles
    ! use donnees,only:x2,w2,x3,w3
    ! use comon,only:auxig,typeof

    ! Implicit none
    
    ! double precision,intent(out)::ss
    ! integer,intent(in)::choix
    
    ! double precision::auxfunca,func6J
    ! external::func6J
    ! integer::j

    ! ss=0.d0
    ! if (typeof.eq.0) then
        ! do j=1,20
            ! if (choix.eq.3) then
                ! auxfunca=func6J(x2(j))
                ! ss = ss+w2(j)*(auxfunca)
            ! endif
        ! end do
    ! else
        ! do j=1,32
            ! if (choix.eq.3) then
                ! auxfunca=func6J(x3(j))
                ! ss = ss+w3(j)*(auxfunca)
            ! endif
        ! end do
    ! endif

    ! return
    
    ! END SUBROUTINE gauherJ

! !=====================================================================

    ! double precision function func6J(frail)

    ! use tailles
    ! use comon,only:auxig,alpha,sig2,res1,res3,aux1,nig,cdc

    ! IMPLICIT NONE

    ! double precision,intent(in)::frail

    ! func6J = frail*(nig(auxig)+alpha*cdc(auxig))- &
    ! dexp(frail)*(res1(auxig)-res3(auxig))- &
    ! dexp(alpha*frail)*aux1(auxig)- &
    ! (frail**2.d0)/(2.d0*sig2)

    ! func6J = dexp(func6J)

    ! return

    ! end function func6J
