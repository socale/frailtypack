!================================================================
! Remarque tres importante: Ne pas nommer les subroutines fortran 
! a appeler dans R en usant des majuscules
!================================================================
subroutine jointsurrogate(nsujet1,ng,ntrials1,maxiter,nst,nparamfrail,indice_a_estime,param_risque_base,nbrevar,&
                          filtre0,donnees,death,p,prop_i,n_sim1,EPS2,kappa0,vect_kappa,logNormal,nsim_node,Param_kendall_boot,&
                          vrai_val_init,param_init,revision_echelle,random_generator0,sujet_equi,prop_trait,paramSimul,&
                          autreParamSim,fichier_kendall,fichier_R2, param_estimes, sizeVect, b, H_hessOut,HIHOut,resOut,&
                          LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out,ni,ier,istop,ziOut, affiche_itter,Varcov,&
                          dataHessian,dataHessianIH,datab,vbetast,vbetastinit)
                          
    ! programme principal permettant le traitement des donnees et l'appel du joint_surogate pour l'estimation des parametres
    
    use sortie
    use Autres_fonctions
    use double_precision
    use var_surrogate, only: graine,aleatoire,nbre_sim,nbre_itter_PGH,nb_procs,random_generator,affiche_itteration, &
         copula_function
    use Autres_fonctions, only:pos_proc_domaine
    !use mpi ! module pour l'environnement MPI
    !$ use OMP_LIB 
                        
    implicit none

    ! =======debut declaration des variables================
    
    ! =====Parametres prises en entree de la subroutine=====
    integer, dimension(3),intent(in)::nbrevar
    integer,dimension(13), intent(inout)::nsim_node
    integer,intent(in)::nsujet1,ng,ntrials1,nst,maxiter,nparamfrail,n_sim1,logNormal,vrai_val_init,random_generator0,sujet_equi,&
                        affiche_itter
                        
    integer,dimension(5),intent(in)::indice_a_estime
    integer,dimension(5),intent(in):: param_risque_base
    integer,dimension(3),intent(in):: Param_kendall_boot
    integer,dimension(16),intent(in):: autreParamSim ! indique si l'on estime (1=oui, 0=non) estime ou pas zeta(1), covST(2), alpha(3), gammaST(4). indique en (5) si on prend en compte l'heterogeneite sur le risque de base
    integer,dimension(nbrevar(3),2), intent(in)::filtre0
    double precision,dimension(nsujet1,5+nbrevar(1)), intent(in):: donnees
    double precision,dimension(ng,5+nbrevar(2)), intent(in):: death
    double precision,dimension(2), intent(in):: kappa0
    double precision,dimension(9), intent(in):: param_init
    double precision,dimension(23+nbrevar(1) + nbrevar(2)-2), intent(in):: paramSimul
    !double precision,dimension(:), intent(in):: paramSimul
    double precision,dimension(3), intent(inout)::EPS2
    character(len=30),dimension(5)::NomFichier
    double precision, intent(in)::prop_trait,revision_echelle
    integer, dimension(5), intent(in)::sizeVect
    double precision, dimension(ntrials1), intent(in)::p,prop_i
    double precision,dimension(n_sim1,2), intent(in):: vect_kappa
    double precision, dimension(nbrevar(3),2), intent(in)::vbetast,vbetastinit

    
    ! ! =====Parametres fournies en sortie par la subroutine=====
    integer, intent(out):: ni, istop, ier
    double precision,dimension(n_sim1,3), intent(out):: fichier_kendall,fichier_R2
    double precision,dimension(n_sim1,nsim_node(13)),intent(out):: param_estimes
    !double precision,dimension(:,:),intent(inout):: param_estimes
    double precision,dimension(sizeVect(1)), intent(out)::b
    double precision,dimension(2), intent(out)::LCV
    double precision,dimension(sizeVect(2)), intent(out)::x1Out
    double precision,dimension(sizeVect(3)), intent(out)::x2Out
    double precision,dimension(sizeVect(1),sizeVect(1)), intent(out)::H_hessOut,HIHOut ! H_hessOut = matrice hesienne (des variance-covariance), HIHOut= matrice hessienne corrigee
    double precision,dimension(sizeVect(1)*n_sim1,sizeVect(1)), intent(out)::dataHessian, dataHessianIH ! sauvegarde des matrices hessiennes des differentes simulations 
    double precision,dimension(n_sim1,sizeVect(1)), intent(out)::datab ! sauvegarde des vecteurs de parametres de toutes les simulations 
    double precision,dimension(sizeVect(2),3), intent(out)::lamOut
    double precision,dimension(sizeVect(3),3), intent(out)::lam2Out
    double precision,dimension(sizeVect(4),3), intent(out)::suOut
    double precision,dimension(sizeVect(5),3), intent(out)::su2Out
    double precision,dimension(param_risque_base(5)+6), intent(out)::ziOut
    double precision,dimension(sizeVect(4)), intent(out)::xSu1
    double precision,dimension(sizeVect(5)), intent(out)::xSu2
    double precision, intent(out)::resOut
    double precision,dimension(3,3), intent(out):: Varcov ! pour la matrice de variance covariance de (sigma_S,sigma_ST_,sigma_T) par la delta-metode 

    ! =====Autres variables utilisees dans la subroutine
    !character(len=30)::donnees
    character(len=10), dimension(nbrevar(3))::nomvarl
    character(len=30)::dateamj,zone,heure1,heure2,param_estime, param_empirique,param_empirique_NC,tableau_rejet
    integer::i,j,effet,ver,nva1,nva2,nva,ag,nz, cpt,cpt_dc,noVar1,noVar2,k,typeJoint,np !ii,jj,ncur
    double precision::ax1,ax2,tp1,tp2 !tempon
    integer, dimension(:),allocatable::vdeces,vsurrogate !contient les dates devenement: deces et progression
    character(len=20),dimension(:),allocatable::nomvart,nomvar2t,nomvar,nomvar2
    double precision,dimension(:),allocatable::tt0dc,tt1dc
    integer,dimension(:),allocatable::icdc,groupe,pourtrial,ic,trials,nigs,cdcs,nigts,cdcts
    integer,dimension(:,:),allocatable::nig_Ts,cdc_Ts
    double precision,dimension(:,:),allocatable::vaxdc,vaxdct
    double precision,dimension(:),allocatable::tt0,tt1,ttU
    double precision,dimension(:,:),allocatable::vax,vaxt
    double precision,dimension(2)::k0,k0_save,k01_save,ckappa
    double precision,dimension(3)::EPS
    integer,dimension(8)::values
    integer, dimension(0:1)::randomisation,deces,surrogate
    double precision::bi,bs,wres !wald
    character(len=30)::kapa !aaa !les fichiers de sortie
    integer,dimension(:),allocatable::filtre,filtre2
!cpm
    integer::mt11,mt12,mt1,mt2,n_sim,ntrials,nsujet
    double precision,dimension(2)::shape_weib,scale_weib
    integer::typeof,nbintervR,nbintervDC,equidistant !nbrecu,nbdeces
    !double precision::Xgamma
!predictor
    double precision,dimension(:,:),allocatable::MartinGales,v_chap_kendall,v_chap_R2,theta_chap_kendall,theta_chap_R2, &
                                                theta_chap_copula, v_chap_copula
    double precision,dimension(:),allocatable::linearpred,vect_kendall_tau,t_chap_kendall,t_chap_R2,vect_kendall_tau_temp,vect_R2
    double precision,dimension(:),allocatable::linearpreddc
    double precision,dimension(:),allocatable::time
    double precision,dimension(:),allocatable::timedc,vbetas,vbetat, vbetas_intit, vbetat_intit
    integer,dimension(4)::mtaille
    integer,dimension(3)::paratps
    double precision,dimension(4)::paraweib
    double precision,dimension(3)::descripSurr,descripDeces
    double precision,dimension(:,:),allocatable:: paGH,matrice_generation ! parametre pour l'adaptative: en ligne les individus, en colone on a respectivement: les ui_cham,
        !     racine carree du determinant de l'inverse de la cholesky,variance des ui_chap,les covariances estimees des fragilites pour chaque individu, sachant que la matrice de variances covariance est bien la cholesky                                                        
    !parametres de simulation
    integer::n_col,mode_cens,n_essai,n_obs,weib,frailty_cor,affiche_stat,s_i,indice_eta,indice_theta&
                ,rangparam,rangparam2,nbre_rejet,ind_temp,seed_,une_donnee,gener_only,kapa_use,&
                i_min,i_max,i_min_t,i_max_t,ind_rech,Rech_kappa,incre_kappa,statut_kappa,statut_kappa1&
                ,indice_varS,indice_varT,indice_covST,rangparam_sigs,rangparam_sigt,rangparam_sigst,&
                np_save,control_kappa,ind_premier_kappa,control,control2,frailt_base,&
                indice_gamma,indice_alpha,rangparam_gamma,nbre_rejet_0,indice_gamma_st,indice_theta_t,&
                indice_theta_st,indice_gamma_t,rangparam_thetat,rangparam_thetast,rangparam_gammat,&
                rangparam_gammast,rangparam_alpha,decoup_simul,incre_decoup,method_int_kendal,N_MC_kendall,&
                param_weibull,donne_reel,indice_seed,npoint1,npoint2,rangparam_eta,nboot_kendal,nparam_kendall,&
                rangparam_theta,erreur_fichier,indicCP,controlgoto,remplnsim,indice_kapa, type_joint_estim,&
                rangparam_copula,pfs
                
                
    double precision::theta,eta,betas,alpha,betat,lambdas,nus,lambdat,nut,temps_cens,&
                        cens0,rsqrt,sigma_s,sigma_t,moy_theta,theta_sim,moy_dec,moy_cens,&
                        moy_pros,moy_theta_est,moy_se_theta,moy_eta,moy_se_eta,&
                        bi2,bs2,n_sim_exact,&
                        moy_trt,taux_couverture_theta,taux_couverture_eta,bi_theta,bs_theta,bi_eta,bs_eta,&
                        gamma1,gamma2,theta2,sigmas_sim,sigmat_sim,rho_sim,varS1,varT1,covST1,varS_es,varT_es,covST_es,&
                        moy_sigmas,moy_sigmat,moy_rho,moy_sigmas_est,moy_sigmat_est,moy_sigmast_est,moy_se_sigmas,moy_se_sigmat,&
                        moy_se_sigmast,bi_sigmas,bi_sigmat,bi_sigmast,bs_sigmas,bs_sigmat,bs_sigmast,taux_couverture_sigmas,&
                        taux_couverture_sigmat,taux_couverture_sigmast,sigmast_vrai,sigma_ss_init,sigma_tt_init,sigma_st_init,&
                        theta_init,betas_init,betat_init,moy_ni, gamma_ui,alpha_ui,gamma_sim,moy_gamma,gamma_init,alpha_init,&
                        moy_gamma_est,moy_se_gamma,bi_gamma,bs_gamma,taux_couverture_gamma,moy_ni_0,moy_trt_0,moy_theta_0,&
                        moy_sigmas_0,moy_sigmat_0,moy_rho_0,moy_gamma_0,tab_var_sigma_0,moy_se_theta_0,taux_couverture_theta_0,&
                        moy_gamma_est_0,moy_se_gamma_0,taux_couverture_gamma_0,moy_sigmas_est_0,moy_se_sigmas_0,&
                        taux_couverture_sigmas_0,moy_sigmat_est_0,moy_se_sigmat_0,taux_couverture_sigmat_0,&
                        moy_sigmast_est_0,moy_se_sigmast_0,taux_couverture_sigmast_0,moy_eta_0,moy_se_eta_0,&
                        taux_couverture_eta_0,moy_betaS_0,moy_betaS_se_0,taux_couvertureS_0,moy_betaT_0,&
                        moy_betaT_se_0,taux_couvertureT_0,se_theta_sim,se_sigmas_sim,se_sigmat_sim,se_rho_sim,&
                        se_gamma_sim,& !se_theta_sim_0,se_sigmas_sim_0,se_sigmat_sim_0,se_rho_sim_0,se_gamma_sim_0
                        n_sim_exact_0,moy_theta_est_0,moy_pros_0,moy_dec_0,theta2_t,rsqrt_theta,gamma_uit,rsqrt_gamma_ui,&
                        thetat_init,thetast_init,gammat_init,gammast_init,theta_simt,rho_sim_wij,gamma_simt,rho_sim_ui,&
                        moy_thetat,moy_rho_wij,moy_gammat,moy_rho_ui,thetast_vrai,gammast_vrai,moy_thetat_est,moy_se_thetat,&
                        taux_couverture_thetat,moy_thetast_est,moy_se_thetast,taux_couverture_thetast,moy_gammat_est,moy_se_gammat,&
                        taux_couverture_gammat,moy_gammast_est,moy_se_gammast,taux_couverture_gammast,thetaS1,thetaT1,thetaST1,&
                        thetaT_es,thetaST_es,gammaS1,gammaT1,gammaST1,gammaS_es,gammaT_es,gammaST_es,moy_alpha,moy_se_alpha,&
                        taux_couverture_alpha,se_theta_simt,se_rho_sim_wij,se_gamma_simt,se_rho_sim_ui,se_theta_est,se_eta_est,&
                        se_beta_s,se_beta_t,se_sigmas_est,se_sigmat_est,se_cov_est,se_gamma_est,se_alpha_est,se_thetat_est,&
                        se_cov_est_wij,se_gamma_estt,se_cov_est_ui,se_theta_est_0,se_eta_est_0,se_beta_s_0,se_beta_t_0,&
                        se_sigmat_est_0,se_cov_est_0,se_gamma_est_0,se_alpha_est_0,se_thetat_est_0,se_gamma_estt_0,thetaS_es,&
                        se_cov_est_ui_0,moy_thetat_0,se_theta_simt_0,moy_thetat_est_0,moy_se_thetat_0,taux_couverture_thetat_0,&
                        moy_rho_wij_0,se_rho_sim_wij_0,moy_thetast_est_0,se_cov_est_wij_0,moy_se_thetast_0,&
                        moy_gammat_0,se_gamma_simt_0,moy_gammat_est_0,moy_se_gammat_0,taux_couverture_gammat_0,&
                        moy_rho_ui_0,se_rho_sim_ui_0,moy_gammast_est_0,moy_se_gammast_0,taux_couverture_gammast_0,&
                        moy_alpha_0,moy_se_alpha_0,taux_couverture_alpha_0,R2_trial,se_R2_trial,moy_R2_trial,&
                        taux_couverture_R2_trial,moy_se_R2_trial,moy_bi_R2_trial,moy_bs_R2_trial,moy_kendal_11,tau_kendal_11,&
                        moy_kendal_10,tau_kendal_10,moy_kendal_01,tau_kendal_01,moy_kendal_00,tau_kendal_00,se_kendal_11,&
                        se_kendal_01,se_kendal_00,moy_tau_boots,IC_Inf,IC_sup,zeta_init,moy_R2_boots,IC_Inf_R2,IC_sup_R2,&
                        CP_R2_boot,CP_ktau_boot,moy_R2_boots_test,se_sigmas_est_0,taux_couverture_thetast_0,se_kendal_10,&
                        bi_R2_trial,bs_R2_trial,thetacopule, thetacopula_init, printnbre, moy_param_cop, moy_se_param_cop,&
                        bi_param_cop, bs_param_cop, pour_ic, taux_couverture_param_cop,taux_couverture_tauk, vrai_tau_copula
                        
    double precision, dimension(:,:),allocatable::don_simul,don_simulS, don_simultamp,don_simulStamp,don_simulS1,&
                        parametre_empirique, parametre_estimes,parametre_empirique_NC,parametre_estimes_MPI,&
                        parametre_estimes_MPI_T,result_bootstrap
    double precision, dimension(:),allocatable::tab_var_theta,tampon,tampon_all, moy_betaS, moy_betaT,moy_betaS_se, moy_betaT_se,&
                        taux_couvertureS, taux_couvertureT
    double precision, dimension(:,:),allocatable::tab_var_sigma
    integer,parameter ::trt1=1,v_s1=2,v_t1=3,trialref1=4,w_ij1=5,timeS1=6,timeT1=7,&
                      timeC1=8,statusS1=9,statusT1=10,initTime1=11,Patienref1=12,u_i1=13,& 
                          w_ijt=14,u_it=15 ! definissent les indices du tableau de donnee simulees
    integer, dimension(:),allocatable::tab_rejet ! pour les rangs des jeux de donnees rejetees
    double precision, dimension(:,:),allocatable::kappa,d_S,d_T !jeu de donnees reelle pour le test
    integer,dimension(:),allocatable::tableEssai,tableNsim ! tableNsim: indique le nombre de simulation a effectuer par chaque processus
    double precision,dimension(:,:),allocatable::donnee_essai,theta_st_2,gamma_st_2,theta_st0_2,gamma_st0_2 !sigma_st_2,sigma_st0_2
    double precision,dimension(2,2)::chol,sigma_st,theta_st,gamma_st,sigma_st0,theta_st0,gamma_st0,Chol_R2,mat_A
    integer, dimension(4)::indice_esti
    integer::nb_processus,rang,code,n_sim_total,suplement,erreur,comm,init_i,max_i,debut_exe,indice_sim_proc,sofeu, &
            rang_proc,init_i_proc,max_i_proc, code_print ! je redefini ces indices car les precedentes sont utilisees autrement: cas OpenMP
    double precision,dimension(10)::t
    double precision,dimension(3,3):: sigmac ! pour la mtrice de variance covariance de Sigma par la delta-metode 
    double precision,dimension(3,3):: hb 
    
    !=====================================================================================
    !*********fin declaration des variables et debut du programme principale**************
    !=====================================================================================
    ! test operators priority
    ! printnbre = 1.d0/2.d0 *5.d0
    ! call dblepr("test 1/2 *5", -1, printnbre, 1)
    ! printnbre = (1.d0/2.d0) *5.d0
    ! call dblepr("test (1/2) *5", -1, printnbre, 1)
    ! call dblepr("voile kappa", -1, vect_kappa, n_sim1)
    ! goto 998
    
    ! for copula model 
    copula_function = nsim_node(12) ! the copula function, can be 1 for clayton or 2 for Gumbel-Hougaard
    type_joint_estim = nsim_node(8) ! type of estimated model
    ! affectation de certains parametres
    nomvarl(1) = "trt"
    NomFichier(1) = "kappa_valid_crois.txt"
    NomFichier(2) = "Parametre_estime.txt"
    NomFichier(3) = "Parametre_empirique.txt"
    NomFichier(4) = "Parametre_empirique_NC.txt"
    NomFichier(5) = "tab_rejet.txt"
    rang_proc = 0
    rangparam_eta = 0
    rangparam_gammast = 0
    rangparam_thetat = 0
    rangparam_thetast = 0
    rangparam_gammat = 0
    statut_kappa = 0
    
    pour_ic = 0.d0
    R2_trial = 0.d0
    se_R2_trial = 0.d0
    tau_kendal_10 = 0.d0
    tau_kendal_01 = 0.d0
    tau_kendal_00 = 0.d0
    tau_kendal_11 = 0.d0
    thetaS1 = 0.d0 
    thetaT1 = 0.d0 
    thetaST1 = 0.d0
    thetaT_es = 0.d0 
    thetaST_es = 0.d0
    thetat_init = 0.d0 
    thetast_init = 0.d0
    varS_es = 0.d0 
    varT_es = 0.d0
    gammaS1 = 0.d0 
    gammaT1 = 0.d0
    gammaST1 = 0.d0
    gammaS_es = 0.d0
    gammaT_es = 0.d0
    gammaST_es = 0.d0
    gammat_init = 0.d0
    gammast_init = 0.d0
    moy_tau_boots = 0.d0
    moy_R2_boots = 0.d0
    
    
    
    affiche_itteration = affiche_itter
    n_sim = n_sim1 ! nombre de simulations
    ntrials = ntrials1
    nsujet = nsujet1
    typeof = param_risque_base(1)    !type de function de risque  0:Splines,  1:Cpm  2:weib
    nbintervR = param_risque_base(2) !Nombre intervalle surrogate
    nbintervDC = param_risque_base(3) !Nombre intervalle deces
    equidistant = param_risque_base(4) !cpm (1:equidistant, 0:percentile)')
    nz = param_risque_base(5) ! nombre de noeud pour les spline
    
    ! indisue si l'on estime (1=oui, 0=non) estime ou pas zeta(1), covST(2), alpha(3), gammaST(4). indique en (5) si on prend en compte l'heterogeneite sur le risque de base
    indice_eta=indice_a_estime(1)    !dit si l'on estime eta (1) ou pas (0)
    indice_covST=indice_a_estime(2) !dit si l'on estime la covariance des frailties essais ou pas
    indice_alpha=indice_a_estime(3)    !dit si l'on estime alpha (1) ou non(0)
    indice_gamma_st=indice_a_estime(4) !dit si l'on estime sigma_us_ut (1) ou non(0)
    frailt_base=indice_a_estime(5) !dit si l'on prend en compte l'heterogeneite sur le risque de base aussi bien dans la generation des donnes que dans l'estimation(1) ou non (0)
    !call intpr("I'm there scl 13:", -1, code_print, 1)
    ver=nbrevar(3) ! nombre total de variables explicatives
    allocate(vbetas(ver),vbetat(ver),vbetas_intit(ver), vbetat_intit(ver))
    
    AG=0 ! andersen-gill approach(1=oui) 
    
    ! Gestion des kappas pour le jeux de donnees reelles
    ax1 = kappa0(1)    ! pour surrogate
    ax2 = kappa0(2)    ! pour true
    kapa= NomFichier(1)
    
    ! Parametres associes au taux de kendall et au bootstrap
    method_int_kendal = Param_kendall_boot(1)
    N_MC_kendall = Param_kendall_boot(2)
    nboot_kendal = Param_kendall_boot(3)
    
    ! Autres noms de fichiers:
    param_estime = NomFichier(2) 
    param_empirique = NomFichier(3)
    param_empirique_NC = NomFichier(4)
    tableau_rejet = NomFichier(5)

    ! Parametres initiaux
    theta_init = param_init(1) ! if we are estimating the joint surrogate model
    thetacopula_init = param_init(1) ! if we are estimating copula modele
    sigma_ss_init = param_init(2)
    sigma_tt_init = param_init(3)
    sigma_st_init = param_init(4)
    gamma_init = param_init(5)
    alpha_init = param_init(6)
    zeta_init = param_init(7)
    betas_init = param_init(8)
    betat_init = param_init(9)
    random_generator = random_generator0

    ! parametres de simulation
    gamma1 = paramSimul(1) 
    gamma2 = paramSimul(2)
    theta2 = paramSimul(3)
    eta = paramSimul(4)
    gamma_ui = paramSimul(5)
    alpha_ui = paramSimul(6)
    theta2_t = paramSimul(7)
    rsqrt_theta = paramSimul(8)
    gamma_uit = paramSimul(9)
    rsqrt_gamma_ui = paramSimul(10)
    betas = paramSimul(11)
    betat = paramSimul(12)
    lambdas = paramSimul(13)
    nus = paramSimul(14)
    lambdat = paramSimul(15)
    nut = paramSimul(16)
    mode_cens = int(paramSimul(17))
    temps_cens = paramSimul(18)
    cens0 = paramSimul(19)
    rsqrt = paramSimul(20)
    sigma_s = paramSimul(21)
    sigma_t = paramSimul(22)
    thetacopule = paramSimul(23)
    ! les les parametres de simulation pour les autres covariables se trouvent a la fin du tableau paramSimul
    
    !call dblepr("paramSimul = ",-1,paramSimul,size(paramSimul))
    if(nsim_node(11) == 3) then ! joint frailty copula, remplissage des vecteurs des variables explicatives
        do i = 1, size(vbetast,1)
            vbetas(i) = vbetast(i,1) ! beta_s
            vbetat(i) = vbetast(i,2) ! beta_t
            vbetas_intit(i) = vbetastinit(i,1) ! beta_s
            vbetat_intit(i) = vbetastinit(i,2) ! beta_t
            !call dblepr("vbetast", -1, vbetast(i,:), size(vbetast,2))
        enddo
    endif
    
    if(nsim_node(11) == 1) then ! joint surrogate, remplissage des vecteurs des variables explicatives
        vbetas(1) = betas ! beta_s
        vbetat(1) = betat ! beta_t
    endif
    
    ! autres parametres de simulation
    weib = autreParamSim(1)
    param_weibull = autreParamSim(2)
    frailty_cor = autreParamSim(3)
    affiche_stat = autreParamSim(4)
    seed_ = autreParamSim(5)
    une_donnee = autreParamSim(6)
    donne_reel = autreParamSim(7)
    gener_only = autreParamSim(8)
    kapa_use = autreParamSim(9)
    decoup_simul = autreParamSim(10)
    aleatoire = autreParamSim(11)
    nbre_sim = autreParamSim(12)
    graine = autreParamSim(13)
    ckappa(1) =dble(autreParamSim(14))
    ckappa(2) =dble(autreParamSim(15))
    pfs = autreParamSim(16)
    ! call intpr("avant appel joint:pfs", -1,pfs, 1)
    
    np=sizeVect(1)
    call date_and_time(dateamj,heure1,zone,values) ! pour la date de debut du programme
    
    rang=0
    nb_processus=1
    nb_procs=1
    controlgoto=0

    ! initialisation de l'environnement de travail pour le parallelisme MPI 
    !call MPI_INIT(code) ! cet environnement est desactive a la fin du programme a l'aide de MPI_FINALIZE
    !call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code) ! commaitre le nombre de processus associe a l'identificateur code
    !call MPI_COMM_RANK(MPI_COMM_WORLD,rang_proc,code)
    !desactivation de l'environnement de travail pour le programme parallele
    !!call MPI_FINALIZE(comm) 

    if(rang_proc==0) then
        !write(*,*)
        !write(*,*)'******************************************'
        !write(*,*)'******* DEBUT PROGRAMME SURROGACY ********'
        !write(*,*)'******************************************'
    endif
    
    typeJoint=8 ! 8 pour les models conjoint surrogate (sans reccurence) et true endpoint

    !Ouverture et lecture des fichiers d'entree sortie 
    ! open(2,file='joint_scl_simul.inf')
    !open(4,file='OutJoint_simul.txt')
    
    
    ! read(2,*)nsujet ! nombre d'observations fichier de donnes surrogate
    allocate(groupe(nsujet),ic(nsujet),tt0(nsujet),tt1(nsujet),ttU(nsujet),linearpred(nsujet),pourtrial(nsujet))

    ! read(2,*)ng ! nombre d'observations fichier de donnes deces (true endpoint)
    ! read(2,*)ntrials ! nombre total d'essai
    allocate(tt0dc(ng),tt1dc(ng),icdc(ng),trials(ntrials),nigs(ntrials),cdcs(ntrials),nigts(ntrials),cdcts(ntrials),&
    nig_Ts(ntrials,2),cdc_Ts(ntrials,2))
    allocate(MartinGales(ng,4),linearpreddc(ng))
    ! read(2,*)maxiter         !nb maximum d'iteration
    ! read(2,*)nst            !nb de fonction de risque (2 pour )
    ! read(2,*)nparamfrail     !nombre de parametres associes aux effets aleatoires: eta+theta+varvs+varvt+covst
    ! read(2,*)indice_eta        !dit si l'on estime eta (1) ou pas (0)
    ! read(2,*)indice_covST     !dit si l'on estime la covariance des frailties essais ou pas
    ! read(2,*)indice_alpha   !dit si l'on estime alpha (1) ou non(0)
    ! read(2,*)indice_gamma_st !dit si l'on estime sigma_us_ut (1) ou non(0)
    ! read(2,*)frailt_base    !dit si l'on prend en compte l'heterogeneite sur le risque de base aussi bien dans la generation des donnes que dans l'estimation(1) ou non (0)

    if(indice_eta==0 .and. (nparamfrail==2 .or. nparamfrail==5)) then
        !!print*,"Attention si indice_eta=0, nparamfrail==1 ou >2"
        !stop
    endif
    
    indice_esti(1)=indice_eta
    indice_esti(2)=frailt_base
    indice_esti(3)=indice_alpha
    indice_esti(4)=indice_gamma_st
    
    !indicateur de presence deffet aleatoire (0=Non, 1=Oui)
    if(nparamfrail.eq.0) then
        effet=0                
    else
        effet=1
        indice_theta=1
        indice_varS=1
        indice_varT=1
        !indice_covST=1
        if(frailt_base==1) then
            indice_gamma=1
            !indice_alpha=1
        endif
        
        if(nsim_node(8)==2) then
            indice_gamma=1
            indice_theta_t=1
            indice_theta_st=1
            indice_gamma_t=1
        endif
    end if
    
    if(type_joint_estim == 3) then ! joint frailty-copula model
        indice_theta = 0
        indice_theta_t = 0
        indice_theta_st = 0
    endif

    !!write(4,*)'**************************************************'
    !!write(4,*)'*****************JOINT MODEL *********************'
    !!write(4,*)'**********SURROGATE AND TRUE ENDPOINTS ***********'

    ! read(2,*)typeof !type de function de risque  0:Splines,  1:Cpm  2:weib
    
    select case(typeof)
        case(0)
            !!write(4,*)'*** Splines ***'
            ! read(2,*)         ! on saute la ligne associee aux parametres de constance par morceaux
        case(1)
            !!write(4,*)'*** Constante par morceaux 2010 ***'
        case(2)
            !!write(4,*)'*** Weibull 2010 ***'
    end select
    if(rang_proc==0) then
        !!write(4,*)'**************************************************'
        !!write(4,*)trim('---> le nb dessais vaut ='),ntrials 
        !write(*,*)trim('---> le nb dessais vaut ='),ntrials 
        !!write(4,*)trim('---> le nb de sujets ='),ng 
        !write(*,*)trim('---> le nb de sujets ='),ng 
        !!write(4,*)trim('---> nombre de fonctions de risque de base ='),nst
        !write(*,*)trim('---> nombre de fonctions de risque de base ='),nst

        
        if (effet.eq.1)then
            !write(*,*)trim('---> Nombre de parametres associes aux effets aleatoires:'), nparamfrail
            !!write(4,*)trim('---> Nombre de parametres associes aux effets aleatoires:'), nparamfrail
        else
            !write(*,*)'-----modele SANS effet aleatoire'
            !!write(4,*)'-----modele SANS effet aleatoire'
        endif    
            

        !!write(4,*)trim('---> Type de fonction de risque (0:Splines 1:cpm 2:Weib)'),typeof
        !write(*,*)trim('---> Type de fonction de risque (0:Splines 1:cpm 2:Weib)'),typeof
    endif
    
    mt11 = 100 !scl  A quoi sert ce parametre?
    mt12 = 100
    
    if(typeof == 1) then
        ! read(2,*)nbintervR,nbintervDC,equidistant
        !!write(4,*)trim('---> Nombre intervalle surrogate'),nbintervR
        !!write(4,*)trim('---> cpm (1:equidistant, 0:percentile)'),equidistant
        mt1 = 3*nbintervR    !scl: pourquoi *3
        mt2 = 3*nbintervDC
    else
        if(typeof == 2) then
            ! read(2,*)nbintervR,nbintervDC
        end if
        mt1 = 100
        mt2 = 100

    end if
    
    ! allocate(suOut(mt11,3),su2Out(mt12,3)) ! matrice des estimates de la survie a baseline(avec IC) surrogate et deces
    ! allocate(x1Out(mt1))     ! temps pour la representation des fonction de risque et de la survies surrogate
    ! allocate(x2Out(mt2))     ! temps pour la representation des fonction de risque et de la survies true endpoint
    ! allocate(lamOut(mt1,3))    ! estimates du risque de base surrogate avec IC
    ! allocate(lam2Out(mt2,3))! estimates du risque de base deces avec IC
    
    !indice_alpha = 1 !pas necessaire pour moi

    allocate(time(nbintervR+1),timedc(nbintervDC+1))

    !read(2,*)ver    ! nombre de variables explicatives
    !if(rang_proc==0)!write(*,*)
    nva1 = 0         ! nb de var expli pour donnees recurrentes
    nva2 = 0         ! nb de var expli pour deces
     
    allocate(filtre(ver),filtre2(ver),nomvart(ver),nomvar2t(ver))!scl filtre indique si la variable est prise en compte pour les reccures et filtre2 dit si elle est prise en compte pour les deces
    allocate(vaxt(nsujet,ver),vaxdct(ng,ver)) ! matrice de toutes les variables explicatives
    
       ! indicatrice de prise en compte de variable
    do i=1,ver
        filtre(i)=filtre0(1,i)
        filtre2(i)=filtre0(2,i)
    enddo
    
    if(rang_proc==0) then
        ! recuperation des noms des variables
        !write(*,*)trim("---> la variable: appartient-elle aux fichiers Surrogate et True endpoint resp.?(1=Oui, 0=Non)")
        !write(*,*)trim("     Nom des variable"),trim("     Surrogate"),trim("     Deces")
    endif

    if(ver.gt.0)then ! aumoins une variable explicative comme c'est le cas
        do j=1,ver
           ! read(2,*)nomvarl,filtre(j),filtre2(j) ! nom de la variable + indicateur d'appartenance aux deux fichiers surrogate et true
            !if(rang_proc==0) !write(*,*)"       ",nomvarl(j)," ",filtre(j)," ",filtre2(j)
            nva1 = nva1 + filtre(j) ! adjustment for recurrent events
            nva2 = nva2 + filtre2(j) ! adjustment for survival
            if(filtre(j).eq.1)then
                nomvart(nva1) = nomvarl(j)
            endif  
            if(filtre2(j).eq.1)then
                nomvar2t(nva2) = nomvarl(j)
            endif   
        end do
    endif    

    !On active ou pas les filtres
    !------> filtre 1
    if (nva1 .ne. 0) then
        noVar1 = 0 !
        
      !  if(rang_proc==0) !write(*,*)'****** Présence de variables explicatives pour surrogate   *******'
    else
        noVar1 = 1
    end if
    
    !------> filtre 2
    if (nva2 .ne. 0) then
        noVar2 = 0
       ! if(rang_proc==0) !write(*,*)'********* Présence de variables explicatives pour décès **********'
    else
        noVar2 = 1
    end if
    
    if(rang_proc==0) then
        !write(*,*)trim('---> noVar1:(1= aucune covariable, 0= aumoins une covariable)'),novar1
        !write(*,*)trim('---> noVar2:(1= aucune covariable, 0= aumoins une covariable)'),novar2
    endif
    
    allocate(vax(nsujet,nva1)) ! matrice des variables explicatives presentes dans le jeux de donnees surrogate
    allocate(nomvar(nva1),nomvar2(nva2)) !Nnom des variables explicatives presentes dans le jeux de donnees surrogate et deces
    
    if(ver.gt.0)then
        nomvar=nomvart(1:nva1)
        nomvar2=nomvar2t(1:nva2)
    end if
    
    nva = nva1+nva2

    allocate(vaxdc(ng,nva2)) ! matrice des variables explicatives presentes dans le jeux de donnees deces
    allocate(moy_betaS(nva1), moy_betaT(nva2),moy_betaS_se(nva1), moy_betaT_se(nva2),&
            taux_couvertureS(nva1), taux_couvertureT(nva2))
    !!write(4,*)trim('---> Nombre de variables explicatives pour les donnees surrrogate:'),nva1
   ! if(rang_proc==0) !write(*,*)trim('---> Nombre de variables explicatives pour les donnees surrrogate:'),nva1
    !!write(4,*)trim('---> Nombre de variables explicatives pour les donnees deces:'),nva2
   ! if(rang_proc==0) !write(*,*)trim('---> Nombre de variables explicatives pour les donnees deces:'),nva2
         
    ! read(2,*)donnees     ! Nom du fichier des donnees surrogacy
    ! read(2,*)death        ! Nom du fichier des donnees deces
    !if(rang_proc==0) !write(*,*)trim('---> Nom des fichiers de donnees = '),donnees,'et' , death
    !!write(4,*)trim('---> Nom des fichiers de donnees = '),donnees,'et' , death  
    
    !read(2,*)n_sim        !nombre de simulation a faire pour evaluer les moyennes empirique
    !if(rang_proc==0) !write(*,*)trim('---> Nombre de simulation a faire pour evaluer les moyennes empirique = '),n_sim
    !!write(4,*)trim('---> Nombre de simulation a faire pour evaluer les moyennes empirique = '),n_sim  
    
    ! read(2,*)AG                ! andersen-gill approach(1=oui)
    ! read(2,*)nz                ! nombre de noeud pour les spline
    ! allocate(ziOut(nz+6))    ! vecteur des noeuds duu risque de base estimer par les splines
    
    ! read(2,*)EPS2    ! critere de convergence du modele
    
    !Deux parametres de lissage :
    ! read(2,*)ax1 ! kappa pour surrrogate
    ! read(2,*)ax2 ! kappa pour deces
    ! read(2,*)kapa ! nom du fichier pour les kappa de la validation croisee

    ! fichiers qui contiendont les estimation des fonctions de risque    
    if(effet.eq.1)then
        ! read(2,*)fich1 !hazard  surrogate
        ! read(2,*)fich2 !Survie  surrogate
        ! read(2,*)fich3 !hazard true
        ! read(2,*)fich4 !Survie true
        ! read(2,*)fich5 !cumul hazard str 1
        ! read(2,*)fich6 !cumul hazard str 2
    end if
     ! read(2,*)logNormal!indique si on a une distribution lognormale des effets aleatoires (1) ou Gamma (0)
    ! read(2,*)nsim_node(1) ! nombre de simulation pour l'integration par Monte carlo, vaut 0 si on ne veut pas faire du MC
    ! read(2,*)nsim_node(2) ! nombre de points de quadrature a utiliser (preference 5 points pour l'adaptatice et 32 poits pour la non adaptatice)
    npoint1=nsim_node(2) ! nombre de point de quadrature a privilegier initialement
    npoint2 = nsim_node(9) ! nombre de point de quadrature a utiliser en cas de non convergence de prefenrence 7 ou 9 pour la pseudo adaptative et 32 pour la non adaptative
    ! read(2,*)nsim_node(3) ! doit-on faire de l'adaptative(1) ou de la non-adaptative(0)
    nbre_itter_PGH = nsim_node(10) !nombre d'itteration aubout desquelles reestimer les effects aleatoires a posteriori pour la pseude adaptative. si 0 pas de resestimation
    ! read(2,*)nsim_node(4) ! indique la methode d'integration 0=Monte carlo,1= quadrature, 2=quadrature individuel+MC essai, 3=Laplace, 4= monte carlo individuel + quadrature essai
    !nsim_node(5)=nparamfrail ! indique le nombre de parametres associes aux effets aleatoires
    
    ! if(nsim_node(4) .ne. 3) then ! cac pour laplace on parallelise dans le calcul integral
        ! !call MPI_INIT(comm) ! cet environnement est desactive a la fin du programme a l'aide de MPI_FINALIZE
        ! !call MPI_COMM_SIZE( MPI_COMM_WORLD,nb_processus,comm) ! commaitre le nombre de processus associe a l'identificateur code
        ! !call MPI_COMM_RANK( MPI_COMM_WORLD,rang,comm)! pour chaque processus associe a l'identificateur code retourne son rang
    ! endif
    
    ! read(2,*)nsim_node(6) ! indique si lon fait de la vectorisation dans le calcul integral (1) ou non (0). rmq: la vectorisation permet de reduire le temps de calcul
    ! read(2,*)nsim_node(7) ! indique le nombre d'effet aleatoire cas quadrature adaptative
    ! read(2,*)nsim_node(8) ! type de modele a estimer: 0=joint classique avec un effet aleatoire partage au niveau individuel,1=joint surrogate avec 1 frailty partagé indiv et 2 frailties correles essai,
                          ! 2=joint surrogate sans effet aleatoire partage donc deux effets aleatoires a chaque fois"
    np_save=nsim_node(2)
    ! read(2,*)method_int_kendal ! methode d'integration pour le taux de kendall: 0= montecarle, 1= quadrature quaussienne classique, 2= approximation de Laplace
    ! read(2,*)N_MC_kendall ! nombre de boucle MC pour le calcul du taux de kendal en approximant l'integrale par montye carlo
    ! read(2,*)nboot_kendal ! nombre d'echantillon bootstrap pour le calcul de l'IC du taux de ke,ndall
    ! read(2,*)fichier_kendall ! fichier dans lequel saugarder les taux de kendall avec les IC par boostrap
    
    nparam_kendall=4 ! on a 4 parametres qui rentrent dans le calcul du tau de kendall: theta, alpha, gamma, zeta
    
    if(method_int_kendal==4) then
        if(indice_alpha==0) nparam_kendall=nparam_kendall-1
        if(indice_eta==0) nparam_kendall=nparam_kendall-1        
    endif
    if(method_int_kendal==5) then
        nparam_kendall=2 ! dans ce cas on a fixe u_i =0, ce qui annule gamma et alpha
        if(indice_eta==0) nparam_kendall=nparam_kendall-1
    endif

    allocate(vect_kendall_tau(nboot_kendal),v_chap_kendall(nparam_kendall,nparam_kendall),vect_R2(nboot_kendal),&
            theta_chap_kendall(1,nparam_kendall),t_chap_kendall(nparam_kendall),v_chap_R2(3,3),t_chap_R2(3),&
            theta_chap_R2(1,3),result_bootstrap(n_sim,6), theta_chap_copula(1,1), v_chap_copula(1,1))

    indice_kapa = 1    
    n_essai=ntrials
    n_obs=nsujet
    ! le jeu de donnee doit contenir 13 ou 17 colonnes
    !read(5,*)n_col    ! je ne lis plus cette variable
    if(nsim_node(8)==2)then
        n_col=15 + nbrevar(3)-1 ! j'ajoute le surplus des covariables, -1 pour le traitement qui est deja pris en compte
    else
        n_col=13 + nbrevar(3)-1 ! j'ajoute le surplus des covariables
    endif
    alpha = eta    ! alpha associe a u_i chez les deces
    
    !generation des donnees par joint failty-copula
    if(nsim_node(11)==3) allocate(don_simultamp(n_obs,n_col),don_simulStamp(n_obs,n_col))
        
    allocate(don_simul(n_obs,n_col),don_simulS1(n_obs,n_col))
    
    if(nsim_node(8)==2)then
        allocate(donnee_essai(n_essai,5))
    else
        allocate(donnee_essai(n_essai,4))
    endif
    
    don_simul=0.d0
    ! read(5,*)weib ! 0= on simule les temps par une loi exponentielle, 1= on simule par une weibull
    ! read(5,*)param_weibull ! parametrisation de la weibull utilisee: 0= parametrisation par defaut dans le programme de Virginie, 1= parametrisation a l'aide de la fonction de weibull donnee dans le cous de Pierre
    ! read(5,*)frailty_cor ! indique si l'on considere pour le modele de simulation deux effets aleatoire correles au niveau essai(=1) ou un effet aleatoire partage(=0) ou encore on simule sans effet aleatoire au niveau essai(=2, model conjoint classique)
    ! read(5,*)affiche_stat ! dit si l'on affiche les statistiques des donnees simulees(1) ou non (0)
    ! read(5,*)seed_  !jeux de donnees a retenir pour la validation croisee
    ! read(5,*)une_donnee ! pour dire si on simule avec un seul jeu de donnees(1) ou pas (0). ceci pour tester le programme d'estimation
    ! read(5,*)donne_reel !dit si 1 a la question precedente dit s'il sagit du jeux de donnees reel (1) ou non (0)
    ! read(5,*)gener_only ! dit si on voudrait seulement generer les donnees(1) ou generer et faire des simulation(0)
    ! read(5,*)kapa_use ! dit si on utilise un kappa a chaque generation de donnee (1) ou le premier kappa pour tous les jeux de donnees(0)
    ! read(5,*)decoup_simul ! dans le cas où l'on a decoupe les simulations en plusieurs paquets, donne le nombre de generation de donnees a ne pas considerer avant d'engager les simulations. ceci empêche de reproduire les meme 
    !          jeux de donnees pour tous les paquets de simulation. vaut 0 si pas de decoupage pevu sinon pour chaque jeux de simulation mettre cette valeur a jour. Exp si 10 paquets de simul pour un total de 100, on affecte 0 
    !          pour le premier paquet, 10 pour le second, 20 pour le 3 ieme, ... 90 pour le 10ieme
    ! read(5,*)aleatoire    ! dit si on reinitialise la generation des nombre aleatoire avec un environnement different a chaque appel (1) ou non(O).En cas de generation differente, on utilise l'horloge (heure) de l'ordinateur comme graine. Dans ce cas, il n'est pas possible de reproduire les donnees simulees
    ! read(5,*)nbre_sim    ! dans le cas ou aleatoire=1, cette variable indique le nombre de generation qui vont etre faites
    ! read(5,*)graine    ! dans le cas ou l'on voudrait avoir la possibilite de reproduire les donnees generees alors on met la variable aleatoire=0 et on donne dans cette variable la graine a utiliser pour la generation

    sigmast_vrai=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
    thetast_vrai=rsqrt_theta*dsqrt(theta2)*dsqrt(theta2_t)
    gammast_vrai=rsqrt_gamma_ui*dsqrt(gamma_ui)*dsqrt(gamma_uit)
    
    if(nsim_node(11) == 3) then ! si joint frailty copula, alors on ajout les covariable aux jeux de donnees
      allocate(d_S(nsujet*n_sim,6 +size(vbetast,1)-1),d_T(ng*n_sim,6+size(vbetast,1)-1)) 
    else
      allocate(d_S(nsujet*n_sim,6),d_T(ng*n_sim,6)) 
    endif
    
    if(une_donnee==1) then
        ! on recupere le jeu de donnees reelles
        d_S=donnees
        d_T=death
    else
        ! on recupere le jeu de donnees reelles
    endif
    
    moy_theta=0.d0
    moy_theta=0.d0
    moy_sigmas=0.d0
    moy_sigmat=0.d0
    moy_rho=0.d0
    moy_gamma=0.d0
    moy_trt=0.d0
    moy_dec=0.d0
    moy_pros=0.d0
    moy_cens=0.d0
    moy_theta_est=0.d0
    moy_se_theta=0.d0
    moy_sigmas_est=0.d0
    moy_sigmat_est=0.d0
    moy_sigmast_est=0.d0
    moy_eta=0.d0
    moy_se_eta=0.d0
    moy_se_sigmas=0.d0
    moy_se_sigmat=0.d0
    moy_se_sigmast=0.d0
    moy_betaS=0.d0
    moy_betaS_se=0.d0
    moy_betaT=0.d0
    moy_betaT_se=0.d0
    moy_thetat=0.d0
    moy_rho_wij=0.d0
    moy_gammat=0.d0
    moy_rho_ui=0.d0
    taux_couvertureS=0.d0
    taux_couvertureT=0.d0
    taux_couverture_theta=0.d0
    taux_couverture_eta=0.d0
    taux_couverture_sigmas=0.d0
    taux_couverture_sigmat=0.d0
    taux_couverture_sigmast=0.d0
    nbre_rejet=0.d0
    moy_ni=0
    moy_gamma_est=0.d0
    moy_se_gamma=0.d0
    taux_couverture_gamma=0.d0
    moy_thetat_est=0.d0
    moy_se_thetat=0.d0
    taux_couverture_thetat=0.d0
    moy_thetast_est=0.d0
    moy_se_thetast=0.d0
    taux_couverture_thetast=0.d0
    moy_gammat_est=0.d0
    moy_se_gammat=0.d0
    taux_couverture_gammat=0.d0
    moy_gammast_est=0.d0
    moy_se_gammast=0.d0
    taux_couverture_gammast=0.d0
    moy_alpha=0.d0
    moy_se_alpha=0.d0
    taux_couverture_alpha=0.d0
    moy_R2_trial=0.d0
    moy_se_R2_trial=0.d0
    moy_bi_R2_trial=0.d0
    moy_bs_R2_trial=0.d0
    taux_couverture_R2_trial=0.d0
    moy_kendal_11=0.d0
    moy_kendal_10=0.d0
    moy_kendal_01=0.d0
    moy_kendal_00=0.d0
    debut_exe=1  ! pour dire on peut commencer l'execution pour ce processus
    indice_sim_proc=1 ! pour indicer le table des parametres estimes MPI
    indicCP=1 ! pour indicer le vecteur des IC par bootstrap
    moy_param_cop = 0.d0
    moy_se_param_cop = 0.d0
    taux_couverture_param_cop = 0.d0
    taux_couverture_tauk = 0.d0
    
    allocate(parametre_empirique(n_sim,12),parametre_empirique_NC(n_sim,12)) !parametre_empirique: contient les parametres empirique pour chaque simul (trt,theta,surrrogate,deces,vars,vart,rhost,gamma) 
    if(nsim_node(8)==0)then !model conjoint surrogate classique
        allocate(parametre_estimes(n_sim,8)) !parametres estimes: contient les parametres estimes(theta_chap+sd,zeta+sd,beta_s+sd,beta8t+sd)
        !allocate(param_estimes(n_sim,24))
    else
        if(nsim_node(8)==1)then ! modele avec les effets aleatoires partages
            allocate(parametre_estimes(n_sim,24)) !parametres estimes: contient les parametres estimes(theta_chap+sd,zeta+sd,beta_s+sd,beta8t+sd,sigma_s+sd,sigma_t+sd,sigmast+sd,gamma_ui+sd,alpha_ui+sd, R2 reduit et sd, taux de kendall)
            allocate(parametre_estimes_MPI(n_sim,24))! contient les parametres estimes par chaque processus dans MPI
            allocate(parametre_estimes_MPI_T(n_sim,24)) ! contient tous les parametres, de tous les processus
        endif
        
        if(nsim_node(8)==2)then ! modele complet avec les effets correles. on a 11 parametre a estimer dans le pire des cas avec leur SE
            allocate(parametre_estimes(n_sim,32)) !parametres estimes: contient les parametres estimes(theta_chap+sd,zeta+sd,beta_s+sd,beta8t+sd,sigma_s+sd,sigma_t+sd,sigmast+sd,gamma_ui+sd,alpha_ui+sd, R2 reduit et sd, taux de kendall)
            allocate(parametre_estimes_MPI(n_sim,32))! contient les parametres estimes par chaque processus dans MPI
            allocate(parametre_estimes_MPI_T(n_sim,32)) ! contient tous les parametres, de tous les processus
        endif
        
        if(nsim_node(8)==3)then ! joint frailty-copula model
            allocate(parametre_estimes(n_sim,25 + nva -2)) !parametres estimes: contient les parametres estimes(theta_copula+sd,zeta+sd,beta_s+sd,beta8t+sd,sigma_s+sd,sigma_t+sd,sigmast+sd,gamma_ui+sd,alpha_ui+sd, R2 reduit et sd, taux de kendall + sd) + variables explicative supplementaires
            allocate(parametre_estimes_MPI(n_sim,25 + nva -2))! contient les parametres estimes par chaque processus dans MPI
            allocate(parametre_estimes_MPI_T(n_sim,25 + nva -2)) ! contient tous les parametres, de tous les processus
        endif
        
        !allocate(param_estimes(n_sim,size(parametre_estimes,2)))
    endif
    
    
    
    allocate(tab_rejet(n_sim),tab_var_theta(n_sim),tab_var_sigma(n_sim,8))
    !allocate(Vect_sim_MC(nsim_node(1),1))
    parametre_empirique=0.d0
    parametre_estimes=0.d0
    parametre_estimes_MPI=0.d0
    parametre_estimes_MPI_T=0.d0
    tab_rejet=0
    allocate(kappa(n_sim,2))
    !if(kapa_use.eq.1) then ! on usilise un nouveau kappa pour chaque jeu
        !on recupere tous les kappa
    !    do i=1,n_sim
    !        read(15,*)kappa(i,1),kappa(i,2)
    !    enddo
    !    close(15)
    !    open(15,file=kapa)
    !endif
    incre_kappa=1 ! pour contenir le nombre de kappa ayant permis la convergence
    statut_kappa1=0
    s_i=1
    !if(une_donnee.ne.1) then !si on doit generer les donnees alors
    !    s_i=1
    !else ! si les donnees proviennent d'un fichier de donnee
    !    s_i=decoup_simul+1
    !endif
    
    !========= gestion de la variable n_sim pour le processus courant===========
    !=========autrement dit, du nombre de simulation par processeur=============
    n_sim_total=n_sim ! nombre de simulations a effectuer au total
    n_sim=INT(n_sim/nb_processus)
    suplement=n_sim_total-n_sim*nb_processus ! donne le nombre de simulation a partager entre les premiers processus seulement
    
    ! remplissage du table du nombre de simulation a effectuer par processus
    !!print*,nb_processus,n_sim,suplement
    allocate(tableNsim(nb_processus))
    tableNsim(1:nb_processus)=n_sim
    !!print*,tableNsim
    tableNsim(1:suplement)=n_sim+1 ! tous les essais jusqu'au rang supplement-1 recoivent une tâche supplementaire a realiser
    n_sim=tableNsim(rang+1) ! nombre de simulations a effectuer par le processus courant
    
    if (rang==0) then
       ! if(rang_proc==0) !print*, "tableNsim=",tableNsim
        init_i=1 ! ce processus commence a la premiere simulation
    else
        init_i=sum(tableNsim(1:rang))+1 ! ce processus commence a la simulation qui respecte son ordre et doit s'arreter au nombre de simultation dont il a le droit d'executer
    endif
    
    max_i=init_i+tableNsim(rang+1)-1!rang maximale de la simulation a executer (-1 car on a deja incrementer init_i de 1)
    !!call MPI_ABORT(MPI_COMM_WORLD,erreur,code)! on stop tous les programmes appartenant au communicateur code, equivalent de l'instruction stop en sequantiel
    incre_decoup=0
    ! !print*,"max_i=",max_i,tableNsim(rang+1)
    ! stop
    
    ! open(445,file='ResultPreliminair_model_complet_simul.txt') ! juste pour voir ce que xa donne
    do while (s_i<=n_sim_total)
        !Write(aaa,'(i3)') s_i ! instruction pour convertir un entier en chaine de caractere (3 caracteres)
        !aaa="sde"//aaa ! instruction pour concatener deux chaines de caracteres
        !!print*,"aaa=",aaa
        !stop
        if(une_donnee.eq.1) then ! dans ce cas on utilise un seul jeu de donnees pour les simul, de preference les donnees reel
            !if(s_i.eq.seed_) then
            !    call simulation(don_simul,don_simulS1,ind_temp,n_col,theta,eta,betas,alpha,betat,p,prop_i,lambdas,nus,&
            !                    lambdat,nut,mode_cens,temps_cens,cens0,n_essai,n_obs,&
            !                    rsqrt,sigma_s,sigma_t,weib,frailty_cor,affiche_stat)    
            !endif
            ! indice des donnees a utiliser pour la simulation encours
            i_min=s_i*nsujet-nsujet+1
            i_max=s_i*nsujet
            i_min_t=s_i*ng-ng+1
            i_max_t=s_i*ng
            if(donne_reel==1) then
                don_simulS1(:,initTime1)=d_S(i_min:i_max,4)
                don_simulS1(:,timeS1)=d_S(i_min:i_max,5)
                don_simulS1(:,statusS1)=d_S(i_min:i_max,6)
                don_simulS1(:,trialref1)=d_S(i_min:i_max,1)
                don_simulS1(:,Patienref1)=d_S(i_min:i_max,2)
                don_simulS1(:,trt1)=d_S(i_min:i_max,3)
                don_simul(:,initTime1)=d_T(i_min_t:i_max_t,4)
                don_simul(:,timeT1)=d_T(i_min_t:i_max_t,5)
                don_simul(:,statusT1)=d_T(i_min_t:i_max_t,6)
                don_simul(:,trialref1)=d_T(i_min_t:i_max_t,1)
                don_simul(:,Patienref1)=d_T(i_min_t:i_max_t,2)
                don_simul(:,trt1)=d_T(i_min_t:i_max_t,3)    
            else ! alors la position des variables n'est plus la meme
                don_simulS1(:,initTime1)=d_S(i_min:i_max,1)
                don_simulS1(:,timeS1)=d_S(i_min:i_max,2)
                don_simulS1(:,statusS1)=d_S(i_min:i_max,3)
                don_simulS1(:,trialref1)=d_S(i_min:i_max,4)
                don_simulS1(:,Patienref1)=d_S(i_min:i_max,5)
                don_simulS1(:,trt1)=d_S(i_min:i_max,6)
                don_simul(:,initTime1)=d_T(i_min_t:i_max_t,1)
                don_simul(:,timeT1)=d_T(i_min_t:i_max_t,2)
                don_simul(:,statusT1)=d_T(i_min_t:i_max_t,3)
                don_simul(:,trialref1)=d_T(i_min_t:i_max_t,4)
                don_simul(:,Patienref1)=d_T(i_min_t:i_max_t,5)
                don_simul(:,trt1)=d_T(i_min_t:i_max_t,6)
            endif
            
            ! j'ajoute les autres variables a la fin
            do i = 2,nbrevar(3)
                if(filtre(i).eq.1)then
                    don_simulS1(:,size(don_simulS1,2) - nbrevar(3) + i - 1)=d_S(i_min:i_max,size(don_simulS1,2) &
                    - nbrevar(3) + i - 1)
                endif
                if(filtre2(i).eq.1)then
                    don_simul(:,size(don_simul,2)- nbrevar(3) + i - 1)=d_T(i_min_t:i_max_t,size(don_simul,2)- nbrevar(3) + i - 1)
                endif
            enddo
            ! on met à jour le nombre d'essais
            20041 continue
            
            ! pour la gestion des paquets de simulation, avance dans le fichier des kappas pour se placer au bon endroit
            if(incre_decoup<decoup_simul) then !incre_decoup<decoup_simul c'est pour gerer le cas des simpulations par paquet
                ax1 = vect_kappa(indice_kapa,1)
                ax2 = vect_kappa(indice_kapa,2)
                indice_kapa = indice_kapa +1
                incre_decoup=incre_decoup+1
                goto 20041 ! pour etre sur qu'on n'utilise pas les kappas des jeux de donnee a ne pas considerer dans ce paquet de simulation
            endif

            
            if((s_i<init_i).or.s_i>max_i) then 
                debut_exe=0 ! pour dire le processus ne considere pas ce jeu de donnee
            else
                debut_exe=1 !jeux de donnees a considerer pour le processus courant
            endif
               
            !!print*,d_T(1:20,1)
            !!print*,don_simul(:,trialref1)
            allocate(tableEssai(ng))
            
            tableEssai=table_essai(nint(don_simul(:,trialref1)))
            
            n_essai=0
            do i=1,ng
                if(tableEssai(i).ne.0)then
                    n_essai=n_essai+1
                endif
            enddo
            ntrials=n_essai
            deallocate(tableEssai)
            ind_temp=ng
            theta=theta2
            !stop
        else            
            theta=0.d0
            !call intpr("nsim_node(8) =", -1, nsim_node(8), 1)
            if(nsim_node(8)==0) then
                indice_seed=0
                call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                10 continue
                call generation_Gamma(don_simul,don_simulS1,ng,n_col,logNormal,affiche_stat,theta,&
                    ng,ver,alpha,cens0,temps_cens,gamma1,gamma2,theta2,lambdas,nus,lambdat,nut,betas,betat)
                indice_seed=indice_seed+1
                if(indice_seed<s_i) goto 10
            else
                !Je suppose que le processus de rang n simule d'abor n-1 jeux de donnees a ne pas utiliser et n'utilise finalement que la donnees de son rang. ceci permet de considerer les donnees comme
                ! si elles etaient toutes simulees par un seul processus. sinon les premiers jeux de donnees auront tendance a avoir les meme distributions. par ailleurs bien faire la repartition des 
                ! donnees a utiliser par chaque processus car pour le moment je suppose un jeu de donnees par processus donc autant de coeur que de replication a faire
                !!print*,"Je suis le processus ",rang," parmi ",nb_processus
                ! do k=1,rang
                ! !jeux de donnees a surssoire car ont deja ete utilisees par les processus de rang inrerieure
                !incre_decoup=0
                !call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                                
                indice_seed=0
                call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                11 continue
                
                if(nsim_node(11)==1) then !modele avec effets aleatoires partages
                   
                    call Generation_surrogate(don_simul,don_simulS1,ng,n_col,logNormal,affiche_stat,theta,&
                        ng,ver,alpha,cens0,temps_cens,gamma1,gamma2,theta2,lambdas,nus,lambdat,nut,vbetas,vbetat,&
                        n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma_ui,alpha_ui,frailt_base,pfs)    
                endif
                
                ! call dblepr("sigma_s", -1, sigma_s, 1)    
                ! call dblepr("sigma_s", -1, sigma_s, 1)    
                ! call dblepr("rsqrt", -1, rsqrt, 1)                    
                
                if(nsim_node(11)==3) then ! joint frailty copula model
                    call Generation_surrogate_copula(don_simultamp,don_simulStamp,ng,n_col,logNormal,affiche_stat,theta,&
                        ng,ver,alpha,cens0,temps_cens,gamma1,gamma2,theta2,lambdas,nus,lambdat,nut,vbetas,vbetat,&
                        n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma_ui,alpha_ui,frailt_base,thetacopule, filtre,&
                        filtre2,pfs)  
                    don_simul(:,1:4) = don_simultamp(:,1:4)
                    don_simul(:,5) = 0.d0 ! on le met a 0 car je ne prends pas en compte les w_ij au moment de generation avec les copule. ducoup matricce avec -1 colone, par rapport a la generation a partir du modele joint surrogate
                    don_simul(:,6:size(don_simul,2)) = don_simultamp(:,5:(size(don_simultamp,2)-1)) ! -1 car la derniere colonne n'est pas rempli
                    don_simulS1(:,1:4) = don_simulStamp(:,1:4)
                    don_simulS1(:,5) = 0.d0
                    don_simulS1(:,6:size(don_simulS1,2)) = don_simulStamp(:,5:(size(don_simulStamp,2)-1)) ! -1 car la derniere colonne n'est pas rempli

                endif
                
                if(nsim_node(11)==2) then ! modele complet avec effets aleatoires correles
                    ! call Generation_surrogate_complet(don_simul,don_simulS1,ng,n_col,logNormal,affiche_stat,theta,&
                        ! ng,ver,alpha,cens0,temps_cens,gamma1,theta2,lambdas,nus,lambdat,nut,betas,betat,&
                        ! n_essai,rsqrt,sigma_s,sigma_t,prop_i,gamma_ui,theta2_t,rsqrt_theta,&
                        ! gamma_uit,rsqrt_gamma_ui,indice_gamma_st)    
                endif
     
        !==================================================================
        !==================22/05/2019==== ================================= 
        !==================================================================
        
         ! je commente ces lignes qui ne sont plus utilisees, car le vecteur des kappas produit correspond exactement au nombre de simulation a faire dans le pacquet ourant.
         ! ces lignes seraient a decommenter si j'utilisais mon programme fortran qui prend en entree le fichier de tous les kappas issus de la falidation croisee sur R         

                ! 2004    continue                
                ! ! pour la gestion des paquets de simulation, on genere sans utiliser les decoup_simul premiers jeux de donnees
                ! if(incre_decoup<decoup_simul) then !incre_decoup<decoup_simul c'est pour gerer le cas des simpulations par paquet
                    ! ! read(15,*)ax1,ax2 ! on incremente les kappa pour etre sur que pour les jeux donnees a utiliser on utilise le kappa correspondant
                    ! ax1 = vect_kappa(indice_kapa,1)
                    ! ax2 = vect_kappa(indice_kapa,2)
                    ! indice_kapa = indice_kapa+1
                    ! incre_decoup=incre_decoup+1
                    ! goto 2004 ! pour etre sur qu'on n'utilise pas le jeu de donnee
                ! endif
        !==================================================================
        !==================22/05/2019 Fin commentaires==== ================================= 
        !==================================================================
                
                ! si je suis avec des paquets de simulation, alors quand je suis la, je genere les premier jeu de donnees a ne pas considere, apres reinitialisation de l'environnement de generation
                indice_seed=indice_seed+1
                if(indice_seed<(s_i+decoup_simul)) goto 11
                
                if((s_i<init_i).or.s_i>max_i) then 
                    debut_exe=0 ! pour dire le processus ne considere pas ce jeu de donnee
                else
                    debut_exe=1 !jeux de donnees a considerer pour le processus courant
                endif
                
                ! jeux de donnees a considerer pour le processus courant
                
                !!print*,"=============suis la s_i===================",s_i,"rang=",rang
                ! call Generation_surrogate(don_simul,don_simulS1,ng,n_col,logNormal,affiche_stat,theta,&
                        ! ng,ver,alpha,cens0,temps_cens,gamma1,gamma2,theta2,lambdas,nus,lambdat,nut,betas,betat,&
                        ! n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma_ui,alpha_ui,frailt_base)                
            endif
            ! ind_temp=ng
            !sauvegarde de tous le jeu de donnees simulees
            do k=1,ng
                ! !write(18,*)don_simul(k,:)
                ! !write(19,*)don_simulS1(k,:)
            enddo
        endif
        ! do i=1,10
            ! !print*,"entête don_simul(1:5)",don_simul(i,1:5)
        ! end do
        ! do i=1,10
            ! !print*,"entête don_simul (6-10)",don_simul(i,6:10)
        ! end do
               
        !!print*,don_simul(:,:)    
        !!print*,don_simulS1(1:5,)        
        !stop    
            
! call dblepr("don_simultamp", -1, dble(don_simultamp(1,:)), size(don_simultamp,2))        
! call dblepr("don_simulStamp", -1, dble(don_simulStamp(1,:)), size(don_simulStamp,2))        
! call dblepr("don_simulS1", -1, dble(don_simulS1(1,:)), size(don_simulS1,2))        

        
        
        ind_temp=ng
        allocate(don_simulS(ind_temp,n_col))
        allocate(tableEssai(n_essai))
        tableEssai=0
        don_simulS=don_simulS1(1:ind_temp,:) ! on recupere les significatives du tableaux don_simulS1
        nsujet=ind_temp ! on met à jour le nombre d'observations sur les surrogates
        !!print*,sum(don_simulS(:,statusS1))
    !do i=1,10
    !!print*,"entête don_simul(1:5)",don_simul(i,1:5)
    !end do
    !do i=1,10
    !!print*,"entête don_simul (6-10)",don_simul(i,6:10)
    !end do
    
    !******************************************************************
    !*************Debut de la lecture des donnees**********************
    !******************************************************************
   ! if(rang_proc==0) !write(*,*)
   ! if(rang_proc==0) !write(*,*)'****** Debut de la lecture des donnees******'
    !if(s_i.eq.seed_) then
    !    open(9,file=donnees)
    !    open(10,file=death)
        !!write(9,*)"initTime","timeS","statusS","trialref","Patienref","trt"    
        !!write(10,*)"initTime","timeT","statusT","trialref","Patienref","trt" 
    !endif
    !!print*,ax1,ax2
    
    !---------------------->  Deces
    trials=0 ! vecteur contenant le nombre de sujet par etude
    vaxdc=0.d0
    Deces=0
    randomisation=0
    surrogate=0
    nigs=0
    cdcs=0
    nigts=0
    cdcts=0
  
    do i= 1,ng !sur les groupes uniquement (= sujet)
        !read(10,*)
        tt0dc(i)=don_simul(i,initTime1)/revision_echelle ! temps initial
        tt1dc(i)=don_simul(i,timeT1)/revision_echelle! temps de deces
        icdc(i)=Int(don_simul(i,statusT1)) ! delta_star
        pourtrial(i)=don_simul(i,trialref1) !indice de l'essai
        groupe(i)=don_simul(i,Patienref1) ! numero de l'individu
        vaxdct(i,1)=don_simul(i,trt1) ! vecteur des variables explicatives
        tableEssai(pourtrial(i))=tableEssai(pourtrial(i))+1
        ! j'ajoute les autres variables a la fin
        
        if(type_joint_estim == 3) then! joint frailty copula
            if(ver > 1) then ! I add the rest of covariates
                vaxdct(i,2:ver) = don_simul(i,(size(don_simul,2) - ver +2):size(don_simul,2))
            endif
        endif
        
        !ecriture des donnees dans le fichier (juste pour le premier jeux de donnee)
        ! on sauvegarde seulement si on est dans la simple generation des donnees
        !!print*,seed_
        if(une_donnee==0 .and. gener_only==1)then 
            if(seed_.eq.0) then
                !!print*,"donnees-death",donnees,death
              !  if(rang_proc==0) !write(16,*)tt0dc(i),tt1dc(i),icdc(i),pourtrial(i),groupe(i),(vaxdct(i,j),j=1,ver)
            endif
                    !stop
            if(s_i.eq.seed_) then
                !!print*,"donnees-death",donnees,death
               ! if(rang_proc==0) !write(16,*)tt0dc(i),tt1dc(i),icdc(i),pourtrial(i),groupe(i),(vaxdct(i,j),j=1,ver)
            endif
        endif
        k=0
        do j=1,ver
            if (filtre2(j).eq.1)then
                k=k+1
                vaxdc(i,k)=vaxdct(i,j) ! on maintient parmi les variables explicatives de deces celles indiquees dans le filtre
            end if
        end do
        ! elements de statistique
        trials(int(pourtrial(i)))=trials(int(pourtrial(i)))+1 ! on incremente le nombre de personnes dans letude
        Deces(int(icdc(i)))=Deces(int(icdc(i)))+1 ! compte le nombre de decedes et de vivant
        randomisation(int(vaxdct(i,1)))=randomisation(int(vaxdct(i,1)))+1 ! compte le nombre de personne sous traitements et non traites
        if(icdc(i).eq.1) then !on incremente le nombre de personnes decedees par essai
            cdcs(int(pourtrial(i)))=cdcs(int(pourtrial(i)))+1 ! nombre de cedes par essai
            cdcts(int(pourtrial(i)))=cdcts(int(pourtrial(i)))+int(vaxdct(i,1))! nombre de deces traites par essai
        endif
    end do
   
    ! !print*,"trials",trials
    ! if(s_i==2)    then
        ! stop
    ! endif
    deallocate(groupe)
    allocate(groupe(nsujet))
    !close(10)
    
!---------------------->  Surrogate
    vax=0.d0
    
! call intpr("suis danc funcpan n_col=", -1, n_col, 1)
! call dblepr("suis danc funcpan don_simulS=", -1, dble(don_simulS(1,:)), size(don_simulS,2))

    do i = 1,nsujet     !sur les observations
        k=0
        !read(9,*)tt0(i),tt1(i),ic(i),pourtrial(i),groupe(i),(vaxt(i,j),j=1,ver)
        tt0(i)=don_simulS(i,initTime1)/revision_echelle ! temps initial
        tt1(i)=don_simulS(i,timeS1)/revision_echelle! temps de deces
        ic(i)=Int(don_simulS(i,statusS1)) ! delta_star
        pourtrial(i)=don_simulS(i,trialref1) !indice de l'essai
        groupe(i)=don_simulS(i,Patienref1) ! numero de l'individu
        vaxt(i,1)=don_simulS(i,trt1) ! vecteur des variables explicatives
        ! j'ajoute les autres variables a la fin
        
        if(type_joint_estim == 3) then! joint frailty copula
            if(ver > 1) then ! I add the rest of covariates
                vaxt(i,2:ver) = don_simulS(i,(size(don_simulS,2) - ver +2):size(don_simulS,2))
            endif
        endif
                
        !ecriture des donnees dans le fichier (juste pour le premier jeux de donnee)
        ! on sauvegarde seulement si on est dans la simple generation des donnees
        
        if(une_donnee==0  .and. gener_only==1)then
            if(seed_==0) then
             !   if(rang_proc==0) !write(9,*)tt0(i),tt1(i),ic(i),pourtrial(i),groupe(i),(vaxt(i,j),j=1,ver)
            endif
            
            if(s_i==seed_) then
             !   if(rang_proc==0) !write(9,*)tt0(i),tt1(i),ic(i),pourtrial(i),groupe(i),(vaxt(i,j),j=1,ver)
            endif
        endif
            
        do j=1,ver
            if (filtre(j).eq.1)then    
                k=k+1
                vax(i,k)=vaxt(i,j) ! on maintient parmi les variables explicatives de surrogacy celles indiquees dans le filtre
            end if
        end do
        surrogate(int(ic(i)))=surrogate(int(ic(i)))+1 ! compte le nombre de personne avec et sans progression
        if(ic(i).eq.1) then !on incremente le nombre de personnes avec une progression par essai
            nigs(int(pourtrial(i)))=nigs(int(pourtrial(i)))+1 ! nbre personnes avec une progression
            nigts(int(pourtrial(i)))=nigts(int(pourtrial(i)))+int(vaxt(i,1)) ! nbre personnes avec une progression traites
        endif
    end do
    !close(9)
    
    ! construction de la matrice du nombre d'observations
    nig_Ts(:,1)=nigs
    nig_Ts(:,2)=nigts
    cdc_Ts(:,1)=cdcs
    cdc_Ts(:,2)=cdcts
    !!print*,"nigs=",nigs
    !!print*,"nigs=",nig_Ts(:,1)
    !!print*,"cdcts=",cdcts
    !!print*,"cdcts=",cdc_Ts(:,2)
    
    allocate(vdeces(deces(1)),vsurrogate(surrogate(1)))
    ! on recupere les dates de deces et de progression: cas des progression sans recidive
    
    k=1
    j=1
    do i=1,nsujet
        if(ic(i).eq.1) then
            vsurrogate(k)=int(tt1(i))
            k=k+1
        end if
        if(icdc(i).eq.1) then
            vdeces(j)=int(tt1dc(i))
            j=j+1
        end if
    end do 
    
    descripSurr(1)=minval(vsurrogate)
    descripSurr(2)=maxval(vsurrogate)
    descripSurr(3)=Median(vsurrogate, surrogate(1))
    
    descripDeces(1)=minval(vdeces)
    descripDeces(2)=maxval(vdeces)
    descripDeces(3)=Median(vdeces, deces(1))
    
    
    !====================================================================
    !                         quelques statistiques
    !====================================================================
    if(rang_proc==0) then
        !write(*,*)"===================================="
        !write(*,*)"*** Quelques statistiques:*********"
        !write(*,*)"===================================="    
        !write(*,*)'*** Nombre de personnes par essai:*********'
    endif
    ! do i=1,ntrials
        ! if(rang_proc==0) !print*,"Essai",i,"=",trials(i)
    ! end do
    
    !!write(4,*)"*** Nombre de personnes par essai:*********"
    ! do i=1,ntrials
        ! !!write(4,*)"Essai",i,"=",trials(i)
    ! end do
    
    if(rang_proc==0) then
        !write(*,*)trim("--->Nombre de personnes traitées globalement: "),randomisation(1),trim("soit(%)"),&
         !             (randomisation(1)/dble(nsujet))*100
        !write(*,*)trim("--->Nombre de personnes Décédées: "),deces(1),trim("soit(%)"),(deces(1)/dble(nsujet))*100
        !write(*,*)trim("--->Nombre de personnes avec une progression: "),surrogate(1),trim("soit(%)"),&
     !                 (surrogate(1)/dble(nsujet))*100
        !stop
        !write(*,*)trim("--->Description des temps de progression:")
        !write(*,*)trim("    Min="),descripSurr(1),trim("Max="),descripSurr(2),trim("Mediane="),descripSurr(3)
        !write(*,*)trim("--->Description des temps de deces:")
        !write(*,*)trim("    Min="),descripDeces(1),trim("Max="),descripDeces(2),trim("Mediane="),descripDeces(3)
        !write(*,*)
    endif
    

    
    allocate(paGH(ng,nsim_node(7)+1+nsim_node(7) + (nsim_node(7)*(nsim_node(7)-1))/2))
    paGH=0.d0
    !k0(1) = ax1 ! plus loin
    !k0(2) = ax2 ! plus loin
    ! call intpr("nz", -1, nz, 1)
    ! call intpr("nva", -1, nva, 1)
    ! call intpr("nparamfrail", -1, nparamfrail, 1)
    ! Nombre de parametres du modele
    select case(typeof)
        case(0)
            np = 2*(nz+2) + nva + nparamfrail !scl
        case(1)
            np = nbintervDC + nbintervR + nva + nparamfrail
        case(2)
            np = 2*nst + nva + effet + nparamfrail
    end select 
    
    if(rang_proc==0) then
        !write(*,*)'--->Methode dintegration:0=Monte carlo,1=quadrature,2=quadrature (frailties individuels)',&
        !'+MC(frailties essais),3=Laplace',nsim_node(4)
        !write(*,*)'--->nombre de simulation pour lapproximation de lintegrale par monte carlo',nsim_node(1)
        !write(*,*)'--->nombre de points de quadrature a utiliser en cas de quadrature',nsim_node(2)
        !write(*,*)'--->nombre de points de quadrature a utiliser en cas de non convergence',npoint2
        !write(*,*)'--->doit-on faire de ladaptative(1) ou de la pseudo-adaptative(0)',nsim_node(3)
        !write(*,*)'--->nombre total de parametres',np
    endif
    
    !allocate(b(np),H_hessOut(np,np),HIHOut(np,np))
    mtaille(1) = mt1 ! vecteur contenant la taille des vecteurs utilises pour ploter lees fonction de risque et de survie
    mtaille(2) = mt2
    mtaille(3) = mt11
    mtaille(4) = mt12
    
    ! paramete de la weibull a prendre en entree
    shape_weib(1) = paraweib(1)  ! pour le surrogate
    shape_weib(2) = paraweib(2)  ! pour le deces
    scale_weib(1) = paraweib(3)  ! pour le surrogate
    scale_weib(2) = paraweib(4)  ! pour le deces

    paratps = 0 ! indicateur de presence de variable temps-dependente

    MartinGales = 0 ! matrice des residus de martingale pour surrogate, deces, la predistion des fragilits, variance des fragilites predites

    ttU = 0.d0        ! en cas de censure par intervalle, borne superieure de l'intervalle
    
    !"=======Initialisation des parametres du modele a laide du modele conjoint classique============="
    ! transformer le main du modele conjoint initiale en une fonction qui retourne les betas puis lappeller ici pour linitialisation des betas
    !!print*,sum(nigs)
    !!print*,sum(cdcs)
    
    if(gener_only.eq.1) then ! ceci voudrait dire qu'on voudrait simplement generer les donnees pas d'estimation
      !if(rang_proc==0) !print*,"suis a la simulation",s_i
      goto 1000
    endif
 
    !==========Gestion des kappa====================
    !===============================================
    ! ici on considere soit le premier kappa pour toutes les simulation, soit un kappa par simul ou alors un kappa par simul avec possibilite en cas de non convergence
    ! de rentrer utiliser un des 4 premiers kappas ayant permis la convergence, dans ce cas, si le premier kappa ne permet pas la convergence, on prend le suivant 
    statut_kappa=0 ! pour dire que c'est un nouveau kappa, s'il permet la convergence on le retient pour les cas de non convergence
    ind_rech=1 ! s'incremente a la sortie du joint si le modele n'a pas converge
    Rech_kappa=0 ! controle la convergence du modele en faisant varie le kappa
    control_kappa=0 ! controle si on a deja divise les kappas dans le cas kappa_use=4
    ind_premier_kappa=0 ! qui est incremente a chaque fois que le kappa considere ne permet pas la convergence du premier jeu de donnee si on teste 5 kappa sur le premier jeu et on n'a pas de resultat, alors on continu avec le jeu de donnee suivant
    control=0 ! pour controler l'acces unique a la modification du kappa en cas de non convergence et lorsque le changement du nbre de point et des kappas issus de la cross validation ne marchent pas
    control2=0! pour s'assurer qu'on ne boucle pas lorsque ni le changement de kappa ni le nbre de point ne permet pas la convergence (en fait on doit relancer le modele une seule fois en changeant simultanement les deux parametre)
2001 continue
    if(kapa_use.ne.0) then ! on usilise un nouveau kappa pour chaque jeu de donnees
        !if(une_donnee.ne.1 .or. donne_reel .ne.1)read(15,*)ax1,ax2 ! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
        if(une_donnee.ne.1 .or. donne_reel .ne.1)then
            ax1 = vect_kappa(indice_kapa,1)
            ax2 = vect_kappa(indice_kapa,2)
            indice_kapa = indice_kapa+1 ! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
        end if
        k0(1)=ax1
        k0(2)=ax2
        k0_save=k0
        if(s_i==1)    k01_save=k0 ! on sauvegarde le premier kappa a utiliser plutard en cas de besoin
        !!print*,"kappa considere vaut:",k01_save
        if(statut_kappa1==1)then ! on utilise le kappa suivant (deuxieme)
            !read(15,*)ax1,ax2
            !k0(1)=ax1
            !k0(2)=ax2
            ! on revient au kappa courant
            ! close(15)
            ! open(15,file=kapa)
            indice_kapa = 1
            !read(15,*)ax1,ax2
            ax1 = vect_kappa(indice_kapa,1)
            ax2 = vect_kappa(indice_kapa,2)
            indice_kapa = indice_kapa + 1
            statut_kappa1=0
        endif
    else
        if((s_i.eq.1).or.((statut_kappa1==0).and.(ind_rech<=n_sim))) then !on considere le premier kappa qui marche pour toute les simul
            !if(une_donnee.ne.1 .or. donne_reel .ne.1) read(15,*)ax1,ax2  ! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
            if(une_donnee.ne.1 .or. donne_reel .ne.1)then
                ax1 = vect_kappa(indice_kapa,1)
                ax2 = vect_kappa(indice_kapa,2) 
                indice_kapa = indice_kapa +1! si les deux vallent un alors on utilise les kappas dournis dans le fichiers des parametres: joint_scl_simul
             end if
            k0(1)=ax1
            k0(2)=ax2
            k0_save=k0
            !statut_kappa1=1
            ind_rech=ind_rech+1
        endif
    endif
    
    
    !!print*,"kappa=",k0
2000 continue
2003 continue
    if (((Rech_kappa==1).and.((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4))) .or.(controlgoto==1)) then ! on recherhce parmi les kappas deja utilise s'il ya un ki permet la convergence
        if(((kapa_use.eq.2).or.(kapa_use.eq.4)).and.(controlgoto==0))then
            k0(1)=kappa(ind_rech-1,1)
            k0(2)=kappa(ind_rech-1,2)
        else
            controlgoto=0 ! ceci me permet de mieux gerer le goto 2003
            if((k0(1)-k0(2))>=1000)then
                k0(1)=k0(1)/10.d0
            else if((k0(1)-k0(2))<=-1000)then
                k0(2)=k0(2)/10.d0
            else
                k0(1)=k0(1)/10.d0
                k0(2)=k0(2)/10.d0
            endif
            !if(rang_proc==0) !print*,"division par 10 de kappa"
            goto 2002
        endif
        ! initialisation du vecteur b des parametres a l'aide des parametres de simulation
        !b=0.5d0
        !b(np-1)=1.0d0 !b_s
        !b(np)=0.7d0   !b_t
        !if(logNormal==1)then
        !    b(np-2)=0.5d0 !theta2
        !else
        !    b(np-2)=dsqrt(0.5d0) !theta
        !endif
    endif
2002 continue
    ! initialisation du vecteur b des parametres a l'aide des parametres de simulation
    b=0.5d0
    
    if(nsim_node(8)==0)then !model conjoint surrogate classique
        !if(logNormal==1)then
            b(np-2-nparamfrail+indice_eta+indice_theta)=theta !theta2
            !b(np-2-nparamfrail+indice_eta+indice_theta)=0.5d0 !theta2
        !else
        !    b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta) !theta
        !endif
    else
        if(nsim_node(8)==1)then !effets aleatoires partages
            b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=dsqrt(sigma_s)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)=dsqrt(sigma_t)
            if(indice_covST==1) then
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST)=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
            endif
            if(frailt_base==1) then
                if(indice_alpha==1) b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha)=&
                  alpha_ui
                b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+indice_gamma)=&
                  dsqrt(gamma_ui)
            endif
        endif
        
        if(nsim_node(8)==2)then !effets aleatoires correles
            b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=dsqrt(sigma_s)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)=dsqrt(sigma_t)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST)=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma)=dsqrt(gamma_ui)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+indice_theta_t)=&
                dsqrt(theta2_t)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+indice_theta_t+&
                indice_theta_st)=rsqrt_theta*dsqrt(theta)*dsqrt(theta2_t)
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+indice_theta_t+&
                indice_theta_st+indice_gamma_t)=dsqrt(gamma_uit)
            if(indice_gamma_st==1)b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                indice_theta_t+indice_theta_st+indice_gamma_t+indice_gamma_st)=rsqrt_gamma_ui*dsqrt(gamma_ui)*dsqrt(gamma_uit)
            
        endif
        
        if(nsim_node(8)==3)then !joint frailty-copula
            b(np-nva-nparamfrail+indice_varS)=dsqrt(sigma_s)
            b(np-nva-nparamfrail+indice_varS+indice_varT)=dsqrt(sigma_t)
            if(indice_covST==1) then
                b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST)=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
            endif
            if(frailt_base==1) then
                if(indice_alpha==1) b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_alpha)=&
                  alpha_ui
                b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_alpha+indice_gamma)=&
                  dsqrt(gamma_ui)
            endif
            if(copula_function == 1) b(np-nva) = dlog(thetacopule) ! claton: exp transform
            if(copula_function == 2) b(np-nva) = dsqrt(thetacopule)  ! Gumbel: choleschy transform
            b((np-nva + 1) : (np - nva + nva1)) = vbetas
            b((np-nva2 + 1) : np) = vbetat
            ! call dblepr(" gamma_ui joint =", -1, gamma_ui, 1)
        endif
    endif
    
    !affectation manuelle des parametres initiaux (choleschy obtenue de R)
    if(vrai_val_init==1)then ! on initialise avec les vrai valeurs
        if(theta>=1.d0)then
            theta_init=0.7d0 ! on fait ceci car pour un theta initialiser a 1 le modèle est difficile a estimer
        else
            theta_init=theta
        endif
        theta_init=theta
        sigma_ss_init=sigma_s
        sigma_tt_init=sigma_t
        gamma_init=gamma_ui
        alpha_init=alpha_ui
        if(indice_covST==1) then
            sigma_st_init=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
        endif
        betas_init=betas
        betat_init=betat
        gamma_init=gamma_ui
        alpha_init=alpha_ui
        zeta_init=eta
        
        ! pour le modele complet
        thetat_init=theta2_t
        thetast_init=rsqrt_theta*dsqrt(theta)*dsqrt(theta2_t)
        gammat_init=gamma_uit
        gammast_init=rsqrt_gamma_ui*dsqrt(gamma_ui)*dsqrt(gamma_uit)
        
        ! frailty-copula 
        thetacopula_init = thetacopule
        vbetas_intit = vbetas
        vbetat_intit = vbetat
        
    endif
    
    if(sigma_ss_init.eq.0.d0) then
        !if(rang_proc==0) !print*,"Attention: revoir la valeur initiale de sigma_SS"
    else
        if((sigma_tt_init-(sigma_st_init**2.d0)/sigma_ss_init)<0)then
        !    if(rang_proc==0) !print*,"Attention: revoir les valeurs initiales de la matrice sigma des effets aleatoire"
        endif
    endif
    
    if(nsim_node(8)==3)then !joint frailty-copula
        b(np-nva-nparamfrail+indice_varS)=dsqrt(sigma_ss_init)
        b(np-nva-nparamfrail+indice_varS+indice_varT)=dsqrt(sigma_tt_init-&
                (sigma_st_init**2.d0)/sigma_ss_init)
        if(indice_covST==1) then
            b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST)=sigma_st_init/dsqrt(sigma_ss_init)
        endif
        if(frailt_base==1) then
            if(indice_alpha==1) b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_alpha)=&
                alpha_init
                b(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_alpha+indice_gamma)=&
                    dsqrt(gamma_init)
            endif
        if(copula_function == 1)  b(np-nva) =  dlog(thetacopula_init)! claton: exp transform
        if(copula_function == 2) b(np-nva) = dsqrt(thetacopula_init)  ! Gumbel: choleschy transform
        b((np-nva + 1) : (np - nva + nva1)) = vbetas_intit
        b((np-nva2 + 1) : np) = vbetat_intit
    else 
        b(np-1)=betas_init !b_s
        b(np)=betat_init
        if(indice_eta==1)then
            b(np-2-nparamfrail+indice_eta)=zeta_init !eta
        endif
    endif
    
    ! si modele complet avec effets aleatoires partages
    if(nsim_node(8)==1)then
        b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=dsqrt(sigma_ss_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)=dsqrt(sigma_tt_init-&
                (sigma_st_init**2.d0)/sigma_ss_init)
        if(indice_covST==1) then
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST)=sigma_st_init/&
                dsqrt(sigma_ss_init)
         endif
        if(frailt_base==1) then
            if(indice_alpha==1) b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+&
                indice_alpha)=alpha_init
            b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+indice_gamma)=&
                dsqrt(gamma_init)
        endif
    endif
    
    ! si modele complet avec effets aleatoires correles
    if(nsim_node(8)==2)then
        b(np-2-nparamfrail+indice_eta+indice_theta)=dsqrt(theta_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS)=dsqrt(sigma_ss_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)=sigma_st_init/dsqrt(sigma_ss_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST)=dsqrt(sigma_tt_init-&
            (sigma_st_init**2.d0)/sigma_ss_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+indice_gamma)=&
            dsqrt(gamma_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+indice_theta_t)=&
                thetast_init/dsqrt(theta_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+indice_theta_t+&
                indice_theta_st)=dsqrt(thetat_init-(thetast_init**2.d0)/theta_init)
        b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+indice_theta_t+&
        indice_theta_st+indice_gamma_t)=gammast_init/dsqrt(gamma_init)
        if(indice_gamma_st==1)b(np-2-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                indice_theta_t+indice_theta_st+indice_gamma_t+indice_gamma_st)=dsqrt(gammat_init-(gammast_init**2.d0)/&
                gamma_init)
    endif
    
    if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3))then
        !if(rang_proc==0) !print*,nsim_node(2),np_save
        if((nsim_node(2).ne.np_save).and.(s_i==1))then
            !if(rang_proc==0) !print*,"le nombre de point de quadrature est passé de",np_save,"a:",nsim_node(2)
        endif
    endif

    ! !print*,"b",b(np-1-nparamfrail : np)
    ! stop
2005 continue

    EPS =EPS2    ! critere de convergence du modele
    !!print*,"b=",b
    !******************Fin lecture et initialisation des parametres***********************
    if(rang_proc==0) then
        !print*,""
        !print*,"======================================================================="
        !print*,"============debut de la maximisation de la vraisemblance==============="
        ! code_print = 0
        ! call intpr("============Begin of the maximization procedure...==============:", -1, code_print, 1)
        !print*,"======================================================================="
    endif
    call cpu_time(tp1)
    !if(s_i<=3) then
    !    !print*,"Je continue a l'itteration 2 de la boucle"
    !    goto 1000
    !endif
    if(debut_exe==1)then
        !if(rang_proc==0) !print*,"Je suis à la simulation",s_i,"kappa=",k0,"EPS=",EPS
    else
        goto 1002 ! le processus courant n'utilise pas ce jeu de donnee
    endif
    
    !if(s_i .ne.12) goto 1000
    !nsim_node(2)=20 ! on fait ceci juste pour debugger le programme
    !goto 101  ! on fait ceci juste pour debugger le programme
    !b=0.5d0
    !if(s_i<=1) goto 1000
    ! call intpr("np", -1, np, 1)
    ! call intpr("size(H_hessOut) avant appel joint:", -1, size(H_hessOut,1), 1)
    ! call intpr("size(H_hessOut) avant appel joint:", -1, size(H_hessOut,2), 1)
    ! call dblepr("H_hessOut(17,17) avant appel joint:", -1, H_hessOut(17,17), 1)
    ! j'incremente le kappa d'une constante en cas de precision pour la gestion des probleme de convergence
    k0(1) = k0(1) + ckappa(1)
    k0(2) = k0(2) + ckappa(2)
    
    if(affiche_itteration == 1) then
        !call dblepr("avant appel joint:ckappa", -1,ckappa , 2)
        ! call dblepr("avant appel joint:k0", -1,k0 , 2)
        ! call intpr("avant appel joint:nsujet", -1, nsujet, 1)
        ! call intpr("avant appel joint:ng", -1, ng, 1)
        ! call intpr("avant appel joint:ntrials", -1, ntrials, 1)
        ! call intpr("avant appel joint:nz", -1,nz , 1)
        ! call intpr("avant appel joint:nst", -1,nst, 1)
        !call intpr("avant appel joint:pourtrial", -1, pourtrial, nsujet)
        ! call intpr("avant appel joint:trials", -1,trials , ntrials)
        ! call intpr("personnes avec progression:nig_Ts(:,1)", -1, nig_Ts(:,1), ntrials)
        ! call intpr("personnes avec progression traitees:nig_Ts(:,2)", -1,nig_Ts(:,2) , ntrials)
        ! call intpr("personnes decedees:cdc_Ts(:,1)", -1, cdc_Ts(:,1), ntrials)
        ! call intpr("personnes decedees traitees:cdc_Ts(:,2)", -1, cdc_Ts(:,2), ntrials)
        ! call intpr("--->Nombre de personnes traitees globalement", -1, randomisation(1), 1)
        ! call dblepr("en %", -1, (randomisation(1)/dble(nsujet))*100, 1)
        ! call intpr("--->Nombre de personnes Decedees globalement", -1, deces(1), 1)
        ! call dblepr("en %", -1, (deces(1)/dble(nsujet))*100, 1)
        ! call intpr("--->Nombre de personnes avec une progression", -1, surrogate(1), 1)
        ! call dblepr("en %", -1, (surrogate(1)/dble(nsujet))*100, 1)
        ! call dblepr("temps de progression Min - Max - Mediane", -1, (/descripSurr(1), descripSurr(2), descripSurr(3)/) , 3)
        ! call dblepr("temps de deces Min - Max - Mediane", -1, (/descripDeces(1), descripDeces(2), descripDeces(3)/) , 3)
            
        ! call intpr("avant appel joint:icdc(1:10)", -1,icdc(1:10) , 10)
        ! call intpr("avant appel joint:nva1", -1,nva1 , 1)
        ! call intpr("avant appel joint:noVar1", -1,noVar1 , 1)
        ! call intpr("avant appel joint:noVar2", -1, noVar2, 1)
        ! call dblepr("avant appel joint:tt0", -1, tt0, size(tt0))
        ! call dblepr("avant appel joint:tt1(1:10)", -1, tt1(1:10), 10)
        ! call intpr("avant appel joint:ic(1:10)", -1,ic(1:10) , 10)
        ! call intpr("avant appel joint:groupe", -1,groupe , ntrials)

        !call dblepr("avant appel joint:tt0dc", -1,tt0dc , nsujet)
        ! call dblepr("avant appel joint:tt1dc(1:10)", -1,tt1dc(1:10) , 10)
        !call dblepr("avant appel joint:vax(:,1)", -1,vax(:,1) , nsujet)
        !call dblepr("avant appel joint:vaxdc(:,1)", -1,vaxdc(:,1) , nsujet)
        ! call dblepr("avant appel joint:paraweib", -1,paraweib , 4)
        ! call dblepr("avant appel joint:ziOut", -1, ziOut, nz+6)
        ! call dblepr("avant appel joint:EPS", -1,EPS , 3)
        ! call intpr("avant appel joint:ag", -1, ag, 1)
        ! call intpr("avant appel joint:maxiter", -1,maxiter , 1)
        ! call intpr("avant appel joint:np", -1,np , 1)
        ! call intpr("avant appel joint:typeof", -1,typeof , 1)
        ! call intpr("avant appel joint:equidistant", -1,equidistant , 1)
        ! call intpr("avant appel joint:nbintervR", -1, nbintervR, 1)
        ! call intpr("avant appel joint:nbintervDC", -1, nbintervDC, 1)
        ! call intpr("avant appel joint:mtaille", -1, mtaille, 4)
        ! call intpr("avant appel joint:logNormal", -1,logNormal , 1)
        ! call intpr("avant appel joint:nsim_node", -1,nsim_node , 13)
        ! call intpr("avant appel joint:indice_esti", -1, indice_esti, 4)
        ! call intpr("avant appel joint:indice_covST", -1,indice_covST , 1)
        ! call intpr("==============avant appel joint:param_weibull======", -1,param_weibull , 1)
    endif
                        
    Call joint_surrogate(nsujet,ng,ntrials,0,nz,nst,k0,tt0,tt1,ic,groupe,trials,pourtrial,nig_Ts,cdc_Ts,0, &
                        tt0dc,tt1dc,icdc,0.d0,0,nva1,vax,nva2,vaxdc,0.d0,noVar1,noVar2,ag,maxiter,np,b,H_hessOut,&
                        HIHOut,resOut,LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out,typeof,equidistant,&
                        nbintervR,nbintervDC,mtaille,ni,cpt,cpt_dc,ier,istop,paraweib,MartinGales,linearpred,&
                        linearpreddc,ziOut,time,timedc,0.d0 , 1 , 0 , ttU , logNormal , paratps , 0 , 0.d0 , 0.d0 , &
                        EPS,nsim_node,indice_esti,indice_covST,0.d0,param_weibull)
    ! call intpr("Nombre itteration:", -1, ni, 1)
    if (istop.eq.1) then
        call dblepr("voila le vecteur b des parametres", -1, b(2*(nz+2)+1:np), nva + nparamfrail)
    endif
!122     continue
    !if(s_i<5) nsim_node(2)=32 ! on fait ceci juste pour debugger le programme
    !nsim_node(2)=32 ! on fait ceci juste pour debugger le programme
    !!print*,"covar========================"
    if (istop.ne.1) then
        !call intpr("je suis la :", -1, ni, 1)
        ! if(nsim_node(3).ne.1) then ! cas ou on ne fait pas de l'adaptative
            if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3) .and. (kapa_use.eq.4))then !cas quadrature
                if ((control_kappa==3) .or. (control_kappa==4)) then ! on emet les vraies valeurs
                    if(np_save==npoint1) nsim_node(2)=npoint2
                    if(np_save==npoint2) nsim_node(2)=npoint1
                endif
                
                if((nsim_node(2).eq.np_save))then ! on change le nombre de point de quadrature avant de continuer
                    if(np_save==npoint1) nsim_node(2)=npoint2
                    if(np_save==npoint2) nsim_node(2)=npoint1
                    goto 2002 ! apres changement du nbre de point on relance l'estimation 
                else
                    if(control_kappa<4)then ! on vient de changer le nombre de point alors on revient
                    
                        if(control_kappa<2)then ! on va aller divise les kappas par 10, en maintenant change le nbre de point
                            !if(rang_proc==0) !print*,"division du kappa par 10 ou par 100"
                            control_kappa=control_kappa+1
                            controlgoto=1
                           ! goto 2003
                        else
                            if(control_kappa==2)then ! on remet les valeurs de kappa pour la suite
                                k0=k0_save
                            endif
                        endif
                        if(controlgoto==0)then
                            nsim_node(2)=np_save
                            if(control_kappa<4)then ! on relance le modele avec les kappas diviss par 10, sans change le nbre de point
                                !if(rang_proc==0) !print*,"division du kappa par 10 ou par 100 avec maintien du nombre de", &
                                 !   "points de quadrature"
                                control_kappa=control_kappa+1
                                controlgoto=1
                                !goto 2003
                            endif 
                        endif
                    else !on a deja changer aussi bien le nbre de point que le kappa et xa ne marche toujours pas, on va donc changer les deux
                        if(control_kappa==4)then !dans ce cas on change les kappa et les points de quadrature
                            if(np_save==npoint1) nsim_node(2)=npoint2
                            if(np_save==npoint2) nsim_node(2)=npoint1
                            
                            if(control2==0) then !on considere dans ce cas le dernier kappa/100 avec changement aussi du nombre de points de quadrature
                               ! if(rang_proc==0) !print*,"changement simultané du nbre de point et kappa division par 100"
                                control2=1 ! pour controler qu'on relance une seule fois le modele
                                goto 2002 !
                            endif
                            
                            if(s_i>1) then
                                if(control2==1)then ! on considere dans ce cas le premier kappa qui marche ou le premier kappa avec les points de quadrature de depart
                                    nsim_node(2)=np_save
                                    control2=2 
                               !     if(rang_proc==0) !print*,"changement du kappa (premier kappa qui marche ou tout premier si", &
                                 !       "aucun n'a marché jusqu'ici)"
                                    if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                                        k0(1)=kappa(1,1)
                                        k0(2)=kappa(1,2)
                                    else ! sinon on recupere le tout premier kappa
                                        k0=k01_save
                                    endif
                                    goto 2002 ! apres changement du nbre de point on relance l'estimation une seule fois
                                endif
                            endif
                            
                            if(s_i>1) then
                                if(control2==2)then ! on considere dans ce cas le premier kappa qui marche ou le premier kappa avec changement aussi du nombre de points de quadrature
                                    control2=3 
                                   ! if(rang_proc==0) !print*,"changement simultane du nbre de point et kappa (premier kappa", &
                                    !        "qui marche ou tout premier si aucun n'a marche jusqu'ici)"
                                    if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                                        k0(1)=kappa(1,1)
                                        k0(2)=kappa(1,2)
                                    else ! sinon on recupere le tout premier kappa
                                        k0=k01_save
                                    endif
                                    goto 2002 ! apres changement du nbre de point on relance l'estimation une seule fois
                                endif
                            endif
                        endif
                        
                    endif
                    if(controlgoto==0)then
                       ! if(rang_proc==0) !print*,"suis en fin là, control_kappa vaut",control_kappa
                    endif
                endif
                
                if(controlgoto==1) goto 2003
                
               ! if(rang_proc==0) !print*,"suis en fin la, aucune des deux combinaisons n'a marche"
                ! si j'arrive ici c'est qu'aucune methode n'a converge alors je change les valeurs initiales avec les dernieres estimees et je reprends le premier kappa qui marche s'il y en a ou simplement le premier kappa

    !101 continue    !pour le debuggaga        
                !!print*,"k0 avant vaut",k0
                if(ind_rech==1)then
                    if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                        k0(1)=kappa(1,1)
                        k0(2)=kappa(1,2)
                        ind_rech=2    ! vu que j'ai utilise le premier kappa qui a permis la convergence alors je passe au suivant
                    else ! sinon on recupere le tout premier kappa
                        k0=k01_save
                    endif
                    !!print*,"k0 pour cette nouvelle relance vaut",k0
                endif
                
                if(control==0) then ! on relance le modèle avec le premier kappa (converge ou pas) en changent les valeurs initiales aux dernieres estimations qui n'ont pas permises la convergence
                    ! control=1
                    ! !print*,"suis dans controle 1"
                    ! goto 2002    ! les parametres sont initialises avec les vraies valeurs
                !    if(rang_proc==0) !print*,"changement simultane du nbre de point et kappa (premier kappa", &
                  !          "qui marche ou tout premier si aucun n'a marche jusqu'ici) ainsi que des valeurs", & 
                  !          "initiales des parametres"
                    control=1
                    goto 2005    ! les parametres sont initialises avec les dernieres estimations n'ayant pas permis la convergence
                 else
                    if(control==1) then ! ici j'utilise le nombre de points de quadrature de base et je change juste le kappa et les valeurs initiales
                       ! if(rang_proc==0) !print*,"changement du kappa (premier kappa qui marche ou tout", &
                         !       "premier si aucun n'a marché jusqu'ici) ainsi que des valeurs initiales des parametres"
                        nsim_node(2)=np_save 
                        control=2
                        goto 2005    ! les parametres sont initialises avec les dernieres estimations n'ayant pas permis la convergence
                     endif
                endif
            endif
        ! endif
        
        if((nsim_node(4).eq.3 .or. nsim_node(4).eq.0 .or. ((nsim_node(4).eq.1 .or. nsim_node(4).eq.2).and.&
            nsim_node(3).eq.1)) .and. (kapa_use.eq.4))then !cas Monte carlo et Laplace ou alors pseudo adaptative
        
                if(control_kappa<2)then ! on vient de changer le nombre de point alors on revient
                
                    if(control_kappa<2)then ! on va aller divise les kappas par 10, en maintenant change le nbre de point
                      !  if(rang_proc==0) !print*,"division du kappa par 10 ou par 100"
                        control_kappa=control_kappa+1
                        controlgoto=1
                        !goto 2003
                    else
                        if(control_kappa==2)then ! on remet les valeurs de kappa pour la suite
                            k0=k0_save
                        endif
                    endif                
                else !on a deja changer le kappa et xa ne marche toujours pas
                    if(control_kappa==2)then !dans ce cas on change les kappa et les points de quadrature                        
                        if(s_i>1) then
                            if(control2==1)then ! on considere dans ce cas le premier kappa qui marche ou le premier kappa avec les points de quadrature de depart
                                control2=2 
                                ! if(rang_proc==0) !print*,"changement du kappa (premier kappa qui marche ou tout premier si ",&
                                    ! "aucun n'a marché jusqu'ici)"
                                if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                                    k0(1)=kappa(1,1)
                                    k0(2)=kappa(1,2)
                                else ! sinon on recupere le tout premier kappa
                                    k0=k01_save
                                endif
                                goto 2002 ! apres changement du kappa on relance l'estimation une seule fois
                            endif
                        endif
                    endif
                endif
            if(controlgoto==1) goto 2003
            ! if(rang_proc==0) !print*,"suis en fin là, control_kappa vaut",control_kappa            
            ! if(rang_proc==0) !print*,"suis en fin la, aucune des kappas n'a marche"
            ! si j'arrive ici c'est qu'aucune methode n'a converge alors je change les valeurs initiales avec les dernieres estimees et je reprends le premier kappa qui marche s'il y en a ou simplement le premier kappa

!101 continue    !pour le debuggaga        
            !!print*,"k0 avant vaut",k0
            if(ind_rech==1)then
                if(incre_kappa>=2) then ! on a aumoins un kappa qui marche on recupere le tout premier
                    k0(1)=kappa(1,1)
                    k0(2)=kappa(1,2)
                    ind_rech=2    ! vu que j'ai utilise le premier kappa qui a permis la convergence alors je passe au suivant
                else ! sinon on recupere le tout premier kappa
                    k0=k01_save
                endif
                !!print*,"k0 pour cette nouvelle relance vaut",k0
            endif
            
            if(control==0) then ! on relance le modèle avec le premier kappa (converge ou pas) en changent les valeurs initiales aux dernieres estimations qui n'ont pas permises la convergence
                ! control=1
                ! !print*,"suis dans controle 1"
                ! goto 2002    ! les parametres sont initialises avec les vraies valeurs
                ! if(rang_proc==0) !print*,"changement simultane du nbre de point et kappa (premier kappa", &
                    ! "qui marche ou tout premier si aucun n'a marche jusqu'ici) ainsi que des valeurs ",&
                    ! "initiales des parametres"
                control=1
                goto 2005    ! les parametres sont initialises avec les dernieres estimations n'ayant pas permis la convergence
            endif
        endif
        
        ! if(rang_proc==0) !print*,"je vais a la recherche des kappas qui on converges pour n'avoir rien trouver"
        ! if(rang_proc==0) !print*,"ind_rech=",ind_rech,"incre_kappa=",incre_kappa
        if((ind_rech<incre_kappa).and.((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4))) then ! on recherhce parmi les kappas deja utilise s'il ya un qui permet la convergence
            ! if(rang_proc==0) !print*,"test du kappa numero",ind_rech-1,"qui a deja eu à converger" 
            Rech_kappa=1
            ind_rech=ind_rech+1            
            statut_kappa=1 ! evite qu'on enregistre le kappa trouve s'il permet la convergence
            
            if(ind_rech==4)then ! on s'arrete et passe a la suivante simulation
                goto 1000
            endif
            goto 2000
        endif
        
        ! si le tout premier kappa ne converge pas
        !if((s_i==1).and.(n_sim.ne.1).and.((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4))) then 
        if((s_i==1).and.(ind_premier_kappa==4))then
            ! if(rang_proc==0) !print*, "Deja 5 kappas testes pour le premier jeu de donnees et pas toujouts", &
                ! "convergence, plus la peine de continuer"
        endif
        
        if((s_i==1).and.(n_sim.ne.1) .and.(ind_premier_kappa<4)) then 
            ! if(rang_proc==0) !print*,"suis la pour recherche kappa: on utilise le suivant"
            ind_premier_kappa=ind_premier_kappa+1
            statut_kappa1=1
            goto 2001 ! on utilise le kappa suivant
        endif

    else 
        if(kapa_use==0)then ! si c'est un seul kappa qu'il faut utiliser alors on le maintient et on continu
            statut_kappa1=1
        endif
    endif
        
1000 continue
    Rech_kappa=0
    ind_rech=1
    
    ! extraction des fragilites niveau essai pour evaluation
    k=1
    do i=1,n_essai
        !if(nsim_node(8)==2 .or. nsim_node(8)==3)then
        if(nsim_node(8)==2)then
            donnee_essai(i,:)=don_simulS(k,(/v_s1,v_t1,trialref1,u_i1,u_it/))
        else
            donnee_essai(i,:)=don_simulS(k,(/v_s1,v_t1,trialref1,u_i1/))
        endif
        k=k+tableEssai(i)
        !!print*,donnee_essai(i,4)
    enddo
    
    nsim_node(2)=np_save ! on remet le nombre de point de quadrature s'il y'a eu modification

    if (istop.ne.1) then
        ! if(rang_proc==0) !write(*,*)"ERREUR : LE MODELE N'A PAS CONVERGE. Istop=",istop
        !if(rang_proc==0) !write(5,*)"ERREUR : LE MODELE N'A PAS CONVERGE. Istop=",istop
        if(rang_proc==0) then
            code_print = 0
           !call intpr("===Model did not converged!!! please try to modified initial values or others parameters===:", -1,&
           !code_print, 1)
            !write(*,*) "===Model did not converged!!! please try to modified initial values or others parameters===:"
        endif
        ! si pas de convergence on simule un nouveau jeux de donnees
        nbre_rejet=nbre_rejet+1
        tab_rejet(nbre_rejet)=s_i
        
        ! controle de l'ecriture dans le fichier
        erreur_fichier=0
        996 continue
        !write(8,*,iostat=erreur_fichier)s_i
        if(erreur_fichier .ne.0) then
            ! if(rang_proc==0) !print*,"ATTENTION!! erreur d'ecriture de l'indice dans le fichier ",&
            ! "de sortie en cas de non convergence on insiste jusqu'a l'ecriture"
            goto 996
        endif
        
        ! parametres empiriques en cas de non convergence
        parametre_empirique_NC(nbre_rejet,1)=sum(don_simul(:,trt1))*100.d0/n_obs! personne traitees
        parametre_empirique_NC(nbre_rejet,2)=variance(don_simul(:,w_ij1))! theta empirique
        parametre_empirique_NC(nbre_rejet,3)=sum(don_simul(:,statusS1))*100.d0/n_obs! progression
        parametre_empirique_NC(nbre_rejet,4)=sum(don_simul(:,statusT1))*100.d0/n_obs! deces
        
        if(nsim_node(8).ne.0)then !model conjoint surrogate
            parametre_empirique_NC(nbre_rejet,5)=variance(donnee_essai(:,1))! sigma_s empirique
            parametre_empirique_NC(nbre_rejet,6)=variance(donnee_essai(:,2))! sigma_t empirique
            parametre_empirique_NC(nbre_rejet,7)=covariance(donnee_essai(:,1),donnee_essai(:,2))/&
                (dsqrt(variance(donnee_essai(:,1)))*dsqrt(variance(donnee_essai(:,2))))! correlation empirique: rho(x,y)=cov(x,y)/(sigma_x*sigma_y)
            parametre_empirique_NC(nbre_rejet,8)=variance(donnee_essai(:,4))! gamma S empirique
            
            if(nsim_node(8)==2)then
                parametre_empirique_NC(nbre_rejet,9)=variance(don_simul(:,w_ijt))! theta t empirique
                parametre_empirique_NC(nbre_rejet,10)=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/&
                (dsqrt(variance(don_simul(:,w_ij1)))*dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
                parametre_empirique_NC(nbre_rejet,11)=variance(donnee_essai(:,5))! gamma T empirique
                parametre_empirique_NC(nbre_rejet,12)=covariance(donnee_essai(:,4),donnee_essai(:,5))/&
                (dsqrt(variance(donnee_essai(:,4)))*dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST
            endif
        endif
            
        !!write(12,*)parametre_empirique_NC(nbre_rejet,:)
        if(gener_only.eq.1) then
            theta_sim=variance(don_simul(:,w_ij1))! variance des effets aleatoires simules
            sigmas_sim=variance(donnee_essai(:,1))! sigma_s empirique
            sigmat_sim=variance(donnee_essai(:,2))! sigma_t empirique
            rho_sim=covariance(donnee_essai(:,1),donnee_essai(:,2))/(dsqrt(variance(donnee_essai(:,1)))*&
                        dsqrt(variance(donnee_essai(:,2))))
            gamma_sim=variance(donnee_essai(:,4)) !gamma_ui empirique
            theta_simt=variance(don_simul(:,w_ijt))! variance des effets aleatoires simules T
            rho_sim_wij=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/(dsqrt(variance(don_simul(:,w_ij1)))*&
                        dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
            gamma_simt=variance(donnee_essai(:,5))! gamma T empirique
            rho_sim_ui=covariance(donnee_essai(:,4),donnee_essai(:,5))/(dsqrt(variance(donnee_essai(:,4)))*&
                        dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST
            !!print*,"rho=",rho_sim
            !stop
            tab_var_theta(s_i)=theta_sim
            tab_var_sigma(s_i,1)=sigmas_sim
            tab_var_sigma(s_i,2)=sigmat_sim
            tab_var_sigma(s_i,3)=rho_sim
            tab_var_sigma(s_i,4)=gamma_sim
            tab_var_sigma(s_i,5)=theta_simt
            tab_var_sigma(s_i,6)=rho_sim_wij
            tab_var_sigma(s_i,7)=gamma_simt
            tab_var_sigma(s_i,8)=rho_sim_ui
            moy_trt=moy_trt+(sum(don_simul(:,trt1))*100.d0/n_obs)! proportion des traites
            moy_theta=moy_theta+(theta_sim)  ! on calcule la moyenne empirique des theta S
            moy_thetat=moy_thetat+(theta_simt)  ! on calcule la moyenne empirique des theta T
            moy_rho_wij=moy_rho_wij+(rho_sim_wij)  ! on calcule la moyenne empirique des rho wij
            moy_sigmas=moy_sigmas+(sigmas_sim)  ! on calcule la moyenne empirique des sigmas
            moy_sigmat=moy_sigmat+(sigmat_sim)  ! on calcule la moyenne empirique des sigmat
            moy_rho=moy_rho+(rho_sim)  ! on calcule la moyenne empirique des rho
            moy_gamma=moy_gamma+gamma_sim ! on calcule la moyenne empirique des gamma_ui
            moy_gammat=moy_gammat+gamma_simt ! on calcule la moyenne empirique des gamma_ui
            moy_rho_ui=moy_rho_ui+(rho_sim_ui)  ! on calcule la moyenne empirique des rho wij
            moy_pros=moy_pros+(sum(don_simulS(:,statusS1))*100.d0/n_obs) !progression
            moy_dec=moy_dec+(sum(don_simul(:,statusT1))*100.d0/n_obs) ! deces
        endif
                  
    else
        !on conserve le kappa
        if(((kapa_use.eq.2).or.(kapa_use.eq.3).or.(kapa_use.eq.4)).and.(statut_kappa==0)) then ! on conserve ce kappa juste s'il est nouveau
        !on recupere tous les kappa qui marchent
            kappa(incre_kappa,1)=k0(1)
            kappa(incre_kappa,2)=k0(2)
            incre_kappa=incre_kappa+1
            ! if(rang_proc==0) !print*,"j'incremente le kappa s_i=",s_i,"incre_kappa=",incre_kappa
        endif
        
        !!print*,"suis la main2 b=",b
        nva=nva1+nva2
        rangparam=np-nva-nparamfrail+indice_eta+indice_theta 
        rangparam_theta=rangparam
        if(indice_eta==1) rangparam_eta=np-nva-nparamfrail+indice_eta
        rangparam_sigs=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS
        rangparam_sigt=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT
        rangparam_sigst=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST
        rangparam_alpha=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha
        rangparam_gamma=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+indice_gamma
        rangparam_copula=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+&
        indice_alpha+indice_gamma + 1
        
        ! call intpr("indice_eta", -1, indice_eta, 1)
        ! call intpr("indice_theta", -1, indice_theta, 1)
        ! call intpr("rangparam_sigs", -1, rangparam_sigs, 1)
        ! call intpr("rangparam_sigt", -1, rangparam_sigt, 1)
        ! call dblepr("H_hessOut", -1, H_hessOut, size(H_hessOut,2)**2)
        
        if(nsim_node(8)==2)then !effets aleatoires correles
            rangparam_thetat=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                indice_gamma+indice_theta_t
            rangparam_thetast=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                indice_gamma+indice_theta_t+indice_theta_st
            rangparam_gammat=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                indice_gamma+indice_theta_t+indice_theta_st+indice_gamma_t
            rangparam_gammast=np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha+&
                indice_gamma+indice_theta_t+indice_theta_st+indice_gamma_t+indice_gamma_st
        endif
        
        rangparam2=rangparam
        ! elements de calcul des estimations
        theta_sim=variance(don_simul(:,w_ij1))! variance des effets aleatoires simules
        tab_var_theta(s_i-nbre_rejet)=theta_sim
        sigmas_sim=variance(donnee_essai(:,1))! sigma_s empirique
        sigmat_sim=variance(donnee_essai(:,2))! sigma_t empirique
        theta_simt=variance(don_simul(:,w_ijt))! variance des effets aleatoires simules T
        rho_sim_wij=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/(dsqrt(variance(don_simul(:,w_ij1)))&
                *dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
        gamma_simt=variance(donnee_essai(:,5))! gamma T empirique
        rho_sim_ui=covariance(donnee_essai(:,4),donnee_essai(:,5))/(dsqrt(variance(donnee_essai(:,4)))*&
                dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST
            
        rho_sim=covariance(donnee_essai(:,1),donnee_essai(:,2))/&
            (dsqrt(variance(donnee_essai(:,1)))*dsqrt(variance(donnee_essai(:,2))))
        gamma_sim=variance(donnee_essai(:,4)) !gamma_ui empirique
        tab_var_sigma(s_i-nbre_rejet,1)=sigmas_sim
        tab_var_sigma(s_i-nbre_rejet,2)=sigmat_sim
        tab_var_sigma(s_i-nbre_rejet,3)=rho_sim
        tab_var_sigma(s_i-nbre_rejet,4)=gamma_sim
        tab_var_sigma(s_i-nbre_rejet,5)=theta_simt
        tab_var_sigma(s_i-nbre_rejet,6)=rho_sim_wij
        tab_var_sigma(s_i-nbre_rejet,7)=gamma_simt
        tab_var_sigma(s_i-nbre_rejet,8)=rho_sim_ui
        moy_trt=moy_trt+(sum(don_simul(:,trt1))*100.d0/n_obs)! proportion des traites
        moy_theta=moy_theta+(theta_sim)  ! on calcule la moyenne empirique des theta
        moy_sigmas=moy_sigmas+(sigmas_sim)  ! on calcule la moyenne empirique des sigmas
        moy_sigmat=moy_sigmat+(sigmat_sim)  ! on calcule la moyenne empirique des sigmat
        moy_rho=moy_rho+(rho_sim)  ! on calcule la moyenne empirique des rho
        moy_gamma=moy_gamma+gamma_sim ! on calcule la moyenne empirique des gamma_ui
        moy_thetat=moy_thetat+(theta_simt)  ! on calcule la moyenne empirique des theta T
        moy_rho_wij=moy_rho_wij+(rho_sim_wij)  ! on calcule la moyenne empirique des rho wij
        moy_gammat=moy_gammat+gamma_simt ! on calcule la moyenne empirique des gamma_ui
        moy_rho_ui=moy_rho_ui+(rho_sim_ui)  ! on calcule la moyenne empirique des rho wij
        !!print*,statusS1,sum(don_simulS(:,statusS1))
        !stop
        moy_pros=moy_pros+(sum(don_simulS(:,statusS1))*100.d0/n_obs) !progression
        moy_dec=moy_dec+(sum(don_simul(:,statusT1))*100.d0/n_obs) ! deces
        
        ! dans ce qui suit, les se de sigma_s, sigma_t et sigma_st sont calcules par la delta methode, suite au changement de variable
        ! theta estime
        ! if(rang_proc==0) !print*,"moy se des theta:",rangparam,H_hessOut(rangparam,rangparam),&
                ! (dsqrt(((2.d0*b(rangparam))**2.d0)*H_hessOut(rangparam,rangparam))),&
                ! b(rangparam),&
                ! moy_se_theta+(dsqrt(((2.d0*b(rangparam))**2.d0)*H_hessOut(rangparam,rangparam)))
                
        moy_theta_est=moy_theta_est+(b(rangparam)**2.d0) ! theta estime
        moy_se_theta=moy_se_theta+(dsqrt(((2.d0*b(rangparam))**2.d0)*H_hessOut(rangparam,rangparam)))
        bi_theta = b(rangparam)**2.d0 - 1.96d0*(dsqrt(((2.d0*b(rangparam))**2.d0)*&
                    H_hessOut(rangparam,rangparam)))
        bs_theta = b(rangparam)**2.d0 + 1.96d0*(dsqrt(((2.d0*b(rangparam))**2.d0)*&
                    H_hessOut(rangparam,rangparam)))
        !taux de couverture
        if(theta>=bi_theta .and. theta<=bs_theta)then ! taux de couverture
            taux_couverture_theta=taux_couverture_theta+1.d0
        endif
        
        if(nsim_node(8)==3) then
            if(copula_function == 1) then  ! claton: exp transform
                moy_param_cop = moy_param_cop+dexp(b(rangparam_copula)) ! param copule estime
                ! par delta methode, SE = sqrt(exp(theta_chap) * H_theta_chap * exp(theta_chap))
                pour_ic = dexp(b(rangparam_copula)) * dsqrt(H_hessOut(rangparam_copula,rangparam_copula))
                bi_param_cop = dexp(b(rangparam_copula)) - 1.96d0 * pour_ic
                bs_param_cop = dexp(b(rangparam_copula)) + 1.96d0 * pour_ic
            endif
            if(copula_function == 2) then  ! Gumbel: choleschy transform
                moy_param_cop = moy_param_cop+(b(rangparam_copula)**2.d0) ! param copule estime
                pour_ic = 2.d0*b(rangparam_copula) * dsqrt(H_hessOut(rangparam_copula,rangparam_copula))
                bi_param_cop = b(rangparam_copula)**2.d0 - 1.96d0 * pour_ic
                bs_param_cop = b(rangparam_copula)**2.d0 + 1.96d0 * pour_ic
            endif
            moy_se_param_cop=moy_se_param_cop + pour_ic
            !taux de couverture
            if(thetacopule>=bi_param_cop .and. thetacopule<=bs_param_cop)then ! taux de couverture
                taux_couverture_param_cop = taux_couverture_param_cop+1.d0
            endif
        
        endif
        if(nsim_node(8)==2) then
            !theta estime. j'utilise les variables au niveau essai juste pour ne pas avoir a declarer de nouvelles et pas pour grand chose
            if(logNormal==1)then
                thetaS1 = b(rangparam)
                thetaT1 = b(rangparam_thetat)
                thetaST1 =b(rangparam_thetast)
                ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
                !Chol: matrice triangulaire inferieur
                Chol=0.d0 
                Chol(1,1)=thetaS1
                Chol(2,1)=thetaST1
                !Chol(1,2)=covST1
                Chol(2,2)=thetaT1
                theta_st=MATMUL(Chol,TRANSPOSE(Chol))
                thetaS_es=theta_st(1,1)
                thetaT_es=theta_st(2,2)
                thetaST_es=theta_st(1,2)
                
                !theta_t
                moy_thetat_est=moy_thetat_est+thetaT_es ! sigma_t estime
                moy_se_thetat=moy_se_thetat+2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                            2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                            thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))
                bi_sigmat = thetaT_es - 1.96d0*2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                            2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                            thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))
                bs_sigmat = thetaT_es + 1.96d0*2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                            2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                            thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))
                !taux de couverture
                
                if(theta2_t>=bi_sigmat .and. theta2_t<=bs_sigmat)then ! taux de couverture
                    taux_couverture_thetat=taux_couverture_thetat+1.d0
                endif
            
                !sigma_st
                moy_thetast_est=moy_thetast_est+thetaST_es ! sigma_t estime
                moy_se_thetast=moy_se_thetast+dsqrt(thetaST1**2.d0*H_hessOut(rangparam,rangparam)+&
                            2.d0*thetaS1*thetaST1*H_hessOut(rangparam,rangparam_thetast)+&
                            thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))
                bi_sigmast = thetaST_es - 1.96d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam,rangparam)+&
                            2.d0*thetaS1*thetaST1*H_hessOut(rangparam,rangparam_thetast)+&
                            thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))
                bs_sigmast = thetaST_es + 1.96d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam,rangparam)+&
                            2.d0*thetaS1*thetaST1*H_hessOut(rangparam,rangparam_thetast)+&
                            thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))
                !taux de couverture
                !sigmast_vrai=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                if(thetast_vrai>=bi_sigmast .and. thetast_vrai<=bs_sigmast)then ! taux de couverture
                    taux_couverture_thetast=taux_couverture_thetast+1.d0
                endif
            endif    
        endif
        
        moy_ni=moy_ni+ni ! Nombre moyen d'itteration pour la convergence
                        
        if(nsim_node(8).ne.0)then !model conjoint surrogate
            ! gamma estime
            if(frailt_base==1) then        
                ! gamma estime S             
                moy_gamma_est=moy_gamma_est+(b(rangparam_gamma)**2.d0) ! theta estime
                moy_se_gamma=moy_se_gamma+(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*H_hessOut(rangparam_gamma,&
                            rangparam_gamma)))
                bi_gamma = b(rangparam_gamma)**2.d0 - 1.96d0*(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*H_hessOut(rangparam_gamma,&
                            rangparam_gamma)))
                bs_gamma = b(rangparam_gamma)**2.d0 + 1.96d0*(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*H_hessOut(rangparam_gamma,&
                            rangparam_gamma)))
                !taux de couverture
                if(gamma_ui>=bi_gamma .and. gamma_ui<=bs_gamma)then ! taux de couverture
                    taux_couverture_gamma=taux_couverture_gamma+1.d0
                endif
                
                if(nsim_node(8)==2)then
                    gammaS1 = b(rangparam_gamma)
                    gammaT1 = b(rangparam_gammat)
                    gammaST1 =b(rangparam_gammast)
                    ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
                    !Chol: matrice triangulaire inferieur
                    Chol=0.d0 
                    Chol(1,1)=gammaS1
                    Chol(2,1)=gammaST1
                    !Chol(1,2)=gammaST1
                    Chol(2,2)=gammaT1
                    gamma_st=MATMUL(Chol,TRANSPOSE(Chol))
                    gammaS_es=gamma_st(1,1)
                    gammaT_es=gamma_st(2,2)
                    gammaST_es=gamma_st(1,2)
                    
                    ! gamma estime T
                    moy_gammat_est=moy_gammat_est+gammaT_es ! sigma_t estime
                    moy_se_gammat=moy_se_gammat+2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                                2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                                gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))
                    bi_sigmat = gammaT_es - 1.96d0*2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                                2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                                gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))
                    bs_sigmat = gammaT_es + 1.96d0*2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                                2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                                gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))
                    !taux de couverture
                    
                    if(gamma_uit>=bi_sigmat .and. gamma_uit<=bs_sigmat)then ! taux de couverture
                        taux_couverture_gammat=taux_couverture_gammat+1.d0
                    endif
                
                    !gamma_st
                    moy_gammast_est=moy_gammast_est+gammaST_es ! sigma_t estime
                    moy_se_gammast=moy_se_gammast+dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                                2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                                gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))
                    bi_sigmast = gammaST_es - 1.96d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                                2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                                gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))
                    bs_sigmast = gammaST_es + 1.96d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                                2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                                gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))
                    !taux de couverture
                    !sigmast_vrai=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
                    if(gammast_vrai>=bi_sigmast .and. gammast_vrai<=bs_sigmast)then ! taux de couverture
                        taux_couverture_gammast=taux_couverture_gammast+1.d0
                    endif
                endif
            endif
            
            !sigma estimes
            covST_es = 0.d0
            if(logNormal==1)then
                varS1 = b(rangparam_sigs)
                varT1 = b(rangparam_sigt)
                covST1 =b(rangparam_sigst)
                ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
                !Chol: matrice triangulaire inferieur
                Chol=0.d0 
                Chol(1,1)=varS1
                Chol(2,1)=covST1
                !Chol(1,2)=covST1
                Chol(2,2)=varT1
                sigma_st=MATMUL(Chol,TRANSPOSE(Chol))
                varS_es=sigma_st(1,1)
                varT_es=sigma_st(2,2)
                covST_es=sigma_st(1,2)
                
                ! =========Delta methode pour varcov des elements de sigme_v===============
                ! recherche de la matrice de variance-covariance de (sigma_S,sigma_ST,sigmaT) par la delta methode:
                ! à partir de la hessienne. voir le raisonnement dans le cahier à la date du 04/01/2019
                hb(1,:) = (/ 2.d0*Chol(1,1), 0.d0, 0.d0 /)
                hb(2,:) = (/ 0.d0, 2.d0*Chol(2,2), 2.d0*Chol(2,1) /)
                hb(3,:) = (/ Chol(2,1), 0.d0, Chol(1,1) /)
                sigmac(1,:) = (/H_hessOut(rangparam_sigs,rangparam_sigs), H_hessOut(rangparam_sigs,rangparam_sigt), &
                                H_hessOut(rangparam_sigs,rangparam_sigst)/)
                sigmac(2,:) = (/H_hessOut(rangparam_sigt,rangparam_sigs), H_hessOut(rangparam_sigt,rangparam_sigt), &
                                H_hessOut(rangparam_sigt,rangparam_sigst)/)
                sigmac(3,:) = (/H_hessOut(rangparam_sigst,rangparam_sigs), H_hessOut(rangparam_sigst,rangparam_sigt), &
                                H_hessOut(rangparam_sigst,rangparam_sigst)/)
                
                hb = TRANSPOSE(hb)
                varcov = MATMUL(TRANSPOSE(hb), sigmac)
                varcov = MATMUL(varcov, hb)
                
                
                ! ========== Fin delta methode ==================
                
                ! ====sauvegarde de la hessienne et du vecteur b des parametres====
                
                do i=1,np
                    dataHessian(np*(s_i-nbre_rejet-1) + i,:) = H_hessOut(i,:)
                    dataHessianIH(np*(s_i-nbre_rejet-1) + i,:) = HIHOut(i,:)
                enddo
                
                datab(s_i-nbre_rejet,:) = b
                
                ! Write(aaa,'(i3)') s_i ! instruction pour convertir un entier en chaine de caractere (3 caracteres)
                ! aaa="H_hessOut"//aaa ! instruction pour concatener deux chaines de caracteres
                ! aaa=aaa//".txt"
                ! open(270,file=aaa)                
                ! write(270,*) dataHessian                
                ! close(270)
                
                !====Fin sauvegarde hessienne et b======
                !calcul du R2(trial) reduit, c'est a dire sans prise en compte des effets aleatoires sur le risque de base
                R2_trial=(covST1**2)/(covST1**2+varT1**2)
                se_R2_trial=2.d0*dsqrt((covST1**2 * varT1**4 * H_hessOut(rangparam_sigst,rangparam_sigst)-2.d0*covST1**3 &
                * varT1**3 *H_hessOut(rangparam_sigst,rangparam_sigt) + varT1**2 * covST1**4 * H_hessOut(rangparam_sigt,&
                rangparam_sigt)))/(covST1**2 + varT1**2)**2
                moy_R2_trial=moy_R2_trial+R2_trial
                moy_se_R2_trial=moy_se_R2_trial+se_R2_trial
                bi_R2_trial = R2_trial - 1.96d0*se_R2_trial
                bs_R2_trial = R2_trial + 1.96d0*se_R2_trial
                moy_bi_R2_trial=moy_bi_R2_trial+bi_R2_trial
                moy_bs_R2_trial=moy_bs_R2_trial+bs_R2_trial
                !taux de couverture
                if((rsqrt**2)>=bi_R2_trial .and. (rsqrt**2)<=bs_R2_trial)then ! taux de couverture
                    taux_couverture_R2_trial=taux_couverture_R2_trial+1.d0
                endif
            endif
            
            ! dans ce qui suit, les se de sigma_s, sigma_t et sigma_st sont calcules par la delta methode, suite au changement de variable
            
            !sigma_s
            moy_sigmas_est=moy_sigmas_est+varS_es ! sigma_s estime
            ! moy_se_sigmas=moy_se_sigmas+(dsqrt(((2.d0*b(rangparam_sigs))**2.d0)*H_hessOut(rangparam_sigs,rangparam_sigs)))
            ! bi_sigmas = varS_es - 1.96d0*(dsqrt(((2.d0*b(rangparam_sigs))**2.d0)*H_hessOut(rangparam_sigs,rangparam_sigs)))
            ! bs_sigmas = varS_es + 1.96d0*(dsqrt(((2.d0*b(rangparam_sigs))**2.d0)*H_hessOut(rangparam_sigs,rangparam_sigs)))
            moy_se_sigmas=moy_se_sigmas+ dsqrt(varcov(1,1))
            bi_sigmas = varS_es - 1.96d0*dsqrt(varcov(1,1))
            bs_sigmas = varS_es + 1.96d0*dsqrt(varcov(1,1))
            !taux de couverture
            if(sigma_s>=bi_sigmas .and. sigma_s<=bs_sigmas)then ! taux de couverture
                taux_couverture_sigmas=taux_couverture_sigmas+1.d0
            endif
        
            !sigma_t
            moy_sigmat_est=moy_sigmat_est+varT_es ! sigma_t estime
            ! moy_se_sigmat=moy_se_sigmat+2.d0*dsqrt(covST1**2.d0*H_hessOut(rangparam_sigst,rangparam_sigst)+&
                        ! 2.d0*varT1*covST1*H_hessOut(rangparam_sigst,rangparam_sigt)+&
                        ! varT1**2.d0*H_hessOut(rangparam_sigt,rangparam_sigt))
            ! bi_sigmat = varT_es - 1.96d0*2.d0*dsqrt(covST1**2.d0*H_hessOut(rangparam_sigst,rangparam_sigst)+&
                        ! 2.d0*varT1*covST1*H_hessOut(rangparam_sigst,rangparam_sigt)+&
                        ! varT1**2.d0*H_hessOut(rangparam_sigt,rangparam_sigt))
            ! bs_sigmat = varT_es + 1.96d0*2.d0*dsqrt(covST1**2.d0*H_hessOut(rangparam_sigst,rangparam_sigst)+&
                        ! 2.d0*varT1*covST1*H_hessOut(rangparam_sigst,rangparam_sigt)+&
                        ! varT1**2.d0*H_hessOut(rangparam_sigt,rangparam_sigt))
            moy_se_sigmat=moy_se_sigmat+dsqrt(varcov(2,2))
            bi_sigmat = varT_es - 1.96d0*dsqrt(varcov(2,2))
            bs_sigmat = varT_es + 1.96d0*dsqrt(varcov(2,2))
            !taux de couverture
            
            if(sigma_t>=bi_sigmat .and. sigma_t<=bs_sigmat)then ! taux de couverture
                taux_couverture_sigmat=taux_couverture_sigmat+1.d0
            endif
        
            !sigma_st
            moy_sigmast_est=moy_sigmast_est+covST_es ! sigma_t estime
            ! moy_se_sigmast=moy_se_sigmast+dsqrt(covST1**2.d0*H_hessOut(rangparam_sigs,rangparam_sigs)+&
                        ! 2.d0*varS1*covST1*H_hessOut(rangparam_sigs,rangparam_sigst)+&
                        ! varS1**2.d0*H_hessOut(rangparam_sigst,rangparam_sigst))
            ! bi_sigmast = covST_es - 1.96d0*dsqrt(covST1**2.d0*H_hessOut(rangparam_sigs,rangparam_sigs)+&
                        ! 2.d0*varS1*covST1*H_hessOut(rangparam_sigs,rangparam_sigst)+&
                        ! varS1**2.d0*H_hessOut(rangparam_sigst,rangparam_sigst))
            ! bs_sigmast = covST_es + 1.96d0*dsqrt(covST1**2.d0*H_hessOut(rangparam_sigs,rangparam_sigs)+&
                        ! 2.d0*varS1*covST1*H_hessOut(rangparam_sigs,rangparam_sigst)+&
                        ! varS1**2.d0*H_hessOut(rangparam_sigst,rangparam_sigst))
            moy_se_sigmast=moy_se_sigmast+dsqrt(varcov(3,3))
            bi_sigmast = covST_es - 1.96d0*dsqrt(varcov(3,3))
            bs_sigmast = covST_es + 1.96d0*dsqrt(varcov(3,3))
            !taux de couverture
            !sigmast_vrai=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
            if(sigmast_vrai>=bi_sigmast .and. sigmast_vrai<=bs_sigmast)then ! taux de couverture
                taux_couverture_sigmast=taux_couverture_sigmast+1.d0
            endif
        endif
        
        !eta
        if(indice_eta==1)then
            moy_eta=moy_eta+(b(np-nva-nparamfrail+indice_eta))
            moy_se_eta=moy_se_eta+(dsqrt(H_hessOut(np-nva-nparamfrail+indice_eta,np-nva-nparamfrail+indice_eta)))
            bi_eta = b(np-nva-nparamfrail+indice_eta) - 1.96d0*(dsqrt(H_hessOut(np-nva-nparamfrail+indice_eta,np-&
                    nva-nparamfrail+indice_eta)))
            bs_eta = b(np-nva-nparamfrail+indice_eta) + 1.96d0*(dsqrt(H_hessOut(np-nva-nparamfrail+indice_eta,np-&
                    nva-nparamfrail+indice_eta)))
            !taux de couverture
            if(eta>=bi_eta .and. eta<=bs_eta)then ! taux de couverture
                taux_couverture_eta=taux_couverture_eta+1.d0
            endif
        endif
        
        !alpha
        if(indice_alpha==1)then            
            moy_alpha=moy_alpha+(b(rangparam_alpha))
            moy_se_alpha=moy_se_alpha+(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha)))
            bi_eta = b(rangparam_alpha) - 1.96d0*(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha)))
            bs_eta = b(rangparam_alpha) + 1.96d0*(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha)))
            !taux de couverture
            if(alpha_ui>=bi_eta .and. alpha_ui<=bs_eta)then ! taux de couverture
                taux_couverture_alpha=taux_couverture_alpha+1.d0
            endif
        endif
        
        if(nsim_node(8).ne.3) then 
            do i=1,nva
                rangparam=np-nva+i
                !Intervalle de confiance
                bi2 = b(rangparam) - 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                bs2 = b(rangparam) + 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                    
                if(i.eq.1)then
                    ! if(rang_proc==0) !write(*,*)'*** For surrogate ***'
                    moy_betaS(1)=moy_betaS(1)+b(rangparam)
                    moy_betaS_se(1)=moy_betaS_se(1)+(dsqrt(H_hessOut(rangparam,rangparam)))
                    parametre_estimes(s_i-nbre_rejet,5)=b(rangparam) !beta_s_chap
                    parametre_estimes(s_i-nbre_rejet,6)=dsqrt(H_hessOut(rangparam,rangparam)) !sd beta_s_chap
                    if(betas>=bi2 .and. betas<=bs2)then ! taux de couverture
                        taux_couvertureS(1)=taux_couvertureS(1)+1.d0
                    endif
                endif
                    
                if(i.eq.nva1+1)then
                    ! if(rang_proc==0) !write(*,*)'*** For true endpoint ***',nva,i
                    moy_betaT(1)=moy_betaT(1)+b(rangparam)
                    moy_betaT_se(1)=moy_betaT_se(1)+(dsqrt(H_hessOut(rangparam,rangparam)))
                    parametre_estimes(s_i-nbre_rejet,7)=b(rangparam) !beta_t_chap
                    parametre_estimes(s_i-nbre_rejet,8)=(dsqrt(H_hessOut(rangparam,rangparam))) !sd beta_t_chap
                    if(betat>=bi2 .and. betat<=bs2)then ! taux de couverture
                        taux_couvertureT(1)=taux_couvertureT(1)+1.d0
                    endif
                endif
            end do 
        else ! copula model
            do i=1,nva1
                rangparam=np-nva+i
                !Intervalle de confiance
                bi2 = b(rangparam) - 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                bs2 = b(rangparam) + 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                moy_betaS(i)=moy_betaS(i)+b(rangparam)
                moy_betaS_se(i)=moy_betaS_se(i)+(dsqrt(H_hessOut(rangparam,rangparam)))
                if(i.eq.1) then ! traitement
                    parametre_estimes(s_i-nbre_rejet,5)=b(rangparam) !beta_s_chap
                    parametre_estimes(s_i-nbre_rejet,6)=dsqrt(H_hessOut(rangparam,rangparam)) !sd beta_s_chap
                else ! sauvegarde des autres covariables
                    parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*(nva2-1)-2*nva1+2*(i-1)+1)=b(rangparam) 
                    parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*(nva2-1)-2*nva1+2*(i-1)+2)=&
                    dsqrt(H_hessOut(rangparam,rangparam)) 
                endif
                if(betas>=bi2 .and. betas<=bs2)then ! taux de couverture
                    taux_couvertureS(i)=taux_couvertureS(i)+1.d0
                endif
            enddo
            do i=1,nva2
                rangparam = np - nva + nva1 + i
                !Intervalle de confiance
                bi2 = b(rangparam) - 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                bs2 = b(rangparam) + 1.96d0*dsqrt(H_hessOut(rangparam,rangparam))
                    
                !if(i.eq.nva1+1)then
                ! if(rang_proc==0) !write(*,*)'*** For true endpoint ***',nva,i
                moy_betaT(i)=moy_betaT(i)+b(rangparam)
                moy_betaT_se(i)=moy_betaT_se(i)+(dsqrt(H_hessOut(rangparam,rangparam)))
                if(i.eq.1) then ! traitement
                    parametre_estimes(s_i-nbre_rejet,7)=b(rangparam) !beta_t_chap
                    parametre_estimes(s_i-nbre_rejet,8)=dsqrt(H_hessOut(rangparam,rangparam)) !sd beta_t_chap
                else ! sauvegarde des autres covariables
                    parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*nva2+2*(i-1)+1)=b(rangparam) 
                    parametre_estimes(s_i-nbre_rejet,size(parametre_estimes,2)-2*nva2+2*(i-1)+2)=&
                    dsqrt(H_hessOut(rangparam,rangparam)) 
                endif    
                if(betat>=bi2 .and. betat<=bs2)then ! taux de couverture
                    taux_couvertureT(i)=taux_couvertureT(i)+1.d0
                endif
                !endif
            end do 
        endif
        
        ! parametres empiriques
        parametre_empirique(s_i-nbre_rejet,1)=sum(don_simul(:,trt1))*100.d0/n_obs ! personne traitees
        parametre_empirique(s_i-nbre_rejet,2)=variance(don_simul(:,w_ij1)) ! theta empirique
        parametre_empirique(s_i-nbre_rejet,3)=sum(don_simul(:,statusS1))*100.d0/n_obs ! progression
        parametre_empirique(s_i-nbre_rejet,4)=sum(don_simul(:,statusT1))*100.d0/n_obs
            
        if(nsim_node(8).ne.0)then !model conjoint surrogate
            parametre_empirique(s_i-nbre_rejet,5)=variance(donnee_essai(:,1))! sigma_s empirique
            parametre_empirique(s_i-nbre_rejet,6)=variance(donnee_essai(:,2))! sigma_t empirique
            parametre_empirique(s_i-nbre_rejet,7)=covariance(donnee_essai(:,1),donnee_essai(:,2))/&
                (dsqrt(variance(donnee_essai(:,1)))*dsqrt(variance(donnee_essai(:,2))))! correlation empirique: rho(x,y)=cov(x,y)/(sigma_x*sigma_y)
            parametre_empirique(s_i-nbre_rejet,8)=variance(donnee_essai(:,4))! gamma empirique
            if(nsim_node(8)==2)then
                parametre_empirique(s_i-nbre_rejet,9)=variance(don_simul(:,w_ijt))! theta t empirique
                parametre_empirique(s_i-nbre_rejet,10)=covariance(don_simul(:,w_ij1),don_simul(:,w_ijt))/&
                    (dsqrt(variance(don_simul(:,w_ij1)))*dsqrt(variance(don_simul(:,w_ijt))))! correlation empirique wij_ST
                parametre_empirique(s_i-nbre_rejet,11)=variance(donnee_essai(:,5))! gamma T empirique
                parametre_empirique(s_i-nbre_rejet,12)=covariance(donnee_essai(:,4),donnee_essai(:,5))/&
                    (dsqrt(variance(donnee_essai(:,4)))*dsqrt(variance(donnee_essai(:,5))))! correlation empirique ui_ST        
            endif
        endif
            !!write(7,*)parametre_empirique(s_i-nbre_rejet,:)
            
            !!print*,"suis là 1 s_i=",s_i
        !Parametres estimes
        if(nsim_node(8).ne.3) then
            parametre_estimes(s_i-nbre_rejet,1)=(b(rangparam2)**2.d0) !theta
            parametre_estimes(s_i-nbre_rejet,2)=(dsqrt(((2.d0*b(rangparam2))**2.d0)*H_hessOut(rangparam2,rangparam2)))! se theta
        else ! copula
            if(copula_function == 1) then ! claton: exp transform
                parametre_estimes(s_i-nbre_rejet,1) = dexp(b(rangparam_copula)) !theta_copula
                parametre_estimes(s_i-nbre_rejet,2) = pour_ic
            endif
            if(copula_function == 2) then ! Gumbel
                parametre_estimes(s_i-nbre_rejet,1) = b(rangparam_copula)**2.d0 !theta_copula
                parametre_estimes(s_i-nbre_rejet,2) = pour_ic
            endif
            
        endif
        
        if(indice_eta==1)then
            parametre_estimes(s_i-nbre_rejet,3)=(b(np-nva-nparamfrail+indice_eta))!zeta
            parametre_estimes(s_i-nbre_rejet,4)=(dsqrt(H_hessOut(np-nva-nparamfrail+indice_eta,np-nva-nparamfrail+indice_eta))) !se zeta
        else
            parametre_estimes(s_i-nbre_rejet,3)=1.d0!zeta
            parametre_estimes(s_i-nbre_rejet,4)=0.d0 !se zeta
        endif
            
        if(nsim_node(8).ne.0)then !model conjoint surrogate (se calcule par la delta methode)    
            parametre_estimes(s_i-nbre_rejet,9)= varS_es!siqma_s
            parametre_estimes(s_i-nbre_rejet,10)=dsqrt(varcov(1,1))! se sigma_s
            parametre_estimes(s_i-nbre_rejet,11)= varT_es!siqma_t
            parametre_estimes(s_i-nbre_rejet,12)=dsqrt(varcov(2,2))! se sigma_t
            parametre_estimes(s_i-nbre_rejet,13)= covST_es !siqma_st
            parametre_estimes(s_i-nbre_rejet,14)=dsqrt(varcov(3,3))! se sigma_st
            if(frailt_base==1)then
                parametre_estimes(s_i-nbre_rejet,15)=(b(rangparam_gamma)**2.d0) !gamma
                parametre_estimes(s_i-nbre_rejet,16)=(dsqrt(((2.d0*b(rangparam_gamma))**2.d0)*H_hessOut(rangparam_gamma,&
                    rangparam_gamma)))! se gamma
                if(indice_alpha==1)then
                    parametre_estimes(s_i-nbre_rejet,17)=(b(rangparam_alpha))!alpha
                    parametre_estimes(s_i-nbre_rejet,18)=(dsqrt(H_hessOut(rangparam_alpha,rangparam_alpha))) !se alpha
                else
                    parametre_estimes(s_i-nbre_rejet,17)=1.d0!zeta
                    parametre_estimes(s_i-nbre_rejet,18)=0.d0 !se zeta
                endif
                
            endif
            
            !!print*,"suis là 2 s_i=",s_i
            ! pour les R2 trial 
            if(nsim_node(8)==1 .or. nsim_node(8)==3)then
                parametre_estimes(s_i-nbre_rejet,19)=R2_trial     !r2_trial reduite
                parametre_estimes(s_i-nbre_rejet,20)=se_R2_trial    !se r2_trial reduite
            endif
            
            if(nsim_node(8)==2)then
                parametre_estimes(s_i-nbre_rejet,19)= thetaT_es!theta_t
                parametre_estimes(s_i-nbre_rejet,20)=2.d0*dsqrt(thetaST1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast)+&
                        2.d0*thetaT1*thetaST1*H_hessOut(rangparam_thetast,rangparam_thetat)+&
                        thetaT1**2.d0*H_hessOut(rangparam_thetat,rangparam_thetat))! se theta_t
                parametre_estimes(s_i-nbre_rejet,21)= thetaST_es !theta_st
                parametre_estimes(s_i-nbre_rejet,22)=dsqrt(thetaST1**2.d0*H_hessOut(rangparam_theta,rangparam_theta)+&
                        2.d0*thetaS1*thetaST1*H_hessOut(rangparam_theta,rangparam_thetast)+&
                        thetaS1**2.d0*H_hessOut(rangparam_thetast,rangparam_thetast))! se theta_st
                !gamma
                parametre_estimes(s_i-nbre_rejet,23)= gammaT_es!gamma_t
                parametre_estimes(s_i-nbre_rejet,24)=2.d0*dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast)+&
                        2.d0*gammaT1*gammaST1*H_hessOut(rangparam_gammast,rangparam_gammat)+&
                        gammaT1**2.d0*H_hessOut(rangparam_gammat,rangparam_gammat))! se gamma_t
                parametre_estimes(s_i-nbre_rejet,25)= gammaST_es !gamma_st
                parametre_estimes(s_i-nbre_rejet,26)=dsqrt(gammaST1**2.d0*H_hessOut(rangparam_gamma,rangparam_gamma)+&
                        2.d0*gammaS1*gammaST1*H_hessOut(rangparam_gamma,rangparam_gammast)+&
                        gammaS1**2.d0*H_hessOut(rangparam_gammast,rangparam_gammast))! se gamma_st
                
                ! pour les R2 trial reduits et complets
                parametre_estimes(s_i-nbre_rejet,27)=R2_trial     !r2_trial reduite
                parametre_estimes(s_i-nbre_rejet,28)=se_R2_trial    !se r2_trial reduite
            endif
            
            ! !print*,"suis là 3 theta_ST=",theta_ST
            ! !print*,"suis là 3 gamma_st=",gamma_st
            ! !print*,"suis là 3 sigma_st=",sigma_st
            
            ! ======calcul du taux de kendall =================
            call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour le calcul du tau de kendall
            if(nsim_node(8)==2)then
                if(method_int_kendal==4 .or.method_int_kendal==5)then ! 1 seul taux de kendall et modele complet
                    tau_kendal_00=tau_kendall(theta_ST,gamma_st,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                    parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des non traites z_11=0,z_21=0
                    moy_kendal_00=moy_kendal_00+tau_kendal_00
                    
                else    ! 4 taux de kendall (suivant les bras de traitement) et modele complet
                    tau_kendal_11=tau_kendall(theta_ST,gamma_st,sigma_st,1,1,method_int_kendal,N_MC_kendall,&
                        parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des traites z_11=1,z_21=1
                    tau_kendal_10=tau_kendall(theta_ST,gamma_st,sigma_st,1,0,method_int_kendal,N_MC_kendall,&
                        parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                    tau_kendal_01=tau_kendall(theta_ST,gamma_st,sigma_st,0,1,method_int_kendal,N_MC_kendall,&
                        parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des traite z_11=0,z_21=1
                    tau_kendal_00=tau_kendall(theta_ST,gamma_st,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                        parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),1)    !tau de kendal des non traites z_11=0,z_21=0
                    moy_kendal_11=moy_kendal_11+tau_kendal_11
                    moy_kendal_10=moy_kendal_10+tau_kendal_10
                    moy_kendal_01=moy_kendal_01+tau_kendal_01
                    moy_kendal_00=moy_kendal_00+tau_kendal_00
                    
                endif
            endif

            if(nsim_node(8)==1 .or. nsim_node(8)==3)then ! 1 seul taux de kendall et modele reduit
                allocate(theta_ST_2(1,1),gamma_st_2(1,1))
                theta_ST_2(1,1)= parametre_estimes(s_i-nbre_rejet,1) !theta estime
                gamma_st_2(1,1)= parametre_estimes(s_i-nbre_rejet,15) ! gamma estime
                if(nsim_node(8)==1) then
                    tau_kendal_00=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                    parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),0)    !tau de kendal des non traites z_11=0,z_21=0
                    
                else !copula
                    if(copula_function == 1) then! claton
                        tau_kendal_00 = parametre_estimes(s_i-nbre_rejet,1) / (parametre_estimes(s_i-nbre_rejet,1) +2.d0) 
                    endif
                    if(copula_function == 2) then ! Gumbel
                        tau_kendal_00 = parametre_estimes(s_i-nbre_rejet,1) / (parametre_estimes(s_i-nbre_rejet,1) +1.d0)
                    endif                        
                endif
                
                moy_kendal_00=moy_kendal_00+tau_kendal_00
                
                !calcul intervalle de confiance du tau de kendall par bootstrap parametrique, a partir de la distribution a posteriorie des parametres
                if(nsim_node(8).ne.3)then ! si pas modele de copule
                    if(indice_alpha==0 .and. indice_eta==0)then ! on estime ni alpha, ni eta
                        v_chap_kendall=0.d0
                        t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma)/) ! parametres necessaire: theta, gamma, zeta, alpha
                        v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),&
                        H_hessOut(rangparam_theta,rangparam_gamma)/)
                        v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),&
                        H_hessOut(rangparam_gamma,rangparam_gamma)/)
                    else ! on estime au moins un des deux
                        if(indice_alpha==1 .and. indice_eta == 1)then !on estime les deux
                            t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma),b(rangparam_eta),b(rangparam_alpha)/) ! parametres necessaire: theta, gamma, zeta, alpha
                            v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),H_hessOut(rangparam_theta,&
                                                rangparam_gamma),H_hessOut(rangparam_theta,rangparam_eta),&
                                                H_hessOut(rangparam_theta,rangparam_alpha)/)
                            v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),H_hessOut(rangparam_gamma,&
                                                rangparam_gamma),H_hessOut(rangparam_eta,rangparam_gamma),&
                                                H_hessOut(rangparam_alpha,rangparam_gamma)/)
                            v_chap_kendall(3,:)=(/H_hessOut(rangparam_theta,rangparam_eta),H_hessOut(rangparam_eta,&
                                                rangparam_gamma),H_hessOut(rangparam_eta,rangparam_eta),&
                                                H_hessOut(rangparam_eta,rangparam_alpha)/)
                            v_chap_kendall(4,:)=(/H_hessOut(rangparam_theta,rangparam_alpha),&
                                                H_hessOut(rangparam_alpha,rangparam_gamma) ,H_hessOut(rangparam_alpha,&
                                                rangparam_eta),H_hessOut(rangparam_alpha,rangparam_alpha)/)
                        else ! on estime seulement un des deux
                            if(indice_alpha==1)then ! c'est alpha on estime
                                v_chap_kendall=0.d0
                                t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma),b(rangparam_alpha)/) ! parametres necessaire: theta, gamma, alpha
                                v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),&
                                                     H_hessOut(rangparam_theta,rangparam_gamma)&
                                                    ,H_hessOut(rangparam_theta,rangparam_alpha)/)
                                v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),&
                                                    H_hessOut(rangparam_gamma,rangparam_gamma)&
                                                    ,H_hessOut(rangparam_alpha,rangparam_gamma)/)    
                                v_chap_kendall(3,:)=(/H_hessOut(rangparam_theta,rangparam_alpha),&
                                                     H_hessOut(rangparam_gamma,rangparam_alpha)&
                                                    ,H_hessOut(rangparam_alpha,rangparam_alpha)/)
                            else ! c'est eta on estime
                                t_chap_kendall=(/b(rangparam_theta),b(rangparam_gamma),b(rangparam_eta)/) ! parametres necessaire: theta, gamma, zeta
                                v_chap_kendall(1,:)=(/H_hessOut(rangparam_theta,rangparam_theta),&
                                                    H_hessOut(rangparam_theta,rangparam_gamma)&
                                                    ,H_hessOut(rangparam_theta,rangparam_eta)/)
                                v_chap_kendall(2,:)=(/H_hessOut(rangparam_theta,rangparam_gamma),&
                                                    H_hessOut(rangparam_gamma,rangparam_gamma),&
                                                    H_hessOut(rangparam_eta,rangparam_gamma)/)
                                v_chap_kendall(3,:)=(/H_hessOut(rangparam_theta,rangparam_eta),&
                                                    H_hessOut(rangparam_eta,rangparam_gamma),&
                                                    H_hessOut(rangparam_eta,rangparam_eta)/)
                            endif
                        endif
                    endif
                endif
                
                ! pour R2
                v_chap_R2=0.d0
                t_chap_R2=(/b(rangparam_sigs),b(rangparam_sigt),b(rangparam_sigst)/) ! parametres necessaire: sigma_s,sigma_t,sigma_st 
                v_chap_R2(1,:)=(/H_hessOut(rangparam_sigs,rangparam_sigs),H_hessOut(rangparam_sigs,&
                                rangparam_sigt), H_hessOut(rangparam_sigs,rangparam_sigst)/)
                v_chap_R2(2,:)=(/H_hessOut(rangparam_sigt,rangparam_sigs),H_hessOut(rangparam_sigt,&
                                rangparam_sigt), H_hessOut(rangparam_sigt,rangparam_sigst)/)
                v_chap_R2(3,:)=(/H_hessOut(rangparam_sigst,rangparam_sigs),H_hessOut(rangparam_sigst,&
                                rangparam_sigt), H_hessOut(rangparam_sigst,rangparam_sigst)/)

                moy_tau_boots=0.d0
                moy_R2_boots=0.d0
                ! moy_R2_boots_test=0.d0
                
                call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                allocate(matrice_generation(1,np))
                v_chap_copula(1,1) = H_hessOut(rangparam_copula,rangparam_copula)
                
                do i=1,nboot_kendal
                    ! call rmvnorm(b,H_hessOut,1,0,matrice_generation)
                    ! theta_ST_2(1,1)= matrice_generation(1,rangparam_theta)**2.d0 !theta simule
                    ! gamma_st_2(1,1)= matrice_generation(1,rangparam_gamma)**2.d0  ! gamma simule
                    
                    if(nsim_node(8).ne.3) then 
                        call rmvnorm(t_chap_kendall,v_chap_kendall,1,0,theta_chap_kendall)
                        theta_ST_2(1,1)= theta_chap_kendall(1,1)**2.d0 !theta simule
                        gamma_st_2(1,1)= theta_chap_kendall(1,2)**2.d0  ! gamma simule
                    else
                        call rmvnorm((/b(rangparam_copula)/),v_chap_copula,1,0,theta_chap_copula)
                    endif
                    
                    if(indice_alpha==0 .and. indice_eta==0)then
                        if(nsim_node(8).ne.3) then
                            vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                 method_int_kendal,N_MC_kendall, 1.d0,1.d0,0)    !tau de kendal des non traites z_11=0,z_21=0
                        else
                            ! if(copula_function == 1) vect_kendall_tau(i) = parametre_estimes(s_i-nbre_rejet,1) / &
                                ! (parametre_estimes(s_i-nbre_rejet,1) +1.d0) ! claton
                            ! if(copula_function == 2) vect_kendall_tau(i) = parametre_estimes(s_i-nbre_rejet,1) / &
                                ! (parametre_estimes(s_i-nbre_rejet,1) +2.d0)  ! Gumbel
                            if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                            if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                        endif
                    else
                        if(indice_alpha==1 .and. indice_eta==1)then
                            if(nsim_node(8).ne.3) then
                                vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                    method_int_kendal,N_MC_kendall,&
                                                    theta_chap_kendall(1,4),theta_chap_kendall(1,3),0)
                            else
                                if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                    (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                                if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                    (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                            endif
                        else
                            if(indice_alpha==1)then
                                if(nsim_node(8).ne.3) then
                                    vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                    method_int_kendal,N_MC_kendall,&
                                                    theta_chap_kendall(1,4),1.d0,0)
                                else
                                    if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                        (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                                    if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                        (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                                endif
                            else
                                if(nsim_node(8).ne.3) then
                                    vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,&
                                                    method_int_kendal,N_MC_kendall,&
                                                    1.d0,theta_chap_kendall(1,3),0)
                                else
                                    if(copula_function == 1) vect_kendall_tau(i) = dexp(theta_chap_copula(1,1)) / &
                                        (dexp(theta_chap_copula(1,1)) +2.d0) ! claton
                                    if(copula_function == 2) vect_kendall_tau(i) = theta_chap_copula(1,1)**2.d0 / &
                                        (theta_chap_copula(1,1)**2.d0 +1.d0)  ! Gumbel
                                endif
                            endif
                        endif
                    endif
                    
                    moy_tau_boots=moy_tau_boots+vect_kendall_tau(i)        
                    !R2
                    !call dblepr("Avant ", -1, v_chap_R2, size(v_chap_R2,2)**2)
                    call rmvnorm(t_chap_R2,v_chap_R2,1,0,theta_chap_R2)
                    Chol_R2=0.d0 
                    Chol_R2(1,1)=theta_chap_R2(1,1) ! associe a sigma S
                    Chol_R2(2,2)=theta_chap_R2(1,2) ! associe a sigma T
                    Chol_R2(2,1)=theta_chap_R2(1,3) ! associe a sigma ST
                    mat_A=MATMUL(Chol_R2,TRANSPOSE(Chol_R2))
                    ! vect_R2(i)=(matrice_generation(1,rangparam_sigst)**2.d0)/&
                                    ! (matrice_generation(1,rangparam_sigst)**2.d0 &
                                    ! + matrice_generation(1,rangparam_sigt)**2.d0)
                    vect_R2(i)=(Chol_R2(2,1)**2.d0)/(Chol_R2(2,1)**2.d0+Chol_R2(2,2)**2.d0)
                    moy_R2_boots=moy_R2_boots+vect_R2(i)
                    !!write(18,*)vect_R2(i)
                enddo
                
                deallocate(matrice_generation)
                ! endif
                call percentile_scl(vect_kendall_tau,nboot_kendal,0.025d0,IC_Inf) !borne inf
                call percentile_scl(vect_kendall_tau,nboot_kendal,0.975d0,IC_sup) !borne sup

                !R2
                call percentile_scl(vect_R2,nboot_kendal,0.025d0,IC_Inf_R2) !borne inf
                call percentile_scl(vect_R2,nboot_kendal,0.975d0,IC_sup_R2) !borne sup

                ! sauvegarde des resultats dans un tableau
                result_bootstrap(s_i,1:3)=(/moy_R2_boots/dble(nboot_kendal),IC_Inf_R2,IC_sup_R2/) ! R2
                result_bootstrap(s_i,4:6)=(/moy_tau_boots/dble(nboot_kendal),IC_Inf,IC_sup/)    ! Tau de kendall
                deallocate(theta_ST_2,gamma_st_2)
            endif
            
            
            ! if(nsim_node(8)==1)then ! 1 seul taux de kendall et modele reduit
                ! allocate(theta_ST_2(1,1),gamma_st_2(1,1))
                ! theta_ST_2(1,1)= parametre_estimes(s_i-nbre_rejet,1) !theta estime
                ! gamma_st_2(1,1)= parametre_estimes(s_i-nbre_rejet,15) ! gamma estime
                ! tau_kendal_00=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                ! parametre_estimes(s_i-nbre_rejet,17),parametre_estimes(s_i-nbre_rejet,3),0)    !tau de kendal des non traites z_11=0,z_21=0
                ! moy_kendal_00=moy_kendal_00+tau_kendal_00
                ! moy_tau_boots=0.d0
                ! moy_R2_boots=0.d0
                ! ! moy_R2_boots_test=0.d0
                
                ! call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation
                ! allocate(matrice_generation(1,np))
                
                ! do i=1,nboot_kendal
                    ! call rmvnorm(b,H_hessOut,1,0,matrice_generation)
                    ! theta_ST_2(1,1)= matrice_generation(1,rangparam_theta)**2.d0 !theta simule
                    ! gamma_st_2(1,1)= matrice_generation(1,rangparam_gamma)**2.d0  ! gamma simule
                    
                    ! ! call rmvnorm(t_chap_kendall,v_chap_kendall,1,0,theta_chap_kendall)
                    ! ! theta_ST_2(1,1)= theta_chap_kendall(1,1)**2.d0 !theta simule
                    ! ! gamma_st_2(1,1)= theta_chap_kendall(1,2)**2.d0  ! gamma simule
                    ! if(indice_alpha==0 .and. indice_eta==0)then
                        ! vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                                              ! 1.d0,1.d0,0)    !tau de kendal des non traites z_11=0,z_21=0
                    ! else
                        ! if(indice_alpha==1)then
                            ! ! vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                                            ! ! theta_chap_kendall(1,4),1.d0,0)
                            ! vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                                            ! matrice_generation(1,rangparam_alpha),1.d0,0)
                        ! else
                            ! ! vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                                            ! ! theta_chap_kendall(1,4),theta_chap_kendall(1,3),0)
                            ! vect_kendall_tau(i)=tau_kendall(theta_ST_2,gamma_st_2,sigma_st,0,0,method_int_kendal,N_MC_kendall,&
                                            ! matrice_generation(1,rangparam_alpha),matrice_generation(1,rangparam_eta),0)
                        ! endif
                    ! endif
                    
                    ! moy_tau_boots=moy_tau_boots+vect_kendall_tau(i)        
                    ! !R2
                    ! ! call rmvnorm(t_chap_R2,v_chap_R2,1,0,theta_chap_R2)
                    ! ! Chol_R2=0.d0 
                    ! ! Chol_R2(1,1)=theta_chap_R2(1,1) ! associe a sigma S
                    ! ! Chol_R2(2,2)=theta_chap_R2(1,2) ! associe a sigma T
                    ! ! Chol_R2(2,1)=theta_chap_R2(1,3) ! associe a sigma ST
                    ! ! mat_A=MATMUL(Chol_R2,TRANSPOSE(Chol_R2))
                    ! vect_R2(i)=(matrice_generation(1,rangparam_sigst)**2.d0)/&
                                    ! (matrice_generation(1,rangparam_sigst)**2.d0 &
                                    ! + matrice_generation(1,rangparam_sigt)**2.d0)
                    ! ! vect_R2(i)=(Chol_R2(2,1)**2.d0)/(Chol_R2(2,1)**2.d0+Chol_R2(2,2)**2.d0)
                    ! moy_R2_boots=moy_R2_boots+vect_R2(i)
                    ! !!write(18,*)vect_R2(i)
                ! enddo
                
                ! deallocate(matrice_generation)
                ! ! endif
                ! call percentile_scl(vect_kendall_tau,nboot_kendal,0.025d0,IC_Inf) !borne inf
                ! call percentile_scl(vect_kendall_tau,nboot_kendal,0.975d0,IC_sup) !borne sup
                if(rang_proc==0) then
                    !write(3,*)moy_tau_boots/dble(nboot_kendal),IC_Inf,IC_sup
                    fichier_kendall(s_i,:)=(/moy_tau_boots/dble(nboot_kendal),IC_Inf,IC_sup/)
                endif
                ! if(rang_proc==0) !print*,"moyenne empirique des taux de kendall par bootstrap",moy_tau_boots/dble(nboot_kendal)
                ! if(rang_proc==0) !print*,"vrai tau estime:",tau_kendal_00    
                ! if(rang_proc==0) !print*,"Intervalle de confiance:","[",IC_Inf,"-",IC_sup,"]"                
                
                ! !R2
                ! call percentile_scl(vect_R2,nboot_kendal,0.025d0,IC_Inf_R2) !borne inf
                ! call percentile_scl(vect_R2,nboot_kendal,0.975d0,IC_sup_R2) !borne sup
                if(rang_proc==0) then
                    !write(18,*) moy_R2_boots/dble(nboot_kendal),IC_Inf_R2,IC_sup_R2
                    fichier_R2(s_i,:)=(/moy_R2_boots/dble(nboot_kendal),IC_Inf_R2,IC_sup_R2/)
                endif
                ! if(rang_proc==0) !print*,"moyenne empirique des taux des R2 par bootstrap",moy_R2_boots/dble(nboot_kendal)
                ! if(rang_proc==0) !print*,"vrai R2 estime:",R2_trial    
                ! if(rang_proc==0) !print*,"Intervalle de confiance:","[",IC_Inf_R2,"-",IC_sup_R2,"]"                
                
                ! ! sauvegarde des resultats dans un tableau
                ! result_bootstrap(s_i,1:3)=(/moy_R2_boots/dble(nboot_kendal),IC_Inf_R2,IC_sup_R2/) ! R2
                ! result_bootstrap(s_i,4:6)=(/moy_tau_boots/dble(nboot_kendal),IC_Inf,IC_sup/)    ! Tau de kendall
                ! deallocate(theta_ST_2,gamma_st_2)
            ! endif
            
            !!print*,"suis là 4 s_i=",s_i
            
            !if(nsim_node(8)==1 .or. nsim_node(8)==3)then
            if(nsim_node(8)==1)then
                parametre_estimes(s_i-nbre_rejet,21)=tau_kendal_11    !tau de kendal des traites z_11=1,z_21=1
                parametre_estimes(s_i-nbre_rejet,22)=tau_kendal_10    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                parametre_estimes(s_i-nbre_rejet,23)=tau_kendal_01    !tau de kendal des traite z_11=0,z_21=1
                parametre_estimes(s_i-nbre_rejet,24)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
            endif
            if(nsim_node(8)==2)then
                parametre_estimes(s_i-nbre_rejet,29)=tau_kendal_11    !tau de kendal des traites z_11=1,z_21=1
                parametre_estimes(s_i-nbre_rejet,30)=tau_kendal_10    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                parametre_estimes(s_i-nbre_rejet,31)=tau_kendal_01    !tau de kendal des traite z_11=0,z_21=1
                parametre_estimes(s_i-nbre_rejet,32)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
            endif
            
            if(nsim_node(8)==3)then
                if(nva1 > 1 .or. nva2 >1) then !si plus d'une variable explicative
                    parametre_estimes(s_i-nbre_rejet,22)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
                else
                    parametre_estimes(s_i-nbre_rejet,24)=tau_kendal_00    !tau de kendal des non traites z_11=0,z_21=0
                endif
                ! SE de kendall_tau par delta method. voir cahier le 18/04/2019 pour demonstration
                if(copula_function == 1) then
                    if(nva1 > 1 .or. nva2 >1) then !si plus d'une variable explicative
                        parametre_estimes(s_i-nbre_rejet,23) =  2.d0 * dexp(b(rangparam_copula))/&
                            (dexp(b(rangparam_copula)) + 2.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                            rangparam_copula))! se_tau_kendal
                    else
                        parametre_estimes(s_i-nbre_rejet,25) =  2.d0 * dexp(b(rangparam_copula))/&
                            (dexp(b(rangparam_copula)) + 2.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                            rangparam_copula))! se_tau_kendal
                    endif
                    bi_sigmas = dexp(b(rangparam_copula)) - 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                    bs_sigmas = dexp(b(rangparam_copula)) + 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                    !taux de couverture
                    vrai_tau_copula = thetacopule/(thetacopule + 2.d0)
                endif
                if(copula_function == 2) then
                    if(nva1 > 1 .or. nva2 >1) then !si plus d'une variable explicative
                        parametre_estimes(s_i-nbre_rejet,23) =   2.d0 * b(rangparam_copula)/&
                            (b(rangparam_copula)**2.d0 + 1.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                            rangparam_copula))! se_tau_kendal
                    else
                        parametre_estimes(s_i-nbre_rejet,25) =   2.d0 * b(rangparam_copula)/&
                            (b(rangparam_copula)**2.d0 + 1.d0)**2.d0 * dsqrt(H_hessOut(rangparam_copula,&
                            rangparam_copula))! se_tau_kendal
                    endif
                    bi_sigmas = b(rangparam_copula)**2.d0 - 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                    bs_sigmas = b(rangparam_copula)**2.d0  + 1.96d0*parametre_estimes(s_i-nbre_rejet,25)
                    !taux de couverture
                    vrai_tau_copula = thetacopule/(thetacopule + 1.d0)
                endif
                if(vrai_tau_copula >= bi_sigmas .and. vrai_tau_copula <= bs_sigmas)then ! taux de couverture
                    taux_couverture_tauk=taux_couverture_tauk+1.d0
                endif    
                
            endif
            
        endif
        
        ! !print*," Attention!!! sauvegarde des resultats, nb_processus=",nb_processus,"rang=",rang
        erreur_fichier=0
        if(nb_processus<=1)then ! on ecrit le resultat dans le fichier que si l'on n'a pas plusieurs processus. si oui l'ecriture est faite apres la sortie de la boucle par le processus de rang 0
            ! !print*,"ce que je veux suvagarder: ",s_i-nbre_rejet
            ! !print*,parametre_estimes(s_i-nbre_rejet,:)
            
            997 continue
            !write(11,*,iostat=erreur_fichier)parametre_estimes(s_i-nbre_rejet,:)
            param_estimes(s_i,:) = parametre_estimes(s_i-nbre_rejet,:)
            if(erreur_fichier .ne.0) then
                ! if(rang_proc==0) !print*,"ATTENTION!! erreur d'ecriture des parametres estimes dans le fichier de sortie",&
                ! "oninsiste jusqu'a l'ecriture"
                goto 997
            endif
            ! if(rang_proc==0) !write(20,*,iostat=erreur_fichier)parametre_estimes(s_i-nbre_rejet,:)
        endif
        ! !print*,"Fin sauvegarde des parametres estimes,erreur_fichier=",erreur_fichier
        ! sauvegarde des donnees a ranger dans le fichier resultat
        
        if(nb_processus>1)then ! on ne le fait que dans on a plus d'un processus
            parametre_estimes_MPI(indice_sim_proc,:)=parametre_estimes(s_i-nbre_rejet,:) ! juste pour le MIP
            indice_sim_proc=indice_sim_proc+1
        endif

    endif
    
    1002 continue
    s_i=s_i+1
    !deallocate(b)
    !deallocate(H_hessOut,HIHOut)
    deallocate(vdeces,vsurrogate,paGH,don_simulS,tableEssai)
    !close(2)
    !close(4)
    !if(s_i.eq.seed_) then
    !    close(9)
    !    close(10)
    !endif
    !!print*,"suis là s_i=",s_i
end do
    EPS2 = EPS ! on restitue les critères de convergence pour le dernier jeux de donnees

    ! ecart-type des parametres simules
    se_theta_sim=dsqrt(variance(tab_var_theta))
    if(nsim_node(8).ne.0)then
        se_sigmas_sim= dsqrt(variance(tab_var_sigma(:,1)))
        se_sigmat_sim= dsqrt(variance(tab_var_sigma(:,2)))
        se_rho_sim= dsqrt(variance(tab_var_sigma(:,3)))
        se_gamma_sim=dsqrt(variance(tab_var_sigma(:,4)))
        if(nsim_node(8)==2)then
            se_theta_simt=dsqrt(variance(tab_var_sigma(:,5)))
            se_rho_sim_wij=dsqrt(variance(tab_var_sigma(:,6)))
            se_gamma_simt=dsqrt(variance(tab_var_sigma(:,7)))
            se_rho_sim_ui=dsqrt(variance(tab_var_sigma(:,8)))
        endif
    endif
    
    ! ecart-type des parametres estime(SD)
    se_theta_est=dsqrt(variance(parametre_estimes(:,1)))
    if(indice_eta==1) se_eta_est=dsqrt(variance(parametre_estimes(:,3)))
    se_beta_s=dsqrt(variance(parametre_estimes(:,5)))
    se_beta_t=dsqrt(variance(parametre_estimes(:,7)))
    if(nsim_node(8).ne.0)then
        se_sigmas_est= dsqrt(variance(parametre_estimes(:,9)))
        se_sigmat_est= dsqrt(variance(parametre_estimes(:,11)))
        se_cov_est= dsqrt(variance(parametre_estimes(:,13)))
        se_gamma_est=dsqrt(variance(parametre_estimes(:,15)))
        if(indice_alpha==1) se_alpha_est=dsqrt(variance(parametre_estimes(:,17)))
        
        if(nsim_node(8)==1)then
            se_R2_trial=dsqrt(variance(parametre_estimes(:,19)))
            se_kendal_11=dsqrt(variance(parametre_estimes(:,21)))
            se_kendal_10=dsqrt(variance(parametre_estimes(:,22)))
            se_kendal_01=dsqrt(variance(parametre_estimes(:,23)))
            se_kendal_00=dsqrt(variance(parametre_estimes(:,24)))
        endif
        
        if(nsim_node(8)==2)then
            se_thetat_est=dsqrt(variance(parametre_estimes(:,19)))
            se_cov_est_wij=dsqrt(variance(parametre_estimes(:,21)))
            se_gamma_estt=dsqrt(variance(parametre_estimes(:,23)))
            se_cov_est_ui=dsqrt(variance(parametre_estimes(:,25)))
            se_R2_trial=dsqrt(variance(parametre_estimes(:,27)))
            se_kendal_11=dsqrt(variance(parametre_estimes(:,29)))
            se_kendal_10=dsqrt(variance(parametre_estimes(:,30)))
            se_kendal_01=dsqrt(variance(parametre_estimes(:,31)))
            se_kendal_00=dsqrt(variance(parametre_estimes(:,32)))
        endif
    endif
    
    if(nb_processus>1)then ! on ne le fait que dans on a plus d'un processus
        n_sim_exact=dble(n_sim-nbre_rejet)
        ! envois des parametres estimes au processus de rang 0
        ! if(rang_proc==0) !print*
        ! if(rang_proc==0) !print*
        ! !call MPI_GATHER(parametre_estimes_MPI(1:tableNsim(rang+1),:),tableNsim(rang+1), MPI_DOUBLE_PRECISION,parametre_estimes_MPI_T,&
                        ! n_sim_total,MPI_DOUBLE_PRECISION,0, comm,code)
        if (rang==0) then
            init_i=1 ! ce processus commence a la premiere simulation
        else
            init_i=sum(tableNsim(1:rang))+1 ! ce processus commence a la simulation qui respecte son ordre et doit s'arreter au nombre de simultation dont il a le droit d'executer
        endif
        max_i=init_i+tableNsim(rang+1)-1
        ! ce qui suit pour l'envoie des element de la matrice vecteur par vecteur
        allocate(tampon(tableNsim(rang+1)*size(parametre_estimes_MPI,2)))
        allocate(tampon_all(n_sim_total*size(parametre_estimes_MPI,2)))
        
        !sauvegarde des element de la matrice ligne par ligne
        k=1
        do i=1,tableNsim(rang+1)
        !do i=1,n_sim_exact
            j=k+size(parametre_estimes_MPI,2)
            tampon(k:j)=parametre_estimes_MPI(i,:)
            k=k+size(parametre_estimes_MPI,2)
        enddo
        
        if(rang==0)then
            n_sim_exact=n_sim_exact_0
            ! if(rang_proc==0) !print*, "taux de convergence=",n_sim_exact*100.d0/n_sim_total,"(",n_sim_exact,"/",n_sim_total,")"
            !reccuperation des element sauvegardes dans la matrice parametre_estimes_MPI_T ligne par ligne
            k=1
            !do i=1,n_sim_total
            i=1
            do while (i<=n_sim_total)
            !do i=1,max_i ! avec max_i le nombre de simulation qui a converge
                j=k+size(parametre_estimes_MPI,2)
                !!print*,"sum(tampon_all(k:j))",sum(tampon_all(k:j))
                !if(sum(tampon_all(k:j)).ne.0.d0) then ! on ne recupere pas les donnees qui n'ont pas converge (et dont les 00000...)
                    parametre_estimes_MPI_T(i,:)=tampon_all(k:j)
                    !sauvegarde tout les resultats dans un seul fichier cree et gere par le processus de rang 0
                    ! if(nb_processus>1) !write(11,*)parametre_estimes_MPI_T(i,:)
                    ! if(rang_proc==0) !print*,"donnees fusionnees",parametre_estimes_MPI_T(i,:)
                    i=i+1
                !endif
                k=k+size(parametre_estimes_MPI,2)
            enddo
            
            ! recherche des lignes correspondantes aux simulations qui n'ont pas convergees
            deallocate(parametre_estimes)
            allocate(parametre_estimes(INT(n_sim_exact),size(parametre_estimes_MPI,2)))
            
            k=1
            do i=1,n_sim_total
                if(sum(parametre_estimes_MPI_T(i,:)).ne.0.d0) then
                    !mise à jour des parametres estimes a partir des valeurs issues de la fusion
                    parametre_estimes(k,:)=parametre_estimes_MPI_T(i,:)
                    k=k+1
                endif
            enddo
            
            moy_theta_est=moy_theta_est_0
            moy_ni=moy_ni_0
            moy_trt=moy_trt_0
            moy_theta=moy_theta_0
            moy_sigmas=moy_sigmas_0
            moy_sigmat=moy_sigmat_0
            moy_rho=moy_rho_0
            moy_gamma=moy_gamma_0
            moy_pros=moy_pros_0
            moy_dec=moy_dec_0
            tab_var_sigma=tab_var_sigma_0
            moy_se_theta=moy_se_theta_0
            taux_couverture_theta=taux_couverture_theta_0
            moy_gamma_est=moy_gamma_est_0
            moy_se_gamma=moy_se_gamma_0
            taux_couverture_gamma=taux_couverture_gamma_0
            moy_sigmas_est=moy_sigmas_est_0
            moy_se_sigmas=moy_se_sigmas_0
            taux_couverture_sigmas=taux_couverture_sigmas_0
            moy_sigmat_est=moy_sigmat_est_0
            moy_se_sigmat=moy_se_sigmat_0
            taux_couverture_sigmat=taux_couverture_sigmat_0
            moy_sigmast_est=moy_sigmast_est_0
            moy_se_sigmast=moy_se_sigmast_0
            taux_couverture_sigmast=taux_couverture_sigmast_0
            moy_eta=moy_eta_0
            moy_se_eta=moy_se_eta_0
            taux_couverture_eta=taux_couverture_eta_0
            moy_betaS=moy_betaS_0
            moy_betaS_se=moy_betaS_se_0
            taux_couvertureS=taux_couvertureS_0
            moy_betaT=moy_betaT_0
            moy_betaT_se=moy_betaT_se_0
            taux_couvertureT=taux_couvertureT_0
            nbre_rejet=nbre_rejet_0
            se_theta_est=se_theta_est_0
            se_eta_est=se_eta_est_0
            se_beta_s=se_beta_s_0
            se_beta_t=se_beta_t_0
            se_sigmas_est=se_sigmas_est_0
            se_sigmat_est=se_sigmat_est_0
            se_cov_est=se_cov_est_0
            se_gamma_est=se_gamma_est_0
            se_alpha_est=se_alpha_est_0
            se_thetat_est=se_thetat_est_0
            se_cov_est_wij=se_cov_est_wij_0
            se_gamma_estt=se_gamma_estt_0
            se_cov_est_ui=se_cov_est_ui_0
            moy_thetat=moy_thetat_0
            se_theta_simt=se_theta_simt_0
            moy_thetat_est=moy_thetat_est_0
            moy_se_thetat=moy_se_thetat_0
            taux_couverture_thetat=taux_couverture_thetat_0
            moy_rho_wij=moy_rho_wij_0
            se_rho_sim_wij=se_rho_sim_wij_0
            moy_thetast_est=moy_thetast_est_0
            se_cov_est=se_cov_est_0
            moy_se_thetast=moy_se_thetast_0
            taux_couverture_thetast=taux_couverture_thetast_0
            moy_gammat=moy_gammat_0
            se_gamma_estt=se_gamma_estt_0
            se_gamma_simt=se_gamma_simt_0
            moy_gammat_est=moy_gammat_est_0
            moy_se_gammat=moy_se_gammat_0
            taux_couverture_gammat=taux_couverture_gammat_0
            moy_rho_ui=moy_rho_ui_0
            se_rho_sim_ui=se_rho_sim_ui_0
            moy_gammast_est=moy_gammast_est_0
            se_cov_est_ui=se_cov_est_ui_0
            moy_se_gammast=moy_se_gammast_0
            taux_couverture_gammast=taux_couverture_gammast_0
            moy_alpha=moy_alpha_0
            se_alpha_est=se_alpha_est_0
            moy_se_alpha=moy_se_alpha_0
            taux_couverture_alpha=taux_couverture_alpha_0
            
        endif
        
        deallocate(tampon,tampon_all)
        
    endif    
    !else ! si l'on a qu'un seul processus alors on suppose que l'on fait que du OpenMP et on fait donc ce qui suit:
    if(rang==0)then
        if(nb_processus<=1) n_sim_exact=dble(n_sim-nbre_rejet)
        if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then
            if(nsim_node(4).eq.1) then
                if(rang_proc==0) then
                    if(logNormal==1)then
                        !print*,"==============================================================="
                        !print*,"Simulation lognormale et extimation par spline et gauss Hermite" 
                        !print*,"==============================================================="
                    else
                        !print*,"==============================================================="
                        !print*,"Simulation Gamma et extimation par spline et gauss Laguerre"
                        !print*,"==============================================================="
                    endif
                endif
            endif
            if(nsim_node(4).eq.2) then
                if(rang_proc==0) then
                    if(logNormal==1)then
                        !print*,"==============================================================="
                        !print*,"Simulation lognormale et extimation par spline et gauss Hermite individuel+Monte-carlo essai" 
                        !print*,"==============================================================="
                    else
                        !print*,"cas impossible pour ce type de modele"
                    endif
                endif
            endif
        else
            if(rang_proc==0) then
                if(logNormal==1)then
                    !print*,"==============================================================="
                    !print*,"Simulation lognormale et extimation par spline et Monte-carlo" 
                    !print*,"==============================================================="
                else
                    !print*,"==============================================================="
                    !print*,"Simulation Gamma et extimation par spline et Monte-carlo"
                    !print*,"==============================================================="
                endif
            endif
        endif
        if(rang_proc==0) then
            !print*,""
            !print*,"Paramètres de simulation"
            !print*,"-----------------------"
            !print*,"N_sujet:   ",nsujet, "n_essai:",n_essai
            !print*,"Theta(S):     ",theta2
            if(nsim_node(8)==2)then 
                !print*,"Theta_T:     ",theta2_t
                !print*,"rho_w_ST:     ",rsqrt_theta
            endif
        
            if(nsim_node(8).ne.0)then !model conjoint surrogate
                !print*,"sigma_s:     ",sigma_s
                !print*,"sigma_t:     ",sigma_t
                !print*,"rho_v_st:     ",rsqrt
                if(frailt_base==1) then
                    !print*,"gamma(S):        ",gamma_ui
                    ! if(indice_alpha==1) !print*,"alpha:        ",alpha_ui
                endif
                if(nsim_node(8)==2)then 
                    !print*,"gamma_T:        ",gamma_uit
                    !print*,"rho_u_st:     ",rsqrt_gamma_ui
                endif
            endif
            !print*,"lambdas:   ",lambdas
            !print*,"nus:       ",nus
            !print*,"lambdat:   ",lambdat
            !print*,"nut:       ",nut
            !print*,"beta_S:    ",betas
            !print*,"beta_T:    ",betat
            !print*,"censA:     ",temps_cens
            ! if(nb_processus<=1) !print*,"N_sim:     ",n_sim
            ! if(nb_processus>1) !print*,"N_sim:     ",n_sim_total
            !print*,""
            !print*,"Paramètres d'estimation"
            !print*,"-----------------------"
            if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then
                if(nsim_node(4).eq.1) then
                    !print*,"nbre_point:",nsim_node(2)
                endif
                if(nsim_node(4).eq.2) then
                    !print*,"nbre_point:",nsim_node(2)
                    !print*,"nbre bouble MC:",nsim_node(1)
                endif
            else
                !print*,"nbre bouble MC:",nsim_node(1)
            endif
            !print*,"kapa_use:  ",kapa_use
            !print*,"noeuds utilise:",nz
            !!print*,"valeurs initiales param de variances:",theta_init,sigma_ss_init,sigma_tt_init,sigma_st_init,gamma_init
            !print*,"Nombre moyen d'itterarion pour la convergence:",moy_ni/n_sim_exact
            ni=int(moy_ni/n_sim_exact)
            ! call intpr("Nombre itteration:", -1, ni, 1)
            !print*,""

            call cpu_time(tp2)
            !!print*,'************temps_CPU en secondes******* =',tp2-tp1
            
            !print*,"==============================================================="
            !print*,"===================Moyenne empirique des estimations==========="
            !print*,"==============================================================="
            
            !print*,"parametres empiriques"
            !print*,"*********************"
        endif
            if(gener_only.eq.1) then
                nbre_rejet=dble(0)    
                if(nb_processus<=1)then 
                    n_sim_exact=dble(n_sim-nbre_rejet)
                else
                    n_sim_exact=dble(n_sim_total-nbre_rejet)
                endif
            if(rang_proc==0) then
                !print*,"moyenne empirique des personnes traitées= ",    moy_trt/n_sim_exact
                !print*,"moyenne empirique des theta(S)_sim= ",moy_theta/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types des theta(S)_sim= ",se_theta_sim
                if(nsim_node(8)==2)then
                    !print*,"moyenne empirique des thetaT_sim= ",moy_thetat/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types des thetaT_sim= ",se_theta_simt
                    !print*,"moyenne empirique des rho_wij_ST_sim= ",moy_rho_wij/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types des rho_wij_ST_sim= ",se_rho_sim_wij
                endif
                
                if(nsim_node(8).ne.0)then !model conjoint surrogate
                    !print*,"moyenne empirique des sigmas_sim= ",moy_sigmas/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types des sigmas_sim= ",se_sigmas_sim
                    !print*,"moyenne empirique des sigmat_sim= ",moy_sigmat/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types des sigmat_sim= ",se_sigmat_sim
                    !print*,"moyenne empirique des rho_sim= ",moy_rho/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types des rho_sim= ",se_rho_sim
                    
                    if(frailt_base==1) then
                        !print*,"moyemme empiriques des gamma(S)=",moy_gamma/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-type empirique des gamma(S)=",se_gamma_sim
                        if(nsim_node(8)==2)then
                            !print*,"moyenne empirique des gammaT_sim= ",moy_gammat/n_sim_exact
                            ! if(n_sim_exact>1) !print*,"ecart-types des gammaT_sim= ",se_gamma_simt
                            !print*,"moyenne empirique des rho_ui_ST_sim= ",moy_rho_ui/n_sim_exact
                            ! if(n_sim_exact>1) !print*,"ecart-types des rho_ui_ST_sim= ",se_rho_sim_ui
                        endif
                    endif
                    
                    ! if(nsim_node(8)==2)then
                        ! !print*,"moyenne empirique des sigmas_sim= ",moy_sigmas/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des sigmas_sim= ",se_sigmas_sim
                        ! !print*,"moyenne empirique des sigmat_sim= ",moy_sigmat/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des sigmat_sim= ",se_sigmat_sim
                        ! !print*,"moyenne empirique des rho_sim= ",moy_rho/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des rho_sim= ",se_rho_sim
                    ! endif
                endif
                !print*,"moyenne empirique des proportion progression=",moy_pros/n_sim_exact
                !print*,"moyenne empirique des proportions décès= ",moy_dec/n_sim_exact
            endif
                goto 1001
            endif
        
        
        if(nb_processus<=1) n_sim_exact=dble(n_sim-nbre_rejet) !nombre de simulations qui ont converges
        
        if(rang_proc==0) then
            !print*,"moyenne empirique des personnes traitées=",    moy_trt/n_sim_exact
            if(une_donnee==0)then !si on ne genere pas les donnees, pas la peine de les afficher
                ! if(rang_proc==0) !print*,"moyenne empirique des theta(S)_sim=",moy_theta/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types des theta_sim=",se_theta_sim
                    if(nsim_node(8)==2)then
                        !print*,"moyenne empirique des thetaT_sim= ",moy_thetat/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des thetaT_sim= ",se_theta_simt
                        !print*,"moyenne empirique des rho_wij_ST_sim= ",moy_rho_wij/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des rho_wij_ST_sim= ",se_rho_sim_wij
                    endif
                        
                    if(nsim_node(8).ne.0)then !model conjoint surrogate
                        !print*,"moyenne empirique des sigmas_sim= ",moy_sigmas/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des sigmas_sim= ",se_sigmas_sim
                        !print*,"moyenne empirique des sigmat_sim= ",moy_sigmat/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des sigmat_sim= ",se_sigmat_sim
                        !print*,"moyenne empirique des rho_sim= ",moy_rho/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types des rho_sim= ",se_rho_sim
                        
                        if(frailt_base==1) then
                            !print*,"moyemme empiriques des gamma(S)=",moy_gamma/n_sim_exact
                            ! if(n_sim_exact>1) !print*,"ecart-type empirique des gamma=",se_gamma_sim
                            if(nsim_node(8)==2)then
                                !print*,"moyenne empirique des gammaT_sim= ",moy_gammat/n_sim_exact
                                ! if(n_sim_exact>1) !print*,"ecart-types des gammaT_sim= ",se_gamma_simt
                                !print*,"moyenne empirique des rho_ui_ST_sim= ",moy_rho_ui/n_sim_exact
                                ! if(n_sim_exact>1) !print*,"ecart-types des rho_ui_ST_sim= ",se_rho_sim_ui
                            endif
                        endif
                    endif
            endif
            !print*,"moyenne empirique des proportion progression=",moy_pros/n_sim_exact
            !print*,"moyenne empirique des proportions décès=",moy_dec/n_sim_exact
            !!print*,"moyenne empirique des proportions censure=",moy_cens
            
            !print*,""
            !print*,"Parametres estimés"
            !print*,"******************"
            !print*,"theta(S)",theta2
            !print*,"moyenne des theta(S) estimés=",moy_theta_est/n_sim_exact
            ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des theta(S) estimés=",se_theta_est
            !print*,"moyenne des ecart-types des theta(S) estimés=",moy_se_theta/n_sim_exact
            !print*,"taux de couverture des theta(S)=",taux_couverture_theta*100.d0/n_sim_exact
        
            if(nsim_node(8)==2)then
                !print*,""
                !print*,"theta_T",theta2_t
                !print*,"moyenne des theta_T estimés=",moy_thetat_est/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des theta_T estimés=",se_thetat_est
                !print*,"moyenne des ecart-types des theta_T estimés=",moy_se_thetat/n_sim_exact
                !print*,"taux de couverture des theta_T=",taux_couverture_thetat*100.d0/n_sim_exact
                
                !print*,""
                !print*,"theta_ST",thetast_vrai
                !print*,"moyenne des theta T estimés=",moy_thetast_est/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des theta_ST estimés=",se_cov_est_wij
                !print*,"moyenne des ecart-types des theta_ST estimés=",moy_se_thetast/n_sim_exact
                !print*,"taux de couverture des theta_ST=",taux_couverture_thetast*100.d0/n_sim_exact
            endif
        
            if(indice_eta==1)then
                !print*,""
                !print*,"eta",eta
                !print*,"moyenne des eta estimés=",moy_eta/n_sim_exact
                ! ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des eta estimés=",se_eta_est
                !print*,"moyenne des ecart-types des eta estimés=",moy_se_eta/n_sim_exact
                !print*,"taux de couverture des eta=    ",taux_couverture_eta*100.d0/n_sim_exact
            endif
        
            if(nsim_node(8).ne.0)then !model conjoint surrogate
                if(frailt_base==1) then
                    !print*,""
                    !print*,"gamma(S)",gamma_ui
                    !print*,"moyenne des gamma(S) estimés=",moy_gamma_est/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des gamma(S) estimés=",se_gamma_est
                    !print*,"moyenne des ecart-types des gamma(S) estimés=",moy_se_gamma/n_sim_exact
                    !print*,"taux de couverture des gamma(S)=",taux_couverture_gamma*100.d0/n_sim_exact
                    
                    if(indice_alpha==1)then
                        !print*,""
                        !print*,"Alpha",alpha_ui
                        !print*,"moyenne des alpha estimés=",moy_alpha/n_sim_exact
                        ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des alpha estimés=",se_alpha_est
                        !print*,"moyenne des ecart-types des alpha estimés=",moy_se_alpha/n_sim_exact
                        !print*,"taux de couverture des alpha=    ",taux_couverture_alpha*100.d0/n_sim_exact
                    endif
                endif
            
                if(nsim_node(8)==2)then
                    !print*,""
                    !print*,"gamma_T",gamma_uit
                    !print*,"moyenne des gamma_T estimés=",moy_gammat_est/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des gamma_T estimés=",se_gamma_estt
                    !print*,"moyenne des ecart-types des gamma_T estimés=",moy_se_gammat/n_sim_exact
                    !print*,"taux de couverture des gamma_T=",taux_couverture_gammat*100.d0/n_sim_exact
                    
                    !print*,""
                    !print*,"gamma_ST",gammast_vrai
                    !print*,"moyenne des gamma_ST estimés=",moy_gammast_est/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des gamma_ST estimés=",se_cov_est_ui
                    !print*,"moyenne des ecart-types des gamma_ST estimés=",moy_se_gammast/n_sim_exact
                    !print*,"taux de couverture des gamma_ST=",taux_couverture_gammast*100.d0/n_sim_exact
                endif
        
                !print*,""
                !print*,"sigma_S",sigma_s
                !print*,"moyenne des sigma_S estimés=",moy_sigmas_est/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des sigma_S estimés=",se_sigmas_est
                !print*,"moyenne des ecart-types des sigma_S estimés=",moy_se_sigmas/n_sim_exact
                !print*,"taux de couverture des sigma_S=",taux_couverture_sigmas*100.d0/n_sim_exact
                
                !print*,""
                !print*,"sigma_T",sigma_t
                !print*,"moyenne des sigma_T estimés=",moy_sigmat_est/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des sigma_T estimés=",se_sigmat_est
                !print*,"moyenne des ecart-types des sigma_T estimés=",moy_se_sigmat/n_sim_exact
                !print*,"taux de couverture des sigma_T=",taux_couverture_sigmat*100.d0/n_sim_exact
                
                !print*,""
                !print*,"sigma_ST",sigmast_vrai
                !print*,"moyenne des sigma_ST estimés=",moy_sigmast_est/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des sigma_ST estimés=",se_cov_est
                !print*,"moyenne des ecart-types des sigma_ST estimés=",moy_se_sigmast/n_sim_exact
                !print*,"taux de couverture des sigma_ST=",taux_couverture_sigmast*100.d0/n_sim_exact
            endif
        
            !print*,""
            !print*,"beta_S",betas
            !print*,"moyenne des beta_S estimés=    ",moy_betaS/n_sim_exact
            ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des beta_S estimés=",se_beta_s
            !print*,"moyenne des ecart-types des beta_S estimés=",moy_betaS_se/n_sim_exact
            !print*,"taux de couverture des beta_S=",taux_couvertureS*100.d0/n_sim_exact
            
            !print*,""
            !print*,"beta_T",betat
            !print*,"moyenne des beta_T estimés=",moy_betaT/n_sim_exact
            ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des beta_T estimés=",se_beta_t
            !print*,"moyenne des ecart-types des beta_T estimés=",moy_betaT_se/n_sim_exact
            !print*,"taux de couverture des beta_T=",taux_couvertureT*100.d0/n_sim_exact
        endif
        ! parametres de validation du surrogate
        if(nsim_node(8).ne.0)then
            if(rang_proc==0) then
                !print*,""
                !print*,"R2_trial=",rsqrt**2
                !print*,"moyenne des R2_trial estimés=",moy_R2_trial/n_sim_exact
                ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des R2_trial =",se_R2_trial
                !print*,"moyenne des ecart-types (SE) des R2_trial estimés=",moy_se_R2_trial /n_sim_exact
                !print*,"taux de couverture des R2_trial=",taux_couverture_R2_trial*100.d0/n_sim_exact
                !print*,"moyenne des IC de R2_trial",moy_bi_R2_trial/n_sim_exact,";",moy_bs_R2_trial/n_sim_exact
            endif
            ! calcul des vrais taux de kendall estimes
            if(nsim_node(8).eq.2)then ! modele complet
                theta_ST0(:,1)= (/theta2,thetast_vrai/)
                theta_ST0(:,2)= (/thetast_vrai,theta2_t/)
                gamma_st0(:,1)= (/gamma_ui,gammast_vrai/)
                gamma_st0(:,2)= (/gammast_vrai,gamma_uit/)
                sigma_st0(:,1)= (/sigma_s,sigmast_vrai/)
                sigma_st0(:,2)= (/sigmast_vrai,sigma_t/)
            endif
            
            if((nsim_node(8).eq.1) .or. (nsim_node(8).eq.3))then ! modele reduit
                allocate(theta_ST0_2(1,1),gamma_st0_2(1,1))
                theta_ST0_2(1,1)= theta2
                gamma_st0_2(1,1)= gamma_ui
                if(indice_eta==0) eta=1.d0
                if(indice_alpha==0) alpha_ui=1.d0
            endif
            
                ! !print*,"theta_ST_2",theta_ST0_2
                ! !print*,"gamma_st_2",gamma_st0_2
                ! !print*,"sigma_st",sigma_st0
                ! !print*,"method_int_kendal",method_int_kendal
                ! !print*,"alpha",alpha_ui
                ! !print*,"zeta",eta
            
            if(nsim_node(8).eq.1)then ! modele reduit
                if(method_int_kendal==4 .or. method_int_kendal==5)then ! 1 seul taux de kendall
                    tau_kendal_00=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,0,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des non traites z_11=0,z_21=0
                else
                    tau_kendal_11=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,1,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des traites z_11=1,z_21=1
                    tau_kendal_10=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,1,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                    tau_kendal_01=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,0,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des traite z_11=0,z_21=1
                    tau_kendal_00=tau_kendall(theta_ST0_2,gamma_st0_2,sigma_st0,0,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,0)!tau de kendal des non traites z_11=0,z_21=0
                endif
            endif
            if(nsim_node(8).eq.3)then ! modele reduit copule
                tau_kendal_00=vrai_tau_copula
            endif
            
            deallocate(theta_ST0_2,gamma_st0_2)
            
            if(nsim_node(8).eq.2)then ! modele complet
                if(method_int_kendal==4)then ! 1 seul taux de kendall
                    tau_kendal_00=tau_kendall(theta_ST0,gamma_st0,sigma_st0,0,0,method_int_kendal,&
                    N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des non traites z_11=0,z_21=0
                else
                    tau_kendal_11=tau_kendall(theta_ST0,gamma_st0,sigma_st0,1,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des traites z_11=1,z_21=1
                    tau_kendal_10=tau_kendall(theta_ST0,gamma_st0,sigma_st0,1,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
                    tau_kendal_01=tau_kendall(theta_ST0,gamma_st0,sigma_st0,0,1,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des traite z_11=0,z_21=1
                    tau_kendal_00=tau_kendall(theta_ST0,gamma_st0,sigma_st0,0,0,method_int_kendal,&
                                  N_MC_kendall,alpha_ui,eta,1)    !tau de kendal des non traites z_11=0,z_21=0
                endif
            endif
            
            if(rang_proc==0) then
                if(method_int_kendal==4 .or. method_int_kendal==5)then ! 1 seul taux de kendall et MC
                    !print*,""
                    !print*,"tau_kendal =",tau_kendal_00
                    !print*,"moyenne des tau_kendal estimés=",moy_kendal_00/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des tau_kendal =",se_kendal_00
                    !print*,"biais=",tau_kendal_00-moy_kendal_00/n_sim_exact
                else 
                    !print*,""
                    !print*,"tau_kendal_11=",tau_kendal_11
                    !print*,"moyenne des tau_kendal_11 estimés=",moy_kendal_11/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des tau_kendal_11 =",se_kendal_11
                    !print*,"biais=",tau_kendal_11-moy_kendal_11/n_sim_exact
                    
                    !print*,""
                    !print*,"tau_kendal_10=",tau_kendal_10
                    !print*,"moyenne des tau_kendal_10 estimés=",moy_kendal_10/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des tau_kendal_10 =",se_kendal_10
                    !print*,"biais=",tau_kendal_10-moy_kendal_10/n_sim_exact
                    
                    !print*,""
                    !print*,"tau_kendal_01=",tau_kendal_01
                    !print*,"moyenne des tau_kendal_01 estimés=",moy_kendal_01/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des tau_kendal_01 =",se_kendal_01
                    !print*,"biais=",tau_kendal_01-moy_kendal_01/n_sim_exact
                    
                    !print*,""
                    !print*,"tau_kendal_00=",tau_kendal_00
                    !print*,"moyenne des tau_kendal_00 estimés=",moy_kendal_00/n_sim_exact
                    ! if(n_sim_exact>1) !print*,"ecart-types empirique (SD) des tau_kendal_00 =",se_kendal_00
                    !print*,"biais=",tau_kendal_00-moy_kendal_00/n_sim_exact
                endif
                ! recherche des taux de couverture par bootsrap
                CP_R2_boot=0.d0
                CP_ktau_boot=0.d0
                remplnsim=INT(n_sim_exact) ! jute pour avoir un INt pour le Do
                
                do i=1,remplnsim
                    if((result_bootstrap(i,2)<=rsqrt**2) .and. (result_bootstrap(i,3)>=rsqrt**2)) CP_R2_boot=&
                        CP_R2_boot+1
                    if((result_bootstrap(i,5)<=tau_kendal_00) .and. (result_bootstrap(i,6)>=tau_kendal_00)) &
                        CP_ktau_boot=CP_ktau_boot+1
                enddo
                if(rang_proc==0) then
                    !print*,"=============Taux de couverture des parametres de validation par bootstrap===="
                    !print*,"CP R2(en %):",CP_R2_boot*100.d0/n_sim_exact
                    !print*,"CP Kendall tau(en %):",CP_ktau_boot*100.d0/n_sim_exact
                endif
            endif
            
        endif
        
        if(rang_proc==0) then
            !print*,""
            !print*,"Kappa utilise surrogate et true=",ax1,ax2
            !print*,"nbre_rejet:",nbre_rejet
            !print*,"Taux de rejet:",nbre_rejet*100.d0/(n_sim_exact+nbre_rejet)
        
        !===============sauvegarde des resultats de simulation dans le fichier=============
        !==================================================================================
            !open(17,file="Resultat_simulation.txt")
                if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then
                    if(nsim_node(4).eq.1) then
                        if(logNormal==1)then
                            !write(17,*)"==============================================================="
                            !write(17,*)"Simulation lognormale et extimation par spline et gauss Hermite" 
                            !write(17,*)"==============================================================="
                        else
                            !write(17,*)"==============================================================="
                            !write(17,*)"Simulation Gamma et extimation par spline et gauss Laguerre"
                            !write(17,*)"==============================================================="
                        endif
                    endif
                    if(nsim_node(4).eq.2) then
                        if(logNormal==1)then
                            !write(17,*)"==============================================================="
                            !write(17,*)"Simulation lognormale et extimation par spline et GH individuel+MC essai" 
                            !write(17,*)"==============================================================="
                        else
                            !write(17,*)"cas impossible pour ce type de modele"
                        endif
                    endif
            else
                if(logNormal==1)then
                    !write(17,*)"==============================================================="
                    !write(17,*)"Simulation lognormale et extimation par spline et Monte-carlo" 
                    !write(17,*)"==============================================================="
                else
                    !write(17,*)"==============================================================="
                    !write(17,*)"Simulation Gamma et extimation par spline et Monte-carlo"
                    !write(17,*)"==============================================================="
                endif
            endif
            
            !write(17,*)""
            !write(17,*)"Paramètres de simulation"
            !write(17,*)"-----------------------"
            !write(17,*)"N_sujet:   ",nsujet, "n_essai:",n_essai
            !write(17,*)"Theta(S):     ",theta2
            if(nsim_node(8)==2)then 
                !write(17,*)"Theta_T:     ",theta2_t
                !write(17,*)"rho_w_ST:     ",rsqrt_theta
            endif
            
            if(nsim_node(8).ne.0)then !model conjoint surrogate
                !write(17,*)"sigma_s:     ",sigma_s
                !write(17,*)"sigma_t:     ",sigma_t
                !write(17,*)"rho_v_st:     ",rsqrt
                if(frailt_base==1) then
                    !write(17,*)"gamma(S):        ",gamma_ui
                    ! if(indice_alpha==1) !write(17,*)"alpha:        ",alpha_ui
                endif
                if(nsim_node(8)==2)then 
                    !print*,"gamma_T:        ",gamma_ui
                    !print*,"rho_u_st:     ",rsqrt_gamma_ui
                endif
            endif
            !write(17,*)"lambdas:   ",lambdas
            !write(17,*)"nus:       ",nus
            !write(17,*)"lambdat:   ",lambdat
            !write(17,*)"nut:       ",nut
            !write(17,*)"beta_S:    ",betas
            !write(17,*)"beta_T:    ",betat
            !write(17,*)"censA:     ",temps_cens
            !write(17,*)"N_sim:     ",n_sim
            !write(17,*)""
            !write(17,*)"Paramètres d'estimation"
            !write(17,*)"-----------------------"
            if((nsim_node(4).ne.0) .and.(nsim_node(4).ne.3)) then
                if(nsim_node(4).eq.1) then
                    !write(17,*)"nbre_point:",nsim_node(2)
                endif
                if(nsim_node(4).eq.2) then
                    !write(17,*)"nbre_point:",nsim_node(2)
                    !write(17,*)"nbre bouble MC:",nsim_node(1)
                endif
            else
                !write(17,*)"nbre bouble MC:",nsim_node(1)
            endif
            !write(17,*)"kapa_use:  ",kapa_use
            !write(17,*)"noeuds utilise:",nz
            !!write(17,*)"valeurs initiales param de variances:",theta_init,sigma_ss_init,sigma_tt_init,sigma_st_init
            !write(17,*)"Nombre moyen d'itterarion pour la convergence:",moy_ni/n_sim_exact
            !write(17,*)""

            call cpu_time(tp2)
            !!write(17,*)'************temps_CPU en secondes******* =',tp2-tp1

            !write(17,*)"==============================================================="
            !write(17,*)"===================Moyenne empirique des estimations==========="
            !write(17,*)"==============================================================="

            !write(17,*)"parametres empiriques"
            !write(17,*)"*********************"
        endif
        if(gener_only.eq.1) then
            nbre_rejet=dble(0)
            !n_sim_exact=dble(n_sim-nbre_rejet)
            if(nb_processus<=1)then 
                n_sim_exact=dble(n_sim-nbre_rejet)
            else
                n_sim_exact=dble(n_sim_total-nbre_rejet)
            endif
            if(rang_proc==0) then
                !write(17,*)"moyenne empirique des personnes traitées=", moy_trt/n_sim_exact
                !write(17,*)"moyenne empirique des theta(S)_sim=",moy_theta/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types des theta_sim=",dsqrt(variance(tab_var_theta))
                if(nsim_node(8)==2)then
                    !write(17,*)"moyenne empirique des thetaT_sim= ",moy_thetat/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des thetaT_sim= ",se_theta_simt
                    !write(17,*)"moyenne empirique des rho_wij_ST_sim= ",moy_rho_wij/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des rho_wij_ST_sim= ",se_rho_sim_wij
                endif
            
                if(nsim_node(8).ne.0)then !model conjoint surrogate
                    !write(17,*)"moyenne empirique des sigmas_sim= ",moy_sigmas/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des sigmas_sim= ",dsqrt(variance(tab_var_sigma(:,1)))
                    !write(17,*)"moyenne empirique des sigmat_sim= ",moy_sigmat/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des sigmat_sim= ",dsqrt(variance(tab_var_sigma(:,2)))
                    !write(17,*)"moyenne empirique des rho_sim= ",moy_rho/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des rho_sim= ",dsqrt(variance(tab_var_sigma(:,3)))
                    if(frailt_base==1) then
                        !write(17,*)"moyemme empiriques des gamma(S)=",moy_gamma/n_sim_exact
                        !write(17,*)"ecart-type empirique des gamma(S)=",dsqrt(variance(tab_var_sigma(:,4)))
                        if(nsim_node(8)==2)then
                            !write(17,*)"moyenne empirique des gammaT_sim= ",moy_gammat/n_sim_exact
                            ! if(n_sim_exact>1) !write(17,*)"ecart-types des gammaT_sim= ",se_gamma_simt
                            !write(17,*)"moyenne empirique des rho_ui_ST_sim= ",moy_rho_ui/n_sim_exact
                            ! if(n_sim_exact>1) !write(17,*)"ecart-types des rho_ui_ST_sim= ",se_rho_sim_ui
                        endif
                    endif
                endif
                !write(17,*)"moyenne empirique des proportion progression=",moy_pros/n_sim_exact
                !write(17,*)"moyenne empirique des proportions décès=",moy_dec/n_sim_exact
            endif
            goto 1001
        endif

        if(nb_processus<=1) n_sim_exact=dble(n_sim-nbre_rejet) !nombre de simulations qui ont converges
        
        if(rang_proc==0) then
            !write(17,*)"moyenne empirique des personnes traitées=", moy_trt/n_sim_exact
            if(une_donnee==0)then
            !write(17,*)"moyenne empirique des theta(S)_sim=",moy_theta/n_sim_exact
            ! if(n_sim_exact>1) !write(17,*)"ecart-types des theta(S)_sim=",se_theta_sim
                if(nsim_node(8)==2)then
                    !write(17,*)"moyenne empirique des thetaT_sim= ",moy_thetat/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des thetaT_sim= ",se_theta_simt
                    !write(17,*)"moyenne empirique des rho_wij_ST_sim= ",moy_rho_wij/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des rho_wij_ST_sim= ",se_rho_sim_wij
                endif
                    
                if(nsim_node(8).ne.0)then !model conjoint surrogate
                    !write(17,*)"moyenne empirique des sigmas_sim= ",moy_sigmas/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des sigmas_sim= ",se_sigmas_sim
                    !write(17,*)"moyenne empirique des sigmat_sim= ",moy_sigmat/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des sigmat_sim= ",se_sigmat_sim
                    !write(17,*)"moyenne empirique des rho_sim= ",moy_rho/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types des rho_sim= ",se_rho_sim
                    if(frailt_base==1) then
                        !write(17,*)"moyemme empiriques des gamma(S)=",moy_gamma/n_sim_exact
                        !write(17,*)"ecart-type empirique des gamma(S)=",se_gamma_sim
                        if(nsim_node(8)==2)then
                            !write(17,*)"moyenne empirique des gammaT_sim= ",moy_gammat/n_sim_exact
                            ! if(n_sim_exact>1) !write(17,*)"ecart-types des gammaT_sim= ",se_gamma_simt
                            !write(17,*)"moyenne empirique des rho_ui_ST_sim= ",moy_rho_ui/n_sim_exact
                            ! if(n_sim_exact>1) !write(17,*)"ecart-types des rho_ui_ST_sim= ",se_rho_sim_ui
                        endif
                    endif
                endif
            endif
            !write(17,*)"moyenne empirique des proportion progression=",moy_pros/n_sim_exact
            !write(17,*)"moyenne empirique des proportions décès=",moy_dec/n_sim_exact
            !!write(17,*)"moyenne empirique des proportions censure=",moy_cens

            !write(17,*)""
            !write(17,*)"Parametres estimés"
            !write(17,*)"******************"
            !write(17,*)"theta(S)",theta2
            !write(17,*)"moyenne des theta(S) estimés=",moy_theta_est/n_sim_exact
            ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des theta(S) estimés=",se_theta_est
            !write(17,*)"moyenne des ecart-types des theta(S) estimés=",moy_se_theta/n_sim_exact
            !write(17,*)"taux de couverture des theta(S)=",taux_couverture_theta*100.d0/n_sim_exact
            
            if(nsim_node(8)==2)then
                !write(17,*)""
                !write(17,*)"theta_T",theta2_t
                !write(17,*)"moyenne des theta_T estimés=",moy_thetat_est/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des theta_T estimés=",se_thetat_est
                !write(17,*)"moyenne des ecart-types des theta_T estimés=",moy_se_thetat/n_sim_exact
                !write(17,*)"taux de couverture des theta_T=",taux_couverture_thetat*100.d0/n_sim_exact
                
                !write(17,*)""
                !write(17,*)"theta_ST",thetast_vrai
                !write(17,*)"moyenne des theta T estimés=",moy_thetast_est/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des theta_ST estimés=",se_cov_est_wij
                !write(17,*)"moyenne des ecart-types des theta_ST estimés=",moy_se_thetast/n_sim_exact
                !write(17,*)"taux de couverture des theta_ST=",taux_couverture_thetast*100.d0/n_sim_exact
            endif
            
            if(indice_eta==1)then
                !write(17,*)""
                !write(17,*)"eta",eta
                !write(17,*)"moyenne des eta estimés=",moy_eta/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des eta estimés=",se_eta_est
                !write(17,*)"moyenne des ecart-types des eta estimés=",moy_se_eta/n_sim_exact
                !write(17,*)"taux de couverture des eta=",taux_couverture_eta*100.d0/n_sim_exact
            endif

            if(nsim_node(8).ne.0)then !model conjoint surrogate
                if(frailt_base==1) then
                    !write(17,*)""
                    !write(17,*)"gamma(S)",gamma_ui
                    !write(17,*)"moyenne des gamma(S) estimés=",moy_gamma_est/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des gamma(S) estimés=",se_gamma_est
                    !write(17,*)"moyenne des ecart-types des gamma(S) estimés=",moy_se_gamma/n_sim_exact
                    !write(17,*)"taux de couverture des gamma(S)=",taux_couverture_gamma*100.d0/n_sim_exact
                    
                    if(indice_alpha==1)then
                        !write(17,*)""
                        !write(17,*)"Alpha",alpha_ui
                        !write(17,*)"moyenne des alpha estimés=",moy_alpha/n_sim_exact
                        ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des alpha estimés=",se_alpha_est
                        !write(17,*)"moyenne des ecart-types des alpha estimés=",moy_se_alpha/n_sim_exact
                        !write(17,*)"taux de couverture des alpha=    ",taux_couverture_alpha*100.d0/n_sim_exact
                    endif
                endif
                
                if(nsim_node(8)==2)then
                    !write(17,*)""
                    !write(17,*)"gamma_T",gamma_uit
                    !write(17,*)"moyenne des gamma_T estimés=",moy_gammat_est/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des gamma_T estimés=",se_gamma_estt
                    !write(17,*)"moyenne des ecart-types des gamma_T estimés=",moy_se_gammat/n_sim_exact
                    !write(17,*)"taux de couverture des gamma_T=",taux_couverture_gammat*100.d0/n_sim_exact
                    
                    !write(17,*)""
                    !write(17,*)"gamma_ST",gammast_vrai
                    !write(17,*)"moyenne des gamma_ST estimés=",moy_gammast_est/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des gamma_ST estimés=",se_cov_est_ui
                    !write(17,*)"moyenne des ecart-types des gamma_ST estimés=",moy_se_gammast/n_sim_exact
                    !write(17,*)"taux de couverture des gamma_ST=",taux_couverture_gammast*100.d0/n_sim_exact
                endif
                
                !write(17,*)""
                !write(17,*)"sigma_S",sigma_s
                !write(17,*)"moyenne des sigma_S estimés=",moy_sigmas_est/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des sigma_S estimés=",se_sigmas_est
                !write(17,*)"moyenne des ecart-types des sigma_S estimés=",moy_se_sigmas/n_sim_exact
                !write(17,*)"taux de couverture des sigma_S=",taux_couverture_sigmas*100.d0/n_sim_exact
                
                !write(17,*)""
                !write(17,*)"sigma_T",sigma_t
                !write(17,*)"moyenne des sigma_T estimés=",moy_sigmat_est/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des sigma_T estimés=",se_sigmat_est
                !write(17,*)"moyenne des ecart-types des sigma_T estimés=",moy_se_sigmat/n_sim_exact
                !write(17,*)"taux de couverture des sigma_T=",taux_couverture_sigmat*100.d0/n_sim_exact
                
                !write(17,*)""
                !write(17,*)"sigma_ST",sigmast_vrai
                !write(17,*)"moyenne des sigma_ST estimés=",moy_sigmast_est/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des sigma_ST estimés=",se_cov_est
                !write(17,*)"moyenne des ecart-types des sigma_ST estimés=",moy_se_sigmast/n_sim_exact
                !write(17,*)"taux de couverture des sigma_ST=",taux_couverture_sigmast*100.d0/n_sim_exact
            endif
        
            !write(17,*)""
            !write(17,*)"beta_S",betas
            !write(17,*)"moyenne des beta_S estimés=",moy_betaS/n_sim_exact
            ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des beta_S estimés=",se_beta_s
            !write(17,*)"moyenne des ecart-types des beta_S estimés=",moy_betaS_se/n_sim_exact
            !write(17,*)"taux de couverture des beta_S=",taux_couvertureS*100.d0/n_sim_exact

                !write(17,*)""
            !write(17,*)"beta_T",betat
            !write(17,*)"moyenne des beta_T estimés=",moy_betaT/n_sim_exact
            ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des beta_T estimés=",se_beta_t
            !write(17,*)"moyenne des ecart-types des beta_T estimés=",moy_betaT_se/n_sim_exact
            !write(17,*)"taux de couverture des beta_T=",taux_couvertureT*100.d0/n_sim_exact

            ! parametres de validation du surrogate
            if(nsim_node(8).ne.0)then
                !write(17,*)""
                !write(17,*)"R2_trial=",rsqrt**2
                !write(17,*)"moyenne des R2_trial estimés=",moy_R2_trial/n_sim_exact
                ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des R2_trial =",se_R2_trial
                !write(17,*)"moyenne des ecart-types (SE) des R2_trial estimés=",moy_se_R2_trial /n_sim_exact
                !write(17,*)"taux de couverture des R2_trial=",taux_couverture_R2_trial*100.d0/n_sim_exact
                !write(17,*)"moyenne des IC de R2_trial",moy_bi_R2_trial/n_sim_exact,";",moy_bs_R2_trial/n_sim_exact
                
                if(method_int_kendal==4 .or. method_int_kendal==5)then ! 1 seul taux de kendall et MC
                    !write(17,*)""
                    !write(17,*)"tau_kendal =",tau_kendal_00
                    !write(17,*)"moyenne des tau_kendal estimés=",moy_kendal_00/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des tau_kendal =",se_kendal_00
                    !write(17,*)"biais=",tau_kendal_00-moy_kendal_00/n_sim_exact
                else 
                    !write(17,*)""
                    !write(17,*)"tau_kendal_11=",tau_kendal_11
                    !write(17,*)"moyenne des tau_kendal_11 estimés=",moy_kendal_11/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des tau_kendal_11 =",se_kendal_11
                    
                    !write(17,*)""
                    !write(17,*)"tau_kendal_10=",tau_kendal_10
                    !write(17,*)"moyenne des tau_kendal_10 estimés=",moy_kendal_10/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des tau_kendal_10 =",se_kendal_10
                    
                    !write(17,*)""
                    !write(17,*)"tau_kendal_01=",tau_kendal_01
                    !write(17,*)"moyenne des tau_kendal_01 estimés=",moy_kendal_01/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des tau_kendal_01 =",se_kendal_01
                    
                    !write(17,*)""
                    !write(17,*)"tau_kendal_00=",tau_kendal_00
                    !write(17,*)"moyenne des tau_kendal_00 estimés=",moy_kendal_00/n_sim_exact
                    ! if(n_sim_exact>1) !write(17,*)"ecart-types empirique (SD) des tau_kendal_00 =",se_kendal_00
                endif
                
                ! taux de couverture par bootstrap
                if(rang_proc==0) then
                    !write(17,*)"=============Taux de couverture des parametres de validation par bootstrap===="
                    !write(17,*)"CP R2(en %):",CP_R2_boot*100.d0/n_sim
                    !write(17,*)"CP Kendall tau(en %):",CP_ktau_boot*100.d0/n_sim
                endif
            endif
            
            !write(17,*)""
            !write(17,*)"Kappa utilise surrogate et true=",ax1,ax2
            !write(17,*)"nbre_rejet",nbre_rejet
            !write(17,*)"Taux de rejet:",nbre_rejet*100.d0/(n_sim_exact+nbre_rejet)
            !close(17)
        endif
    endif
    
    1001 continue
    
    !desactivation de l'environnement de travail pour le programme parallele
    ! if(nsim_node(4).ne.3) then
        ! !call MPI_FINALIZE(comm) 
    ! endif
    !call MPI_FINALIZE(code) 
    !!print*,"suis là================="
    ! fermeture des fichiers des parametres
    !close(3)
    !close(11)
    !close(20)
    ! close(7)
    !close(8)
    ! close(9)
    ! close(16)
    ! close(12)
    !close(13)
    !!print*,"suis là=================1"
    !close(14)
    !!print*,"suis là=================2"
    ! close(15)
    !close(18)
    !!print*,"suis là=================3"
    !!print*,"suis là=================4"
    ! close(18)
    ! !print*,"suis là=================5"
    ! close(19)
    ! close(445)
    !!print*,"suis là=================7"   
   goto 998 
    if (istop.ne.1) then
        !write(*,*)"ERREUR : LE MODELE N'A PAS CONVERGE"
    else
        if (effet.eq.1) then
            !write(*,*)'** nombre d observations',nsujet
            !write(*,*)'** nombre de donnees recurrentes ',cpt
            !write(*,*)'** nombre de donnees censureees ',cptcens
            !write(*,*)'** nombre de deces ',cpt_dc
            !write(*,*)'** nombre de donnees tronquees gauche',cptaux
            !!write(4,*)'** nombre d observations',nsujet
            !!write(4,*)'** nombre de donnees recurrentes ',cpt
            !!write(4,*)'** nombre de donnees censureees ',cptcens
            !!write(4,*)'** nombre de deces ',cpt_dc
            !!write(4,*)'** nombre de donnees tronquees gauche',cptaux
        endif    

        !!write(4,*)'nombre de noeuds :',nz

        !!write(4,*)'** nombre de sujet SANS donnees recurrentes ',nb0recu
        !!write(4,*)'** nombre de donnees recurrentes ',moyrecu    
        moyrecu=moyrecu/ng
        !!write(4,*)'** nombre moyen de donnees recurrentes ',moyrecu
        !write(*,*)'** nombre de sujet SANS donnees recurrentes ',nb0recu
        !write(*,*)'** nombre de donnees recurrentes ',moyrecu
        moyrecu=moyrecu/ng
        !write(*,*)'** nombre moyen de donnees recurrentes ',moyrecu    
            
        !!write(4,*)'nombre total de parametres',np

        !!write(4,*)'kappa1',ax1
        !!write(4,*)'kappa2',ax2
        !!write(4,*)'nombre de noeuds:',nz

        
        !!write(4,*)'valeur de ni',ni 
        !write(*,*)'valeur de ni',ni
        
        
        if (effet.eq.1)then
            if(AG.eq.1)then
                !!write(4,*)'*************************************** ' 
                !!write(4,*)'**** ANDERSEN-GILL APPROACH *********** ' 
                !write(*,*)'**** ANDERSEN-GILL APPROACH *********** ' 
                !!write(4,*)'*************************************** ' 
            endif         
            !!write(4,*)'**************** ' 
            !!write(4,*)'THETA = variance de Z dans Z.exp(bX)'
            !!write(4,*)'(Z suit une GAMMA)'
            !!write(4,*)b(np-nva-indice_alpha)*b(np-nva-indice_alpha)
            
            !write(*,*)'**************** ' 
            
            !write(*,*)'THETA = ',b(np-nva-indice_alpha)*b(np-nva-indice_alpha)
    !     pour tenir compte du changt de var : delta methode
            !write(*,*)'SE theta (=H)',dsqrt(((2.d0*b(np-nva-indice_alpha))**2)* &
            ! H_hessOut(np-nva-indice_alpha,np-nva-indice_alpha))
            !!write(4,*)'SE theta(=H)',dsqrt(((2.d0*b(np-nva-indice_alpha))**2)* &
            ! H_hessOut(np-nva-indice_alpha,np-nva-indice_alpha))
            !write(*,*)'SE theta (=HIH)',dsqrt(((2.d0*b(np-nva-indice_alpha))**2)* &
            ! HIHOut(np-nva-indice_alpha,np-nva-indice_alpha))
            !!write(4,*)'SE theta(=HIH)',dsqrt(((2.d0*b(np-nva-indice_alpha))**2)* &
            ! HIHOut(np-nva-indice_alpha,np-nva-indice_alpha))
            
            !!write(4,*)'**************** ' 
            !write(*,*)'**************** '           
            !!write(4,*)'**************** '     

            !!write(4,*)'ALPHA dans z_i ^ alpha, pour le deces'
            !!write(4,*)b(np-nva)
            !write(*,*)'**************** ' 
            !write(*,*)'ALPHA = ',b(np-nva)
            
    !     pour tenir compte du changt de var : delta methode
        
            !write(*,*)'SE ALPHA (=H)',dsqrt(H_hessOut(np-nva,np-nva))
            !!write(4,*)'SE ALPHA (=H)',dsqrt(H_hessOut(np-nva,np-nva))
            !write(*,*)'SE alpha (=HIH)',dsqrt(HIHOut(np-nva,np-nva))
            !!write(4,*)'SE alpha (=HIH)',dsqrt(HIHOut(np-nva,np-nva))        
            !!write(4,*)'**************** ' 
            !write(*,*)'**************** ' 
            
        endif    
        
        if(nva.gt.0)then
        
            do i=1,nva
                j=(np-nva+i)*(np-nva+i+1)/2
                bi = b(np-nva+i) - 1.96*dsqrt(H_hessOut(np-nva+i,np-nva+i))
                bs = b(np-nva+i) + 1.96*dsqrt(H_hessOut(np-nva+i,np-nva+i))
                
                if(i.eq.1)then
                    !!write(4,*)'*** FOR RECURRENT EVENTS ***'
                endif
                
                if(i.eq.nva1+1)then
                    !!write(4,*)'*** FOR DEATH ***',nva,i
                endif
                
                !!write(4,*)'**************** '
                if(i.gt.nva1)then
                    !!write(4,*)'Variable : ',nomvar2(i-nva1)
                else
                    !!write(4,*)'Variable : ',nomvar(i)
                endif
                
                !!write(4,*)i,')','beta=',b(np-nva+i)
                !!write(4,*)' '
                !!write(4,*)i,')',' SE (=H)',dsqrt(H_hessOut(np-nva+i,np-nva+i))
                wres=(b(np-nva+i))/dsqrt(H_hessOut(np-nva+i,np-nva+i))
                !!write(4,*)'---> WALD',wres
                !!write(4,*)' '
                !!write(4,*)i,')',' SE (=HIH)',dsqrt(HIHOut(np-nva+i,np-nva+i))
                !!write(4,*)'---> WALD',(b(np-nva+i))/dsqrt(HIHOut(np-nva+i,np-nva+i))
                !!write(4,*)' '
                !!write(4,*)'RR : ',dexp(b(np-nva+i)),'  IC',dexp(bi),dexp(bs) 
                !!write(4,*)'**************** '        
            end do 
                
            !!write(4,*)'**************** ' 
            !!write(4,*)'---> log vraisemb marginale complete penalisee',resOut
            !!write(4,*)'**************** ' 
            !!write(4,*)'**************** ' 
            !!write(4,*)'---> The approximate likehood cross-valiadation criterion' 
            !!write(4,*)'---- in the semi parametric case LCV = ',LCV(1)
            !!write(4,*)'**************** '     

            if (typeof .ne.0) then
                !!write(4,*)'**************** ' 
                !!write(4,*)'---> AIC = 2k-2l =',LCV(2)
                !!write(4,*)'**************** ' 
            end if    
                
        endif
    endif

    ! if (nst == 2) then
        ! open(21,file=fich1) !hazard strate 1
        ! open(22,file=fich2) !surv strate 1
        ! open(23,file=fich3) !hazard strate 2
        ! open(24,file=fich4) !surv strate 2
        ! open(25,file=fich5) !Cumul hazard strate 1
        ! open(26,file=fich6) !Cumul hazard strate 2        
        
    ! else
        ! open(21,file=fich1)
        ! open(22,file=fich2)
        ! open(25,file=fich5) !Cumul hazard strate 1        
    ! end if
    
    ! select case(typeof)
    ! case(0)    
    
    ! do i=1,100
        ! !write(21,*)real(x1Out(i)),real(lamOut(i,1)),real(lamOut(i,2)),real(lamOut(i,3))
        ! !write(22,*)real(x1Out(i)),real(suOut(i,1)),real(suOut(i,2)),real(suOut(i,3))
        ! !write(25,*)real(x1Out(i)),real(-dlog(suOut(i,1))),real(-dlog(suOut(i,2))),real(-dlog(suOut(i,3)))    
    ! end do    
    
    ! if (nst.eq.2) then 
        ! do i=1,100
            ! !write(23,*)real(x2Out(i)),real(lam2Out(i,1)),real(lam2Out(i,2)),real(lam2Out(i,3))
            ! !write(24,*)real(x2Out(i)),real(su2Out(i,1)),real(su2Out(i,2)),real(su2Out(i,3))
            ! !write(26,*)real(x2Out(i)),real(-dlog(su2Out(i,1))),real(-dlog(su2Out(i,2))),real(-dlog(su2Out(i,3)))
        ! end do
    ! end if
!------------------------------------ Cpm    
    ! case(1)    
        ! do i=1,mt1
            ! !write(21,*)real(x1Out(i)),real(lamOut(i,1)),real(lamOut(i,2)),real(lamOut(i,3))
        ! end do
        ! do i=1,100
            ! !write(22,*)real(xSu1(i)),real(SuOut(i,1)),real(SuOut(i,2)),real(SuOut(i,3))
            ! !write(25,*)real(xSu1(i)),real(-dlog(suOut(i,1))),real(-dlog(suOut(i,2))),real(-dlog(suOut(i,3)))
        ! end do

    
    ! if (nst.eq.2) then 
        ! do i=1,mt2
            ! !write(23,*)real(x2Out(i)),real(lam2Out(i,1)),real(lam2Out(i,2)),real(lam2Out(i,3))
            
        ! end do
        ! do i=1,100
            ! !write(24,*)real(xSu2(i)),real(Su2Out(i,1)),real(Su2Out(i,2)),real(Su2Out(i,3))
            ! !write(26,*)real(xSu2(i)),real(-dlog(Su2Out(i,1))),real(-dlog(Su2Out(i,2))),real(-dlog(Su2Out(i,3)))
        ! end do    

    ! end if    
!------------------------------------ Weibull    
    ! case(2)    
        ! do i=1,100
            ! !write(21,*)real(x1Out(i)),real(lamOut(i,1)),real(lamOut(i,2)),real(lamOut(i,3))
            ! !write(22,*)real(xSu1(i)),real(suOut(i,1)),real(suOut(i,2)),real(suOut(i,3))
            ! !write(25,*)real(xSu1(i)),real(-dlog(suOut(i,1))),real(-dlog(suOut(i,2))),real(-dlog(suOut(i,3)))    
        ! end do        

    ! if (nst.eq.2) then 
        ! do i=1,100
            ! !write(23,*)real(x2Out(i)),real(lam2Out(i,1)),real(lam2Out(i,2)),real(lam2Out(i,3))
            ! !write(24,*)real(xSu2(i)),real(su2Out(i,1)),real(su2Out(i,2)),real(su2Out(i,3))
             ! !write(26,*)real(xSu2(i)),real(-dlog(su2Out(i,1))),real(-dlog(su2Out(i,2))),real(-dlog(su2Out(i,3)))
        ! end do
    ! end if        
    

    ! end select
            
    998 continue
    ! if (nst ==2) then
        ! close(21)
        ! close(22)
        ! close(23)
        ! close(24)
        ! close(25)
        ! close(26)
    ! else
        ! close(21)
        ! close(22)
        ! close(25)        
    ! end if

    call date_and_time(dateamj,heure2,zone,values)
    !!write(4,*) '***************************************************'   
    !!write(4,*) '*** starting time: ***', dateamj,heure1
    !!write(4,*) '*** Ending time:(hhmmss.sss) ***', dateamj,heure2        
    !close(4)
    ! close(2)    
    
    ! deallocate(groupe,ic,tt0,tt1,ttU,tt0dc,tt1dc,icdc,&
    ! !,suOut,su2Out,x1Out,x2Out,lamOut,lam2Out,ziOut,&
    ! filtre,filtre2,nomvart,nomvar2t, &
    ! vaxt,vaxdct,vax,nomvar,nomvar2,vaxdc)
    ! deallocate(MartinGales,linearpred,linearpreddc,time,timedc)    
    ! deallocate(nigs,cdcs,pourtrial,nigts,cdcts,nig_Ts,cdc_Ts)
    ! deallocate(don_simul,don_simulS1,p,prop_i,parametre_empirique,parametre_estimes,tab_rejet,parametre_empirique_NC)
    ! deallocate(kappa,tab_var_theta,donnee_essai,tableNsim,parametre_estimes_MPI,parametre_estimes_MPI_T)
    ! deallocate(vect_kendall_tau,v_chap_kendall,theta_chap_kendall,t_chap_kendall,v_chap_R2,theta_chap_R2,t_chap_R2,result_bootstrap)
    ! !deallocate(Vect_sim_MC)
    
    !generation des donnees par joint failty-copula 
    !call intpr("position test 1:", -1, ni, 1)
    if(nsim_node(11)==3) deallocate(don_simultamp,don_simulStamp)
    deallocate(moy_betaS, moy_betaT,moy_betaS_se, moy_betaT_se,taux_couvertureS, taux_couvertureT, &
               theta_chap_copula, v_chap_copula)
    deallocate(d_S,d_T,vbetas,vbetat)
    !call intpr("position test 2:", -1, ni, 1)
    endsubroutine jointsurrogate
    !complilation:
    !mpif90 -fopenmp -O3 -o exe_joint_surr_MPI_OMP  Adonnees.f90 Aparameters.f90 autres_fonctions.f90 Integrant_scl.f90 aaOptim_New_scl.f90 aaOptim_New_scl2.f90 funcpa_laplace.f90 aaOptim.f90 aaOptim_SCL_0.f90 aaOptimres.f90 funcpa_adaptative.f90 Integrale_mult_scl.f90 Pour_Adaptative.f90 aaUseFunction.f90 funcpajsplines_surrogate_scl_1.f90 funcpajsplines_surrogate_scl_2.f90 afuncpasres.f90 aresidusMartingale.f90 distance.f90 joint_surrogate.f90 main_Surr_simulation.f90
    !execution
    !time ./exe_joint_surr_MPI_OMP
    !time mpirun -n 1 ./exe_joint_surr_MPI_OMP