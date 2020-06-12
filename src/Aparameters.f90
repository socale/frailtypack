    module tailles
    integer :: npmax,NSUJETMAX,nvarmax,nsujetymax  
    integer :: ngmax                    !AD:,maxiter
    integer :: ndatemax,ndatemaxdc,nzmax
    integer :: nssgbyg,nssgmax
    integer :: nboumax,NSIMAX,NOBSMAX
    integer :: nsujetBmax ! add TwoPart
    save
    end module tailles

    !module genparamsShort
    !integer , parameter :: nboumax=1000
    !integer , parameter :: NSIMAX=5000
    !integer , parameter :: NOBSMAX=15000
    !save
    !end module genparamsShort

    !module genparams
    !integer , parameter :: npmax=70
    !integer , parameter :: NSUJETMAX=15000
    !integer , parameter :: nvarmax=45
    !integer , parameter :: ngmax=5000
    !integer , parameter :: nboumax=1000
    !integer , parameter :: NSIMAX=5000
    !integer , parameter :: ndatemax=30000
    !integer , parameter :: NOBSMAX=15000
    !save
    !end module genparams

    module parameters
        double precision,save::epsa,epsb,epsd
        double precision,save::epsa_res,epsb_res,epsd_res
        integer,save::maxiter
    end module parameters

    module commun
    implicit none
    integer,save::ngexact,nssgexact
    integer,dimension(:,:),allocatable,save::ssg
    integer,dimension(:),allocatable,save:: mid
    integer,dimension(:,:),allocatable,save::mij,mij2
    integer,dimension(:),allocatable,save::mij_ind
    integer,save::nbpara
    double precision,dimension(:,:),allocatable,save::aux1,aux2
    end module commun

! time dependant janvier 2013
    module betatttps
    integer,save::npbetatps,npbetatps1,npbetatps2,npbetatps3,nbrecuTPS,nbrecumetaTPS,nbdecesTPS &
    ,entTPS,entmetaTPS,entdcTPS,npbetatpscross,nbinnerknots,qorder,equidistantTPS
        integer,dimension(:),allocatable,save::filtretps,filtre2tps,filtre3tps
    double precision,dimension(:),allocatable,save::knotsTPS,knotsdcTPS,knotsmetaTPS,the1TPS,the2TPS &
    ,betatps,betatps2,betatps3,betatpsX,theTPS,betatpsminX,betatpsmaxX,varBetatps,varBetatpsHIH &
    ,innerknots,innerknotsdc,innerknotsmeta,BasisSinhaTPS
    double precision,save::censtps,boundaryknots(2)
    end module betatttps

!=================================================
! pour le modele avec des donnees longitudinales

!AK
    module donnees_indiv
      implicit none
      double precision,dimension(:,:),allocatable,save::bb1,bb2,grandb
          double precision,save :: sigmav,range
      double precision,dimension(:,:),allocatable,save::mu,mu1
      double precision,dimension(:,:),allocatable,save::Z1,Z2
      integer,save::numpat,nmescur,nmescur2,nmescurr,nmescurr1,nmes,it_cur
      integer,parameter ::nf=1
      double precision,dimension(:),allocatable,save::b1
                double precision,dimension(:),allocatable,save :: ycurrent,ycurrent2,current_mean
                double precision,dimension(:),allocatable,save :: xeacurrent,part
                double precision,dimension(:,:),allocatable,save:: x2,x22,z22,z11,x2cur,z1cur
                integer,dimension(:),allocatable,save:: nii,nii2,nmes_o,nmes_o2
                integer,save:: it_rec
                double precision:: ut2cur,frailpol,frailpol2,frailpol3,frailpol4
                double precision,dimension(:),allocatable,save :: res1cur,res3cur,res2cur
    double precision,dimension(:,:),allocatable,save::Z1B, muB,XB,mu1B,x2Bcur,z1Bcur ! add TwoPart
    double precision,dimension(:),allocatable,save :: Bcurrent, current_meanRaw ! add TwoPart
    integer,save::nmescurB, it_curB, interceptBin !add TwoPart
    double precision,save::fixed_Binary
    integer, save::GLMloglink0,MTP0 ! glm log lionk + marginal two part                                                           
    end module donnees_indiv


        module random_effect
                double precision,dimension(:),allocatable,save::re

        end module random_effect

                module prediction
        double precision,dimension(2)::survDC,predtime_cm
        double precision::XbetapredDCi,XbetapredRi
        double precision,dimension(:),allocatable :: survRi, hazRi
        end module prediction

        module choix_epoce
                integer:: choix_e
        end module choix_epoce    
    
!=====================================================================================
    module comon
    implicit none
!*****************************************************************
    double precision,save :: K_G0, K_D0, lambda,y0
    double precision,dimension(:,:),allocatable:: nodes,weights
      integer,dimension(2),save::genz
      integer,save::npp,ni_cur, which_random
        double precision,save::vals
        integer ,save:: all, nnodes_all
        double precision,dimension(:),allocatable,save ::range
        integer,dimension(:),allocatable,save::RE_which
        !double precision,parameter::pi = 3.141592653589793238462643d0
        double precision,dimension(:),allocatable,save:: vi
        integer,dimension(:),allocatable,save:: nmesrec,nmesy,nmesrec1
        integer,save :: maxmesy, maxmesrec
        integer,dimension(:),allocatable,save:: groupee,groupeey
!*****dace1
    double precision,dimension(:),allocatable,save::date,datedc
    double precision,dimension(:),allocatable,save::zi,zidc

!*****dace2
    double precision,dimension(:),allocatable,save::t0dc,t1dc
    double precision,dimension(:),allocatable,save::t0,t1,t2,t3,yy,tU
    integer,dimension(:),allocatable,save:: c, cdc
    integer,dimension(:),allocatable,save:: nt0,nt1,ntU
    integer,dimension(:),allocatable,save:: nt0dc,nt1dc
    integer,save::nsujet,nva,nva1,nva2,nva3,nva4,ndate,ndatedc,nst,nstRec,nsujety,nobs,np_e
!*****dace4
    integer,dimension(:),allocatable,save::stra
!*****family
    integer,save::nfam
       integer,dimension(:),allocatable,save:: fam, fsize
!*****ve1
    double precision,dimension(:,:),allocatable,save::ve
    double precision,dimension(:,:),allocatable,save::vedc
    double precision,dimension(:,:),allocatable,save::vey
     !*** IJ: vector of weights
     double precision,dimension(:),allocatable,save::wtsvec

!*** donnees longitudinales
    double precision,save::vet3
!***** random effects
        double precision,dimension(:,:),allocatable,save :: ziy,ziyd,ziyr!random effects for y
         double precision,save :: sigmae                         !sigma of epsilon
        double precision,save :: etaydc1, etaydc2,etayr1,etayr2                 !reg coef for link functions
        double precision,dimension(:),allocatable,save:: etaydc,etayr              !reg coef for link functions
        integer,save :: nb_re,netar,netadc
        integer,save :: linkidyr,linkidyd,link
        double precision,dimension(:,:),allocatable,save::Ut,Utt,varcov_marg,sum_mat
         !****** censure a gauche
        double precision,save :: s_cag, box_cox_par
        integer,save :: s_cag_id, box_cox1
!*****dace3
    double precision,save::pe
    integer,save::effet,nz1,nz2,nzloco,nzdc
!*****dace7
    double precision,dimension(:,:),allocatable,save::I_hess,H_hess
    double precision,dimension(:,:),allocatable,save::Hspl_hess!(npmax,npmax)
    double precision,dimension(:,:),allocatable,save::PEN_deri!(npmax,1)
    double precision,dimension(:,:),allocatable,save::hess!(npmax,npmax)
    double precision,dimension(:,:),allocatable,save::H_hess_GH,I_hess_Gh,b_paGh
!*****contrib
    integer,save::ng          !nb de gpes
!*****groupe
    integer,dimension(:),allocatable,save::g!(nsujetmax)
    integer,dimension(:),allocatable,save::nig!(ngmax)  ! nb d events recurrents par sujet
!*****mem1
    double precision,dimension(:),allocatable,save::mm3,mm2!(ndatemax)
    double precision,dimension(:),allocatable,save::mm1,mm!(ndatemax)
!*****mem2
    double precision,dimension(:),allocatable,save::im3,im2,im1,im
!*****pen1
    double precision,dimension(:),allocatable,save::m3m3,m2m2,m1m1,mmm,m3m2
!*****pen2
    double precision,dimension(:),allocatable,save:: m3m1,m3m,m2m1,m2m,m1m,mi

!************************************************************
!AD: add for death
    double precision, dimension(:),allocatable,save::mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc
!AD:end
!************************************************************
! %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%%
    integer,save::AG
! %%%%%%%%%%%%% indic ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%
    integer,save::indic_ALPHA, indic_xi
!****  theta/alpha
    double precision,save::theta,alpha,eta,xi !en exposant pour la frailty deces
!****** indicateur de troncature
    integer,save:: indictronq,indictronqdc ! =0 si donnees non tronquees reellement
!*****auxig
    integer,save :: auxig, auxif
!******  aux1 aux2
    double precision,dimension(:),allocatable,save::res1,res3,res4,res5
    !double precision,dimension(:),allocatable,save::res1dc,res2dc,res3dc,res4dc,res5dc
    double precision,dimension(:),allocatable,save::aux1,aux2
    !double precision,dimension(:),allocatable,save::integrale4,integrale3gap
    double precision,save::resnonpen
    double precision,dimension(2)::kkapa
    double precision,dimension(:),allocatable,save::k0T
!****** Type du modele
    integer,save::model
    double precision,dimension(:),allocatable::vvv
!cpm
    double precision :: cens,mint ! rajout de mint
    integer,save:: nbrecu,nbdeces,nbintervR,nbintervDC
    integer,save::indic_eta
!double precision,save::eta !en exposant pour la frailty deces

    integer,save::typeof,typeof2
        double precision,dimension(:),allocatable,save::ttt,tttdc
        double precision,dimension(:),allocatable,save::betacoef
!Weib
    double precision,save::etaR,etaD,betaR,betaD
    integer,save::indic_tronc,typeJoint
    double precision,dimension(:),allocatable,save::etaT,betaT
!censure par intervalle
    integer,dimension(:),allocatable,save::d
    integer,save::dmax
    integer::intcens
    double precision,dimension(:),allocatable,save::resL,resU
! distribution des frailty par une log-normale
    integer::logNormal,timedep
    double precision,save::sig2, det
    
     double precision,dimension(:),allocatable,save::b_e
        integer,save::nea,nb1
        double precision,dimension(:),allocatable,save::timecur,timecur2
        double precision,dimension(:),allocatable,save::the1_e
        !parametres pour GH pseudo-adaptative
        double precision,dimension(:),allocatable,save::invBi_cholDet,vet22
        double precision,dimension(:,:),allocatable,save:: invBi_chol,b_lme,mat,matb_chol
        double precision,dimension(:),allocatable,save::v_jf,varv_jf
        integer::methodGH
        integer,save :: method_GH,nodes_number,initGH
        integer :: res_ind,it,n_wezly
        integer :: nb_gh,nb_gl
        !add current-level association - interaction with time
        integer,save :: numInter, numInterB
        integer,dimension(:),allocatable,save::positionVarT
        
    ! add TwoPart
        integer,save :: TwoPart, nsujetB, nbB,nby, maxmesB, nvaB
        double precision,dimension(:),allocatable,save::bb ! add TwoPart
        double precision,dimension(:,:),allocatable,save :: ziB,varcov_margB, sum_matB
    double precision,dimension(:,:),allocatable,save::veB
        integer,dimension(:),allocatable,save:: nmesB,nmes_oB,groupeeB !add TwoPart
        integer :: itB             
    end module comon
!=====================================================================================

    module comongroup
!=== add:18/04/2012
    integer,save::lignedc
    double precision,save::vet,vet2
    double precision,dimension(:),allocatable,save::the1,the2
    integer,dimension(:),allocatable,save::gsuj!attention gpe pour un sujet
    integer,dimension(:),allocatable,save::nigdc  ! nb de recurr ou dc par gpe
    integer,save::indic_joint
    double precision,save::expb1,expb2
    integer,dimension(:),allocatable,save::ictemp
    double precision,dimension(:),allocatable,save::temps0dc,temps1dc
    double precision,dimension(:,:),allocatable,save::variable
    double precision,dimension(:),allocatable,save::Binit
    double precision,dimension(:,:),allocatable,save::ve1,ve2,ve3
    double precision,dimension(:,:),allocatable,save::ve4 ! add TwoPart
    end module comongroup

    module jointmods
    !double precision,dimension(:),allocatable,save:: ut1,ut2,dut1,dut2
    !double precision,dimension(:),allocatable,save::b,bh
    double precision,dimension(:),allocatable,save::dut1,dut2
    double precision,dimension(:),allocatable,save::ut1,ut2
    double precision,dimension(:),allocatable,save::res1dc,res2dc,res3dc,res4dc,res5dc
    double precision,dimension(:),allocatable,save::res2,integrale4
    double precision,dimension(:),allocatable,save::integrale1,integrale2,integrale3,integrale3gap
    end module jointmods

    module additiv
    implicit none
        integer,save::correlini,correl
!*****contrib
        integer,save::ngexact       !nb EXACT de gpes
!*****mij
        integer,dimension(:),allocatable,save::mid ! nb de dc dans gpe i
!******indicateur du nb de parametres
        integer,save::nbpara
        double precision,dimension(:,:),allocatable,save::ve2
!*****inversion
        integer,save::indic_sousv,sousm
!*******sigma2tau2rho
        double precision,save::sigma2,tau2,rho,cov
!*****ut1ut2
        double precision,dimension(:),allocatable,save::dut1,dut2
        double precision,dimension(:),allocatable,save::ut1,ut2
!**** betaaux
        double precision,dimension(:),allocatable,save::betaaux
!*****invD
        double precision,dimension(:,:),allocatable,save::invD
        double precision,dimension(:),allocatable,save::aux1,aux2
        double precision,dimension(:,:),allocatable,save::Xbeta

    end module additiv


    module residusM
        double precision,dimension(:),allocatable,save::Residus, &
            varResidus,cumulhaz,vecuiRes,post_esp,post_SD,som_Xbeta
        double precision,dimension(:),allocatable,save::ResidusRec,&
            Residusdc,Rrec,Nrec,Rdc,Ndc,vecviRes,RisqCumul,ResidusLongi,&
            Rdc_res ,Nrec_fam,Ndc_fam,Nrec_ind
        double precision,save::cares,cbres,ddres
        double precision,dimension(:),allocatable,save:: vres
        integer , save :: ierres,nires,istopres,effetres,indg,it_res,it_res_rec
        integer, save :: it_resB ! add TwoPart
        double precision,save::rlres,varuiR,moyuiR,varviR,moyviR,corruiviR
        double precision,dimension(:),allocatable::vuu,b_temp
        integer,save::indic_cumul
        integer,dimension(:),allocatable::n_ssgbygrp
        double precision,dimension(:,:),allocatable,save::vecuiRes2,cumulhaz1,cumulhaz0,&
                cumulhazdc,invsigma,zet,zetd,zetr,XbetaY_res,Pred_y
        double precision,dimension(:,:),allocatable,save::XbetaB_res !add TwoPart
        double precision,save::detSigma
           integer,save :: nig_mc,np_mc
    double precision,save :: sig2_mc,res1_mc
      double precision,dimension(:),allocatable,save::mu1_res
     double precision,dimension(:,:),allocatable,save::ZetB ! add TwoPart
      double precision,dimension(:),allocatable,save::mu1_resB
    end module residusM

    module splines
    double precision,dimension(:),allocatable,save::aux,auxmeta
    double precision,dimension(:),allocatable,save:: v
    double precision,dimension(:,:),allocatable,save::I1_hess,H1_hess
    double precision,dimension(:,:),allocatable,save::I2_hess,H2_hess
    double precision,dimension(:,:),allocatable,save::HI2
    double precision,dimension(:,:),allocatable,save::HIH,IH,HI
    double precision,dimension(:,:),allocatable,save::BIAIS
    double precision,dimension(:),allocatable,save:: vax,vaxdc,vaxmeta,vaxy
    double precision,dimension(:),allocatable,save:: vaxB ! add TwoPart
    integer,dimension(:),allocatable,save::filtre,filtre2,filtre3, filtre4
    integer,dimension(:),allocatable,save::filtreB ! add TwoPart                                                  
    integer,save::ver

    end module splines

    module var_surrogate
    
    implicit none
    
! scl nouvelles variables pour surrogacy
        integer, save::Nmax, nsim,ntrials,posind_i,cpteu,position_i, affiche_itteration !Nmax= nombre total de sujet, nsim= nombre de simulation pour le montecarlo aumoins 10000,posind_i=position de l'individu courant pour l'integrant
        integer, save::methodInt ! methode d'integration,0= MC,1=MC+quadrature,2=quadrature guaussienne
        integer, dimension(:),allocatable,save ::nsujeti,pourtrial ! nombe de sujet par cluster, la taille vaut ng; indice essai des individu: la taille vaut nsujet
        integer,dimension(:),allocatable,save ::delta,deltastar ! Indicateur d'evenement surrogate, indocateur d'evenement true endpoint, la taille vaut Nmax
        integer,dimension(:),allocatable,save ::nigs,cdcs,nigts,cdcts! nombre de personnes avec progression et decede par essai, nigts,cdcts avec traitement
        integer, save:: indice_alpha,indice_sigma,indice_eta,indice_theta,indice_varS, indice_varT,indice_covST,&
                        indice_gamma,indice_alpha_ui,indice_gamma_t,indice_gamma_st,indice_theta_t,indice_theta_st! ces variables sont initialisees a 1 si le parametre correspondant est pris en compte dans le modele
        integer, save::type_mod_surr,type_joint ! model considere: 1= frailti individel et essais gaussiens, 0= individuel gamma et essai guassien
        double precision, save::sigma2,theta2,varS, varT,covST,alpha_ui,gamma_ui,theta2_t,theta_st,gamma_ui_t,gamma_ui_st !variances et covariances frailty trial specific(V_S,V_T)
        double precision,dimension(:),allocatable,save::const_res4,const_res5,res2s,res2_dcs,res2s_sujet,res2_dcs_sujet ! utiliser dans le calsul integrale, fourni dans funcpac,  la taille vaut nsujet. ,res2s,res2_dcs la taille vaut ntrial
        double precision,dimension(:),allocatable,save::const_res1,const_aux1 ! utiliser pour le calcul integrale dans funcpa
        double precision,dimension(2,2), save::varcov,varcovinv !matrice des variances covariance des fragilites (Vsi,Vti) et son inverse
        double precision,save::determinant ! contiendra le determinantde la matrice de la variance covariance pour le calcul de la vraisemblance
        !double precision,dimension(:,:),allocatable,save::ve ! matrice des variable explicatives avec le traitement comme premiere variable
        integer,save::ndim,npoint! ndim: dimension de l'integrale ou encore le nombre d'integration, npoint=nbre de point d'integration
        integer,save::zeta,nparamfrail! nombre de parametres associes aux effets aleatoires: eta+theta+varvs+varvt+covst=5
        logical, save::adaptative ! adaptative un boulien qui dit si on veux approximer par adaptative(.true.) ou pas (.false.)
        double precision,dimension(:),allocatable,save:: mu,mui,xx1,ww1 !mu= vecteur des mu pour monte carlo et mui=vecteur des bi chapeau pour l'adaptative, xx1,ww1 points et poids de quadrature
        double precision,dimension(:,:),allocatable,save:: vc,vcinv !vc= matrice de var-cov des enffet aleatoires MC, invBi_chol=inverse de la matrice de cholesky pour l'adaptative,vcinv=matrice inverse des vc
        double precision,dimension(:),allocatable,save:: invBi_chol_Essai,invBi_chol_Individuel! contient les element de la cholesky du determinant de la hessienne, niveau essai et individuel
        double precision,dimension(:,:),allocatable,save:: ui_chap,ui_chap_Essai, IhessLaplace,H_hess_laplace,hess_laplace ! ui_chap=matrice des effets aleatoires estimes pour l'adaptative, ui_chap_Essai pour les estimation au niveau essai
        double precision,dimension(:),allocatable,save:: invBi_cholDet_Essai,vvv_laplace,b_i_laplace,v_i_laplace ! invBi_cholDet racine carree du determinant de l'inverse de  la matrice de choloesky au niveau essa
        !double precision,dimension(:),allocatable,save:: invBi_cholDet ! invBi_cholDet racine carree du determinant de l'inverse de  la matrice de choloesky: utiliser la subroutine dmfsd pour le calcul de la matrice de cholesky, niveau individuel
        integer,save::a_deja_simul! dit si on a deja simule une fois les donnees pour le calcul integrale par MC
        integer,save::sujet_essai_max,vectorisation,individu_j ! taille du plus grand essai, vectorisation= indique si l'on fait de la vectorisation dans le calcul integral pour reduire les temps de calcul ou pas
        double precision,dimension(:,:),allocatable,save::Vect_sim_MC ! vecteur des nombre simules pour l'estimation de l'integrale par Monte carlo, wij_chap: pour contenir solution de k(w_ij) laplaca
        double precision,parameter::pi=3.141592653589793d0
        double precision,dimension(:,:),allocatable,save::Chol,wij_chap ! la cholesky de la matrice de variance covariance des effets aleatoires au niveau essai
        double precision,save::vs_i,vt_i,u_i,penalisation! pour l'adaptative
        integer,save::estim_wij_chap ! dit si l'on a deja estimer les wij pour la pseudo adaptative
        !integer,save::indice_B_essai ! compte le nombe d'element du vecteur invBi_chol_Essai des elements de la matrice B
        integer,save::Int_essai,control_param ! dit si l'on est entrain d'integrer sur les w_ij (0) ou alors sur les (vs_i,vt_i)(1)
        integer,save::essai_courant,control_adaptative ! donne la position du cluster courant, control_adaptative dit di on fait de l'adaptative (1) ou de la pseudo adaptative(0)
        integer, save::frailt_base    !dit si l'on prend en compte l'heterogeneite sur le risque de base aussi bien dans la generation des donnes que dans l'estimation(1) ou non (0)
        integer, save::indicej,comm2,nb_procs! utiliser dans la pseudo adaptative pour gerer la somme cumulee des sujet par essai
        integer,save::param_weibull,Test! parametrisation de la weibull utilisee: 0= parametrisation par defaut dans le programme de Virginie, 1= parametrisation a l'aide de la fonction de weibull donnee dans le cous de Pierre
        double precision, save::rho ! pour le calcul integrale par Laplace
        integer, save::graine,aleatoire,nbre_sim ! Pour la gestion du seed lors de la generation des nombres aleatoires
        integer,dimension(:),allocatable,save:: control_wij_chap,table_par_pro
        integer, save::switch_adaptative ! variable qui prend la valeur 0 si on veut ommetre la pseudo adaptative (cas typique lorsqu'echou l'estimation des effets aleatoires a posteriori) ou 1 si tout se passe bien
        integer, save::nbre_itter_PGH ! nombre d'itteration aubout desquelles reestimer les effects aleatoires a posteriori pour la pseude adaptative. si 0 pas de resestimation
        integer,save::random_generator ! generateur des nombre aleatoire, (1) si Random_number() et (2) si uniran(). Random_number() me permet de gerer le seed
        ! Add for the joint frailty-copula model -scl - 05-04-2019
        integer,save:: copula_function, control_affichage, control_adaptative_laplace ! the copula function, can be 1 for clayton or 2 for Gumbel-Hougaard
        double precision, save:: theta_copule ! copula parameters
    end module var_surrogate
    
    !gestion de la double precision
    module double_precision
    implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    end module double_precision