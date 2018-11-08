    module tailles
    integer :: npmax,NSUJETMAX,nvarmax,nsujetymax  
    integer :: ngmax                    !AD:,maxiter
    integer :: ndatemax,ndatemaxdc,nzmax
    integer :: nssgbyg,nssgmax
    integer :: nboumax,NSIMAX,NOBSMAX
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
         !****** censure à gauche
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
        integer :: methodGH,nodes_number,initGH
        integer :: res_ind,it,n_wezly
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
        double precision,save::rlres,varuiR,moyuiR,varviR,moyviR,corruiviR
        double precision,dimension(:),allocatable::vuu,b_temp
        integer,save::indic_cumul
        integer,dimension(:),allocatable::n_ssgbygrp
        double precision,dimension(:,:),allocatable,save::vecuiRes2,cumulhaz1,cumulhaz0,&
                cumulhazdc,invsigma,zet,zetd,zetr,XbetaY_res,Pred_y
        double precision,save::detSigma
           integer,save :: nig_mc,np_mc
    double precision,save :: sig2_mc,res1_mc
      double precision,dimension(:),allocatable,save::mu1_res
              

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
    integer,dimension(:),allocatable,save::filtre,filtre2,filtre3, filtre4
    integer,save::ver

    end module splines

