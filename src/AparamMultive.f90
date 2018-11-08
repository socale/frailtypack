    module taillesmultiv
    integer,save:: npmax,NSUJETMAX,nvarmax,nsujetmetamax
    integer,save:: ngmax                    !AD:,maxiter
    integer,save:: ndatemax,ndatemaxdc,nzmax,ndatemetamax
    integer,save::nssgbyg,nssgmax
    end module taillesmultiv

    module parametersmultiv
        double precision,save::epsa,epsb,epsd
        integer,save::maxiter
    end module parametersmultiv

    module sortiemultive
    integer,save::cptaux,cptauxmeta,cptcens,cptcensmeta,nb0recu
    double precision,save::moyrecu
    end module sortiemultive

    module comonmultiv
    implicit none
    integer,save::indiv
    double precision,dimension(:),allocatable,save::zero
!*****************************************************************
!*****dace1
    double precision,dimension(:),allocatable,save::date,datedc,datemeta
    double precision,dimension(:),allocatable,save::zi,zimeta,zidc

!*****dace2
    double precision,dimension(:),allocatable,save::t0dc,t1dc
    double precision,dimension(:),allocatable,save::t0,t1,t2,t3,t4, &
    t0meta,t1meta
    integer,dimension(:),allocatable,save:: c, cdc,cmeta
    integer,dimension(:),allocatable,save:: nt0,nt1
    integer,dimension(:),allocatable,save:: nt0meta,nt1meta
    integer,dimension(:),allocatable,save:: nt0dc,nt1dc
    integer,save::nsujet,nva,nva1,nva2,nva3,ndate,ndatedc,nst,ndatemeta,&
    nsujetmeta


!*****dace4
    integer,dimension(:),allocatable,save::stra
!*****ve1
    double precision,dimension(:,:),allocatable,save::ve,vemeta
    double precision,dimension(:,:),allocatable,save::vedc
!*****dace3
    double precision,save::pe
    integer,save::effet,nz1,nz2,nz3
    integer,save::nzloco,nzdc,nzmeta
!*****dace7
    double precision,dimension(:,:),allocatable,save::I_hess,H_hess
    double precision,dimension(:,:),allocatable,save::Hspl_hess!(npmax,npmax)
    double precision,dimension(:,:),allocatable,save::PEN_deri!(npmax,1)
    double precision,dimension(:,:),allocatable,save::hess!(npmax,npmax)
!*****contrib
    integer,save::ng          !nb de gpes
    integer,dimension(:),allocatable,save::vectn
!*****groupe
    integer,dimension(:),allocatable,save::g,gmeta!(nsujetmax)

    integer,dimension(:),allocatable,save::nig,nigmeta!(ngmax)  ! nb d events recurrents par sujet
!*****mem1
    double precision,dimension(:),allocatable,save::mm3,mm2!(ndatemax)
    double precision,dimension(:),allocatable,save::mm1,mm!(ndatemax)
!*****mem2
    double precision,dimension(:),allocatable,save::im3,im2,im1,im
!*****pen1
    double precision,dimension(:),allocatable,save::m3m3,m2m2,m1m1,mmm,m3m2
!*****pen2
    double precision,dimension(:),allocatable,save:: m3m1,m3m,m2m1,m2m,m1m,mi

    double precision, dimension(:),allocatable,save::mm3meta,mm2meta,mm1meta,mmmeta &
    ,im3meta,im2meta,im1meta,immeta

    double precision,dimension(:),allocatable,save::m3m3b,m2m2b,m1m1b,mmmb,m3m2b

    double precision,dimension(:),allocatable,save:: m3m1b,m3mb,m2m1b,m2mb,m1mb,mib

    double precision,dimension(:),allocatable,save::m3m3c,m2m2c,m1m1c,mmmc,m3m2c

    double precision,dimension(:),allocatable,save:: m3m1c,m3mc,m2m1c,m2mc,m1mc,mic

!************************************************************
!AD: add for death
    double precision, dimension(:),allocatable,save::mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc

!AD:end
!************************************************************
! %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%%
    integer,save::AG
! %%%%%%%%%%%%% indi! ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%
    integer,save::indic_ALPHA
!****  theta/alpha
    double precision,save::theta,alpha,eta,alpha1,alpha2 !en exposant pour la frailty deces
!****** indicateur de troncature
    integer,save:: indictronq,indictronqdc,indictronqmeta ! =0 si donnees non tronquées reellement
!*****auxig
    integer,save :: auxig
!******  aux1 aux2
    double precision,dimension(:),allocatable,save::res1,res3,res1meta,res3meta,aux1,aux2
    double precision,save::resnonpen
    double precision,dimension(:),allocatable::kkapa
!****** Type du modele
    integer,save::model
    double precision,dimension(:),allocatable::vvv
!cpm
    double precision ::cens
    integer,save:: nbrecu,nbdeces,nbrecumeta,nbintervR,nbintervDC,nbintervM


    integer,save::indic_eta,indic_rho,indic_a1,indic_a2
!    double precision,save::eta !en exposant pour la frailty deces
    double precision,dimension(:),allocatable,save::res4
    integer,save::typeof,typeof2
        double precision,dimension(:),allocatable,save::ttt,tttdc,tttmeta
        double precision,dimension(:),allocatable,save::betacoef
!Weib
    double precision,save::etaR,etaD,betaR,betaD,etaM,betaM
    integer,save::indic_tronc,indic_frail


    end module comonmultiv


    module residusMmultiv
!add multiv
        double precision,dimension(:),allocatable,save::ResidusRec2,Rrec2,Nrec2
!
        double precision,dimension(:),allocatable,save::Residus &
        ,varResidus,cumulhaz,vecuiRes,post_esp,post_SD,som_Xbeta
        double precision,dimension(:),allocatable,save::ResidusRec,&
        Residusdc,Rrec,Nrec,Rdc,Ndc,vecviRes,RisqCumul
        double precision,save::cares,cbres,ddres
        double precision,dimension(:),allocatable,save:: vres
        integer , save :: ierres,nires,istopres,effetres,indg
        double precision,save::rlres,varuiR,moyuiR,varviR,moyviR,corruiviR
        double precision,dimension(:),allocatable::vuu,b_temp
        integer,save::indic_cumul
        integer,dimension(:),allocatable::n_ssgbygrp
        double precision,dimension(:,:),allocatable,save::cumulhaz1,invsigma
        double precision,save::detSigma


    end module residusMmultiv
