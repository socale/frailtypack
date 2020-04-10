


!========================          FUNCPAJ_SPLINES         ====================
    double precision function funcpajsplines_surrogate(b,np,id,thi,jd,thj,k0)
    
    use tailles
    use comon
    use residusM
    use var_surrogate ! fichier Aparameters_scl.f90
    use fonction_A_integrer ! integrant (fichier Integrant_scl.f90)
    use GaussHermi_mult ! pour la fonction d'integration (fichier Integrale_mult_scl.f90)
    use monteCarlosMult_Gaus ! pour integration par monte carlo (fichier GuassHermi_mult_scl.f90)
    use InverseMatrix !se trouve dans le fichier autres_fonctions.f90, pour l'inverse de la matrice des effets aleatoires et le calcul du determinant
    use func_adaptative
    use parameters, only: maxiter
    use optim_scl2, only:marq98j_scl2  ! pour faire appel a marquard
    use Autres_fonctions,only:Determinant_2,cholesky_factorisation
    !use func_laplace  ! pour tout ce qui est de l'approximation de laplace: fichier funcpa_laplace.f90
    use Laplace_contribution ! pour tout ce qui est de l'approximation de laplace: fichier Integrale_mult_scl.f90
    ! !$ use OMP_LIB
    !use mpi ! module pour l'environnement MPI
    use Autres_fonctions, only:init_random_seed
    
    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :
    
    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    
    integer::n,i,j,k,vj,ig,choix,l,vcdiag,nsujet_trial,dimint,dimint_Ind,init_i,max_i !erreur,code,compteur
    integer,dimension(ngmax)::cpt
    double precision::pe1,pe2,inv,som1,som2,res,vet,vet2,h1,varS1,varT1,covST1,som_cont_0 !inc,som
    double precision,dimension(3):: resultatInt
    
    double precision,dimension(-2:npmax):: the1,the2
    double precision,dimension(np)::bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc &
    ,res3dc,integrale1,integrale2
!AD: for death,change dimension 
    double precision,dimension(ndatemax)::dut1
    double precision,dimension(ndatemaxdc)::dut2
!AD:end
    double precision,dimension(0:ndatemax)::ut1
    double precision,dimension(0:ndatemaxdc)::ut2
    !double precision,dimension(:),allocatable::frail
    !double precision::int,logGammaJ,c3,c4,pourgam
    double precision,dimension(ntrials)::integrale3
    double precision,dimension(:,:),allocatable:: mat_A
    !double precision,dimension(ng,1)::wij_chap1
    
    ! var utilisees en vue de l'estimation des fragilites a posteriorie
    integer::ier,istop,ss,sss,ni,model_save,nparamfrail_save,maxiter_save,nmax_2,& !nb_pro2,comm,rang2,
             np_2,indice_B_essai,non_conv,rang,n_par_pro !indice_ind_util_essai,control_est
    integer::lm !frail_essai_deja_est !variable qui dit si pour un essai donne l'on a deja estimes les vsi et vti (1) ou non (0)
    integer,parameter::effet2=0,np_1=1
    double precision::ca,cb,dd,som_cont,usim,x22,SX
    !double precision::res
    double precision, dimension(2)::k0_2
    !double precision, dimension(1)::v,b_2
    double precision, dimension(:),allocatable::b_i      ! pour les 2 parametres des effets aleatoires a predire niveau essai
    double precision, dimension(:),allocatable::v_i    ! pour les 2 parametres des effets aleatoires a predire niveau essai
    double precision, allocatable, dimension(:,:)::H_hessOut,HIH,HIHOut,IH,invBi_chol_2,H_hess_scl,I_hess_scl
    double precision,dimension(:,:), allocatable::hess_scl
    double precision,dimension(:), allocatable::vvv_scl

!    !print*,'debut funcpa'
    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0
    ut1=0.d0
    ut2=0.d0
    dut2=0.d0
    dut1=0.d0
    varS1=0.d0
    varT1=0.d0
    covST1=0.d0 
    rang = 0    
    do i=1,np
        bh(i)=b(i)
    end do 

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj    
    n = (np-nva-nparamfrail)/nst
    !n = (np-nva-effet-indic_ALPHA)/nst
    !!print*,"nparamfrail=",nparamfrail,"effet=",effet,"indic_ALPHA=",indic_ALPHA
    !!print*,"(np-nva-nparamfrail)/nst",(np-nva-nparamfrail)/nst,"n=",n
    
    ! reparametrisation des parametres de la spline (>=0) pour être sur d'avoir une fonction des risquer de base positive
    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i 
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do
    
    !!print*,indice_eta+indice_theta+indice_varS+indice_varT+indice_covST
    !stop
    ! !print*,"effet=",effet
    if(effet.eq.1) then
        if(logNormal==1)then
            theta2 = bh(np-nva-nparamfrail+indice_eta+indice_theta)**2.d0 ! scl on recupere theta du vecteur des parametre, au carree car c'est bien la variance
            !!print*,"theta2=",theta2
            varS1 = bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS)
            varT1 = bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)
            sig2=theta2 ! je fais appel a sig2 car c'est la variable utilisee dans la suite des procedures pour le joint classique
            !!print*,"vs,vt,covst=",varS1,varT1,covST1
            if(frailt_base==1) then
                if(indice_alpha_ui==1)then
                    alpha_ui=bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha_ui)
                else
                    alpha_ui=1.d0
                endif
                gamma_ui=bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha_ui+&
                         indice_gamma)
            endif
        else
            !theta = bh(np-nva-nparamfrail+indice_eta+indice_theta)**2.d0
            !sig2=theta ! je fais appel a sig2 car c'est la variable utilisee dans la suite des procedures pour le joint classique
            !print*,"Loi gamma non disponible pour le modele complet de surrogacy; probleme de la listribution gamma bivariée"
        endif
        
        if(indice_eta==0)then
            eta=1.d0! on fixe eta a 1
        else
            eta = bh(np-nva-nparamfrail+indice_eta)
        endif
        !!print*,"theta2=",theta2
        if(indice_covST==1)then
            covST1 =bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST)
        else
            covST1=0.d0
        endif
        
        if(type_joint==2) then !on ajoute les parametre associes au modele complet avec effets aleatoires correles
            theta2=bh(np-nva-nparamfrail+indice_eta+indice_theta)
            theta2_t=bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        indice_theta_t)
            theta_st=bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        indice_theta_t+indice_theta_st)
            !!print*,"indice_theta_st=",indice_theta_st
            gamma_ui_t=bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        indice_theta_t+indice_theta_st+indice_gamma_t)
            if(indice_gamma_st==0) then !si on impose ne correleation nulle des frailties associes au risque de base
                gamma_ui_st=0.d0    !dans ce cas on ne l'estime pas
            else
                gamma_ui_st=bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                            indice_theta_t+indice_theta_st+indice_gamma_t+indice_gamma_st)
            endif
            ! !print*,"c=",bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        ! indice_theta_t),np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        ! indice_theta_t,"theta2_t=",theta2_t
            ! !print*,"indice=",np,nva,nparamfrail,indice_eta,indice_theta,indice_varS,indice_varT,indice_covST,indice_gamma,&
                        ! indice_theta_t
            ! !print*,"param=",theta2,theta2_t,theta_st,gamma_ui,gamma_ui_t,gamma_ui_st,varS1,covST1,varT1
            
            allocate(mat_A(6,6))
            Chol=0.d0 
            Chol(1,1)=theta2
            Chol(2,2)=theta2_t
            Chol(2,1)=theta_st
            Chol(3,3)=gamma_ui
            Chol(4,4)=gamma_ui_t
            Chol(4,3)=gamma_ui_st
            Chol(5,5)=varS1
            Chol(6,6)=covST1
            Chol(6,5)=varT1
            
            mat_A=MATMUL(Chol,TRANSPOSE(Chol))
            theta2=mat_A(1,1)
            theta2_t=mat_A(2,2)
            theta_st=mat_A(2,1)
            gamma_ui=mat_A(3,3)
            gamma_ui_t=mat_A(4,4)
            gamma_ui_st=mat_A(4,3)
            varS=mat_A(5,5)
            varT=mat_A(6,6)
            covST=mat_A(6,5)
            
            if((covST/(sqrt(varS)*sqrt(varT))>=1) .or. ((covST/(sqrt(varS)*sqrt(varT))<=-1)))then
                !print*, "la correlation sur vi est au dela de 1 ou < -1",covST/(sqrt(varS)*sqrt(varT))
                !stop
            endif
            
            if((gamma_ui_st/(sqrt(gamma_ui)*sqrt(gamma_ui_t))>=1) .or. (gamma_ui_st/(sqrt(gamma_ui)*sqrt(gamma_ui_t))>=1))then
                !print*, "la correlation sur ui est au dela de 1 ou < -1",gamma_ui_st/(sqrt(gamma_ui)*sqrt(gamma_ui_t))
                !stop
            endif
            
            if((theta_st/(sqrt(theta2)*sqrt(theta2_t))>=1) .or. (theta_st/(sqrt(theta2)*sqrt(theta2_t))>=1))then
                !print*, "la correlation sur wi est au dela de 1 ou < -1",theta_st/(sqrt(theta2)*sqrt(theta2_t))
                !stop
            endif
        endif
        ! !print*,"Chol=",Chol
        ! !print*,"mat_A=",mat_A
        
                ! varcov(1,1)=theta2
                ! varcov(1,2)=theta_st
                ! varcov(2,1)=theta_st
                ! varcov(2,2)=theta2_t
                
    endif
      ! !print*,"gamma_ui=",gamma_ui,"alpha_ui=",alpha_ui
      ! !print*,np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha_ui,&
            ! np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_alpha_ui+indice_gamma
            ! !print*,b(21),b(22)
      ! !print*,"size(b)=",size(b),np
      ! !print*,b(np-nva-nparamfrail+1:np)
      ! stop
    
    !!print*,np,nva,nparamfrail,indice_eta,indice_theta,indice_varS,indice_varT,indice_covST
    ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
    !Chol: matrice triangulaire inferieur. pour eviter de refaire la factorisation de cholesky pour l'algo MC, j'utilise directement cette matrice de cholesky a la place de la matrice de variance-covariance
    if(type_joint==1) then !cas modele a fragilites partages
        if(frailt_base==0)then
            allocate(mat_A(2,2))
            Chol=0.d0 
            Chol(1,1)=varS1
            Chol(2,1)=covST1
            !Chol(1,2)=covST1
            Chol(2,2)=varT1
        else
            allocate(mat_A(3,3))
            Chol=0.d0 
            Chol(1,1)=varS1
            Chol(2,1)=covST1
            !Chol(1,2)=covST1
            Chol(2,2)=varT1
            Chol(3,3)=gamma_ui
        endif
        mat_A=MATMUL(Chol,TRANSPOSE(Chol))
        varS=mat_A(1,1)
        varT=mat_A(2,2)
        covST=mat_A(1,2)
        if(frailt_base==1)    gamma_ui=mat_A(3,3)
    endif
    
    ! !print*,"chol",chol
    ! !print*,"mat_A",mat_A
    ! !print*,"gamma_ui",gamma_ui
    !stop
    !!print*,"chol vs,vt,covst=",varS1,varT1,covST1
    !!print*,"mat_A vs,vt,covst=",mat_A(1,1),mat_A(2,2),mat_A(1,2)
!    theta=theta**2 ! pour rester dans la logique de la cholesky
!!print*,"suis la dans funcpa============2"
!----------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!AD:modify
    dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
    dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
    
    ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
    ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
    
    ut1(0) = 0.d0
    ut2(0) = 0.d0
     
!//// NEW AMADOU vvv :
!--- strate1
    som1 = 0.d0
    vj = 0
    do i=2,ndate-1
        do k = 2,n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som1 = som1 + the1(j-4)
                vj  = j
                endif   
            endif
        end do 
    
        ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
        +(the1(j-1)*im1(i))+(the1(j)*im(i))
        
        dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
        +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

    end do
                   
!--- strate2
    vj = 0
    som2 = 0.d0

    do i=2,ndatedc-1
        do k = 2,n-2
            if (((datedc(i)).ge.(zi(k-1))).and.(datedc(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                som2 = som2 + the2(j-4)
                vj  = j
                endif   
            endif
        end do

        if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
            +(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))
            dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
            +(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))
        endif
            
    end do

!-------------fin strate2  
    i = n-2
    h1 = (zi(i)-zi(i-1))
    
    ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
    dut1(ndate) = (4.d0*the1(i-1)/h1)

        
    ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
    dut2(ndatedc) = (4.d0*the2(i-1)/h1)
         
!    !print*,ndatemaxdc
!    !print*,'ut2',ut2
    
    
!//// fin NEW AMADOU-vvv
!AD:end
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
        !integrale3_(k) = 0.d0
        aux1(k)=0.d0
        aux2(k)=0.d0
    end do
    
    integrale3 = 0.d0
    const_res1=0.d0
    const_aux1=0.d0
    const_res4=0.d0
    const_res5=0.d0
    res2s=0.d0
    res2_dcs=0.d0
    res2s_sujet=0.d0
    res2_dcs_sujet=0.d0

!*******************************************         
!-----avec un effet aleatoire dans le modele
!*********************************************

    inv = 1.d0/sigma2

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

    !!print*,"g=",g

    do i=1,nsujet 
        cpt(g(i))=cpt(g(i))+1  
        if(nva1.gt.0)then
            vet = 0.d0   
            do j=1,nva1
                vet =vet + bh(np-nva+j)*dble(ve(i,j))
!                  !write(*,*)'*** funcpaj_splines vet',vet,ve(i,j),i,j
            end do
        else
            vet=1.d0
        endif
         
    res2s_sujet(i)=dlog(dut1(nt1(i)))+vet
         
        if((c(i).eq.1))then
            res2s(pourtrial(i)) = res2s(pourtrial(i))+dlog(dut1(nt1(i)))+vet
        endif  
        if ((res2s(pourtrial(i)).ne.res2s(pourtrial(i))).or.(abs(res2s(pourtrial(i))).ge. 1.d30)) then
            funcpajsplines_surrogate=-1.d9
            goto 123
        end if    
        
        ! scl pour calcul integrale
        !const_res1(pourtrial(i))=const_res1(pourtrial(i))+ut1(nt1(i))*dexp(vet)
        !if ((const_res1(pourtrial(i)).ne.const_res1(pourtrial(i))).or.(abs(const_res1(pourtrial(i))).ge. 1.d30)) then
        !    funcpajsplines=-1.d9
        !    goto 123
        !end if
        
        const_res4(g(i)) = const_res4(g(i))+ut1(nt1(i))* dexp(vet)
        if ((const_res4(g(i)).ne.const_res4(g(i))).or.(abs(const_res4(g(i))).ge. 1.d30)) then !scl si risque cumule >10^30
            funcpajsplines_surrogate=-1.d9
            goto 123
        end if
    end do
!    !print*,'const_res1',const_res1
!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces 
!ccccccccccccccccccccccccccccccccccccccccc 

    do k=1,ng  
        if(nva2.gt.0)then
            vet2 = 0.d0   
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
            end do
            !vet2 = dexp(vet2)!scl
        else
            vet2=1.d0
        endif
        !!print*,"cdc=",cdc
        !stop
        
        res2_dcs_sujet(k)=dlog(dut2(nt1dc(k)))+vet2
        
        if(cdc(k).eq.1)then
            !!print*,"k=",k,"size(nt1dc)",size(nt1dc),"nt1dc(k)=",nt1dc(k),"dut2(nt1dc(k))=",dut2(nt1dc(k))
            !!print*,"dut2(1:10)",dut2(1:10)
            res2_dcs(pourtrial(k)) =res2_dcs(pourtrial(k))+dlog(dut2(nt1dc(k)))+vet2
            !!print*,"res2_dcs(pourtrial(k))=",res2_dcs(pourtrial(k))
            !stop
            if ((res2_dcs(pourtrial(k)).ne.res2_dcs(pourtrial(k))).or.(abs(res2_dcs(pourtrial(k))).ge. 1.d30)) then
                funcpajsplines_surrogate=-1.d9
                goto 123
            end if    
        endif 
      
        ! scl pour calcul integrale
        !const_aux1(pourtrial(k))=const_aux1(pourtrial(k))+ut2(nt1dc(k))*dexp(vet2)
        !if ((const_aux1(pourtrial(k)).ne.const_aux1(pourtrial(k))).or.(abs(const_aux1(pourtrial(k))).ge. 1.d30)) then
        !    funcpajsplines=-1.d9
        !    goto 123
        !end if
        const_res5(g(k)) = const_res5(g(k))+ut2(nt1dc(k))* dexp(vet2) 
!        !print*,"je suis la bien funcpajsplines_surr ligne 260, k=",k,"const_res5(k)=",const_res5(k)
        if ((const_res5(g(k)).ne.const_res5(g(k))).or.(abs(const_res5(g(k))).ge. 1.d30)) then
            funcpajsplines_surrogate=-1.d9
            goto 123
        end if
    end do
    
! !print*,const_res5(1),const_res4(1)
!**************INTEGRALES ****************************
!!print*,"ng=",ng,"ntrials=",ntrials

    !================================================================================
    !==========distribution lognormale des effects aleatoires==============================
    !================================================================================
!!print*,"suis la dans funcpa============3"    
    if (logNormal==1) then 
        select case(methodInt)
        case(3) ! estimation par Approximation de Laplace
            
            model_save=model
            nparamfrail_save=nparamfrail
            maxiter_save=maxiter
            model = 9 !scl pour le model effet aleatoires
            maxiter=20
            
            position_i=1
            !call cpu_time(c3)
            if(type_joint==1) then !cas modeles a effets aleatoires partages
                ! matrice des variances-covariance sur sigma_v
                varcov(1,1)=varS
                varcov(1,2)=covST
                varcov(2,1)=covST
                varcov(2,2)=varT
                call matinv(varcov,varcovinv,determinant) ! scl calcul de l'inverse et du determinant de la matrice de variance-covariance 
                rho=varcov(1,2)/dsqrt(varcov(1,1)*varcov(2,2))
                ! calcul des contribution des essais a la logvraisemblance
                som_cont=0.d0
                som_cont_0=0.d0
                ! !print*,"varcov",varcov
                
                ! !call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_pro2,code)
                 !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
                ! !print*,"nombre de processus pour comm2",nb_pro2, "mon rang",rang
                
                ! !call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
                ! !call MPI_COMM_RANK(MPI_COMM_WORLD,rang2,code)! pour chaque processus associe a l'identificateur code retourne son rang
                !!print*,"nombre de processus pour comm",nb_procs, "mon rang",rang
                
                ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
                !!call MPI_ABORT(comm2,erreur,code)
                n_par_pro=table_par_pro(rang+1) ! nombre de simulations a effectuer par le processus courant
                
                ! indice des calculs a effectuer par le processus courant
                if (rang==0) then
                    ! !print*, "table_par_pro=",table_par_pro
                    init_i=1 ! ce processus commence a la premiere simulation
                else
                    init_i=sum(table_par_pro(1:rang))+1 ! ce processus commence a la simulation qui respecte son ordre et doit s'arreter au nombre de simultation dont il a le droit d'executer
                endif
                
                max_i=init_i+table_par_pro(rang+1)-1!rang maximale de la simulation a executer (-1 car on a deja incrementer init_i de 1)
    
                ! !print*,"processus de rang",rang
                ! !print*,"table_par_pro",table_par_pro
                ! !print*,"rang=",rang,"min",init_i,"max=",max_i
                ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)! on stop tous les programmes appartenant au communicateur code, equivalent de l'instruction stop en sequantiel
                
                !compteur=compteur+1
                !!print*,"suis dans funcpa",size(wij_chap1),size(wij_chap1,1),size(wij_chap1,2)
                ! !$OMP PARALLEL DO default(none) firstprivate (model) PRIVATE (ig,resultatInt)& 
                ! !$OMP shared(ntrials,nsujeti,integrale3,determinant) REDUCTION(+:som_cont)
                    do ig=1,ntrials 
                        if((ig<init_i).or.ig>max_i) then 
                            goto 1000 ! pour dire le processus ne considere pas ce jeu de donnee
                        endif
                        
                        ! !print*,"trial",ig,"processus",rang
                        if(ig==1)then
                            position_i=1
                        else
                            position_i=sum(nsujeti(1:(ig-1)))+1
                        endif
                        ! !print*,"position_i",position_i
                        essai_courant=ig
                        resultatInt=0.d0
                        ! if(ig>0) then
                        ! if(ig==5) then
                            resultatInt=Cont_Laplace_Essai(determinant)                    
                        ! endif
                        ! som_cont=som_cont+resultatInt(1) 
                        ! !print*,"resultatInt(1)=",resultatInt(1)-(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant)&
                                ! +dlog(2.d0*pi*gamma_ui))
                        ! stop
                        if(resultatInt(1) .ne. -1.d9) then
                             som_cont=som_cont+resultatInt(1) + res2s(ig) &
                                + res2_dcs(ig)&
                                -(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant)&
                                +dlog(2.d0*pi*gamma_ui))
                        endif
                        integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                        !!print*,"contribution essai",ig,":",resultatInt
                        1000 continue
                    end do
                     ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
                ! !$OMP END PARALLEL DO
                ! !print*,"voila ma somme des contribution",som_cont,"processus",rang
                ! synthese des contributions
                ! !call MPI_REDUCE(som_cont,som_cont_0,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
                ! !print*,rang,"pour reduce",som_cont_0
                !call MPI_ALLREDUCE(som_cont,som_cont_0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
                ! !print*,rang,"som_cont avant",som_cont,som_cont_0
                ! call sleep(1)
                som_cont=som_cont_0
                ! !print*,rang,"som_cont aprest",som_cont
                ! !print*,"voila la somme totale des contribution",som_cont,"processus",rang
                !!call MPI_ABORT(MPI_COMM_WORLD,erreur,comm)
                
                ! !print*,"log-likelihood:",som_cont
                ! stop
                !!call MPI_Barrier(MPI_COMM_WORLD,comm2) ! pour la synchronisation globale avant 
            endif
            
            ! restitution des parametres
            model=model_save
            nparamfrail=nparamfrail_save
            maxiter=maxiter_save
            !desactivation de l'environnement de travail pour le programme parallele
            
            
        case(0) ! estimation par monte carlo
            posind_i=1
            !call cpu_time(c3)
            if(type_joint==1) then !cas modeles a effets aleatoires partages
                do ig=1,ntrials 
                    !auxig=ig pas utiliser
                    
                    allocate(mu(nsujeti(ig))) !initialisation du vecteur des moyennes des effects alatoires
                    allocate(vc(nsujeti(ig),nsujeti(ig))) !initialisation de la matrice de variance covariance pour le MC
                    !mu=0.d0
                    !vc=0.d0
                    !do l=1,nsujeti(ig)
                    !    vc(l,l)=theta2    
                    !end do        
                    vcdiag=0 ! la matrice n'est plus diagonale
                    resultatInt=0.d0
                    nsujet_trial=nsujeti(ig)
                    !calcul de l'integrale par monte carlo pour l'integrale multiple et quadrature adaptative ou pas pour l'integrale su vsi et vti
                    if(frailt_base==1)then
                        dimint=3 ! deux integrations au niveau essai correspondant aux effets aleatoires correles et un effet aleatoire associé au risque de base
                    else
                        dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles 
                    endif
                    resultatInt=MC_MultInd_Essai(Integrale_Individuel_MC,MC_Multiple_surr,dimint,nsujet_trial,ig,mat_A)
    !                !print*,"funcpa: resultatInt(ig)=",resultatInt(ig),"ig=",ig
                    !stop
                    
                    posind_i=posind_i+nsujeti(ig)
                    if(resultatInt(1).eq.0.d0) then
                        integrale3(ig)=0.1d-300
                        !!print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                    else
                        integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                    end if
                    !!print*,"funcpajsplines_surr ligne 312 nsujeti(ig)=",nsujeti(ig),"resultatInt=",integrale3(ig)  
                    deallocate(mu,vc)
                end do
            endif
            
            ! cas modele a effet aleatoires correles
            if(type_joint==2) then 
                        
                do ig=1,ntrials 
                    allocate(mu(nsujeti(ig))) !initialisation du vecteur des moyennes des effects alatoires
                    !allocate(vc(nsujeti(ig),nsujeti(ig))) !initialisation de la matrice de variance covariance pour le MC
                    
                    vcdiag=0 ! la matrice n'est plus diagonale
                    resultatInt=0.d0
                    nsujet_trial=nsujeti(ig)
                    
                    !calcul de l'integrale par monte carlo au niveau essai et quadrature gaussienne au niveau individuel
                    dimint_Ind=2
                    if(frailt_base==1)then
                        dimint=4 ! 4 integrations au niveau essai correspondant aux effets aleatoires correles et un effet aleatoire associé au risque de base
                    else
                        dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles 
                    endif
                    resultatInt=MC_MultInd_Essai_Cor(Integrale_Individuel_MC_cor,MC_Multiple_surr_cor,dimint_Ind,dimint,&
                                                     nsujet_trial,ig)
        !            resultatInt=MC_MultInd_Essai(Integrale_Individuel_MC,MC_Multiple_surr,dimint,nsujet_trial,ig,mat_A)
                     ! !print*,"funcpa: resultatInt(ig)=",dlog(resultatInt(ig)),"ig=",ig
                     ! stop
                    
                    posind_i=posind_i+nsujeti(ig)
                    if(resultatInt(1).eq.0.d0) then
                        integrale3(ig)=0.1d-300
                        !print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                    else
                        integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                    end if
                    !!print*,"funcpajsplines_surr ligne 534 nsujeti(ig)=",nsujeti(ig),"resultatInt=",integrale3(ig)  
                    !deallocate(mu,vc)
                    deallocate(mu)
                end do    
            endif
        case(2)! quadature classique (non-adaptative) ou pseudo-adaptative selon le contenu de la variable adaptative au niveau individuel et monte-carlo au niveau essai
            posind_i=1
            !cas modele avec effets aleatoires partages
            if(type_joint==1) then 
                if(rang==0)then
                    if(adaptative .and. control_adaptative==1) then ! on effectue le changement de variable                                    
                        ! !print*,"==============================================="
                        ! !print*,"Recherche des effets aleatoires  à postériorie"
                        ! !print*,"==============================================="
                    endif
                endif
                
                do ig=1,ntrials                 
                    allocate(mu(nsujeti(ig))) !initialisation du vecteur des moyennes des effects alatoires
                    allocate(vc(nsujeti(ig),nsujeti(ig))) !initialisation de la matrice de variance covariance pour le MC
                    
                    vcdiag=0 ! la matrice n'est plus diagonale
                    resultatInt=0.d0
                    nsujet_trial=nsujeti(ig)
                    
                    !calcul de l'integrale par monte carlo pour l'integrale multiple et quadrature adaptative ou pas pour l'integrale sur vsi et vti
                    if(frailt_base==1)then
                        dimint=3 ! deux integrations au niveau essai correspondant aux effets aleatoires correles et un effet aleatoire associé au risque de base
                    else
                        dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles 
                    endif
                    resultatInt=MC_Gauss_MultInd_Essai(Integrale_Individuel,gauss_HermMultA_surr,dimint,nsujet_trial,ig,npoint)
                    
                    deallocate(mu,vc)
                    
                    ! if(adaptative .and. resultatInt(1).eq.-1.d9) then
                        ! funcpajsplines_surrogate=-1.d9
                        ! goto 123
                    ! endif
                    
                    posind_i=posind_i+nsujeti(ig)
                    if(resultatInt(1).eq.0.d0) then
                        integrale3(ig)=0.1d-300
                        !!print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                    else
                        integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                    end if
                    !!print*,"funcpajsplines_surr ligne 312 nsujeti(ig)=",nsujeti(ig),"resultatInt=",integrale3(ig)  
                end do 
				!call intpr("control_affichage", -1, control_affichage, 1)
				! if (control_affichage < 2) then
					! call dblepr("valeur des integrales itteration 1", -1, integrale3, ntrials)
					! control_affichage = control_affichage + 1
				! endif
				
            endif
            
            !cas effets aleatoires correles
            if(type_joint==2) then 
                ! matrice des variances-covariance sur sigma_w
                varcov(1,1)=theta2
                varcov(1,2)=theta_st
                varcov(2,1)=theta_st
                varcov(2,2)=theta2_t
                call matinv(varcov,varcovinv,determinant) ! scl calcul de l'inverse de la matrice de variance-covariance et du determinant
                ! if(determinant.eq.0.d0) then ! mais ce cas n'est plus suppose arrive grace a la cholesky
                    ! !print*,"Attention determinant vaut 0"
                    ! !stop
                    ! determinant=0.d-10 ! ceci permet d'eviter les division par 0 si le determinant est =0
                ! end if
                
                do ig=1,ntrials                 
                    allocate(mu(nsujeti(ig))) !initialisation du vecteur des moyennes des effects alatoires
                    !allocate(vc(nsujeti(ig),nsujeti(ig))) !initialisation de la matrice de variance covariance pour le MC
                    
                    vcdiag=0 ! la matrice n'est plus diagonale
                    resultatInt=0.d0
                    nsujet_trial=nsujeti(ig)
                    
                    !calcul de l'integrale par monte carlo au niveau essai et quadrature gaussienne au niveau individuel
                    dimint_Ind=2
                    if(frailt_base==1)then
                        dimint=4 ! 4 integrations au niveau essai correspondant aux effets aleatoires correles et un effet aleatoire associé au risque de base
                    else
                        dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles 
                    endif
                    resultatInt=MC_Gauss_MultInd_Essai_Cor(Integrale_Individuel_cor,gauss_HermMultInd_cor,dimint_Ind,dimint,&
                                                           nsujet_trial,ig,npoint)
                    ! !print*,"funcpa: resultatInt(ig)=",dlog(resultatInt(ig))-nsujet_trial*(log(2*pi)+log(determinant)/2.d0)&
                            ! ,"ig=",ig
                    ! stop
                    
                    posind_i=posind_i+nsujeti(ig)
                    if(resultatInt(1).eq.0.d0) then
                        integrale3(ig)=0.1d-300
                        !!print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                    else
                        integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                    end if
                    !!print*,"funcpajsplines_surr ligne 534 nsujeti(ig)=",nsujeti(ig),"resultatInt=",integrale3(ig)  
                    !deallocate(mu,vc)
                    deallocate(mu)
                end do    
            endif
        case(1)! quadature classique (non-adaptative) ou pseudo-adaptative selon le contenu de la variable adaptative
!            !print*,"funcpajsplines_surrogate_scl.f90 lige 307: quadrature pas encore implémenté"
            
            !!print*, "suis la dans funcpa adaptative debut",adaptative,control_adaptative
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code) ! recherche du rang du processus
            l=1
            res = 0.d0
            ! matrice des variances-covariance
            varcov(1,1)=varS
            varcov(1,2)=covST
            varcov(2,1)=covST
            varcov(2,2)=varT
            call matinv(varcov,varcovinv,determinant) ! scl calcul de l'inverse de la matrice de variance-covariance et du determinant
            if(determinant.eq.0.d0) then ! mais ce cas n'est plus suppose arrive a grace a la cholesky
                !print*,"Attention determinant vaut 0"
                !stop
                determinant=0.d-10 ! ceci permet d'eviter les division par 0 si le determinant est =0
            end if
            
            
            !================================================================================
            !estimation des fragilites a posteriori, a utiliser dans le calcul integral
            !================================================================================
         
            if(adaptative .and. control_adaptative==1) then ! on effectue le changement de variable
                !if(estim_wij_chap.eq.0) then ! on n'a pas encore estime les wij_chap 
                    ! !print*,""
                    if(rang==0)then
                         ! !print*,"==============================================="
                         ! !print*,"Recherche des effets aleatoires  à postériorie"
                         ! !print*,"==============================================="
                    endif
                    ! !print*,""
                    k0_2=k0 
                                        
                    !initialisation des variables de module
                    individu_j=1
                    ni=0
                    ca=0.d0
                    cb=0.d0
                    dd=0.d0
                    model_save=model
                    nparamfrail_save=nparamfrail
                    maxiter_save=maxiter
                    model = 9 !scl pour le model effet aleatoires
                    maxiter=10
                    non_conv=0
                    ui_chap=0.d0
                    !posind_i_save=posind_i
                    
                    indice_B_essai=1 ! compte le nombe d'element du vecteur invBi_chol_Essai des elements de la matrice B
                    !indice_ind_util_essai=0 ! indice de l'individu utilise pour l'estimation des frailties dans l'essai
                    i=1
                    nmax_2=0 ! pour la somme cumulee du nombre de sujet par essai
                    posind_i=1 
                    do k=1,ntrials !k permet d'indicer les essais. ce changement empeche l'ambiguite avec le J passe en parametre
                        indicej=i
                        nmax_2=nmax_2+nsujeti(k)
                        essai_courant=k
                        ! ====================================================================================================
                        ! estimation des vs_i_chapeau et vt_i_chapeau
                        ! ====================================================================================================
                        !!print*
                        !!print*, "suis la dans funcpa adaptative"
                        
                        if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de l'heterogeneite sur les risque des bas
                            np_2=2
                            nparamfrail=2
                        else
                            np_2=3
                            nparamfrail=3
                        endif
                        !deallocate(H_hess_scl,I_hess_scl,H_hessOut,HIH,HIHOut,IH,invBi_chol_2,hess_scl,vvv_scl)
                        allocate(I_hess_scl(np_2,np_2),H_hess_scl(np_2,np_2),invBi_chol_2(np_2,np_2),H_hessOut(np_2,np_2),&
                                 b_i(np_2),v_i(np_2*(np_2+3)/2),HIH(np_2,np_2),HIHOut(np_2,np_2),IH(np_2,np_2),&
                                hess_scl(np_2,np_2),vvv_scl(np_2*(np_2+1)/2))
                        b_i=0.5d0
                        v_i=0.d0
                        
                        10 continue
                        call marq98J_scl2(k0_2,b_i,np_2,ni,v_i,res,ier,istop,effet2,ca,cb,dd,funcpafrailtyPred_Essai,&
                                         I_hess_scl,H_hess_scl,hess_scl,vvv_scl)
                                
                        ! !print*,"suis la2"
                        ! !print*,"b_i=",b_i
                        ! !print*,"H_hess_scl=",H_hess_scl
                        ! !print*,"I_hess_scl=",I_hess_scl
                        ! stop
                        
                        if (istop.ne.1 .and. non_conv<=10) then ! on passe à l'individu suivant, juste pour le test
                            b_i=-0.5*non_conv
                            non_conv=non_conv+1 !compte le nombre de fois qu'on n'a pas pu estime les frailties niveau essai sur certains individus
                            goto 10
                        endif
                                
                        if(non_conv==11 .and. istop .ne. 1)then
                            !print*,"le nombre de tentative sans convergence vaut:",non_conv
                            !print*,"istop=",istop,"essai k=",k
                            non_conv=0
                            funcpajsplines_surrogate=-1.d9
                            goto 123
                        endif
                                
                        if(non_conv>0 .and. non_conv<=10) then
                            !print*,"le nombre de tentative pour la convergence vaut:",non_conv
                            !print*,"istop=",istop,"essai k=",k
                            non_conv=0
                            !stop
                        endif
                        ! !print*,"k=",k,"istop=",istop
                        
                        ui_chap_Essai(k,1)=b_i(1) ! ui_chap_Essai: contient uniquement les vs_i_chapeau et vt_i_chapeau et u_i_chapeau
                        ui_chap_Essai(k,2)=b_i(2)
                        if(frailt_base==1) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de l'heterogeneite sur les risque des bas
                            ui_chap_Essai(k,3)=b_i(3)
                        endif
                        
                        ! !print*,"invBi_chol_2=",invBi_chol_2
                        ! !print*,"size(invBi_chol_2)",size(invBi_chol_2),"np_2=",np_2
                        do ss=1,np_2
                            do sss=1,np_2
                                !HIHOut(ss,sss) = HIH(ss,sss)
                                H_hessOut(ss,sss)= I_hess_scl(ss,sss)
                                invBi_chol_2(sss,ss)=H_hess_scl(sss,ss) ! je fais ceci juste pour le calcul de la cholesky I_hess_scl est que l'inverse de la hessienne
                                !!print*,invBi_chol_2(sss,ss)
                            end do
                        end do
                        ! !print*,"invBi_chol_2=",invBi_chol_2
                        ! !print*,"H_hessOut=",H_hessOut
                        call Cholesky_Factorisation(invBi_chol_2)! calcul de la cholesky de l'inverse de la hessienne 
                                
                        ! je sauvegarde ls element de la matrice B pour l'essai k
                        do ss=1,np_2     
                            do sss=1,np_2
                                invBi_chol_Essai(indice_B_essai)=invBi_chol_2(sss,ss)
                                indice_B_essai=indice_B_essai+1
                            enddo
                        enddo
                        ! je sauvegarde ls element de la matrice B pour l'individu i
                                
                        !calcul du determinant de la cholesky de l'inverse de la hessienne
                        ! invBi_chol_2(1,1)=1.d0
                        ! invBi_chol_2(2,1)=2.d0
                        ! invBi_chol_2(1,2)=3.d0
                        ! invBi_chol_2(2,2)=4.d0
                        ! !print*,"Determinant_2((/1,2,3,4/),np_2)",Determinant_2(invBi_chol_2,np_2)
                        ! stop
                        invBi_cholDet_Essai(k)=Determinant_2(invBi_chol_2,np_2) ! essai
                        !invBi_cholDet_Essai(k)=1.d0/dsqrt(Determinant_2(H_hess_scl,np_2)) ! essai
                        ! !print*,"invBi_cholDet_Essai(k)=",Determinant_2(invBi_chol_2,np_2)
                        ! ! on calcule plutôt le determinant de la hessienne inverse et on le passe sous la racine
                        ! !print*,"dsqrt(Determinant_2(I_hess_scl,np_2))",dsqrt(Determinant_2(I_hess_scl,np_2))
                        ! !print*,"1/dsqrt(Determinant_2(H_hess_scl,np_2))",1.d0/dsqrt(Determinant_2(H_hess_scl,np_2))
                        !stop
                        !invBi_cholDet(i)=invBi_chol_Individuel(i) !individuel    
                        ! !print*,"suis la1"        
                        deallocate(H_hessOut,HIH,HIHOut,IH,invBi_chol_2,I_hess_scl,H_hess_scl,hess_scl,vvv_scl,b_i,v_i)
                        ! allocate(I_hess_scl(np_1,np_1),H_hess_scl(np_1,np_1),invBi_chol_2(np_1,np_1),H_hessOut(np_1,np_1))
                        ! allocate(HIH(np_1,np_1),HIHOut(np_1,np_1),IH(np_1,np_1),hess_scl(np_1,np_1),vvv_scl(np_1*(np_1+1)/2))
                        
                        posind_i=posind_i+nsujeti(k) ! a utiliser dans funcpafrailtyPred_Essai
                        i=nmax_2+1 ! on continu avec le premier sujet du prochain cluster
                    enddo ! fin estimation des vs_i_chapeau et vt_i_chapeau
                    ! sauvegarde des resultats dans les fichiers
                    
                    ! open(20,file='Prediction_wij_chapeau.txt')
                    ! open(21,file='Prediction_vsi_vti_chapeau.txt')
                    
                    ! !!print*,"impression des frailties au individuel dans le fichier Prediction_wij_chapeau"
                    ! !write(20,*)"sujet"," ","w_ij_chapeau"
                    ! do ss=1,(i-1)
                        ! !write(20,*)ss,ui_chap(ss,1)
                    ! enddo            
                    
                    ! !!print*,"impression des frailties au niveau essai dans le fichier Prediction_vsi_vti_chapeau"
                    ! !write(21,*)"sujet"," ","vs_i_chapeau"," ","vt_i_chapeau"
                    ! do ss=1,ntrials
                        ! !write(21,*)ss,ui_chap_Essai(ss,1),ui_chap_Essai(ss,2)
                    ! enddo
                    
                    ! close(20)
                    ! close(21)
                    
                    !stop
                    ! deallocate(H_hessOut,HIH,HIHOut,IH,invBi_chol_2,hess_scl,vvv_scl)
                    estim_wij_chap=1 ! pour eviter de faire les estimations dans l'integrale    
                    ! deallocate(I_hess_scl,H_hess_scl) 
                    
                    model=model_save
                    nparamfrail=nparamfrail_save
                    maxiter=maxiter_save
                    !!print*,"suis la======================="
                    ! je remets les position i et jcol
                    individu_j=1
                    
                    ! sauvegarde des resultats dans les fichiers
                    
                    ! open(20,file='Prediction_wij_chapeau.txt')
                    ! open(21,file='Prediction_vsi_vti_chapeau.txt')
                    
                    !!print*,"impression des frailties au individuel dans le fichier Prediction_wij_chapeau"
                    !write(20,*)"sujet"," ","w_ij_chapeau"
                    do ss=1,(i-1)
                        !write(20,*)ss,ui_chap(ss,1)
                    enddo            
                    
                    !!print*,"impression des frailties au niveau essai dans le fichier Prediction_vsi_vti_chapeau"
                    !write(21,*)"sujet"," ","vs_i_chapeau"," ","vt_i_chapeau"
                    do ss=1,ntrials
                        !write(21,*)ss,ui_chap_Essai(ss,1),ui_chap_Essai(ss,2),ui_chap_Essai(ss,3)
                    enddo
                    
                    ! close(20)
                    ! close(21)
                    control_adaptative=0
                    !!print*,"suis la====="
                    !stop
                !endif
                if(rang==0)then
                    ! !print*,"================================================================================"
                    ! !print*,"Fin estimation des fragilites a posteriorie"
                    ! !print*,"================================================================================"
                endif
            endif
            !!print*,"suis sorti====","estim_wij_chap=",estim_wij_chap
            
    
            
            posind_i=1
            do ig=1,ntrials 
                nsujet_trial=nsujeti(ig)        
                resultatInt=0.d0
                if(frailt_base==0) dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles
                if(frailt_base==1) dimint=3 ! deux integrations au niveau essai correspondant aux effets aleatoires correles
                resultatInt=gauss_HermMultInd_Essai(Integrale_Individuel,gauss_HermMultA_surr,npoint,dimint,nsujet_trial,ig)
                !!print*,"12"
                !!print*,"funcpa resultatInt(1)=",-(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi)+dlog(theta2**nsujeti(ig)))+resultatInt(1)
                !stop
                !resultatInt=gauss_HermMult(integrant_indiv_1,npoint,nsujet_trial)
                !!print*,resultatInt(1),npoint,posind_i,nsujeti(ig)
                posind_i=posind_i+nsujeti(ig)
                !!print*,"funcpajsplines_surr ligne 285: resultatInt=",resultatInt(1),"ni=",nsujeti(ig),"essai",ig  
                if(resultatInt(1).eq.0.d0) then
                    integrale3(ig)=0.1d-300
                    !!print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                else
                    integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                endif 
                !deallocate(vc,frail,vcinv)
            end do
        !determinant=theta2**nsujeti(ig)
        !!print*,"cas quadrature"
        !!print*,-(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi)+dlog(determinant))+dlog(integrale3)
        !stop
        case(4)! quadature classique (non-adaptative) ou pseudo-adaptative (selon le contenu de la variable adaptative) niveau essai et monte-carlo niveau individuel
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code) ! recherche du rang du processus
            l=1
            res = 0.d0
            ! matrice des variances-covariance
            varcov(1,1)=varS
            varcov(1,2)=covST
            varcov(2,1)=covST
            varcov(2,2)=varT
            call matinv(varcov,varcovinv,determinant) ! scl calcul de l'inverse de la matrice de variance-covariance et du determinant
            if(determinant.eq.0.d0) then ! mais ce cas n'est plus suppose arrive a grace a la cholesky
                !print*,"Attention determinant vaut 0"
                !stop
                determinant=0.d-10 ! ceci permet d'eviter les division par 0 si le determinant est =0
            end if
            
            !===============================================================================
            ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
            !===============================================================================
            if(a_deja_simul.eq.0) then
                call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour lagraine
                Vect_sim_MC=0.d0 
                x22=0.d0
                lm=1
                do while(lm.le.nsim)
                    ! pour integrer sur un seul effet aleatoire au niveau individuel, on genere seulement des normales centree reduites, la transformation se fait dans le calcul de l'integran (fichier integrant.f90)
                    usim=0.d0
                    SX=1.d0
                    call bgos(SX,0,Vect_sim_MC(lm,1),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                    lm=lm+1
                end do                            
                a_deja_simul=1 ! pour dire qu'on ne simule plus
            endif
                        
            !================================================================================
            !estimation des fragilites a posteriori, a utiliser dans le calcul integral
            !================================================================================
         
            if(adaptative .and. control_adaptative==1) then ! on effectue le changement de variable
                !if(estim_wij_chap.eq.0) then ! on n'a pas encore estime les wij_chap 
                    ! !print*,""
                    if(rang==0)then
                         ! !print*,"==============================================="
                         ! !print*,"Recherche des effets aleatoires  à postériorie"
                         ! !print*,"==============================================="
                    endif
                    ! !print*,""
                    k0_2=k0 
                                        
                    !initialisation des variables de module
                    individu_j=1
                    ni=0
                    ca=0.d0
                    cb=0.d0
                    dd=0.d0
                    model_save=model
                    nparamfrail_save=nparamfrail
                    maxiter_save=maxiter
                    model = 9 !scl pour le model effet aleatoires
                    maxiter=10
                    non_conv=0
                    ui_chap=0.d0
                    !posind_i_save=posind_i
                    
                    indice_B_essai=1 ! compte le nombe d'element du vecteur invBi_chol_Essai des elements de la matrice B
                    !indice_ind_util_essai=0 ! indice de l'individu utilise pour l'estimation des frailties dans l'essai
                    i=1
                    nmax_2=0 ! pour la somme cumulee du nombre de sujet par essai
                    posind_i=1 
                    do k=1,ntrials !k permet d'indicer les essais. ce changement empeche l'ambiguite avec le J passe en parametre
                        indicej=i
                        nmax_2=nmax_2+nsujeti(k)
                        essai_courant=k
                        ! ====================================================================================================
                        ! estimation des vs_i_chapeau et vt_i_chapeau
                        ! ====================================================================================================
                        !!print*
                        !!print*, "suis la dans funcpa adaptative"
                        
                        if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de l'heterogeneite sur les risque des bas
                            np_2=2
                            nparamfrail=2
                        else
                            np_2=3
                            nparamfrail=3
                        endif
                        !deallocate(H_hess_scl,I_hess_scl,H_hessOut,HIH,HIHOut,IH,invBi_chol_2,hess_scl,vvv_scl)
                        allocate(I_hess_scl(np_2,np_2),H_hess_scl(np_2,np_2),invBi_chol_2(np_2,np_2),H_hessOut(np_2,np_2),&
                                 b_i(np_2),v_i(np_2*(np_2+3)/2),HIH(np_2,np_2),HIHOut(np_2,np_2),IH(np_2,np_2),&
                                hess_scl(np_2,np_2),vvv_scl(np_2*(np_2+1)/2))
                        b_i=0.5d0
                        v_i=0.d0
                        
                        11 continue
                        call marq98J_scl2(k0_2,b_i,np_2,ni,v_i,res,ier,istop,effet2,ca,cb,dd,funcpafrailtyPred_Essai,&
                                         I_hess_scl,H_hess_scl,hess_scl,vvv_scl)
                                
                        ! !print*,"suis la2"
                        ! !print*,"b_i=",b_i
                        ! !print*,"H_hess_scl=",H_hess_scl
                        ! !print*,"I_hess_scl=",I_hess_scl
                        ! stop
                        
                        if (istop.ne.1 .and. non_conv<=10) then ! on passe à l'individu suivant, juste pour le test
                            b_i=-0.5*non_conv
                            non_conv=non_conv+1 !compte le nombre de fois qu'on n'a pas pu estime les frailties niveau essai sur certains individus
                            goto 11
                        endif
                                
                        if(non_conv==11 .and. istop .ne. 1)then
                            !print*,"le nombre de tentative sans convergence vaut:",non_conv
                            !print*,"istop=",istop,"essai k=",k
                            non_conv=0
                            funcpajsplines_surrogate=-1.d9
                            goto 123
                        endif
                                
                        if(non_conv>0 .and. non_conv<=10) then
                            !print*,"le nombre de tentative pour la convergence vaut:",non_conv
                            !print*,"istop=",istop,"essai k=",k
                            non_conv=0
                            !stop
                        endif
                        ! !print*,"k=",k,"istop=",istop
                        
                        ui_chap_Essai(k,1)=b_i(1) ! ui_chap_Essai: contient uniquement les vs_i_chapeau et vt_i_chapeau et u_i_chapeau
                        ui_chap_Essai(k,2)=b_i(2)
                        if(frailt_base==1) then! on ignore simplement le terme avec ui si on ne doit pas tenir compte de l'heterogeneite sur les risque des bas
                            ui_chap_Essai(k,3)=b_i(3)
                        endif
                        
                        ! !print*,"invBi_chol_2=",invBi_chol_2
                        ! !print*,"size(invBi_chol_2)",size(invBi_chol_2),"np_2=",np_2
                        do ss=1,np_2
                            do sss=1,np_2
                                !HIHOut(ss,sss) = HIH(ss,sss)
                                H_hessOut(ss,sss)= I_hess_scl(ss,sss)
                                invBi_chol_2(sss,ss)=H_hess_scl(sss,ss) ! je fais ceci juste pour le calcul de la cholesky I_hess_scl est que l'inverse de la hessienne
                                !!print*,invBi_chol_2(sss,ss)
                            end do
                        end do
                        ! !print*,"invBi_chol_2=",invBi_chol_2
                        ! !print*,"H_hessOut=",H_hessOut
                        call Cholesky_Factorisation(invBi_chol_2)! calcul de la cholesky de l'inverse de la hessienne 
                                
                        ! je sauvegarde ls element de la matrice B pour l'essai k
                        do ss=1,np_2     
                            do sss=1,np_2
                                invBi_chol_Essai(indice_B_essai)=invBi_chol_2(sss,ss)
                                indice_B_essai=indice_B_essai+1
                            enddo
                        enddo
                        ! je sauvegarde ls element de la matrice B pour l'individu i
                                
                        !calcul du determinant de la cholesky de l'inverse de la hessienne
                        ! invBi_chol_2(1,1)=1.d0
                        ! invBi_chol_2(2,1)=2.d0
                        ! invBi_chol_2(1,2)=3.d0
                        ! invBi_chol_2(2,2)=4.d0
                        ! !print*,"Determinant_2((/1,2,3,4/),np_2)",Determinant_2(invBi_chol_2,np_2)
                        ! stop
                        invBi_cholDet_Essai(k)=Determinant_2(invBi_chol_2,np_2) ! essai
                        !invBi_cholDet_Essai(k)=1.d0/dsqrt(Determinant_2(H_hess_scl,np_2)) ! essai
                        ! !print*,"invBi_cholDet_Essai(k)=",Determinant_2(invBi_chol_2,np_2)
                        ! ! on calcule plutôt le determinant de la hessienne inverse et on le passe sous la racine
                        ! !print*,"dsqrt(Determinant_2(I_hess_scl,np_2))",dsqrt(Determinant_2(I_hess_scl,np_2))
                        ! !print*,"1/dsqrt(Determinant_2(H_hess_scl,np_2))",1.d0/dsqrt(Determinant_2(H_hess_scl,np_2))
                        !stop
                        !invBi_cholDet(i)=invBi_chol_Individuel(i) !individuel    
                        ! !print*,"suis la1"        
                        deallocate(H_hessOut,HIH,HIHOut,IH,invBi_chol_2,I_hess_scl,H_hess_scl,hess_scl,vvv_scl,b_i,v_i)
                        ! allocate(I_hess_scl(np_1,np_1),H_hess_scl(np_1,np_1),invBi_chol_2(np_1,np_1),H_hessOut(np_1,np_1))
                        ! allocate(HIH(np_1,np_1),HIHOut(np_1,np_1),IH(np_1,np_1),hess_scl(np_1,np_1),vvv_scl(np_1*(np_1+1)/2))
                        
                        posind_i=posind_i+nsujeti(k) ! a utiliser dans funcpafrailtyPred_Essai
                        i=nmax_2+1 ! on continu avec le premier sujet du prochain cluster
                    enddo ! fin estimation des vs_i_chapeau et vt_i_chapeau
                    ! sauvegarde des resultats dans les fichiers
                    
                    ! open(20,file='Prediction_wij_chapeau.txt')
                    ! open(21,file='Prediction_vsi_vti_chapeau.txt')
                    
                    ! !!print*,"impression des frailties au individuel dans le fichier Prediction_wij_chapeau"
                    ! !write(20,*)"sujet"," ","w_ij_chapeau"
                    ! do ss=1,(i-1)
                        ! !write(20,*)ss,ui_chap(ss,1)
                    ! enddo            
                    
                    ! !!print*,"impression des frailties au niveau essai dans le fichier Prediction_vsi_vti_chapeau"
                    ! !write(21,*)"sujet"," ","vs_i_chapeau"," ","vt_i_chapeau"
                    ! do ss=1,ntrials
                        ! !write(21,*)ss,ui_chap_Essai(ss,1),ui_chap_Essai(ss,2)
                    ! enddo
                    
                    ! close(20)
                    ! close(21)
                    
                    !stop
                    ! deallocate(H_hessOut,HIH,HIHOut,IH,invBi_chol_2,hess_scl,vvv_scl)
                    estim_wij_chap=1 ! pour eviter de faire les estimations dans l'integrale    
                    ! deallocate(I_hess_scl,H_hess_scl) 
                    
                    model=model_save
                    nparamfrail=nparamfrail_save
                    maxiter=maxiter_save
                    !!print*,"suis la======================="
                    ! je remets les position i et jcol
                    individu_j=1
                    
                    ! sauvegarde des resultats dans les fichiers
                    
                    ! open(20,file='Prediction_wij_chapeau.txt')
                    ! if(rang==0)    open(21,file='Prediction_vsi_vti_chapeau.txt')
                    
                    !!print*,"impression des frailties au individuel dans le fichier Prediction_wij_chapeau"
                    ! !write(20,*)"sujet"," ","w_ij_chapeau"
                    ! do ss=1,(i-1)
                        ! !write(20,*)ss,ui_chap(ss,1)
                    ! enddo            
                    
                    !!print*,"impression des frailties au niveau essai dans le fichier Prediction_vsi_vti_chapeau"
                    if(rang==0)then
                        !write(21,*)"sujet"," ","vs_i_chapeau"," ","vt_i_chapeau"
                        do ss=1,ntrials
                            !write(21,*)ss,ui_chap_Essai(ss,1),ui_chap_Essai(ss,2),ui_chap_Essai(ss,3)
                        enddo
                    endif
                    ! close(20)
                   ! if(rang==0)close(21)
                    control_adaptative=0
                    !!print*,"suis la====="
                    !stop
                !endif
                if(rang==0)then
                    ! !print*,"================================================================================"
                    ! !print*,"Fin estimation des fragilites a posteriorie"
                    ! !print*,"================================================================================"
                endif
            endif
            !!print*,"suis sorti====","estim_wij_chap=",estim_wij_chap
            
    
            
            posind_i=1
            do ig=1,ntrials 
                nsujet_trial=nsujeti(ig)        
                resultatInt=0.d0
                if(frailt_base==0) dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles
                if(frailt_base==1) dimint=3 ! deux integrations au niveau essai correspondant aux effets aleatoires correles
                resultatInt=gauss_HermMultInd_Essai_MC(Integrale_Individuel_MC,gauss_HermMultA_surr_MC,npoint,dimint,&
                                                       nsujet_trial,ig)
                !!print*,"12" 
                !!print*,"funcpa resultatInt(1)=",-(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi)+dlog(theta2**nsujeti(ig)))+resultatInt(1)
                !stop
                !resultatInt=gauss_HermMult(integrant_indiv_1,npoint,nsujet_trial)
                !!print*,resultatInt(1),npoint,posind_i,nsujeti(ig)
                posind_i=posind_i+nsujeti(ig)
                !!print*,"funcpajsplines_surr ligne 285: resultatInt=",resultatInt(1),"ni=",nsujeti(ig),"essai",ig  
                if(resultatInt(1).eq.0.d0) then
                    integrale3(ig)=0.1d-300
                    !!print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                else
                    integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                endif 
                !deallocate(vc,frail,vcinv)
            end do
        !determinant=theta2**nsujeti(ig)
        !!print*,"cas quadrature"
        !!print*,-(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi)+dlog(determinant))+dlog(integrale3)
        !stop
        end select
!!print*,"suis la dans funcpa============4"
!!print*,"dlog(integrale3)",integrale3
    !!print*,"res=",res
    
    !************* FIN INTEGRALES **************************
                   
        res=0.d0
        select case(methodInt)
        case(3) ! estimation par approximation de laplace
            ! je teste si l'une des estimations dans le calcul intégrale n'a pas marché 
            do ig=1,ntrials
                if(integrale3(ig).eq.-1.d9) then
                    funcpajsplines_surrogate=-1.d9
                    Rrec = 0.d0
                    Nrec = 0.d0
                    Rdc = 0.d0
                    Ndc = 0.d0
                    goto 123
                endif
            end do
                
            if(sigma2.gt.(1.d-8)) then      
                res=som_cont
            else
                res=som_cont
            endif
            
            
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajsplines_surrogate=-1.d9
                goto 123
            end if
        case(0) ! estimation par monte carlo
            do k=1,ntrials!ng  
                if(cpt(k).gt.0)then
                    if(sigma2.gt.(1.d-8)) then      
                        res= res + res2s(k) &
                        + res2_dcs(k)&
                        !+ integrale3(k)
                        + dlog(integrale3(k))
                    else
                        res= res + res2s(k) &
                        + res2_dcs(k)  &
                        !+ integrale3(k)
                        + dlog(integrale3(k))
                        ! !print*,"dlog(integrale3(k))=",dlog(integrale3(k))
                        ! stop
                    endif
                    
                    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        funcpajsplines_surrogate=-1.d9
                        goto 123
                    end if    
                endif
    !            !print*,'k',k
                !!print*,'res',res2s(k) + res2_dcs(k),'dlog(integrale3(k))', dlog(integrale3(k))     
            end do
        case(2) ! estimation par monte carlo niveau essai et quadrature niveau individuel
            ! effets aleatoires partages
            if(type_joint==1) then 
                do k=1,ntrials!ng  
                    if(cpt(k).gt.0)then
                        !!print*,"nsujeti(k)=",k,nsujeti(k)
                        if(sigma2.gt.(1.d-8)) then      
                            res= res + res2s(k)-(1.d0/2.d0)*nsujeti(k)*dlog(2.d0*pi*theta2) &
                            + res2_dcs(k)&
                            !+ integrale3(k)
                            + dlog(integrale3(k))
                        else
                            res= res + res2s(k)-(1.d0/2.d0)*nsujeti(k)*dlog(2.d0*pi*theta2) &
                            + res2_dcs(k)  &
                            !+ integrale3(k)
                            + dlog(integrale3(k))
                            ! !print*,"dlog(integrale3(k))=",dlog(integrale3(k))&
                            ! -(1.d0/2.d0)*nsujeti(k)*dlog(2.d0*pi*theta2)
                            ! stop
                        endif
                        
                        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                            funcpajsplines_surrogate=-1.d9
                            goto 123
                        end if    
                    endif
        !            !print*,'k',k
                    !!print*,'res',res2s(k) + res2_dcs(k),'dlog(integrale3(k))', dlog(integrale3(k))     
                end do
            endif
            
            ! effets aleatoirres correles
            if(type_joint==2) then 
                do k=1,ntrials!ng  
                    if(cpt(k).gt.0)then
                        !!print*,"nsujeti(k)=",k,nsujeti(k)
                        if(sigma2.gt.(1.d-8)) then      
                            res= res + res2s(k)-nsujeti(k)*(dlog(2.d0*pi) + (1.d0/2.d0)*dlog(determinant)) &
                            + res2_dcs(k)&
                            !+ integrale3(k)
                            + dlog(integrale3(k))
                        else
                            res= res + res2s(k)-nsujeti(k)*(dlog(2.d0*pi) + (1.d0/2.d0)*dlog(determinant))  &
                            + res2_dcs(k)  &
                            !+ integrale3(k)
                            + dlog(integrale3(k))
                        endif
                        
                        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                            funcpajsplines_surrogate=-1.d9
                            goto 123
                        end if    
                    endif
        !            !print*,'k',k
                    !!print*,'res',res2s(k) + res2_dcs(k),'dlog(integrale3(k))', dlog(integrale3(k))     
                end do
            endif
            
        case(1) !quadrature
            !!print*,"sum(res2s)=",SUM(res2s(1:ntrials)),"size(res2s)=",size(res2s)
            !!print*,"sum(res2_dcs)=",SUM(res2_dcs(1:ntrials)),"size(res2_dcs)=",size(res2_dcs)
            do k=1,ntrials!ng  
                if(cpt(k).gt.0)then
                    if(frailt_base==1) then
                        if(sigma2.gt.(1.d-8)) then  
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant)&
                            !+ integrale3(k)
                            +dlog(2.d0*pi*gamma_ui))+ dlog(integrale3(k))
                            
                        else
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant)&
                            !+ integrale3(k)
                            +dlog(2.d0*pi*gamma_ui))+ dlog(integrale3(k))
                            
                            ! !print*,"dlog(integrale3(k))=",dlog(integrale3(k))&
                            ! -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant)&
                            ! !+ integrale3(k)
                            ! +dlog(2.d0*pi*gamma_ui))
                            ! stop
                        endif
                    else
                        if(sigma2.gt.(1.d-8)) then  
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant))&
                            !+ integrale3(k)
                            + dlog(integrale3(k))
                        else
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant))&
                            !+ integrale3(k)
                            + dlog(integrale3(k))
                        endif
                    endif
     
                    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        funcpajsplines_surrogate=-1.d9
                        goto 123
                    end if    
                endif
    !            !print*,'k',k
                !!print*,'res',res2s(k) + res2_dcs(k),'dlog(integrale3(k))', -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))&
                !                                        +dlog(determinant))+dlog(integrale3(k))     
        end do
    case(4) !quadrature essai et monte-carlo individuel
            !!print*,"sum(res2s)=",SUM(res2s(1:ntrials)),"size(res2s)=",size(res2s)
            !!print*,"sum(res2_dcs)=",SUM(res2_dcs(1:ntrials)),"size(res2_dcs)=",size(res2_dcs)
            do k=1,ntrials!ng  
                if(cpt(k).gt.0)then
                    if(frailt_base==1) then
                        if(sigma2.gt.(1.d-8)) then  
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*((2.d0*dlog(2.d0*pi))+dlog(determinant)&
                            +dlog(2.d0*pi*gamma_ui))+ dlog(integrale3(k))
                            
                        else
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*((2.d0*dlog(2.d0*pi))+dlog(determinant)&
                            +dlog(2.d0*pi*gamma_ui))+ dlog(integrale3(k))
                            
                            ! !print*,"dlog(integrale3(k))=",dlog(integrale3(k))&
                            ! -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant)&
                            ! !+ integrale3(k)
                            ! +dlog(2.d0*pi*gamma_ui))
                            ! stop
                        endif
                    else
                        if(sigma2.gt.(1.d-8)) then  
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*((2.d0*dlog(2.d0*pi))+dlog(determinant))&
                            + dlog(integrale3(k))
                        else
                            !determinant=theta2**nsujeti(k) ! produit des elements de la diagonale
                            res= res + res2s(k) &
                            + res2_dcs(k)&
                            -(1/2.d0)*((2.d0*dlog(2.d0*pi))+dlog(determinant))&
                            + dlog(integrale3(k))
                        endif
                    endif
     
                    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        funcpajsplines_surrogate=-1.d9
                        goto 123
                    end if    
                endif
    !            !print*,'k',k
                !!print*,'res',res2s(k) + res2_dcs(k),'dlog(integrale3(k))', -(1/2.d0)*(nsujeti(k)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))&
                !                                        +dlog(determinant))+dlog(integrale3(k))     
        end do
        
        end select
        
 !!print*,"integrale3=",integrale3        
 ! !print*,"suis la dans funcpa============5, res=",res  
 ! stop

    endif
!!print*,res2s,res2_dcs,sum(res2s)+sum(res2_dcs),-pourgam,integrale3(k)
! !print*,"res =",res
! stop

    !================================================================================
    !=========================calcul de la penalisation==============================
    !================================================================================

    pe1 = 0.d0
    pe2 = 0.d0
    
    do i=1,n-3
        pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
        *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
        the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
        m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
        the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
        m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
        *the1(i)*m1m(i))
        if(nst.eq.1)then
            pe2=0.d0
        else
            pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
            *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
            the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
            m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
            the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
            m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
            *the2(i)*m1m(i))
        endif
    end do
    
    pe = k0(1)*pe1 + k0(2)*pe2 
    resnonpen = res
    res = res - pe
    ! !print*,"funcpajsplines ligne 434, vraisemblance penalisee res=",res,"resnonpen",resnonpen
    !cpteu=cpteu+1
    !if(cpteu.eq.1000) then
        !stop
    !endif
    deallocate(mat_A)
    if ((res.ne.res).or.(abs(res).ge. 1.d30).or.(res .ge. 0.d0)) then
        funcpajsplines_surrogate=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123

    else
        funcpajsplines_surrogate = res 
        ! section encore a definir en fonction de la suite
        do k=1,ng
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
        end do
    end if
!Ad:
123     continue
 !!print*,"suis la dans funcpa============6, res=",res  
 !stop 
    return
    
    end function funcpajsplines_surrogate



