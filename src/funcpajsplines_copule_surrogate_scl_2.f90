


!========================          FUNCPAJ_SPLINES         ====================
    double precision function funcpajsplines_copule_surrogate(b,np,id,thi,jd,thj,k0)
    
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
    use optim_scl, only:marq98j_scl  ! pour faire appel a marquard
    use Autres_fonctions,only:Determinant_2,cholesky_factorisation
    !use func_laplace  ! pour tout ce qui est de l'approximation de laplace: fichier funcpa_laplace.f90
    use Laplace_contribution ! pour tout ce qui est de l'approximation de laplace: fichier Integrale_mult_scl.f90
    ! !$ use OMP_LIB
    !use mpi ! module pour l'environnement MPI
    use Autres_fonctions, only:init_random_seed
    use func_laplace, only: funcpaLaplace_copula
    
    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :
    
    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    
    integer::n,i,j,k,vj,ig,choix,l,vcdiag,nsujet_trial,dimint,dimint_Ind,erreur,code,compteur,init_i,max_i
    integer,dimension(ngmax)::cpt
    double precision::pe1,pe2,som,inv,som1,som2,res,vet,vet2,h1,inc,varS1,varT1,covST1,som_cont_0, v_si, v_ti, ui
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
    double precision,dimension(:),allocatable::frail
    double precision::int,gammaJ,c3,c4,pourgam
    double precision,dimension(ntrials)::integrale3
    double precision,dimension(:,:),allocatable:: mat_A
    !double precision,dimension(ng,1)::wij_chap1
    
    ! var utilisees en vue de l'estimation des fragilites a posteriorie
    integer::ier,istop,ss,sss,ni,model_save,nparamfrail_save,maxiter_save,nmax_2,nb_pro2,comm,rang2,&
             np_2,indice_B_essai,indice_ind_util_essai,non_conv,control_est,rang,n_par_pro
    integer::frail_essai_deja_est,lm ! variable qui dit si pour un essai donne l'on a deja estimes les vsi et vti (1) ou non (0)
    integer,parameter::effet2=0,np_1=1
    double precision::ca,cb,dd,som_cont,usim,x22,SX, jacobien, f_vi
    !double precision::res
    double precision, dimension(2)::k0_2
    !double precision, dimension(1)::v,b_2
    double precision, dimension(:),allocatable::b_i      ! pour les 2 parametres des effets aleatoires a predire niveau essai
    double precision, dimension(:),allocatable::v_i    ! pour les 2 parametres des effets aleatoires a predire niveau essai
    double precision, allocatable, dimension(:,:)::H_hessOut,HIH,HIHOut,IH,invBi_chol_2,H_hess_scl,I_hess_scl
    double precision,dimension(:,:), allocatable::hess_scl
    double precision,dimension(:), allocatable::vvv_scl
    double precision, dimension(:,:),allocatable::m1,m3  
    double precision, dimension(:,:),allocatable::m    

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
    do i=1,np
        bh(i)=b(i)
    end do 

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj   
    !call intpr("nva", -1, nva, 1)    
    n = (np-nva-nparamfrail)/nst
    ! reparametrisation des parametres de la spline (>=0) pour être sur d'avoir une fonction des risquer de base positive
    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i 
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do
    
    if(effet.eq.1) then
        if(logNormal==1)then
            !== 17/05/2019== introduction of some corrections on the value of theta, to avoid very high values==
            if(copula_function == 1) then 
                theta_copule = dexp(bh(np-nva)) ! clayton: exp transform
                !theta_copule = dexp(minval((/6.d0,bh(np-nva)/))) ! clayton: exp transform
            endif
            if(copula_function == 2)then 
                theta_copule = (bh(np-nva))**2.d0  ! Gumbel: choleschy transform
                !theta_copule = minval((/bh(np-nva),15.d0/))**2.d0
            endif
            !theta_copule = bh(np-nva) ! sans transformation
            !call dblepr("theta_copule = ", -1, theta_copule, 1)
            varS1 = bh(np-nva-nparamfrail+indice_varS)
            varT1 = bh(np-nva-nparamfrail+indice_varS+indice_varT)
            !sig2=theta2 ! je fais appel a sig2 car c'est la variable utilisee dans la suite des procedures pour le joint classique

            if(frailt_base==1) then
                if(indice_alpha_ui==1)then
                    alpha_ui=bh(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_alpha_ui)
                else
                    alpha_ui=1.d0
                endif
                gamma_ui=bh(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_alpha_ui+&
                         indice_gamma)
            endif
        endif
        
        if(indice_covST==1)then
            covST1 =bh(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST)
        else
            covST1=0.d0
        endif
        
        if(type_joint==2) then !on ajoute les parametre associes au modele complet avec effets aleatoires correles
            gamma_ui_t=bh(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        1 +indice_gamma_t)
            if(indice_gamma_st==0) then !si on impose ne correleation nulle des frailties associes au risque de base
                gamma_ui_st=0.d0    !dans ce cas on ne l'estime pas
            else
                gamma_ui_st=bh(np-nva-nparamfrail+indice_varS+indice_varT+indice_covST+indice_gamma+&
                            1 +indice_gamma_t+indice_gamma_st)
            endif
            ! !print*,"c=",bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        ! indice_theta_t),np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT+indice_covST+indice_gamma+&
                        ! indice_theta_t,"theta2_t=",theta2_t
            ! !print*,"indice=",np,nva,nparamfrail,indice_eta,indice_theta,indice_varS,indice_varT,indice_covST,indice_gamma,&
                        ! indice_theta_t
            ! !print*,"param=",theta2,theta2_t,theta_st,gamma_ui,gamma_ui_t,gamma_ui_st,varS1,covST1,varT1
            
            allocate(mat_A(6,6)) ! je laisse telqu'el pour ne pas commettre des erreurs en modifiant
            Chol=0.d0 
            Chol(3,3)=gamma_ui
            Chol(4,4)=gamma_ui_t
            Chol(4,3)=gamma_ui_st
            Chol(5,5)=varS1
            Chol(6,6)=covST1
            Chol(6,5)=varT1
            
            mat_A=MATMUL(Chol,TRANSPOSE(Chol))
            gamma_ui=mat_A(3,3)
            gamma_ui_t=mat_A(4,4)
            gamma_ui_st=mat_A(4,3)
            varS=mat_A(5,5)
            varT=mat_A(6,6)
            covST=mat_A(6,5)
           
        endif
                
    endif
    
    !!print*,np,nva,nparamfrail,indice_eta,indice_theta,indice_varS,indice_varT,indice_covST
    ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
    !Chol: matrice triangulaire inferieur. pour eviter de refaire la factorisation de cholesky pour l'algo MC, j'utilise directement cette matrice de cholesky a la place de la matrice de variance-covariance
    if(type_joint==1 .or. type_joint==3) then !cas modele a fragilites partages
        allocate(mat_A(3,3))
        if(frailt_base==0)then
            Chol=0.d0 
            Chol(1,1)=varS1
            Chol(2,1)=covST1
            Chol(2,2)=varT1
        else
            Chol=0.d0 
            Chol(1,1)=varS1
            Chol(2,1)=covST1
            Chol(2,2)=varT1
            Chol(3,3)=gamma_ui
        endif
        mat_A = MATMUL(Chol,TRANSPOSE(Chol))
        varS=mat_A(1,1)
        varT=mat_A(2,2)
        covST=mat_A(1,2)
        if(frailt_base==1)    gamma_ui = mat_A(3,3)
    
        ! if(control_affichage == 0) then
            ! call dblepr("bh = ", -1, bh(np-nva-nparamfrail +1 :np),nva+nparamfrail)
            ! call dblepr(" gamma_ui funcpa =", -1, gamma_ui, 1)
            ! call dblepr(" Chol =", -1, Chol, 9)
            ! call dblepr(" TRANSPOSE(Chol) =", -1, TRANSPOSE(Chol), 9)
            ! call dblepr(" mat_A =", -1, mat_A, size(mat_A,1)*size(mat_A,2))
            ! control_affichage = 1
        ! endif
    endif
    ! call intpr(" size(Chol, 2)=", -1, size(Chol, 2), 1)
    
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
!     pour le surrogate
!ccccccccccccccccccccccccccccccccccccccccc

    
    do i=1,nsujet 
        cpt(g(i))=cpt(g(i))+1  
        if(nva1.gt.0)then
            vet = 0.d0   
            do j=1,nva1
                vet =vet + bh(np-nva+j)*dble(ve(i,j))
            end do
        else
            vet=1.d0
        endif
        ! call intpr("nva1", -1, nva1, 1)
        ! call dblepr(" ve(i,j)=", -1, ve(i,:), 2)
         
    !res2s_sujet(i)=dlog(dut1(nt1(i)))+vet
    res2s_sujet(i)=dut1(nt1(i)) * dexp(vet) ! baseline hazard for subject i
         
        ! if((c(i).eq.1))then
            ! res2s(pourtrial(i)) = res2s(pourtrial(i))+dlog(dut1(nt1(i)))+vet
        ! endif  
        ! if ((res2s(pourtrial(i)).ne.res2s(pourtrial(i))).or.(abs(res2s(pourtrial(i))).ge. 1.d30)) then
            ! funcpajsplines_copule_surrogate=-1.d9
            ! goto 123
        ! end if    
        
        const_res4(g(i)) = const_res4(g(i))+ut1(nt1(i))* dexp(vet) ! cumulative baseline hazard for subject i
        if ((const_res4(g(i)).ne.const_res4(g(i))).or.(abs(const_res4(g(i))).ge. 1.d30)) then !scl si risque cumule >10^30
            funcpajsplines_copule_surrogate=-1.d9
            goto 123
        end if
    end do
    
!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces 
!ccccccccccccccccccccccccccccccccccccccccc 

    do k=1,ng  
        !call dblepr("suis danc funcpan vedc=", -1, dble(vedc(k,:)), size(vedc,2))
        if(nva2.gt.0)then
            vet2 = 0.d0   
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
            end do
        else
            vet2=1.d0
        endif
        
        !res2_dcs_sujet(k)=dlog(dut2(nt1dc(k)))+vet2
        res2_dcs_sujet(k)=dut2(nt1dc(k))* dexp(vet2)
        
        ! if(cdc(k).eq.1)then
            ! res2_dcs(pourtrial(k)) =res2_dcs(pourtrial(k))+dlog(dut2(nt1dc(k)))+vet2
            ! if ((res2_dcs(pourtrial(k)).ne.res2_dcs(pourtrial(k))).or.(abs(res2_dcs(pourtrial(k))).ge. 1.d30)) then
                ! funcpajsplines_copule_surrogate=-1.d9
                ! goto 123
            ! end if    
        ! endif 
      
        const_res5(g(k)) = const_res5(g(k))+ut2(nt1dc(k))* dexp(vet2) 
        if ((const_res5(g(k)).ne.const_res5(g(k))).or.(abs(const_res5(g(k))).ge. 1.d30)) then
            funcpajsplines_copule_surrogate=-1.d9
            goto 123
        end if
    end do
    ! call dblepr("const_res4=", -1, const_res4, nsujet)
    ! call dblepr("const_res5=", -1, const_res5, nsujet)
    ! call dblepr("res2s_sujet=", -1, res2s_sujet, nsujet)
    ! call dblepr("res2_dcs_sujet=", -1, res2_dcs_sujet, nsujet)
!**************INTEGRALES ****************************

    !================================================================================
    !==========distribution lognormale des effects aleatoires==============================
    !================================================================================
    ! call intpr(" dans methodInt=", -1, methodInt, 1)
     ! call intpr("nsujeti=", -1, nsujeti, ntrials)
    if (logNormal==1) then 
        select case(methodInt)
        case(0) ! estimation par monte carlo
            posind_i=1
            ! varcov(1,1)=varS
            ! varcov(1,2)=covST
            ! varcov(2,1)=covST
            ! varcov(2,2)=varT
            ! call matinv(varcov,varcovinv,determinant) ! scl calcul de l'inverse de la matrice de variance-covariance et du determinant
            !call cpu_time(c3)
            if(type_joint==3) then !cas modeles a effets aleatoires partages
                do ig=1,ntrials 
                    allocate(mu(nsujeti(ig))) !initialisation du vecteur des moyennes des effects alatoires
                    allocate(vc(nsujeti(ig),nsujeti(ig))) !initialisation de la matrice de variance covariance pour le MC      
                    vcdiag=0 ! la matrice n'est plus diagonale
                    resultatInt=0.d0
                    nsujet_trial=nsujeti(ig)
                    !calcul de l'integrale par monte carlo pour l'integrale multiple et quadrature adaptative ou pas pour l'integrale su vsi et vti
                    if(frailt_base==1)then
                        dimint=3 ! deux integrations au niveau essai correspondant aux effets aleatoires correles et un effet aleatoire associé au risque de base
                    else
                        dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles 
                    endif
                    resultatInt=MC_Copula_Essai(Integrant_Copula,dimint,nsujet_trial,ig)
                    
                    posind_i=posind_i+nsujeti(ig)
                    if(resultatInt(1).eq.0.d0) then
                        integrale3(ig)=0.1d-300
                    else
                        integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                    end if
                    deallocate(mu,vc)
                end do
            endif
            
            ! ========= End for now===============
            !call dblepr("integrale3=", -1, integrale3, ntrials)
            !call dblepr("log integrale3=", -1, dlog(integrale3), ntrials)
            ! call dblepr(" dans sum integrale3=", -1, sum(integrale3), 1)
            ! call dblepr(" dans log sum integrale3=", -1, dlog(sum(integrale3)), 1)
        
        case(1)! quadature classique (non-adaptative) ou pseudo-adaptative selon le contenu de la variable adaptative
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
                        !call dblepr("Recherche des effets aleatoires  à postériorie", -1, integrale3(1), 1)
                    endif
                    ! !print*,""
                    k0_2=k0 
                                        
                    !initialisation des variables de module
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
                    do k=1,ntrials
                        essai_courant=k
                        ! ====================================================================================================
                        ! estimation des ui_chapeau, vs_i_chapeau et vt_i_chapeau
                        ! ====================================================================================================
                        
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
                        call marq98J_scl2(k0_2,b_i,np_2,ni,v_i,res,ier,istop,effet2,ca,cb,dd,funcpafrailtyPred_copula,&
                                         I_hess_scl,H_hess_scl,hess_scl,vvv_scl)
                            
                        if (istop.ne.1 .and. non_conv<=10) then ! on passe à l'individu suivant, juste pour le test
                            b_i=-0.5*non_conv
                            non_conv=non_conv+1 !compte le nombre de fois qu'on n'a pas pu estime les frailties niveau essai sur certains individus
                            goto 10
                        endif
                                
                        if(non_conv==11 .and. istop .ne. 1)then
                            !print*,"le nombre de tentative sans convergence vaut:",non_conv
                            !print*,"istop=",istop,"essai k=",k
                            non_conv=0
                            funcpajsplines_copule_surrogate=-1.d9
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
                        invBi_cholDet_Essai(k)=Determinant_2(invBi_chol_2,np_2) ! essai     
                        deallocate(H_hessOut,HIH,HIHOut,IH,invBi_chol_2,I_hess_scl,H_hess_scl,hess_scl,vvv_scl,b_i,v_i)
                        posind_i=posind_i+nsujeti(k) ! a utiliser dans funcpafrailtyPred_Essai
                        i=nmax_2+1 ! on continu avec le premier sujet du prochain cluster
                    enddo ! fin estimation des vs_i_chapeau et vt_i_chapeau    
                    model=model_save
                    nparamfrail=nparamfrail_save
                    maxiter=maxiter_save                   
                    control_adaptative=0

                if(rang==0)then
                    !call dblepr("Fin estimation des fragilites a posteriorie", -1, integrale3(1), 1)
                endif
            endif

            posind_i=1
            do ig=1,ntrials 
                nsujet_trial=nsujeti(ig)        
                resultatInt=0.d0
                if(frailt_base==0) dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles
                if(frailt_base==1) dimint=3 ! deux integrations au niveau essai correspondant aux effets aleatoires correles
                resultatInt=gauss_Herm_copula_Int(Integrant_Copula,npoint,dimint,nsujet_trial,ig)
                posind_i=posind_i+nsujeti(ig)
                if(resultatInt(1).eq.0.d0) then
                    integrale3(ig)=0.1d-300
                    !!print*,"integrale nulle et affectation de la valeur 0.1d-300, trials:",ig
                else
                    integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                endif 
            end do
   
            ! call dblepr("integrale3=", -1, integrale3, ntrials)
            ! call dblepr("log integrale3=", -1, dlog(integrale3), ntrials)
            
            
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
        case(3) ! estimation par Approximation de Laplace
            res = 0.d0
            ! matrice des variances-covariance
            varcov(1,1)=varS
            varcov(1,2)=covST
            varcov(2,1)=covST
            varcov(2,2)=varT
            call matinv(varcov,varcovinv,determinant) ! scl calcul de l'inverse de la matrice de variance-covariance et du determinant
            if(determinant.eq.0.d0) then ! mais ce cas n'est plus suppose arrive a grace a la cholesky
                determinant=0.d-10 ! ceci permet d'eviter les division par 0 si le determinant est =0
            end if
            
            
            !================================================================================
            !estimation des fragilites a posteriori, a utiliser dans le calcul integral
            !================================================================================
            if(rang==0)then
                !call dblepr("Recherche des effets aleatoires  à postériorie", -1, integrale3(1), 1)
            endif
            k0_2=k0                     
            !initialisation des variables de module
            ni=0
            ca=0.d0
            cb=0.d0
            dd=0.d0
            model_save=model
            nparamfrail_save=nparamfrail
            maxiter_save=maxiter
            model = 10 !scl pour le model effet aleatoires
            maxiter=30
            non_conv=0
            !i=1
            !nmax_2=0 ! pour la somme cumulee du nombre de sujet par essai
            
            ! call intpr("je vais pour le calcul integral=", -1, posind_i, 1)
            if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de l'heterogeneite sur les risque des bas
                np_2=2
                nparamfrail=2
            else
                np_2=3
                nparamfrail=3
            endif
            if(control_adaptative_laplace == 0) then
                b_i_laplace = 0.5d0
                v_i_laplace = 0.d0
            endif
                
            posind_i=1 
            do k=1,ntrials
                essai_courant=k
                ! ====================================================================================================
                ! estimation des ui_chapeau, vs_i_chapeau et vt_i_chapeau
                ! ====================================================================================================
                
                if(control_adaptative_laplace == 0) then    ! ici on voudrait estimer les une seule fois les v_i             
                    100 continue
                    call marq98J_scl2(k0_2,b_i_laplace,np_2,ni,v_i_laplace,res,ier,istop,&
                        effet2,ca,cb,dd,funcpaLaplace_copula,IhessLaplace,H_hess_laplace,&
                        hess_laplace,vvv_laplace)
                    
                    ! if(control_affichage == 0) then
                        ! control_affichage = 1
                        ! call intpr("istop=", -1, istop, 1)        
                        ! call dblepr("b_i_laplace=", -1, b_i_laplace, np_2)    
                    ! endif
                    if (istop.ne.1 .and. non_conv<=10) then ! pas de convergence, on modifie la valeur initiale et recommence l'optimisation
                        b_i_laplace=-0.5*non_conv
                        non_conv=non_conv+1 !compte le nombre de fois qu'on n'a pas pu estime les frailties niveau essai sur certains individus
                        goto 100
                    endif
                                    
                    if(non_conv==11 .and. istop .ne. 1)then
                        !print*,"le nombre de tentative sans convergence vaut:",non_conv
                        !print*,"istop=",istop,"essai k=",k
                        non_conv=0
                        funcpajsplines_copule_surrogate=-1.d9
                        goto 123
                    endif
                                    
                    if(non_conv>0 .and. non_conv<=10) then
                        non_conv=0
                        ! il y'a eu concergence
                        if(adaptative)control_adaptative_laplace = 1
                    endif 
                endif    

                
                jacobien = Determinant_2(IhessLaplace,np_2) ! determinant de la hesienne
                v_si = b_i_laplace(1)
                v_ti = b_i_laplace(2)
                if(frailt_base==1) then
                    ui = b_i_laplace(3)
                else
                    ui = 0.d0
                endif
                
                ! if(control_affichage == 0) then
                    ! control_affichage = 1
                    ! call dblepr("jacobien=", -1, jacobien, 1)        
                    ! call dblepr("b_i_laplace=", -1, b_i_laplace, np_2)    
                ! endif
                
                ! allocate(m(1,1),m1(1,2),m3(1,2))
                ! m1(1,1)= v_si
                ! m1(1,2)= v_ti
                ! m3=MATMUL(m1,varcovinv)
                ! m=MATMUL(m3,TRANSPOSE(m1))
                ! f_vi = 1.d0/(2.d0 * pi *  dsqrt(2.d0 * pi * gamma_ui * determinant)) * &
                       ! dexp(- 1.d0/2.d0 * m(1,1) - 1.d0/2.d0 * ui**2.d0 / gamma_ui)
                
                ! deallocate(m,m1,m3) 
                ! integrale3(k) = f_vi * (2.d0 * pi)**(np_2/2.d0) * Integrant_Copula(v_si,v_ti,ui,essai_courant,nsujeti(essai_courant))*&
                                    ! jacobien**(-1.d0/2.d0)
                integrale3(k) = (2.d0 * pi)**(np_2/2.d0) * Integrant_Copula(v_si,v_ti,ui,essai_courant,nsujeti(essai_courant))*&
                                    jacobien**(-1.d0/2.d0)
                posind_i=posind_i+nsujeti(k)
               !i=nmax_2+1 ! on continu avec le premier sujet du prochain cluster
            enddo ! fin calcul integral    

            ! call dblepr("integrale3=", -1, integrale3, ntrials)
            
            model=model_save
            nparamfrail=nparamfrail_save
            maxiter=maxiter_save                   

            if(rang==0)then
                !call dblepr("Fin estimation des fragilites a posteriorie", -1, integrale3(1), 1)
            endif
 
                ! !call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_pro2,code)
                 !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
                
                ! !call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
                ! !call MPI_COMM_RANK(MPI_COMM_WORLD,rang2,code)! pour chaque processus associe a l'identificateur code retourne son rang
                
                ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
                !!call MPI_ABORT(comm2,erreur,code)
                ! n_par_pro=table_par_pro(rang+1) ! nombre de simulations a effectuer par le processus courant
                
                ! ! indice des calculs a effectuer par le processus courant
                ! if (rang==0) then
                    ! ! !print*, "table_par_pro=",table_par_pro
                    ! init_i=1 ! ce processus commence a la premiere simulation
                ! else
                    ! init_i=sum(table_par_pro(1:rang))+1 ! ce processus commence a la simulation qui respecte son ordre et doit s'arreter au nombre de simultation dont il a le droit d'executer
                ! endif
                
                ! max_i=init_i+table_par_pro(rang+1)-1!rang maximale de la simulation a executer (-1 car on a deja incrementer init_i de 1)
    
                ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)! on stop tous les programmes appartenant au communicateur code, equivalent de l'instruction stop en sequantiel
                
                !compteur=compteur+1
                ! !$OMP PARALLEL DO default(none) firstprivate (model) PRIVATE (ig,resultatInt)& 
                ! !$OMP shared(ntrials,nsujeti,integrale3,determinant) REDUCTION(+:som_cont)
                    ! do ig=1,ntrials 
                        ! if((ig<init_i).or.ig>max_i) then 
                            ! goto 1000 ! pour dire le processus ne considere pas ce jeu de donnee
                        ! endif
                        
                        ! if(ig==1)then
                            ! position_i=1
                        ! else
                            ! position_i=sum(nsujeti(1:(ig-1)))+1
                        ! endif

                        ! essai_courant=ig
                        ! resultatInt=0.d0
                        ! resultatInt=Cont_Laplace_Essai(determinant)                    
                        ! if(resultatInt(1) .ne. -1.d9) then
                             ! som_cont=som_cont+resultatInt(1) + res2s(ig) &
                                ! + res2_dcs(ig)&
                                ! -(1/2.d0)*(nsujeti(ig)*dlog(2.d0*pi*theta2)+(2.d0*dlog(2.d0*pi))+dlog(determinant)&
                                ! +dlog(2.d0*pi*gamma_ui))
                        ! endif
                        ! integrale3(ig) = resultatInt(1) ! on recupere le premier element car les deux autres sont supposes etre la precision et la variance 
                        ! 1000 continue
                    ! end do
                     ! ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
                ! ! !$OMP END PARALLEL DO
                ! ! synthese des contributions
                ! ! !call MPI_REDUCE(som_cont,som_cont_0,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code)
                ! !call MPI_ALLREDUCE(som_cont,som_cont_0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
                ! ! call sleep(1)
                ! som_cont=som_cont_0
                ! !!call MPI_ABORT(MPI_COMM_WORLD,erreur,comm)
                ! !!call MPI_Barrier(MPI_COMM_WORLD,comm2) ! pour la synchronisation globale avant     
        end select
    
    !************* FIN INTEGRALES **************************
                   
        res=0.d0
        select case(methodInt)
        case(0) ! estimation par monte carlo
            ! call dblepr("integrale3=", -1, integrale3, ntrials)
            ! call dblepr("log(integrale3)=", -1, dlog(integrale3), ntrials)
            res = sum(dlog(integrale3))
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajsplines_copule_surrogate=-1.d9
                goto 123
            end if 
        case(1) ! estimation par monte carlo
            res = sum(dlog(integrale3))
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajsplines_copule_surrogate=-1.d9
                goto 123
            end if 
        case(3) ! estimation par approximation de laplace
            res = sum(dlog(integrale3))
            ! if(control_affichage == 0) then
                ! control_affichage = 1
                ! call dblepr("integrale3=", -1, integrale3, ntrials)        
                ! call dblepr("dlog(integrale3)=", -1, dlog(integrale3), ntrials)    
                ! call dblepr("res=", -1, res, 1)
            ! endif
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajsplines_copule_surrogate=-1.d9
                goto 123
            end if
        end select
    endif

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
    ! if(control_affichage == 0) then
        ! control_affichage = 1
        ! call dblepr("resnonpen=", -1, resnonpen, 1)        
        ! call dblepr("k0=", -1, k0, 2)    
        ! call dblepr("pe2=", -1, pe2, 1)
        ! call dblepr("pe1=", -1, pe1, 1)
        ! call dblepr("pe=", -1, pe, 1)
        ! call dblepr("res=", -1, res, 1)    
    ! endif
    ! call dblepr("k0 ", -1, k0, 2)
    ! call dblepr("res ", -1, res, 1)
    deallocate(mat_A)
    if ((res.ne.res).or.(abs(res).ge. 1.d30).or.(res .ge. 0.d0)) then
        funcpajsplines_copule_surrogate=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123

    else
        funcpajsplines_copule_surrogate = res 
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
    return
    
    end function funcpajsplines_copule_surrogate



