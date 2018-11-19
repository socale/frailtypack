subroutine jointsurrokendall(theta, gamma, alpha, eta, adaptative, npg, ndim,&
                        ui_chap_Essai, invBi_chol_Essai_k, xx1, ww1,&
                        sigma_v, methodeInt, N_MC_kendall, method_int_kendal,&
                        random_gen,aleatoire,nbre_sim,seed, &
                        tau_kendal_00, tau_kendal_01, tau_kendal_10, tau_kendal_11, ss)
    
    use fonction_A_integrer, only:funcJointSurroKendall
    use Autres_fonctions, only:tau_kendall, init_random_seed
    use var_surrogate, only: random_generator
    ! !$ use OMP_LIB
    
    implicit none
    
    integer ::ii,jj,kk,ll
    integer,intent(in):: ndim,npg,adaptative,methodeInt, method_int_kendal,&
                         N_MC_kendall,random_gen,seed,aleatoire,nbre_sim
    double precision, intent(out):: ss,tau_kendal_00, tau_kendal_01, tau_kendal_10,&
                        tau_kendal_11
    double precision,intent(in)::theta, gamma, alpha, eta
    double precision,dimension(2,2),intent(in)::sigma_v
    double precision,dimension(ndim,ndim),intent(in)::invBi_chol_Essai_k ! pour recuperer les matrice B_k dans le vecteur des matrices B des essais K
    double precision,dimension(ndim),intent(in)::ui_chap_Essai 
    double precision,dimension(npg),intent(in)::xx1, ww1
    double precision,dimension(ndim)::xxl, m
    double precision,dimension(:,:),allocatable::theta_0,gamma_0
    double precision::ss1,ss2,ss3,auxfunca
    
    double precision, external::gauss_HermMultA_surr    
    
    ss = 0.d0
    ! call dblepr("ss",-1,ss,1)
    if(methodeInt==1) then ! integration par gauss hermite
        if(ndim == 4) then ! model complet avec prise en compte de l'heterogeneite sur les risque de base
            ! ! $OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl,kk,ss2,ll,ss3) firstprivate(auxfunca) &
            ! ! $OMP SHARED(npg,xx1,ww1,invBi_chol_Essai_k,ui_chap_Essai,adaptative, theta, gamma, alpha, eta)&
            ! ! $OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do ll=1,npg
                ss3=0.d0    
                do kk=1,npg
                    ss2=0.d0
                    do jj=1,npg
                        ss1=0.d0
                        do ii=1,npg    
                            xxl(1)=xx1(ll)        
                            xxl(2)=xx1(kk)                    
                            xxl(3)=xx1(jj)
                            xxl(4)=xx1(ii)
                            ! call dblepr("premier vecteur de points",-1,xxl,1)
                            ! changement de variable en cas de quadrature adaptative
                            if(adaptative==1) then ! on effectue le changement de variable
                                m=matmul(invBi_chol_Essai_k,xxl)
                                xxl=ui_chap_Essai(:)+dsqrt(2.d0)*m        
                            end if
                            auxfunca=funcJointSurroKendall(xxl(1), xxl(2), xxl(3), xxl(4), theta, gamma, alpha, eta, 1.d0)
                            ! call dblepr("auxfunca",-1,auxfunca,1)
                            ss1 = ss1 + ww1(ii) * auxfunca
                        end do
                        ss2 = ss2+ww1(jj)*ss1
                    end do
                    ss3 = ss3+ww1(kk)*ss2
                end do
                ss = ss+ww1(ll)*ss3
            end do
            !!$OMP END PARALLEL DO
        
        
        else
            ! !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl,kk,ss2) firstprivate(auxfunca) &
            ! !$OMP SHARED(npg,xx1,ww1,invBi_chol_Essai_k,ui_chap_Essai,adaptative, theta, gamma, alpha, eta)&
            ! !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)   
            do jj=1,npg
                ss1=0.d0
                do ii=1,npg                                
                    xxl(1)=xx1(jj)
                    xxl(2)=xx1(ii)
                    !changement de variable en cas de quadrature adaptative
                    if(adaptative==1) then ! on effectue le changement de variable
                        m=matmul(invBi_chol_Essai_k,xxl)
                        xxl=ui_chap_Essai(:)+dsqrt(2.d0)*m        
                    end if
                    auxfunca=funcJointSurroKendall(0.d0, xxl(1), 0.d0, xxl(2), theta, gamma, alpha, eta, 0.d0)
                    ss1 = ss1+ww1(ii)*(auxfunca)
                end do
                ss = ss+ww1(jj)*ss1
            end do
            ! !$OMP END PARALLEL DO
        endif
        ss = 2.d0*ss -1.d0
    else ! dans ce cas a O, integration pat monte carlo
        random_generator = random_gen
        if(random_generator==1) call init_random_seed(seed,aleatoire,nbre_sim)! initialisation de l'environnement de generation
        
        allocate(theta_0(1,1),gamma_0(1,1))
        theta_0(1,1)= theta
        gamma_0(1,1)= gamma
        
        if(method_int_kendal==4 .or. method_int_kendal==5)then ! 1 seul taux de kendall
            tau_kendal_00=tau_kendall(theta_0,gamma_0,sigma_v,0,0,method_int_kendal,N_MC_kendall,alpha,eta,0)!tau de kendal des non traites z_11=0,z_21=0
            ss = tau_kendal_00
        else
            tau_kendal_11=tau_kendall(theta_0,gamma_0,sigma_v,1,1,method_int_kendal,N_MC_kendall,alpha,eta,0)!tau de kendal des traites z_11=1,z_21=1
            tau_kendal_10=tau_kendall(theta_0,gamma_0,sigma_v,1,0,method_int_kendal,N_MC_kendall,alpha,eta,0)!tau de kendal des 1 traite et l'autre non traite z_11=1,z_21=0
            tau_kendal_01=tau_kendall(theta_0,gamma_0,sigma_v,0,1,method_int_kendal,N_MC_kendall,alpha,eta,0)!tau de kendal des traite z_11=0,z_21=1
            tau_kendal_00=tau_kendall(theta_0,gamma_0,sigma_v,0,0,method_int_kendal,N_MC_kendall,alpha,eta,0)!tau de kendal des non traites z_11=0,z_21=0
        endif
    endif
            
end subroutine jointsurrokendall