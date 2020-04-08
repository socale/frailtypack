module fonction_A_integrer
! ici je definis les fonctions a integrer (integrants)
implicit none
contains
    
!==================================================================
!=====================functions a integrer surrogacy===============
!==================================================================

     ! 1 effet aleatoire partage niveau individeul +2 effet aleatoire correles en interaction avec le traitement:
     ! copula model
 double precision function Integrant_Copula(vsi,vti,ui,ig,nsujet_trial)
    ! vsi= frailtie niveau essai associe a s
    ! vti= frailtie niveau essai associe a t
    ! ui = random effect associated xith the baseline hazard
    ! ig = current cluster
    ! nsujet_trial = number of subjects in the current trial
    
    use var_surrogate, only: posind_i, alpha_ui, const_res4, const_res5, res2_dcs_sujet,res2s_sujet, &
        theta_copule, delta, deltastar, copula_function, methodInt, pi, gamma_ui, determinant, &
        varcovinv, adaptative, control_affichage
    use comon, only: ve
    
    IMPLICIT NONE
    integer,intent(in):: ig, nsujet_trial
    double precision,intent(in)::vsi,vti,ui
    integer::j,n
    double precision::integrant, f_Sij, f_Tij, fbar_Sij, fbar_Tij, C_theta, phimun_S, phimun_T,phiprim_ST,&
                      sumphimun_ST, phisecond_ST, phiprimphimun_S, derivphi_ij, contri_indiv,phiprimphimun_T,&
                      f_V, tempon, logfbar_Sij, logfbar_Tij
    double precision, dimension(:,:),allocatable::m1,m3  
    double precision, dimension(:,:),allocatable::m    
     
    C_theta = 0.d0
    phimun_S = 0.d0
    phimun_T = 0.d0
    phiprim_ST = 0.d0
    sumphimun_ST = 0.d0
    phisecond_ST = 0.d0
    phiprimphimun_S = 0.d0
    phiprimphimun_T = 0.d0
    logfbar_Sij = 0.d0
    logfbar_Tij = 0.d0
    integrant = 1.d0
    do j = 1, nsujet_trial
        ! Expression in the log-vraisamblance
        f_Sij = res2s_sujet(posind_i-1+j) * dexp(ui + vsi*dble(ve(posind_i-1+j,1))) &
                * dexp(- const_res4(posind_i-1+j) * dexp(ui + vsi*dble(ve(posind_i-1+j,1))))
        f_Tij = res2_dcs_sujet(posind_i-1+j) * dexp(alpha_ui * ui + vti*dble(ve(posind_i-1+j,1)))&
                * dexp(- const_res5(posind_i-1+j) * dexp(alpha_ui * ui + vti*dble(ve(posind_i-1+j,1))))
        fbar_Sij = dexp(- const_res4(posind_i-1+j) * dexp(ui + vsi*dble(ve(posind_i-1+j,1))))
        fbar_Tij = dexp(- const_res5(posind_i-1+j) * dexp(alpha_ui * ui + vti*dble(ve(posind_i-1+j,1))))
        select case(copula_function)
            case(1) ! clayton copula 
                C_theta = (fbar_Sij**(-theta_copule) + fbar_Tij**(-theta_copule) - 1.d0)&
                        **(-1.d0/theta_copule)
                phimun_S = (fbar_Sij**(-theta_copule) - 1.d0) / theta_copule
                phimun_T = (fbar_Tij**(-theta_copule) - 1.d0) / theta_copule
                phiprim_ST = - (fbar_Sij**(-theta_copule) + fbar_Tij**(-theta_copule) - 1.d0)**&
                          (-(1.d0+theta_copule)/theta_copule)
                sumphimun_ST = phimun_S + phimun_T
                phisecond_ST = (1.d0 + theta_copule) * (1.d0 + theta_copule * sumphimun_ST)**(- (1.d0 + &
                             2.d0 * theta_copule)/theta_copule)
                phiprimphimun_S = - fbar_Sij**(1.d0 + theta_copule)
                phiprimphimun_T = - fbar_Tij**(1.d0 + theta_copule)
            case(2) ! Gumbel copula
                logfbar_Sij = - const_res4(posind_i-1+j) * dexp(ui + vsi*dble(ve(posind_i-1+j,1)))
                logfbar_Tij = - const_res5(posind_i-1+j) * dexp(alpha_ui * ui + vti*dble(ve(posind_i-1+j,1)))
                C_theta = dexp(-((-logfbar_Sij)**(theta_copule + 1.d0) + (- logfbar_Tij)**&
                        (theta_copule + 1.d0))**(1.d0/(theta_copule + 1.d0)))
                phimun_S = (-logfbar_Sij)**(1.d0 + theta_copule)
                phimun_T = (-logfbar_Tij)**(1.d0 + theta_copule)
                phiprim_ST = - 1.d0/(1.d0+theta_copule) * ((- logfbar_Sij)**(1.d0+theta_copule) + &
                          (-logfbar_Tij)**(1.d0+theta_copule))**(-theta_copule/(1.d0+theta_copule)) &
                          * C_theta
                sumphimun_ST = phimun_S + phimun_T
                phisecond_ST = (1.d0 / (1.d0 + theta_copule)**2.d0) * (theta_copule * sumphimun_ST**&
                             (-(2.d0 * theta_copule + 1.d0)/(theta_copule + 1.d0)) + sumphimun_ST**&
                             (-(2.d0 * theta_copule)/(theta_copule + 1.d0))) * dexp(- sumphimun_ST**&
                             (1.d0/(1.d0 + theta_copule)))
                phiprimphimun_S = (-fbar_Sij/(1.d0 + theta_copule)) * (- logfbar_Sij)**(- theta_copule)
                phiprimphimun_T = (-fbar_Tij/(1.d0 + theta_copule)) * (- logfbar_Tij)**(- theta_copule)
        end select
        
        ! expression with derrivatives
        derivphi_ij = dble(delta(posind_i-1+j)) * dble(deltastar(posind_i-1+j)) * phisecond_ST + (dble(delta(posind_i-1+j))&
                    * (1.d0 - dble(deltastar(posind_i-1+j))) + (1.d0 - dble(delta(posind_i-1+j))) * &
                    dble(deltastar(posind_i-1+j))) *  phiprim_ST + (1.d0 - dble(delta(posind_i-1+j))) *&
                    (1.d0 - dble(deltastar(posind_i-1+j))) * C_theta
        
        ! individual contributions
        if(phiprimphimun_S > -1.d-299) phiprimphimun_S = -1.d-299
        if(phiprimphimun_T > -1.d-299) phiprimphimun_T = -1.d-299
        
        contri_indiv = derivphi_ij * (f_Sij / phiprimphimun_S)**dble(delta(posind_i-1+j)) * &
                       (f_Tij / phiprimphimun_T)&
                       **dble(deltastar(posind_i-1+j))
        
        tempon = integrant
        integrant = integrant * contri_indiv
        
        if(adaptative .and. integrant==0) then
            ! call intpr("posind_i-1+j ", -1, posind_i-1+j, 1)
            ! call dblepr("contri_indiv = ", -1, contri_indiv, 1)
            ! call dblepr("integrant = ", -1, integrant, 1)
        endif
        !if((integrant .ne. integrant) .and. (control_affichage == 0)) then
        ! ! if(control_affichage == 0) then
            ! control_affichage = 1
            ! call dblepr("vsi = ", -1, vsi, 1)
            ! call dblepr("vti = ", -1, vti, 1)
            ! call dblepr("ui = ", -1, ui, 1)
            ! call intpr("ig = ", -1, ig, 1)
            ! call intpr("nsujet_trial = ", -1, nsujet_trial, 1)
            ! call intpr("posind_i = ", -1, posind_i, 1)
            ! call dblepr("alpha_ui = ", -1, alpha_ui, 1)
            ! call dblepr("const_res4 = ", -1, const_res4, 1)
            ! call dblepr("const_res5 = ", -1, const_res5, 1)
            ! call dblepr("res2_dcs_sujet = ", -1, res2_dcs_sujet, 1)
            ! call dblepr("res2s_sujet = ", -1, res2s_sujet, 1)
            ! call dblepr("theta_copule = ", -1, theta_copule, 1)
            ! call intpr("delta(posind_i-1+1) = ", -1, delta(posind_i-1+1), 1)
            ! call intpr("deltastar(posind_i-1+1) = ", -1, deltastar(posind_i-1+1), 1)
            ! call intpr("delta(posind_i-1+2) = ", -1, delta(posind_i-1+2), 1)
            ! call intpr("deltastar(posind_i-1+2) = ", -1, deltastar(posind_i-1+2), 1)
            ! call intpr("copula_function = ", -1, copula_function, 1)
            ! call intpr("methodInt = ", -1, methodInt, 1)
            ! call dblepr("pi = ", -1, pi, 1)
            ! call dblepr("gamma_ui = ", -1, gamma_ui, 1)
            ! call dblepr("determinant = ", -1, determinant, 1)
            ! if (methodInt == 1) call dblepr("varcovinv = ", -1, varcovinv, 9)
            ! call intpr("adaptative = ", -1, adaptative, 1)
            ! call dblepr("ve(posind_i-1+1,1) = ", -1, ve(posind_i-1+1,1), 1)
            ! call dblepr("ve(posind_i-1+2,1) = ", -1, ve(posind_i-1+2,1), 1)
            ! call dblepr("f_Sij = ", -1, f_Sij, 1)
            ! call dblepr("f_Tij = ", -1, f_Tij, 1)
            ! call dblepr("fbar_Sij = ", -1, fbar_Sij, 1)
            ! call dblepr("fbar_Tij = ", -1, fbar_Tij, 1)
            ! call dblepr("C_theta = ", -1, C_theta, 1)
            ! call dblepr("phimun_S = ", -1, phimun_S, 1)
            ! call dblepr("phimun_T = ", -1, phimun_T, 1)
            ! call dblepr("phiprim_ST = ", -1, phiprim_ST, 1)
            ! call dblepr("sumphimun_ST = ", -1, sumphimun_ST, 1)
            ! call dblepr("phisecond_ST = ", -1, phisecond_ST, 1)
            ! call dblepr("phiprimphimun_S = ", -1, phiprimphimun_S, 1)
            ! call dblepr("phiprimphimun_T = ", -1, phiprimphimun_T, 1) 
            ! call dblepr("f_Sij / phiprimphimun_S = ", -1, f_Sij / phiprimphimun_S, 1)
            ! call dblepr("dble(delta(posind_i-1+j)) = ", -1, dble(delta(posind_i-1+j)), 1)
            ! call dblepr("f_Tij /phiprimphimun_T = ", -1, f_Tij / phiprimphimun_T, 1)
            ! call dblepr("dble(deltastar(posind_i-1+j)) = ", -1, dble(deltastar(posind_i-1+j)), 1)
            ! call dblepr("derivphi_ij = ", -1, derivphi_ij, 1)
            ! call dblepr("contri_indiv = ", -1, contri_indiv, 1)
            ! call dblepr("tempon = ", -1, tempon, 1)
            ! call dblepr("integrant = ", -1, integrant, 1)
        ! endif
        
        ! if((integrant < 0.d0) .and. (control_affichage == 0)) then
            ! control_affichage = 1
            ! call dblepr("vsi = ", -1, vsi, 1)
            ! call dblepr("vti = ", -1, vti, 1)
            ! call dblepr("ui = ", -1, ui, 1)
            ! call intpr("ig = ", -1, ig, 1)
            ! call intpr("nsujet_trial = ", -1, nsujet_trial, 1)
            ! call intpr("posind_i = ", -1, posind_i, 1)
            ! call dblepr("alpha_ui = ", -1, alpha_ui, 1)
            ! call dblepr("const_res4 = ", -1, const_res4, 1)
            ! call dblepr("const_res5 = ", -1, const_res5, 1)
            ! call dblepr("res2_dcs_sujet = ", -1, res2_dcs_sujet, 1)
            ! call dblepr("res2s_sujet = ", -1, res2s_sujet, 1)
            ! call dblepr("theta_copule = ", -1, theta_copule, 1)
            ! call intpr("delta(posind_i-1+1) = ", -1, delta(posind_i-1+1), 1)
            ! call intpr("deltastar(posind_i-1+1) = ", -1, deltastar(posind_i-1+1), 1)
            ! call intpr("delta(posind_i-1+2) = ", -1, delta(posind_i-1+2), 1)
            ! call intpr("deltastar(posind_i-1+2) = ", -1, deltastar(posind_i-1+2), 1)
            ! call intpr("copula_function = ", -1, copula_function, 1)
            ! call intpr("methodInt = ", -1, methodInt, 1)
            ! call dblepr("pi = ", -1, pi, 1)
            ! call dblepr("gamma_ui = ", -1, gamma_ui, 1)
            ! call dblepr("determinant = ", -1, determinant, 1)
            ! if (methodInt == 1) call dblepr("varcovinv = ", -1, varcovinv, 9)
            ! call intpr("adaptative = ", -1, adaptative, 1)
            ! call dblepr("ve(posind_i-1+1,1) = ", -1, ve(posind_i-1+1,1), 1)
            ! call dblepr("ve(posind_i-1+2,1) = ", -1, ve(posind_i-1+2,1), 1)
            ! call dblepr("f_Sij = ", -1, f_Sij, 1)
            ! call dblepr("f_Tij = ", -1, f_Tij, 1)
            ! call dblepr("fbar_Sij = ", -1, fbar_Sij, 1)
            ! call dblepr("fbar_Tij = ", -1, fbar_Tij, 1)
            ! call dblepr("C_theta = ", -1, C_theta, 1)
            ! call dblepr("phimun_S = ", -1, phimun_S, 1)
            ! call dblepr("phimun_T = ", -1, phimun_T, 1)
            ! call dblepr("phiprim_ST = ", -1, phiprim_ST, 1)
            ! call dblepr("sumphimun_ST = ", -1, sumphimun_ST, 1)
            ! call dblepr("phisecond_ST = ", -1, phisecond_ST, 1)
            ! call dblepr("phiprimphimun_S = ", -1, phiprimphimun_S, 1)
            ! call dblepr("phiprimphimun_T = ", -1, phiprimphimun_T, 1)
            ! call dblepr("f_Sij / phiprimphimun_S = ", -1, f_Sij / phiprimphimun_S, 1)
            ! call dblepr("dble(delta(posind_i-1+j)) = ", -1, dble(delta(posind_i-1+j)), 1)
            ! call dblepr("f_Tij /phiprimphimun_T = ", -1, f_Tij / phiprimphimun_T, 1)
            ! call dblepr("dble(deltastar(posind_i-1+j)) = ", -1, dble(deltastar(posind_i-1+j)), 1)
            ! call dblepr("derivphi_ij = ", -1, derivphi_ij, 1)
            ! call dblepr("contri_indiv = ", -1, contri_indiv, 1)
            ! call dblepr("tempon = ", -1, tempon, 1)
            ! call dblepr("integrant = ", -1, integrant, 1)
        ! endif
    enddo
    
    if (methodInt == 1 .or. methodInt == 3)then
        allocate(m(1,1),m1(1,2),m3(1,2))
        m1(1,1)=vsi
        m1(1,2)=vti
        m3=MATMUL(m1,varcovinv)
        m=MATMUL(m3,TRANSPOSE(m1))
        f_V = 1.d0/(2.d0 * pi *  dsqrt(2.d0 * pi * gamma_ui * determinant)) * &
        dexp(- 1.d0/2.d0 * m(1,1) - 1.d0/2.d0 * ui**2.d0 / gamma_ui)
        Integrant_Copula = integrant * f_V
        deallocate(m,m1,m3)
    endif
    
    if (methodInt == 0) Integrant_Copula = integrant

    ! if(control_affichage == 0)then
        ! control_affichage = 1
        ! call dblepr("vsi = ", -1, vsi, 1)
        ! call dblepr("vti = ", -1, vti, 1)
        ! call dblepr("ui = ", -1, ui, 1)
        ! call intpr("ig = ", -1, ig, 1)
        ! call intpr("nsujet_trial = ", -1, nsujet_trial, 1)
        ! call intpr("posind_i = ", -1, posind_i, 1)
        ! call dblepr("alpha_ui = ", -1, alpha_ui, 1)
        ! call dblepr("const_res4 = ", -1, const_res4, 1)
        ! call dblepr("const_res5 = ", -1, const_res5, 1)
        ! call dblepr("res2_dcs_sujet = ", -1, res2_dcs_sujet, 1)
        ! call dblepr("res2s_sujet = ", -1, res2s_sujet, 1)
        ! call dblepr("theta_copule = ", -1, theta_copule, 1)
        ! call intpr("delta(posind_i-1+1) = ", -1, delta(posind_i-1+1), 1)
        ! call intpr("deltastar(posind_i-1+1) = ", -1, deltastar(posind_i-1+1), 1)
        ! call intpr("delta(posind_i-1+2) = ", -1, delta(posind_i-1+2), 1)
        ! call intpr("deltastar(posind_i-1+2) = ", -1, deltastar(posind_i-1+2), 1)
        ! call intpr("copula_function = ", -1, copula_function, 1)
        ! call intpr("methodInt = ", -1, methodInt, 1)
        ! call dblepr("pi = ", -1, pi, 1)
        ! call dblepr("gamma_ui = ", -1, gamma_ui, 1)
        ! call dblepr("determinant = ", -1, determinant, 1)
        ! if (methodInt == 1) call dblepr("varcovinv = ", -1, varcovinv, 9)
        ! call intpr("adaptative = ", -1, adaptative, 1)
        ! call dblepr("ve(posind_i-1+1,1) = ", -1, ve(posind_i-1+1,1), 1)
        ! call dblepr("ve(posind_i-1+2,1) = ", -1, ve(posind_i-1+2,1), 1)
        ! call dblepr("f_Sij = ", -1, f_Sij, 1)
        ! call dblepr("f_Tij = ", -1, f_Tij, 1)
        ! call dblepr("fbar_Sij = ", -1, fbar_Sij, 1)
        ! call dblepr("fbar_Tij = ", -1, fbar_Tij, 1)
        ! call dblepr("C_theta = ", -1, C_theta, 1)
        ! call dblepr("phimun_S = ", -1, phimun_S, 1)
        ! call dblepr("phimun_T = ", -1, phimun_T, 1)
        ! call dblepr("phiprim_ST = ", -1, phiprim_ST, 1)
        ! call dblepr("sumphimun_ST = ", -1, sumphimun_ST, 1)
        ! call dblepr("phisecond_ST = ", -1, phisecond_ST, 1)
        ! call dblepr("phiprimphimun_T = ", -1, phiprimphimun_T, 1)
        ! call dblepr("derivphi_ij = ", -1, derivphi_ij, 1)
        ! call dblepr("integrant = ", -1, integrant, 1)        
    ! endif
    return
    end function Integrant_Copula
    
    
    double precision function funcSurrNN_MC(frail,n,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! n= taille du vecteur frail qui vaut nombre de sujet+2 pour vs et vt
    ! i = position du premier individu du cluster, sera increentee de la taille du cluster + 1 a la sortie de la fonction
     use var_surrogate
    use comon, only: ve,vedc,eta

    IMPLICIT NONE
    integer ::j,k !n=taille du tableau frail
    integer,intent(in):: i
    integer,intent(in)::n
    double precision:: s1,C1,c2,c3,c4,c5,vs,vt
    double precision,intent(in),dimension(n)::frail

    ! fin declaration et debut de la fonction
    s1=0.d0
    C1=0.d0
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c5=0.d0
    vs=0.d0
    vt=0.d0
    s1=0.d0
    j=0
    k=0
    vs=frail(n-1)
    vt=frail(n)

    k=1
    !!print*,"i=",i,"n=",n,"i+n-2-1=",i+n-2-1
    do j=i,(i+n-2-1) ! nombre d'individu du cluster
      c3=frail(k)*(delta(j)+deltastar(j)*eta)
      c2=(vs*delta(j)+vt*deltastar(j))*ve(j,1)
      c4=const_res4(j)*dexp(vs*ve(j,1))*dexp(frail(k))
      c5=const_res5(j)*dexp(vt*vedc(j,1))*dexp(frail(k)*eta)
      s1=s1+c3+c2-c4-c5
      k=k+1
!      !print*,"c3=",c3,"c2=",c2,"c4=",c4,"c=5",c5,"s1=",s1
    end do
    !!print*,"k=",k
    !stop
    !!print*,"integrant ligne 58"
!   end if
    funcSurrNN_MC=dexp(s1)
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    return
    end function funcSurrNN_MC
    
    ! cas un effets aleatoires partage
    double precision function funcSurrNN_MC_Essai_t1(frail,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! i = position du cluster sur lequel on se trouve
     use var_surrogate
    use comon, only: ve,vedc,alpha

    IMPLICIT NONE
    integer ::n 
    integer,intent(in):: i
    !integer,intent(in)::n
    double precision:: c2,c3,c4,vs
    double precision,intent(in),dimension(1)::frail

    ! fin declaration et debut de la fonction
    !s1=0.d0
    !c1=0.d0
    c2=0.d0
    c3=0.d0
    c4=0.d0
    !c5=0.d0
    !vs=0.d0
    !j=0
    !k=0
    !vs=frail(1)

    !k=1
    !do j=i,(i+n-1) ! nombre d'individu du cluster
    !  c2=(delta(j)+alpha*deltastar(j))*vs*ve(j,1)
    !  c4=const_res4(j)*dexp(vs*ve(j,1))
    !  c5=const_res5(j)*dexp(vs*alpha*vedc(j,1))
    !  s1=s1+c2-c4-c5
    !end do
    
    n=nsujeti(i)
    vs=frail(1)
    !c1=-vs**2/(2.d0*sigma2)
    c2=(nigts(i)+alpha*cdcts(i))*vs
    c3=SUM(const_res4(posind_i:(posind_i+n-1))*dexp(vs*dble(ve(posind_i:(posind_i+n-1),1))))
    c4=SUM(const_res5(posind_i:(posind_i+n-1))*dexp(alpha*vs*dble(vedc(posind_i:(posind_i+n-1),1))))

    funcSurrNN_MC_Essai_t1=dexp(c2-c3-c4)
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    return
    end function funcSurrNN_MC_Essai_t1
    
    ! cas un effet aleatoires partages niveau individuel
    double precision function funcSurrNN_MC_Essai_indiv_1(frail,i)
    ! fonction a integrer: cas des effets aleatoires individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! i = position du cluster sur lequel on se trouve
     use var_surrogate
    use comon, only: eta

    IMPLICIT NONE
    integer ::n_i 
    integer,intent(in):: i
    !integer,intent(in)::n
    double precision:: C1,c2,c3,c4
    double precision,intent(in),dimension(nsujeti(i))::frail

    ! fin declaration et debut de la fonction
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c1=0.d0
    !dimInt=size(frail)
    n_i=nsujeti(i)
    !vs=frail(dimInt)
    ! calculs
    c1=SUM(frail(1:n_i)*(delta(posind_i:(posind_i+n_i-1))+deltastar(posind_i:(posind_i+n_i-1))*eta))
    !c2=(nigts(i)+alpha*cdcts(i))*vs
    c3=SUM(const_res4(posind_i:(posind_i+n_i-1))*dexp(frail(1:n_i)))
    c4=SUM(const_res5(posind_i:(posind_i+n_i-1))*dexp(eta*frail(1:n_i)))

    funcSurrNN_MC_Essai_indiv_1=dexp(c1-c3-c4)
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    return
    end function funcSurrNN_MC_Essai_indiv_1
    
    ! cas un effet aleatoires partages niveau individuel: integrale par quadrature non adaptative
    double precision function funcSurrNN_MC_Essai_indiv_1QNA(frail,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! i = position du cluster sur lequel on se trouve
     use var_surrogate
    use comon, only: eta

    IMPLICIT NONE
    integer ::n_i
    integer,intent(in):: i
    !integer,intent(in)::n
    double precision:: c1,c3,c4
    double precision,intent(in),dimension(:)::frail
    !double precision,dimension(1,nsujeti(i))::c5,frail_mat
    double precision,dimension(1,1)::c2

    ! fin declaration et debut de la fonction
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c1=0.d0
    !dimInt=size(frail)
    n_i=nsujeti(i)
    !frail_mat(1,:)=frail
    !c5=MATMUL(frail_mat,vcinv)
    !c2=-(1/2.d0)*MATMUL(c5,TRANSPOSE(frail_mat))
    c2=SUM(-frail(1:n_i)**2/(2.d0*theta2))
    !!print*,"==============================="
    !!print*,"c2=",c2
    !stop
    ! calculs
    c1=SUM(frail(1:n_i)*(delta(posind_i:(posind_i+n_i-1))+deltastar(posind_i:(posind_i+n_i-1))*eta))
    c3=SUM(const_res4(posind_i:(posind_i+n_i-1))*dexp(frail(1:n_i)))
    c4=SUM(const_res5(posind_i:(posind_i+n_i-1))*dexp(eta*frail(1:n_i)))
    funcSurrNN_MC_Essai_indiv_1QNA=dexp(c2(1,1)+c1-c3-c4)
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    return
    end function funcSurrNN_MC_Essai_indiv_1QNA
    
    ! cas un effet aleatoires partages niveau individuel: integrant pour un individu
    double precision function integrant_indiv_1(frailij,j)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! j = individu j du cluster i
     use var_surrogate
    use comon, only: eta,lognormal,theta

    IMPLICIT NONE
    integer,intent(in):: j
    double precision,intent(in)::frailij
    double precision:: c1,c2,c3,c4

    ! fin declaration et debut de la fonction
    !!print*,"j=",j,"posind_i=",posind_i-1,"posind_i+j=",posind_i-1+j
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c1=0.d0
    if(lognormal==1)then ! integrant cas lognormal
        c1=-(frailij**2.d0)/(2.d0*theta2)
        c2=frailij*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)
        c3=const_res4(posind_i-1+j)*dexp(frailij)
        c4=const_res5(posind_i-1+j)*dexp(eta*frailij)
        integrant_indiv_1=dexp(c1+c2-c3-c4)
    else ! integrant cas gamma
        c1=dlog(frailij)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta+1./theta-1.)
        c2=const_res4(posind_i-1+j)*frailij
        c3=const_res5(posind_i-1+j)*frailij**eta
        c4=frailij/theta
        integrant_indiv_1=dexp(c1-c2-c3-c4)
    endif
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    end function integrant_indiv_1
        
    ! cas un effet aleatoires partages niveau individuel: integrant pour un individu
    double precision function integrant_indiv_1A(j,npoint1)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! j = individu j du cluster i
    ! npoint = nombre de point de quadrature
     use var_surrogate
    use comon, only: eta,lognormal,theta,invBi_chol
    IMPLICIT NONE
    integer,intent(in):: j,npoint1
    double precision,dimension(npoint1):: xx

    ! fin declaration et debut de la fonction
    !!print*,"j=",j,"posind_i=",posind_i-1,"posind_i+j=",posind_i-1+j
    if(adaptative) then ! on effectue le changement de variable
        xx=ui_chap(posind_i-1+j,1)+dsqrt(2.d0)*invBi_chol(posind_i-1+j,posind_i-1+j)**(1.d0/2.d0)    
    else
        xx=xx1
    end if
        
    ! Je calcule pour chaque individu du cluster par vectorisation son integrale par gaussHermite non adaptative
    if(lognormal==1)then ! integrant cas lognormal
        integrant_indiv_1A=SUM(ww1(1:npoint1)*dexp(&
            - (xx(1:npoint1)**2.d0)/(2.d0*theta2)&
            + xx(1:npoint1)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)&
            - const_res4(posind_i-1+j)*dexp(xx(1:npoint1))&
            - const_res5(posind_i-1+j)*dexp(eta*xx(1:npoint1))&
        ))
    else ! integrant cas gamma
        integrant_indiv_1A=SUM(ww1(1:npoint1)*&
            dexp(&
            dlog(xx(1:npoint1))*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta+1./theta-1.)&
                - const_res4(posind_i-1+j)*xx(1:npoint1)&
                - const_res5(posind_i-1+j)*xx(1:npoint1)**eta&
                - xx(1:npoint1)/theta))
                
                !!print*,xx(1:npoint1)
                !!print*,dlog(xx(1:npoint1))
    endif
    !!print*,"ww=",ww1(1:npoint1)
    !!print*,"xx=",xx1(1:npoint1)
    !!print*,"integrant_indiv_1A: je vaut",integrant_indiv_1A
    end function integrant_indiv_1A
    
    !==============Integration par gauss hermite non adaptative=======================
    
     ! 1 effet aleatoire partage (ou correle) niveau individeul +2 effet aleatoire correles niveau essai en interaction avec le traitement+ 2 effets aleatoires associes au risque de base:
     ! integrale au niveau individuel
 double precision function Integrale_Individuel(vsi,vti,ui,j,npoint1)
    ! vsi= frailtie niveau essai associe a s
    ! vti= frailtie niveau essai associe a t
    ! j = individu j du cluster i
    ! npoint = nombre de point de quadrature
     use var_surrogate, only:posind_i,ui_chap,theta2,const_res5,const_res4,& !vs_i,vt_i,individu_j
        deltastar,delta,adaptative,xx1,ww1,& !estim_wij_chap,nparamfrail,ntrials,nsujeti
        invBi_chol_Individuel,& !invBi_chol_Essai,invBi_cholDet_Essai,essai_courant,ui_chap_Essai
        alpha_ui,frailt_base,switch_adaptative
    use comon, only: eta,lognormal,ve,eta,invBi_cholDet !I_hess,H_hess,model,nsujet,theta
    !use parameters, only: maxiter
    !use mod_Adaptative ! pour lintegrant a maximiser pour l'adaptative
!    use func_adaptative
    use Autres_fonctions
    !use tailles, only:npmax 
    !use optim_SCL_0
    IMPLICIT NONE
    integer,intent(in):: j,npoint1
    double precision,intent(in)::vsi,vti,ui
    double precision,dimension(npoint1):: xx
    !integer::i,ier,istop,ss,sss,ni,model_save,nparamfrail_save,maxiter_save,k,nmax,indicej, & 
    !         np,indice_B_essai,indice_ind_util_essai,non_conv
    !integer::frail_essai_deja_est  ! variable qui dit si pour un essai donne l'on a deja estimes les vsi et vti (1) ou non (0)
    integer,parameter::effet=0,np_1=1
    !double precision::ca,cb,dd
    !double precision::res
    !double precision, dimension(2)::k0
!    double precision, dimension(1)::v,b
    !double precision, dimension(2)::b_i      ! pour les 2 parametres des effets aleatoires a predire niveau essai
    !double precision, dimension(2,2)::v_i    ! pour les 2 parametres des effets aleatoires a predire niveau essai
    !double precision, allocatable, dimension(:,:)::H_hessOut,HIH,HIHOut,IH,invBi_chol

    if(adaptative .and. switch_adaptative==1) then ! on effectue le changement de variable 
        xx=ui_chap(posind_i-1+j,1)+dsqrt(2.d0)*xx1*invBi_chol_Individuel(posind_i-1+j)            
        !!print*,"posind_i-1+j",posind_i-1+j
        !!print*,"ui_chap(posind_i-1+j,1)",ui_chap(posind_i-1+j,1)
    else
        xx=xx1
    end if
    
    ! if(frailt_base==1) xx=xx1
    
    !!print*,"xx=", xx,"xx1=",xx1
    ! Je calcule pour chaque individu du cluster par vectorisation son integrale 
    !!print*,"posind_i",posind_i,"j=",j
    !!print*,"delta(posind_i-1+j)",delta(posind_i-1+j)
    ! !print*,"ww1(1:npoint1)",ww1(1:npoint1)
    ! !print*,"xx(1:npoint1)",xx(1:npoint1)
    ! !print*,"xx(1:npoint1)",xx(1:npoint1)
    ! !print*,"delta(posind_i-1+j)",delta(posind_i-1+j)
    ! !print*,"deltastar(posind_i-1+j)",deltastar(posind_i-1+j)
    ! !print*,"const_res4(posind_i-1+j)",const_res4(posind_i-1+j)
    ! !print*,"ve(posind_i-1+j,1)",ve(posind_i-1+j,1)
    ! !print*,"const_res5(posind_i-1+j)",const_res5(posind_i-1+j)
    ! !print*,"vsi",vsi
    ! !print*,"vti",vti
        
    if(lognormal==1)then ! integrant cas lognormal
        if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de l'heterogeneite sur les risque des bas
            Integrale_Individuel=SUM(ww1(1:npoint1)*dexp((vsi*delta(posind_i-1+j)+vti*deltastar(posind_i-1+j))&
                *dble(ve(posind_i-1+j,1))&
                -(xx(1:npoint1)**2.d0)/(2.d0*theta2)&
                + xx(1:npoint1)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)&
                - const_res4(posind_i-1+j)*dexp(xx(1:npoint1)+vsi*dble(ve(posind_i-1+j,1)))&
                - const_res5(posind_i-1+j)*dexp(eta*xx(1:npoint1)+vti*dble(ve(posind_i-1+j,1)))&
            ))
        else
            Integrale_Individuel=SUM(ww1(1:npoint1)*dexp(ui*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*alpha_ui)&
                +(vsi*delta(posind_i-1+j)+vti*deltastar(posind_i-1+j))*dble(ve(posind_i-1+j,1))&
                -(xx(1:npoint1)**2.d0)/(2.d0*theta2)&
                + xx(1:npoint1)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)&
                - const_res4(posind_i-1+j)*dexp(xx(1:npoint1)+ui+vsi*dble(ve(posind_i-1+j,1)))&
                - const_res5(posind_i-1+j)*dexp(eta*xx(1:npoint1)+alpha_ui*ui+vti*dble(ve(posind_i-1+j,1)))&
            ))
        endif
    else ! integrant cas gamma
        !print*,"quadrature par la loi gamma non disponible pour le modele complet de surrogacy;",&
         !       "probleme de la listribution gamma bivariée"
    endif
    
    !!print*,"Integrale_Individuel avant=",Integrale_Individuel
    
    if(adaptative .and. switch_adaptative==1) then
        ! if(frailt_base==0) Integrale_Individuel=Integrale_Individuel*2.d0**(1.d0/2.d0)*invBi_cholDet(posind_i-1+j)  !2^(q/2), même determinant pour tous les individus
        Integrale_Individuel=Integrale_Individuel*2.d0**(1.d0/2.d0)*invBi_cholDet(posind_i-1+j)  !2^(q/2), même determinant pour tous les individus
    end if
   
    !!print*,"12"
    !!print*,"Integrale_Individuel=",Integrale_Individuel,"theta=",theta2
    !stop
   ! 100 continue
    return
    end function Integrale_Individuel
    
     ! 2 effets aleatoires correle niveau individeul +2 effet aleatoire correles niveau essai en interaction avec le traitement+ 2 effets aleatoires associes au risque de base:
     
    double precision function Integrale_Individuel_cor(vsi,vti,ui,uti,nnodes,ndim,j)
    ! vsi= frailtie niveau essai associe a s
    ! vti= frailtie niveau essai associe a t
    ! ui=  fragilites niveau essai base a S
    ! uti= fragilites niveau essai base a T
    ! nnodes: nombre de point de quadrature
    ! ndim= dimension de l'integrale 2 ou 1 integrations?
    ! j = individu j du cluster i
    ! npoint = nombre de point de quadrature
     use var_surrogate, only:posind_i,const_res5,const_res4,&
        deltastar,delta,varcovinv,varcov
    use comon, only: ve
    use var_surrogate, only: adaptative,xx1,ww1,posind_i,frailt_base
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    !$ use OMP_LIB
    !use mod_Adaptative ! pour lintegrant a maximiser pour l'adaptative
!    use func_adaptative
    !use Autres_fonctions
    IMPLICIT NONE
    integer,intent(in):: j,ndim,nnodes
    double precision,intent(in)::vsi,vti,ui,uti
    double precision, dimension(:,:),allocatable::m1,m3
    double precision, dimension(:,:),allocatable::m
    integer ::ii,jj,npg
    double precision::ss1,auxfunca,ss,wsij,wtij,test,rho
    double precision, dimension(ndim)::xxl !vecteur qui contiendra à chaque fois les points de quadrature
    double precision,dimension(ndim,ndim)::invBi_chol_Essai_k ! pour recuperer les matrice B_k dans le vecteur des matrices B des essais K

    !!print*,"xx=", xx,"xx1=",xx1
    ! Je calcule pour chaque individu du cluster par vectorisation son integrale 
    !!print*,"posind_i",posind_i,"j=",j
    !!print*,"delta(posind_i-1+j)",delta(posind_i-1+j)
    ! !print*,"ww1(1:npoint1)",ww1(1:npoint1)
    ! !print*,"xx(1:npoint1)",xx(1:npoint1)
    ! !print*,"xx(1:npoint1)",xx(1:npoint1)
    ! !print*,"delta(posind_i-1+j)",delta(posind_i-1+j)
    ! !print*,"deltastar(posind_i-1+j)",deltastar(posind_i-1+j)
    ! !print*,"const_res4(posind_i-1+j)",const_res4(posind_i-1+j)
    ! !print*,"ve(posind_i-1+j,1)",ve(posind_i-1+j,1)
    ! !print*,"const_res5(posind_i-1+j)",const_res5(posind_i-1+j)
    ! !print*,"vsi",vsi
    ! !print*,"vti",vti
    allocate(m(1,1))    
    allocate(m1(1,2),m3(1,2))
    ! !print*,"matrice varcovinv=",varcovinv
    ! !print*,"matrice m1=",m1
    ! !print*,"matrice m3=",m3
    ! !print*,"matrice m=",m
    
    !Partie a adapter au cas des individus
    if(frailt_base==0) then
        if(adaptative) then
        ! je recupere la matrice B_i de l'individu i dans le vecteur des B_i, dans ce cas il s'agit de la matrice associe aux individus et plus aux essais
            ! cpt=((i-1)*ndim**2)+1 !ceci permet de parcourir le vecteur invBi_chol_Essai en fonction du cluster sur lequel on se trouve
            do jj=1,ndim     
                do ii=1,ndim
                    ! invBi_chol_Essai_k(ii,jj)=invBi_chol_Essai(cpt) ! en effet le vecteur "invBi_chol_Essai" a ete rempli par des matrices colonne apres colonne
                    ! cpt=cpt+1

                enddo
            enddo        
            
            !!print*,"invBi_chol_Essai=",invBi_chol_Essai
            !!print*,"invBi_chol_Essai_k=",invBi_chol_Essai_k
            
        endif
    endif
    
    npg=nnodes
    auxfunca=0.d0
    ss=0.d0
    ss1=0.d0
    !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,xxl,wsij,wtij,m1,m3,m) firstprivate(ss1,auxfunca,vsi,vti,ui,uti)& 
    !$OMP SHARED(npg,varcovinv,nnodes,xx1,ww1,invBi_chol_Essai_k,ndim,adaptative,posind_i,test,rho,varcov,&
    !$OMP delta,j,deltastar,const_res4,ve,const_res5)&
    !$OMPREDUCTION(+:ss) SCHEDULE(Dynamic,1)
        do ii=1,npg
            ss1=0.d0
            do jj=1,npg
                xxl(1)=xx1(ii)
                xxl(2)=xx1(jj)
                !changement de variable en cas de quadrature adaptative
                if(adaptative) then ! on effectue le changement de variable
                    ! m2=matmul(invBi_chol_Essai_k,xxl)                            
                    ! xxl=ui_chap_Essai(i,1:2)+dsqrt(2.d0)*m2    
                !!print*,"xxl",xxl
                end if
                wsij=xxl(1)
                wtij=xxl(2)
                m1(1,1)=wsij
                m1(1,2)=wtij
                m3 = matmul(m1,varcovinv) !produit matriciel, resultat dans m3
                m = matmul(m3,transpose(m1)) !produit matriciel, resultat dans m
                
                !pour verification:
                ! rho=varcov(1,2)/(sqrt(varcov(1,1)*varcov(2,2)))
                ! test=-1.d0/(2.d0*(1.d0-rho**2))*(wsij**2/varcov(1,1)+wtij**2/varcov(2,2)&
                     ! -2.d0*wsij*wtij*rho/(sqrt(varcov(1,1)*varcov(2,2))))
                ! !print*,"expression calculee vaut:", test
                ! !print*,"expression par fonction vaut:", -m(1,1)/2.d0
                ! stop
                
                !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                auxfunca=dexp(wsij*delta(posind_i-1+j)+wtij*deltastar(posind_i-1+j)-m(1,1)/2.d0&
                            - const_res4(posind_i-1+j)*dexp(wsij+ui+vsi*dble(ve(posind_i-1+j,1)))&
                            - const_res5(posind_i-1+j)*dexp(wtij+uti+vti*dble(ve(posind_i-1+j,1)))&
                )
                !!print*,"12"
                ss1 = ss1+ww1(jj)*(auxfunca)
                ! !print*,vsi,vti,ui,uti,xxl(1),xxl(2)
                ! !print*,"Integrant niveau Individuel=",auxfunca
                ! stop
            end do
            ss = ss+ww1(ii)*ss1
        end do
    !$OMP END PARALLEL DO
    
    if(adaptative) then
        ! if(frailt_base==0) ss=ss*2.d0**(dble(ndim)/2.d0)*invBi_cholDet_Essai(i)  ! en admettant que "invBi_cholDet_Essai" c'est bien pour les individus
    end if
                
    Integrale_Individuel_cor=ss
    deallocate(m,m1,m3)
    !!print*,"12"
    !!print*,"Integrale_Individuel=",Integrale_Individuel,"theta=",theta2
    !stop
    return
    end function Integrale_Individuel_cor
    
    
    ! cas un effet aleatoires partages niveau individuel: integrant pour un individu, Monte carlo
    double precision function integrant_indiv_1MC(frailij,j)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! j = individu j du cluster i
     use var_surrogate
    use comon, only: eta

    IMPLICIT NONE
    integer,intent(in):: j
    double precision,intent(in)::frailij
    double precision:: c1,c2,c3,c4

    ! fin declaration et debut de la fonction
    !!print*,"j=",j,"posind_i=",posind_i-1,"posind_i+j=",posind_i-1+j
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c1=0.d0
    c2=frailij*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)
    c3=const_res4(posind_i-1+j)*dexp(frailij)
    c4=const_res5(posind_i-1+j)*dexp(eta*frailij)
    integrant_indiv_1MC=dexp(c2-c3-c4)
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    end function integrant_indiv_1MC
    
    ! cas un effet aleatoires partages niveau individuel: integrant pour toutes les realisations simulees de Monte carlo d'un individu
    double precision function integrant_indiv_1MCA(nsimu,j,mu1,vc1)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! nsimu= nombre de simulation
    ! j = individu j du cluster i
    ! mu1= moyenne du frailty
    ! vc1= variance du frailty
     use var_surrogate
    use comon, only: eta

    IMPLICIT NONE
    integer,intent(in):: nsimu,j
    double precision, intent(in)::mu1,vc1
    double precision,dimension(1:nsimu)::frailij,usim
    double precision:: c1,c2,c3,c4
    integer::n

    ! fin declaration et debut de la fonction
    !!print*,"j=",j,"posind_i=",posind_i-1,"posind_i+j=",posind_i-1+j
    n=size(frailij) ! qui correspond au nombre de simulation
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c1=0.d0
    ! on utilise la vectorisation pour reduire le temps de calcul
!    call cpu_time(c3)
    usim=Vect_sim_MC(1:nsimu,1) ! on lui affecte a chaque fois le premier vecteur des donnees simulees
    frailij=mu1+vc1*usim
    integrant_indiv_1MCA=SUM(dexp(&
        frailij(1:n)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)&
        -const_res4(posind_i-1+j)*dexp(frailij(1:n))&
        -const_res5(posind_i-1+j)*dexp(eta*frailij(1:n))&
    ))/dble(nsimu)
!    call cpu_time(c4)
    
!    !print*,"t1=",c4-c3
    
!    call cpu_time(c3)
!    do ttt=1,nsim
!      c2=c2+dexp(&
!        frailij(ttt)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)&
!        -const_res4(posind_i-1+j)*dexp(frailij(ttt))&
!        -const_res5(posind_i-1+j)*dexp(eta*frailij(ttt))&
!        )
!    end do
!    call cpu_time(c4)
    
!    !print*,"t2=",c4-c3
!    stop
    
    end function integrant_indiv_1MCA
    
    !==============Integration par Monte carlo=======================
    
     ! 1 effet aleatoire partage niveau individeul +2 effet aleatoire correles en interaction avec le traitement:
     ! integrale au niveau individuel
 double precision function Integrale_Individuel_MC(vsi,vti,ui,j,nsimu,mu1,vc1)
    ! vsi= frailtie niveau essai associe a s
    ! vti= frailtie niveau essai associe a t
    ! j = individu j du cluster i
    ! nsimu = nombre de boucle MC
    ! mu1= moyenne du frailty
    ! vc1= variance du frailty
     use var_surrogate
    use comon, only: eta,ve
    IMPLICIT NONE
    integer,intent(in):: j,nsimu
    double precision,intent(in)::vsi,vti,ui
    double precision, intent(in)::mu1,vc1
    double precision,dimension(1:nsimu)::frailij,usim
    integer::n
    
    ! Je calcule pour chaque individu du cluster par vectorisation son integrale par gaussHermite non adaptative
    n=size(frailij) ! qui correspond au nombre de simulation
    ! on utilise la vectorisation pour reduire le temps de calcul
!    call cpu_time(c3)
    usim=Vect_sim_MC(1:nsimu,1) ! on lui affecte a chaque fois le premier vecteur des donnees simulees
    frailij=mu1+vc1*usim
    if(frailt_base==0)then
        Integrale_Individuel_MC=SUM(dexp(&
            (vsi*delta(posind_i-1+j)+vti*deltastar(posind_i-1+j))*dble(ve(posind_i-1+j,1))+&
            frailij(1:n)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)&
            -const_res4(posind_i-1+j)*dexp(frailij(1:n)+vsi*dble(ve(posind_i-1+j,1)))&
            -const_res5(posind_i-1+j)*dexp(eta*frailij(1:n)+vti*dble(ve(posind_i-1+j,1)))&
        ))/dble(nsimu)
    else
        Integrale_Individuel_MC=SUM(dexp(&
            ui*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*alpha_ui)&
            +(vsi*delta(posind_i-1+j)+vti*deltastar(posind_i-1+j))*dble(ve(posind_i-1+j,1))+&
            frailij(1:n)*(delta(posind_i-1+j)+deltastar(posind_i-1+j)*eta)&
            -const_res4(posind_i-1+j)*dexp(frailij(1:n)+ui+vsi*dble(ve(posind_i-1+j,1)))&
            -const_res5(posind_i-1+j)*dexp(eta*frailij(1:n)+alpha_ui*ui+vti*dble(ve(posind_i-1+j,1)))&
        ))/dble(nsimu)
    endif
    return
    end function Integrale_Individuel_MC
    
    ! 2 effet aleatoire correlee au niveau individeul +2 effet aleatoire correles en interaction avec le traitement + 2 effets aleatoires risque de base:
     ! integrale au niveau individuel
    double precision function Integrale_Individuel_MC_cor(vsi,vti,ui,uti,j,nsimu,ndim,mu1,frailij)
        ! vsi= frailtie niveau essai associe a s
        ! vti= frailtie niveau essai associe a t
        ! j = individu j du cluster i
        ! nsimu = nombre de boucle MC
        ! ndim= dimension de l'integrale 2 ou 1 integrations?
        ! mu1= moyenne du frailty
        ! ui=  fragilites niveau essai base a S
        ! uti= fragilites niveau essai base a T
        ! frailij= frailty au niveau individuel generes suivant la gaussienne multivariee centree
        use var_surrogate, only: delta,deltastar,const_res4,const_res5,frailt_base,posind_i
        use comon, only: ve
        IMPLICIT NONE
        integer,intent(in):: j,nsimu,ndim
        double precision,intent(in)::vsi,vti,ui,uti
        double precision, intent(in)::mu1
        double precision,dimension(nsimu,ndim),intent(in)::frailij
        !double precision,dimension(:,:),allocatable::usim
        integer::n !l
        
        
        ! Je calcule pour chaque individu du cluster par vectorisation son integrale par gaussHermite non adaptative
        if(ndim==2)then
            !allocate(usim(1:nsimu,ndim))
            
            n=nsimu ! qui correspond au nombre de simulation
            ! on utilise la vectorisation pour reduire le temps de calcul
            ! call cpu_time(c3)
            ! usim=Vect_sim_MC(1:nsimu,1:2) ! on lui affecte a chaque fois le premier vecteur des donnees simulees
            ! frailij=mu1+vc1*usim
            ! l=1
            ! do while(l.le.nsimu)
                ! frailij(l,:)=0.d0+MATMUL(vc1,Vect_sim_MC(l,1:2)) 
                ! l=l+1
            ! end do
            
            if(frailt_base==0)then
                Integrale_Individuel_MC_cor=SUM(dexp(&
                    frailij(1:n,1)*delta(posind_i-1+j)+frailij(1:n,2)*deltastar(posind_i-1+j)&
                    -const_res4(posind_i-1+j)*dexp(frailij(1:n,1)+vsi*dble(ve(posind_i-1+j,1)))&
                    -const_res5(posind_i-1+j)*dexp(frailij(1:n,2)+vti*dble(ve(posind_i-1+j,1)))&
                ))/dble(nsimu)
            else
                Integrale_Individuel_MC_cor=SUM(dexp(&
                    frailij(1:n,1)*delta(posind_i-1+j)+frailij(1:n,2)*deltastar(posind_i-1+j)&
                    -const_res4(posind_i-1+j)*dexp(frailij(1:n,1)+ui+vsi*dble(ve(posind_i-1+j,1)))&
                    -const_res5(posind_i-1+j)*dexp(frailij(1:n,2)+uti+vti*dble(ve(posind_i-1+j,1)))&
                ))/dble(nsimu)
            endif
        endif
        !deallocate(usim)
        return
    end function Integrale_Individuel_MC_cor

    
     !gauss hermite pour une dimension: viens de Agnieszca
 
    SUBROUTINE gauherJ_scl(func,ss,nnodes,position_i)
        ! ss = resultat de l'integrale
        ! nnodes = nombre de noeudds d'integration
        ! position_i = la position j de l'individu sur lequel on se trouve dans le cluster i
        use donnees
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::nnodes,position_i
        double precision::auxfunca
        integer::j,methodGH
        double precision,dimension(nnodes):: xx1,ww1
        
        interface
          double precision function func30(arg,ndim)
          integer,intent(in)::ndim
          double precision,dimension(ndim), intent(in):: arg
          end function func30
          
          double precision function func(arg,i)
            integer, intent(in)::i! i est la position de l'individu sur lequel on se trouve
            double precision, intent(in):: arg
          end function func
        end interface
        
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
    
        ss=0.d0
        !!print*,"position_i quadrarure=",position_i
            auxfunca = 0.d0
            methodGH=0
            if(methodGH.eq.0) then !non-adaptative
                do j=1,nnodes  
                    auxfunca=func(xx1(j),position_i)
                    ss = ss+ww1(j)*(auxfunca)                     
                end do
            endif    
        return
    
    END SUBROUTINE gauherJ_scl
    
    !gauss Laguerre pour une dimension: viens de FRAILTYPACK
 
    SUBROUTINE gaulagJ_scl(func,ss,nnodes,position_i)
        ! ss = resultat de l'integrale
        ! nnodes = nombre de noeudds d'integration
        ! position_i = la position j de l'individu sur lequel on se trouve dans le cluster i
        use donnees
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::nnodes,position_i
        double precision::auxfunca
        integer::j
        double precision,dimension(nnodes):: xx1,ww1
        
        interface
          double precision function func30(arg,ndim)
          integer,intent(in)::ndim
          double precision,dimension(ndim), intent(in):: arg
          end function func30
          
          double precision function func(arg,i)
            integer, intent(in)::i! i est la position de l'individu sur lequel on se trouve
            double precision, intent(in):: arg
          end function func
        end interface
        
            if(nnodes.eq.20) then
                xx1(1:nnodes) = x(1:nnodes)
                ww1(1:nnodes) = w(1:nnodes)
            else if (nnodes.eq.32) then
                xx1(1:nnodes) = x1(1:nnodes)
                ww1(1:nnodes) = w1(1:nnodes)
                else 
                    !print*,"nous ne considErons que 20 ou 32 points de quadrature pour l'integration pas", &
                 !   "par quadrature guass Laguerre"
                    ! stop
            end if
    
        ss=0.d0
        !!print*,"position_i quadrarure=",position_i
            auxfunca = 0.d0
            do j=1,nnodes  
                auxfunca=func(xx1(j),position_i)
                ss = ss+ww1(j)*(auxfunca)                 
            end do   
        return
    
    END SUBROUTINE gaulagJ_scl
 
    
    ! cas deux effets aleatoires partages niveau individuel et essai
    double precision function funcSurrNN_MC_Essai_indiv(frail,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! i = position du cluster sur lequel on se trouve
     use var_surrogate
    use comon, only: ve,vedc,eta,alpha

    IMPLICIT NONE
    integer ::n,dimInt 
    integer,intent(in):: i
    !integer,intent(in)::n
    double precision:: C1,c2,c3,c4,vs
    double precision,intent(in),dimension(:)::frail

    ! fin declaration et debut de la fonction
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c1=0.d0
    dimInt=size(frail)
    n=nsujeti(i)
    vs=frail(dimInt)
    ! calculs
    c1=SUM(frail(1:(dimInt-1))*(delta(posind_i:(posind_i+n-1))+deltastar(posind_i:(posind_i+n-1))*eta))
    c2=(nigts(i)+alpha*cdcts(i))*vs
    c3=SUM(const_res4(posind_i:(posind_i+n-1))*dexp(frail(1:(dimInt-1))+vs*dble(ve(posind_i:(posind_i+n-1),1))))
    c4=SUM(const_res5(posind_i:(posind_i+n-1))*dexp(eta*frail(1:(dimInt-1))+alpha*vs*dble(vedc(posind_i:(posind_i+n-1),1))))

    funcSurrNN_MC_Essai_indiv=dexp(c1+c2-c3-c4)
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    return
    end function funcSurrNN_MC_Essai_indiv
    
    ! cas deux effets aleatoires correles
    double precision function funcSurrNN_MC_Essai_t2(frail,n,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! n= taille du vecteur frail
    ! i = position du premier individu du cluster, sera increentee de la taille du cluster + 1 a la sortie de la fonction
     use var_surrogate
    use comon, only: ve,vedc

    IMPLICIT NONE
    integer ::j,k !n=taille du tableau frail
    integer,intent(in):: i
    integer,intent(in)::n
    double precision:: s1,C1,c2,c3,c4,c5,vs,vt
    double precision,intent(in),dimension(n)::frail

    ! fin declaration et debut de la fonction
    s1=0.d0
    C1=0.d0
    c2=0.d0
    c3=0.d0
    c4=0.d0
    c5=0.d0
    vs=0.d0
    vt=0.d0
    s1=0.d0
    j=0
    k=0
    vs=frail(1)
    vt=frail(2)

    k=1
    !!print*,"i=",i,"n=",n,"i+n-2-1=",i+n-2-1
    do j=i,(i+n-1) ! nombre d'individu du cluster
      !c3=frail(k)*(delta(j)+deltastar(j)*eta)
      c2=(vs*delta(j)+vt*deltastar(j))*ve(j,1)
      c4=const_res4(j)*dexp(vs*ve(j,1))
      c5=const_res5(j)*dexp(vt*vedc(j,1))
      !c4=const_res4(j)*dexp(vs*ve(j,1))*dexp(frail(k))
      !c5=const_res5(j)*dexp(vt*vedc(j,1))*dexp(frail(k)*eta)
      s1=s1+c2-c4-c5
      !k=k+1
!      !print*,"c3=",c3,"c2=",c2,"c4=",c4,"c=5",c5,"s1=",s1
    end do
    !!print*,"k=",k
    !stop
    !!print*,"integrant ligne 58"
!   end if
    funcSurrNN_MC_Essai_t2=dexp(s1)
!    !print*,"funcSurrNN_MC ligne 47: je vaut",funcSurrNN_MC
    return
    end function funcSurrNN_MC_Essai_t2
    
    !====Integrant cas MC+quadrature
    
    double precision function funcSurrNN1(frail,frailst,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! frailst= vecteur des fragilites de taille 2 vs et vt
    ! i = position du premier individu du cluster, sera increentee de la taille du cluster + 1 a la sortie de la fonction
    !use tailles
    !use comon,only:nig,auxig,alpha,theta,res1,res3,aux1,cdc
    use var_surrogate
    use comon, only: ve,eta

    IMPLICIT NONE
    integer ::j,n,k !n=taille du tableau frail
    integer,intent(inout):: i
    double precision:: s1,C1,c2,c3,c4,c5,vs,vt
    double precision,intent(in),dimension(:)::frail
    double precision,intent(in),dimension(2)::frailst !pour les effets aleatoire vs et vt
    double precision, dimension(1,2)::m1,m3
    double precision, dimension(1,1)::m

    ! fin declaration et debut de la fonction
    
    n=size(frail)
    s1=0.d0
    !calcul de c1
    vs=frailst(1)
    vt=frailst(2)
    m1(1,1)=vs
    m1(1,2)=vt
    call multiJ(m1,varcovinv,1,2,2,m3) !produit matriciel, resultat dans m3
    call multiJ(m3,transpose(m1),1,2,1,m) !produit matriciel, resultat dans m
    c1=(1.d0/(2.d0*pi*dsqrt(determinant)))*dexp((-1.d0/2.d0) * m(1,1)) ! attention qd deteminant vaut 0
    !!print*,"vs=",vs,"vt=",vt,"c1=",c1,"varcovinv=",varcovinv,"determinant=",determinant,"pi=",pi
    !!print*,"m3=",m3,"m=",m
    !stop
    !!print*,"******je suis dans la subroutine funcSurrNN1 ligne 103: c1=",c1
    
    ! calcul des elements dans la somme
    k=1
    do j=i,(i+n) ! n=nombre d'individu du cluster
      c3=frail(k)*(delta(j)+deltastar(j)*eta)
      c2=(vs*delta(j)+vt*deltastar(j))*ve(j,1)
      c4=const_res4(j)*dexp(vs*ve(j,1))*dexp(frail(k))
      c5=const_res5(j)*dexp(vt*ve(j,1))*dexp(frail(k)*eta)
      s1=s1+c3+c2-c4-c5
      k=k+1
    end do
    
    !i=j+1 ! prochain cluster (on le met à jour dans funcpajspline_scl)
    
    funcSurrNN1=c1*dexp(s1)
    !!print*,"funcSurrNN1 ligne 57: je vaut",funcSurrNN1
    return
    end function funcSurrNN1
    
    !==================================
    ! func pour monte carlo
    !==================================
    
    double precision function funcSurrNNMC(frail,ndim,i)
    ! dans cette fonction qui sera utilisee pour le monte carlo, on doit faire une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! frail= vecteur des effets aleatoire! 
    ! ndim= dimension de l'integrale ou nombre d'individu du cluster
    ! i position du premier individu du cluster dans lequel on integre
    
    use var_surrogate, only: npoint
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    
    implicit none
    integer ::nnodes,k
    integer,intent(in):: ndim
    integer, intent(inout)::i !position de l'individu
    double precision,dimension(ndim), intent(in):: frail
    double precision,dimension(2):: frailst
    double precision,dimension(:),allocatable ::xx1,ww1
    double precision::resultat,inc
    !double precision, external::gaussHermMultMC    
    
    allocate(xx1(npoint))
    allocate(ww1(npoint))
    
    nnodes=npoint
    !!print*,"npoint=",npoint
    
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
    end if
    
    !integration sur vsi et vti
    !adaptative=.false.
    frailst=xx1(1)
    k=2
    inc=0.d0
    !!print*,"je suis dans la subroutine funcSurrNNMC l=112: frail(1:2)=",frail(1:2),"frailst",frailst,"i=",i,"k=",k,"xx1(1:2)=",xx1(1:2)
    !!print*,"ww1(1:2)=",ww1(1:2),"inc=",inc
    !resultat=gaussHermMultMC(funcSurrNN1,frail,frailst,i,k,xx1,ww1,inc)
    !!print*,"******je suis dans la subroutine funcSurrNNMC ligne 117: pour i=",i
    resultat=gaussHermMultMC(frail,frailst,i,k,xx1,ww1,inc)
    !!print*,"******je suis dans la subroutine funcSurrNNMC ligne 119:, resultat=",resultat,"pour i=",i
    funcSurrNNMC=resultat
    
    deallocate(xx1)
    deallocate(ww1)
    return
  end function funcSurrNNMC
  
!==================================================================
!=======fonction gauss hermite multidimension=======
!==================================================================
  
  ! funcSurrNN est l'integrant, herm le resultat de l integrale sur -infty , +infty
! l'appel de cette fonction se fait en initialisant le vecteur frail a x(1)
! ie do j=1,size(frail)
!        frail[j]=x(1)
!     end do
! initialiser k a size(frail)   

recursive function gaussHermMultMC(frail1,frail,i,k,x,w,inc) result(herm)
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! frail1: vecteur donc la taille correspond au nombre d'integrale pour l'integrale multiple
   ! frail: vecteur de taille 2 pour les integrale sur vsi et vti
   ! i: position du premier individu du cluster dans lequel on integre
   ! k permet de gerer la recursivite et vaut initialement la taille de frail
   ! x et w: vecteur des poids et points de quadrature definis dans le fichier Adonnees.f90
   ! inc un increment pour le controle, vaut 0 initialement
   
   use var_surrogate, only:adaptative
   !use comon, only:invBi_cholDet
   
   implicit none
   
   integer, intent(in)::k
   integer,intent(inout):: i
   double precision, intent(inout)::inc
   integer :: l1,lk,n,npoint
   double precision,dimension(:),intent(in) :: frail1
   double precision,dimension(:),intent(inout) :: frail
   double precision,dimension(:),intent(in) ::x,w
   double precision ::herm,s
   
   ! fin declaration et debut du programme
   !!print*,"******je suis dans la subroutine gaussHermMultMC fichier integrant_scl ligne 157:"
   
   n=size(frail)
   npoint=size(x)
   if (k.eq.1) then
     s=0.d0
     do l1=1,npoint
       frail(n)=x(l1)
       s=s+w(l1)*funcSurrNN1(frail1,frail,i)
       !!print*,"s=",s
       inc=inc+1.d0
     end do
     herm=s
     !!print*,"ligne 175: herm=",s
   else
     s=0.d0
     do lk=1,npoint
       frail(n-k+1)=x(lk)
       s=s+w(lk)*gaussHermMultMC(frail1,frail,i,k-1,x,w,inc)
       !!print*,"ligne 181: w(lk)=",w(lk),"i=",i,"s=",s
     end do
     herm=s
     !!print*,"ligne 184: herm=",s,"i=",i,"npoint=",npoint
   end if
   ! ici on complete avec le produit matrice de cholesky de 2^r/2 de l'adaptative
   if(adaptative .and. inc .eq.(npoint**n)) then
     !print*," sub routine gaussHermMultMC, je suis dans le if fichier integrant_scl.f90"
     ! stop
   end if
   !!print*,"ligne 191: herm=",herm
 end function gaussHermMultMC  
 
 ! integrant joint classique: 2 effet aleatoire correles en interaction avec le traitement
 double precision function funcSurrNN_Essai_2t(frail,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai normalement distribues
    ! frail= vecteur des frailties de taille 2
    ! i = position du cluster sur lequel on se trouve

    use var_surrogate
    use comon, only: ve,vedc

    IMPLICIT NONE
    integer ::n
    integer,intent(in):: i
    double precision:: C1,c2,c3,c4,vs,vt
    double precision,dimension(2),intent(in)::frail
    double precision, dimension(1,2)::m1,m3
    double precision, dimension(1,1)::m
    

    vs=frail(1)
    vt=frail(2)
    m1(1,1)=vs
    m1(1,2)=vt
    m3=MATMUL(m1,varcovinv)
    m=MATMUL(m3,TRANSPOSE(m1))
    !!print*,varcov
    !!print*,varcovinv
    !!print*,determinant
    
    c1=-m(1,1)/2.d0
    !!print*,m(1,1)
    !!print*,m

    n=nsujeti(i)
    c2=nigts(i)*vs+cdcts(i)*vt
    c3=SUM(const_res4(posind_i:(posind_i+n-1))*dexp(vs*dble(ve(posind_i:(posind_i+n-1),1))))
    c4=SUM(const_res5(posind_i:(posind_i+n-1))*dexp(vt*dble(vedc(posind_i:(posind_i+n-1),1))))
    !!print*,"posind_i",posind_i
    funcSurrNN_Essai_2t=dexp(c1+c2-c3-c4)
    return
    end function funcSurrNN_Essai_2t
    
 ! integrant joint classique: 1 effet aleatoire partage en interaction avec le traitement
 double precision function funcSurrNN_Essai_t(frail,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai normalement distribues
    ! frail= vecteur des frailties de taille 2
    ! i = position du cluster sur lequel on se trouve

    use var_surrogate
    use comon, only: ve,vedc,alpha !res1

    IMPLICIT NONE
    integer ::n
    integer,intent(in):: i
    double precision:: C1,c2,c3,c4,vs
    double precision,intent(in)::frail
    !double precision, dimension(1,2)::m1,m3
    
    n=nsujeti(i)
    vs=frail
    c1=-vs**2/(2.d0*sigma2)
    c2=(nigts(i)+alpha*cdcts(i))*vs
    c3=SUM(const_res4(posind_i:(posind_i+n-1))*dexp(vs*dble(ve(posind_i:(posind_i+n-1),1))))
    c4=SUM(const_res5(posind_i:(posind_i+n-1))*dexp(alpha*vs*dble(vedc(posind_i:(posind_i+n-1),1))))
    
    funcSurrNN_Essai_t=dexp(c1+c2-c3-c4)
    return
    end function funcSurrNN_Essai_t
    
    ! sans interaction avec le traitement
     double precision function funcSurrNN_Essai(frail,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai normalement distribues
    ! frail= vecteur des frailties de taille 2
    ! i = position du cluster sur lequel on se trouve
    !use tailles
    !use comon,only:nig,auxig,alpha,theta,res1,res3,aux1,cdc
    use var_surrogate
    use comon, only: alpha !vedc

    IMPLICIT NONE
    !integer ::j !n=taille du tableau frail
    integer,intent(in):: i
    double precision:: C1,c2,c3,c4,vs
    double precision,intent(in)::frail
    !double precision, dimension(1,2)::m1,m3

    ! fin declaration et debut de la fonction
    
    !n=size(frail)
    vs=frail
    c1=-vs**2/(2.d0*sigma2)
    c2=(nigs(i)+alpha*cdcs(i))*vs
    c3=const_res1(i)*dexp(vs)
    c4=const_aux1(i)*dexp(alpha*vs)
    
    funcSurrNN_Essai=dexp(c1+c2-c3-c4)
    return
    end function funcSurrNN_Essai

    
    !sauvegarde
    double precision function funcSurrNN(frail,i)
    ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues
    ! frail= vecteur des frailties de taille nombre d'individus dans le cluster ou on se trouve
    ! i = position du premier individu du cluster, sera increentee de la taille du cluster + 1 a la sortie de la fonction
    !use tailles
    !use comon,only:nig,auxig,alpha,theta,res1,res3,aux1,cdc
    use var_surrogate
    use comon, only: ve,eta,theta

    IMPLICIT NONE
    integer ::j,n !n=taille du tableau frail
    integer,intent(inout):: i
    double precision:: s1,C1,c2,c3,c4,c5,vs,vt
    double precision,intent(in),dimension(:)::frail
    double precision, dimension(1,2)::m1,m3
    double precision, dimension(:,:), allocatable::m

    ! fin declaration et debut de la fonction
    
    allocate(m(1,1))
    n=size(frail)
    s1=0
    !calcul de c1
    vs=frail(n-1)
    vt=frail(n)
    m1(1,1)=vs
    m1(1,2)=vt
    call multiJ(m1,varcovinv,1,2,2,m3) !produit matriciel, resultat dans m3
    call multiJ(m3,transpose(m1),1,2,1,m) !produit matriciel, resultat dans m
    c1=(-1.d0/2.d0) * m(1,1)
    
    ! calcul des elements dans la somme
    do j=i,(n-2) ! on s'arrete a n-2 car les w_ij sarrete a n-2. a n-1 on a vs et a n o na vt
      c3=delta(j)+deltastar(j)*eta
      c2=vs*delta(j)+vt*deltastar(j)*ve(j,1)
      c4=const_res4(j)*dexp(vs*ve(j,1))
      c5=const_res5(j)*dexp(vt*ve(j,1))
      s1=s1+frail(j)*c3+c2-c4*exp(frail(j))-c5*exp(frail(j))-(1.d0/(2.d0*theta))*frail(j)**2.d0
    end do
    
    !i=j+1 ! prochain cluster
    
    funcSurrNN=dexp(c1+s1)
    return
    end function funcSurrNN
  
  !==========fonction a integrer exp(-x^2-y^2-z^2...)
    
  double precision function func(arg,ndim)
    use var_surrogate, only:adaptative,mui 
    use comon, only:invBi_chol
    
    implicit none
    integer ::i
    integer, intent(in)::ndim
    double precision ::s,s2
    double precision,dimension(ndim), intent(in):: arg
    double precision,dimension(ndim):: frail
    double precision,parameter::pi=3.141592653589793d0
    
    if(adaptative) then
      frail = mui+sqrt(2.d0)*MATMUL(invBi_chol,arg)
    else
      frail = arg
    end if
            
    s=0.d0
    s2=0.d0
    do i=1,ndim
      s=s-(frail(i)**2.d0)
      s2=S2-frail(i)**2.d0
    end do
    func=dexp(s2)*((1.d0/(2.d0*pi))**(ndim/2.d0))*dexp(s/2.d0)
  end function func
  
  !1 integrale
  double precision function func1(arg,ndim)
    implicit none
    !integer ::i
    integer, intent(in)::ndim
    double precision ::s
    double precision,dimension(ndim), intent(in):: arg
    double precision,parameter::pi=3.141592653589793d0
    
    s=dexp(-arg(1)**2.d0)*1.d0/(2.d0*pi)**(1.d0/2.d0)*dexp(-arg(1)**2.d0/2.d0)
    func1=s
  end function func1
  
  !func_x2y2
  double precision function func_x2y2(arg,ndim)
    implicit none
    integer ::i
    integer,intent(in)::ndim
    double precision ::s,s2
    double precision,dimension(2), intent(in):: arg
    double precision,parameter::pi=3.141592653589793d0
    s=0.d0
    s2=0.d0
    do i=1,ndim
      s=s-(arg(i)**2.d0)
      s2=S2-arg(i)
    end do
    func_x2y2=dexp(s2)*((1.d0/(2.d0*pi))**(ndim/2.d0))*dexp(s/2.d0)
  end function func_x2y2
  
  !func 2
  double precision function func2(arg,ndim)
    implicit none
    integer ::i
    integer, intent(in)::ndim
    double precision ::s,s2
    double precision,dimension(ndim), intent(in):: arg
    double precision,parameter::pi=3.141592653589793d0
    s=0.d0
    s2=1.d0
    do i=1,ndim
      s=s-arg(i)**2.d0
      s2=S2*(arg(i))
    end do
    func2=s2*(1.d0/(2.d0*pi)**(ndim/2.d0))*dexp((1.d0/2.d0)*s)
  end function func2
  
    !==========fonction a integrer gaussien exp(-x^2-y^2-z^2...)
    
  double precision function funcMC(arg,ndim)
    implicit none
    integer ::i
    integer, intent(in)::ndim
    double precision ::s
    double precision,dimension(ndim), intent(in):: arg
    s=0.d0
    do i=1,ndim
      s=s-arg(i)**2.d0
    end do
    funcMC=dexp(s)
  end function funcMC
  
  !deuxieme
  double precision function funcMC2(arg,ndim)
    implicit none
    integer ::i
    integer, intent(in)::ndim
    double precision ::s
    double precision,dimension(ndim), intent(in):: arg
    s=1.d0
    do i=1,ndim
      s=s*arg(i)
    end do
    funcMC2=s
  end function funcMC2
  
  !================== multiplication de matrice  ==================

! multiplie A par B avec le resultat dans C

    subroutine multiJ(A,B,IrowA,JcolA,JcolB,C)
!     remarque :  jcolA=IrowB
    use tailles
    
    IMPLICIT NONE
    
    integer,intent(in)::IrowA,JcolA,JcolB
    double precision,dimension(IrowA,JcolA),intent(in):: A
    double precision,dimension(JcolA,JcolB),intent(in):: B
    double precision,dimension(IrowA,JcolB),intent(out)::C       
    integer::i,j,k
    double precision::sum
    
    !!print*,"A=",A
    !!print*,"B=",B
    do I=1,IrowA
        do J=1,JcolB
            sum=0
            do K=1,JcolA
                sum=sum+A(I,K)*B(k,J)
            end do
            C(I,J)=sum
            !!print*,"sum=",sum
        end do
    end do
    
    return
    
    end subroutine multiJ
    
    ! integrant pour le calcul du tau de kendall par quadrature de Gauss Hermite 
    
    double precision function funcJointSurroKendall(u, w, up, wp, theta, gamma, alpha, eta, ui)
        
        use var_surrogate, only:pi
        implicit none
        
        double precision, intent(in):: u, w, up, wp, theta, gamma, alpha, eta, ui
        double precision ::num,den
        
        if(ui==1) then
            num = (dexp(w + u + eta * w + alpha * u) + dexp(wp + up + eta * wp + alpha * up)) *&
                dexp(-(wp**2.d0)/(2.d0 * theta)) * dexp(-(up**2.d0)/(2.d0 * gamma)) * &
                dexp(-(w**2.d0)/(2.d0 * theta)) * dexp((-u**2.d0)/(2.d0 * gamma))
          
            den = (dexp(wp + up) + dexp(w + u)) * (dexp(eta * wp + alpha * up) + &
                dexp(eta * w + alpha * u)) * 4.d0 * pi**2.d0 * theta * gamma
        else
            num = (dexp(w + eta * w) + dexp(wp + eta * wp)) * dexp(-(wp**2.d0)/(2.d0 * theta))* &
                dexp(-(w**2.d0)/(2.d0 * theta))
          
            den = (dexp(wp) + dexp(w)) * (dexp(eta * wp ) + &
                dexp(eta * w)) * 2.d0 * pi * theta
        
        endif
        funcJointSurroKendall = num/den
  
  end function funcJointSurroKendall
  
end module fonction_A_integrer