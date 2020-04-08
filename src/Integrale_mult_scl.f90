!========================================================================
!=======module pour la fonction monteCarlosMult : calcul integral =======
!========================================================================

module Laplace_contribution

implicit none

contains
! ================Calcul de la contribution d'un essai a la log-vraisemblance par approximation de Laplace=====================
    
    
    double precision function Cont_Laplace_Essai(determin)
        !determin : le determinant de la matrice de variance-covariance effects aleatoires niveau essai
        !essaicourant: essai courant
        !posindi : poisition du sujet dans le jeu de donnee
        use var_surrogate, only:vs_i,vt_i,u_i, & !theta2,const_res5,const_res4
            pi,nsujeti,essai_courant,position_i,& !deltastar,delta,varcovinv
            nparamfrail,rho,varcov,gamma_ui,wij_chap,& !alpha_ui,res2s_sujet,res2_dcs_sujet
            Test
            
        use comon, only: eta !ve
        use Autres_fonctions,only:Determinant
        use optim_scl2, only:marq98j_scl2  ! pour faire appel a marquard 
        use func_laplace, only:funcpaXi_chapeau ! se traouve dans le fichier funcpa_laplace.f90 pour les autres fonction necessaires a laplace
        use func_laplace, only:Int_Laplace_ind
        !$ use OMP_LIB
        
        implicit none
         
        double precision, intent(in)::determin    
        !integer, intent(in)::essaicourant,posindi
        integer,parameter::effet2=0
        double precision::ca,cb,dd,zeta,& !k_second,h_second_ui,h_ui_vsi,h_ui_vti,h_second_vsi,h_vsi_vti,h_second_vti
                            jacobien,h2,h,ui,vsi,vti,res,B_Lap,control !h1,tp1,tp2
        double precision, dimension(2)::k0_2
        double precision, allocatable, dimension(:,:)::H_hess_scl,I_hess_scl,hess_scl
        double precision,dimension(:), allocatable::vvv_scl,v,b_2
        integer::ier,istop,ni,np_2,nparamfrail_save,i,non_conv !individu_j
        !double precision,dimension(3,3)::mat_J ! matrice jacobienne
        
        zeta=eta        
        !====================================================================================================
        ! ==============================estimation des X_i_chapeau==========================================
        !====================================================================================================    
        k0_2=100.d0 ! valeur arbitraire, ne sert a rien
        ni=0
        ca=0.d0
        cb=0.d0
        dd=0.d0
        Test=0                        
        
        
        !=================pour le test=============================
        if(Test==1)then ! je fais ceci pour evaluer le calcul integral par laplace
            np_2=2
            allocate(v(np_2*(np_2+3)/2))
            allocate(b_2(np_2))
            allocate(I_hess_scl(np_2,np_2))
            allocate(H_hess_scl(np_2,np_2),hess_scl(np_2,np_2),vvv_scl(np_2*(np_2+1)/2))
                            
            b_2=0.5d0
            v=0.d0 ! matrice des derrivees premieres et seconde
            nparamfrail_save=nparamfrail
            nparamfrail=2
            non_conv=0
            
            !!print*, model
            ! call marq98o2(b_2,np_2,ni,v,res,ier,istop,funcpaXi_chapeau)
            call marq98J_scl2(k0_2,b_2,np_2,ni,v,res,ier,istop,effet2,ca,cb,dd,funcpaXi_chapeau,I_hess_scl,H_hess_scl,&
                             hess_scl,vvv_scl)
                         
            ! calcul integrale J(x,y,z): posons k=10
            if(istop==1)then
                !print*,"I_hess_scl",I_hess_scl
                !print*,"H_hess_scl",H_hess_scl
                !print*,"hess_scl",hess_scl
                !print*,"vvv_scl",vvv_scl
                !print*,"v=",v
                !print*,"determinant hessienne",determinant(I_hess_scl,2)
                h=dexp(-wij_chap(1,1) + 5.d0*dlog(wij_chap(1,1)))*dsqrt(2*pi*wij_chap(1,1)**2.d0/5.d0)
                !print*,"x_0 pour vaut:",wij_chap(1,1)
                !print*,"la valeur de l'integrale pour  k*n!, avec k=10 et n=5 vaut",h*10
                !print*,"la valeur de l'integrale pour n! vaut",h
                !print*,"y_0 et z_0 valent:",b_2(1),b_2(2)
                !print*,"la valeur de l'integrale pour  apres appel de la function pour l'integrationn! vaut",&
                 !   Int_Laplace_ind(position_i,i,vs_i,vt_i,u_i)
                h=2.d0*pi*dexp(b_2(1)**2.d0 + 2.d0*b_2(2)+dlog(h))*1.d0/dsqrt(8.d0)
                !print*,"voila la valeur de int(exp(y^2 + 2z^2 + n!) pour n=5",h
                !print*,"===========Pour annuler ce calcul, attribuer la valeur 0 à la variable Test dans Cont_Laplace_Essai() "
                stop
            else
                h=dexp(-wij_chap(1,1) + 5.d0*dlog(wij_chap(1,1)))*dsqrt(2*pi*wij_chap(1,1)**2.d0/5.d0)
                !print*,"x_0 pour vaut:",wij_chap(1,1)
                !print*,"la valeur de l'integrale pour  k*n!, avec k=10 et n=5 vaut",h*10
                !print*,"la valeur de l'integrale pour n! vaut",h
                !print*,"la valeur de l'integrale pour  apres appel de la function pour l'integrationn! vaut",&
                 !   Int_Laplace_ind(position_i,i,vs_i,vt_i,u_i)
                !print*,"probleme de convergence essai"
                !print*,"===========Pour annuler ce calcul, attribuer la valeur 0 à la variable Test dans Cont_Laplace_Essai() "
                stop
            endif
            !goto 125
            
            ! stop
        endif
        !=================fin pour le test==========================
        
        ! !print*,"suis la 0",essai_courant,nsujeti(essai_courant),size(wij_chap),size(wij_chap,1),size(wij_chap,2)
        !stop
        np_2=3
        allocate(v(np_2*(np_2+3)/2))
        allocate(b_2(np_2))
        allocate(I_hess_scl(np_2,np_2))
        allocate(H_hess_scl(np_2,np_2),hess_scl(np_2,np_2),vvv_scl(np_2*(np_2+1)/2))
                        
        b_2=0.5d0
        v=0.d0 ! matrice des derrivees premieres et seconde
        nparamfrail_save=nparamfrail
        nparamfrail=3
        non_conv=0
        
        10 continue
        !!print*, model
        ! call marq98o2(b_2,np_2,ni,v,res,ier,istop,funcpaXi_chapeau)
        call marq98J_scl2(k0_2,b_2,np_2,ni,v,res,ier,istop,effet2,ca,cb,dd,funcpaXi_chapeau,I_hess_scl,H_hess_scl,&
                         hess_scl,vvv_scl)
        
             
        if (istop.ne.1 .and. non_conv<=10) then ! on passe à l'individu suivant, juste pour le test
            !!print*,"Fin estimation_b essai bloquer",b_2,essai_courant
            b_2=-0.5*non_conv
            non_conv=non_conv+1 !compte le nombre de fois qu'on n'a pas pu estime les frailties niveau essai sur certains individus
            !stop
            goto 10
        endif
                                
        if(non_conv==11 .and. istop .ne. 1)then
            ! !print*,"le nombre de tentative sans convergence vaut:",non_conv,"essai ig_=",essai_courant
            ! !print*,"istop=",istop,"essai k=",essai_courant
            non_conv=0
            Cont_Laplace_Essai =-1.d9
            goto 126
        endif
                                
        if(non_conv>0 .and. non_conv<=10) then
            !!print*,"le nombre de tentative pour la convergence vaut:",non_conv
            ! !print*,"istop=",istop,"essai k=",essai_courant
            non_conv=0
            !stop
        endif        
        
        nparamfrail=nparamfrail_save ! on restitu sa valeur avant de continuer
        ! model = model_save
        
        ! if(essai_courant==1)!print*,"Fin estimation_b essai",b_2,"essai",essai_courant
        ! stop
        !if(essai_courant==5) stop
        u_i=b_2(1)        
        vs_i=b_2(2)
        vt_i=b_2(3)

        ! ===========================fin estimation des X_i_chapeau===========================================================
        jacobien=determinant(I_hess_scl,3) ! determinant de la hesienne

        ! ===============calculer h(X_chapeau)===============================
        
        ! calcul des integrales au niveau individuel pour tous les sujets du cluster
        ui=u_i
        vsi=vs_i 
        vti=vt_i
        
        ! calcul de l'integrale au niveau individuel
            
            h2=0.d0
            control=0
    !!print*,"suis la 2"        
            !$OMP PARALLEL DO default(none) PRIVATE (i,h)& 
            !$OMP  shared(control,nsujeti,essai_courant,position_i,vs_i,vt_i,u_i) REDUCTION(+:h2)
                do i=1,nsujeti(essai_courant)
                    ! individu_j=position_i-1+i
                    h=Int_Laplace_ind(position_i,i,vs_i,vt_i,u_i)
                    if(h==-1.d9) then
                        control=1
                        !funcpaXi_chapeau =-1.d9
                        !goto 124
                    endif
                    h2=h2 + h
                enddo
            
            !$OMP END PARALLEL DO
    ! !print*,"suis la 3"        
    ! stop
            ! en cas de non convergence sur un individu
            if(control==1) then
                Cont_Laplace_Essai =-1.d9
                goto 126
            endif
    
        ! B_Lap=dlog(2.d0*pi)+(1.d0/2.d0)*(dlog(determin)+dlog(2.d0*pi*gamma_ui))
        B_Lap=0.d0
        
        h=  B_Lap  &
            + (u_i**2.d0)/(2.d0*gamma_ui) +1.d0/(2.d0*(1.d0-(rho**2.d0))) &
            *((vs_i**2.d0)/varcov(1,1) + (vt_i**2.d0)/varcov(2,2)&
            - (2.d0*vs_i*vt_i*rho)/dsqrt(varcov(1,1)*varcov(2,2)))&
            - h2
        
        if(jacobien<0) then 
            ! !print*,essai_courant,"jacobien",jacobien,"res",res
            ! !print*,essai_courant,"jacobien negatif",jacobien
            ! !print*,"jacobien 1:",mat_J(1,:)
            ! !print*,"jacobien 2:",mat_J(2,:)
            ! !print*,"jacobien 3:",mat_J(3,:)
            ! !print*,"determinant",Determinant(mat_J, 3),Determinant(mat_J, 3),Determinant(mat_J, 3)
            ! stop
            ! Cont_Laplace_Essai =-1.d9
            ! goto 126
        endif
        
        res=3.d0/2.d0 * dlog(2*pi) - h - (1.d0/2.d0)*dlog(abs(jacobien))
        
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        Cont_Laplace_Essai =-1.d9
        goto 126
    else
        Cont_Laplace_Essai = res
    end if
    
    
    126    continue
    
    ! if(essai_courant .ne. 1)then
        ! deallocate(wij_chap)
    ! endif
    deallocate(v,b_2)        
    deallocate(H_hess_scl,I_hess_scl,hess_scl,vvv_scl)
    
    return
    
    endfunction Cont_Laplace_Essai
    

end module laplace_contribution




module monteCarlosMult_Gaus

    implicit none
 
    contains
!========================================================================
!
!SUBROUTINES simulation pour le calcul integrale multiple par monte carlo
!
!========================================================================
double precision function MC_Copula_Essai(func,ndim,nsujet_trial,i)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! func: fonction a integrer au niveau individuel
    ! ndim= dimension de l'integrale 2 ou 3 integrations?
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    use Autres_fonctions, only:init_random_seed
    use var_surrogate, only: Vect_sim_MC,a_deja_simul,nsim,chol,frailt_base,&
                             graine,aleatoire,nbre_sim,nb_procs, control_affichage
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use Autres_fonctions, only:pos_proc_domaine
    !use mpi
    !$ use OMP_LIB
    
    implicit none
    integer,intent(in):: ndim,nsujet_trial,i
    integer ::ii,jj,l,m,maxmes,nsimu,init_i,max_i,code,erreur,rang
    double precision:: ss,SX,x22 
    double precision,dimension(:,:),allocatable::vc, fraili
    double precision,dimension(:),allocatable::usim,vi

    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,ig,nsujet_trial)
            ! vsi= frailtie niveau essai associe a s
            ! vti= frailtie niveau essai associe a t
            ! ui = random effect associated xith the baseline hazard
            ! ig = current cluster
            ! nsujet_trial = number of subjects in the current trial
            integer,intent(in):: ig, nsujet_trial
            double precision,intent(in)::vsi,vti,ui
        end function func
    end interface
    
    allocate(vc(ndim,ndim),fraili(nsim,ndim))
    !vc=ABS(chol)
    !vc = chol
    if(frailt_base==0)then
        vc = 0.d0 
        vc(1,1) = Chol(1,1)
        vc(2,1) = Chol(2,1)
        vc(2,2) = Chol(2,2)
    else
        vc = 0.d0 
        vc(1,1) = Chol(1,1)
        vc(2,1) = Chol(2,1)
        vc(2,2) = Chol(2,2)
        vc(3,3) = Chol(3,3)
    endif
    nsimu=nsim
    x22=0.d0

    maxmes=size(vc,2)
     allocate(vi(maxmes*(maxmes+1)/2))
     allocate(usim(maxmes))    
    
  ! --------------------- boucle du MC ------------------------
    l=1
    !stemp=0
    !===============================================================================
    ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
    !===============================================================================
    
    if(a_deja_simul.eq.0) then
        call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour lagraine
        Vect_sim_MC=0.d0
        do while(l.le.nsimu)
            ! on genere suivant des normales centrees reduites pour les variables aleatoires correlees
             usim=0.d0
             do m=1,maxmes
                 SX=1.d0
                 call bgos(SX,0,Vect_sim_MC(l,m),x22,0.d0)
             end do
            l=l+1
        end do    
        
        a_deja_simul=1 ! pour dire qu'on ne simule plus
    endif
    
    ! on utilise les generations precedentes pour obtenir deux variables correlees suivant une multinormale centree de covariance vc
    l=1
    do while(l.le.nsimu)
        if(frailt_base==0)then
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,1:2)) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
        else
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,1:3))
        endif
        l=l+1
    end do
    
    ! call intpr(" dans nb_procs=", -1, nb_procs, 1)
    ! call intpr(" dans ndim=", -1, ndim, 1)
    !integration sur vsi et vti
    ss=0.d0
    ! call OMP_SET_NUM_THREADS(1)
    if(nb_procs==1) then !on fait du open MP car un seul processus
        rang=0
        if(ndim.eq.2) then
            !$OMP PARALLEL DO default(none) PRIVATE (ii) SHARED(nsimu,nsujet_trial,i,fraili)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
                do ii=1,nsimu
                    ss=ss+func(fraili(ii,1),fraili(ii,2),0.d0,i,nsujet_trial)
                    !!print*,"ss",ss
                end do
            !$OMP END PARALLEL DO
        else ! cas de 3 points
           !$OMP PARALLEL DO default(none) PRIVATE (ii) SHARED(nsimu,nsujet_trial,i,fraili)&
           !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
                do ii=1,nsimu
                    ss=ss+func(fraili(ii,1),fraili(ii,2),fraili(ii,3),i,nsujet_trial)
                    ! call dblepr(" dans ss=", -1, ss, 1)
                end do
           !$OMP END PARALLEL DO
        end if
        ! call intpr("cluster i ", -1, i, 1)
        ! call dblepr("integrant ss ", -1, ss, 1)
    else ! dans ce cas on va faire du MPI
        ! rang du processus courang
        !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
        ! on cherche les position initiale et finale pour le processus courant
        call pos_proc_domaine(nsimu,nb_procs,rang,init_i,max_i)
        if(ndim.eq.2) then
            do ii=1,nsimu
                if((ii<init_i).or.ii>max_i) then 
                    goto 1003 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                ss=ss+func(fraili(ii,1),fraili(ii,2),0.d0,i,nsujet_trial)
                ! !print*,"ss",ss
                1003 continue
            end do
        else ! cas de 3 points
            do ii=1,nsimu
                if((ii<init_i).or.ii>max_i) then 
                    goto 1004 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                ss=ss+func(fraili(ii,1),fraili(ii,2),fraili(ii,3),i,nsujet_trial)
                ! !print*,"ss",ss
                1004 continue
            end do
        end if
        ! !print*,"rang",rang, "mon ss vaut",ss
        ! on fait la reduction et redistribu le resultat a tous les procesus
        !call MPI_ALLREDUCE(ss,ss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
        ! !print*,"rang",rang, "voila ss general",ss
        ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
    endif
    MC_Copula_Essai=ss/dble(nsimu)
    ! if(control_affichage == 0)then
        ! control_affichage = 1
        ! call intpr("ss=", -1, ss, 1)
        ! call intpr("MC_Copula_Essai=", -1, MC_Copula_Essai, 1)
    ! endif

    deallocate(vi,usim,vc,fraili)
    return
  end function MC_Copula_Essai

    subroutine monteCarlosMult(funcMC,mu,vc,nsim,vcdiag,posind_i,result)
    ! mu: l'esperance de mes variables
    ! VC: matrice de variance-covariance
    ! func: fonction a moyenner ou encore la fonction dont il faut calculer l'experance
    ! nsim: nombre de simultions
    ! vcdiag: un entier(1=oui, 0=non) qui dit si la matrice de variance covariance est diagonale ou pas. pour eviter la transformation de cholesky
    ! posind_i: position du cluster courant
    ! result: vecteur contenant le resltats de l'integrale, la variance et la precision
    use Autres_fonctions, only:init_random_seed
    use var_surrogate,only:Vect_sim_MC,a_deja_simul,graine,aleatoire,nbre_sim
    !$ use OMP_LIB
        
    implicit none
    integer :: jj,j,k,ier,l,maxmes,stemp,tid1 !maxmes= nombre de dimension ou encore dimension de X
    integer, intent(in)::nsim,vcdiag
    integer, intent(in)::posind_i
    double precision::eps,ymarg,SX,x22,somp ! ymarg contient le resultat de l'integrale
    double precision, intent(out),dimension(3):: result
    double precision,dimension(:),allocatable::usim,ysim
    double precision,dimension(:),allocatable::vi
    double precision, intent(in),dimension(:)::mu
    double precision,dimension(:,:),intent(inout)::vc
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function funcMC(nsimu,j,mu1,vc1)
            integer, intent(in)::j,nsimu! j est la position de l'individu dans le cluster i
            double precision, intent(in)::mu1,vc1
        end function funcMC
    
    ! fonction avec un seul argument
     double precision function funcSurrNN(arg)
      double precision,dimension(2), intent(in):: arg
     end function funcSurrNN
     
     ! fonction avec deux parametres
     double precision function func(arg,ndim)
      integer, intent(in)::ndim
      double precision,dimension(ndim), intent(in):: arg
     end function func
    end interface
    
    !=============debut de la fonction=============================
    x22=0.d0
    somp=0.d0
    result=0.d0
    ! VC en vecteur

    maxmes=size(vc,2)
    !!print*,"maxmes=",maxmes
    !!print*,'voila la matrive de vac-cov avant vc:',vc
     allocate(vi(maxmes*(maxmes+1)/2))
     jj=0
     Vi=0.d0
     do j=1,maxmes
        do k=j,maxmes
           jj=j+k*(k-1)/2
           Vi(jj)=VC(j,k)
        end do
     end do
     
     !!print*,'voila la matrive de vi:',vi
         
     EPS=10.d-10
     
     !!print*,'voila la matrive de vac-cov avant:',vi(1:6)
     !!print*,'dim=',size(vi,dim=1)
     
     !!print*,"vcdiag=",vcdiag
     if(vcdiag.eq.0) then
     !!print*,"suis la"
       CALL DMFSD(Vi,maxmes,eps,ier) ! si matice diagonale on na pas besoin de ceci
     end if

     !!print*,'voila la matrive de vac-cov apres:',vi
     if (ier.eq.-1) then
        ymarg=9999.d0
        result(1)=-1
        result(2)=0.d0
        result(3)=0.d0
        goto 654
     end if
     
     VC=0.d0
     do j=1,maxmes
        do k=1,j
           VC(j,k)=Vi(k+j*(j-1)/2)
        end do
     end do
    
    !!print*,'voila la cholesky de vc:',vc
    
    
  ! --------------------- boucle du MC ------------------------
    ymarg=0.d0
      allocate(usim(maxmes))
      allocate(ysim(maxmes))
    
      l=1
      stemp=0
      !===============================================================================
      ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
      !===============================================================================
!    if(a_deja_simul.eq.0.d0) then
!        do while(l.le.nsim)
!            usim=0.d0
!            do m=1,sujet_essai_max
!                SX=1.d0
!                call bgos(SX,0,Vect_sim_MC(l,m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
!            end do
            !Vect_sim_MC(l,1:sujet_essai_max)=usim
!            l=l+1
!        end do
!        a_deja_simul=1 ! pour dire qu'on ne simule plus
!    endif
    
    if(a_deja_simul.eq.0.d0) then
        call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour lagraine
        ! do m=1,sujet_essai_max
            ! do while(l.le.nsim)
                ! usim=0.d0
                ! SX=1.d0
                ! call bgos(SX,0,Vect_sim_MC(l,m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                ! l=l+1
            ! end do    
        ! end do
        ! je considere ici 1 seul jeu de simulation a utiliser pour tus les individus
        !do m=1,sujet_essai_max
            call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour le seed
            do while(l.le.nsim)
                usim=0.d0
                SX=1.d0
                call bgos(SX,0,Vect_sim_MC(l,1),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                l=l+1
            end do    
        
        !call rmvnorm(0.d0,sigma,n_essai,0,Vect_sim_MC(:,2:3))
        !end do
        a_deja_simul=1 ! pour dire qu'on ne simule plus
    endif
    
    ! debut de la boucle MC
!    l=1
!    stemp=0
!    do while(l.le.nsim)
!        usim=0.d0
!        ysim=0.d0
!        usim=Vect_sim_MC(l,1:maxmes) ! on lui affecte ce qui avait deja ete simule     
!        ysim=mu+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
!        somp=funcMC(ysim,posind_i)
!        ymarg=ymarg+somp
!        l=l+1
!        stemp=stemp+1
!    end do
    ymarg=1.d0
    !call OMP_SET_NUM_THREADS(50)
    !$OMP PARALLEL DO PRIVATE (l,tid1) firstprivate(nsim,mu,vc) REDUCTION(*:ymarg)
        do l=1,maxmes
            !tid1 = OMP_GET_THREAD_NUM ()
            !!print*,"calcul integrale bouble=",l,"thread num=",tid1,OMP_GET_NUM_THREADS()
            !ymarg=ymarg+dlog(funcMC(nsim,l,mu(l),vc(l,l)))
            ymarg=ymarg*funcMC(nsim,l,mu(l),vc(l,l))
        end do
    !$OMP END PARALLEL DO

    !ymarg=ymarg/dble(nsim)**maxmes ! la division est deja faite au fil des itterations dans l'evaluation de l'integrant
    result(1)=ymarg
    result(2)=0.d0
    result(3)=0.d0
    !!print*,"prog monteCarlosMult ligne 656: resultat MC pour trial posind_i=",posind_i,": ymarg=",ymarg
    deallocate(usim,ysim,vi)
654    continue 
    return
    end subroutine monteCarlosMult
    
    
    
    ! monte carlo cas d'un effet aleatoire individuel partage, ici le calcul se fait individu par individu
    
    subroutine monteCarlosMult_ind(funcMC,mu,vc,nsim,vcdiag,posind_ind,result)
    ! mu: l'esperance de mes variables
    ! VC: matrice de variance-covariance
    ! func: fonction a moyenner ou encore la fonction dont il faut calculer l'experance
    ! nsim: nombre de simultions
    ! vcdiag: un entier(1=oui, 0=non) qui dit si la matrice de variance covariance est diagonale ou pas. pour eviter la transformation de cholesky
    ! posind_ind: position de l'individu courant dans le cluster
    ! result: vecteur contenant le resltats de l'integrale, la variance et la precision
    use Autres_fonctions, only:init_random_seed
    use var_surrogate,only:Vect_sim_MC,a_deja_simul,sujet_essai_max,graine,aleatoire,nbre_sim
        
    implicit none
    integer :: jj,j,k,ier,l,m,maxmes,stemp !maxmes= nombre de dimension ou encore dimension de X
    integer, intent(in)::nsim,vcdiag
    integer, intent(in)::posind_ind
    double precision::eps,ymarg,SX,x22,somp ! ymarg contient le resultat de l'integrale
    double precision, intent(out),dimension(3):: result
    double precision,dimension(:),allocatable::usim,ysim
    double precision,dimension(:),allocatable::vi
    double precision, intent(in),dimension(:)::mu
    double precision,dimension(:,:),intent(inout)::vc
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function funcMC(arg,j)
            integer, intent(in)::j! j est la position de l'individu dans le cluster i
            double precision,intent(in):: arg ! fragilite wij
        end function funcMC
    
    ! fonction avec un seul argument
     double precision function funcSurrNN(arg)
      double precision,dimension(2), intent(in):: arg
     end function funcSurrNN
     
     ! fonction avec deux parametres
     double precision function func(arg,ndim)
      integer, intent(in)::ndim
      double precision,dimension(ndim), intent(in):: arg
     end function func
    end interface
    
    !=============debut de la fonction=============================
    x22=0.d0
    somp=0.d0
    result=0.d0
    ! VC en vecteur

    maxmes=size(vc,2)
    !!print*,"maxmes=",maxmes
    !!print*,'voila la matrive de vac-cov avant vc:',vc
     allocate(vi(maxmes*(maxmes+1)/2))
     jj=0
     Vi=0.d0
     do j=1,maxmes
        do k=j,maxmes
           jj=j+k*(k-1)/2
           Vi(jj)=VC(j,k)
        end do
     end do
     
     !!print*,'voila la matrive de vi:',vi
         
     EPS=10.d-10
     
     !!print*,'voila la matrive de vac-cov avant:',vi(1:6)
     !!print*,'dim=',size(vi,dim=1)
     
     !!print*,"vcdiag=",vcdiag
     if(vcdiag.eq.0) then
     !!print*,"suis la"
       CALL DMFSD(Vi,maxmes,eps,ier) ! si matice diagonale on na pas besoin de ceci
     end if

     !!print*,'voila la matrive de vac-cov apres:',vi
     if (ier.eq.-1) then
        ymarg=9999.d0
        result(1)=-1
        result(2)=0.d0
        result(3)=0.d0
        goto 654
     end if
     
     VC=0.d0
     do j=1,maxmes
        do k=1,j
           VC(j,k)=Vi(k+j*(j-1)/2)
        end do
     end do
    
    !!print*,'voila la cholesky de vc:',vc
    
    
  ! --------------------- boucle du MC ------------------------
    ymarg=0.d0
      allocate(usim(maxmes))
      allocate(ysim(maxmes))
    
      l=1
      stemp=0
      !===============================================================================
      ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
      !===============================================================================
!    if(a_deja_simul.eq.0.d0) then
!        do while(l.le.nsim)
!            usim=0.d0
!            do m=1,sujet_essai_max
!                SX=1.d0
!                call bgos(SX,0,Vect_sim_MC(l,m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
!            end do
!            l=l+1
!        end do
!        a_deja_simul=1 ! pour dire qu'on ne simule plus
!    endif

    if(a_deja_simul.eq.0.d0) then
        call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour le
        do m=1,sujet_essai_max
            do while(l.le.nsim)
                usim=0.d0
                SX=1.d0
                call bgos(SX,0,Vect_sim_MC(l,m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                l=l+1
            end do    
        end do
        a_deja_simul=1 ! pour dire qu'on ne simule plus
    endif
    
    ! ici on simule un eul vecteur et l'utilise pour tout. pas bon

!    if(a_deja_simul.eq.0.d0) then
!        do while(l.le.nsim)
!            usim=0.d0
!            SX=1.d0
!            call bgos(SX,0,Vect_sim_MC(l,1),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
!            l=l+1
!        end do
!        a_deja_simul=1 ! pour dire qu'on ne simule plus
!    endif
    
    ! debut de la boucle MC
      l=1
      stemp=0
      do while(l.le.nsim)
         usim=0.d0
         ysim=0.d0
        ! if(a_deja_simul.eq.0.d0) then ! a revoir
        !    do m=1,maxmes
                !!print*,"suis la m=",m
        !        SX=1.d0
        !        call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
        !    end do
        !    Vect_sim_MC(l,1:maxmes)=usim
        !else ! on a deja simule les donnees une fois 
            usim=Vect_sim_MC(l,1) ! on lui affecte ce qui avait deja ete simule
            !usim=Vect_sim_MC(l,posind_ind) ! on lui affecte ce qui avait deja ete simule
        !endif
         
         ysim=mu+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
         ! pour lexperance:
         !do j=1,maxmes
         !   ymarg(j)=ymarg(j)+ysim(j)/dble(nsim) ! ymarg contient l'esperance
         !end do
         somp=funcMC(ysim(1),posind_ind)
         ymarg=ymarg+somp
         l=l+1
         stemp=stemp+1
      end do
    !  !print*,"ymarg=",ymarg
      ymarg=ymarg/dble(nsim)
      result(1)=ymarg
      result(2)=0.d0
      result(3)=0.d0
      !!print*,"prog monteCarlosMult ligne 656: resultat MC pour trial posind_i=",posind_i,": ymarg=",ymarg
      deallocate(usim,ysim,vi)
654      continue 
      return
    end subroutine monteCarlosMult_ind
    
! --------------------- fin MC ------------------------      


!C ******************** BGOS ********************************
! pour la simulation des X_i suivant une gaussienne centree reduite

    SUBROUTINE BGOS(SX,ID,X1,X2,RO)
      
!C     ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)
      use var_surrogate, only: random_generator
      
      implicit none
      double precision ::RO,SX
      integer ::ID
      double precision ::F,V1,V2,S,DLS,RO2
      double precision ::X1,X2!,UNIRAN
!C     !write(*,*)'dans bgos'


 5    CONTINUE

!C     !write(*,*)'avant rand :'

!C     X1=RAND()
!C     X2=RAND()
      ! scl 27/03/2018: remplacement de uniran() par random_number(), pour pouvoir gerer le seed
      if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
          X1=UNIRAN()
          X2=UNIRAN()
      else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
          CALL RANDOM_NUMBER(X1)
          CALL RANDOM_NUMBER(X2)
      endif

      IF(ID.NE.1) GO TO 10
      F=2.d0*dSQRT(3.d0)
      X1=(X1-0.5)*F
      X2=(X2-0.5)*F
      GO TO 20
10    CONTINUE
      V1=2.d0*X1-1
      V2=2.d0*X2-1
      S=V1*V1+V2*V2
      IF(S.GE.1.) GO TO 5
      DLS=dSQRT(-2.d0*dLOG(S)/S)
      X1=V1*DLS
      X2=V2*DLS
20    CONTINUE
      RO2=RO*RO
      IF(ABS(RO).GT.1.E-10) X2=(X1+X2*dSQRT(1.d0/RO2-1.d0))*RO
      X1=X1*SX
      X2=X2*SX

!C      !write(*,*) 'X1 ',X1,' X2 ',X2
!C OK, X1 et X2 sont créés

!C      !write(*,*)'fin bgos'

      RETURN
    END subroutine bgos
!C ------------------- FIN SUBROUTINE BGOS -----------------

! =====================subroutine uniran=====================

    DOUBLE PRECISION FUNCTION UNIRAN()
!C
!C     Random number generator(RCARRY), adapted from F. James
!C     "A Review of Random Number Generators"
!C      Comp. Phys. Comm. 60(1990), pp. 329-344.
!C
      implicit none
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY, ONE
      PARAMETER ( ONE = 1, TWOM24 = ONE/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS /      &
    0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,      &
    0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
    0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
    0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      UNIRAN = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNIRAN .LT. 0 ) THEN
         UNIRAN = UNIRAN + 1
         CARRY = TWOM24
      ELSE
         CARRY = 0
      ENDIF
      SEEDS(I) = UNIRAN
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )

      return
      END function uniran
      

!C ******************** DMFSD ********************************


    subroutine dmfsd(a,n,eps,ier)
!
!   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
!   MATRICE = TRANSPOSEE(T)*T
!   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
!            PAR COLONNE DE LA METRICE A FACTORISER
!   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
!
!   SUBROUTINE APPELE PAR DSINV
!
!   N : DIM. MATRICE
!   EPS : SEUIL DE TOLERANCE
!   IER = 0 PAS D'ERREUR
!   IER = -1 ERREUR
!   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
!
      implicit none

      integer,intent(in)::n
      integer,intent(out)::ier
      double precision,intent(in)::eps
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision :: dpiv,dsum,tol
      integer::i,k,l,kpiv,ind,lend,lanf,lind

!
!   TEST ON WRONG INPUT PARAMETER N
!
      dpiv=0.d0
      if (n-1.lt.0) goto 12
      if (n-1.ge.0) ier=0
!
!   INITIALIZE DIAGONAL-LOOP
!
      kpiv=0
      do k=1,n
          kpiv=kpiv+k
          ind=kpiv
          lend=k-1
!
!   CALCULATE TOLERANCE
!
          tol=dabs(eps*sngl(A(kpiv)))
!
!   START FACTORIZATION-LOOP OVER K-TH ROW
!
         do i=k,n
            dsum=0.d0
            if (lend.lt.0) goto 2
            if (lend.eq.0) goto 4
            if (lend.gt.0) goto 2
!
!   START INNER LOOP
!
2           do l=1,lend
               lanf=kpiv-l
               lind=ind-l
               dsum=dsum+A(lanf)*A(lind)
            end do

!
!   END OF INNEF LOOP
!
!   TRANSFORM ELEMENT A(IND)
!
4           dsum=A(ind)-dsum
            if (i-k.ne.0) goto 10
            if (i-k.eq.0) goto 5
!   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
!


5           if (sngl(dsum)-tol.le.0) goto 6
            if (sngl(dsum)-tol.gt.0) goto 9
6           if (dsum.le.0) goto 12
            if (dsum.gt.0) goto 7
7           if (ier.le.0) goto 8
            if (ier.gt.0) goto 9
8           ier=k-1
!
!   COMPUTE PIVOT ELEMENT
!
9           dpiv=dsqrt(dsum)
            A(kpiv)=dpiv
            dpiv=1.D0/dpiv
            goto 11
!
!   CALCULATE TERMS IN ROW
!
10          A(ind)=dsum*dpiv
11          ind=ind+i
         end do
      end do

!
!   END OF DIAGONAL-LOOP
!
      if(ier.eq.-1) then 
        !print*,'Erreur dans le calcul de la cholesky, subroutine dmfsd: ier1=',ier
      end if
      
      return
12    ier=-1
      !print*,'Erreur dans le calcul de la cholesky, subroutine dmfsd: ier1=',ier
      return

    end subroutine dmfsd

!C ------------------- FIN SUBROUTINE DMFSD -----------------

 double precision function MC_Multiple_surr(func,vsi,vti,ui,nsimu,mu1,vc1,n,i) 
   !Monte carlo a l'aide du produit des integrales de chaque individu du cluster
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! nsimu = nombre de boucle MC
   ! mu1= moyenne du frailty
   ! vc1= variance du frailty
   ! n: nombre de sujet dans le cluster courant
   ! vsi= frailtie niveau essai associe a s
   ! vti= frailtie niveau essai associe a t
   ! ui fragilite associe au risque de base
   ! i= cluster courant
   
   use var_surrogate, only:nb_procs !&
                             !nigts,nigs
   !use comon, only:invBi_cholDet
   use comon, only: lognormal
   !$ use OMP_LIB
   
   implicit none
   
   integer ::k2
   integer, intent(in)::n,i,nsimu
   double precision,intent(in)::vsi,vti,mu1,vc1,ui
   double precision ::herm,I1
   
   ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,j,nsimu,mu1,vc1)
            ! vsi= frailtie niveau essai associe a s
            ! vti= frailtie niveau essai associe a t
            ! j = individu j du cluster i
            ! nsimu = nombre de boucle MC
            ! mu1= moyenne du frailty
            ! vc1= variance du frailty
            use var_surrogate
            use comon, only: eta,ve
            integer,intent(in):: j,nsimu
            double precision,intent(in)::vsi,vti,ui
            double precision, intent(in)::mu1,vc1
        end function func
    end interface
    
   ! fin declaration et debut du programme
    herm = 0.d0
    I1=1.d0
    if(lognormal==1)then
        herm =1.d0
        if(nb_procs <= 1) then ! on ne fait du openMp que si on a un seul processus car dans le cas MPI on parallelise en amont
            !$OMP PARALLEL DO default(none) PRIVATE (k2,I1) shared(n,vsi,vti,nsimu,mu1,vc1,ui) REDUCTION(*:herm)
                do k2=1,n
                    I1=func(vsi,vti,ui,k2,nsimu,mu1,vc1)
                    herm=herm*I1
                end do
            !$OMP END PARALLEL DO
        else
            do k2=1,n
                I1=func(vsi,vti,ui,k2,nsimu,mu1,vc1)
                herm=herm*I1
            end do
        endif
    else
        !print*,"quadrature par la loi gamma non disponible pour le modele complet de surrogacy",&
         !       "probleme de la listribution gamma bivariée"
    endif
   
   ! calcul du terme I1(vsi,vti)
    ! if(frailt_base==0)then
        ! I1=dexp(nigts(i)*vsi+cdcts(i)*vti)
    ! else
        ! I1=dexp(ui*(nigs(i)+alpha_ui*cdcs(i))+nigts(i)*vsi+cdcts(i)*vti)
    ! endif
    I1=1.d0
    MC_Multiple_surr=I1*herm
    return
 end function MC_Multiple_surr
 
    ! produit des integrales au niveau individuel MC pour modele complet
    
    double precision function MC_Multiple_surr_cor(func,vsi,vti,ui,uti,nsimu,mu1,frailij,ndim,n,i) 
       ! Monte carlo a l'aide du produit des integrales de chaque individu du cluster
       ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
       ! nsimu = nombre de boucle MC
       ! mu1= moyenne du frailty
       ! frailij= frailty au niveau individuel generes suivant la gaussienne multivariee centree
       ! n: nombre de sujet dans le cluster courant
       ! vsi= frailtie niveau essai associe a s
       ! vti= frailtie niveau essai associe a t
       ! ui fragilite associe au risque de base
       ! uti= fragilites niveau essai base a T
       ! i= cluster courant
       ! ndim= dimension de l'integrale 2 ou 1 integrations?
       
       use var_surrogate, only:cdcts,nigts,frailt_base,nigs,cdcs !,&
                                !alpha_ui
       use comon, only: lognormal
       !$ use OMP_LIB
       
       implicit none
       
       integer ::k2
       integer, intent(in)::n,i,nsimu,ndim
       double precision,intent(in)::vsi,vti,mu1,ui,uti
       double precision,dimension(nsimu,ndim),intent(in)::frailij
       double precision ::herm,I1
       
       ! bloc interface pour la definition de la fonction func
        interface
            double precision function func(vsi,vti,ui,uti,j,nsimu,ndim,mu1,frailij)
                use var_surrogate, only: delta,deltastar,const_res4,const_res5,Vect_sim_MC,frailt_base,posind_i
                use comon, only: eta,ve
                integer,intent(in):: j,nsimu,ndim
                double precision,intent(in)::vsi,vti,ui,uti
                double precision, intent(in)::mu1
                double precision,dimension(nsimu,ndim),intent(in)::frailij
            end function func
        end interface
        
       ! fin declaration et debut du programme
        herm = 0.d0
        I1=1.d0
        !!print*,"suisi laaa1"
        if(lognormal==1)then
            herm =1.d0
            !$OMP PARALLEL DO default(none) PRIVATE (k2,I1) shared(n,vsi,vti,nsimu,mu1,frailij,ui,uti,ndim) REDUCTION(*:herm)
                do k2=1,n
                    I1=func(vsi,vti,ui,uti,k2,nsimu,ndim,mu1,frailij)
                    herm=herm*I1
                end do
            !$OMP END PARALLEL DO
        else
            !print*,"quadrature par la loi gamma non disponible pour le modele complet de surrogacy",&
             !       "probleme de la listribution gamma bivariée"
        endif
       
       !!print*,"suisi laaa 2"
       !stop
       ! calcul du terme I1(vsi,vti)
        if(frailt_base==0)then
            I1=dexp(nigts(i)*vsi+cdcts(i)*vti)
        else
            I1=dexp(ui*nigs(i)+uti*cdcs(i)+nigts(i)*vsi+cdcts(i)*vti)
        endif
        
        ! !print*,"trial=",i,herm,I1
        ! stop
        MC_Multiple_surr_cor=I1*herm
        return
    end function MC_Multiple_surr_cor
    
      ! calcul de l'integral au niveau essai par Monte-carlo et par quadrature au niveau individuel
 double precision function MC_Gauss_MultInd_Essai(func,func2,ndim,nsujet_trial,i,npoint)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! func: fonction a integrer au niveau individuel
    ! mu1= moyenne du frailty
    ! vc1= variance du frailty
    ! ndim= dimension de l'integrale 2 ou 3 integrations?
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    ! npoint= nombre de point de quadrature
    use Autres_fonctions, only:init_random_seed
    use parameters, only: maxiter
    use comon, only: model
    use func_adaptative, only: funcpafrailtyPred_ind
    use optim_scl, only:marq98j_scl  ! pour faire appel a marquard 
    use var_surrogate, only: Vect_sim_MC,a_deja_simul,nsim,chol,frailt_base,& !sujet_essai_max,theta2
                             graine,aleatoire,nbre_sim,nsujeti,essai_courant,indicej,& !gamma_ui,alpha_ui
                             vs_i,vt_i,u_i,invBi_chol_Individuel,ui_chap,adaptative,control_adaptative,&
                             nparamfrail,ntrials,switch_adaptative,nb_procs
    use Autres_fonctions, only:pos_proc_domaine
    use comon, only:invBi_cholDet
    !use mpi
                             
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    !$ use OMP_LIB
    
    implicit none
    integer ::ii,jj,kk,l,m,maxmes,nsimu,ig,ind_frail,init_i,max_i,rang  !maxmes= nombre de dimension ou encore dimension de X !npg,j,k,nbrejet,stemp,tid1,code,erreur
            
    integer,intent(in):: ndim,nsujet_trial,i
    !double precision,dimension(1:nnodes) ::xx1,ww1
    double precision::auxfunca,ss,ca,cb,dd,res
    double precision,dimension(:,:),allocatable::vc
    double precision,dimension(:,:),allocatable::fraili
    double precision,dimension(:),allocatable::usim
    double precision::ymarg,SX,x22,somp ! ymarg contient le resultat de l'integrale
    !double precision,dimension(:),allocatable::ysim
    double precision,dimension(:),allocatable::vi
    integer,intent(in):: npoint
    double precision, dimension(2)::k0_2
    double precision, dimension(:),allocatable::v,b_2
    double precision, allocatable, dimension(:,:)::H_hessOut,H_hess_scl,I_hess_scl
    double precision,dimension(:,:), allocatable::hess_scl
    double precision,dimension(:), allocatable::vvv_scl
    integer::ier,istop,sss,ni,nmax_2,np_2,nparamfrail_save,model_save,maxiter_save,individu_j,np_1,effet2,&
            non_conv

    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,j,npoint1)
        ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
        ! j = individu j du cluster i
        ! npoint1 = nombre de point de quadrature
            integer,intent(in):: j,npoint1
            double precision,intent(in)::vsi,vti,ui
        end function func
        
        double precision function func2(func,vsi,vti,ui,npoint1,n,i)
            integer, intent(in)::n,npoint1,i
            double precision,intent(in)::vsi,vti,ui
            
            interface
                double precision function func(vsi,vti,ui,j,npoint1)
                    integer,intent(in):: j,npoint1
                    double precision,intent(in)::vsi,vti,ui
                end function func
            end interface
        end function func2
    end interface
    
    allocate(vc(size(chol,2),size(chol,2)),fraili(nsim,size(chol,2)))
    !vc=ABS(chol)
    vc=chol
    nsimu=nsim
    x22=0.d0
    somp=0.d0

    maxmes=size(vc,2)
     allocate(vi(maxmes*(maxmes+1)/2))
     allocate(usim(maxmes))
    
    !!print*,'voila la cholesky de vc:',vc
    !!print*,"mat_A vs,vt,covst=",mat_A(1,1),mat_A(2,2),mat_A(1,2)
    !!print*,"chol vs,vt,covst=",vc(1,1),vc(2,1),vc(1,2),vc(2,2)
    
    
  ! --------------------- boucle du MC ------------------------
    ymarg=0.d0
    l=1
    !stemp=0
    !===============================================================================
    ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
    !===============================================================================
    
    if(a_deja_simul.eq.0) then
        call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour le
        !!print*,"nsimu=",nsimu,size(Vect_sim_MC,1),size(Vect_sim_MC,2)
        !stop
        Vect_sim_MC=0.d0
        do while(l.le.nsimu)
            ! pour integrer sur un seul effet aleatoire au niveau individuel, on genere seulement des normales centree reduites, la transformation se fait dans le calcul de l'integran (fichier integrant.f90)
            usim=0.d0
            SX=1.d0
            call bgos(SX,0,Vect_sim_MC(l,1),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
            
            ! on genere suivant des normales centrees reduites pour les variables aleatoires correlees
             usim=0.d0
             do m=1,maxmes
                 SX=1.d0
                 !call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                 call bgos(SX,0,Vect_sim_MC(l,m+1),x22,0.d0)
             end do
            !Vect_sim_MC(l,2:3)=0.d0+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
            !!print*,"l=",l,"Vect_sim_MC(l,2:3)=",Vect_sim_MC(l,2:3),maxmes
            l=l+1
        end do    
        
        a_deja_simul=1 ! pour dire qu'on ne simule plus
    endif
    
    ! on utilise les generations precedentes pour obtenir deux variables correlees suivant une multinormale centree de covariance vc
    l=1
    do while(l.le.nsimu)
        if(frailt_base==0)then
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,2:3)) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
        else
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,2:4))
        endif
        l=l+1
    end do
    
    !prediction des effects aleatoires a posteori au niveau individuel en cas de quadrature pseudo-adaptative
    if(adaptative .and. control_adaptative==1) then ! on effectue le changement de variable
        nmax_2=0
        kk=1
        do ig=1,ntrials                 
            indicej=kk
            nmax_2=nmax_2+nsujeti(ig)
            essai_courant=ig
            
            do ii=indicej,nmax_2
                individu_j=ii
                                    
                !====================================================================================================
                ! estimation des w_ij
                !====================================================================================================    
                ! !print*,"==============================================="
                ! !print*,"Recherche des effets aleatoires  à postériorie"
                ! !print*,"==============================================="
                
                k0_2=0.d0
                ni=0
                ca=0.d0
                cb=0.d0
                dd=0.d0
                non_conv=0        
                ind_frail=1                
                np_2=1
                np_1=1
                effet2=0
                !call intpr("je suis la pour pseudo-adpdative 1136", -1, adaptative, 1)
                allocate(I_hess_scl(np_2,np_2),v(np_2*(np_2+3)/2),b_2(1))
                allocate(H_hess_scl(np_2,np_2),hess_scl(np_2,np_2),vvv_scl(np_2*(np_2+1)/2))
                allocate(H_hessOut(np_2,np_2))
                !allocate(HIH(np_2,np_2),HIHOut(np_2,np_2),IH(np_2,np_2),invBi_chol_2(np_2,np_2))
                b_2(1)=0.5d0
                v=0.d0
                
                1241 continue
                
                nparamfrail_save=nparamfrail
                model_save=model
                maxiter_save=maxiter
                nparamfrail=1
                model = 9 !scl pour le model effet aleatoires
                maxiter=10
                vs_i=fraili(ind_frail,1)
                vt_i=fraili(ind_frail,2)
                u_i=fraili(ind_frail,3)
                
                call marq98J_scl(k0_2,b_2,np_1,ni,v,res,ier,istop,effet2,ca,cb,dd,funcpafrailtyPred_ind,I_hess_scl,H_hess_scl,&
                                hess_scl,vvv_scl,individu_j)
                                    
                nparamfrail=nparamfrail_save ! on restitu sa valeur avant de continuer
                model=model_save
                maxiter=maxiter_save
                
                
                !call dblepr("b_2 pseudo-adpd 1165", -1, b_2, 1)
                if (istop.ne.1 .and. ind_frail.eq.5) then
                    ind_frail=ind_frail+1 ! on prend un autre jeux d'effets aleatoire: le suivant
                    goto 1241
                endif
                
                !call intpr("istop pour pseudo-adpdative 1171", -1, istop, 1)
                if(istop .ne.1) then
                    non_conv=1
                    ! !print*,"2-individu",ii,"wij=",b_2,"istop=",istop,"ier=",ier,"v=",v
                    ! !print*,"le modele d'estimation des w_ij_chapeau n'a pas converge, i=",ii
                    switch_adaptative=0 ! pour cette itteration, calculer l'intergralle avec la pseudo adaptative
                    goto 1242
                endif
                
                switch_adaptative=1 ! Cas de convergence, car goto non execute
                ui_chap(ii,1)=b_2(1)        
                do jj=1,np_1
                    do sss=1,np_1
                        H_hessOut(jj,sss)= I_hess_scl(jj,sss)
                    end do
                end do
                                    
                invBi_chol_Individuel(ii)=dsqrt(H_hess_scl(1,1))    
                !calcul du determinant de la cholesky de l'inverse de la hessienne                    
                invBi_cholDet(ii)=invBi_chol_Individuel(ii) !individuel
                
                
                deallocate(H_hessOut)
                !deallocate(HIH,HIHOut,IH,invBi_chol_2)
                deallocate(H_hess_scl)
                
                deallocate(I_hess_scl)
                deallocate(hess_scl)
                deallocate(vvv_scl)
                deallocate(v)
                deallocate(b_2)
            enddo ! fin estimation des w_ij_chapeau
            kk=nmax_2+1 ! on continu avec le premier sujet du prochain cluster
        enddo
        
        ! !print*,"Fin Recherche des effets aleatoires  à postériorie"
        ! !print*,"==============================================="
        control_adaptative=0
        1242 continue 
        ! if(non_conv.eq.1) MC_Gauss_MultInd_Essai=-1.d9
    endif
            
    
                    
                    ! ====================================================================================================
                    ! Fin estimation des ws_ij_chapeau et wt_ij_chapeau
                    ! ====================================================================================================
    
    !integration sur vsi et vti
    auxfunca=0.d0
    ss=0.d0
    if(nb_procs==1) then !on fait du open MP car un seul processus
        rang=0
        if(ndim.eq.2) then
          !!print*,'je suis la'
            !ii=0
            !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,nsujet_trial,i,fraili,npoint)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
                do ii=1,nsimu
                    auxfunca=func2(func,fraili(ii,1),fraili(ii,2),0.d0,npoint,nsujet_trial,i)
                    ss=ss+auxfunca
                    !!print*,"ss",ss
                end do
            !$OMP END PARALLEL DO
        else ! cas de 3 points
            !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,nsujet_trial,i,fraili,npoint)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
                do ii=1,nsimu
                    auxfunca=func2(func,fraili(ii,1),fraili(ii,2),fraili(ii,3),npoint,nsujet_trial,i)
                    ss=ss+auxfunca
                    !!print*,"ss",ss
                end do
            !$OMP END PARALLEL DO
        end if
    else ! dans ce cas on va faire du MPI
        !rang du processus courang
        !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
        !on cherche les position initiale et finale pour le processus courant
        call pos_proc_domaine(nsimu,nb_procs,rang,init_i,max_i)
        if(ndim.eq.2) then
            do ii=1,nsimu
                if((ii<init_i).or.ii>max_i) then 
                    goto 1005 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                auxfunca=func2(func,fraili(ii,1),fraili(ii,2),0.d0,npoint,nsujet_trial,i)
                ss=ss+auxfunca
                1005 continue
            end do
        else ! cas de 3 points
            do ii=1,nsimu
                if((ii<init_i).or.ii>max_i) then 
                    goto 1006 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                auxfunca=func2(func,fraili(ii,1),fraili(ii,2),fraili(ii,3),npoint,nsujet_trial,i)
                ss=ss+auxfunca
                1006 continue
            end do
        end if
        !on fait la reduction et redistribu le resultat a tous les procesus
        !call MPI_ALLREDUCE(ss,ss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
    endif
    
    ! !print*,"ss=",ss
    ! stop
    MC_Gauss_MultInd_Essai=ss/dble(nsimu)
    !!print*,"ss",ss,"dble(nsimu)",dble(nsimu),"MC_MultInd_Essai",MC_MultInd_Essaii
    !stop
    deallocate(vi,usim,vc,fraili)
    ! 1242 continue 
    ! if(non_conv.eq.1) MC_Gauss_MultInd_Essai=-1.d9
    return
  end function MC_Gauss_MultInd_Essai
  
! calcul de l'integral au niveau essai pour le modele surrogate final par Monte-carlo
 double precision function MC_MultInd_Essai(func,func2,ndim,nsujet_trial,i,mat_A)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! func: fonction a integrer au niveau individuel
    ! mu1= moyenne du frailty
    ! vc1= variance du frailty
    ! ndim= dimension de l'integrale 2 ou 3 integrations?
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    use Autres_fonctions, only:init_random_seed
    use var_surrogate, only: Vect_sim_MC,a_deja_simul,theta2,nsim,chol,frailt_base,&
                             graine,aleatoire,nbre_sim,nb_procs !alpha_ui,gamma_ui,sujet_essai_max
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use Autres_fonctions, only:pos_proc_domaine
    !use mpi
    !$ use OMP_LIB
    
    implicit none
    integer ::ii,l,m,maxmes,nsimu,init_i,max_i,rang !maxmes= nombre de dimension ou encore dimension de X !nbrejet,stemp,tid1,jj,npg,kk,j,k,ier,code,erreur
    integer,intent(in):: ndim,nsujet_trial,i
    !double precision,dimension(1:nnodes) ::xx1,ww1
    double precision::auxfunca,ss
    double precision,dimension(:,:),allocatable::vc
    double precision,dimension(:,:),allocatable::fraili
    double precision,dimension(:),allocatable::usim
    double precision::ymarg,SX,x22,somp ! ymarg contient le resultat de l'integrale
    !double precision,dimension(:),allocatable::ysim
    double precision,dimension(:),allocatable::vi
    double precision,dimension(2,2),intent(in):: mat_A

    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,j,nsimu,mu1,vc1)
            ! vsi= frailtie niveau essai associe a s
            ! vti= frailtie niveau essai associe a t
            ! j = individu j du cluster i
            ! nsimu = nombre de boucle MC
            ! mu1= moyenne du frailty
            ! vc1= variance du frailty
            integer,intent(in):: j,nsimu
            double precision,intent(in)::vsi,vti,ui
            double precision, intent(in)::mu1,vc1
        end function func
        
        double precision function func2(func,vsi,vti,ui,nsimu,mu1,vc1,n,i)
            integer, intent(in)::n,nsimu,i
            double precision,intent(in)::vsi,vti,ui
            double precision, intent(in)::mu1,vc1
            
            interface
                double precision function func(vsi,vti,ui,j,nsimu,mu1,vc1)
                    integer,intent(in):: j,nsimu
                    double precision,intent(in)::vsi,vti,ui
                    double precision, intent(in)::mu1,vc1
                end function func
            end interface
        end function func2
    end interface
    
    allocate(vc(size(chol,2),size(chol,2)),fraili(nsim,size(chol,2)))
    !vc=ABS(chol)
    vc=chol
    nsimu=nsim
    x22=0.d0
    somp=0.d0

    maxmes=size(vc,2)
     allocate(vi(maxmes*(maxmes+1)/2))
     allocate(usim(maxmes))
    
    !!print*,'voila la cholesky de vc:',vc
    !!print*,"mat_A vs,vt,covst=",mat_A(1,1),mat_A(2,2),mat_A(1,2)
    !!print*,"chol vs,vt,covst=",vc(1,1),vc(2,1),vc(1,2),vc(2,2)
    
    
  ! --------------------- boucle du MC ------------------------
    ymarg=0.d0
    l=1
    !stemp=0
    !===============================================================================
    ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
    !===============================================================================
    
    if(a_deja_simul.eq.0) then
        call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour lagraine
        !!print*,"nsimu=",nsimu,size(Vect_sim_MC,1),size(Vect_sim_MC,2)
        ! !print*,"debut generation des donnes pour le MC"
        !stop
        Vect_sim_MC=0.d0
        do while(l.le.nsimu)
            ! pour integrer sur un seul effet aleatoire au niveau individuel, on genere seulement des normales centree reduites, la transformation se fait dans le calcul de l'integran (fichier integrant.f90)
            usim=0.d0
            SX=1.d0
            call bgos(SX,0,Vect_sim_MC(l,1),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
            
            ! on genere suivant des normales centrees reduites pour les variables aleatoires correlees
             usim=0.d0
             do m=1,maxmes
                 SX=1.d0
                 !call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                 call bgos(SX,0,Vect_sim_MC(l,m+1),x22,0.d0)
             end do
            !Vect_sim_MC(l,2:3)=0.d0+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
            !!print*,"l=",l,"Vect_sim_MC(l,2:3)=",Vect_sim_MC(l,2:3),maxmes
            l=l+1
        end do    
        
        a_deja_simul=1 ! pour dire qu'on ne simule plus
    endif
    
    ! on utilise les generations precedentes pour obtenir deux variables correlees suivant une multinormale centree de covariance vc
    l=1
    do while(l.le.nsimu)
        if(frailt_base==0)then
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,2:3)) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
        else
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,2:4))
        endif
        l=l+1
    end do
    
    !usim=Vect_sim_MC(1:nsimu,1) ! on lui affecte a chaque fois le premier vecteur des donnees simulees
    !fraili=Vect_sim_MC(1:nsimu,2:3)
    ! do l=1,nsim
      ! !print*,fraili(l,:)
    ! enddo
    ! stop
    
    !integration sur vsi et vti
    auxfunca=0.d0
    ss=0.d0
    if(nb_procs==1) then !on fait du open MP car un seul processus
        rang=0
        if(ndim.eq.2) then
          !!print*,'je suis la'
            !ii=0
            !$OMP PARALLEL DO default(none) PRIVATE (ii) SHARED(nsimu,nsujet_trial,i,fraili,theta2)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
                do ii=1,nsimu
                    ss=ss+func2(func,fraili(ii,1),fraili(ii,2),0.d0,nsimu,0.d0,theta2,nsujet_trial,i)
                    !!print*,"ss",ss
                end do
            !$OMP END PARALLEL DO
        else ! cas de 3 points
            !$OMP PARALLEL DO default(none) PRIVATE (ii) SHARED(nsimu,nsujet_trial,i,fraili,theta2)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
                do ii=1,nsimu
                    ss=ss+func2(func,fraili(ii,1),fraili(ii,2),fraili(ii,3),nsimu,0.d0,theta2,nsujet_trial,i)
                    ! !print*,"ss",ss
                end do
            !$OMP END PARALLEL DO
        end if
    else ! dans ce cas on va faire du MPI
        ! rang du processus courang
        !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
        ! on cherche les position initiale et finale pour le processus courant
        call pos_proc_domaine(nsimu,nb_procs,rang,init_i,max_i)
        if(ndim.eq.2) then
          !!print*,'je suis la'
            !ii=0
            do ii=1,nsimu
                if((ii<init_i).or.ii>max_i) then 
                    goto 1003 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                ss=ss+func2(func,fraili(ii,1),fraili(ii,2),0.d0,nsimu,0.d0,theta2,nsujet_trial,i)
                ! !print*,"ss",ss
                1003 continue
            end do
        else ! cas de 3 points
            do ii=1,nsimu
                if((ii<init_i).or.ii>max_i) then 
                    goto 1004 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                ss=ss+func2(func,fraili(ii,1),fraili(ii,2),fraili(ii,3),nsimu,0.d0,theta2,nsujet_trial,i)
                ! !print*,"ss",ss
                1004 continue
            end do
        end if
        ! !print*,"rang",rang, "mon ss vaut",ss
        ! on fait la reduction et redistribu le resultat a tous les procesus
        !call MPI_ALLREDUCE(ss,ss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
        ! !print*,"rang",rang, "voila ss general",ss
        ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
    endif
    
    ! !print*,"rang",rang, "voila ss general",ss
    ! stop
    
    MC_MultInd_Essai=ss/dble(nsimu)
    !!print*,"ss",ss,"dble(nsimu)",dble(nsimu),"MC_MultInd_Essai",MC_MultInd_Essaii
    !stop
    deallocate(vi,usim,vc,fraili)
    return
  end function MC_MultInd_Essai

  
    ! calcul de l'integral au niveau essai par Monte-carlo et par quadrature au niveau individuel
 double precision function MC_Gauss_MultInd_Essai_Cor(func,func2,ndim_Ind,ndim,nsujet_trial,i,npoint)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire wsi et wti
    ! func: fonction a integrer au niveau individuel
    ! mu1= moyenne du frailty
    ! vc1= variance du frailty
    ! ndim_Ind= dimension de l'integrale niveau individuel
    ! ndim= dimension de l'integrale niveau essai
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    ! npoint= nombre de point de quadrature
    use Autres_fonctions, only:init_random_seed
    use var_surrogate, only: Vect_sim_MC,a_deja_simul,nsim,chol,frailt_base,graine,aleatoire,nbre_sim,nb_procs
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use Autres_fonctions, only:pos_proc_domaine
    !use mpi
    !$ use OMP_LIB
    
    implicit none
    integer ::ii,jj,l,m,maxmes,nsimu,init_i,max_i,rang !maxmes= nombre de dimension ou encore dimension de X !npg,kk,j,k,ier,nbrejet,stemp,tid1,code,erreur
    integer,intent(in):: ndim,nsujet_trial,i,ndim_Ind
    !double precision,dimension(1:nnodes) ::xx1,ww1
    double precision::auxfunca,ss
    double precision,dimension(:,:),allocatable::vc
    double precision,dimension(:,:),allocatable::fraili
    double precision,dimension(:),allocatable::usim
    double precision::ymarg,SX,x22,somp ! ymarg contient le resultat de l'integrale
    !double precision,dimension(:),allocatable::ysim
    double precision,dimension(:),allocatable::vi
    integer,intent(in):: npoint

    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,uti,nnodes,ndim,j)
            integer,intent(in):: j,ndim,nnodes
            double precision,intent(in)::vsi,vti,ui,uti
        end function func
        
        double precision function func2(func,vsi,vti,ui,uti,nnodes,ndim,nsujet_trial,i)
            double precision,intent(in)::vsi,vti,ui,uti
            integer,intent(in):: ndim,nnodes,nsujet_trial,i
            
            interface
                double precision function func(vsi,vti,ui,uti,nnodes,ndim,j)
                    integer,intent(in):: j,ndim,nnodes
                    double precision,intent(in)::vsi,vti,ui,uti
                end function func
            end interface
        end function func2
    end interface
    
    allocate(vc(ndim,ndim),fraili(nsim,ndim))
    !vc=ABS(chol)
    ! do jj=1,ndim
        ! do ii=1,ndim
            ! vc(ii,jj)=chol(ii,jj)
        ! enddo
    ! enddo
    
    if(ndim_Ind==2) then    
        !cholesky niveau essai
        do jj=(ndim_Ind+1),(ndim_Ind+ndim)
            do ii=(ndim_Ind+1),(ndim_Ind+ndim)
                !vc(ii,jj)=chol(ii,jj)
                vc(ii-ndim_Ind,jj-ndim_Ind)=chol(ii,jj)
            enddo
        enddo
    else
        !print*,"bien vouloir gerer le cas avec une fragilité partagé au niveau individuel avant de continuer"
    endif    
    
    nsimu=nsim
    x22=0.d0
    somp=0.d0

    maxmes=size(vc,2)
     allocate(vi(maxmes*(maxmes+1)/2))
     allocate(usim(maxmes))
    
    !!print*,'voila la cholesky de vc:',vc
    !!print*,"mat_A vs,vt,covst=",mat_A(1,1),mat_A(2,2),mat_A(1,2)
    !!print*,"chol vs,vt,covst=",vc(1,1),vc(2,1),vc(1,2),vc(2,2)
    
    
  ! --------------------- boucle du MC ------------------------
    ymarg=0.d0
    l=1
    !stemp=0
    !===============================================================================
    ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
    !===============================================================================
    
    if(a_deja_simul.eq.0) then
        call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour lagraine
        !!print*,"nsimu=",nsimu,size(Vect_sim_MC,1),size(Vect_sim_MC,2)
        !stop
        Vect_sim_MC=0.d0
        do while(l.le.nsimu)
            ! pour integrer sur un seul effet aleatoire au niveau individuel, on genere seulement des normales centree reduites, la transformation se fait dans le calcul de l'integran (fichier integrant.f90)
            usim=0.d0
            SX=1.d0
            call bgos(SX,0,Vect_sim_MC(l,1),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite pour ws_ij
            call bgos(SX,0,Vect_sim_MC(l,2),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite pour wt_ij
            
            ! on genere suivant des normales centrees reduites pour les variables aleatoires correlees
             usim=0.d0
             do m=1,maxmes
                 SX=1.d0
                 !call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                 call bgos(SX,0,Vect_sim_MC(l,m+2),x22,0.d0)
             end do
            !Vect_sim_MC(l,2:3)=0.d0+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
            !!print*,"l=",l,"Vect_sim_MC(l,2:3)=",Vect_sim_MC(l,2:3),maxmes
            l=l+1
        end do    
        
        a_deja_simul=1 ! pour dire qu'on ne simule plus
    endif
    
    ! on utilise les generations precedentes pour obtenir deux variables correlees suivant une multinormale centree de covariance vc
    l=1
    do while(l.le.nsimu)
        if(frailt_base==0)then
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,3:4)) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
        else
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,3:6))
        endif
        l=l+1
    end do
    
    !integration sur usi, uti, vsi et vti
    auxfunca=0.d0
    ss=0.d0
    select case(ndim)
    case(2)
      !!print*,'je suis la'
        !ii=0
        !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,nsujet_trial,i,fraili,npoint,ndim_Ind)&
        !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do ii=1,nsimu
                auxfunca=func2(func,fraili(ii,1),fraili(ii,2),0.d0,0.d0,npoint,ndim_Ind,nsujet_trial,i)
                ss=ss+auxfunca
                !!print*,"ss",ss
            end do
        !$OMP END PARALLEL DO

    case(3) ! cas de 3 points
        if(nb_procs==1) then !on fait du open MP car un seul processus
            rang=0
            !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,nsujet_trial,i,fraili,npoint,ndim_Ind)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
                do ii=1,nsimu
                    auxfunca=func2(func,fraili(ii,2),fraili(ii,3),fraili(ii,1),0.d0,npoint,ndim_Ind,nsujet_trial,i)
                    ss=ss+auxfunca
                    !!print*,"ss",ss
                end do
            !$OMP END PARALLEL DO
        else ! dans ce cas on va faire du MPI
            ! rang du processus courang
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
            ! on cherche les position initiale et finale pour le processus courant
            call pos_proc_domaine(nsimu,nb_procs,rang,init_i,max_i)
            do ii=1,nsimu
                if((ii<init_i).or.ii>max_i) then 
                    goto 1002 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                auxfunca=func2(func,fraili(ii,2),fraili(ii,3),fraili(ii,1),0.d0,npoint,ndim_Ind,nsujet_trial,i)
                ss=ss+auxfunca
                1002 continue
            end do
            ! on fait la reduction et redistribu le resultat a tous les procesus
            !call MPI_ALLREDUCE(ss,ss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
        endif
        
    case(4)! vsi vti usi uti, 
        !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,nsujet_trial,i,fraili,npoint,ndim_Ind)&
        !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do ii=1,nsimu
                auxfunca=func2(func,fraili(ii,3),fraili(ii,4),fraili(ii,1),fraili(ii,2),npoint,ndim_Ind,nsujet_trial,i)
                ss=ss+auxfunca
                !!print*,"ss",ss
            end do
        !$OMP END PARALLEL DO
    end select
    
    MC_Gauss_MultInd_Essai_Cor=ss/dble(nsimu)
    !!print*,"ss",ss,"dble(nsimu)",dble(nsimu),"MC_MultInd_Essai",MC_MultInd_Essaii
    !stop
    deallocate(vi,usim,vc,fraili)
    return
  end function MC_Gauss_MultInd_Essai_Cor
  
      ! calcul de l'integral au niveau essai et individuel par Monte-carlo 
 double precision function MC_MultInd_Essai_Cor(func,func2,ndim_Ind,ndim,nsujet_trial,i)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire wsi et wti
    ! func: fonction a integrer au niveau individuel
    ! mu1= moyenne du frailty
    ! vc1= variance du frailty
    ! ndim_Ind= dimension de l'integrale niveau individuel
    ! ndim= dimension de l'integrale niveau essai
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    use Autres_fonctions, only:init_random_seed
    use var_surrogate, only: Vect_sim_MC,a_deja_simul,nsim,chol,frailt_base,graine,aleatoire,nbre_sim
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    !$ use OMP_LIB
    
    implicit none
    integer ::ii,jj,l,m,maxmes,nsimu,nfrail2 !maxmes= nombre de dimension ou encore dimension de X !npg,kk,j,k,ier,nbrejet,stemp,tid1
    integer,intent(in):: ndim,nsujet_trial,i,ndim_Ind
    !double precision,dimension(1:nnodes) ::xx1,ww1
    double precision::auxfunca,ss !ss1,ss2,mu1
    double precision,dimension(:,:),allocatable::vc,vc1!,usim_
    double precision,dimension(:,:),allocatable::fraili,frailij
    double precision,dimension(:),allocatable::usim
    double precision::ymarg,SX,x22,somp ! ymarg contient le resultat de l'integrale
    !double precision,dimension(:),allocatable::ysim
    double precision,dimension(:),allocatable::vi

    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,uti,j,nsimu,ndim,mu1,frailij)
        use Autres_fonctions, only:init_random_seed
            use var_surrogate, only: delta,deltastar,const_res4,const_res5,Vect_sim_MC,frailt_base,posind_i,&
                graine,aleatoire,nbre_sim
            
            integer,intent(in):: j,nsimu,ndim
            double precision,intent(in)::vsi,vti,ui,uti
            double precision, intent(in)::mu1
            double precision,dimension(nsimu,ndim),intent(in)::frailij
        end function func
               
        double precision function func2(func,vsi,vti,ui,uti,nsimu,mu1,frailij,ndim,n,i)
            integer, intent(in)::n,i,nsimu,ndim
            double precision,intent(in)::vsi,vti,mu1,ui,uti
            double precision,dimension(nsimu,ndim),intent(in)::frailij
            
            interface
                double precision function func(vsi,vti,ui,uti,j,nsimu,ndim,mu1,frailij)
                    use var_surrogate, only: delta,deltastar,const_res4,const_res5,Vect_sim_MC,frailt_base,posind_i
                    integer,intent(in):: j,nsimu,ndim
                    double precision,intent(in)::vsi,vti,ui,uti
                    double precision, intent(in)::mu1
                    double precision,dimension(nsimu,ndim),intent(in)::frailij
                end function func
            end interface
        end function func2
    end interface
    
    allocate(vc(ndim,ndim),fraili(nsim,ndim),frailij(nsim,ndim_Ind),vc1(ndim_Ind,ndim_Ind))
    
    !!print*,"size(chol)",size(chol)
    !stop
    if(ndim_Ind==2) then
        ! cholesky niveau individuel
        do jj=1,ndim_Ind
            do ii=1,ndim_Ind
                vc1(ii,jj)=chol(ii,jj)
                !!print*,ii,jj,vc1(ii,jj),chol(ii,jj)
            enddo
        enddo
        
        !cholesky niveau essai
        do jj=(ndim_Ind+1),(ndim_Ind+ndim)
            do ii=(ndim_Ind+1),(ndim_Ind+ndim)
                vc(ii-ndim_Ind,jj-ndim_Ind)=chol(ii,jj)
                !!print*,ii,jj,vc(ii-ndim_Ind,jj-ndim_Ind),chol(ii,jj)
            enddo
        enddo
    else
      !print*,"bien vouloir gerer le cas avec une fragilité partagé au niveau individuel avant de continuer"
    endif
    
    !!print*,"suis lààààà 1"
    !stop
    
    nsimu=nsim
    x22=0.d0
    somp=0.d0

    maxmes=size(vc,2)
    allocate(vi(maxmes*(maxmes+1)/2))
    allocate(usim(maxmes))
    
    !!print*,'voila la cholesky de vc:',vc
    !!print*,"mat_A vs,vt,covst=",mat_A(1,1),mat_A(2,2),mat_A(1,2)
    !!print*,"chol vs,vt,covst=",vc(1,1),vc(2,1),vc(1,2),vc(2,2)
    
    
  ! --------------------- boucle du MC ------------------------
    ymarg=0.d0
    l=1
    !stemp=0
    !===============================================================================
    ! initialisation de la matrice des donnees generees pour l'estimation de l'integrale 
    !===============================================================================
    
    if(a_deja_simul.eq.0) then
        call init_random_seed(graine,aleatoire,nbre_sim)! initialisation de l'environnement de generation pour lagraine
        !!print*,"nsimu=",nsimu,size(Vect_sim_MC,1),size(Vect_sim_MC,2)
        !print*,"debut generation des donnees pour le MC"
        !stop
        Vect_sim_MC=0.d0
        do while(l.le.nsimu)
            ! pour integrer sur un seul effet aleatoire au niveau individuel, on genere seulement des normales centree reduites, la transformation se fait dans le calcul de l'integran (fichier integrant.f90)
            usim=0.d0
            SX=1.d0
            call bgos(SX,0,Vect_sim_MC(l,1),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite pour ws_ij
            call bgos(SX,0,Vect_sim_MC(l,2),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite pour wt_ij
            l=l+1
        enddo    
        if(frailt_base==0)then
            l=1
            do while(l.le.nsimu)
                ! on genere suivant des normales centrees reduites pour les variables aleatoires correlees au niveau essai
                 usim=0.d0
                 do m=1,maxmes
                     SX=1.d0
                     !call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                     call bgos(SX,0,Vect_sim_MC(l,m+2),x22,0.d0)
                 end do
                !Vect_sim_MC(l,2:3)=0.d0+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
                !!print*,"l=",l,"Vect_sim_MC(l,2:3)=",Vect_sim_MC(l,2:3),maxmes
                l=l+1
            end do    
        else
            ! donnes pour les effets aleatoires associes au risque de base
            nfrail2=maxmes-2
            l=1
            do while(l.le.nsimu)
                ! on genere suivant des normales centrees reduites pour les variables aleatoires correlees au niveau essai
                 usim=0.d0
                 do m=1,nfrail2
                     SX=1.d0
                     !call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                     call bgos(SX,0,Vect_sim_MC(l,m+2),x22,0.d0)
                 end do
                !Vect_sim_MC(l,2:3)=0.d0+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
                !!print*,"l=",l,"Vect_sim_MC(l,2:3)=",Vect_sim_MC(l,2:3),maxmes
                l=l+1
            end do    
            
            ! donnes pour les effets aleatoires associes en interaction avec le traitement
            l=1
            do while(l.le.nsimu)
                ! on genere suivant des normales centrees reduites pour les variables aleatoires correlees au niveau essai
                 usim=0.d0
                 do m=1,nfrail2
                     SX=1.d0
                     !call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
                     call bgos(SX,0,Vect_sim_MC(l,m+4),x22,0.d0)
                 end do
                !Vect_sim_MC(l,2:3)=0.d0+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
                !!print*,"l=",l,"Vect_sim_MC(l,2:3)=",Vect_sim_MC(l,2:3),maxmes
                l=l+1
            end do    
        endif
                    
        a_deja_simul=1 ! pour dire qu'on ne simule plus
    endif
    
    ! on utilise les generations precedentes pour obtenir deux variables correlees suivant une multinormale centree de covariance vc
    l=1
    !!print*,"suis lààààà 1"
    
    ! allocate(usim_(2,1))
    do while(l.le.nsimu)
        if(frailt_base==0)then
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,3:4)) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
        else
            fraili(l,:)=0.d0+MATMUL(vc,Vect_sim_MC(l,3:6))
        endif        
        ! frailties gaussiens correlees au niveau individuel
        ! usim_(1,1)=Vect_sim_MC(l,1)
        ! usim_(2,1)=Vect_sim_MC(l,2)
        !!print*,"MATMUL(vc1,usim_)=",MATMUL(vc1,Vect_sim_MC(l,2))
        frailij(l,:)=0.d0+MATMUL(vc1,Vect_sim_MC(l,1:2))
        !!print*,"frailij(l,:)",frailij(l,:)
        l=l+1
    end do
    
    ! !print*,"suis lààààà 2"
    !stop
    
    !integration sur usi, uti, vsi et vti
    auxfunca=0.d0
    ss=0.d0
    select case(ndim)    
    case(4)! usi uti, vsi vti
        !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,nsujet_trial,i,fraili,ndim_Ind,frailij)&
        !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do ii=1,nsimu
                auxfunca=func2(func,fraili(ii,3),fraili(ii,4),fraili(ii,1),fraili(ii,2),nsimu,0.d0,frailij,ndim_Ind,nsujet_trial,i)
                ss=ss+auxfunca
                !!print*,"ss",ss (func,vsi,vti,ui,uti,nsimu,mu1,vc1,ndim,n,i)
            end do
        !$OMP END PARALLEL DO
    end select
    
    !  !print*,"suis lààààà 3"
    !  stop
    MC_MultInd_Essai_Cor=ss/dble(nsimu)
    !!print*,"ss",ss,"dble(nsimu)",dble(nsimu),"MC_MultInd_Essai_Cor",MC_MultInd_Essai_Cor
    !stop
    deallocate(vi,usim,vc,fraili,frailij,vc1)
    return
  end function MC_MultInd_Essai_Cor
  
end module monteCarlosMult_Gaus

!==================================================================
!=======module pour la fonction gauss hermite multidimension=======
!================================================================== 

module GaussHermi_mult
 ! function fortran pour le calcul integral
 implicit none
 
 contains
 
! funcSurrNN est l'integrant, herm le resultat de l integrale sur -infty , +infty
! l'appel de cette fonction se fait en initialisant le vecteur frail a x(1)
! ie do j=1,size(frail)
!        frail[j]=x(1)
!     end do
! initialiser k a size(frail)   

recursive function gaussHermMult(func,frail1,frail,i,k,x,w,inc) result(herm)
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
   
   ! bloc interface pour la definition de la fonction func
    interface
    ! func pour surrogate cas de deux variables normalement distribuees
     double precision function func(frail,frailst,i)
     integer,intent(inout):: i
     double precision,intent(in),dimension(:)::frail
     double precision,intent(in),dimension(2)::frailst !pour les effets aleatoire vs et vt
     end function func 
    end interface

   ! fin declaration et debut du programme
   !print*,"******je suis dans la subroutine gaussHermMult ligne 51:"
   
   n=size(frail)
   npoint=size(x)
   if (k.eq.1) then
     s=0
     do l1=1,npoint
       frail(n)=x(l1)
       s=s+w(l1)*func(frail1,frail,i)
       inc=inc+1.d0
     end do
     herm=s
   else
     s=0.d0
     do lk=1,npoint
       frail(n-k+1)=x(lk)
       s=s+w(lk)*gaussHermMult(func,frail1,frail,i,k-1,x,w,inc)
     end do
     herm=s
   end if
   ! ici on complete avec le produit matrice de cholesky de 2^r/2 de l'adaptative
   if(adaptative .and. inc .eq.(npoint**n)) then
     !print*," sub routine gaussHermMult, je suis dans le if fichier integrale_mult_scl.f90"
     ! stop
   end if
 end function gaussHermMult  

recursive function gaussHermMultGen(func,frail,k,x,w,inc,i) result(herm) 
   ! quadrature cas general avec la fonction func
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! frail: vecteur donc la taille correspond au nombre d'integrale
   ! k permet de gerer la recursivite et vaut initialement la taille de frail
   ! x et w: vecteur des poids et points de quadrature definis dans le fichier Adonnees.f90
   ! inc: un increment pour le controle, vaut 0 initialementt
   ! i: trial courant dans le quel on effectue le calcul integrale
   
   use var_surrogate, only:adaptative
   !use comon, only:invBi_cholDet
   
   implicit none
   
   integer, intent(in)::i
   integer, intent(inout)::k
   double precision, intent(inout)::inc
   integer :: l1,lk,n,npoint,k2
   double precision,dimension(:),intent(inout) :: frail
   double precision,dimension(:),intent(in) ::x,w
   double precision ::herm,s
   
   ! bloc interface pour la definition de la fonction func
    interface
    ! double precision function funcSurrNN(arg)
    !  double precision,dimension(2), intent(in):: arg
    ! end function funcSurrNN
     
    ! fonction avec deux parametres
     double precision function func(arg,i)
      integer, intent(in)::i! i est le cluster courant pour le calcul integral
      double precision,dimension(:), intent(in):: arg
     end function func
    end interface
    
    
   ! fin declaration et debut du programme
   n=size(frail)
   npoint=size(x)
   if (k.eq.1) then
     s=0
     do l1=1,npoint
       frail(n)=x(l1)
       s=s+w(l1)*func(frail,i)
       inc=inc+1.d0
     end do
     herm=s
   else
     s=0.d0
     do lk=1,npoint
       frail(n-k+1)=x(lk)
       k2=k-1
       s=s+w(lk)*gaussHermMultGen(func,frail,k2,x,w,inc,i)
     end do
     herm=s
   end if
   ! ici on complete avec le produit matrice de cholesky de 2^r/2 de l'adaptative
   if(adaptative .and. inc .eq.(npoint**n)) then
     !print*," sub routine gaussHermMultGen, je suis dans le if fichier integrale_mult_scl.f90"
     ! stop
   end if
 end function gaussHermMultGen
 
 
 double precision function gauss_HermMult(func,func2,npoint,n) 
   ! quadrature a l'aide du produit des integrales de chaque individu du cluster
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! npoint: nombre de point de quadratures
   ! n: nombre de suijet dans le cluster courant
   
   use var_surrogate, only:adaptative
   use comon, only: lognormal
   
   implicit none
   
   integer ::k2
   integer, intent(in)::n,npoint
   double precision ::herm,I1!,gauss_HermMult ! en commentaire car le type du resultat est place avant la declaration de la fonction
    external gauherJ1_scl
   
   ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(arg,j)
            integer, intent(in)::j! j est la position de l'individu dans le cluster i
            double precision,intent(in):: arg ! fragilite wij
        end function func
     
        SUBROUTINE func2(func,ss,nnodes,position_i)
        ! ss = resultat de l'integrale
        ! nnodes = nombre de noeudds d'integration
        ! position_i = la position j de l'individu sur lequel on se trouve dans le cluster i  
            double precision,intent(out)::ss
            integer,intent(in)::nnodes,position_i
            interface
                double precision function func(arg,i)
                    integer, intent(in)::i! i est la position de l'individu sur lequel on se trouve
                    double precision, intent(in):: arg
                end function func
            end interface
        end subroutine func2

    end interface
    
   ! fin declaration et debut du programme
    I1=0.d0
    if(lognormal==1)then
        herm =0.d0
    !        !print*,n
    !    stop
        do k2=1,n
            !call gauherJ1_scl(func,I1,npoint,k2)
            call func2(func,I1,npoint,k2)
            !!print*,"I1=",I1
            herm=herm+dlog(I1)
        end do
    else
        herm =0.d0
        do k2=1,n
            call func2(func,I1,npoint,k2)
            !!print*,"I1=",I1
            herm=herm+dlog(I1)
        end do
    endif
    !!print*,herm !premier sujet
    !stop

    ! *******************************A voir********************************************
   ! ici on complete avec le produit matrice de cholesky de 2^r/2 de l'adaptative
   if(adaptative) then
          !print*," sub routine gaussHermMult, je suis dans le if fichier integrale_mult_scl.f90"
     ! stop
   end if
    gauss_HermMult=herm
    return
 end function gauss_HermMult
 
  double precision function gauss_HermMultA(func,npoint1,n) 
   ! quadrature a l'aide du produit des integrales de chaque individu du cluster
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! npoint: nombre de point de quadratures
   ! n: nombre de suijet dans le cluster courant
   
   use var_surrogate, only:adaptative
   use comon, only: lognormal,invBi_cholDet
   !$ use OMP_LIB
   
   implicit none
   
   integer ::k2
   integer, intent(in)::n,npoint1
   double precision ::herm,I1!,gauss_HermMultA
   
   ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(j,npoint1)
        ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
        ! j = individu j du cluster i
        ! npoint = nombre de point de quadrature
            integer,intent(in):: j,npoint1
        end function func
    end interface
    
   ! fin declaration et debut du programme
    I1=0.d0
    if(lognormal==1)then
        herm =1.d0
        !$OMP PARALLEL DO PRIVATE (k2,I1) firstprivate(n,npoint1) REDUCTION(*:herm)
            do k2=1,n
                I1=func(k2,npoint1)
                !herm=herm+dlog(I1)
                herm=herm*I1
            end do
        !$OMP END PARALLEL DO
    else
        herm =0.d0
        do k2=1,n
            I1=func(k2,npoint1)
            herm=herm+dlog(I1)
            !!print*,herm !premier sujet
            !stop
        end do
    endif

   ! ici on complete avec le produit matrice de cholesky de 2^r/2 de l'adaptative
   if(adaptative) then
     herm=herm*dsqrt(2.d0)**n*invBi_cholDet(1)  !2^(q/2), même determinant pour tous les individus
   end if
    gauss_HermMultA=herm
    return
 end function gauss_HermMultA
 
 double precision function gauss_HermMultA_surr(func,vsi,vti,ui,npoint1,n,i) 
   ! quadrature a l'aide du produit des integrales de chaque individu du cluster
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! npoint: nombre de point de quadratures
   ! n: nombre de suijet dans le cluster courant
   ! vsi= frailtie niveau essai associe a s
   ! vti= frailtie niveau essai associe a t
   ! i= cluster courant
   
   use var_surrogate, only:varcovinv,gamma_ui,& !adaptative,cdcts,nigts,estim_wij_chap
                           frailt_base,methodInt,nb_procs,methodInt& !alpha_ui,nigs,cdcs,nsim,theta2
                           ,nb_procs
   use comon, only: lognormal
   use Autres_fonctions, only:pos_proc_domaine
   !use mpi
   !$ use OMP_LIB
   
   implicit none
   
   integer ::k2,init_i,max_i,rang
   integer, intent(in)::n,npoint1,i
   double precision,intent(in)::vsi,vti,ui
   double precision ::herm,I1,c1,c2
   double precision, dimension(:,:),allocatable::m1,m3  
   double precision, dimension(:,:),allocatable::m
   
   ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,j,npoint1)
        ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
        ! j = individu j du cluster i
        ! npoint = nombre de point de quadrature
            integer,intent(in):: j,npoint1
            double precision,intent(in)::vsi,vti,ui
        end function func
    end interface
    
    ! estimation des effets aleatoire pour la pseudo adaptative
    ! if(adaptative .and. estim_wij_chap.eq.0) then ! on n'a pas encore estime les wij_chap
        ! I1=func(vsi,vti,k2,npoint1) ! ceci permet d'estimer les effets aleatoires au niveau individuelle mais on utilise pas le resultat du calcul integral
        ! goto 100
    ! endif

   ! fin declaration et debut du programme
    !!print*,"suis la=======================1"
    herm = 0.d0
    I1=1.d0
    if(lognormal==1)then
        herm =1.d0
        if(nb_procs==1) then !on fait du open MP car un seul processus
            ! rang=0
            !$OMP PARALLEL DO default(none) PRIVATE (k2,I1) shared(n,npoint1,vsi,vti,ui,methodInt) REDUCTION(*:herm)
                do k2=1,n
                    !!print*,"1"
                    ! if(methodInt==4)then !integration au niveau individuel par monte-carlo
                        ! I1=func(vsi,vti,ui,k2,nsim,0.d0,theta2)
                    ! else
                        I1=func(vsi,vti,ui,k2,npoint1)
                    ! endif
                    !!print*,"vsi",vsi
                    !!print*,"vti",vti
                    !!print*,"12"
                    !herm=herm+dlog(I1)
                    herm=herm*I1
                    !!print*,"I1=",I1
                    !pos_ui_chap=pos_ui_chap+1 !permet de compter le nombre cumule des sujets deja parcouru 
                end do
            !$OMP END PARALLEL DO
        else ! dans ce cas on va faire du MPI donc plus de openMP
            ! rang du processus courang
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
            ! on cherche les position initiale et finale pour le processus courant
            call pos_proc_domaine(n,nb_procs,rang,init_i,max_i)
            do k2=1,n
                ! if((k2<init_i).or.k2>max_i) then 
                    ! goto 1001 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                ! endif
                ! if(methodInt==4)then !integration au niveau individuel par monte-carlo
                    ! I1=func(vsi,vti,ui,k2,nsim,0.d0,theta2)
                ! else
                    I1=func(vsi,vti,ui,k2,npoint1)
                ! endif
                herm=herm*I1
                1001 continue
            end do
            ! on fait la reduction et redistribu le resultat a tous les procesus
            ! !call MPI_ALLREDUCE(herm,herm,1,MPI_DOUBLE_PRECISION,MPI_PROD,MPI_COMM_WORLD,code)
        endif
    else
        !print*,"quadrature par la loi gamma non disponible pour le modele complet de surrogacy",&
         !       "probleme de la listribution gamma bivariée"
    endif

   ! ici on complete avec le produit matrice de cholesky de 2^r/2 de l'adaptative
   ! if(adaptative) then
    ! !!print*,"suis là dans gauss_HermMultA_surr pour adaptative"
     ! herm=herm*dsqrt(2.d0)**n*invBi_cholDet(1)  !2^(q/2), même determinant pour tous les individus
     ! !print*,"invBi_cholDet(1)=", invBi_cholDet(1)
     ! goto 100
   ! end if
   
   ! calcul du terme I1(vsi,vti)
   if(methodInt==1)then ! integration par quadrature
       allocate(m(1,1),m1(1,2),m3(1,2))
        m1(1,1)=vsi
        m1(1,2)=vti
        !call multiJ(m1,varcovinv,1,2,2,m3) !produit matriciel, resultat dans m3
        !call multiJ(m3,transpose(m1),1,2,1,m) !produit matriciel, resultat dans m
        m3 = matmul(m1,varcovinv) !produit matriciel, resultat dans m3
        m = matmul(m3,transpose(m1)) !produit matriciel, resultat dans m
        if(frailt_base==1) then
            !!print*,"suis la gamma_ui",gamma_ui,"i=",i
            c1=-((ui**2.d0)/(2.d0*gamma_ui))-(1.d0/2.d0) * m(1,1)
            ! c2=ui*(nigs(i)+alpha_ui*cdcs(i))+nigts(i)*vsi+cdcts(i)*vti
        else
            c1=-(1.d0/2.d0) * m(1,1)
            !c2=nigts(i)*vsi+cdcts(i)*vti ! deja pris en compte dans le calcul integrale: "Integrale_Individuel()"
        endif
        
        c2=0.d0 ! car deja pris en compte dans "Integrale_Individuel()"
        I1=dexp(c1+c2)
        !!print*,"I1=",I1
        !!print*,"herm=",herm,"I1=",I1
        deallocate(m,m1,m3)
    endif
    
    ! cas integration par monte carlo niveau essai et quadrature au niveau individuel
    if(methodInt==2)then ! integration par monte carlo niveau essai et quadrature au niveau individuel
       ! allocate(m(1,1))
        ! m1(1,1)=vsi
        ! m1(1,2)=vti
        !call multiJ(m1,varcovinv,1,2,2,m3) !produit matriciel, resultat dans m3
        !call multiJ(m3,transpose(m1),1,2,1,m) !produit matriciel, resultat dans m
        ! m3 = matmul(m1,varcovinv) !produit matriciel, resultat dans m3
        ! m = matmul(m3,transpose(m1)) !produit matriciel, resultat dans m
        if(frailt_base==1) then
            !!print*,"suis la gamma_ui",gamma_ui,"i=",i
            ! c1=-((ui**2.d0)/(2.d0*gamma_ui))-(1.d0/2.d0) * m(1,1)
            ! c2=ui*(nigs(i)+alpha_ui*cdcs(i))+nigts(i)*vsi+cdcts(i)*vti
            ! !print*,sum(nigs)
            ! !print*,sum(cdcs)
            !stop
        else
            ! c1=-(1.d0/2.d0) * m(1,1)
            ! c2=nigts(i)*vsi+cdcts(i)*vti
        endif
        I1=1
        !I1=dexp(c2)
        !!print*,"I1=",I1
        !!print*,"herm=",herm,"I1=",I1
        ! deallocate(m)
    endif
    
    gauss_HermMultA_surr=I1*herm
    100 continue
    return
 end function gauss_HermMultA_surr
 
 !integration au niveau individuel par monte-carlo  et essai par guass hermite
  double precision function gauss_HermMultA_surr_MC(func,vsi,vti,ui,npoint1,n,i) 
   ! quadrature a l'aide du produit des integrales de chaque individu du cluster
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! npoint: nombre de point de quadratures
   ! n: nombre de suijet dans le cluster courant
   ! vsi= frailtie niveau essai associe a s
   ! vti= frailtie niveau essai associe a t
   ! i= cluster courant
   
   use var_surrogate, only:varcovinv,gamma_ui,& !adaptative,cdcts,nigts,estim_wij_chap
                           frailt_base,methodInt,nb_procs,nsim,theta2,methodInt& !alpha_ui,nigs,cdcs
                           ,nb_procs
   use comon, only: lognormal
   use Autres_fonctions, only:pos_proc_domaine
   !use mpi
   !$ use OMP_LIB
   
   implicit none
   
   integer ::k2,init_i,max_i,rang
   integer, intent(in)::n,npoint1,i
   double precision,intent(in)::vsi,vti,ui
   double precision ::herm,I1,c1,c2
   double precision, dimension(:,:),allocatable::m1,m3  
   double precision, dimension(:,:),allocatable::m
   
   ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,j,npoint1,mu1,vc1)
        ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
        ! j = individu j du cluster i
        ! npoint = nombre de point de quadrature
            integer,intent(in):: j,npoint1
            double precision,intent(in)::vsi,vti,ui,mu1,vc1
        end function func
    end interface
    
    ! estimation des effets aleatoire pour la pseudo adaptative
    ! if(adaptative .and. estim_wij_chap.eq.0) then ! on n'a pas encore estime les wij_chap
        ! I1=func(vsi,vti,k2,npoint1) ! ceci permet d'estimer les effets aleatoires au niveau individuelle mais on utilise pas le resultat du calcul integral
        ! goto 100
    ! endif

   ! fin declaration et debut du programme
    !!print*,"suis la=======================1"
    I1=1.d0
    if(lognormal==1)then
        herm =1.d0
        if(nb_procs==1) then !on fait du open MP car un seul processus
            ! rang=0
            !$OMP PARALLEL DO default(none) PRIVATE (k2,I1) shared(n,npoint1,vsi,vti,ui,methodInt,theta2,nsim) REDUCTION(*:herm)
                do k2=1,n
                    !!print*,"1"
                    I1=func(vsi,vti,ui,k2,nsim,0.d0,theta2)

                    !!print*,"vsi",vsi
                    !!print*,"vti",vti
                    !!print*,"12"
                    !herm=herm+dlog(I1)
                    herm=herm*I1
                    !!print*,"I1=",I1
                    !pos_ui_chap=pos_ui_chap+1 !permet de compter le nombre cumule des sujets deja parcouru 
                end do
            !$OMP END PARALLEL DO
        else ! dans ce cas on va faire du MPI donc plus de openMP
            ! rang du processus courang
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
            ! on cherche les position initiale et finale pour le processus courant
            call pos_proc_domaine(n,nb_procs,rang,init_i,max_i)
            do k2=1,n
                if((k2<init_i).or.k2>max_i) then 
                    goto 1007 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                endif
                
                I1=func(vsi,vti,ui,k2,nsim,0.d0,theta2)
                herm=herm*I1
                1007 continue
            end do
            ! on fait la reduction et redistribu le resultat a tous les procesus
            !call MPI_ALLREDUCE(herm,herm,1,MPI_DOUBLE_PRECISION,MPI_PROD,MPI_COMM_WORLD,code)
        endif
    else
        !print*,"quadrature par la loi gamma non disponible pour le modele complet de surrogacy",&
        !        "probleme de la listribution gamma bivariée"
    endif

   ! ici on complete avec le produit matrice de cholesky de 2^r/2 de l'adaptative
   ! if(adaptative) then
    ! !!print*,"suis là dans gauss_HermMultA_surr pour adaptative"
     ! herm=herm*dsqrt(2.d0)**n*invBi_cholDet(1)  !2^(q/2), même determinant pour tous les individus
     ! !print*,"invBi_cholDet(1)=", invBi_cholDet(1)
     ! goto 110
   ! end if
   
   ! calcul du terme I1(vsi,vti)
   if(methodInt==1)then ! integration par quadrature
       allocate(m(1,1),m1(1,2),m3(1,2))
        m1(1,1)=vsi
        m1(1,2)=vti
        !call multiJ(m1,varcovinv,1,2,2,m3) !produit matriciel, resultat dans m3
        !call multiJ(m3,transpose(m1),1,2,1,m) !produit matriciel, resultat dans m
        m3 = matmul(m1,varcovinv) !produit matriciel, resultat dans m3
        m = matmul(m3,transpose(m1)) !produit matriciel, resultat dans m
        if(frailt_base==1) then
            !!print*,"suis la gamma_ui",gamma_ui,"i=",i
            c1=-((ui**2.d0)/(2.d0*gamma_ui))-(1.d0/2.d0) * m(1,1)
            ! c2=ui*(nigs(i)+alpha_ui*cdcs(i))+nigts(i)*vsi+cdcts(i)*vti
        else
            c1=-(1.d0/2.d0) * m(1,1)
            !c2=nigts(i)*vsi+cdcts(i)*vti ! deja pris en compte dans le calcul integrale: "Integrale_Individuel()"
        endif
        
        c2=0.d0 ! car deja pris en compte dans "Integrale_Individuel()"
        I1=dexp(c1+c2)
        !!print*,"I1=",I1
        !!print*,"herm=",herm,"I1=",I1
        deallocate(m,m1,m3)
    endif
    
    ! cas integration par monte carlo niveau essai et quadrature au niveau individuel
    if(methodInt==2)then ! integration par monte carlo niveau essai et quadrature au niveau individuel
       ! allocate(m(1,1))
        ! m1(1,1)=vsi
        ! m1(1,2)=vti
        !call multiJ(m1,varcovinv,1,2,2,m3) !produit matriciel, resultat dans m3
        !call multiJ(m3,transpose(m1),1,2,1,m) !produit matriciel, resultat dans m
        ! m3 = matmul(m1,varcovinv) !produit matriciel, resultat dans m3
        ! m = matmul(m3,transpose(m1)) !produit matriciel, resultat dans m
        if(frailt_base==1) then
            !!print*,"suis la gamma_ui",gamma_ui,"i=",i
            ! c1=-((ui**2.d0)/(2.d0*gamma_ui))-(1.d0/2.d0) * m(1,1)
            ! c2=ui*(nigs(i)+alpha_ui*cdcs(i))+nigts(i)*vsi+cdcts(i)*vti
            ! !print*,sum(nigs)
            ! !print*,sum(cdcs)
            !stop
        else
            ! c1=-(1.d0/2.d0) * m(1,1)
            ! c2=nigts(i)*vsi+cdcts(i)*vti
        endif
        I1=1
        !I1=dexp(c2)
        !!print*,"I1=",I1
        !!print*,"herm=",herm,"I1=",I1
        ! deallocate(m)
    endif
    
    gauss_HermMultA_surr_MC=I1*herm
    110 continue
    return
 end function gauss_HermMultA_surr_MC
 
 ! computation of the integrale using gaussian-Hermite quadrature for the joint frailty-copula model
 double precision function gauss_Herm_copula_Int(func,nnodes,ndim,nsujet_trial,i)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! func: fonction a integrer au niveau individuel
    ! nnodes: nombre de point de quadrature
    ! ndim= dimension de l'integrale 2 ou 3 integrations?
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    
    use var_surrogate, only: adaptative,xx1,ww1,invBi_chol_Essai,ui_chap_Essai,&
        invBi_cholDet_Essai,nb_procs
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use fonction_A_integrer, only:multiJ
    use Autres_fonctions, only:pos_proc_domaine
    !use mpi
    !$ use OMP_LIB
    
    
    implicit none
    integer ::ii,jj,npg,kk,cpt,init_i,max_i,code,erreur,rang
    integer,intent(in):: ndim,nnodes,nsujet_trial,i
    double precision::ss1,ss2,auxfunca,ss
    double precision, dimension(ndim)::xxl,m 
    double precision,dimension(ndim,ndim)::invBi_chol_Essai_k ! pour recuperer les matrice B_k dans le vecteur des matrices B des essais K 
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,ig,nsujet_trial)
            ! vsi= frailtie niveau essai associe a s
            ! vti= frailtie niveau essai associe a t
            ! ui = random effect associated xith the baseline hazard
            ! ig = current cluster
            ! nsujet_trial = number of subjects in the current trial    
            IMPLICIT NONE
            integer,intent(in):: ig, nsujet_trial
            double precision,intent(in)::vsi,vti,ui
        end function func
    end interface
    
    npg=nnodes    
    auxfunca=0.d0
    ss=0.d0
    ss1=0.d0
    ss2=0.d0
        if(adaptative) then
        ! je recupere la matrice B_i de l'essai i dans le vecteur des B_i
            cpt=((i-1)*ndim**2)+1 !ceci permet de parcourir le vecteur invBi_chol_Essai en fonction du cluster sur lequel on se trouve
            do jj=1,ndim     
                do ii=1,ndim
                    invBi_chol_Essai_k(ii,jj)=invBi_chol_Essai(cpt) ! en effet le vecteur "invBi_chol_Essai" a ete rempli par des matrices colonne apres colonne
                    cpt=cpt+1
                enddo
            enddo         
        endif
    
    if(ndim.eq.2) then
        !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl) firstprivate(auxfunca) SHARED(npg,nsujet_trial,i,xx1,ww1,&
        !$OMP invBi_chol_Essai_k,ndim,ui_chap_Essai,adaptative)&
        !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do ii=1,npg
                ss1=0.d0
                do jj=1,npg
                    xxl(1)=xx1(ii)
                    xxl(2)=xx1(jj)
                    !changement de variable en cas de quadrature adaptative
                    if(adaptative) then ! on effectue le changement de variable
                        m=matmul(invBi_chol_Essai_k,xxl)
                        xxl=ui_chap_Essai(i,1:2)+dsqrt(2.d0)*m        
                        !!print*,"xxl",xxl
                    end if
                    auxfunca=func(xxl(1),xxl(2),0.d0,i,nsujet_trial)
                    ss1 = ss1+ww1(jj)*(auxfunca)
                end do
                ss = ss+ww1(ii)*ss1
            end do
      !$OMP END PARALLEL DO
    else ! cas de 3 points
        if(nb_procs==1) then !on fait du open MP car un seul processus
            rang=0
            !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl,kk,ss2) firstprivate(auxfunca) &
            !$OMP SHARED(npg,nsujet_trial,i,xx1,ww1,&
            !$OMP invBi_chol_Essai_k,ndim,ui_chap_Essai,adaptative)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do kk=1,npg
                ss2=0.d0
                do ii=1,npg
                    ss1=0.d0
                    do jj=1,npg    
                        xxl(2)=xx1(jj)        
                        xxl(1)=xx1(ii)                    
                        xxl(3)=xx1(kk)
                        !changement de variable en cas de quadrature adaptative
                        if(adaptative) then ! on effectue le changement de variable
                            m=matmul(invBi_chol_Essai_k,xxl)
                            xxl=ui_chap_Essai(i,1:3)+dsqrt(2.d0)*m        
                        end if                      
                        auxfunca=func(xxl(1),xxl(2),xxl(3),i,nsujet_trial)
                        ss1 = ss1+ww1(jj)*(auxfunca)
                    end do
                    ss2 = ss2+ww1(ii)*ss1
                end do
                ss = ss+ww1(kk)*ss2
            end do
            !$OMP END PARALLEL DO
        else ! dans ce cas on va faire du MPI
            ! rang du processus courang
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
            ! on cherche les position initiale et finale pour le processus courant
            call pos_proc_domaine(npg,nb_procs,rang,init_i,max_i)
            ! !print*,"nb_procs=",nb_procs,"rang=",rang
            ! !print*,"init_i,max_i=",init_i,max_i
            
            do kk=1,npg
                ! if((kk<init_i).or.kk>max_i) then 
                    ! goto 1000 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                ! endif
                ss2=0.d0
                do ii=1,npg
                    ss1=0.d0
                    do jj=1,npg    
                        xxl(2)=xx1(jj)        
                        xxl(1)=xx1(ii)                    
                        xxl(3)=xx1(kk)
                        
                        !changement de variable en cas de quadrature adaptative
                        if(adaptative) then ! on effectue le changement de variable
                            m=matmul(invBi_chol_Essai_k,xxl)
                            xxl=ui_chap_Essai(i,1:3)+dsqrt(2.d0)*m        
                        end if
                        
                        !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                        auxfunca=func(xxl(1),xxl(2),xxl(3),i,nsujet_trial)
                        !!print*,"12"
                        ss1 = ss1+ww1(jj)*(auxfunca)
                        ! !print*,ww1(jj)
                        ! stop
                        !!print*,"13"
                        !!print*,"Integrant niveau essai=",auxfunca
                    end do
                    ss2 = ss2+ww1(ii)*ss1
                    ! !print*,ww1(ii)
                    ! stop
                end do
                ss = ss+ww1(kk)*ss2
                ! !print*,ww1(kk)
                ! stop
                1000 continue
            end do
            ! on fait la reduction et redistribu le resultat a tous les procesus
            ! !call MPI_ALLREDUCE(ss,ss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
        endif
    end if
    
    ! !print*,"rang=",rang,"voile la valeur de ton calcul integral",ss
    ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
    !!print*,"ss=",ss
    if(adaptative) then
        ss=ss*2.d0**(dble(ndim)/2.d0)*invBi_cholDet_Essai(i) 
        ! if(frailt_base==0) ss=ss*2.d0**(dble(ndim)/2.d0)*invBi_cholDet_Essai(i)  
    end if
    !!print*,"ss_transformé=",ss
    
    gauss_Herm_copula_Int=ss
    101 continue
    return
  end function gauss_Herm_copula_Int
  
 ! calcul de l'integral au niveau essai pour le modele surrogate final par quadrature gaussienne
 double precision function gauss_HermMultInd_Essai(func,func2,nnodes,ndim,nsujet_trial,i)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! func: fonction a integrer au niveau individuel
    ! nnodes: nombre de point de quadrature
    ! ndim= dimension de l'integrale 2 ou 3 integrations?
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    
    use var_surrogate, only: adaptative,xx1,ww1,posind_i,invBi_chol_Essai,ui_chap_Essai,&
        invBi_cholDet_Essai,frailt_base,nb_procs
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use fonction_A_integrer, only:multiJ
    use Autres_fonctions, only:pos_proc_domaine
    !use mpi
    !$ use OMP_LIB
    
    
    implicit none
    integer ::ii,jj,npg,kk,cpt,init_i,max_i,rang
    integer,intent(in):: ndim,nnodes,nsujet_trial,i
    !double precision,dimension(1:nnodes) ::xx1,ww1
    double precision::ss1,ss2,auxfunca,ss
    double precision, dimension(ndim)::xxl,m !vecteur qui contiendra à chaque fois les points de quadrature
    double precision,dimension(ndim,ndim)::invBi_chol_Essai_k ! pour recuperer les matrice B_k dans le vecteur des matrices B des essais K
    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,j,npoint1)
        ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
        ! j = individu j du cluster i
        ! npoint1 = nombre de point de quadrature
            integer,intent(in):: j,npoint1
            double precision,intent(in)::vsi,vti,ui
        end function func
        
        double precision function func2(func,vsi,vti,ui,npoint1,n,i)
            integer, intent(in)::n,npoint1,i
            double precision,intent(in)::vsi,vti,ui
            
            interface
                double precision function func(vsi,vti,ui,j,npoint1)
                    integer,intent(in):: j,npoint1
                    double precision,intent(in)::vsi,vti,ui
                end function func
            end interface
        end function func2
    end interface
    !!print*,"suis la=======================1",adaptative, estim_wij_chap
    ! if(adaptative .and. estim_wij_chap.eq.0) then ! on n'a pas encore estime les wij_chap
        ! auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i) ! ceci permet d'estimer les effets aleatoires au niveau individuelle mais on utilise pas le resultat du calcul integral
        ! goto 101
    ! endif
    
    npg=nnodes    
    !integration sur vsi et vti
    auxfunca=0.d0
    ss=0.d0
    ss1=0.d0
    ss2=0.d0
    ! if(frailt_base==0) then
        if(adaptative) then
        ! je recupere la matrice B_i de l'essai i dans le vecteur des B_i
            cpt=((i-1)*ndim**2)+1 !ceci permet de parcourir le vecteur invBi_chol_Essai en fonction du cluster sur lequel on se trouve
            do jj=1,ndim     
                do ii=1,ndim
                    invBi_chol_Essai_k(ii,jj)=invBi_chol_Essai(cpt) ! en effet le vecteur "invBi_chol_Essai" a ete rempli par des matrices colonne apres colonne
                    cpt=cpt+1
                    !!print*,"invBi_chol_Essai_k(cpt)=",invBi_chol_Essai_k(cpt)
                enddo
            enddo        
            
            !!print*,"invBi_chol_Essai=",invBi_chol_Essai
            !!print*,"invBi_chol_Essai_k=",invBi_chol_Essai_k
            
        endif
    ! endif
    
    !!print*,"ndim=",ndim
    
    if(ndim.eq.2) then
      !!print*,'je suis la'
        !ii=0
        !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl) firstprivate(auxfunca) SHARED(npg,nnodes,nsujet_trial,i,xx1,ww1,&
        !$OMP invBi_chol_Essai_k,ndim,ui_chap_Essai,adaptative,posind_i)&
        !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do ii=1,npg
                ss1=0.d0
                do jj=1,npg
                    xxl(1)=xx1(ii)
                    xxl(2)=xx1(jj)
                    !changement de variable en cas de quadrature adaptative
                    if(adaptative) then ! on effectue le changement de variable
                        m=matmul(invBi_chol_Essai_k,xxl)
                        xxl=ui_chap_Essai(i,1:2)+dsqrt(2.d0)*m        
                        !!print*,"xxl",xxl
                    end if
                    !stop
                    !xxl(1)=xx1(ii)
                    !xxl(2)=xx1(jj)
                    
                    !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                    auxfunca=func2(func,xxl(1),xxl(2),0.d0,nnodes,nsujet_trial,i)
                    !!print*,"12"
                    ss1 = ss1+ww1(jj)*(auxfunca)
                    !!print*,"13"
                    !!print*,"Integrant niveau essai=",auxfunca
                end do
                ss = ss+ww1(ii)*ss1
            end do
      !$OMP END PARALLEL DO
    else ! cas de 3 points
        if(nb_procs==1) then !on fait du open MP car un seul processus
            rang=0
            !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl,kk,ss2) firstprivate(auxfunca) &
            !$OMP SHARED(npg,nnodes,nsujet_trial,i,xx1,ww1,&
            !$OMP invBi_chol_Essai_k,ndim,ui_chap_Essai,adaptative,posind_i,frailt_base)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do kk=1,npg
                ss2=0.d0
                do ii=1,npg
                    ss1=0.d0
                    do jj=1,npg    
                        xxl(2)=xx1(jj)        
                        xxl(1)=xx1(ii)                    
                        xxl(3)=xx1(kk)
                        
                        !changement de variable en cas de quadrature adaptative
                        if(adaptative) then ! on effectue le changement de variable
                            m=matmul(invBi_chol_Essai_k,xxl)
                            xxl=ui_chap_Essai(i,1:3)+dsqrt(2.d0)*m        
                            !!print*,"xxl",xxl
                        end if
                        
                        !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                        auxfunca=func2(func,xxl(1),xxl(2),xxl(3),nnodes,nsujet_trial,i)
                        !!print*,"12"
                        ss1 = ss1+ww1(jj)*(auxfunca)
                        ! !print*,ww1(jj)
                        ! stop
                        !!print*,"13"
                        !!print*,"Integrant niveau essai=",auxfunca
                    end do
                    ss2 = ss2+ww1(ii)*ss1
                    ! !print*,ww1(ii)
                    ! stop
                end do
                ss = ss+ww1(kk)*ss2
                ! !print*,ww1(kk)
                ! stop
            end do
            !$OMP END PARALLEL DO
        else ! dans ce cas on va faire du MPI
            ! rang du processus courang
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
            ! on cherche les position initiale et finale pour le processus courant
            call pos_proc_domaine(npg,nb_procs,rang,init_i,max_i)
            ! !print*,"nb_procs=",nb_procs,"rang=",rang
            ! !print*,"init_i,max_i=",init_i,max_i
            
            do kk=1,npg
                ! if((kk<init_i).or.kk>max_i) then 
                    ! goto 1000 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                ! endif
                ss2=0.d0
                do ii=1,npg
                    ss1=0.d0
                    do jj=1,npg    
                        xxl(2)=xx1(jj)        
                        xxl(1)=xx1(ii)                    
                        xxl(3)=xx1(kk)
                        
                        !changement de variable en cas de quadrature adaptative
                        if(adaptative) then ! on effectue le changement de variable
                            m=matmul(invBi_chol_Essai_k,xxl)
                            xxl=ui_chap_Essai(i,1:3)+dsqrt(2.d0)*m        
                        end if
                        
                        !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                        auxfunca=func2(func,xxl(1),xxl(2),xxl(3),nnodes,nsujet_trial,i)
                        !!print*,"12"
                        ss1 = ss1+ww1(jj)*(auxfunca)
                        ! !print*,ww1(jj)
                        ! stop
                        !!print*,"13"
                        !!print*,"Integrant niveau essai=",auxfunca
                    end do
                    ss2 = ss2+ww1(ii)*ss1
                    ! !print*,ww1(ii)
                    ! stop
                end do
                ss = ss+ww1(kk)*ss2
                ! !print*,ww1(kk)
                ! stop
                1000 continue
            end do
            ! on fait la reduction et redistribu le resultat a tous les procesus
            ! !call MPI_ALLREDUCE(ss,ss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
        endif
    end if
    
    ! !print*,"rang=",rang,"voile la valeur de ton calcul integral",ss
    ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
    !!print*,"ss=",ss
    if(adaptative) then
        ss=ss*2.d0**(dble(ndim)/2.d0)*invBi_cholDet_Essai(i) 
        ! if(frailt_base==0) ss=ss*2.d0**(dble(ndim)/2.d0)*invBi_cholDet_Essai(i)  
    end if
    !!print*,"ss_transformé=",ss
    
    gauss_HermMultInd_Essai=ss
    101 continue
    return
  end function gauss_HermMultInd_Essai
  
 ! calcul de l'integral au niveau essai pour le modele surrogate final par quadrature gaussienne sachant qu'au niveau individuel on va integrer oar Monte-carlo
 double precision function gauss_HermMultInd_Essai_MC(func,func2,nnodes,ndim,nsujet_trial,i)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! func: fonction a integrer au niveau individuel
    ! nnodes: nombre de point de quadrature
    ! ndim= dimension de l'integrale 2 ou 3 integrations?
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    
    use var_surrogate, only: adaptative,xx1,ww1,posind_i,invBi_chol_Essai,ui_chap_Essai,&
        invBi_cholDet_Essai,frailt_base,nb_procs
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use fonction_A_integrer, only:multiJ
    use Autres_fonctions, only:pos_proc_domaine
    !use mpi
    !$ use OMP_LIB
    
    
    implicit none
    integer ::ii,jj,npg,kk,cpt,init_i,max_i,rang
    integer,intent(in):: ndim,nnodes,nsujet_trial,i
    !double precision,dimension(1:nnodes) ::xx1,ww1
    double precision::ss1,ss2,auxfunca,ss
    double precision, dimension(ndim)::xxl,m !vecteur qui contiendra à chaque fois les points de quadrature
    double precision,dimension(ndim,ndim)::invBi_chol_Essai_k ! pour recuperer les matrice B_k dans le vecteur des matrices B des essais K
    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,j,npoint1,mu1,vc1)
        ! fonction a integrer: cas des effets aleatoires niveau essai et individuel normalement distribues: ici on integre seulement par MC
        ! j = individu j du cluster i
        ! npoint1 = nombre de point de quadrature
        ! mu1= moyenne du frailty
        ! vc1= variance du frailty
            integer,intent(in):: j,npoint1
            double precision,intent(in)::vsi,vti,ui,mu1,vc1
        end function func
        
        double precision function func2(func,vsi,vti,ui,npoint1,n,i)
            integer, intent(in)::n,npoint1,i
            double precision,intent(in)::vsi,vti,ui
            
            interface
                double precision function func(vsi,vti,ui,j,npoint1,mu1,vc1)
                    integer,intent(in):: j,npoint1
                    double precision,intent(in)::vsi,vti,ui,mu1,vc1
                end function func
            end interface
        end function func2
    end interface
    !!print*,"suis la=======================1",adaptative, estim_wij_chap
    ! if(adaptative .and. estim_wij_chap.eq.0) then ! on n'a pas encore estime les wij_chap
        ! auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i) ! ceci permet d'estimer les effets aleatoires au niveau individuelle mais on utilise pas le resultat du calcul integral
        ! goto 101
    ! endif
    
    npg=nnodes    
    !integration sur vsi et vti
    auxfunca=0.d0
    ss=0.d0
    ss1=0.d0
    ss2=0.d0
    ! if(frailt_base==0) then
        if(adaptative) then
        ! je recupere la matrice B_i de l'essai i dans le vecteur des B_i
            cpt=((i-1)*ndim**2)+1 !ceci permet de parcourir le vecteur invBi_chol_Essai en fonction du cluster sur lequel on se trouve
            do jj=1,ndim     
                do ii=1,ndim
                    invBi_chol_Essai_k(ii,jj)=invBi_chol_Essai(cpt) ! en effet le vecteur "invBi_chol_Essai" a ete rempli par des matrices colonne apres colonne
                    cpt=cpt+1
                    !!print*,"invBi_chol_Essai_k(cpt)=",invBi_chol_Essai_k(cpt)
                enddo
            enddo        
            
            !!print*,"invBi_chol_Essai=",invBi_chol_Essai
            !!print*,"invBi_chol_Essai_k=",invBi_chol_Essai_k
            
        endif
    ! endif
    
    !!print*,"ndim=",ndim
    
    if(ndim.eq.2) then
      !!print*,'je suis la'
        !ii=0
        !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl) firstprivate(auxfunca) SHARED(npg,nnodes,nsujet_trial,i,xx1,ww1,&
        !$OMP invBi_chol_Essai_k,ndim,ui_chap_Essai,adaptative,posind_i)&
        !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do ii=1,npg
                ss1=0.d0
                do jj=1,npg
                    xxl(1)=xx1(ii)
                    xxl(2)=xx1(jj)
                    !changement de variable en cas de quadrature adaptative
                    if(adaptative) then ! on effectue le changement de variable
                        m=matmul(invBi_chol_Essai_k,xxl)
                        xxl=ui_chap_Essai(i,1:2)+dsqrt(2.d0)*m        
                        !!print*,"xxl",xxl
                    end if
                    !stop
                    !xxl(1)=xx1(ii)
                    !xxl(2)=xx1(jj)
                    
                    !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                    auxfunca=func2(func,xxl(1),xxl(2),0.d0,nnodes,nsujet_trial,i)
                    !!print*,"12"
                    ss1 = ss1+ww1(jj)*(auxfunca)
                    !!print*,"13"
                    !!print*,"Integrant niveau essai=",auxfunca
                end do
                ss = ss+ww1(ii)*ss1
            end do
      !$OMP END PARALLEL DO
    else ! cas de 3 points
        if(nb_procs==1) then !on fait du open MP car un seul processus
            rang=0
            !$OMP PARALLEL DO default(none) PRIVATE (ii,jj,ss1,m,xxl,kk,ss2) firstprivate(auxfunca) &
            !$OMP SHARED(npg,nnodes,nsujet_trial,i,xx1,ww1,&
            !$OMP invBi_chol_Essai_k,ndim,ui_chap_Essai,adaptative,posind_i,frailt_base)&
            !$OMP    REDUCTION(+:ss) SCHEDULE(Dynamic,1)
            do kk=1,npg
                ss2=0.d0
                do ii=1,npg
                    ss1=0.d0
                    do jj=1,npg    
                        xxl(2)=xx1(jj)        
                        xxl(1)=xx1(ii)                    
                        xxl(3)=xx1(kk)
                        
                        !changement de variable en cas de quadrature adaptative
                        if(adaptative) then ! on effectue le changement de variable
                            m=matmul(invBi_chol_Essai_k,xxl)
                            xxl=ui_chap_Essai(i,1:3)+dsqrt(2.d0)*m        
                            !!print*,"xxl",xxl
                        end if
                        
                        !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                        auxfunca=func2(func,xxl(1),xxl(2),xxl(3),nnodes,nsujet_trial,i)
                        !!print*,"12"
                        ss1 = ss1+ww1(jj)*(auxfunca)
                        ! !print*,ww1(jj)
                        ! stop
                        !!print*,"13"
                        !!print*,"Integrant niveau essai=",auxfunca
                    end do
                    ss2 = ss2+ww1(ii)*ss1
                    ! !print*,ww1(ii)
                    ! stop
                end do
                ss = ss+ww1(kk)*ss2
                ! !print*,ww1(kk)
                ! stop
            end do
            !$OMP END PARALLEL DO
        else ! dans ce cas on va faire du MPI
            ! rang du processus courang
            !call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
            ! on cherche les position initiale et finale pour le processus courant
            call pos_proc_domaine(npg,nb_procs,rang,init_i,max_i)
            ! !print*,"nb_procs=",nb_procs,"rang=",rang
            ! !print*,"init_i,max_i=",init_i,max_i
            
            do kk=1,npg
                ! if((kk<init_i).or.kk>max_i) then 
                    ! goto 1000 ! pour dire le processus ne considere pas cet itteration car n'appartient pas a son domaine
                ! endif
                ss2=0.d0
                do ii=1,npg
                    ss1=0.d0
                    do jj=1,npg    
                        xxl(2)=xx1(jj)        
                        xxl(1)=xx1(ii)                    
                        xxl(3)=xx1(kk)
                        
                        !changement de variable en cas de quadrature adaptative
                        if(adaptative) then ! on effectue le changement de variable
                            m=matmul(invBi_chol_Essai_k,xxl)
                            xxl=ui_chap_Essai(i,1:3)+dsqrt(2.d0)*m        
                        end if
                        
                        !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                        auxfunca=func2(func,xxl(1),xxl(2),xxl(3),nnodes,nsujet_trial,i)
                        !!print*,"12"
                        ss1 = ss1+ww1(jj)*(auxfunca)
                        ! !print*,ww1(jj)
                        ! stop
                        !!print*,"13"
                        !!print*,"Integrant niveau essai=",auxfunca
                    end do
                    ss2 = ss2+ww1(ii)*ss1
                    ! !print*,ww1(ii)
                    ! stop
                end do
                ss = ss+ww1(kk)*ss2
                ! !print*,ww1(kk)
                ! stop
                1000 continue
            end do
            ! on fait la reduction et redistribu le resultat a tous les procesus
            ! !call MPI_ALLREDUCE(ss,ss,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)
        endif
    end if
    
    ! !print*,"rang=",rang,"voile la valeur de ton calcul integral",ss
    ! !call MPI_ABORT(MPI_COMM_WORLD,erreur,code)
    !!print*,"ss=",ss
    if(adaptative) then
        ss=ss*2.d0**(dble(ndim)/2.d0)*invBi_cholDet_Essai(i) 
        ! if(frailt_base==0) ss=ss*2.d0**(dble(ndim)/2.d0)*invBi_cholDet_Essai(i)  
    end if
    !!print*,"ss_transformé=",ss
    
    gauss_HermMultInd_Essai_MC=ss
    101 continue
    return
  end function gauss_HermMultInd_Essai_MC
  
   ! calcul de l'integral au niveau individuel pour le modele surrogate final par quadrature gaussienne
 double precision function gauss_HermMultInd_cor(func,vsi,vti,ui,uti,nnodes,ndim,nsujet_trial,i)
    ! dans cette fonction on fait une quadrature adaptative ou non pour les deux effets aleatoire vsi et vti
    ! func: fonction a integrer au niveau individuel
    ! nnodes: nombre de point de quadrature
    ! ndim= dimension de l'integrale 2 ou 3 integrations?
    ! nsujet_trial= nombre de sujets dans le cluster courant
    ! i= cluster courant
    
    use var_surrogate, only: frailt_base,&
        nigs,cdcs,nigts,cdcts,pi !adaptative,xx1,ww1,estim_wij_chap,posind_i,invBi_chol_Essai,ui_chap_Essai,determinant,invBi_cholDet_Essai
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use fonction_A_integrer, only:multiJ
    !$ use OMP_LIB
    
    implicit none
    integer :: k2
    double precision,intent(in)::vsi,vti,ui,uti
    integer,intent(in):: ndim,nnodes,nsujet_trial,i
    !double precision,dimension(1:nnodes) ::xx1,ww1
    double precision::ss,herm,I1,c2
    !double precision, dimension(ndim)::xxl,m,xx !vecteur qui contiendra à chaque fois les points de quadrature
    !double precision,dimension(ndim,ndim)::invBi_chol_Essai_k ! pour recuperer les matrice B_k dans le vecteur des matrices B des essais K
    !double precision, external::gauss_HermMultA_surr    
    
    ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(vsi,vti,ui,uti,nnodes,ndim,j)
            integer,intent(in):: j,ndim,nnodes
            double precision,intent(in)::vsi,vti,ui,uti
        end function func
        
    end interface
    !!print*,"suis la=======================1",adaptative, estim_wij_chap
    ! if(adaptative .and. estim_wij_chap.eq.0) then ! on n'a pas encore estime les wij_chap
        ! auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i) ! ceci permet d'estimer les effets aleatoires au niveau individuelle mais on utilise pas le resultat du calcul integral
        ! goto 101
    ! endif
    
    ! npg=nnodes    
    ! !integration sur vsi et vti
    ! auxfunca=0.d0
    ! ss=0.d0
    ! ss1=0.d0
    ! ss2=0.d0
    
    ! !pas encore implementer la quadrature pseudoadaptative dans ce cas
    ! if(frailt_base==0) then
        ! if(adaptative) then
        ! ! je recupere la matrice B_i de l'individu i dans le vecteur des B_i, dans ce cas il s'agit de la matrice associe aux individus et plus aux essais
            ! cpt=((i-1)*ndim**2)+1 !ceci permet de parcourir le vecteur invBi_chol_Essai en fonction du cluster sur lequel on se trouve
            ! do jj=1,ndim     
                ! do ii=1,ndim
                    ! invBi_chol_Essai_k(ii,jj)=invBi_chol_Essai(cpt) ! en effet le vecteur "invBi_chol_Essai" a ete rempli par des matrices colonne apres colonne
                    ! cpt=cpt+1
                    ! !!print*,"invBi_chol_Essai_k(cpt)=",invBi_chol_Essai_k(cpt)
                ! enddo
            ! enddo        
            
            ! !!print*,"invBi_chol_Essai=",invBi_chol_Essai
            ! !!print*,"invBi_chol_Essai_k=",invBi_chol_Essai_k
            
        ! endif
    ! endif
    
    !!print*,"ndim=",ndim
    herm =1.d0
    if(ndim.eq.2) then
     ! !print*,'je suis la'
        !ii=0
        !$OMP PARALLEL DO default(none) PRIVATE (k2,ss) SHARED(nnodes,nsujet_trial,vsi,vti,ui,uti,ndim)&
        !$OMP    REDUCTION(*:herm) SCHEDULE(Dynamic,1)
            do k2=1,nsujet_trial
                ! do ii=1,npg
                    ! ss1=0.d0
                    ! do jj=1,npg
                        ! xxl(1)=xx1(ii)
                        ! xxl(2)=xx1(jj)
                        ! !changement de variable en cas de quadrature adaptative
                        ! if(adaptative) then ! on effectue le changement de variable
                            ! m=matmul(invBi_chol_Essai_k,xxl)                            
                            ! xxl=ui_chap_Essai(i,1:2)+dsqrt(2.d0)*m        
                            ! !!print*,"xxl",xxl
                        ! end if                        
                        ! !auxfunca=func2(func,xx1(ii),xx1(jj),nnodes,nsujet_trial,i)
                        ! auxfunca=func(vsi,vti,ui,uti,xxl(1),xxl(2),k2)
                        ! !!print*,"12"
                        ! ss1 = ss1+ww1(jj)*(auxfunca)
                        ! ! !print*,vsi,vti,ui,uti,xxl(1),xxl(2)
                        ! ! !print*,"Integrant niveau Individuel=",auxfunca
                        ! ! stop
                    ! end do
                    ! ss = ss+ww1(ii)*ss1
                ! end do
                ss=func(vsi,vti,ui,uti,nnodes,ndim,k2)
                herm=herm*ss
            enddo
      !$OMP END PARALLEL DO
    end if
    
    !calcul de I1
    if(frailt_base==1) then
        c2=ui*nigs(i)+uti*cdcs(i)+nigts(i)*vsi+cdcts(i)*vti
    else
        c2=nigts(i)*vsi+cdcts(i)*vti
    endif
    I1=dexp(c2)
     ! !print*,"trial=",i,"I=",herm*exp(-nsujet_trial*(log(2*pi)+log(determinant)/2.d0)),I1
     ! stop
    
    gauss_HermMultInd_cor=I1*herm
    101 continue
    return
  end function gauss_HermMultInd_cor
    
  double precision function monteCarlos_ind(func,mu,vc,nsim,vcdiag) 
   ! Monte carlo a l'aide du produit des integrales de chaque individu du cluster, compte tenu de l'independance des effets aleatoires
   ! func est la fonction a integrer, definie dans le module fonction_A_integrer (Integrant_scl.f90)
   ! mu: vecteur des moyenne des effets aleatoires=0
   ! vc: matrice diagonale des variances des effecs aleatoires
   ! nsim: nombre de simulation
   ! vcdiag: dit si la matrice vc est diagonale
   
    !use var_surrogate, only:adaptative
    !use comon, only:invBi_cholDet
    use monteCarlosMult_Gaus ! pour l'integrale par monte carlo   
   
    implicit none
   
    integer ::k2
    integer::n
    double precision ::mc
    double precision,intent(in),dimension(:,:)::vc
    double precision,intent(in),dimension(:)::mu
    double precision,dimension(1)::mu1
    double precision,dimension(1,1)::vc1
    double precision,dimension(3)::I1
    integer,intent(in)::nsim,vcdiag
   
   ! bloc interface pour la definition de la fonction func
    interface
        double precision function func(arg,j)
            integer, intent(in)::j! j est la position de l'individu dans le cluster i
            double precision,intent(in):: arg ! fragilite wij
        end function func
    end interface
    
   ! fin declaration et debut du programme
    I1=0.d0
    mc =1.d0
    n=size(mu)
    !!print*,"n=",n
    do k2=1,n
        mu1=mu(k2)
        vc1=vc(k2,k2)
        call monteCarlosMult_ind(func,mu1,vc1,nsim,vcdiag,k2,I1)
!        !print*,"I1=",I1
        mc=mc*I1(1)
    end do

    monteCarlos_ind=mc
    return
 end function monteCarlos_ind
 
 !===========================================
 !gauss hermite pour une dimension: viens de Agnieszca
 
 SUBROUTINE gauherJ1_scl(func,ss,nnodes,position_i)
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
            auxfunca = 0.d0
            methodGH=0
            if(methodGH.eq.0) then !non-adaptative
                do j=1,nnodes  
                    auxfunca=func(xx1(j),position_i)
                    ss = ss+ww1(j)*(auxfunca)                     
                end do
            endif    
        return
    
        END SUBROUTINE gauherJ1_scl
 
 !========Guass hermite pour deux et 3 dimensions====================================
 !==============================================================================
 
 subroutine gausshermite_2_3(func30,ss,npg,npoint,trial_i)
        !npoint nombre d'integrale: vaut 2 initialement
        !trial_i ordre du trial

        double precision,intent(out)::ss
        double precision::ss1,ss2
        double precision::auxfunca
        INTEGER::ii,jj,npg,kk
        integer,intent(in)::npoint,trial_i
        double precision,dimension(npg)::X,W
        double precision,dimension(npoint)::tab
        interface
          !double precision function func30(arg,ndim)
          !integer,intent(in)::ndim
          !double precision,dimension(ndim), intent(in):: arg
          !end function func30
          double precision function func30(arg,trial_i)
            integer, intent(in)::trial_i! trial_i est le cluster courant pour le calcul integral
            double precision,dimension(2), intent(in):: arg
          end function func30
        end interface

    if(npg.eq.10) then
    
        x(1)=-3.43615911884d0
        x(2)=-2.53273167423d0
        x(3)=-1.7566836493d0
        x(4)=-1.03661082979d0
        x(5)=-0.342901327224d0
        x(6)=0.342901327224d0
        x(7)=1.03661082979d0
        x(8)=1.7566836493d0
        x(9)=2.53273167423d0
        x(10)=3.43615911884d0
        
        w(1)=1.02545169137d0
        w(2)=0.820666126405d0
        w(3)=0.741441931944d0
        w(4)=0.703296323105d0
        w(5)=0.687081853951d0
        w(6)=0.687081853951d0
        w(7)=0.703296323105d0
        w(8)=0.741441931944d0
        w(9)=0.820666126405d0
        w(10)=1.02545169137d0
    endif


    if(npg.eq.15) then

        x(1)=-4.49999070731d0
        x(2)=-3.6699503734d0
        x(3)=-2.96716692791d0
        x(4)=-2.32573248617d0
        x(5)=-1.71999257519d0
        x(6)=-1.13611558521d0
        x(7)=-0.565069583256d0
        x(8)=0.d0
        x(9)=0.565069583256d0
        x(10)=1.13611558521d0
        x(11)=1.71999257519d0
        x(12)=2.32573248617d0
        x(13)=2.96716692791d0
        x(14)=3.6699503734d0
        x(15)=4.49999070731d0
    
        w(1)=0.948368970828d0
        w(2)=0.748607366017d0
        w(3)=0.666166005109d0
        w(4)=0.620662603527d0
        w(5)=0.593027449764d0
        w(6)=0.576193350284d0
        w(7)=0.567021153447d0
        w(8)=0.564100308726d0
        w(9)=0.567021153447d0
        w(10)=0.576193350284d0
        w(11)=0.593027449764d0
        w(12)=0.620662603527d0
        w(13)=0.666166005109d0
        w(14)=0.748607366017d0
        w(15)=0.948368970828d0
    endif

    if(npg.eq.20) then
        x(1)=-5.38748089001d0
        x(2)=-4.60368244955d0
        x(3)=-3.94476404012d0
        x(4)=-3.34785456738d0
        x(5)=-2.78880605843d0
        x(6)=-2.25497400209d0
        x(7)=-1.73853771212d0
        x(8)=-1.2340762154d0
        x(9)=-0.737473728545d0
        x(10)=-0.245340708301d0
        x(11)=0.245340708301d0
        x(12)=0.737473728545d0
        x(13)=1.2340762154d0
        x(14)=1.73853771212d0
        x(15)=2.25497400209d0
        x(16)=2.78880605843d0
        x(17)=3.34785456738d0
        x(18)=3.94476404012d0
        x(19)=4.60368244955d0
        x(20)=5.38748089001d0

        w(1)=0.898591961453d0
        w(2)=0.704332961176d0
        w(3)=0.62227869619d0
        w(4)=0.575262442852d0
        w(5)=0.544851742366d0
        w(6)=0.524080350949d0
        w(7)=0.509679027117d0
        w(8)=0.499920871336d0
        w(9)=0.493843385272d0
        w(10)=0.490921500667d0
        w(11)=0.490921500667d0
        w(12)=0.493843385272d0
        w(13)=0.499920871336d0
        w(14)=0.509679027117d0
        w(15)=0.524080350949d0
        w(16)=0.544851742366d0
        w(17)=0.575262442852d0
        w(18)=0.62227869619d0
        w(19)=0.704332961176d0
        w(20)=0.898591961453d0
    endif


    if(npg.eq.25) then

        x(1)=-6.16427243405d0
        x(2)=-5.41363635528d0
        x(3)=-4.78532036736d0
        x(4)=-4.21860944438d0
        x(5)=-3.69028287701d0
        x(6)=-3.18829492442d0
        x(7)=-2.70532023717d0
        x(8)=-2.23642013027d0
        x(9)=-1.77800112434d0
        x(10)=-1.32728070207d0
        x(11)=-0.881982756214d0
        x(12)=-0.440147298645d0
        x(13)=0.d0
        x(14)=0.440147298645d0
        x(15)=0.881982756214d0
        x(16)=1.32728070207d0
        x(17)=1.77800112434d0
        x(18)=2.23642013027d0
        x(19)=2.70532023717d0
        x(20)=3.18829492442d0
        x(21)=3.69028287701d0
        x(22)=4.21860944438d0
        x(23)=4.78532036736d0
        x(24)=5.41363635528d0
        x(25)=6.16427243405d0
    
        w(1)=0.862401988724d0
        w(2)=0.673022290239d0
        w(3)=0.592081693052d0
        w(4)=0.54491777224d0
        w(5)=0.513655789745d0
        w(6)=0.49150688189d0
        w(7)=0.475249738004d0
        w(8)=0.463141046576d0
        w(9)=0.454155885528d0
        w(10)=0.447661256587d0
        w(11)=0.443259189252d0
        w(12)=0.440705828912d0
        w(13)=0.439868722169d0
        w(14)=0.440705828912d0
        w(15)=0.443259189252d0
        w(16)=0.447661256587d0
        w(17)=0.454155885528d0
        w(18)=0.463141046576d0
        w(19)=0.475249738004d0
        w(20)=0.49150688189d0
        w(21)=0.513655789745d0
        w(22)=0.54491777224d0
        w(23)=0.592081693052d0
        w(24)=0.673022290239d0
        w(25)=0.862401988724d0
    endif

    if(npg.eq.30) then

        x(1)=-6.863345294d0
        x(2)=-6.13827922d0
        x(3)=-5.533147152d0
        x(4)=-4.988918969d0
        x(5)=-4.483055357d0
        x(6)=-4.003908604d0
        x(7)=-3.544443873d0
        x(8)=-3.09997053d0
        x(9)=-2.667132125d0
        x(10)=-2.243391468d0
        x(11)=-1.826741144d0
        x(12)=-1.4155278d0
        x(13)=-1.008338271d0
        x(14)=-0.603921059d0
        x(15)=-0.201128577d0
        x(16)=0.201128577d0
        x(17)=0.603921059d0
        x(18)=1.008338271d0
        x(19)=1.4155278d0
        x(20)=1.826741144d0
        x(21)=2.243391468d0
        x(22)=2.667132125d0
        x(23)=3.09997053d0
        x(24)=3.544443873d0
        x(25)=4.003908604d0
        x(26)=4.483055357d0
        x(27)=4.988918969d0
        x(28)=5.533147152d0
        x(29)=6.13827922d0
        x(30)=6.863345294d0

        w(1)=0.834247471d0
        w(2)=0.649097981d0
        w(3)=0.569402692d0
        w(4)=0.522525689d0
        w(5)=0.491057996d0
        w(6)=0.468374813d0
        w(7)=0.451321036d0
        w(8)=0.438177023d0
        w(9)=0.427918063d0
        w(10)=0.419895004d0
        w(11)=0.413679364d0
        w(12)=0.408981575d0
        w(13)=0.405605123d0
        w(14)=0.403419817d0
        w(15)=0.402346067d0
        w(16)=0.402346067d0
        w(17)=0.403419817d0
        w(18)=0.405605123d0
        w(19)=0.408981575d0
        w(20)=0.413679364d0
        w(21)=0.419895004d0
        w(22)=0.427918063d0
        w(23)=0.438177023d0
        w(24)=0.451321036d0
        w(25)=0.468374813d0
        w(26)=0.491057996d0
        w(27)=0.522525689d0
        w(28)=0.569402692d0
        w(29)=0.649097981d0
        w(30)=0.834247471d0

    endif
    auxfunca=0.d0
    ss=0.d0
    ss1=0.d0
    ss2=0.d0
    if(npoint.eq.2) then
      !!print*,'je suis la'
      ii=0
      do ii=1,npg
          ss1=0.d0
          do jj=1,npg
              tab(1)=x(ii)
              tab(2)=x(jj)
              auxfunca=func30(tab,trial_i)
              ss1 = ss1+w(jj)*(auxfunca)
          end do
          ss = ss+w(ii)*ss1
      end do
    else ! cas de 3 points
    !!print*,'non je suis plutot la', npoint
      do kk=1,npg
        ss2=0.d0
        do ii=1,npg
          ss1=0.d0
          do jj=1,npg
              tab(1)=x(kk)
              tab(2)=x(ii)
              tab(3)=x(jj)
              auxfunca=func30(tab,trial_i)
              ss1 = ss1+w(jj)*(auxfunca)
          end do
          ss2 = ss2+w(ii)*ss1
        end do
        ss = ss+w(kk)*ss2
      end do
    end if

    end subroutine gausshermite_2_3
end module GaussHermi_mult

