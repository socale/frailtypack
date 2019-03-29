subroutine Generation_surrogate_copula(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
            ng,ver,truealpha,propC,cens_A,gamma1,gamma2,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,&
            betat,n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma,alpha,frailt_base,thetacopule)
    ! lognormal: dit si la distribution des effets aleatoires est lognormal pour le modele complet (1) ou lognormal pour le joint classique de 2007 (2) ou gamma pour le joint classique de 2007(0)
    ! use Autres_fonctions
    ! theta2: variance des frailties gaussiens associe a S
    ! gamma1,gamma2: parametres de la gamma
    ! cens_A:censure administrative
    ! propC:proportion de personnes censurees
    ! n_essai:  nombre d'essai a generer
    ! rsqrt: niveau de correlation souhaite entre les fragilites sepecifiques aux traitement au niveau essai S et T
    ! sigma_s: variance des frailties au niveau essai associee au surrogate
    ! sigma_t: variance des frailties au niveau essai associee au true
    ! p: proportion des personnes traitees par essai
    ! prop_i: proportion des sujets par essai
    ! gamma: variance de l'effet aleatoire u_i associe au risque de base chez S
    ! alpha: parametre de puissance (zeta) associe a u_i pour les deces
    ! frailt_base: dit si l'on prend en compte l'heterogeneite sur le risque de base aussi bien dans la generation des donnes que dans l'estimation(1) ou non (0)
	! thetacopule : parametre de la copule de clayton
     use var_surrogate, only: random_generator

      integer, intent(in)::n_essai,frailt_base,affiche_stat,n_obs,n_col,lognormal,ng,ver
      double precision, intent(in)::truealpha,propC,cens_A,gamma1,gamma2,theta2,gamma,alpha,&
                                    lambda_S,nu_S,lambda_T,nu_T,rsqrt,sigma_s,sigma_t,&
									thetacopule
	  double precision, dimension(ver), intent(in)::betas,betat ! vecteur des coefficients du model. contiont le meme nombre d'element, et donc doit etre bien rempli
      double precision,dimension(n_essai),intent(in)::prop_i,p      
      double precision,intent(out)::vrai_theta
      double precision,dimension(n_obs,n_col),intent(out)::don_simulS1,don_simul
      
	  integer, parameter::npmax=70,NOBSMAX=15000,nvarmax=45,ngmax=5000,nboumax=1000,NSIMAX=5000,&
                          ndatemax=30000
      integer  j,k,n,ii,iii,iii2,i,nbou
      integer filtre(nvarmax), filtre2(nvarmax), ind(nboumax), values(8)
      double precision maxtemps
      character*18 nomvarl,nomvar(nvarmax),nomvar2(nvarmax)
      character*20 dateamj, zone, heure1

!************declarations pour donnees generees **********
      integer :: nb_recur,nb_dc,nb_cens,delta,deltadc,jj,ig,nrecurr,nobs,max_recu
      real , dimension(:), allocatable:: v1
      real :: piece,demi
      double precision :: ui,temps1,temps1_S,gapx,gapdc,moy_idnum,x,xdc,cens,cbeta1,cbeta2,cbeta3,&
      auxbeta1,auxbeta2,uniran,moyui, cbeta2, cbeta4
      double precision, dimension(2):: bg1,bw1,bw2
      integer, dimension(ngmax):: idnum
      double precision, dimension(ngmax):: vecui

      integer :: nmax,nva,nva1,nva2, ibou
      common /nmax/nmax
      double precision date(ndatemax), zi(-2:npmax)
      common /dace1/date,zi
      double precision t0(NOBSMAX),t1(NOBSMAX),t1_S(NOBSMAX), ve(n_obs,ver),ve2(n_obs,ver)
      integer c(NOBSMAX), cdc(NOBSMAX), g(NOBSMAX), nig(ngmax)

      ! Ajout SCL
      double precision::x22,sigma_st,u_i,tempon
      double precision,dimension(:),allocatable::tempsD,mu
      integer::Aretenir
      integer,parameter ::trt1=1,v_s1=2,v_t1=3,trialref1=4,w_ij1=5,timeS1=6,timeT1=7,&
                      timeC1=8,statusS1=9,statusT1=10,initTime1=11,Patienref1=12,u_i1=13 ! definissent les indices du tableau de donnee simulees
      double precision,dimension(n_essai)::n_i
      double precision,dimension(:,:),allocatable::sigma,x_ 
    
      
!CCCCCCCCCCCCCCCCChosur9.f CCCCCCCCCCCCCCCCCCCCCCCC
      allocate(v1(ver))
	  don_simulS1 = 0.d0
	  don_simul = 0.d0
      
      call date_and_time(dateamj,heure1,zone,values)                  
      nrecurr=1
      nbou=1
      
      allocate(tempsD(ng))
      !nst=2

         nva1 = 0 ! nb de var expli pour donnees recurrentes
         nva2 = 0 ! nb de var expli pour deces

         if(ver.gt.0)then
            do 44 j=1,ver
               nomvarl="trt"
               filtre(j)=1
               filtre2(j)=1
               nva1 = nva1 + filtre(j) ! adjustment for recurrent events
               nva2 = nva2 + filtre2(j) ! adjustment for survival
               if(filtre(j).eq.1)then
                  nomvar(nva1) = nomvarl
               endif 
               if(filtre2(j).eq.1)then
                  nomvar2(nva2) = nomvarl
               endif    
 44         continue
         endif

 !        nva = nva1+nva2
    
    Aretenir=1
        
!c*******************************************
!c***** DEBUT SIMULATIONS *******************
!c*******************************************
    maxtemps = 0.d0
    don_simulS1=0.d0
    don_simul=0.d0
    
!==============initialisation des parametres======================
    ! call dblepr("prop_i", -1, prop_i(:), size(prop_i))
    n_i=NINT(n_obs*prop_i) ! nombre de sujet par essai

    if(sum(n_i)<n_obs) then
        n_i(minloc(n_i,mask=n_i==minval(n_i)))=n_i(minloc(n_i,mask=n_i==minval(n_i)))+(n_obs-sum(n_i)) ! on ajoute au premier essai de plus petite taille le nombre d'individu non encore affecte (1 generalement) a cause des problemes d'arrondi
    endif
    if(sum(n_i)>n_obs) then 
        n_i(maxloc(n_i,mask=n_i==maxval(n_i)))=n_i(maxloc(n_i,mask=n_i==maxval(n_i)))-(sum(n_i)-n_obs) ! on soustrait au premier essai de plus grande taille le nombre d'individu affecte en trop (1 generalement) a cause des problemes d'arrondi
    endif
    
    ! simulation des effets aleatoires specifiques aux essais
    
    ! matrice des covariances des frailties au niveau essai, en supposant une correlation de 0.95 (erreurs liees aux effets du traitement sur S et T niveau essai) 
    sigma_st=rsqrt*sqrt(sigma_s)*sqrt(sigma_t)
    
    if(frailt_base==1) then
        allocate(sigma(3,3),mu(3),x_(n_essai,3))
        sigma(1,1)=sigma_s
        sigma(1,2)=sigma_st
        sigma(2,1)=sigma_st
        sigma(2,2)=sigma_t
        !pour u_i
        sigma(3,1:2)=0.d0
        sigma(1:2,3)=0.d0
        sigma(3,3)=gamma
        !!print*,sigma
        !stop
    else
        allocate(sigma(2,2),mu(2),x_(n_essai,2))
        sigma(1,1)=sigma_s
        sigma(1,2)=sigma_st
        sigma(2,1)=sigma_st
        sigma(2,2)=sigma_t
    endif
    mu=0.d0
    
    !generation de (vs_i, vt_i) suivant une multinormale
    call rmvnorm(mu,sigma,n_essai,0,x_)    
    
    k=1
    do i=1,n_essai
        ! effet aleatoires specifiques aux essais
        don_simul(k:((k+NINT(n_i(i)))-1),v_s1)=x_(i,1)
        don_simul(k:((k+NINT(n_i(i)))-1),v_t1)=x_(i,2)
        !don_simul(k:((k+NINT(n_i(i)))-1),u_i1)=x_(i,3)
        ! simulation des u_i associee aux risques de base suivant une normale
        if(frailt_base==1) then
            ! call bgos(sqrt(gamma),0,u_i,x22,0.d0)
            u_i=x_(i,3)
        else
            u_i=0.d0
        endif
        don_simul(k:((k+NINT(n_i(i)))-1),u_i1)=u_i
        ! variable trialref
        don_simul(k:((k+NINT(n_i(i)))-1),trialref1)=i
        !!print*,don_simul(k:((k+NINT(n_i(i)))-1),trialref1)
        k=k+n_i(i)
    enddo
    
    don_simulS1=don_simul
    do 1000 ibou=1,nbou 
        ind(ibou)=0
        nig = 0
        g=0
        maxtemps = 0.d0
        auxbeta1=0.d0
        auxbeta2=0.d0

!c********************************************************
!c******** DEBUT generation des donnees *******************
!c********************************************************
            bg1(1)=gamma1
            bg1(2)=gamma1
           if(ibou.eq.1)then
                if (lognormal==0)then ! Gamma
                    vrai_theta=bg1(1)/(bg1(2)*bg1(2))
                else !lognormale
                    vrai_theta=theta2
                endif
           endif

!!c-- parametres de la WEibull for recurrent events 
         bw1(1)=lambda_S 
         bw1(2)=nu_S
!!c-- parametres de la WEibull for death
         bw2(1)=lambda_T
         bw2(2)=nu_T

         demi = 0.5             ! pour var expli
         cbeta1=betas          
         cbeta3=betat        
         nb_recur = 0           ! nb de temps 
         nb_dc = 0              ! nb de dc
         nb_cens = 0            ! nb de cens
         nobs=0
         max_recu= 1
         moyui = 0.d0
        !close(10)
    do 30 ig=1,ng ! sur les groupes
            
        if (lognormal==0)then ! Gamma    
            call gamgui(bg1(1),ui) !genere gamma pour zi
            ui = ui/bg1(2)
!c     verification des ui
            vecui(ig) = ui
            moyui = moyui + ui
        else !lognormale
            x22=0.d0
            !bg1(1)= variance des wij
            call bgos(sqrt(theta2),0,ui,x22,0.d0) !on simule les w_ij suivant une normale
            !verification des ui
            vecui(ig) = ui
            moyui = moyui + ui
        endif
        
!!c---  variables explicatives par sujet
         do 111 j=1,ver
                if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
                    tempon= uniran()
                else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
                    CALL RANDOM_NUMBER(tempon)
                endif
               piece=real(tempon)
               if (piece.le.demi) then
                  v1(j) = 0.
               else
                  v1(j) = 1.
               endif
         111        continue		   
         x=0.d0
         xdc=0.d0
         cens=0.d0

         do 10 k=1,nrecurr ! observations max / sujet
			 if(k.gt.max_recu)then
				max_recu=k 
			 endif
             nobs=nobs+1   ! indice l ensemble des observations
             idnum(ig) = k
!c-----------Surrogate --------------------------------------------
!c---  genere temps de progression a partir 2 var explic :
			if (lognormal==0)then ! Gamma

			else ! lognormale
				if (lognormal==1)then !joint surrogate
					if(frailt_base==1) then ! on tient compte des u_i
						auxbeta1=don_simul(ig,u_i1)+don_simul(ig,v_s1)*dble(v1(1))+betas(1)*dble(v1(1))
						auxbeta2=alpha*don_simul(ig,u_i1)+don_simul(ig,v_t1)*dble(v1(1))+betat(1)*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
					    do ii = 2,ver ! on considere a partir de la 2 ieme variable car le traitement est prise en compte deja
						  if(filtre(ii).eq.1)then
						    auxbeta1 = auxbeta1 + betas(ii)*dble(v1(ii))
						  endif
						  if(filtre2(ii).eq.1)then
							auxbeta2 = auxbeta2 + betat(ii)*dble(v1(ii))
						  endif
                        enddo
					else ! on ne tient pas compte des u_i dans la generation des temps de survie
						auxbeta1=don_simul(ig,v_s1)*dble(v1(1))+betas(1)*dble(v1(1))
						auxbeta2=don_simul(ig,v_t1)*dble(v1(1))+betat(1)*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
					    do ii = 2,ver ! on considere a partir de la 2 ieme variable car le traitement est prise en compte deja
						  if(filtre(ii).eq.1)then
						    auxbeta1 = auxbeta1 + betas(ii)*dble(v1(ii))
						  endif
						  if(filtre2(ii).eq.1)then
							auxbeta2 = auxbeta2 + betat(ii)*dble(v1(ii))
						  endif
                        enddo
					endif
				endif
			endif

			call weiguicopule(bw1(1),bw2(1),bw1(2),bw2(2),auxbeta1,auxbeta2,thetacopule,gapx,gapdc)
			x=gapx ! temps de progression
			xdc=gapdc ! temps de deces
			tempsD(ig)=xdc

		! scl============censure====================
				cens=cens_A

				if(xdc.le.cens)then ! patient decede
				  deltadc=1.d0
				  temps1 = xdc
				  nb_dc =nb_dc + 1
			    else    !patient censuree administrativement
				  deltadc=0.d0
				  temps1 = cens
				  nb_cens =nb_cens + 1
			    endif

				!on construit les temps de progression
				if(x < temps1)then ! evenement avant la censure
					delta=1.d0
					temps1_S = x
					nb_recur =nb_recur + 1
					nig(ig) = nig(ig)+1 !nb events recurrents
				else
					if((x.eq.cens).and.(deltadc==0.d0)) then !evenement a la date de censure et patient vivant
						delta=1.d0
						temps1_S = x
						nb_recur =nb_recur + 1
						nig(ig) = nig(ig)+1 !nb events recurrents
					else ! progression le meme jour que le deces ou sans progression
						delta=0.d0             ! on suppose pas d'evenement si le meme jour que le deces
						temps1_S=temps1! et on censure a la date de deces(ou censure)
					endif
				endif
		!c****** for gap time :         
					   t0(nobs) = 0.d0
		!c fin gap
					   t1(nobs) = temps1
					   t1_S(nobs) = temps1_S
					   c(nobs) = delta
					   cdc(nobs) = deltadc
					   g(nobs)= ig
					   iii = 0
					   iii2 = 0
					   do 110 ii = 1,ver
						  if(filtre(ii).eq.1)then
							 iii = iii + 1
							 ve(nobs,iii) = dble(v1(ii))
						  endif
						  if(filtre2(ii).eq.1)then
							 iii2 = iii2 + 1
							 ve2(nobs,iii2) = dble(v1(ii))
						  endif
		 110           continue
	           !c*** pour le tester sur un autre programme: on complete les nouveaux parametres simules dans le jeux de donnees
				don_simulS1(ig,initTime1)=t0(nobs)
				don_simulS1(ig,timeS1)=t1_S(nobs)
				don_simulS1(ig,statusS1)=c(nobs)
				don_simulS1(ig,Patienref1)=g(nobs)
                don_simulS1(ig,trt1)=ve2(nobs,1)
				don_simulS1(ig,w_ij1)=ui		
				   
                don_simul(ig,initTime1)=t0(nobs)
				don_simul(ig,timeT1)=t1(nobs)
				don_simul(ig,statusT1)=cdc(nobs)
				don_simul(ig,Patienref1)=g(nobs)
				don_simul(ig,trt1)=ve2(nobs,1)
				don_simul(ig,w_ij1)=ui       
                
               do ii = 2,ver
					if(filtre(ii).eq.1)then
					  don_simulS1(ig,6+ii-1)=ve2(nobs,ii)
					endif
					if(filtre2(ii).eq.1)then
					  don_simul(ig,6+ii-1)=ve2(nobs,ii)
					endif
		        enddo				

				!c****************************************************
				!c******** FIN generation des donnees****************
				!c****************************************************
          
				if (maxtemps.lt.t1(nobs))then
				   maxtemps = t1(nobs)
				endif  

				if (delta.eq.0.d0)then ! deces ou censure
				   goto 30          ! on change de sujet
				endif				   
		   10      continue                   !observations par sujet
		   30   continue                  !sujets=groupes		   
			  if(ibou.eq.Aretenir)then
				  moy_idnum=0
				  do 444 jj=1,ng
					 moy_idnum =moy_idnum + idnum(jj)
				  444     continue
					  moy_idnum =moy_idnum / ng 
			  endif		   		   
	!c=========     on retourne a l'iteration suivante
	 1000 continue
     ! scl============censure conseillee pour la proportion souhaitee de censure==
     call percentile_scl(tempsD,ng,1.d0-propC,cens)
    deallocate(tempsD,sigma,mu,x_,v1)
 
endsubroutine Generation_surrogate_copula




















	subroutine weiguicopule(a,at,b,bt,betau,betaut,theta,Sij,Tij)
	!theta : parametre de copula
    ! fonction de densité de la loi de weibull = f(x)=b**a . a . x**(a-1) . exp(-(bx)**a) (voir cours de Piere Jolie page 41)
    ! generation du temps de deces conditionnellement au surrogate: T_ij|S_ij
        use var_surrogate, only:param_weibull
        use var_surrogate, only: random_generator
        double precision ::a,b,at,bt,Sij,Tij,u,ut,v,betau,betaut,vt,utij,vtij,theta
        double precision ::uniran
        real ::ran2
        if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
            u = uniran()
			ut = uniran()
			utij = ((ut**(-theta/(1.d0 + theta)) - 1.d0) * u**(- theta) + 1.d0)**(-1.d0 / theta)
        else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)    
            CALL RANDOM_NUMBER(u)
        endif
        v = (1.d0-u)
		vt = (1.d0-ut)
		vtij = ((vt**(-theta/(1.d0 + theta)) - 1.d0) * v**(- theta) + 1.d0)**(-1.d0 / theta)
		
        if(param_weibull==0)then !parametrisation de weibull par defaut dans le programme de Virginie: fonction de densite differente de celle donnee ci-dessus
            Sij = (1.d0/b)*((-dexp(-betau)*dlog(v))**(1.d0/a))
			Tij = (1.d0/bt)*((-dexp(-betaut)*dlog(vtij))**(1.d0/at))
        else 
            ! parametrisation de la weibull donnee par la fonction de densite ci-dessus
            ! F^-1(t)=1/b*(-log(1-t))**(1/a)
            Sij=    (1.d0/b)*((-dlog(1+ dlog(u)*dexp(-betau)))**(1.d0/a))
			Tij=    (1.d0/bt)*((-dlog(1+ dlog(utij)*dexp(-betaut)))**(1.d0/at))
        endif

    return
    endsubroutine weiguicopule
