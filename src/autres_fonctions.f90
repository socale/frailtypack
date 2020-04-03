module InverseMatrix
!fonction qui retourne l'inverse et le determinant d'une matrice carree

    implicit none
    
    contains
    
    subroutine matinv(A,B,det)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    implicit none
    
    double precision, intent(in):: A(2,2)   !! Matrix
    double precision,intent(out):: B(2,2)   !! Inverse matrix
    double precision:: detinv
    double precision, intent(out):: det

    ! Calculate the inverse determinant of the matrix
    det=A(1,1)*A(2,2) - A(1,2)*A(2,1)
    detinv = 1.d0/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
    ! !print*,"detinv=",detinv,"det=",det
    ! !print*,"A=",A
    ! !print*,"B=",B
    ! stop
  end subroutine matinv
  end module InverseMatrix
  
   Module Autres_fonctions
  
    implicit none
    
    contains
    
    ! fonction pour le calcul du taux de kendall
    double precision function tau_kendall(sigma_w,sigma_u,sigma_v,z_11,z_21,method_int,N_MC,alpha,zeta,model_complet)
        !sigma_w : matrice des covariances frailties niveau individuel
        !sigma_u : matrice des covariances frailties associes aux risque de base
        !sigma_v: matrice des covariances frailties niveau essai en interaction avec le traitement
        !z_11 : indicatrice de traitement individu 1
        !z_21 : indicatice de tritement individu 2
        !method_int : methode d'integration: 0= montecarle, 1= quadrature quaussienne classique, 2= approximation de Laplace, 
        !             4=integration par monte-carlo, 1 seul taux de kendall, 5=integration par monte-carlo, 1 seul taux de kendall et 
        !             pas de stratification sur les risques de base
        !N_MC: nombre de boucle MC ou nombre de points de quadrature si method_int=1
        !alpha: fragilite associe a w_ij
        !zeta: fragilite associe a u_i
        !model_complet: dit si on utilise le model(1) complet ou pas (0)
        
        use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
        
        implicit none
        
        double precision,dimension(:,:),intent(in)::sigma_w,sigma_u,sigma_v
        double precision,intent(in)::alpha,zeta
        integer, intent(in)::z_11,z_21, method_int,N_MC,model_complet
        double precision,dimension(:),allocatable::mu,xx1,ww1
        double precision,dimension(:,:),allocatable::u,up,v,vp,w11,wp,sigma,x_
        integer:: i,nnodes
        double precision::integral,rho_wst,rho_ust,rho_vst,somm
        double precision,parameter::pi=3.141592653589793d0
        
        ! generation des variables aleatoires suivant des multinormales
    
        integral=0.d0
        !!print*,"sigma_w=",sigma_w
        !!print*,"sigma_u=",sigma_u
        !!print*,"sigma_v=",sigma_v

        
        select case(method_int)
        case(0) ! integration par monte-carlo
            allocate(mu(2),u(1,2),up(1,2),v(1,2),vp(1,2),w11(1,2),wp(1,2))
            mu=0.d0
            do i=1,N_MC
                !generation des donnees
                call rmvnorm(mu,sigma_w,1,0,w11)    ! (ws_ij,wt_ij)
                call rmvnorm(mu,sigma_w,1,0,wp)    ! (ws_i'j',wt_i'j')
                call rmvnorm(mu,sigma_u,1,0,u)    ! (us_i,ut_i)
                call rmvnorm(mu,sigma_u,1,0,up)    ! (us_i',ut_i')
                call rmvnorm(mu,sigma_v,1,0,v)    ! (vs_i,vt_i)
                call rmvnorm(mu,sigma_v,1,0,vp)    ! (vs_i',vt_i')
                
                !evaluation de l'integrant
                integral=integral+(dexp(w11(1,1)+u(1,1)+v(1,1)*z_11+w11(1,2)+u(1,2)+v(1,2)*z_21)+&
                     dexp(wp(1,1)+up(1,1)+vp(1,1)*z_11+wp(1,2)+up(1,2)+vp(1,2)*z_21))/&
                    ((dexp(wp(1,1)+up(1,1)+vp(1,1)*z_11)+dexp(w11(1,1)+u(1,1)+v(1,1)*z_11))*&
                    (dexp(wp(1,2)+up(1,2)+vp(1,2)*z_21)+dexp(w11(1,2)+u(1,2)+v(1,2)*z_21)))
            enddo
            deallocate(mu,u,up,v,vp,w11,wp)
            tau_kendall=2.d0*integral/N_MC -1.d0
            
        case(2) ! integration par quadrature approximation de laplace
        
        case(3) ! integration par monte-carlo: je suppose une multinormale avec une seule matrice
            allocate(mu(6),u(1,2),up(1,2),v(1,2),vp(1,2),w11(1,2),wp(1,2),sigma(12,12),x_(1,12))
            mu=0.d0
            sigma=0.d0
            sigma(1:2,1)=(/sigma_w(1:2,1)/)
            sigma(1:2,2)=(/sigma_w(1:2,2)/)
            sigma(3:4,3)=(/sigma_w(1:2,1)/)
            sigma(3:4,4)=(/sigma_w(1:2,2)/)
            sigma(5:6,5)=(/sigma_u(1:2,1)/)
            sigma(5:6,6)=(/sigma_u(1:2,2)/)
            sigma(7:8,7)=(/sigma_u(1:2,1)/)
            sigma(7:8,8)=(/sigma_u(1:2,2)/)
            sigma(9:10,9)=(/sigma_v(1:2,1)/)
            sigma(9:10,10)=(/sigma_v(1:2,2)/)
            sigma(11:12,11)=(/sigma_v(1:2,1)/)
            sigma(11:12,12)=(/sigma_v(1:2,2)/)
            
            do i=1,N_MC
                !generation des donnees
                call rmvnorm(mu,sigma,1,0,x_)    ! (ws_ij,wt_ij)
                w11(1,1:2)=(/x_(1,1:2)/)
                wp(1,1:2)=(/x_(1,3:4)/)
                u(1,1:2)=(/x_(1,5:6)/)
                u(1,1:2)=(/x_(1,7:8)/)
                v(1,1:2)=(/x_(1,9:10)/)
                v(1,1:2)=(/x_(1,11:12)/)
                
                !evaluation de l'integrant
                integral=integral+(dexp(w11(1,1)+u(1,1)+v(1,1)*z_11+w11(1,2)+u(1,2)+v(1,2)*z_21)+&
                     dexp(wp(1,1)+up(1,1)+vp(1,1)*z_11+wp(1,2)+up(1,2)+vp(1,2)*z_21))/&
                    ((dexp(wp(1,1)+up(1,1)+vp(1,1)*z_11)+dexp(w11(1,1)+u(1,1)+v(1,1)*z_11))*&
                    (dexp(wp(1,2)+up(1,2)+vp(1,2)*z_21)+dexp(w11(1,2)+u(1,2)+v(1,2)*z_21)))
            enddo
            deallocate(mu,u,up,v,vp,w11,wp,sigma,x_)
            tau_kendall=2.d0*integral/N_MC -1.d0
            
        case(4) ! integration par monte-carlo, 1 seul taux de kendall
            if(model_complet==0) then ! alors modele reduit avec effect aleatoires partages
                allocate(mu(1),u(1,1),up(1,1),w11(1,1),wp(1,1))
                mu=0.d0
                integral=0.d0
                somm=0.d0
                do i=1,N_MC
                    !generation des donnees
                    call rmvnorm(mu,sigma_w,1,0,w11)    ! (ws_ij,wt_ij)
                    call rmvnorm(mu,sigma_w,1,0,wp)    ! (ws_i'j',wt_i'j')
                    call rmvnorm(mu,sigma_u,1,0,u)    ! (us_i,ut_i)
                    call rmvnorm(mu,sigma_u,1,0,up)    ! (us_i',ut_i')
                    ! call rmvnorm(mu,sigma_v,1,0,v)    ! (vs_i,vt_i)
                    ! call rmvnorm(mu,sigma_v,1,0,vp)    ! (vs_i',vt_i')
                    
                    !evaluation de l'integrant
                    somm=((dexp(w11(1,1)+u(1,1)+w11(1,1)*zeta+u(1,1)*alpha)+&
                         dexp(wp(1,1)+up(1,1)+wp(1,1)*zeta+up(1,1)*alpha))/&
                        ((dexp(wp(1,1)+up(1,1))+dexp(w11(1,1)+u(1,1)))*&
                        (dexp(wp(1,1)*zeta+up(1,1)*alpha)+dexp(w11(1,1)*zeta+u(1,1)*alpha))))
                    integral=integral+somm
                    ! !print*,i,somm
                enddo
                deallocate(mu,u,up,w11,wp)
                tau_kendall=2.d0*integral/N_MC -1.d0
                ! !print*,"tau_kendall",tau_kendall
                ! stop
            else ! modele complet
                allocate(mu(2),u(1,2),up(1,2),w11(1,2),wp(1,2))
                mu=0.d0
                do i=1,N_MC
                    !generation des donnees
                    call rmvnorm(mu,sigma_w,1,0,w11)    ! (ws_ij,wt_ij)
                    call rmvnorm(mu,sigma_w,1,0,wp)    ! (ws_i'j',wt_i'j')
                    call rmvnorm(mu,sigma_u,1,0,u)    ! (us_i,ut_i)
                    call rmvnorm(mu,sigma_u,1,0,up)    ! (us_i',ut_i')
                    ! call rmvnorm(mu,sigma_v,1,0,v)    ! (vs_i,vt_i)
                    ! call rmvnorm(mu,sigma_v,1,0,vp)    ! (vs_i',vt_i')
                    
                    !evaluation de l'integrant
                    integral=integral+((dexp(w11(1,1)+u(1,1)+w11(1,2)+u(1,2))+&
                         dexp(wp(1,1)+up(1,1)+wp(1,2)+up(1,2)))/&
                        ((dexp(wp(1,1)+up(1,1))+dexp(w11(1,1)+u(1,1)))*&
                        (dexp(wp(1,2)+up(1,2))+dexp(w11(1,2)+u(1,2)))))
                enddo
                deallocate(mu,u,up,w11,wp)
                tau_kendall=2.d0*integral/N_MC -1.d0
            endif
            case(5) ! integration par monte-carlo, 1 seul taux de kendall oon fixe u_i a 0
            if(model_complet==0) then ! alors modele reduit avec effect aleatoires partages
                allocate(mu(1),u(1,1),up(1,1),w11(1,1),wp(1,1))
                mu=0.d0
                do i=1,N_MC
                    !generation des donnees
                    call rmvnorm(mu,sigma_w,1,0,w11)    ! (ws_ij,wt_ij)
                    call rmvnorm(mu,sigma_w,1,0,wp)    ! (ws_i'j',wt_i'j')
                    ! call rmvnorm(mu,sigma_u,1,0,u)    ! (us_i,ut_i)
                    ! call rmvnorm(mu,sigma_u,1,0,up)    ! (us_i',ut_i')
                    ! call rmvnorm(mu,sigma_v,1,0,v)    ! (vs_i,vt_i)
                    ! call rmvnorm(mu,sigma_v,1,0,vp)    ! (vs_i',vt_i')
                    
                    !evaluation de l'integrant
                    integral=integral+((dexp(w11(1,1)+w11(1,1)*zeta)+&
                         dexp(wp(1,1)+wp(1,1)*zeta))/&
                        ((dexp(wp(1,1))+dexp(w11(1,1)))*&
                        (dexp(wp(1,1)*zeta)+dexp(w11(1,1)*zeta))))
                enddo
                deallocate(mu,u,up,w11,wp)
                tau_kendall=2.d0*integral/N_MC -1.d0
            else ! modele complet
                allocate(mu(2),u(1,2),up(1,2),w11(1,2),wp(1,2))
                mu=0.d0
                do i=1,N_MC
                    !generation des donnees
                    call rmvnorm(mu,sigma_w,1,0,w11)    ! (ws_ij,wt_ij)
                    call rmvnorm(mu,sigma_w,1,0,wp)    ! (ws_i'j',wt_i'j')
                    ! call rmvnorm(mu,sigma_u,1,0,u)    ! (us_i,ut_i)
                    ! call rmvnorm(mu,sigma_u,1,0,up)    ! (us_i',ut_i')
                    ! call rmvnorm(mu,sigma_v,1,0,v)    ! (vs_i,vt_i)
                    ! call rmvnorm(mu,sigma_v,1,0,vp)    ! (vs_i',vt_i')
                    
                    !evaluation de l'integrant
                    integral=integral+((dexp(w11(1,1)+w11(1,2))+&
                         dexp(wp(1,1)+wp(1,2)))/&
                        ((dexp(wp(1,1))+dexp(w11(1,1)))*&
                        (dexp(wp(1,2))+dexp(w11(1,2)))))
                enddo
                deallocate(mu,u,up,w11,wp)
                tau_kendall=2.d0*integral/N_MC -1.d0
            endif
                case(1) ! integration par quadrature gaussienne
            nnodes=N_MC
            allocate(xx1(nnodes))
            allocate(ww1(nnodes))
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
            
            rho_wst = sigma_w(1,2)/dsqrt(sigma_w(1,1)*sigma_w(2,2))
            rho_ust = sigma_u(1,2)/dsqrt(sigma_u(1,1)*sigma_u(2,2))
            rho_vst = sigma_v(1,2)/dsqrt(sigma_v(1,1)*sigma_v(2,2))
            
            integral=0.d0
            tau_kendall=integral-1.d0
            deallocate(xx1,ww1)
        endselect
        
        
        return
        
    end function tau_kendall
    
    
    !fin fonction tau de kendall
    
    
    ! --------------------------------------------------------------------
    !      REAL FUNCTION  Median() :
    !    This function receives an array X of N entries, copies its value
    !      to a local array Temp(), sorts Temp() and computes the median.
    !    The returned value is of REAL type.
    ! --------------------------------------------------------------------

   REAL FUNCTION  Median(X, N)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(IN) :: X
      INTEGER, INTENT(IN)                :: N
      INTEGER, DIMENSION(1:N)            :: Temp
      INTEGER                            :: i

      DO i = 1, N                       ! make a copy
         Temp(i) = X(i)
      END DO
      CALL  Sort(Temp, N)               ! sort the copy
      IF (MOD(N,2) == 0) THEN           ! compute the median
         Median = (Temp(N/2) + Temp(N/2+1)) / 2.0
      ELSE
         Median = Temp(N/2+1)
      END IF
   END FUNCTION  Median
   
   ! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      INTEGER                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)        ! assume the first is the min
      Location = Start            ! record its position
      DO i = Start+1, End        ! start with next elements
         IF (x(i) < Minimum) THEN    !   if x(i) less than the min?
            Minimum  = x(i)        !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      INTEGER, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1            ! except for the last
         Location = FindMinimum(x, i, Size)    ! find min from this to last
         CALL  Swap(x(i), x(Location))    ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
   
   subroutine percentile_scl(t1,n,q,tq)
    ! q= quantile rechercher EXple: 0.10d0
    ! n= taille du vecteur t
    ! Adapter de la fonction du fichier aaUseFunction.f90

    implicit none
    ! Based on discussion in this paper http://www.haiweb.org/medicineprices/manual/quartiles_iTSS.pdf

    integer::n,ib,i,indd
    double precision::a,b,c,temp
    double precision,dimension(n),intent(in)::t1
    double precision,dimension(n)::t
    double precision,intent(in):: q !quantile recherche(pourcentatge)
    double precision,intent(out)::tq 
    double precision::t25,t975
    
    !!print*,"suis dans la fonction"
    t=t1
    ! n=size(t,1)
    
    ! !print*,"data=",t
    ! tri des temps
    indd=1
    do while (indd.eq.1)
        indd=0
        do i=1,(n-1)
            if (t(i).gt.t(i+1)) then
                temp=t(i)
                t(i)=t(i+1)
                t(i+1)=temp
                indd=1
            end if
        end do
    end do

    ! quantile d'ordre 2.5%
    a=(n-1)*0.025d0
    ! b=mod(a,1.0d0) ! bad formula # 19/11/2018
    b = a - int(a)
    c=a-b
    ib=int(c)
    t25= (1-b)*t(ib+1)+b*t(ib+2)

    ! quantile d'ordre 97.5%
    a=(n-1)*0.975d0
    ! b=mod(a,1.0d0) ! bad formula # 19/11/2018
    b = a - int(a)
    c=a-b
    ib=int(c)
    t975= (1-b)*t(ib+1)+b*t(ib+2)
    
    ! quantile d'ordre q%
    a=(n-1)*dble(q)
    ! b=mod(a,1.0d0) ! bad formula # 19/11/2018
    b = a - int(a)
    c=a-b
    ib=int(c) ! pb: si q = 1, ib = n-1 et donc ib + 2 =n + 1 > n pour la dimension de t
    if(ib <= n-2)then
        tq = (1-b)*t(ib+1)+b*t(ib+2)
    else
        tq = t(n) ! l'on suppose ici qu'on cherche le 100th percentile de la serie, ce qui est normale car dans ce cas, q est tres proche de 1
    endif

    end subroutine percentile_scl
   
   ! generation d'une uniforme dans l'intervalle [a,b]
     ! extrait de la fonction C "runif" du package "stats": fichier runif.c des sources de R.
    
    subroutine runif(a, b, rgener)
        ! a,b : borne de l'interval 
        ! rgener : nombre aleatoire genere
        use var_surrogate, only: random_generator
        implicit none
        double precision, intent(in)::a,b
        double precision, intent(out)::rgener
        double precision::u
        
        if(b < a .or. a < 0 .or. b <0) then
            rgener = -1
        else
            if(a == b) then 
                rgener = a
            else
                if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
                    u = UNIRAN()
                else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
                    CALL RANDOM_NUMBER(u)
                endif
                rgener = a + (b - a) * u
            endif
        endif 
        
        return
        
    end subroutine runif
    
   
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
   
    double precision function uniran()
!
!     Random number generator(RCARRY), adapted from F. James
!     "A Review of Random Number Generators"
!      Comp. Phys. Comm. 60(1990), pp. 329-344.
    
    double precision,save::carry
    double precision,dimension(24),save::seeds 
    double precision,parameter::one=1 
    double precision,parameter::twom24 = ONE/16777216
    integer,save::i,j
    data i, j, carry / 24, 10, 0.0 /
    data seeds / &
    0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635, &
    0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
    0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
    0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/

    uniran = seeds(i) - seeds(j) - carry
    
    if (uniran .lt. 0) then
        uniran = uniran + 1
        carry = twom24
    else
        carry = 0
    end if
    
    seeds(I) = uniran
    I = 24 - MOD( 25-I, 24 )
    J = 24 - MOD( 25-J, 24 )
    
    end function uniran
   
 ! --------------------------------------------------------------------
! SUBROUTINE  simulation():
!    permet de simuler un jeu de donnee a partir d'un modele conjoint.
! --------------------------------------------------------------------
   
   SUBROUTINE simulation(donnee,donneeS,ind_temp,n_col,theta,ksi,betas,alpha,betat,p,prop_i,lambdas,nus,&
                         lambdat,nut,mode_cens,temps_cens,cens0,n_essai,n_obs,&
                         rsqrt,sigma_s,sigma_t,weib,frailty_cor,affiche_stat)
  ! donnee: donnee simulee a retourner pour les deces
  ! donneeS: donnee simulee a retourner pour les surrogate
  ! ind_temp: donne la taille du tableau final complete pour les cas de progression sans deces 
  ! n_col: nombre de colonne du jeu de donnee simulee, vaut 10 si une seule variables explicative: traitement
  ! theta:  variance de la fragilite lognormale au niveau individuel
  ! ksi: parametre de fragilite devant w_ij
  ! betas: effect fixe du traitement sur S
  ! alpha: parametre de fragilite devant v_i, cas des effets aleatoire partage niveau essai
  ! betat:effet fixe du traitement sur T
  ! p: proportion des personnes traitees par essai
  ! prop_i: proportion des sujets par essai
  ! l'idee c'est de retenir les parametres pour lesquelles on obtient une proportion de censure pour la PFS egale a la propotion du jeux de donnes de depart
  ! lambdas,nus,lambdat,nut: Parametres de la weibull pour S et T: 
  ! mode_cens: on utilise les percentiles, soit le taux de censure obtenu des donnees(1),censure fixe(2) et dans ce cas on renseigne le temps de censure
  ! temps_cens: temps de censure considere si censure fixe
  ! cens: proportion des personnes censuree
  ! n_essai:  nombre d'essai a generer
  ! n_obs: nombre total d'observations (ou taille des donnees a simuler)
  ! rsqrt: niveau de correlation souhaite entre les fragilites sepecifiques aux traitement au niveau essai S et T
  ! sigma_s: variance des frailties au niveau essai associee au surrogate
  ! sigma_t: variance des frailties au niveau essai associee au true
  ! weib: indique si on simule par une weibull(1) ou par la loi exponentiel(0)
  ! frailty_cor: indique si l'on considere pour le modele de simulation deux effets aleatoire correles au niveau essai(=1) ou un effet aleatoire partage(=0) ou encore on simule sans effet aleatoire au niveau essai(=2, model conjoint classique)
  ! affiche_stat: dit si l'on affiche les statistiques des donnees simulees(1) ou non (0)
  use var_surrogate, only: random_generator
  
  Implicit none
  
  integer, intent(in)::mode_cens,n_essai,n_obs,weib,frailty_cor,n_col,affiche_stat
  double precision,intent(in)::theta,ksi,betas,alpha,betat,lambdas,nus,lambdat,nut,temps_cens,cens0,rsqrt,sigma_s,sigma_t
  double precision,dimension(n_essai),intent(in)::prop_i,p
  double precision, dimension(n_obs,n_col),intent(out)::donnee
  double precision, dimension(n_obs,n_col),intent(out)::donneeS ! pour les donnees completees surrogate
  integer,intent(out) ::ind_temp
  integer ::k,i,n1,cpte
  integer,parameter ::trt1=1,v_s1=2,v_t1=3,trialref1=4,w_ij1=5,timeS1=6,timeT1=7,&
                      timeC1=8,statusS1=9,statusT1=10,initTime1=11,Patienref1=12 ! definissent les indices du tableau de donnee
  double precision::n_rand,x22,quant_c,t25,t75,t50,cens
  double precision,dimension(n_essai)::n_i
  double precision,dimension(n_obs)::u
  integer,dimension(n_essai,2)::tab ! pour la table de contingence des essais
  
  !==============initialisation des parametres======================
  n_i=NINT(n_obs*prop_i) ! nombre de sujet par essai

  if(sum(n_i)<n_obs) then
    n_i(minloc(n_i,mask=n_i==minval(n_i)))=n_i(minloc(n_i,mask=n_i==minval(n_i)))+(n_obs-sum(n_i)) ! on ajoute au premier essai de plus petite taille le nombre d'individu non encore affecte (1 generalement) a cause des problemes d'arrondi
  endif
  if(sum(n_i)>n_obs) then 
    n_i(maxloc(n_i,mask=n_i==maxval(n_i)))=n_i(maxloc(n_i,mask=n_i==maxval(n_i)))-(sum(n_i)-n_obs) ! on soustrait au premier essai de plus grande taille le nombre d'individu affecte en trop (1 generalement) a cause des problemes d'arrondi
  endif

  !=======================Generation des donnees=============================
  donnee(:,trt1)=2.d0 ! trt
  !donnee=data.frame(trt)
  donnee(:,v_s1)=0.d0 !v_s
  donnee(:,v_t1)=0.d0 !v_t
  donnee(:,trialref1)=0.d0 ! trialref
  donnee(:,initTime1)=0.d0 !initTime1
  
  if(frailty_cor==2)then ! modele conjoint classique: un seul effet aleatoire partage au niveau individuel

! ============simulation randomisation par essai======  
    k=1
    do i=1,n_essai
      !variable traitement
      n1=k+n_i(i)-1
!      do l=k,n1
!        n_rand=uniran()
!        if(n_rand<=p(i)) then
!            donnee(l,trt1)=1.d0 ! attention au sens ici car si on est <p alors on est traite et pas le contraire.
!        else
!            donnee(l,trt1)=0.d0
!        endif
!      end do
      ! reference des essais
      donnee(k:n1,trialref1)=i
      k=k+n_i(i)
    end do
! ============fin randomisation par essai======
    ! ici on randomise sans tenir compte des essai
    k=1
    do i=1,n_obs
        if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
            n_rand=uniran()
        else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
            CALL RANDOM_NUMBER(n_rand)
        endif
        if(n_rand<=p(1)) then  !on suppose que p(1) contient la prportion des traitees
            donnee(i,trt1)=1.d0 ! attention au sens ici car si on est <p alors on est traite et pas le contraire.
        else
            donnee(i,trt1)=0.d0
        endif
        ! reference des essais
        !n1=k+n_i(i)-1
        !donnee(k:n1,trialref1)=i
        !k=k+n_i(i)
    end do

    
    ! fragilites specifiques aux sujets
    x22=0.d0
    !print*,"sqrt(theta)=",sqrt(theta)
    do i=1,n_obs
        !call bgos(theta,0,donnee(i,w_ij1),x22,0.d0) !on simule les w_ij suivant une normale
        call bgos(sqrt(theta),0,donnee(i,w_ij1),x22,0.d0) !on simule les w_ij suivant une normale
    end do
    !!print*,donnee(:,w_ij1)
    !!print*,sum(donnee(:,w_ij1))
    !stop
    ! Generation des temps de suivi (S et T) (voir Austin P.C., statist. Med., 2012, page 3)
    ! nous simulons le risque de base par une loi de weibull (lambda > 0 et nu > 0)
    do i=1,n_obs
        if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
            u(i)=uniran()
        else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
            call RANDOM_NUMBER(u(i))
        endif
    end do
    donnee(:,timeS1)=0.d0 !timeS
    donnee(:,timeT1)=0.d0 !timeT
    if(weib==1)then
      donnee(:,timeS1)=(-dlog(u)/(lambdas*dexp(donnee(:,w_ij1)+betas*donnee(:,trt1))))**(1/nus)
      donnee(:,timeT1)=(-dlog(u)/(lambdat*dexp(ksi*donnee(:,w_ij1)+betat*donnee(:,trt1))))**(1/nut)
    else ! alors on simule suivant la loi exponentielle
      donnee(:,timeS1)=-dlog(u)/(lambdas*dexp(donnee(:,w_ij1)+betas*donnee(:,trt1)))
      donnee(:,timeT1)=-dlog(u)/(lambdat*dexp(ksi*donnee(:,w_ij1)+betat*donnee(:,trt1)))
      !!print*,donnee(1:5,timeS1)
      !!print*,donnee(1:5,timeT1)
      !!print*,"eta=",ksi,"betas=",betas,"betat=",betat,"theta=",theta
      !stop
    endif
    
  endif
    if(affiche_stat==1) then ! on affiche les statistique des donnees simulees
        tab=table(donnee(:,trialref1),n_essai)
        do i=1, n_essai
            !print*,tab(i,1),tab(i,2)
        end do
    endif

  ! calcul du temps de censure: la censure ne depend pas du type d'evenement, c'est par rapport au deces
  if(mode_cens==1)then ! cenure a l'aide des percentiles obtenus des donnees
    !cens= prop.table(table(don$statusT=="0"))(2) !proportion des personnes censurees
    cens= cens0 !proportion des personnes censurees
    !quant_c=min(quantile(donnee$timeT,probs = 1-cens),temps_cens)
    call percentile_scl(donnee(:,timeT1),n_obs,1.d0-cens,quant_c)
    !print*,"1.d0-cens",1.d0-cens,"cens",cens
    !!print*,"quant_c",quant_c
    !stop
  else ! censure fixe
    quant_c=temps_cens
  endif
  
  donnee(:,timeC1)=quant_c
  donnee(:,statusS1)=-1.d0
  donnee(:,statusT1)=-1.d0
  !!print*,"min tempsS",minval(donnee(:,timeS1))
  !!print*,"min tempsT",minval(donnee(:,timeT1))
  !!print*,"max tempsS",maxval(donnee(:,timeS1))
  !!print*,"max tempsT",maxval(donnee(:,timeT1))
  
  !=============Description des caracteristiques des temps de suivi simules============
    if(affiche_stat==1) then ! on affiche les statistique des donnees simulees
        call percentile_scl(donnee(:,timeS1),n_obs,0.25d0,t25)
        call percentile_scl(donnee(:,timeS1),n_obs,0.50d0,t50)
        call percentile_scl(donnee(:,timeS1),n_obs,0.75d0,t75)
  
        !print*,"summary:    Minimum        25%        50%        75%        Maximum"
        !print*,"summary tempsS",minval(donnee(:,timeS1)),t25,t50,t75,maxval(donnee(:,timeS1)) 
  
        call percentile_scl(donnee(:,timeT1),n_obs,0.25d0,t25)
        call percentile_scl(donnee(:,timeT1),n_obs,0.50d0,t50)
        call percentile_scl(donnee(:,timeT1),n_obs,0.75d0,t75)
  
        !print*,"summary tempsT",minval(donnee(:,timeT1)),t25,t50,t75,maxval(donnee(:,timeT1))

        !print*,"le temps de censure vaut:",quant_c
    endif
    !!print*,nint(donnee(:,statusS1))
    !!print*,nint(donnee(:,statusT1))
    ! construction de la variable statut (Surrogate ou true)
    ind_temp=1 ! pour avoir le nombre d'observations dans le jeux de donnees donneeS
    cpte=0
    do i=1,n_obs
        donnee(i,Patienref1)=i
        !on construit les temps de deces et statut deces
        if(donnee(i,timeT1)<=donnee(i,timeC1))then ! patient decede
            donnee(i,statusT1)=1.d0
        else    !patient censuree administrativement
            donnee(i,statusT1)=0.d0
            donnee(i,timeT1)=donnee(i,timeC1)
        endif
        
        !on construit les temps de progression
        if(donnee(i,timeS1)<donnee(i,timeT1))then
            cpte=cpte+1
            donnee(i,statusS1)=1.d0
            !on complete une ligne pour la censure 
            donneeS(ind_temp,:)=donnee(i,:)
!            donneeS(ind_temp+1,:)=donnee(i,:)
!            donneeS(ind_temp+1,initTime1)=0.d0 !initTime vaut la date d'evenement
!            donneeS(ind_temp+1,timeS1)=donnee(i,timeT1)-donnee(i,timeS1) ! temps de survie avant deces depuis la progression
!            donneeS(ind_temp+1,statusS1)=0.d0 ! on ne considere que le premier evenement
!            donneeS(ind_temp+1,Patienref1)=donnee(i,Patienref1) !
!            ind_temp=ind_temp+2 ! +2 a cause de la nouvelle ligne augmentee
            ind_temp=ind_temp+1 ! +2 a cause de la nouvelle ligne augmentee
        else
            if((donnee(i,timeS1)==donnee(i,timeT1)).and.(donnee(i,statusT1)==0.d0)) then !evenement a la date de censure
                donnee(i,statusS1)=1.d0
                donneeS(ind_temp,:)=donnee(i,:)
                ind_temp=ind_temp+1
                cpte=cpte+1
            else ! progression le meme jour que le deces ou sans progression
                donnee(i,statusS1)=0.d0             ! on suppose pas d'evenement si le meme jour que le deces
                donnee(i,timeS1)=donnee(i,timeT1)! et on censure a la date de deces(ou censure)
                donneeS(ind_temp,:)=donnee(i,:)
                ind_temp=ind_temp+1
                cpte=cpte+1
            endif
        endif
    enddo
    !!print*,"cpte=",cpte
    !stop
    !!print*,nint(donnee(:,statusT1))
    !stop
     !affichage de quelques lignes
    if(affiche_stat==1) then
        !print*,"surrogate"
        do i=1,10
            !print*,(donneeS(i,6:n_col))
        enddo
        !print*,"Deces"
        do i=1,10
            !print*,(donnee(i,6:n_col))
        enddo
    endif
    !stop
    ind_temp=ind_temp-1 !pour rester vec la taille souhaitee
    !on ajoute 0.5 aux temps avant d'arrondir pour eviter le min des temps a 0, ie pas suivi
    !donnee(:,timeS1)=nint(donnee(:,timeS1)+0.5d0) ! arrondi
    !donnee(:,timeT1)=nint(donnee(:,timeT1)+0.5d0) !arrondi    
    !donneeS(:,timeS1)=nint(donneeS(:,timeS1)+0.5d0) ! arrondi
    !donneeS(:,timeT1)=nint(donneeS(:,timeT1)+0.5d0) !arrondi    
!    do i=1,n_obs
!    if(donnee(i,timeC1)<min(donnee(i,timeS1),donnee(i,timeT1)))then ! censure
!      donnee(i,statusS1)=0.d0
!      donnee(i,statusT1)=0.d0
!      donnee(i,timeS1)=donnee(i,timeC1)
!      donnee(i,timeT1)=donnee(i,timeC1)
!   goto 1
!    else
!        if((donnee(i,timeS1)<=donnee(i,timeC1)).and.(donnee(i,timeT1)>donnee(i,timeC1))) then! progression + censure
!            donnee(i,statusS1)=1.d0
!            donnee(i,statusT1)=0.d0
!            donnee(i,timeT1)=donnee(i,timeC1) ! on censure l'individu
!            goto 1
!        else !
!            if((donnee(i,timeS1)>donnee(i,timeC1)).and.(donnee(i,timeT1)<=donnee(i,timeC1)))then ! deces sans progression
!                donnee(i,statusS1)=0.d0
!                donnee(i,statusT1)=1.d0
!                donnee(i,timeS1)=donnee(i,timeT1) ! on censure l'individu a la date de deces
!                goto 1
!            else ! (donnee(i,timeS<=donnee(i,timeC)&<=(donnee(i,timeTdonnee(i,timeC), les deux evenements se sont produits
!                if(donnee(i,timeS1)<=donnee(i,timeT1))then ! progression + deces  
!                    donnee(i,statusS1)=1.d0
!                    donnee(i,statusT1)=1.d0
!                    goto 1
!                else ! deces sans progression
!                    if(donnee(i,timeS1)>donnee(i,timeT1))then
!                        donnee(i,statusS1)=0.d0
!                        donnee(i,statusT1)=1.d0
!                        donnee(i,timeS1)=donnee(i,timeT1) ! on censure l'individu a la date de deces
!                    endif
!                endif
!            endif
!        endif
!    endif
!    1    continue
!  end do
  return
end subroutine simulation

! fonction permettant de calculer la variance d'un echantillon: on divise par n-1 pour avoir un estimateur non biaise

function variance(array)
    implicit none
    double precision::variance,x_bar
    double precision,dimension(:),intent(in)::array
    x_bar=sum(array)/dble(size(array))
    ! !print*,"suis dans la variance"
    ! if(size(array)==1)then
        ! !print*,"Attention: variance d'un seul element recherché"
        ! variance=0.d0
        ! return
    ! else
        variance=sum((array-x_bar)**2.d0)/dble(size(array)-1)
        ! return
    ! endif
end function variance

! fonction permettant de calculer la covariance d'un echantillon: on divise par n-1 pour avoir un estimateur non biaise
function covariance(x,y)
    implicit none
    integer ::m
    double precision::covariance,xmn,ymn
    double precision,dimension(:),intent(in)::x,y
    double precision,dimension(size(x))::xdev,ydev
    
    m=size(x)
    xmn = SUM(x) / dble(m)
    ymn = SUM(y) / dble(m)
    xdev = x - xmn
    ydev = y - ymn    
    covariance = SUM( xdev * ydev) / dble(m-1)
endfunction

! fonction pour la table de contingence appliquee a un vecteur
function table(tab,n)
    !n: nombre de modalites
    implicit none
    integer,intent(in)::n
    double precision,dimension(:),intent(in)::tab
    integer,dimension(n,2)::table
    integer,dimension(size(tab))::t1
    integer::i,j,n_obs
    
    n_obs=size(tab)
    t1=0
    do i=1,n_obs
        j=nint(tab(i))
        t1(j)=t1(j)+1
    end do
!    n=count(t1==0)
    j=1
    do i=1,n_obs
        if(t1(i).ne.0)then
            table(j,1)=i
            table(j,2)=t1(i)
            j=j+1
        endif
    end do
    
end function table

! fonction pour la table de contingence appliquee a un vecteur pour les essais
function table_essai(tab)
    !n: les case avec 0 correspondent aux modalites sans effectif
    implicit none
    !integer::n
    integer,dimension(:),intent(in)::tab
    integer,dimension(size(tab))::t2,table_essai
    integer::i
        
    t2=0
    do i=1,size(tab)
        t2(tab(i))=t2(tab(i))+1
    enddo
    table_essai=t2    
        
end function table_essai

!c===================================   GAMGUI    ============================

    subroutine gamgui(a,x)
        use var_surrogate, only: random_generator
        double precision :: a,b,c,u,v,w,x,y,z
        double precision ::uniran
        !real ::ran2
        integer ::accept

        accept = 0
        b = a-1.
        c = 3.*a - 0.75
 1      if(accept.eq.0) then
            if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
                u = uniran()!dble(rand())
                v = uniran()!dble(rand())
            else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
                CALL RANDOM_NUMBER(u)
                CALL RANDOM_NUMBER(v)
            endif
           w = u*(1.-u)
           y = sqrt(c/w)*(u-0.5)
           x = b+y
           if (x.ge.0.) then
              z = 64.*w*w*w*v*v
              if(z.le.(1.-2.*y*y/x).or.log(z)&
           .le.2.*(b*log(x/b)-y)) accept = 1
              goto 1
           else
              goto 1
           endif   
        else
           return
        endif   
    return
    endsubroutine gamgui
    
    !c===================================   WEIGUI2    ============================

    subroutine weigui2(a,b,betau,x)
!c fonction de densité de la loi de weibull = f(x)=b**a . a . x**(a-1) . exp(-(bx)**a) (voir cours de Piere Jolie page 41)
        use var_surrogate, only:param_weibull
        use var_surrogate, only: random_generator
        double precision ::a,b,x,u,v,betau
        double precision ::uniran
        !real ::ran2
        if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
            u = uniran()!dble(rand())
        else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)    
            CALL RANDOM_NUMBER(u)
        endif
        v = (1.d0-u)
        if(param_weibull==0)then !parametrisation de weibull par defaut dans le programme de Virginie: fonction de densite differente de celle donnee ci-dessus
            x = (1.d0/b)*((-dexp(-betau)*dlog(v))**(1.d0/a))
        else 
            ! parametrisation de la weibull donnee par la fonction de densite ci-dessus
            ! F^-1(t)=1/b*(-log(1-t))**(1/a)
            x=    (1.d0/b)*((-dlog(1+ dlog(u)*dexp(-betau)))**(1.d0/a))
        endif

    return
    endsubroutine weigui2
    
    subroutine weiguicopule(a,at,b,bt,betau,betaut,theta,Sij,Tij)
    !theta : parametre de copula
    ! fonction de densité de la loi de weibull = f(x)=b.a. x**(a-1)
    ! a, at parametres d'echelle, b, bt : parametres de forme
    ! generation du temps de deces conditionnellement au surrogate: T_ij|S_ij
        use var_surrogate, only:param_weibull
        use var_surrogate, only: random_generator
        double precision ::a,b,at,bt,Sij,Tij,u,ut,v,betau,betaut,vt,utij,vtij,theta,wij
        double precision ::uniran
        real ::ran2
        if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
            u = uniran()
            ut = uniran()
        else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)    
            CALL RANDOM_NUMBER(u)
            CALL RANDOM_NUMBER(ut)
        endif
        v = (1.d0-u)
        vt = (1.d0-ut)
        wij = (1.d0 - u)**(- theta)
        vtij = 1.d0 - wij + wij * vt**(-theta/(1.d0 + theta))
        
       !parametrisation de weibull par defaut dans le programme de Virginie: fonction de densite differente de celle donnee ci-dessus, idem a Takeshi
        Sij = ((1.d0/b)*(-dexp(-betau)*dlog(v)))**(1.d0/a)
        Tij = ((1.d0/(theta *bt))*dexp(-betaut)*dlog(vtij))**(1.d0/at)
        
        ! call dblepr("theta", -1, theta, 1)
        ! call dblepr("bt", -1, bt, 1)
        ! call dblepr("at", -1, at, 1)
        ! call dblepr("vtij", -1, vtij, 1)
        ! call dblepr("betaut", -1, betaut, 1)
        ! call dblepr("Sij", -1, Sij, 1)
        ! call dblepr("Tij", -1, Tij, 1)
    return
    endsubroutine weiguicopule
    
    !==============================================================================================================================
    ! subroutine pour la generation  des donnees pour une distribution gamma des effects aleatoires et weibull des risques de bases
    !==============================================================================================================================
    subroutine generation_Gamma(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
            ng,ver,truealpha,propC,cens_A,gamma1,gamma2,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,betat)
      !lognormal: dit si la distribution des effets aleatoires est lognormal (1) ou gamma (0)
      !use Autres_fonctions
      use var_surrogate, only: random_generator
      
      integer, intent(in)::n_obs,n_col,lognormal,ng,ver
      double precision, intent(in)::truealpha,propC,cens_A,gamma1,gamma2,theta2 !theta2: variance des frailties gaussiens,gamma1,gamma2: parametres de la gamma,cens_A:censure administrative,propC:proportion de personnes censurees
      double precision, intent(in)::lambda_S,nu_S,lambda_T,nu_T,betas,betat
      double precision,dimension(n_obs,n_col),intent(out)::don_simulS1,don_simul
      integer, parameter::npmax=70,NOBSMAX=15000,nvarmax=45,ngmax=5000
      integer,parameter::nboumax=1000,NSIMAX=5000,ndatemax=30000
      
      integer  j,k,nz,cpt,cpt_dc,ii,iii,iii2
      integer  cptstr1,cptstr2
      !integer  i,ic,ic2,ni,ier,istop,ef
      integer  cptni,cptni1,cptni2,nb_echec,nb_echecor
      integer  nbou2,cptbiais
      integer  nbou
      integer filtre(nvarmax), filtre2(nvarmax) 
      !integer cpt1(nboumax) 
      !integer cpt2(nboumax) 
      !integer cpt3(nboumax) 
      integer ind(nboumax) , cptaux 
      
      !real vax(nvarmax)
      !double precision tt0
      !double precision tt1,tt2
      
      !double precision h
      !double precision ro,wres,csi,csi1
      double precision maxtemps
      !double precision bi,bs,wald,str
      !double precision moyvar1,moyse1,moyse_cor1!theta
      !double precision moyvar2,moyse2,moyse_cor2!alpha
      !double precision moybeta1, moysebeta1,moysebeta_cor1
      !double precision moybeta2, moysebeta2,moysebeta_cor2
      !double precision moybeta3, moysebeta3,moysebeta_cor3
      !double precision moybweib1,moybweib2,moybweib3,moybweib4!weibull
      !double precision moysebweib1,moysebweib2,moysebweib3,moysebweib4
!c     double precision  tt0,tt1
!c     double precision, pointer  tt0
!c     double precision,pointer  tt1
      !double precision varxij,eca,varsmarg,smoy,smoyxij
      double precision lrs
      double precision BIAIS_moy
      
      !double precision aux(2*NOBSMAX)
      !double precision v((npmax*(npmax+3)/2))
      !double precision k0(2)
      !double precision b(npmax)
      !double precision se1(nboumax),se_cor1(nboumax)!theta
      !double precision se2(nboumax),se_cor2(nboumax)!alpha
      !double precision tvars1(nboumax),tvars2(nboumax)
      !double precision beta1(nboumax),beta2(nboumax),beta3(nboumax)
      !double precision , dimension(nboumax):: sebeta1,sebeta_cor1
      !double precision , dimension(nboumax):: sebeta2,sebeta_cor2
      !double precision , dimension(nboumax):: sebeta3,sebeta_cor3

      !double precision bweib1(nboumax),bweib2(nboumax)!weibull
      !double precision bweib3(nboumax),bweib4(nboumax)
      !double precision , dimension(nboumax):: sebweib1,sebweib2
      !double precision , dimension(nboumax):: sebweib3,sebweib4

      !double precision biais_theta(nboumax)
      !double precision I1_hess(npmax,npmax),H1_hess(npmax,npmax)
      !double precision I2_hess(npmax,npmax),H2_hess(npmax,npmax)
      !double precision HI1,HI2(npmax,npmax)
      !double precision HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax)
      !double precision BIAIS(npmax,1)
      
      character(18) :: nomvarl
      character(18) :: nomvar(nvarmax),nomvar2(nvarmax)
      !character(18) :: donnees
      character(24) :: ficpar
      !character(14) :: fich1
      !character(14) :: fich2
      !character(14) :: fich3
      !character(14) :: fich4
      !character(14) :: fich1b
      !character(14) :: fich2b
      !character(14) :: fich3b
      !character(14) :: fich4b
      character(20) :: dateamj
      character(20) :: zone
      character(20) :: heure1
      !character(20) :: heure2
      integer values(8)

!c************declarations pour donnees generees **********
      integer :: nb_recur,nb_dc,nb_cens,delta,deltadc,jj
      integer :: ig,nrecurr,nobs,max_recu
      real , dimension(2):: v1
      real :: piece,rien,demi
      double precision :: ui,temps1,temps1_S !random effect
      double precision :: gapx,gapdc,moy_idnum
      double precision :: x,xdc,cens,cbeta1,cbeta3
      double precision :: auxbeta1,auxbeta2 ! for recurr and death
      double precision :: uniran
      double precision, dimension(2):: bg1,bw1,bw2
      integer, dimension(ngmax):: idnum
      double precision, dimension(ngmax):: vecui
      double precision :: moyui
!c*****************************************************************
      
!c*****************************************************************
!c***** nmax
         integer :: nmax
         common /nmax/nmax
!c*****dace1 
      double precision date(ndatemax)
      double precision zi(-2:npmax)
      common /dace1/date,zi
      
!c*****dace2
      double precision t0(NOBSMAX),t1(NOBSMAX),t1_S(NOBSMAX)
      integer c(NOBSMAX), cdc(NOBSMAX)
      !integer nt0(NOBSMAX),nt1(NOBSMAX)
      integer  nva1,nva2,nst!,nobs
      !common /dace2/t0,t1,c,cdc,nt0,nt1,nobs,nva,nva1,nva2,ndate,nst
!c*****dace4
      !integer  stra(NOBSMAX)
      !common /dace4/stra
!c*****ve1
      double precision ve(NOBSMAX,nvarmax),ve2(NOBSMAX,nvarmax)
      !common /ve1/ve,ve2
!c*****dace3
      !double precision  pe
      integer  effet,nz1,nz2
      !common /dace3/pe,effet,nz1,nz2
!c*****dace7
      !double precision I_hess(npmax,npmax),H_hess(npmax,npmax)
      !double precision Hspl_hess(npmax,npmax)
      !double precision PEN_deri(npmax,1)
      !double precision hess(npmax,npmax)
      !common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
!c*****contrib
      !common /contrib/ng       
!c*****groupe
      integer g(NOBSMAX)
      integer nig(ngmax)        ! nb d events recurrents , different de mi()!
      !common /gpe/g,nig
      
!c*****mem1
      !double precision mm3(ndatemax),mm2(ndatemax)
      !double precision mm1(ndatemax),mm(ndatemax)
      !common /mem1/mm3,mm2,mm1,mm
!c     %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      !integer AG
      !common /andersengill/AG
!c     %%%%%%%%%%%%% indic ALPHA %%%%%%%%%%%%%%%%%%%%%%%%% 
      integer indic_ALPHA
      !common /alpha/indic_ALPHA ! pour preciser un para en plus 
!c**** theta/alpha
      !double precision  theta,alpha !en exposant pour la frailty deces 
      !common /thetaalpha/alpha,theta
!c******indicateur de troncature
      !integer :: indictronq     ! =0 si donnees non tronquées reellement
      !common /troncature/indictronq
      
!c******indicateur iteration
      integer  ibou
      !common /boucle/ibou
!c************FIN COMMON ***********************************
      
!c     nst: deux finctions de risque a estimer (meme bases de splines)
!c     ist: appartenance aux strates
!c icen=1: censure
!c     ib: matrice de variance de beta chapeau
!c     I_hess : -hessienne non inversee sur vraisemblance non penalisee
!c     H_hess : inverse de -hessienne  sur vraisemblance penalisee
      !Ajout SCL
      double precision::x22,tempon!,theta2
      double precision,dimension(:),allocatable::tempsD
      integer::Aretenir
      integer,parameter ::trt1=1,v_s1=2,v_t1=3,trialref1=4,w_ij1=5,timeS1=6,timeT1=7,&
                      timeC1=8,statusS1=9,statusT1=10,initTime1=11,Patienref1=12 ! definissent les indices du tableau de donnee simulees
      integer,intent(in)::affiche_stat
      double precision,intent(out)::vrai_theta
!CCCCCCCCCCCCCCCCChosur9.f CCCCCCCCCCCCCCCCCCCCCCCC
    if(affiche_stat==1)then
      !write(*,*)'    ******************************************'
      !write(*,*)'  ****** DEBUT PROGRAMME FRAILTY.F**********'
      !write(*,*)'******************************************'
    endif
      nmax = 300 !nb iterations max dans marquard
      
      indic_alpha=1 ! on precise que l on a un parametre en plus estimer
!c      allocate (tt0)
!c      allocate (tt1)
      
      call date_and_time(dateamj,heure1,zone,values)                  
!c     !write(*,*)'Starting time: ', dateamj,heure1,zone,values
      
      lrs=0.d0      
      nb_echec=0
      nb_echecor=0
      nbou2=0
      ficpar='joint2.inf'
      !open(2,file=ficpar)
      !open(4,file='outjoint')
      
      !read(2,*)ng
      !read(2,*)nrecurr !nb de dobservations recurrentes max par sujet
      nrecurr=1
      !read(2,*)nbou
      nbou=1
    if(affiche_stat==1)then
      !write(4,*)'** nb de simulations',nbou
      !write(*,*)'** nb de simulations',nbou
    endif
      
      allocate(tempsD(ng))
      nst=2
    
    if(affiche_stat==1)then
      !write(4,*)'**************************************************'
      !write(4,*)'************ JOINT MODEL *************************'
      !write(4,*)'*** RECURRENT EVENTS and TERMINATING EVENT *******'
      !write(4,*)'**************************************************'
      !write(4,*)'** nb de groupes=sujets =',ng 
      !write(*,*)'** nb de groupes=sujets =',ng 
      !write(4,*)'** deux fonctions de risque de base  = ',nst
      !write(*,*)'** fonction de risque de base  = ',nst
    endif

!c     !write(4,*)'** nb de simulations = ',nbou
!c     !write(4,*)'** indicateur de censure (1 censure existante)',icen
        cptni=0
        cptni1=0
        cptni2=0
        biais_moy=0.d0
        cptbiais=0
        cptaux=0

!c**************************************************
!c********************* prog spline****************
         !read(2,*) !scl pour faire passer a la ligne suivante
         !read(2,*)ver

         nva1 = 0 ! nb de var expli pour donnees recurrentes
         nva2 = 0 ! nb de var expli pour deces
         !!print*,"ver=",ver

         if(ver.gt.0)then
            do 44 j=1,ver
               !read(2,*)nomvarl,filtre(j) ,filtre2(j)
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
    if(affiche_stat==1)then
      !write(4,*)'** explanatory variables for recurrent events:',nva1
      !write(*,*)'** explanatory variables for recurrent events:',nva1
      !write(4,*)'** explanatory variables for deaths:',nva2
      !write(*,*)'** explanatory variables for deaths:',nva2
    endif
    
    !  read(2,*)truealpha 
      
    if(affiche_stat==1)then
      !write(4,*)'** Vraie valeur de Alpha',truealpha 
      !write(*,*)'** Vraie valeur de Alpha',truealpha 
    endif

!c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
!      read(2,*)AG
!c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      read(2,*)nz
        nz=6
      
      nz1=nz
      nz2=nz
      if(nz.gt.20)then
         nz = 20
      endif
      if(nz.lt.4)then
         nz = 4
      endif
      
    if(affiche_stat==1)then
     !write(4,*)'nombre de noeuds :',nz
    endif
    
!C deux parametres de lissage :
         !read(2,*)ax1
        ! ax1=94000
         !if(nst.eq.2)then
            !read(2,*)ax2
        !    ax2=135000
        ! endif
         !   !write(4,*)'kappa1',ax1
         !   if(nst.eq.2)then
        !        if(affiche_stat==1)then
        !            !write(4,*)'kappa2',ax2
        !            !write(4,*)'nombre de noeuds:',nz
        !        endif
         !   endif
         !read(2,*)fich1
        ! if(nst.eq.2)then
            !read(2,*)fich2
         !endif
        !read(2,*)fich3 
        ! if(nst.eq.2)then
            !read(2,*)fich4
        ! endif
         !read(2,*)fich1b
        ! if(nst.eq.2)then
            !read(2,*)fich2b
        ! endif
         !read(2,*)fich3b 
        ! if(nst.eq.2)then
           ! read(2,*)fich4b
         !endif
!c         close(2)
        ! k0(1) = ax1
         !if(nst.eq.2)then
         !   k0(2) = ax2
        ! endif
        ! read(2,*)propC ! scl= proportion de personnes censurees
        ! read(2,*)cens_A ! censure administrative
        ! read(2,*)Aretenir! Jeu de donnee a retenir dans les simulations pour le test
        Aretenir=1
        
        if(Aretenir>nbou) then
            !print*,"l'indice du jeu de données à retenitr doit etre inferieure aunombre de simulation"
            !stop
        endif
!c         !write(*,*),'lissage',k0(1)
!c         !write(*,*),'lissage2',k0(2)

!c*******************************************
!c***** DEBUT SIMULATIONS *******************
!c*******************************************
         maxtemps = 0.d0
         cpt = 0
         cpt_dc = 0
         cptstr1 = 0
         cptstr2 = 0
                 !!print*,"=1suislà"
        !do jj=1,ng
         don_simulS1=0.d0
         don_simul=0.d0
        !enddo
                 !!print*,"=suislà"
         do 1000 ibou=1,nbou 
            if(affiche_stat==1)then
                !write(*,*)'************ Generation NUMERO :      ',ibou
            endif
           ind(ibou)=0
           effet = 1
           nig = 0
           g=0
           maxtemps = 0.d0
           auxbeta1=0.d0
           auxbeta2=0.d0

!c********************************************************
!c******** DEBUT generation des donnees *******************
!c********************************************************

        !open(10,file="parametre_2007.inf")
!!c------------------ UI = FRAILTY ---------------------
!c 2 parametres de la gamma tq E(ui)=bg1(1)/bg1(2)  var(ui)=bg1(1)/(bg1(2)*bg1(2))
           !read(10,*)bg1(1)
            bg1(1)=gamma1
           !read(10,*)bg1(2)
            bg1(2)=gamma1
           !read(10,*)theta2 ! variance des frailties gaussien wij
           !!print*,"=====================",bg1(1),bg1(2)
           if(ibou.eq.1)then
                !write(4,*)'** '
                if (lognormal==0)then ! Gamma
                    if(affiche_stat==1)then    
                        !write(4,*)'** vraie valeur de theta = **'&
                         !   ,bg1(1)/(bg1(2)*bg1(2))
                        !write(*,*)'** vraie valeur de theta = **'&
                         !   ,bg1(1)/(bg1(2)*bg1(2))
                    endif
                    vrai_theta=bg1(1)/(bg1(2)*bg1(2))
                else !lognormale
                    if(affiche_stat==1)then    
                        !write(4,*)'** vraie valeur de theta = **'&
                         !   ,theta2
                        !write(*,*)'** vraie valeur de theta = **'&
                         !   ,theta2
                    endif
                    vrai_theta=theta2
                endif
           endif

!!c---------------- X --------------------

!!c-- parametres de la WEibull for recurrent events 
         bw1(1)=lambda_S 
         bw1(2)=nu_S
!!c-- parametres de la WEibull for death
         bw2(1)=lambda_T
         bw2(2)=nu_T

         demi = 0.5             ! pour var expli
         cbeta1=betas          
         cbeta3=betas        
         
         if(affiche_stat==1)then
            if(ibou.eq.1)then
                !write(4,*)'**vraie valeur de beta1 (recurrent) = **',cbeta1
                !write(4,*)'**vraie valeur de beta2  (recurrent) = **',cbeta2
                !write(4,*)'** vraie valeur de beta3 (deces) = **',cbeta3
            endif
        endif
         nb_recur = 0           ! nb de temps 
         nb_dc = 0              ! nb de dc
         nb_cens = 0            ! nb de cens
         nobs=0
         max_recu= 1
         moyui = 0.d0
        !close(10)
!!c--------------------------------------------------------------
!!c--------------------------------------------------------------
!!c--------------------------------------------------------------

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
            !call gamgui(bg1(1),ui) !genere gamma pour zi
            !ui = ui/bg1(2)
!c     verification des ui
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
               piece=real(tempon) !rand()
                !piece=real(uniran()) !rand()
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
!c----------- CENSORING --------------------------------------------
                  !cens =  1.d0 + 250.d0*uniran() !scl voir plubas
              
!c-----------RECURRENT --------------------------------------------
!c---  genere temps recurrents a partir 2 var explic :
    if (lognormal==0)then ! Gamma
        auxbeta1=dlog(ui)+cbeta1*dble(v1(1))!+cbeta2*dble(v1(2)) ! scl je considere uniquement le traitement
        auxbeta2=truealpha*dlog(ui)+cbeta3*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
    else ! lognormale
        auxbeta1=ui+cbeta1*dble(v1(1))!+cbeta2*dble(v1(2)) ! scl je considere uniquement le traitement
        auxbeta2=truealpha*ui+cbeta3*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
    endif

        call weigui2(bw1(1),bw1(2),auxbeta1,gapx)
!c**************gap time: 
        x=gapx
                 
!c-----------DECES --------------------------------------------
!c---     genere temps de dc a partir d une var explic :
        call weigui2(bw2(1),bw2(2),auxbeta2,gapdc) 
!c************** calendar time: 
!c               xdc=xdc+gapdc
!c**************  gap time: 
        xdc=gapdc
        tempsD(ig)=xdc
                 !!print*,"xdc=",xdc,tempsD(i)
                 
         
! scl============censure====================
        cens=cens_A
!c------------------------------------------------bilan:
                ! if ((xdc.le.x).and.(xdc.le.cens)) then !deces
                !    deltadc = 1
                !    delta = 0
                !    temps1 = xdc
                !    nb_dc =nb_dc + 1
                ! endif
                
                if(xdc.le.cens)then ! patient decede
                    deltadc=1.d0
                    temps1 = xdc
                    nb_dc =nb_dc + 1
                else    !patient censuree administrativement
                    deltadc=0.d0
                    temps1 = cens
                    nb_cens =nb_cens + 1
                endif


                ! if ((cens.le.x).and.(cens.le.xdc)) then ! censure
                !    deltadc = 0
                !    delta = 0
                 !   temps1 = cens
                 !   nb_cens =nb_cens + 1
                ! endif

                ! if ((x.le.cens).and.(x.le.xdc)) then ! recurrent event
                !    deltadc = 0
                !    delta = 1
                !   temps1 = x
                !    nb_recur =nb_recur + 1
                !    nig(ig) = nig(ig)+1 !nb events recurrents
                ! endif  
                 
        
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

!c*** pour le tester sur un autre programme 
            don_simulS1(ig,initTime1)=t0(nobs)
            don_simulS1(ig,timeS1)=t1_S(nobs)
            don_simulS1(ig,statusS1)=c(nobs)
            don_simulS1(ig,trialref1)=1
            don_simulS1(ig,Patienref1)=g(nobs)
            don_simulS1(ig,trt1)=ve2(nobs,1)
            don_simulS1(ig,w_ij1)=ui
            
            don_simul(ig,initTime1)=t0(nobs)
            don_simul(ig,timeT1)=t1(nobs)
            don_simul(ig,statusT1)=cdc(nobs)
            don_simul(ig,trialref1)=1
            don_simul(ig,Patienref1)=g(nobs)
            don_simul(ig,trt1)=ve2(nobs,1)
            don_simul(ig,w_ij1)=ui
            
            ! if(ibou.eq.Aretenir)then
            !    open(9,file="parametres.txt")
            !    !write(9,112)(t0(nobs)),t1_S(nobs),c(nobs),t1(nobs),cdc(nobs)&
            !       ,g(nobs),int(ve(nobs,1))&
            !        ,int(ve(nobs,2)),int(ve2(nobs,1))
! 112                 format(f7.1,f7.1,I2,f7.1,I2,I5,I2,I2,I2)
               !!print*,"size(g)",size(g)
            !    open(13,file="gastadv_T.txt")
            !    !write(13,113)1,g(nobs),int(ve2(nobs,1)),int(t0(nobs)),t1(nobs),cdc(nobs)
                    
! 113            format(I5,I5,I2,I2,f7.1,I2)
            !    open(8,file="gastadv_S.txt")
            !    !write(8,1131)1,g(nobs),int(ve(nobs,1)),int(t0(nobs)),t1_S(nobs),c(nobs)!&
                    !,int(ve(nobs,2))
 !1131           format(I5,I5,I2,I2,f7.1,I2)
            ! endif
!c******************************************************************
!c            
            if (maxtemps.lt.t1(nobs))then
               maxtemps = t1(nobs)
            endif  
!c         !write(*,*)'*données**',t0(ig),t1(ig),c(ig),ve(ig,1),ve(ig,2),ig

            if (delta.eq.0.d0)then ! deces ou censure
               goto 30          ! on change de sujet
            endif
 10      continue                   !observations par sujet
 30   continue                  !sujets=groupes
 !                         !print*,"suis la"
      if(ibou.eq.Aretenir)then
        if(affiche_stat==1)then
            !write(*,*)'** nombre total d observations',nobs
            !write(*,*)"** nb d'evenements surrogate",nb_recur
            !write(*,*)'** nb donnees censurees ',nb_cens
            !write(*,*)'** nb deces ',nb_dc

            !write(*,*)'** proportion de deces (en %) ',nb_dc*100.d0/nobs
            !write(*,*)'** proportion de surrogate (en %) ',nb_recur*100.d0/nobs     
            !write(4,*)'** nombre total d observations',nobs
            !write(4,*)"** nb d'evenements surrogate",nb_recur
            !write(4,*)'** nb donnees censurees ',nb_cens
            !write(4,*)'** nb deces ',nb_dc
            !write(4,*)'** proportion de surrogate (en %) ',nb_recur*100.d0/nobs
            !write(4,*)'** proportion de deces (en %) ',nb_dc*100.d0/nobs
        endif
      moy_idnum=0
      do 444 jj=1,ng
         moy_idnum =moy_idnum + idnum(jj)
 444     continue
          moy_idnum =moy_idnum / ng
    !  !write(*,*)'** nb moyen d observations par sujet ',moy_idnum,ng
    !  !write(4,*)'** nb moyen d observations par sujet ',moy_idnum
      
      endif
!c****************************************************
!c******** FIN generation des donnees****************
!c****************************************************

!c=========     on retourne a l'iteration suivante
 1000 continue
    ! scl============censure conseillee pour la proportion souhaitee de censure==
        !!print*,tempsD
        call percentile_scl(tempsD,ng,1.d0-propC,cens)
        if(affiche_stat==1)then
            !print*,"la censure conseillée pour avoir:",propC*100,"% de personnes censurée vaut:",cens
            !print*,"la max de temps de deces vaux:",maxval(tempsD)
        endif
    deallocate(tempsD)
    !close(2)
    !close(4)
endsubroutine generation_Gamma !FIN prog principal



!==============================================================================================================================
!     subroutine pour la generation  des donnees du modèle surrogate complet avec un effet aleatoire partage au niveau individuel 
!    et 2 effets correles au niveau essai en interaction avec le traitement
!==============================================================================================================================

subroutine Generation_surrogate(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
            ng,ver,truealpha,propC,cens_A,gamma1,gamma2,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,&
            betat,n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma,alpha,frailt_base, pfs)
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
    ! pfs : used to specified if the time to progression should be censored by the death time (0) or not (1). The default is 0 as in sofeu et al. (2019). In this case, death is not included in the surrogate endpoint. 
     use var_surrogate, only: random_generator

      integer, intent(in)::n_essai,frailt_base,affiche_stat,n_obs,n_col,lognormal,ng,ver, pfs
      double precision, intent(in)::truealpha,propC,cens_A,gamma1,gamma2,theta2,gamma,alpha,&
                                    lambda_S,nu_S,lambda_T,nu_T,rsqrt,sigma_s,sigma_t
      double precision, dimension(ver), intent(in)::betas,betat ! vecteur des coefficients associes aux effets fixe du model. contiont le meme nombre d'element, et donc doit etre bien rempli
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
      real , dimension(2):: v1
      real :: piece,demi
      double precision :: ui,temps1,temps1_S,gapx,gapdc,moy_idnum,x,xdc,cens,cbeta1,cbeta2,cbeta3,&
      auxbeta1,auxbeta2,uniran,moyui
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
         cbeta1=betas(1)         
         cbeta3=betat(1)     
         nb_recur = 0           ! nb de temps 
         nb_dc = 0              ! nb de dc
         nb_cens = 0            ! nb de cens
         nobs=0
         max_recu= 1
         moyui = 0.d0
        !close(10)
    do ig=1,ng ! sur les groupes
            
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
!c-----------RECURRENT --------------------------------------------
!c---  genere temps recurrents a partir 2 var explic :
            if (lognormal==0)then ! Gamma
                !print*,"generation gamma avec 2 effets aleatoires correles au niveau essai non encore", &
                !      "implementee suis au probleme de la loi gamma multivariee"
            else ! lognormale
                !ui represente les w_ij dans cette expression
                if (lognormal==1)then !joint surrogate
                    if(frailt_base==1) then ! on tient compte des u_i
                        auxbeta1=ui+don_simul(ig,u_i1)+don_simul(ig,v_s1)*dble(v1(1))+cbeta1*dble(v1(1))!+cbeta2*dble(v1(2)) ! scl je considere uniquement le traitement
                        auxbeta2=truealpha*ui+alpha*don_simul(ig,u_i1)+don_simul(ig,v_t1)*dble(v1(1))+cbeta3*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
                    else ! on ne tient pas compte des u_i dans la generation des temps de survie
                        auxbeta1=ui+don_simul(ig,v_s1)*dble(v1(1))+cbeta1*dble(v1(1))!+cbeta2*dble(v1(2)) ! scl je considere uniquement le traitement
                        auxbeta2=truealpha*ui+don_simul(ig,v_t1)*dble(v1(1))+cbeta3*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
                    endif
                else !(2)joint classique, 2007
                    ! on ne tient pas compte des u_i dans la generation des temps de survie
                    auxbeta1=ui+cbeta1*dble(v1(1))!+cbeta2*dble(v1(2)) ! scl je considere uniquement le traitement
                    auxbeta2=truealpha*ui+cbeta3*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
                endif
            endif
            
            call weigui2(bw1(1),bw1(2),auxbeta1,gapx)
                x=gapx
                         
        !c-----------DECES --------------------------------------------
        !c---     genere temps de dc a partir d une var explic :
                call weigui2(bw2(1),bw2(2),auxbeta2,gapdc) 
        !c**************  gap time: 
                xdc=gapdc
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
                        if(deltadc == 0.d0) then ! si le patient est vivant alors pas de progression
                            delta=0.d0           
                            temps1_S=temps1
                        else ! le patient fait la progression le meme jour que le deces
                            if(pfs == 0) then ! le deces censure la progression et donc on considere qu'il n'ya pas eu de progression
                                delta=0.d0             ! on suppose pas d'evenement si le meme jour que le deces
                                temps1_S=temps1! et on censure a la date de deces(ou censure)
                            else ! dans ce cas la progression inclue le deces: cas de la PFS ou DFS
                                delta=1.d0
                                temps1_S = temps1
                                nb_recur =nb_recur + 1
                                nig(ig) = nig(ig)+1 !nb events recurrents
                            endif
                        endif
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
                       do ii = 1,ver
                          if(filtre(ii).eq.1)then
                             iii = iii + 1
                             ve(nobs,iii) = dble(v1(ii))
                          endif
                          if(filtre2(ii).eq.1)then
                             iii2 = iii2 + 1
                             ve2(nobs,iii2) = dble(v1(ii))
                          endif
                        enddo
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

                !c****************************************************
                !c******** FIN generation des donnees****************
                !c****************************************************
          
                if (maxtemps.lt.t1(nobs))then
                   maxtemps = t1(nobs)
                endif  

                if (delta.eq.0.d0)then ! deces ou censure
                   goto 30          ! on change de sujet
                endif                   
            10      continue          !observations par sujet
           30   continue                  !sujets=groupes    
            enddo        
            
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
    deallocate(tempsD,sigma,mu,x_)
 
endsubroutine Generation_surrogate

!==============================================================================================================================
!     subroutine pour la generation  des donnees du modèle surrogate complet avec un effet aleatoire partage au niveau individuel 
!    et 2 effets correles au niveau essai en interaction avec le traitement
!==============================================================================================================================

subroutine Generation_surrogate_copula(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
            ng,ver,truealpha,propC,cens_A,gamma1,gamma2,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,&
            betat,n_essai,rsqrt,sigma_s,sigma_t,p,prop_i,gamma,alpha,frailt_base,thetacopule,filtre, filtre2,&
            pfs)
    ! lognormal: dit si la distribution des effets aleatoires est lognormal pour le modele complet (1) ou lognormal pour le joint classique de 2007 (2) ou gamma pour le joint classique de 2007(0)
    ! use Autres_fonctions
    ! theta2: variance des frailties gaussiens associe a S
    ! gamma1,gamma2: parametres de la gamma
    ! cens_A:censure administrative. 
    ! propC:proportion minimale de personnes censurees aleatoirement. si 0, alors censure fixe
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
    ! filtre: vecteur qui dit si une variable est prise en compte pour le surrogate
    ! filtre2: vecteur qui dit si une variable est prise en compte pour le true endpoint
    ! pfs : used to specified if the time to progression should be censored by the death time (0) or not (1). The default is 0 as in sofeu et al. (2019). In this case, death is not included in the surrogate endpoint. 
    
     use var_surrogate, only: random_generator

      integer, intent(in)::n_essai,frailt_base,affiche_stat,n_obs,n_col,lognormal,ng,ver, pfs
      double precision, intent(in)::truealpha,propC,cens_A,gamma1,gamma2,theta2,gamma,alpha,&
                                    lambda_S,nu_S,lambda_T,nu_T,rsqrt,sigma_s,sigma_t,&
                                    thetacopule
      double precision, dimension(ver), intent(in)::betas,betat ! vecteur des coefficients associes aux effets fixe du model. contiont le meme nombre d'element, et donc doit etre bien rempli
      double precision,dimension(n_essai),intent(in)::prop_i,p      
      double precision,intent(out)::vrai_theta
      double precision,dimension(n_obs,n_col),intent(out)::don_simulS1,don_simul
      
      integer, parameter::npmax=70,NOBSMAX=15000,nvarmax=45,ngmax=5000,nboumax=1000,NSIMAX=5000,&
                          ndatemax=30000
      integer  j,k,n,ii,iii,iii2,i,nbou
      integer, dimension(ver), intent(in):: filtre, filtre2
      integer ind(nboumax), values(8)
      double precision maxtemps
      character*18 nomvar(nvarmax),nomvar2(nvarmax)
      character*20 dateamj, zone, heure1

!************declarations pour donnees generees **********
      integer :: nb_recur,nb_dc,nb_cens,delta,deltadc,jj,ig,nrecurr,nobs,max_recu, nobs_save, nobs_temp
      real , dimension(:), allocatable:: v1
      real :: piece,demi
      double precision, dimension(n_obs):: x, xdc ! permet le stockage pour utilisation par la suite
      double precision :: ui,temps1,temps1_S,gapx,gapdc,moy_idnum,cens,cbeta2,&
      auxbeta1,auxbeta2,uniran,moyui, cbeta4 
      double precision, dimension(n_essai)::qi
      double precision, dimension(2):: bg1,bw1,bw2
      integer, dimension(ngmax):: idnum
      double precision, dimension(ngmax):: vecui

      integer :: nmax, ibou
      common /nmax/nmax
      double precision date(ndatemax), zi(-2:npmax)
      common /dace1/date,zi
      double precision t0(NOBSMAX),t1(NOBSMAX),t1_S(NOBSMAX), ve(n_obs,ver),ve2(n_obs,ver)
      integer c(NOBSMAX), cdc(NOBSMAX), g(NOBSMAX), nig(ngmax)

      ! Ajout SCL
      double precision::x22,sigma_st,u_i,tempon
      double precision,dimension(:),allocatable::tempsD,mu
      integer::Aretenir
      integer,parameter ::trt1=1,v_s1=2,v_t1=3,trialref1=4,timeS1=5,timeT1=6,&
                      timeC1=7,statusS1=8,statusT1=9,initTime1=10,Patienref1=11,u_i1=12 ! definissent les indices du tableau de donnee simulees
      double precision,dimension(n_essai)::n_i
      double precision,dimension(:,:),allocatable::sigma,x_ 
    
     ! ! some print
      ! call intpr("n_obs", -1,n_obs , 1)    
      ! call intpr("n_col", -1,n_col , 1)
      ! call intpr("lognormal", -1, lognormal, 1)    
      ! call dblepr("vrai_theta", -1,vrai_theta , 1)
      ! call intpr("ng", -1,ng , 1)    
      ! call intpr("ver", -1,ver , 1)
      ! call dblepr("truealpha", -1, truealpha, 1)    
      ! call dblepr("propC", -1,propC , 1)
      ! call dblepr("cens_A", -1, cens_A, 1)    
      ! call dblepr("gamma1", -1,gamma1 , 1)
      ! call dblepr("gamma2", -1,gamma2 , 1)    
      ! call dblepr("theta2", -1, theta2, 1)
      ! call dblepr("lambda_S", -1, lambda_S, 1)    
      ! call dblepr("nu_S", -1,nu_S , 1)
      ! call dblepr("lambda_T", -1, lambda_T, 1)    
      ! call dblepr("nu_T", -1,nu_T , 1)
      ! call dblepr("betas", -1,betas , ver)
      ! call dblepr("betat", -1,betat, ver)    
      ! call intpr("n_essai", -1, n_essai, 1)
      ! call dblepr("rsqrt", -1, rsqrt, 1)
      ! call dblepr("sigma_s", -1,sigma_s , 1)    
      ! call dblepr("sigma_t", -1,sigma_t , 1)
      ! call dblepr("p", -1,p, size(p))
      ! call dblepr("prop_i", -1, prop_i, size(prop_i))    
      ! call dblepr("gamma", -1,gamma , 1)
      ! call dblepr("alpha", -1,alpha , 1)
      ! call intpr("frailt_base", -1, frailt_base, 1)    
      ! call dblepr("thetacopule", -1,thetacopule , 1)
      ! call intpr("filtre", -1,filtre,size(filtre))    
      ! call intpr("filtre2", -1,filtre, size(filtre2))
      ! call intpr("frailt_base", -1, pfs, 1)
!CCCCCCCCCCCCCCCCChosur9.f CCCCCCCCCCCCCCCCCCCCCCCC
      allocate(v1(ver))
      don_simulS1 = 0.d0
      don_simul = 0.d0
      
      call date_and_time(dateamj,heure1,zone,values)                  
      nrecurr=1
      nbou=1
      
      allocate(tempsD(ng))

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
    !call dblepr("sigma =", -1, sigma, size(sigma)**2)
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
        k=k+NINT(n_i(i))
    enddo
    
    don_simulS1=don_simul
        ibou=1 
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
         nb_recur = 0           ! nb de temps 
         nb_dc = 0              ! nb de dc
         nb_cens = 0            ! nb de cens
         nobs=0
         max_recu= 1
         moyui = 0.d0
         x = 0.d0
         xdc = 0.d0
        !close(10)
    do ig=1,ng ! sur les groupes                 
!!c---  variables explicatives par sujet
        do 111 j=1,ver
        if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
            tempon= uniran()
        else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
            CALL RANDOM_NUMBER(tempon)
        endif
            
        piece = real(tempon)
        if (piece.le.demi) then
            v1(j) = 0.
        else
            v1(j) = 1.
        endif
        111        continue
         
        cens = 0.d0
        k = 1
        if(k.gt.max_recu)then
            max_recu = k 
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
                    if(ver > 1) then
                        do ii = 2,ver ! on considere a partir de la 2 ieme variable car le traitement est prise en compte deja
                            if(filtre(ii).eq.1)then
                                auxbeta1 = auxbeta1 + betas(ii)*dble(v1(ii))
                            endif
                            if(filtre2(ii).eq.1)then
                                auxbeta2 = auxbeta2 + betat(ii)*dble(v1(ii))
                            endif
                        enddo
                    endif
                else ! on ne tient pas compte des u_i dans la generation des temps de survie
                    auxbeta1=don_simul(ig,v_s1)*dble(v1(1))+betas(1)*dble(v1(1))
                    auxbeta2=don_simul(ig,v_t1)*dble(v1(1))+betat(1)*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
                    if(ver > 1) then
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
        endif

        call weiguicopule(bw1(1),bw2(1),bw1(2),bw2(2),auxbeta1,auxbeta2,thetacopule,gapx,gapdc)
        x(ig)=gapx ! temps de progression
        xdc(ig)=gapdc ! temps de deces
        tempsD(ig)=xdc(ig)
            
        !c****** for gap time :  
        iii = 0
        iii2 = 0        
        do ii = 1,ver
            if(filtre(ii).eq.1)then
                iii = iii + 1
                ve(nobs,iii) = dble(v1(ii))
            endif
            if(filtre2(ii).eq.1)then
                iii2 = iii2 + 1
                ve2(nobs,iii2) = dble(v1(ii))   
            endif
        enddo                     
    enddo
    
    nobs = 0
    ! recherche des quantiles par essai sur lequel on s'appuie pour la generation uniforme des temps de censures. Exp: 75ieme percentile = 0.75 dans propC
    if(propC > 0.d0) then! censure aleatoire 
        k = 1
        do i = 1, n_essai
            call percentile_scl(xdc(k:(k+NINT(n_i(i))-1)),NINT(n_i(i)),propC,qi(i))
            k = k + NINT(n_i(i))
        enddo
    endif
    
    !call dblepr("temps de deces pour les deux premier essais: xdc(1:40)", -1, xdc(1:40), 40)
    !call intpr("NINT(n_i(1))", -1, NINT(n_i(1)), 1)
    !call dblepr("Voila les quantiles qi", -1, qi, n_essai)
    
    do ig=1,ng ! sur les groupes
        
        nobs = nobs + 1   ! indice l ensemble des observations
        ! scl============censure====================
        if(propC == 0.d0) then! censure fixe ou administrative
            cens = cens_A
        else ! censure aleatoire: generation uniforme entre 1 et la quantile  de l'essai i calcule precedemment
            ! sortie d'etude
            call runif(1.d0, qi(NINT(don_simul(ig,trialref1))), cens)
            ! je prends le min entre la censure administrative et la date de sortie d'etude.
            cens = min(cens, cens_A)
        endif
        
        if(xdc(ig).le.cens)then ! patient decede
            deltadc=1.d0
            temps1 = xdc(ig)
            nb_dc =nb_dc + 1
        else    !patient censuree administrativement
            deltadc=0.d0
            temps1 = cens
            nb_cens =nb_cens + 1
        endif

                !on construit les temps de progression
        if(x(ig) < temps1)then ! evenement avant la censure
            delta=1.d0
            temps1_S = x(ig)
            nb_recur =nb_recur + 1
            nig(ig) = nig(ig)+1 !nb events recurrents
        else
            if((x(ig).eq.cens).and.(deltadc==0.d0)) then !evenement a la date de censure et patient vivant
                delta=1.d0
                temps1_S = x(ig)
                nb_recur =nb_recur + 1
                nig(ig) = nig(ig)+1 !nb events recurrents
            else ! progression le meme jour que le deces ou sans progression
                 ! delta=0.d0           
                ! temps1_S=temps1
                if(deltadc == 0.d0) then ! si le patient est vivant alors pas de progression
                    delta=0.d0           
                    temps1_S=temps1
                else ! le patient fait la progression le meme jour que le deces
                    if(pfs == 0) then ! le deces censure la progression et donc on considere qu'il n'ya pas eu de progression
                        delta=0.d0             ! on suppose pas d'evenement si le meme jour que le deces
                        temps1_S=temps1! et on censure a la date de deces(ou censure)
                    else ! dans ce cas la progression inclue le deces: cas de la PFS ou DFS
                        delta=1.d0
                        temps1_S = temps1
                        nb_recur =nb_recur + 1
                        nig(ig) = nig(ig)+1 !nb events recurrents
                    endif
                endif
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

        !c*** pour le tester sur un autre programme: on complete les nouveaux parametres simules dans le jeux de donnees
        don_simulS1(ig,initTime1)=t0(nobs)
        don_simulS1(ig,timeS1)=t1_S(nobs)
        don_simulS1(ig,statusS1)=c(nobs)
        don_simulS1(ig,Patienref1)=g(nobs)
        don_simulS1(ig,trt1)=ve2(nobs,1)
        !don_simulS1(ig,w_ij1)=ui        
                   
        don_simul(ig,initTime1)=t0(nobs)
        don_simul(ig,timeT1)=t1(nobs)
        don_simul(ig,statusT1)=cdc(nobs)
        don_simul(ig,Patienref1)=g(nobs)
        don_simul(ig,trt1)=ve2(nobs,1)
        !don_simul(ig,w_ij1)=ui                   
                
        ! j'ajoute les autres variables a la fin
        if(ver > 1) then
            do ii = 2,ver
                if(filtre(ii).eq.1)then
                    don_simulS1(ig,size(don_simulS1,2) - ver + ii - 1)=ve(nobs,ii)
                endif
                if(filtre2(ii).eq.1)then
                    don_simul(ig,size(don_simul,2)- ver + ii - 1)=ve2(nobs,ii)
                endif
            enddo    
        endif                    

        !c****************************************************
        !c******** FIN generation des donnees****************
        !c****************************************************
          
        if (maxtemps.lt.t1(nobs))then
            maxtemps = t1(nobs)
        endif                    
    enddo                  !observations par sujet            
     ! scl============censure conseillee pour la proportion souhaitee de censure==
    call percentile_scl(tempsD,ng,1.d0-propC,cens)
    deallocate(tempsD,sigma,mu,x_,v1)
    ! call dblepr("betas", -1, betas, ver)
    ! call dblepr("betat", -1, betat, ver)
    ! call dblepr("voile don_simul", -1, don_simul(1,:), size(don_simul,2))
    ! call dblepr("voile don_simul", -1, don_simul(2,:), size(don_simul,2))
endsubroutine Generation_surrogate_copula

!==============================================================================================================================
!     subroutine pour la generation  des donnees du modèle surrogate complet avec des effets aleatoires correles aussi bien au niveau individuel 
!    qu'au niveau essai avec et sans interaction avec le traitement
!==============================================================================================================================

! subroutine Generation_surrogate_complet(don_simul,don_simulS1,n_obs,n_col,lognormal,affiche_stat,vrai_theta,&
            ! ng,ver,truealpha,propC,cens_A,gamma1,theta2,lambda_S,nu_S,lambda_T,nu_T,betas,&
            ! betat,n_essai,rsqrt,sigma_s,sigma_t,prop_i,gamma,theta2_t,rsqrt_theta,&
            ! gamma_t,rsqrt_gamma,use_gamma_st)
    ! ! lognormal: dit si la distribution des effets aleatoires est lognormal (1) ou gamma (0)
    ! ! use Autres_fonctions
    ! ! theta2: variance des frailties gaussiens associe a S
    ! ! gamma1,gamma2: parametres de la gamma
    ! ! cens_A:censure administrative
    ! ! propC:proportion de personnes censurees
    ! ! n_essai:  nombre d'essai a generer
    ! ! rsqrt: niveau de correlation souhaite entre les fragilites sepecifiques aux traitement au niveau essai S et T
    ! ! sigma_s: variance des frailties au niveau essai associee au surrogate
    ! ! sigma_t: variance des frailties au niveau essai associee au true
    ! ! p: proportion des personnes traitees par essai
    ! ! prop_i: proportion des sujets par essai
    ! ! gamma: variance de l'effet aleatoire u_i associe au risque de base chez S
    ! ! truealpha: fragilile associe a w_ij pour les deces
    ! ! alpha: fragilile associe a u_i pour les deces
    ! ! frailt_base: dit si l'on prend en compte l'heterogeneite sur le risque de base aussi bien dans la generation des donnes que dans l'estimation(1) ou non (0)
    ! ! theta2_t: variance des frailties gaussiens associe a T
    ! ! rsqrt_theta: niveau de correlation entre wij_s et wij_t
    ! ! gamma_t: variance de l'effet aleatoire u_i associe au risque de base chez T
    ! ! rsqrt_gamma: niveau de correlation entreus_i et ut_i
    ! ! use_gamma_st: dit si l'on va estimer(1) ou pas(0) la correlation entre us_i et ut_i 
      ! use var_surrogate, only: random_generator
      
      ! integer, intent(in)::n_essai,use_gamma_st!,frailt_base
      ! double precision,dimension(n_essai),intent(in)::prop_i!,p
      ! integer, intent(in)::n_obs,n_col,lognormal,ng,ver
      ! double precision, intent(in)::truealpha,propC,cens_A,gamma1,theta2,gamma,alpha,theta2_t,gamma_t,rsqrt_gamma!,gamma2
      ! double precision, intent(in)::lambda_S,nu_S,lambda_T,nu_T,betas,betat,rsqrt,sigma_s,sigma_t,rsqrt_theta
      ! double precision,dimension(n_obs,n_col),intent(out)::don_simulS1,don_simul
      ! integer, parameter::npmax=70,NOBSMAX=15000,nvarmax=45,ngmax=5000
      ! integer,parameter::nboumax=1000,NSIMAX=5000,ndatemax=30000
      
      ! integer  j,k,nz,cpt,cpt_dc,ii,iii,iii2!,kk,ij,groupe,n,np
      ! integer  cptstr1,cptstr2,trace,trace1,trace2
      ! integer  i!,ni,ef,ic,ic2,ier,istop
      ! integer  cptni,cptni1,cptni2,nb_echec,nb_echecor
      ! integer  nbou2,cptbiais!,id,l
      ! integer  nbou!,m,icen,idum
      ! integer filtre(nvarmax), filtre2(nvarmax) 
      ! ! integer cpt1(nboumax) 
      ! ! integer cpt2(nboumax) 
      ! ! integer cpt3(nboumax) 
      ! integer ind(nboumax) , cptaux 
      
      ! real vax(nvarmax)
      ! double precision tt0
      ! double precision tt1,tt2
      
      ! ! double precision h
      ! double precision ro,wres!,csi,csi1
      ! double precision res,maxtemps !ax2,ax1,max,min
      ! double precision wald,str !, bs,bi
      ! ! double precision moyvar1!theta,moyse1,moyse_cor1
      ! ! double precision moyvar2!,moysebeta3,moysebeta1,alpha,moyse2,moyse_cor2
      ! ! double precision moysebeta_cor1!,moysebeta2,moybeta3, moybeta2,moybeta1 
      ! ! double precision moysebeta_cor2
      ! ! double precision moysebeta_cor3
      ! ! double precision moybweib1,moybweib2,moybweib3,moybweib4!weibull
      ! ! double precision moysebweib1,moysebweib2,moysebweib3,moysebweib4
! !c     double precision  tt0,tt1
! !c     double precision, pointer  tt0
! !c     double precision,pointer  tt1
      ! double precision varxij,varsmarg,smoy,smoyxij!,eca
      ! double precision lrs!,pe2,pe1,
      ! double precision BIAIS_moy
      
      ! ! double precision aux(2*NOBSMAX)
      ! double precision v((npmax*(npmax+3)/2))
      ! !double precision k0(2)
      ! ! double precision b(npmax)
      ! double precision se1(nboumax),se_cor1(nboumax)!theta
      ! double precision se2(nboumax),se_cor2(nboumax)!alpha
      ! double precision tvars1(nboumax),tvars2(nboumax)
      ! ! double precision beta1(nboumax),beta2(nboumax),beta3(nboumax)
      ! double precision , dimension(nboumax):: sebeta1,sebeta_cor1
      ! double precision , dimension(nboumax):: sebeta2,sebeta_cor2
      ! double precision , dimension(nboumax):: sebeta3,sebeta_cor3

      ! ! double precision bweib1(nboumax),bweib2(nboumax)!weibull
      ! ! double precision bweib3(nboumax),bweib4(nboumax)
      ! double precision , dimension(nboumax):: sebweib1,sebweib2
      ! double precision , dimension(nboumax):: sebweib3,sebweib4

      ! ! double precision biais_theta(nboumax)
      ! ! double precision I1_hess(npmax,npmax)!,H1_hess(npmax,npmax)
      ! !double precision I2_hess(npmax,npmax)!,H2_hess(npmax,npmax)
      ! ! double precision HI1,HI2(npmax,npmax)
      ! !double precision IH(npmax,npmax)!,HIH(npmax,npmax),HI(npmax,npmax)
      ! !double precision BIAIS(npmax,1)
      
      ! character*18 nomvarl
      ! character*18 nomvar(nvarmax),nomvar2(nvarmax)
      ! ! character*18 donnees
      ! character*24 ficpar
      ! !character*14 fich1
      ! ! character*14 fich2
      ! ! character*14 fich3
      ! ! character*14 fich4
      ! ! character*14 fich1b
      ! ! character*14 fich2b
      ! ! character*14 fich3b
      ! ! character*14 fich4b
      ! character*20 dateamj
      ! character*20 zone
      ! character*20 heure1
      ! ! character*20 heure2
      ! integer values(8)

! !c************declarations pour donnees generees **********
      ! integer :: nb_recur,nb_dc,nb_cens,delta,deltadc,jj!,no
      ! integer :: ig,sg,nrecurr,nobs,max_recu
      ! real , dimension(2):: v1
      ! real :: piece,rien,demi
      ! double precision :: ui,temps1,temps1_S !random effect
      ! double precision :: gapx,gapdc,moy_idnum!,gapcens
      ! double precision :: x,xdc,tronc,cens,cbeta1,cbeta3!,cbeta2
      ! double precision :: auxbeta1,auxbeta2 ! for recurr and death
      ! double precision :: zbqlexp,uniran
      ! double precision, dimension(2):: bg1,bw1,bw2
      ! integer, dimension(ngmax):: idnum
      ! double precision, dimension(ngmax):: vecui
      ! double precision :: moyui
! !c*****************************************************************
      
! !c*****************************************************************
! !c***** nmax
         ! integer :: nmax
         ! common /nmax/nmax
! !c*****dace1 
      ! double precision date(ndatemax)
      ! double precision zi(-2:npmax)
      ! common /dace1/date,zi
      
! !c*****dace2
      ! double precision t0(NOBSMAX),t1(NOBSMAX),t1_S(NOBSMAX)
      ! integer c(NOBSMAX), cdc(NOBSMAX)
      ! ! integer nt0(NOBSMAX),nt1(NOBSMAX)
      ! integer  nva1,nva2,nst!,nva,nobs,ndate
      ! !common /dace2/t0,t1,c,cdc,nt0,nt1,nobs,nva,nva1,nva2,ndate,nst
! !c*****dace4
      ! integer  stra(NOBSMAX)
      ! !common /dace4/stra
! !c*****ve1
      ! double precision ve(NOBSMAX,nvarmax),ve2(NOBSMAX,nvarmax)
      ! !common /ve1/ve,ve2
! !c*****dace3
      ! ! double precision  pe
      ! integer  effet,nz1,nz2
      ! !common /dace3/pe,effet,nz1,nz2
! !c*****dace7
      ! !double precision I_hess(npmax,npmax)!,H_hess(npmax,npmax)
      ! ! double precision Hspl_hess(npmax,npmax)
      ! ! double precision PEN_deri(npmax,1)
      ! ! double precision hess(npmax,npmax)
      ! !common /dace7/PEN_deri,I_hess,H_hess,Hspl_hess,hess
! !c*****contrib
      ! !common /contrib/ng       
! !c*****groupe
      ! integer g(NOBSMAX)
      ! integer nig(ngmax)        ! nb d events recurrents , different de mi()!
      ! !common /gpe/g,nig
      
! !c*****mem1
      ! ! double precision mm3(ndatemax),mm2(ndatemax)
      ! ! double precision mm1(ndatemax),mm(ndatemax)
      ! !common /mem1/mm3,mm2,mm1,mm
! !c     %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
      ! ! integer AG
      ! !common /andersengill/AG
! !c     %%%%%%%%%%%%% indic ALPHA %%%%%%%%%%%%%%%%%%%%%%%%% 
      ! integer indic_ALPHA
      ! !common /alpha/indic_ALPHA ! pour preciser un para en plus 
! !c**** theta/alpha
      ! double precision  theta !en exposant pour la frailty deces 
      ! !common /thetaalpha/alpha,theta
! !c******indicateur de troncature
      ! !integer :: indictronq     ! =0 si donnees non tronquées reellement
      ! !common /troncature/indictronq
      
! !c******indicateur iteration
      ! integer  ibou
      ! !common /boucle/ibou
! !c************FIN COMMON ***********************************
      
! !c     nst: deux finctions de risque a estimer (meme bases de splines)
! !c     ist: appartenance aux strates
! !c icen=1: censure
! !c     ib: matrice de variance de beta chapeau
! !c     I_hess : -hessienne non inversee sur vraisemblance non penalisee
! !c     H_hess : inverse de -hessienne  sur vraisemblance penalisee
      ! !Ajout SCL
      ! double precision::x22,sigma_st,u_i,gamma_st,theta_st,tempon !, cens_C
      ! double precision,dimension(:),allocatable::tempsD
      ! integer::Aretenir
      ! ! definissent les indices du tableau de donnee simulees
      ! integer,parameter ::trt1=1,v_s1=2,v_t1=3,trialref1=4,w_ij1=5,timeS1=6,timeT1=7,&
                          ! timeC1=8,statusS1=9,statusT1=10,initTime1=11,Patienref1=12,u_i1=13,& 
                          ! w_ijt=14,u_it=15
      ! integer,intent(in)::affiche_stat
      ! double precision,intent(out)::vrai_theta
      ! double precision,dimension(n_essai)::n_i
      ! double precision,dimension(:,:),allocatable::sigma,sigma_wij
      ! double precision,dimension(:),allocatable::mu
      ! double precision,dimension(:,:),allocatable::x_      
      
! !CCCCCCCCCCCCCCCCChosur9.f CCCCCCCCCCCCCCCCCCCCCCCC
    ! if(affiche_stat==1)then
      ! !write(*,*)'    ******************************************'
      ! !write(*,*)'  ****** DEBUT PROGRAMME FRAILTY.F**********'
      ! !write(*,*)'******************************************'
    ! endif
      ! nmax = 300 !nb iterations max dans marquard
      
      ! indic_alpha=1 ! on precise que l on a un parametre en plus estimer
! !c      allocate (tt0)
! !c      allocate (tt1)
      
      ! call date_and_time(dateamj,heure1,zone,values)                  
! !c     !write(*,*)'Starting time: ', dateamj,heure1,zone,values
      
      ! lrs=0.d0      
      ! nb_echec=0
      ! nb_echecor=0
      ! nbou2=0
      ! ficpar='joint2.inf'
      ! !open(4,file='outjoint')
      
      ! nrecurr=1
      ! nbou=1
    ! if(affiche_stat==1)then
      ! !write(4,*)'** nb de simulations',nbou
      ! !write(*,*)'** nb de simulations',nbou
    ! endif
      
      ! allocate(tempsD(ng))
      ! nst=2
    
    ! if(affiche_stat==1)then
      ! !write(4,*)'**************************************************'
      ! !write(4,*)'************ JOINT MODEL *************************'
      ! !write(4,*)'*** RECURRENT EVENTS and TERMINATING EVENT *******'
      ! !write(4,*)'**************************************************'
      ! !write(4,*)'** nb de groupes=sujets =',ng 
      ! !write(*,*)'** nb de groupes=sujets =',ng 
      ! !write(4,*)'** deux fonctions de risque de base  = ',nst
      ! !write(*,*)'** fonction de risque de base  = ',nst
    ! endif

! !c     !write(4,*)'** nb de simulations = ',nbou
! !c     !write(4,*)'** indicateur de censure (1 censure existante)',icen
        ! cptni=0
        ! cptni1=0
        ! cptni2=0
        ! biais_moy=0.d0
        ! cptbiais=0
        ! cptaux=0

! !c**************************************************
! !c********************* prog spline****************
         ! !read(2,*) !scl pour faire passer a la ligne suivante
         ! !read(2,*)ver

         ! nva1 = 0 ! nb de var expli pour donnees recurrentes
         ! nva2 = 0 ! nb de var expli pour deces
         ! !!print*,"ver=",ver

         ! if(ver.gt.0)then
            ! do 44 j=1,ver
               ! !read(2,*)nomvarl,filtre(j) ,filtre2(j)
               ! nomvarl="trt"
               ! filtre(j)=1
               ! filtre2(j)=1
               ! nva1 = nva1 + filtre(j) ! adjustment for recurrent events
               ! nva2 = nva2 + filtre2(j) ! adjustment for survival
               ! if(filtre(j).eq.1)then
                  ! nomvar(nva1) = nomvarl
               ! endif 
               ! if(filtre2(j).eq.1)then
                  ! nomvar2(nva2) = nomvarl
               ! endif    
 ! 44         continue
         ! endif

 ! !        nva = nva1+nva2
    ! if(affiche_stat==1)then
      ! !write(4,*)'** explanatory variables for recurrent events:',nva1
      ! !write(*,*)'** explanatory variables for recurrent events:',nva1
      ! !write(4,*)'** explanatory variables for deaths:',nva2
      ! !write(*,*)'** explanatory variables for deaths:',nva2
    ! endif
    
    ! !  read(2,*)truealpha 
      
    ! if(affiche_stat==1)then
      ! !write(4,*)'** Vraie valeur de Alpha',truealpha 
      ! !write(*,*)'** Vraie valeur de Alpha',truealpha 
    ! endif

! !c %%%%%%%%%%%%% ANDERSEN-GILL %%%%%%%%%%%%%%%%%%%%%%%%% 
! !      read(2,*)AG
! !c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! !      read(2,*)nz
        ! nz=6
      
      ! nz1=nz
      ! nz2=nz
      ! if(nz.gt.20)then
         ! nz = 20
      ! endif
      ! if(nz.lt.4)then
         ! nz = 4
      ! endif
      
    ! if(affiche_stat==1)then
     ! !write(4,*)'nombre de noeuds :',nz
    ! endif
    
    ! Aretenir=1
        
    ! if(Aretenir>nbou) then
        ! !print*,"l'indice du jeu de données à retenitr doit etre inferieure au nombre de simulation"
        ! stop
    ! endif


! !c*******************************************
! !c***** DEBUT SIMULATIONS *******************
! !c*******************************************
    ! maxtemps = 0.d0
    ! cpt = 0
    ! cpt_dc = 0
    ! cptstr1 = 0
    ! cptstr2 = 0
    ! don_simulS1=0.d0
    ! don_simul=0.d0
    
! !==============initialisation des parametres======================
    ! n_i=NINT(n_obs*prop_i) ! nombre de sujet par essai

    ! if(sum(n_i)<n_obs) then
        ! n_i(minloc(n_i,mask=n_i==minval(n_i)))=n_i(minloc(n_i,mask=n_i==minval(n_i)))+(n_obs-sum(n_i)) ! on ajoute au premier essai de plus petite taille le nombre d'individu non encore affecte (1 generalement) a cause des problemes d'arrondi
    ! endif
    ! if(sum(n_i)>n_obs) then 
        ! n_i(maxloc(n_i,mask=n_i==maxval(n_i)))=n_i(maxloc(n_i,mask=n_i==maxval(n_i)))-(sum(n_i)-n_obs) ! on soustrait au premier essai de plus grande taille le nombre d'individu affecte en trop (1 generalement) a cause des problemes d'arrondi
    ! endif
    
    ! ! simulation des effets aleatoires specifiques aux essais
    
    ! ! matrice des covariances des frailties au niveau essai, en supposant une correlation de 0.95 (erreurs liees aux effets du traitement sur S et T niveau essai) 
    ! sigma_st=rsqrt*dsqrt(sigma_s)*dsqrt(sigma_t)
    ! gamma_st=rsqrt_gamma*dsqrt(gamma)*dsqrt(gamma_t)
    ! theta_st=rsqrt_theta*dsqrt(theta2)*dsqrt(theta2_t)
    
    ! !=======================================sigma============
    ! allocate(sigma(2,2),mu(2),x_(n_essai,2))
    ! sigma=0.d0 ! initialisation
    ! x_=0.d0
    ! !pour (vs_i,vt_i)
    ! sigma(1,1)=sigma_s
    ! sigma(1,2)=sigma_st
    ! sigma(2,1)=sigma_st
    ! sigma(2,2)=sigma_t
    
   ! ! if(affiche_stat==1) !print*,"Matrice de var-cov frailties essai sigma:",sigma
    
    ! call rmvnorm(mu,sigma,n_essai,0,x_)    
    ! ! !print*,sum(x_(1:n_essai,1))/n_essai,sum(x_(1:n_essai,2))/n_essai
    
     ! ! do i=1,n_essai
        ! ! ! call rmvnorm(mu,sigma,1,0,x_)
         ! ! !print*,i,x_(i,:)
     ! ! enddo
    ! ! stop
    
    
    ! k=1
    ! do i=1,n_essai
        ! ! effet aleatoires specifiques aux essais
        ! don_simul(k:((k+NINT(n_i(i)))-1),v_s1)=x_(i,1)
        ! don_simul(k:((k+NINT(n_i(i)))-1),v_t1)=x_(i,2)
        ! don_simulS1(k:((k+NINT(n_i(i)))-1),v_s1)=x_(i,1)
        ! don_simulS1(k:((k+NINT(n_i(i)))-1),v_t1)=x_(i,2)
        ! k=k+n_i(i)
    ! enddo
    
    ! !======gamma============
    ! sigma=0.d0 ! initialisation
    ! x_=0.d0
    ! !pour (us_i,ut_i)
    ! sigma(1,1)=gamma
    ! if(use_gamma_st==1) then ! juste si on veut tenir compte de la correlation entre (us_i,ut_i)
        ! sigma(1,2)=gamma_st
        ! sigma(2,1)=gamma_st
    ! endif
    ! sigma(2,2)=gamma_t
    
    ! mu=0.d0
    
    ! !if(affiche_stat==1)    !print*,"Matrice de var-cov frailties essai gamma:",sigma
    
    ! !generation de (us_i, ut_i) suivant une multinormale
    ! call rmvnorm(mu,sigma,n_essai,0,x_)    
      ! ! do i=1,10
        ! ! ! call rmvnorm(mu,sigma,1,0,x_)
         ! ! !print*,i,x_(i,:)
     ! ! enddo
    
    ! !===========================================
        
    ! ! allocate(sigma(4,4),mu(4),x_(n_essai,4))
    ! ! sigma=0.d0 ! initialisation
    ! ! !pour (vs_i,vt_i)
    ! ! sigma(1,1)=sigma_s
    ! ! sigma(1,2)=sigma_st
    ! ! sigma(2,1)=sigma_st
    ! ! sigma(2,2)=sigma_t
    
    ! ! if(affiche_stat==1) !print*,"Matrice de var-cov frailties essai sigma:",sigma
    
    ! ! !pour (us_i,ut_i)
    ! ! sigma(3,3)=gamma
    ! ! if(use_gamma_st==1) then ! juste si on veut tenir compte de la correlation entre (us_i,ut_i)
        ! ! sigma(3,4)=gamma_st
        ! ! sigma(4,3)=gamma_st
    ! ! endif
    ! ! sigma(4,4)=gamma_t
    
    ! ! mu=0.d0
    
    ! ! if(affiche_stat==1)    !print*,"Matrice de var-cov frailties essai sigma+gamma:",sigma
    
    ! ! !generation de (vs_i, vt_i) suivant une multinormale
    ! ! call rmvnorm(mu,sigma,n_essai,0,x_)    
      ! ! do i=1,10
        ! ! ! call rmvnorm(mu,sigma,1,0,x_)
         ! ! !print*,i,x_(i,:)
     ! ! enddo
    ! ! stop
    
    ! k=1
    ! do i=1,n_essai
        ! ! effet aleatoires  associee aux risques de base
        ! don_simul(k:((k+NINT(n_i(i)))-1),u_i1)=x_(i,1)
        ! don_simul(k:((k+NINT(n_i(i)))-1),u_it)=x_(i,2)
        ! don_simulS1(k:((k+NINT(n_i(i)))-1),u_i1)=x_(i,1)
        ! don_simulS1(k:((k+NINT(n_i(i)))-1),u_it)=x_(i,2)
        
        ! ! variable trialref
        ! don_simul(k:((k+NINT(n_i(i)))-1),trialref1)=i
        ! don_simulS1(k:((k+NINT(n_i(i)))-1),trialref1)=i
        ! !!print*,don_simul(k:((k+NINT(n_i(i)))-1),trialref1)
        ! k=k+n_i(i)
    ! enddo
    
    ! !!print*,"covar=",covariance(don_simul(:,v_s1),don_simul(:,v_s1))
    ! !!print*,"vs=",variance(don_simul(:,v_s1))
    ! !!print*,"vt=",variance(don_simul(:,v_t1))
    ! !!print*,"r=",covariance(don_simul(:,v_s1),don_simul(:,v_s1))/(dsqrt(variance(don_simul(:,v_s1)))*dsqrt(variance(don_simul(:,v_t1))))
    ! !don_simulS1=don_simul
    
    ! do 1000 ibou=1,nbou 
        ! if(affiche_stat==1)then
            ! !write(*,*)'************ Generation NUMERO :      ',ibou
        ! endif
        ! ind(ibou)=0
        ! effet = 1
        ! nig = 0
        ! g=0
        ! maxtemps = 0.d0
        ! auxbeta1=0.d0
        ! auxbeta2=0.d0

! !c********************************************************
! !c******** DEBUT generation des donnees *******************
! !c********************************************************

        ! !open(10,file="parametre_2007.inf")
! !!c------------------ UI = FRAILTY ---------------------
! !c 2 parametres de la gamma tq E(ui)=bg1(1)/bg1(2)  var(ui)=bg1(1)/(bg1(2)*bg1(2))
           ! !read(10,*)bg1(1)
            ! bg1(1)=gamma1
           ! !read(10,*)bg1(2)
            ! bg1(2)=gamma1
           ! !read(10,*)theta2 ! variance des frailties gaussien wij
           ! !!print*,"=====================",bg1(1),bg1(2)
           ! if(ibou.eq.1)then
                ! !write(4,*)'** '
                ! if (lognormal==0)then ! Gamma
                    ! if(affiche_stat==1)then    
                        ! !write(4,*)'** vraie valeur de theta = **'&
                         ! !   ,bg1(1)/(bg1(2)*bg1(2))
                        ! !write(*,*)'** vraie valeur de theta = **'&
                         ! !   ,bg1(1)/(bg1(2)*bg1(2))
                    ! endif
                    ! vrai_theta=bg1(1)/(bg1(2)*bg1(2))
                ! else !lognormale
                    ! if(affiche_stat==1)then    
                        ! !write(4,*)'** vraie valeur de theta = **'&
                        ! !    ,theta2
                        ! !write(*,*)'** vraie valeur de theta = **'&
                          ! !  ,theta2
                    ! endif
                    ! vrai_theta=theta2
                ! endif
           ! endif

! !!c---------------- X --------------------

! !!c-- parametres de la WEibull for recurrent events 
         ! bw1(1)=lambda_S 
         ! bw1(2)=nu_S
! !!c-- parametres de la WEibull for death
         ! bw2(1)=lambda_T
         ! bw2(2)=nu_T

         ! demi = 0.5             ! pour var expli
         ! cbeta1=betas          
         ! cbeta3=betat        
         
         ! if(affiche_stat==1)then
            ! if(ibou.eq.1)then
                ! !write(4,*)'**vraie valeur de beta1 (recurrent) = **',cbeta1
                ! !write(4,*)'**vraie valeur de beta2  (recurrent) = **',cbeta2
                ! !write(4,*)'** vraie valeur de beta3 (deces) = **',cbeta3
                ! !write(*,*)'**vraie valeur de beta1 (recurrent) = **',cbeta1
                ! !write(*,*)'**vraie valeur de beta2  (recurrent) = **',cbeta2
                ! !write(*,*)'** vraie valeur de beta3 (deces) = **',cbeta3
            ! endif
        ! endif
        ! !stop
         ! nb_recur = 0           ! nb de temps 
         ! nb_dc = 0              ! nb de dc
         ! nb_cens = 0            ! nb de cens
         ! nobs=0
         ! max_recu= 1
         ! moyui = 0.d0
        ! !close(10)
! !!c--------------------------------------------------------------
! !!c--------------------------------------------------------------
! !!c--------------------------------------------------------------
    
    ! ! on reinitialise les matrices pour la generation au niveau individuel
     ! deallocate(sigma,mu,x_)
     ! allocate(sigma(2,2),mu(2),x_(1,2))
    ! ! allocate(sigma(2,2),mu(2),x_(ng,2))
    ! sigma=0.d0
    ! mu=0.d0
    
    ! !pour (ws_ij,wt_ij)
    ! sigma(1,1)=theta2
    ! sigma(1,2)=theta_st
    ! sigma(2,1)=theta_st
    ! sigma(2,2)=theta2_t    
    ! ! !print*,sigma
    
    ! ! stop
    
   ! ! if(affiche_stat==1)    !print*,"Matrice de var-cov frailties theta:",sigma
    ! !!print*,"mu",mu
    
    ! !generation de (ws_i, wt_i) suivant une multinormale
    ! ! call rmvnorm(mu,sigma,ng,0,x_)
      ! ! do i=1,10
        ! ! ! call rmvnorm(mu,sigma,1,0,x_)
         ! ! !print*,i,x_(i,:)
     ! ! enddo
    ! ! stop
    ! ! ! ! effet aleatoires specifiques aux individus
    ! ! don_simul(:,w_ij1)=x_(:,1)
    ! ! don_simul(:,w_ijt)=x_(:,2)
    ! ! don_simulS1=don_simul
    
    
    ! do 30 ig=1,ng ! sur les groupes
            
        ! if (lognormal==0)then ! Gamma    
            ! call gamgui(bg1(1),ui) !genere gamma pour zi
            ! ui = ui/bg1(2)
! !c     verification des ui
            ! vecui(ig) = ui
            ! moyui = moyui + ui
        ! else !lognormale
            ! ! x22=0.d0
            ! ! !bg1(1)= variance des wij
            ! ! call bgos(sqrt(theta2),0,ui,x22,0.d0) !on simule les w_ij suivant une normale
            ! ! !verification des ui
            ! ! vecui(ig) = ui
            ! ! moyui = moyui + ui
            
            ! ! simulation des ws_ij et wt_ij
            ! call rmvnorm(mu,sigma,1,0,x_)    
            ! !if(ig .le.10) !print*,ig,x_
    ! ! effet aleatoires specifiques aux individus
            ! don_simul(ig,w_ij1)=x_(1,1)
            ! don_simul(ig,w_ijt)=x_(1,2)
            ! don_simulS1(ig,w_ij1)=x_(1,1)
            ! don_simulS1(ig,w_ijt)=x_(1,2)
        ! endif
        
        
! !!c---  variables explicatives par sujet
            ! do 111 j=1,ver
                ! if(random_generator==2)then ! on generer avec uniran(mais gestion du seed pas garanti)
                    ! tempon= uniran()
                ! else !on generer avec RANDOM_NUMBER(avec gestion du seed garanti)
                    ! CALL RANDOM_NUMBER(tempon)
                ! endif
               ! piece=real(tempon)
               ! !piece=real(uniran()) !rand()
               ! if (piece.le.demi) then
                  ! v1(j) = 0.
               ! else
                  ! v1(j) = 1.
               ! endif
 ! 111        continue
            
               ! x=0.d0
               ! xdc=0.d0
               ! cens=0.d0

               ! do 10 k=1,nrecurr ! observations max / sujet
                  ! if(k.gt.max_recu)then
                     ! max_recu=k 
                  ! endif

                  ! nobs=nobs+1   ! indice l ensemble des observations
                  ! idnum(ig) = k
! !c----------- CENSORING --------------------------------------------
                  ! !cens =  1.d0 + 250.d0*uniran() !scl voir plubas
              
! !c-----------RECURRENT --------------------------------------------
! !c---  genere temps recurrents a partir 2 var explic :
    ! if (lognormal==0)then ! Gamma
        ! !print*,"generation gamma avec 2 effets aleatoires correles au niveau essai non encore", &
        ! !        "implementee suis au probleme de la loi gamma multivariee"
        ! stop
    ! else ! lognormale
        ! ! if(frailt_base==1) then ! on tient en compte les u_i
        ! !ui represente les w_ij dans cette expression
        ! auxbeta1=don_simul(ig,w_ij1)+don_simul(ig,u_i1)+don_simul(ig,v_s1)*dble(v1(1))+cbeta1*dble(v1(1))!+cbeta2*dble(v1(2)) ! scl je considere uniquement le traitement
        ! auxbeta2=don_simul(ig,w_ijt)+don_simul(ig,u_it)+don_simul(ig,v_t1)*dble(v1(1))+cbeta3*dble(v1(1)) ! on utilise le log pour pouvour mettre l'expression dans l'exponentiel
    ! endif

        ! call weigui2(bw1(1),bw1(2),auxbeta1,gapx)
! !c**************gap time: 
        ! x=gapx
                 
! !c-----------DECES --------------------------------------------
! !c---     genere temps de dc a partir d une var explic :
        ! call weigui2(bw2(1),bw2(2),auxbeta2,gapdc) 
! !c************** calendar time: 
! !c               xdc=xdc+gapdc
! !c**************  gap time: 
        ! xdc=gapdc
        ! tempsD(ig)=xdc
                 ! ! !print*,"xdc=",xdc,"x=",x,"cens_A",cens_A
                ! ! stop
                 
         
! ! scl============censure====================
        ! cens=cens_A
! !c------------------------------------------------bilan:
                ! ! if ((xdc.le.x).and.(xdc.le.cens)) then !deces
                ! !    deltadc = 1
                ! !    delta = 0
                ! !    temps1 = xdc
                ! !    nb_dc =nb_dc + 1
                ! ! endif
                
                ! if(xdc.le.cens)then ! patient decede
                    ! deltadc=1.d0
                    ! temps1 = xdc
                    ! nb_dc =nb_dc + 1
                ! else    !patient censuree administrativement
                    ! deltadc=0.d0
                    ! temps1 = cens
                    ! nb_cens =nb_cens + 1
                ! endif


                ! ! if ((cens.le.x).and.(cens.le.xdc)) then ! censure
                ! !    deltadc = 0
                ! !    delta = 0
                 ! !   temps1 = cens
                 ! !   nb_cens =nb_cens + 1
                ! ! endif

                ! ! if ((x.le.cens).and.(x.le.xdc)) then ! recurrent event
                ! !    deltadc = 0
                ! !    delta = 1
                ! !   temps1 = x
                ! !    nb_recur =nb_recur + 1
                ! !    nig(ig) = nig(ig)+1 !nb events recurrents
                ! ! endif  
                 
        
        ! !on construit les temps de progression
        ! if(x < temps1)then ! evenement avant la censure
            ! delta=1.d0
            ! temps1_S = x
            ! nb_recur =nb_recur + 1
            ! nig(ig) = nig(ig)+1 !nb events recurrents
        ! else
            ! if((x.eq.cens).and.(deltadc==0.d0)) then !evenement a la date de censure et patient vivant
                ! delta=1.d0
                ! temps1_S = x
                ! nb_recur =nb_recur + 1
                ! nig(ig) = nig(ig)+1 !nb events recurrents
            ! else ! progression le meme jour que le deces ou sans progression
                ! delta=0.d0             ! on suppose pas d'evenement si le meme jour que le deces
                ! temps1_S=temps1! et on censure a la date de deces(ou censure)
            ! endif
        ! endif
        
        ! ! !print*,"temps1_S=",temps1_S,"delta",delta,"temps1_T=",temps1,"deltadc",deltadc

! !c****** for gap time :         
               ! t0(nobs) = 0.d0
! !c fin gap
               ! t1(nobs) = temps1
               ! t1_S(nobs) = temps1_S
               ! c(nobs) = delta
               ! cdc(nobs) = deltadc
               ! g(nobs)= ig
               ! iii = 0
               ! iii2 = 0
               ! do 110 ii = 1,ver
                  ! if(filtre(ii).eq.1)then
                     ! iii = iii + 1
                     ! ve(nobs,iii) = dble(v1(ii))
                  ! endif
                  ! if(filtre2(ii).eq.1)then
                     ! iii2 = iii2 + 1
                     ! ve2(nobs,iii2) = dble(v1(ii))
                  ! endif
 ! 110           continue

! !c*** pour le tester sur un autre programme: on complete les nouveaux parametres simules dans le jeux de donnees
            ! don_simulS1(ig,initTime1)=t0(nobs)
            ! don_simulS1(ig,timeS1)=t1_S(nobs)
            ! don_simulS1(ig,statusS1)=c(nobs)
            ! don_simulS1(ig,Patienref1)=g(nobs)
            ! don_simulS1(ig,trt1)=ve2(nobs,1)
            ! ! don_simulS1(ig,w_ij1)=ui
            ! !!print*,"don_simulS1(ig,statusS1)",ig,don_simulS1(ig,statusS1)
            ! !!print*,"don_simulS1(:,statusS1)",don_simulS1(:,statusS1)
            ! !!print*,"sum(don_simulS1(:,statusS1))",sum(don_simulS1(:,statusS1))
            
            ! don_simul(ig,initTime1)=t0(nobs)
            ! don_simul(ig,timeT1)=t1(nobs)
            ! don_simul(ig,statusT1)=cdc(nobs)
            ! don_simul(ig,Patienref1)=g(nobs)
            ! don_simul(ig,trt1)=ve2(nobs,1)
            ! ! don_simul(ig,w_ij1)=ui
            
            ! ! if(ibou.eq.Aretenir)then
            ! !    open(9,file="parametres.txt")
            ! !    !write(9,112)(t0(nobs)),t1_S(nobs),c(nobs),t1(nobs),cdc(nobs)&
            ! !       ,g(nobs),int(ve(nobs,1))&
            ! !        ,int(ve(nobs,2)),int(ve2(nobs,1))
! ! 112                 format(f7.1,f7.1,I2,f7.1,I2,I5,I2,I2,I2)
               ! !!print*,"size(g)",size(g)
            ! !    open(13,file="gastadv_T.txt")
            ! !    !write(13,113)1,g(nobs),int(ve2(nobs,1)),int(t0(nobs)),t1(nobs),cdc(nobs)
                    
! ! 113            format(I5,I5,I2,I2,f7.1,I2)
            ! !    open(8,file="gastadv_S.txt")
            ! !    !write(8,1131)1,g(nobs),int(ve(nobs,1)),int(t0(nobs)),t1_S(nobs),c(nobs)!&
                    ! !,int(ve(nobs,2))
 ! !1131           format(I5,I5,I2,I2,f7.1,I2)
            ! ! endif
! !c******************************************************************
! !c            
            ! if (maxtemps.lt.t1(nobs))then
               ! maxtemps = t1(nobs)
            ! endif  
! !c         !write(*,*)'*données**',t0(ig),t1(ig),c(ig),ve(ig,1),ve(ig,2),ig

            ! if (delta.eq.0.d0)then ! deces ou censure
               ! goto 30          ! on change de sujet
            ! endif
 ! 10      continue                   !observations par sujet
 ! 30   continue                  !sujets=groupes
 ! !                         !print*,"suis la"
      ! if(ibou.eq.Aretenir)then
        ! if(affiche_stat==1)then
            ! !write(*,*)'** nombre total d observations',nobs
            ! !write(*,*)"** nb d'evenements surrogate",nb_recur
            ! !write(*,*)'** nb donnees censurees ',nb_cens
            ! !write(*,*)'** nb deces ',nb_dc

            ! !write(*,*)'** proportion de deces (en %) ',nb_dc*100.d0/nobs
            ! !write(*,*)'** proportion de surrogate (en %) ',nb_recur*100.d0/nobs     
            ! !write(4,*)'** nombre total d observations',nobs
            ! !write(4,*)"** nb d'evenements surrogate",nb_recur
            ! !write(4,*)'** nb donnees censurees ',nb_cens
            ! !write(4,*)'** nb deces ',nb_dc
            ! !write(4,*)'** proportion de surrogate (en %) ',nb_recur*100.d0/nobs
            ! !write(4,*)'** proportion de deces (en %) ',nb_dc*100.d0/nobs
        ! endif
      ! moy_idnum=0
      ! do 444 jj=1,ng
         ! moy_idnum =moy_idnum + idnum(jj)
 ! 444     continue
          ! moy_idnum =moy_idnum / ng
    ! !  !write(*,*)'** nb moyen d observations par sujet ',moy_idnum,ng
    ! !  !write(4,*)'** nb moyen d observations par sujet ',moy_idnum
      
      ! endif
! !c****************************************************
! !c******** FIN generation des donnees****************
! !c****************************************************

! !c=========     on retourne a l'iteration suivante
 ! 1000 continue
    ! ! scl============censure conseillee pour la proportion souhaitee de censure==
        ! !!print*,tempsD
        ! call percentile_scl(tempsD,ng,1.d0-propC,cens)
        ! if(affiche_stat==1)then
            ! !print*,"la censure conseillée pour avoir:",propC*100,"% de personnes censurée vaut:",cens
            ! !print*,"la max de temps de deces vaux:",maxval(tempsD)
        ! endif
    ! deallocate(tempsD,sigma,mu,x_)
    ! !close(2)
    ! !close(4)
! endsubroutine Generation_surrogate_complet !FIN prog principal

! =========== subroutine pour la generation des donnees suivant une multinormal et vecteur de moyenne et da matrice de variance covariance donnees

subroutine rmvnorm(mu,vc1,nsim,vcdiag,ysim)
    ! mu: l'esperance de mes variables
    ! VC1: matrice de variance-covariance
    ! nsim: nombre de generations a faire
    ! vcdiag: un entier(1=oui, 0=non) qui dit si la matrice de variance covariance est diagonale ou pas. pour eviter la transformation de cholesky
    ! ysim: vecteur des realisations d'ne normale de moyenne mu et de matrice de covariance vc
        
    implicit none
    integer :: jj,j,k,ier,l,m,maxmes !maxmes= nombre de dimension ou encore dimension de X
    integer, intent(in)::nsim,vcdiag
    double precision::eps,ymarg,SX,x22 ! ymarg contient le resultat de l'integrale
    double precision, intent(in),dimension(:)::mu
    double precision,dimension(:,:),intent(in)::vc1
    double precision,dimension(:,:),allocatable::vc
    double precision,dimension(:),allocatable::usim
    !double precision,dimension(nsim,size(vc,2)),intent(out)::ysim
    double precision,dimension(:,:),intent(out)::ysim
    double precision,dimension(:),allocatable::vi
    
    !=============debut de la fonction=============================
    !!print*,vc
    !stop
    x22=0.d0
    maxmes=size(vc1,2)
    allocate(vi(maxmes*(maxmes+1)/2),usim((size(vc1,2))),vc(size(vc1,1),size(vc1,2)))
    vc=vc1
    jj=0
    Vi=0.d0
    do j=1,maxmes
        do k=j,maxmes
           jj=j+k*(k-1)/2
           Vi(jj)=VC(j,k)
        end do
    end do
    ! !print*,vi
    EPS=10.d-10
    !call dblepr("Vi =", -1, Vi, size(Vi))
    if(vcdiag.eq.0) then
        CALL DMFSD(Vi,maxmes,eps,ier) ! si matice diagonale on na pas besoin de ceci
    end if
    !!print*,vi
    if (ier.eq.-1) then
        call intpr("Problem with the cholesky transformation in the program", -1, ier, 1)
    else ! ysim sera un vecteur de 0
     
        VC=0.d0
        do j=1,maxmes
            do k=1,j
                VC(j,k)=Vi(k+j*(j-1)/2)
            end do
        end do    
        
        ! --------------------- Generation des donnees ------------------------
        ymarg=0.d0
        !!print*,vc
        !stop
        l=1
        do while(l.le.nsim)
            usim=0.d0
            do m=1,maxmes
                SX=1.d0
                call bgos(SX,0,usim(m),x22,0.d0) !usim contient des valeurs simulees d'une Normale centre reduite
            end do
            ysim(l,:)=mu+MATMUL(vc,usim) ! ysim contient des realisations d'une Normale de moyenne mu et de matrice de variance VC telle que chVC'chVC = VC
            l=l+1
        end do
    endif
            
    deallocate(vi,usim,vc)
    return
end subroutine rmvnorm

!subroutine pour la factorisation de cholesky 

subroutine Cholesky_Factorisation(vc)
    ! VC: matrice de variance-covariance (symetrique) a factoriser
        
    implicit none
    integer :: jj,j,k,ier,maxmes !maxmes= nombre de dimension ou encore dimension de X
    double precision::eps
    double precision,dimension(:,:),intent(inout)::vc
    double precision,dimension(:),allocatable::vi
    
    !=============debut de la fonction=============================

    maxmes=size(vc,2)
    allocate(vi(maxmes*(maxmes+1)/2))
    jj=0
    Vi=0.d0
    do j=1,maxmes
        do k=j,maxmes
           jj=j+k*(k-1)/2
           Vi(jj)=VC(j,k)
        end do
    end do
         
    EPS=10.d-10
    CALL DMFSD(Vi,maxmes,eps,ier)! fonction qui fait la factorisation de cholesky
    VC=0.d0
    if (ier.eq.-1) then
        !print*,"Probleme dans la transformation de cholesky pour la generation multinormale"
        ! stop
    else ! on retourne un vecteur de 0 car pas possible de transformer
        do j=1,maxmes
            do k=1,j
                VC(j,k)=Vi(k+j*(j-1)/2)
            end do
        end do    
    end if
    
end subroutine Cholesky_Factorisation


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

!===================Determinant of a Matrix (Fortran 90)=====================================

!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
double precision FUNCTION Determinant(matrix,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    double precision, DIMENSION(n,n) :: matrix,matrix2
    double precision :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    matrix2=matrix
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0.d0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0.d0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                Determinant = 0.d0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    Determinant = l
    DO i = 1, n
        Determinant = Determinant * matrix(i,i)
    END DO
    matrix=matrix2
END FUNCTION Determinant

! nouvelle definition avec un nouveau nombre
double precision FUNCTION Determinant_2(matrix, n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    double precision, DIMENSION(n,n) :: matrix,matrix2
    double precision :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    matrix2=matrix
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0.d0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0.d0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                Determinant_2 = 0.d0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    Determinant_2 = l
    DO i = 1, n
        Determinant_2 = Determinant_2 * matrix(i,i)
    END DO
    matrix=matrix2
END FUNCTION Determinant_2

    !fonction permettant de calculer le determinant d'une matrice carrée d'ordre 3
    double precision function determinant_scl_3(mat)

    implicit none

    double precision, dimension(3,3),intent(in)::mat
    double precision::resul
    if (size(mat,2).ne.3) then
        !print*,"desole ce programme n'est utilisé que pour des matrices carrées d'ordre 3"
        resul=-1.d9
        !stop
    else
        resul=mat(1,1)*(mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3))&
              -mat(2,1)*(mat(1,2)*mat(3,3)-mat(3,2)*mat(1,3))&
              +mat(3,1)*(mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3))
    endif
    
    determinant_scl_3= resul

    end function determinant_scl_3
    
    SUBROUTINE init_random_seed(graine,aleatoire,nbre_sim)
    ! Cette fonction permet de reinitialiser l'environnement de generation des nombre aleatoires lorsque l'on utilise la fontion RANDOM_NUMBER() pour la generation des donnees
    ! Elle est intéressante lorsqu'aucours de la procedure d'estimation, on fait d'autres appels a la fonction RANDOM_NUMBER(). pour pouvoir reproduire les donnees genere il est important 
    ! de reinitialiser l'espace de generation en se donnant une graine. ceci est l'equivalent de la fonction set.seed() de R 
    ! =======Version amelioree par casimir Ledoux SOFEU du programme de Matilde : 27/03/2018=======
    
        ! aleatoire: dit si on reinitialise la generation des nombre aleatoire avec un environnement different a chaque appel (1) ou non(O).
        ! En cas de generation differente, on utilise l'horloge (heure) de l'ordinateur comme graine. Dans ce cas, il n'est pas possible de reproduire les donnees simulees
        ! nbre_sim: dans le cas ou aleatoire=1, cette variable indique le nombre de generation a faire
        ! graine: dans le cas ou l'on voudrait avoir la possibilite de reproduire les donnees generees alors on met la variable aleatoire=0 et on donne dans cette variable la graine a utiliser pour la generation
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        integer, intent(in)::graine,nbre_sim,aleatoire
          
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
          
        CALL SYSTEM_CLOCK(COUNT=clock)
         
        if(aleatoire==1)then
            seed = clock + nbre_sim * (/ (i - 1, i = 1, n) /)
        else
            seed=graine
        endif
        CALL RANDOM_SEED(PUT = seed)
          
        DEALLOCATE(seed)
    END SUBROUTINE init_random_seed
    
    subroutine pos_proc_domaine(taille_domaine,nb_procs,rang,init_i,max_i)
        ! cette procedure permet pour chaque processus de calculer ses positions initiale et finale dans un domaine a nb_procs processus
        ! nb_procs: nombre de processus qui se partagent le domaine
        ! rang: rang du processus
        ! init_i: position initiale
        ! max_i : position finale
        ! taille_domaine : taille du sous domaine
        IMPLICIT NONE
        integer, intent(in)::rang,taille_domaine,nb_procs
        integer, intent(out)::init_i,max_i
        integer::suplement,n_par_pro
        integer,dimension(:),allocatable::table_par_pro
        
        n_par_pro=INT(taille_domaine/nb_procs)
        suplement=taille_domaine-n_par_pro*nb_procs ! donne le nombre de simulation a partager entre les premiers processus seulement
                    
        ! remplissage du table du nombre de simulation a effectuer par processus
        allocate(table_par_pro(nb_procs))
        table_par_pro(1:nb_procs)=n_par_pro
        table_par_pro(1:suplement)=n_par_pro+1
        
        n_par_pro=table_par_pro(rang+1) ! tous les processus jusqu'au rang supplement-1 recoivent une tâche supplementaire a realiser
        
        ! indice des calculs a effectuer par le processus courant
        if (rang==0) then
            init_i=1 ! ce processus commence a la premiere simulation
        else
            init_i=sum(table_par_pro(1:rang))+1 ! ce processus commence a la simulation qui respecte son ordre et doit s'arreter au nombre de simultation dont il a le droit d'executer
        endif
                
        max_i=init_i+table_par_pro(rang+1)-1!rang maximale de la simulation a executer (-1 car on a deja incrementer init_i de 1)
        deallocate(table_par_pro)
        
    endsubroutine pos_proc_domaine

end module Autres_fonctions
 