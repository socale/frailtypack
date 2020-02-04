! **** Author : Myriam L.
! **** Subject : Prediction of a new recurrent event, with a shared model
! **** Created on the 2016, 17 october
! **** Modified : 2016, 21 october by Myriam, L.

    subroutine predict_recurr_sha(LogN,npred0,surv_s,surv_t,surv_r,betapred,var,&
                 predAll,nreci,ntimeAll,icproba,nsample,varalea,&
                 surv_smc,surv_tmc,surv_rmc,betapredmc,predAlllow,predAllhigh)
    
      implicit none
      
      integer::i,iii,j
      integer,intent(in)::npred0,ntimeAll,icproba,nsample,LogN
      double precision,dimension(npred0,ntimeall),intent(in)::surv_t,surv_s
      double precision,dimension(npred0),intent(in)::surv_r,betapred
      double precision,intent(in)::var
      double precision,dimension(npred0)::predProba
      double precision,dimension(3)::survR,survRalea
      integer, dimension(npred0),intent(in)::nreci
      double precision,dimension(npred0,ntimeAll),intent(out):: predAll      
      double precision::ss1,ss2, nrecj, Xbeta,Xbetamc
      double precision,dimension(nsample*npred0,ntimeAll)::surv_tmc,surv_smc      
      double precision,dimension(nsample,npred0)::surv_rmc
      double precision,dimension(nsample),intent(in) ::varalea
      double precision,dimension(nsample,npred0)::predProbaalea
      double precision,dimension(npred0,ntimeAll) :: predAlllow,predAllhigh,betapredmc
            
      do iii=1,ntimeAll
          do i=1,npred0
              survR(1) = surv_s(i,iii) ! Recuperation des valeurs de survie calculees dans R
              survR(2) = surv_t(i,iii)
              survR(3) = surv_r(i)
              nrecj = nreci(i)
              Xbeta = betapred(i)
              if (LogN.eq.1) then
                  call gauher_LogNsha(ss1,ss2,var,survR,nrecj,Xbeta)
              else 
                  call gaulag_gammasha(ss1,ss2,var,survR,nrecj,Xbeta)
              endif
              predProba(i) = ss1/ss2       
          end do
          
          predAll(:,iii) = predProba    

      !================= Intervalle de confiance==
          if (icproba.eq.1) then !(Calcul uniquement si demande)
              do j=1,nsample !(MC.sample)
                  ss1 = 0.d0 
                  ss2 = 0.d0
                  do i=1,npred0  !(nb patients)             
                      survRalea(1) = surv_smc(npred0*(j-1)+i,iii)
                      survRalea(2) = surv_tmc(npred0*(j-1)+i,iii)
                      survRalea(3) = surv_rmc(j,i)
                      nrecj = nreci(i)
                      
                      Xbetamc = betapredmc(i,j)

                      if (LogN.eq.1) then
                          call gauher_LogNsha(ss1,ss2,varalea(j),survRalea,nrecj,Xbetamc)
                      else
                          call gaulag_gammasha(ss1,ss2,varalea(j),survRalea,nrecj,Xbetamc)
                      endif
                      predProbaalea(j,i) = ss1/ss2
                  end do
              end do
    
            ! utilisation de la fonction percentile2 de aaUseFunction
              do i=1,npred0
                  call percentile2(predProbaalea(:,i),nsample,predAlllow(i,iii),predAllhigh(i,iii))
              end do            
          endif

      end do 

    end subroutine predict_recurr_sha    

!=========================
! Prediction : numerator
!=========================
    double precision function func1predShaRec(frail,theta,survR,nrecj,Xbeta)
      ! calcul des integrandes (numerateur de la fonction de prediction)
      implicit none

      double precision,intent(in)::frail
      double precision,dimension(3)::survR
      double precision::theta,nrecj,logGammaJ,Xbeta

      func1predShaRec = ((survR(1)**(frail*Xbeta))-(survR(2)**(frail*Xbeta))) & 
                      *(frail**nrecj)*(survR(3)**(frail*Xbeta)) &
                      *(frail**(1.d0/theta -1.d0) * exp(-frail/theta))&
                      /(theta**(1.d0/theta) * dexp(logGammaJ(1.d0/theta)))
      return
    end function func1predShaRec

!=========================
! Prediction : denominator
!=========================
    double precision function func2predShaRec(frail,theta,survR,nrecj,Xbeta)
    ! calcul des integrandes (denominateur de la fonction de prediction)
      implicit none

      double precision,intent(in)::frail
      double precision,dimension(3)::survR
      double precision::theta,nrecj,logGammaJ,Xbeta

      func2predShaRec = (survR(1)**(frail*Xbeta))*(frail**nrecj)&
             *(survR(3)**(frail*Xbeta))*(frail**(1.d0/theta -1.d0) &
       * exp(-frail/theta)) / (theta**(1.d0/theta) *dexp( logGammaJ(1.d0/theta)))
      return
    end function func2predShaRec


!=========================
! Calcul des intégrales
!=========================
    subroutine gaulag_gammasha(ss1,ss2,theta,survR,nrecj,Xbeta)
    
      use donnees,only:w1,x1     
      implicit none
      
      double precision,intent(out)::ss1,ss2
      double precision::var1,Xbeta
      double precision::auxfunca1,auxfunca2
      double precision,external :: func1predShaRec,func2predShaRec
      double precision,dimension(3)::survR
      double precision::theta,nrecj
      integer:: j     
      
      ss1=0.d0
      ss2=0.d0
      
      do j=1,32
          var1 = x1(j)
          auxfunca1 = func1predShaRec(var1,theta,survR,nrecj,Xbeta)
          ss1 = ss1 + w1(j)*(auxfunca1)
          auxfunca2 = func2predShaRec(var1,theta,survR,nrecj,Xbeta)
          ss2 = ss2 + w1(j)*(auxfunca2)             
      end do

    end subroutine gaulag_gammasha

!***************************************
!***************************************
!******* Distribution LogN *************
!***************************************
!***************************************

!=========================
! Prediction : numerator
!=========================
    double precision function func1predShaRecLogN(frail,psig2,survR,nrecj,Xbeta)
      ! calcul des integrandes (numerateur de la fonction de prediction)
      implicit none
      
      double precision,intent(in)::frail
      double precision,dimension(3)::survR
      double precision::psig2,nrecj,Xbeta
      double precision,parameter::pi=3.141592653589793d0

      func1predShaRecLogN = &
          ((survR(1)**(dexp(frail)*Xbeta))-(survR(2)**(dexp(frail)*Xbeta)))&
            *(dexp(frail)**nrecj) * (survR(3)**(dexp(frail)*Xbeta)) &
            *dexp(-(frail**2.d0/(2.d0*psig2)))&
            *(1.d0/dsqrt(2.d0*pi*psig2))
      
      return    
    end function func1predShaRecLogN

!=========================
! Prediction : denominator
!=========================
    double precision function func2predShaRecLogN(frail,psig2,survR,nrecj,Xbeta)
    ! calcul des integrandes (denominateur de la fonction de prediction)
      implicit none
      
      double precision,intent(in)::frail
      double precision,dimension(3)::survR
      double precision::psig2,nrecj,Xbeta
      double precision,parameter::pi=3.141592653589793d0
      
      func2predShaRecLogN = (survR(1)**(dexp(frail)*Xbeta))&
             *(dexp(frail)**nrecj)*(survR(3)**(dexp(frail)*Xbeta))&
             *dexp(-(frail**2.d0/(2.d0*psig2)))&
             *(1.d0/dsqrt(2.d0*pi*psig2))
      return    
    end function func2predShaRecLogN

!=========================
! Calcul des intégrales 
!=========================
    subroutine gauher_LogNsha(ss1,ss2,psig2,survR,nrecj,Xbeta)
    
      use donnees,only:x3,w3      
      implicit none
      
      double precision,intent(out)::ss1,ss2
      double precision::var1,Xbeta
      double precision::auxfunca1,auxfunca2
      double precision,external :: func1predShaRecLogN,func2predShaRecLogN
      double precision,dimension(3)::survR
      double precision::psig2,nrecj
      integer:: j     
      
      ss1=0.d0
      ss2=0.d0
      
      do j=1,32
          var1 = x3(j)
          auxfunca1 = func1predShaRecLogN(var1,psig2,survR,nrecj,Xbeta)
          ss1 = ss1 + w3(j)*(auxfunca1)
          auxfunca2 = func2predShaRecLogN(var1,psig2,survR,nrecj,Xbeta)
          ss2 = ss2 + w3(j)*(auxfunca2)             
      end do

    end subroutine gauher_LogNsha