
!AK 19/02/2015 shared log normal predictions and MC confidence bands
! ============================================== prediction Shared - logNormale
    
    subroutine predict_logn_sha(npred0,surv_s,surv_t,betapred,sigma2, &
        predAll,icproba,ntimeAll,nsample,sig2alea,surv_smc,surv_tmc, &
        betapredmc,predAlllow,predAllhigh)
    
    implicit none
    
    integer::i,iii,j
    integer,intent(in)::npred0,icproba,ntimeAll,nsample
    double precision,dimension(npred0,ntimeall),intent(in)::surv_t,surv_s
    double precision,dimension(npred0),intent(in)::betapred
    double precision,intent(in)::sigma2
    double precision,dimension(npred0)::predProba
    double precision,dimension(2)::survDC,survDCalea
    double precision,dimension(nsample*npred0,ntimeAll)::surv_tmc,surv_smc
    double precision,dimension(npred0,ntimeAll) :: predAlllow,predAllhigh,betapredmc
    double precision,dimension(nsample),intent(in) :: sig2alea 
    double precision,dimension(npred0,ntimeAll),intent(out):: predAll
    double precision,dimension(nsample,npred0)::predProbaalea
    double precision::ss1,ss2,Xbeta,Xbetamc
    
    do iii=1,ntimeAll
        do i=1,npred0
            survDC(1) = surv_s(i,iii)
            survDC(2) = surv_t(i,iii)
            Xbeta = betapred(i) 
    
            call gauher_shapred(ss1,ss2,sigma2,survDC,Xbeta)
            predProba(i) = ss1/ss2
         
        end do
        
       predAll(:,iii) = predProba
        !=============================================
        ! Variabilite des proba predites
        if (icproba.eq.1) then ! calcul de l'intervalle de confiance seulement si demande
            do j=1,nsample
                ss1 = 0.d0
                ss2 = 0.d0
                do i=1,npred0 
              
                    survDCalea(1) = surv_smc(npred0*(j-1)+i,iii)
                    survDCalea(2) = surv_tmc(npred0*(j-1)+i,iii)
                    Xbetamc = betapredmc(i,j)
              
                    call gauher_shapred(ss1,ss2,sig2alea(j),survDCalea,Xbetamc)
                    predProbaalea(j,i) = ss1/ss2
             
                end do
            end do

            ! utilisation de la fonction percentile2 de aaUseFunction
            do i=1,npred0
                call percentile2(predProbaalea(:,i),nsample,predAlllow(i,iii),predAllhigh(i,iii))
            end do
            
        endif ! calcul de l'intervalle de confiance seulement si demande

    end do
 
    end subroutine predict_logn_sha
    
    
!=========================
! Prediction : numerator
!=========================
    double precision function func1predLogN(frail,psig2,survDC,Xbeta)
    ! calcul de l integrant (numerateur de la fonction de prediction)
    implicit none
    
    double precision,intent(in)::frail
    double precision,dimension(2)::survDC
    double precision::psig2,Xbeta
    double precision,parameter::pi=3.141592653589793d0
    func1predLogN =  &
      ((survDC(1)**(dexp(frail)*Xbeta))-(survDC(2)**(dexp(frail)*Xbeta))) &
      *dexp(-(frail**2.d0/(2.d0*psig2)))*(1.d0/dsqrt(2.d0*pi*psig2))
     
    return
    
    end function func1predLogN

!=========================
! Prediction : denominator
!=========================

    double precision function func2predLogN(frail,psig2,survDC,Xbeta)
    ! calcul de l integrant (denominateur de la fonction de prediction)
    implicit none
    
    double precision,intent(in)::frail
    double precision,dimension(2)::survDC
    double precision::psig2,Xbeta
    double precision,parameter::pi=3.141592653589793d0
    
    func2predLogN = (survDC(1)**(dexp(frail)*Xbeta))*dexp(-(frail**2.d0/(2.d0*psig2)))&
        *(1.d0/dsqrt(2.d0*pi*psig2))
    
    return
    
    end function func2predLogN


!=========================
! Calcul des intégrales
!=========================
    subroutine gauher_shapred(ss1,ss2,psig2,survDC,Xbeta) !! 
    
!    use tailles  
    use donnees,only:w3,x3
    
    implicit none

    double precision,intent(out)::ss1,ss2
    double precision::var1
    double precision::auxfunca1,auxfunca2
    double precision,external :: func1predLogN,func2predLogN
    double precision,dimension(2)::survDC
    double precision::psig2,Xbeta
    integer:: j

! gauss hermite
  
    ss1=0.d0
    ss2=0.d0
    
    do j=1,32
        var1 = x3(j)
        auxfunca1 = func1predLogN(var1,psig2,survDC,Xbeta)
        ss1 = ss1 + w3(j)*(auxfunca1)
        auxfunca2 = func2predLogN(var1,psig2,survDC,Xbeta)
        ss2 = ss2 + w3(j)*(auxfunca2)           
    end do    
    end subroutine gauher_shapred
    
    
    
    !==========================================
    !=================== AK frailty.mc ==============
    !======== pour les predictions conditionnelles - shared model
    !========= distribution normal ====================
    
    
    subroutine frailpred_sha_nor_mc(np0,frailtypred,sig20,res10,nig0)
    
   
    use optimres
    use residusM

    implicit none
    
    integer,intent(in)::np0,nig0
    double precision,external::funcpasres_mc
    double precision,intent(out)::frailtypred
    double precision,intent(in)::res10,sig20

        sig2_mc = sig20
        res1_mc = res10
        nig_mc = nig0
        np_mc = np0
    
    allocate(vuu(2),vecuiRes(1),vres((1*(1+3)/2)))
 
    vecuiRes=0.d0
    moyuiR=0.d0
    varuiR=0.d0
    
    cares=0.d0
    cbres=0.d0
    ddres=0.d0
    frailtypred = 0.d0
  
    
    
            vuu=0.9d0
            
            call marq98res(vuu,1,nires,vres,rlres,ierres,istopres,cares,cbres,ddres,funcpasres_mc)

            
            if (istopres.eq.1) then
                frailtypred = vuu(1)*vuu(1)
                       else
                ! non convergence ou erreur de calcul de la fonction a maximiser
            
                frailtypred= 0.d0
            
            endif
     
    deallocate(vuu,vecuiRes,vres)
    end subroutine frailpred_sha_nor_mc

      
!========================    FUNCPARES_MC FRAILPRED DENSITE A POSTERIORI       ====================

!!!!
!!!! Calcul frailtypred shared log-normal mc
!!!!
    double precision function funcpasres_mc(uu,np,id,thi,jd,thj)

   use residusM

    implicit none

    integer,intent(in)::id,jd,np
    double precision,dimension(np)::bh
    double precision,dimension(np),intent(in)::uu
    double precision,intent(in)::thi,thj
    double precision::frail

    
    !np = np_mc
    bh=uu
    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    frail = bh(1)*bh(1)

    funcpasres_mc = dexp(nig_mc*frail - &
    dexp(frail)*res1_mc - (frail**2.d0)/(2.d0*sig2_mc))
  
    return

    end function funcpasres_mc
