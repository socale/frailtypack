module mod_Adaptative
    implicit none
    
    contains
    !========================          fonction qui une fois lancee permet l'estimation a posteriorie des effets aleatoires        ====================
    double precision function funcpa_estimatin_fragilite(b,np,id,thi,jd,thj,k0)
    
    use tailles
    use comon
    use residusM
    use var_surrogate ! fichier Aparameters_scl.f90
    use fonction_A_integrer ! integrant (fichier Integrant_scl.f90)
    use GaussHermi_mult ! pour la fonction d'integration (fichier Integrale_mult_scl.f90)
    use monteCarlosMult_Gaus ! pour integration par monte carlo (fichier GuassHermi_mult_scl.f90)
    use InverseMatrix !se trouve dans le fichier autres_fonctions.f90, pour l'inverse de la matrice des effets aleatoires et le calcul du determinant
        
    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :
    
    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    
    integer::n,i,j,k,vj,ig,choix,l,nsujet_trial,dimint !vcdiag
    integer,dimension(ngmax)::cpt
    double precision::inv,som1,som2,res,vet,vet2,h1,varS1,varT1,covST1 !pe1,pe2,som,inc
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
    double precision,dimension(2,2):: mat_A

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
    n = (np-nva-nparamfrail)/nst
    !n = (np-nva-effet-indic_ALPHA)/nst
    !!print*,"nparamfrail=",nparamfrail,"effet=",effet,"indic_ALPHA=",indic_ALPHA
    !!print*,"(np-nva-nparamfrail)/nst",(np-nva-nparamfrail)/nst,"n=",n

    do i=1,n
        the1(i-3)=(bh(i))*(bh(i))
        j = n+i 
        if (nst.eq.2) then
            the2(i-3)=(bh(j))*(bh(j))
        endif
    end do
    
    !!print*,indice_eta+indice_theta+indice_varS+indice_varT+indice_covST
    !stop
!!print*,"bh=",bh
    
    if(effet.eq.1) then
        if(logNormal==1)then
            theta2 = bh(np-nva-nparamfrail+indice_eta+indice_theta)**2.d0 ! scl on recupere theta du vecteur des parametre, au carree car c'est bien la variance
            !!print*,"theta2=",theta2
            varS1 = bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS)
            varT1 = bh(np-nva-nparamfrail+indice_eta+indice_theta+indice_varS+indice_varT)
            sig2=theta2 ! je fais appel a sig2 car c'est la variable utilisee dans la suite des procedures pour le joint classique
            !!print*,"vs,vt,covst=",varS1,varT1,covST1
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
    endif
    !!print*,np,nva,nparamfrail,indice_eta,indice_theta,indice_varS,indice_varT,indice_covST
    ! pour eviter d'avoir des matrices de variances-covariances non defini positive, je suppose que c'est la cholesky qui est generee. par consequent sigma=Chol*Chol^T
    !Chol: matrice triangulaire inferieur. pour eviter de refaire la factorisation de cholesky pour l'algo MC, j'utilise directement cette matrice de cholesky a la place de la matrice de variance-covariance
    Chol=0.d0 
    Chol(1,1)=varS1
    Chol(2,1)=covST1
    !Chol(1,2)=covST1
    Chol(2,2)=varT1
    mat_A=MATMUL(Chol,TRANSPOSE(Chol))
    varS=mat_A(1,1)
    varT=mat_A(2,2)
    covST=mat_A(1,2)
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
            
        if((c(i).eq.1))then
            res2s(pourtrial(i)) = res2s(pourtrial(i))+dlog(dut1(nt1(i)))+vet
        endif  
        if ((res2s(pourtrial(i)).ne.res2s(pourtrial(i))).or.(abs(res2s(pourtrial(i))).ge. 1.d30)) then
            funcpa_estimatin_fragilite=-1.d9
            goto 123
        end if    
        
        ! scl pour calcul integrale
        
        const_res4(g(i)) = const_res4(g(i))+ut1(nt1(i))* dexp(vet)
        if ((const_res4(g(i)).ne.const_res4(g(i))).or.(abs(const_res4(g(i))).ge. 1.d30)) then !scl si risque cumule >10^30
            funcpa_estimatin_fragilite=-1.d9
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
        if(cdc(k).eq.1)then
            !!print*,"k=",k,"size(nt1dc)",size(nt1dc),"nt1dc(k)=",nt1dc(k),"dut2(nt1dc(k))=",dut2(nt1dc(k))
            !!print*,"dut2(1:10)",dut2(1:10)
            res2_dcs(pourtrial(k)) =res2_dcs(pourtrial(k))+dlog(dut2(nt1dc(k)))+vet2
            !!print*,"res2_dcs(pourtrial(k))=",res2_dcs(pourtrial(k))
            !stop
            if ((res2_dcs(pourtrial(k)).ne.res2_dcs(pourtrial(k))).or.(abs(res2_dcs(pourtrial(k))).ge. 1.d30)) then
                funcpa_estimatin_fragilite=-1.d9
                goto 123
            end if    
        endif 
      
        const_res5(g(k)) = const_res5(g(k))+ut2(nt1dc(k))* dexp(vet2) 
!        !print*,"je suis la bien funcpajsplines_surr ligne 260, k=",k,"const_res5(k)=",const_res5(k)
        if ((const_res5(g(k)).ne.const_res5(g(k))).or.(abs(const_res5(g(k))).ge. 1.d30)) then
            funcpa_estimatin_fragilite=-1.d9
            goto 123
        end if
    end do

!!print*,"size(const_res5)",size(const_res5)
!**************INTEGRALES ****************************
!!print*,"ng=",ng,"ntrials=",ntrials

    !================================================================================
    !==========distribution lognormale des effects aleatoires==============================
    !================================================================================
    if (logNormal==1) then 
        posind_i=1
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
            
        ! on lance ceci juste pour besoin de la prédiction des effests aleatoires            
        ig=1 
        nsujet_trial=nsujeti(ig)        
        resultatInt=0.d0
        dimint=2 ! deux integrations au niveau essai correspondant aux effets aleatoires correles
        resultatInt=gauss_HermMultInd_Essai(Integrale_Individuel,gauss_HermMultA_surr,npoint,dimint,nsujet_trial,ig)                
    endif
    
    123 continue
    return
    endfunction funcpa_estimatin_fragilite

endmodule mod_Adaptative 