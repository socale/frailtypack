!mtaille =c(mt1,mt2,mt11,mt12)
!paraweib =c(shapeweib(1),shapeweib(2),scaleweib(1),scaleweib(2))
!kendall: 1ere colonne ss0, 2eme colonne tau
!paratps = c(timedep0,nbinnerknots,qorder0)
!counts = c(ni,cpt,cpt_dc)
!noVar = c(noVar1,noVar2,noVar3)

    
    
    
    !--entC*te pour fortran
     subroutine longiuninl(nsujety0,ng0,yy0,groupey0,nb0,which_random0,box_cox0,matzy0,cag0,&
     nva30,nva40,vaxy0,noVar,maxit0,np,b,H_hessOut,HIHOut,resOut,LCV,&
     counts,ier_istop,EPS,GH,paGH,b_pred,weights0,nodes0,nnodes_all0,initialGH,axT)
    
    !AD: add for new marq
        use donnees_indiv
        use parameters
        use comon
        use tailles
        use lois_normales
        use optim
   !   use ParametresPourParallelisation
    !AD:pour fortran
        use sortie
        use residusM
        use comongroup
        use splines
    !AD:
        implicit none
    
        integer::maxit0,initialGH
        integer,intent(in)::nsujety0,ng0,nva30,nva40,nb0,nnodes_all0,which_random0
        integer::np
        integer,dimension(nsujety0),intent(in) :: groupey0
        double precision,dimension(2),intent(in) :: cag0,box_cox0
         integer,dimension(2),intent(in):: GH
         double precision,dimension(ng0,nb0+1+nb0 + (nb0*(nb0-1))/2)::b_pred
        double precision,dimension(nnodes_all0,nb0),intent(in):: nodes0,weights0 
        double precision,dimension(2)::axT
    !    integer::typeJoint0
    
        double precision,dimension(nsujety0,nva30+nva40),intent(in):: vaxy0
        double precision,dimension(nsujety0,2) :: matzy0
        double precision,dimension(nsujety0) :: yy0
        double precision,dimension(np,np)::H_hessOut,HIHOut
        double precision::resOut
        integer::ss,sss
        double precision,dimension(np):: b
        double precision,dimension(2),intent(out)::LCV
       
        integer,dimension(2),intent(in)::noVar
        integer::noVar3,noVar4!! rajout
        integer::cpt,cpt_dc,ni
        integer,dimension(2),intent(out)::ier_istop
        integer,dimension(3),intent(out)::counts
        integer::groupey,j,ii,iii,   &
        i, &
        ier,istop
        double precision::res,max, &
        rl_temp !! rajout
    !AD: add for new marq
        double precision::ca,cb,dd
        double precision,external::funcpalongi_uni
        double precision,dimension(2)::k0
    
  !  double precision,dimension(ng0,nb0+1+nb0 + (nb0*(nb0-1))/2)::b_paGH0
    double precision,dimension(1,nva30+nva40)::coefBetaY
     
    double precision,external::funcpares_uni,funcprep
    
        integer::ngtemp
    
         double precision,dimension(3),intent(inout)::EPS ! seuils de convergence : on recupere les valeurs obtenues lors de l'algorithme a la fin
       
            double precision,dimension(ng0,nb0+1+nb0 + (nb0*(nb0-1))/2):: paGH
                 !double precision,dimension(:,:),allocatable :: mat_sigma,varcov_marg_inv
            !double precision,dimension(:),allocatable :: matv
            !double precision,dimension(nva30,nva30) :: element
        
 

    allocate(nodes(nnodes_all0,nb0),weights(nnodes_all0,nb0))!,b_lme_nnodes(GH(2)**(nb0+effet0),ng0,nb0))
        
    nodes = nodes0

    weights= weights0
    
    nnodes_all = nnodes_all0 !GH(2)**(nb0+effet0)
    
    which_random = which_random0

     k0(1) = axT(1)
     k0(2) = axT(2)
    
            ier = ier_istop(1)
            istop = ier_istop(2)
 
      
    
        typeJoint = 0
    
   
        s_cag_id = int(cag0(1))
        s_cag = cag0(2)
        
        box_cox1 = int(box_cox0(1))
        if(box_cox1.eq.1) then 
        box_cox_par = box_cox0(2)
        else 
        box_cox_par = int(box_cox0(2))
        end if
    
        maxiter = maxit0
    !AD:add for new marq
    !write(*,*)'eps',eps
        epsa = EPS(1) !1.d-4
        epsb = EPS(2) !1.d-4
        epsd = EPS(3) !1.d-4
    
        epsa_res = 1.d-4
        epsb_res = 1.d-4
        epsd_res = 1.d-4
    
    
        genz(1) = 30
        genz(2) = 500
    
            
    
        ngmax=ng0
        ng=ng0
    
    
    
            ngtemp=ng
    
     
        allocate(vaxy(nva30+nva40),ResidusLongi(nsujety0),Pred_y(nsujety0,2))
    
  
    
     
        nsujetymax = nsujety0
        nsujety = nsujety0
    

    
        allocate(yy(nsujetymax)) !! chgt dimension aux
    
        !** number of random effects ***
        nea = nb0

        nb1 = nb0
        !********* longitudinal data *************

        yy = yy0
        nb_re = nb0 !+ (nb0*(nb0-1))/2.d0
       
   
  
  
    
    
            allocate (Ut(nea,nea),Utt(nea,nea),ziy(nsujety0,2),sum_mat(nva30,nva30),matb_chol(nea,nea),mat(nea,nea))!,&
    !    mat_all(nnodes_all*nea,nnodes_all*nea))
    
        ziy = matzy0
  
    
        allocate(vuu(nea))
        !*** find the number of recurrent measures per patient
        allocate(nmesy(ng),nmes_o(ng))
        allocate(groupeey(nsujety))
   
        nmesy = 1
        
        groupeey = groupey0
     
        
        nmes_o = 0
    
        
        
       
    
        i = 1
        do j=2,nsujety
            if(groupeey(j).eq.groupeey(j-1))then
                nmesy(i)=nmesy(i)+1
                if((s_cag_id.eq.1).and.(yy(j).gt.s_cag))nmes_o(i) = nmes_o(i)+1
                                       
            else
            if((s_cag_id.eq.1).and.(yy(j).gt.s_cag))nmes_o(i) = nmes_o(i)+1
                i = i+1
            end if
        end do
        if((s_cag_id.eq.1).and.(yy(nsujety).gt.s_cag))nmes_o(i) = nmes_o(i)+1
        maxmesy=0
    
        do i=1,ng
            if (nmesy(i).gt.maxmesy) then
                maxmesy=nmesy(i)
            end if
        end do
    
      if(s_cag_id.eq.0)nmes_o = nmesy
    

    
        allocate(mu1_res(maxmesy),mu1(maxmesy,1))
        allocate(varcov_marg(nsujety,maxmesy))
    
    
    
        nst=1
      
   

    !------------  entre non fichier et nombre sujet -----
    
    
        nva3=nva30
        nva4=nva40    
        

        noVar3 = noVar(1)
        noVar4 = noVar(2)
        
      nva = nva3 + nva4
        nvarmax=nva
        allocate(vey(nsujetymax,nvarmax))
        allocate(ve3(nsujetymax,nva3+nva4))
        allocate(filtre3(max(1,nva30)),filtre4(max(1,nva40)))
    
  
    ! AK: longitudinal growth
        if (noVar3.eq.1) then
    !        write(*,*)'filtre 3 desactive'
            filtre3=0
            nva3=0
        else
            filtre3=1
        end if
        
        
        ! AK: longitudinal growth
        if (noVar4.eq.1) then
    !        write(*,*)'filtre 3 desactive'
            filtre4=0
            nva4=0
        else
            filtre4=1
        end if
    
    
        nva = nva3 + nva4

    !AD:end


    
        !cccccccccccccccccccccccccccccccccc
    ! pour les donnees longitudinales
    !cccccccccccccccccccccccccccccccccc
        vaxy = 0
        do i = 1,nsujety     !sur les observations
            groupey=groupey0(i)
            do j=1,nva3
                vaxy(j)=vaxy0(i,j)
            end do
            iii = 0
            do ii = 1,nva3
                if(filtre3(ii).eq.1)then
                    iii = iii + 1
                    vey(i,iii) = dble(vaxy(ii)) !ici sur les observations
                endif
            end do
             do j=1,nva4
                vaxy(j)=vaxy0(i,nva3+j)
            end do
             do ii = 1,nva4
                if(filtre4(ii).eq.1)then
                    iii = iii + 1
                    vey(i,iii) = dble(vaxy(ii)) !ici sur les observations
                endif
            end do
        end do
  
        deallocate(filtre3,filtre4)

   
        npmax=np
    
    
        allocate(hess(npmax,npmax),I1_hess(npmax,npmax) &
        ,H1_hess(npmax,npmax),I2_hess(npmax,npmax),H2_hess(npmax,npmax) &
        ,HI2(npmax,npmax),HIH(npmax,npmax),IH(npmax,npmax),HI(npmax,npmax),BIAIS(npmax,1))
        !Hspl_hess(npmax,npmax),PEN_deri(npmax,1),
    !,I3_hess(npmax,npmax),H3_hess(npmax,npmax),HI3(npmax,npmax)&
        
    !------- initialisation des parametres
    
        ! savoir si l'utilisateur a entre des parametres initiaux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        if (sum(b).eq.0.d0) then
            do i=1,npmax
                b(i)=5.d-1
            end do
        endif
    
       !     b(np-nva-nb_re) = 1.d0
    
      
    
    !############################
        ca=0.d0
        cb=0.d0
        dd=0.d0
    
    
    
        allocate(Z1(maxmesy,nb0),Zet(nsujetymax,nb0))
        allocate(mu(maxmesy,2),ycurrent(maxmesy),b1(np))!,mu_all(nsujety,nnodes_all)
    !    allocate(ycurrent_all(nnodes_all,nsujety),param_alnorm(nnodes_all,22),s_cag_all(nnodes_all,1),&
    !    upper_all(nnodes_all),aux_all(nnodes_all,ng),cdc_all(nnodes_all,ng),weights_all(nnodes_all))
      
        allocate(x2(maxmesy,nva3),x2cur(1,nva3),z1cur(1,nb1),current_mean(1))
        
        allocate(I_hess(np,np),H_hess(np,np),v((np*(np+3)/2)))
    
       allocate(vvv((np*(np+1)/2)))
        
        
        
         if(initialGH.eq.1) then
     allocate(H_hess_GH(nea,nea),I_hess_GH(nea,nea),b_paGH(ng0,nea+1+nea + (nea*(nea-1))/2))
    
   
                allocate(vecuiRes2(ng,nea),&
                        vres(nea*(nea+3)/2))
          
    coefBetaY = 0.d0
    coefBetaY(1,:) = b((np-nva3-nva4+1):np)
         
                b1 = b
                npp = np
           rl_temp=funcprep(b1,npp)
         
     Call Residus_uni(b,np,funcpares_uni)
    
        
 paGH(1:ng,1:nea+1+nea + (nea*(nea-1))/2) = b_paGH(1:ng,1:nea+1+nea + (nea*(nea-1))/2)
    
                deallocate(vecuiRes2,vres,b_paGH,H_hess_GH,I_hess_GH)
    
   end if
   
            ! Parametres pour GH pseudo-adaptative
            allocate(b_lme(ng,nea),invBi_cholDet(ng),invBi_chol(ng,nea + (nea*(nea-1))/2))
                    do i=1,ng
                            b_lme(i,1:nea) = paGH(i,1:nea)
                            invBi_cholDet(i) = paGH(i,nea+1)
                            
                            invBi_chol(i,1:nea + (nea*(nea-1))/2) = paGH(1,(nea+2):(nea+1+nea + (nea*(nea-1))/2))
                            
                    end do
    
            methodGH = GH(1)
            nodes_number = GH(2)
        
        
          allocate(H_hess_GH(nea,nea), vres(nea*(nea+3)/2))
        
        model = 1
        v = 0.d0
        res = 0.d0
        
        call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpalongi_uni)
   
    
    deallocate(H_hess_GH, vres)
        
        resOut=res
    !Al:
        EPS(1) = ca
        EPS(2) = cb
        EPS(3) = dd
    
        ier_istop(1) = ier
        ier_istop(2) = istop
    !Al:
        if (istop.ne.1) goto 1000
    
        call multiJ(I_hess,H_hess,np,np,np,IH)
        call multiJ(H_hess,IH,np,np,np,HIH)
    
    
       
    
    
    
        do ss=1,npmax
            do sss=1,npmax
                HIHOut(ss,sss) = HIH(ss,sss)
                H_hessOut(ss,sss)= H_hess(ss,sss)
            end do
        end do
    
    !AD:add LCV
    !LCV(1) The approximate like cross-validation Criterion
    !LCV(2) Akaike information Criterion
    !     calcul de la trace, pour le LCV (likelihood cross validation)
        LCV=0.d0
     
    !        write(*,*)'=========> Akaike information Criterion <========='
          
                LCV(2) = (1.d0 / (nsujety)) *(np - resOut)
          
    
        1000 continue
    !AD:end

                
      
    !    write(*,*)'======== Calcul des residus de martingale ========'
    
        deallocate(I_hess,H_hess)
    
    
    allocate(H_hess_GH(nea,nea),I_hess_GH(nea,nea),b_paGH(ng0,nea+1+nea + (nea*(nea-1))/2))
    
   
                allocate(vecuiRes2(ng,nea),&
                        vres(nea*(nea+3)/2))
          
    coefBetaY = 0.d0
    coefBetaY(1,:) = b((np-nva3-nva4+1):np)
       
                b1 = b
                npp = np
             b1 = b
                npp = np
           rl_temp=funcprep(b1,npp)
             
     Call Residus_uni(b,np,funcpares_uni)
        
  b_pred(1:ng,1:nea+1+nea + (nea*(nea-1))/2) = b_paGH(1:ng,1:nea+1+nea + (nea*(nea-1))/2)
    
                deallocate(vecuiRes2,vres,b_paGH)
    
    
    
    
    
    
    
       

        counts(1) = ni
        counts(2) = cpt
        counts(3) = cpt_dc
    
 
        deallocate( vey,vaxy)
    
     deallocate(H_hess_GH,I_hess_GH)
        deallocate(hess,v,I1_hess,H1_hess,I2_hess,H2_hess,HI2,HIH,IH,HI,BIAIS)
    
        deallocate(ResidusLongi,Pred_y)
        deallocate(ve3,Zet)
    
    
        deallocate(Ut,Utt,varcov_marg,sum_mat,mat)!,mat_all)

    
        deallocate(ziy,b1,yy)
        deallocate(nmesy,groupeey,nmes_o,mu1_res)
    
    !   deallocate(I3_hess,H3_hess,HI3)
    
        deallocate(Z1,mu,ycurrent)!,mu_all
       
    deallocate(x2,x2cur,z1cur,current_mean,mu1)
    
      
        deallocate(b_lme,invBi_chol,invBi_cholDet,matb_chol)
    
        deallocate(nodes,weights,vuu,vvv)!,b_lme_nnodes)
    
        return
    
        end subroutine longiuninl
    
    

    !=========================================================
        
        
        
         SUBROUTINE gauherJ1_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ1_uni
        
        
         SUBROUTINE gauherJ2_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&
        (dexp(-dexp(lambda)*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ2_uni
        
        
        !************************************************
        
             SUBROUTINE gauherJ3_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(1)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ3_uni
    
    
    
    !********************************************
    
         SUBROUTINE gauherJ4_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0-lambda-Xea22(1)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(1))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ4_uni
    
    !*******************************************
    
         SUBROUTINE gauherJ5_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ5_uni
    
    !*****************************************
    
         SUBROUTINE gauherJ6_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ6_uni
        
        !************************************************
        
             SUBROUTINE gauherJ7_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ7_uni
    
    
    !**********************************************
    
         SUBROUTINE gauherJ8_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ8_uni
    
    
    !****************************************************
    
         SUBROUTINE gauherJ9_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ9_uni
        
        
        !**********************************************
        
        
             SUBROUTINE gauherJ10_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(1)-lambda-Xea22(2)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(2))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ10_uni
    
    
    !********************************************
    
         SUBROUTINE gauherJ11_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(3)-lambda+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda)*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ11_uni
        
        !*******************************************
        
             SUBROUTINE gauherJ12_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ12_uni
    
    
    !***********************************************
    
         SUBROUTINE gauherJ13_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ13_uni
    
    
    !***************************************
    
         SUBROUTINE gauherJ14_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(3))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ14_uni
    
    
    !********************************************
    
         SUBROUTINE gauherJ15_uni(ss)
    
    use optim
    
        use tailles
        use donnees
        use comon,only:sigmae,nb1,& !nmesy,npp,nva3,nea
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt
            nodes,weights,yy,it,ziy,mat,det, K_G0, K_D0, lambda, ziy,y0,invBi_cholDet, &
            box_cox1,box_cox_par
        use donnees_indiv
         Implicit none
    
        double precision,intent(out)::ss
    
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision :: func72J
        double precision,parameter::pi=3.141592653589793d0
         
        upper = .false.
        ss=0.d0
            
    !        if(numpat.eq.2)write(*,*)'stad',nnodes_all
            
            do n=1,nnodes_all
       
            mu1 = 0.d0
           i = numpat
           
             matb_chol = 0.d0
          
          do k=1,nb1 
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do

    

        if(methodGH.eq.1) then
            Xea(1:nb1) = nodes(n,1:nb1)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nb1) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            else
            Xea2(1:nb1,1) = nodes(n,1:nb1)*sqrt(2.d0)
        Xea22(1:nb1) = nodes(n,1:nb1)*sqrt(2.d0)
            end if

        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

        !************ODE - analytic solution *********************

    
                    mu1(1:nmescur,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
            ziy((it+1):(it+nmescur),2)*dexp(K_D0+Xea22(3)-lambda-Xea22(4)+mu(1:nmescur,2))*&!-Xea22(3)+Xea22(2)
        (dexp(-dexp(lambda+Xea22(4))*ziy((it+1):(it+nmescur),1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(yy(it+k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            
                    else
                yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
            
        yscalar = yscalar + (yy(it+k)-mu1(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)


        func72J =dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)
                       
    
    
        func72J = dexp(func72J)

                    ss = ss+product(weights(n,1:nb1))*func72J
                    
    

            end do
    if(methodGH.eq.1) ss = ss*invBi_cholDet(i)*2.d0**(nb1/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nb1/2.d0)

        return

        END SUBROUTINE gauherJ15_uni
    