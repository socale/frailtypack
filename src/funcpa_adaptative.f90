module func_adaptative

    implicit none
    
    contains
   ! fonction a maximiser pour l'estimation des effes aleatoires pour un sujet
    double precision function funcpafrailtyPred_ind(b,np,id,thi,jd,thj,k0,individu_j)
        
        use var_surrogate, only:vs_i,vt_i,u_i,theta2,const_res5,const_res4,&
            deltastar,delta,pi,& !varcovinv,penalisation,essai_courant,cdcts,nigts
            alpha_ui,frailt_base
        use comon, only: eta,ve !lognormal,resnonpen
        
        implicit none
         
        integer,intent(in)::id,jd,np,individu_j
        !integer::i
        double precision,dimension(np),intent(in)::b 
        double precision,dimension(2),intent(in)::k0
        double precision,intent(in)::thi,thj
        double precision::vsi,vti,res,ui !test
        double precision,dimension(:),allocatable::bh
        double precision::wij
        
        allocate(bh(np))
        !call dblepr("b(1) funcpafrailtyPred_ind 41", -1, b(1), 1)             
        bh(1)=b(1)
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
        
        wij=bh(1)
        ! on fixe vsi et vti à leurs valeurs initiles donnees par le premier point de quadrature
        vsi=vs_i
        vti=vt_i    
        ui=u_i
        
        !call intpr("funcpafrailtyPred_ind 60", -1, np, 1)
        if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de l'heterogeneite sur les risque des bas
            res=dexp((vsi*delta(individu_j)+vti*deltastar(individu_j))*dble(ve(individu_j,1))&
            -(wij**2.d0)/(2.d0*theta2)&
                + wij*(delta(individu_j)+deltastar(individu_j)*eta)&
                - const_res4(individu_j)*dexp(wij+vsi*dble(ve(individu_j,1)))&
                - const_res5(individu_j)*dexp(eta*wij+vti*dble(ve(individu_j,1)))&
            )
        else
            res=dexp(ui*(delta(individu_j)+deltastar(individu_j)*alpha_ui)&
                +(vsi*delta(individu_j)+vti*deltastar(individu_j))*dble(ve(individu_j,1))&
                -(wij**2.d0)/(2.d0*theta2)&
                + wij*(delta(individu_j)+deltastar(individu_j)*eta)&
                - const_res4(individu_j)*dexp(wij+ui+vsi*dble(ve(individu_j,1)))&
                - const_res5(individu_j)*dexp(eta*wij+alpha_ui*ui+vti*dble(ve(individu_j,1)))&
            )
        endif
        res=dlog(res)
        
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpafrailtyPred_ind =-1.d9
        goto 123
    else
        funcpafrailtyPred_ind = res
    end if
    
    
    123     continue
    deallocate(bh)
    !call dblepr("funcpa95", -1, funcpafrailtyPred_ind, 1)
    return
    
    endfunction funcpafrailtyPred_ind
    
    ! fonction a maximiser pour l'estimation des effes aleatoires pour un essai
    double precision function funcpafrailtyPred_Essai(b,np,id,thi,jd,thj,k0)
        
        use var_surrogate, only:vs_i,vt_i,u_i,theta2,& !const_res5,const_res4,
            pi,varcovinv,cdcts,nigts,essai_courant,ui_chap,& !deltastar,delta
            nsujeti,npoint,estim_wij_chap,& !adaptative,penalisation,posind_i
            indicej,invBi_chol_Individuel,individu_j,nparamfrail,&
            gamma_ui,frailt_base,methodInt,nsim
        use comon, only: invBi_cholDet !eta,lognormal,ve,resnonpen
        use fonction_A_integrer,only:Integrale_Individuel,Integrale_Individuel_MC
        use optim_scl, only:marq98j_scl  ! pour faire appel a marquard 
        
        implicit none
         
        integer,intent(in)::id,jd,np
        integer::i
        double precision,dimension(np),intent(in)::b
        double precision,dimension(2),intent(in)::k0
        double precision,intent(in)::thi,thj
        double precision::vsi,vti,res,ui
        double precision,dimension(np)::bh
        !double precision::wij
        double precision ::I1,c1,c2,integral
        double precision, dimension(:,:),allocatable::m1,m3  
        double precision, dimension(:,:),allocatable::m
        
        ! var utilisees en vue de l'estimation des fragilites a posteriorie
        integer,parameter::effet2=0,np_1=1
        double precision::ca,cb,dd
        !double precision::res
        double precision, dimension(2)::k0_2
        double precision, dimension(:),allocatable::v,b_2
        double precision, allocatable, dimension(:,:)::H_hessOut,HIH,HIHOut,IH,invBi_chol_2,H_hess_scl,I_hess_scl
        double precision,dimension(:,:), allocatable::hess_scl
        double precision,dimension(:), allocatable::vvv_scl
        integer::ier,istop,ss,sss,ni,nmax_2,np_2,nparamfrail_save
        
                        
        do i=1,np
            bh(i)=b(i)
        end do
        
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
        
        vsi=bh(1)
        vti=bh(2)
        if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de
            ui=0.d0
        else
            ui=bh(3)
        endif

        if(methodInt==4) goto 125 ! dans ce cas on n'a plus besoin destimer les w_ij_chapeau car on integre au niveau individuel par MC
        !===========estimation des wij_chapeaux associes au vsi=======================
        !=============================================================================
        !            i=1
                    ! posind_i=1 
                    ! do k=1,ntrials !k permet d'indicer les essais. ce changement empeche l'ambiguite avec le J passe en parametre
        !                indicej=i
                        nmax_2=sum(nsujeti(1:essai_courant)) !pour la somme cumulee du nombre de sujet par essai
                        !frail_essai_deja_est=0 ! l'on n'a pas encore estime les vsi et vti pour cet essai
        !                essai_courant=k
        ! !print*,nmax_2,indicej,i
        ! stop
                        do i=indicej,nmax_2
                            !b(1)=-20.d0
                            individu_j=i
                            
                            !====================================================================================================
                            ! estimation des w_ij
                            !====================================================================================================    
                            k0_2=k0
                            vs_i=vsi
                            vt_i=vti
                            u_i=ui
                            ni=0
                            ca=0.d0
                            cb=0.d0
                            dd=0.d0
                            
                            np_2=1
                            allocate(I_hess_scl(np_2,np_2),v(np_2*(np_2+3)/2),b_2(1))
                            allocate(H_hess_scl(np_2,np_2),invBi_chol_2(np_2,np_2),hess_scl(np_2,np_2),vvv_scl(np_2*(np_2+1)/2))
                            allocate(H_hessOut(np_2,np_2),HIH(np_2,np_2),HIHOut(np_2,np_2),IH(np_2,np_2))
                    
                            b_2(1)=0.5d0
                            v=0.d0
                            nparamfrail_save=nparamfrail
                            nparamfrail=1
                            !!print*,"ok pour le premier"
                            call marq98J_scl(k0_2,b_2,np_1,ni,v,res,ier,istop,effet2,ca,cb,dd,funcpafrailtyPred_ind,&
                                             I_hess_scl,H_hess_scl,hess_scl,vvv_scl,individu_j)
                            nparamfrail=nparamfrail_save ! on restitu sa valeur avant de continuer
                            
                            !!print*,"Fin estimation de w_ij_chapeau pour l'individu",i,"istop=",b_2
                            ! if(i==2)then
                                ! !print*,"Fin estimation",i,"istop=",b_2,vs_i,vt_i
                            ! endif
                            ! stop
                            !!print*,"b_2(1)=",b_2(1),"ui_chap(i,1)=",ui_chap(i,1)
                            !stop                            
                            ui_chap(i,1)=b_2(1)
                            !!print*,"ok pour le premier"
                            if (istop.ne.1) then
                                !print*,"2-individu",i,"wij=",b_2,"istop=",istop,"ier=",ier,"v=",v
                                !print*,"le modele d'estimation des w_ij_chapeau n'a pas converge, i=",i
                                funcpafrailtyPred_Essai =-1.d9
                                goto 124
                            endif
                            
                            do ss=1,np_1
                                do sss=1,np_1
                                    !HIHOut(ss,sss) = HIH(ss,sss)
                                    H_hessOut(ss,sss)= I_hess_scl(ss,sss)
                                end do
                            end do
                            
                            ! !print*,H_hess_scl(1,1),I_hess_scl(1,1)
                            ! stop
                                
                            invBi_chol_Individuel(i)=dsqrt(H_hess_scl(1,1))
                            !calcul du determinant de la cholesky de l'inverse de la hessienne                    
                            invBi_cholDet(i)=invBi_chol_Individuel(i) !individuel
                            deallocate(H_hess_scl,I_hess_scl,H_hessOut,HIH,HIHOut,IH,invBi_chol_2,hess_scl,vvv_scl,v,b_2)
                        enddo ! fin estimation des w_ij_chapeau
                        
                        
                        
                        ! !print*,"impression des frailties au individuel dans le fichier Prediction_wij_chapeau"
                        ! !print*,"sujet"," ","w_ij_chapeau"
                        ! do ss=1,(i-1)
                            ! !print*,ss,ui_chap(ss,1)
                        ! enddo    
                    
                        ! ====================================================================================================
                        ! Fin estimation des ws_ij_chapeau et wt_ij_chapeau
                        ! ====================================================================================================
        125    continue
        
        ! calcul du terme I1(vsi,vti)
        !!print*,"essai_courant=",essai_courant
        allocate(m1(1,2),m3(1,2),m(1,1))
        m1(1,1)=vsi
        m1(1,2)=vti
        m3 = matmul(m1,varcovinv) !produit matriciel, resultat dans m3
        m = matmul(m3,transpose(m1)) !produit matriciel, resultat dans m
        c1=(-1.d0/2.d0) * m(1,1)
        c2=nigts(essai_courant)*vsi+cdcts(essai_courant)*vti ! car deja pris en compte dans "Integrale_Individuel()"
        c2=0.d0 ! car deja pris en compte dans "Integrale_Individuel()"
        
        if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de
            I1=dexp(c1+c2)
        else
            I1=dexp(c1+c2-0.5d0*ui**2/gamma_ui)
        endif
        ! !print*,"I1=",I1
        ! !print*,"m1=",m1
        ! !print*,"varcovinv=",varcovinv
        ! !print*,"m3=",m3
        ! !print*,"m=",m
        ! !print*,"c1=",c1
        ! !print*,"c2=",c2
        res=1.d0
        estim_wij_chap=1    
        !adaptative=.false.
        if(methodInt==4)then !integration au niveau individuel par monte-carlo
            do individu_j=1,nsujeti(essai_courant)
                res=res*Integrale_Individuel_MC(vsi,vti,ui,individu_j,nsim,0.d0,theta2)
            enddo        
        else
            do individu_j=1,nsujeti(essai_courant)
                ! res=res*Integrale_Individuel(vsi,vti,ui,individu_j,npoint)
                res=res*Integrale_Individuel(vsi,vti,ui,individu_j,npoint)
            enddo
        endif
        
        !adaptative=.true.
        res=I1*res
        res=dlog(res)
        
        !resnonpen = res
        !res = res - penalisation
        deallocate(m1,m3,m)
        ! !print*,"res=",res
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpafrailtyPred_Essai =-1.d9
        goto 124
    else
        funcpafrailtyPred_Essai = res
    end if
    
    
    124    continue
    
    return
    
    endfunction funcpafrailtyPred_Essai
    
    ! joint frailty_copula
    double precision function funcpafrailtyPred_copula(b,np,id,thi,jd,thj,k0)
        use var_surrogate, only:essai_courant,nsujeti,frailt_base
        use fonction_A_integrer,only:Integrant_Copula
        
        implicit none
         
        integer,intent(in)::id,jd,np
        integer::i
        double precision,dimension(np),intent(in)::b
        double precision,dimension(2),intent(in)::k0
        double precision,intent(in)::thi,thj
        double precision::vsi,vti,res,ui
        double precision,dimension(np)::bh
         
        do i=1,np
            bh(i)=b(i)
        end do
        
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
        
        vsi=bh(1)
        vti=bh(2)
        if(frailt_base==0) then! on annule simplement le terme avec ui si on ne doit pas tenir compte de
            ui=0.d0
        else
            ui=bh(3)
        endif
        
        res = Integrant_Copula(vsi,vti,ui,essai_courant,nsujeti(essai_courant))
        !call dblepr("res funcpaadaptativ = ", -1, res, 1)
        res = dlog(res)
        !call dblepr("log res funcpaadaptativ = ", -1, res, 1)
        
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then ! pour test infini et NaN
            funcpafrailtyPred_copula =-1.d9
            !call dblepr("log res funcpaadaptativ = ", -1, res, 1)
            goto 124
        else
            funcpafrailtyPred_copula = res
        end if
    
    
        124    continue
        
        return
    
    endfunction funcpafrailtyPred_copula

endmodule func_adaptative