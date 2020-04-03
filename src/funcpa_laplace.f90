module func_laplace

    implicit none
    
    contains
   ! fonction a maximiser pour recherche les solution de k(w_ij) niveau individuel
    double precision function funcpaw_ij_chapeau(b,np,id,thi,jd,thj,k0,individu_j)
        
        use var_surrogate, only:vs_i,vt_i,u_i,theta2,const_res5,const_res4,&
            deltastar,delta,pi,res2s_sujet,res2_dcs_sujet,alpha_ui,Test
        use comon, only: eta,ve
        
        implicit none
         
        integer,intent(in)::id,jd,np,individu_j
        double precision,dimension(np),intent(in)::b
        double precision,dimension(2),intent(in)::k0
        double precision,intent(in)::thi,thj
        double precision::vsi,vti,res,ui
        double precision,dimension(np)::bh
        double precision::wij,zeta,A_Lap
                        
        bh(1)=b(1)
        
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
        
        ! !print*,id,jd
        ! stop
        
        wij=bh(1)
        ! on recupere les effets aleatoires au niveau essai, a modifier a chaque itteration
        vsi=vs_i
        vti=vt_i    
        ui=u_i
        zeta=eta
        
        
        ! calcul de A
        ! A_Lap= delta(individu_j)*res2s_sujet(individu_j)+deltastar(individu_j)*res2_dcs_sujet(individu_j)&
            ! +ui*(delta(individu_j)+deltastar(individu_j)*alpha_ui)&
            ! +(vsi*delta(individu_j)+vti*deltastar(individu_j))*dble(ve(individu_j,1))&
            ! -dlog(2.d0*pi*theta2)/2.d0
        A_Lap= ui*(delta(individu_j)+deltastar(individu_j)*alpha_ui)&
            +(vsi*delta(individu_j)+vti*deltastar(individu_j))*dble(ve(individu_j,1))
            
        ! calcul de M
!        M_Lap=(A_Lap+1.d0)**2.d0
        
        ! calcul de K(w_ij)
        
        res = A_Lap + wij*(delta(individu_j)+deltastar(individu_j)*zeta) - wij**2.d0/(2.d0*theta2)&
            - const_res4(individu_j)*dexp(wij + ui + vsi*dble(ve(individu_j,1)))&
            - const_res5(individu_j)*dexp(zeta*wij + alpha_ui*ui + vti*dble(ve(individu_j,1)))
            
        if(Test==1)then ! je fais ceci pour evaluer le calcul integral par laplace
            res=-bh(1) + 5.d0*dlog(bh(1))
        endif
        ! !print*,"res",res
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaw_ij_chapeau =-1.d9
        goto 123
    else
        !funcpaw_ij_chapeau = -bh(1) + 5.d0*dlog(bh(1)) ! pour tester la maximisation avec Marquard. marche tres bien
        funcpaw_ij_chapeau = res
    end if
    
    123     continue
    
    return
    
    endfunction funcpaw_ij_chapeau

    ! Calcul integrale au niveau individuel par approximation de Laplace
    
    double precision function Int_Laplace_ind(position_i,individu_essai,vsi,vti,ui)
        !individu_essai: individu courant dans l'essai
        !vsi,vti,ui : effets aleatoires au niveau essai
        !wij_chap : contient des w_ij_chapeau
        use var_surrogate, only:theta2,const_res5,const_res4,vs_i,vt_i,u_i,wij_chap,& !nsujeti
            deltastar,delta,pi,alpha_ui,Test!,individu_j !res2s_sujet,res2_dcs_sujet,control_wij_chap
        use comon, only: eta,ve,model
        use optim_scl, only:marq98j_scl  ! pour faire appel a marquard 
        
        implicit none
         
        integer,intent(in)::position_i,individu_essai
        integer::model_save,individu !i

        double precision,intent(in)::vsi,vti,ui
        double precision::res
        
        ! var utilisees en vue de l'estimation des fragilites a posteriorie
        integer,parameter::effet2=0
        double precision::ca,cb,dd,k_second,zeta,A_Lap,integr
        !double precision::res
        !double precision,dimension(:,:),intent(inout)::wij_chap
        double precision, dimension(2)::k0_2
        double precision, dimension(:),allocatable::v,b_2
        double precision, allocatable, dimension(:,:)::H_hess_scl,I_hess_scl
        double precision,dimension(:,:), allocatable::hess_scl
        double precision,dimension(:), allocatable::vvv_scl
        integer::ier,istop,ni,np_2,non_conv
        
        zeta=eta
        individu=position_i-1+individu_essai
        !individu_j=individu
        !!print*,"individu_j=",individu_j
    ! if(control_wij_chap(individu_j)==0)then ! on estime une seule fois wij_chap        
        !====================================================================================================
        ! ==============================estimation des w_ij_chapeau==========================================
        !====================================================================================================    
        k0_2=100
        ni=0
        ca=0.d0
        cb=0.d0
        dd=0.d0
        vs_i=vsi
        vt_i=vti
        u_i=ui
                                
        np_2=1
        allocate(v(np_2*(np_2+3)/2))
        allocate(b_2(np_2))
        allocate(I_hess_scl(np_2,np_2))
        allocate(H_hess_scl(np_2,np_2),hess_scl(np_2,np_2),vvv_scl(np_2*(np_2+1)/2))

                        
        b_2(1)=0.5d0
        v=0.d0 ! matrice des derrivees premieres et seconde
        
        ! !print*,"b avant",b_2
        model_save=model
        model = 9
        !nparamfrail_save=nparamfrail !jenleve car pas manipuler etant donné le type de modele
        !nparamfrail=1
        non_conv=0
        10 continue
        ! call marq98o(b_2,np_2,ni,v,res,ier,istop,funcpaw_ij_chapeau)
        call marq98J_scl(k0_2,b_2,np_2,ni,v,res,ier,istop,effet2,ca,cb,dd,funcpaw_ij_chapeau,I_hess_scl,H_hess_scl,&
                         hess_scl,vvv_scl,individu)
        
        !=============juste pour test=================
        if(Test==1)then ! je fais ceci pour evaluer le calcul integral par laplace
            if(istop==1)then
                ! !print*,"istop",istop,"x0 vaut:",b_2
                integr=dexp(-b_2(1) + 5.d0*dlog(b_2(1)))*dsqrt(2*pi*b_2(1)**2.d0/5.d0)
                ! !print*,"la valeur de l'integrale vaut",integr
                Int_Laplace_ind=integr
                wij_chap(1,1)=b_2(1)
                goto 124
                ! stop
            else
                !print*,"probleme de convergence individuel"
                ! stop
            endif
        endif
        !==============fin pour test=================
        
        if (istop.ne.1 .and. non_conv<=10) then ! on passe à l'individu suivant, juste pour le test
            !!print*,"Fin estimation_b essai bloquer",b_2,essai_courant
            b_2=-0.5*non_conv
            non_conv=non_conv+1 !compte le nombre de fois qu'on n'a pas pu estime les frailties niveau essai sur certains individus
            !stop
            goto 10
        endif
                                
        if(non_conv==11 .and. istop .ne. 1)then
            ! !print*,"le nombre de tentative sans convergence vaut:",non_conv,"sujet k=",individu
            non_conv=0
            Int_Laplace_ind =-1.d9
            goto 124
        endif
                                
        if(non_conv>0 .and. non_conv<=10) then
            !!print*,"le nombre de tentative pour la convergence vaut:",non_conv
            ! !print*,"istop=",istop,"essai k=",essai_courant
            non_conv=0
            !stop
        endif    
        
        model=model_save
        ! !print*,    "sizeof(wij_chap)=",size(wij_chap,1)    
        ! !print*,"wij_chap",wij_chap
        if(istop .eq. 1) then
            wij_chap(individu_essai,1)=b_2(1)
        endif
        
        ! ===========================fin estimation des w_ij_chapeau===========================================================
    ! endif
! !print*,"suis la 1"
    ! !print*,size(ve,1),size(ve,2),sum(ve)
    ! !print*,ve
    ! stop
        
        ! A_Lap= delta(individu)*res2s_sujet(individu)+deltastar(individu)*res2_dcs_sujet(individu)&
            ! +ui*(delta(individu)+deltastar(individu)*alpha_ui)&
            ! +(vsi*delta(individu)+vti*deltastar(individu))*dble(ve(individu,1))&
            ! -dlog(2.d0*pi*theta2)/2.d0
        
        A_Lap= ui*(delta(individu)+deltastar(individu)*alpha_ui)&
            +(vsi*delta(individu)+vti*deltastar(individu))*dble(ve(individu,1))
        
        ! calcul de la derivee seconde
        ! k_second=-(1.d0/theta2 + const_res4(individu)*dexp(wij_chap(individu_essai,1)+ui+vsi*dble(ve(individu,1)))&
            ! + zeta**2.d0 * const_res5(individu)*dexp(zeta*wij_chap(individu_essai,1)+alpha_ui*ui+vti*dble(ve(individu,1))))
        
        k_second=I_hess_scl(1,1) ! il s'agit de la hessienne
        ! !print*,"k_second",k_second,"hessienn",I_hess_scl
        ! stop
        
        ! calcul de la log integrale au niveau individuel
        res= (1.d0/2.d0)*dlog(2.d0*pi) + A_Lap + wij_chap(individu_essai,1)*(delta(individu)+deltastar(individu)*zeta)&
            - wij_chap(individu_essai,1)**2.d0/(2.d0*theta2) - const_res4(individu)*dexp(wij_chap(individu_essai,1)+ui&
            + vsi*dble(ve(individu,1)))&
            - const_res5(individu)*dexp(zeta*wij_chap(individu_essai,1)+alpha_ui*ui+vti*dble(ve(individu,1)))&
            -(1.d0/2.d0) * dlog(abs(k_second)) ! on fait * M pour annuler le 1/M utiliser dans le calcul de k_second
        
        ! !print*,"A=",A_Lap
        ! !print*,"M=",M_Lap
         !!print*,"k_second=",k_second,v,M_Lap
         !stop
        ! !print*,"Mk_second",k_second,"v",v
        ! !print*,"model",model
        ! !print*,"Integrale sur le sujet",individu_j,":",res
         !stop
        
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        Int_Laplace_ind =-1.d9
        goto 124
    else
        Int_Laplace_ind = res
    end if
    
    
    124    continue
    deallocate(v,b_2)
    deallocate(H_hess_scl,I_hess_scl,hess_scl,vvv_scl)
    return
    
    endfunction Int_Laplace_ind

    ! fonction h(X_i) a maximiser pour le calcul integrale au niveau essai
    
    double precision function funcpaXi_chapeau(b,np,id,thi,jd,thj,k0)
        !wij_chap: contient les valeur estimees de w_ij_chapeau
        use var_surrogate, only:pi,essai_courant,position_i,nsujeti,&
            gamma_ui,varcov,rho,Test !determinant,wij_chap
        !use comon, only: ve
        !use optim_scl, only:marq98j_scl  ! pour faire appel a marquard 
        !use fonction_A_integrer, only:Int_Laplace_ind
        
        implicit none
         
        integer,intent(in)::id,jd,np
        integer::i
        !double precision,dimension(:,:),intent(inout)::wij_chap
        double precision,dimension(np),intent(in)::b
        double precision,dimension(2),intent(in)::k0
        double precision,intent(in)::thi,thj
        double precision::res,vs_i,vt_i,u_i,h,h2,B_Lap,control !h1
        double precision,dimension(np)::bh

        
        !!print*,"vs_i=",vs_i
        !!print*,"vt_i=",vt_i
        !!print*,"individu_j=",individu_j
        !!print*,"position_i=",position_i
        
        !!print*,"b_i=",b
        !!print*,"np=",np
                        
        ! !print*,"suis la 3333",size(wij_chap),size(wij_chap,1),size(wij_chap,2)
        ! stop
        do i=1,np
            bh(i)=b(i)
        end do
        
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
        
        ! wij=bh(1)
        ! if(np==1) then
            ! vsi=vs_i
            ! vti=vt_i
        ! else
        u_i=bh(1)
        vs_i=bh(2)
        vt_i=bh(3)
        
        !=================pour le test=============================
        if(Test==1)then ! je fais ceci pour evaluer le calcul integral par laplace
            ! calcul integrale k*I(x): posons k=10
            
            h=10*Int_Laplace_ind(position_i,i,vs_i,vt_i,u_i)
            ! !print*,"voila la valeur de 10 *I(x)",h
            funcpaXi_chapeau=(-bh(1)**2.d0 - 2.d0*bh(2)**2.d0 -dlog(h/10))
            
            goto 125
            
            ! stop
        endif
        !=================fin pour le test==========================
        
        ! calcul des integrales au niveau individuel pour tous les sujets du cluster
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
            funcpaXi_chapeau =-1.d9
            goto 125
        endif
    
        ! calcul des parametres
        ! B_Lap=dlog(2.d0*pi)+(1.d0/2.d0)*(dlog(determinant)+dlog(2.d0*pi*gamma_ui))
        B_Lap=0.d0
        
        h=B_Lap &
            + (u_i**2.d0)/(2.d0*gamma_ui) + 1.d0/(2.d0*(1.d0-(rho**2.d0))) &
            * ((vs_i**2.d0)/varcov(1,1) + (vt_i**2.d0)/varcov(2,2)&
            - (2.d0*vs_i*vt_i*rho)/dsqrt(varcov(1,1)*varcov(2,2)))&
            - h2
        
        res=h    
        ! res=B_Lap+h
        ! !print*,"fin estimation des wij_chapea"
        ! !print*,wij_chap(1:nsujeti(essai_courant),1)
        !!print*,"vraisemblance=",res
        ! stop
    
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpaXi_chapeau =-1.d9
        goto 125
    else
        funcpaXi_chapeau = -res
    end if
    
    
    125    continue
    
    return
    endfunction funcpaXi_chapeau
    
    ! joint frailty_copula
    double precision function funcpaLaplace_copula(b,np,id,thi,jd,thj,k0)
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
        if(res == 0.d0) res = 1.d-299
        ! call dblepr("res funcpaadaptativ = ", -1, res, 1)
        res = -(- dlog(res))
        ! call dblepr("log res funcpaadaptativ = ", -1, res, 1)
        
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then ! pour test infini et NaN
            funcpaLaplace_copula =-1.d9
            !call dblepr("log res funcpaadaptativ = ", -1, res, 1)
            goto 126
        else
            funcpaLaplace_copula = res
        end if
    
        126    continue
        
        !funcpaLaplace_copula = -((bh(1)+1.d0) ** 2.d0 + (bh(2)+2.d0) ** 2.d0 + (bh(3)-6.d0) ** 2.d0 )
        return
    
    endfunction funcpaLaplace_copula
    
endmodule func_laplace
    