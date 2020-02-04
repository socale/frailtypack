    
    
    
    !========================          FUNCPA NEW         ====================
        double precision function funcpajlongiweib_nl(b,np,id,thi,jd,thj,k0)
    
        use lois_normales
        use donnees_indiv
        use tailles
        use comon
        !use ParametresPourParallelisation
        use comongroup,only:vet,vet2
        use residusM
            use optim
        implicit none
    
    ! *** NOUVELLLE DECLARATION F90 :
    
        integer,intent(in)::id,jd,np
        double precision,dimension(np),intent(in)::b
        double precision,dimension(2)::k0
        double precision,intent(in)::thi,thj
    
        integer::n,i,j,k,vj,ig,choix,ier,jj
        integer,dimension(ngmax)::cpt
        double precision::sum,res,eps
    
    double precision,dimension(np)::bh
        double precision,dimension(ngmax)::res2,res1dc,res2dc &
        ,res3dc,integrale1,integrale2,integrale3
        double precision::int
            ! for the numerical integral hrmsym
        integer :: restar,nf2
        double precision:: epsabs,epsrel
        !double precision,dimension(2):: result, abserr2
        !double precision,dimension(1000) :: work
        external :: vraistot,vraistot_weib
        double precision,dimension(nea) :: xea
        !integer ::neval,ifail
    
        !integer::l
    double precision :: eps_w,finddet
        double precision,dimension(ngmax)::integrale4
        double precision,parameter::pi=3.141592653589793d0
            !double precision,dimension(:,:),allocatable :: mat_sigma,varcov_marg_inv
                    !double precision,dimension(nva3,nva3) :: element
             double precision,dimension(nea*(nea+1)/2)::matv
    
  
        npp = np
        eps_w=1.d-7
        kkapa=k0
        choix=0
        ig=0
        k=0
        vj=0
        n=0
        j=0
        do i=1,np
            bh(i)=b(i)
        end do
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
        b1 =bh
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar-2)/(effet+1)!/nst
    
        if(typeJoint.eq.3) then
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= bh(3)**2
        etaD= bh(4)**2
        else
        betaD= bh(1)**2
        etaD= bh(2)**2
        end if
    
            K_G0 =  bh(np-nva-nb_re-1-netadc - netar-2*effet-1)
        K_D0 =  bh(np-nva-nb_re-1-netadc - netar-2*effet)
    
        if(typeJoint.eq.3) then
            sigmav = bh(np-nva-nb_re-1-netadc - netar-1)
            alpha = bh(np-nva-nb_re-1-netadc - netar)
        else
        sigmav = -1.d0
        end if
    
        if(nea.ge.1) then
                etaydc(1:netadc) = bh(np-nva-nb_re-netadc:np-nva-nb_re-1)
                 etayr(1:netar) =bh(np-nva-nb_re-netadc - netar:np-nva-nb_re-netadc - 1)
           sigmae = bh(np-nva-nb_re)*bh(np-nva-nb_re)

            Ut = 0.d0
            Utt = 0.d0

                  do j=1,nb1
                do k=1,j
                 
                       Ut(j,k)=sqrt(bh(np-nva-nb_re+k)**2.d0)
                  Utt(k,j)=sqrt(bh(np-nva-nb_re+k)**2.d0)
            
                end do
            end do
    
            if(typeJoint.eq.3) then
                Ut(nea,nea) = sqrt(sigmav**2.d0)
                Utt(nea,nea) = sqrt(sigmav**2.d0)
            end if
    
        end if
    
        
        
        mat = matmul(ut,utt)
    

        jj=0
    ! jjj = 0
        do j=1,nea
        do k=j,nea
        jj=j+k*(k-1)/2
    !    jjj = jjj +1
    
        matv(jj)=mat(j,k)
        !  bb2vv(jjj)=bb2(j,k)
    
        end do
        end do
        ier = 0
        eps = 1.d-10
    
 
        call dsinvj(matv,nea,eps,ier)
    
        mat=0.d0
        do j=1,nea
                do k=1,nea
                            if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                end do
    det = finddet(matmul(ut,utt),nea)
    
   
    !---- avec ou sans variable explicative  ------
        funcpajlongiweib_nl = 0.d0
        do k=1,ng
            res1(k) = 0.d0
            res3(k) = 0.d0
            res2(k) = 0.d0
            res1dc(k) = 0.d0
            res2dc(k) = 0.d0
            res3dc(k) = 0.d0
            cpt(k) = 0
            integrale1(k) = 0.d0
            integrale2(k) = 0.d0
            integrale3(k) = 0.d0
        !    integrale4(k) = 1.d0
            aux1(k)=0.d0
            aux2(k)=0.d0
        end do
    
    
    !ccccccccccccccccccccccccccccccccccccccccc
    !     pour les donnees recurrentes
    !ccccccccccccccccccccccccccccccccccccccccc
    
        if(typeJoint.eq.1.or.typeJoint.eq.3) then
    
        do i=1,nsujet
        !  cpt(g(i))=cpt(g(i))+1
            if(nva1.gt.0)then
                vet = 0.d0
                do j=1,nva1
                    vet =vet + bh(np-nva3-nva2-nva1+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif
    
    
            if((c(i).eq.1))then
    
    
                res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
    
                if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                !    funcpajlongiweib=-1.d9
                !    write(*,*)'ok 1',bh
                    goto 123
                end if
            endif
    
    !     nouvelle version
            res1(g(i)) = res1(g(i))+((t1(i)/etaR)**betaR)*vet
    
    
            res3(g(i)) = res3(g(i))+((t0(i)/etaR)**betaR)*vet
    
    
        end do
    
        end if
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    vet22(1:ng) = MATMUL(dble(vedc(1:ng,1:nva2)),bh(np-nva4-nva3-nva2+1:(np-nva3-nva4)))
    vet22 = dexp(vet22)
   
        do k=1,ng
        !   if(nva2.gt.0)then
        !       vet2 = 0.d0
        !       do j=1,nva2
        !
        !           vet2 =vet2 + bh(np-nva3-nva2+j)*dble(vedc(k,j))
        !       end do
        !       vet2 = dexp(vet2)
        !   else
        !       vet2=1.d0
        !   endif
            if(cdc(k).eq.1)then
                res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet22(k))
                !write(*,*)res2dc(k)
                if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                    funcpajLongiweib_nl=-1.d9
                            !       write(*,*)'ble',t1dc(k)
        !        print*,'ok 4'
                    goto 123
                end if
            endif
    
    !pour le calcul des integrales / pour la survie, pas les donnes recurrentes:
        aux1(k)=((t1dc(k)/etaD)**betaD)*vet22(k)
    
    
        RisqCumul(k) = aux1(k)
            Rdc_res(k) = aux1(k)/vet2
        end do
    
    !**************INTEGRALES ****************************
        
        res = 0.d0
    
            it = 0
            it_rec = 1
            epsabs = 1.d-100
            epsrel = 1.d-100
            restar = 0
            nf2 = nf
       !     nmes_o=0
            integrale4 = 0.d0
    
            do ig=1,ng
    
            !    ycurrent  = 0.d0
                auxig=ig
                choix = 4
                numpat = ig
                nmescur =nmesy(ig)
                nmescurr = nmesrec(ig)
                nmescurr1 = nmesrec1(ig)
                it_cur = it
    
                x2 = 0.d0
                x2cur = 0.d0
                z1cur = 0.d0
                current_mean = 0.d0
                          
         
         
    
                res1cur = 0.d0
                res2cur = 0.d0
                res3cur = 0.d0
    
       
    
        
              mu = 0.d0
            mu(1:nmescur,1) = matmul(vey((it+1):(it+nmescur),1:(nva3)),bh((np-nva4-nva3+1):(np-nva4)))
    
            mu(1:nmescur,2) = matmul(vey((it+1):(it+nmescur),(nva3+1):(nva3+nva4)),bh((np-nva4+1):np))
            xea = 0.d0
    
    
              
                    select case(which_random)
        case(1)
            call gauherJ1_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(2)
            call gauherJ2_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(3)
            call gauherJ3_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(4)
            call gauherJ4_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(5)
            call gauherJ5_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(6)
            call gauherJ6_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(7)
            call gauherJ7_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(8)
            call gauherJ8_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(9)
            call gauherJ9_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(10)
            call gauherJ10_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(11)
            call gauherJ11_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(12)
            call gauherJ12_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(13)
            call gauherJ13_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(14)
            call gauherJ14_nl(int,aux1(ig),res1(ig)-res3(ig))
        case(15)
            call gauherJ15_nl(int,aux1(ig),res1(ig)-res3(ig))        
    end select    
              
                integrale4(ig) =int !result(1) !
        
  
            it_rec = it_rec + nmescurr
            it = it + nmescur
    
            if(integrale4(ig).gt.1.E+30) then
                integrale4(ig) = 1.E+30
            end if
    
      
    
            sum=0.d0
          
    
            !********* vraisemblance ***************
                if (integrale4(ig).le.0.d0) then
                        res= res +res2(ig) &
                            + res2dc(ig)-0.5d0*dlog(2.d0*pi)*nmes_o(ig) -0.5d0*dlog(sigmae)*nmes_o(ig)  &
                            -112.d0
            else
                        res= res + res2(ig) &
                            + res2dc(ig)- 0.5d0*dlog(2.d0*pi)*nmes_o(ig) -0.5d0*dlog(sigmae)*nmes_o(ig)  &
                            +dlog(integrale4(ig))
            endif
    
    
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpajLongiweib_nl=-1.d9
                    goto 123
                end if
    
    
        end do
    !    write(*,*)res
    !--------------------------------------------------------
    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpajLongiweib_nl =-1.d9
        if(typeJoint.eq.3)   Rrec = 0.d0
        if(typeJoint.eq.3)   Nrec = 0.d0
            Rdc = 0.d0
            Ndc = 0.d0
            !print*,'ok 7'
            goto 123
        else
            funcpajLongiweib_nl = res
    
            do k=1,ng
                if(typeJoint.eq.3)Rrec(k)=res1(k)
                if(typeJoint.eq.3)Nrec(k)=nig(k)
                Rdc(k)=RisqCumul(k)
                Ndc(k)=cdc(k)
            end do
        end if
    
    !Ad:
        123     continue
    
        return
    
        end function funcpajlongiweib_nl
    
    !=================================================================================================
   
            
    !========================          FUNCPA NEW         ====================
        double precision function funcweib_nl(b,np,id,thi,jd,thj,k0)
    
        use lois_normales
        use donnees_indiv
        use tailles
        use comon
        !use ParametresPourParallelisation
        use comongroup,only:vet,vet2
        use residusM
            use optim
        implicit none
    
    ! *** NOUVELLLE DECLARATION F90 :
    
        integer,intent(in)::id,jd,np
        double precision,dimension(np),intent(in)::b
        double precision,dimension(2)::k0
        double precision,intent(in)::thi,thj
    
        integer::n,i,j,k,vj,ig,choix,ier,jj
        integer,dimension(ngmax)::cpt
        double precision::res,eps
    
    double precision,dimension(np)::bh
        double precision,dimension(ngmax)::res2,res1dc,res2dc &
        ,res3dc,integrale1,integrale2,integrale3
        !double precision::int
            ! for the numerical integral hrmsym
        !integer :: restar,nf2
        !double precision:: epsabs,epsrel
        !double precision,dimension(2):: result, abserr2
        !double precision,dimension(1000) :: work
        external :: vraistot,vraistot_weib
        !double precision,dimension(nea) :: xea
        !integer ::neval,ifail
    
        !integer::l
    double precision :: eps_w,finddet
        !double precision,dimension(ngmax)::integrale4
        double precision,parameter::pi=3.141592653589793d0
            !double precision,dimension(:,:),allocatable :: mat_sigma,varcov_marg_inv
                    !double precision,dimension(nva3,nva3) :: element
             double precision,dimension(nea*(nea+1)/2)::matv
    
  
        npp = np
        eps_w=1.d-7
        kkapa=k0
        choix=0
        ig=0
        k=0
        vj=0
        n=0
        j=0
        do i=1,np
            bh(i)=b(i)
        end do
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
        b1 =bh
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar-2)/(effet+1)!/nst
    
        if(typeJoint.eq.3) then
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= bh(3)**2
        etaD= bh(4)**2
        else
        betaD= bh(1)**2
        etaD= bh(2)**2
        end if
    
            K_G0 =  bh(np-nva-nb_re-1-netadc - netar-2*effet-1)
        K_D0 =  bh(np-nva-nb_re-1-netadc - netar-2*effet)
    
        if(typeJoint.eq.3) then
            sigmav = bh(np-nva-nb_re-1-netadc - netar-1)
            alpha = bh(np-nva-nb_re-1-netadc - netar)
        else
        sigmav = -1.d0
        end if
    
        if(nea.ge.1) then
                etaydc(1:netadc) = bh(np-nva-nb_re-netadc:np-nva-nb_re-1)
                 etayr(1:netar) =bh(np-nva-nb_re-netadc - netar:np-nva-nb_re-netadc - 1)
           sigmae = bh(np-nva-nb_re)*bh(np-nva-nb_re)

            Ut = 0.d0
            Utt = 0.d0

                  do j=1,nb1
                do k=1,j
                    if(j.eq.k) then
                    
                       Ut(j,k)=sqrt(bh(np-nva-nb_re+k)**2.d0)
                  Utt(k,j)=sqrt(bh(np-nva-nb_re+k)**2.d0)
            !      Ut(j,k)=sqrt(bh(np-nva-nb_re+k+j*(j-1)/2)**2.d0)
            !     Utt(k,j)=sqrt(bh(np-nva-nb_re+k+j*(j-1)/2)**2.d0)
                  else 
             !           Ut(j,k)=bh(np-nva-nb_re+k+j*(j-1)/2)
              !      Utt(k,j)=bh(np-nva-nb_re+k+j*(j-1)/2)
                end if
                end do
            end do
    
            if(typeJoint.eq.3) then
                Ut(nea,nea) = sqrt(sigmav**2.d0)
                Utt(nea,nea) = sqrt(sigmav**2.d0)
            end if
    
        end if
    
        
        
        mat = matmul(ut,utt)
    

        jj=0
    ! jjj = 0
        do j=1,nea
        do k=j,nea
        jj=j+k*(k-1)/2
    !    jjj = jjj +1
    
        matv(jj)=mat(j,k)
        !  bb2vv(jjj)=bb2(j,k)
    
        end do
        end do
        ier = 0
        eps = 1.d-10
    
 
        call dsinvj(matv,nea,eps,ier)
    
        mat=0.d0
        do j=1,nea
                do k=1,nea
                            if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                end do
    det = finddet(matmul(ut,utt),nea)
    
   
    !---- avec ou sans variable explicative  ------

        do k=1,ng
            res1(k) = 0.d0
            res3(k) = 0.d0
            res2(k) = 0.d0
            res1dc(k) = 0.d0
            res2dc(k) = 0.d0
            res3dc(k) = 0.d0
            cpt(k) = 0
            integrale1(k) = 0.d0
            integrale2(k) = 0.d0
            integrale3(k) = 0.d0
        !    integrale4(k) = 1.d0
            aux1(k)=0.d0
            aux2(k)=0.d0
        end do
    
    
    !ccccccccccccccccccccccccccccccccccccccccc
    !     pour les donnees recurrentes
    !ccccccccccccccccccccccccccccccccccccccccc
    
        if(typeJoint.eq.1.or.typeJoint.eq.3) then
    
        do i=1,nsujet
        !  cpt(g(i))=cpt(g(i))+1
            if(nva1.gt.0)then
                vet = 0.d0
                do j=1,nva1
                    vet =vet + bh(np-nva3-nva2-nva1+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif
    
    
            if((c(i).eq.1))then
    
    
                res2(g(i)) = res2(g(i))+(betaR-1.d0)*dlog(t1(i))+dlog(betaR)-betaR*dlog(etaR)+dlog(vet)
    
                if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                !    funcpajlongiweib=-1.d9
                !    write(*,*)'ok 1',bh
                    goto 123
                end if
            endif
    
    !     nouvelle version
            res1(g(i)) = res1(g(i))+((t1(i)/etaR)**betaR)*vet
    
    
            res3(g(i)) = res3(g(i))+((t0(i)/etaR)**betaR)*vet
    
    
        end do
    
        end if
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    vet22(1:ng) = MATMUL(dble(vedc(1:ng,1:nva2)),bh(np-nva4-nva3-nva2+1:(np-nva3-nva4)))
    vet22 = dexp(vet22)
   
        do k=1,ng
   
    
    !pour le calcul des integrales / pour la survie, pas les donnes recurrentes:
        aux1(k)=((t1dc(k)/etaD)**betaD)*vet22(k)
    
    
        RisqCumul(k) = aux1(k)*vet22(k)
            Rdc_res(k) = aux1(k)/vet2
        end do
    
        
        res = 0.d0
    
    
    !--------------------------------------------------------
    
   
            funcweib_nl = res
    
            do k=1,ng
                if(typeJoint.eq.3)Rrec(k)=res1(k)-res3(k)
                if(typeJoint.eq.3)Nrec(k)=nig(k)
                Rdc(k)=RisqCumul(k)
                Ndc(k)=cdc(k)
            end do
       
    
    !Ad:
        123     continue
    
        return
    
        end function funcweib_nl
    