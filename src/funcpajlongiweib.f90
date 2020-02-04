    
    
    
    !========================          FUNCPA NEW         ====================
        double precision function funcpajlongiweib(b,np,id,thi,jd,thj,k0)
    
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
        double precision,dimension(2):: result, abserr2
        double precision,dimension(1000) :: work
        external :: vraistot,vraistot_weib
        double precision,dimension(nea) :: xea
        integer ::neval,ifail
    
        integer::l!,it
    double precision :: eps_w,finddet
        double precision,dimension(ngmax)::integrale4
        double precision,parameter::pi=3.141592653589793d0
            double precision,dimension(:,:),allocatable :: mat_sigma,varcov_marg_inv
            double precision,dimension(:),allocatable :: matv
            double precision,dimension(nva3,nva3) :: element
    
    
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
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar)/(effet+1)!/nst
    
        if(typeJoint.eq.3) then
        betaR= bh(1)**2
        etaR= bh(2)**2
        betaD= bh(3)**2
        etaD= bh(4)**2
        else
        betaD= bh(1)**2
        etaD= bh(2)**2
        end if
    
    
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
                    Ut(j,k)=sqrt(bh(np-nva-nb_re+k+j*(j-1)/2)**2.d0)
                    Utt(k,j)=sqrt(bh(np-nva-nb_re+k+j*(j-1)/2)**2.d0)
                    else 
                        Ut(j,k)=bh(np-nva-nb_re+k+j*(j-1)/2)
                    Utt(k,j)=bh(np-nva-nb_re+k+j*(j-1)/2)
                end if
                    
                end do
        
                end do
    
            if(typeJoint.eq.3) then
                Ut(nea,nea) = sigmav!sqrt(sigmav**2.d0)
                Utt(nea,nea) = sigmav!sqrt(sigmav**2.d0)
            end if
    
        end if
    
    
    !---- avec ou sans variable explicative  ------
        funcpajlongiweib = 0.d0
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
    
        do k=1,ng
            if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
    
                    vet2 =vet2 + bh(np-nva3-nva2+j)*dble(vedc(k,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
            if(cdc(k).eq.1)then
                res2dc(k) = (betaD-1.d0)*dlog(t1dc(k))+dlog(betaD)-betaD*dlog(etaD)+dlog(vet2)
                !write(*,*)res2dc(k)
                if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                    funcpajLongiweib=-1.d9
                            !       write(*,*)'ble',t1dc(k)
        !        print*,'ok 4'
                    goto 123
                end if
            endif
    
    !pour le calcul des integrales / pour la survie, pas les donnes recurrentes:
        aux1(k)=((t1dc(k)/etaD)**betaD)*vet2
    
    
        RisqCumul(k) = aux1(k)
            Rdc_res(k) = aux1(k)/vet2
        end do
    
    !**************INTEGRALES ****************************
    
    
        if(nea.ge.1) then
    
        it = 0
        it_rec = 1
        nmes_o=0
        EPSABS=1.d-50
        EPSREL=1.d-50
        restar = 0
        nf2 = nf
        integrale4 = 0.d0
    
            sum_mat= 0.d0
    
        do ig=1,ng
                    element = 0.d0
            ycurrent  = 0.d0
            auxig=ig
            choix = 4
            numpat = ig
            nmescur =nmesy(ig)
            nmescurr = nmesrec(ig)
            nmescurr1 = nmesrec1(ig)
            it_cur = it
            !call gaulagJ(int,choix)
    
    
                    allocate(mat_sigma(nmescur,nmescur))
    
            x2 = 0.d0
            x2cur = 0.d0
            z1cur = 0.d0
    
            current_mean = 0.d0
                    mat_sigma = 0.d0
                if(nmescur.gt.0) then
                    do i= 1,nmescur
                        ycurrent(i) = yy(it+i)
                                            mat_sigma(i,i) = sigmae**2.d0
                        if(s_cag_id.eq.1)then
                            if(ycurrent(i).gt.s_cag) then
                                nmes_o(ig) = nmes_o(ig)+1
                            end if
                        else
                            nmes_o(ig) = nmescur
                        end if
                    end do
    
    
    ! creation de Zi
    
            Z1=0.d0
                l=0
    
    
                if(nmescur.gt.0) then
                    do k=1,nb1
                        l=l+1
                        do i=1,nmescur
                            Z1(i,l)=dble(ziy(it+i,k))
                        end do
                    end do
                else
                    do i=1,nmescur
                        Z1(i,1)=0.d0
                    end do
                end if
    
                l = 0
                X2 = 0.d0
    
                do k=1,nva3
                    l = l + 1
                    do j=1,nmescur
                        X2(j,l) = dble(vey(it+j,k))
                    end do
                end do
    
            varcov_marg((it+1):(it+nmescur),1:nmescur) =Matmul( MATMUL(ziy((it+1):(it+nmescur),1:nb1), &
                    MATMUL(Ut(1:nb1,1:nb1),Utt(1:nb1,1:nb1))),transpose(ziy((it+1):(it+nmescur),1:nb1)))+ &
                    mat_sigma
    
            allocate(matv(nmescur*(nmescur+1)/2),varcov_marg_inv(nmescur,nmescur))
            matv = 0.d0
                            do j=1,nmescur
        do k=j,nmescur
        jj=j+k*(k-1)/2
        matv(jj)=varcov_marg(it+j,k)
    
        end do
        end do
        ier = 0
        eps = 1.d-10
    
    
    
            call dsinvj(matv,nmescur,eps,ier)
    
        varcov_marg_inv=0.d0
        do j=1,nmescur
                do k=1,nmescur
                            if (k.ge.j) then
                varcov_marg_inv(j,k)=matv(j+k*(k-1)/2)
                else
                varcov_marg_inv(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                end do
    
    
            element =  Matmul(Matmul(Transpose(vey(it+1:it+nmescur,1:nva3)), &
                                                            varcov_marg_inv(1:nmescur,1:nmescur)), vey(it+1:it+nmescur,1:nva3))
    
                    do j=1,nva3
                            do k=1,nva3
                                            sum_mat(j,k) = sum_mat(j,k) +  element(j,k)
                            end do
                    end do
    
    
                            deallocate(matv,varcov_marg_inv)
    
            mu = 0.d0
            mu(1:nmescur,1) = matmul(X2(1:nmescur,1:(nva3)),bh((np-nva3+1):np))
        
            xea = 0.d0
    
        choix = 3
    
            if(methodGH.le.1) then
        if(nmesy(numpat).gt.0) then
            allocate(mu1(nmesy(numpat),1))
        else
            allocate(mu1(1,1))
        end if
                if(typeJoint.eq.2.and.nb1.eq.1) then
                    call gauherJ21(int,choix,nodes_number)
                else if(typeJoint.eq.2.and.nb1.eq.2) then
                    call gauherJ22(int,choix,nodes_number)
                else if(typeJoint.eq.2.and.nb1.eq.3) then 
                    call gauherJ23(int, choix, nodes_number)
                else if(typeJoint.eq.3.and.nb1.eq.1) then
                    call gauherJ31(int,choix,nodes_number)
                else if(typeJoint.eq.3.and.nb1.eq.2) then
                    call gauherJ32(int,choix,nodes_number)
                else if(typeJoint.eq.3.and.nb1.eq.3) then
                
                    call gauherJ33(int,choix,nodes_number)
                end if
        deallocate(mu1) 
                integrale4(ig) =int 
            else
    
                call  hrmsym(nea, nf2,genz(1),genz(2),vraistot_weib, epsabs, &
                    epsrel, restar, result, abserr2, neval, ifail, work)
                integrale4(ig) =result(1)
            end if
   
    
        it_rec = it_rec + nmescurr
    
            it = it + nmescur
    
        if(integrale4(ig).gt.1.E+30) then
            integrale4(ig) = 1.E+30
        end if
    
        else
        do k=1,nb1
        Z1 (k,1)=0.d0
            end do
    
    
    
        mu = 0.d0
        !call gauherJ(int,choix)
            call  hrmsym( nea, nf2, genz(1), genz(2) ,vraistot_weib, epsabs, &
            epsrel, restar, result, abserr2, neval, ifail, work)
        it_rec = it_rec + nmescur
    
            integrale4(ig) =result(1)
        !    integrale4(ig) =int
    
    
    
            it = it + nmescur
    
        if(integrale4(ig).gt.1.E+30) then
        integrale4(ig) = 1.E+30
        end if
        end if
        !    end if
    
        deallocate(mat_sigma)
        end do
    !stop
    !************* FIN INTEGRALES **************************
        else
    
        sigmae = 1.d0
    
        end if
        res = 0.d0
    
                det = finddet(matmul(ut,utt),nea)
    
        do k=1,ng
    
        if(nb1.eq.0) then
            nmescur = 0
        else
        nmescur = nmesy(k)
    
        end if
            sum=0.d0
            if(effet.eq.0) then
            res2(k) = 0.d0
            end if
    
            !********* vraisemblance ***************
                if (integrale4(k).le.0.d0) then
                        res= res +res2(k) &
                            + res2dc(k)-0.5d0*dlog(2.d0*pi)*nmes_o(k) -0.5d0*dlog(sigmae)*nmes_o(k)  &
                            -112.d0
            else
                        res= res + res2(k) &
                            + res2dc(k)- 0.5d0*dlog(2.d0*pi)*nmes_o(k) -0.5d0*dlog(sigmae)*nmes_o(k)  &
                            +dlog(integrale4(k))
            endif
    
   ! write(*,*)k,res,integrale4(k)
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpajLongiweib=-1.d9
                    goto 123
                end if
    
    
        end do
 !   stop
    !--------------------------------------------------------
    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpajLongiweib =-1.d9
        if(typeJoint.eq.3)   Rrec = 0.d0
        if(typeJoint.eq.3)   Nrec = 0.d0
            Rdc = 0.d0
            Ndc = 0.d0
            !print*,'ok 7'
            goto 123
        else
            funcpajLongiweib = res
    
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
    
        end function funcpajlongiweib
    
    !=================================================================================================
    
    
    
        subroutine vraistot_weib(nea2,xea,nf2,funvls)
    
        use donnees_indiv
        use lois_normales
        use comon
        use optim
    !   use ParametresPourParallelisation
        use comongroup,only:vet,vet2
        use random_effect
    
        implicit none
        double precision,dimension(nea):: xea,xea2
        double precision :: funvls,vraisind,vraisindtemp,alnorm,prod_cag
        double precision,dimension(:),allocatable::ui
        integer :: nea2 ,j,i,k,nf2,nrec,ii
        double precision :: yscalar,resultdc,abserr,resabs,resasc
        double precision,external::risqindivdcCM,survRCM,survdcCM
        double precision:: res,res22
        logical :: upper
        double precision,dimension(1):: current_meanR
        double precision :: resultR
    
        upper = .false.
    
        nea2 = nea
        nf2 = nf
    
        i = numpat
    
        nrec = nig(i)
    
        Xea2=0.d0
        do j=1,nea
            Xea2(j)=Xea(j)
        end do
    
    
    
        if(nmesy(i).gt.0) then
            allocate(mu1(nmesy(i),1))
        else
            allocate(mu1(1,1))
        end if
    
        allocate(re(nea),ui(nea))
    
        ui = MATMUL(Utt,Xea2)
        mu1 = 0.d0
    
        if(nmescur.gt.0) then
            mu1(1:nmescur,1) = mu(1:nmescur,1) + MATMUL(Z1(1:nmescur,1:nb1),ui(1:nb1))
        else
            mu1(1:nmescur,1)  = mu(1:nmescur,1)
        end if
    
    
            res22 = 0.d0
    
        res = 0.d0
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le recurrences
    !ccccccccccccccccccccccccccccccccccccccccc
            
        current_meanR = 0.d0
            if(typeJoint.eq.3) then
                    if(link.eq.1) then!********* Link Random effects ************
                            do ii=1,nmescurr
    
                                    if(nva1.gt.0)then
                                            vet = 0.d0
                                            do j=1,nva1
                                                    vet =vet + b1(npp-nva+j)*dble(ve(it_rec+ii-1,j))
                                            end do
                                            vet = dexp(vet)
                                    else
                                            vet=1.d0
                                    endif
    
    
    
                                    res1(i) = res1(i)+(t1(ii)/etaR)**betaR
                                        res3(i) = res3(i)+(t0(ii)/etaR)**betaR
    
   
    
                    end do
    
       
                    res = res1(i) - res3(i)
    
        else !*********** Current Mean ******************
                    current_meanR = 0.d0
                            do ii=it_rec,it_rec+nmescurr-1
                                    resultR = 0.d0
                                    call integrationdc(survRCM,t0(ii),t1(ii),resultR,abserr,resabs,resasc,i,b1,npp,ui)
    
                                    res = res + resultR     !c'est deja res1-res3
    
                            if((c(ii).eq.1))then
                                    X2cur(1,1) = 1.d0
                                    X2cur(1,2) = t1(ii)
                                    if(nva3.gt.0) then
                                            do k=3,nva3
                                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                                            end do
                                    end if
    
                                    Z1cur(1,1) = 1.d0
                            if(nb1.eq.2)Z1cur(1,2) = t1(ii)
                                    current_meanR = current_meanR + MATMUL(X2cur,b1((npp-nva3+1):npp))&
                                            +MATMUL(Z1cur,ui(1:nb1))
                            end if
                    end do
    
            end if
            end if
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    
    
            if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
            vet2 =vet2 + b1(npp-nva3-nva2+j)*dble(vedc(i,j))
    
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
            if(link.eq.1) then
    ! pour le calcul des integrales / pour la survie, pas les donnes recurrentes:
    
    
                aux1(i)=((t1dc(numpat)/etaD)**betaD)*vet2
        
    
        else !********* Current Mean ****************
    
    call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,ui)
    
            aux1(i) = resultdc
    
            X2cur(1,1) = 1.d0
            X2cur(1,2) = t1dc(numpat)
            if(nva3.gt.0) then
                do k=3,nva3
                    X2cur(1,k) = dble(vey(it_cur+1,k))
                end do
            end if
    
        !     write(*,*)numpat,'x2cur',x2cur
            Z1cur(1,1) = 1.d0
            if(nb1.eq.2) then
                Z1cur(1,2) =t1dc(numpat)
            end if
    
            current_mean = 1.d0
    
            if(nb1.eq.1) then
                current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Z1cur(1,1)*ui(1)
            else
            current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+MATMUL(Z1cur,ui(1:nb1))
            end if
    
    
    
        end if
    
    
    
    
    
        vraisindtemp = 1.d0
        yscalar = 0.d0
        prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
            do k = 1,nmescur
                if(ycurrent(k).le.s_cag) then
                    prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
                else
                    yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                end if
            end do
        else
            do k=1,nmescur
                yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
        end if
    
        yscalar = dsqrt(yscalar)
    
    if(link.eq.1)then
            if(typeJoint.eq.3) then
    
            vraisindtemp = dlog(prod_cag) -(yscalar**2.d0)/(sigmae*2.d0)&
                                    -dexp(ui(nb1+1)+dot_product(etayr,ui(1:nb1)))&
                                    *(res1(i)-res3(i)) +nig(i)*dot_product(etayr,ui(1:nb1)) &
                                    -dexp(ui(nb1+1)*alpha+dot_product(etaydc,ui(1:nb1)))&
                                    *aux1(i) & !(res1(auxig)) &
                                    +cdc(i)*dot_product(etaydc,ui(1:nb1)) &
                                    +ui(nb1+1)*(nig(i)+alpha*cdc(i)) !&
    
           else  if(typeJoint.eq.2) then
    
            vraisindtemp = dlog(prod_cag) -(yscalar**2.d0)/(sigmae*2.d0)&
                                    -aux1(i) *dexp(dot_product(etaydc,ui(1:nb1)))  & !(res1(auxig)) &
                                    +cdc(i)*dot_product(etaydc,ui(1:nb1))
    
        
            end if
    
        else if(link.eq.2) then
            if(typeJoint.eq.3) then
    
                vraisindtemp = dlog(prod_cag)-(yscalar**2.d0)/(sigmae*2.d0)&
                    -dexp(ui(nb1+1))*res +etayr(1)*current_meanR(1) &
                    -dexp((ui(nb1+1)*alpha))*aux1(i) & !(res1(auxig)) &
                    +cdc(i)*(etaydc(1)*current_mean(1)) &
                    +ui(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
         else  if(typeJoint.eq.2) then
    
                vraisindtemp = dlog(prod_cag) -(yscalar**2.d0)/(sigmae*2.d0)&
                    -aux1(i) & !(res1(auxig)) &
                    +cdc(i)*(etaydc(1)*current_mean(1))
            end if
        end if
        vraisind = dexp(vraisindtemp)
    
        funvls = vraisind
    
    ! xeacurrent = xea
        deallocate(mu1,re)
    
        end subroutine vraistot_weib
    
        
        
        