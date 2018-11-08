    
    
    
    !========================          FUNCPAJ_SPLINES         ====================
        double precision function funcpajlongisplines(b,np,id,thi,jd,thj,k0)
    
        use donnees_indiv
        use lois_normales
        use tailles
        use comon
        !use ParametresPourParallelisation
            use residusM
            use optim
        IMPLICIT NONE
    
    ! *** NOUVELLLE DECLARATION F90 :
    
        integer,intent(in)::id,jd,np
        double precision,dimension(np),intent(in)::b
        double precision,dimension(2),intent(in)::k0
        double precision,intent(in)::thi,thj
    
        ! for the numerical integral hrmsym
        integer :: restar,nf2,jj,ier
        double precision:: epsabs,epsrel
        double precision,dimension(2):: result, abserr2
        double precision,dimension(1000) :: work
        external :: vraistot,vraistot_splines,vraistot_weib
        double precision,dimension(nea) :: xea
        integer ::neval,ifail
        !
        integer::n,i,j,k,vj,ig,choix,l!,it
        integer,dimension(ngmax)::cpt
        double precision::pe1,pe2,sum,som1,som2,res,vet,vet2,h1
        double precision :: eps_s
    
        double precision,dimension(-2:npmax):: the1,the2
        double precision,dimension(np)::bh
        double precision,dimension(ngmax)::res2,res1dc,res2dc &
        ,res3dc,integrale1,integrale2,integrale3,integrale4
    !AD: for death,change dimension
        double precision,dimension(ndatemax)::dut1
        double precision,dimension(ndatemaxdc)::dut2
    !AD:end
        double precision,dimension(0:ndatemax)::ut1
        double precision,dimension(0:ndatemaxdc)::ut2
        double precision::int,eps
        double precision,parameter::pi=3.141592653589793d0
            double precision,dimension(:,:),allocatable :: mat_sigma,varcov_marg_inv
            double precision,dimension(:),allocatable :: matv
            double precision,dimension(nva3,nva3) :: element
    
    
            npp = np
        eps_s = 1.d-7
    !    print*,'debut funcpa'
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
        do i=1,np
        bh(i)=b(i)
        end do
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
            b1 = bh
    
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar)/(effet+1) !nst        !to znaczy ze dzielimy lliczbe wezlow na 2
    
        the1 = 0.d0
            the2 = 0.d0
        if(typeJoint.ne.2) then
        do i=1,n
                the1(i-3)=(bh(i))*(bh(i))
            j = n+i
            if (nst.eq.2) then
                the2(i-3)=(bh(j))*(bh(j))
            endif
        end do
        else
        do i=1,n
                the2(i-3)=(bh(i))*(bh(i))
        end do
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
                Ut(nea,nea) = sqrt(sigmav**2.d0)
                Utt(nea,nea) = sqrt(sigmav**2.d0)
            end if
    
        end if

    !----------  calcul de ut1(ti) et ut2(ti) ---------------------------
    !    attention the(1)  sont en nz=1
    !        donc en ti on a the(i)
        if(effet.eq.1) then
            dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
            ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))
        end if
        dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
    
        ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
    
        ut1(0) = 0.d0
        ut2(0) = 0.d0
    
    
    !//// NEW AMADOU vvv :
    !--- strate1
        som1 = 0.d0
        vj = 0
    
        
        if(effet.eq.1) then
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
        end if
    
    
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
    
    
                ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
                +(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))
                dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
                +(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))
    
    
        end do
    
    
        i = n-2
        h1 = (zi(i)-zi(i-1))
        if(effet.eq.1) then
        ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
            dut1(ndate) = (4.d0*the1(i-1)/h1)
        end if
    
        ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
    
        dut2(ndatedc) = (4.d0*the2(i-1)/h1)
    
    
    !-------------fin strate2
    
    
    !    print*,ndatemaxdc
    !    print*,'ut2',ut2
    
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
            integrale3(k) = 0.d0
            !integrale4(k) = 1.d0
            aux1(k)=0.d0
            aux2(k)=0.d0
        end do
    
    
    !ccccccccccccccccccccccccccccccccccccccccc
    !     pour les donnees recurrentes
    !ccccccccccccccccccccccccccccccccccccccccc
    
        
        if(effet.eq.1) then
    
        do i=1,nsujet
            cpt(g(i))=cpt(g(i))+1
            if(nva1.gt.0)then
                vet = 0.d0
                do j=1,nva1
                    vet =vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif
    
    
            if((c(i).eq.1))then
                res2(g(i)) = res2(g(i))+dlog(dut1(nt1(i))*vet)
            endif
    
            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
    
                funcpajLongisplines=-1.d9
                goto 123
            end if
    !     nouvelle version
    
            Rrec(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
    
    
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
    
                res2dc(k) = dlog(dut2(nt1dc(k))*vet2)
    
                if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
    
                    funcpajLongisplines=-1.d9
                ! print*,'gt 1'
                    goto 123
                end if
            endif
    
                RisqCumul(k) = ut2(nt1dc(k))*vet2
                Rdc_res(k) = ut2(nt1dc(k))!*vet2
    
        end do
    
    
    
    !**************INTEGRALES ****************************
    
        if(nea.ge.1) then
            it = 0
            it_rec = 1
            epsabs = 1.d-100
            epsrel = 1.d-100
            restar = 0
            nf2 = nf
            nmes_o=0
            integrale4 = 0.d0
                            sum_mat= 0.d0
    
            do ig=1,ng
    
                ycurrent  = 0.d0
                auxig=ig
                choix = 4
                numpat = ig
                nmescur =nmesy(ig)
                nmescurr = nmesrec(ig)
                nmescurr1 = nmesrec1(ig)
                it_cur = it
    
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
    
                res1cur = 0.d0
                res2cur = 0.d0
                res3cur = 0.d0
    
                if(typeJoint.eq.3) then
                do i= 0,nmescurr-1
    
                res1cur(i+1) = ut1(nt1(it_rec+i))
                    res3cur(i+1)  = ut1(nt0(it_rec+i))
                res2cur(i+1) = dut1(nt1(it_rec+i))
                end do
            end if
    
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
    
            ut2cur = ut2(nt1dc(ig))
    
                    choix = 3
            if(methodGH.le.1) then
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
                integrale4(ig) =int !result(1) !
            else
                call  hrmsym(nea, nf2,genz(1),genz(2),vraistot_splines, epsabs, &
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
            call  hrmsym( nea, nf2, genz(1), genz(2) ,vraistot_splines, epsabs, &
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
    
    !************* FIN INTEGRALES **************************
        else
            sigmae = 1.d0
        end if
    
        res = 0.d0
        do k=1,ng
            sum=0.d0
            if(nb1.ge.1) then
                nmescur = nmesy(k)
            else
                nmescur = 0
            end if
    !*************************************************************************
    !     vraisemblnce
    !*************************************************************************

            if(effet.eq.0.or.link.eq.2) then
                res2(k) = 0.d0
    
            end if
    
            if (integrale4(k).le.0.d0) then
    
                        res= res +res2(k) &
                            + res2dc(k)-0.5d0*dlog(2.d0*pi)*nmes_o(k) -0.5d0*dlog(sigmae)*nmes_o(k)  &
                            -112.d0
            else
                        res= res + res2(k) &
                            + res2dc(k)- 0.5d0*dlog(2.d0*pi)*nmes_o(k) -0.5d0*dlog(sigmae)*nmes_o(k)  &
                            +dlog(integrale4(k))
                                endif
    
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajLongisplines=-1.d9
                goto 123
            end if
  !  write(*,*)k,res,integrale4(k)
        end do
   ! stop
    !---------- calcul de la penalisation -------------------
    
        pe1 = 0.d0
        pe2 = 0.d0
    
        do i=1,n-3
        if(effet.eq.1) then
            pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
            *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
            the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
            m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
            the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
            m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
            *the1(i)*m1m(i))
            else
            pe1 = 0.d0
            end if
            if(nst.eq.1)then
                pe2=0.d0
            else
    
                pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
                *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
                the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
                m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
                the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
                m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
                *the2(i)*m1m(i))
    
            endif
        end do
    
        pe = k0(1)*pe1 + k0(2)*pe2
        resnonpen = res
    
        res = res - pe
    
    
    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpajLongisplines=-1.d9
        if(typeJoint.eq.3)    Rrec = 0.d0
        if(typeJoint.eq.3)    Nrec = 0.d0
            Rdc = 0.d0
            Ndc = 0.d0
            goto 123
    
        else
    
            funcpajLongisplines = res
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
    
    
        end function funcpajlongisplines
    
    !****************** for GENZ algorithm ********************
        subroutine vraistot_splines(nea2,xea,nf2,funvls)
    
        use donnees_indiv
        use lois_normales
        use comon
    !   use ParametresPourParallelisation
        use optim
        use random_effect
        use comongroup,only:vet,vet2
    
        implicit none
        double precision,dimension(nea):: xea,xea2
        double precision :: funvls,vraisind,vraisindtemp,alnorm,prod_cag
        double precision,dimension(:),allocatable::ui
        integer :: nea2 ,j,i,k,nf2,nrec,ii
        double precision :: yscalar
        double precision::res,res22
        logical :: upper
        double precision,external::survdcCM, survRCM
        double precision :: resultdc,abserr,resabs,resasc
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
    
    
        ui = MATMUL(Ut,Xea2)
    
    
        if(nmescur.gt.0) then
            mu1(1:nmescur,1) = mu(1:nmescur,1) + MATMUL(Z1(1:nmescur,1:nb1),ui(1:nb1))
        else
            mu1(1:nmescur,1)  = mu(1:nmescur,1)
        end if
    
        res22 = 0.d0
        res1(i) = 0.d0
        res3(i) = 0.d0
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
    
    
                                            res1(i) = res1(i)+res1cur(ii)*vet!*dexp(etayr1*ui(1)+etayr2*ui(2))
                                            res3(i) = res3(i)+res3cur(ii)*vet!*dexp(etayr1*ui(1)+etayr2*ui(2))
    

    
                    end do
    
                   res = res1(i) - res3(i)
    
        else !*********** Current Mean ******************
                    
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
                 aux1(i)=ut2cur*vet2!*dexp(etaydc1*ui(1))
      
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
    
            if ((aux1(i).ne.aux1(i)) ) then!.or.(abs(aux1(i)).ge. 1.d30)) then
    
              end if
    
    
        vraisind = 1.d0
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
    
        vraisindtemp = 0.d0
        if(link.eq.1)then
            if(typeJoint.eq.3) then
    
            vraisindtemp = dlog(prod_cag) -(yscalar**2.d0)/(sigmae*2.d0)&
                                    -dexp(ui(nb1+1)+dot_product(etayr,ui(1:nb1)))&
                                    *(res1(i)-res3(i)) +nig(i)*(dot_product(etayr,ui(1:nb1))) &
                                    -dexp(ui(nb1+1)*alpha+dot_product(etaydc,ui(1:nb1)))&
                                    *aux1(i) & !(res1(auxig)) &
                                    +cdc(i)*(dot_product(etaydc,ui(1:nb1))) &
                                    +ui(nb1+1)*(nig(i)+alpha*cdc(i)) !&
    
          
            else  if(typeJoint.eq.2) then
    
            vraisindtemp = dlog(prod_cag) -(yscalar**2.d0)/(sigmae*2.d0)&
                               -aux1(i) *dexp(dot_product(etaydc,ui(1:nb1)))  & !(res1(auxig)) &
                                    +cdc(i)*(dot_product(etaydc,ui(1:nb1)))
    
            end if
    
        else if(link.eq.2) then
            if(typeJoint.eq.3) then
    
                vraisindtemp = dlog(prod_cag)  -(yscalar**2.d0)/(sigmae*2.d0)&
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
    
        deallocate(mu1,re)
    
        end subroutine vraistot_splines
    
    
    
    