    
    
    
    
    
    !------------------------------------------------------------------------------------
    !                   Function CVPL for Joint Model (recur/surv, longi/surv, recur/longi/surv)
    !------------------------------------------------------------------------------------
    
    
    
            subroutine cvpl_long(ng0,nsujet0,nsujety0,groupe0,groupey0,c0,cdc0,Y0,nva10,nva20,nva30,&
                            nb10,netar0,netadc0,link0,ve0,vedc0,velong0,matzy0,&
                            s_cag0,s_cag_id0,typeof0,typeJoint0,nz0,zi0,np0,b0,H_1, &
                            t00,t10,t0dc0,t1dc0,nt,valT,rl_cond,epoir,contribt,atrisk)
    
        use tailles,only:npmax
        ! use comon,only:nbintervR,nbintervDC,ttt,tttdc,ut,utt
        use comon,only:typeJoint,typeof,nst,nva,nva1,nva2,nva3,ve,vey,vedc, &
        nz1,nz2,zi,yy,c,cdc,date,ndate,datedc,ndatedc,t0,t1,t0dc,t1dc, &
            ng,nsujet,nsujety,npp,nea,nb_re,nb1,ziyd,ziyr,g,b_e,vals,effet,&
            s_cag,s_cag_id,link,netar,netadc,indic_alpha,res_ind,etaydc,etayr!,groupeey
            use lois_normales
            use donnees_indiv
    
        implicit none
    
        integer,intent(in)::np0,ng0,nsujet0,nsujety0,nt,nva10,nva20,nva30,typeof0,nz0,nb10,&
                                                                    typeJoint0,link0,netar0,netadc0
            double precision,dimension(np0),intent(in)::b0
            integer,dimension(ng0),intent(in)::cdc0
        integer,dimension(nsujet0),intent(in)::groupe0,c0
            integer,dimension(nsujety0),intent(in)::groupey0
            double precision,dimension(nsujety0),intent(in)::Y0
            double precision,dimension(nsujet0,nva10),intent(in)::ve0
        double precision,dimension(ng0,nva20),intent(in)::vedc0
        double precision,dimension(nsujety0,nva30),intent(in)::velong0
            double precision,dimension(nsujety0,nb10),intent(in)::matzy0
        double precision,dimension(-2:nz0+3),intent(in)::zi0
        double precision,dimension(np0,np0),intent(in)::H_1
        double precision,dimension(nsujet0),intent(in)::t00,t10
        double precision,dimension(ng0),intent(in)::t0dc0,t1dc0
        double precision,dimension(nt),intent(in)::valT
            double precision,intent(in)::s_cag0
            integer, intent(in)::s_cag_id0
    
        integer::i,k,t,nsujet_t,j,maxmesy
        double precision::rl_condt,trace3,min,mint,max,mindc,maxdc,maxt,mintdc,maxtdc
        double precision,dimension(0:(2*nsujety0))::aux
        double precision,dimension(ng0)::rlindiv
        double precision,dimension(np0,np0)::J_condt
            double precision,dimension(np0,np0)::mat3
        integer,dimension(ng0)::indT
    
        double precision,dimension(nt),intent(out)::rl_cond,epoir
            double precision,dimension(nt*ng0),intent(out)::contribt
            double precision,dimension(nt),intent(out)::atrisk
    
        ! redefinition des variables du module
        nst = 1

        nva1 = nva10
            nva2 = nva20
            nva3 = nva30
            if(typeJoint0.eq.3) nva = nva1+nva2+nva3
            if(typeJoint0.eq.2) nva = nva2 +nva3
            nb1 = nb10      !number of random effects in longitudinal part
            nea = nb1
    
            if(typeJoint0.eq.3) nea = nea + 1
            nb_re  =  nb1 + INT((nb1*(nb1-1))/2.d0) !number of elements to estimate from matrix B1
    
            res_ind = 0  !for integral for the current level
    ! effet = 1
        typeof = typeof0
            ng = ng0
            nz1 = nz0
        nz2 = nz0
            npp = np0
            npmax = np0
            !nz = nz0
            nsujet = nsujet0
            nsujety = nsujety0
    typeJoint = typeJoint0
    
    link = link0
    if(typeJoint.eq.2) then
    netar = 0
    netadc = netadc0
    else if(typejoint.eq.3) then
    netar = netar0
    netadc = 0
    allocate(etayr(netar))
    end if
            if(typeJoint.ne.2) then
            effet = 1
            else
            effet = 0
            end if
    
            allocate(nmes_o(ng),nmes_o2(ng))
        allocate(yy(nsujety))
            allocate(cdc(ng),c(nsujet))
            allocate(vedc(ng,nva2),ve(nsujet,nva1),vey(nsujety,nva3))
            allocate(zi(-2:nz1+3))
            allocate(ziyd(nsujety,nb1 ),ziyr(nsujety,nb1 ))
            allocate(etaydc(netadc))
            
            allocate(nii(ng),g(nsujet),b_e(npp))
    
            zi(-2:nz1+3) = zi0(-2:nz1+3)
            ziyd = matzy0
            ziyr = matzy0
        yy = y0
            c = c0
        cdc = cdc0
            ve = ve0
        vedc = vedc0
        vey = velong0
    
            if(link.eq.2)then
            allocate(x2cur(1,nva3),z1cur(1,nb1),current_mean(1))
            end if
            netadc = netadc0
            netar = netar0
    
    
            s_cag = s_cag0
            s_cag_id = s_cag_id0
    
            if(typeJoint.ne.2) then
            indic_alpha = 1
            else
            indic_alpha = 0
            end if
    
            zi(-2:nz1+3) = zi0(-2:nz1+3)
            ziyd = matzy0
            ziyr = matzy0
        yy = y0
            c = c0
        cdc = cdc0
            ve = ve0
        vedc = vedc0
        vey = velong0
    
    
            netadc = netadc0
            netar = netar0
    
    
            s_cag = s_cag0
            s_cag_id = s_cag_id0
    
            if(typeJoint.ne.2) then
            indic_alpha = 1
            else
            indic_alpha = 0
            end if
    
    
            nii = 1
            i = 1
        do j=2,nsujety
            if(groupey0(j).eq.groupey0(j-1))then
                nii(i)=nii(i)+1
            else
                i = i+1
            end if
        end do
    
        maxmesy=0
    
        do i=1,ng
            if (nii(i).gt.maxmesy) then
                maxmesy=nii(i)
            end if
        end do
    
            allocate(t0(nsujet))
            allocate(t1(nsujet))
            allocate(t0dc(ng))
            allocate(t1dc(ng))
    !       allocate(groupeey(nsujety))
            allocate(z1(maxmesy,nb1),z11(maxmesy,nb1),z2(maxmesy,nea),z22(maxmesy,nb1))
            allocate(ycurrent(maxmesy),x2(maxmesy,nva3),x22(maxmesy,nva3),ycurrent2(maxmesy))
            allocate(nii2(ng))
    
                    t0 = t00
            t1 = t10
        t0dc = t0dc0
        t1dc = t1dc0
    !  groupeey = groupey0
            g = groupe0
    
    allocate(date(0:(2*nsujet)),datedc(0:(2*nsujety)))
        date = 0.d0
    
        maxt = 0.d0
        mint = 0.d0
    
            do i=1,nsujet
            if (i.eq.1) then
                mint = t0(i) ! affectation du min juste une fois
            endif
            if (maxt.lt.t1(i)) then
                maxt = t1(i)
            endif

            if (mint.gt.t0(i)) then
                mint = t0(i)
            endif
        end do
    
        maxtdc = 0.d0
        mintdc = 0.d0
        contribt = 0.d0
    
    
        do i=1,ng
            if (i.eq.1) then
                mintdc = t0dc(i) ! affectation du min juste une fois
            endif
            if (maxtdc.lt.t1dc(i)) then
                maxtdc = t1dc(i)
            endif
            if (mintdc.gt.t0dc(i)) then
                mintdc = t0dc(i)
            endif
        end do
    
        min = 1.d-10
        max = maxt
        aux = 0.d0
    
    
                            if(typeJoint.ne.2) then
        do i = 1,(2*nsujet)
            do k = 1,nsujet
                if (t0(k).ge.min) then
                    if(t0(k).lt.max)then
                        max = t0(k)
                    endif
                endif
                if (t1(k).ge.min) then
                    if(t1(k).lt.max)then
                        max = t1(k)
                    endif
                endif
    
            end do
            aux(i) = max
            min = max + 1.d-12 ! pour virer les doublons
            max = maxt
        end do
    
        date(1) = aux(1)
        k = 1
        do i=2,(2*nsujet)
            if (aux(i).gt.aux(i-1)) then
                k = k+1
                date(k) = aux(i)
            endif
        end do
        ndate = k
    
            end if
    
        mindc = 0.d0
        maxdc = maxtdc
        aux = 0.d0
    
            do i = 1,(2*ng)
            do k = 1,ng
                if (t0dc(k).ge.mindc) then
                    if (t0dc(k).lt.maxdc) then
                        maxdc = t0dc(k)
                    endif
                endif
                if (t1dc(k).ge.mindc) then
                    if (t1dc(k).lt.maxdc) then
                        maxdc = t1dc(k)
                    endif
                endif
            end do
            aux(i) = maxdc
            mindc = maxdc + 1.d-12
            maxdc = maxtdc
        end do
    
    
            datedc(1) = aux(1)
        k = 1
        do i=2,(2*ng)
            if (aux(i).gt.aux(i-1)) then
                k = k+1
                datedc(k) = aux(i)
            endif
        end do
        ndatedc = k
            ! fin de la redefinition des variables du module
    
        b_e = b0
    
    atrisk = 0.d0
    contribt = 0.d0
    
    
        do t=1,nt ! boucle sur les temps de validation
            vals = valT(t)
            indT = 0
            nsujet_t = 0
            do i=1,ng ! boucle sur les individus
                if (t1dc(i).ge.valT(t)) then
    
                    indT(i) = 1
                    nsujet_t = nsujet_t + 1
                end if
            end do
            atrisk(t) = nsujet_t
            j=1
    
            nii2 = 0
    
            do i=1,nsujety
    
            if(vey(i,2).lt.valT(t)) then
    
            nii2(j)= nii2(j) + 1
            end if
            if(i.lt.nsujety) then
            if(groupey0(i+1).ne.groupey0(i)) then
    
            j = j + 1
            end if
            end if
            end do
    
    
                    J_condt = 0.d0
                rlindiv = 0.d0
    
    
            call derivc_condT_long(b_e,npp,J_condt,rlindiv,ng,nsujet,indT)
    
   
            rl_condt = 0.d0
            do i=1,ng
    
    
                    contribt(ng*(t-1)+i)=rlindiv(i)
    
                if (rlindiv(i).eq.-1.d9) then
                     rl_cond(t) = -1.d9
                    epoir(t) = -1.d9
                    goto 5289
                end if
                rl_condt = rl_condt + rlindiv(i)
    
           end do
    
            mat3 = MATMUL(H_1,J_condt)
              trace3 = 0.d0
    
            do k=1,npp
               trace3 = trace3 + mat3(k,k)
            end do
              epoir(t) = -rl_condt/dble(nsujet_t)+(trace3*dble(ng)/(dble(nsujet_t)*dble(ng-1))) !cvpl
            rl_cond(t) = -rl_condt/dble(nsujet_t) !mpl
    
            if (epoir(t).ne.epoir(t)) then
                epoir(t) = -1.d9
            end if
            if (rl_cond(t).ne.rl_cond(t)) then
                rl_cond(t) = -1.d9
            end if
       5289    continue
    
        end do
    
    
            deallocate(yy,cdc,vedc,vey,ve,zi)
            deallocate(nmes_o,nmes_o2)
            deallocate(ziyd,ziyr,c)
        deallocate(t0,t1,t0dc,g)
             deallocate(t1dc)
                    deallocate(nii)
            deallocate(nii2)
                    deallocate(b_e)
                    deallocate(z1)
                    deallocate(z11)
            deallocate(z2)
        deallocate(etaydc)
        if(typeJoint.eq.3)deallocate(etayr)
            deallocate(z22)
    
            deallocate(ycurrent,date,datedc)
    
            deallocate(X2)
    
            deallocate(x22)
            if(link.eq.2)deallocate(x2cur,z1cur,current_mean)
            deallocate(ycurrent2)
        end subroutine cvpl_long
    
    
    !-----------------------------------------------------------
    !                        derivc_condt
    !------------------------------------------------------------
        subroutine derivc_condt_long(b,m,V,rlindiv,nobs,nsujet,indT)
        !use comon,only:g,groupeey,nea,t1,ve,vedc,
            use comon,only  : nva3,yy, vey,ziyd,ziyr,all,b_e,&
            typeJoint,s_cag,s_cag_id,nb1
            use lois_normales
            use donnees_indiv
            use choix_epoce
        implicit none

        integer::m,i,k,id,nsujet,nobs,it,l
        double precision::funcpi_long,funcpi2_long,thn,th,z,temp1,temp2
        double precision,dimension(m,1)::Uscore,Uscore2
        double precision,dimension(m)::b
        double precision,dimension(m,m)::V
        double precision,dimension(nsujet)::rlindiv
        
        integer,dimension(nsujet)::indT
    
    
           allocate(b1(m))
            b1 = b_e
    V= 0.d0
    
         z = 0.d0
        id = 0
            it = 0
    
        do i=1,nobs
                     z1 = 0.d0
            z11 = 0.d0
            Uscore = 0.d0
            Uscore2 = 0.d0
                    nmes_o = 0
                    nmes_o2 = 0
                    nmescur = nii(i)
                    nmescur2 = nii2(i)
                    it_cur = it
    
                     ycurrent = 0.d0
                            x2 = 0.d0
   
            z22 = 0.d0
            z2 = 0.d0
  
                    do l=1,nmescur
    
    
                            if(s_cag_id.eq.1)then
                                                    if(ycurrent(l).gt.s_cag) then
                                                    nmes_o(i) = nmes_o(i)+1
                                                    end if
                                            else
                                            nmes_o(i) = nmescur
                                            end if
    
                    ycurrent(l) = yy(l+it)
                    do k=1,nva3
                   x2(l,k) =  vey(it+l,k)
                    end do
    
    
                    do k = 1,nb1
                     z2(l,k) = ziyd(it+l,k)
                      end do
    
    
                    if(typeJoint.ne.2) then
                    do k = 1,nb1
                    z1(l,k) = ziyr(it+l,k)
                   end do
                    end if
                    end do
    
    
                    all = 0
                    ycurrent2= 0.d0
                    x22 =0.d0
           
                    do l=1,nmescur2
           ycurrent2(l) = yy(l+it)
    
                    if(s_cag_id.eq.1)then
                                                    if(ycurrent2(l).gt.s_cag) then
                                                    nmes_o2(i) = nmes_o2(i)+1
                                                    end if
                                            else
                                            nmes_o2(i) = nmescur2
                                            end if
    
                    do k=1,nva3
                    x22(l,k) =  vey(it+l,k)
                    end do
    
    
                    do k = 1,nb1
                     z22(l,k) = ziyd(it+l,k)
                     end do
    
                    if(typeJoint.ne.2) then
                    do k = 1,nb1
                    z11(l,k) = ziyr(it+l,k)
                     end do
                    end if
    
                    end do
    
            if (indT(i).eq.1) then ! contribution de i sachant T
            choix_e = 1
            all = 0
    
             rlindiv(i) =funcpi_long(b,m,id,z,id,z,i)- funcpi2_long(b,m,id,z,id,z,i)
    
                 if (rlindiv(i).eq.-1.d9) then
                    V = 0.d0
                    rlindiv = -1.d9
                    goto 777
                end if
            end if
            do k=1,m
                     th =  1.d-6!DMAX1(1.d-3, (1.d-4)*DABS(b(k)))!
            thn = -1.D0*th
            if (indT(i).eq.1) then ! calcul des derivees sachant T
                                    all = 0
                                    choix_e = 1
                temp1 =  funcpi_long(b,m,k,th,id,z,i)-funcpi2_long(b,m,k,th,id,z,i)
                                    choix_e = 1
                    temp2 = funcpi_long(b,m,k,thn,id,z,i)-funcpi2_long(b,m,k,thn,id,z,i)
     
                    if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then
                        V = 0.d0
                       goto 777
                    end if
                     Uscore(k,1) = -(temp1-temp2)/(2.d0*th)
       
                            end if
                ! calcul des derivees
                            all = 1
                            choix_e = 2
                    temp1 = funcpi_long(b,m,k,th,id,z,i)!-funcpi2(b,m,k,th,id,z,i)
                temp2 = funcpi_long(b,m,k,thn,id,z,i)!-funcpi2(b,m,k,th,id,z,i)
            if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then
                    V = 0.d0
                    goto 777
                end if
                Uscore2(k,1) = -(temp1-temp2)/(2.d0*th)
    
            end do
    
                V = V + MATMUL(Uscore,transpose(Uscore2))
    
                    it = it + nmescur
    
    
    
        end do
        
    777 continue
    
            deallocate(b1)
    
        return
    
        end subroutine derivc_condt_long
    
    
    
    !-----------------------------------------------------------
    !                        FUNCPI
    !------------------------------------------------------------
    
    
        double precision function funcpi_long(b,npp,id,thi,jd,thj,i)
    
            use lois_normales
        !use comon,only:t1,g,vedc,betacoef,cdc,date,datedc,ndate,typeof,&
        !nst,vey,nz1,zidc,stra,ndatedc,nva1,nva2,t0,t1dc,tttdc,nva3,nz2,&
        !the1_e
        use comon,only:nva,&
            nea,sigmae,netar,etaydc,netadc,&
            alpha,all,effet,typeJoint,&
            etayr,betaR,etaR,betaD,etaD,ut,utt,nb1,nb_re
        use donnees_indiv
        use choix_epoce
        implicit none
    
        integer::n,npp,id,jd,i,j,k
    double precision::thi,thj
          double precision,dimension(npp)::b,bh
       double precision::vrais,int!,temp
            integer::ndim,mintps,maxtps,restar,nf2
            double precision::epsabs,epsrel, integrale4
            double precision,dimension(nea) :: xea
            external :: vraistot_e
            double precision,dimension(2):: result, abserr2
            integer ::neval,ifail,choix!,nmes
                    double precision,dimension(1000) :: work
            double precision,parameter::pi=3.141592653589793d0
    
    
    
        n = 0
        betaR = 0.d0
        etaR = 0.d0
        betaD = 0.d0
        etaD = 0.d0
    
        bh = b
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
            b1 = bh
    
            if(effet.eq.1) then
                    sigmav = bh(npp-nva-1-nb_re-effet-netadc - netar)
   
                    alpha = bh(npp-nva-nb_re-effet-netadc - netar)
    
            end if
    
    
            if(all.eq.1) then
            nmes = nmescur
            else
            nmes = nmescur2
            end if
    
            etaydc = bh(npp-nva-nb_re-netadc:npp-nva-nb_re-1)
            if(typeJoint.eq.3)etayr = bh(npp-nva-nb_re-netadc - netar:npp-nva-nb_re-netadc - 1)
    
    
            sigmae = bh(npp-nva-nb_re)*bh(npp-nva-nb_re)!
    
          
    
            allocate(Ut(nea,nea),Utt(nea,nea))
                Ut = 0.d0
            Utt = 0.d0
    
                    do j=1,nb1
                do k=1,j
                Ut(j,k)=bh(npp-nva-nb_re+k+j*(j-1)/2)
                Utt(k,j)=bh(npp-nva-nb_re+k+j*(j-1)/2)
    
                end do
            end do
    
            if(typeJoint.eq.3) then
                Ut(nea,nea) = sigmav
                Utt(nea,nea) = sigmav
            end if
    
        vrais = 0.d0
    
                    !**************INTEGRALES ****************************
    
    
            ndim = nea
    
    mintps =800
    maxtps = 1000
    epsabs = 1.d-100
    epsrel = 1.d-100
    restar = 0
    
    nf2 = nf
    integrale4 = 1.d0
            numpat =i
    
    if(nmes.gt.0) then
    
                    xea = 0.d0
    
                    choix = 1
    
            if(nb1.lt.3) then
    
            if(typeJoint.eq.2.and.nb1.eq.1) then
                    call gauherJcvpl(int,choix_e)
            else if(typeJoint.eq.2.and.nb1.eq.2) then
                    call gauherJ2cvpl(int,choix_e)
            else if(typeJoint.eq.3.and.nb1.eq.1) then
                    call gauherJ3cvpl(int,choix_e)
            else if(typeJoint.eq.3.and.nb1.eq.2) then
    
                    call gauherJ4cvpl(int,choix_e)
            end if
    
            integrale4 =int !result(1) !
    
            else
    
    
            call  hrmsym( ndim, nf2,30,1500,vraistot_e, epsabs, &
            epsrel, restar, result, abserr2, neval, ifail, work)
            integrale4 =result(1)
                    end if
    
    
            else
            integrale4= 1.d0
            end if
    
            if(integrale4.gt.1.E+30) then
    
            integrale4 = 1.E+30
            end if
    
    
    !               write(*,*)"funcpi",result(1)
            if(all.eq.1) then
                    if(integrale4.le.0) then
    
                    vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o(i)  &
                    - 0.5d0*dlog(sigmae)*nmes_o(i)  -720.d0
                    else
            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o(i) &
                    - 0.5d0*dlog(sigmae)*nmes_o(i) +dlog(integrale4)
            end if
            if ((vrais.ne.vrais).or.(abs(vrais).gt.1.d30)) then
                funcpi_long = -1.d9
        !      print*,"funcpi all1 4",vrais,integrale4
                goto 1000
            end if
            else
    
                            if(integrale4.le.0.d0) then
    
                    vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o2(i)  &
                    - 0.5d0*dlog(sigmae)*nmes_o2(i)  -720.d0
                    else
            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o2(i) &
                    - 0.5d0*dlog(sigmae)*nmes_o2(i) +dlog(integrale4)
            end if
            if ((vrais.ne.vrais).or.(abs(vrais).gt.1.d30)) then
                funcpi_long = -1.d9
        !     print*,"funcpi all0 4",vrais,integrale4,(integrale4.lt.0.d0)
                goto 1000
            end if
    
    
            end if
    
        funcpi_long= vrais
    
        if ((funcpi_long.ne.funcpi_long).or.(abs(funcpi_long).gt.1.d30)) then
            funcpi_long = -1.d9
        end if
    
    1000 continue
    
            deallocate(ut,utt)
    
        return
    
        end function funcpi_long
    
    !------------------------------------------------------------
    !               VRAISTOT_E
    !------------------------------------------------------------
    subroutine vraistot_e(ndim2,xea,nf2,funvls)

        use optim    
        use tailles
        use comongroup,only:vet,vet2
        !use comon,only:res1,res3,aux1,nig,utt
        use comon,only:alpha,cdc,sigmae,&
            nva2,nva1,npp,nva3,vedc,ve,netadc,netar,betaD,etaD,t1dc,nb1,&
            t0,t1,betaR,etaR,effet,etaydc,etayr,typeof,link,nva,vey,c,s_cag_id,s_cag,&
            all,zi,ndatedc,ndate,nea,nz1,nz2,indic_ALPHA,&
            date,datedc,vals,nsujet,g,typeJoint,ut
        use donnees_indiv
        use choix_epoce
    
    
    implicit none
    double precision,dimension(nea):: xea,xea2
    double precision :: funvls,vraisind,alnorm,prod_cag
    double precision,dimension(nea)::ui
    integer :: ndim2,jj,nf2,i,k,j,n,ndim
    double precision :: yscalar
            double precision,dimension(nmes)::mu11
        double precision,dimension(npp)::bh
        double precision::sudc
        double precision::lamdc,temp
            double precision,dimension(-2:npp)::the1,the2
            double precision,dimension(2)::su,sut1,sut0
        double precision::lam
            double precision::T, tempscl
            logical :: upper
            double precision,parameter::pi=3.141592653589793d0
            upper = .false.
            !ndim2 = nea
            bh = b1
    nf2 = nf
    ndim = ndim2
            if(choix_e.eq.2) all = 0
            if(all.eq.1) then
                    nmes = nmescur
            else
                    nmes = nmescur2
            end if
            n=0
            mu11 = 0.d0
            sut0  = 0.d0
            sut1 = 0.d0
    
            i = numpat
    Xea2=0.d0
    
            do jj=1,nea
            Xea2(jj)=Xea(jj)
        end do
    
    
    
    
                            select case(typeof)
            case(0)
    
    
                            n = (npp-nva-nea-effet-indic_ALPHA-1-netadc-netar)/(effet+1)
    
            if(effet.eq.1) then
                    do k=1,n
                            the1(k-3)=(bh(k))*(bh(k))
                    j = n+k
                            the2(k-3)=(bh(j))*(bh(j))
    
            end do
            else
            do k=1,n
                            the2(k-3)=(bh(k))*(bh(k))
            end do
            end if
    
    
            case(2)
            betaR = b1(1)**2
                            etaR = b1(2)**2
                            betaD = b1(3)**2
                            etaD = b1(4)**2
    
        end select
    
            vraisind = 1.d0
    
            ui = MATMUL(Ut,Xea2)
    
    
    
                    !ccccccccccccccccccccccccccccccccccccccccc
            ! pour les recurrences
            !ccccccccccccccccccccccccccccccccccccccccc
    if(typeJoint.ne.2) then
            if(link.eq.1) then
    
                    do k=1,nsujet ! toutes les recurrences pour l'individu i
                            if (g(k).eq.i) then
    
    
                                    if ((t1(k).gt.vals).and.(all.eq.0)) then ! si on veut filtrer les recurrences (all = 0)
                                            goto 7                               ! et si la recurrence est apres t, on sort de la boucle
                                    endif
    
                                    if(nva1.gt.0)then
                                            vet = 0.d0
                                            do j=1,nva1
                                                    vet =vet + b1(npp-nva3-2-nva2-nva1+j)*dble(ve(k,j))
                                            end do
                                            vet = dexp(vet)
                                    else
                                            vet=1.d0
                                    endif
    
                                    ! CALCUL DE LA SURVIE
                                    select case(typeof)
                                            case(0)
                                                    call susps(t1(k),the1,nz1,temp,lam,zi)
                                                    sut1(1) = temp
                                                    call susps(t0(k),the1,nz1,temp,lam,zi)
                                                    sut0(1) = temp
                    case(2)
                                                    sut1(1) = dexp(-(t1(k)/etaR)**betaR)
                                                    sut0(1) = dexp(-(t0(k)/etaR)**betaR)
                                    end select
    
                                vraisind = vraisind * (sut1(1)/sut0(1))&
                                           **(dexp(ui(nea))*vet*dexp(dot_product(etayr,ui(1:nb1))))
    
    
                                    if ((vraisind.ne.vraisind).or.(abs(vraisind).gt.1.d30)) then
                                            vraisind = -1.d9
                                !           print*,"1"
                end if
    
                                    ! CALCUL DU RISQUE
                                    if (c(k).eq.1) then
                                            select case(typeof)
                                                    case(0)
                                                            call susps(t1(k),the1,nz1,tempscl,lam,zi)
                                                            su = tempscl
                                                            if (t1(k).eq.date(ndate)) then
                                                                    lam = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
                                                            endif
                            case(2)
                                                            if (t1(k).eq.0.d0) t1(k) = 1d-12 ! utile car log(0) => -Inf
                            ! ecriture en exp(log) pour virer l'exposant
                                                            lam = (betaR*dexp((betaR-1.d0)*dlog(t1(k)))/(etaR**betaR))
                                                    end select
    
                    vraisind = vraisind * dexp(ui(nea))*lam*vet*dexp(dot_product(etayr,ui(1:nb1)))
    
    
                    if ((vraisind.ne.vraisind).or.(abs(vraisind).gt.1.d30)) then
                        vraisind = -1.d9
                !      print*,"2",vraisind,lam,vet
                                            !        goto 1000
                    end if
                            endif
                    end if
    7   continue
    
        end do
    
    
            else !*********** Current Mean ******************
        !   write(*,*)'trzeba zrobic dla current mean'
        !   stop
    
            end if
    
                    end if
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    
        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
                    vet2 =vet2 + b1(npp-nva3-nva2+j)*dble(vedc(numpat,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
                    if(choix_e.eq.1) then
                    T = t1dc(i)
                    else
                    T = vals
                    end if
    
    
                                ! CALCUL DE LA SURVIE
        select case(typeof)
            case(0)
                call susps(T,the2,nz2,temp,lamdc,zi)
                sudc = temp
        case(2)
                sudc = dexp(-(T/etaD)**betaD)
        end select
    
                    if(link.eq.1) then !******** Random Effects **********************
    
                    if(typeJoint.ne.2) then
    vraisind = vraisind*sudc**(vet2&
                            *dexp(ui(nea)*alpha+dot_product(etaydc,ui(1:nb1))))
                    else
                    vraisind = vraisind*sudc**(vet2&
                            *dexp(dot_product(etaydc,ui(1:nb1))))
                    end if
    
            !end if
    
            else !********** Current Mean ****************
    
                    X2cur(1,1) = 1.d0
                    X2cur(1,2) = T
                    if(nva3.gt.0) then
                            do k=3,nva3
                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                            end do
                    end if
    
    
                    z1cur(1,1) = 1.d0
                    z1cur(1,2) =T
    
    
                            current_mean = 1.d0
                            current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Matmul(z1cur,ui(1:2))
    
            vraisind = vraisind*sudc**(vet2&
                            *dexp(ui(nea)*alpha+etaydc(1)*current_mean(1)) )
    
            end if
    
            if(choix_e.eq.1) then
            ! CALCUL DU RISQUE
        if (cdc(i).eq.1) then
            select case(typeof)
                case(0)
                    call susps(t1dc(i),the2,nz2,sudc,lamdc,zi)
                    if (t1dc(i).eq.datedc(ndatedc)) then
                        lamdc = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
                    endif
    
                case(2)
                    if (t1dc(i).eq.0.d0) t1dc(i) = 1d-12
                    lamdc = (betaD*dexp((betaD-1.d0)*dlog(t1dc(i)))/(etaD**betaD))
            end select
    
            if(link.eq.1) then
                    if(typeJoint.ne.2) then
                            vraisind = vraisind*lamdc*vet2*dexp(ui(nea)*alpha+dot_product(etaydc,ui(1:nb1)) )
                    else
                            vraisind = vraisind*lamdc*vet2*dexp(dot_product(etaydc,ui(1:nb1)))
                    end if
            else
                            vraisind =vraisind*lamdc*vet2*dexp(ui(nea)*alpha+etaydc(1)*current_mean(1) )
            end if
    
            end if
    end if
    
    
    if(all.eq.1) then
    if(nmescur.gt.0) then
    
            mu11 =MATMUL(X2(1:nmescur,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z2(1:nmescur,1:max(netadc,netar)),ui(1:nb1))
            end if
            else
            if(nmescur2.gt.0) then
    
            mu11 =MATMUL(X22(1:nmescur2,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z22(1:nmescur2,1:max(netadc,netar)),ui(1:nb1))
            end if
            end if
    
    prod_cag = 1.d0
            !----- censure ? gauche-----------
    
    
            if(s_cag_id.eq.1)then
                    if(all.eq.1) then
    
                            do k = 1,nmescur
    
                                    if(ycurrent(k).le.s_cag) then
    
                                    prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
    
                                            mu11(k) = ycurrent(k)
    
    
                                    end if
                            end do
            else
            do k = 1,nmescur2
            if(ycurrent2(k).le.s_cag) then
                    prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
    
            mu11(k) = ycurrent2(k)
    
    
                    end if
            end do
            end if
            endif
    
            yscalar = 0.d0
            if(all.eq.1) then
            do k=1,nmescur
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
            else
                    do k=1,nmescur2
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
            end if
                    yscalar = dsqrt(yscalar)
    
    
            vraisind= vraisind *prod_cag*dexp( -(yscalar**2.d0)/(sigmae*2.d0))
    
    
            funvls = vraisind
    
    
    end subroutine vraistot_e
    
    
    
    !-----------------------------------------------------------
    !                        FUNCPI2
    !------------------------------------------------------------
    
    
        double precision function funcpi2_long(b,npp,id,thi,jd,thj,i)
    
            use lois_normales
    !use comon,only:g,t1,t0,vedc,betacoef,cdc,date,datedc,ndate,typeof,nst,&
    !vey,nz1,zidc,stra,ndatedc,nva1,nva2,t1dc,tttdc,nva3,nz2,the1_e,
    use comon,only:nva, nea,sigmae,netar,etaydc,netadc,&
            alpha,all,effet,typeJoint,&
            etayr,betaR,etaR,betaD,etaD,ut,utt,nb1,nb_re
            use donnees_indiv
            use choix_epoce
        implicit none
    
        integer::n,npp,id,jd,i,j,k
        double precision::thi,thj
      double precision,dimension(npp)::b,bh
        double precision::vrais,int!,temp
            integer::ndim,mintps,maxtps,restar,nf2
            double precision::epsabs,epsrel, integrale4
            double precision,dimension(nea):: xea
            external :: vraistot_e
            double precision,dimension(2):: result, abserr2
            integer ::neval,ifail,choix!nmes,
                    double precision,dimension(1000) :: work
            double precision,parameter::pi=3.141592653589793d0
    
        n = 0
        betaR = 0.d0
        etaR = 0.d0
        betaD = 0.d0
        etaD = 0.d0
    
        bh = b
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
            b1 = bh
            if(effet.eq.1) then
                sigmav = bh(npp-nva-1-nb_re-effet-netadc - netar)
                alpha = bh(npp-nva-nb_re-effet-netadc - netar)
    
            end if
    
    
            if(all.eq.1) then
            nmes = nmescur
            else
            nmes = nmescur2
            end if
    
               etaydc = bh(npp-nva-nb_re-netadc:npp-nva-nb_re-1)
            if(typeJoint.eq.3)etayr = bh(npp-nva-nb_re-netadc - netar:npp-nva-nb_re-netadc - 1)
    
    
            sigmae = bh(npp-nva-nb_re)*bh(npp-nva-nb_re)!
    
    
            allocate(Ut(nea,nea),Utt(nea,nea))
                Ut = 0.d0
            Utt = 0.d0
    
                    do j=1,nb1
                do k=1,j
                Ut(j,k)=bh(npp-nva-nb_re+k+j*(j-1)/2)
                Utt(k,j)=bh(npp-nva-nb_re+k+j*(j-1)/2)
    
                end do
            end do
    
            if(typeJoint.eq.3) then
                Ut(nea,nea) = sigmav
                Utt(nea,nea) = sigmav
            end if
    
    
    
        vrais = 0.d0
    
    
                    !**************INTEGRALES ****************************
    
    
    ndim = nea
    mintps = 800
    maxtps = 1000
    epsabs = 1.d-100
    epsrel = 1.d-100
    restar = 0
    nf2 = nf
    integrale4 = 1
            numpat =i
    if(nmes.gt.0) then
    
                    xea = 0.d0
         
    
                    choix = 2
                    choix_e = 2
            if(nb1.lt.3) then
    if(nea.eq.1.and.typeJoint.eq.2) then
                    call gauherJcvpl(int,choix_e)
            else if(nea.eq.2.and.typeJoint.eq.2) then
                    call gauherJ2cvpl(int,choix_e)
            else if(nea.eq.2.and.typeJoint.eq.3) then
                    call gauherJ3cvpl(int,choix_e)
            else if(nea.eq.3.and.typeJoint.eq.3) then
                    call gauherJ4cvpl(int,choix_e)
            end if
    
            integrale4 =int !result(1) !
    
            else
    
            call  hrmsym( ndim, nf2,30,1500,vraistot_e, epsabs, &
            epsrel, restar, result, abserr2, neval, ifail, work)
            integrale4 =result(1)
                    end if
    
            else
            integrale4= 1.d0
            end if
    
            if(integrale4.gt.1.E+30) then
    
            integrale4 = 1.E+30
            end if
    !               stop
    !               write(*,*)"funcpi2",result(1)
                    if(integrale4.le.0) then
    
                            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o2(i) &
                                    - 0.5d0*dlog(sigmae)*nmes_o2(i) -720.d0
                    else
                            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o2(i)&
                            - 0.5d0*dlog(sigmae)*nmes_o2(i)+dlog(integrale4)
                    end if
            if ((vrais.ne.vrais).or.(abs(vrais).gt.1.d30)) then
                funcpi2_long = -1.d9
            !   print*,"funcpi 2 4",integrale4
                goto 1000
            end if
    
    
            !end if
    
        funcpi2_long = vrais
    
        if ((funcpi2_long.ne.funcpi2_long).or.(abs(funcpi2_long).gt.1.d30)) then
    !    write(*,*)integrale4
            funcpi2_long= -1.d9
        end if
    
    
    
    1000 continue
    
    deallocate(ut,utt)
    
        return
    
        end function funcpi2_long
    
    
    
            !------------------------------------------------------------
    !               VRAISTOTLONG
    !------------------------------------------------------------
    subroutine vraistotlong(ndim2,xea,nf2,funvls)
    
    use lois_normales
    !use comon,only:cdc,datedc,vey,zidc,stra,ndatedc,t1dc,tttdc,all,etayr2,utt,the1_e
    use comon,only:typeof,nst,nbintervDC,nva,nz1,nb1, &
        date,ndate,nva1,nva2,vedc,t0,t1,&
            nb_re,nva3,nz2,sigmae,netar,&
            alpha,effet,typeJoint,netadc,&
            etayr,etaydc,nbintervR,etaR,etaD,betaR,betaD,nsujet,&
            vals,ve,g,zi,ttt,c,npp,s_cag,s_cag_id,ut,nea
    use donnees_indiv
    
    implicit none
    double precision,dimension(nea):: xea,xea2
    double precision ::funvls,vraisind,prod_cag,alnorm
    double precision,dimension(nea)::ui
    integer :: ndim2,ii,jj,nf2,k,n,j
    double precision :: yscalar
    double precision,dimension(nmes)::mu11
        integer::gg
        double precision::vet,vet2
         double precision,dimension(npp)::bh
        double precision::sudc
        double precision::temp,lamdc
            double precision,dimension(-2:npp)::the1,the2
            double precision,dimension(npp)::betacoef
            double precision,dimension(2)::su,sut1,sut0
        double precision::lam, tempscl
            logical::upper
            double precision,parameter::pi=3.141592653589793d0
    
            upper=.false.
    
            sudc = 0.d0
    bh = b1
                    ndim2 = nea
    n=0
    nf2 = nf
    ii = numpat
           nmes = nmescur2
    
    mu11 = 0.d0
     
    Xea2=0.d0
    
            do jj=1,nea
            Xea2(jj)=Xea(jj)
        end do
    
           select case(typeof)
            case(0)
    
                            if(typeJoint.ne.2) then
                            n = (npp-nva-nb_re-effet-2-netar-netadc)/2
                do k=1,n
                    the1(k-3) = (bh(k))*(bh(k))
                    j = n+k
                    the2(k-3) = (bh(j))*(bh(j))
                end do
                            else
    
                            n = (npp-nva-nb_re-1-netar-netadc)
                do k=1,n
                    the2(k-3) = (bh(k))*(bh(k))
                end do
                            end if
            case(1)
                betacoef = 0.d0
                do k=1,(nbintervR+nbintervDC)
                    betacoef(k) = bh(k)**2
                end do
            case(2)
            if(typeJoint.ne.2) then
                            betaR = bh(1)**2
                            etaR = bh(2)**2
                            betaD = bh(3)**2
                            etaD = bh(4)**2
                    else
                            betaD = bh(1)**2
                            etaD = bh(2)**2
                    end if
        end select
    
            vraisind = 1
    
            ui = MATMUL(Ut,Xea2)
    
                    ! ------------ Reccurent ------------- !
            if(typeJoint.ne.2) then
    
        do k=1,nsujet ! toutes les recurrences pour l'individu i
            if ((g(k).eq.ii).and.(t1(k).le.vals)) then
    
                if(nva1.gt.0)then
                    vet = 0.d0
                    do j=1,nva1
                        vet = vet + bh(npp-nva+j)*dble(ve(k,j))
                    end do
                    vet = dexp(vet)
                else
                    vet = 1.d0
                endif
    
                ! CALCUL DE LA SURVIE
                select case(typeof)
                    case(0)
                        call susps(t1(k),the1,nz1,temp,lam,zi)
                        sut1(1) = temp
                        call susps(t0(k),the1,nz1,temp,lam,zi)
                        sut0(1) = temp
                    case(1)
                        call survival_cpm(t1(k),bh,nst,nbintervR,ttt,sut1)
                        call survival_cpm(t0(k),bh,nst,nbintervR,ttt,sut0)
                case(2)
                        sut1(1) = dexp(-(t1(k)/etaR)**betaR)
                        sut0(1) = dexp(-(t0(k)/etaR)**betaR)
                end select
    
                vraisind  = vraisind  * (sut1(1)/sut0(1))&
                            **(dexp(ui(nea))*vet*dexp(dot_product(etayr,ui(1:nb1))))
    
                if ((vraisind .ne.vraisind ).or.(abs(vraisind ).gt.1.d30)) then
    !                 print*,"1",func2E
                    vraisind  = -1.d9
        !            goto 1000
                end if
    
                ! CALCUL DU RISQUE
                if (c(k).eq.1) then
                    select case(typeof)
                        case(0)
                            call susps(t1(k),the1,nz1,tempscl,lam,zi)
                            su = tempscl
                            if (t1(k).eq.date(ndate)) then
                                lam = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
                            endif
                        case(1)
                            do gg=1,nbintervR
                                if ((t1(k).ge.ttt(gg-1)).and.(t1(k).lt.ttt(gg))) then
                                    lam = betacoef(gg)
                                end if
                            end do
                            if (t1(k).ge.ttt(nbintervR)) then
                                lam = betacoef(nbintervR)
                            end if
                        case(2)
                            if (t1(k).eq.0.d0) t1(k) = 1d-12 ! utile car log(0) => -Inf
                            ! ecriture en exp(log) pour virer l'exposant
                            lam = (betaR*dexp((betaR-1.d0)*dlog(t1(k)))/(etaR**betaR))
                    end select
    
                vraisind  = vraisind*dexp(ui(nea))*lam*vet*dexp(dot_product(etayr,ui(1:nb1)))
    
                    if ((vraisind.ne.vraisind).or.(abs(vraisind).gt.1.d30)) then
                        vraisind = -1.d9
            !          print*,"2vraistotlong",vraisind , betaR,etaR,dlog(t1(k))
                                    end if
                endif
            endif
    
            end do
    
    
            end if
    
            ! ------------ Death ------------- !
    
        if(nva2.gt.0)then
            vet2 = 0.d0
            do jj=1,nva2
                vet2 = vet2 + b1(npp-nva3-nva2+jj)*dble(vedc(ii,jj))
            end do
    
    
            vet2 = dexp(vet2)
    
        else
            vet2 = 1.d0
        endif
    
        ! CALCUL DE LA SURVIE
    
            !    sudc = dexp(-(vals/(etaD))**(betaD))
                                ! CALCUL DE LA SURVIE
        select case(typeof)
            case(0)
                call susps(vals,the2,nz2,temp,lamdc,zi)
                sudc = temp
        case(2)
                sudc = dexp(-(vals/etaD)**betaD)
        end select
    
    
    
            if(typeJoint.ne.2) then
    vraisind = vraisind*sudc**(dexp(alpha*ui(nea))*vet2&
                            *dexp(dot_product(etaydc,ui(1:nb1))))
            else
            vraisind = vraisind*sudc**(dexp(vet2)&
                            *dexp(dot_product(etaydc,ui(1:nb1))))
    
            end if
    
        if (( vraisind .ne. vraisind ).or.(abs( vraisind ).gt.1.d30)) then
        vraisind= 1E+30
        !    print*,"3",vet2
        ! goto 745
        end if
    
    
            if(nmescur2.gt.0) then
    
            mu11 =MATMUL(X22(1:nmescur2,1:nva3),b1((npp-nva3):npp))+MATMUL(Z11(1:nmescur2,1:nb1),ui(1:nb1))&
                  +MATMUL(Z22(1:nmescur2,1:nb1),ui(1:nb1))
            end if
    
    
            prod_cag =1.d0
            !----- censure ? gauche-----------
    
            if(s_cag_id.eq.1)then
  
            do k = 1,nmescur2
            if(ycurrent2(k).le.s_cag) then
                    prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
    
            mu11(k) = ycurrent2(k)
    
    
                    end if
            end do
            end if
    
            if(nmescur2.gt.0) then
            yscalar = 0.d0
            do k=1,nmescur2
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
                    yscalar = dsqrt(yscalar)
    
    
            else
            yscalar = 0.d0
            end if
    !       end if
    
    
    
            vraisind= vraisind *prod_cag*dexp( -(yscalar**2.d0)/(sigmae*2.d0))!&
    
    funvls = vraisind
    
    end subroutine vraistotlong    
    
            !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJcvpl(ss,choix)
    
        use tailles
        !use donnees,only:w2,
        use donnees,only:x2,x3,w3,x9,w9
        !use comon,only:auxig,netadc,effet,netar,
        use comon,only:typeof,typeJoint,nea
        use donnees_indiv,only : frailpol,frailpol2
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca,func6Jcvpl,func7Jcvpl,func8Jcvpl,func9Jcvpl
        external::func6Jcvpl,func7Jcvpl,func8Jcvpl,func9Jcvpl
        integer::j
    
    
            auxfunca = 0.d0
        ss=0.d0
        if (typeof.eq.0) then
            do j=1,9
        !      if (choix.eq.3) then
    
                    if(nea.eq.1.and.typeJoint.eq.2) then
    
                                    auxfunca=func6Jcvpl(x2(j),choix)
                                    else if(nea.eq.2.and.typeJoint.eq.2) then
                                                    auxfunca = func7Jcvpl(frailpol,x9(j),choix)
                            else if(nea.eq.2.and.typeJoint.eq.3) then
                                                    auxfunca = func8Jcvpl(frailpol,x2(j),choix)
                                            else if(nea.eq.3.and.typeJoint.eq.3) then
    
                                                    auxfunca = func9Jcvpl(frailpol2,frailpol,x9(j),choix)
                                    endif
                    ss = ss+w9(j)*(auxfunca)
            !    endif
            end do
        else
            !     if (choix.eq.3) then
                            if(nea.eq.1.and.typeJoint.eq.2) then
                                    do j=1,32
                                    auxfunca=func6Jcvpl(x3(j),choix)
                                                    ss = ss+w3(j)*(auxfunca)
                                    end do
                            else if(nea.eq.2.and.typeJoint.eq.2) then
                                    do j=1,32
                                                    auxfunca = func7Jcvpl(frailpol,x3(j),choix)
                                                    ss = ss+w3(j)*(auxfunca)
                                    end do
                            else if(nea.eq.2.and.typeJoint.eq.3) then
                                    do j=1,32!15!
                                                    auxfunca = func8Jcvpl(frailpol,x3(j),choix)
                                                    ss = ss+w3(j)*(auxfunca)
                                    end do
                            else if(nea.eq.3.and.typeJoint.eq.3) then
                                    do j=1,32!15!
                                                    auxfunca = func9Jcvpl(frailpol2,frailpol,x3(j),choix)
                                                    ss = ss+w3(j)*(auxfunca)
                                    end do
                            endif
    
                    !    end if
    
        endif
    
        return
    
        END SUBROUTINE gauherJcvpl
    
            !***********************************
            !********* Gauss-Hermit pour la dimension 2*
            !*************************************
    
        SUBROUTINE gauherJ2cvpl(ss,choix)
    
        use tailles
        !use donnees,only:x2,w2,x3,w3,
        use donnees,only:x9,w9
        !use comon,only:typeof,netadc,netar,effet
        use donnees_indiv,only : frailpol!,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca,func6Jcvpl
        external::func6Jcvpl,gauherJcvpl
        integer::j
    
        ss=0.d0
    
        do j=1,9
    
            frailpol = x9(j)
    
            call gauherJcvpl(auxfunca,choix)
            ss = ss+w9(j)*(auxfunca)
    
        end do

        return
    

        END SUBROUTINE gauherJ2cvpl
    
    
                    !***********************************
            !********* Gauss-Hermit pour la dimension 2 - mod?le trviarie b_10, v*
            !*************************************
    
        SUBROUTINE gauherJ3cvpl(ss,choix)
    
        use tailles
        !use donnees,only:x2,w2,x3,w3
        use donnees,only:x9,w9
        !use comon,only:auxig,typeof,netadc
        !use donnees_indiv,only :numpat
        use donnees_indiv,only : frailpol
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca,func6Jcvpl
        external::func6Jcvpl,gauherJcvpl
        integer::j
    
        ss=0.d0
    
            do j=1,9
        !      if (choix.eq.3) then
                    frailpol = x9(j)
    
                            call gauherJcvpl(auxfunca,choix)
                    ss = ss+w9(j)*(auxfunca)
    
    
            !   endif
            end do
        
    
        return
    
        END SUBROUTINE gauherJ3cvpl
    
    
                            !***********************************
            !********* Gauss-Hermit pour la dimension 3 - mod?le trviarie b_10,b_11, v
            !*************************************
    
            SUBROUTINE gauherJ4cvpl(ss,choix)
    
        use tailles
        !use donnees,only:x2,w2,x3,w3
        use donnees,only:x9,w9
        !use comon,only:auxig,typeof,netadc
        !use donnees_indiv,only :numpat
        use donnees_indiv,only : frailpol2
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        double precision::auxfunca,func6Jcvpl
        external::func6Jcvpl,gauherJcvpl
        integer::j
    
        ss=0.d0
    
            do j=1,9
            !  if (choix.eq.3) then
                frailpol2 = x9(j)
                            call gauherJ2cvpl(auxfunca,choix)
                    ss = ss+w9(j)*(auxfunca)
            ! endif
            !         write(*,*)ss
            end do
    
    
        return
            !stop
    
        END SUBROUTINE gauherJ4cvpl
    
    
    !=====================================================================
    ! pour la loi log-normale
        double precision function func6Jcvpl(frail,choix)
    
        use tailles
        use comongroup,only:vet2!,vet
        !use comon,only:nea,date,auxig,alpha,sig2,res1,res3,aux1,nig,netar,utt,
        use comon,only:sigmae,nea,&
            nva2,npp,nva3,vedc,netadc,betaD,etaD,t1dc,etaydc,link,&
            vey,typeof,s_cag_id,s_cag,cdc,all,zi,ndatedc,nva,nz2,&
            datedc,ut,nb_re,t0dc,vals,nzdc
        use donnees_indiv
        IMPLICIT NONE
    
        double precision,dimension(nea),intent(in)::frail
            integer,intent(in)::choix
            double precision :: yscalar,alnorm,prod_cag,vraisind
            integer :: j,i,k,n
            logical :: upper
            double precision,dimension(nmescur)::mu11
            double precision,dimension(-2:npp)::the2
            double precision::sudc,T
        double precision::lamdc,temp
            double precision,parameter::pi=3.141592653589793d0
            double precision,external::survdcCM
        double precision :: abserr,resabs,resasc
    
    
            upper = .false.
            i = numpat
            if(all.eq.1) then
                    nmes = nmescur
            else
                    nmes = nmescur2
            end if
            n=0
            if(choix.eq.2) nmes = nmescur2
            mu11 = 0.d0
    
            select case(typeof)
                    case(0)
                            n = (npp-nva-nb_re-1-netadc)
    
                do k=1,n
                    the2(k-3) = (b1(k))*(b1(k))
                end do
                    case(2)
                    betaD = b1(1)**2
                            etaD = b1(2)**2
            end select
            vraisind = 1.d0
    
            if(choix.eq.1) then
                    T = t1dc(i)
                    else
                    T = vals
                    end if
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    
        if(nva2.gt.0)then
                    vet2 = 0.d0
            do j=1,nva2
                            vet2 =vet2 + b1(npp-nva3-nva2+j)*dble(vedc(numpat,j))
                    end do
                vet2 = dexp(vet2)
        else
                    vet2=1.d0
        endif
    
            ! CALCUL DE LA SURVIE
    
            if(link.eq.1) then !******** Random Effects **********************
                    select case(typeof)
                    case(0)
                            call susps(T,the2,nz2,temp,lamdc,zi)
                sudc = temp
                    case(2)
                sudc = dexp(-(T/etaD)**betaD)
        end select
    
                    vraisind = vraisind*sudc**(vet2&
                            *dexp(etaydc(1)*frail(1)))
    
            else !********** Current Mean ****************
    
                            nzdc = nz2
                call integrationdc(survdcCM,t0dc(i),T,sudc,abserr,resabs,resasc,i,b1,npp,frail)
    
    
    
                    X2cur(1,1) = 1.d0
            X2cur(1,2) =T
            if((nva3-2).gt.0) then
                do k=3,nva3
                        X2cur(1,k) = dble(vey(it_cur+1,k))
                    end do
            end if
    
            Z1cur(1,1) = 1.d0
            current_mean = 0.d0
    
            current_mean(1) =dot_product(X2cur(1,1:nva3),b1((npp-nva3+1):npp))+Z1cur(1,1)*frail(1)
    
    
                    vraisind = vraisind*dexp(-sudc)!**(vet2&
                            !*dexp(etaydc1*current_mean(1)) )
    
            end if
    
            if(choix.eq.1) then
            ! CALCUL DU RISQUE
                    if (cdc(i).eq.1) then
                            select case(typeof)
                                    case(0)
                                            call susps(T,the2,nz2,sudc,lamdc,zi)
                    if (T.eq.datedc(ndatedc)) then
                        lamdc = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
                    endif
    
                                    case(2)
                                            if (T.eq.0.d0)T = 1d-12
                                            lamdc = (betaD*dexp((betaD-1.d0)*dlog(T))/(etaD**betaD))
    
                            end select
    
                            if(link.eq.1) then
                                    vraisind = vraisind*lamdc*vet2*dexp(etaydc(1)*frail(1) )
                            else
                                    vraisind =vraisind*lamdc*vet2*dexp(etaydc(1)*current_mean(1) )
                            end if
                    end if
                    !if(frail.eq.-7.12581396102905.and.lamdc.ge.1.d0)write(*,*)i,lamdc,t1dc(i),betaD,etaD
                    !if(frail.eq.-7.12581396102905)write(*,*)i,lamdc
    
            end if
    
            if(all.eq.1) then
                    if(nmescur.gt.0) then
                            mu11 =MATMUL(X2(1:nmescur,1:nva3),b1((npp-nva3+1):npp))&
                                    +Z2(1:nmescur,1)*frail(1)
                    end if
            else
                    if(nmescur2.gt.0) then
                            mu11 =MATMUL(X22(1:nmescur2,1:nva3),b1((npp-nva3+1):npp))&
                                    +Z22(1:nmescur2,1)*frail(1)
                    end if
            end if
    
            prod_cag = 1.d0
            !----- censure ? gauche-----------
    
            if(s_cag_id.eq.1)then
                    if(all.eq.1) then
                            do k = 1,nmescur
                                    if(ycurrent2(k).le.s_cag) then
                                            prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
                                            mu11(k) = ycurrent2(k)
                                    end if
                            end do
                    else
                            do k = 1,nmescur2
                                    if(ycurrent2(k).le.s_cag) then
                                            prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
                                            mu11(k) = ycurrent2(k)
                                    end if
                            end do
                    end if
            endif
    
            yscalar = 0.d0
            if(all.eq.1) then
                    do k=1,nmescur
                            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
                    end do
            else
                    do k=1,nmescur2
                            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
                    end do
            end if
            yscalar = dsqrt(yscalar)
    
            vraisind = vraisind*prod_cag*dexp( -(yscalar**2.d0)/(sigmae*2.d0)&
                                    - (frail(1)**2.d0)/(2.d0*ut(1,1)**2))*&
                                            (1/ut(1,1))*(2.d0*pi)**(1.d0/2.d0)
    
    
        func6Jcvpl = vraisind
        return
    
        end function func6Jcvpl
    
    
            !***********************************************
            !****************** func7J  - bivariate dim 2**********************
            !***********************************************
            double precision function func7Jcvpl(frail,frail2,choix)
            
            use optim    
            use tailles
            !use comongroup,only:vet
            use comongroup,only:vet2
            !use comon,only:nea,date,auxig,alpha,sig2,res1,res3,aux1,nig,netar,&
            !res_ind
            use comon,only:cdc,sigmae,&
                nva2,npp,nva3,vedc,netadc,betaD,etaD,t1dc,etaydc,link,&
                vey,typeof,s_cag_id,s_cag,all,zi,ndatedc,nva,nz2,&
                datedc,vals,ut,utt,nb_re,t0dc,nb1,nzdc
            use donnees_indiv
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2
            integer,intent(in)::choix
            double precision :: yscalar,eps,finddet,det,alnorm,prod_cag,vraisind
            integer :: j,i,jj,k,ier,n
    double precision,dimension(nb1*(nb1+1)/2)::matv
    double precision,dimension(2,1)::  Xea2
    double precision,dimension(2):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,dimension(2,2)::mat
            double precision,dimension(nmescur)::mu11
            double precision,dimension(-2:npp)::the2
            double precision::sudc
        double precision::lamdc,temp
            double precision:: T
            logical :: upper
            double precision,parameter::pi=3.141592653589793d0
            double precision,external::survdcCM
        double precision :: abserr,resabs,resasc
    
    
    
            upper = .false.
            func7Jcvpl  = 0.d0
            i = numpat
            if(choix.eq.2) all = 0
            if(all.eq.1) then
                    nmes = nmescur
            else
                    nmes = nmescur2
            end if
            n=0
            mu11 = 0.d0
    
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea22(1) = frail
            Xea22(2) = frail2
    !       ui = MATMUL(Ut,Xea2)
    
    !       uii = matmul(Xea22,Utt)
            mat = matmul(ut,utt)
    jj=0
    ! jjj = 0
    do j=1,2
        do k=j,2
        jj=j+k*(k-1)/2
    !       jjj = jjj +1
        matv(jj)=mat(j,k)
            !  bb2vv(jjj)=bb2(j,k)
    
        end do
    end do
            ier = 0
            eps = 1.d-10
    
    
                    call dsinvj(matv,nb1,eps,ier)
    
            mat=0.d0
        do j=1,2
                    do k=1,2
                                        if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                        end do
    
    
            uii = matmul(Xea22,mat)
                    det = finddet(matmul(ut,utt),2)
    
    
                    uiiui=matmul(uii,Xea2)
    
    
    
            select case(typeof)
            case(0)
    
    
                            n = (npp-nva-nb_re-1-netadc)
                do k=1,n
                    the2(k-3) = (b1(k))*(b1(k))
    
                end do
    
    
            case(2)
    
    
                            betaD = b1(1)**2
                            etaD = b1(2)**2
    
        end select
    
    
            vraisind = 1.d0
            if(choix.eq.1) then
                    T = t1dc(i)
    
                    else
                    T = vals
                    end if
    
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    
        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
                    vet2 =vet2 + b1(npp-nva3-nva2+j)*dble(vedc(numpat,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
    
    
                                ! CALCUL DE LA SURVIE
    
    
                    if(link.eq.1) then !******** Random Effects **********************
                    select case(typeof)
            case(0)
                call susps(T,the2,nz2,temp,lamdc,zi)
                sudc = temp
        case(2)
                sudc = dexp(-(T/etaD)**betaD)
        end select
   
    vraisind = vraisind*sudc**(vet2&
                            *dexp(etaydc(1)*xea22(1)+etaydc(2)*xea22(2)))
    
            !end if
    
            else !********** Current Mean ****************
    
    !write(*,*)t1dc(i),T
    
            nzdc = nz2
    !       if(xea22(1).eq.-5.38748073577881.and.xea22(2).eq.-5.38748073577881) then
    !if((i.eq.72.or.i.eq.71).and.choix.eq.2) then
    !if(i.eq.72.or.i.eq.71) then
            call integrationdc(survdcCM,t0dc(i),T,sudc,abserr,resabs,resasc,i,b1,npp,Xea22)
    !write(*,*)i,sudc ,choix
    !end if
    !end if
            !                       write(*,*)'ble'
    !if(Xea22(1).eq.-5.38748073577881.and.Xea22(2).eq.-5.38748073577881)write(*,*)sudc
                    X2cur(1,1) = 1.d0
            X2cur(1,2) =t1dc(i)
            if((nva3-2).gt.0) then
                do k=3,nva3
                        X2cur(1,k) = dble(vey(it_cur+1,k))
                    end do
            end if
    
            Z1cur(1,1) = 1.d0
                    Z1cur(1,2) = T
            current_mean = 0.d0
    
            current_mean(1) =dot_product(X2cur(1,1:nva3),b1((npp-nva3+1):npp))&
                                    +dot_product(Z1cur(1,1:nb1),Xea22(1:nb1))
    
    ! if(T.eq.1.d0.and.Xea22(1).eq.-5.38748073577881.and.Xea22(2).eq.-5.38748073577881)write(*,*)i,&
            !                       all,choix,'cl',dexp(-sudc),sudc,res_ind
                    vraisind = vraisind*dexp(-sudc)!**(vet2&
                            !*dexp(etaydc1*current_mean(1)) )
    
            end if
    
            if(choix.eq.1) then
            ! CALCUL DU RISQUE
        if (cdc(i).eq.1) then
            select case(typeof)
                case(0)
                call susps(t1dc(i),the2,nz2,sudc,lamdc,zi)
                    if (t1dc(i).eq.datedc(ndatedc)) then
                        lamdc = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
                    endif
    
                case(2)
                    if (t1dc(i).eq.0.d0) t1dc(i) = 1d-12
                    lamdc = (betaD*dexp((betaD-1.d0)*dlog(t1dc(i)))/(etaD**betaD))
            end select
    
            if(link.eq.1) then
            
                    vraisind = vraisind*lamdc*vet2*dexp(etaydc(1)*xea22(1)+etaydc(2)*xea22(2) )
            else
                            vraisind =vraisind*lamdc*vet2*dexp(etaydc(1)*current_mean(1) )
            end if
    
            end if
    end if

    if(all.eq.1) then
    if(nmescur.gt.0) then
    
            mu11 =MATMUL(X2(1:nmescur,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
            end if
            else
            if(nmescur2.gt.0) then
    
            mu11 =MATMUL(X22(1:nmescur2,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z22(1:nmescur2,1:nb1),Xea22(1:nb1))
            end if
            end if
    
    prod_cag = 1.d0
            !----- censure ? gauche-----------
    
            if(s_cag_id.eq.1)then
                    if(all.eq.1) then
    
                            do k = 1,nmescur
    
                                    if(ycurrent(k).le.s_cag) then
                                            prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
                                            mu11(k) = ycurrent(k)
    
    
                                    end if
                            end do
            else
            do k = 1,nmescur2
            if(ycurrent2(k).le.s_cag) then
                    prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
            mu11(k) = ycurrent2(k)
    
    
                    end if
            end do
            end if
            endif
    !               if(Xea22(1).eq.-5.38748073577881.and.Xea22(2).eq.-5.38748073577881)write(*,*)mu11
            yscalar = 0.d0
            if(all.eq.1) then
            do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu11(k))**2
            end do
            else
                    do k=1,nmescur2
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
            end if
                    yscalar = dsqrt(yscalar)
    
    
    
                    vraisind = vraisind*prod_cag*dexp( -(yscalar**2.d0)/(sigmae*2.d0)&
                                    -uiiui(1)/2.d0)*det**(-0.5d0)*&
                                            (2.d0*pi)
    
    
        func7Jcvpl = vraisind
    
    
    !       deallocate(mu11)
    
        return
    
        end function func7Jcvpl
    
                    !***********************************************
            !****************** func8J **********************
            !***********************************************
            double precision function func8Jcvpl(frail,frail2,choix)
    
            use optim
    
        use tailles
        use comongroup,only:vet,vet2
        !use comon,only:res1,res3,aux1,nig,etaydc2,etayr2,indic_ALPHA,
        use comon,only:alpha,cdc,sigmae,&
            nva2,nva1,npp,nva3,vedc,ve,netadc,netar,betaD,etaD,t1dc,etaydc,etayr,&
            t0,t1,betaR,etaR,effet,typeof,link,nva,vey,c,s_cag_id,s_cag,&
            all,zi,ndatedc,ndate,nb_re,nz1,nz2,&
            date,datedc,vals,nsujet,g,ut,utt,nb1,t0dc,nzdc
        use donnees_indiv
    
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2
            integer,intent(in)::choix
            double precision :: yscalar,eps,finddet,det,alnorm,prod_cag,vraisind
            integer :: j,i,jj,k,ier,n
            double precision,dimension((netar+effet)*(netar+effet+1)/2)::matv
            double precision,dimension((netar+effet),1)::  Xea2
            double precision,dimension((netar+effet)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,dimension((netar+effet),(netar+effet))::mat
            double precision,dimension(nmescur)::mu11
            double precision,dimension(2)::su,sut1,sut0
                    double precision,dimension(-2:npp)::the2,the1
            logical :: upper
            double precision,parameter::pi=3.141592653589793d0
            double precision,external::survdcCM,survRCM
            double precision :: abserr,resabs,resasc
            double precision :: resultR
            double precision,dimension(1):: current_meanR
                    double precision::sudc
        double precision::lamdc,temp,lam, tempscl
            double precision:: T
    
    
            func8Jcvpl = 0.d0
            upper = .false.
            if(choix.eq.2) all = 0
            if(all.eq.1) then
                    nmes = nmescur
            else
                    nmes = nmescur2
            end if
            n=0
            mu11 = 0.d0
        sut1 = 0.d0 
        sut0 = 0.d0
            i = numpat
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea22(1) = frail
            Xea22(2) = frail2
    
            mat = matmul(ut,utt)
    
            jj=0
    do j=1,2
        do k=j,2
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
            end do
            end do
            ier = 0
            eps = 1.d-10
    
            call dsinvj(matv,(netar+effet),eps,ier)
    
            mat=0.d0
        do j=1,2
                    do k=1,2
                                        if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                        end do
    
                    uii = matmul(Xea22,mat)
                    det = finddet(matmul(ut,utt),2)
    
    
                    uiiui=matmul(uii,Xea2)
    
    
                    select case(typeof)
            case(0)
    
                    the1 = 0.d0
                    the2 = 0.d0
    
                            n = (npp-nva-nb_re-1-netadc-netar-2)/2
    
                    do k=1,nz1
                    the1(k-3) = (b1(k))*(b1(k))
                    j = nz1+k
                    the2(k-3) = (b1(j))*(b1(j))
    
                end do
    
            case(2)
    
    
                            betaR = b1(1)**2
                            etaR = b1(2)**2
                    betaD = b1(3)**2
                            etaD = b1(4)**2
        end select
    
    
            vraisind = 1.d0
    
            if(choix.eq.1) then
                    T = t1dc(i)
    
                    else
                    T = vals
                    end if
    
            !ccccccccccccccccccccccccccccccccccccccccc
            ! pour le risque
            !ccccccccccccccccccccccccccccccccccccccccc
            do k=1,nsujet ! toutes les recurrences pour l'individu i
            if ((g(k).eq.i).and.(t1(k).le.vals)) then
                            if(nva1.gt.0)then
                                            vet = 0.d0
                                            do j=1,nva1
                                                    vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(k,j))
                                            end do
                                            vet = dexp(vet)
                                    else
                                            vet=1.d0
                                    endif
    
                                    ! CALCUL DE LA SURVIE
    
                                    if(link.eq.1) then
                                            select case(typeof)
                                            case(0)
                                                    call susps(t1(k),the1,nz1,temp,lam,zi)
                                                    sut1(1) = temp
                                                    call susps(t0(k),the1,nz1,temp,lam,zi)
                                                    sut0(1) = temp
                    case(2)
                                                    sut1(1) = dexp(-(t1(k)/etaR)**betaR)
                                                    sut0(1) = dexp(-(t0(k)/etaR)**betaR)
                                    end select
    
                                            vraisind = vraisind * (sut1(1)/sut0(1))&
                                                            **(dexp(Xea22(2))*vet*dexp(etayr(1)*Xea22(1)))
                                    else
                                            resultR = 0.d0
                                            call integrationdc(survRCM,t0(k),t1(k),resultR,abserr,resabs,resasc,k,b1,npp,xea22(1))
                                            vraisind = vraisind * dexp(-resultR*dexp(Xea22(2)))
                                    end if
    
                                    if(sut1(1).ge.sut0(1))then
                        !         write(*,*)'ble',sut1(1),sut0(1),t1(k),t0(k),'zi',zi ,'the1',the1,'nz1',nz1
                        !          write(*,*)the1(6)
                            !       call susps2(t1(k),the1,nz1,temp,lam,zi)
                            vraisind = -1.d9
                            !       stop
                                    end if
    
            !       write(*,*)'tutaj przed', sut1(1),sut0(1),&
            !                                               dexp(Xea22(2))*vet*dexp(etayr1*Xea22(1))
                                    if ((vraisind.ne.vraisind).or.(abs(vraisind).gt.1.d30)) then
                        !           write(*,*)'tutaj', sut1(1),sut0(1),&
                        !                                    dexp(Xea22(2))*vet*dexp(etayr1*Xea22(1))
                            !               stop
                                            vraisind = -1.d9
    
                end if
    
                                    ! CALCUL DU RISQUE
                                    if (c(k).eq.1) then
                                            select case(typeof)
                                                    case(0)
                                                            call susps(t1(k),the1,nz1,tempscl,lam,zi)
                                                            su = tempscl
                                                            if (t1(k).eq.date(ndate)) then
                                                                    lam = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
                                                            endif
                            case(2)
                                                            if (t1(k).eq.0.d0) t1(k) = 1d-12 ! utile car log(0) => -Inf
                            ! ecriture en exp(log) pour virer l'exposant
                                                            lam = (betaR*dexp((betaR-1.d0)*dlog(t1(k)))/(etaR**betaR))
                                                    end select
    
                                    if(link.eq.1) then
                                            vraisind = vraisind * dexp(Xea22(2))*lam*vet*dexp(etayr(1)*Xea22(1))
                                    else
                                            X2cur(1,1) = 1.d0
                                            X2cur(1,2) = t1(k)
                                            if(nva3.gt.0) then
                                                    do j=3,nva3
                                                            X2cur(1,j) = dble(vey(it_cur+1,j))
                                                    end do
                                            end if
    
                                            Z1cur(1,1) = 1.d0
                                                    current_meanR = MATMUL(X2cur,b1((npp-nva3+1):npp))+Z1cur(1,1)*Xea22(1)
                                                    vraisind = vraisind * dexp(Xea22(2))*lam*vet*dexp(current_meanR(1))
                                    end if
    
                    if ((vraisind.ne.vraisind).or.(abs(vraisind).gt.1.d30)) then
                        vraisind = -1.d9
                      !  print*,"2",vraisind,lam,vet
                                            !        goto 1000
                    end if
                            endif
                    end if
    
    
        end do
    
            !ccccccccccccccccccccccccccccccccccccccccc
            ! pour le deces
            !ccccccccccccccccccccccccccccccccccccccccc
    
        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
                    vet2 =vet2 + b1(npp-nva3-nva2+j)*dble(vedc(numpat,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
    
                    if(link.eq.1) then !******** Random Effects **********************
                                                ! CALCUL DE LA SURVIE
        select case(typeof)
            case(0)
                call susps(T,the2,nz2,temp,lamdc,zi)
                sudc = temp
        case(2)
                sudc = dexp(-(T/etaD)**betaD)
        end select
    
    vraisind = vraisind*sudc**(vet2&
                            *dexp(alpha*Xea22(2)+etaydc(1)*xea22(1)))
    
    
            !end if
    
            else !********** Current Mean ****************
    
            nzdc = nz2
                        call integrationdc(survdcCM,t0dc(i),T,sudc,abserr,resabs,resasc,i,b1,npp,Xea22(1))
    
    
                    X2cur(1,1) = 1.d0
            X2cur(1,2) =t1dc(i)
            if((nva3-2).gt.0) then
                do k=3,nva3
                        X2cur(1,k) = dble(vey(it_cur+1,k))
                    end do
            end if
    
            Z1cur(1,1) = 1.d0
            !       Z1cur(1,2) = t1dc(i)
            current_mean = 0.d0
    
            current_mean(1) =dot_product(X2cur(1,1:nva3),b1((npp-nva3+1):npp))&
                                    +dot_product(Z1cur(1,1:nb1),Xea22(1:nb1))
    
    
                    vraisind = vraisind*dexp(-sudc*dexp(alpha*xea22(2)))!**(vet2&
                            !*dexp(etaydc1*current_mean(1)) )
            end if
    
            if(choix.eq.1) then
            ! CALCUL DU RISQUE
        if (cdc(i).eq.1) then
            select case(typeof)
                case(0)
                    call susps(t1dc(i),the2,nz2,sudc,lamdc,zi)
                    if (t1dc(i).eq.datedc(ndatedc)) then
                        lamdc = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
                    endif
    
                case(2)
                    if (t1dc(i).eq.0.d0) t1dc(i) = 1d-12
                    lamdc = (betaD*dexp((betaD-1.d0)*dlog(t1dc(i)))/(etaD**betaD))
            end select
    
            if(link.eq.1) then
            vraisind = vraisind*lamdc*vet2*dexp(alpha*Xea22(2)+etaydc(1)*xea22(1))
            else
                            vraisind =vraisind*lamdc*vet2*dexp(alpha*Xea22(2)+etaydc(1)*current_mean(1) )
            end if
    
            end if
    end if
    
            if(all.eq.1) then
                    if(nmescur.gt.0) then
    
                    mu11 =MATMUL(X2(1:nmescur,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
                    end if
            else
                    if(nmescur2.gt.0) then
    
                            mu11 =MATMUL(X22(1:nmescur2,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z22(1:nmescur2,1:nb1),Xea22(1:nb1))
                    end if
            end if
    
            prod_cag = 1.d0
            !----- censure ? gauche-----------
    
            if(s_cag_id.eq.1)then
                    if(all.eq.1) then
    
                            do k = 1,nmescur
    
                                    if(ycurrent(k).le.s_cag) then
                                            prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
                                            mu11(k) = ycurrent(k)
    
    
                                    end if
                            end do
            else
            do k = 1,nmescur2
            if(ycurrent2(k).le.s_cag) then
                    prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
            mu11(k) = ycurrent2(k)
    
    
                    end if
            end do
            end if
            endif
    
            yscalar = 0.d0
            if(all.eq.1) then
            do k=1,nmescur
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
            else
                    do k=1,nmescur2
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
            end if
                    yscalar = dsqrt(yscalar)
    
    
                            vraisind = vraisind*prod_cag*dexp( -(yscalar**2.d0)/(sigmae*2.d0)&
                                    -uiiui(1)/2.d0)*det**(-0.5d0)*&
                                            (2.d0*pi)**(-1.d0)
    
    
    
        func8Jcvpl = vraisind
    
    
        return
    
        end function func8Jcvpl
    
    
                            !***********************************************
            !****************** func9J **********************
            !***********************************************
            double precision function func9Jcvpl(frail,frail2,frail3,choix)
    use optim
    
        use tailles
        use comongroup,only:vet,vet2
        !use comon,only:res1,res3,aux1,nig,indic_ALPHA
        use comon,only:alpha,cdc,sigmae,&
            nva2,nva1,npp,nva3,vedc,ve,netadc,netar,betaD,etaD,t1dc,etaydc,etayr,&
            t0,t1,betaR,etaR,effet,typeof,link,nva,vey,c,s_cag_id,s_cag,&
            all,zi,ndatedc,ndate,nb_re,nz1,nz2,&
            date,datedc,vals,nsujet,g,ut,utt,nb1,t0dc,nzdc
        use donnees_indiv
    
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2,frail3
            integer,intent(in)::choix
            double precision :: yscalar,eps,finddet,det,alnorm,prod_cag,vraisind
            integer :: j,i,jj,k,ier,n
            double precision,dimension((netar+effet)*(netar+effet+1)/2)::matv
            !double precision,dimension(2) :: vec,vec2
            !double precision,dimension(2):: p2,p1
            double precision,dimension((netar+effet),1)::  Xea2
            double precision,dimension((netar+effet)):: uii, Xea22
            double precision,dimension(1)::uiiui
            double precision,dimension((netar+effet),(netar+effet))::mat
            double precision,dimension(nmescur)::mu11
            double precision,dimension(2)::su,sut1,sut0
            double precision,dimension(-2:npp)::the1,the2
            double precision::sudc
        double precision::lamdc,temp,lam, tempscl
            double precision:: T
            logical :: upper
            double precision,parameter::pi=3.141592653589793d0
            double precision,external::survdcCM,survRCM
            double precision :: abserr,resabs,resasc
            double precision :: resultR
            double precision,dimension(1):: current_meanR
    
    
            func9Jcvpl = 0.d0
            upper = .false.
            if(choix.eq.2) all = 0
            if(all.eq.1) then
                    nmes = nmescur
            else
                    nmes = nmescur2
            end if
    
            mu11 = 0.d0
            n = 0
            sut0 = 0.d0
            sut1 = 0.d0
            i = numpat
    
            Xea2(1,1) = frail
            Xea2(2,1) = frail2
            Xea2(3,1) = frail3
            Xea22(1) = frail
            Xea22(2) = frail2
            Xea22(3) = frail3
    
                    mat = matmul(ut,utt)
    
            jj=0
    do j=1,3
        do k=j,3
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
            end do
            end do
            ier = 0
            eps = 1.d-10
    
            call dsinvj(matv,(netar+effet),eps,ier)
    
            mat=0.d0
        do j=1,3
                    do k=1,3
                                        if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                        end do
    
                    uii = matmul(Xea22,mat)
                    det = finddet(matmul(ut,utt),3)
    
    
                    uiiui=matmul(uii,Xea2)
    
    
                    select case(typeof)
            case(0)
    
                    the1 = 0.d0
                    the2 = 0.d0
    
    
                            n = (npp-nva-nb_re-1-netadc-netar-2)/2
    
                    do k=1,n
                    the1(k-3) = (b1(k))*(b1(k))
                    j = n+k
                    the2(k-3) = (b1(j))*(b1(j))
    
                end do
    
            case(2)
    
    
                            betaR = b1(1)**2
                            etaR = b1(2)**2
                    betaD = b1(3)**2
                            etaD = b1(4)**2
        end select
    
    
            vraisind = 1.d0
    
            if(choix.eq.1) then
                    T = t1dc(i)
    
                    else
                    T = vals
                    end if
    
            !ccccccccccccccccccccccccccccccccccccccccc
            ! pour le risque
            !ccccccccccccccccccccccccccccccccccccccccc
            do k=1,nsujet ! toutes les recurrences pour l'individu i
            if ((g(k).eq.i).and.(t1(k).le.vals)) then
                            if(nva1.gt.0)then
                                            vet = 0.d0
                                            do j=1,nva1
                                                    vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(k,j))
                                            end do
                                            vet = dexp(vet)
                                    else
                                            vet=1.d0
                                    endif
    
                                    ! CALCUL DE LA SURVIE
    
                                    if(link.eq.1) then
    
                                    select case(typeof)
                                            case(0)
                                                    call susps(t1(k),the1,nz1,temp,lam,zi)
                                                    sut1(1) = temp
                                                    call susps(t0(k),the1,nz1,temp,lam,zi)
                                                    sut0(1) = temp
                    case(2)
                                                    sut1(1) = dexp(-(t1(k)/etaR)**betaR)
                                                    sut0(1) = dexp(-(t0(k)/etaR)**betaR)
                                    end select
    
                                            vraisind = vraisind * (sut1(1)/sut0(1))&
                                                            **(dexp(Xea22(3))*vet*dexp(etayr(1)*Xea22(1)+etayr(2)*Xea22(2)))
                                    else
                                            resultR = 0.d0
         call integrationdc(survRCM,t0(k),t1(k),resultR,abserr,resabs,resasc,k,b1,npp,xea22(1:nb1))
                                            vraisind = vraisind * dexp(-resultR*dexp(Xea22(3)))
                                    end if
    
    
                                    if ((vraisind.ne.vraisind).or.(abs(vraisind).gt.1.d30)) then
                                            vraisind = -1.d9
                                       !     print*,"1"
                end if
    
                                    ! CALCUL DU RISQUE
                                    if (c(k).eq.1) then
                                            select case(typeof)
                                                    case(0)
                                                            call susps(t1(k),the1,nz1,tempscl,lam,zi)
                                                            su = tempscl
                                                            if (t1(k).eq.date(ndate)) then
                                                                    lam = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
                                                            endif
                            case(2)
                                                            if (t1(k).eq.0.d0) t1(k) = 1d-12 ! utile car log(0) => -Inf
                            ! ecriture en exp(log) pour virer l'exposant
                                                            lam = (betaR*dexp((betaR-1.d0)*dlog(t1(k)))/(etaR**betaR))
                                                    end select
    
                                    if(link.eq.1) then
                                            vraisind = vraisind * dexp(Xea22(3))*lam*vet*dexp(etayr(1)*Xea22(1)+etayr(2)*Xea22(2))
                                    else
                                            X2cur(1,1) = 1.d0
                                            X2cur(1,2) = t1(k)
                                            if(nva3.gt.0) then
                                                    do j=3,nva3
                                                            X2cur(1,j) = dble(vey(it_cur+1,j))
                                                    end do
                                            end if
    
                                            Z1cur(1,1) = 1.d0
                                            Z1cur(1,2) = t1(k)
                                                    current_meanR = MATMUL(X2cur,b1((npp-nva3+1):npp))+MATMUL(Z1cur,Xea22(1:2))
                                                    vraisind = vraisind * dexp(Xea22(3))*lam*vet*dexp(current_meanR(1))
                                    end if
    
                    if ((vraisind.ne.vraisind).or.(abs(vraisind).gt.1.d30)) then
                        vraisind = -1.d9
                      !  print*,"2",vraisind,lam,vet
                                            !        goto 1000
                    end if
                            endif
                    end if
    
    
        end do
    
            !ccccccccccccccccccccccccccccccccccccccccc
            ! pour le deces
            !ccccccccccccccccccccccccccccccccccccccccc
    
        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
                    vet2 =vet2 + b1(npp-nva3-nva2+j)*dble(vedc(numpat,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
                                ! CALCUL DE LA SURVIE
    
    
                    if(link.eq.1) then !******** Random Effects **********************
                    select case(typeof)
            case(0)
                call susps(T,the2,nz2,temp,lamdc,zi)
                sudc = temp
        case(2)
                sudc = dexp(-(T/etaD)**betaD)
        end select
    
    vraisind = vraisind*sudc**(vet2&
                            *dexp(alpha*Xea22(3)+etaydc(1)*xea22(1)+etaydc(2)*xea22(2)))
    
    
            !end if
    
            else !********** Current Mean ****************
    
            nzdc = nz2
                        call integrationdc(survdcCM,t0dc(i),T,sudc,abserr,resabs,resasc,i,b1,npp,Xea22(1:2))
    
    
                    X2cur(1,1) = 1.d0
            X2cur(1,2) =t1dc(i)
            if((nva3-2).gt.0) then
                do k=3,nva3
                        X2cur(1,k) = dble(vey(it_cur+1,k))
                    end do
            end if
    
            Z1cur(1,1) = 1.d0
                    Z1cur(1,2) = t1dc(i)
            current_mean = 0.d0
    
            current_mean(1) =dot_product(X2cur(1,1:nva3),b1((npp-nva3+1):npp))&
                                    +dot_product(Z1cur(1,1:nb1),Xea22(1:nb1))
    
    
                    vraisind = vraisind*dexp(-sudc*dexp(alpha*xea22(2)))!**(vet2&
                            !*dexp(etaydc1*current_mean(1)) )
            !       write(*,*)sudc
            !              stop
            end if
    
            if(choix.eq.1) then
            ! CALCUL DU RISQUE
        if (cdc(i).eq.1) then
            select case(typeof)
                case(0)
                    call susps(t1dc(i),the2,nz2,sudc,lamdc,zi)
                    if (t1dc(i).eq.datedc(ndatedc)) then
                        lamdc = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
                    endif
    
                case(2)
                    if (t1dc(i).eq.0.d0) t1dc(i) = 1d-12
                    lamdc = (betaD*dexp((betaD-1.d0)*dlog(t1dc(i)))/(etaD**betaD))
            end select
    
            if(link.eq.1) then
            vraisind = vraisind*lamdc*vet2*dexp(alpha*Xea22(3)+etaydc(1)*xea22(1)+etaydc(2)*xea22(2))
            else
                            vraisind =vraisind*lamdc*vet2*dexp(alpha*Xea22(3)+etaydc(1)*current_mean(1) )
            end if
    
            end if
    end if
    
            if(all.eq.1) then
                    if(nmescur.gt.0) then
    
                    mu11 =MATMUL(X2(1:nmescur,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z2(1:nmescur,1:nb1),Xea22(1:nb1))
                    end if
            else
                    if(nmescur2.gt.0) then
    
                            mu11 =MATMUL(X22(1:nmescur2,1:nva3),b1((npp-nva3+1):npp))&
                                    +MATMUL(Z22(1:nmescur2,1:nb1),Xea22(1:nb1))
                    end if
            end if
    
            prod_cag = 1.d0
            !----- censure ? gauche-----------
    
            if(s_cag_id.eq.1)then
                    if(all.eq.1) then
    
                            do k = 1,nmescur
    
                                    if(ycurrent(k).le.s_cag) then
                                            prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
                                            mu11(k) = ycurrent(k)
    
    
                                    end if
                            end do
            else
            do k = 1,nmescur2
            if(ycurrent2(k).le.s_cag) then
                    prod_cag = prod_cag*(1.d0-alnorm((mu11(k)-s_cag)/sqrt(sigmae),upper))
            mu11(k) = ycurrent2(k)
    
    
                    end if
            end do
            end if
            endif
    
            yscalar = 0.d0
            if(all.eq.1) then
            do k=1,nmescur
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
            else
                    do k=1,nmescur2
            yscalar = yscalar + (ycurrent2(k)-mu11(k))**2
            end do
            end if
                    yscalar = dsqrt(yscalar)
    
            vraisind = vraisind*prod_cag*dexp( -(yscalar**2.d0)/(sigmae*2.d0)&
                                    -uiiui(1)/2.d0)*det**(-0.5d0)*&
                                            (2.d0*pi)**(-3.d0/2.d0)
    
    
                            func9Jcvpl = vraisind
            return
    
        end function func9Jcvpl
    
    