    
    
    
    
    
    !------------------------------------------------------------------------------------
    !                   Function CVPL for Joint Model NL (recur/longi/surv)
    !------------------------------------------------------------------------------------
    
    
    
            subroutine cvplNL(ng0,nsujet0,nsujety0,groupe0,groupey0,c0,cdc0,Y0,nva10,nva20,nva30,nva40,&
                            nb10,which_random0,box_cox0,netar0,netadc0,link0,ve0,vedc0,velong0,matzy0,&
                            s_cag0,s_cag_id0,typeof0,nz0,zi0,np0,b0,H_1, &
                            t00,t10,t0dc0,t1dc0,nt,valT,rl_cond,epoir,contribt,atrisk,&
                            GH,paGH0,weights0,nodes0,nnodes_all0)
    
            use tailles,only:npmax
        use comon,only:typeof,nst,nva,nva1,nva2,& !nbintervR,nbintervDC
        nva3,nva4,ve,vey,vedc, &
        nz1,nz2,zi,yy,c,cdc,date,ndate,datedc,ndatedc,t0,t1,t0dc,t1dc,g, & !ttt,tttdc
            ng,nsujet,nsujety,npp,nea,nb_re,nb1,ziy,g,b_e,vals,netadc,netar,&
            s_cag,s_cag_id,link,netar,netadc,indic_alpha,res_ind,& !ut,utt,
            nodes,weights,invBi_chol,b_lme,invBi_cholDet,methodGH,nnodes_all,nodes_number,&
            etaydc,etayr,mat,matb_chol,box_cox1,box_cox_par,which_random
            use lois_normales
            use donnees_indiv
    
        implicit none
    
        integer,intent(in)::np0,ng0,nsujet0,nsujety0,nt,nva10,nva20,nva30,nva40,typeof0,nz0,nb10,&
                           link0,netar0,netadc0,nnodes_all0,which_random0
            double precision,dimension(np0),intent(in)::b0
            integer,dimension(ng0),intent(in)::cdc0
        integer,dimension(nsujet0),intent(in)::groupe0,c0
        double precision,dimension(2),intent(in) :: box_cox0
            integer,dimension(nsujety0),intent(in)::groupey0
            double precision,dimension(nsujety0),intent(in)::Y0
            double precision,dimension(nsujet0,nva10),intent(in)::ve0
        double precision,dimension(ng0,nva20),intent(in)::vedc0
        double precision,dimension(nsujety0,nva30+nva40),intent(in)::velong0
            double precision,dimension(nsujety0,2),intent(in)::matzy0
        double precision,dimension(-2:nz0+3),intent(in)::zi0
        double precision,dimension(np0,np0),intent(in)::H_1
        double precision,dimension(nsujet0),intent(in)::t00,t10
        double precision,dimension(ng0),intent(in)::t0dc0,t1dc0
        double precision,dimension(nt),intent(in)::valT
            double precision,intent(in)::s_cag0
            integer, intent(in)::s_cag_id0
        integer,dimension(2),intent(in):: GH
        double precision,dimension(ng0,nb10+1+1+nb10 +1+ ((nb10+1)*(nb10+1-1))/2),intent(in):: paGH0
        double precision,dimension(nnodes_all0,(nb10+1)),intent(in):: nodes0,weights0 
    
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
    
    
    allocate(nodes(nnodes_all0,nb10+1),weights(nnodes_all0,nb10+1))
    
    nodes = nodes0
    weights= weights0
    nnodes_all = nnodes_all0
    

        ! redefinition des variables du module
        nst = 1

        nva1 = nva10
            nva2 = nva20
            nva3 = nva30
            nva4 = nva40
          nva = nva1+nva2+nva3+nva4
            nb1 = nb10      !number of random effects in longitudinal part
            nea = nb1
    
             nea = nea + 1
            nb_re  =  nb1 !+ (nb1*(nb1-1))/2.d0 !number of elements to estimate from matrix B1
    
            res_ind = 0  !for integral for the current level

        typeof = typeof0
            ng = ng0
            nz1 = nz0
        nz2 = nz0
            npp = np0
            npmax = np0
            !nz = nz0
            nsujet = nsujet0
            nsujety = nsujety0
    
    link = link0
   
    netar = netar0
    netadc = 0

    
    
            allocate(nmes_o(ng),nmes_o2(ng))
        allocate(yy(nsujety))
            allocate(cdc(ng),c(nsujet))
            allocate(vedc(ng,nva2),ve(nsujet,nva1),vey(nsujety,nva3+nva4))
            allocate(zi(-2:nz1+3))
            allocate(ziy(nsujety,2))
    
            allocate(nii(ng),g(nsujet),b_e(npp))
    
            zi(-2:nz1+3) = zi0(-2:nz1+3)
            ziy= matzy0
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
        allocate(etaydc(netadc),etayr(netar),mat(nea,nea),matb_chol(nea,nea))
    
            s_cag = s_cag0
            s_cag_id = s_cag_id0

            
            box_cox1 = int(box_cox0(1))
        if(box_cox1.eq.1)then 
        box_cox_par = box_cox0(2)
        else 
        box_cox_par = int(box_cox0(2))
        end if
        
        which_random = which_random0
              indic_alpha = 1

       
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
          allocate(z1(maxmesy,2),z11(maxmesy,2))!,z2(maxmesy,nea),z22(maxmesy,nb1))
            allocate(ycurrent(maxmesy),x2(maxmesy,nva3+nva4),x22(maxmesy,nva3+nva4),ycurrent2(maxmesy))
            allocate(nii2(ng), mu(maxmesy,2))
    
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
    !         if ((maxt.lt.tU(i)).and.(tU(i).ne.t1(i))) then
    !             maxt = tU(i)
    !         endif
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
  
    
    allocate(b_lme(ng,nea),invBi_cholDet(ng),invBi_chol(ng,nea + (nea*(nea-1))/2))
     do i=1,ng
        b_lme(i,1:nea) = paGH0(i,1:nea)
        invBi_cholDet(i) = paGH0(i,nea+1)
        invBi_chol(i,1:nea + (nea*(nea-1))/2) = paGH0(i,(nea+2):(nea+1+nea + (nea*(nea-1))/2))
    end do
    
            methodGH = GH(1)
            nodes_number = GH(2)
    
    
    
    
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
                    if(ziy(i,1).lt.valT(t)) then
    
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
    
            call derivc_condT_NL(b_e,npp,J_condt,rlindiv,ng,nsujet,indT)
    
    
            rl_condt = 0.d0
            do i=1,ng

    
                    contribt(ng*(t-1)+i)=rlindiv(i)
    
                if (rlindiv(i).eq.-1.d9) then
                !    print*,"oups",i,rlindiv(i)
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
            deallocate(ziy,c)
        deallocate(t0,t1,t0dc,g)
            deallocate(t1dc)
                    deallocate(nii)
            deallocate(nii2)
                    deallocate(b_e)
                    deallocate(z1)
                    deallocate(z11,mu)
           
    
            deallocate(ycurrent,date,datedc)
    
            deallocate(X2)
    
            deallocate(x22)
            if(link.eq.2)deallocate(x2cur,z1cur,current_mean)
            deallocate(ycurrent2)
            deallocate(nodes,weights)
            deallocate(b_lme,invBi_cholDet,invBi_chol)
            deallocate(etaydc,etayr,mat,matb_chol)
            
        end subroutine cvplNL
    
    
    !-----------------------------------------------------------
    !                        derivc_condt
    !------------------------------------------------------------
        subroutine derivc_condt_NL(b,m,V,rlindiv,nobs,nsujet,indT)
        !use comon,only:g
            use comon,only  : nva3,nva4,yy, &
                    vey,ziy,all,b_e,&
                    s_cag,s_cag_id !nb1,g,vedc,ve,groupeey,t1 !nea
                    use lois_normales
                    use donnees_indiv
                    use choix_epoce
        implicit none
    
    integer::m,i,k,id,nsujet,nobs,it,l
        double precision::funcpi_NL,funcpi2_NL,thn,th,z,temp1,temp2
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
 
    
     do l=1,nmescur
        if(s_cag_id.eq.1)then
            if(ycurrent(l).gt.s_cag) then
                nmes_o(i) = nmes_o(i)+1
            end if
         else
             nmes_o(i) = nmescur
         end if

         ycurrent(l) = yy(l+it)
        do k=1,nva3+nva4
             x2(l,k) =  vey(it+l,k)
        end do
      
        z1(l,1:2) = ziy(it+l,1:2)
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
    
        do k=1,nva3+nva4
            x22(l,k) =  vey(it+l,k)
        end do
    
        z11(l,1:2) = ziy(it+l,1:2)
    end do
    
            if (indT(i).eq.1) then ! contribution de i sachant T
            choix_e = 1
            all = 0
    

  rlindiv(i) =funcpi_NL(b,m,id,z,id,z,i)- funcpi2_NL(b,m,id,z,id,z,i)
             if (rlindiv(i).eq.-1.d9) then
                    V = 0.d0
                    rlindiv = -1.d9
                    goto 777
                end if
            else 
            
            end if
            do k=1,m
                    th = 1.d-3! DMAX1(1.d-3, (1.d-4)*DABS(b(k)))!1.d-6
            thn = -1.D0*th
            if (indT(i).eq.1) then ! calcul des derivees sachant T
                                    all = 0
                                    choix_e = 1                
               temp1 =  funcpi_NL(b,m,k,th,id,z,i)-funcpi2_NL(b,m,k,th,id,z,i)
                                   choix_e = 1
                    temp2 = funcpi_NL(b,m,k,thn,id,z,i)-funcpi2_NL(b,m,k,thn,id,z,i)
   
                    if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then
                        V = 0.d0
                        goto 777
                    end if
                    Uscore(k,1) = -(temp1-temp2)/(2.d0*th)
                             end if
                ! calcul des derivees
                            all = 1
                            choix_e = 2
                   temp1 = funcpi_NL(b,m,k,th,id,z,i)
               temp2 = funcpi_NL(b,m,k,thn,id,z,i)
                if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then
                    V = 0.d0
                    !rlindiv = -1.d9
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
    
        end subroutine derivc_condt_NL
    
    
    
    !-----------------------------------------------------------
    !                        FUNCPI
    !------------------------------------------------------------
    
    
        double precision function funcpi_NL(b,npp,id,thi,jd,thj,i)
    use optim
            use lois_normales
        use comon,only:nva, & !typeof,nst,vey,nz1,
        !zidc,stra,date,datedc,ndate,ndatedc,nva1,nva2,vedc,t0,t1,t1dc,tttdc,g,cdc,&
            nea,sigmae,netar,etaydc,netadc,netar ,& !nva3,nz2,the1_e,betacoef,g
            alpha,all,netadc,netar,&
            etayr,alpha,betaR,etaR,betaD,etaD,ut,utt,nb1,nb_re,&
            K_G0,K_D0,lambda,y0,det,mat,indic_alpha
            use donnees_indiv
            use choix_epoce
        implicit none
    
        integer::n,npp,id,jd,i,j,k,jj,ier
    double precision::thi,thj,eps,finddet
        double precision,dimension(npp)::b,bh
        double precision::vrais,int!,temp
            integer::ndim,mintps,maxtps,restar,nf2
            double precision::epsabs,epsrel, integrale4
            double precision,dimension(nea) :: xea
            integer ::choix!,nmes
            double precision,parameter::pi=3.141592653589793d0
     double precision,dimension(nea*(nea+1)/2)::matv
    
        
        n = 0
        betaR = 0.d0
        etaR = 0.d0
        betaD = 0.d0
        etaD = 0.d0
    
        bh = b
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
   
            b1 = bh
    
       
            sigmav = bh(npp-nva-nb_re-1-netadc - netar-indic_alpha)
            alpha = bh(npp-nva-nb_re-1-netadc - netar)
    
    
        K_G0 =  bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha-3)
        K_D0 =  bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha-2)
        lambda =  bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha-1)
        y0 = bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha)
    
            if(all.eq.1) then
            nmes = nmescur
            else
            nmes = nmescur2
            end if
    
            etaydc(1:nb1) = bh(npp-nva-nb_re-netadc:npp-nva-nb_re-1)
            sigmae = bh(npp-nva-nb_re)*bh(npp-nva-nb_re)!
            etayr(1:nb1) = bh(npp-nva-nb_re-netadc-netar:npp-nva-nb_re-netadc-1)
    

            allocate(Ut(nea,nea),Utt(nea,nea))
                Ut = 0.d0
            Utt = 0.d0
    
                    do j=1,nb1
                do k=1,j
                 if(j.eq.k) then
                Ut(j,k)=sqrt(bh(npp-nva-nb_re+k)**2.d0)
                  Utt(k,j)=sqrt(bh(npp-nva-nb_re+k)**2.d0)
                  end if
                end do
            end do
    
            
                Ut(nea,nea) = sigmav
                Utt(nea,nea) = sigmav
          
            
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
    
             call gauherJ3cvplNL(int,choix_e)

    
        integrale4 =int
    
     else
        integrale4= 1.d0
     end if
            
    if(integrale4.gt.1.E+30) then
        integrale4 = 1.E+30
    end if
    
    if(all.eq.1) then
        if(integrale4.le.0) then
            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o(i)  &
                        - 0.5d0*dlog(sigmae)*nmes_o(i)  -720.d0
        else
            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o(i) &
                        - 0.5d0*dlog(sigmae)*nmes_o(i) +dlog(integrale4)
        end if
        if ((vrais.ne.vrais).or.(abs(vrais).gt.1.d30)) then
                funcpi_NL = -1.d9
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
                funcpi_NL = -1.d9
        !     print*,"funcpi all0 4",vrais,integrale4,(integrale4.lt.0.d0)
                goto 1000
        end if
    end if
    
        funcpi_NL= vrais
    
        if ((funcpi_NL.ne.funcpi_NL).or.(abs(funcpi_NL).gt.1.d30)) then
            funcpi_NL = -1.d9
        end if
    
    1000 continue
    
            deallocate(ut,utt)
     
        return
    
    end function funcpi_NL
    
 
    
    
    !-----------------------------------------------------------
    !                        FUNCPI2
    !------------------------------------------------------------
    
    
        double precision function funcpi2_NL(b,npp,id,thi,jd,thj,i)
    
    use optim
            use lois_normales
    use comon,only:nva, & !typeof,nst,vey,nz1,
        !zidc,stra,date,datedc,ndate,ndatedc,nva1,nva2,vedc,t0,t1,t1dc,tttdc,g,cdc,&
            nea,sigmae,netar,etaydc,netadc,netar ,& !nva3,nz2,the1_e,betacoef,g
            alpha,all,netadc,netar,&
            etayr,alpha,betaR,etaR,betaD,etaD,ut,utt,nb1,nb_re,&
            K_G0,K_D0,lambda,y0,det,mat,indic_alpha        
            use donnees_indiv
            use choix_epoce
        implicit none
    
        integer::n,npp,id,jd,i,j,k,ier,jj
        double precision::thi,thj,eps,finddet
        double precision,dimension(npp)::b,bh
        double precision::vrais,int!,temp
            integer::ndim,mintps,maxtps,restar,nf2
            double precision::epsabs,epsrel, integrale4
            double precision,dimension(nea):: xea
            integer ::choix!nmes,
            double precision,parameter::pi=3.141592653589793d0
    double precision,dimension(nea*(nea+1)/2)::matv
    
        n = 0
        betaR = 0.d0
        etaR = 0.d0
        betaD = 0.d0
        etaD = 0.d0
    
        bh = b
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
            b1 = bh
         
            sigmav = bh(npp-nva-nb_re-1-netadc - netar-indic_alpha)
            alpha = bh(npp-nva-nb_re-1-netadc - netar)
    
    
        K_G0 =  bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha-3)
        K_D0 =  bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha-2)
        lambda =  bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha-1)
        y0 = bh(npp-nva-nb_re-1-netadc - netar-1-indic_alpha)
    
    
            if(all.eq.1) then
            nmes = nmescur
            else
            nmes = nmescur2
            end if
    
            etaydc(1:nb1) = bh((npp-nva-nb_re-netadc):(npp-nva-nb_re-1))
    
            sigmae = bh(npp-nva-nb_re)*bh(npp-nva-nb_re)!
    
            etayr(1:nb1) = bh((npp-nva-nb_re-netadc-netar):(npp-nva-nb_re-netadc-1))
        
        
            allocate(Ut(nea,nea),Utt(nea,nea))
                Ut = 0.d0
            Utt = 0.d0
    
                       do j=1,nb1
                do k=1,j
                 if(j.eq.k) then
                Ut(j,k)=sqrt(bh(npp-nva-nb_re+k)**2.d0)
                  Utt(k,j)=sqrt(bh(npp-nva-nb_re+k)**2.d0)
                  end if
                end do
            end do
    
          
                Ut(nea,nea) = sigmav
                Utt(nea,nea) = sigmav
           
    
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
            call gauherJ3cvplNL(int,choix_e)
            
         
      
            integrale4 =int !result(1) !
    
           
    
            else
            integrale4= 1.d0
            end if
    
            if(integrale4.gt.1.E+30) then
    
            integrale4 = 1.E+30
            end if
    !               stop
  
  if(integrale4.le.0) then
    
                            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o2(i) &
                                    - 0.5d0*dlog(sigmae)*nmes_o2(i) -720.d0
                    else
                            vrais =   -0.5d0*dlog(2.d0*pi)*nmes_o2(i)&
                            - 0.5d0*dlog(sigmae)*nmes_o2(i)+dlog(integrale4)
                    end if
            if ((vrais.ne.vrais).or.(abs(vrais).gt.1.d30)) then
                funcpi2_NL = -1.d9
            !   print*,"funcpi 2 4",integrale4
                goto 1000
            end if
    
    
            !end if
    
        funcpi2_NL = vrais
    
        if ((funcpi2_NL.ne.funcpi2_NL).or.(abs(funcpi2_NL).gt.1.d30)) then
   
   funcpi2_NL= -1.d9
        end if
    
    
    
    1000 continue
    
    deallocate(ut,utt)
    
        return
    
        end function funcpi2_NL
    
    
    
    
    
    
    
        
    
    
            !***********************************
            !********* Gauss-Hermit - modele trviarie 
            !*************************************
    
            SUBROUTINE gauherJ3cvplNL(ss,choix)
    
        
           use optim
        use tailles
        use comongroup,only:vet,vet2
        
        use comon,only:cdc,sigmae,& !auxig,alpha,sig2,res1,res3,aux1,nig,nmesy
            nva2,npp,nva3,nva4,vedc,nea,nb1,betaD,etaD,t1dc,etaydc,link,& !t0dc
            typeof,s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,ut,utt,
            nodes,weights,mat,det, K_G0, K_D0, lambda, y0,invBi_cholDet,netadc,& !yy,ziy,initGh,it
            all,vals,nva,ndatedc,nz2,zi,datedc,t1,betaR,etaR,etayr,indic_alpha,&
            ndate,netar,nsujet,nva1,nz1,g,ve,t0,c,date,box_cox_par,box_cox1,which_random
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix
        !double precision::aux
       
        !double precision::frail,frail2,frail3
        double precision :: yscalar,alnorm,prod_cag,T,vraisind !eps
        integer :: j,i,k,nnn,nn !jj,ier
    !    double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
       ! double precision,dimension(nb1,nb1)::mat
        logical :: upper
        double precision,external::survdcCM,survRCM
            double precision :: func10J !abserr,resabs,resasc
            !double precision :: resultR
            !double precision,dimension(1):: current_meanR
        double precision,parameter::pi=3.141592653589793d0
            double precision,dimension(nmescur,1)::mu11
          double precision,dimension(2)::su,sut1,sut0
            double precision,dimension(-2:npp)::the1,the2
            double precision::sudc
        double precision::lamdc,temp,lam, tempscl
        
        upper = .false.
        ss=0.d0
            
        if(choix.eq.2) all = 0    
        i = numpat
        
        if(all.eq.1) then
            nmes = nmescur
         else
           nmes = nmescur2
         end if
    mu = 0.d0
    mu(1:nmes,1) = matmul(x2(1:nmes,1:(nva3)),b1((npp-nva4-nva3+1):(npp-nva4)))
    mu(1:nmes,2) = matmul(x2(1:nmes,(nva3+1):(nva3+nva4)),b1((npp-nva4+1):npp))

    do nnn=1,nnodes_all
           func10J = 0.d0
            mu11 = 0.d0
           
             matb_chol = 0.d0
          
         if(methodGH.eq.1) then
         do k=1,nea
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(nnn,1:nea)
            Xea22(1:nea) = b_lme(i,1:nea) +  Matmul(matb_chol(1:nea,1:nea),Xea(1:nea))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(nnn,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
    !        Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)! v_jf(i)+sqrt(2.d0*varv_jf(i))*Xea(nb1+1)!
              Xea2(1:nea,1) = Xea22(1:nea)
        
            else if(methodGH.eq.3) then 
             do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(nnn,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
            Xea22(nb1+1) = Xea(nb1+1)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            
            else if(methodGH.eq.0) then
            Xea2(1:nea,1) = nodes(nnn,:)
        Xea22(1:nea) = nodes(nnn,:)
            end if

       uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)

    select case(typeof)
        case(0)
            nn = (npp-nva-1-indic_ALPHA-1-nb1 - netadc - netar-4)/2
            do k=1,nn
                the1(k-3) = (b1(k))*(b1(k))
                j = nn+k
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
                     vet =vet + b1(npp-nva3-nva4-nva2-nva1+j)*dble(ve(k,j))
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
                  **(dexp(dot_product(etayr,Xea22(1:netar))))
         else

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
                         lam = 4.d0*the1(nn-2-1)/(zi(nn-2)-zi(nn-2-1))
                     endif
                case(2)
                    if (t1(k).eq.0.d0) t1(k) = 1d-12 ! utile car log(0) => -Inf
                            ! ecriture en exp(log) pour virer l'exposant
                    lam = (betaR*dexp((betaR-1.d0)*dlog(t1(k)))/(etaR**betaR))
                end select
    
         if(link.eq.1) then
            vraisind = vraisind * dexp(Xea22(nea))*lam*vet*&
                    dexp(dot_product(etayr,Xea22(1:netar)))
                 
         else
                                          
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
             vet2 =vet2 + b1(npp-nva4-nva3-nva2+j)*dble(vedc(numpat,j))
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
    
        vraisind = vraisind*sudc**(vet2*dexp(dot_product(etaydc,Xea22(1:netadc))))
    else !********** Current Mean ****************
       
    end if
   
    if(choix.eq.1) then
            ! CALCUL DU RISQUE
        if (cdc(i).eq.1) then
        
            select case(typeof)
            case(0)
                call susps(t1dc(i),the2,nz2,sudc,lamdc,zi)
                 if (t1dc(i).eq.datedc(ndatedc)) then
                    lamdc = 4.d0*the2(nn-2-1)/(zi(nn-2)-zi(nn-2-1))
                    
                 endif
            case(2)
                if (t1dc(i).eq.0.d0) t1dc(i) = 1d-12
                lamdc = (betaD*dexp((betaD-1.d0)*dlog(t1dc(i)))/(etaD**betaD))
            end select
    
            if(link.eq.1) then
    
                    vraisind = vraisind*lamdc*vet2*dexp(dot_product(etaydc,Xea22(1:netadc)))
            
            else
            
            end if
    
        end if
    end if
    
    
    !************ODE - analytic solution *********************

    
        if(nmes.gt.0) then
    
        
        
             select case(which_random)
        case(1)
        mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0-lambda+mu(1:nmes,2))*&
        (dexp(-dexp(lambda)*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(2)
            mu11(1:nmes,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0-lambda+mu(1:nmes,2))*&
        (dexp(-dexp(lambda)*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(3)
            mu11(1:nmes,1)  =((dexp(y0+dexp(K_G0+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(1)-lambda+mu(1:nmes,2))*&
        (dexp(-dexp(lambda)*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(4)
            mu11(1:nmes,1)  =((dexp(y0+dexp(K_G0+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0-lambda-Xea22(1)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(1))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(5)
        mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0-lambda+mu(1:nmes,2))*&
        (dexp(-dexp(lambda)*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(6)
            mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmes,2))*&
        (dexp(-dexp(lambda)*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(7)
            mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(2))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(8)
            mu11(1:nmes,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(2)-lambda+mu(1:nmes,2))*&
        (dexp(-dexp(lambda)*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(9)
            mu11(1:nmes,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0-lambda-Xea22(2)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(2))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(10)
            mu11(1:nmes,1)  =((dexp(y0+dexp(K_G0+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(1)-lambda-Xea22(2)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(2))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(11)
        mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(3)-lambda+mu(1:nmes,2))*&
        (dexp(-dexp(lambda)*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(12)
            mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0-lambda-Xea22(3)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(3))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(13)
            mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(3))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(14)
            mu11(1:nmes,1)  =((dexp(y0+dexp(K_G0+Xea22(1)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(2)-lambda-Xea22(3)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(3))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        case(15)
            mu11(1:nmes,1)  =((dexp(y0+Xea22(1)+dexp(K_G0+Xea22(2)+mu(1:nmes,1))*Z1(1:nmes,1)+&
            Z1(1:nmes,2)*dexp(K_D0+Xea22(3)-lambda-Xea22(4)+mu(1:nmes,2))*&
        (dexp(-dexp(lambda+Xea22(4))*Z1(1:nmes,1))-1)))**box_cox_par-box_cox1)/box_cox_par
        
        end select
        
        
        
        
        end if



        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmes
    
                if(ycurrent(k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu11(k,1)-s_cag)/sqrt(sigmae),upper))
       
                else
                yscalar = yscalar + (ycurrent(k)-mu11(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmes
            
        yscalar = yscalar + (ycurrent(k)-mu11(k,1))**2.d0
        end do
        end if

    
        yscalar = dsqrt(yscalar)
        if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321

         vraisind = vraisind*prod_cag*dexp( -(yscalar**2.d0)/(sigmae*2.d0)&
                                    -uiiui(1)/2.d0)*det**(-0.5d0)*&
                                            (2.d0*pi)**(-nb1/2.d0)


        func10J =vraisind 

        ss = ss+product(weights(nnn,:))*func10J
    

            end do
    
     if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)

        return
    
        END SUBROUTINE gauherJ3cvplNL
    
    
    
    