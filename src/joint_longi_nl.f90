!mtaille =c(mt1,mt2,mt11,mt12)
!paraweib =c(shapeweib(1),shapeweib(2),scaleweib(1),scaleweib(2))
!kendall: 1ere colonne ss0, 2eme colonne tau
!paratps = c(timedep0,nbinnerknots,qorder0)
!counts = c(ni,cpt,cpt_dc)
!noVar = c(noVar1,noVar2,noVar3)

    
    
    
    !--entC*te pour fortran
        subroutine jointlonginl(nsujet0,nsujety0,ng0,nz0,k0,tt00,tt10,ic0,groupe0      &
        ,tt0dc0,tt1dc0,icdc0,link0,yy0,groupey0,nb0,which_random0,box_cox0,matzy0,cag0&
        ,nva10,vax0,nva20,vaxdc0,nva30,nva40,vaxy0,noVar,ag0,maxit0   &
        ,np,neta0,b,H_hessOut,HIHOut,resOut,LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out &
        ,typeof0,equidistant,mtaille &
        ,counts,ier_istop,paraweib &
       ! ,ResLongi&!Pred_y0 &
       ,ziOut ,EPS,&
    GH,paGH0, b_pred0,effet0,indic_alpha0,weights0,nodes0,nnodes_all0, RE_which0)
    

    
    !AD: add for new marq
        use donnees_indiv
        use parameters
        use splines
        use comon
        use tailles
        use lois_normales
        use optim
    !     use ParametresPourParallelisation
    !AD:pour fortran
        use sortie
        use residusM
        use comongroup
        use betatttps
    !AD:
        implicit none
    
        integer::maxit0,npinit,nvatmp,indic_alphatmp,mt1,mt2,mt11,mt12,effet0,indic_ALPHA0    !nn
        integer,dimension(4),intent(in)::mtaille
        integer,intent(in)::nsujet0,nsujety0,ng0,nz0,nva10,nva20,nva30,nva40,ag0,nb0,nnodes_all0
        integer,dimension(2),intent(in)::link0
    double precision,dimension(nz0+6),intent(out)::ziOut
        integer::np,equidistant,which_random0
        integer,dimension(2),intent(in) :: neta0
        integer,dimension(nsujet0),intent(in)::groupe0,ic0
        integer,dimension(nsujety0),intent(in) :: groupey0
        integer,dimension(ng0),intent(in)::icdc0
        double precision,dimension(2),intent(in) :: cag0,box_cox0
         integer,dimension(3),intent(in):: GH
         double precision,dimension(ng0,nb0+1+effet0+nb0+effet0 +((nb0+effet0)*(nb0+effet0-1))/2),intent(out)::b_pred0
        double precision,dimension(nnodes_all0,(nb0+effet0)),intent(in):: nodes0,weights0 
    !    integer::typeJoint0
    
        double precision,dimension(ng0)::tt0dc0,tt1dc0
        double precision,dimension(nsujet0)::tt00,tt10 !! rajout
        double precision,dimension(2)::k0
        double precision,dimension(nsujet0,nva10),intent(in):: vax0
        double precision,dimension(ng0,nva20),intent(in):: vaxdc0
        double precision,dimension(nsujety0,nva30+nva40),intent(in):: vaxy0
        double precision,dimension(nsujety0,2) :: matzy0
        double precision,dimension(nsujety0) :: yy0
        double precision,dimension(np,np)::H_hessOut,HIHOut
        double precision::resOut
        double precision,dimension(mtaille(1))::x1Out
        double precision,dimension(mtaille(2))::x2Out
        double precision,dimension(mtaille(1),3)::lamOut
        double precision,dimension(mtaille(3),3)::suOut
        double precision,dimension(mtaille(2),3)::lam2Out
        double precision,dimension(mtaille(4),3)::su2Out
        integer::ss,sss,itj
        double precision,dimension(np):: b
        double precision,dimension(2),intent(out)::LCV
        double precision,dimension(2)::shapeweib,scaleweib
        double precision,dimension(4),intent(out)::paraweib
    
        integer,dimension(4),intent(in)::noVar
        integer::noVar1,noVar2,noVar3,noVar4!! rajout
        integer::cpt,cpt_dc,ni
        integer,dimension(2),intent(out)::ier_istop
        integer,dimension(3),intent(out)::counts
        integer::groupe,groupey,ij,kk,j,k,nz,n,ii,iii,iii2,cptstr1,cptstr2   &
        ,i,ic,icdc,cptni,cptni1,cptni2,nb_echec,nb_echecor,id,cptbiais &
        ,cptauxdc,ier,istop
        double precision::tt0,tt0dc,tt1,tt1dc,h,hdc,res,min,mindc,max,pord, &
        maxdc,maxt,maxtdc,moy_peh0,moy_peh1,lrs,BIAIS_moy,mintdc !! rajout
        double precision,dimension(2)::res01
    !AD: add for new marq
        double precision::ca,cb,dd
        double precision,external::funcpajLongisplines_nl,funcpajLongiweib_nl
        double precision,external::funcpaj_tps,funcpaG_tps
        double precision,dimension(100)::xSu1,xSu2
    
        integer::typeof0
    !predictor
        !double precision,dimension(ng0)::Resmartingale,Resmartingaledc!,frailtypred!,frailtyvar
    
        !double precision,dimension(ng0,nb0+effet0)::re_pred
        !double precision,dimension(ng0,nb0+effet0+2),intent(out)::MartinGales
        !double precision,dimension(nsujety0,2):: Pred_y0
      !      double precision,dimension(nsujety0,4),intent(out):: ResLongi
      !      double precision,dimension(nsujety0) :: ResLongi_cond0,ResLongi_marg0,&
       !                    ResLongi_chol0,ResLongi_cond_st0
    
        double precision,external::funcpajres,funcpajres_log,funcpajres_triNL
        double precision,external::funcsplines_nl,funcweib_nl
        !double precision,dimension(nsujet0)::linearpred
        !double precision,dimension(ng0)::linearpreddc
        double precision,dimension(1,nva10)::coefBeta
        double precision,dimension(1,nva20)::coefBetadc
        double precision,dimension(1,nva30)::coefBetaY
        !double precision::rl_temp
        double precision,dimension(1,nsujet0)::XBeta
        !double precision,dimension(1,nsujety0)::XBetaY
        double precision,dimension(1,ng0)::XBetadc
        integer,dimension(nb0)::RE_which0
    
        integer::ngtemp
    

        double precision,dimension(3),intent(inout)::EPS ! seuils de convergence : on recupere les valeurs obtenues lors de l'algorithme a la fin
       
           double precision,dimension(ng0,nb0+1+nb0+ (nb0*(nb0-1))/2),intent(in):: paGH0
               double precision,dimension(ng0,nb0+1+nb0 + (nb0*(nb0-1))/2)::paGH
                 !double precision,dimension(:,:),allocatable :: varcov_marg_inv
            !double precision,dimension(:),allocatable :: matv
            !double precision,dimension(nva30,nva30) :: element
        
        mt1=mtaille(1)
        mt2=mtaille(2)
        mt11=mtaille(3)
        mt12=mtaille(4)

    allocate(nodes(nnodes_all0,nb0+effet0),weights(nnodes_all0,nb0+effet0))!,b_lme_nnodes(GH(2)**(nb0+effet0),ng0,nb0))
    
    nodes = nodes0
    paGH = paGH0
    weights= weights0
    nnodes_all = nnodes_all0 !GH(2)**(nb0+effet0)
    
    which_random = which_random0
   !     ResLongi_marg0 = 0.d0
    !        ResLongi_chol0 = 0.d0
     !       ResLongi_cond0 = 0.d0
      !      ResLongi_cond_st0 = 0.d0
       ! Pred_y0  =0.d0
    
            ier = ier_istop(1)
            istop = ier_istop(2)
 
        allocate(vaxdc(nva20),vaxy(nva30+nva40))
            if(nsujet0.gt.1)allocate(vax(nva10))
        
    
     
    allocate(RE_which(nb0)) 
    
    RE_which = RE_which0
    
    
        if(nsujet0.gt.1.and.nsujety0.gt.0) typeJoint = 3
        if(nsujet0.eq.1.and.nsujety0.gt.0) typeJoint = 2
    
        ag = ag0
        typeof = typeof0
        model = 1
    
    
        s_cag_id = int(cag0(1))
        s_cag = cag0(2)
        
        box_cox1 = int(box_cox0(1))
        if(box_cox1.eq.1)then 
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
    
    !     if (intcens0.eq.1) then
    !         epsa = 1.d-2
    !         epsb = 1.d-2
    !         epsd = 1.d-2
    !     endif
    !AD:end
    
        genz(1) = 30
        genz(2) = 500
    
        lrs = 0.d0
        moy_peh0 = 0.d0
        moy_peh1 = 0.d0
    
        nb_echec = 0
        nb_echecor = 0
        nb0recu = 0
        moyrecu =0.d0
    
        ngmax=ng0
        ng=ng0
    
    

            ngtemp=ng
    
    
        
                    res_ind = 0
    
    !Prise en compte des temps dc
    
    
        allocate(variable(ngtemp,nva20),ictemp(ngtemp))
    
        allocate(temps0dc(ngtemp),temps1dc(ngtemp))
    
    
            temps0dc=tt0dc0
            temps1dc=tt1dc0
            ictemp=icdc0
            variable=vaxdc0
    
    !end
        allocate(RisqCumul(ng))
        allocate(ResidusRec(ngtemp),Residusdc(ngtemp),ResidusLongi(nsujety0))!,Pred_y(nsujety0,2))
        allocate(Rrec(ngtemp),Nrec(ngtemp),Rdc(ngtemp),Ndc(ngtemp),Rdc_res(ngtemp))
    
        allocate(cdc(ngtemp))
        allocate(t0dc(ngtemp))
        allocate(t1dc(ngtemp))
        allocate(aux1(ngtemp))
        allocate(aux2(ngtemp))
        allocate(res1(ngtemp))
        !allocate(res4(ngtemp))
        allocate(res3(ngtemp))!,mi(ngtemp))
    
        allocate(nig(ng),nigdc(ng))
        shapeweib = 0.d0
        scaleweib = 0.d0
        nsujetmax=nsujet0
        if(nsujetmax.eq.1) nsujetmax = nsujety0
        nsujet=nsujet0
        nsujetymax = nsujety0
        nsujety = nsujety0

    !Al
        if(typeJoint.eq.1.or.typeJoint.eq.3) allocate(t0(nsujetmax),t1(nsujetmax),c(nsujetmax),stra(nsujetmax),g(nsujetmax))
    
        allocate(yy(nsujetymax), aux(3*nsujetmax)) !! chgt dimension aux
    
        !** number of random effects ***
         if(typeJoint.eq.0)nea = nb0
        if(typeJoint.eq.2)nea = nb0
        if(typeJoint.eq.3) nea = nb0 + 1
        nb1 = nb0
        !********* longitudinal data *************
    !    if(linkidyd0+linkidyr0.ge.1) then
        link = link0(1)
        linkidyd = link0(1)
        if(typeJoint.eq.3)linkidyr = link0(2)
        yy = yy0
        nb_re = nb0 !+ (nb0*(nb0-1))/2.d0
        netadc = neta0(1)
        netar = neta0(2)
   
        
  
        if(typeJoint.eq.3) then
            effet = 1
        else
                effet = 0
        end if
    !        if(effet.eq.1) then
    !            indic_alpha = 1
   !                 else
                            indic_ALPHA = indic_alpha0
   !         end if
    
            allocate (Ut(nea,nea),Utt(nea,nea),ziy(nsujety0,2),sum_mat(nva30,nva30),matb_chol(nea,nea),mat(nea,nea))!,&
    !    mat_all(nnodes_all*nea,nnodes_all*nea))
    
        ziy = matzy0
  
    
        allocate(vuu(nea))
        !*** find the number of recurrent measures per patient
        allocate(nmesrec(ng),nmesrec1(ng),nmesy(ng),nmes_o(ng))
        allocate(groupee(nsujet),groupeey(nsujety))
        nmesrec =1
        nmesy = 1
        
        groupeey = groupey0
        groupee = groupe0
        
        nmes_o = 0
    
        
        
        if(typeJoint.eq.1.or.typeJoint.eq.3) then
            if(groupee(1).eq.1) then
                nmesrec(1) = 1
                if(ic0(1).eq.1) nmesrec1(1) = 1
            else
                nmesrec(1) = 0
            end if
            i = 1
            do j=2,nsujet
                if(groupee(j).eq.groupee(j-1)) then
                if(ic0(j).eq.1) nmesrec1(i) = nmesrec1(i)+1
                    nmesrec(i) = nmesrec(i) + 1
                else
                    i = i+1
                end if
            end do
        end if
    
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
    
        maxmesrec=0
    
        do i=1,ng
            if (nmesrec(i).gt.maxmesrec) then
                maxmesrec=nmesrec(i)
            end if
        end do
    
    if(s_cag_id.eq.0)nmes_o = nmesy
    
        allocate(mu1_res(maxmesy),mu1(maxmesy,1))
        allocate(varcov_marg(nsujety,maxmesy))
    
    
        ndatemaxdc=2*ng0
        if (typeof == 0) then
            allocate(nt0dc(ngtemp),nt1dc(ngtemp),nt0(nsujetmax),nt1(nsujetmax))!! rajout
            allocate(mm3dc(ndatemaxdc),mm2dc(ndatemaxdc),mm1dc(ndatemaxdc),mmdc(ndatemaxdc) &
            ,im3dc(ndatemaxdc),im2dc(ndatemaxdc),im1dc(ndatemaxdc),imdc(ndatemaxdc))
        end if
    
        if(typeJoint.eq.2.or.typeJoint.eq.0)nst=1
        if(typeJoint.eq.3)nst = 2
        ni=0
    !---debut des iterations de simulations
        id=1
        cptni=0
        cptni1=0
        cptni2=0
        biais_moy=0.d0
        cptbiais=0
        cptaux=0
        cptauxdc=0
        ij=0
        kk=0
        groupe=0
        n=0
        nz=0
    !**************************************************
    !********************* prog spline****************
    
    ! effet=1
        res01(1)=0.d0
        res01(2)=0.d0
    !------------  entre non fichier et nombre sujet -----
    
        if(typeJoint.eq.3) nva1=nva10
            if(typeJoint.eq.2.or.typeJoint.eq.0) nva1 = 0
        nva2=nva20
        nva3=nva30
        nva4=nva40    
        
        noVar1 = noVar(1)
        noVar2 = novar(2)
        noVar3 = noVar(3)
        noVar4 = noVar(4)
        
        nva = nva1 +  nva2 + nva3 +nva4
        
        nvarmax=nva
        allocate(ve(nsujetmax,nvarmax),vedc(ngtemp,nvarmax),vey(nsujetymax,nvarmax))
        allocate(ve1(nsujetmax,nva10),ve2(ngtemp,nva2),ve3(nsujetymax,nva3+nva4),vet22(ngtemp))
        allocate(filtre(nva10),filtre2(nva20),filtre3(nva30),filtre4(nva40))
    
    ! AD: recurrent
        if (noVar1.eq.1) then
    !        write(*,*)'filtre 1 desactive'
            filtre=0
            nva1=0
        else
            filtre=1
        end if
    !AD:death
        if (noVar2.eq.1) then
    !        write(*,*)'filtre 2 desactive'
            filtre2=0
            nva2=0
        else
            filtre2=1
        end if
    
    
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
    
    
            nva = nva1+nva2+nva3 + nva4
    !AD:end
 
    !------------  lecture fichier -----------------------
        maxt = 0.d0
        mint = 0.d0
    
        maxtdc = 0.d0
        mintdc = 0.d0
    
        cpt = 0
        cptcens = 0
        cpt_dc = 0
        k = 0
        cptstr1 = 0
        cptstr2 = 0
        nigdc = 0
        nig=0
    
    !ccccccccccccccccccccc
    ! pour le deces
    !cccccccccccccccccccc

         
    do k = 1,ngtemp
    
            if (k.eq.1) then
                mintdc = temps0dc(k) ! affectation du min juste une fois
            endif
    
            tt0dc=temps0dc(k)
            tt1dc=temps1dc(k)
            icdc=ictemp(k)
            if(typeJoint.eq.1.or.typeJoint.eq.3)then
                groupe=groupe0(k) !!!! inutile et surtout faux car groupe0 contient les groupes
                !!!! pour les evenements et non les deces (les lignes ne sont pas les memes)
            end if
    
            do j=1,nva20
                vaxdc(j)=variable(k,j)
            enddo
    
            if(tt0dc.gt.0.d0)then
                cptauxdc=cptauxdc+1
            endif
    !------------------   deces c=1 pour donnees de survie
            if(icdc.eq.1)then
                cpt_dc = cpt_dc + 1
                cdc(k)=1
                t0dc(k) = tt0dc      !/100.d0
                t1dc(k) = tt1dc      !+0.000000001
    
                iii = 0
                iii2 = 0
                do ii = 1,nva20
                    if(filtre2(ii).eq.1)then
                        iii2 = iii2 + 1
    
                        vedc(k,iii2) = dble(vaxdc(ii))
                    endif
                end do
    
            else
    !------------------   censure a droite ou event recurr  c=0
                if(icdc.eq.0)then
                    cdc(k) = 0
                    iii = 0
                    iii2 = 0
                    do ii = 1,nva20
                        if(filtre2(ii).eq.1)then
                        iii2 = iii2 + 1
                        vedc(k,iii2) = dble(vaxdc(ii))
                        endif
                    end do
                    t0dc(k) = tt0dc
                    t1dc(k) = tt1dc
                endif
            endif
    
            if (maxtdc.lt.t1dc(k))then
                maxtdc = t1dc(k)
            endif
    
            if (mintdc.gt.t0dc(k)) then
                mintdc = t0dc(k)
            endif
        end do
    
        deallocate(temps0dc,temps1dc,ictemp,variable)
    !AD:
        if (typeof .ne. 0) then
            cens = maxtdc
        end if
    !Ad
        k = 0
        cptstr1 = 0
        cptstr2 = 0
    

    !cccccccccccccccccccccccccccccccccc
    ! pour les donnees recurrentes
    !cccccccccccccccccccccccccccccccccc
            if(typeJoint.eq.3) then
        do i = 1,nsujet     !sur les observations
    
            if (i.eq.1) then
                mint = tt00(i) ! affectation du min juste une fois
            endif
    
            tt0=tt00(i)
            tt1=tt10(i)
            ic=ic0(i)
    
                    groupe=groupe0(i)
    
    
    !------------------
            do j=1,nva10
    
                vax(j)=vax0(i,j)
            enddo
    
    
    !--------------
            if(tt0.gt.0.d0)then
                cptaux=cptaux+1
            endif
    !-----------------------------------------------------
    
    !     essai sans troncature
    !     tt0=0
    !------------------   observation c=1 pour donnees recurrentes
    
    
            if(ic.eq.1)then
                cpt = cpt + 1
                c(i)=1
                t0(i) = tt0
                t1(i) = tt1
                t1(i) = t1(i)
    
                g(i) = groupe
    
                nig(groupe) = nig(groupe)+1 ! nb d event recurr dans un groupe
                iii = 0
                iii2 = 0
    !                  do ii = 1,ver
                do ii = 1,nva10
                    if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) = dble(vax(ii)) !ici sur les observations
                    endif
                end do
    
    
            else
    !------------------   censure a droite  c=0 pour donnees recurrentes
    
            if(ic.eq.0)then
                cptcens=cptcens+1
    
                    c(i) = 0
                    iii = 0
                    iii2 = 0
            !                     do ii = 1,ver
                    do ii = 1,nva10
                        if(filtre(ii).eq.1)then
                        iii = iii + 1
                        ve(i,iii) = dble(vax(ii))
                        endif
                    end do
                    t0(i) =  tt0
                    t1(i) = tt1
                    t1(i) = t1(i)
                    g(i) = groupe
    
                endif
            endif
            if (maxt.lt.t1(i))then
                maxt = t1(i)
            endif
    
    
    
            if (mint.gt.t0(i)) then
                mint = t0(i)
            endif
        end do
    
            end if
    
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
  

        deallocate(filtre,filtre2,filtre3,filtre4)
    ! nsujet=i-1
 
        if (typeof == 0) then
            nz=nz0
            nz1=nz
            nz2=nz
    
            nzloco=nz1
            nzdc=nz2
    
            if(nz.gt.20)then
                nz = 20
            endif
            if(nz.lt.4)then
                nz = 4
            endif
        end if
    
        if(typeJoint.eq.2.or.typeJoint.eq.0) then
            ndatemax = 2*nsujety
        else
            ndatemax=2*nsujet+sum(ic0) !! rajout
        end if
    
        allocate(date(ndatemax),datedc(ndatemax))
    
        if(typeof == 0) then
            allocate(mm3(ndatemax),mm2(ndatemax) &
            ,mm1(ndatemax),mm(ndatemax),im3(ndatemax),im2(ndatemax),im1(ndatemax),im(ndatemax))
        end if
    
     
    !!!  DONNEES DECES
        aux = 0.d0
        mindc = 0.d0
        maxdc = maxtdc
    
        do i = 1,2*ngtemp
            do k = 1,ngtemp
    
                if((t0dc(k).ge.mindc))then
                    if(t0dc(k).lt.maxdc)then
                        maxdc = t0dc(k)
                    endif
                endif
                if((t1dc(k).ge.mindc))then
                    if(t1dc(k).lt.maxdc)then
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
        do i=2,2*ngtemp
            if(aux(i).gt.aux(i-1))then
                k = k+1
                datedc(k) = aux(i)
            endif
    
        end do
    
        if(typeof == 0) then
            ndatedc = k
        end if
    
    
    !--------------- zi- ----------------------------------
    
    !      construire vecteur zi (des noeuds)
    !!! DONNEES RECURRENTES
            if(typeJoint.eq.1.or.typeJoint.eq.3) then
                    min = 0.d0
                    aux = 0.d0
                    max = maxt
    
                    do i = 1,(2*nsujet+sum(ic0)) !! rajout
                            do k = 1,nsujet
                            
                                    if (t0(k).gt.min) then
                                            if (t0(k).le.max) then
                                                    max = t0(k)
                                            endif
                                    endif
                                    if (t1(k).gt.min) then
                                            if (t1(k).le.max) then
                                                    max = t1(k)
                                        endif
                                    endif
    
                            end do
            
                    
                            aux(i) = max
                            min = max + 1.d-17
                            max = maxt
                    end do
    
                    date(1) = aux(1)
                    k = 1
                    do i=2,(2*nsujet+sum(ic0))
            
                            if(aux(i).gt.aux(i-1))then
                                    k = k+1
                                    date(k) = aux(i)
                            endif
    
                    end do
            end if

     

            if(typeof == 0) then
                    if(typeJoint.eq.3) then
    ! Al:10/03/2014 emplacement des noeuds splines en percentile (sans censure par intervalle)
                            if(equidistant.eq.0) then ! percentile
                                    i=0
                                    j=0
    !----------> taille - nb de recu
                                    do i=1,nsujet
                                            if(t1(i).ne.(0.d0).and.c(i).eq.1) then
                                                    j=j+1
                                            endif
                                    end do
                                    nbrecu=j
    
    !----------> allocation des vecteur temps
                                    allocate(t2(nbrecu))
    
    !----------> remplissage du vecteur de temps
                                    j=0
                                    do i=1,nsujet
                                            if (t1(i).ne.(0.d0).and.c(i).eq.1) then
                                                    j=j+1
                                                    t2(j)=t1(i)
                                            endif
                                    end do
    
                                    nzmax=nz+3
                                    allocate(zi(-2:nzmax))
                                    
                                    ndate = k
                                    zi(-2) = mint
                                    zi(-1) = mint !date(1)
                                    zi(0) = mint !date(1)
                                    zi(1) = mint !date(1)
                                    j=0
                                    do j=1,nz-2
                                            pord = dble(j)/(dble(nz)-1.d0)
                                call percentile3(t2,nbrecu,pord,zi(j+1))
                                    end do
                                    zi(nz) = maxt !date(ndate)
                                    zi(nz+1) = maxt !zi(nz)
                                    zi(nz+2) = maxt !zi(nz)
                                    zi(nz+3) = maxt !zi(nz)
                                    ziOut = zi
                                    deallocate(t2)
                            else ! equidistant
                                    nzmax=nz+3
                                    allocate(zi(-2:nzmax))
                                
                                    ndate = k
    
                                    zi(-2) = date(1)
                                    zi(-1) = date(1)
                                    zi(0) = date(1)
                                    zi(1) = date(1)
                                    h = (date(ndate)-date(1))/dble(nz-1)
                                     do i=2,nz-1
                                            zi(i) =zi(i-1) + h
                                    end do
                                    zi(nz) = date(ndate)
                                    zi(nz+1)=zi(nz)
                                    zi(nz+2)=zi(nz)
                                    zi(nz+3)=zi(nz)
                                    ziOut = zi
                            endif
    
                    end if
    
    ! ajout : noeuds des deces
            if(equidistant.eq.0) then ! percentile
                i=0
                j=0
    !----------> taille - nb de deces
                do i=1,ngtemp !nsujet
                    if(t1dc(i).ne.(0.d0).and.cdc(i).eq.1) then
                        j=j+1
                    endif
                end do
                nbdeces=j
    
    !----------> allocation des vecteur temps
                allocate(t3(nbdeces))
    
    !----------> remplissage du vecteur de temps
                j=0
                do i=1,ngtemp !nsujet
                    if (t1dc(i).ne.(0.d0).and.cdc(i).eq.1) then
                        j=j+1
                        t3(j)=t1dc(i)
                    endif
                end do
    
                allocate(zidc(-2:(nzdc+3)))
    
                zidc(-2) = mintdc
                zidc(-1) = mintdc
                zidc(0) = mintdc
                zidc(1) = mintdc
                j=0
                do j=1,nzdc-2
                    pord = dble(j)/(dble(nzdc)-1.d0)
                    call percentile3(t3,nbdeces,pord,zidc(j+1))
                end do
                zidc(nzdc) = maxtdc
                zidc(nzdc+1) = maxtdc
                zidc(nzdc+2) = maxtdc
                zidc(nzdc+3) = maxtdc
                deallocate(t3)
            else ! equidistant
    
                allocate(zidc(-2:(nzdc+3)))
    
                zidc(-2) = datedc(1)
                zidc(-1) = datedc(1)
                zidc(0) = datedc(1)
                zidc(1) = datedc(1)
                hdc = (datedc(ndatedc)-datedc(1))/dble(nzdc-1)
    
                do i=2,nzdc-1
                    zidc(i) =zidc(i-1) + hdc
                end do
                zidc(nzdc) = datedc(ndatedc)
                zidc(nzdc+1)=zidc(nzdc)
                zidc(nzdc+2)=zidc(nzdc)
                zidc(nzdc+3)=zidc(nzdc)
            endif
            if(typeJoint.eq.2.or.typeJoint.eq.0) ziOut = zidc
    ! fin ajout
    
        end if
    
    
    !---------- affectation nt0dc,nt1dc DECES ----------------------------
    
        indictronqdc=0
        do k=1,ngtemp
            if (typeof == 0) then
                if(nig(k).eq.0)then
                    nb0recu = nb0recu + 1 !donne nb sujet sans event recu
                endif
                moyrecu =  moyrecu + dble(nig(k))
    
                if(t0dc(k).eq.0.d0)then
                    nt0dc(k) = 0
                endif
            end if
    
            if(t0dc(k).ne.0.d0)then
                indictronqdc=1
            endif
    
            if (typeof == 0) then
                do j=1,ndatedc
                    if(datedc(j).eq.t0dc(k))then
                        nt0dc(k)=j
                    endif
    
                    if(datedc(j).eq.t1dc(k))then
                        nt1dc(k)=j
                    endif
                end do
            end if
        end do
    
    !---------- affectation nt0,nt1,ntU RECURRENTS----------------------------
    
            if(typeJoint.eq.1.or.typeJoint.eq.3) then
        indictronq=0
        do i=1,nsujet
            if (typeof == 0) then
                if(t0(i).eq.0.d0)then
                    nt0(i) = 0
                endif
            end if
    
            if(t0(i).ne.0.d0)then
                indictronq=1
            endif
            if (typeof == 0) then
                do j=1,ndate
                    if(date(j).eq.t0(i))then
                        nt0(i)=j
                    endif
                
                    if(date(j).eq.t1(i))then
                        nt1(i)=j
                    endif
                end do
            end if
        end do
            end if
    
        if (typeof == 0) then
    !---------- affectation des vecteurs de splines -----------------
            n  = nz+2
    
    
    !AD:add argument:ndatedc
    
        if(typeJoint.eq.1.or.typeJoint.eq.3) then
        call vecspliJ(n,ndate,ndatedc)
        else
    
        call vecspliJNoRecur(n,ndatedc)
        end if
    
            if(typeJoint.eq.2.or.typeJoint.eq.0) then
                    nzmax = nzdc+3
                    allocate(zi(-2:nzmax))
    
                    zi(-2:nzmax) = zidc(-2:nzmax)
                    end if
    !AD:end
            allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax) &
            ,m3m1(nzmax),m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax))
    
        m3m3 = 0.d0 
    m2m2 = 0.d0 
    m1m1 = 0.d0 
    mmm = 0.d0 
    m3m2 = 0.d0 
    m3m1 = 0.d0 
    m3m = 0.d0 
    m2m1 = 0.d0 
    m2m = 0.d0 
    m1m = 0.d0
    
            call vecpenJ(n)
        end if
  
        npmax=np
    npp = np
        allocate(the1(-2:npmax),the2(-2:npmax))!,the3(-2:npmax))
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
    
      !      b(np-nva-nb_re) = 1.d0
    
        if (typeof == 2) then
            if(typeJoint.eq.2)    b(1:2)=0.8d0
            if(typeJoint.eq.3)     b(1:4)=0.8d0
        end if
    
          
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     INITIALISATION SUR SHARED Frailty MODEL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(typeof.eq.0.and.typeJoint.eq.4) then
            indic_joint=0
            indic_alphatmp = indic_alpha
            indic_alpha=0
            nvatmp=nva
            nva=nva1
            
                    npinit = nz+2+nva1+effet
            
    
            allocate(Binit(npinit))
        allocate(vvv((npinit*(npinit+1)/2)))
      
    
            nst=1
            stra=1
    !        select case(initialisation)
    !            case(1)
                    !=======> initialisation par shared
                    Binit=1.d-1 !5.d-1
                    
                        Binit((nz+2+1):(nz+2+nva1))=1.d-1
                        Binit(nz+2+nva1+effet)=1.d0
                
                
    !    write(*,*),'===================================',effet,npinit
    !    write(*,*),'== sur un SHARED FRAILTY model ====',(nz+2),nva1,nva2
    !    write(*,*),'=== donnees recurrentes uniquement ='
    !    write(*,*),'==================================='
    !     write(*,*),(Binit(i),i=1,npinit)
   
            allocate(I_hess(npinit,npinit),H_hess(npinit,npinit),v((npinit*(npinit+3)/2)))
    
     

   
                !    if (timedep.eq.0) then
                !            call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaGsplines)
                !    else
                        call marq98J(k0,Binit,npinit,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaG_tps)
                !    endif
            !
            
            deallocate(I_hess,H_hess,v)
    
    !    write(*,*)'======= Fin Shared ======='
    
                b(1:(nz+2)) = Binit(1:(nz+2))
            deallocate(Binit)
    
            
                ! b((np-nva1+1):np)=b((nz+4):(nz+3+nva1))
        
    
            
    
    !=====initialisation des splines
    
        b((nz+3):(2*(nz+2)))=b(1:(nz+2))
        
        indic_alpha = indic_alphatmp
        nst = 2
        nva = nvatmp
        deallocate(vvv)
            end if
    !############################
        ca=0.d0
        cb=0.d0
        dd=0.d0
    
    RisqCumul = 0.d0

        allocate(Z1(maxmesy,nb0),Zet(nsujetymax,2))
        allocate(mu(maxmesy,2),ycurrent(maxmesy),b1(np))!,mu_all(nsujety,nnodes_all)
    !    allocate(ycurrent_all(nnodes_all,nsujety),param_alnorm(nnodes_all,22),s_cag_all(nnodes_all,1),&
    !    upper_all(nnodes_all),aux_all(nnodes_all,ng),cdc_all(nnodes_all,ng),weights_all(nnodes_all))
        if(typeof == 0) then
            allocate(res1cur(maxmesrec),res2cur(maxmesrec),res3cur(maxmesrec))
        else
            allocate(res1cur(1),res2cur(1),res3cur(1))
        end if
        allocate(x2(maxmesy,nva3),x2cur(1,nva3),z1cur(1,nb1),current_mean(1))
        allocate(part(ngtemp))
    
        allocate(I_hess(np,np),H_hess(np,np),v((np*(np+3)/2)))
    
        if (typeof .ne. 0)allocate(vvv((np*(np+1)/2)))
        
        
              allocate(etaydc(netadc),etayr(netar),b_lme(ng,nb1),invBi_cholDet(ng),invBi_chol(ng,nea + (nea*(nea-1))/2))
            
        initGH = GH(3)
        methodGH = GH(1)

!        if(initGH.eq.1.and.GH(1).eq.1) then
!            allocate(H_hess_GH(nea,nea),I_hess_GH(nea,nea),b_paGH(ng0,nea+1+nea + (nea*(nea-1))/2))
!    
!            if(typeJoint.ge.2) then
!                allocate(vecuiRes2(ng,nea),vres(nea*(nea+3)/2),XbetaY_res(1,nsujety))
!                effetres = effet
!                XbetaY_res = XbetaY
!    
!                b1 = b
!                npp = np
!                etaydc(1:netadc) = b1(np-nva-nb_re-netadc:np-nva-nb_re-1)
!                etayr(1:netar) =b1(np-nva-nb_re-netadc - netar:np-nva-nb_re-netadc - 1)
!                
!                if(typeof.eq.0)then
!                    rl_temp=funcsplines_nl(b1,np,1,0.d0,1,0.d0,k0)
!                else if(typeof.eq.2) then 
!                    rl_temp=funcweib_nl(b1,np,1,0.d0,1,0.d0,k0)
!                end if
!    
!                if (typeJoint.eq.2) then
!                    Resmartingale = 0.d0
!                    Resmartingaledc = 0.d0
!               
!                
!       Call Residusj_biv_nl(b,np,funcpajres_biv,Resmartingaledc,ResLongi_cond0,ResLongi_cond_st0,ResLongi_marg0,&
!                                                                   ResLongi_chol0,Pred_y0,re_pred)
!                else
!
!   !     Call Residusj_tri_nl(b,np,funcpajres_tri,Resmartingale,Resmartingaledc,ResLongi_cond0,ResLongi_cond_st0,&
!    !                           ResLongi_marg0,ResLongi_chol0,Pred_y0,re_pred)
!                                                                
!                endif
!  
!    pagH(1:ng,1:nb1) = b_paGH(1:ng,1:nb1)
!    pagH(1:ng,nb1+1) = b_paGH(1:ng,nea + 1)
!    pagH(1:ng,(nb1+1):nb1+1+nb1 + (nb1*(nb1-1))/2) = b_paGH(1:ng,(nea+1):nea+1+nb1 + (nb1*(nb1-1))/2)
!    
!               
!    
!                deallocate(vecuiRes2,XbetaY_res,vres,H_hess_GH,I_hess_GH,b_paGH)
!    
!    
!    
!            end if
!    
!        
!        end if
!        end if
        
            ! Parametres pour GH pseudo-adaptative

                    do i=1,ng
                            if(methodGH.le.2) then 
                            b_lme(i,1:nb1) = paGH(i,1:nb1)
                            invBi_cholDet(i) = paGH(i,nb1+1)
                                
                            invBi_chol(i,1:nb1 + (nb1*(nb1-1))/2) = paGH(i,(nb1+2):(nb1+1+nb1 + (nb1*(nb1-1))/2))
                                  end if 
                              end do
                    
          !  methodGH = GH(1)
            nodes_number = GH(2)

     allocate( vres(nea*(nea+3)/2),H_hess_GH(nea,nea))

!    
!    open(24,file='paGH1.txt')!,status='old')
!    do i=1,ng
!    write(24,*)paGh(i,:)
!    
!    end do
!    
!    close(24)
!    stop
!
    
        !    ress= 0
        

!        end if
        !### vectorielle
    !    do ig=1,ng
    !    do i=1,nnodes_all
    !            b_lme_nnodes(i,:) = b_lme(ig,:)
    !            yy_all(
    !    end do
    !    end do
           
    
        !if (istop .ne. 1)goto 1000 ! si l'initialisation ne marche pas, ne pas faire le modele
    
   
  
            select case(typeof)
                case(0)
    !                 if (timedep.eq.0) then
    !                     if (logNormal.eq.0) then
    !                         if (intcens.eq.1) then
    !                             call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajsplines_intcens)
    !                         else
    
        
                                call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajlongisplines_nl)
    !                         endif
    !                     else
    !                         call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajsplines_log)
    !                     endif
    !                 else
    !                     call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaj_tps)
    !                 endif
                case(2)
    !                 if (timedep.eq.0) then
    !                     if (logNormal.eq.0) then
    !                         if (intcens.eq.1) then
    !                             call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajweib_intcens)
    !                         else
    
    
                                call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajlongiweib_nl)
    
    !                         endif
    !                     else
    !                         call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajweib_log)
    !                     endif
    !                 else
    !                     call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaj_tps)
    !                 endif
            end select
        
      
    deallocate( vres,H_hess_GH)
        
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
    
    
        if(effet.eq.1.and.ier.eq.-1)then
            v((np-nva-indic_alpha)*(np-nva-indic_alpha+1)/2)=10.d10
        endif
    
        res01(effet+1)=res
    
            if (typeJoint.eq.2) then
            scaleweib(1) = 0.d0
            shapeweib(1) = 0.d0
            scaleweib(2) = b(2)!
            shapeweib(2) = b(1)
        else
            scaleweib(1) = b(2) !etaR
            shapeweib(1) = b(1) !betaR
    
            scaleweib(2) = b(4) !etaD
            shapeweib(2) = b(3) !betaD
        end if
        paraweib(1) = shapeweib(1)
        paraweib(2) = shapeweib(2)
        paraweib(3) = scaleweib(1)
        paraweib(4) = scaleweib(2)
    ! --------------  Lambda and survival estimates JRG January 05
    
            xsu1 = 0.d0
            xsu2 = 0.d0
            nstRec = 1
    
            if(typeJoint.eq.2) then
            select case(typeof)
            case(0)
                call distancelongisplines(nz2,b,mt2,x2Out,lam2Out,su2Out)
            case(2)
                Call distancelongiweib(b,np,mt2,x2Out,lam2Out,xSu2,su2Out)
        end select
            else if(typeJoint.eq.3) then
                    allocate(etaT(nstRec),betaT(nstRec))
        select case(typeof)
            case(0)
                call distanceJsplines(nz1,nz2,b,mt1,mt2,x1Out,lamOut,suOut,x2Out,lam2Out,su2Out)
            case(2)
    
                Call distanceJweib(b,np,mt1,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out)
        end select
        deallocate(etaT,betaT)
                    
            end if
    
    
    
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
        if (typeof == 0) then
            !write(*,*)'The approximate like cross-validation Criterion in the non parametric case'
            call multiJ(H_hess,I_hess,np,np,np,HI)
            do i =1,np
                LCV(1) = LCV(1) + HI(i,i)
            end do
                if(typeJoint.eq.1) then
                    LCV(1) = (LCV(1) - resnonpen) / nsujet
                else if(typeJoint.eq.2) then
                    LCV(1) = (LCV(1) - resnonpen) /(ng+nsujety)
                else
                    LCV(1) = (LCV(1) - resnonpen) /(nsujet+nsujety)
                end if
            else
    !        write(*,*)'=========> Akaike information Criterion <========='
            if(typeJoint.eq.1) then
                LCV(2) = (1.d0 / nsujet) *(np - resOut)
            else if(typeJoint.eq.2) then
                LCV(2) = (1.d0 /( nsujety+ng)) *(np - resOut)
            else
                LCV(2) = (1.d0 / (nsujet+nsujety)) *(np - resOut)
            !        write(*,*)'======== AIC :',LCV(2)
            end if
        end if
    
        1000 continue
    !AD:end

                
 
    !    write(*,*)'======== Calcul des residus de martingale ========'
    
        deallocate(I_hess,H_hess)
    
        coefBeta = 0.d0
        coefBetadc = 0.d0
        coefBetaY = 0.d0
        Xbeta = 0.d0
        Xbetadc = 0.d0
        Zet = 0.d0
            if(typeJoint.eq.3) then
        do i=1,nsujet
            do j=1,nva1
                ve1(i,j)=ve(i,j)
            end do
        end do
                    end if
        do i=1,ngtemp
            do j=1,nva2
                ve2(i,j)=vedc(i,j)
            end do
        end do
    
        if(typeJoint.eq.2) then
            do i=1,nsujety
                do j=1,nva3
                    ve3(i,j) = vey(i,j)
                end do
                do j=1,2
                    Zet(i,j) = ziy(i,j)
                end do
            end do
        else
            do i=1,nsujety
                do j=1,nva3
                    ve3(i,j) = vey(i,j)
                end do
                do j=1,2
                    Zet(i,j) = ziy(i,j)
                end do
            end do
        end if
    
     
    
                    if(typeJoint.eq.2) then
                            coefBetadc(1,:) = b((np-nva+1):(np-nva+nva2))
               !             coefBetaY(1,:) = b((np-nva+nva2+1):np)
                            Xbetadc = matmul(coefBetadc,transpose(ve2))
                    !        XbetaY = matmul(coefBetaY,transpose(ve3))
                    else
                            coefBeta(1,:) = b((np-nva+1):(np-nva+nva1))
                            coefBetadc(1,:) = b((np-nva+nva1+1):(np-nva+nva1+nva2))
                !            coefBetaY(1,:) = b((np-nva+nva1+nva2+1):np)
                            Xbeta = matmul(coefBeta,transpose(ve1))
                            Xbetadc = matmul(coefBetadc,transpose(ve2))
              !              XbetaY = matmul(coefBetaY,transpose(ve3))
                    end if
    
    
            
        !deallocate(I_hess,H_hess)
        
        itj = 0
        
!        do ig=1,ng
!    
!            nmescur =nmesy(ig)
!            allocate(matv(nmescur*(nmescur+1)/2),varcov_marg_inv(nmescur,nmescur),&
!                        mat_sigma(nmescur,nmescur))
!    
!            
!                    do i= 1,nmescur
!                        mat_sigma(i,i) = sigmae**2.d0
!                     end do
!    
!     
!     varcov_marg((itj+1):(itj+nmescur),1:nmescur) =Matmul( MATMUL(ziy((itj+1):(itj+nmescur),1:nb1), &
!        MATMUL(Ut(1:nb1,1:nb1),Utt(1:nb1,1:nb1))),transpose(ziy((itj+1):(itj+nmescur),1:nb1)))+ &
!                    mat_sigma
!    
!            
!            matv = 0.d0
!            do j=1,nmescur
!                do k=j,nmescur
!                    jj=j+k*(k-1)/2
!                    matv(jj)=varcov_marg(itj+j,k)
!                end do
!            end do
!        ier = 0
!        epsj = 1.d-10
!    
!        call dsinvj(matv,nmescur,epsj,ier)
!    
!        varcov_marg_inv=0.d0
!        do j=1,nmescur
!                do k=1,nmescur
!                    if (k.ge.j) then
!                        varcov_marg_inv(j,k)=matv(j+k*(k-1)/2)
!                    else
!                        varcov_marg_inv(j,k)=matv(k+j*(j-1)/2)
!                    end if
!            end do
!        end do
!    
!        element =  Matmul(Matmul(Transpose(vey(itj+1:itj+nmescur,1:nva3)), &
!                        varcov_marg_inv(1:nmescur,1:nmescur)), vey(itj+1:itj+nmescur,1:nva3))
!    
!        do j=1,nva3
!            do k=1,nva3
!                 sum_mat(j,k) = sum_mat(j,k) +  element(j,k)
!            end do
!        end do
!    
!    
!            itj = itj + nmescur
!    
!        
!            deallocate(matv,varcov_marg_inv,mat_sigma)
!
!        end do
        
        
        
        
    
    ! if((istop == 1) ) then
    
    
    allocate(H_hess_GH(nea,nea),I_hess_GH(nea,nea))!,b_paGH(ng0,nb1+1+nb1 + (nb1*(nb1-1))/2))
    
        if((istop == 1) ) then
    
    
     allocate(vecuiRes2(ng,nea), vres(nea*(nea+3)/2))
     
    
        b1 = b
        npp = np

        Call PostRE_triNL(b,np,funcpajres_triNL, b_pred0)                        
              
    deallocate(vecuiRes2,vres)
    
    
    
            end if
    
 
    !    if(typeJoint.ge.2) then
    !            allocate(vecuiRes2(ng,nea),&
    !                    vres(nea*(nea+3)/2),&
    !                    XbetaY_res(1,nsujety))
    !            !I_hess = 0.d0
    !            !H_hess = 0.d0
    !
    !
    !            effetres = effet
    !            XbetaY_res = XbetaY
    !
    !    b1 = b
    !            npp = np
    !
    !
    !            if(typeof.eq.0)then
    !            rl_temp=funcsplines_nl(b1,np,1,0.d0,1,0.d0,k0)
    !            else if(typeof.eq.2) then 
    !            
    !            rl_temp=funcweib_nl(b1,np,1,0.d0,1,0.d0,k0)
    !            end if
    !
    !            if (typeJoint.eq.2) then
    !                Resmartingale = 0.d0
    !            Resmartingaledc = 0.d0
    !                      write(*,*)'ok3'

                

    !                if(typeof.eq.0)then
    !            rl_temp=funcsplines(b1,np,1,0.d0,1,0.d0,k0)
    !            else if(typeof.eq.2) then 
    !            
    !            rl_temp=funcweib(b1,np,1,0.d0,1,0.d0,k0)
    !            end if
                
    !   Call Residusj_biv_nl(b,np,funcpajres_biv,Resmartingaledc,ResLongi_cond0,ResLongi_cond_st0,ResLongi_marg0,&
     !                                                              ResLongi_chol0,Pred_y0,re_pred)
                                    
        
     !           else
                  
     !           Call Residusj_tri_nl(b,np,funcpajres_tri,Resmartingale,Resmartingaledc,ResLongi_cond0,ResLongi_cond_st0,&
      !                                                                             ResLongi_marg0,ResLongi_chol0,Pred_y0,re_pred)
     !           endif
  
    !b_pred0(1:ng,1:nea+1+nea + (nea*(nea-1))/2) = b_paGH(1:ng,1:nea+1+nea + (nea*(nea-1))/2)
    
    !            if (istopres.eq.1) then
    !            if(link.eq.1) then
    !                if(typeJoint.eq.3)then
    !                    do i=1,nsujet
    !                        linearpred(i)=Xbeta(1,i)+re_pred(g(i),nb1+1)+dot_product(etayr,re_pred(g(i),1:netar))
    !                    end do
    !
    !                    do i=1,ng
    !                        linearpreddc(i)=Xbetadc(1,i)+alpha*re_pred(i,nb1+1)+dot_product(etaydc,re_pred(i,1:netadc))
    !                    end do
    !                else
    !                    do i=1,ng
    !                        linearpreddc(i)=Xbetadc(1,i)+dot_product(etaydc(1:netadc),re_pred(i,1:netadc))
    !                    end do
    !                end if
    !            else 
    !                  if(typeJoint.eq.3)then
    !                
    !                    do i=1,nsujet
    !                        linearpred(i)=Xbeta(1,i)+re_pred(g(i),nb1+1)+XBetaY(1,g(i))
    !                    end do
    !
    !                    do i=1,ng
    !                        linearpreddc(i)=Xbetadc(1,i)+alpha*re_pred(i,nb1+1)+XBetaY(1,g(i))
    !                    end do
    !                else
    !                    do i=1,ng
    !                        linearpreddc(i)=Xbetadc(1,i)+XBetaY(1,g(i))
    !                    end do
    !                end if
    !            
    !            end if
    !            endif
   
    !        MartinGales(:,1)=Resmartingale
    !            MartinGales(:,2)=Resmartingaledc
    !            !Residuals(1:nsujety) = ResLongi_marg
    !    
    !            MartinGales(:,3:(2+nea))=re_pred
    !                        ResLongi(:,1) = ResLongi_cond0
    !                        ResLongi(:,2) = ResLongi_cond_st0
    !                        ResLongi(:,3) = ResLongi_marg0
    !                        ResLongi(:,4) = ResLongi_chol0
        
    
    !                linearpred=0.d0
    !                linearpreddc=0.d0
    !                        ResidusRec = 0.d0
    !                        ResidusDC = 0.d0
    !                        ResidusLongi = 0.d0
    !                        Pred_y = 0.d0
    
    !            deallocate(vecuiRes2)
    !            deallocate(XbetaY_res)
    !            deallocate(vres)
    !
    !
    !        end if
    
     
        counts(1) = ni
        counts(2) = cpt
        counts(3) = cpt_dc

    deallocate(nig,risqCumul)
    
        deallocate(cdc)
    
    deallocate(H_hess_GH,I_hess_GH)!,b_paGH)
    deallocate(t0dc)
    
    deallocate(t1dc)
    
    deallocate(aux1)!,aux2)
             
     
            deallocate(res1,res3,nigdc,etayr,etaydc)
    if(typeJoint.ne.2) deallocate(t0,t1,c,stra,g,vax)
    
        deallocate( vey,vedc,vaxy,vaxdc,aux)
    
            deallocate(ve)
        deallocate(hess,v,I1_hess,H1_hess,I2_hess,H2_hess,HI2,HIH,IH,HI,BIAIS,date,datedc)
        deallocate(ResidusRec,Residusdc,ResidusLongi,Rrec,Nrec,Rdc,Ndc,Rdc_res)!Pred_y,
        deallocate(vuu,ve1,ve2,ve3,Zet,vet22)
    
        deallocate(the1,the2)
    
    deallocate(RE_which)
       deallocate(Ut,Utt,varcov_marg,sum_mat,mat)!,mat_all)

    
        deallocate(ziy,b1,yy)
        deallocate(nmesrec,nmesrec1,nmesy,groupee,groupeey,nmes_o,mu1_res)
    !   deallocate(I3_hess,H3_hess,HI3)
    
        deallocate(Z1,mu,ycurrent,part)!,mu_all
        deallocate(res1cur,res2cur,res3cur)
    deallocate(x2,x2cur,z1cur,current_mean,mu1)
        deallocate(aux2)!,res1,res4,res3)
    !     deallocate(knotsTPS,knotsdcTPS,innerknots,innerknotsdc)
   ! deallocate(ycurrent_all,param_alnorm,s_cag_all,upper_all,cdc_all,aux_all,weights_all)
    
        if (typeof == 0) then
            deallocate(nt0dc,nt1dc,nt0,nt1)
    
        deallocate(mm3dc,mm2dc,mm1dc,&
            mmdc,im3dc,im2dc,im1dc,imdc)
    
            deallocate(mm3,mm2,mm1,mm,im3,im2,im1,im,zi,zidc,m3m3,m2m2,m1m1,mmm)
    
                deallocate(m3m2,m3m1,m3m,m2m1,m2m,m1m)
        end if
        if (typeof .ne. 0)deallocate(vvv) !,t1dc,kkapa)
        deallocate(b_lme,invBi_chol,invBi_cholDet,matb_chol)
        deallocate(nodes,weights)!,b_lme_nnodes)

        return
    
        end subroutine jointlonginl
    
    
  
        
    
    
    
        
        
            !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ1_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,betaD,etaD,t0dc,t1dc,link
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt,
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox1,box_cox_par
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ1_nl
        
        
        
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ2_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !nmesy,auxig
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,betaD,etaD,t0dc,t1dc,link
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey,typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)


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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ2_nl
        
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ3_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,betaD,etaD,t0dc,t1dc,link
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
                if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ3_nl
        
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ4_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,betaD,etaD,t0dc,t1dc,link
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ4_nl
        
    
        
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ5_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy 
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,betaD,etaD,t0dc,t1dc,link
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0+Xea22(1) + dexp(K_G0+Xea22(2) +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ5_nl
        
    
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ6_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0+Xea22(1) + dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ6_nl
        
    
                    !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ7_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0+Xea22(1) + dexp(K_G0+mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ7_nl
        
    
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ8_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0 + dexp(K_G0+Xea22(1) +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ8_nl
        
    
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ9_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0 + dexp(K_G0+Xea22(1) +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ9_nl
        
    
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ10_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0 + dexp(K_G0 +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ10_nl
        
    
        
                    !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ11_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0+Xea22(1) + dexp(K_G0+Xea22(2) +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ11_nl
        
    
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ12_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0+Xea22(1) + dexp(K_G0+Xea22(2) +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ12_nl
        
    
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ13_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0+Xea22(1) + dexp(K_G0 +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ13_nl
        
    
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ14_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0 + dexp(K_G0+Xea22(1) +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ14_nl
        
    
        
                !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ15_nl(ss,aux,res)
    
    use optim
        use tailles
        use donnees
        use comon,only:alpha,nig,cdc,sigmae,& !auxig,nmesy
            nea,nb1,etaydc,etayr,& !nva2,npp,nva3,vedc,link,betaD,etaD,t0dc,t1dc
            s_cag_id,s_cag,methodGH,b_lme,invBi_chol,matb_chol,nnodes_all,& !vey, typeof,ut,utt
            nodes,weights,yy,it,ziy,mat,det,invBi_cholDet,K_D0,K_G0,lambda,y0,&
            netar,netadc,RE_which,box_cox_par,box_cox1
        use donnees_indiv
        Implicit none
    
        double precision,intent(out)::ss
        double precision::aux,res
       
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k,n,licznik
        double precision,dimension(nea,1)::  Xea2
        double precision,dimension(nea):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        logical :: upper
        double precision,external::survdcCM
        double precision :: func10J
        double precision,parameter::pi=3.141592653589793d0
        double precision:: b_y0, b_K_G0,b_K_D0,b_lambda
         
        upper = .false.
        ss=0.d0
            
        
            do n=1,nnodes_all
       
            mu1 = 0.d0 
           i = numpat
             matb_chol = 0.d0
        
        if(methodGH.eq.1) then
         do k=1,nb1
            do j=1,k
                 matb_chol(k,j) = invBi_chol(i,j+k*(k-1)/2)
        end do
        end do
            Xea(1:nea) = nodes(n,1:nea)
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol(1:nb1,1:nb1),Xea(1:nb1))*sqrt(2.d0)
                    if(b_lme(i,1).eq.0.d0)then 
             Xea22(1:nea) = nodes(n,:)*sqrt(2.d0)
             invBi_cholDet(i)  = 1.d0
            end if
            Xea22(nb1+1) =Xea(nb1+1)*sqrt(2.d0)
              Xea2(1:nea,1) = Xea22(1:nea)
        
            end if

        b_y0 = 0.d0
        b_K_G0 = 0.d0 
        b_K_D0 = 0.d0 
        b_lambda  = 0.d0 
        
        licznik  = 1
        do k=1,nb1
        if(RE_which(k).eq.1) then 
            b_y0 = Xea22(licznik)
        else if(RE_which(k).eq.2) then 
            b_K_G0= Xea22(licznik)
        else if(RE_which(k).eq.3) then 
            b_K_D0= Xea22(licznik)
        else if(RE_which(k).eq.4) then 
            b_lambda= Xea22(licznik)
        end if
        licznik = licznik + 1
        end do
    
        uii = matmul(Xea22,mat)
    
            uiiui=matmul(uii,Xea2)
    

                mu1(1:nmescur,1)  =((dexp(y0+Xea22(1) + dexp(K_G0+Xea22(2) +mu(1:nmescur,1))*ziy((it+1):(it+nmescur),1)+&
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
       if(prod_cag.le.(1d-320))prod_cag = 1d-320

    
  
                        
       func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                            -uiiui(1)/2.d0-0.5d0*dlog(det)-dble(nea)/2.d0*dlog(2.d0*pi)&       !
                            -dexp(Xea22(nb1+1)*alpha+dot_product(etaydc,Xea22(1:netadc)))*aux&
                            + cdc(i)*dot_product(etaydc,Xea22(1:netadc))&
                            -dexp(Xea22(nb1+1)+dot_product(etayr,Xea22(1:netar)))*res &
                            +nig(i)*dot_product(etayr,Xea22(1:netar))&
                            +Xea22(nb1+1)*(nig(i)+alpha*cdc(i))! &
    
    
    

        func10J = dexp(func10J)
   


                    ss = ss+product(weights(n,:))*func10J
            end do
        
            if(methodGH.eq.1) ss = ss*invBi_cholDet(i) *2.d0**(nea/2.d0)
    if(methodGH.eq.0) ss = ss*2.d0**(nea/2.d0)
   
        return
    
        END SUBROUTINE gauherJ15_nl
        
    
        
                                    
            
        