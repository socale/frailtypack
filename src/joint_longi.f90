!mtaille =c(mt1,mt2,mt11,mt12)
!paraweib =c(shapeweib(1),shapeweib(2),scaleweib(1),scaleweib(2))
!kendall: 1ere colonne ss0, 2eme colonne tau
!paratps = c(timedep0,nbinnerknots,qorder0)
!counts = c(ni,cpt,cpt_dc)
!noVar = c(noVar1,noVar2,noVar3)

    
    
    
    !--entC*te pour fortran
        subroutine joint_longi(VectNsujet,ngnzag,k0,tt00,tt10,ic0,groupe0      &
        ,tt0dc0,tt1dc0,icdc0,link0,yy0,bb0,groupey0,groupeB0,Vectnb0,matzy0,matzB0 &
        ,cag0,VectNvar,vax0,vaxdc0,vaxy0,vaxB0,noVar,maxit0   &
        ,np,neta0,b,H_hessOut,HIHOut,resOut,LCV,x1Out,lamOut,xSu1,suOut,x2Out,lam2Out,xSu2,su2Out &
        ,typeof0,equidistant,mtaille &
        ,counts,ier_istop,paraweib &
        ,MartinGales,ResLongi,Pred_y0,GLMlog0 &
        ,positionVarTime,numInterac &
        ,linearpred,linearpreddc,ziOut &
    ,paratps,filtretps0,BetaTpsMat,BetaTpsMatDc,BetaTpsMatY,EPS,GH,paGH)
    
    !AD: add for new marq
        use donnees_indiv
        use parameters
        use splines
        use comon
        use tailles
        use lois_normales
        use optim
     use var_surrogate, only:Chol,a_deja_simul,Vect_sim_MC,graine ! add Monte-carlo
    !     use ParametresPourParallelisation
    !AD:pour fortran
        use sortie
        use residusM
        use comongroup
        use betatttps
    !AD:
        implicit none
    
    integer,dimension(3), intent(in)::VectNsujet
    integer,dimension(4), intent(in)::VectNvar

        integer::maxit0,npinit,nvatmp,indic_alphatmp,mt1,mt2,mt11,mt12 !nn
        integer,dimension(4),intent(in)::mtaille
        integer,dimension(4),intent(in)::ngnzag
        integer,dimension(2),intent(in)::Vectnb0
        integer,dimension(2),intent(in)::link0
    double precision,dimension(ngnzag(2)+6),intent(out)::ziOut
    integer, intent(in)::equidistant
        integer::np,nva10,nva20,nva30
        integer,dimension(2),intent(in) :: neta0
        integer,dimension(VectNsujet(1)),intent(in)::groupe0,ic0
        integer,dimension(VectNsujet(2)),intent(in) :: groupey0
        integer,dimension(VectNsujet(3)),intent(in) :: groupeB0 ! add TwoPaart
        integer,dimension(ngnzag(1)),intent(in)::icdc0
        double precision,dimension(2),intent(in) :: cag0
        
        !add for interaction terms in current-level association
        integer,dimension(2),intent(in):: numInterac
        integer,dimension((numInterac(1)+numInterac(2))*4),intent(in) :: positionVarTime

        integer::ng0,nz0,ag0
        
        double precision,dimension(ngnzag(1))::tt0dc0,tt1dc0
        double precision,dimension(VectNsujet(1))::tt00,tt10 !! rajout
        double precision,dimension(2)::k0
        double precision,dimension(VectNsujet(1),VectNvar(1)),intent(in):: vax0
        double precision,dimension(ngnzag(1),VectNvar(2)),intent(in):: vaxdc0
        double precision,dimension(VectNsujet(2),VectNvar(3)),intent(in):: vaxy0
        double precision,dimension(VectNsujet(3),VectNvar(4)),intent(in):: vaxB0 ! add TwoPart
        double precision,dimension(VectNsujet(2),(Vectnb0(1)-Vectnb0(2))) :: matzy0
        double precision,dimension(VectNsujet(3),Vectnb0(2)) :: matzB0
        double precision,dimension(VectNsujet(2)) :: yy0
        double precision,dimension(VectNsujet(3)) :: bb0 !add TwoPart
        double precision,dimension(np,np)::H_hessOut,HIHOut
        double precision::resOut
        double precision,dimension(mtaille(1))::x1Out
        double precision,dimension(mtaille(2))::x2Out
        double precision,dimension(mtaille(1),3)::lamOut
        double precision,dimension(mtaille(3),3)::suOut
        double precision,dimension(mtaille(2),3)::lam2Out
        double precision,dimension(mtaille(4),3)::su2Out
        integer::ss,sss,nsujet0,nsujety0
        double precision,dimension(np):: b
        double precision,dimension(2),intent(out)::LCV
        double precision,dimension(2)::shapeweib,scaleweib
        double precision,dimension(4),intent(out)::paraweib
    integer, dimension(2),intent(in)::GLMlog0 ! for glm with log link + marginal two-part
    
        integer,dimension(4),intent(in)::noVar ! modif TwoPart
        integer::noVar1,noVar2,noVar3!! rajout
        integer::cpt,cpt_dc,ni
        integer,dimension(2),intent(out)::ier_istop
        integer,dimension(3),intent(out)::counts
        integer::groupe,groupey,ij,kk,j,k,nz,n,ii,iii,iii2,cptstr1,cptstr2   &
        ,i,ic,icdc,cptni,cptni1,cptni2,nb_echec,nb_echecor,id,cptbiais &
        ,cptauxdc,p,ier,istop
        double precision::tt0,tt0dc,tt1,tt1dc,h,hdc,res,min,mindc,max,pord, &
        maxdc,maxt,maxtdc,moy_peh0,moy_peh1,lrs,BIAIS_moy,mintdc !! rajout
        double precision,dimension(2)::res01
    !AD: add for new marq
        double precision::ca,cb,dd
        double precision,external::funcpajLongisplines,funcpajLongiweib

        double precision,external::funcpaj_tps,funcpaG_tps
        double precision,dimension(100)::xSu1,xSu2
    
        integer::typeof0
    !predictor
        double precision,dimension(ngnzag(1))::Resmartingale,Resmartingaledc!,frailtypred!,frailtyvar
    
        double precision,dimension(ngnzag(1),(Vectnb0(1)+1))::re_pred
       double precision,dimension(ngnzag(1),(2+Vectnb0(1)+1)),intent(out)::MartinGales
        double precision,dimension(VectNsujet(2),2),intent(out):: Pred_y0
       double precision,dimension(VectNsujet(2),4),intent(out):: ResLongi
        double precision,dimension(VectNsujet(2)) :: ResLongi_cond0,ResLongi_marg0,&
                            ResLongi_chol0,ResLongi_cond_st0
    
        double precision,external::funcpajres,funcpajres_log,funcpajres_biv,funcpajres_tri
        double precision,dimension(VectNsujet(1)),intent(out)::linearpred
        double precision,dimension(ngnzag(1)),intent(out)::linearpreddc
        double precision,dimension(1,VectNvar(1))::coefBeta
        double precision,dimension(1,VectNvar(2))::coefBetadc
        double precision,dimension(1,VectNvar(3))::coefBetaY
        double precision,dimension(1,VectNvar(4))::coefBetaB ! add TwoPart
        double precision::coefBeta2
        double precision,dimension(1,VectNsujet(1))::XBeta
        double precision,dimension(1,VectNsujet(2))::XBetaY
        double precision,dimension(1,VectNsujet(3))::XBetaB !add TwoPart
        double precision,dimension(1,ngnzag(1))::XBetadc    
    
        integer::ngtemp
    
        integer,dimension(3),intent(in)::paratps
        integer,dimension(VectNvar(1)+VectNvar(2)+VectNvar(3)),intent(in)::filtretps0
        double precision,dimension(0:100,0:4*sum(filtretps0(1:VectNvar(1))))::BetaTpsMat !!! a refaire
        double precision,dimension(0:100,0:4*sum(filtretps0(VectNvar(1)+1:VectNvar(1)+VectNvar(2))))::BetaTpsMatDc
        double precision,dimension(0:100,0:4*sum( &
        filtretps0(VectNvar(1)+VectNvar(2)+1:VectNvar(1)+VectNvar(2)+VectNvar(3))))::BetaTpsMatY
        double precision,dimension(paratps(2)+paratps(3))::basis
        double precision,dimension(3),intent(inout)::EPS ! seuils de convergence : 
        ! on recupere les valeurs obtenues lors de l'algorithme a la fin
        integer,dimension(2),intent(in):: GH
        double precision,dimension(ngnzag(1),Vectnb0(1)+1+Vectnb0(1) + (Vectnb0(1)*(Vectnb0(1)-1))/2),intent(in):: paGH
            
        character(len=100)::bar
   
   !add TwoPart
   integer::nvaB0,groupeB,nbB0,noVarB,nsujetB0, nb0

   a_deja_simul=0 ! add Monte-carlo
   ng0=ngnzag(1)
   nz0=ngnzag(2)
   ag0=ngnzag(3)
   graine=ngnzag(4) ! seed for Monte-carlo method

   GLMloglink0=0
   GLMloglink0=GLMlog0(1) ! lambda box-cox transformation
   MTP0=0
   MTP0=GLMlog0(2)
   
    nsujet0=VectNsujet(1)
    nsujety0=VectNsujet(2)
    nsujetB0=VectNsujet(3)
            
    nb0=Vectnb0(1)
    nbB0=Vectnb0(2)
        

    nva10=VectNvar(1)        
    nva20=VectNvar(2)
    nva30=VectNvar(3)
    nvaB0=VectNvar(4)
    
    TwoPart=0
    if(nsujetB0.ne.0) TwoPart=1


        mt1=mtaille(1)
        mt2=mtaille(2)
        mt11=mtaille(3)
        mt12=mtaille(4)
    
        ResLongi_marg0 = 0.d0
            ResLongi_chol0 = 0.d0
            ResLongi_cond0 = 0.d0
            ResLongi_cond_st0 = 0.d0
        Pred_y0  =0.d0
    
            ier = ier_istop(1)
            istop = ier_istop(2)
            
                
    ! add for current-level interaction with time
    numInter = numInterac(1)
    numInterB = numInterac(2)
    !if(link0(1).ge.2) then
    if(positionVarTime(1).lt.400) then
    allocate(positionVarT((numInterac(1)+numInterac(2))*4))
     positionVarT = positionVarTime
    else
     numInter = 0
    end if
    !end if

        allocate(vaxdc(nva20),vaxy(nva30))
        if(TwoPart.eq.1) allocate(vaxB(nvaB0)) ! add TwoPart
        
            if(nsujet0.gt.1)allocate(vax(nva10))
        timedep = paratps(1)
        nbinnerknots = paratps(2)
        qorder = paratps(3)
    
        allocate(filtretps(nva10),filtre2tps(nva20),filtre3tps(nva30))
        allocate(betatps(nva10),betatps2(nva20),betatps3(nva30))
        filtretps = filtretps0(1:nva10)
        filtre2tps = filtretps0((nva10+1):(nva10+nva20))
    
        filtre3tps = filtretps0((nva10+nva20+1):(nva10+nva20+nva30))
    
        npbetatps1 = (nbinnerknots+qorder-1)*sum(filtretps)
        npbetatps2 = (nbinnerknots+qorder-1)*sum(filtre2tps)
        npbetatps3 = (nbinnerknots+qorder-1)*sum(filtre3tps)
        npbetatps = npbetatps1 + npbetatps2 + npbetatps3
    
    
    
        if(nsujet0.gt.1.and.nsujety0.gt.0) typeJoint = 3
        if(nsujet0.eq.1.and.nsujety0.gt.0) typeJoint = 2
    
        ag = ag0
        typeof = typeof0
        model = 7
    
    
        s_cag_id = int(cag0(1))
        s_cag = cag0(2)
    
        maxiter = maxit0
    !AD:add for new marq
    !write(*,*)'eps',eps
        epsa = EPS(1) !1.d-4
        epsb = EPS(2) !1.d-4
        epsd = EPS(3) !1.d-4
    
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
        allocate(ResidusRec(ngtemp),Residusdc(ngtemp),ResidusLongi(nsujety0),Pred_y(nsujety0,2))
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
        if(TwoPart.eq.1) nsujetB = nsujetB0 ! ADD TwoPart
        if(TwoPart.eq.1) nsujetBmax = nsujetB0 ! ADD TwoPart

    !Al
        if(typeJoint.ne.2) allocate(t0(nsujetmax),t1(nsujetmax),c(nsujetmax),stra(nsujetmax),g(nsujetmax))
    
        allocate(yy(nsujetymax), aux(3*nsujetmax)) !! chgt dimension aux
    
        !** number of random effects ***
        if(typeJoint.eq.2)nea = nb0
        if(typeJoint.eq.3) nea = nb0 + 1
        nb1 = nb0
        !********* longitudinal data *************
    !    if(linkidyd0+linkidyr0.ge.1) then
        link = link0(1)
        linkidyd = link0(1)
        if(typeJoint.eq.3)linkidyr = link0(2)
        yy = yy0
        if(TwoPart.eq.1) allocate(bb(nsujetBmax))!TwoPart
        if(TwoPart.eq.1) bb = bb0 ! binary values for TwoPart
        nbB = nbB0 ! number of random effects in binary part
        nb_re = nb0 + INT((nb0*(nb0-1))/2.d0)
        netadc = neta0(1)
        netar = neta0(2)
        nby = nb0-nbB !add TwoPart

            ! Parametres pour GH pseudo-adaptative
            allocate(etaydc(netadc),etayr(netar),b_lme(ng,nb1),invBi_cholDet(ng),invBi_chol(ng,nb_re))
                    do i=1,ng
                            b_lme(i,1:nb1) = paGH(i,1:nb1)
                            invBi_cholDet(i) = paGH(i,nb1+1)
                            invBi_chol(i,1:nb_re) = paGH(1,(nb1+2):(nb1+1+nb_re))
                    end do
    
            method_GH = GH(1)
            nodes_number = GH(2)
    
        if(typeJoint.eq.3) then
            effet = 1
        else
                effet = 0
        end if
            if(effet.eq.1) then
                indic_alpha = 1
                    else
                            indic_ALPHA = 0
            end if
    
                                                                                  
    
        if(TwoPart.eq.0) allocate (Ut(nea,nea),Utt(nea,nea),ziy(nsujety0,nb1),sum_mat(nva30,nva30))
        if(TwoPart.eq.1) allocate (Ut(nea,nea),Utt(nea,nea),sum_mat(nva30,nva30),ziy(nsujety0,nby))!modif TwoPart
        if(TwoPart.eq.1) allocate (ziB(nsujetB0, nbB),sum_matB(nvaB0,nvaB0)) ! add TwoPart (+modif ziy)
        if(TwoPart.eq.1) ziB = matzB0 ! binary random effects covariates ! add TwoPart 

        ziy = matzy0 ! random effects covariates matrix
    
    
        allocate(vuu(nea))
        !*** find the number of recurrent measures per patient
        allocate(nmesrec(ng),nmesrec1(ng),nmesy(ng),nmes_o(ng))
        allocate(groupee(nsujet),groupeey(nsujety))
        nmesrec =1
        nmesy = 0
        nmes_o = 0
        groupeey = groupey0
        groupee = groupe0
        allocate(nmesB(ng))
        nmesB = 0

        if(TwoPart.eq.1) allocate(groupeeB(nsujetB),nmes_oB(ng)) !add TwoPart // nmes_oB useless ?
        if(TwoPart.eq.1) groupeeB = groupeB0 !add TwoPart
    
    
        if(typeJoint.ne.2) then
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
        do j=1,nsujety
            if(groupeey(j).eq.i) then
                    nmesy(i)=nmesy(i)+1
                else
                    i = i+1
                    if(groupeey(j).eq.i) then
                        nmesy(i)=nmesy(i)+1
                    end if
                end if
        end do
    
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
    
    if(TwoPart.eq.1)then
        nmesB = 0
        i = 1
        do j=1,nsujetB
        if(groupeeB(j).eq.i) then
            nmesB(i)=nmesB(i)+1
        else
        i=i+1
            if(groupeeB(j).eq.i) then
            nmesB(i)=nmesB(i)+1
            end if
        end if
        end do

        
        
        maxmesB=0
        do i=1,ng
            if (nmesB(i).gt.maxmesB) then
                maxmesB=nmesB(i) ! looking for the maximum numnber of repeated measurement
            end if
        end do
    end if
    
        if(TwoPart.eq.1) allocate(mu1_resB(maxmesB), varcov_margB(nsujetB,maxmesB))


        allocate(mu1_res(maxmesy))
        allocate(varcov_marg(nsujety,maxmesy))
    
    
        ndatemaxdc=2*ng0
        if (typeof == 0) then
            allocate(nt0dc(ngtemp),nt1dc(ngtemp),nt0(nsujetmax),nt1(nsujetmax))!! rajout
            allocate(mm3dc(ndatemaxdc),mm2dc(ndatemaxdc),mm1dc(ndatemaxdc),mmdc(ndatemaxdc) &
            ,im3dc(ndatemaxdc),im2dc(ndatemaxdc),im1dc(ndatemaxdc),imdc(ndatemaxdc))
        end if
    
        if(typeJoint.eq.2)nst=1
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
            if(typeJoint.eq.2) nva1 = 0
        nva2=nva20
        nva3=nva30
        nvaB=nvaB0 ! number of fixed effects in binary part (TwoPart)
    
        noVar1 = noVar(1)
        noVar2 = novar(2)
        noVar3 = noVar(3)
        noVarB = noVar(4)

        nva = nva1 +  nva2 + nva3 + nvaB
        nvarmax=nva
        allocate(ve(nsujetmax,nvarmax),vedc(ngtemp,nvarmax),vey(nsujetymax,nvarmax))
        allocate(ve1(nsujetmax,nva10),ve2(ngtemp,nva2),ve3(nsujetymax,nva3))
        allocate(filtre(nva10),filtre2(nva20),filtre3(nva30))
        if(TwoPart.eq.1) allocate(ve4(nsujetBmax,nvaB), filtreB(nvaB0),veB(nsujetBmax,nvarmax))
        if(TwoPart.eq.0) allocate(filtreB(nvaB0+1))
 
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
    
    
    ! AK: longitudinal
        if (noVar3.eq.1) then
    !        write(*,*)'filtre 3 desactive'
            filtre3=0
            nva3=0
        else
            filtre3=1
        end if
    
        if (noVarB.eq.1) then ! ADD TwoPart
            filtreB=0
            nvaB=0
        else
            filtreB=1   
        end if
    
            nva = nva1+nva2+nva3+nvaB ! add TwoPart

    
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
            if(typeJoint.ne.2)then
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
        end do
    
    !cccccccccccccccccccccccccccccccccc
    ! binary data ! add TwoPart
    !cccccccccccccccccccccccccccccccccc
    if(TwoPart.eq.1) then
    vaxB = 0
        do i = 1,nsujetB     !for each longi observation
            groupeB=groupeB0(i) ! id (useful ?)
            do j=1,nvaB ! for each fixed covariate
                vaxB(j)=vaxB0(i,j) ! set values of covariate j for observation i
            end do
            iii = 0
            do ii = 1,nvaB
                if(filtreB(ii).eq.1)then
                    iii = iii + 1
                    veB(i,iii) = dble(vaxB(ii)) !ici sur les observations
                endif
            end do
        end do
    end if
    

        deallocate(filtre,filtre2,filtre3)

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
    
        if(typeJoint.eq.2) then
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
            if(typeJoint.ne.2) then
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

                        zi(-2) = mint
                        zi(-1) = mint !date(1)
                        zi(0) = mint !date(1)
                        zi(1) = mint !date(1)
                        j=0
                        do j=1,nz-2
                            pord = dble(j)/(dble(nz)-1.d0)
                !     call percentile3(t2,nbrecu,pord,zi(j+1))
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
                if(typeJoint.eq.2) ziOut = zidc
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
    
    
            if(typeJoint.ne.2) then
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
        if(typeJoint.ne.2) then
        call vecspliJL(n,ndate,ndatedc)
        else
    
        call vecspliJNoRecur(n,ndatedc)
        end if
    
            if(typeJoint.eq.2) then
                    nzmax = nzdc+3
                    allocate(zi(-2:nzmax))
    
                    zi(-2:nzmax) = zidc(-2:nzmax)
                    end if
    !AD:end
            allocate(m3m3(nzmax),m2m2(nzmax),m1m1(nzmax),mmm(nzmax),m3m2(nzmax) &
            ,m3m1(nzmax),m3m(nzmax),m2m1(nzmax),m2m(nzmax),m1m(nzmax))
    
    
            call vecpenJ(n)
        end if
    
        npmax=np
    
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
    
            b(np-nva-nb_re) = 1.d0
    
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
    
    
        allocate(Z1(maxmesy,nb0),Zet(nsujetymax,nb0))
        allocate(mu(maxmesy,1),ycurrent(maxmesy),b1(np))
        if(typeof == 0) then
            allocate(res1cur(maxmesrec),res2cur(maxmesrec),res3cur(maxmesrec))
        else
            allocate(res1cur(1),res2cur(1),res3cur(1))
        end if
        allocate(x2(maxmesy,nva3),x2cur(1,nva3),z1cur(1,nb1),current_mean(1), current_meanRaw(1))!add TwoPart
        allocate(part(ngtemp))
    
        allocate(I_hess(np,np),H_hess(np,np),v((np*(np+3)/2)))
    
        if (typeof .ne. 0)allocate(vvv((np*(np+1)/2)))
    
        !if (istop .ne. 1)goto 1000 ! si l'initialisation ne marche pas, ne pas faire le modele
    
    ! add TwoPart
    if(TwoPart.eq.1) allocate(Z1B(maxmesB,nbB),ZetB(nsujetBmax,nbB))
    if(TwoPart.eq.1) allocate(muB(maxmesB,1),Bcurrent(maxmesB))
    if(TwoPart.eq.1) allocate(XB(maxmesB,nvaB)) 

        allocate(chol(nb1,nb1)) ! Monte-carlo

     !       a_deja_simul=0 ! add Monte-carlo

        if(typeJoint.ge.2) then
            select case(typeof)
                case(0)
    !                 if (timedep.eq.0) then
    !                     if (logNormal.eq.0) then
    !                         if (intcens.eq.1) then
    !                             call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajsplines_intcens)
    !                         else
                           call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajlongisplines)
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
    
    
                                call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajlongiweib)
    
    !                         endif
    !                     else
    !                         call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpajweib_log)
    !                     endif
    !                 else
    !                     call marq98J(k0,b,np,ni,v,res,ier,istop,effet,ca,cb,dd,funcpaj_tps)
    !                 endif
            end select
! last iteration print removed for now because printed multiple times when using MPI parallel computing.
if(1.eq.0) then
            if(maxiter < 200) then
                bar(1:3) = "0%|"    
                do k=1, maxiter
                    bar(3+k:3+k)="*"
                enddo
                !bar(3+maxiter+1:3+maxiter+5) = "|100%"
                bar(maxiter+4:maxiter+4) = "|"
                do k=maxiter+5, 100
                    bar(k:k+1)=" "
                enddo
            else                
                bar(1:3) = "0%|"    
                do k=1, 70
                    bar(3+k:3+k)="*"
                enddo
                !bar(3+maxiter+1:3+maxiter+5) = "|100%"                
                bar(73:100) = "|                          "                
            endif
            
            if (ni > 300) then !Au dela de 300 iteration la bar ne peut pas s'afficher (inpr limite a 255 caracteres), 
            ! on fait le choix de laisser apparaître les iteration malgre tout
                call intpr('Iteration:', -1, ni, 1)
            else
                call intpr(bar, -1, ni, 0)
                call intpr('Iteration:', -1, ni, 1)
            endif     
        end if
end if

        resOut=res
    !Al:
        EPS(1) = ca
        EPS(2) = cb
        EPS(3) = dd
    
        ier_istop(1) = ier
        ier_istop(2) = istop
    !Al:
!    istop = 1
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
                    if(TwoPart.ne.1) then
                        LCV(1) = (LCV(1) - resnonpen) /(ng+nsujety)
                    else if(TwoPart.eq.1) then
                        LCV(1) = (LCV(1) - resnonpen) /(ng+nsujetB) ! add TwoPart
                    end if
                else
                    LCV(1) = (LCV(1) - resnonpen) /(nsujet+nsujety)
                end if
            else
    !        write(*,*)'=========> Akaike information Criterion <========='
            if(typeJoint.eq.1) then
                LCV(2) = (1.d0 / nsujet) *(np - resOut)
            else if(typeJoint.eq.2) then
                    if(TwoPart.ne.1) then
                        LCV(2) = (1.d0 /( nsujety+ng)) *(np - resOut)
                    else if(TwoPart.eq.1) then
                        LCV(2) = (1.d0 /( nsujetB+ng)) *(np - resOut)
                    end if
            else
                LCV(2) = (1.d0 / (nsujet+nsujety)) *(np - resOut)
            !        write(*,*)'======== AIC :',LCV(2)
            end if
        end if
    
        1000 continue
    !AD:end
    
                
        if (timedep.eq.1) then
        BetaTpsMat = 0.d0
            BetaTpsMatDc = 0.d0
                    BetaTpsMatY = 0.d0
    
    
    ! CALCUL DES ESTIMATIONS DES beta DEPENDANTS DU TEMPS
    
            if (nva1.gt.0) then
                call drawTimeCoef(np,b,nva1,filtretps,BetaTpsMat)
            end if
    
            if (nva2.gt.0) then
                call drawTimeCoefdc(np,b,nva2,filtre2tps,BetaTpsMatDc)
            end if
    
                    if (nva3.gt.0) then
                            call drawTimeCoefdc(np,b,nva3,filtre3tps,BetaTpsMatY)  !a faire
            end if
    
        endif
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
                do j=1,nby ! modif TwoPart
                    Zet(i,j) = ziy(i,j)
                end do
            end do
        else
            do i=1,nsujety
                do j=1,nva3
                    ve3(i,j) = vey(i,j)
                end do
                do j=1,nby ! modif TwoPart
                    Zet(i,j) = ziy(i,j)
                end do
            end do
        end if

        if(TwoPart.eq.1) then ! add TwoPart
            do i=1,nsujetB
                do j=1,nvaB
                    ve4(i,j) = veB(i,j)
                end do
                do j=1,nbB
                    ZetB(i,j) = ziB(i,j)
                end do
            end do
        end if

        if (timedep.eq.0) then
    
                    if(typeJoint.eq.2) then
                            coefBetadc(1,:) = b((np-nva+1):(np-nva+nva2))
                            coefBetaY(1,:) = b((np-nva+nva2+1):(np-nva+nva2+nva3))
                            if(TwoPart.eq.1) coefBetaB(1,:) = b((np-nva+nva2+nva3+1):(np-nva+nva2+nva3+nvaB)) ! add TwoPart
                            Xbetadc = matmul(coefBetadc,transpose(ve2))
                            XbetaY = matmul(coefBetaY,transpose(ve3))
                            if(TwoPart.eq.1) XbetaB = matmul(coefBetaB,transpose(ve4))
                    else
                            coefBeta(1,:) = b((np-nva+1):(np-nva+nva1))
                            coefBetadc(1,:) = b((np-nva+nva1+1):(np-nva+nva1+nva2))
                            coefBetaY(1,:) = b((np-nva+nva1+nva2+1):np)
                            Xbeta = matmul(coefBeta,transpose(ve1))
                            Xbetadc = matmul(coefBetadc,transpose(ve2))
                            XbetaY = matmul(coefBetaY,transpose(ve3))
                    end if
    
        else
                    if(typeJoint.eq.3)then
            do j=1,nsujet
                p=1
                call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,nbinnerknots+qorder,t1(j), &
                innerknots,boundaryknots,basis)
                do i=1,nva1
                    coefBeta2 = 0.d0
                    if (filtretps(i).eq.1) then
                        do k=-qorder+1,nbinnerknots
                            coefBeta2 = coefBeta2 + b(np-(nva+npbetatps)+p-1+k+qorder)*basis(k+qorder)
                        end do
                    else
                        coefBeta2 = b(np-(nva+npbetatps)+p)
                    endif
                    Xbeta(1,j) = Xbeta(1,j) + coefBeta2*dble(ve(j,i))
                    p=p+filtretps(i)*(nbinnerknots+qorder-1)+1
                end do
            enddo
                    end if
            do j=1,ngtemp
                p=1
                call splinebasisIndiv(qorder-1,nbinnerknots+2*qorder,nbinnerknots,nbinnerknots+qorder,t1dc(j), &
                innerknots,boundaryknots,basis)
                do i=1,nva2
                    coefBeta2 = 0.d0
                    if (filtre2tps(i).eq.1) then
                        do k=-qorder+1,nbinnerknots
                            coefBeta2 = coefBeta2 + b(np-(nva+npbetatps)+nva1+p-1+k+qorder)*basis(k+qorder)
                        end do
                    else
                        coefBeta2 = b(np-(nva+npbetatps)+nva1+p)
                    endif
                Xbetadc(1,j) = Xbetadc(1,j) + coefBeta2*dble(vedc(j,i))
                    p=p+filtre2tps(i)*(nbinnerknots+qorder-1)+1
                end do
            enddo
        endif
    
        !deallocate(I_hess,H_hess)
    
    ! if((istop == 1) ) then
    
    
    
    
        if(typeJoint.ge.2.and.TwoPart.eq.0) then ! value should be 2, modified for TwoPart models
                allocate(vecuiRes2(ng,nb1+1),&
                        vres(nea*(nea+3)/2),&
                        XbetaY_res(1,nsujety))
                
                if(TwoPart.eq.1) allocate(XbetaB_res(1,nsujetB)) !add TwoPart
                !I_hess = 0.d0
                !H_hess = 0.d0
    
    
                effetres = effet
                XbetaY_res = XbetaY
                if(TwoPart.eq.1) XbetaB_res = XbetaB

                if (typeJoint.eq.2) then
                    Resmartingale = 0.d0
                Resmartingaledc = 0.d0
    !                      write(*,*)'ok3'
    if(TwoPart.eq.0) then
        Call Residusj_biv(b,np,funcpajres_biv,Resmartingaledc,ResLongi_cond0,ResLongi_cond_st0,ResLongi_marg0,&
                                                                    ResLongi_chol0,Pred_y0,re_pred)
    else if(TwoPart.eq.1) then
            Call Residusj_biv(b,np,funcpajres_biv,Resmartingaledc,ResLongi_cond0,ResLongi_cond_st0,ResLongi_marg0,&
                                                                    ResLongi_chol0,Pred_y0,re_pred)
    end if
    
!open(2,file='C:/Users/dr/Documents/Docs pro/Docs/1_DOC TRAVAIL/2_TPJM/GIT_2019/debug.txt')  
!       write(2,*)'ping'
!     close(2)
        !    write(*,*)'ok4'
                ! re_pred = 0.d0
                re_pred(:,nb1+1) = 0.d0 ! modified (binary part not included yet)
                else
                Call Residusj_tri(b,np,funcpajres_tri,Resmartingale,Resmartingaledc,ResLongi_cond0,ResLongi_cond_st0,&
                                                                                    ResLongi_marg0,ResLongi_chol0,Pred_y0,re_pred)
                                                                                    
                endif
        linearpred=0.d0
                    linearpreddc=0.d0
                if (istopres.eq.1) then
                    if(typeJoint.eq.3)then
                        do i=1,nsujet
                            linearpred(i)=Xbeta(1,i)+re_pred(g(i),nb1+1)+dot_product(etayr,re_pred(g(i),1:nb1))
                        end do
    
                        do i=1,ng
                            linearpreddc(i)=Xbetadc(1,i)+alpha*re_pred(i,nb1+1)+dot_product(etaydc,re_pred(i,1:nb1))
                            end do
                    else
                        do i=1,ng
                            linearpreddc(i)=Xbetadc(1,i)+dot_product(etaydc,re_pred(i,1:nb1)) ! invalid read of size 8 ! removed binary part for now
                        end do
                    end if
                endif

                MartinGales(:,1)=Resmartingale
               MartinGales(:,2)=Resmartingaledc
        
                MartinGales(:,3:(3+nb1))=re_pred ! modified (binary not included yet)
                         ResLongi(1:nsujety,1) = ResLongi_cond0(1:nsujety)
                         ResLongi(1:nsujety,2) = ResLongi_cond_st0(1:nsujety)
                         ResLongi(1:nsujety,3) = ResLongi_marg0(1:nsujety)
                         ResLongi(1:nsujety,4) = ResLongi_chol0(1:nsujety)

                        
        !    open(5,file="residuals_FFCDcag.txt")
        !    do i =1,nsujety
        !    write(5,*)i,ResLongi_marg(i),Pred_ymarg(i)
        !!    write(*,*)i,ResLongi_marg(i),Pred_ymarg(i)
        !    end do
        !    close(5)
        !    open(6,file="martingales_FFCDcag.txt")
        !    do i=1,ng
        !        write(6,*)i,Resmartingaledc(i)
        !        end do
        !    close(6)
    
                
                            ResidusRec = 0.d0
                            ResidusDC = 0.d0
                            ResidusLongi = 0.d0
       !                     Pred_y = 0.d0
    
                deallocate(vecuiRes2,XbetaY_res,vres)
    if(TwoPart.eq.1) deallocate(XbetaB_res)
    
    
            end if
    !  else
    !              !deallocate(I_hess,H_hess)
    !
    !      ! les 4 variables suivantes sont maintenant dans MartinGales
    !      !         Resmartingale=0.d0
    !      !         Resmartingaledc=0.d0
    !      !         frailtypred=0.d0
    !      !         frailtyvar=0.d0
    !
    !
    !              MartinGales=0.d0
    !              linearpred=0.d0
    !              linearpreddc=0.d0
    !
    !  end if

        counts(1) = ni
        counts(2) = cpt
        counts(3) = cpt_dc
    
    deallocate(nig,risqCumul)
    
        deallocate(cdc)
    
    
    deallocate(t0dc)
    
    deallocate(t1dc)
    
    deallocate(aux1)!,aux2)
    
            deallocate(res1,res3,nigdc)
    
    !,aux2,res1,res4,res3,mi,
    if(typeJoint.ne.2) deallocate(t0,t1,c,stra,g,vax)

        deallocate( vey,vedc,vaxy,vaxdc,aux)
    deallocate(nmesB)
        if(TwoPart.eq.1) deallocate(groupeeB,vaxB, bb, ziB) ! add TwoPart

       if(TwoPart.eq.1) deallocate( nmes_oB, veB, ve4) ! add TwoPart
        if(TwoPart.eq.1) deallocate(filtreB,Z1B, ZetB, muB, Bcurrent, XB)!, z1Bcur)

        if(TwoPart.eq.1) deallocate(sum_matB, varcov_margB)
        if(TwoPart.eq.0) deallocate(filtreB)
        
        if(method_GH.eq.3) deallocate(Vect_sim_MC)

            deallocate(ve)
        deallocate(hess,v,I1_hess,H1_hess,I2_hess,H2_hess,HI2,HIH,IH,HI,BIAIS,date,datedc)
    
        deallocate(ResidusRec,Residusdc,ResidusLongi,Pred_y,Rrec,Nrec,Rdc,Ndc,Rdc_res)
        deallocate(vuu,ve1,ve2,ve3,Zet)
    
        deallocate(the1,the2)
    
        deallocate(filtretps,filtre2tps,filtre3tps)
        deallocate(betatps,betatps2,betatps3)
    
        deallocate(Ut,Utt,varcov_marg,sum_mat)
    
    
        deallocate(ziy,b1,yy)
        deallocate(nmesrec,nmesrec1,nmesy,groupee,groupeey,nmes_o,mu1_res)
        if(TwoPart.eq.1) deallocate(mu1_resB)
    
    !   deallocate(I3_hess,H3_hess,HI3)
    
        deallocate(Z1,mu,ycurrent,part)
        deallocate(res1cur,res2cur,res3cur)
    deallocate(x2,x2cur,z1cur,current_mean, current_meanRaw) !add TwoPart
    
        deallocate(aux2)!,res1,res4,res3)
    !     deallocate(knotsTPS,knotsdcTPS,innerknots,innerknotsdc)
    
    
        if (typeof == 0) then
            deallocate(nt0dc,nt1dc,nt0,nt1)
    
        deallocate(mm3dc,mm2dc,mm1dc,&
            mmdc,im3dc,im2dc,im1dc,imdc)
    
            deallocate(mm3,mm2,mm1,mm,im3,im2,im1,im,zi,zidc,m3m3,m2m2,m1m1,mmm)
    
                deallocate(m3m2,m3m1,m3m,m2m1,m2m,m1m)
        end if
    
        if (typeof .ne. 0)deallocate(vvv) !,t1dc,kkapa)
        deallocate(etaydc,etayr,b_lme,invBi_chol,invBi_cholDet)
    deallocate(chol) !Monte-carlo
    if(link.ge.2) then
    if(positionVarTime(1).lt.400) deallocate(positionVarT)
end if
        return
    
        end subroutine joint_longi
    
    
    !========================== VECSPLI ==============================
    !AD:add argument:ndatedc
        subroutine vecspliJL(n,ndate,ndatedc)
    !AD:end
        use tailles
    !AD:
        use comon,only:date,datedc,zi,mm3,mm2,mm1,mm,im3,im2,im1,im &
        ,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc
    !AD:end
        IMPLICIT NONE
    
        integer,intent(in)::n,ndate,ndatedc
        integer::i,j,k
        double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
    
    !----------  calcul de u(ti) :  STRATE1 ---------------------------
    !    attention the(1)  sont en nz=1
    !        donc en ti on a the(i)
        j=0
        do i=1,ndate-1
            do k = 2,n-2
                if ((date(i).ge.zi(k-1)).and.(date(i).lt.zi(k)))then
                    j = k-1
                endif
            end do
            ht = date(i)-zi(j)
            htm= date(i)-zi(j-1)
            h2t= date(i)-zi(j+2)
            ht2 = zi(j+1)-date(i)
            ht3 = zi(j+3)-date(i)
            hht = date(i)-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3(i) = (0.25d0*(date(i)-zi(j-3))*mm3(i))+(0.25d0*hh2 &
            *mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
            im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0) &
                +(h4*mm(i)*0.25d0)
            im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
            im(i)  = ht*mm(i)*0.25d0
    
        end do
    
    !AD: add for death
    !----------  calcul de u(ti) :  STRATE2 ---------------------------
    !    attention the(1)  sont en nz=1
    !        donc en ti on a the(i)
    
        do i=1,ndatedc-1
            do k = 2,n-2
                if ((datedc(i).ge.zi(k-1)).and.(datedc(i).lt.zi(k)))then
                    j = k-1
                endif
            end do
            ht = datedc(i)-zi(j)
            htm= datedc(i)-zi(j-1)
            h2t= datedc(i)-zi(j+2)
            ht2 = zi(j+1)-datedc(i)
            ht3 = zi(j+3)-datedc(i)
            hht = datedc(i)-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3dc(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2dc(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1dc(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mmdc(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3dc(i) = (0.25d0*(datedc(i)-zi(j-3))*mm3dc(i))+(0.25d0*hh2 &
            *mm2dc(i))+(0.25d0*h3m*mm1dc(i))+(0.25d0*h4*mmdc(i))
            im2dc(i) = (0.25d0*hht*mm2dc(i))+(h3m*mm1dc(i)*0.25d0) &
                +(h4*mmdc(i)*0.25d0)
            im1dc(i) = (htm*mm1dc(i)*0.25d0)+(h4*mmdc(i)*0.25d0)
            imdc(i)  = ht*mmdc(i)*0.25d0
    
        end do
    !AD:end
        end subroutine vecspliJL
    
            !========================== VECSPLI NO RECUR ==============================
    !AD:add argument:ndatedc
        subroutine vecspliJNoRecur(n,ndatedc)
    !AD:end
        use tailles
    !AD:
        use comon,only:datedc,zidc,mm3dc,mm2dc,&
         mm1dc,mmdc,im3dc,im2dc,im1dc,imdc
        !,date,im,im1,im3,im2,mm3,mm2,mm1,mm
    !AD:end
        IMPLICIT NONE
    
        integer,intent(in)::n,ndatedc
        integer::i,j,k
        double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
    
    !AD: add for death
    !----------  calcul de u(ti) :  STRATE2 ---------------------------
    !    attention the(1)  sont en nz=1
    !        donc en ti on a the(i)
    
        j=0
        do i=1,ndatedc-1
            do k = 2,n-2
    
                if ((datedc(i).ge.zidc(k-1)).and.(datedc(i).lt.zidc(k)))then
                    j = k-1
                endif
            end do
    
            ht = datedc(i)-zidc(j)
            htm= datedc(i)-zidc(j-1)
            h2t= datedc(i)-zidc(j+2)
            ht2 = zidc(j+1)-datedc(i)
            ht3 = zidc(j+3)-datedc(i)
            hht = datedc(i)-zidc(j-2)
            h = zidc(j+1)-zidc(j)
            hh= zidc(j+1)-zidc(j-1)
            h2= zidc(j+2)-zidc(j)
            h3= zidc(j+3)-zidc(j)
            h4= zidc(j+4)-zidc(j)
            h3m= zidc(j+3)-zidc(j-1)
            h2n=zidc(j+2)-zidc(j-1)
            hn= zidc(j+1)-zidc(j-2)
            hh3 = zidc(j+1)-zidc(j-3)
            hh2 = zidc(j+2)-zidc(j-2)
    
            mm3dc(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
                    mm2dc(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1dc(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mmdc(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3dc(i) = (0.25d0*(datedc(i)-zidc(j-3))*mm3dc(i))+(0.25d0*hh2 &
            *mm2dc(i))+(0.25d0*h3m*mm1dc(i))+(0.25d0*h4*mmdc(i))
            im2dc(i) = (0.25d0*hht*mm2dc(i))+(h3m*mm1dc(i)*0.25d0) &
                +(h4*mmdc(i)*0.25d0)
            im1dc(i) = (htm*mm1dc(i)*0.25d0)+(h4*mmdc(i)*0.25d0)
            imdc(i)  = ht*mmdc(i)*0.25d0
    
        end do
    
    !AD:end
        end subroutine vecspliJNoRecur
    
    
    !=============================================================
    ! gauss hermite
    ! func est l integrant, ss le resultat de l integrale sur -infty , +infty
    SUBROUTINE gauherJ21(ss,choix,nnodes)
    
        use tailles
        use donnees
        use comon,only:nb1,typeJoint,nea,method_GH,invBi_cholDet!auxig,typeof
        use donnees_indiv,only : frailpol,frailpol2,frailpol3,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix,nnodes
        double precision::auxfunca,func8J,func9J,func10J,func11J, funcTP4J,funcG ! add TwoPart
        external::func8J,func9J,func10J,func11J, funcTP4J,funcG ! add TwoPart
        integer::j
            double precision,dimension(nnodes):: xx1,ww1
    
            if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
            else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
            else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
            else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
            else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
            else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
            else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
            end if

        ss=0.d0
            auxfunca = 0.d0
            if(method_GH.eq.0) then

            do j=1,nnodes
                if (choix.eq.3) then
                        if(typeJoint.eq.2.and.nb1.eq.1) then
    
                            auxfunca=funcG(0.d0,0.d0,0.d0,0.d0,xx1(j))
                    else if(typeJoint.eq.2.and.nb1.eq.2) then
                            auxfunca = funcG(0.d0,0.d0,0.d0,frailpol,xx1(j))
                    else if(typeJoint.eq.2.and.nb1.eq.3) then
                            auxfunca = funcG(0.d0,0.d0,frailpol2,frailpol,xx1(j))
                    else if(typeJoint.eq.2.and.nb1.eq.4) then
                            auxfunca = funcTP4J(frailpol3,frailpol2,frailpol,xx1(j)) ! add TwoPart
                    else if(typeJoint.eq.3.and.nea.eq.2) then
                            auxfunca = func8J(frailpol,xx1(j))
    
    
                    else if(typeJoint.eq.3.and.nea.eq.3) then
                            auxfunca = func9J(frailpol2,frailpol,xx1(j))
                    else if(typeJoint.eq.3.and.nea.eq.4) then
                            auxfunca = func11J(frailpol3,frailpol2,frailpol,xx1(j))
                  
                    endif
                    ss = ss+ww1(j)*(auxfunca)
                endif
    
            end do
 else
            
                
                        do j=1,nnodes
                if (choix.eq.3) then
                        if(typeJoint.eq.2.and.nb1.eq.1) then
    
                            auxfunca=funcG(0.d0,0.d0,0.d0,0.d0,xx1(j))
                    else if(typeJoint.eq.2.and.nb1.eq.2) then
                            auxfunca = funcG(0.d0,0.d0,0.d0,frailpol,xx1(j))
                    else if(typeJoint.eq.2.and.nb1.eq.3) then
                            auxfunca = funcG(0.d0,0.d0,frailpol2,frailpol,xx1(j))
                    else if(typeJoint.eq.3.and.nea.eq.2) then
                            auxfunca = func8J(frailpol,xx1(j))
                    else if(typeJoint.eq.3.and.nea.eq.3) then
                            auxfunca = func9J(frailpol2,frailpol,xx1(j))
                    else if(typeJoint.eq.3.and.nea.eq.4) then
                            auxfunca = func11J(frailpol3,frailpol2,frailpol,xx1(j))
                  
                    endif
                    ss = ss+ww1(j)*(auxfunca)
                endif
    
            end do
    
    
              if(typeJoint.eq.2.and.nb1.eq.1) ss = ss*invBi_cholDet(numpat)*2.d0**(1.d0/2.d0)
            endif
    
        return
    
        END SUBROUTINE gauherJ21
    
    
        !***********************************
        !********* Gauss-Hermit pour la dimension 2*
        !*************************************
    
        SUBROUTINE gauherJ22(ss,choix,nnodes)
    
        use tailles
        use donnees
        use comon,only:method_GH,invBi_cholDet,nea,typeJoint!auxig,typeof,nb1,nea
        use donnees_indiv,only : frailpol,numpat
   !      !$ use OMP_LIB
   Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix,nnodes
    !  double precision :: frail2
        double precision::auxfunca
        integer::j
            double precision,dimension(nnodes):: xx1,ww1
    
                    if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
            else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
            else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
            else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
            else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
            else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
            else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
            end if
    
        ss=0.d0
            if(method_GH.eq.0) then
!  !$OMP PARALLEL DO default(none) PRIVATE (j,auxfunca,frailpol)& 
!  !$OMP  shared(nnodes,xx1,ww1,choix,typeJoint)&
!   !$OMP  REDUCTION(+:ss)  SCHEDULE(Dynamic,1)
            do j=1,nnodes
            !  if (choix.eq.3) then
                frailpol = xx1(j)
                    call gauherJ21(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
            !    endif
            end do
!   !$OMP END PARALLEL DO          

    
            ss = ss
            else
!  !$OMP PARALLEL DO default(none) PRIVATE (j,auxfunca,frailpol)& 
!  !$OMP  shared(nnodes,xx1,ww1,choix,typeJoint)&
!   !$OMP  REDUCTION(+:ss)  SCHEDULE(Dynamic,1)
                    do j=1,nnodes
            !  if (choix.eq.3) then
                    frailpol = xx1(j)
                    call gauherJ21(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
            end do
!  !$OMP END PARALLEL DO          

            if(typeJoint.eq.2.and.nea.eq.2)ss = ss*invBi_cholDet(numpat) *2.d0**(nea/2.d0)
    
            end if
                        
        return
    
        END SUBROUTINE gauherJ22
    
    
       !***********************************
        !********* Gauss-Hermit pour la dimension 3 (bivarie)*
        !*************************************
    
        SUBROUTINE gauherJ23(ss,choix,nnodes)
    
        use tailles
        use donnees
        use comon,only:method_GH,invBi_cholDet,nea!auxig,typeof,nb1,typeJoint,nea
        use donnees_indiv,only : frailpol2,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix,nnodes
    !  double precision :: frail2
        double precision::auxfunca
        integer::j
            double precision,dimension(nnodes):: xx1,ww1
    
                    if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
            else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
            else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
            else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
            else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
            else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
            else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
            end if
    
        ss=0.d0
            if(method_GH.eq.0) then
    
            do j=1,nnodes
            !  if (choix.eq.3) then
                frailpol2 = xx1(j)
                    call gauherJ22(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
            !    endif
            end do
    
    
            ss = ss
            else
                    do j=1,nnodes
            !  if (choix.eq.3) then
                    frailpol2 = xx1(j)
                    call gauherJ22(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
            end do
    
            ss = ss*invBi_cholDet(numpat)*2.d0**(nea/2.d0)
    
            end if
        return
    
        END SUBROUTINE gauherJ23
    
    
    
        !***********************************
        !********* Gauss-Hermit pour la dimension 4 (bivarie)*
        !*************************************
    
        SUBROUTINE gauherJ24(ss,choix,nnodes)
    
        use tailles
        use donnees
        use comon,only:method_GH,invBi_cholDet,nea!auxig,typeof,nb1,typeJoint,nea
        use donnees_indiv,only : frailpol3,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix,nnodes
    !  double precision :: frail2
        double precision::auxfunca
        integer::j
            double precision,dimension(nnodes):: xx1,ww1
    
                    if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
            else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
            else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
            else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
            else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
            else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
            else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
            end if
    
        ss=0.d0
            if(method_GH.eq.0) then
    
            do j=1,nnodes
            !  if (choix.eq.3) then
                frailpol3 = xx1(j)
                    call gauherJ23(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
            !    endif
            end do
    
    
            ss = ss
            else
            do j=1,nnodes
            !  if (choix.eq.3) then
                    frailpol3 = xx1(j)
                    call gauherJ23(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
            end do
    
            ss = ss*invBi_cholDet(numpat)*2.d0**(nea/2.d0)
    
            end if
        return
    
        END SUBROUTINE gauherJ24
    
    
            !***********************************
        !********* Gauss-Hermit pour la dimension 2 - modC(le trviarie b_10, v*
        !*************************************
    
        SUBROUTINE gauherJ31(ss,choix,nnodes)
    
        use tailles
        use donnees
        use comon,only:method_GH,invBi_cholDet,nea!auxig,typeof,nb1,typeJoint,nea
        use donnees_indiv,only : frailpol,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix,nnodes
    !    double precision :: frail2
        double precision::auxfunca
        integer::j
            double precision,dimension(nnodes):: xx1,ww1
    
                            if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
            else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
            else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
            else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
            else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
            else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
            else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
            end if
    
    
        ss=0.d0
            if(method_GH.eq.0) then
    
            do j=1,nnodes
                if (choix.eq.3) then
                    frailpol = xx1(j)
                call gauherJ21(auxfunca,choix,nnodes)
    
                    ss = ss+ww1(j)*(auxfunca)
                endif
            end do
    
            else
            do j=1,nnodes
                if (choix.eq.3) then
                frailpol = xx1(j)
                call gauherJ21(auxfunca,choix,nnodes)
                ss = ss+ww1(j)*(auxfunca)
            endif
            end do
            ss = ss*invBi_cholDet(numpat)*2.d0**(nea/2.d0)
    
            end if
        return
    
        END SUBROUTINE gauherJ31
    
    
                !***********************************
        !********* Gauss-Hermit pour la dimension 3 - modC(le trviarie b_10,b_11, v
        !*************************************
    
        SUBROUTINE gauherJ32(ss,choix,nnodes)
    
        use tailles
        use donnees
        use comon,only:method_GH,invBi_cholDet,nea!auxig,typeof,nb1,typeJoint,nea
        use donnees_indiv,only : frailpol2,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix,nnodes
    !    double precision :: frail2
        double precision::auxfunca
        integer::j
            double precision,dimension(nnodes):: xx1,ww1

                            if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
            else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
            else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
            else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
            else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
            else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
            else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
            end if
    
        ss=0.d0
            if(method_GH.eq.0) then
    
            do j=1,nnodes
                if (choix.eq.3) then
                frailpol2 = xx1(j)
                call gauherJ22(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
                endif
            end do
    
            else
    
                do j=1,nnodes
                frailpol2 = xx1(j)
                call gauherJ22(auxfunca,choix,nnodes)
    
                ss = ss+ww1(j)*(auxfunca)
                
                end do
                    ss = ss*invBi_cholDet(numpat) *2.d0**(nea/2.d0)
    
            end if
        return
    
        END SUBROUTINE gauherJ32
    
           !***********************************
        !********* Gauss-Hermit pour la dimension 4 - modC(le trviarie b_10,b_11,b_12, v
        !*************************************
    
        SUBROUTINE gauherJ33(ss,choix,nnodes)
    
        use tailles
        use donnees
        use comon,only:method_GH,invBi_cholDet,nea!auxig,typeof,nb1,typeJoint,nea
        use donnees_indiv,only : frailpol3,numpat
        Implicit none
    
        double precision,intent(out)::ss
        integer,intent(in)::choix,nnodes
    !    double precision :: frail2
        double precision::auxfunca
        integer::j
            double precision,dimension(nnodes):: xx1,ww1
    
                            if(nnodes.eq.5) then
                    xx1(1:nnodes) = x5(1:nnodes)
                    ww1(1:nnodes) = w5(1:nnodes)
            else if (nnodes.eq.7) then
                    xx1(1:nnodes) = x7(1:nnodes)
                    ww1(1:nnodes) = w7(1:nnodes)
            else if (nnodes.eq.9) then
                    xx1(1:nnodes) = x9(1:nnodes)
                    ww1(1:nnodes) = w9(1:nnodes)
            else if (nnodes.eq.12) then
                    xx1(1:nnodes) = x12(1:nnodes)
                    ww1(1:nnodes) = w12(1:nnodes)
            else if (nnodes.eq.15) then
                    xx1(1:nnodes) = x15(1:nnodes)
                    ww1(1:nnodes) = w15(1:nnodes)
            else if (nnodes.eq.20) then
                    xx1(1:nnodes) = x2(1:nnodes)
                    ww1(1:nnodes) = w2(1:nnodes)
            else if (nnodes.eq.32) then
                    xx1(1:nnodes) = x3(1:nnodes)
                    ww1(1:nnodes) = w3(1:nnodes)
            end if
    
    
        ss=0.d0
            if(method_GH.eq.0) then
    
            do j=1,nnodes
                if (choix.eq.3) then
                frailpol3 = xx1(j)
                call gauherJ32(auxfunca,choix,nnodes)
                    ss = ss+ww1(j)*(auxfunca)
                endif
            end do
    
            else
    
                do j=1,nnodes
                frailpol3 = xx1(j)
                call gauherJ32(auxfunca,choix,nnodes)
    
                ss = ss+ww1(j)*(auxfunca)
                
                end do
                    ss = ss*invBi_cholDet(numpat) *2.d0**(nea/2.d0)
    
            end if
        return
    
        END SUBROUTINE gauherJ33
    
    
    !=====================================================================
    ! pour la loi log-normale
        double precision function func6JL(frail)
    
        use tailles
        use comongroup,only:vet2!vet
        use optim
        use comon,only:aux1,cdc,sigmae,nmesy,nb1,&
            nva2,npp,nva3,vedc,betaD,etaD,t1dc,etaydc,link,t0dc,&
            vey,typeof,s_cag_id,s_cag,ut,method_GH,b_lme,invBi_chol
            !auxig,alpha,sig2,res1,res3,nb1,nea,nig,utt,
        use donnees_indiv
        IMPLICIT NONE
    
        double precision,intent(in)::frail
        double precision :: yscalar,alnorm,prod_cag
        integer :: j,i,k
        logical :: upper
        double precision,external::survdcCM
        double precision :: resultdc,abserr,resabs,resasc
        double precision, dimension(nb1)::Xea
        double precision,parameter::pi=3.141592653589793d0
    
        upper = .false.
        i = numpat
    
            if(method_GH.eq.1) then
            Xea(1) = b_lme(i,1) +invBi_chol(i,1)*frail*sqrt(2.d0)
            else
            Xea(1) = frail!*sqrt(2.d0)*ut(1,1)
            !Xea(1) = frail*sqrt(2.d0)
            end if
    
    
!        if(nmesy(numpat).gt.0) then
!            allocate(mu1(nmesy(numpat),1))
!        else
!            allocate(mu1(1,1))
!        end if
    
        if(nmescur.gt.0) then
            mu1(1:nmescur,1) = mu(1:nmescur,1) +Xea(1)*Z1(1:nmescur,1)
        else
            mu1(1:nmescur,1)  = mu(1:nmescur,1)
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
    
            if(link.eq.1) then !******** Random Effects **********************
    ! pour le calcul des integrales / pour la survie, pas les donno?=es recurrentes:
    
    
    
            if(typeof.eq.2) then
                aux1(numpat)=((t1dc(numpat)/etaD)**betaD)*vet2!*dexp(etaydc1*frail)
            else if(typeof.eq.0) then
                aux1(numpat)=ut2cur*vet2!*dexp(etaydc1*frail)
    
            end if
        !end if
        else !********** Current Mean ****************
    
    
    
    
            !if(typeof.eq.2) then
            !    aux1(i)=((t1dc(i)/etaD)**betaD)*vet2*dexp(etaydc1*current_mean(1))
        call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,Xea)
    
        aux1(i) = resultdc
    !         if(aux1(i).ge.1.d0) write(*,*)i,aux1(i),Xea(1)
    
        !    else if(typeof.eq.0) then
        !        aux1(i)=ut2cur*vet2*dexp(etaydc1*current_mean(1))
        !    end if
        X2cur(1,1) = 1.d0
            X2cur(1,2) =t1dc(numpat)
            if((nva3-2).gt.0) then
                do k=3,nva3
                        X2cur(1,k) = dble(vey(it_cur+1,k))
                    end do
            end if
    
            Z1cur(1,1) = 1.d0
            current_mean = 0.d0
    
            current_mean(1) =dot_product(X2cur(1,1:nva3),b1((npp-nva3+1):npp))+Z1cur(1,1)*Xea(1)
    
    
        end if
    
            if ((aux1(numpat).ne.aux1(numpat)) ) then!.or.(abs(aux1(i)).ge. 1.d30)) then
    !          print*,'ok 5'
    !           write(*,*)aux1(numpat),vet2
            end if
    
    
        yscalar = 0.d0
        !********* Left-censoring ***********
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(ycurrent(k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
                !(0.5d0*(1.d0-erf((mu1(k)-s_cag)/(sigmae*dsqrt(2.d0)))))
    
                    else
                yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                end if
            end do
        else
            do k=1,nmescur
        yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
        end do
        end if
    
        if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321
    
        yscalar = dsqrt(yscalar)
    
            if(link.eq.1) then
        func6JL =   dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    - (Xea(1)**2.d0)/(2.d0*ut(1,1)**2) &
                    - dlog(ut(1,1))-dlog(2.d0*pi)/2.d0&
                        -aux1(numpat)*dexp(etaydc(1)*Xea(1))  + cdc(numpat)*etaydc(1)*Xea(1)
        else
        func6JL =   dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                        - (Xea(1)**2.d0)/(2.d0*ut(1,1)**2)&
                    - dlog(ut(1,1))-dlog(2.d0*pi)/2.d0&
                        -aux1(numpat) &
                        + cdc(numpat)*etaydc(1)*current_mean(1)
    
        
        end if
    
    func6JL = dexp(func6JL)
    
    
!        deallocate(mu1)
        return
    
        end function func6JL
    
    
        !***********************************************
        !****************** func7J **********************
        !***********************************************
        double precision function func7J(frail2,frail)
        use optim
    
        use tailles
        use comongroup,only:vet2!vet
        use comon,only:aux1,cdc,sigmae,nmesy,&
            nva2,npp,nva3,vedc,nb1,betaD,etaD,t0dc,t1dc,etaydc,link,&
            vey, typeof,s_cag_id,s_cag,ut,utt,method_GH,b_lme,invBi_chol
            !auxig,alpha,sig2,res1,res3,nig,nea
        use donnees_indiv
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2
        double precision :: yscalar,eps,finddet,det,alnorm,prod_cag
        integer :: j,i,jj,k,ier
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(2,1)::  Xea2
        double precision,dimension(2):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        double precision,dimension(2,2)::mat,matb_chol
        logical :: upper
        double precision,external::survdcCM
        double precision :: resultdc,abserr,resabs,resasc
        double precision,parameter::pi=3.141592653589793d0
    
        upper = .false.
!        if(nmesy(numpat).gt.0) then
!            allocate(mu1(nmesy(numpat),1))
!        else
!            allocate(mu1(1,1))
!        end if
    
        i = numpat
        matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
        if(method_GH.eq.1) then
            Xea(1) = frail
            Xea(2) = frail2
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            else
            Xea2(1,1) = frail!
        Xea2(2,1) = frail2!
        Xea22(1) = frail!
        Xea22(2) = frail2!
            end if
        mat = matmul(ut,utt)
    
    
    
    
        jj=0
    ! jjj = 0
        do j=1,2
        do k=j,2
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
    
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
    
    
        if(nmescur.gt.0) then
            mu1(1:nmescur,1) = mu(1:nmescur,1) +MATMUL(Z1(1:nmescur,1:nb1),Xea2(1:2,1))
        else
            mu1(1:nmescur,1)  = mu(1:nmescur,1)
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
    
                if(link.eq.1) then !************ Random Effects ****************
    
            if(typeof.eq.2) then
                aux1(numpat)=((t1dc(numpat)/etaD)**betaD)*vet2
            else if(typeof.eq.0) then
                aux1(numpat)=ut2cur*vet2
            end if
        !end if
    
        else !******* Current Level ***************
    
    
    
            call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,xea22)
    
            aux1(i) = resultdc
    
                X2cur(1,1) = 1.d0
            X2cur(1,2) = t1dc(numpat)
            if(nva3.gt.2) then
                do k=3,nva3
                    X2cur(1,k) = dble(vey(it_cur+1,k))
                end do
            end if
    
    
            Z1cur(1,1) = 1.d0
            Z1cur(1,2) = t1dc(numpat)
    
    
                current_mean = 0.d0
                current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Matmul(Z1cur,Xea22)
    
        end if
    
        
    
        !********* Left-censoring ***********
    
    
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(ycurrent(k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
         
                    else
                yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
        yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2.d0
        end do
        end if
    
        !yscalar = norm2(ycurrent(1:nmescur) - mu1(1:nmescur))
    
        yscalar = dsqrt(yscalar)
        if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321
        !if(prod_cag.le.(1d-320))prod_cag = 1d-320

        if(link.eq.1) then
        func7J = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dlog(2.d0*pi)&
                        -aux1(numpat)*dexp(etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))&
                        + cdc(numpat)*(etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))
    
        else
        func7J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                        -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dlog(2.d0*pi)&
                        -aux1(numpat)&
                        + cdc(numpat)*(etaydc(1)*current_mean(1))
        end if
        
        func7J = dexp(func7J)
!        deallocate(mu1)
    
        return
    
        end function func7J
    
    
         !***********************************************
        !****************** func10J **********************
        !***********************************************
        double precision function func10J(frail3,frail2,frail)
        use optim
    
        use tailles
        use comongroup,only:vet2!vet
        use comon,only:aux1,cdc,sigmae,nmesy,&
            nva2,npp,nva3,vedc,nb1,betaD,etaD,t0dc,t1dc,etaydc,link,&
            vey, typeof,s_cag_id,s_cag,ut,utt,method_GH,b_lme,invBi_chol,&
            numInter,numInterB,positionVarT,&
            nbB, nby, nvaB, nmesB,TwoPart,veB! add TwoPart
        use donnees_indiv
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2,frail3
        double precision :: yscalar,eps,finddet,det,alnorm,prod_cag
        integer :: j,i,jj,k,ier
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(3,1)::  Xea2
        double precision,dimension(3):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        double precision,dimension(3,3)::mat,matb_chol
        logical :: upper
        double precision,external::survdcCM
        double precision :: resultdc,abserr,resabs,resasc
        double precision,parameter::pi=3.141592653589793d0
        double precision :: Bscalar ! add TwoPart
        double precision,dimension(1) :: Bcv,Bcurrentvalue, cmY
        double precision, dimension(1,nb1)::z1YcurG
        double precision, dimension(1,nb1)::z1BcurG
                double precision, dimension(1,nvaB)::x2BcurG

        integer::counter, counter2
        
        upper = .false.
        
        !add TwoPart
    if(TwoPart.eq.1) then
        if(nmesB(numpat).gt.0) then
            allocate(mu1B(nmesB(numpat),1))
        else
            allocate(mu1B(1,1))
        end if
    end if
    
!        if(nmesy(numpat).gt.0) then
!            allocate(mu1(nmesy(numpat),1))
!        else
!            allocate(mu1(1,1))
!        end if
    
        i = numpat
        matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
            matb_chol(3,1) =  invBi_chol(i,4)
            matb_chol(3,2) =  invBi_chol(i,5)
            matb_chol(3,3) =  invBi_chol(i,6)
        if(method_GH.eq.1) then
            Xea(1) = frail
            Xea(2) = frail2
            Xea(3) = frail3
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            else
            Xea2(1,1) = frail!
        Xea2(2,1) = frail2!
        Xea2(3,1) = frail3!
        Xea22(1) = frail!
        Xea22(2) = frail2
        Xea22(3) = frail3
            end if
        mat = matmul(ut,utt)
    
 
    
    
    
        jj=0
    ! jjj = 0
        do j=1,3
        do k=j,3
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
    
        end do
        end do
        ier = 0
        eps = 1.d-10
    
    
            call dsinvj(matv,nb1,eps,ier)
     
    
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
    
    
        if(nmescur.gt.0) then
            mu1(1:nmescur,1) = mu(1:nmescur,1) +MATMUL(Z1(1:nmescur,1:nb1),Xea2(1:3,1))
        else
            mu1(1:nmescur,1)  = mu(1:nmescur,1)
        end if
    
            ! add TwoPart
        if(TwoPart.eq.1) then
            if(nmescurB.gt.0) then 
                mu1B(1:nmescurB,1) = muB(1:nmescurB,1) +MATMUL(Z1B(1:nmescurB,1:nbB),Xea2(nby+1:nb1,1))
            else
                mu1B(1:nmescurB,1)  = muB(1:nmescurB,1)
            end if
        end if
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    
        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
            vet2 =vet2 + b1(npp-nva3-nva2-nvaB+j)*dble(vedc(numpat,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
                if(link.eq.1) then !************ Random Effects ****************
    
            if(typeof.eq.2) then
                aux1(numpat)=((t1dc(numpat)/etaD)**betaD)*vet2
            else if(typeof.eq.0) then
                aux1(numpat)=ut2cur*vet2
            end if
      
    
        else !******* Current Level ***************
    
               counter=0
    
        counter2=1
    
            call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,xea22)
    
            aux1(i) = resultdc
    
if((nva3-1).gt.0) then
    X2cur(1,1) = 1.d0
    do k=2,nva3
    X2cur(1,k) = dble(vey(it_cur+1,k))
    end do
    if(numInter.ge.1) then
    do counter = 1,numInter ! compute time and interactions at t1dc
    X2cur(1,positionVarT(counter2+1)) =t1dc(numpat) ! time effect
    X2cur(1,positionVarT(counter2+2)) =t1dc(numpat)*dble(vey(it_cur+1,positionVarT(counter2)))! interaction
    counter2=counter2+3
    end do
    end if
end if
    
    if(TwoPart.eq.1) then
    if((nvaB-1).gt.0) then
    X2BcurG(1,1) = 1.d0
    do k=2,nvaB
    X2BcurG(1,k) = dble(veB(it_curB+1,k))
    end do
    if(numInterB.ge.1) then
    do counter = 1,numInter ! compute time and interactions at t1dc
    X2BcurG(1,positionVarT(counter2+1)) =t1dc(numpat)! time effect
    X2BcurG(1,positionVarT(counter2+2)) =t1dc(numpat)*dble(veB(it_curB+1,positionVarT(counter2)))! interaction
    counter2=counter2+3
    end do
    end if
end if
end if
    
          !  Z1cur(1,1) = 1.d0
          !  Z1cur(1,2) = t1dc(numpat)
          !  Z1cur(1,3) = 1.d0
    
                current_mean = 0.d0
                
                
                
                                   if(TwoPart.eq.1) then

                        z1YcurG(1,1) = 1.d0 ! random intercept only for now
                        z1YcurG(1,2) = t1dc(numpat)!z1YcurG(1,2) = t1dc(numpat)
                        z1YcurG(1,3) = 0.d0
                        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
                        z1BcurG(1,2) = 0.d0
                        z1BcurG(1,3) = 1.d0
                        Bcurrentvalue=0.d0
                        Bcv=0.d0

                        Bcv=MATMUL(X2BcurG,b1((npp-nvaB+1):npp))+Matmul(z1BcurG,Xea22)
                        Bcurrentvalue=dexp(Bcv)/(1+dexp(Bcv))
                    

        cmY = (MATMUL(x2cur,b1((npp-nva3-nvaB+1):(npp-nvaB)))+Matmul(z1YcurG,Xea22))

        current_mean = cmY*Bcurrentvalue
        
                    else if(TwoPart.eq.0) then
                            current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Matmul(Z1cur,Xea22)

                    end if    
        end if
    
        !********* Left-censoring ***********
    
    
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(ycurrent(k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
      
                    else
                yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
        yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2.d0
        end do
        end if
  
          !add TwoPart
    Bscalar=0.d0
    if(TwoPart.eq.1) then
        do k=1,nmescurB
            Bscalar = Bscalar + (Bcurrent(k)*mu1B(k,1)+dlog(1-(dexp(mu1B(k,1))/(1+dexp(mu1B(k,1))))))
        end do
    end if
    
        yscalar = dsqrt(yscalar)
        if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321
        

      
        if(link.eq.1) then
        func10J = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                  -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -3.d0/2.d0*dlog(2.d0*pi)&   !-(nb1/2.d0)*dlog(det*2.d0*pi)&
                     +Bscalar& ! add TwoPart
                   -aux1(numpat)*dexp(dot_product(etaydc,Xea22(1:nb1)))&
                    + cdc(numpat)*dot_product(etaydc,Xea22(1:nb1))
    
        else
        func10J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                        -uiiui(1)/2.d0-(nb1/2.d0)*dlog(det*2.d0*pi)&
                        +Bscalar& ! add TwoPart
                        -aux1(numpat)&
                        + cdc(numpat)*(etaydc(1)*current_mean(1))
        end if
      
        func10J = dexp(func10J)
   
!        deallocate(mu1)
 if(TwoPart.eq.1) deallocate(mu1B) ! add TwoPart

        return
    
        end function func10J                    
    
    
    
            !#################################################################################### add TwoPart (4 correlated random effects)
            double precision function funcTP4J(frail4,frail3,frail2,frail)
        use optim
    
        use tailles
        use comongroup,only:vet2!vet
        use comon,only:aux1,cdc,sigmae,nmesy,&
            nva2,npp,nva3,vedc,nb1,betaD,etaD,t0dc,t1dc,etaydc,link,&
            vey, typeof,s_cag_id,s_cag,ut,utt,method_GH,b_lme,invBi_chol,&
            nbB, nby, nvaB, nmesB,TwoPart,veB! add TwoPart
        use donnees_indiv
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2,frail3,frail4
        double precision :: yscalar,eps,finddet,det,alnorm,prod_cag
        integer :: j,i,jj,k,ier
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(4,1)::  Xea2
        double precision,dimension(4):: uii, Xea22,Xea
        double precision,dimension(1)::uiiui
        double precision,dimension(4,4)::mat,matb_chol
        logical :: upper
        double precision,external::survdcCM
        double precision :: resultdc,abserr,resabs,resasc
        double precision,parameter::pi=3.141592653589793d0
        double precision :: Bscalar ! add TwoPart
        double precision,dimension(1) :: Bcv,Bcurrentvalue, cmY
                double precision, dimension(1,nb1)::z1YcurG
        double precision, dimension(1,nb1)::z1BcurG
                double precision, dimension(1,nvaB)::x2BcurG

    ! for each node on quadrature, we compute the integrand (it is then multiplied by weight in gauher subroutine)

        upper = .false.
 !       if(nmesy(numpat).gt.0) then
 !           allocate(mu1(nmesy(numpat),1))
 !       else
 !           allocate(mu1(1,1))
 !       end if
    
    !add TwoPart
    if(TwoPart.eq.1) then
        if(nmesB(numpat).gt.0) then
            allocate(mu1B(nmesB(numpat),1))
        else
            allocate(mu1B(1,1))
        end if
    end if
    
        i = numpat
        matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
            matb_chol(3,1) =  invBi_chol(i,4)
            matb_chol(3,2) =  invBi_chol(i,5)
            matb_chol(3,3) =  invBi_chol(i,6)
            matb_chol(4,1) =  invBi_chol(i,7)
            matb_chol(4,2) =  invBi_chol(i,8)
            matb_chol(4,3) =  invBi_chol(i,9)
            matb_chol(4,4) =  invBi_chol(i,10)
           




           
        if(method_GH.eq.1) then ! if Pseudo-adaptive
            Xea(1) = frail
            Xea(2) = frail2
            Xea(3) = frail3
            Xea(4) = frail4
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            else ! Standard
            Xea2(1,1) = frail! values defined by quadrature
        Xea2(2,1) = frail2!
        Xea2(3,1) = frail3!
        Xea2(4,1) = frail4
        Xea22(1) = frail!
        Xea22(2) = frail2
        Xea22(3) = frail3
        Xea22(4) = frail4
            end if
        mat = matmul(ut,utt)
    
 
    
        jj=0
    ! jjj = 0
        do j=1,4
        do k=j,4
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
    
        end do
        end do
        ier = 0
        eps = 1.d-10
    
    
            call dsinvj(matv,nb1,eps,ier)
     
    
        mat=0.d0
        do j=1,4
                do k=1,4
                            if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                end do
    
    
        uii = matmul(Xea22,mat)
            det = finddet(matmul(ut,utt),4)
    
            uiiui=matmul(uii,Xea2)
    
    ! mu : covariates values*parameters for each observations
    ! Z1 : covariates values for random effects
    ! mu1 : mu + random effects
        if(nmescur.gt.0) then
            mu1(1:nmescur,1) = mu(1:nmescur,1) +MATMUL(Z1(1:nmescur,1:nby),Xea2(1:nby,1))
        else
            mu1(1:nmescur,1)  = mu(1:nmescur,1)
        end if
    
        ! add TwoPart
        if(TwoPart.eq.1) then
            if(nmescurB.gt.0) then 
                mu1B(1:nmescurB,1) = muB(1:nmescurB,1) +MATMUL(Z1B(1:nmescurB,1:nbB),Xea2(nby+1:nb1,1))
            else
                mu1B(1:nmescurB,1)  = muB(1:nmescurB,1)
            end if
        end if

       
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
    
        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
            vet2 =vet2 + b1(npp-nva3-nva2-nvaB+j)*dble(vedc(numpat,j)) ! modif TwoPart
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
        if(link.eq.1) then !************ Random Effects ****************
    
            if(typeof.eq.2) then
                aux1(numpat)=((t1dc(numpat)/etaD)**betaD)*vet2
            else if(typeof.eq.0) then
                aux1(numpat)=ut2cur*vet2
            end if
      
    
        else !******* Current Level *************** # TO DO
    
    
    
            call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,xea22)
    
            aux1(i) = resultdc
    
                x2cur(1,1) = 1.d0
            x2cur(1,2) = t1dc(numpat)
            if(nva3.gt.2) then ! modif TwoPart
                do k=3,nva3
                    x2cur(1,k) = dble(vey(it_cur+1,k))
                end do
            end if

            if(TwoPart.eq.1) then
                            X2BcurG(1,1) = 1.d0
                    X2BcurG(1,2) = t1dc(numpat)
                    if(nvaB.gt.2) then
                        do k=3,nvaB
                            X2BcurG(1,k) = dble(veB(it_curB+1,k))
                        end do
                    end if
            end if
            

            Z1cur(1,1) = 1.d0
            Z1cur(1,2) = t1dc(numpat)
            
                current_mean = 0.d0
                    if(TwoPart.eq.1) then

                        z1YcurG(1,1) = 1.d0
                        z1YcurG(1,2) = t1dc(numpat)
                        z1YcurG(1,3) = 0.d0
                        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
                        z1BcurG(1,2) = 0.d0
                        z1BcurG(1,3) = 1.d0
                        Bcurrentvalue=0.d0
                        Bcv=0.d0
                        Bcv=MATMUL(X2BcurG,b1((npp-nvaB+1):npp))+Matmul(z1BcurG,Xea22)
                        Bcurrentvalue=dexp(Bcv)/(1+dexp(Bcv))
                    
        cmY = (MATMUL(x2cur,b1((npp-nva3-nvaB+1):(npp-nvaB)))+Matmul(z1YcurG,Xea22))

        current_mean = cmY*Bcurrentvalue


        
                    else if(TwoPart.eq.0) then
                            current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Matmul(Z1cur,Xea22)

                    end if
        end if
    
            if ((aux1(numpat).ne.aux1(numpat)) ) then
    
  
            end if
    
        !********* Left-censoring ***********
    
    
    
        yscalar = 0.d0
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(ycurrent(k).le.s_cag) then
                prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
      
                    else
                yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2.d0
                end if
            end do
        else
            do k=1,nmescur
        yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2.d0
        end do
        end if
        ! yscalar=(Y*-Y)^2
  
    !add TwoPart
    Bscalar=0.d0
    if(TwoPart.eq.1) then
        do k=1,nmescurB
            Bscalar = Bscalar + (Bcurrent(k)*mu1B(k,1)+dlog(1-(dexp(mu1B(k,1))/(1+dexp(mu1B(k,1))))))
        end do
    end if
 
     ! yscalar=yscalar+Bscalar
      
        yscalar = dsqrt(yscalar)
        if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321

        if(link.eq.1) then 
        funcTP4J = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                  -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -3.d0/2.d0*dlog(2.d0*pi)&   !-(nb1/2.d0)*dlog(det*2.d0*pi)&
                   +Bscalar& ! add TwoPart
                   -aux1(numpat)*dexp(dot_product(etaydc,Xea22(1:nb1)))&
                    + cdc(numpat)*dot_product(etaydc,Xea22(1:nb1))
    
        else
        funcTP4J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                        -uiiui(1)/2.d0-(nb1/2.d0)*dlog(det*2.d0*pi)&
                   +Bscalar& ! add TwoPart
                        -aux1(numpat)&
                        + cdc(numpat)*(etaydc(1)*current_mean(1))
        end if

        funcTP4J = dexp(funcTP4J)
   
!        deallocate(mu1)
        if(TwoPart.eq.1) deallocate(mu1B) ! add TwoPart
        return
    
        end function funcTP4J
    
    
    
    
                    !***********************************************
            !****************** func8J **********************
            !***********************************************
            double precision function func8J(frail,frail2)
    use optim
    
        use tailles
        use comongroup,only:vet,vet2
        use comon,only:alpha,res1,res3,aux1,nig,cdc,sigmae,nmesy,&
            nva2,nva1,npp,nva3,vedc,ve,nb1,nea,betaD,etaD,t1dc,etaydc,etayr,&
            t0,t1,betaR,etaR,typeof,vey,c,link,s_cag_id,s_cag,ut,utt,&
            t0dc,method_GH,b_lme,invBi_chol!effet,nb_re,nva
        use donnees_indiv
    
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2
            double precision :: yscalar,eps,finddet,det,res22,alnorm,prod_cag
            integer :: j,i,jj,k,ier,ii
            double precision,dimension(nea*(nea+1)/2)::matv
            double precision,dimension(nea,1)::  Xea2
            double precision,dimension(nea):: uii, Xea22,Xea
            double precision,dimension(1)::uiiui
            double precision,dimension(nea,nea)::mat
            logical :: upper
        double precision,external::survdcCM,survRCM
            double precision :: resultdc,abserr,resabs,resasc
            double precision :: resultR
            double precision,dimension(1):: current_meanR
            double precision,parameter::pi=3.141592653589793d0
    
    
            upper = .false.
                    i = numpat
                    if(method_GH.eq.1) then
            Xea(1) = b_lme(i,1) +invBi_chol(i,1)*frail2*sqrt(2.d0)
            Xea(2) = frail
            else
            Xea(1) = frail!2!*sqrt(2.d0)*ut(1,1)
            Xea(2) = frail2
            end if
    
!            if(nmesy(numpat).gt.0) then
!                    allocate(mu1(nmesy(numpat),1))
!            else
!                    allocate(mu1(1,1))
!            end if
    
    
    
            Xea2(1:2,1) = Xea(1:2)
    
            Xea22(1:2) = Xea(1:2)
    
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
    
            call dsinvj(matv,nea,eps,ier)
    
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
    
    
    
    
            if(nmescur.gt.0) then
                    mu1(1:nmescur,1) = mu(1:nmescur,1) +Z1(1:nmescur,1)*Xea22(1)
            else
                    mu1(1:nmescur,1)  = mu(1:nmescur,1)
            end if
            res1(i) = 0.d0
            res3(i) = 0.d0
            res22 = 0.d0
    
            !ccccccccccccccccccccccccccccccccccccccccc
            ! pour les recurrences
            !ccccccccccccccccccccccccccccccccccccccccc
    current_meanR = 0.d0
                if(link.eq.1) then
    
            if(typeof.eq.2) then
                    do ii=it_rec,it_rec+nmescurr-1
    
                            if(nva1.gt.0)then
                                    vet = 0.d0
                                    do j=1,nva1
                                            vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(ii,j))
                                    end do
                                    vet = dexp(vet)
                            else
                                    vet=1.d0
                            endif
    
    
                            res1(i) = res1(i)+((t1(ii)/etaR)**betaR)*vet!*dexp(etayr1*Xea22(1))
                            res3(i) = res3(i)+((t0(ii)/etaR)**betaR)*vet!*dexp(etayr1*Xea22(1))
    
    
                    end do
                    else if(typeof.eq.0) then
    
                            do ii=1,nmescurr
    
                            if(nva1.gt.0)then
                                    vet = 0.d0
                                    do j=1,nva1
                                            vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(it_rec+ii-1,j))
                                    end do
                                    vet = dexp(vet)
                            else
                                    vet=1.d0
                            endif
    
                            res1(i) = res1(i)+res1cur(ii)*vet!*dexp(etayr1*Xea22(1))
                            res3(i) = res3(i)+res3cur(ii)*vet!*dexp(etayr1*Xea22(1))
    
      
                    end do
                    end if
    
            else if(link.eq.2) then
    
            
                    do ii=it_rec,it_rec+nmescurr-1
            resultR = 0.d0
                            call integrationdc(survRCM,t0(ii),t1(ii),resultR,abserr,resabs,resasc,ii,b1,npp,xea22)
    
                    res1(i) = res1(i) + resultR     !c'est deja res1-res3
    
    
                            if((c(ii).eq.1))then
                                    X2cur(1,1) = 1.d0
                                    X2cur(1,2) = t1(ii)
                                    if(nva3.gt.2) then
                                            do k=3,nva3
                                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                                            end do
                                    end if
    
                                    Z1cur(1,1) = 1.d0
                                    current_meanR = current_meanR + MATMUL(X2cur,b1((npp-nva3+1):npp))+Z1cur(1,1)*Xea22(1)
                            end if
                    end do
    
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
    
                            if(link.eq.1) then !************ Random Effects ****************
    
                    if(typeof.eq.2) then
                            aux1(numpat)=((t1dc(numpat)/etaD)**betaD)*vet2!*dexp(etaydc1*Xea22(1)+etaydc2*Xea22(2))
                    else if(typeof.eq.0) then
                            aux1(numpat)=ut2cur*vet2!*dexp(etaydc1*Xea22(1)+etaydc2*Xea22(2))
                    end if
            !end if
    
    
            else !******* Current Level ***************
    
    
                    resultdc = 0.d0
                    call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,xea22)
    
                    aux1(i) = resultdc
    
                            X2cur(1,1) = 1.d0
                    X2cur(1,2) = t1dc(numpat)
                    if(nva3.gt.2) then
                            do k=3,nva3
                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                            end do
                    end if
    
    
                    Z1cur(1,1) = 1.d0
    
    
                            current_mean = 0.d0
                            current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Matmul(Z1cur,Xea22(1:nb1))
    
            end if
    
            if ((aux1(numpat).ne.aux1(numpat)) ) then!.or.(abs(aux1(i)).ge. 1.d30)) then
    
    
        !    write(*,*)aux1(numpat),vet2
            end if
    
    
            yscalar = 0.d0
                    prod_cag = 1.d0
            if(s_cag_id.eq.1)then
                    do k = 1,nmescur
                            if(ycurrent(k).le.s_cag) then
    
         
    
                            prod_cag = prod_cag*(1.d0-alnorm((mu1(k,1)-s_cag)/sqrt(sigmae),upper))
            !               mu1(k) = ycurrent(k)
                            else
                            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
                            end if
                    end do
            else
                    do k=1,nmescur
            yscalar = yscalar + (ycurrent(k)-mu1(k,1))**2
            end do
            end if
    
            !yscalar = norm2(ycurrent(1:nmescur) - mu1(1:nmescur))
    
    
            func8J= 0.d0
    
            yscalar = dsqrt(yscalar)
    
            if(link.eq.1) then
            func8J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                       -uiiui(1)/2.d0-0.5d0*dlog(det)&
                       -dlog(2.d0*pi)&
                      -dexp(Xea22(2)*alpha+etaydc(1)*Xea22(1))*aux1(numpat) + cdc(numpat)*(etaydc(1)*Xea22(1))&
                      -dexp(Xea22(2)+etayr(1)*Xea22(1))*(res1(i)-res3(i)) +nig(i)*etayr(1)*Xea22(1) &
                                            +Xea22(2)*(nig(i)+alpha*cdc(i))
            else if(link.eq.2) then
                    func8J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                             -uiiui(1)/2.d0-0.5d0*dlog(det)&
                             -dlog(2.d0*pi)&
                            -dexp(Xea22(2)*alpha)*aux1(numpat)&
                             + cdc(numpat)*(etaydc(1)*current_mean(1))&
                             -dexp(Xea22(2))*res1(i)+etayr(1)*current_meanR(1)  &
                                            +Xea22(2)*(nig(i)+alpha*cdc(i))
            end if
    
        func8J = dexp(func8J)
        
!    deallocate(mu1)
    
        return
    
        end function func8J
    
    
                            !***********************************************
            !****************** func9J **********************
            !***********************************************
            double precision function func9J(frail,frail2,frail3)
    use optim
    
        use tailles
        use comongroup,only:vet,vet2
        use comon,only:alpha,res1,res3,aux1,nig,cdc,sigmae,nmesy,&
            nva2,nva1,npp,nva3,vedc,ve,nea,nb1,betaD,etaD,t1dc,etaydc,etayr,&
            t0,t1,betaR,etaR,typeof,link,vey,c,s_cag_id,s_cag,&
            ut,utt,t0dc,t1dc,method_GH,b_lme,invBi_chol!effet,nva
        use donnees_indiv
            IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2,frail3
            double precision :: yscalar,eps,finddet,det,res22,alnorm,prod_cag
            integer :: j,i,jj,k,ier,ii
            double precision,dimension(nea*(nea+1)/2)::matv
            !double precision,dimension(2) :: vec,vec2
            !double precision,dimension(2):: p2,p1
            double precision,dimension(nea,1)::  Xea2
            double precision,dimension(nea):: uii, Xea22,Xea
            double precision,dimension(1)::uiiui
            double precision,dimension(nea,nea)::mat
            logical :: upper
            double precision,external::survdcCM,survRCM
            double precision :: resultdc,abserr,resabs,resasc
            double precision :: resultR
            double precision,dimension(1):: current_meanR
            double precision,parameter::pi=3.141592653589793d0
            double precision,dimension(2,2)::matb_chol
    
            upper = .false.
 
!            if(nmesy(numpat).gt.0) then
!                    allocate(mu1(nmesy(numpat),1))
!            else
!                    allocate(mu1(1,1))
!            end if
  
            i = numpat
            matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) = invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
    
        if(method_GH.eq.1) then
            Xea(1) = frail2
            Xea(2) = frail3
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea(1:2))*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            Xea22(3) = frail
            Xea2(3,1) = frail
            else
                    Xea2(1,1) = frail2
            Xea2(2,1) = frail3
            Xea2(3,1) = frail
            Xea22(1) = frail2
            Xea22(2) = frail3
            Xea22(3) = frail
            end if
    
    
    
            mat = matmul(ut,utt)
    
            jj=0
    do j=1,nea
        do k=j,nea
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
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
    
                    uii = matmul(Xea22,mat)
                    det = finddet(matmul(ut,utt),nea)
    
    
                    uiiui=matmul(uii,Xea2)
    
            if(nmescur.gt.0) then
                    mu1(1:nmescur,1) = mu(1:nmescur,1) +MATMUL(Z1(1:nmescur,1:nb1),Xea22(1:nb1))
            else
                    mu1(1:nmescur,1)  = mu(1:nmescur,1)
            end if
    
                    res1(i) = 0.d0
            res3(i) = 0.d0
            res22 = 0.d0
            !ccccccccccccccccccccccccccccccccccccccccc
            ! pour les recurrences
            !ccccccccccccccccccccccccccccccccccccccccc
            current_meanR = 0.d0
            if(link.eq.1) then
    
            if(typeof.eq.2) then
                    do ii=it_rec,it_rec+nmescurr-1
    
                            if(nva1.gt.0)then
                                    vet = 0.d0
                                    do j=1,nva1
                                            vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(ii,j))
                                    end do
                                    vet = dexp(vet)
                            else
                                    vet=1.d0
                            endif
    
    
                            res1(i) = res1(i)+((t1(ii)/etaR)**betaR)*vet!*dexp(etayr1*Xea22(1))
                            res3(i) = res3(i)+((t0(ii)/etaR)**betaR)*vet!*dexp(etayr1*Xea22(1))
    
    
                    end do
                    else if(typeof.eq.0) then
    
                            do ii=1,nmescurr
    
                            if(nva1.gt.0)then
                                    vet = 0.d0
                                    do j=1,nva1
                                            vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(it_rec+ii-1,j))
                                    end do
                                    vet = dexp(vet)
                            else
                                    vet=1.d0
                            endif
    
                            res1(i) = res1(i)+res1cur(ii)*vet!*dexp(etayr1*Xea22(1))
                            res3(i) = res3(i)+res3cur(ii)*vet!*dexp(etayr1*Xea22(1))
    
                    end do
                    end if
    
            else if(link.eq.2) then
                
                    do ii=it_rec,it_rec+nmescurr-1
            resultR = 0.d0
    
                            call integrationdc(survRCM,t0(ii),t1(ii),resultR,abserr,resabs,resasc,ii,b1,npp,xea22)
    
                    res1(i) = res1(i) + resultR     !c'est déjà res1-res3
    
    
    
                            if((c(ii).eq.1))then
                                    X2cur(1,1) = 1.d0
                                    X2cur(1,2) = t1(ii)
                                    if(nva3.gt.2) then
                                            do k=3,nva3
                                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                                            end do
                                    end if
    
                                    Z1cur(1,1) = 1.d0
                                    Z1cur(1,2) = t1(ii)
                                    current_meanR = current_meanR + MATMUL(X2cur,b1((npp-nva3+1):npp))&
                                            +MATMUL(Z1cur,Xea22(1:2))
                            end if
                    end do
    
    
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
    
            if(link.eq.1) then
    
                    if(typeof.eq.2) then
                            aux1(numpat)=((t1dc(numpat)/etaD)**betaD)*vet2!*dexp(etaydc1*Xea22(1)+etaydc2*Xea22(2))
                    else if(typeof.eq.0) then
                            aux1(numpat)=ut2cur*vet2!*dexp(etaydc1*Xea22(1)+etaydc2*Xea22(2))
                    end if
    
            else !******* Current Level ***************
    
    
                    resultdc = 0.d0
                    call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,xea22)
    
                    aux1(i) = resultdc
    
                            X2cur(1,1) = 1.d0
                    X2cur(1,2) = t1dc(numpat)
                    if(nva3.gt.2) then
                            do k=3,nva3
                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                            end do
                    end if
    
    
                    Z1cur(1,1) = 1.d0
                    Z1cur(1,2) = t1dc(numpat)
    
                            current_mean = 0.d0
                            current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Matmul(Z1cur,Xea22(1:nb1))
    
            end if
    
    
    
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
    
            !yscalar = norm2(ycurrent(1:nmescur) - mu1(1:nmescur))
    
    
            yscalar = dsqrt(yscalar)
            func9J = 0.d0
            if(link.eq.1) then
            func9J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                        -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -3.d0/2.d0*dlog(2.d0*pi)&       !
                        -dexp(Xea22(3)*alpha+etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))*aux1(numpat)&
                        + cdc(numpat)*(etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))&
                        -dexp(Xea22(3)+etayr(1)*Xea22(1)+etayr(2)*Xea22(2))*(res1(i)-res3(i)) &
                         +nig(i)*(etayr(1)*Xea22(1)+etayr(2)*Xea22(2) )&
                          +Xea22(3)*(nig(i)+alpha*cdc(i))! &
            else
            func9J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                                            -uiiui(1)/2.d0-0.5d0*dlog(det)&
                                            -3.d0/2.d0*dlog(2.d0*pi)&
                                            -dexp(Xea22(3)*alpha)*aux1(numpat)&
                                            + cdc(numpat)*(etaydc(1)*current_mean(1))&
                                            -dexp(Xea22(3))*res1(i)+etayr(1)*current_meanR(1)    &
                                            +Xea22(3)*(nig(i)+alpha*cdc(i))
            end if
    
                            func9J = dexp(func9J)
    
!            deallocate(mu1)
    
        return
    
        end function func9J
    
              !***********************************************
            !****************** func11J **********************
            !***********************************************
            double precision function func11J(frail,frail2,frail3,frail4)
    use optim
    
        use tailles
        use comongroup,only:vet,vet2
        use comon,only:alpha,res1,res3,aux1,nig,cdc,sigmae,nmesy,&
            nva2,nva1,npp,nva3,vedc,ve,nea,nb1,betaD,etaD,t1dc,etaydc,etayr,&
            t0,t1,betaR,etaR,typeof,link,vey,c,s_cag_id,s_cag,&
            ut,utt,t0dc,t1dc,method_GH,b_lme,invBi_chol!effet,nva
        use donnees_indiv
            IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2,frail3,frail4
            double precision :: yscalar,eps,finddet,det,res22,alnorm,prod_cag
            integer :: j,i,jj,k,ier,ii
            double precision,dimension(nea*(nea+1)/2)::matv
            !double precision,dimension(2) :: vec,vec2
            !double precision,dimension(2):: p2,p1
            double precision,dimension(nea,1)::  Xea2
            double precision,dimension(nea):: uii, Xea22,Xea
            double precision,dimension(1)::uiiui
            double precision,dimension(nea,nea)::mat
            logical :: upper
            double precision,external::survdcCM,survRCM
            double precision :: resultdc,abserr,resabs,resasc
            double precision :: resultR
            double precision,dimension(1):: current_meanR
            double precision,parameter::pi=3.141592653589793d0
            double precision,dimension(3,3)::matb_chol
    
            upper = .false.
!            if(nmesy(numpat).gt.0) then
!                    allocate(mu1(nmesy(numpat),1))
!            else
!                    allocate(mu1(1,1))
!            end if
    
            i = numpat
            matb_chol = 0.d0
             matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
            matb_chol(3,1) =  invBi_chol(i,4)
            matb_chol(3,2) =  invBi_chol(i,5)
            matb_chol(3,3) =  invBi_chol(i,6)
    
        if(method_GH.eq.1) then
            Xea(1) = frail2
            Xea(2) = frail3
            Xea(3) = frail4
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,Xea(1:nb1))*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            Xea22(4) = frail
            Xea2(4,1) = frail
            else
                    Xea2(1,1) = frail2
            Xea2(2,1) = frail3
            Xea2(3,1) = frail4
            Xea2(4,1) = frail
            Xea22(1) = frail2
            Xea22(2) = frail3
            Xea22(3) = frail4
            Xea22(4) = frail
            end if
    
    
    
            mat = matmul(ut,utt)
    
            jj=0
    do j=1,nea
        do k=j,nea
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
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
    
                    uii = matmul(Xea22,mat)
                    det = finddet(matmul(ut,utt),nea)
    
    
                    uiiui=matmul(uii,Xea2)
    
            if(nmescur.gt.0) then
                    mu1(1:nmescur,1) = mu(1:nmescur,1) +MATMUL(Z1(1:nmescur,1:nb1),Xea22(1:nb1))
            else
                    mu1(1:nmescur,1)  = mu(1:nmescur,1)
            end if
    
                    res1(i) = 0.d0
            res3(i) = 0.d0
            res22 = 0.d0
            !ccccccccccccccccccccccccccccccccccccccccc
            ! pour les recurrences
            !ccccccccccccccccccccccccccccccccccccccccc
            current_meanR = 0.d0
            if(link.eq.1) then
    
            if(typeof.eq.2) then
                    do ii=it_rec,it_rec+nmescurr-1
    
                            if(nva1.gt.0)then
                                    vet = 0.d0
                                    do j=1,nva1
                                            vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(ii,j))
                                    end do
                                    vet = dexp(vet)
                            else
                                    vet=1.d0
                            endif
    
    
                            res1(i) = res1(i)+((t1(ii)/etaR)**betaR)*vet!*dexp(etayr1*Xea22(1))
                            res3(i) = res3(i)+((t0(ii)/etaR)**betaR)*vet!*dexp(etayr1*Xea22(1))
    
    
                    end do
                    else if(typeof.eq.0) then
    
                            do ii=1,nmescurr
    
                            if(nva1.gt.0)then
                                    vet = 0.d0
                                    do j=1,nva1
                                            vet =vet + b1(npp-nva3-nva2-nva1+j)*dble(ve(it_rec+ii-1,j))
                                    end do
                                    vet = dexp(vet)
                            else
                                    vet=1.d0
                            endif
    
                            res1(i) = res1(i)+res1cur(ii)*vet!*dexp(etayr1*Xea22(1))
                            res3(i) = res3(i)+res3cur(ii)*vet!*dexp(etayr1*Xea22(1))
    
                    end do
                    end if
    
            else if(link.eq.2) then
                
                    do ii=it_rec,it_rec+nmescurr-1
            resultR = 0.d0
    
                            call integrationdc(survRCM,t0(ii),t1(ii),resultR,abserr,resabs,resasc,ii,b1,npp,xea22)
    
                    res1(i) = res1(i) + resultR     !c'est déjà res1-res3
    
    
    
                            if((c(ii).eq.1))then
                                    X2cur(1,1) = 1.d0
                                    X2cur(1,2) = t1(ii)
                                    if(nva3.gt.2) then
                                            do k=3,nva3
                                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                                            end do
                                    end if
    
                                    Z1cur(1,1) = 1.d0
                                    Z1cur(1,2) = t1(ii)
                                    current_meanR = current_meanR + MATMUL(X2cur,b1((npp-nva3+1):npp))&
                                            +MATMUL(Z1cur,Xea22(1:2))
                            end if
                    end do
    
    
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
    
            if(link.eq.1) then
    
                    if(typeof.eq.2) then
                            aux1(numpat)=((t1dc(numpat)/etaD)**betaD)*vet2!*dexp(etaydc1*Xea22(1)+etaydc2*Xea22(2))
                    else if(typeof.eq.0) then
                            aux1(numpat)=ut2cur*vet2!*dexp(etaydc1*Xea22(1)+etaydc2*Xea22(2))
                    end if
    
            else !******* Current Level ***************
    
    
                    resultdc = 0.d0
                    call integrationdc(survdcCM,t0dc(numpat),t1dc(numpat),resultdc,abserr,resabs,resasc,numpat,b1,npp,xea22)
    
                    aux1(i) = resultdc
    
                            X2cur(1,1) = 1.d0
                    X2cur(1,2) = t1dc(numpat)
                    if(nva3.gt.2) then
                            do k=3,nva3
                                    X2cur(1,k) = dble(vey(it_cur+1,k))
                            end do
                    end if
    
    
                    Z1cur(1,1) = 1.d0
                    Z1cur(1,2) = t1dc(numpat)
    
                            current_mean = 0.d0
                            current_mean = MATMUL(X2cur,b1((npp-nva3+1):npp))+Matmul(Z1cur,Xea22(1:nb1))
    
            end if
    
    
    
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
    
            !yscalar = norm2(ycurrent(1:nmescur) - mu1(1:nmescur))
    
    
            yscalar = dsqrt(yscalar)
            func11J = 0.d0
            if(link.eq.1) then
            func11J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                        -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -2.d0*dlog(2.d0*pi)&       !
                        -dexp(Xea22(nea)*alpha+dot_product(etaydc,Xea22(1:nb1)))*aux1(numpat)&
                        + cdc(numpat)*dot_product(etaydc,Xea22(1:nb1))&
                        -dexp(Xea22(nea)+dot_product(etayr,Xea22(1:nb1)))*(res1(i)-res3(i)) &
                         +nig(i)*dot_product(etayr,Xea22(1:nb1) )&
                          +Xea22(nea)*(nig(i)+alpha*cdc(i))! &
            else
            func11J =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)&
                                            -uiiui(1)/2.d0-0.5d0*dlog(det)&
                                            -2.d0*dlog(2.d0*pi)&
                                            -dexp(Xea22(nea)*alpha)*aux1(numpat)&
                                            + cdc(numpat)*(etaydc(1)*current_mean(1))&
                                            -dexp(Xea22(nea))*res1(i)+etayr(1)*current_meanR(1)    &
                                            +Xea22(nea)*(nig(i)+alpha*cdc(i))
            end if
 
                            func11J = dexp(func11J)
    
!            deallocate(mu1)
    
        return
    
        end function func11J
    
        !====================================================================
        double precision function risqindivdcCM2(tps,i,bh,np)
    
        use tailles
        use comon
        use betatttps
        use donnees_indiv
        use random_effect
    
        integer::j,i,np,k,n
        double precision::vet2,tps
        double precision,dimension(-2:npmax)::the2
        double precision::bbb,su
        double precision,dimension(np)::bh
    !  double precision::BasisSinhaT1(nbinnerknots+qorder)
        double precision,dimension(nea)::frail
    
    
        k=0
        j=0
        su=0.d0
        bbb=0.d0
    
        frail(1:nea) = re(1:nea)
  
    
        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
            vet2 =vet2 + bh(np-nva3-nva2+j)*dble(vedc(i,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
            X2cur(1,1) = 1.d0
    
        X2cur(1,2) = tps
        if(nva3.gt.0) then
            do k=3,nva3
                        do j=1,nmescur
                    X2cur(1,k) = dble(vey(it_cur+j,k+2))
                    end do
    
                end do
        end if
    
        Z1cur(1,1) = 1.d0
        if(netadc.eq.2) then
        Z1cur(1,2) =tps
        end if
    
    
            current_mean = 1.d0
    
            if(netar.eq.1) then
        current_mean = MATMUL(X2cur,bh((np-nva3+1):np))+Z1cur(1,1)*frail(1)
        else
    !        current_mean = MATMUL(X2cur,bh((np-nva3-2+1):np))+MATMUL(Z2cur,frail(1:netar))
            end if
    
    
    
        select case(typeof)
            case(0) ! calcul du risque splines
    
        if(netar+netadc.ge.1) then
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar)/(effet+1) !nst        !to znaczy ze dzielimy lliczbe wezlow na 2
        else
        n = (np-nva-effet-indic_ALPHA-nb_re - netadc - netar)/(effet+1)
        endif
    
    
            do k=1,n
            !   the1(k-3)=(bh(k))**2.d0
                the2(k-3)=(bh(k))**2.d0
            end do
    
            call susps(tps,the2,nzdc,su,bbb,zi)
    
            if (tps.eq.datedc(ndatedc)) then
                bbb = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
            endif
    
    
            case(2) ! calcul du risque weibull
    
        ! betaD = bh(3)**2
        !  etaD = bh(4)**2
    
            if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf
    
            bbb = (betaD*dexp((betaD-1.d0)*dlog(tps))/(etaD**betaD))
    
        end select
    
        risqindivdcCM2 = bbb*vet2*dexp(etaydc(1)*current_mean(1)+alpha*re(2))!+cdc(i)*etaydc1*current_mean(1)
    
        return
    
        end function risqindivdcCM2
    
    !====================================================================
    !====================================================================
    
    
    
    !====================================================================
        double precision function survRCM(tps,it_surv,bh,np,frail)
    
        use tailles
        use comon
        use betatttps
            use donnees_indiv
            use random_effect
    
        integer::j,np,k,n,it_surv
        double precision::vet2
            double precision::tps
        double precision,dimension(-2:npmax)::the1
        double precision::bbb,su
        double precision,dimension(np)::bh
            double precision,dimension(nea)::frail
    
    
        k=0
        j=0
        su=0.d0
        bbb=0.d0
    !       frail(1:nea) = re(1:nea)
    
    
        if(nva1.gt.0)then
                vet2 = 0.d0
                do j=1,nva1
                    vet2 =vet2 + bh(npp-nva3-nva2-nva1+j)*dble(ve(it_surv,j))
                            end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
    
                            X2cur(1,1) = 1.d0
           X2cur(1,2) = tps
            if(nva3.gt.2) then
                    do k=3,nva3
                                    X2cur(1,k) = dble(vey(it_cur+1,k))
    
                            end do
            end if
    
        Z1cur(1,1) = 1.d0
            if(nb1.eq.2)  Z1cur(1,2) =tps
    
    
            !       write(*,*)'frail',frail
                    current_mean = 1.d0
                    if(nea.gt.1) then
                            current_mean =dot_product(X2cur(1,1:nva3),bh((np-nva3+1):np))&
                                            +dot_product(Z1cur(1,1:nb1),frail(1:nb1))
                    else
                            current_mean = dot_product(X2cur(1,1:nva3),bh((np-nva3+1):np))+Z1cur(1,1:nb1)*frail(1:nb1)
                    end if
    
    
    
        select case(typeof)
            case(0) ! calcul du risque splines
    
        if(netar+netadc.ge.1) then
            n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar)/(effet+1) !nst          !to znaczy ze dzielimy lliczbe wezlow na 2
            else
            n = (np-nva-effet-indic_ALPHA-nb_re - netadc - netar)/(effet+1)
            endif
    
    
            do k=1,n
            !   the1(k-3)=(bh(k))**2.d0
                the1(k-3)=(bh(k))**2.d0
            end do
    
            call susps(tps,the1,nz1,su,bbb,zi)
    
            if (tps.eq.date(ndate)) then
                bbb = 4.d0*the1(n-2-1)/(zi(n-2)-zi(n-2-1))
            endif
    
    
    
            case(2) ! calcul du risque weibull
    
        ! betaD = bh(3)**2
        !  etaD = bh(4)**2
    
                    if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf
    !write(*,*)'betaD',betaD,etaD
            bbb = (betaR*dexp((betaR-1.d0)*dlog(tps))/(etaR**betaR))
    
        end select
    
        survRCM = bbb*vet2*dexp(etayr(1)*current_mean(1))!+re(2))!+cdc(i)*etaydc1*current_mean(1)
    
    !write(*,*)bbb,vet2,dexp(etayr1*current_mean(1))
        return
    
        end function survRCM
    
    
    
    
    
    !====================================================================
         double precision function survdcCM(tps,i,bh,np,frail)
    
        use tailles
        use comon
        use betatttps
        use donnees_indiv
            use residusM,only:Rdc_res
    
        double precision, dimension(1,nva3)::x2curG
        double precision, dimension(1,nb1)::z1curG
    !    double precision, dimension(1,nb1)::z1YcurG
        double precision, dimension(1,nb1)::z1BcurG
        double precision, dimension(1)::current_meanG
        integer::j,i,np,k,n
        double precision::vet2,tps
        double precision,dimension(-2:npmax)::the2
        double precision::bbb,su
        double precision,dimension(np)::bh
    double precision,dimension(nb1)::frail
            double precision,dimension(1) :: Bcv,Bcurrentvalue, cmY,cmGtemp ! add TwoPart
    integer::counter, counter2
            double precision, dimension(1,nvaB)::x2BcurG
    double precision::resultf1, resultf2,f1,f2 !add TwoPart

resultf1=0.d0
resultf2=0.d0
the2=0.d0
        k=0
        j=0
        su=0.d0
        bbb=0.d0
        cmGtemp=0.d0
    counter=0
    counter2=1

        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
            vet2 =vet2 + bh(np-nva3-nva2-nvaB+j)*dble(vedc(i,j))
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif
x2curG=0.d0
    if((nva3-1).gt.0) then ! set the value of covariates at time to event! (interaction must be computed accordingly)
        x2curG(1,1) = 1.d0
        do k=2,nva3
            x2curG(1,k) = dble(vey(it_cur+1,k))
        end do
        
            resultf1=f1(tps) ! maybe remove this calculus if not needed (non-linear slopes)
            resultf2=f2(tps)
if(numInter.ge.1)then
            do counter = 1,numInter !in case of multiple interactions
            if(positionVarT(counter2+3).eq.0) then ! linear
                x2curG(1,positionVarT(counter2+1)) =tps
                if(positionVarT(counter2+2).ne.0) then
                    if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                    x2curG(1,positionVarT(counter2+2)) =tps*dble(vey(it_cur+1,positionVarT(counter2)))
                    else if(positionVarT(counter2).gt.100) then
                    x2curG(1,positionVarT(counter2+2)) =tps*dble(vedc(i,positionVarT(counter2)-100))
                    end if
                end if
            else if(positionVarT(counter2+3).eq.1) then ! f1
                x2curG(1,positionVarT(counter2+1)) =resultf1
                if(positionVarT(counter2+2).ne.0) then
                    if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                    x2curG(1,positionVarT(counter2+2)) =resultf1*dble(vey(it_cur+1,positionVarT(counter2)))
                    else if(positionVarT(counter2).gt.100) then
                    x2curG(1,positionVarT(counter2+2)) =resultf1*dble(vedc(i,positionVarT(counter2)-100))
                    end if   
                end if                    
            else if(positionVarT(counter2+3).eq.2) then ! f2
                x2curG(1,positionVarT(counter2+1)) =resultf2
                if(positionVarT(counter2+2).ne.0) then
                    if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                    x2curG(1,positionVarT(counter2+2)) =resultf2*dble(vey(it_cur+1,positionVarT(counter2)))
                    else if(positionVarT(counter2).gt.100) then
                    x2curG(1,positionVarT(counter2+2)) =resultf2*dble(vedc(i,positionVarT(counter2)-100))
                    end if     
                end if
            end if
                counter2=counter2+4
            end do
        end if
    end if

! open(2,file='C:/Users/dr/Documents/Docs pro/Docs/1_DOC TRAVAIL/2_TPJM/GIT_2019/debug.txt')  
!         write(2,*)'X2BcurG',X2BcurG
!     write(2,*)'nvaB',nvaB
!        write(2,*)'veB',veB
!          write(2,*)'positionVarT',positionVarT
!         write(2,*)'TwoPart',TwoPart
!         write(2,*)'resultf1',resultf1
!          write(2,*)'resultf2',resultf2
!           write(2,*)'counter2',counter2
!          close(2)
    
    if(TwoPart.eq.1) then
    X2BcurG=0.d0
    if((nvaB-1).gt.0) then
    X2BcurG(1,1) = 1.d0
    do k=2,nvaB
    X2BcurG(1,k) = dble(veB(it_curB+1,k))
    end do
    if(numInterB.ge.1) then
        do counter = 1,numInterB ! compute time and interactions at tps
        if(positionVarT(counter2+3).eq.0) then !linear
        X2BcurG(1,positionVarT(counter2+1)) =tps ! time effect
            if(positionVarT(counter2+2).ne.0) then
                if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                X2BcurG(1,positionVarT(counter2+2)) =tps*dble(veB(it_curB+1,positionVarT(counter2)))! interaction
                else if(positionVarT(counter2).gt.100) then
                X2BcurG(1,positionVarT(counter2+2)) =tps*dble(vedc(i,positionVarT(counter2)-100))
                end if
            end if
        else if(positionVarT(counter2+3).eq.1) then !f1
        X2BcurG(1,positionVarT(counter2+1)) =resultf1 ! time effect
            if(positionVarT(counter2+2).ne.0) then
                if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                X2BcurG(1,positionVarT(counter2+2)) =resultf1*dble(veB(it_curB+1,positionVarT(counter2)))! interaction
                else if(positionVarT(counter2).gt.100) then
                X2BcurG(1,positionVarT(counter2+2)) =resultf1*dble(vedc(i,positionVarT(counter2)-100))
                end if
            end if
        else if(positionVarT(counter2+3).eq.2) then !f2
        X2BcurG(1,positionVarT(counter2+1)) =resultf2 ! time effect
            if(positionVarT(counter2+2).ne.0) then
                if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                X2BcurG(1,positionVarT(counter2+2)) =resultf2*dble(veB(it_curB+1,positionVarT(counter2)))! interaction
                else if(positionVarT(counter2).gt.100) then
                X2BcurG(1,positionVarT(counter2+2)) =resultf2*dble(vedc(i,positionVarT(counter2)-100))
                end if
            end if
        end if
        counter2=counter2+4  
        end do
    end if
end if
end if


        z1curG(1,1) = 1.d0
        if(nb1.eq.2)  z1curG(1,2) =tps
     
            current_meanG = 0.d0
if(TwoPart.eq.1) then
    if(nb1.eq.2) then
        z1curG(1,1) = 1.d0 ! random intercept only for now
        z1curG(1,2) = 0.d0!z1Ycur(1,2) = tps
        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
        z1BcurG(1,2) = 1.d0
    else if(nb1.eq.3) then
        z1curG(1,1) = 1.d0 !
        z1curG(1,2) = tps
        z1curG(1,3) = 0.d0
        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
        z1BcurG(1,2) = 0.d0
        z1BcurG(1,3) = 1.d0
    else if(nb1.eq.4) then
        if(nbY.eq.2) then ! intercept + linear slope in each model
        z1curG(1,1) = 1.d0 !
        z1curG(1,2) = tps
        z1curG(1,3) = 0.d0
        z1curG(1,4) = 0.d0
        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
        z1BcurG(1,2) = 0.d0
        z1BcurG(1,3) = 1.d0
        z1BcurG(1,4) = tps
        else ! 2 parametric functions of time + random intercept
            z1curG(1,1) = 1.d0 !
            z1curG(1,2) = resultf1
            z1curG(1,3) = resultf2
            z1curG(1,4) = 0.d0
            z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
            z1BcurG(1,2) = 0.d0
            z1BcurG(1,3) = 0.d0
            z1BcurG(1,4) = 1.d0
        end if
        else if(nb1.eq.5) then
                resultf1=f1(tps) ! need to compute function of time at each point of gauss-kronrod approx.
                resultf2=f2(tps)
        z1curG(1,1) = 1.d0 !intercept continuous WATCHOUT ORDER OF RE!!
        z1curG(1,2) = resultf1
        z1curG(1,3) = resultf2
        z1curG(1,4) = 0.d0 
        z1curG(1,5) = 0.d0
        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
        z1BcurG(1,2) = 0.d0
        z1BcurG(1,3) = 0.d0
        z1BcurG(1,4) = 1.d0
        z1BcurG(1,5) = tps

    end if
                            
    Bcurrentvalue=0.d0
    Bcv=0.d0
    Bcv=dot_product(X2BcurG(1,1:nvaB),bh((np-nvaB+1):np))+dot_product(z1BcurG(1,1:nb1),frail(1:nb1))
    Bcurrentvalue=dexp(Bcv)/(1+dexp(Bcv))
                        
    cmY = (dot_product(x2curG(1,1:nva3),bh((np-nva3-nvaB+1):(np-nvaB)))+dot_product(z1curG(1, 1:nb1),frail(1:nb1)))

    if(GLMloglink0.eq.0) then
        current_meanG = cmY*Bcurrentvalue
    else if(GLMloglink0.eq.1) then
        if(MTP0.eq.0) then
            current_meanG = dexp(cmY)*Bcurrentvalue
        else if(MTP0.eq.1) then
            current_meanG = dexp(cmY)
        end if
    end if

            
else if(TwoPart.eq.0) then
    if(nb1.eq.2) then
        z1curG(1,1) = 1.d0 ! random intercept only for now
        z1curG(1,2) = tps!z1Ycur(1,2) = t1dc(i)
    else if(nb1.eq.3) then
        resultf1=f1(tps) ! need to compute function of time at each point of gauss-kronrod approx.
        resultf2=f2(tps)
        z1curG(1,1) = 1.d0 !
        z1curG(1,2) = resultf1
        z1curG(1,3) = resultf2
    end if
              
     if(GLMloglink0.eq.0) then
                if(nea.gt.1) then
                    current_meanG =dot_product(x2curG(1,1:nva3),bh((np-nva3+1):np))&
 +dot_product(z1curG(1,1:nb1),frail(1:nb1))
                else
                    current_meanG = dot_product(x2curG(1,1:nva3),bh((np-nva3+1):np))+z1curG(1,1:nb1)*frail(1:nb1)
                end if
    else  if(GLMloglink0.eq.1) then
                if(nea.gt.1) then
                    cmY =dot_product(x2curG(1,1:nva3),bh((np-nva3+1):np))&
                    +dot_product(z1curG(1,1:nb1),frail(1:nb1))
                    current_meanG=dexp(cmY)
                else
                    cmY = dot_product(x2curG(1,1:nva3),bh((np-nva3+1):np))+z1curG(1,1:nb1)*frail(1:nb1)
                    current_meanG=dexp(cmY)
                end if
    end if
end if
                    

        select case(typeof)
            case(0) ! calcul du risque splines
    
        if(netar+netadc.ge.1) then
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar)/(effet+1) !nst        !to znaczy ze dzielimy lliczbe wezlow na 2
        else
        n = (np-nva-effet-indic_ALPHA-nb_re - netadc - netar)/(effet+1)
        endif
    
    
            do k=1,n
                if(typeJoint.eq.2)the2(k-3)=(bh(k))**2.d0
                if(typeJoint.eq.3)the2(k-3)=(bh(k+n))**2.d0
            end do
    
            call susps(tps,the2,nzdc,su,bbb,zi)
    
            if (tps.eq.datedc(ndatedc)) then
                bbb = 4.d0*the2(n-2-1)/(zi(n-2)-zi(n-2-1))
            endif
    
            case(2) ! calcul du risque weibull
    
       
            if (tps.eq.0.d0) tps = 1d-12 ! utile car log(0) => -Inf
    
            bbb =     (betaD*dexp((betaD-1.d0)*dlog(tps))/(etaD**betaD)) !((tps/etaD)**betaD)!
    
        end select

            if(res_ind.eq.1) bbb = Rdc_res(i)
    
                                                                                              
    
    if(link.eq.2) then
            survdcCM =bbb*vet2*dexp(etaydc(1)*current_meanG(1))!+cdc(i)*etaydc1*current_meanG(1)  
else if(link.eq.3) then
            survdcCM =bbb*vet2*dexp(etaydc(1)*Bcurrentvalue(1)+etaydc(2)*cmY(1))!+cdc(i)*etaydc1*current_meanG(1) 
            end if
            
        return
    
        end function survdcCM
    
    !====================================================================
    !====================================================================
    
    
    
    
    
    
    
    
    
    
    
    
            double precision function funcG(frail5,frail4,frail3,frail2,frail)
    
        use tailles
        !use comongroup,only:vet2!vet
        use optim
        use comon,only:cdc,sigmae,nmesy,&
            nva2,npp,nva3,vedc,nb1,betaD,etaD,t0dc,t1dc,etaydc,link,&
            vey, typeof,s_cag_id,s_cag,ut,utt,method_GH,b_lme,invBi_chol,&
            nbB, nby, nvaB, nmesB,TwoPart,veB,numInter,numInterB,positionVarT
            !auxig,alpha,sig2,res1,res3,nb1,nea,nig,utt,
        use donnees_indiv
        IMPLICIT NONE
    
        double precision,intent(in)::frail,frail2,frail3,frail4,frail5
        double precision, dimension(numpat)::auxG
        double precision, dimension(1,nva3)::x2curG
        double precision, dimension(1,nvaB)::x2BcurG
    !    double precision, dimension(1,nb1)::z1curG
        double precision, dimension(1,nb1)::z1YcurG
        double precision, dimension(1,nb1)::z1BcurG
double precision, dimension(nmesy(numpat),1):: mu1G
double precision, dimension(nmesB(numpat),1):: mu1BG
        double precision,dimension(nb1*(nb1+1)/2)::matv
        double precision,dimension(nb1,1)::  Xea2
        double precision,dimension(nb1):: uii, Xea22,XeaG,Xea
        double precision,dimension(1)::uiiui
        double precision,dimension(nb1,nb1)::mat,matb_chol
        double precision, dimension(1)::current_meanG
        double precision :: yscalar,yscalarlog,eps,finddet,det,alnorm,prod_cag,yscal1,yscal2
        integer :: j,jj,i,k,ier
        logical :: upper
        double precision,external::survdcCM
        double precision :: resultdc,abserr,resabs,resasc,vet2
        double precision,parameter::pi=3.141592653589793d0
        double precision :: resultf1, resultf2, f1, f2 ! add TwoPart
        double precision,dimension(1) :: Bscalar,Bscal,Bcv,Bcurrentvalue, cmY,cmGtemp
        integer :: counter, counter2 ! add for current-level interaction
        upper = .false.
        i = numpat

 uiiui=0.d0
 funcG=0.d0
 current_meanG = 0.d0
 cmGtemp=0.d0
 det=0.d0
Xea(1)=0.d0
Xea2=0.d0
Xea22=0.d0
XeaG=0.d0
uii=0.d0
mu1G=0.d0
mu1BG=0.d0
auxG=0.d0
resultdc=0.d0

resultf1=0.d0
resultf2=0.d0
!open(2,file='C:/Users/dr/Documents/Docs pro/Docs/1_DOC TRAVAIL/2_TPJM/GIT_2019/debug.txt')  
!       write(2,*)'ping'
!    close(2)
!    stop

    if(nb1.eq.1) then
            if(method_GH.eq.1) then
            Xea(1) = b_lme(i,1) +invBi_chol(i,1)*frail*sqrt(2.d0)
            else
            Xea(1) = frail!*sqrt(2.d0)*ut(1,1)
            end if
else if(nb1.gt.1) then
if(nb1.eq.2) then
        matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
        if(method_GH.eq.1) then
            XeaG(1) = frail
            XeaG(2) = frail2
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,XeaG)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
            else
            Xea2(1,1) = frail!
        Xea2(2,1) = frail2!
        Xea22(1) = frail!
        Xea22(2) = frail2!
            end if
else if(nb1.eq.3) then

        matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
            matb_chol(3,1) =  invBi_chol(i,4)
            matb_chol(3,2) =  invBi_chol(i,5)
            matb_chol(3,3) =  invBi_chol(i,6)
        if(method_GH.eq.1) then
            XeaG(1) = frail
            XeaG(2) = frail2
            XeaG(3) = frail3
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,XeaG)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            else
            Xea2(1,1) = frail!
        Xea2(2,1) = frail2!
        Xea2(3,1) = frail3!
        Xea22(1) = frail!
        Xea22(2) = frail2
        Xea22(3) = frail3
            end if
else if(nb1.eq.4) then

        matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
            matb_chol(3,1) =  invBi_chol(i,4)
            matb_chol(3,2) =  invBi_chol(i,5)
            matb_chol(3,3) =  invBi_chol(i,6)
            matb_chol(4,1) =  invBi_chol(i,7)
            matb_chol(4,2) =  invBi_chol(i,8)
            matb_chol(4,3) =  invBi_chol(i,9)
            matb_chol(4,4) =  invBi_chol(i,10)
        if(method_GH.eq.1) then
            XeaG(1) = frail
            XeaG(2) = frail2
            XeaG(3) = frail3
            XeaG(4) = frail4
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,XeaG)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            else
            Xea2(1,1) = frail!
        Xea2(2,1) = frail2!
        Xea2(3,1) = frail3!
        Xea2(4,1) = frail4!
        Xea22(1) = frail!
        Xea22(2) = frail2
        Xea22(3) = frail3
        Xea22(4) = frail4
            end if
else if(nb1.eq.5) then

        matb_chol = 0.d0
            matb_chol(1,1) = invBi_chol(i,1)
            matb_chol(2,1) =  invBi_chol(i,2)
            matb_chol(2,2) =  invBi_chol(i,3)
            matb_chol(3,1) =  invBi_chol(i,4)
            matb_chol(3,2) =  invBi_chol(i,5)
            matb_chol(3,3) =  invBi_chol(i,6)
            matb_chol(4,1) =  invBi_chol(i,7)
            matb_chol(4,2) =  invBi_chol(i,8)
            matb_chol(4,3) =  invBi_chol(i,9)
            matb_chol(4,4) =  invBi_chol(i,10)
            matb_chol(5,1) =  invBi_chol(i,11)
            matb_chol(5,2) =  invBi_chol(i,12)
            matb_chol(5,3) =  invBi_chol(i,13)
            matb_chol(5,4) =  invBi_chol(i,14)
            matb_chol(5,5) =  invBi_chol(i,15)

        if(method_GH.eq.1) then
            XeaG(1) = frail
            XeaG(2) = frail2
            XeaG(3) = frail3
            XeaG(4) = frail4
            XeaG(5) = frail5
            Xea22(1:nb1) = b_lme(i,1:nb1) +  Matmul(matb_chol,XeaG)*sqrt(2.d0)
            Xea2(1:nb1,1) = Xea22(1:nb1)
    
            else
            Xea2(1,1) = frail!
        Xea2(2,1) = frail2!
        Xea2(3,1) = frail3!
        Xea2(4,1) = frail4!
        Xea2(5,1) = frail5!
        Xea22(1) = frail!
        Xea22(2) = frail2
        Xea22(3) = frail3
        Xea22(4) = frail4
        Xea22(5) = frail5
            end if
            end if            

        mat=0.d0
        matv=0.d0
        mat = matmul(ut,utt)

        jj=0
    ! jjj = 0
        do j=1,nb1
        do k=j,nb1
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
        end do
        end do
        ier = 0
        eps = 1.d-10
    

            call dsinvj(matv,nb1,eps,ier)
    mat=0.d0
        do j=1,nb1
                do k=1,nb1
            if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                end do
    
    
        uii = matmul(Xea22,mat)
            det = finddet(matmul(ut,utt),nb1)
    
            uiiui=matmul(uii,Xea2)
 end if


        if(nmescur.gt.0) then
        if(nb1.eq.1) then
            mu1G(1:nmescur,1) = mu(1:nmescur,1) +Xea(1)*Z1(1:nmescur,1)
            else if(nb1.gt.1) then
            mu1G(1:nmescur,1) = mu(1:nmescur,1) +MATMUL(Z1(1:nmescur,1:nby),Xea2(1:nby,1))
            end if
        else
            mu1G(1:nmescur,1)  = mu(1:nmescur,1)
        end if

        
                ! add TwoPart
        if(TwoPart.eq.1) then
            if(nmescurB.gt.0) then 
                mu1BG(1:nmescurB,1) = muB(1:nmescurB,1) +MATMUL(Z1B(1:nmescurB,1:nbB),Xea2(nby+1:nb1,1))
            else
                mu1BG(1:nmescurB,1)  = muB(1:nmescurB,1)
            end if
        end if


    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc

        if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
    
            vet2 =vet2 + b1(npp-nva3-nva2-nvaB+j)*dble(vedc(i,j))
    
                end do
                vet2 = dexp(vet2)
            else
                vet2=1.d0
            endif

            if(link.eq.1) then !******** Random Effects **********************
    ! pour le calcul des integrales / pour la survie, pas les donno?=es recurrentes:
        
            if(typeof.eq.2) then
                auxG(i)=((t1dc(i)/etaD)**betaD)*vet2!*dexp(etaydc1*frail)
            else if(typeof.eq.0) then
                auxG(i)=ut2cur*vet2!*dexp(etaydc1*frail)
            end if
        !end if
        else if(link.ge.2) then !********** Current Mean ****************
        counter=0
        counter2=1
      if(nb1.eq.1) then
        call integrationdc(survdcCM,t0dc(i),t1dc(i),resultdc,abserr,resabs,resasc,i,b1,npp,Xea)
    else if(nb1.gt.1) then
        call integrationdc(survdcCM,t0dc(i),t1dc(i),resultdc,abserr,resabs,resasc,i,b1,npp,xea22)
    end if
        auxG(i) = resultdc
        
x2curG=0.d0
    if((nva3-1).gt.0) then ! set the value of covariates at time to event! (interaction must be computed accordingly)
        x2curG(1,1) = 1.d0
        do k=2,nva3
            x2curG(1,k) = dble(vey(it_cur+1,k))
        end do
        
            resultf1=f1(t1dc(i)) ! maybe remove this calculus if not needed (non-linear slopes)
            resultf2=f2(t1dc(i))
if(numInter.ge.1)then

            do counter = 1,numInter !in case of multiple interactions
            if(positionVarT(counter2+3).eq.0) then ! linear
                x2curG(1,positionVarT(counter2+1)) =t1dc(i)
                if(positionVarT(counter2+2).ne.0) then
                    if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                    x2curG(1,positionVarT(counter2+2)) =t1dc(i)*dble(vey(it_cur+1,positionVarT(counter2)))
                    else if(positionVarT(counter2).gt.100) then
                    x2curG(1,positionVarT(counter2+2)) =t1dc(i)*dble(vedc(i,positionVarT(counter2)-100))
                    end if
                end if
            else if(positionVarT(counter2+3).eq.1) then ! f1
                x2curG(1,positionVarT(counter2+1)) =resultf1
                if(positionVarT(counter2+2).ne.0) then
                    if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                    x2curG(1,positionVarT(counter2+2)) =resultf1*dble(vey(it_cur+1,positionVarT(counter2)))
                    else if(positionVarT(counter2).gt.100) then
                    x2curG(1,positionVarT(counter2+2)) =resultf1*dble(vedc(i,positionVarT(counter2)-100))
                    end if
                end if
            else if(positionVarT(counter2+3).eq.2) then ! f2
                x2curG(1,positionVarT(counter2+1)) =resultf2
                if(positionVarT(counter2+2).ne.0) then
                    if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                    x2curG(1,positionVarT(counter2+2)) =resultf2*dble(vey(it_cur+1,positionVarT(counter2)))
                    else if(positionVarT(counter2).gt.100) then
                    x2curG(1,positionVarT(counter2+2)) =resultf2*dble(vedc(i,positionVarT(counter2)-100))
                    end if
                end if
            end if
                counter2=counter2+4
            end do
        end if
    end if

    
    
    if(TwoPart.eq.1) then
    X2BcurG=0.d0
    if((nvaB-1).gt.0) then
    X2BcurG(1,1) = 1.d0
    do k=2,nvaB
    X2BcurG(1,k) = dble(veB(it_curB+1,k))
    end do
    if(numInterB.ge.1) then
    do counter = 1,numInterB ! compute time and interactions at t1dc(i)
        if(positionVarT(counter2+3).eq.0) then !linear
            X2BcurG(1,positionVarT(counter2+1)) =t1dc(i) ! time effect
            if(positionVarT(counter2+2).ne.0) then
                if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                X2BcurG(1,positionVarT(counter2+2)) =t1dc(i)*dble(veB(it_curB+1,positionVarT(counter2)))! interaction
                else if(positionVarT(counter2).gt.100) then
                X2BcurG(1,positionVarT(counter2+2)) =t1dc(i)*dble(vedc(i,positionVarT(counter2)-100))
                end if
            end if
        else if(positionVarT(counter2+3).eq.1) then !f1
            X2BcurG(1,positionVarT(counter2+1)) =resultf1 ! time effect
            if(positionVarT(counter2+2).ne.0) then
                if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                X2BcurG(1,positionVarT(counter2+2)) =resultf1*dble(veB(it_curB+1,positionVarT(counter2)))! interaction
                else if(positionVarT(counter2).gt.100) then
                X2BcurG(1,positionVarT(counter2+2)) =resultf1*dble(vedc(i,positionVarT(counter2)-100))
                end if
            end if
        else if(positionVarT(counter2+3).eq.2) then !f2
            X2BcurG(1,positionVarT(counter2+1)) =resultf2 ! time effect
            if(positionVarT(counter2+2).ne.0) then
                if(positionVarT(counter2).le.100) then ! if interaction terms not included 
                X2BcurG(1,positionVarT(counter2+2)) =resultf2*dble(veB(it_curB+1,positionVarT(counter2)))! interaction
                else if(positionVarT(counter2).gt.100) then
                X2BcurG(1,positionVarT(counter2+2)) =resultf2*dble(vedc(i,positionVarT(counter2)-100))
                end if
            end if
        end if
        counter2=counter2+4  
    end do
    end if
end if
end if    

z1YcurG=0.d0
z1BcurG=0.d0

             ! no ping here             
             
            z1YcurG(1,1) = 1.d0
            if(nb1.eq.2) then
                z1YcurG(1,2) = t1dc(i)
            end if  
            

    if(nb1.eq.1) then
            current_meanG(1) =dot_product(x2curG(1,1:nva3),b1((npp-nva3+1):npp))+z1YcurG(1,1)*Xea(1)
    else if(nb1.gt.1) then
    
    if(TwoPart.eq.1) then
if(nb1.eq.2) then
    z1YcurG(1,1) = 1.d0 ! random intercept only for now
    z1YcurG(1,2) = 0.d0!z1Ycur(1,2) = t1dc(i)
    z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
    z1BcurG(1,2) = 1.d0
else if(nb1.eq.3) then
    z1YcurG(1,1) = 1.d0 !
    z1YcurG(1,2) = t1dc(i)
    z1YcurG(1,3) = 0.d0
    z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
    z1BcurG(1,2) = 0.d0
    z1BcurG(1,3) = 1.d0
else if(nb1.eq.4) then
    if(nbY.eq.2) then ! random intercept and slope in each model
        z1YcurG(1,1) = 1.d0 !
        z1YcurG(1,2) = t1dc(i)
        z1YcurG(1,3) = 0.d0
        z1YcurG(1,4) = 0.d0
        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
        z1BcurG(1,2) = 0.d0
        z1BcurG(1,3) = 1.d0
        z1BcurG(1,4) = t1dc(i)
    else
        resultf1=f1(t1dc(i)) ! need to compute function of time at each point of gauss-kronrod approx.
        resultf2=f2(t1dc(i))
        
        z1YcurG(1,1) = 1.d0 !
        z1YcurG(1,2) = resultf1
        z1YcurG(1,3) = resultf2
        z1YcurG(1,4) = 0.d0
        z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
        z1BcurG(1,2) = 0.d0
        z1BcurG(1,3) = 0.d0
        z1BcurG(1,4) = 1.d0
    end if
else if(nb1.eq.5) then
    resultf1=f1(t1dc(i)) ! need to compute function of time at each point of gauss-kronrod approx.
    resultf2=f2(t1dc(i))
    z1YcurG(1,1) = 1.d0 !intercept continuous WATCHOUT ORDER OF RE!!
    z1YcurG(1,2) = resultf1
    z1YcurG(1,3) = resultf2
    z1YcurG(1,4) = 0.d0 
    z1YcurG(1,5) = 0.d0
    z1BcurG(1,1) = 0.d0 ! need to decide intercept / time here !
    z1BcurG(1,2) = 0.d0
    z1BcurG(1,3) = 0.d0
    z1BcurG(1,4) = 1.d0
    z1BcurG(1,5) = t1dc(i)

end if


                        Bcurrentvalue=0.d0
                        Bcv=0.d0

                        Bcv=MATMUL(X2BcurG,b1((npp-nvaB+1):npp))+Matmul(z1BcurG,Xea22)
                        Bcurrentvalue=dexp(Bcv)/(1+dexp(Bcv))
                    
cmY=0.d0
        cmY = (MATMUL(x2curG,b1((npp-nva3-nvaB+1):(npp-nvaB)))+Matmul(z1YcurG,Xea22))
!open(2,file='C:/Users/dr/Documents/Docs pro/Docs/1_DOC TRAVAIL/2_TPJM/GIT_2019/debug.txt')  
!       write(2,*)'boxcoxlambda',boxcoxlambda
!     close(2)
     
    if(GLMloglink0.eq.0) then
            current_meanG = cmY*Bcurrentvalue
    else if(GLMloglink0.eq.1) then
        if(MTP0.eq.0) then
            current_meanG = dexp(cmY)*Bcurrentvalue
        else if(MTP0.eq.1) then
            current_meanG = dexp(cmY) 
        end if
    end if
 else if(TwoPart.eq.0) then
 
 if(nb1.eq.2) then
    z1YcurG(1,1) = 1.d0 ! random intercept only for now
    z1YcurG(1,2) = t1dc(i)!z1Ycur(1,2) = t1dc(i)
else if(nb1.eq.3) then
    resultf1=f1(t1dc(i)) ! need to compute function of time at each point of gauss-kronrod approx.
    resultf2=f2(t1dc(i))
    z1YcurG(1,1) = 1.d0 !
    z1YcurG(1,2) = resultf1
    z1YcurG(1,3) = resultf2

end if
 
 
 if(GLMloglink0.eq.0) then
current_meanG = MATMUL(X2curG,b1((npp-nva3+1):npp))+Matmul(z1YcurG,Xea22)
else if(GLMloglink0.eq.1) then
cmY = MATMUL(X2curG,b1((npp-nva3+1):npp))+Matmul(z1YcurG,Xea22)
current_meanG = dexp(cmY)
end if
                    end if    
            end if
    
        end if
    


        yscalar = 0.d0
        yscal1 = 0.d0
        yscal2 = 0.d0
        yscalarlog = 0.d0
        !********* Left-censoring ***********
            prod_cag = 1.d0
    
        if(s_cag_id.eq.1)then
    
            do k = 1,nmescur
    
                if(ycurrent(k).le.s_cag) then
                    prod_cag = prod_cag*(1.d0-alnorm((mu1G(k,1)-s_cag)/sqrt(sigmae),upper))
                    !(0.5d0*(1.d0-erf((mu1G(k)-s_cag)/(sigmae*dsqrt(2.d0)))))
                else
                    if(GLMloglink0.eq.0) then
                        yscalar = yscalar + (ycurrent(k)-mu1G(k,1))**2
                    else if(GLMloglink0.eq.1) then
                        if(TwoPart.eq.0) then
                            yscalar = yscalar + (dlog(ycurrent(k))-(mu1G(k,1)-(sigmae/2)))**2
                            yscalarlog = yscalarlog - dlog(ycurrent(k))
                        else if(TwoPart.eq.1) then
                            if(ycurrent(k).ne.0) then
                                if(MTP0.eq.0) then
                                    yscalar = yscalar + (dlog(ycurrent(k))-(mu1G(k,1)-(sigmae/2)))**2
                                    yscalarlog = yscalarlog - dlog(ycurrent(k))
                                else if (MTP0.eq.1) then

                                end if
                            end if
                        end if
                    end if
                end if
            end do
        else ! no left censoring
            do k=1,nmescur
                if(GLMloglink0.eq.0) then ! original scale (no transformation of the longi outcome)
                        yscalar = yscalar + (ycurrent(k)-mu1G(k,1))**2
                else if(GLMloglink0.eq.1) then ! lognormal density
                    if(TwoPart.eq.0) then
                        yscalar = yscalar + (dlog(ycurrent(k))-(mu1G(k,1)-(sigmae/2)))**2
                        yscalarlog = yscalarlog - dlog(ycurrent(k)) 
                    else if(TwoPart.eq.1) then ! two-part model for the longitudinal outcome
                        if(ycurrent(k).ne.0) then
                            if(MTP0.eq.0) then
                                yscalar = yscalar + (dlog(ycurrent(k))-(mu1G(k,1)-(sigmae/2)))**2
                                yscalarlog = yscalarlog - dlog(ycurrent(k))
                            else if (MTP0.eq.1) then ! marginal two-part model
yscalar = yscalar + (dlog(ycurrent(k))-(mu1G(k,1)-dlog(Bcurrentvalue(1))-(sigmae/2)))**2
    yscalarlog = yscalarlog - dlog(ycurrent(k))
                            end if
                        end if
                    end if
               end if
            end do
        end if

        yscalar = dsqrt(yscalar)    
        if(prod_cag.lt.0.1d-321)prod_cag= 0.1d-321

Bscal(1)=0.d0
    Bscalar(1)=0.d0
    if(TwoPart.eq.1) then ! binary part contribution
        do k=1,nmescurB
            if(MTP0.eq.0) then
            Bscal(1)=(Bcurrent(k)*mu1BG(k,1))+dlog(1.d0-(dexp(mu1BG(k,1))/(1+dexp(mu1BG(k,1)))))
            Bscalar(1) = Bscalar(1) + Bscal(1)       
            else if(MTP0.eq.1) then
Bscalar(1) = Bscalar(1) - dlog(1.d0+dexp(mu1BG(k,1)))
            end if
        end do
    end if
     
funcG=0.d0
       if (method_GH.ne.3) then ! loglikelihood
    if(nb1.eq.1) then
            if(link.eq.1) then
        funcG =   dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog& !longi part
                    - (Xea(1)**2.d0)/(2.d0*ut(1,1)**2) & ! random effects density
                    - dlog(ut(1,1))-dlog(2.d0*pi)/2.d0& ! random effects density
                        -auxG(i)*dexp(etaydc(1)*Xea(1))  + cdc(i)*etaydc(1)*Xea(1) ! survival part
        else
        funcG =   dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        - (Xea(1)**2.d0)/(2.d0*ut(1,1)**2)&
                    - dlog(ut(1,1))-dlog(2.d0*pi)/2.d0&
                        -auxG(i) &
                        + cdc(i)*etaydc(1)*current_meanG(1)
    
        
        end if
        else if(nb1.eq.2) then
        if(link.eq.1) then
        funcG = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                    -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dlog(2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)*dexp(etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))&
                        + cdc(i)*(etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))
    
        else if(link.eq.2) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dlog(2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*current_meanG(1))
        else if(link.eq.3) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dlog(2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*Bcurrentvalue(1)+etaydc(2)*cmY(1))
        end if
        else 
        if(link.eq.1) then
        funcG = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                  -uiiui(1)/2.d0-0.5d0*dlog(det)&
                        -dble(nb1)/2.d0*dlog(2.d0*pi)&   !-(nb1/2.d0)*dlog(det*2.d0*pi)&
                     +Bscalar(1)& ! add TwoPart
                   -auxG(i)*dexp(dot_product(etaydc,Xea22(1:nb1)))&
                    + cdc(i)*dot_product(etaydc,Xea22(1:nb1))
        else if(link.eq.2) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        -uiiui(1)/2.d0-(dble(nb1)/2.d0)*dlog(det*2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*current_meanG(1))
        else if(link.eq.3) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        -uiiui(1)/2.d0-(dble(nb1)/2.d0)*dlog(det*2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*Bcurrentvalue(1)+etaydc(2)*cmY(1))
        end if
        end if
        else ! Monte-carlo
            if(nb1.eq.1) then
            if(link.eq.1) then
        funcG =   dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog& !longi part
                   ! - (Xea(1)**2.d0)/(2.d0*ut(1,1)**2) & ! random effects density
                   ! - dlog(ut(1,1))-dlog(2.d0*pi)/2.d0& ! random effects density
                        -auxG(i)*dexp(etaydc(1)*Xea(1))  + cdc(i)*etaydc(1)*Xea(1) ! survival part
        else
        funcG =   dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                    !    - (Xea(1)**2.d0)/(2.d0*ut(1,1)**2)&
                    !- dlog(ut(1,1))-dlog(2.d0*pi)/2.d0&
                        -auxG(i) &
                        + cdc(i)*etaydc(1)*current_meanG(1)
        end if
        else if(nb1.eq.2) then
        if(link.eq.1) then
        funcG = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                    !-uiiui(1)/2.d0-0.5d0*dlog(det)&
                    !    -dlog(2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)*dexp(etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))&
                        + cdc(i)*(etaydc(1)*Xea22(1)+etaydc(2)*Xea22(2))
    
        else if(link.eq.2) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        !-uiiui(1)/2.d0-0.5d0*dlog(det)&
                        !-dlog(2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*current_meanG(1))
        else if(link.eq.3) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        !-uiiui(1)/2.d0-0.5d0*dlog(det)&
                        !-dlog(2.d0*pi)&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*Bcurrentvalue(1)+etaydc(2)*cmY(1))
        end if
        else if(nb1.gt.2) then
        if(link.eq.1) then
        funcG = dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                     +Bscalar(1)& ! add TwoPart
                   -auxG(i)*dexp(dot_product(etaydc,Xea22(1:nb1)))&
                    + cdc(i)*dot_product(etaydc,Xea22(1:nb1))
        else if(link.eq.2) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*current_meanG(1))
        else if(link.eq.3) then
        funcG =  dlog(prod_cag)-(yscalar**2.d0)/(2.d0*sigmae)+yscalarlog&
                        +Bscalar(1)& ! add TwoPart
                        -auxG(i)&
                        + cdc(i)*(etaydc(1)*Bcurrentvalue(1)+etaydc(2)*cmY(1))
        end if
        end if
        end if
        
    funcG = dexp(funcG)

        return
    
        end function funcG
    
    
    
        !====================================================================
    !====================================================================
    
    
    
    
    
    
    
    
    
    !! Monte-carlo
    subroutine MC_JointModels(ss,func2,ndim,intpoints)
    use var_surrogate, only: nbre_sim
    use donnees ! pour les points et poids de quadrature (fichier Adonnees.f90)
    use comon, only:nb1,nodes_number
    use donnees_indiv
    !use mpi
    !$ use OMP_LIB
    
    implicit none
    integer ::ii,nsimu !code,erreur,nbrejet,stemp,tid1,npg,kk,j,k,ier !maxmes= nombre de dimension ou encore dimension de X
    integer,intent(in):: ndim
    double precision,intent(out)::ss
    double precision,dimension(nodes_number,nb1),intent(in)::intpoints
    double precision::auxfunca !,mu1,vc1,ss2
    double precision,dimension(nodes_number)::ss1
    double precision::x2222,somp !eps !ymarg contient le resultat de l'integrale
    double precision::func2
        
    nsimu=nbre_sim
    x2222=0.d0
    somp=0.d0

  ! --------------------- boucle du MC ------------------------
    auxfunca=0.d0
    ss=0.d0
    ss1=0.d0
   select case(ndim)
     case(1)
        ii=0
         !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,intpoints,ss1)&
         !$OMP SCHEDULE(Static)
            do ii=1,nsimu
                auxfunca=func2(0.d0,0.d0,0.d0,0.d0,intpoints(ii,1))
                ss1(ii)=auxfunca
            end do
         !$OMP END PARALLEL DO
   case(2)

        ii=0
         !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,intpoints,ss1)&
         !$OMP SCHEDULE(Static)
            do ii=1,nsimu
                auxfunca=func2(0.d0,0.d0,0.d0,intpoints(ii,2),intpoints(ii,1))
                ss1(ii)=auxfunca
            end do
         !$OMP END PARALLEL DO
         
       case(3)

        ii=0
         !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,intpoints)&
         !$OMP REDUCTION(+:ss1) SCHEDULE(Dynamic,1)
            do ii=1,nsimu
                auxfunca=func2(0.d0,0.d0,intpoints(ii,3),intpoints(ii,2),intpoints(ii,1))
                ss1(ii)=auxfunca
            end do
         !$OMP END PARALLEL DO
         
               case(4)

        ii=0
         !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,intpoints,ss1)&
         !$OMP SCHEDULE(Static)
            do ii=1,nsimu
                auxfunca=func2(0.d0,intpoints(ii,4),intpoints(ii,3),intpoints(ii,2),intpoints(ii,1))
                ss1(ii)=auxfunca
            end do
         !$OMP END PARALLEL DO
          
                 case(5)

        ii=0
         !$OMP PARALLEL DO default(none) PRIVATE (ii,auxfunca) SHARED(nsimu,intpoints,ss1)&
         !$OMP SCHEDULE(Static)
            do ii=1,nsimu
                auxfunca=func2(intpoints(ii,5),intpoints(ii,4),intpoints(ii,3),intpoints(ii,2),intpoints(ii,1))
                ss1(ii)=auxfunca
            end do
         !$OMP END PARALLEL DO
                  
    end select

    ss=dble(SUM(ss1))/dble(nsimu) 

    return 
  end subroutine MC_JointModels
  
  

  ! convert these fonctions to make them arguments in the R function so the user can define specific functions.
  ! conversion can be achieved by using a C-wrapper to call the R function directly from FORTRAN 
  !f1
    function f1(t) result(res)
       double precision, intent(in) :: t !input
       double precision :: res! output
       res = dexp(-3.5d0*t)
    end function f1

!f2
     function f2(t) result(res)
       double precision, intent(in) :: t !input
       double precision :: res! output
       res = t
    end function f2

  
  