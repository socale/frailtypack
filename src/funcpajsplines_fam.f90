


!========================          FUNCPAJ_SPLINES  FOR JOINT NESTED MODEL        ====================
    double precision function funcpajsplines_fam(b,np,id,thi,jd,thj,k0)

    use tailles
    !use comon,only:AG,auxig,auxif,fam,fsize,t0,t1,t0dc,t1dc,nt0dc,res4
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m,mm3,mm2,mm1,mm,&
    im3,im2,im1,im,mm3dc,mm2dc,mm1dc,mmdc,im3dc,im2dc,im1dc,imdc,date,datedc,zi,&
    c,cdc,nt0,nt1,nt1dc,nsujet,nva,nva1,nva2,ndate,ndatedc,nst,& 
    effet,stra,ve,vedc,pe,ng,g,nig,indic_ALPHA,ALPHA,theta,nstRec,k0T, & 
    fam,nfam,fsize,indic_xi, xi, eta, & !for family 
    aux1,aux2,res1,res3,kkapa,resnonpen,nb_gl 
    use residusM
    use comongroup,only:vet,vet2,the2!,the1

    IMPLICIT NONE

! *** NOUVELLLE DECLARATION F90 :

    integer,intent(in)::id,jd,np
    double precision,dimension(np),intent(in)::b
    double precision,dimension(2),intent(in)::k0
    double precision,intent(in)::thi,thj
    double precision,dimension(-2:npmax,nstRec)::the1T
    integer::n,i,j,k,vj,ig,choix,jj
    integer,dimension(ngmax)::cpt
    double precision::pe2,sum,som2,res,h1
!deleted unused inv argument %myriam
!added inv2 for 1/eta
!deleted unused inv2 argument %myriam
    double precision,dimension(nstRec)::pe1T,som1T
    double precision,dimension(np)::bh
    double precision,dimension(ngmax)::res2,res1dc,res2dc &
    ,res3dc,integrale1,integrale2,integrale3
    double precision,dimension(nfam)::integrale3f
!AD: for death,change dimension 
    double precision,dimension(ndatemax,nstRec)::dut1T
    double precision,dimension(ndatemaxdc)::dut2
!AD:end
    double precision,dimension(0:ndatemax,nstRec)::ut1T
    double precision,dimension(0:ndatemaxdc)::ut2
    double precision::int
    !deleted unused logGammaJ argument %myriam
        
    kkapa=k0
    choix=0
    ig=0
    k=0
    vj=0
    n=0
    j=0
    ut1T=0.d0
    ut2=0.d0
    dut2=0.d0
    dut1T=0.d0
    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    n = (np-nva-effet-indic_ALPHA-indic_xi)/(nstRec+1)

    do jj=1,nstRec+1 !en plus strates A.Lafourcade 07/2014
        if (jj .ne. nstRec+1) then
            do i=1,n
                the1T(i-3,jj)=(bh((jj-1)*n+i))*(bh((jj-1)*n+i))
            end do
        else
            do i=1,n
                the2(i-3)=(bh(nstRec*n+i))*(bh(nstRec*n+i))
            end do
        end if
    end do
    
    if(effet.ge.1) then
        theta = bh(np-nva-indic_ALPHA-indic_xi-1)*bh(np-nva-indic_ALPHA-indic_xi-1)
        eta = bh(np-nva-indic_ALPHA-indic_xi)*bh(np-nva-indic_ALPHA-indic_xi)
        if (indic_alpha.eq.1) then ! new : joint more flexible alpha = 1 
            alpha = bh(np-nva-indic_xi)
        else
            alpha = 1.d0
        endif
        if (indic_xi.eq.1) then ! new : joint more flexible xi = 1 
            xi = bh(np-nva)
        else
            xi = 1.d0
        endif

    endif

!----------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!AD:modify
    do jj=1,nstRec !en plus strates A.Lafourcade 07/2014
        dut1T(1,jj) = (the1T(-2,jj)*4.d0/(zi(2)-zi(1)))
        ut1T(1,jj) = the1T(-2,jj)*dut1T(1,jj)*0.25d0*(zi(1)-zi(-2))
        ut1T(0,jj) = 0.d0
    end do

    dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
    ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
    ut2(0) = 0.d0

!//// NEW AMADOU vvv :
!--- strate1
    som1T = 0.d0
    vj = 0
    do i=2,ndate-1
        do k = 2,n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                    do jj=1,nstRec !en plus strates A.Lafourcade 07/2014
                        som1T(jj) = som1T(jj)+the1T(j-4,jj)
                    end do
                    vj  = j
                endif
            endif
        end do

        do jj=1,nstRec !en plus strates A.Lafourcade 07/2014
            ut1T(i,jj) = som1T(jj) +(the1T(j-3,jj)*im3(i))+(the1T(j-2,jj)*im2(i)) &
            +(the1T(j-1,jj)*im1(i))+(the1T(j,jj)*im(i))

            dut1T(i,jj) = (the1T(j-3,jj)*mm3(i))+(the1T(j-2,jj)*mm2(i)) &
            +(the1T(j-1,jj)*mm1(i))+(the1T(j,jj)*mm(i))
        end do

    end do

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

        if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3dc(i))+(the2(j-2)*im2dc(i)) &
            +(the2(j-1)*im1dc(i))+(the2(j)*imdc(i))
            dut2(i) = (the2(j-3)*mm3dc(i))+(the2(j-2)*mm2dc(i)) &
            +(the2(j-1)*mm1dc(i))+(the2(j)*mmdc(i))
        endif

    end do

!-------------fin strate2
    i = n-2
    h1 = (zi(i)-zi(i-1))

    do jj=1,nstRec !en plus strates A.Lafourcade 07/2014
        ut1T(ndate,jj)=som1T(jj)+the1T(i-4,jj)+the1T(i-3,jj)+the1T(i-2,jj)+the1T(i-1,jj)
        dut1T(ndate,jj) = (4.d0*the1T(i-1,jj)/h1)
    end do

    ut2(ndatedc)=som2+the2(i-4)+the2(i-3)+the2(i-2)+the2(i-1)!am the1(i-4)
    dut2(ndatedc) = (4.d0*the2(i-1)/h1)

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
        aux1(k)=0.d0
        aux2(k)=0.d0
    end do

    Nrec_fam = 0.d0 
    Nrec_ind = 0.d0
    Ndc_fam = 0.d0

    do i=1,nfam
        integrale3f(i) = 0.d0
    end do

!*********************************************
!-----avec un effet aleatoire dans le modele
!*********************************************

!    inv = 1.d0/theta
!    inv2 = 1.d0/eta

!ccccccccccccccccccccccccccccccccccccccccc
!     pour les donnees recurrentes
!ccccccccccccccccccccccccccccccccccccccccc

    do i=1,nsujet
        cpt(g(i))=cpt(g(i))+1
        Nrec_fam(fam(g(i))) = Nrec_fam(fam(g(i))) +1
        Nrec_ind(g(i)) = Nrec_ind(g(i)) +1
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
            res2(g(i)) = res2(g(i))+dlog(dut1T(nt1(i),stra(i))*vet)
        endif
        if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
!             print*,"here1"
            funcpajsplines_fam=-1.d9
            goto 123
        end if

!     nouvelle version
        res1(g(i)) = res1(g(i)) + ut1T(nt1(i),stra(i))*vet

        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
!             print*,"here2"
            funcpajsplines_fam=-1.d9
            goto 123
        end if

!     modification pour nouvelle vraisemblance / troncature:
        res3(g(i)) = res3(g(i)) + ut1T(nt0(i),stra(i))*vet
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
!             print*,"here3"
            funcpajsplines_fam=-1.d9
            goto 123
        end if    
    end do

!ccccccccccccccccccccccccccccccccccccccccc
! pour le deces
!ccccccccccccccccccccccccccccccccccccccccc

    do k=1,ng ! dans Joint ng=nb individus
        if(nva2.gt.0)then
            vet2 = 0.d0
            do j=1,nva2
                vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))
            end do
            vet2 = dexp(vet2)
        else
            vet2=1.d0
        endif
        if(cdc(k).eq.1)then
            Ndc_fam(fam(k)) = Ndc_fam(fam(k)) + 1
            res2dc(k) = dlog(dut2(nt1dc(k))*vet2)
            if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
                funcpajsplines_fam=-1.d9
!                 print*,'gt 1'
                goto 123
            end if
        endif

! pour le calcul des integrales / pour la survie, pas les donnees recurrentes:
        aux1(k)=ut2(nt1dc(k))*vet2
        aux2(k)=aux2(k)+ut2(nt0(k))*vet2 !vraie troncature

        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            funcpajsplines_fam=-1.d9
!             print*,'gt 2'
            goto 123
        end if
        if ((aux2(k).ne.aux2(k)).or.(abs(aux2(k)).ge. 1.d30)) then
            funcpajsplines_fam=-1.d9
            goto 123
        end if
    end do

!**************INTEGRALES ****************************
!    do ig=1,ng
!        auxig=ig
!        choix = 3
!        call gaulagJ(int,choix)
!        integrale3(ig) = int !moins bon
!    end do

!**************INTEGRALES ****************************

        choix = 3
        call gaulagJf(int,nb_gl)
        res=int
 !       write(*,*) 'funcpajslines: sum(loglik_family)', res


!************* FIN INTEGRALES **************************

 !   res = 0.d0
    do k=1,ng
        sum=0.d0
        if(cpt(k).gt.0)then
            if(theta.gt.(1.d-8)) then
!cccc ancienne vraisemblance : pour calendar sans vrai troncature cccccccc

                res= res + res2(k) &
                    + res2dc(k)  !--      pour le deces:
            else
!*************************************************************************
!     developpement de taylor d ordre 3
!*************************************************************************
!                   write(*,*)'************** TAYLOR *************'                   
                res= res + res2(k) &
                    + res2dc(k)  

            endif
            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
!                 print*,"here",k,res2(k),res2dc(k),logGammaJ(1./theta),dlog(theta),dlog(integrale3(k))
                funcpajsplines_fam=-1.d9
                goto 123
            end if
        endif

    end do

!    write(*,*) 'funcpaslines, final loglik', res
    
!res=0.d0
!do k=1,ng
!    res=res + res2(k) +res2dc(k)
!end do
!write(*,*) 'funcpaj.fam: res2', res

!write(*,*) 'funcpaj.fam: integrale',int
    
!write(*,*) 'funcpaj.fam: b', b
!write(*,*) 'funcpaj: res', res
!stop 



!---------- calcul de la penalisation -------------------

    pe1T=0.d0
    pe2 = 0.d0

    do jj=1,nstRec !en plus strates A.Lafourcade 07/2014
        do i=1,n-3
            pe1T(jj) = pe1T(jj)+(the1T(i-3,jj)*the1T(i-3,jj)*m3m3(i))+(the1T(i-2,jj) &
            *the1T(i-2,jj)*m2m2(i))+(the1T(i-1,jj)*the1T(i-1,jj)*m1m1(i))+( &
            the1T(i,jj)*the1T(i,jj)*mmm(i))+(2.d0*the1T(i-3,jj)*the1T(i-2,jj)* &
            m3m2(i))+(2.d0*the1T(i-3,jj)*the1T(i-1,jj)*m3m1(i))+(2.d0* &
            the1T(i-3,jj)*the1T(i,jj)*m3m(i))+(2.d0*the1T(i-2,jj)*the1T(i-1,jj)* &
            m2m1(i))+(2.d0*the1T(i-2,jj)*the1T(i,jj)*m2m(i))+(2.d0*the1T(i-1,jj) &
            *the1T(i,jj)*m1m(i))
        end do
    end do

    do i=1,n-3
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

    pe=0.d0
    do jj=1,nstRec !en plus strates A.Lafourcade 05/2014
        pe=pe+pe1T(jj)*k0T(jj)
    end do
    pe=pe+pe2*k0T(nstRec+1)

    resnonpen = res

    res = res - pe

    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
        funcpajsplines_fam=-1.d9
        Rrec = 0.d0
        Nrec = 0.d0
        Rdc = 0.d0
        Ndc = 0.d0
        goto 123

    else
    funcpajsplines_fam = res
        do k=1,ng
        
            Rrec(k)=res1(k)
            Nrec(k)=nig(k)
            Rdc(k)=aux1(k)
            Ndc(k)=cdc(k)
        end do
        k=0
        do i= 1,nfam
        do j=1,fsize(i)
        !     if(i.eq.1) then 
        !     write(*,*)i,fam(i),fsize(i),fsize(1)+j,aux1(fsize(1)+j)
        !    else  
        !    write(*,*)i,fam(i),fsize(i),aux1(k+j)
        !    end if
        cumulhaz1(i,j) = res1(k+j)
                cumulhaz0(i,j) = res3(k+j)
                cumulhazdc(i,j) = aux1(k+j)
                
            end do
            
           ! if(i.lt.ng.and.fam(i).ne.fam(i+1)) then 
                k = k +fsize(i)
            !end if
        end do 
    !    stop
    end if
!Ad:
123     continue

! write(*,*) 'funcpajsplines: loglik after penalty', res
    return

    end function funcpajsplines_fam




