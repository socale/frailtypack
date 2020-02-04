    
    
    
    !========================          FUNCPAJ_SPLINES         ====================
        double precision function funcpajlongisplines_nl(b,np,id,thi,jd,thj,k0)
    
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
        INTEGER ::MINPTS, MAXPTS

        double precision:: epsabs,epsrel
        !double precision,dimension(2):: result, abserr2
        !double precision,dimension(1000) :: work
        external :: vraistot,vraistot_splines,vraistot_weib
        double precision,dimension(nea) :: xea
        !integer ::neval,ifail
        !
        integer::n,i,j,k,vj,ig,choix
        integer,dimension(ngmax)::cpt
        double precision::pe1,pe2,sum,som1,som2,res,vet,h1
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
       
        external:: func7j
    double precision::func7j,finddet
       double precision,dimension(nea*(nea+1)/2)::matv
       !double precision,dimension(nsujety,1)::vey_bh
        
 
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
    
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar-4)/(effet+1) !nst        !to znaczy ze dzielimy lliczbe wezlow na 2


        n_wezly = n
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
        


        K_G0 =  bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)
        K_D0 =  bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)
        lambda =  bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)
        y0 = bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha)
        
    
        
        
        if(typeJoint.eq.3) then
            sigmav = bh(np-nva-nb_re-1-netadc - netar-indic_alpha)
            alpha = bh(np-nva-nb_re-1-netadc - netar)
            if(indic_alpha.eq.0)alpha=1.d0
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
    
                funcpajLongisplines_nl=-1.d9
                goto 123
            end if
   !     nouvelle version
   
        res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet  
        if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
             funcpajLongisplines_nl=-1.d9
            goto 123
        end if             
!     modification pour nouvelle vraisemblance / troncature:
        res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet 
        if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge. 1.d30)) then
             funcpajLongisplines_nl=-1.d9
            goto 123
        end if    
    
           ! Rrec(g(i)) = Rrec(g(i)) + ut1(nt1(i))*vet
    
    
        end do
        end if
    
     
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc

    vet22(1:ng) = MATMUL(dble(vedc(1:ng,1:nva2)),bh(np-nva4-nva3-nva2+1:(np-nva3-nva4)))
    vet22 = dexp(vet22)

        do k=1,ng
    

    
            if(cdc(k).eq.1)then
    
                res2dc(k) = dlog(dut2(nt1dc(k))*vet22(k))
        
                if ((res2dc(k).ne.res2dc(k)).or.(abs(res2dc(k)).ge. 1.d30)) then
    
                    funcpajLongisplines_nl=-1.d9
                ! print*,'gt 1'
                    goto 123
                end if
            endif
    
    
            aux1(k)=ut2(nt1dc(k))*vet22(k)
            
    !    aux2(k)=aux2(k)+ut2(nt0(k))*vet2 !vraie troncature
        
        if ((aux1(k).ne.aux1(k)).or.(abs(aux1(k)).ge. 1.d30)) then
            
                    funcpajLongisplines_nl=-1.d9
!            print*,'gt 2'
            goto 123
        end if    
        if ((aux2(k).ne.aux2(k)).or.(abs(aux2(k)).ge. 1.d30)) then
        
                    funcpajLongisplines_nl=-1.d9
!            print*,'gt 3'
!            print*,k,'(aux2(k)',aux2(k),'nt0(k)',nt0(k),'vet2',vet2,'ut2(nt0(k))',ut2(nt0(k))
            goto 123
        end if    
    
    
                RisqCumul(k) = ut2(nt1dc(k))*vet22(k)

                Rdc_res(k) = ut2(nt1dc(k))!*vet2
    
        end do

        if(effet.eq.0.or.link.eq.2) then
                res2= 0.d0
    
            end if
    
    !**************INTEGRALES ****************************


    
        res = 0.d0
            it = 0
            it_rec = 1
            MINPTS=30
         MAXPTS=500               

            epsabs = 1.d-100
            epsrel = 1.d-100
            restar = 0
        !    nmes_o = 0
            nf2 = nf
            integrale4 = 0.d0
    
            do ig=1,ng
    
           !     ycurrent  = 0.d0
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
    
                if(typeJoint.eq.3) then
                do i= 0,nmescurr-1
    
                res1cur(i+1) = ut1(nt1(it_rec+i))
                    res3cur(i+1)  = ut1(nt0(it_rec+i))
                res2cur(i+1) = dut1(nt1(it_rec+i))
                end do
            end if
    
  
        
    ut2cur = ut2(nt1dc(ig))
 
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
    
      
     
    
  
    !************* FIN INTEGRALES **************************
   
    
            sum=0.d0
   
    !*************************************************************************
    !     vraisemblnce
    !*************************************************************************

            if (integrale4(ig).le.0.d0) then
    
                        res= res +res2(ig) &
                            + res2dc(ig)-0.5d0*dlog(2.d0*pi)*nmes_o(ig) -0.5d0*dlog(sigmae)*nmes_o(ig)  &
                            -112.d0
            else
                        res= res + res2(ig) &
                            + res2dc(ig)- 0.5d0*dlog(2.d0*pi)*nmes_o(ig) -0.5d0*dlog(sigmae)*nmes_o(ig)  &
                            +dlog(integrale4(ig))
            end if            
       

            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpajLongisplines_nl=-1.d9
                goto 123
            end if



        end do
    

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
  
                pe2 = pe2+(the2(i-3)*the2(i-3)*m3m3(i))+(the2(i-2) &
                *the2(i-2)*m2m2(i))+(the2(i-1)*the2(i-1)*m1m1(i))+( &
                the2(i)*the2(i)*mmm(i))+(2.d0*the2(i-3)*the2(i-2)* &
                m3m2(i))+(2.d0*the2(i-3)*the2(i-1)*m3m1(i))+(2.d0* &
                the2(i-3)*the2(i)*m3m(i))+(2.d0*the2(i-2)*the2(i-1)* &
                m2m1(i))+(2.d0*the2(i-2)*the2(i)*m2m(i))+(2.d0*the2(i-1) &
                *the2(i)*m1m(i))
                
      
        end do
 
        pe = k0(1)*pe1 + k0(2)*pe2
        resnonpen = res
   
        res = res - pe
  
      do k=1,ng
                            if(typeJoint.eq.3) then 
                            Rrec(k)=res1(k)-res3(k)
                            Nrec(k)=nig(k)
                            end if
                Rdc(k)=RisqCumul(k)
                Ndc(k)=cdc(k)
            end do
            

        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpajLongisplines_nl=-1.d9
   !      

            goto 123
    
        else
    
            funcpajLongisplines_nl = res
          
        end if
    !Ad:
        123     continue
    
        return
    
    
        end function funcpajlongisplines_nl
    
  
    
        
    !========================          FUNCPAJ_SPLINES         ====================
        double precision function funcsplines_nl(b,np,id,thi,jd,thj,k0)
    
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
        integer :: jj,ier
        !double precision:: epsabs,epsrel
        !double precision,dimension(2):: result, abserr2
        !double precision,dimension(1000) :: work
        external :: vraistot,vraistot_splines,vraistot_weib
        !double precision,dimension(nea) :: xea
        !integer ::neval,ifail
        !
        integer::n,i,j,k,vj,ig,choix
        integer,dimension(ngmax)::cpt
        double precision::som1,som2,res,vet,h1
        double precision :: eps_s
    
        double precision,dimension(-2:npmax):: the1,the2
        double precision,dimension(np)::bh
        double precision,dimension(ngmax)::res2,res1dc,res2dc &
        ,res3dc,integrale1,integrale2,integrale3
    !AD: for death,change dimension
        double precision,dimension(ndatemax)::dut1
        double precision,dimension(ndatemaxdc)::dut2
    !AD:end
        double precision,dimension(0:ndatemax)::ut1
        double precision,dimension(0:ndatemaxdc)::ut2
        double precision::eps
        double precision,parameter::pi=3.141592653589793d0
       
        external:: func7j
    double precision::func7j,finddet
       double precision,dimension(nea*(nea+1)/2)::matv
       !double precision,dimension(nsujety,1)::vey_bh
        

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
    
  

            b1 = bh
    
        n = (np-nva-effet-indic_ALPHA-1-nb_re - netadc - netar-4)/(effet+1) !nst        !to znaczy ze dzielimy lliczbe wezlow na 2

        n_wezly = n
        
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
    
    
    K_G0 =  bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-3)
        K_D0 =  bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-2)
        lambda =  bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha-1)
        y0 = bh(np-nva-nb_re-1-netadc - netar-effet-indic_alpha)
        
        
        
        if(typeJoint.eq.3) then
            sigmav = bh(np-nva-nb_re-1-netadc - netar-indic_alpha)
            alpha = bh(np-nva-nb_re-1-netadc - netar)
                if(indic_alpha.eq.0)alpha=1.d0
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
    
    
        
    
            res1(g(i)) = res1(g(i)) + ut1(nt1(i))*vet
    res3(g(i)) = res3(g(i)) + ut1(nt0(i))*vet 
    
        end do
        end if
    
    
    !ccccccccccccccccccccccccccccccccccccccccc
    ! pour le deces
    !ccccccccccccccccccccccccccccccccccccccccc
 
    vet22(1:ng) = MATMUL(dble(vedc(1:ng,1:nva2)),bh(np-nva4-nva3-nva2+1:(np-nva3-nva4)))
    vet22 = dexp(vet22)

        do k=1,ng
    
  
    
                RisqCumul(k) = ut2(nt1dc(k))*vet22(k)
                Rdc_res(k) = ut2(nt1dc(k))!*vet2
    
        end do

     

     
        res = 0.d0


    
            funcsplines_nl = res
            do k=1,ng
                            if(typeJoint.eq.3) then 
                            Rrec(k)=res1(k)!-res3(k)
                            Nrec(k)=nig(k)
                            end if
                Rdc(k)=RisqCumul(k)
                Ndc(k)=cdc(k)
            end do

    !Ad:
        123     continue
    
        return
    
    
        end function funcsplines_nl