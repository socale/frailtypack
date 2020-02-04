    
    
    
    !========================          FUNCPAJ_SPLINES         ====================
        double precision function funcpalongi_uni(b,np,id,thi,jd,thj,k0)
    
        use donnees_indiv
        use lois_normales
        use tailles
        use comon
   !     use ParametresPourParallelisation
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
        !double precision,dimension(2):: result, abserr2
        !double precision,dimension(1000) :: work
        external :: vraistot,vraistot_splines,vraistot_weib
        double precision,dimension(nea) :: xea
        !integer ::neval,ifail
        !
        integer::n,i,j,k,vj,ig,choix
        !integer,dimension(ngmax)::cpt
        double precision::sum,res
        double precision :: eps_s
    
        !double precision,dimension(-2:npmax):: the1,the2
        double precision,dimension(np)::bh
        double precision,dimension(ngmax):: &
        integrale4
    !AD: for death,change dimension
        !double precision,dimension(ndatemax)::dut1
        !double precision,dimension(ndatemaxdc)::dut2
    !AD:end
        !double precision,dimension(0:ndatemax)::ut1
        !double precision,dimension(0:ndatemaxdc)::ut2
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
 
        do i=1,np
        bh(i)=b(i)
        end do
    
        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj
    
            b1 = bh
    
 
        K_G0 =  bh(np-nva-nb_re-1-3)
        K_D0 =  bh(np-nva-nb_re-1-2)
        lambda =  bh(np-nva-nb_re-1-1)
    
        y0 = bh(np-nva-nb_re-1)
          
      
            sigmae = bh(np-nva-nb_re)*bh(np-nva-nb_re)
    

            Ut = 0.d0
            Utt = 0.d0
    
                    do j=1,nb1
                do k=1,j
                
                       Ut(j,k)=bh(np-nva-nb_re+k)!sqrt(bh(np-nva-nb_re+k)**2.d0)
                  Utt(k,j)=bh(np-nva-nb_re+k)!sqrt(bh(np-nva-nb_re+k)**2.d0)
                  
          
                end do
            end do
    
  
        
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
    

    

    !**************INTEGRALES ****************************

    
        res = 0.d0
            it = 0
            it_rec = 1
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

                it_cur = it
    
                x2 = 0.d0
                x2cur = 0.d0
                z1cur = 0.d0
                current_mean = 0.d0
                          
        
   mu = 0.d0

            mu(1:nmescur,1) = matmul(vey((it+1):(it+nmescur),1:(nva3)),bh((np-nva4-nva3+1):(np-nva4)))
    
            mu(1:nmescur,2) = matmul(vey((it+1):(it+nmescur),(nva3+1):(nva3+nva4)),bh((np-nva4+1):np))
            xea = 0.d0
 

    select case(which_random)
        case(1)
            call gauherJ1_uni(int)
        case(2)
            call gauherJ2_uni(int)
        case(3)
            call gauherJ3_uni(int)
        case(4)
            call gauherJ4_uni(int)
        case(5)
            call gauherJ5_uni(int)
        case(6)
            call gauherJ6_uni(int)
        case(7)
            call gauherJ7_uni(int)
        case(8)
            call gauherJ8_uni(int)
        case(9)
            call gauherJ9_uni(int)
        case(10)
            call gauherJ10_uni(int)
        case(11)
            call gauherJ11_uni(int)
        case(12)
            call gauherJ12_uni(int)
        case(13)
            call gauherJ13_uni(int)
        case(14)
            call gauherJ14_uni(int)
        case(15)
            call gauherJ15_uni(int)        
    end select    
                integrale4(ig) =int !result(1) !

 
        
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

                        res= res -0.5d0*dlog(2.d0*pi)*nmes_o(ig) -0.5d0*dlog(sigmae)*nmes_o(ig)  &
                            -112.d0
            else
                        res= res - 0.5d0*dlog(2.d0*pi)*nmes_o(ig) -0.5d0*dlog(sigmae)*nmes_o(ig)  &
                            +dlog(integrale4(ig))
                                endif

            if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                funcpalongi_uni=-1.d9
                goto 123
            end if


        end do


    !---------- calcul de la penalisation -------------------
    
    
        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
            funcpalongi_uni=-1.d9
         
            
            goto 123
    
        else
    
            funcpalongi_uni = res
          
        end if
    !Ad:
        123     continue
    
        return
    
    
        end function funcpalongi_uni
    

    
    
    
     !========================          FUNCPAJ_SPLINES         ====================
        double precision function funcprep(b,np)
    
        use donnees_indiv
        use comon
 !  use lois_normales
  !    use tailles
     !      use residusM
            use optim
    !        use alnorm
        IMPLICIT NONE
    
    ! *** NOUVELLLE DECLARATION F90 :
    
        integer,intent(in)::np
        double precision,dimension(np),intent(in)::b
      
        ! for the numerical integral hrmsym
        integer :: ier
        !
        integer::i,j,k,jj
            double precision,dimension(np)::bh
          double precision::eps
        double precision,parameter::pi=3.141592653589793d0
       

    double precision::finddet
       double precision,dimension(nea*(nea+1)/2)::matv
        
        
    !    print*,'debut funcpa'
       
        do i=1,np
        bh(i)=b(i)
        end do
    npp = np
        b1 = bh
    
        K_G0 =  bh(np-nva-nb_re-1-3)
        K_D0 =  bh(np-nva-nb_re-1-2)
        lambda =  bh(np-nva-nb_re-1-1)
    
        y0 = bh(np-nva-nb_re-1)
          
      
         
            sigmae = bh(np-nva-nb_re)*bh(np-nva-nb_re)
    

            Ut = 0.d0
            Utt = 0.d0
    
                    do j=1,nb1
                do k=1,j
               
                       Ut(j,k)=bh(np-nva-nb_re+k)!sqrt(bh(np-nva-nb_re+k)**2.d0)
                  Utt(k,j)=bh(np-nva-nb_re+k)!sqrt(bh(np-nva-nb_re+k)**2.d0)
            
                end do
            end do
    
   
        
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
    
    funcprep = 0.d0
    
    
        end function funcprep
    