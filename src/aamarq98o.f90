

!===============================    MARQ98AUX =========================
!

     
    subroutine marq98o(b,m,ni,v,rl,ier,istop,fctnames)
    
    use tailles
    use optim
    use parameters,only:maxiter
    !use comon,only:t0,t1,c,nt0,nt1,nsujet,nva,ndate,nst,auxig
    
    Implicit none
    
    integer::m,ni,nql,i,ii,nfmax,idpos,ier,istop,igrad,ncount,id,jd,i0
    double precision::da,dm,ga,tr,ca,cb,epsa,epsb,rl &
    ,step,eps,epsd,vw,fi,z,rl1,th,ep,dd,auxmax     
    double precision,dimension(m*(m+3)/2)::V,fu
    double precision,dimension(m)::b,delta,bh,b1
    double precision,dimension(2)::zero
    double precision::maxtadd
    double precision,external::fctnames
    
    zero(1)=0.d0
    zero(2)=0.d0
    id=0
    jd=0
    z=0.d0
    th=1.d-5
    eps=1.d-7
    epsa=1.d-6
    epsb=1.d-6
    epsd=1.d-6
    
    nfmax=m*(m+1)/2
    ca=epsa+1.d0
    cb=epsb+1.d0
    rl1=-1.d+10
    ni=0
    istop=0
    da=0.01d0
    dm=5.d0
    
    nql=1
    
10   continue

        
    z=0.d0
    i0=0     
    rl=fctnames(b,m,i0,z,i0,z)
!    write(*,*)'iter',ni,'vrais',rl1
    rl1=rl 
    
    call derivao(b,m,v,rl,fctnames)
    
    dd = 0.d0
    do i=m*(m+1)/2+1,m*(m+3)/2
        dd = dd + v(i)*v(i)
    end do
    
    dd=dd/dabs(RL)

    if (ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) goto 100
    
        tr=0.d0
        
        do i=1,m
            ii=i*(i+1)/2
            tr=tr+dabs(v(ii))
        end do 
        
        tr=tr/dble(m)
        ncount=0
        ga=0.01d0
        
400          do i=1,nfmax+m
            fu(i)=v(i)
        end do
        
        do i=1,m
            ii=i*(i+1)/2
            if (v(ii).ne.0) then
                fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
            else
                fu(ii)=da*ga*tr
            end if
        end do
        
        call dcholej(fu,m,nql,idpos)
    
        if (idpos.ne.0) then  
    
            if(b(1).lt.-1.d0.or.b(1).gt.1.d0.or.b(2).lt.-1.d0.or.b(2).gt.1.d0) then
                b(1)=0.5d0
                b(2)=0.5d0
                goto 110
            endif
                
            ncount=ncount+1
            
            if (ncount.le.3.or.ga.ge.1.d0) then
                da=da*dm
            else
                ga=ga*dm
                if (ga.gt.1.d0) ga=1.d0
            end if
            
            if (ncount > 10) then
                fu=0.d0
                do i=1,m
                    ii=i*(i+1)/2
                    fu(ii)=1
                end do
                return
            end if
            
            goto 400
            
        else
            do i=1,m
                delta(i)=fu(nfmax+i)
                b1(i)=b(i)+delta(i)
            end do
    
            rl=fctnames(b1,m,id,z,jd,z)
            
            igrad=1
            
            if (rl1.lt.rl) then
                if(da.lt.eps) then
                    da=eps
                else
                    da=da/(dm+2.d0)
                end if
                goto 800
            end if
        end if
        
        call dmaxt(maxtadd,delta,m)
        auxmax=maxtadd
        
        if(auxmax.eq.0.d0) then
            vw=1.D-5
        else
        !    call dmaxt(maxtadd,delta,m)
            vw=th/auxmax
        end if

        step=dlog(1.5d0)
        
        call searpaso(vw,step,b,bh,m,delta,fi,fctnames) 
        
        rl=-fi
        
        IF(rl1.gt.rl) then
        end if
        
        do i=1,m
            delta(i)=vw*delta(i)
        end do
        
        da=(dm-3.d0)*da  
    
800       cb=dabs(rl1-rl)

        ca=0.d0
        
        do i=1,m
            ca=ca+delta(i)*delta(i)
        end do

        do i=1,m
            b(i)=b(i)+delta(i)
        end do

        ni=ni+1
        
        if (ni.gt.maxiter) then
            istop=2    
            goto 110
        end if
        goto 10
    
100     continue

    istop=1
    
    ep=10.d-10
    
    call dsinvj(v,m,ep,ier)    
    if (ier.eq.-1) then
        istop=3
    end if    


110       continue

    return
    
    end subroutine marq98o
         


    subroutine derivao(b,m,v,rl,fctnames)

    implicit none
    
    integer,intent(in)::m
    double precision,intent(inout)::rl
    double precision,dimension(m),intent(in)::b
    double precision,dimension((m*(m+3)/2)),intent(out)::v     
    double precision,dimension(m)::fcith
    integer ::i0,m1,ll,i,k,j,iun
    double precision::fctnames,thn,th,z,vl,th2,vaux  
    external::fctnames

    th=5.d-3
    
    thn=-th
    th2=th*th
    z=0.d0
    i0=0
    iun =1

    rl=fctnames(b,m,iun,z,iun,z)
        if(rl.eq.-1.d9) then
            rl=-1.d9
            goto 123
        end if    
    do i=1,m
        fcith(i)=fctnames(b,m,i,th,i0,z)
                if(fcith(i).eq.-1.d9) then
                    rl=-1.d9
                    goto 123
                end if
    end do

    k=0
    m1=m*(m+1)/2
    ll=m1
!    
    do i=1,m
        ll=ll+1
        vaux=fctnames(b,m,i,thn,i0,z)
                if(vaux.eq.-1.d9) then
                    rl=-1.d9
                    goto 123
                end if    
        vl=(fcith(i)-vaux)/(2.d0*th)
        v(ll)=vl
        do j=1,i
            k=k+1
            v(k)=-(fctnames(b,m,i,th,j,th)-fcith(j)-fcith(i)+rl)/th2
        end do
    end do
      
123   continue    
    return
    
    end subroutine derivao

!================================  SEARPAS joly    ==============================

    subroutine searpaso(vw,step,b,bh,m,delta,fim,fctnames)
    
    implicit none
    
    integer::m,i      
    double precision,dimension(m)::b
    double precision::vw
    double precision,dimension(m)::bh,delta    
    double precision::fim,step   
    double precision::vlw,vlw1,vlw2,vlw3,vm,fi1,fi2,fi3    
    double precision,external::fctnames

    
    vlw1=dlog(vw)
    vlw2=vlw1+step
    call valfpao(vlw1,fi1,b,bh,m,delta,fctnames)
    call valfpao(vlw2,fi2,b,bh,m,delta,fctnames)       
    
    if(fi2.ge.fi1) then
        vlw3=vlw2
        vlw2=vlw1
        fi3=fi2
        fi2=fi1
        step=-step
    
        vlw1=vlw2+step
        call valfpao(vlw1,fi1,b,bh,m,delta,fctnames)   
        if(fi1.gt.fi2) goto 50
    else 
        vlw=vlw1
        vlw1=vlw2
        vlw2=vlw
        fim=fi1
        fi1=fi2
        fi2=fim
    end if
    
    do i=1,40
        vlw3=vlw2
        vlw2=vlw1
        fi3=fi2
        fi2=fi1
    
        vlw1=vlw2+step
        call valfpao(vlw1,fi1,b,bh,m,delta,fctnames)
        if(fi1.gt.fi2) goto 50
        if(fi1.eq.fi2) then
        fim=fi2
        vm=vlw2 
        goto 100
        end if
    end do
    
50     continue
    
    vm=vlw2-step*(fi1-fi3)/(2.d0*(fi1-2.d0*fi2+fi3))   
    call valfpao(vm,fim,b,bh,m,delta,fctnames)    
    if(fim.le.fi2) goto 100
    vm=vlw2
    fim=fi2
100   continue
    vw=dexp(vm)
    
    return
    
    end subroutine searpaso

   
!========================   VALFPAAUX   ============================== 
       
    subroutine valfpao(vw,fi,b,bk,m,delta,fctnames)
    
    implicit none
    
    integer::m  
    double precision,dimension(m)::b,delta  
    double precision,dimension(m)::bk 
    double precision::fi 
    double precision::vw,z    
    integer::i0,i
    double precision,external::fctnames
    
    z=0.d0
    i0=1
    do i=1,m
    bk(i)=b(i)+dexp(vw)*delta(i)
    end do
    fi=-fctnames(bk,m,i0,z,i0,z)
    
    return
        
    end subroutine valfpao

!========================    DEBUT FUNCPAAUX       ====================
    double precision function funcpao(b,np,id,thi,jd,thj)
    
    use tailles
    !use comon,only:t0,t1,nt0,nst,ndate
    use comon,only:g,c,nt1,nsujet,nva,stra,ve,auxig
    !use comon,only:nig
    !use additiv,only:dut1,dut2,ve2,rho
    use additiv,only:ut1,ut2,betaaux,sigma2,tau2,cov
    implicit none
    
    integer::np,id,jd,i,k,ip
    double precision,dimension(np)::bhaux,b
    double precision::thi,thj,res,vet
!****** u_tilde
    double precision  :: u_tilde,v_tilde
!****** derivanal
    double precision , dimension(ngmax)::res3,res4
    double precision , dimension(ngmax)::res5,res6,res8

!==============================================
!================POUR UN GROUPE AUXIG donn� !
!==============================================
    
    res3=0.d0
    res4=0.d0
    res5=0.d0
    res6=0.d0
    res8=0.d0
    
    do i=1,np
        bhaux(i)=b(i)
    end do 
    
    
    if (id.ne.0) bhaux(id)=bhaux(id)+thi 
    if (jd.ne.0) bhaux(jd)=bhaux(jd)+thj    
        
    u_tilde=bhaux(1)
    v_tilde=bhaux(2)
!--------------------------------------------------------
!----------calcul de la vraisemblance ------------------
!---------------------------------------------------------
    
    do k=1,nsujet
        if(nva.gt.0.and.g(k).eq.auxig)then
            vet = 0.d0 
            do ip=1,nva
                vet =vet + betaaux(ip)*ve(k,ip)
            end do
            vet = dexp(vet)
        else
            vet=1.d0
        endif
        
        if(g(k).eq.auxig)then
        
            if(c(k).eq.1)then
                res3(auxig) = res3(auxig)+u_tilde+v_tilde*ve(k,1)
                res8(auxig) = res8(auxig)+ve(k,1)
            endif
            if(stra(k).eq.1)then
                res4(auxig) =res4(auxig)+ut1(nt1(k))*vet*dexp(u_tilde+v_tilde*ve(k,1))
                
                res5(auxig) = res5(auxig)+ut1(nt1(k))*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)
                
                res6(auxig) = res6(auxig)+ut1(nt1(k))*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
            endif
            if(stra(k).eq.2)then
                res4(auxig) = res4(auxig) &
                +ut2(nt1(k))*vet*dexp(u_tilde+v_tilde*ve(k,1))
                
                res5(auxig) = res5(auxig)+ut2(nt1(k))*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*ve(k,1)
                
                res6(auxig) = res6(auxig)+ut2(nt1(k))*vet* &
                dexp(u_tilde+v_tilde*ve(k,1))*(ve(k,1))**2
            endif
        
        endif 
    end do 
        
    res=-res3(auxig) + res4(auxig)+0.5d0*(((u_tilde)**2)/sigma2+((v_tilde)**2)/tau2 &
    -2.d0*u_tilde*v_tilde*cov/(sigma2*tau2))/(1.d0-(cov**2)/(sigma2*tau2))         !-ka
    
    funcpao= -res
    
    return
    
    end function funcpao


    logical function isnann(x)

    implicit none
    
    double precision,intent(in)::x
    
    if (x .ne. x) then
        isnann=.true.
    else
        isnann=.false.
    end if

    end function isnann
