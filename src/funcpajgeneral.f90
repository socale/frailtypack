double precision function funcpajgeneral(b,np,id,thi,jd,thj,k0)

    !use comon,only:AG,d,dmax,datedc,ndatedc,tU,nt0dc,nt1dc,nt1,ntU,stra
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
    mm3,mm2,mm1,mm,im3,im2,im1,im,date,zi,t0,t1,c,nt0,nsujet,nva,ndate, &
    nst,ve,pe,effet,nz1,nz2,ng,g,nig,resnonpen,theta,eta, &
    nva1,nva2,t0dc,t1dc,cdc,res1,res3,res4,res5,&
    vedc,aux1,aux2, auxig, indic_eta,nb_gl

    use tailles
    use comongroup

    implicit none
      integer  n,np,id,jd,i,j,k,vj,ig,gap,choix,kk
      double precision  thi,thj,pe1,pe2,sum,inv,som1,som2,logGammaJ
      double precision,dimension(2)::k0
      integer,dimension(ng)::cpt
      double precision,dimension(np)::b,bh
! --------------------------------------------------------------------
      double precision,dimension(0:ndate)::dut1,dut2
      double precision,dimension(0:ndate)::ut1,ut2
      double precision,dimension(ng)::res1dc,res2dc,res3dc,res4dc,res5dc
      double precision ,dimension(ng)::res2,integrale4
      double precision ,dimension(ng)::suLIST, lamLIST
      double precision,dimension(ng)::integrale1,integrale2,integrale3,integrale3gap

      double precision  res, h1
      double precision lam,gl
      double precision int,int4,int3gap

! ***** indicateur de troncature
      integer :: indictronq,indictronqdc ! =0 si donnees non tronquées reellement
      common /troncature/indictronq,indictronqdc

! ***********************************************************
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                     initializing values
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        indic_eta=1  !!! FB comment
        vet=0.d0
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
        bh=b
        nst=2

        lam = 0.d0      !!! FB comment
        gl = 0.d0      !!! FB comment

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !print *,"             .... ON EST DANS FUNGPAJGENERAL ici ...."
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do i=1,np
            bh(i)=b(i)
         end do

         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj

         n = (np-nva-effet-indic_eta)/nst !n=12
        !print*,"np = ", np          !=32
        !print*,"nva = ", nva        !=6
        !print*,"nst = ", nst        !=2
        !print*,"effet = ", effet    !=1
        !stop
        do kk=1, n+3           !!! FB comment DO to delete ! !!!!!!!!!!!!!!!!
            the1(kk)=0.d0      !!! FB comment DO to delete ! !!!!!!!!!!!!!!!!
            the2(kk)=0.d0      !!! FB comment DO to delete ! !!!!!!!!!!!!!!!!
        end do                 !!! FB comment DO to delete ! !!!!!!!!!!!!!!!!


        do i=1,n
            the1(i-3)=(bh(i))*(bh(i))
            j = n+i
            if (nst.eq.2) then
                the2(i-3)=(bh(j))*(bh(j))
            endif
        end do

        if(effet.eq.1) then
            theta = bh(np-nva-indic_eta)*bh(np-nva-indic_eta)
            eta = bh(np-nva)*bh(np-nva)
        endif

!  --------  calcul de ut1(ti) et ut2(ti) ---------------------------
!       attention the(1)  sont en nz=1
!       donc en ti on a the(i)

         vj = 0
         som1 = 0.d0
         som2 = 0.d0
         dut1(1) = (the1(-2)*4.d0/(zi(2)-zi(1)))
         dut2(1) = (the2(-2)*4.d0/(zi(2)-zi(1)))
         ut1(1) = the1(-2)*dut1(1)*0.25d0*(zi(1)-zi(-2))


         ut2(1) = the2(-2)*dut2(1)*0.25d0*(zi(1)-zi(-2))
         ut1(0) = 0.d0
         ut2(0) = 0.d0

         do i=2,ndate-1
            do k = 2,n-2
               if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                  j = k-1
                  if ((j.gt.1).and.(j.gt.vj))then
                     som1 = som1 + the1(j-4)
                     som2 = som2 + the2(j-4)
                     vj  = j
                  endif
               endif
            end do

            ut1(i) = som1 +(the1(j-3)*im3(i))+(the1(j-2)*im2(i)) &
                +(the1(j-1)*im1(i))+(the1(j)*im(i))
            dut1(i) = (the1(j-3)*mm3(i))+(the1(j-2)*mm2(i)) &
                +(the1(j-1)*mm1(i))+(the1(j)*mm(i))

            if(nst.eq.2)then
            ut2(i) = som2 +(the2(j-3)*im3(i))+(the2(j-2)*im2(i)) &
                +(the2(j-1)*im1(i))+(the2(j)*im(i))
            dut2(i) = (the2(j-3)*mm3(i))+(the2(j-2)*mm2(i)) &
                +(the2(j-1)*mm1(i))+(the2(j)*mm(i))
            endif

        end do
        i = n-2
        h1 = (zi(i)-zi(i-1))
        ut1(ndate)=som1+the1(i-4)+the1(i-3)+the1(i-2)+the1(i-1)
        ut2(ndate)=som2+the1(i-4)+the2(i-3)+the2(i-2)+the2(i-1)
        dut1(ndate) = (4.d0*the1(i-1)/h1)
        dut2(ndate) = (4.d0*the2(i-1)/h1)


!  -------------------------------------------------------
!  --------- calcul de la vraisemblance ------------------
!  -------------------------------------------------------
!  --------- avec ou sans variable explicative  ----------

        do k=1,ng

            res1(k) = 0.d0
            res2(k) = 0.d0
            res3(k) = 0.d0
            res4(k) = 0.d0
            res5(k) = 0.d0
            res1dc(k) = 0.d0
            res2dc(k) = 0.d0
            res3dc(k) = 0.d0
            res4dc(k) = 0.d0
            res5dc(k) = 0.d0
            cpt(k) = 0
            integrale1(k) = 0.d0
            integrale2(k) = 0.d0
            integrale3(k) = 0.d0
            integrale4(k) = 0.d0
            integrale3gap(k) = 0.d0
            aux1(k)=0.d0
            aux2(k)=0.d0
            suLIST(k)=0.d0
            lamLIST(k)=0.d0
        end do

        !print*,"fin calcul des variables expliquatives"
! ******************************************
!  ---avec un effet aleatoire dans le modele
! ********************************************

         inv = 1.d0/theta
        !nsujet = 694

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     pour les donnees recurrentes          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do i=1,nsujet
            !print *, "i = ", i
            !print* ,"  cpt(g(i)) = ", cpt(g(i))
            cpt(g(i))=cpt(g(i))+1
            if(nva1.gt.0)then
               vet = 0.d0

               do j=1,nva1
                  vet =vet + bh(np-nva+j)*dble(ve(i,j))
               end do

               vet = exp(vet)
            else
               vet=1.d0
            endif

            if((c(i).eq.1))then

                call suspJ(t1(i),the1,nz1,lam,gl,zi)
                 res2(g(i)) =res2(g(i))+ log(lam*vet)
            endif
!     nouvelle version
            call suspJ(t1(i),the1,nz1,lam,gl,zi)
            res1(g(i)) = res1(g(i)) + gl*vet
            call suspJ(t0(i),the1,nz1,lam,gl,zi)
            res3(g(i)) = res3(g(i)) + gl*vet
         end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pour le deces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do k=1,ng
             if(nva2.gt.0)then
                vet2 = 0.d0
                do j=1,nva2
                   vet2 =vet2 + bh(np-nva2+j)*dble(vedc(k,j))

                end do
                vet2 = exp(vet2)
             else
                vet2=1.d0
             endif
             if(cdc(k).eq.1)then
!                res2dc(k) = log(dut2(nt1dc(k))*vet2)
                call suspJ(t1dc(k),the2,nz2,lam,gl,zi)
                res2dc(k) = log(lam*vet2)
             endif

! pour le calcul des integrales / pour la survie, pas les données recurrentes:
!             aux1(k)=ut2(nt1dc(k))*vet2
             call suspJ(t1dc(k),the2,nz2,lam,gl,zi)
             aux1(k) = gl*vet2
            if (i.le.(size(g))) then
                aux2(g(i))=aux2(g(i))+ut2(nt0(i))*vet2 !vraie troncature
            end if
! yas RES4 CORRESPOND AU TERME : int_0^(L_i)\lambda_i(t)dt
!             res4(k) = ut2(nt0dc(k))*vet2theta
             call suspJ(t0dc(k),the2,nz2,lam,gl,zi)
             res4(k) = gl*vet2
             call suspJ(t0dc(k),the1,nz1,lam,gl,zi)
             res5(k) = gl*vet
        end do


! *************INTEGRALES ****************************
          do ig=1,ng
                 auxig=ig
                 choix = 3
                 call gaulagJ(int,choix,nb_gl)
                 integrale3(ig) = int !moins bon

                 choix = 4
                 call gaulagJ(int4,choix,nb_gl)
                 integrale4(ig) = int4

                 choix = 5
                 call gaulagJ(int3gap,choix,nb_gl)
                 integrale3gap(ig) = int3gap
           end do
! ************ FIN INTEGRALES **************************

          gap =0

          res = 0.d0
          do k=1,ng
             sum=0.d0
             if(cpt(k).gt.0)then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!     JOINT FRAILTY MODEL 2010 SANS eta     !!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(gap.eq.0) then

                       if(t0dc(k).ne.0) then
                             res= res + res2(k) &
                                + res2dc(k)+ logGammaJ(1./theta+nig(k)+cdc(k)) &
                                - logGammaJ(1./theta) &
                                + (nig(k)+cdc(k))*log(theta) &
                                + log(integrale3(k)) &
                                - log(integrale4(k))
                        else
                              res= res + res2(k)+ res2dc(k) &
                                + logGammaJ(1./theta+nig(k)+cdc(k)) &
                                - 1./eta*log(eta) &
                                - logGammaJ(1./eta)- logGammaJ(1./theta) &
                                + (nig(k)+cdc(k))*log(theta) &
                                + log(integrale3(k))
                        endif

                else
                      res= res + res2(k) &
                        + res2dc(k) &
                        - logGammaJ(1./eta) &
                        + logGammaJ(1./theta+nig(k)+cdc(k)) &
                        - logGammaJ(1./theta)-log(eta)/eta &
                        + (nig(k)+cdc(k))*log(theta) &
                        + log(integrale3gap(k))
                endif
            endif
          end do


!  -------- calcul de la penalisation -------------------

         pe1 = 0.d0
         pe2 = 0.d0
         do i=1,n-3
            pe1 = pe1+(the1(i-3)*the1(i-3)*m3m3(i))+(the1(i-2) &
                *the1(i-2)*m2m2(i))+(the1(i-1)*the1(i-1)*m1m1(i))+( &
                the1(i)*the1(i)*mmm(i))+(2.d0*the1(i-3)*the1(i-2)* &
                m3m2(i))+(2.d0*the1(i-3)*the1(i-1)*m3m1(i))+(2.d0* &
                the1(i-3)*the1(i)*m3m(i))+(2.d0*the1(i-2)*the1(i-1)* &
                m2m1(i))+(2.d0*the1(i-2)*the1(i)*m2m(i))+(2.d0*the1(i-1) &
                *the1(i)*m1m(i))
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

          resnonpen=res

          res = res - pe


          funcpajgeneral = res

          return
          end function funcpajgeneral
