ccccccc
        module rmatmod
                use parameter
                use channels
                use coulfunc
                use deltaf
                use coulvar
                use wtkvar
                use rmatvar
                use mesh
                use system
                use potential
                use potvar
                use clebsg
                use deform
                implicit none
        contains
ccccccc
        subroutine rmat_int()
                implicit none
ccccccc
                integer::i,k       
                integer::li      
ccccccc
                if(allocated(WTK)) deallocate(WTK)
                if(allocated(WTKP)) deallocate(WTKP)
                if(allocated(xle)) deallocate(xle)
                if(allocated(wle)) deallocate(wle)
                if(allocated(phia)) deallocate(phia)
                if(allocated(Cmat)) deallocate(Cmat)
                if(allocated(Vcouple)) deallocate(Vcouple)
                if(allocated(T)) deallocate(T)
                if(allocated(B_i)) deallocate(B_i)
                if(allocated(Ech)) deallocate(Ech)
                if(allocated(C)) deallocate(C)
                if(allocated(Rmat)) deallocate(Rmat)
                if(allocated(Z_O)) deallocate(Z_O)
                if(allocated(Z_I)) deallocate(Z_I)
                if(allocated(Smat)) deallocate(Smat)
                if(allocated(uij)) deallocate(uij)
                if(allocated(uijp)) deallocate(uijp)
                if(allocated(f))  deallocate(f)
                if(allocated(wf_int)) deallocate(wf_int)
                if(allocated(Gamma)) deallocate(Gamma)
ccccccc
                if(allocated(O)) deallocate(O)
ccccccc
                if(allocated(w1)) deallocate(w1)
                if(allocated(vl1)) deallocate(vl1)
                if(allocated(vr1)) deallocate(vr1)
                if(allocated(index)) deallocate(index)
ccccccc
                allocate(WTK(1:beta%nchmax+1),WTKP(1:beta%nchmax+1))
                allocate(xle(1:nr),wle(1:nr))
                allocate(phia(1:nr))
                allocate(Cmat(1:nr,1:nr,1:beta%nchmax,1:beta%nchmax))
                allocate(Vcouple(1:nr,1:nr,1:beta%nchmax,1:beta%nchmax))
                allocate(T(1:nr,1:nr,1:beta%nchmax))
                allocate(B_i(1:beta%nchmax))
                allocate(Ech(1:nr,1:nr,1:beta%nchmax))
                allocate(C(1:nr*beta%nchmax,1:nr*beta%nchmax))
                allocate(Rmat(1:beta%nchmax,1:beta%nchmax),Smat(1:beta%nchmax,1:beta%nchmax))
                allocate(Z_O(1:beta%nchmax,1:beta%nchmax),Z_I(1:beta%nchmax,1:beta%nchmax))
                allocate(uij(1:beta%nchmax,1:beta%nchmax))
                allocate(uijp(1:beta%nchmax,1:beta%nchmax))
                allocate(f(1:nr,1:beta%nchmax,1:beta%nchmax))
                allocate(wf_int(1:nr,1:beta%nchmax,1:beta%nchmax))
                allocate(Gamma(1:beta%nchmax))
ccccccc
                allocate(O(1:beta%nchmax,1:beta%nchmax))
ccccccc
                allocate(w1(1:nr*beta%nchmax))
                allocate(vl1(1:nr*beta%nchmax,1:nr*beta%nchmax))
                allocate(vr1(1:nr*beta%nchmax,1:nr*beta%nchmax))
                allocate(index(1:nr*beta%nchmax))
ccccccc
!lagrange legendre mesh
                call LEGZO(nr,xle,wle)
ccccccc
                do k=1,nr
                    phia(k)=(-1)**(nr+k)*sqrt(1d0/rmax/xle(k)/(1d0-xle(k)))
                end do
ccccccc              
                do i=1,beta%nchmax
                      if(E>Ec(i)) then
                        k_i=sqrt(2d0*mu*abs(Ec(i)-E)/hbarc**2)
ccccccc
                        li=int(lc(i),4)
ccccccc
                        allocate(FC_i(0:li),GC_i(0:li),FCP_i(0:li),GCP_i(0:li))
ccccccc
                         call COUL90(2*k_i*rmax,z_d*z_alpha*e2*mu/hbarc**2/k_i,0d0,li,FC_i,GC_i,FCP_i,GCP_i,KFN,IFAIL)
ccccccc H^{+}=G+iF, H^{-}=G-iF, only give H_i^{+} and H_i^{-}, i channel is enough                       
                         hlp_i=cmplx(GC_i(li),FC_i(li),kind=8)
                         dhlp_i=cmplx(GCP_i(li),FCP_i(li),kind=8)
                                B_i(i)=0!2*k_i*rmax*dhlp_i/hlp_i
                         deallocate(FC_i,GC_i,FCP_i,GCP_i)
!!B=2ki Rmax*H^{+}'/H^{+}
                        else
                        ki=sqrt(2d0*mu*abs(E-Ec(i))/hbarc**2)
                        eta=z_d*z_alpha*e2*mu/hbarc**2/ki
                        call WHIT(eta,rmax,ki,E,int(lc(i),4),WTK,WTKP,0)
                        B_i(i)=2*ki*rmax*WTKP(i)/WTK(i)
                      end if
                end do
ccccccc
        end subroutine rmat_int
ccccccc
!this subroutine gives the potential
        subroutine getpot(str)
            implicit none
            integer::i,j,k,info,lwork,s
            character(len=*)::str
ccccccc
            real*8::xx,zz,vcen,vtens,vls,rn,an,r
            complex*16::cvn
            real*8,allocatable::wr(:),wi(:),vr(:,:),work(:)
            real*8::x,y1,y2
ccccccc
        if(allocated(Vc)) deallocate(Vc)
        if(allocated(VijN)) deallocate(VijN)
        if(allocated(VijC)) deallocate(VijC)
        if(allocated(wr)) deallocate(wr)
        if(allocated(wi)) deallocate(wi)
        if(allocated(vr)) deallocate(vr)
        if(allocated(work)) deallocate(work)
ccccccc
        allocate(Vc(1:nr,1:beta%nchmax,1:beta%nchmax))
        allocate(VijN(1:nr,1:beta%nchmax,1:beta%nchmax))
        allocate(VijC(1:nr,1:beta%nchmax,1:beta%nchmax))   
        allocate(wr(1:beta%nchmax),wi(1:beta%nchmax))
        allocate(vr(1:beta%nchmax,1:beta%nchmax),work(1:3*beta%nchmax))
ccccccc
        select case(str)
ccccccc
!t: tensor force term included in the coupled pot for neutron-proton scattering

            case('t') 
! c Reid neutron-proton potential (T=1, soft core)(only for 2 channels l=0,2)
!         do i=1,nr
!             xx=0.7d0*xle(i)*rmax
!             zz=exp(-xx)
!             vcen=(-10.463d0*zz+105.468d0*zz**2-3187.8d0*zz**4+9924.3d0*zz**6)/xx
!             vtens=-10.463d0*((1+3/xx+3/xx**2)*zz-(12/xx+3/xx**2)*zz**4)/xx+351.77d0*zz**4/xx-1673.5d0*zz**6/xx
!             vls=708.91d0*zz**4/xx-2713.1d0*zz**6/xx
!             Vc(i,1,1)=vcen-2*(beta%j_tot-1)*vtens/(2*beta%j_tot+1)+(beta%j_tot-1)*vls
!             Vc(i,1,2)=6*vtens*sqrt(beta%j_tot*(beta%j_tot+1.0d0))/(2*beta%j_tot+1)
!             Vc(i,2,1)=Vc(i,1,2)
!             Vc(i,2,2)=vcen-2*(beta%j_tot+2)*vtens/(2*beta%j_tot+1)-(beta%j_tot+2)*vls
!         end do
ccccccc
            open(10,file='cpot.txt')
            do i=1,2
                do j=1,2
                    do k=1,nr
                        read(10,*) x,y1,y2
                        Vc(k,i,j)=cmplx(y1,y2)
                    end do
                end do
            end do    
            close(10)

ccccccc

!!coupled channel
!4th order deformation term and 16th order deformation term
!\sum_{l=2,4}\sqrt{\frac{(2l+1)(2I+1)}{4\pi(2J+1)}}\beta_{l}R_{d}\times[\bra{I,0,l,0}\ket{J,0}]^2
!even-even case: I=2(i-1),J=2(j-1)
        case('c')
            call factorialgen(1000) 
ccccccc
            do i=1,beta%nchmax
                do j=1,beta%nchmax
                    O(i,j)=sqrt((2d0*2+1d0)*(4d0*i-3)/(4d0*pi*(4d0*j-3)))*beta_2*R_d*cleb(2*(2*i-2),0,2*2,0,2*(2*j-2),0)**2
     &          +sqrt((2d0*4+1d0)*(4d0*i-3)/(4d0*pi*(4d0*j-3)))*beta_4*R_d*cleb(2*(2*i-2),0,4*2,0,2*(2*j-2),0)**2
                end do
            end do
ccccccc
            call dgeev('N','V',beta%nchmax,O,beta%nchmax,wr,wi,1,1,vr,beta%nchmax,work,-1,info)
ccccccc
            lwork=work(1)
            deallocate(work)
            allocate(work(lwork))
ccccccc
            call dgeev('N','V',beta%nchmax,O,beta%nchmax,wr,wi,1,1,vr,beta%nchmax,work,lwork,info)
ccccccc
!V_{ij}^{N}(r)
            do i=1,beta%nchmax
                do j=1,beta%nchmax
                    do k=1,nr
                        do s=1,beta%nchmax
                        VijN(k,i,j)=VijN(k,i,j) + vr(i,s) * vr(j,s) * WSpot(v_0,aa,r_0,xle(k)*rmax-wr(s))
                        end do
                    end do
                end do
            end do
ccccccc
!V_{ij}^{C}(r), separated into 2 parts: diagonal and off-diagonal of channel index ij
            do i=1,beta%nchmax
                do j=1,beta%nchmax
ccccccc part of diagonal matrix elements, non-deformed term only in diagonal part
                    if(i==j) then
                        do k=1,nr
                            if(xle(k)*rmax<=R_C) then 
                                VijC(k,i,j) = VijC(k,i,j) + z_alpha*z_d*e2/2d0/R_C*(3d0-(xle(k)*rmax)**2/R_C**2)
                            else
                                VijC(k,i,j) = VijC(k,i,j) + z_alpha*z_d*e2/xle(k)/rmax
                            end if
                        end do
                    end if
ccccccc
!for coulomb case , deformation l=2,4, and deformed coulomb part in both diagonal and off-diagonal
ccccccc
                    do k=1,nr
                        if(xle(k)*rmax<=R_C) then 
                    VijC(k,i,j)=VijC(k,i,j)+3d0*z_alpha*z_d*e2/(2d0*2+1)*(xle(k)*rmax)**2*R_C**(-3d0)
     &             *beta_2*sqrt((2d0*2+1)*(4d0*i-3)/(4d0*pi*(4d0*j-3)))*cleb(2*(2*i-2),0,2*2,0,2*(2*j-2),0)**2 
     &             +3d0*z_alpha*z_d*e2/(2d0*4+1)*(xle(k)*rmax)**4*R_C**(-5d0)
     &             *beta_4*sqrt((2d0*4+1)*(4d0*i-3)/(4d0*pi*(4d0*j-3)))*cleb(2*(2*i-2),0,4*2,0,2*(2*j-2),0)**2     
                    else
                    VijC(k,i,j)=VijC(k,i,j)+3d0*z_alpha*z_d*e2/(2d0*2+1)*R_C**2d0*(xle(k)*rmax)**(-3d0)
     &             *beta_2*sqrt((2d0*2+1)*(4d0*i-3)/(4d0*pi*(4d0*j-3)))*cleb(2*(2*i-2),0,2*2,0,2*(2*j-2),0)**2
     &             +3d0*z_alpha*z_d*e2/(2d0*4+1)*R_C**4d0*(xle(k)*rmax)**(-5d0)
     &             *beta_4*sqrt((2d0*4+1)*(4d0*i-3)/(4d0*pi*(4d0*j-3)))*cleb(2*(2*i-2),0,4*2,0,2*(2*j-2),0)**2               
                        end if 
                    end do
ccccccc
                end do 
            end do
ccccccc
!V_{ij}(r)=V_{ij}^{C}(r)+V_{ij}^{N}(r)
            do i=1,beta%nchmax
                do j=1,beta%nchmax
                    do k=1,nr
                        Vc(k,i,j)=VijC(k,i,j)+VijN(k,i,j)
                    end do
                end do
            end do
ccccccc          
        end select

       
        end subroutine


ccccccc
        subroutine rmatrix()
ccccccc
        implicit none
        integer::i,j,k,mm,nn,kk !sum variables
        integer::li,lj     !lc(i),lc(j)
ccccccc
        do mm=1,nr
            do i=1,beta%nchmax
                Ech(mm,mm,i)=Ec(i)-E
            end do
        end do        
ccccccc coupled potential matrix elements Vcouple_{im,jn}
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                do mm=1,nr
                    Vcouple(mm,mm,i,j)=Vc(mm,i,j)
                end do
            end do
        end do
ccccccc
!T+L(B), kinetic energy and Bloch term
        do mm=1,nr 
            do nn=1,nr
                do i=1,beta%nchmax
                    if(mm==nn) then 
                        T(mm,nn,i)=hbarc**2/2/mu/rmax*phia(nn)**2
     &         *(((4*nr**2+4*nr+3d0)*xle(nn)*(1-xle(nn))-
     &            6d0*xle(nn)+1d0)/3d0/xle(nn)/(1-xle(nn))-B_i(i))
                    else
                        T(mm,nn,i)=hbarc**2/2/mu/rmax
     &         *phia(mm)*phia(nn)*(nr**2+nr+1d0+(xle(nn)+xle(mm)-2d0*xle(mm)*xle(nn))/
     &           (xle(nn)-xle(mm))**2-1d0/(1-xle(nn))-1d0/(1-xle(mm))-B_i(i))        
                    end if
                end do
            end do
        end do
!then we add the several matrix elements together to get Cmatrix (without centrifugal term)
        do mm=1,nr
            do nn=1,nr
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                            if(i==j) then 
                            Cmat(mm,nn,i,j)=Cmat(mm,nn,i,j)+Ech(mm,nn,i)+T(mm,nn,i)+Vcouple(mm,nn,i,j)
                            else
                            Cmat(mm,nn,i,j)=Cmat(mm,nn,i,j)+Vcouple(mm,nn,i,j)
                            end if                 
                    end do
                end do
            end do
        end do     
!centrifugal term added to Cmatrix, now Cmatrix is complete
        do i=1,beta%nchmax
            do mm=1,nr
                Cmat(mm,mm,i,i)=Cmat(mm,mm,i,i)+hbarc**2/2d0/mu*lc(i)*(lc(i)+1d0)/xle(mm)**2/rmax**2
            end do
        end do
ccccccc
!reconstruct the Cmatrix for inversion, the size of C is (nr*beta%nchmax,nr*beta%nchmax)
ccccccc
        do mm=1,nr
            do nn=1,nr
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        C((i-1)*nr+mm,(j-1)*nr+nn)=Cmat(mm,nn,i,j)
                    end do
                end do
            end do
        end do    
ccccccc
!get the inversion of Cmatrix, and the inversion is stored just in C.
        call mat_inv(C,nr*beta%nchmax,nr*beta%nchmax)
!Rmatrix , R_{ij}=hbar^2/(2mu a)*\sum_{mn}φ_n(a)(C^{-1})_{in,jm}φ_m(a)
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                do mm=1,nr
                    do nn=1,nr
                        Rmat(i,j) =Rmat(i,j)+hbarc**2/2/mu/rmax*phia(mm)*C((i-1)*nr+mm,(j-1)*nr+nn)*phia(nn)        
                    end do
                end do
            end do
        end do
ccccccc
!Zmatrix: Z_O,Z_I, Smatrix: S=(Z_O)^{-1}Z_I
ccccccc
        KFN=0
        do i=1,beta%nchmax
            do j=1,beta%nchmax
                k_i=sqrt(2d0*mu*abs(Ec(i)-E)/hbarc**2)
                k_j=sqrt(2d0*mu*abs(Ec(j)-E)/hbarc**2)
ccccccc
                li=int(lc(i),4)
                lj=int(lc(j),4)
ccccccc
                allocate(FC_i(0:li),GC_i(0:li),FCP_i(0:li),GCP_i(0:li))
                allocate(FC_j(0:lj),GC_j(0:lj),FCP_j(0:lj),GCP_j(0:lj))
ccccccc
                call COUL90(k_i*rmax,z_d*z_alpha*e2*mu/hbarc**2/k_i,0d0,li,FC_i,GC_i,FCP_i,GCP_i,KFN,IFAIL)
                call COUL90(k_j*rmax,z_d*z_alpha*e2*mu/hbarc**2/k_j,0d0,lj,FC_j,GC_j,FCP_j,GCP_j,KFN,IFAIL)
ccccccc H^{+}=G+iF, H^{-}=G-iF, and their derivatives
                hlp_i=cmplx(GC_i(li),FC_i(li),kind=8)
                hln_i=cmplx(GC_i(li),-FC_i(li),kind=8)
                dhlp_i=cmplx(GCP_i(li),FCP_i(li),kind=8)
                dhln_i=cmplx(GCP_i(li),-FCP_i(li),kind=8)
                hlp_j=cmplx(GC_j(lj),FC_j(lj),kind=8)
                hln_j=cmplx(GC_j(lj),-FC_j(lj),kind=8)
                dhlp_j=cmplx(GCP_j(lj),FCP_j(lj),kind=8)
                dhln_j=cmplx(GCP_j(lj),-FCP_j(lj),kind=8)
ccccccc
                Z_O(i,j)=(k_j*rmax)**(-0.5d0)*(hlp_i*delta(i,j)-k_j*rmax*Rmat(i,j)*dhlp_j)
                Z_I(i,j)=(k_j*rmax)**(-0.5d0)*(hln_i*delta(i,j)-k_j*rmax*Rmat(i,j)*dhln_j)
ccccccc
                deallocate(FC_i,GC_i,FCP_i,GCP_i)
                deallocate(FC_j,GC_j,FCP_j,GCP_j)
            end do
        end do
ccccccc
! Smatrix: S=(Z_O)^{-1}Z_I,first we get the inverse of Z_O
                call mat_inv(Z_O,beta%nchmax,beta%nchmax)
                Smat=matmul(Z_O,Z_I)
ccccccc
!this loop give uij: to give channel wf value at the boundary  u_{ij}(R=rmax)
ccccccc
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        k_i=sqrt(2d0*mu*abs(Ec(i)-E)/hbarc**2)
ccccccc
                        li=int(lc(i),4)
ccccccc
                        allocate(FC_i(0:li),GC_i(0:li),FCP_i(0:li),GCP_i(0:li))
ccccccc
                        call COUL90(k_i*rmax,z_d*z_alpha*e2*mu/hbarc**2/k_i,0d0,li,FC_i,GC_i,FCP_i,GCP_i,KFN,IFAIL)
ccccccc H^{+}=G+iF, H^{-}=G-iF, only give H_i^{+} and H_i^{-}, i channel is enough                       
                        hlp_i=cmplx(GC_i(li),FC_i(li),kind=8)
                        dhlp_i=cmplx(GCP_i(li),FCP_i(li),kind=8)
                        hln_i=cmplx(GC_i(li),-FC_i(li),kind=8)
                        dhln_i=cmplx(GCP_i(li),-FCP_i(li),kind=8)
ccccccc
!-SH^{+},no hln_i*delta(i,j) term
                    uij(i,j)=(hln_i*delta(i,j)-Smat(i,j)*hlp_i)*(0,0.5d0)
                    uijp(i,j)=(dhln_i*delta(i,j)*k_i-Smat(i,j)*dhlp_i*k_i)*(0,0.5d0)
ccccccc
                        deallocate(FC_i,GC_i,FCP_i,GCP_i)
                    end do
                end do
ccccccc
!expansion coefficients for wf_int
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        do nn=1,nr
                            f(nn,i,j)=0d0
                            do kk=1,beta%nchmax
                                do mm=1,nr
                                    f(nn,i,j)=f(nn,i,j)+C((i-1)*nr+nn,(kk-1)*nr+mm)
     &                              *phia(mm)*hbarc**2/2d0/mu*(uijp(j,kk)-B_i(kk)/rmax*uij(j,kk))
                                end do
                            end do
                        end do
                    end do
                end do
ccccccc
!internal wf,wf_int
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        do nn=1,nr
                           wf_int(nn,i,j)=(rmax*wle(nn))**(-0.5d0)*f(nn,i,j)
                        end do
                    end do
                end do
ccccccc
         do i=1,beta%nchmax
            do j=1,beta%nchmax
                do nn=1,nr
                    write(333,*) xle(nn)*rmax,aimag(wf_int(nn,i,j))
                    write(444,*) xle(nn)*rmax,real(wf_int(nn,i,j))
                    write(555,*) xle(nn)*rmax,real(Vc(nn,i,j))
                end do
                write(444,*) '&'
                write(555,*) '&'
            end do 
        end do

ccccccc
! give the (partial) widths , S(1,i) contributes to the ith channel
! here Gamma has the unit of MeV
                do i=1,beta%nchmax
                    li=int(lc(i),4)
                    k_i=sqrt(2d0*mu*abs(Ec(i)-E)/hbarc**2)
ccccccc
                    allocate(FC_i(0:li),GC_i(0:li),FCP_i(0:li),GCP_i(0:li))
                    call COUL90(k_i*rmax,z_d*z_alpha*e2*mu/hbarc**2/k_i,0d0,li,FC_i,GC_i,FCP_i,GCP_i,KFN,IFAIL)
ccccccc
                    Gamma(i)=hbarc**2d0*k_i/mu * abs(uij(1,i))**2/(GC_i(li)**2+FC_i(li)**2)
ccccccc
                    deallocate(FC_i,GC_i,FCP_i,GCP_i)
                end do 

ccccccc
199     format(25('-'),'S matrix',25('-'))
ccccccc
200     format('Channel',13X, 'Re(S)', 13X, 'Im(S)')
201     format('S(',I2,'->',I2,')', 5X, F15.10, 5X, F15.10)
ccccccc
202     format(25('-'),'Modulus of S matrix',25('-'))
ccccccc
203     format('Channel',13X, 'Abs(S)')
204     format('S(', I2, '->', I2, ')',5X,F15.10)
ccccccc
210     format('Channel',13X,'Decay Width')
211     format(I2,13X,F15.10,1X,'MeV')
ccccccc
212     format(30('-'))

ccccccc
                write(*,199)
ccccccc
                write(*,200)
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        write(*,201) ,i,j,real(Smat(i,j)),aimag(Smat(i,j))
                    end do
                end do
ccccccc
                write(*,202)
ccccccc
                write(*,203)
                do i=1,beta%nchmax
                    do j=1,beta%nchmax
                        write(*,204) i, j, abs(Smat(i,j))
                    end do
                end do
ccccccc
                write(*,212)
ccccccc


ccccccc
300     format('Generated files:')
301     format('fort.1: input parameters')
302     format('fort.333: Im(u) of each channel')
303     format('fort.444: Re(u) of each channel')
304     format('fort.555: Re(V) of each channel')
305     format(50('-'))
ccccccc
        write(*,300)
        write(*,301)
        write(*,302)
        write(*,303)
        write(*,304)
        write(*,305)
ccccccc
        end subroutine
ccccccc

        end module