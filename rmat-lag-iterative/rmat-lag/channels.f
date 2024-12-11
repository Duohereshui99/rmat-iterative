ccccccc
        module channels
        implicit none
!channel index for decay type(alpha), including        
!daughter nucleus spin jd, daughter nucleus excitation Ed,
!orbit angular momenta L, total(d+alpha) parity Pi_tot
!total j(parent nucleus j_tot)
!minimum and maximum orbit angular momenta lmin,lmax
!Pi_d is the parity of the daughter nucleus
!where lmin=|j-s|,lmax=j+s
!channel number variable nch
!number of channels nchmax  
!allocatable Sc stores the spin of daughter nucleus of different channels 
!jc stores the total angular momentum of the 2b system J 
!lc stores the orbit L between the 2 bodies
!with size being the number of the channel
!daughter excitation Ec
        type channel
        real*8::jd
        real*8::Ed
        real*8::Pi_d
        real*8::j_tot
        real*8::Pi_tot
        real*8::Pi_alpha
        real*8::j_alpha
        integer::lmin,lmax
        integer::nch
        integer::nchmax
        end type
ccccccc
        real*8,allocatable::Sc(:)
        real*8,allocatable::jc(:)
        real*8,allocatable::lc(:)
        real*8,allocatable::Ec(:)
!channel notation β=|l jd j>,2b system only requires one index β (1 configuration)
        type(channel)::beta 
!total energy E
        complex*16::E,E_0
ccccccc
        contains
!this subroutine calculates the number of channels for the 2b decay of alpha type
!and gives each channel index beta, then stores them in alpha type allocatable vector
!where J,jd are already known
!output: nchmax,S,j,l,E
        subroutine getchannelalphaD()
        implicit none
        integer::l
        real*8::tol
ccccccc
        tol=1e-6
ccccccc         
        beta%nchmax=0
ccccccc       
!         beta%lmin=min(nint(abs(beta%j_tot-abs(beta%jd-beta%j_alpha))),nint(beta%j_tot+abs(beta%jd-beta%j_alpha))
!      &   ,nint(abs(beta%j_tot-(beta%jd+beta%j_alpha))),nint(beta%j_tot+beta%jd+beta%j_alpha))
!         beta%lmax=max(nint(abs(beta%j_tot-abs(beta%jd-beta%j_alpha))),nint(beta%j_tot+abs(beta%jd-beta%j_alpha))
!      &   ,nint(abs(beta%j_tot-(beta%jd+beta%j_alpha))),nint(beta%j_tot+beta%jd+beta%j_alpha))
        beta%lmin=0d0
        beta%lmax=6d0
ccccccc
        do l=beta%lmin,beta%lmax
             if(abs((-1d0)**l*beta%Pi_alpha*beta%Pi_d-beta%Pi_tot)<tol) then 
                beta%nchmax=beta%nchmax+1 
             end if
        end do
!allocate different channel index vector
!with size being the number of channels nchmax
        if (allocated(Sc)) deallocate(Sc)
        if (allocated(jc)) deallocate(jc)
        if (allocated(lc)) deallocate(lc)
        if (allocated(Ec)) deallocate(Ec)
ccccccc
        allocate(Sc(1:beta%nchmax))
        allocate(jc(1:beta%nchmax))
        allocate(lc(1:beta%nchmax))
        allocate(Ec(1:beta%nchmax))
ccccccc
        beta%nch=1
        do l=beta%lmin,beta%lmax
             if(abs((-1d0)**l*beta%Pi_alpha*beta%Pi_d-beta%Pi_tot)<tol) then 
                Sc(beta%nch)=beta%jd
                jc(beta%nch)=beta%j_tot
                lc(beta%nch)=l
                Ec(beta%nch)=beta%Ed
ccccccc
                beta%nch=beta%nch+1
             end if
        end do
        Ec(1)=0; Ec(2)=0.005d0; Ec(3)=0.126d0; Ec(4)=0.266d0
ccccccc
        end subroutine getchannelalphaD


        end module



