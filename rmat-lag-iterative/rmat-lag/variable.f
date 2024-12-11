        module rmatvar
        implicit none
ccccccc
        complex*16,allocatable::Cmat(:,:,:,:)
        complex*16,allocatable::VijN(:,:,:) !nuclear pot part
        complex*16,allocatable::VijC(:,:,:) !coulomb pot part
        complex*16,allocatable::Vc(:,:,:)
        complex*16,allocatable::Vcouple(:,:,:,:)
        complex*16,allocatable::T(:,:,:)
        complex*16,allocatable::B_i(:)
        complex*16,allocatable::Ech(:,:,:)
        complex*16,allocatable::C(:,:)
        complex*16,allocatable::Rmat(:,:)
        complex*16,allocatable::Z_O(:,:)
        complex*16,allocatable::Z_I(:,:)
        complex*16,allocatable::Smat(:,:)
        complex*16,allocatable::u_exta(:,:)    !size: (nch,nch), u_{ij}^{ext}(a)
        complex*16,allocatable::u_extap(:,:)   !size: (nch,nch), u_{ij}^{ext}'(a)
        complex*16,allocatable::f(:,:,:)       !size: (nr,nch,nch), f_{i,n}^{j}=>f_{ij,n}
        complex*16,allocatable::wf_int(:,:,:)  !size: (nr,nch,nch), internal wf u_{ij} (r)
        complex*16,allocatable::vl1(:,:),vr1(:,:) !size: (nch*nr,nch*nr), left and right eigenvectors
        complex*16,allocatable::w1(:) !size: (nch*nr), eigenvalues
ccccccc
        integer,allocatable::index(:)
        end module      
ccccccc
        module mesh 
            implicit none
            real*8::rmax
            integer::nr
            real*8,allocatable::xle(:),wle(:)
            real*8,allocatable::phia(:)
        end module
ccccccc
!2b sysmtem variables
        module system
            implicit none
            real*8::mass_d      !daughter nucleus mass
            real*8::mass_alpha   !alpha mass
            real*8::mu
            real*8::z_d
            real*8::z_alpha
            real*8::z12
            real*8::I_d         !degree of asymmetry of daughter nucleus
            real*8::R_d         !radius of daughter nucleus 
            real*8::R_C         !charge radius of d
        end module
ccccccc
        module deform
            implicit none
            real*8,allocatable::O(:,:)          !deformation operator matrix element under channel states
            real*8::beta_2 !4th order deformation parameter
            real*8::beta_4 !16th order deformation parameter
        end module
ccccccc
!COUL90's variables
        module coulvar
        implicit none
ccccccc
        integer::KFN,IFAIL      !coul90 variables
        real*8::k_i,k_j !channel wave numbers
        real*8,allocatable::FC_i(:),GC_i(:),FCP_i(:),GCP_i(:) !coul90 variables for channel i
        real*8,allocatable::FC_j(:),GC_j(:),FCP_j(:),GCP_j(:) !coul90 variables for channel j
        complex*16::hlp_i,hln_i,dhln_i,dhlp_i   !hankel H^{\pm} and its derivative for channel i
        complex*16::hlp_j,hln_j,dhln_j,dhlp_j   !hankel H^{\pm} and its derivative for channel j
ccccccc
        end module
ccccccc
!whittaker function's variables
!ki: channel wave numbers k_i, eta=Z1Z2e^2mu/hbar^2/k_i,Sommerfeld parameter
!WTK,WTKP: Whittaker functions and their derivatives
        module wtkvar
        implicit none
        real*8::ki,eta
        real*8,allocatable::WTK(:),WTKP(:)
        end module
ccccccc
!potvar: variables in the potential function
!str: type of potential, char type var
        module potvar
        implicit none
        character(len=20)::str
        !coupling potential parameter
        real*8::aa      !decay rate of the WS potential
        real*8::r_0,v_0    !WSpot parameter
        end module 
ccccccc
        module parameter
            implicit none
            real*8,parameter :: hbarc=1d0!197.3269718d0  !hbar      ! NIST Ref 02.12.2014   ! MeV.fm           
            real*8,parameter :: finec=137.03599d0
            real*8,parameter :: amu=931.49432d0      !MeV
            real*8,parameter :: e2=1.43997d0         !MeV.fm
            real*8,parameter :: PI=acos(-1.0)
            complex*16,parameter :: ii=(0.0d0,1.0d0)
        end module
