ccccccc
      module input
        use rmatmod
        implicit none
ccccccc
        contains
ccccccc
        subroutine readinput()
ccccccc
            namelist /channelbeta/ beta
            namelist /meshs/ nr,rmax
            namelist /systems/  mass_d,mass_alpha,z_d,z_alpha
            namelist /potentials/ str
            namelist /deformation/ beta_2,beta_4,v_0
ccccccc
            read(5,nml=channelbeta)
            read(5,nml=meshs)
            read(5,nml=systems)
            read(5,nml=potentials)
            read(5,nml=deformation)
ccccccc
            write(1,nml=channelbeta)
            write(1,nml=meshs)
            write(1,nml=systems)
            write(1,nml=potentials)
ccccccc
            E=(1d0,1d0)
ccccccc
            z12=z_d*z_alpha
            mu=amu*(mass_d*mass_alpha/(mass_d+mass_alpha))
            write(*,*) 'mu=',mu
ccccccc
!m_d=A_d,N_d=A_d-Z_d,N_d-Z_d=m_d-Z_d-Z_d=m_d-2Z_d
            I_d=(mass_d-z_d-z_d)/mass_d
ccccccc
            R_d=(1d0+0.39*I_d)*mass_d**(1.d0/3.d0)
            r_0=R_d+1.17d0
ccccccc
        write(*,*) 'R_d=',R_d
        write(*,*) 'r_0=',r_0
            aa=0.5d0+0.33d0*I_d
            write(*,*) 'aa=',aa
ccccccc
            R_C=1.2d0*(mass_d**(1.d0/3.d0)+mass_alpha**(1.d0/3.d0))
            write(*,*) 'hbarc**2/amu=',hbarc**2/amu/2
        end subroutine
ccccccc
!-----------------------------------------------------------------------
ccccccc
        subroutine get_info()
#ifdef BASE
        print *, 'Base directory: ', BASE
#endif

#ifdef VERDATE
        print *, 'Version date: ', VERDATE
#endif

#ifdef VERREV
        print *, 'Version revision: ', VERREV
#endif

#ifdef COMPDATE
        print *, 'Compilation date: ', COMPDATE
#endif
            end subroutine  
ccccccc
      end module