		module clebsg
			real*8:: dlfac(0:10000)
			real*8:: dl2fac(0:10000)
			real*8:: FACT(0:10000),FFAKINV(0:10000),WFAK(0:10000)
			real*8 :: WFAKINV(0:10000)
		contains
c-----------------------------------------------------------------------
      subroutine factorialgen(n)           !!!!call factorialgen(2*lmax)
! n>0
! ln(i!) and ln((i)!!) for i is odd number 
      implicit real*8 (a-h,o-z)
      dlfac(0)=0.
      fact(0)=1.
      FFAKINV(0)=1.0d0
      WFAK(0)=1.0d0
      WFAKINV(0)=1.0d0


      do 1 i=1,n
      a=i
1      dlfac(i)=dlfac(i-1)+log(a)
C      fact(i)=I*fact(I-1)
C      FFAKINV(I)=1.0/fact(I)
C      WFAK(I)=SQRT(fact(I))
C1     WFAKINV(I)=1.0/WFAK(I)

      dl2fac(1)=0.d0
      do 2 i=3,n+20,2
      a=i
2     dl2fac(i)=dl2fac(i-2)+log(a)


      dl2fac(0)=0.d0
      do 3 i=2,n+20,2
      a=i
3     dl2fac(i)=dl2fac(i-2)+log(a)
      continue
      end subroutine
c---------------------------------------------

c---------------------------------------------
c******
      function flog(i)
        implicit none
        integer::i
        real*8:: flog
!        real*8 :: logfac
!        flog=logfac(i-1)
        if(i==0) then
           flog=0.0d0
        else
           flog=dlfac(i-1)
        end if
      end function flog

c *** ---------------------------------------------
c Factorial LOG
c *** --------------------------------------------
      real*8 function logfac(n) ! FL(N)
      implicit real*8(a-h,o-z),integer*4(i-n)
       fl=0
       if(n>1) then
       FN = 1.
       DO 10 I = 2,N
       FN = FN + 1.
   10  FL = FL +  LOG(FN)
      endif
      logfac=fl
      END FUNCTION

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Clebsch-Gordan coefficient  <l',m',l,m|JM>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function cleb(ia,id,ib,ie,ic,if)
      implicit real*8(a-h,o-z)
C      real*8:: ria,rid,rib,rie,ric,rif
!      common/clebma/faclog(500)
!      COMMON /PRAHA/ FLOG(100), GM(100), DG(25)
C      ia=2.0d0*(ria+.0001d0)
C      ib=2.0d0*(rib+.0001d0)
C      ic=2.0d0*(ric+.0001d0)
C      id=int(sign(1.0_dpreal,rid)*2.0_dpreal*(abs(rid)+.0001d0))
C      ie=int(sign(1.0_dpreal,rie)*2.0_dpreal*(abs(rie)+.0001d0))
C      if=int(sign(1.0_dpreal,rif)*2.0_dpreal*(abs(rif)+.0001d0))
      wwww=-1.0d0
      cleb=0.0d0
      if(id+ie-if) 7000,105,7000
  105 k1=ia+ib+ic
      if((-1)**k1) 7000,107,107
  107 if(.not.((id.eq.0).and.(ie.eq.0))) go to 110
      k1=k1/2
      if((-1)**k1) 7000,110,110
  110 k1=ia+ib-ic
      k2=ic-iabs(ia-ib)
      k3=min0(k1,k2)
      if(k3) 7000,130,130
  130 if((-1)**(ib+ie)) 7000,7000,140
  140 if((-1)**(ic+if)) 7000,7000,150
  150 if(ia-iabs (id)) 7000,152,152
  152 if(ib-iabs (ie)) 7000,154,154
  154 if(ic-iabs (if)) 7000,160,160
  160 if(ia) 7000,175,165
  165 if(ib) 7000,175,170
  170 if(ic) 7000,180,250
  175 cleb=1.0d0
      go to 7000
  180 fb=float(ib+1)
      cleb=((wwww)**((ia-id)/2))/sqrt(fb)
      go to 7000
  250 fc2=ic+1
      iabcp=(ia+ib+ic)/2+1
      iabc=iabcp-ic
      icab=iabcp-ib
      ibca=iabcp-ia
      iapd=(ia+id)/2+1
      iamd=iapd-id
      ibpe=(ib+ie)/2+1
      ibme=ibpe-ie
      icpf=(ic+if)/2+1
      icmf=icpf-if
      vvv=0.5d0
      sqfclg=vvv*(log(fc2)-flog(iabcp+1)
     1      +flog(iabc)+flog(icab)+flog(ibca)
     2      +flog(iapd)+flog(iamd)+flog(ibpe)
     3      +flog(ibme)+flog(icpf)+flog(icmf))
      nzmic2=(ib-ic-id)/2
      nzmic3=(ia-ic+ie)/2
      nzmi= max0(0,nzmic2,nzmic3)+1
      nzmx= min0(iabc,iamd,ibpe)
      if(nzmx.lt.nzmi) go to 7000
      s1=(wwww)**(nzmi-1)
      do 400 nz=nzmi,nzmx
      nzm1=nz-1
      nzt1=iabc-nzm1
      nzt2=iamd-nzm1
      nzt3=ibpe-nzm1
      nzt4=nz-nzmic2
      nzt5=nz-nzmic3
      termlg=sqfclg-flog(nz)-flog(nzt1)-flog(nzt2)
     1           -flog(nzt3)-flog(nzt4)-flog(nzt5)
      ssterm=s1*exp (termlg)
      cleb=cleb+ssterm
  400 s1=-s1
 7000 return
      end function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		end module clebsg