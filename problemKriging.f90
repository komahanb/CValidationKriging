program CrossValidationKriging
  use dimkrig,only:id_proc

  implicit none  
  include 'mpif.h'

  integer,parameter::ndim=4
  integer::fct,initpts,ncyc,npts,stat,nseed,ierr
  real*8::DAT,cverrout
  !  real*8::sample(ndim,500)  
  real*8::Tsample(ndim,500)
  real*8::testpoint(ndim)
  real*8::CVE(500)
  real*8::MCVE,MAXCVE
  integer::j,k,i,knot,nidx,nfunc
  double precision, allocatable, dimension(:,:) :: sample
  integer::nmodel 

  call MPI_START


  if (id_proc.eq.0) then
     print*,''
     print*,'=========================='
     print*,' CROSS VALIDATION KRIGING '
     print*,'=========================='
     print*,''
  end if


  DAT=77 !screen or file

  stat=0

  do nfunc=2,2

     if (nfunc.eq.1) fct=0
     if (nfunc.eq.2) fct=2
     if (nfunc.eq.3) fct=3
     if (id_proc.eq.0) then
        print*,''
        print*,'=========================='
        print*,'Function number ',fct
        print*,'=========================='
        print*,''
     end if

     !0=exp,2=runge,3=rosen

     if (id_proc.eq.0) then
!!$
!!$        if (fct.eq.0) open(unit=37,file='CVKrigfct00dim2.his',form='formatted',status='replace')
!!$        if (fct.eq.2) open(unit=37,file='CVKrigfct02dim2.his',form='formatted',status='replace')
!!$        if (fct.eq.3) open(unit=37,file='CVKrigfct03dim2.his',form='formatted',status='replace')
     end if

 !    if (id_proc.eq.0) write(37,'(a)') 'npts MeanCVE MaxCVE Nmodels'
  !   nmodel=0
   
     do npts=300,300

        allocate(sample(ndim,npts))
!        allocate(samplewithIdx(ndim,npts))
        !        print*,"allocate"

        ncyc=npts
        initpts=npts

        
        
        if (id_proc.eq.0) then

           !Reading the data file

           open(unit=54,file='Training.csv',form='formatted')
           do i=1,300
              read(54,*)sample(1,i),sample(2,i),sample(3,i),sample(4,i) 
          end do
           close(54)


           !           call get_seed(nseed)
           !           call latin_random(ndim,npts-1,nseed,sample(:,2:npts))       

        end if

        call MPI_BCAST(sample,npts*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !        print*,"bcast"

        ! We now have a set of npts samples (first one is the center of domain)

        CVE(:)=0.0
        MaxCVE=0.0

        !!do k=1,npts !loop over all points to find the CVE and finally MCVE
k=1
           testpoint(:)=sample(:,k) !k-th point is the test point
           knot=k

           nidx=0

   !        do i=1,npts

              ! Check for all points and see if they dont match the test point and store it in the training vector to be passed to kriging

        !      if (i.ne.knot) then
       !          nidx=nidx+1
      !           Tsample(:,nidx)=sample(:,i)
     !         end if

    !       end do

           ! Now we have the test and training points

           !  call Krigingestimate(ndimin,ndimint,fctin,DATIN,initpts,ncyc,nptsin,statin,trainingdataptsin,testptin,fmeanout)
           !        print*,"before kriging routine"


           !       do i=1,npts-1
           !          print*,i,Tsample(:,i)
           !       end do

           !      print*,npts-1

           cverrout=0.0
!           nmodel=nmodel+1
           call Krigingestimate(ndim,ndim,fct,DAT,initpts,ncyc,npts,stat,Tsample(:,1:npts),testpoint,cverrout)
           CVE(k)=cverrout

           !stop
           !        print*,"after k"
           !           print*,'CVE',k,cve(k)

!           if (CVE(K).gt.MaxCVE) MaxCVE=CVE(k)

!        end do

!        MCVE=0.0
!        do j=1,npts
!           MCVE=MCVE+CVE(j)
!        end do
!        MCVE=MCVE/dble(npts)


!       if (id_proc.eq.0) write(37,'(i8,2e15.5,i8)'), npts, MCVE,MaxCVE,nmodel
!        if (id_proc.eq.0) print*,'NPTS:',npts,'MCVE:',MCVE,MaxCVE,nmodel


        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        deallocate(Sample)
        !        print*,"deallocate"
     end do  ! loop npts

     if (id_proc.eq.0) close(37)

  end do  !  loop function

  call stop_all
end program CrossValidationKriging

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine find_st(pst,dx,st)
  implicit none
  ! find st which ensures the probability within [-dx:dx] is pst
  double precision, intent(in)  :: pst,dx
  double precision, intent(out) :: st
  integer :: iter
  double precision :: pout,s,ds,vv,vp,dv,sini

  if(pst.le.0.d0.or.pst.ge.1.d0)stop'pst in find_st'
  pout = (1.d0-pst)/2.d0
  sini = 1.d0
190 continue
  s    = sini
  ds   = sini*1.d-3
  iter = 0
200 continue
  iter = iter + 1
  call CDF(-1.d0*dx,0.d0,s,   vv)
  if(dabs(vv-pout).le.1.d-10)go to 210
  call CDF(-1.d0*dx,0.d0,s+ds,vp)
  dv = (vp-vv)/ds
  !       write(*,'(5e15.5)')s,vv,pout,vv-pout,dv
  if(dv.eq.0.d0)stop'dv = 0.d0 in find_st'
  if(iter.ge.100)stop'iter>100 in find_st'
  s = s - (vv-pout)/dv
  if(s.le.0.d0)then
     sini = sini * 0.1d0
     go to 190
  end if
  go to 200
210 continue
  !       write(*,'(4e15.5)')s,vv,pout,vv-pout
  st = s

end subroutine find_st
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine find_x(mode,ran,xc,st,xout)
  implicit none
  ! find xout which ensures the CDF at xout is ran
  ! mode=0 : analytical CDF (fast but less robust?)
  ! mode=1 : numerical  CDF (time comsuming, but robust)
  integer, intent(in) :: mode
  double precision, intent(in)  :: ran,xc,st
  double precision, intent(out) :: xout
  integer :: iter 
  double precision :: x,vv,dv,vp,dx

  x    = xc
  dx   = 1.d-4
  iter = 0
200 continue
  iter = iter + 1
  if(mode.eq.0)then
     call CDF(x,xc,st,vv)
  else
     call CDF_Numerical(x,xc,st,vv)
  end if
  if(dabs(vv-ran).le.1.d-10)go to 210
  if(mode.eq.0)then
     call CDF(x+dx,xc,st,vp)
  else
     call CDF_Numerical(x+dx,xc,st,vp)
  end if
  dv = (vp-vv)/dx
  !       write(*,'(4e15.5)')x,vv,vv-ran,dv
  if(dv.eq.0.d0)stop'dv=0 in find_x'
  if(iter.ge.100)stop'iter>100 in find_x'
  x = x - (vv-ran)/dv
  go to 200
210 continue
  !       write(*,'(3e15.5)')x,vv,vv-ran
  xout = x

end subroutine find_x
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine CDF(xin,xc,st,vout)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: vout
  double precision :: vtmp
  !       vout = 0.5d0 * (1.d0 + erf( (xin-xc)/(st*dsqrt(2.d0)) ))
  call ERF_MINE1( (xin-xc)/(st*dsqrt(2.d0)), vtmp )
  vout = 0.5d0 * (1.d0 + vtmp)
end subroutine CDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine DCDF(xin,xc,st,dvout)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: dvout
  double precision :: dvtmp
  call DERF_MINE( (xin-xc)/(st*dsqrt(2.d0)), dvtmp )
  dvout = 0.5d0*dvtmp/(st*dsqrt(2.d0))
end subroutine DCDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ERF_MINE1(xin,yout)
  implicit none
  double precision, intent(in)  :: xin
  double precision, intent(out) :: yout
  integer :: i,k,n
  double precision :: vsum,kai
  ! n is the order of Taylor
  ! Maybe accurate within the range of [-4:4] with n=100
  n = 100
  vsum = 0.d0
  do i=0,n
     kai = 1.d0
     do k=1,i
        kai = kai * (-1.d0) * xin**2 / dble(k)
     end do
     vsum = vsum + kai*xin/dble(2*i+1)
  end do
  yout = vsum*2.d0/(dsqrt(3.141592653589793238d0))

  if(yout.gt.1.d0)write(*,'(a,2e15.5)')'*ERF>1 ',xin,yout-1.d0
end subroutine ERF_MINE1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine DERF_MINE(xin,dyout)
  implicit none
  double precision, intent(in)  :: xin
  double precision, intent(out) :: dyout
  double precision :: vsum

  vsum  = exp(-1.d0*xin**2)
  dyout = vsum*2.d0/(dsqrt(3.141592653589793238d0))

end subroutine DERF_MINE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine CDF_Numerical(xin,xc,st,cdf)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: cdf
  double precision :: vtmp
  integer :: i,num
  double precision :: xs,xe,dx,x1,x2,pdf1,pdf2
  if(xin.lt.xc)then
     cdf = 0.d0
     xs  = xin -2.d0
     xe  = xin
  else if(xin.ge.xc)then
     cdf = 0.5d0
     xs  = xc
     xe  = xin
  end if
  num = 1001
  dx  = (xe-xs)/dble(num-1)
  do i=1,num-1
     x1 = xs + dble(i-1)*dx
     x2 = xs + dble(i  )*dx
     call normal_dist(x1,xc,st,pdf1)
     call normal_dist(x2,xc,st,pdf2)
     cdf = cdf + (pdf1+pdf2)*dx*0.5d0
  end do
end subroutine CDF_Numerical
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine normal_dist(xin,xc,st,y)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: y
  double precision :: pi

  pi = 4.d0*datan(1.d0)
  y = exp(-1.d0*(xin-xc)**2/2.d0/st/st)/dsqrt(2.d0*pi*st**2)

end subroutine normal_dist

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function dinvnorm(p)

  implicit none

  real*8 :: dinvnorm,p,p_low,p_high
  real*8 :: a1,a2,a3,a4,a5,a6
  real*8 :: b1,b2,b3,b4,b5
  real*8 :: c1,c2,c3,c4,c5,c6
  real*8 :: d1,d2,d3,d4
  real*8 :: z,q,r
  a1=-39.6968302866538
  a2=220.946098424521
  a3=-275.928510446969
  a4=138.357751867269
  a5=-30.6647980661472
  a6=2.50662827745924
  b1=-54.4760987982241
  b2=161.585836858041
  b3=-155.698979859887
  b4=66.8013118877197
  b5=-13.2806815528857
  c1=-0.00778489400243029
  c2=-0.322396458041136
  c3=-2.40075827716184
  c4=-2.54973253934373
  c5=4.37466414146497
  c6=2.93816398269878
  d1=0.00778469570904146
  d2=0.32246712907004
  d3=2.445134137143
  d4=3.75440866190742
  p_low=0.02425
  p_high=1-p_low
  if(p.lt.p_low) goto 201
  if(p.ge.p_low) goto 301
201 q=dsqrt(-2*dlog(p))
  z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
  goto 204
301 if((p.ge.p_low).and.(p.le.p_high)) goto 202
  if(p.gt.p_high) goto 302
202 q=p-0.5
  r=q*q
  z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
  goto 204
302 if((p.gt.p_high).and.(p.lt.1)) goto 203
203 q=dsqrt(-2*dlog(1-p))
  z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204 dinvnorm=z

  return

end function dinvnorm
