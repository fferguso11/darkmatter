!     **********************************************************************************

!     A simulation of WIMP direct detection for the LUX experiment specification

!     A target dataset is created first, before samples are created to populate a 100x100
!     point phase space in logarithm WIMP mass and WIMP-nucleon cross section
!     50 samples are generated for each point, and a K.S test used for each to quantify
!     the probability that the samples generated arise from the same underlying distribution
      
!     **********************************************************************************      

!     Preliminary Functions

!     Function to generate random variable r using system clock to randomise
!     Can vary the factor value to optimise the seed distribution - may not be optimal currently
      function randtest() result(r)

      implicit none
      
      integer size,count,count_rate,count_max,factor,i
      integer, allocatable :: seed(:)
      real ::  r


      call random_seed(size = size)
      allocate(seed(size))

!     This factor helps to increase spread of seeds
      factor = size*1000

!     Call system clock to select random number
      call system_clock(count,count_rate,count_max)

      do i = 1,size
         seed(i) = mod(count,factor)
      end do

      call random_seed(put=seed)

!     r is our random number
      call random_number(r)

      end function randtest

      
!     ***********************************************************
      
!     Function to sort an array in ascending order

!     ***********************************************************

      subroutine sort(n, fa)

      implicit none

      integer :: n, i, j
      real :: fa(n), x
 
      do i = 2, n
        x = fa(i)
        j = i - 1
        do while (j >= 1)
            if (fa(j) <= x) exit
            fa(j + 1) = fa(j)
            j = j - 1
        end do
        fa(j + 1) = x
      end do
      
      end subroutine

!     ***********************************************************
      
!     Main program for the simulation
!     Have avoided using functions/subroutines where possible as
!     found that many compilers have problems with values being
!     corrupted when being passed between them

!     ***********************************************************      

      program simcode

      implicit none

      integer nsample,nevents,i,j,bins,nplot
      double precision pi,cross_section
      real Mw,earth_mass,sigma,sunv,time,tp,b
      real yr,Mx,mass_ratio,c,sconversion,a,true_cal
      real func(2,-1:1000),sfunc(2,-1:1000)
      double precision yfunc(2,-1:3001),nfunc(2,-1:3001),yreference
      real mlist(188),crlist(188),ycalib(188),yweight(188)
      integer yplot,ystep, ival,step,sstep
      real yield,cross,norm,r,randtest
      real numerator,denominator,exponential
      real numer,denomin,expon,max_value,smax_value
      real xcum(-3000:3000),ycum(-3000:3000),zcum(-3000:3000)
      real sxcum(-3000:3000),sycum(-3000:3000),szcum(-3000:3000)
      real speed_hist(0:50000),magvel_hist(0:50000),theta_hist(-180:180)
      real phi_hist(-180:180),ctheta_hist(-180:180)
      real sspeed_hist(0:50000),smagvel_hist(0:50000)
      real stheta_hist(-180:180),senergy_hist(0:50000)
      real sphi_hist(-180:180),sctheta_hist(-180:180)
      real scatter_hist(0:50000),recoil_hist(0:50000),energy_hist(0:50000)
      real sscatter_hist(0:50000),ssrecoil_hist(0:50000)
      real tan,cmvel,length,scmvel,sfactor,slength,stan
      real,dimension (:),allocatable:: sp,sspeed,smagvel
      real,dimension (:,:),allocatable :: sthreed,svel
      real,dimension (:,:),allocatable :: stwod
      real,dimension (:),allocatable :: stheta_cm
      real,dimension (:),allocatable :: stheta,sphi
      real,dimension (:),allocatable :: sscatter_vel
      real,dimension (:),allocatable :: srecoil_vel
      real,dimension (:),allocatable :: smtransfer,su
      real,dimension (:),allocatable :: snew_array,sscat_energy
      real,dimension (:),allocatable :: srec_energy,sorted1
      real u(10000),new_array(10000),rec_energy(10000),scat_energy(10000)
      real p(10000),speed(10000),magvel(10000),theta_cm(10000),theta(10000)
      real threed(3,10000),vel(3,10000),twod(2,10000),mtransfer(10000),sorted2(10000)
      real a2,probks,fac,D,phi(10000),scatter_vel(10000),recoil_vel(10000)
      double precision eps1, eps2,fn1,fo2,fo1,dt1,dt2,dt,te,value
      double precision tebf,sc,fc,cexponents(100),cvalues(100)
      integer n_events_passed1,nruns,ival2,runs,mvals,cvals
      integer print,n_events_passed2,J1,J2,n_o2
      real sm,fm,mexponents(100),mvalues(100),fn2
      double precision maxprob,prob,probkssum,pvalue,yieldprob,yieldsum,yieldmean
      real kev(123),ff(123)
     
!     ***********************************************************************

!     Preprocessing

!     ***********************************************************************

!     Define pi
      pi = 4.0*DATAN(1.D0)
      
!     preset basic parameters
!     vary nplot to increase/decrease range of calculated functions
      probkssum = 0.
      nplot = 1000
      yfunc(2,3001) = 1.
      
!     Uncomment the following to allow the target mass and cross section to be read from the command line 
      
!     Set target cross section and mass (could be given from command line in future)
      write(6,*) "Target cross section in units of cm^2: "
      read(*,*) cross_section
      write(6,*) "Target WIMP mass in units of GeV/c^2: "
      read(*,*)  Mw                ! WIMP mass in GeV/c^2

      ! Alternatively set target mass and cross section within program
      !Mw = 60.
      !cross_section = 4D-44

!     Needed for correct outputs
      print = 10

!     rename file as desired to store p-values of simulation
!     is organised for three columns:  WIMP mass,WIMP-nucleon cross section , mean p-value
      open(unit=97,file='simulationgrid.dat',status='replace')

!     Read form factor data
      open(unit=300,file='formfactor.dat',status='old',action='read')

      ! Assign form factor data to arrays
      do i = 1,123
         read(300,*) kev(i) , ff(i)
      end do
      
      nevents = 10000

!     Set preliminary values
      earth_mass = 1.78E-25 ! mass of the Earth
      sigma = 220/SQRT(2.)  ! parameter for Maxwell-Boltzmann Distribution    
      sunv = 220. ! Velocity of sun
      time = 0    ! Time of year
      tp = 153    ! Reference day in year
      yr = 365    ! Length of year
      bins = 1000     ! Bins for histogram 
      Mx = 131.29       ! Xenon mass (GeV/c^2) 
      mass_ratio = Mw/Mx             ! mass ratio
      c = 3e8                   ! Speed of light
      sconversion = 1.11875e-5  ! To convert energies to keV later
      sfactor = Mw/(Mw+Mx)      ! recurring variable in kinematics
      b = 2.2827                ! harmonic oscillator length for 128-Xe (form factor)


!     Clear arrays (failsafe)
      do i = -1,nplot
         func(1,i) = 0.
         func(2,i) = 0.
      end do

      do i = -3000,3000
         xcum(i) = 0.
         ycum(i) = 0.
         zcum(i) = 0.
      end do

      do i = 0,50000
         speed_hist(i) = 0.
         magvel_hist(i) = 0.
      end do

      do i = -180,180
         theta_hist(i) = 0.
         phi_hist(i) = 0.
         ctheta_hist(i) = 0.
      end do

      do i = 0,50000
         scatter_hist(i) = 0.
         recoil_hist(i) = 0.
         energy_hist(i) = 0.
      end do

      
!     Set up Maxwell-Boltzmann Distribution
!     (at present only resolved to 1 km/s, can be improved)
      do i = 0,nplot
         numer = -(i**2.)
         denomin = 2*(sigma**2.)
         expon = exp(numer/denomin)

         func(1,i) = (i**2.)*expon
         func(2,i) = func(2,i-1) + func(1,i)
      end do

!     Normalise manually by reading in largest cumulative value
      max_value = func(2,nplot)

      do i = 0,nplot
         func(2,i) = func(2,i)/max_value
      end do

!     *******************************************************

!     Target simulation loop

!     Creates the target array of 10,000 particles and simulates
!     interaction with the LUX detector

!     *******************************************************

      do j = 1,nevents

         
!     Randomly sample from Maxwell-Boltzmann
!     This uses a numerical analog of the inverse cumulative approach
!     (hard to find cumulative function explicitly)
!     Takes some random numnber [0,1] then sets sample as the velocity that
!     returns that value via the cdf
!     Could be optimised using exact form for cumulative function of Maxwell-Boltzmann
         p(j) = rand(0)
         step = 1
         do while(p(j) > func(2,step))
            step = step + 1
         end do
         speed(j) = real(step)

!     Extend to 3D
100      continue

!     To find unit vectors, sample randomly from unit cube
!     then sample out to unit sphere and normalise vector length
         
      ! Populate a box [-1,1] in 3D
         threed(1,j) = 2*rand(0) - 1
         threed(2,j) = 2*rand(0) - 1
         threed(3,j) = 2*rand(0) - 1

         !     Sample out to a sphere
         if (threed(1,j)**2 + threed(2,j)**2 + threed(3,j)**2 > 1) then
            go to 100
         end if

         ! Find length of vectors and scale to unit vector
         length = threed(1,j)**2 + threed(2,j)**2 + threed(3,j)**2
         threed(1,j) = threed(1,j)/SQRT(length)
         threed(2,j) = threed(2,j)/SQRT(length)
         threed(3,j) = threed(3,j)/SQRT(length)

         
         ! Combine speeds with unit vectors for velocities
         vel(1,j) = threed(1,j)*speed(j)
         vel(2,j) = threed(2,j)*speed(j)
         vel(3,j) = threed(3,j)*speed(j)

!     To shift to Earth frame used Galilean Transformation
!     laid out in Copi & Krauss (see report)

         ! Shift from Galactic to Earth Frame
         vel(1,j) = vel(1,j) - 0.13*sunv*sin(2*pi*(time-tp)/yr)
         vel(2,j) = vel(2,j) + 0.11*sunv*cos(2*pi*(time-tp)/yr)
         vel(3,j) = vel(3,j) - sunv*(1.05+0.06*cos(2*pi*(time-tp)/yr))
         magvel(j) = SQRT((vel(1,j)**2)+(vel(2,j)**2)+(vel(3,j)**2))

         ! Collect velocity distributions
         ival = nint(vel(1,j))
         xcum(ival) = xcum(ival) + 1

         ival = nint(vel(2,j))
         ycum(ival) = ycum(ival) + 1

         ival = nint(vel(3,j))
         zcum(ival) = zcum(ival) + 1

!     **************************************************

!     Calculate randomized centre of momentum scattering angles

!     Scattering angle taken to be random in CoM frame (all equally probable)
         
!     **************************************************

 8       continue

!     This merely defined in 2D plane
!     (conservation of momentum)
!     uses same method as above to sample evenly
         
         twod(1,j) = 2*rand(0) - 1
         twod(2,j) = 2*rand(0) - 1

         ! Sample out to circle
         if ( (twod(1,j)*twod(1,j) + twod(2,j)*twod(2,j))  > 1) then
            go to 8
         end if

         ! Find length and scale to unit vector
         length = twod(1,j)*twod(1,j) + twod(2,j)*twod(2,j)
         twod(1,j) = twod(1,j)/SQRT(length)
         twod(2,j) = twod(2,j)/SQRT(length)

         ! Find angles from this
         tan = twod(2,j)/twod(1,j)
         ! Ensure correct quadrant for each angle
         if (twod(1,j) > 0) then
               
               theta_cm(j) = ATAN(tan)

            else
               
               if (twod(2,j) > 0) then
                  
                  theta_cm(j) = ATAN(tan) + pi
               else

                  theta_cm(j) = ATAN(tan) - pi

               end if

            end if

            if( theta_cm(j).LT.-pi) then
               theta_cm(j) = theta_cm(j) + 2*pi
            else if (theta_cm(j).GT.pi) then
               theta_cm(j) = theta_cm(j) - 2*pi
            end if
            

!     ************************************************

!     Find lab frame scattering angle and recoil angle

!     Again, Galilean Transformation based on ratio of masses

!     ************************************************

         theta(j) = sin(theta_cm(j)) / (cos(theta_cm(j)) + a)
         theta(j) = ATAN(theta(j))

         if (theta(j) > pi) then

            theta(j) = theta(j) - 2.*pi

         end if

         if (theta(j) < -pi) then

            theta(j) = theta(j) + 2.*pi
            
         end if

         phi(j) = (pi - theta_cm(j))/2.

!     **************************************************

!     Find resultant velocities

!     Uses simple kinematics of two body elastic collisions

!     **************************************************

         cmvel = magvel(j)/(1+a)
         scatter_vel(j) = SQRT(a*a + 2.*a*cos(theta_cm(j)) + 1)
         scatter_vel(j) = scatter_vel(j) * cmvel

         recoil_vel(j)=sfactor*magvel(j)*SQRT(2.-2.*cos(theta_cm(j)))

         ! Find momentum transfer and dimensionless variable u (alternative means to get form factor)
         !mtransfer(j) = mw*(magvel(j) - scatter_vel(j))
         !u(j) = 0.5*b*b*mtransfer(j)*mtransfer(j)

!     *****************************************************

!     Now the energies

!     ****************************************************

         new_array(j) = 0.5*mw*magvel(j)*magvel(j)

         scat_energy(j) = 0.5*mw*scatter_vel(j)*scatter_vel(j)
         rec_energy(j) = 0.5*mx*recoil_vel(j)*recoil_vel(j)

!        Convert to keV
         new_array(j) = new_array(j)*sconversion
         scat_energy(j) = scat_energy(j)*sconversion
         rec_energy(j) = rec_energy(j)*sconversion

!     *************************************************************

!    Convert to degrees and bin data for analysis

!     *************************************************************
         
         
         ! Convert to degrees
         theta_cm(j) = theta_cm(j) *180./pi
         theta(j) = theta(j) *180./pi
         phi(j) = phi(j) *180./pi

         ! Bin angles
         ival = nint(theta_cm(j))
         ctheta_hist(ival) = ctheta_hist(ival) + 1

         ival = nint(theta(j))
         theta_hist(ival) = theta_hist(ival) + 1

         ival = nint(phi(j))
         phi_hist(ival) = phi_hist(ival) + 1

         ! Bin speeds
         ival = nint(speed(j))
         speed_hist(ival) = speed_hist(ival) + 1

         ival = nint(magvel(j))
         magvel_hist(ival) = magvel_hist(ival) + 1
         
         ! Bin energy
         ival = nint(new_array(j))
         energy_hist(ival) = energy_hist(ival) + 1
       
         ival = nint(rec_energy(j))
         recoil_hist(ival) = recoil_hist(ival) + 1

         ival = nint(scat_energy(j))
         scatter_hist(ival) = scatter_hist(ival) + 1

      end do

      ! Write velocities to histograms
      open(unit=23,file="velx.dat",status='replace')
      open(unit=24,file="vely.dat",status='replace')
      open(unit=25,file="velz.dat",status='replace')
      
      do i = -1000,1000
         write(23,*) i, xcum(i)
         write(24,*) i, ycum(i)
         write(25,*) i, zcum(i)
      end do

      close(23)
      close(24)
      close(25)

      ! Write angles to histograms

      open(unit=11,file="theta.dat",status='replace')
      open(unit=13,file="phi.dat",status='replace')
      open(unit=32,file="ctheta.dat",status='replace')

      do i = -180,180
         
         write(11,*) i, theta_hist(i)
         write(13,*) i, phi_hist(i)
         write(32,*) i, ctheta_hist(i)
      end do

      close(11)
      close(13)
      close(32)

      open(unit=26,file="speed.dat",status='replace')
      open(unit=27,file="magvel.dat",status='replace')

      do i = 0,2000
         write(26,*) speed_hist(i)
         write(27,*) magvel_hist(i)
      end do

      close(26)
      close(27)

      open(unit=30,file="energy.dat",status='replace')
      open(unit=28,file="recoil.dat",status='replace')
      open(unit=29,file="scatter.dat",status='replace')

!     ******************************************************

!     Uncomment relevant lines to add in form factor
!     Simply finds nearest form factor value and multiplies recoil spectrum by it
!     form factor given normalised to [0,1]
      
      do i = 0,1000
!     Modulate recoil energies by form factor
         
          !step = 0
          !do while(kev(step) < recoil_hist(i))
          !   step = step + 1
          !end do
          !recoil_hist(i) = recoil_hist(i)*ff(step)

         
         recoil_hist(i) = recoil_hist(i)
         write(30,*) i,energy_hist(i)
         write(28,*) i,recoil_hist(i)
         write(29,*) i,scatter_hist(i)

      end do
      
      close(30)
      close(28)
      close(29)

      ! Sort energy spectrum for later
      call sort(nevents,rec_energy)

      do i = 1,nevents
         sorted2(i) = rec_energy(i)
      end do

      n_events_passed2 = nevents

      n_o2 = nevents

      ! Basic parameters of code
      yplot = 3000
      yreference = 1D44
      nplot = 1000              ! Points to plot distributions


      
      ! Read in digitized plot for calibration of mean yield values
      open(unit = 1,file = "yield.dat",status = "old",action = "read")

      do i=1,188
         read(1,*) mlist(i),crlist(i)

!     Multiply by yreference here to scale values so they can be manipulated effectively
         yweight(i) = (crlist(i)*yreference)/mlist(i)
         ycalib(i) = 2.3/yweight(i)
      end do

      close(1)

      
      ! Zero functions (failsafe)
      do i = -1,3000
         nfunc(1,i) = 0.
         nfunc(2,i) = 0.
      end do
      
      ystep = 0

!     Cycle through masslist to find nearest values
      if(Mw < mlist(188)) then
         ystep = 1
         do while (Mw > mlist(ystep))
            ystep = ystep + 1
         end do
         
!     Have the corresponding index, so can allocate the required calibration     

         ! a is the interpolation ratio
         a = (Mw-mlist(ystep-1))/(mlist(ystep)-mlist(ystep-1))
         true_cal = ycalib(ystep-1) + a*(ycalib(ystep)-ycalib(ystep-1))

!     Have transformed calibration by reference, so do so again here to cancel out 
         cross = cross_section*yreference

         yield = cross/Mw
         yield = yield*true_cal

         !write(6,*) "This is the yield"
         !write(6,*) yield

      else
         
         ! Can expand on this to allow extrapolation to consider masses beyond range

         !write(6,*) "Mass too large for now"

      end if

!     Set up Poisson distribution from which to sample yield
!     To reduce computation, use Normal Distribution if over 20
!     (Central Limit Theorem)
      
      ! Use Poisson for low numbers
      if(yield<20) then

         maxprob = 0.

         do i = 0,yplot

            nfunc(1,i) = exp(-yield)
            do j = 1,i
               nfunc(1,i) = nfunc(1,i) * yield
            end do
           

            ! Perform factorial
            do j = 1,i
               nfunc(1,i) = nfunc(1,i)/j
            end do

            ! Use a stopgap to prevent NaNs when values get very small
            if(isnan(nfunc(1,i))) then
               nfunc(1,i) = 0.0
            end if


            nfunc(2,i) = nfunc(2,i-1) + nfunc(1,i)

!     Find max probability for rescaling of probabilities later            
            if( nfunc(1,i) > maxprob) then
               maxprob = nfunc(1,i)
            end if

            ! Limit to stop function return infinities
            if (nfunc(1,i) < 1E-35) then
               go to 61
            end if
            

            !write(6,*) i , nfunc(1,i) , nfunc(2,i)

         end do

         

 61      continue

      else
         
!     Use Normal Distribution if the yield is large to reduce computation
!     Though Normal Dist is continuous, have here simply used bins of width 1 keV
!     This is ultimately an approximation so be limited in its accuracy
         
         norm = 1/SQRT(2*pi*yield)

         maxprob = 0.

         do i = 0,yplot
            numerator = -((i-yield)**2)
            denominator = 2*yield
            exponential = exp(numerator/denominator)
            nfunc(1,i) = norm * exponential

            ! Stopgap to prevent NaNs when values get very small
            if(isnan(nfunc(1,i))) then
               nfunc(1,i) = 0.0
            end if

              
            nfunc(2,i) = nfunc(2,i-1) + nfunc(1,i)

            !write(6,*) "nfunc"
            !write(6,*) nfunc(1,i)
            
!     Find max probability for rescaling of probabilities later
            if( nfunc(1,i) > maxprob) then
               maxprob = nfunc(1,i)
            end if

            ! Limit to stop function returning infinities
            if (nfunc(1,i) < 1E-35) then
               go to 21
            end if
            

            !write(6,*) nfunc(1,i) , nfunc(2,i)

         end do

 21      continue

      end if
     
!     ******************************************************

!     Start of loop that generates the samples for each point in the phase space

!     ******************************************************

!     Set up log space to vary mass and cross section
      sm = ALOG10(10.)
      fm = ALOG10(4800.)
      sc = DLOG10(2.D-45)
      fc = DLOG10(1.D-43)

      nruns = 50

      do i = 1,100
         mexponents(i) = sm + i*(fm - sm)/100.
         cexponents(i) = sc + i*(fc - sc)/100.
         mvalues(i) = 10.**mexponents(i)
         cvalues(i) = 10.**cexponents(i)
         !write(6,*) mvalues(i) , cvalues(i)
      end do

      ! Loop through cross sections
      do cvals = 1,100

         ! Serves to give a percentage of simulation complete at the end of every mass loop
         write(6,*) cvals

         ! Loop through masses
         do mvals = 1,100

            Mw = mvalues(mvals)
            cross_section = cvalues(cvals)
            
!           Reassign relevant variables
            mass_ratio = Mw/Mx             ! mass ratio
            sfactor = Mw/(Mw+Mx)      ! recurring variable in kinematics
            !Mw = 100
            !cross_section = 9D-45

            write(6,*) "Running sim with: "
            write(6,*) Mw , cross_section
            
!           Ensure sums are zeroed between phase space locations
            probkssum = 0.
            yieldsum = 0.
            
            ! Perform 50 iterations at each point
            do runs = 1,nruns

 52            continue

!     Set the number of datapoints for the sample and parent distribution

              

!     ***********************************************************
      
!     Yield section of code (subroutines couldn't pass values correctly in and out)

!     ***********************************************************       

      ! Basic processing
      do i = -1,3000
         yfunc(1,i) = 0.
         yfunc(2,i) = 0.
      end do
      
      ystep = 0

      ! Find calibration for mass and cross section to estimate mean yield
!     Cycle through masslist to find nearest values
      if(Mw < mlist(188)) then
         ystep = 1
         do while (Mw > mlist(ystep))
            ystep = ystep + 1
         end do
         
!     Have the corresponding index, so can allocate the required calibration     

         ! a is the interpolation ratio
         a = (Mw-mlist(ystep-1))/(mlist(ystep)-mlist(ystep-1))
         true_cal = ycalib(ystep-1) + a*(ycalib(ystep)-ycalib(ystep-1))

!     Have transformed calibration by reference, so do so again here to cancel out 
         cross = cross_section*yreference

         yield = cross/Mw
         yield = yield*true_cal

         !write(6,*) "This is the yield"
         !write(6,*) yield

      else
         
         ! Can expand on this to allow extrapolation to consider masses beyond range

         !write(6,*) "Mass too large for now"

      end if

      ! Use Poisson for low numbers
      if(yield<20) then

         do i = 0,yplot

            yfunc(1,i) = exp(-yield)
            do j = 1,i
               yfunc(1,i) = yfunc(1,i) * yield
            end do
           

            ! Perform factorial
            do j = 1,i
               yfunc(1,i) = yfunc(1,i)/j
            end do
            

            ! Use a stopgap to prevent NaNs when values get very small
            if(isnan(yfunc(1,i))) then
               yfunc(1,i) = 0.0
            end if
            

            yfunc(2,i) = yfunc(2,i-1) + yfunc(1,i)

            ! To reduce computation end function when suitably close to 1 cumulatively
            if(abs((yfunc(2,i+1)-1))<0.000000000001) then
               if (yfunc(2,i+1).GE.1) then
                  go to 10
               end if
               
            end if

            !write(6,*) i , yfunc(1,i) , yfunc(2,i)

         end do

 10      continue

      else
         
! Use Normal Distribution if the yield is large to reduce computation
         norm = 1/SQRT(2*pi*yield)

         do i = 0,yplot
            numerator = -((i-yield)**2)
            denominator = 2*yield
            exponential = exp(numerator/denominator)
            yfunc(1,i) = norm * exponential

            ! Stopgap to prevent NaNs when values get very small
            if(isnan(yfunc(1,i))) then
               yfunc(1,i) = 0.0
            end if

              
            yfunc(2,i) = yfunc(2,i-1) + yfunc(1,i)
            
                        
            ! To reduce computation end function when suitably close to 1 cumulatively
            if(abs((yfunc(2,i+1)-1))<0.000000000001) then
               go to 20
            end if

            !write(6,*) yfunc(1,i) , yfunc(2,i)

         end do

 20      continue

      end if

      ! Now sample from this distribution
      nsample = 0

      ! Uniform random number in interval [0,1] using system clock method
      r = randtest()
      !write(6,*) "Random number is: "
      !write(6,*) r

      ! Loop through cdf until cdf exceeds r
      do while(r > yfunc(2,nsample))
         nsample = nsample + 1
      end do
      
!     Don't want to attempt analysis on empty datasets, so repeat process
!     until non-zero dataset found
!     Could potentially be made more efficient by rescaling functions
!     so that zero is no longer an option, and other value probabilities
!     rescaled to give same probabilities relative to each other
      
      if (nsample == 0) then
         go to 52
      end if

      !write(6,*) "This is maxprob"
      !write(6,*) maxprob

      ! Record probability of detecting this or fewer events
      yieldprob = nfunc(1,nsample)/maxprob
      !write(6,*) nsample
      
      !write(6,*) yieldprob
      
      nevents = 10000

      write(6,*) "Number of events simulated is: "
      write(6,*) nsample

      !write(6,*) "The parent distribution has this many events: "
      !write(6,*) nevents


!     *************************************************************

!     Yield complete, move on to sample simulation

!     *************************************************************

!     Set lengths of arrays based on yield (dynamic arrays so can vary in length)

      allocate ( sp(nsample))
      allocate ( sthreed(3,nsample))
      allocate ( sspeed(nsample))
      allocate ( svel(3,nsample))
      allocate ( smagvel(nsample))
      allocate ( stwod(2,nsample))
      allocate ( stheta_cm(nsample))
      allocate ( stheta(nsample))
      allocate ( sphi(nsample))
      allocate ( sscatter_vel(nsample))
      allocate ( srecoil_vel(nsample))
      allocate ( smtransfer(nsample))
      allocate ( su(nsample))
      allocate ( snew_array(nsample))
      allocate ( sscat_energy(nsample))
      allocate ( srec_energy(nsample))
      allocate ( sorted1(nsample))
      
!     Clear memory in arrays
      
      do i = -1,nplot
         sfunc(1,i) = 0.
         sfunc(2,i) = 0.
      end do

      do i = -3000,3000
         sxcum(i) = 0.
         sycum(i) = 0.
         szcum(i) = 0.
      end do

      do i = 0,50000
         sspeed_hist(i) = 0.
         smagvel_hist(i) = 0.
      end do

      do i = -180,180
         stheta_hist(i) = 0.
         sphi_hist(i) = 0.
         sctheta_hist(i) = 0.
      end do

      do i = 0,50000
         sscatter_hist(i) = 0.
         ssrecoil_hist(i) = 0.
         senergy_hist(i) = 0.
      end do



!     Set up Maxwell-Boltzmann Distribution
      do i = 0,nplot

!     Set up velocity distribution for simulated sample

         ! Constant for now
         numer = -(i**2.)
         denomin = 2*(sigma**2.)
         expon = exp(numer/denomin)

         sfunc(1,i) = (i**2.)*expon
         sfunc(2,i) = sfunc(2,i-1) + sfunc(1,i)

      end do

!     Normalise
      smax_value = sfunc(2,nplot)
     
      do i = 0,nplot
         sfunc(2,i) = sfunc(2,i)/smax_value
      end do

!     *******************************************************************

!     Sample loop

!     *******************************************************************

!     Perform required changes

!     For approximated Necib distribution, set sunv = 190
!     For extension, can define new velocity profile in full
      sunv = 190.

      
      ! Repeat sampling process as done for the target
      do j = 1,nsample

         sp(j) = rand(0)
         sstep = 1
         do while(p(j) > sfunc(2,sstep))
            sstep = sstep + 1
         end do
         sspeed(j) = real(sstep)

 101     continue

      ! Populate 3D box
         sthreed(1,j) = 2*rand(0) - 1
         sthreed(2,j) = 2*rand(0) - 1
         sthreed(3,j) = 2*rand(0) - 1

      !     Sample out to a sphere
         if(sthreed(1,j)**2+sthreed(2,j)**2+sthreed(3,j)**2 > 1) then
            go to 101
         end if

       ! Find length of vectors and scale to unit vector
         slength = sthreed(1,j)**2 + sthreed(2,j)**2 + sthreed(3,j)**2
         sthreed(1,j) = sthreed(1,j)/SQRT(slength)
         sthreed(2,j) = sthreed(2,j)/SQRT(slength)
         sthreed(3,j) = sthreed(3,j)/SQRT(slength)

         ! Combine speeds with unit vectors for velocities
         svel(1,j) = sthreed(1,j)*sspeed(j)
         svel(2,j) = sthreed(2,j)*sspeed(j)
         svel(3,j) = sthreed(3,j)*sspeed(j)

!     ******************************************************

!     Alter sun velocity to see what happens

!     ******************************************************
         
         ! Shift from Galactic to Earth Frame
         svel(1,j) = svel(1,j) - 0.13*sunv*sin(2*pi*(time-tp)/yr)
         svel(2,j) = svel(2,j) + 0.11*sunv*cos(2*pi*(time-tp)/yr)
         svel(3,j) = svel(3,j) - sunv*(1.05+0.06*cos(2*pi*(time-tp)/yr))
         smagvel(j) = SQRT((svel(1,j)**2)+(svel(2,j)**2)+(svel(3,j)**2))

         ! Collect velocity distributions
         ival = nint(svel(1,j))
         sxcum(ival) = sxcum(ival) + 1

         ival = nint(svel(2,j))
         sycum(ival) = sycum(ival) + 1

         ival = nint(svel(3,j))
         szcum(ival) = szcum(ival) + 1

         ival = nint(sspeed(j))
         sspeed_hist(ival) = sspeed_hist(ival) + 1

         ival = nint(smagvel(j))
         smagvel_hist(ival) = smagvel_hist(ival) + 1

!     **************************************************

!     Calculate randomized centre of momentum scattering angles

!     **************************************************

18        continue
         
         stwod(1,j) = 2*rand(0) - 1
         stwod(2,j) = 2*rand(0) - 1

         ! Sample out to circle
         if ( (stwod(1,j)*stwod(1,j) + stwod(2,j)*stwod(2,j))  > 1) then
            go to 18
         end if

         ! Find length and scale to unit vector
         slength = stwod(1,j)*stwod(1,j) + stwod(2,j)*stwod(2,j)
         stwod(1,j) = stwod(1,j)/SQRT(slength)
         stwod(2,j) = stwod(2,j)/SQRT(slength)

         ! Find angles from this
         stan = stwod(2,j)/stwod(1,j)
         ! Ensure correct quadrant for each angle
         if (stwod(1,j) > 0) then
               
               stheta_cm(j) = ATAN(stan)

            else
               
               if (stwod(2,j) > 0) then
                  
                  stheta_cm(j) = ATAN(stan) + pi
               else

                  stheta_cm(j) = ATAN(stan) - pi

               end if

            end if

            if( stheta_cm(j).LT.-pi) then
               stheta_cm(j) = stheta_cm(j) + 2*pi
            else if (stheta_cm(j).GT.pi) then
               stheta_cm(j) = stheta_cm(j) - 2*pi
            end if

!     ************************************************

!     Find lab fram scattering angle and recoil angle

!     ************************************************

         stheta(j) = sin(stheta_cm(j)) / (cos(stheta_cm(j)) + a)
         stheta(j) = ATAN(stheta(j))

         if (stheta(j) > pi) then

            stheta(j) = stheta(j) - 2.*pi

         end if

         if (stheta(j) < -pi) then

            stheta(j) = stheta(j) + 2.*pi
            
         end if

         sphi(j) = (pi - stheta_cm(j))/2.

!     **************************************************

!     Find resultant velocities

!     **************************************************

         scmvel = smagvel(j)/(1+a)
         sscatter_vel(j) = SQRT(a*a + 2.*a*cos(stheta_cm(j)) + 1)
         sscatter_vel(j) = sscatter_vel(j) * scmvel

         srecoil_vel(j)=sfactor*smagvel(j)*SQRT(2.-2.*cos(stheta_cm(j)))


         ! Find momentum transfer and dimensionless variable u (for form factors)
         smtransfer(j) = mw*(smagvel(j) - sscatter_vel(j))
         su(j) = 0.5*b*b*smtransfer(j)*smtransfer(j)

!     *****************************************************

!     Now the energies

!     ****************************************************

         snew_array(j) = 0.5*mw*smagvel(j)*smagvel(j)

         sscat_energy(j) = 0.5*mw*sscatter_vel(j)*sscatter_vel(j)
         srec_energy(j) = 0.5*mx*srecoil_vel(j)*srecoil_vel(j)

!        Convert to keV
         snew_array(j) = snew_array(j)*sconversion
         sscat_energy(j) = sscat_energy(j)*sconversion
         srec_energy(j) = srec_energy(j)*sconversion

!     *************************************************************

!    Convert to degrees and bin data for analysis

!     *************************************************************
         
         
         ! Convert to degrees
         stheta_cm(j) = stheta_cm(j) *180./pi
         stheta(j) = stheta(j) *180./pi
         sphi(j) = sphi(j) *180./pi

         ! Bin angles
         ival = nint(stheta_cm(j))
         sctheta_hist(ival) = sctheta_hist(ival) + 1

         ival = nint(stheta(j))
         stheta_hist(ival) = stheta_hist(ival) + 1

         ival = nint(sphi(j))
         sphi_hist(ival) = sphi_hist(ival) + 1
         
         ! Bin energy
         ival = nint(snew_array(j))
         senergy_hist(ival) = senergy_hist(ival) + 1
       
         ival2 = nint(srec_energy(j))
         ssrecoil_hist(ival2) = ssrecoil_hist(ival2) + 1

         ival = nint(sscat_energy(j))
         sscatter_hist(ival) = sscatter_hist(ival) + 1



      end do

!     ***********************************************

!     Now perform the two-sample K.S Test

!     ***********************************************


!     Compute observed cumulative function for each set of recoil data

      ! Sort the sample dataset (target already done)
      call sort(nsample, srec_energy)

      do i = 1,nsample
         sorted1(i) = srec_energy(i)
      end do

      n_events_passed1 = nsample

!     **** Perform a Kolmogorov-Smirnov Test, returning a p-value
!     **** code here standard k-s fortran

      J1=1                  !  //Next value of DATA1 to be processed.
      J2=1                  !  //Ditto, DATA2.
      Fo1=0.                !   //value of c.d.f. before the next step.
      Fo2=0.                !  //Ditto, for DATA2.
      D=0.

      !write(6,*) j1,j2,fo1,fo2,n_events_passed1,n_events_passed2
          
 1    IF(J1.LE.n_events_passed1.AND.J2.LE.n_o2) THEN !   //If we are not done...

         IF(sorted1(J1).LT.sorted2(J2))THEN !    //Next step is in DATA1.
              FN1=dble(J1)/dble(n_events_passed1)
                
!             DT=AMAX1(DABS(FN1-Fo2),DABS(Fo1-Fo2))

              dt1=dabs(fn1-fo2)
              dt2=dabs(fo1-fo2)
              if(dt1.gt.dt2)dt=dt1
              if(dt2.gt.dt1)dt=dt2
                

              IF(DT.GT.D)D=DT
                
                Fo1=FN1
                J1=J1+1
                
              ELSE               ! //Next step is in DATA2.
                FN2=dble(J2)/dble(n_o2)

!                DT=AMAX1(DABS(FN2-Fo1),DABS(Fo2-Fo1))

                dt1=dabs(fn2-fo1)
                dt2=dabs(fo2-fo1)
                if(dt1.gt.dt2)dt=dt1
                if(dt2.gt.dt1)dt=dt2

                IF (DT.GT.D)D=DT
                Fo2=FN2
                J2=J2+1
             ENDIF
             Go To 1
          ENDIF


          value=DSQRT(dble(n_events_passed1)*dble(n_o2)/(dble(n_events_passed1)+dble(n_o2)))*D

          EPS1=0.001
          EPS2=1.d-7
          A2=-2.*value**2
          FAC=2.
          PRoBKS=0.
          TEBF=0.             ! //Previous term in sum.
          Do 11 J=1,100
             TE=FAC*EXP(A2*dble(J)**2)
             PRoBKS=PRoBKS+TE
             IF(DABS(TE).LT.EPS1*TEBF.oR.DABS(TE).LT.EPS2*PRoBKS)goto 12
             FAC=-FAC           !   //Alternating signs in sum.
             TEBF=DABS(TE)
 11       CoNTINUE
          write(6,*)'KS did not minimise'
          PROBKS=0.             ! //Get here only by failing to converge.
          

 12       continue

                                !D is the maximum difference in the cumulative distr
                                !probks is the p-value for match of the two samples
                                !as usual, p=1 is the chance that the null hypothesis is true
                                !here, the null hypothesis is that the two distributions
                                !are statistically the same.
    
          !write(6,*) probarray(runs)
          !write(6,*) "Probks"
          !write(6,*) probks
          !write(6,*) yieldprob,maxprob,nfunc(1,nsample)

!     Multiply two contributing probabilities for overall p-value
          prob = probks*yieldprob
          write(6,*) prob,probks,yieldprob
!     Sum p-values to take average
          probkssum  = probkssum + prob

          
!     Must deallocate sample arrays so they can be reinitialised in next iteration
      deallocate ( sthreed)
      deallocate ( sspeed)
      deallocate ( svel)
      deallocate ( smagvel)
      deallocate ( stwod)
      deallocate ( stheta_cm)
      deallocate ( stheta)
      deallocate ( sphi)
      deallocate ( sscatter_vel)
      deallocate ( srecoil_vel)
      deallocate ( smtransfer)
      deallocate ( su)
      deallocate ( sp)
      deallocate ( snew_array)
      deallocate ( sscat_energy)
      deallocate ( srec_energy)
      deallocate ( sorted1)

     
      
      end do                  ! Loop over runs

!     Divide by number of runs to get average
      pvalue = probkssum/nruns
      !write(6,*) pvalue

      write(97,*) Mw,cross_section,pvalue
      write(6,*) Mw,cross_section,pvalue


      end do                 ! Loop over mvals

      end do                 ! Loop over cvals
      

      end program simcode
