!     This is an attempt to generate a Gaussian distribution, sample values from
!     it and use these to form a histogram
!     As yet, subroutines and functions have not been used for the main code to aid 
!     efficiency
!     This model will later be used to simulate dark matter velocity distributions

      program secondtry

      implicit none

!     declare variables
      real p(1000000), r1(10), r2(1000000), r3(500), h(1000)
      integer i,s,step, ival, j,nevents,q, i2,i3,i4,o
      real x(2,1000),y(2,1000),f(2,1000)
      real m,k,T
      real pi,sigma,N,mean,max_value,exponent,numerator,denominator
      

!     Define pi internally
      pi = 4.0*DATAN(1.D0)

!     Failsafe
      do i2=1,1000
         h(i2)=0
      enddo
      
!     Set preliminary values for the distribution
      m = 1.78E-25
      k = 1.38E-23
      T = 300
      sigma = 220/SQRT(2.)
      N = SQRT(m/(2*pi*k*T))
      mean = 220.
      nevents = 10000

!     Create our list of randomly generated numbers [1-1000] held in p
      do i3 = 1,nevents
         p(i3) = rand(0)
         !write(6,*) p(i3)
      end do

!     Set up our arrays holding f(x) and the cumulative function
!     x holds the function f(x)
!     y holds F(x), both having x values in (1,i)    
      y(1,1) = 0.0
      y(2,1) = 0.0
      do i = 1,1000
         x(1,i) = (i-1)
         y(1,i) = x(1,i)
         numerator = -((x(1,i))**2)
         denominator = 2*(sigma**2.)
         exponent = exp(numerator/denominator)
         x(2,i) = ((x(1,i))**2)*exponent
         y(2,i+1) = y(2,i) + (x(2,i))
      enddo

!     normalise y set up so it goes from 0 to 1
!     Build f as the cumulative function scaled to 1000

      max_value = y(2,1000)
      do i = 1,1000
         y(2,i) = y(2,i)/max_value
      end do 
      
!     Now build up a distribution of 'nevents' points plucked from f(x)
!     s labels the each random number
!     step is an  arbritrary integer used to count within the index
      do s = 1,nevents
         step = 0
         do while (p(s) > y(2,step))
            step = step + 1
         end do

         r2(s) = y(1,step)
      end do

!     Bin data in order to generate histogram
      do j=1,nevents
         ival=nint(r2(j))
         h(ival)=h(ival)+1
      enddo

!     Print out the histogram data held in h to be graphed
      do q = 1,1000
         write(6,*) h(q) 
      end do
     

      end program secondtry
      

