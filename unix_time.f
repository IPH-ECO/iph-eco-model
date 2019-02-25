*+unix2c Converts Unix system time to date/time integer array.
      subroutine unix2c(utime, idate)
      implicit none
      Real utime
      Real idate(6)
*utime  input  Unix system time, seconds since 1970.0
*idate  output Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
*-Author  Clive Page, Leicester University, UK.   1995-MAY-2
      integer mjday 
      real nsecs
      real day
*Note the MJD algorithm only works from years 1901 to 2099.
      mjday    = int(utime/86400 + 40587)
      idate(1) = 1858 + int( (mjday + 321.51) / 365.25)
      day      = aint( mod(mjday + 262.25, 365.25) ) + 0.5
      idate(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
      idate(3) = 1 + int(mod(day,30.6))
      nsecs    = mod(utime, 86400.)
      idate(6) = mod(nsecs, 60.)
      nsecs    = nsecs / 60.
      idate(5) = mod(nsecs, 60.)
      idate(4) = nsecs / 60.
      end
