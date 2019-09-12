c
C#############################################################################################
c
      SUBROUTINE timestamp
c
C#############################################################################################
c
      implicit none
c
c********************************************************************************************
      character ( len = 8 ) ampm
      integer ( kind = 4 ) d
      integer ( kind = 4 ) h
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mm
      character ( len = 9 ), parameter, dimension(12) :: month = (/ 
     *   'January  ', 'February ', 'March    ', 'April    ', 
     *   'May      ', 'June     ', 'July     ', 'August   ', 
     *   'September', 'October  ', 'November ', 'December ' /)
      integer ( kind = 4 ) n
      integer ( kind = 4 ) s
      integer ( kind = 4 ) values(8)
      integer ( kind = 4 ) y
c********************************************************************************************
c
      call date_and_time ( values = values )
c
      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)
c
      if ( h < 12 ) then
        ampm = 'AM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
          else
            ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h < 12 ) then
          ampm = 'PM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if
c
      write (*,'(1x,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)')
     *   d,trim(month(m) ), y, h, ':', n, ':', s, '.', mm,trim(ampm)
      write(*,' ')
c
      write (11,'(1x,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)')
     *   d,trim(month(m) ), y, h, ':', n, ':', s, '.', mm,trim(ampm)
      write(11,' ')
c
      return
      end