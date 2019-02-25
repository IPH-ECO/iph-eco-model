SUBROUTINE JULIANDAY(I,J,K,JD)

    !   Program for finding the Julian date for a given year (I), month (J) and
    !   day of the month (K). It works even if you give the day of the year.
    !   In this case press <Return> at the request of the month, and it will
    !   automatically compute the day and month of the year.
    !   The algorithm is valid for any Gregorian date producing a Julian date
    !   greater than zero.
    !   (From Fliegel H.F. and van Flandern T.C. 1968. Comm. ACM, 11, p.657)

    !   Author: Mauro Orlandini - Code 666.0
    !           NASA/Goddard Space Flight Center
    !           Laboratory for High Energy Astrophysics
    !           Greenbelt, MD 20771 - USA
    !           E-mail:   6224::ORLANDINI  (SPAN)
    !                     ORLANDINI@LHEAVX.GSFC.NASA.GOV  (Internet)
    !                     ORLANDINI@LHEAVX.SPAN.NASA.GOV  (Bitnet)

    !   Creation Date: Feb 28, 1991

    !   Revisions: Apr 15, 1991   Modified in order that J=0, K=1 gives the same result of J=1, K=1.
    !              Apr 23, 1991   Added the computation of the day of the week.

    ! Called in routines 0-MAIN, PHYTOPLANCTON, LEVENTO, WRITETIMESERIES
    ! Call the following routines: MONTH and DAYWEEK (both internal)

	CHARACTER*3 DW,M
	LOGICAL DAYEAR
	DAYEAR = .FALSE.

	IF (I.LT.100.AND.I.GT.0) WRITE(*,2) I,I

	IF (J.EQ.0) THEN
	   DAYEAR = .TRUE.
	   K1=K
	   K=1
	   J=1
	ELSE 
	   CALL MONTH(J,M)
	ENDIF

	JD1 = K-32075+1461*(I+4800+(J-14)/12)/4
	JD2 = 367*(J-2-(J-14)/12*12)/12
	JD3 = 3*((I+4900+(J-14)/12)/100)/4

	IF (DAYEAR) THEN
	   JD = JD1+JD2-JD3+K1-1721059
	   L = JD + 68569
	   N = 4*L/146097
	   L = L - (146097*N + 3)/4
	   II = 4000*(L + 1)/1461001
	   L = L - 1461*II/4 + 31
	   J = 80*L/2447
	   K = L - 2447*J/80
	   L = J/11
	   J = J + 2 - 12*L
	   II = 100*(N - 49) + II + L
	   CALL MONTH (J,M)
	   CALL DAYWEEK (JD,DW)
	ELSE
	   JD = JD1+JD2-JD3-1721059
	   CALL DAYWEEK (JD,DW)
	ENDIF

1	FORMAT(2X,'Year (I4)             : ',$)
2	FORMAT(2X,'Warning! I will intend ',I2,' and not 19',I2,' !')
3	FORMAT(2X,'Month (I2 or <Ret>)   : ',$)
4	FORMAT(2X,'Day (of month or year): ',$)
5	FORMAT(2X,'JD corresponding to ',I4,'/',I3,': ',I8, ' (',A,', ',A,' ',I2,', ',I4,')')
6	FORMAT(2X,'JD corresponding to ',A,' ',I2,', ',I4,' (',A,'): ',I8)

100	FORMAT(I4)

END SUBROUTINE JULIANDAY

! ----------------------------------------------------
!                   SUBROUTINES
! ----------------------------------------------------

SUBROUTINE MONTH (J,M)
	CHARACTER*3 M
	IF (J.EQ.1) THEN
	   M='Jan'
	ELSE IF (J.EQ.2) THEN
	   M='Feb'
	ELSE IF (J.EQ.3) THEN
           M='Mar'
	ELSE IF (J.EQ.4) THEN
	   M='Apr'
	ELSE IF (J.EQ.5) THEN
	   M='May'
	ELSE IF (J.EQ.6) THEN
	   M='Jun'
	ELSE IF (J.EQ.7) THEN
	   M='Jul'
	ELSE IF (J.EQ.8) THEN
	   M='Aug'
	ELSE IF (J.EQ.9) THEN
	   M='Sep'
	ELSE IF (J.EQ.10) THEN
	   M='Oct'
	ELSE IF (J.EQ.11) THEN
	   M='Nov'
	ELSE
	   M='Dec'
	ENDIF
	RETURN
END SUBROUTINE MONTH

! ----------------------------------------------------

SUBROUTINE DAYWEEK (JD,DW)
	CHARACTER*3 DW
	D=(JD+1.)/7.
	NDW=INT((D-INT(D))*7+0.5)
	IF (NDW.EQ.0) THEN
	   DW='Sun'
	ELSE IF (NDW.EQ.1) THEN
	   DW='Mon'
	ELSE IF (NDW.EQ.2) THEN
	   DW='Tue'
	ELSE IF (NDW.EQ.3) THEN
	   DW='Wed'
	ELSE IF (NDW.EQ.4) THEN
	   DW='Thu'
	ELSE IF (NDW.EQ.5) THEN
	   DW='Fri'
	ELSE
   	   DW='Sat'
	ENDIF
	RETURN
END SUBROUTINE DAYWEEK

