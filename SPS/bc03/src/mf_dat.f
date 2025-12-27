	SUBROUTINE MF_DAT(J,W,IDXLBL)

!	Returns 3 bandpasses that define a UV spectral index according to M. Fanelli et al. (1987, 1990)
!	J = index number in Maraston et al. (2009) Table 1.

!	w(1), w(2): define blue pseudo continuum
!	w(3), w(4): define central index band
!	w(5), w(6): define red  pseudo continuum

!	Array declaration
	character*40 idxlbl
	real w(6)

	IF     (J ==  1) THEN
!		Si III, Si II, O I
		idxlbl = 'Index 01: BL 1302 (A)'
		w(1) = 1270.000
		w(2) = 1290.000
		w(3) = 1292.000
		w(4) = 1312.000
		w(5) = 1345.000
		w(6) = 1365.000
	ELSEIF (J ==  2) THEN
!		Si IV 1393.8; 1402.8
		idxlbl = 'Index 02: Si IV (A)'
		w(1) = 1345.000
		w(2) = 1365.000
		w(3) = 1387.000
		w(4) = 1407.000
		w(5) = 1475.000
		w(6) = 1495.000
	ELSEIF (J ==  3) THEN
!		C II 1429, Si III 1417, Fe IV, Fe V
		idxlbl = 'Index 03: BL 1425 (A)'
		w(1) = 1345.000
		w(2) = 1365.000
		w(3) = 1415.000
		w(4) = 1435.000
		w(5) = 1475.000
		w(6) = 1495.000
	ELSEIF (J ==  4) THEN
!		Fe V + 20 additional lines
		idxlbl = 'Index 04: Fe 1453 (A)'
		w(1) = 1345.000
		w(2) = 1365.000
		w(3) = 1440.000
		w(4) = 1466.000
		w(5) = 1475.000
		w(6) = 1495.000
	ELSEIF (J ==  5) THEN
!		C IV 1548, in absorption
		idxlbl = 'Index 05: C IV 1548 absorption (A)'
		w(1) = 1500.000
		w(2) = 1520.000
		w(3) = 1530.000
		w(4) = 1550.000
		w(5) = 1577.000
		w(6) = 1597.000
	ELSEIF (J ==  6) THEN
!		C IV 1548, central band
		idxlbl = 'Index 06: C IV 1548 central band (A)'
		w(1) = 1500.000
		w(2) = 1520.000
		w(3) = 1540.000
		w(4) = 1560.000
		w(5) = 1577.000
		w(6) = 1597.000
	ELSEIF (J ==  7) THEN
!		C IV 1548, in emission
		idxlbl = 'Index 07: C IV 1548 emission (A)'
		w(1) = 1500.000
		w(2) = 1520.000
		w(3) = 1550.000
		w(4) = 1570.000
		w(5) = 1577.000
		w(6) = 1597.000
!	ELSEIF (J ==   ) THEN
!		C IV 1548 (P-Cygni)
!		idxlbl = 'Index   : C IV P-Cygni (A)'
!		w(1) = 1530.000
!		w(2) = 1550.000
!		w(3) = 1550.000
!		w(4) = 1550.000
!		ojo ojo ojo
!		w(4) = 1555.000
!		w(5) = 1550.000
!		w(6) = 1570.000
	ELSEIF (J ==  8) THEN
!		Fe IV
		idxlbl = 'Index 08: BL 1617 (A)'
		w(1) = 1577.000
		w(2) = 1597.000
		w(3) = 1604.000
		w(4) = 1630.000
		w(5) = 1685.000
		w(6) = 1705.000
	ELSEIF (J ==  9) THEN
!		C I 1656.9, Al II 1670.8
		idxlbl = 'Index 09: BL 1664 (A)'
		w(1) = 1577.000
		w(2) = 1597.000
		w(3) = 1651.000
		w(4) = 1677.000
		w(5) = 1685.000
		w(6) = 1705.000
	ELSEIF (J == 10) THEN
!		N IV 1718.6, Si IV 1722.5,1727.4, Al II
		idxlbl = 'Index 10: BL 1719 (A)'
		w(1) = 1685.000
		w(2) = 1705.000
		w(3) = 1709.000
		w(4) = 1729.000
		w(5) = 1803.000
		w(6) = 1823.000
	ELSEIF (J == 11) THEN
!		Al II, Al III, Fe II, Fe III
		idxlbl = 'Index 11: BL 1853 (A)'
		w(1) = 1803.000
		w(2) = 1823.000
		w(3) = 1838.000
		w(4) = 1868.000
		w(5) = 1885.000
		w(6) = 1915.000
	ELSEIF (J == 12) THEN
!		Fe II 2402
		idxlbl = 'Index 12: Fe II 2402 (A)'
		w(1) = 2285.000
		w(2) = 2325.000
		w(3) = 2382.000
		w(4) = 2422.000
		w(5) = 2432.000
		w(6) = 2458.000
	ELSEIF (J == 13) THEN
!		Perhaps Fe I
		idxlbl = 'Index 13: BL 2538 (A)'
		w(1) = 2432.000
		w(2) = 2458.000
		w(3) = 2520.000
		w(4) = 2556.000
		w(5) = 2562.000
		w(6) = 2588.000
	ELSEIF (J == 14) THEN
!		Fe II 2609
		idxlbl = 'Index 14: Fe II 2609 (A)'
		w(1) = 2562.000
		w(2) = 2588.000
		w(3) = 2596.000
		w(4) = 2622.000
		w(5) = 2647.000
		w(6) = 2673.000
	ELSEIF (J == 15) THEN
		idxlbl = 'Index 15: Mg II (A)'
		w(1) = 2762.000
		w(2) = 2782.000
		w(3) = 2784.000
		w(4) = 2814.000
		w(5) = 2818.000
		w(6) = 2838.000
	ELSEIF (J == 16) THEN
		idxlbl = 'Index 16: Mg I (A)'
		w(1) = 2818.000
		w(2) = 2838.000
		w(3) = 2839.000
		w(4) = 2865.000
		w(5) = 2906.000
		w(6) = 2936.000
	ELSEIF (J == 17) THEN
		idxlbl = 'Index 17: Mg wide (A)'
		w(1) = 2470.000
		w(2) = 2670.000
		w(3) = 2670.000
		w(4) = 2870.000
		w(5) = 2930.000
		w(6) = 3130.000
	ELSEIF (J == 18) THEN
		idxlbl = 'Index 18: Fe I (A)'
		w(1) = 2906.000
		w(2) = 2936.000
		w(3) = 2965.000
		w(4) = 3025.000
		w(5) = 3031.000
		w(6) = 3051.000
	ELSEIF (J == 19) THEN
!		Al I 3092, Fe I 3091.6
		idxlbl = 'Index 19: BL 3096 (A)'
		w(1) = 3031.000
		w(2) = 3051.000
		w(3) = 3086.000
		w(4) = 3106.000
		w(5) = 3115.000
		w(6) = 3155.000
	ELSEIF (J == 20) THEN
!		CIV stellar absorption:
!			defining the continuum as the median of the flux values in [1510, 1525] and [1560, 1575],
!			and directly-integrating the absorption in [1525, 1548] 
		idxlbl = 'CIV absorption (A)'
		w(1) = 1510.000
		w(2) = 1525.000
		w(3) = 1525.000
		w(4) = 1548.000
		w(5) = 1560.000
		w(6) = 1575.000
	ELSEIF (J == 21) THEN
!		HeII stellar emission:
!			define the continuum as the median in [1600, 1700]
!			and directly-integrate the line flux in [1635, 1645].
		idxlbl = 'HeII emission (A)'
		w(1) = 1600.000
!		w(2) = 1700.000
		w(2) = 1650.000
		w(3) = 1635.000
		w(4) = 1645.000
		w(5) = 1650.000
!		w(5) = 1600.000
		w(6) = 1700.000
	ELSE
		write (6,*) 'Unknown Fanelli index:',j
		stop
	ENDIF
	return
	end
