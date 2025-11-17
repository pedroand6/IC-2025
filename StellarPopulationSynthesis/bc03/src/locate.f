      SUBROUTINE LOCATE(XX,N,X,J)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END

      SUBROUTINE LOCATE8_SUB(XX,N,X,J)
      REAL*8 X,XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END

      INTEGER FUNCTION ILOCATE(XX,N,X)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      ILOCATE=JL
      RETURN
      END

	REAL FUNCTION XCLOSEST(X,XX,N,J)
	DIMENSION XX(N)
	integer ilocate
	! Returns closest value to X in array XX
	j = ilocate(xx,n,x)
	if (j == n) then
		return
	elseif (abs(x-xx(j)) < abs(x-xx(j+1))) then
		xclosest = xx(j)
	else
		xclosest = xx(j+1)
		j=j+1
	endif
	return
	end

      INTEGER FUNCTION ALOCATE(XX,N,X)
      ! Alphanumeric (adapted by GBA from Numerical Recipes SORT2)
      CHARACTER*(*) XX(N),X
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      ALOCATE=JL
      RETURN
      END
