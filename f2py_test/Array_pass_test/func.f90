MODULE PassingArray
  implicit none

  contains
  REAL FUNCTION ArrManip(a1,a2,S)

    INTEGER :: a1,a2
    REAL    :: S(a1,a2)

    ArrManip = SUM(S)

  END FUNCTION ArrManip
END MODULE PassingArray
