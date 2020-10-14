* ID: l1fit_BR.f, last updated 2020-08-16, F.Osorio

************************************************************************
      subroutine l1fit(z, y, n, p, n2, p2, coef, resid, minimum, iter,
     *                 tol, rank, info, work)
*
*     wrapper to 'l1_BR'
*
      integer          n, p, n2, p2, iter, rank, info
      double precision minimum, tol
*
      double precision z(n2,*), y(*), coef(*), resid(*)
      integer          work(*)
*
      call L1(n, p, n2, p2, z, y, tol, coef, resid, work)
*
      minimum = z(n+1,p+1)
      rank = int(z(n+1,p+2))
      info = int(z(n+2,p+1))
      iter = int(z(n+2,p+2))

*
*     End of L1FIT
*
      end

************************************************************************
      subroutine l1(m, n, m2, n2, a, b, toler, x, e, s)
*
*     ALGORITHM 478 Collected Algorithms from ACM
*     Comm. ACM, Vol. 17, No. 06, p. 319.
*
*     scalar arguments
      integer          m, n, m2, n2
      double precision toler
*
*     array arguments
      double precision a(m2,*), b(*), x(*), e(*)
      integer          s(*)
*
*     Purpose:
*
*     l1_BR uses a modification of the simplex method of linear programming
*     to calculate an L1 solution to an over-determined system of linear
*     equations.
*
*     Arguments:
*
*     m     (input) INTEGER
*           The number of equations.
*     n     (input) INTEGER
*           The number of unknowns (m.ge.n).
*     m2    (input) INTEGER
*           set equal to m+2 for adjustable dimensions.
*     n2    (input) INTEGER
*           set equal to n+2 for adjustable dimensions.
*     a     (input) DOUBLE PRECISION array, dimension (m2,n2)
*           on entry, the coefficients of the matrix must be
*           stored in the first m rows and n columns of a.
*           these values are destroyed by the subroutine.
*           on exit from the subroutine, the array a contains
*           the following information.
*           a(m+1,n+1) the minimum sum of the absolute values
*                      of the residuals.
*           a(m+1,n+2) the rank of the matrix of coefficients.
*           a(m+2,n+1) exit code with values.
*                      0 - optimal solution which is probably
*                          non-unique.
*                      1 - unique optimal solution.
*                      2 - calculations terminated prematurely
*                          due to rounding errors.
*           a(m+2,n+2) number of simplex iterations performed.
*     b     (input/output) DOUBLE PRECISION array, dimension (m)
*           on entry, b must contain the right hand side of the
*           equations. these values are destroyed by the subroutine.
*     toler (input) DOUBLE PRECISION
*           a small positive tolerance. empirical evidence suggests
*           toler=10**(-d*2/3) where d represents the number of
*           decimal digits of accuracy avalable (see description).
*     x     (output) DOUBLE PRECISION array, dimension (n)
*           on exit, this array contains a solution to the L1 problem.
*     e     (output) DOUBLE PRECISION array, dimension (m)
*           on exit, this array contains the residuals in the equations.
*     s     (workspace) INTEGER array, dimension(m).
*
*     parameters
      double precision big
      parameter       (big = 1.d75)
*
*     local scalars
      double precision d, sum, pivot
      integer          i, j, k, kount, kr, kl, in, out
      logical          stage, test
*
*     intrinsic functions
      double precision min, max
*
*     Executable Statements
*
*     initialization
*
      m1 = m + 1
      n1 = n + 1
      do 10 j=1,n
        a(m2,j) = j
        x(j) = 0.d0
   10 continue
      do 40 i=1,m
        a(i,n2) = n + i
        a(i,n1) = b(i)
        if (b(i).ge.0.d0) go to 30
        do 20 j=1,n2
          a(i,j) = -a(i,j)
   20   continue
   30   e(i) = 0.d0
   40 continue
*
*     compute the marginal costs
*
      do 60 j=1,n1
        sum = 0.d0
        do 50 i=1,m
          sum = sum + a(i,j)
   50   continue
        a(m1,j) = sum
   60 continue
*
*     STAGE I.
*     determine the vector to enter the basis
*
      stage = .true.
      in = 0
      out = 0
      kount = 0
      kr = 1
      kl = 1
   70 max = -1.
      do 80 j=kr,n
        if (abs(a(m2,j)).gt.n) go to 80
        d = abs(a(m1,j))
        if (d.le.max) go to 80
        max = d
        in = j
   80 continue
      if (a(m1,in).ge.0.d0) go to 100
      do 90 i=1,m2
        a(i,in) = -a(i,in)
   90 continue
*
*     determine the vector to leave the basis
*
  100 k = 0
      do 110 i=kl,m
        d = a(i,in)
        if (d.le.toler) go to 110
        k = k + 1
        b(k) = a(i,n1)/d
        s(k) = i
        test = .true.
  110 continue
  120 if (k.gt.0) go to 130
      test = .false.
      go to 150
  130 min = big
      do 140 i=1,k
        if (b(i).ge.min) go to 140
        j = i
        min = b(i)
        out = s(i)
  140 continue
      b(j) = b(k)
      s(j) = s(k)
      k = k - 1
*
*     check for linear dependence in stage I
*
  150 if (test .or. .not.stage) go to 170
      do 160 i=1,m2
        d = a(i,kr)
        a(i,kr) = a(i,in)
        a(i,in) = d
  160 continue
      kr = kr + 1
      go to 260
  170 if (test) go to 180
      a(m2,n1) = 2.d0
      go to 350
  180 pivot = a(out,in)
      if (a(m1,in)-pivot-pivot.le.toler) go to 200
      do 190 j=kr,n1
        d = a(out,j)
        a(m1,j) = a(m1,j) - d - d
        a(out,j) = -d
  190 continue
      a(out,n2) = -a(out,n2)
      go to 120
*
*     pivot on a(out,in)
*
  200 do 210 j=kr,n1
        if (j.eq.in) go to 210
        a(out,j) = a(out,j)/pivot
  210 continue
      do 230 i=1,m1
        if (i.eq.out) go to 230
        d = a(i,in)
        do 220 j=kr,n1
          if (j.eq.in) go to 220
          a(i,j) = a(i,j) - d*a(out,j)
  220   continue
  230 continue
      do 240 i=1,m1
        if (i.eq.out) go to 240
        a(i,in) = -a(i,in)/pivot
  240 continue
      a(out,in) = 1./pivot
      d = a(out,n2)
      a(out,n2) = a(m2,in)
      a(m2,in) = d
      kount = kount + 1
      if (.not.stage) go to 270
*
*     interchange rows in stage I
*
      kl = kl + 1
      do 250 j=kr,n2
        d = a(out,j)
        a(out,j) = a(kount,j)
        a(kount,j) = d
  250 continue
  260 if (kount+kr.ne.n1) go to 70
*
*     STAGE II.
*
      stage = .false.
*
*     determine the vector to enter the basis
*
  270 max = -big
      do 290 j=kr,n
        d = a(m1,j)
        if (d.ge.0.d0) go to 280
        if (d.gt.(-2.d0)) go to 290
        d = -d - 2.d0
  280   if (d.le.max) go to 290
        max = d
        in = j
  290 continue
      if (max.le.toler) go to 310
      if (a(m1,in).gt.0.d0) go to 100
      do 300 i=1,m2
        a(i,in) = -a(i,in)
  300 continue
      a(m1,in) = a(m1,in) - 2.d0
      go to 100
*
*     prepare output
*
  310 l = kl - 1
      do 330 i=1,l
        if (a(i,n1).ge.0.d0) go to 330
        do 320 j=kr,n2
          a(i,j) = -a(i,j)
  320   continue
  330 continue
      a(m2,n1) = 0.d0
      if (kr.ne.1) go to 350
      do 340 j=1,n
        d = abs(a(m1,j))
        if (d.le.toler .or. 2.d0-d.le.toler) go to 350
  340 continue
      a(m2,n1) = 1.d0
  350 do 380 i=1,m
        k = int(a(i,n2))
        d = a(i,n1)
        if (k.gt.0) go to 360
        k = -k
        d = -d
  360   if (i.ge.kl) go to 370
        x(k) = d
        go to 380
  370   k = k - n
        e(k) = d
  380 continue
      a(m2,n2) = kount
      a(m1,n2) = n1 - kr
      sum = 0.d0
      do 390 i=kl,m
        sum = sum + a(i,n1)
  390 continue
      a(m1,n1) = sum
*
      return
*
*     End of L1
*
      end
