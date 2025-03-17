c ID: CALGO478.f, last updated 2022-10-14, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE L1BR(z, y, n, p, n2, p2, coef, resid, minimum, iter,
     *                tol, rank, info, work)
      INTEGER          n, p, n2, p2, iter, rank, info
      DOUBLE PRECISION minimum, tol
      INTEGER          work(*)
      DOUBLE PRECISION z(n2,*), y(*), coef(*), resid(*)
c
c     wrapper to 'L1'
c
      call L1(n, p, n2, p2, z, y, tol, coef, resid, work)
c
      minimum = z(n+1,p+1)
      rank = INT(z(n+1,p+2))
      info = INT(z(n+2,p+1))
      iter = INT(z(n+2,p+2))
c
c     end of L1BR
c
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE L1(m, n, m2, n2, a, b, toler, x, e, s)
      INTEGER          m, n, m2, n2, s(*)
      DOUBLE PRECISION a(m2,*), b(*), x(*), e(*), toler
c
c     L1FIT uses a modification of the simplex method of linear programming
c     to calculate an L1 solution to an over-determined system of linear
c     equations.
c     Algorithm 478: Commun. ACM 17, 1974, 319-320. doi: 10.1145/355616.361024
c
c     parameters:
c     m       (input) INTEGER
c             The number of equations.
c     n       (input) INTEGER
c             The number of unknowns (m >= n).
c     m2      (input) INTEGER
c             set equal to m + 2 for adjustable dimensions.
c     n2      (input) INTEGER
c             set equal to n + 2 for adjustable dimensions.
c     a       (input) DOUBLE PRECISION array, dimension(m2, n2)
c             on entry, the coefficients of the matrix must be
c             stored in the first m rows and n columns of a.
c             these values are destroyed by the subroutine.
c             on exit from the subroutine, the array a contains
c             the following information.
c             a(m+1,n+1) the minimum sum of the absolute values of
c             the residuals.
c             a(m+1,n+2) the rank of the matrix of coefficients.
c             a(m+2,n+1) exit code with values.
c             = 0: optimal solution which is probably non-unique.
c             = 1: unique optimal solution.
c             = 2: calculations terminated prematurely due to
c             rounding errors.
c             a(m+2,n+2) number of simplex iterations performed.
c     b       (input/output) DOUBLE PRECISION array, dimension(m)
c             on entry, b must contain the right hand side of the
c             equations. these values are destroyed by the subroutine.
c     toler   (input) DOUBLE PRECISION
c             a small positive tolerance. empirical evidence suggests
c             toler = 10**(-d*2/3) where d represents the number of
c             decimal digits of accuracy avalable.
c     x       (output) DOUBLE PRECISION array, dimension(n)
c             on exit, this array contains a solution to the L1 problem.
c     e       (output) DOUBLE PRECISION array, dimension(m)
c             on exit, this array contains the residuals in the equations.
c     s       (workspace) INTEGER array, dimension(m).
c
c     .. parameters ..
      DOUBLE PRECISION big
      PARAMETER       (big = 1.d75)
c     .. local scalars ..
      DOUBLE PRECISION d, sum, pivot
      INTEGER          i, j, k, kount, kr, kl, in, out
      LOGICAL          stage, test
c     .. intrinsic functions ..
      DOUBLE PRECISION min, max
c
c     executable statements
c
c     initialization
c
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
c
c     compute the marginal costs
c
      do 60 j=1,n1
        sum = 0.d0
        do 50 i=1,m
          sum = sum + a(i,j)
   50   continue
        a(m1,j) = sum
   60 continue
c
c     STAGE I.
c     determine the vector to enter the basis
c
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
c
c     determine the vector to leave the basis
c
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
c
c     check for linear dependence in stage I
c
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
c
c     pivot on a(out,in)
c
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
c
c     interchange rows in stage I
c
      kl = kl + 1
      do 250 j=kr,n2
        d = a(out,j)
        a(out,j) = a(kount,j)
        a(kount,j) = d
  250 continue
  260 if (kount+kr.ne.n1) go to 70
c
c     STAGE II.
c
      stage = .false.
c
c     determine the vector to enter the basis
c
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
c
c     prepare output
c
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
c
c     end of L1
c
      return
      END
