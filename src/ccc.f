c ID: ccc.f, last updated 2025-06-23, F.Osorio

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE rho1_ustat(x, y, n, p1, p2)
      INTEGER          n
      DOUBLE PRECISION x(*), y(*), p1(*), p2(*)
c
c     contructs kernels of the U-statistics
c
c     parameters:
c     x       (input) DOUBLE PRECISION array, dimension (n)
c             observations of 1st measurement instrument
c     y       (input) DOUBLE PRECISION array, dimension (n)
c             observations of 2nd measurement instrument
c     n       (input) INTEGER
c             length of vectors x and y. n > 0.
c     p1      (output) DOUBLE PRECISION array, dimension (n)
c             kernel of the 1st component of U-statistic
c     p2      (output) DOUBLE PRECISION array, dimension (n)
c             kernel of the 2nd component of U-statistic
c
c     .. local scalars ..
      INTEGER          i, j
      DOUBLE PRECISION acc1, acc2, d1, d2, d3, d4
c
c     quick return if possible
c
      if (n .LE. 0) return
c
c     computing kernel of U-statistics
c
      do i = 1, n
        acc1 = 0.d0
        acc2 = 0.d0
        do j = 1, n
          if (i .NE. j) then
            d1 = dabs(x(i) - y(i))
            d2 = dabs(x(j) - y(j))
            d3 = dabs(x(i) - y(j))
            d4 = dabs(x(j) - y(i))
            acc1 = acc1 + (d1 + d2) / 2
            acc2 = acc2 + (d3 + d4) / 2
          end if
        end do
        p1(i) = acc1 / (n - 1)
        p2(i) = acc2 / (n - 1)
      end do
c
      return
      END
