library("anticlust")

# test that function for greatest common denominator works correctly
x <- c(20, 20, 40, 40, 60, 60, 80, 80)
expect_equal(anticlust:::gcd_set(x), 20)

x <- 1:5
expect_equal(anticlust:::gcd_set(x), 1)

x <- c(9, 9, 9, 5)
expect_equal(anticlust:::gcd_set(x), 1)

x <- c(2, 8, 10, 10, 10, 10, 20, 20, 30, 30, 30, 1000, 1000, 10000000000)
expect_equal(anticlust:::gcd_set(x), 2)

x <- seq(3, 999, by = 3)
expect_equal(anticlust:::gcd_set(x), 3)

x <- c(1, seq(3, 999, by = 3))
expect_equal(anticlust:::gcd_set(x), 1)
