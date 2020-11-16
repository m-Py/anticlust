## Start R with valgrind:
# R -d valgrind --vanilla

library("anticlust")

anticlusters <- anticlustering(
  iris[, -5],
  K = 3,
  objective = "variance",
  method = "exchange",
  categories = iris[, 5]
)

anticlusters <- anticlustering(
  iris[, -5],
  K = 3,
  objective = "variance",
  method = "exchange"
)

q()

## Valgrind Summary from version 0.5.3 (no C code)

# ==16236== 
# ==16236== HEAP SUMMARY:
# ==16236==     in use at exit: 194,070,518 bytes in 41,252 blocks
# ==16236==   total heap usage: 278,822 allocs, 237,570 frees, 663,645,940 bytes allocated
# ==16236== 
# ==16236== LEAK SUMMARY:
# ==16236==    definitely lost: 0 bytes in 0 blocks
# ==16236==    indirectly lost: 0 bytes in 0 blocks
# ==16236==      possibly lost: 0 bytes in 0 blocks
# ==16236==    still reachable: 194,070,518 bytes in 41,252 blocks
# ==16236==                       of which reachable via heuristic:
# ==16236==                         newarray           : 4,264 bytes in 1 blocks
# ==16236==         suppressed: 0 bytes in 0 blocks
# ==16236== Rerun with --leak-check=full to see details of leaked memory
# ==16236== 
# ==16236== For lists of detected and suppressed errors, rerun with: -s
# ==16236== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

## Valgrind Summary from version > 0.5.5 (including C code)

# ==16463== 
# ==16463== HEAP SUMMARY:
# ==16463==     in use at exit: 176,690,322 bytes in 34,744 blocks
# ==16463==   total heap usage: 66,483 allocs, 31,739 frees, 282,474,156 bytes allocated
# ==16463== 
# ==16463== LEAK SUMMARY:
# ==16463==    definitely lost: 0 bytes in 0 blocks
# ==16463==    indirectly lost: 0 bytes in 0 blocks
# ==16463==      possibly lost: 0 bytes in 0 blocks
# ==16463==    still reachable: 176,690,322 bytes in 34,744 blocks
# ==16463==                       of which reachable via heuristic:
# ==16463==                         newarray           : 4,264 bytes in 1 blocks
# ==16463==         suppressed: 0 bytes in 0 blocks
# ==16463== Rerun with --leak-check=full to see details of leaked memory
# ==16463== 
# ==16463== For lists of detected and suppressed errors, rerun with: -s
# ==16463== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)


## (not more "still reachable" under C code?)
