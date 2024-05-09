*ROWS:         2
*COLUMNS:      4
*INTEGER:      0
*NONZERO:      6
*
*
NAME          relaxed_problem_intersection_cut.mps
ROWS
 N  OBJ
 E  mb_ineq
 E  mb_ineq_2
COLUMNS
    x[1]      OBJ       -1.00000000000000
    x[1]      mb_ineq   -1.00000000000000
    x[2]      OBJ       -2.00000000000000
    x[2]      mb_ineq   1.00000000000000
    x[2]      mb_ineq_2 -1.00000000000000
    x[4]      OBJ       -1.00000000000000
    x[4]      mb_ineq   1.00000000000000
    x[4]      mb_ineq_2 -1.00000000000000
    x[5]      OBJ       -1.00000000000000
    x[5]      mb_ineq_2 1.00000000000000
RHS
BOUNDS
 UP BND       x[1]      10.0000000000000
 LO BND       x[2]      -10.0000000000000
 UP BND       x[2]      10.0000000000000
 LO BND       x[4]      -10.0000000000000
 UP BND       x[4]      10.0000000000000
 UP BND       x[5]      10.0000000000000
ENDATA
