\ SCIP STATISTICS
\   Problem name     : 
\   Variables        : 27 (10 binary, 0 integer, 0 implicit integer, 17 continuous)
\   Constraints      : 43
Maximize
 Obj: +1 x10 +1 x11 +1 x12 +1 x14 +1 x16 +1 x18 +1 x20
Subject to
 c0: +1 x12 -1 x14 -1 x18 -1 x10 = +0
 c1: +1 x14 -1 x16 = +0
 c2: +1 x16 +1 x18 -1 x20 +1 x11 = +0
 c3: +1 x10 -1 x11 = +0
 c4: -1 x12 <= +0
 c5: -1 x14 <= +10
 c6: -1 x16 <= +10
 c7: -1 x18 <= +10
 c8: -1 x20 <= +0
 c9: -1 x10 <= +0
 c10: -1 x11 <= +0
 c11: +1 x12 <= +20
 c12: +1 x14 <= +30
 c13: +1 x16 <= +30
 c14: +1 x18 <= +30
 c15: +1 x20 <= +20
 c16: +1 x10 <= +10
 c17: +1 x11 <= +10
 c18: +1 x0 +1 x5 = +1
 c20: x5 = 1 ->  +1 x14 <= -1e-07
 c21: +1 x1 +1 x6 = +1
 c23: x6 = 1 ->  +1 x16 <= -1e-07
 c24: +1 x2 +1 x7 = +1
 c26: x7 = 1 ->  +1 x18 <= -1e-07
 c27: +1 x3 +1 x8 = +1
 c29: x8 = 1 ->  +1 x10 <= -1e-07
 c30: +1 x4 +1 x9 = +1
 c32: x9 = 1 ->  +1 x11 <= -1e-07
 c34: +1 x14 +1 x22 >= +1e-07
 c36: +1 x16 +1 x23 >= +1e-07
 c38: +1 x18 +1 x24 >= +1e-07
 c40: +1 x10 +1 x25 >= +1e-07
 c42: +1 x11 +1 x26 >= +1e-07
Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
 0 <= x8 <= 1
 0 <= x9 <= 1
 x10 free
 x11 free
 x12 free
 x14 free
 x16 free
 x18 free
 x20 free
 x22 free
 x23 free
 x24 free
 x25 free
 x26 free
Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
SOS
 c33: S1:: x22:0.4 x0:0.6
 c35: S1:: x23:0.4 x1:0.6
 c37: S1:: x24:0.4 x2:0.6
 c39: S1:: x25:0.4 x3:0.6
 c41: S1:: x26:0.4 x4:0.6
End
