\ Model Model_sub
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - C0 - 2 C1 - C2 - 2 C3 + 3 C4 + 4 C5 + 3 C6 + 4 C7
Subject To
 R0: C0 + C2 - C4 - C6 >= -4
 R1: C1 + C3 - C5 - C7 >= 19
Bounds
 C0 <= 10
 C1 <= 10
 C2 <= 10
 C3 <= 10
 C4 <= 10
 C5 <= 10
 C6 <= 10
 C7 <= 10
End
