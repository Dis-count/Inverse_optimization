\ Model Model_sub
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0.2 C0 - 0.2 C1 - 0.8 C2 - 1.7 C3 + 1.8 C4 + 2.2 C5 + 2.8 C6 + 3.7 C7
Subject To
 R0: 0.4 C0 + 0.9 C2 - 0.4 C4 - 0.9 C6 >= 0
 R1: 0.4 C1 + 0.9 C3 - 0.4 C5 - 0.9 C7 >= 0
Bounds
 C0 <= 1000
 C1 <= 1000
 C2 <= 1000
 C3 <= 1000
 C4 <= 1000
 C5 <= 1000
 C6 <= 1000
 C7 <= 1000
End
