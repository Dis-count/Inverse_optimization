\ Model Model_sub
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + C1 + C2 + C3 + C4 + C5
Subject To
 R0: 2 C0 + 1.5 C1 <= 1
 R1: 2 C2 + 1.5 C3 <= 24
 R2: 2 C4 + 1.5 C5 <= 9
 R3: - 0.8 C2 - 0.2 C4 <= -5
 R4: - 0.8 C3 - 0.2 C5 <= -4
 R5: = 0
 R6: 1.6 C2 + 1.2 C3 = 19.2
 R7: 0.4 C4 + 0.3 C5 = 1.8
Bounds
 C0 >= -1
 C1 >= 1
 C2 >= 6
 C3 >= 4
 C4 >= 1
 C5 >= 4
End
