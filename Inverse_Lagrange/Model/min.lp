\ Model Min_sub
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + 24 C1
Subject To
 R0: - C0 + 6 C1 >= 5
 R1: C0 + 4 C1 >= 4
Bounds
 C0 <= 1
 C1 <= 1
End
