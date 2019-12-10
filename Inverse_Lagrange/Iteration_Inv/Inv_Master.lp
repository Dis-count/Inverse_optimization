\ Model Inv_Master
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + C1 + C2 + C3 + C4 + C5 + C6 + C7
Subject To
 R0: 0.6666666666666666 C0 + 4 C1 - 0.6666666666666666 C4 - 4 C5 = 5
 R1: 0.6666666666666666 C2 + 4 C3 - 0.6666666666666666 C6 - 4 C7
   = -2.666666666666666
Bounds
End
