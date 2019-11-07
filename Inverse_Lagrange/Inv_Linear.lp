\ Model Inv_Linear_Production
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + C1 + C2 + C3 + C4 + C5 + C6 + C7
Subject To
 R0: 2 C0 + 2.5 C1 - 2 C4 - 2.5 C5 = 0.5
 R1: 2 C2 + 2.5 C3 - 2 C6 - 2.5 C7 = 2
Bounds
End
