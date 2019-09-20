\ Model LB_Inv_UFL
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 + C11 + C12 + C13
   + C14 + C15 + C16 + C17 + C18 + C19 + C20 + C21 + C22 + C23
Subject To
 R0: C0 + C1 + C3 + C5 + C7 - C12 - C13 - C15 - C17 - C19 = -9
 R1: <= 0
 R2: - C2 + C4 - C11 + C14 - C16 + C23 <= 6
 R3: C0 - C3 - C5 - C12 + C15 + C17 <= -9
Bounds
End
