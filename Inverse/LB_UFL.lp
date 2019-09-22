\ Model LB_Inv_UFL
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 + C11 + C12 + C13
   + C14 + C15 + C16 + C17 + C18 + C19 + C20 + C21 + C22 + C23
Subject To
 R0: C0 + C1 + C3 + C5 + C7 - C12 - C13 - C15 - C17 - C19 = -9
 R1: C3 - C4 - C6 + C7 - C15 + C16 + C18 - C19 <= -9
 R2: - C4 + C5 + C7 - C8 + C16 - C17 - C19 + C20 <= -1
 R3: C1 - C2 + C7 - C10 - C13 + C14 - C19 + C22 <= 0
 R4: - C2 + C3 - C9 + C14 - C15 + C21 <= 15
 R5: - C2 + C5 - C11 + C14 - C17 + C23 <= 2
 R6: C1 - C4 + C7 - C13 + C16 - C19 <= -9
 R7: C3 - C6 - C15 + C18 <= -6
 R8: C5 - C8 - C17 + C20 <= 2
 R9: C0 + C3 + C5 - C6 - C8 - C12 - C15 - C17 + C18 + C20 <= -9
Bounds
End
