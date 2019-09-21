\ Model LB_Inv_UFL
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 + C11 + C12 + C13
   + C14 + C15 + C16 + C17 + C18 + C19 + C20 + C21 + C22 + C23 + C24 + C25
   + C26 + C27 + C28 + C29 + C30 + C31 + C32 + C33 + C34 + C35 + C36 + C37
   + C38 + C39
Subject To
 R0: C1 + C3 + C8 + C9 + C10 + C19 - C21 - C23 - C28 - C29 - C30 - C39
   = -64
 R1: C10 - C11 - C30 + C31 <= 1
 R2: - C17 + C19 + C37 - C39 <= 0
 R3: - C0 + C3 - C7 + C19 + C20 - C23 + C27 - C39 <= -11
 R4: - C2 + C3 - C15 + C19 + C22 - C23 + C35 - C39 <= 31
 R5: - C0 - C4 + C8 + C20 + C24 - C28 <= 18
 R6: - C0 - C5 + C9 + C20 + C25 - C29 <= 6
 R7: - C0 - C6 + C10 + C20 + C26 - C30 <= 13
 R8: - C2 + C8 - C12 + C22 - C28 + C32 <= 59
 R9: - C2 + C9 - C13 + C22 - C29 + C33 <= 47
 R10: - C2 + C10 - C14 + C22 - C30 + C34 <= 43
 R11: C3 - C19 - C23 + C39 <= -27
 R12: - C8 + C28 <= 0
 R13: C9 - C17 - C29 + C37 <= -9
 R14: C10 - C18 - C30 + C38 <= -2
 R15: C1 + C8 + C9 + C10 - C16 - C17 - C18 - C21 - C28 - C29 - C30 + C36
   + C37 + C38 <= -39
Bounds
End
