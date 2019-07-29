\ Model I_UFL
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  9 C0 + 2 C1 + 16 C2 + 9 C3 + C4 + 5 C5 + 2 C6 + 3 C7 + 2 C8 + 8 C9 + 4 C10
   + 4 C11 + C12 + 2 C13 + 5 C14 + 7 C15 + 8 C16 + C17 + 4 C19 + 8 C20
   + 4 C21 + 7 C23
Subject To
 R0: C0 - C4 >= 0
 R1: C0 - C5 >= 0
 R2: C0 - C6 >= 0
 R3: C0 - C7 >= 0
 R4: C0 - C8 >= 0
 R5: C1 - C9 >= 0
 R6: C1 - C10 >= 0
 R7: C1 - C11 >= 0
 R8: C1 - C12 >= 0
 R9: C1 - C13 >= 0
 R10: C2 - C14 >= 0
 R11: C2 - C15 >= 0
 R12: C2 - C16 >= 0
 R13: C2 - C17 >= 0
 R14: C2 - C18 >= 0
 R15: C3 - C19 >= 0
 R16: C3 - C20 >= 0
 R17: C3 - C21 >= 0
 R18: C3 - C22 >= 0
 R19: C3 - C23 >= 0
 R20: C4 + C9 + C14 + C19 = 1
 R21: C5 + C10 + C15 + C20 = 1
 R22: C6 + C11 + C16 + C21 = 1
 R23: C7 + C12 + C17 + C22 = 1
 R24: C8 + C13 + C18 + C23 = 1
Bounds
Binaries
 C0 C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20
 C21 C22 C23
End
