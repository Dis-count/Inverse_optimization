\ Model UFL
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C33 + C34 + C35 + C36 + C37 + C38 + C39 + C40 + C41 + C42 + C43 + C44
   + C45 + C46 + C47 + C48 + C49 + C50 + C51 + C52 + C53 + C54 + C55 + C56
Subject To
 R0: C0 + C1 + C2 >= 28
 R1: C3 + C4 + C5 - C12 = 0
 R2: C6 + C7 + C8 - C13 = 0
 R3: C9 + C10 + C11 - C14 = 0
 R4: C0 - C3 + C15 - C24 = 0
 R5: C1 - C4 + C16 - C25 = 0
 R6: C2 - C5 + C17 - C26 = 0
 R7: C0 - C6 + C18 - C27 = 0
 R8: C1 - C7 + C19 - C28 = 0
 R9: C2 - C8 + C20 - C29 = 0
 R10: C0 - C9 + C21 - C30 = 0
 R11: C1 - C10 + C22 - C31 = 0
 R12: C2 - C11 + C23 - C32 = 0
 R13: C12 + C13 + C24 + C26 + C28 = 28
 R14: C12 - C33 + C36 = 5
 R15: C13 - C34 + C37 = 6
 R16: C14 - C35 + C38 = 7
 R17: C24 - C39 + C48 = 11
 R18: C25 - C40 + C49 = 4
 R19: C26 - C41 + C50 = 8
 R20: C27 - C42 + C51 = 5
 R21: C28 - C43 + C52 = 7
 R22: C29 - C44 + C53 = 10
 R23: C30 - C45 + C54 = 19
 R24: C31 - C46 + C55 = 6
 R25: C32 - C47 + C56 = 3
Bounds
End
