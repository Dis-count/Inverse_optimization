\ Model Inv_facility
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C0 + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 + C11 + C12 + C13
   + C14 + C15 + C16 + C17 + C18 + C19 + C20 + C21 + C22 + C23
Subject To
 R0: C0 + C1 + C3 + C5 + C7 - C12 - C13 - C15 - C17 - C19 = -9
 R1: C0 + C3 + C4 + C5 - C12 - C15 - C16 - C17 >= 0
 R2: C0 + C1 + C3 + C4 + C8 - C12 - C13 - C15 - C16 - C20 >= -8
 R3: C0 + C2 + C3 + C4 + C11 - C12 - C14 - C15 - C16 - C23 >= -2
 R4: C0 + C1 + C3 + C5 + C7 - C12 - C13 - C15 - C17 - C19 >= -9
 R5: C0 + C1 + C3 + C7 + C8 - C12 - C13 - C15 - C19 - C20 >= -11
 R6: C0 + C1 + C2 + C3 + C7 + C11 - C12 - C13 - C14 - C15 - C19 - C23
   >= -11
 R7: C0 + C2 + C3 + C5 + C10 - C12 - C14 - C15 - C17 - C22 >= -9
 R8: C0 + C1 + C2 + C3 + C8 + C10 - C12 - C13 - C14 - C15 - C20 - C22
   >= -17
 R9: C0 + C2 + C3 + C10 + C11 - C12 - C14 - C15 - C22 - C23 >= -4
 R10: C0 + C1 + C4 + C5 + C6 - C12 - C13 - C16 - C17 - C18 >= 0
 R11: C0 + C1 + C4 + C6 + C8 - C12 - C13 - C16 - C18 - C20 >= -2
 R12: C0 + C1 + C2 + C4 + C6 + C11 - C12 - C13 - C14 - C16 - C18 - C23
   >= -2
 R13: C0 + C1 + C5 + C6 + C7 - C12 - C13 - C17 - C18 - C19 >= -3
 R14: C1 + C6 + C7 + C8 - C13 - C18 - C19 - C20 >= 0
 R15: C1 + C2 + C6 + C7 + C11 - C13 - C14 - C18 - C19 - C23 >= 0
 R16: C0 + C1 + C2 + C5 + C6 + C10 - C12 - C13 - C14 - C17 - C18 - C22
   >= -9
 R17: C1 + C2 + C6 + C8 + C10 - C13 - C14 - C18 - C20 - C22 >= -6
 R18: C1 + C2 + C6 + C10 + C11 - C13 - C14 - C18 - C22 - C23 >= 1
 R19: C0 + C2 + C4 + C5 + C9 - C12 - C14 - C16 - C17 - C21 >= -15
 R20: C0 + C1 + C2 + C4 + C8 + C9 - C12 - C13 - C14 - C16 - C20 - C21
   >= -23
 R21: C0 + C2 + C4 + C9 + C11 - C12 - C14 - C16 - C21 - C23 >= -10
 R22: C0 + C1 + C2 + C5 + C7 + C9 - C12 - C13 - C14 - C17 - C19 - C21
   >= -24
 R23: C1 + C2 + C7 + C8 + C9 - C13 - C14 - C19 - C20 - C21 >= -21
 R24: C1 + C2 + C7 + C9 + C11 - C13 - C14 - C19 - C21 - C23 >= -14
 R25: C0 + C2 + C5 + C9 + C10 - C12 - C14 - C17 - C21 - C22 >= -17
 R26: C1 + C2 + C8 + C9 + C10 - C13 - C14 - C20 - C21 - C22 >= -20
 R27: C2 + C9 + C10 + C11 - C14 - C21 - C22 - C23 >= -7
Bounds
End
