#!/usr/bin/env python3.7

# Copyright 2020, Gurobi Optimization, LLC

# This example formulates and solves the following simple bilinear model:
#  maximize    x
#  subject to  x + y + z <= 10
#              x * y <= 2         (bilinear inequality)
#              x * z + y * z = 1  (bilinear equality)
#              x, y, z non-negative (x integral in second version)

import gurobipy as gp
from gurobipy import GRB

# Create a new model
m = gp.Model("bilinear")

# Create variables
x1 = m.addVar(name="x1")
x2 = m.addVar(name="x2")
x3 = m.addVar(name="x3")
x4 = m.addVar(name="x4")
x5 = m.addVar(name="x5")
x6 = m.addVar(name="x6")
x7 = m.addVar(name="x7")
x8 = m.addVar(name="x8")
y1 = m.addVar(name="y1")
y2 = m.addVar(name="y2")


# Set objective: maximize x
m.setObjective(1.0*x1+1.0*x2+1.0*x3+1.0*x4+1.0*x5+1.0*x6+1.0*x7+1.0*x8, GRB.MINIMIZE)

# Add linear constraint: x + y + z <= 10
m.addConstr(-2*x1 -2.5*x2 + 2*x5 + 2.5*x6 >= -1, "c0")
m.addConstr(-2*x3 -2.5*x4 + 2*x7 + 2.5*x8 >= -4, "c1")
m.addConstr(y1 + 24* y2 == 20, "c2")

# Add bilinear inequality constraint: x * y <= 2
m.addConstr(x1*y1 + x3*y2 - x5*y1 - x7*y2 - y1 + 6*y2 >= 5, "bilinear0")
m.addConstr(x2*y1 + x4*y2 - x6*y1 - x8*y2 + y1 + 4*y2 >= 4, "bilinear1")

# Add bilinear equality constraint: x * z + y * z == 1

# First optimize() call will fail - need to set NonConvex to 2
try:
    m.optimize()
except gp.GurobiError:
    print("Optimize failed due to non-convexity")

# Solve bilinear model
m.params.NonConvex = 2
m.optimize()

# m.printAttr('x1')

# Constrain 'x' to be integral and solve again
# x.vType = GRB.INTEGER
# m.optimize()
#
# m.printAttr('x')
