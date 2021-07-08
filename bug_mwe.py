#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimum working example (MWE) in the context of electrostatics demonstrating
erroneous solutions obtained on Dirichlet boundaries when using
`assemble_system` with a source term on an internal boundary.

Created on Thu Jul  8 09:22:38 2021

@author: Connor D. Pierce

This software is made available under the terms of the MIT License (a copy of
which is provided with this software).
"""

from fenics import *
from mshr import Rectangle, Circle, generate_mesh
import numpy as np

# Create mesh: square with circle in the middle
R = 0.25
domain      = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5))
inclusion   = Circle(Point(0, 0), R)
domain.set_subdomain(1, inclusion)
mesh        = generate_mesh(domain, 64)

# Define function space
V = FunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)

# Dirichlet boundary conditions for top and bottom surfaces. (Neumann
# conditions are applied automatically for left and right surfaces.)
u_T = Constant(0.5)
u_B = Constant(0.1)
def on_top(x, on_boundary):
    return near(x[1], 0.5) and on_boundary
def on_bottom(x, on_boundary):
    return near(x[1], -0.5) and on_boundary

# Dirichlet BCs on top and bottom surfaces
bc_T = DirichletBC(V, u_T, on_top)
bc_B = DirichletBC(V, u_B, on_bottom)

# Create a MeshFunction on mesh entities of dimension 1 to identify the
# internal boundary. A default value of 0 is assigned to all facets to start.
markers = MeshFunction('size_t', mesh, 1, 0)

# InternalBoundary defines a geometric test to identify facets on the internal
# boundary
class InternalBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(np.sqrt(x[0]**2 + x[1]**2),R,0.1*R)
internalBound = InternalBoundary()

# Mark all facets on the internal boundary as belonging to subdomain 1:
internalBound.mark(markers, 1)

# Create an integration measure over the internal boundary. Note the use of
# 'dS' and not 'ds'--lowercase 'ds' integrates only over external facets!
dS = Measure('dS', domain=mesh, subdomain_data=markers)#, subdomain_id=1)

# Define variational problem
eps1    = Constant(1)
eps2    = Constant(10)
epsilon = Expression('(x[0]*x[0]+x[1]*x[1]<R*R)?eps1:eps2',
                      degree=2,
                      R=R, eps1=eps1, eps2=eps2)
# epsilon = Constant(1)
sigma = Constant(1)  # Surface charge
# rho   = Expression('(x[0]*x[0]+x[1]*x[1]<R*R)?1:0',degree=2,R=R)
rho   = Constant(0)  # Volume charge
a       = dot(epsilon*grad(u), grad(v))*dx
L       = rho*v*dx + sigma*v("-")*dS(1)

# Solution Method 1: assemble system matrix and load vector using
# assemble_system, then solve the linear system
sol1     = Function(V)
A, b = assemble_system(a, L, bcs=[bc_T, bc_B]) 
solve(A, sol1.vector(), b)

# Solution Method 2: solve the linear variational problem using
# solve(a == L, ...)
sol2 = Function(V)
solve(a == L, sol2, [bc_T, bc_B])

# Solution Method 3: assemble system matrix and load vector using
# assemble, apply bcs, then solve the linear system
sol3 = Function(V)
A = assemble(a)
b = assemble(L)
for bc in [bc_T, bc_B]:
    bc.apply(A)
    bc.apply(b)
solve(A, sol3.vector(), b)

# Project the electric field to a continuous basis and plot
from matplotlib import pyplot as plt
plt.figure()
h0 = plot(sol1)
plt.colorbar(h0)
plt.title("Solution Method 1")
plt.figure()
h0 = plot(sol2)
plt.colorbar(h0)
plt.title("Solution Method 2")
plt.figure()
h0 = plot(sol3)
plt.colorbar(h0)
plt.title("Solution Method 3")
plt.show()