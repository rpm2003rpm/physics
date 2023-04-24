## @package timeIndependentWaveFunction
#  Program based on the appendix G of the book Quantum Physics and Atoms.
# 
#  @author  Rodrigo Pedroso Mendes
#  @version V1.0
#  @date    31/03/23 02:38:30
#
#  #LICENSE# 
#    
#  Copyright (c) 2023 Rodrigo Pedroso Mendes
#
#  Permission is hereby granted, free of charge, to any  person   obtaining  a 
#  copy of this software and associated  documentation files (the "Software"), 
#  to deal in the Software without restriction, including  without  limitation 
#  the rights to use, copy, modify,  merge,  publish,  distribute, sublicense, 
#  and/or sell copies of the Software, and  to  permit  persons  to  whom  the 
#  Software is furnished to do so, subject to the following conditions:        
#   
#  The above copyright notice and this permission notice shall be included  in 
#  all copies or substantial portions of the Software.                         
#   
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE  WARRANTIES  OF  MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
#  AUTHORS OR COPYRIGHT HOLDERS BE  LIABLE FOR ANY  CLAIM,  DAMAGES  OR  OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT  OR  OTHERWISE,  ARISING 
#  FROM, OUT OF OR IN CONNECTION  WITH  THE  SOFTWARE  OR  THE  USE  OR  OTHER  
#  DEALINGS IN THE SOFTWARE. 
#    
#  Solve the unidimentional time independent wave function for a given 
#  potencial curve.
#
################################################################################
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.constants import *

#-------------------------------------------------------------------------------
# Wave function
#-------------------------------------------------------------------------------
def waveFunc(Pot, totE, dx, m, xmax):
    x0    = 1e-18
    psi0 = 1
    psi  = [psi0]
    x    = [x0]
    ndD  = dx*m/(hbar**2)*(Pot(x0) - totE)*psi0
    while x0 < xmax:
        psi0 = psi0 + ndD*dx
        x0   = x0   + dx
        ndD  = ndD  + dx*2*m/(hbar**2)*(Pot(x0) - totE)*psi0
        x.append(x0)
        psi.append(psi0)
    return x, psi


#-------------------------------------------------------------------------------
# Potential
#-------------------------------------------------------------------------------
def pot(x):
    if x > 0.5 or x < -0.5:
        return 1.0
    elif x == 0.5 or x == -0.5:
        return 0.5
    else:
        return 0.0


def potHydrogen(x): 
    return -1/(4*pi*epsilon_0)*e**2/x

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    #m = 32*hbar**2
    #m = 100*m
    #func = lambda totE : waveFunc(pot, totE, 1e-5, m, 2)[1][-1]
    #totE = fsolve(func, 0.09, xtol=1e-8)
    #totE = fsolve(func, 0.4, xtol=1e-8)
    #totE = fsolve(func, 0.8, xtol=1e-8)
    #print(totE)
    #x, psi = waveFunc(pot, totE[0], 1e-5, m, 1.5)
    #plt.plot(x, psi)
    #plt.ylabel('Psi')
    #plt.xlabel('x')
    #plt.show()      
    #Calculate hydrogen radius for the first energy level
    x, psi = waveFunc(potHydrogen, -13.6*e, 1e-14, m_e, 3.5e-10)
    plt.plot(x, psi)
    plt.ylabel('Psi')
    plt.xlabel('x')
    plt.show()      
    #Calculate hydrogen radius for the second energy level
    x, psi = waveFunc(potHydrogen, -3.4*e, 1e-14, m_e, 10e-10)
    plt.plot(x, psi)
    plt.ylabel('Psi')
    plt.xlabel('x')
    plt.show()      
    #Calculate hydrogen radius for the third energy level
    x, psi = waveFunc(potHydrogen, -1.51*e, 1e-14, m_e, 1.6e-9)
    plt.plot(x, psi)
    plt.ylabel('Psi')
    plt.xlabel('x')
    plt.show()      

 

