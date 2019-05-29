from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import pandas as pd
import numpy as np
import time
import sympy
from sympy import *
sampler = EmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi', token='', solver='DW_2000Q_2_1'))
orig_stdout = sys.stdout
f = open('mapColoring4ColorsDwaveResults.txt', 'w')
sys.stdout = f
print("\n# MAP COLOURING PROBLEM WITH 4 COLOURS ON D-WAVE #\n")

print("\nSymbolic Computing\n")

#### Symbolic Computing

h = 0.0000005 #small number

n = 4 #four colors

def alpha(n,h):
  return(h+h)

def beta(n,alfa,h):
  return(((n^3+n^2+1)*alfa)+h)

A = alpha(n,h)
B = beta(n,alpha(n,h),h)

exp1 = ''
for i in range(1,n):
  for k in range(1,n):
    for kl in range(1,n):
      if (k!=kl):
        exp1 = exp1+str(B)+"*"+str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(i)+str(kl)))+"+"

print("1st expression:")
print(exp1)

exp2 = ''
for i in range(1,n):
  for k in range(1,n):
      exp2 = exp2+str(A)+"*(1-"+str(var("vc"+str(i)+str(k)))+")"+"+"

print("2nd expression:")
print(exp2)

exp3 = ''
for i in range(1,n):
  for il in range(1,n):
    if(il!=i):
      for k in range(1,n):
        exp3 = exp3+str(A)+"*"+str(var("vc"+str(i)+str(k)))+"*"+str(var("vc"+str(i)+str(il)))+"*"+str(var("neigh"+str(i)+str(il)))+"+"

print("3rd expression:")
print(exp3)

exp4 = ''
for i in range(1,n):
  for k in range(1,n):
    exp4 = exp4+str(A)+"*"+str(var("vc"+str(i)+str(k)))+"*(1-"+str(var("c"+str(k)))+")"+"+"

print("4th expression:")
print(exp4)

exp5 = ''
for k in range(1,n):
    exp5 = exp5+"1*"+str(var("c"+str(k)))+"+"

print("5th expression:")
print(exp5)

#### FINAL EXPRESSION

exp = exp1 + exp2 + exp3 + exp4 + exp5
exp = exp[:-1]
print(exp)

print("Symbolic expression:")
print(exp)

# symbolic expression
exp = sympy.simplify(exp)
print("Symbolic expression (calculated):")
print(exp)

# with numeric constantes
print("Symbolic expression (calculated, replace ^2 since is a binary problem):")
exp = str(exp)
exp = exp.replace('**2', '')
exp = sympy.simplify(exp)
exp = str(exp)
print(exp)

# split and replace

n = len(exp.split(" + ")) #length(unlist(strsplit(exp, " [\\s+] ")))
exp = exp.split(" + ")
print(exp)
expn = ''
for i in range(1,n):
   e = exp[i]
   if (e.count("*") == 3):
    v = e.split("*")
    if (float(v[0]) > 0):
        e = str(var(v[0]))+"*(w*("+str(var(v[1]))+"+"+str(var(v[2]))+"+"+str(var(v[3]))+"-1)+("+str(var(v[1]))+"*"+str(var(v[2]))+"+"+str(var(v[2]))+"*"+str(var(v[3]))+"+"+str(var(v[3]))+"*"+str(var(v[3]))+")-("+str(var(v[3]))+"+"+str(var(v[2]))+"+"+str(var(v[3]))+")+1)"
    else:
        e = str(var(v[0]))+"*w*("+str(var(v[1]))+"+"+str(var(v[2]))+"+"+str(var(v[3]))+"-2)"

   expn = expn + "+" + e

#### Polynomial Reductions
print("Polynomial Reductions: Reduction by Minimum Selection")
print(expn)
exp = sympy.simplify(expn)
pol = str(exp)

print("\nPOLYNOMIAL:\n")
print(pol)
plos = pol.replace(' + ',' ').replace('*',' ').replace(' - ',' ').split()
l = set() 
plosv = [x for x in plos if x not in l and l.add(x) is None]
var = []
coef = []
for value in plosv:
    try:
        coef.append(float(value))
    except ValueError:
        var.append(str(value))

print("\nCOEFFICIENTS:\n")
print(coef)
print("\nVARIABLES:\n")
print(var)

matrix = np.zeros((len(var),len(var)))
df = pd.DataFrame(matrix, columns=var, index=var)
plos = pol.replace(' + ',' ').replace(' - ',' ').split()

# add coef 1 to non coef variables
lp = []
for c in plos:
    if c[0].isalpha():
       c = '1*' + c
    lp.append(c)

# create Q dictionary
Q = {}
for idx, val in enumerate(lp):
	if val.count('*') == 1:
		sp = val.split("*")
		df[sp[1]][sp[1]] = float(sp[0])
		Q[(sp[1],sp[1])] = float(sp[0])
	if val.count('*') == 2:
		sp = val.split("*")
		df[sp[1]][sp[2]] = float(sp[0])
		Q[(sp[1],sp[2])] = float(sp[0])

print("\nMATRIX:\n")
print(df)
print("\nQUBO DICTIONARY:\n")
print(Q)

print("\nD-WAVE OUTPUT:\n")
start_time = time.time()
response = sampler.sample_qubo(Q, num_reads=5000)
minE = 999999
maxO = 0
# create dataframe if we want to store all values
for datum in response.data(['sample', 'energy', 'num_occurrences']):
    if (datum.energy <= minE and datum.num_occurrences >= maxO):
       minE = datum.energy
       maxO = datum.num_occurrences
       sample = datum.sample
    print(datum.sample, "Energy: ", datum.energy, "Occurrences: ", datum.num_occurrences)

print("\nSAMPLE WITH MINIMUM ENERGY AND MAXIMUM OCCURRENCES:\n")
print(sample, "Energy: ", minE, "Occurrences: ", maxO)

elapsed_time = time.time() - start_time

print("\nTIME (sec):\n")
print(elapsed_time)

sys.stdout = orig_stdout
f.close()
