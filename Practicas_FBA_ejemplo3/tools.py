import cobra
from cobra import Model, Reaction, Metabolite
import numpy as np

def desdeMatriz(nombre,nombreModelo):
  S=np.loadtxt(nombre)
  m=len(S)
  n=len(S[0])
  model=Model(nombreModelo)
  metabolitos=[]
  reacciones=[]
  for i in range(m):
    metabolito=Metabolite("M"+str(i),"Metabolito "+str(i))
    metabolito.compartment="c"
    metabolitos.append(metabolito)
  model.add_metabolites(metabolitos)
  for j in range(n):
     reaccion=Reaction("R"+str(j),"Reaccion "+str(j))
     reaccion.lower_bound=0
     reaccion.upper_bound=1000
     misMetabolitos=dict()
     for i in range(m):
       if not S[i][j]==0:
         misMetabolitos[model.metabolites[i]]=S[i][j]
     reaccion.add_metabolites(misMetabolitos)
     reacciones.append(reaccion)
  model.add_reactions(reacciones)
  return(model)
  
  
def tiposReacciones(modelo):
  bloqueadas=[]
  irreversibles=[]
  irreversibles_i=[]
  reversibles=[]
  for r in modelo.reactions:
      modelo.objective=r
      solucion1=modelo.optimize(objective_sense="minimize")
      solucion2=modelo.optimize(objective_sense="maximize")
      minimo=solucion1.objective_value
      maximo=solucion2.objective_value
      if minimo==0 and maximo==0:
          bloqueadas.append(r)
      if minimo<0 and maximo>0:
          reversibles.append(r)
      if minimo>=0 and maximo>0:
          irreversibles.append(r)
      if maximo<=0 and minimo<0:
          irreversibles_i.append(r)
  print("Irrev",len(irreversibles), "Irrev inversas",len(irreversibles_i),"Reversibles",len(reversibles),"Bloqueadas",len(bloqueadas))
  
  modelo.remove_reactions(bloqueadas)
  
  irreversibles=[]
  irreversibles_i=[]
  reversibles=[]
  for i in range(len(modelo.reactions)):
      r=modelo.reactions[i]
      modelo.objective=r
      solucion1=modelo.optimize(objective_sense="minimize")
      solucion2=modelo.optimize(objective_sense="maximize")
      minimo=solucion1.objective_value
      maximo=solucion2.objective_value
      if minimo==0 and maximo==0:
          bloqueadas.append(i)
      if minimo<0 and maximo>0:
          reversibles.append(i)
      if minimo>=0 and maximo>0:
          irreversibles.append(i)
      if maximo<=0 and minimo<0:
          irreversibles_i.append(i) 
  resultado=[]
  resultado.append(modelo)
  resultado.append(irreversibles)
  resultado.append(irreversibles_i)
  resultado.append(reversibles)
  return(resultado) 

def tiposReaccionesD(modelo):
  reversibles2=[]
  irreversibles2=[]
  for i in range(len(modelo.reactions)):
      r=modelo.reactions[i]
      
      texto=r.id
      
      if texto[len(texto)-3:]=="rev":
        r1=modelo.reactions[i-1]
        texto1=r1.id
        if texto==texto1+"_rev":
          reversibles2.append(i-1)
          reversibles2.append(i)
  for i in range(len(modelo.reactions)):
      if not i in reversibles2:
          irreversibles2.append(i)
  print("Irrev",len(irreversibles2),"Reversibles",len(reversibles2))
  return([reversibles2,irreversibles2])
  
  
def removeDeadEnds(modelo,S,reversibles,irreversibles,irreversibles_i):
  metabolitos_deadEnd=[]
  for i in range(len(S)):
    pos=False
    neg=False
    fila=S[i]
    for j in range(len(fila)):
      if not fila[j]==0:
        if j in reversibles:
          pos=True
          neg=True
        elif j in irreversibles:
          if fila[j]>0:
            pos=True
          else:
            neg=True
        else:
          if fila[j]>0:
            neg=True
          else:
            pos=True
    if (not pos) or (not neg):
        metabolitos_deadEnd.append(modelo.metabolites[i])
  modelo.remove_metabolites(metabolitos_deadEnd)
  return(modelo)
  
def desdoblar(model,nombre):
    error = 10**-12
    import cobra
    from cobra import Model, Reaction, Metabolite 
    desdoble=Model("Desdoblado")
    metabolitos=[]
    for m in model.metabolites:
        mm = Metabolite(m.id,m.name)
        mm.compartment='c'
        metabolitos.append(mm)
    desdoble.add_metabolites(metabolitos)
    reacciones=[]
    for r in model.reactions:
        model.objective=r
        solucion1=model.optimize(objective_sense="minimize")
        solucion2=model.optimize(objective_sense="maximize")
        minimo=solucion1.objective_value
        maximo=solucion2.objective_value
        
        if abs(minimo)>error or abs(maximo)>error:
          if minimo<0 and maximo>0:
            reaccion1=Reaction(r.id,r.name)
            reaccion1.lower_bound=0
            reaccion1.upper_bound=r.upper_bound
            for key in r.metabolites.keys():
                reaccion1.add_metabolites({
                    key: r.metabolites[key],
                })
            
            reaccion2=Reaction(r.id+"_rev",r.name+"_rev")
            reaccion2.lower_bound=0
            reaccion2.upper_bound=-r.lower_bound
            for key in r.metabolites.keys():
                reaccion2.add_metabolites({
                    key: -r.metabolites[key],
                })
            reacciones.append(reaccion1)
            reacciones.append(reaccion2)
        
          if minimo>=0 and maximo>0:
            r.lower_bound=0
            reacciones.append(r)
        
        
          if minimo<0 and maximo<=0:
            reaccion=Reaction(r.id+"_i",r.name+"_i")
            reaccion.upper_bound=-r.lower_bound
            reaccion.lower_bound=0
            for key in r.metabolites.keys():
                reaccion.add_metabolites({
                    key: -r.metabolites[key],
                })
            reacciones.append(reaccion)
            
    desdoble.add_reactions(reacciones)
    print("OK")
    return(desdoble)

def listaEnteros(lista):
  lista2=[]
  for i in lista:
    lista2.append(int(i))
  return lista2


def check(valores,S):
  soporte=[]
  noSoporte=[]
  for i in range(len(valores)):
    if valores[i]==0:
      noSoporte.append(i)
    else:
      soporte.append(i)
  reaccionesActivas=len(soporte)
  S1=np.delete(S,noSoporte,1)
  numeroEcuaciones=np.linalg.matrix_rank(S1)
  return(reaccionesActivas == numeroEcuaciones+1)
  
def calculoEFM_test(model, irreversibles,S,n):
   EFMS=[]
   for i in range(n):
       coeficientes=dict()
       irreversibles_select= np.random.choice(irreversibles,10,replace=False)
       for i in irreversibles_select:
           r=model.reactions[i]
           coeficientes[r.forward_variable]=np.random.rand()
       constraint=model.problem.Constraint(0,lb=1,ub=1)
       model.add_cons_vars(constraint)
       model.solver.update()
       constraint.set_linear_coefficients(coefficients=coeficientes)
       lista= np.random.rand(len(model.reactions))
       valores=dict()
       for i in range(len(lista)):
           r=model.reactions[i]
           valores[r]=lista[i]
       model.objective=valores
       solucion=model.optimize(objective_sense="minimize")
       valores=solucion.fluxes.values
       soporte=[]
       noSoporte=[]
       for i in range(len(valores)):
           if valores[i]==0:
               noSoporte.append(i)
           else:
               soporte.append(i)
       reaccionesActivas=len(soporte)
       S1=np.delete(S,noSoporte,1)
       numeroEcuaciones=np.linalg.matrix_rank(S1)
       if (reaccionesActivas == numeroEcuaciones+1):
         if soporte not in EFMS:
           EFMS.append(soporte)
       model.remove_cons_vars(constraint)
   return EFMS