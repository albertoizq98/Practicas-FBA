#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import numpy as np
import scipy as sc
import cobra.test
model=cobra.test.create_test_model("textbook")


# In[2]:


import tools
model=tools.desdoblar(model,"desdoblado")


# In[3]:


tipos=tools.tiposReaccionesD(model)
reversibles=tipos[0]
irreversibles=tipos[1]


# In[4]:


S=cobra.util.create_stoichiometric_matrix(model,array_type="dense")
irreversibles_i=[]
model=tools.removeDeadEnds(model,S,reversibles,irreversibles,irreversibles_i)


# In[5]:


#Restricciones que le pongo

coeficientes=dict()
for i in irreversibles:
    r=model.reactions[i]
    coeficientes[r.forward_variable]=1
constraint=model.problem.Constraint(0,lb=1,ub=1)
model.add_cons_vars(constraint)
model.solver.update()
constraint.set_linear_coefficients(coefficients=coeficientes)


# In[6]:


#Coger la función

lista= np.ones(len(model.reactions))
valores=dict()
for i in range(len(lista)):
    r=model.reactions[i]
    valores[r]=lista[i]
model.objective=valores


# In[7]:


#Calcular la solucion del EFM
solucion=model.optimize(objective_sense="minimize")
valores=solucion.fluxes.values


# In[8]:


#Chequear si la EFM es una EFM
tools.check(valores,S)


# In[9]:


EFMS_utiles=tools.calculoEFM_test(model, irreversibles,S,1000)


# In[10]:


lista_estadistica=[]
for i in EFMS_utiles:
    lista_estadistica.append(len(i))


# In[11]:


cantidad=len(lista_estadistica)
print("La cantidad de EFMs es:", cantidad)


# In[12]:


#Media
media=np.mean(lista_estadistica)
print("La media es:", media)





# In[13]:


#Varianza
varianza=np.var(lista_estadistica)
print("La varianza es:", varianza)


# In[14]:


#Desviación típica
desviacion=np.std(lista_estadistica)
print("La desviacion tipica es:", desviacion)


# In[15]:


#Mediana
mediana=np.median(lista_estadistica)
print("La mediana es:", mediana)


# In[16]:


import matplotlib.pyplot as plt

#Histograma
plt.hist(lista_estadistica)

plt.savefig('Histograma_basico.png')
plt.close()


# In[17]:


#Boxplot
plt.boxplot(lista_estadistica)

plt.savefig('Boxplot_basico.png')
plt.close()

# In[18]:


cuenta=[]
lista_incrementado=[]
for i in lista_estadistica:
    cuenta.append(i)
    lista_incrementado.append(np.mean(cuenta))
plt.plot(lista_incrementado)

plt.savefig('Plotmedias_basico.png')
plt.close()
