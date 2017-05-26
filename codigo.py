import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc("text", usetex = True)

E_minima = -1.0

E_maxima = -0.01

deltaE = 0.01

x_minimo = 0.001

x_maximo = 20

deltax = 0.001

x = np.linspace(x_minimo, x_maximo, (x_maximo - x_minimo)/deltax + 1)

energias_posibles = np.linspace(E_minima, E_maxima, (E_maxima - E_minima)/deltaE + 1)

error_minimo = 0.000000000000001

class atomo:
    
    def __init__(self, Z, n_array, l_array):
        self.nucleo = nucleo(Z)
        self.electrones = [electron(n, l) for (n, l) in zip(n_array, l_array)]
        self.electrones[0].potencial = self.nucleo.potencial(x)
        self.electrones[0].actualizar()
        for i in range(1,len(self.electrones)):
            self.electrones[i].potencial = self.nucleo.potencial(x)
            for j in range(0,i):
                self.electrones[i].potencial += self.electrones[j].potencial
            self.electrones[i].actualizar()
    
    def actualizar(self):
        for e in self.electrones:
            e.potencial = self.nucleo.potencial(x)
            for eprima in self.electrones:
                if eprima != e:
                    potencialprima = np.zeros(len(x))
                    potencialprima[len(potencialprima)-1] = 0
                    for i in range(1, len(potencialprima)):
                        potencialprima[len(potencialprima)-1-i] = potencialprima[len(potencialprima)-1-(i-1)] + 2*eprima.q[len(potencialprima)-1-i]/(x[len(potencialprima)-1-i]*x[len(potencialprima)-1-i]*self.nucleo.Z)
                    e.potencial += potencialprima
            e.actualizar()
    
class nucleo:
    
    def __init__(self, Z):
        nucleo.Z = Z
        
    def potencial(self, x):
        return -2.0/x
    
class electron:
    
    def __init__(self, n, l):
        self.n = n
        self.l = l
        self.potencial = np.zeros(len(x))
        self.q = np.zeros(len(x))
        self.psi = np.zeros(len(x))
        self.Dpsi = np.zeros(len(x))
        if n == 1:
            self.psi[0] = 1
            self.Dpsi[0] = -1
        else:
            self.psi[0] = 0
            self.Dpsi[0] = 1
        self.Psi = self.psi[0]
        self.DPsi = self.Dpsi[0]
        self.E = 0
    
    def actualizar_paso(self, E, posicion):
        self.DPsi += -deltax*(2.0*self.DPsi/x[posicion] + (E - self.potencial[posicion] - self.l*(self.l + 1)/(x[posicion]*x[posicion]))*self.Psi) 
        self.Psi += deltax*self.DPsi
        
    def integrador_sin_memoria(self, E):
        self.Psi = self.psi[0]
        self.DPsi = self.Dpsi[0]
        for j in range(0, len(x) - 1):
            self.actualizar_paso(E, j)
                                                                                                   
    def actualizar_energia(self):
        energias_intervalos = []
        self.integrador_sin_memoria(energias_posibles[0])
        anterior = np.sign(self.DPsi)
        for i in range(1, len(energias_posibles)):
            self.integrador_sin_memoria(energias_posibles[i])
            presente = np.sign(self.DPsi)
            if presente != anterior:
                energias_intervalos.append([energias_posibles[i-1], energias_posibles[i]])
            anterior = presente
        energia_intervalo = energias_intervalos[self.n-1]
        error = energia_intervalo[1] - energia_intervalo[0]
        while error > error_minimo:
            self.integrador_sin_memoria(energia_intervalo[0])
            izquierda = np.sign(self.DPsi)
            self.integrador_sin_memoria((energia_intervalo[1]+energia_intervalo[0])/2)
            central = np.sign(self.DPsi)
            if izquierda != central:
                energia_intervalo[1] = (energia_intervalo[1]+energia_intervalo[0])/2
            else:
                energia_intervalo[0] = (energia_intervalo[1]+energia_intervalo[0])/2
            error = energia_intervalo[1] - energia_intervalo[0]
        self.E = (energia_intervalo[1]+energia_intervalo[0])/2

    def integrador(self):
        self.Psi = self.psi[0]
        self.DPsi = self.Dpsi[0]
        for j in range(0, len(x) - 1):
            self.actualizar_paso(self.E, j)
            self.psi[j+1] = self.Psi
            self.Dpsi[j+1] = self.DPsi
            
    def actualizar_carga(self):
        P = 4.0*np.pi*x*x*self.psi*self.psi
        P = P/(np.sum(P)*deltax)
        for i in range(1, len(self.q)):
            self.q[i] = self.q[i-1] - P[i]*deltax

    def actualizar(self):
        self.actualizar_energia()
        self.integrador()
        self.actualizar_carga()
            
n_iteraciones = 30
n_electrones = 3

litio = atomo(n_electrones, [1, 1, 2], [0, 0, 1])
energias = np.zeros((n_iteraciones, n_electrones))
for i in range(0, n_iteraciones):
    for j in range(0, n_electrones):
        energias[i, j] = litio.electrones[j].E
    litio.actualizar()
color = ["k", "b", "g"]
fig, ax = plt.subplots(1, 2)
for i in range(0, n_electrones):
    ax[0].scatter(range(0, n_iteraciones), energias[:, i], c = color[i], alpha = 0.6, label = r"Electr\'on %d: $%fE_0$" % (i + 1, energias[n_iteraciones-1,i]))
ax[0].set_title(r"Energ\'ia ($%fE_0$)" % np.sum(energias[n_iteraciones - 1, :]))  
ax[0].set_xlabel(r"Iteraci\'on")
ax[0].set_ylabel(r"Energ\'ia ($E_0$)")
for i in range(0,n_electrones):
	ax[1].plot(x, litio.electrones[i].psi, c = color[i], label = r"Electr\'on %d" % (i+1))
ax[1].set_title(r"Funciones de onda radiales")
ax[1].set_xlabel(r"$r(a_0/3)$")
ax[1].set_ylabel(r"$\R(r)$") 
ax[0].legend()
ax[1].legend()
plt.tight_layout()
fig.savefig("1s_1s_2p.pdf")
plt.close(fig)
