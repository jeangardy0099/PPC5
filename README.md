
import numpy as np
import matplotlib.pyplot as plt

# Função que define o sistema de EDOs da equação de Blasius
def deriv(y):
    y1, y2, y3 = y
    return np.array([y2, y3, -0.5*y1*y3], dtype=float)

# Realiza um passo do método de Runge-Kutta de 4ª ordem
def rk4_step(y, h):
    k1 = deriv(y)
    k2 = deriv(y + 0.5*h*k1)
    k3 = deriv(y + 0.5*h*k2)
    k4 = deriv(y + h*k3)
    return y + h*(k1 + 2*k2 + 2*k3 + k4)/6

# Integra o sistema de EDOs para um valor de s = f''(0)
def integrar(s, h, etamax):
    n = int(etamax/h)+1
    eta = np.linspace(0, etamax, n)
    y = np.zeros((n,3))
    
    # Condições iniciais
    y[0] = [0.0,0.0,s]

    # Integração usando RK4
    for i in range(n-1):
        y[i+1] = rk4_step(y[i], h)
    return eta, y

# Calcula o erro do método do tiro
def erro(s,h,etamax):
    _, y = integrar(s,h,etamax)
    return y[-1,1]-1.0

# Método da Secante para determinar o valor correto de f''(0)
def secante(s0,s1,h,etamax,tol,maxit):
    e0=erro(s0,h,etamax)
    e1=erro(s1,h,etamax)
    for k in range(maxit):
        s2 = s1 - e1*(s1-s0)/(e1-e0)
        e2 = erro(s2,h,etamax)

        # Verifica convergência
        if abs(e2)<tol:
            return s2,k+1,e2

        # Atualiza os chutes    
        s0,e0=s1,e1
        s1,e1=s2,e2
    return s2,maxit,e2

# Programa principal
def main():
    print("===  Equação de Blasius ===")

    # Entrada de dados
    s0=float(input("Primeiro chute f''(0): "))
    s1=float(input("Segundo chute f''(0): "))
    h=float(input("Passo Δη: "))
    etamax=float(input("η máximo: "))
    tol=float(input("Tolerância: "))
    maxit=int(input("Máximo de iterações: "))

    # Método do tiro
    s,it,err=secante(s0,s1,h,etamax,tol,maxit)

    # Integração final
    eta,y=integrar(s,h,etamax)


     # Salva os resultados em arquivo CSV
    np.savetxt("blasius_resultados.csv",
               np.column_stack((eta,y)),
               delimiter=",",
               header="eta,f,f_linha,f_duas_linhas",
               comments="")

    # Determina η99
    eta99 = np.interp(0.99,y[:,1],eta)

    # Exibe resultados
    print("\n===== RESULTADOS =====")
    print(f"f''(0) = {s:.6f}")
    print(f"Valor clássico = 0.332057")
    print(f"Erro final = {err:.3e}")
    print(f"Iterações = {it}")
    print(f"f'(etamax) = {y[-1,1]:.6f}")
    print(f"eta99 = {eta99:.6f}")
    print(f"Cδ ≈ {eta99:.6f}")

    # Cálculo do coeficiente de atrito
    Rex=float(input("\nInforme Rex para calcular Cf: "))
    Cf=2*s/np.sqrt(Rex)
    print(f"Cf = {Cf:.8f}")

    # Gráfico de f(η)
    plt.figure()
    plt.plot(eta,y[:,0])
    plt.grid()
    plt.xlabel("η")
    plt.ylabel("f(η)")
    plt.savefig("f.png")

    # Gráfico de f'(η)
    plt.figure()
    plt.plot(eta,y[:,1])
    plt.grid()
    plt.xlabel("η")
    plt.ylabel("f'(η)")
    plt.savefig("f_linha.png")

    # Gráfico de f''(η)
    plt.figure()
    plt.plot(eta,y[:,2])
    plt.grid()
    plt.xlabel("η")
    plt.ylabel("f''(η)")
    plt.savefig("f_duas_linhas.png")

    # Exibe os gráficos na tela
    plt.show()
# Executa o programa
if __name__=="__main__":
    main()
