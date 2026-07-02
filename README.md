# PPC5 – Equação de Blasius



## Descrição

Este trabalho implementa a solução numérica da Equação de Blasius utilizando o Método do Tiro em conjunto com o Método de Runge-Kutta de 4ª Ordem (RK4).

A Equação de Blasius descreve o perfil de velocidade da camada limite laminar sobre uma placa plana.

A equação diferencial é

\[
f''' + \frac{1}{2}ff'' = 0
\]

sujeita às condições de contorno

\[
f(0)=0
\]

\[
f'(0)=0
\]

\[
f'(\infty)=1
\]

O problema é transformado em um sistema de três equações diferenciais ordinárias de primeira ordem e resolvido pelo Método do Tiro.

---

## Métodos Numéricos Utilizados

- Método do Tiro;
- Método da Secante;
- Método de Runge-Kutta de 4ª Ordem (RK4).

---

## Funcionalidades

O programa é capaz de:

- Ler os parâmetros informados pelo usuário;
- Resolver numericamente a Equação de Blasius;
- Determinar o valor convergido de \(f''(0)\);
- Calcular o erro final do Método do Tiro;
- Determinar o número de iterações;
- Calcular o coeficiente local de atrito

\[
C_f=\frac{2f''(0)}{\sqrt{Re_x}}
\]

- Determinar o valor de

\[
\eta_{99}
\]

- Salvar os resultados em arquivo CSV;
- Gerar os gráficos de

  - \(f(\eta)\)
  - \(f'(\eta)\)
  - \(f''(\eta)\)

---

## Arquivos

```
PPC5_Blasius.py
```

Código principal.

```
blasius_resultados.csv
```

Resultados numéricos.

```
f.png
```

Gráfico de \(f(\eta)\).

```
f_linha.png
```

Gráfico de \(f'(\eta)\).

```
f_duas_linhas.png
```

Gráfico de \(f''(\eta)\).

---

## Bibliotecas utilizadas

- NumPy
- Matplotlib

Instalação:

```bash
pip install numpy matplotlib
```

---

## Como executar

Execute o programa:

```bash
python PPC5_Blasius.py
```

Informe os parâmetros solicitados pelo programa.

Exemplo:

```
Primeiro chute: 0.30

Segundo chute: 0.35

Passo Δη: 0.01

η máximo: 10

Tolerância: 1e-8

Máximo de iterações: 50

Rex: 100000
```

---

## Resultados esperados

O programa deve obter aproximadamente

```
f''(0) ≈ 0.332057
```

Além disso,

```
f'(ηmax) ≈ 1
```

e

```
η99 ≈ 4.92


Universidade de Brasília (UnB)

Faculdade de Tecnologia
