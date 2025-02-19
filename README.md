# Flujo de Potencia DC y Análisis de Contingencias

## Marco Teórico

### Flujo de Potencia DC
El flujo de potencia DC es una simplificación del flujo de potencia AC que considera:

1. Magnitudes de voltaje cercanas a 1.0 p.u.
2. Diferencias angulares pequeñas
3. Resistencias de línea despreciables frente a las reactancias
4. Sin pérdidas reactivas

La ecuación fundamental es:
```math
P_{km} = \frac{\theta_k - \theta_m}{x_{km}}
```
donde:
- $P_{km}$: Flujo de potencia activa entre los nodos k y m
- $\theta_k, \theta_m$: Ángulos de voltaje en los nodos k y m
- $x_{km}$: Reactancia de la línea k-m

### Análisis de Contingencias N-1
Evalúa el impacto de la pérdida de una línea en el sistema mediante:
```math
PI = \sqrt{\left(\frac{P_{ij}^{post}}{P_{ij}^{pre}}\right)^2}
```
donde PI es el índice de severidad de la contingencia.

## Funciones Implementadas

### 1. `calcular_ybus(lines, nodes)`
Construye la matriz de admitancia nodal.

**Librerías utilizadas:**
- `DataFrames`: 
  - `nrow()`: Obtiene dimensiones
  - Acceso a columnas: `lines.FROM`, `lines.TO`
- `LinearAlgebra`: Operaciones matriciales

**Entradas (Inputs):**
- `lines`: DataFrame con información de las líneas
  - FROM: Nodo de inicio
  - TO: Nodo final
  - R: Resistencia
  - X: Reactancia
  - B: Susceptancia
- `nodes`: DataFrame con información de los nodos

**Salidas (Outputs):**
- `ybus`: Matriz de admitancia nodal compleja

### 2. `create_Ykm(lines)`
Crea la matriz de admitancias para flujo DC.

**Librerías:**
- `zeros()`: Inicialización de matrices
- `DataFrames`: Acceso a datos de líneas

**Entradas (Inputs):**
- `lines`: DataFrame con información de las líneas
  - FROM: Nodo de inicio
  - TO: Nodo final
  - X: Reactancia

**Salidas (Outputs):**
- `Ykm`: Matriz de admitancias DC


### 3. `create_P_vector(nodes)`
Genera el vector de potencias netas.

**Librerías:**
- `DataFrames`: Acceso a `PLOAD` y `PGEN`
- `zeros()`: Inicialización de vectores

**Entradas (Inputs):**
- `nodes`: DataFrame con información de los nodos
  - PLOAD: Potencia de carga
  - PGEN: Potencia de generación

**Salidas (Outputs):**
- `P`: Vector de potencias netas en p.u.

### 4. `DC_power_flow(Ykm, P, nlines, nnodes, lines)`
Resuelve el flujo de potencia DC.

**Librerías:**
- `LinearAlgebra`: 
  - Operador `\`: Solución de sistemas lineales
- `Base`:
  - `setdiff()`: Eliminación del nodo slack
  - `pushfirst!()`: Manipulación de vectores

**Entradas (Inputs):**
- `Ykm`: Matriz de admitancias DC
- `P`: Vector de potencias netas
- `nlines`: Número de líneas
- `nnodes`: Número de nodos
- `lines`: DataFrame con datos de líneas

**Salidas (Outputs):**
- `d`: Vector de ángulos de voltaje
- `pf`: Vector de flujos de potencia en las líneas

### 5. `Contingency(nlines, Ykm, P, nnodes, lines)`
Realiza análisis de contingencias N-1.

**Librerías:**
- `Base`:
  - `push!()`: Agregar elementos a vectores
  - `sqrt()`: Cálculo de índices
- `DataFrames`: Acceso a datos de líneas

**Entradas (Inputs):**
- `nlines`: Número de líneas
- `Ykm`: Matriz de admitancias original
- `P`: Vector de potencias netas
- `nnodes`: Número de nodos
- `lines`: DataFrame con datos de líneas

**Salidas (Outputs):**
- `almacenamiento`: Matriz de flujos post-contingencia
- `almrank`: Matriz de índices de severidad

### 6. `graficar_contingencias(rank, num_conting)`
Visualiza resultados del análisis de contingencias.

**Librerías:**
- `Plots`:
  - `heatmap()`: Mapa de calor
  - Configuración de gráficos
- `StatsPlots`: Extensiones estadísticas

**Entradas (Inputs):**
- `rank`: Matriz de índices de contingencia
- `num_conting`: Número de contingencias

**Salidas (Outputs):**
- `p`: Objeto plot con el mapa de calor

## Resultados
![\[Aquí puedes incluir gráficos o resultados relevantes\]](plot_1.svg)

## Licencia
Hecho por: Jhon Edward Bedoya Olarte
j.bedoya3@utp.edu.co

<p xmlns:cc="http://creativecommons.org/ns#">Esta obra está licenciada bajo <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1">CC BY 4.0</a></p>
