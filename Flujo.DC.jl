using LinearAlgebra
using DataFrames
using CSV
using StatsPlots
using SparseArrays
using Plots
# Calcular la matriz de admitancia nodal

function calcular_ybus(lines, nodes)
    """
    Entradas: lines: DataFrame
              nodes: DataFrame
    Salidas:  Ybus: Matriz
    """
    num_nodes = nrow(nodes) # Número de nodos  
    num_lines = nrow(lines) # Número de líneas
    ybus = zeros(num_nodes, num_nodes) * 1im # Matriz de admitancias nodales
    for k in 1:num_lines # Se recorre cada línea
        n1 = lines.FROM[k] # Nodo de inicio
        n2 = lines.TO[k] # Nodo final
        yL = 1 / (lines.R[k] + lines.X[k] * 1im) # Admitancia de la línea
        Bs = lines.B[k] * 1im / 2 # Admitancia shunt

        ybus[n1, n1] += yL + Bs # Diagonal
        ybus[n1, n2] -= yL # Fuera
        ybus[n2, n1] -= yL # Fuera
        ybus[n2, n2] += yL + Bs # Diagonal
    end
    return ybus
end

## Funcion principal que
lines = DataFrame(CSV.File("lines.csv")) # Datos de las líneas
nodes = DataFrame(CSV.File("nodes.csv")) # Datos de los nodos
ybus = calcular_ybus(lines, nodes) # Matriz de admitancias nodales
##ybus = sparse(ybus)
Z_bus = inv(ybus) ## Inversa de la matriz Ybus

#Se crean las funciones  flujo de carga DC""" 

function carga_datos() # Función para cargar los datos
    # Importando los datos y creando los DataFrames
    lines = DataFrame(CSV.File("lines.csv")) # Datos de las líneas
    nodes = DataFrame(CSV.File("nodes.csv")) # Datos de los nodos
    num_nodes = nrow(nodes) # Número de nodos
    num_lines = nrow(lines) # Número de líneas
    
    return lines, nodes, num_lines, num_nodes # Se retornan los DataFrames y el número de nodos y líneas
end

# Función para crear la matriz de admitancias
function create_Ykm(lines) 
    num_nodes = nrow(nodes) # Número de nodos 
    num_lines = nrow(lines) # Número de líneas
    
    Ykm = zeros(num_nodes, num_nodes) # Matriz de admitancias
    
    for i in 1:num_lines # Se recorre cada línea
        k = lines.FROM[i] # Nodo de inicio
        m = lines.TO[i] # Nodo final
        Y_km =  1 / (lines.X[i]) # Admitancia de la línea
        Ykm[k, m] = Ykm[k, m] - Y_km # Fuera
        Ykm[m, k] = Ykm[m, k] - Y_km # Fuera
        Ykm[k, k] = Ykm[k, k] + Y_km # Diagonal
        Ykm[m, m] = Ykm[m, m] + Y_km # Diagonal
    end
    return Ykm # Se retorna la matriz de admitancias
end
# Función para crear vector de potencias
function create_P_vector( nodes) 
    num_nodes = nrow(nodes) # Número de nodos
    P = zeros(num_nodes) # Vector de potencias
    
    for i in 1:num_nodes # Se recorre cada nodo
        Pd = nodes.PLOAD[i] # Carga
        Pg = nodes.PGEN[i] # Generación
        P[i] = Pg - Pd # Potencia 
    end
    return P
end

# Función flujo de potencia DC
function DC_power_flow(Ykm, P, nlines, nnodes, lines) 
    slack = 1 # Nodo slack
    dslack = 0 # Voltaje en el nodo slack
    num_lines = nlines # Número de líneas
    num_nodes = nnodes # Número de nodos
    nodos = setdiff(1:num_nodes, slack) # Nodos sin slack
    Ykm1 = Ykm[nodos, nodos] # Matriz de admitancias reducida
    P = P[nodos] # Vector de potencias
    d = zeros(num_nodes) # Vector de voltajes
    d = Ykm1 \ P # Cálculo de los voltajes
    pf = zeros(num_lines) # Vector de potencia de flujo en las líneas
    d = pushfirst!(d, dslack) # Se añade el slack al vector de voltajes
    for i in 1:num_lines # Cálculo de la potencia de flujo en las líneas
        k = lines.FROM[i] # Nodo de inicio
        m = lines.TO[i] # Nodo final
        pf[i] = (d[k] - d[m]) / lines.X[i] # Potencia de flujo en la línea
    end
    return d, pf # Se retornan los voltajes y las potencias de flujo
end


# Función para análisis de contingencias
function Contingency(nlines, Ykm, P, nnodes, lines)
    num_conting = nlines # Número de contingencias
    Ykm1 = create_Ykm(lines) # Matriz de admitancias
    almacenamiento = zeros(num_conting, nlines) # Almacenamiento de las potencias de flujo
    
    dref, pfref = DC_power_flow(Ykm, P, nlines, nnodes, lines) # Flujo de potencia en operación normal
    almacenamiento = [] # Almacenamiento de las potencias de flujo

    for j in 1:num_conting # Se recorre cada contingencia
        k = lines.FROM[j]  # Nodo de inicio
        m = lines.TO[j] # Nodo final
        Ykm1[k, m] = 0 # Se elimina la línea
        Ykm1[m, k] = 0 # Se elimina la línea
        df, pf = DC_power_flow(Ykm1, P, nlines, nnodes, lines) # Flujo de potencia en contingencia
        push!(almacenamiento, pf) # Se almacena la potencia de flujo
    end
    
    almrank = []
    for k in 1:num_conting # Se recorre cada contingencia
        for i in 1:num_conting # Se recorre cada contingencia
            rank = sqrt((almacenamiento[k][i] / pfref[k])^2) # Índice de contingencia
            push!(almrank, rank) # Se almacena el índice de contingencia
        end
    end

    #Se crea el ciclo for para descomponer el vector de almacenamiento
    
    for i in 1:num_conting # Se recorre cada contingencia
        almrank[i] = almrank[(i - 1) * num_conting + 1:i * num_conting] # Se descompone el vector
    end

    return almacenamiento, almrank 
end



# Función para graficar los análisis de contingencia

function graficar_contingencias(rank, num_conting)
    # Crear y llenar la matriz para el mapa de calor
    matriz_calor = zeros(num_conting, num_conting)
    for i in 1:num_conting
        matriz_calor[i, :] = rank[i]
    end
    
    # Calcular valores para la escala de color
    val_min = minimum(filter(x -> x > 0, matriz_calor))
    val_max = maximum(matriz_calor)
    levels = range(val_min, val_max/2, length=5)
    
    # Configurar etiquetas para los ejes
    xticks_pos = 1:1:num_conting
    xticks_labels = string.(xticks_pos)
    yticks_pos = 1:1:num_conting
    yticks_labels = string.(yticks_pos)
    
    # Crear el mapa de calor con mejor visualización
    p = heatmap(
        matriz_calor,
        xlabel="Línea afectada",
        ylabel="Línea en falla",
        title="Análisis de Contingencias N-1",
        size=(1400, 1200),  # Aumentar tamaño para mejor visualización
        color=:magma,
        colorbar_title="Índice de Contingencia",
        aspect_ratio=:equal,
        xticks=(xticks_pos, xticks_labels),
        yticks=(yticks_pos, yticks_labels),
        xrotation=45,  # Rotación diagonal para mejor lectura
        xtickfont=font(7),  # Ajustar tamaño de fuente
        ytickfont=font(7),  # Ajustar tamaño de fuente
        clims=(val_min, val_max/2),
        colorbar_ticks=levels,
        colorbar_tickfontsize=10,
        framestyle=:box,
        margins=20Plots.mm,  # Aumentar márgenes
        bottom_margin=50Plots.px,  # Margen inferior adicional para las etiquetas
        right_margin=20Plots.px   # Margen derecho para la barra de color
    )
    
    # Mostrar el gráfico
    display(p)
    
    return p
end

function main()
    lines, nodes, nlines, nnodes = carga_datos() # Cargar los datos
    Ykm = create_Ykm(lines) # Crear la matriz de admitancias
    P = create_P_vector(lines, nodes) # Crear el vector de potencias
    dref, pfref = DC_power_flow(Ykm, P, nlines, nnodes, lines) # Flujo de potencia en operación normal
    println("  ") # Se imprime un espacio en blanco
    print("El flujo de potencia en operación normal en las líneas es: ", pfref) # Se imprime el flujo de potencia en operación normal
    println("  ") # Se imprime un espacio en blanco
    println("Los ángulos de voltaje en los nodos son: ", dref)
    println("  ")

"""Para ver el resultado del flujo de potencia (angulo de tension y flujo) se debe comentar desde la linea 162 hasta la linea 194"""


    pfconting, rank = Contingency(nlines, Ykm, P, nnodes, lines) # Flujo de potencia en contingencia
    println("  ") # Se imprime un espacio en blanco
    for i in 1:nlines # Se recorre cada línea
        k = lines.FROM[i] # Nodo de inicio
        m = lines.TO[i] # Nodo final
        println("  ") # Se imprime un espacio en blanco 
        println("El flujo de potencia ante contingencia en la línea $i del nodo $k al $m es: ", pfconting[i]) # Se imprime el flujo de potencia ante contingencia
    end
    println("  ") # Se imprime un espacio en blanco
    
    
    
    x = 0 # Variable para almacenar el índice de contingencia
    # Se crea el ciclo para clasificar las líneas más críticas según el índice correspondiente a cada contingencia
    for j in 1:nlines # Se recorre cada línea
        k = lines.FROM[j] # Nodo de inicio
        m = lines.TO[j] # Nodo final
        indice_ordenado = sortperm(rank[j]) # Índice de contingencia ordenado
        num_lineas_criticas = 5 # Número de líneas más críticas
        lineas_criticas = indice_ordenado[end - num_lineas_criticas + 1:end] # Líneas más críticas
        println("  ") # Se imprime un espacio en blanco
        println("Las $num_lineas_criticas líneas más críticas ante contingencia en la línea $j del nodo $k al $m son:") # Se imprime las líneas más críticas
        println("  ") # Se imprime un espacio en blanco
       
        for i in lineas_criticas # Se recorre cada línea
            k = lines.FROM[i] # Nodo de inicio
            m = lines.TO[i] # Nodo final
            println("Línea $i del nodo $k al $m - Índice de Contingencia: ", rank[j][i]) # Se muestran los datos
        end
    end
    
    # Llamar a la función para graficar los análisis de contingencia
    graficar_contingencias(rank, nlines) 
    
    return nothing
end

# Llamada a la función principal
main()