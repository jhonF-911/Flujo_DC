using LinearAlgebra
using DataFrames
using CSV
using Plots
using SparseArrays

# Calcular la matriz de admitancia nodal
function calcular_ybus(lines, nodes)
num_nodes = nrow(nodes)
num_lines= nrow(lines)
Ybus = zeros(num_nodes, num_nodes)*1im
for k =1:num_lines
    n1 = lines.FROM[k]
    n2 = lines.TO[k]
    yL = 1/(lines.R[k] + lines.X[k]*1im)
    Bs = lines.B[k]*1im/2
    Ybus[n1,n1] += yL + Bs # diagonal
    Ybus[n1,n2] -= yL # fuera
    Ybus[n2,n1] -= yL # fuera
    Ybus[n2,n2] += yL + Bs # diagonal

end  
return Ybus
end

"""
    Entradas : lineas: DataFrames con la información de las líneas
               nodos: DataFrames con la información de los nodos
    Salidas : Ybus: Matriz de admitancia nodal
    -
"""

## Funcion principla
lines = DataFrame(CSV.File("lines.csv"))
nodes = DataFrame(CSV.File("nodes.csv"))
Ybus = calcular_ybus(lines, nodes)
#Ybus  = sparse(Ybus)

# inv(Ybus)= Zbus

