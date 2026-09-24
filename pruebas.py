from simulator import *

# PROVA RÀPIDA A CLASSE

rho_10k = atmosfera_isa(10000) # Densitat a 10.000 ft
print(f"Densitat a 10.000 ft: {rho_10k:.4f} kg/m^3")

v_prova = calcular_v_minrod(aviones['A320-212'], 64500, rho_10k, 5000)
print(f"Velocitat V_minROD per l'A320: {v_prova:.2f} m/s ({v_prova * 1.94384:.2f} nuls)")

