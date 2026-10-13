from simulator import lista_aviones, simular_cdo

mapeo_modelos = {
    "B767": "B767-300ER",
    "B737": "B737",
    "B777": "B777-300",
    "A319": "A319-131",
    "A320": "A320-212",
}

# ==============================================================================
# DADES D'ENTRADA (Scenario 1)
# ==============================================================================
tabla_inicial = [
    {"Aircraft model": "B767", "STAR": "ALBER1Z", "%MLW IAF": 80,  "Dist. STAR entry WP to IAF [km]": 129.0844},
    {"Aircraft model": "B737", "STAR": "PUMAL1Z", "%MLW IAF": 100, "Dist. STAR entry WP to IAF [km]": 94.4520},
    {"Aircraft model": "B777", "STAR": "MARTA3Z", "%MLW IAF": 100, "Dist. STAR entry WP to IAF [km]": 178.3476},
    {"Aircraft model": "B767", "STAR": "MATEX3Z", "%MLW IAF": 80,  "Dist. STAR entry WP to IAF [km]": 190.0152},
    {"Aircraft model": "A319", "STAR": "LOBAR2W", "%MLW IAF": 80,  "Dist. STAR entry WP to IAF [km]": 151.3084},
    {"Aircraft model": "A320", "STAR": "CASPE2W", "%MLW IAF": 100, "Dist. STAR entry WP to IAF [km]": 164.0872},
]


# ==============================================================================
# FUNCIONS: UNA PER A CADA COLUMNA
# ==============================================================================

def obtener_valor_en_distancia(x_sim, y_sim, dist_km):
    x_obj = -dist_km * 1000.0

    mejor_idx = 0
    menor_error = float("inf")

    for i in range(len(x_sim)):
        error = abs(x_sim[i] - x_obj)
        if error < menor_error:
            menor_error = error
            mejor_idx = i

    return y_sim[mejor_idx]


def calcular_alt_at_star_entry_wp_m(x_sim, h_sim, dist_km):
    """Columna E: Altitud exacta en el punt d'entrada de la STAR."""
    return round(obtener_valor_en_distancia(x_sim, h_sim, dist_km), 1)


def calcular_time_from_star_entry_wp_to_iaf_s(x_sim, t_sim, dist_km):
    """Columna F: Temps de vol fins a l'IAF."""
    return round(obtener_valor_en_distancia(x_sim, t_sim, dist_km), 1)


def calcular_toa_iaf(tiempo_vuelo_s):
    """Columna G: Sumar els segons de vol a les 11:45:00."""
    segundos_base = 11 * 3600 + 45 * 60
    segundos_totales = int(segundos_base + tiempo_vuelo_s)

    horas = (segundos_totales // 3600) % 24
    minutos = (segundos_totales % 3600) // 60
    segundos = segundos_totales % 60

    return f"{horas:02d}:{minutos:02d}:{segundos:02d}"


def ordenar_vuelos(lista):
    """Columna H: Assigna l'ordre de seqüència segons el temps d'arribada a l'IAF."""
    n = len(lista)

    # Bucle tradicional para ordenar de menor a mayor según '_tiempo_s'
    for i in range(n):
        for j in range(0, n - i - 1):
            if lista[j]["_tiempo_s"] > lista[j + 1]["_tiempo_s"]:
                temp = lista[j]
                lista[j] = lista[j + 1]
                lista[j + 1] = temp

    # Asignamos el orden secuencial (1, 2, 3...)
    for i in range(n):
        lista[i]["Seq. Order"] = i + 1

    return lista


def calcular_separaciones(lista_ordenada):
    """Columna I: Separació temporal respecte al vol precedent a l'IAF."""
    for i in range(len(lista_ordenada)):
        if i == 0:
            lista_ordenada[i]["Sep. to previous at IAF [s]"] = "-"
        else:
            diff = lista_ordenada[i]["_tiempo_s"] - lista_ordenada[i - 1]["_tiempo_s"]
            lista_ordenada[i]["Sep. to previous at IAF [s]"] = round(diff, 1)
    return lista_ordenada


# ==============================================================================
# CÀLCULS
# ==============================================================================
filas_calculadas = []

for fila in tabla_inicial:
    nombre_real = mapeo_modelos[fila["Aircraft model"]]
    avion_actual = next(av for av in lista_aviones if av.nombre == nombre_real)

    pct = fila["%MLW IAF"]
    d_km = fila["Dist. STAR entry WP to IAF [km]"]

    # Simulació CDO
    x_sim, h_sim, t_sim, _ = simular_cdo(avion_actual, pct)

    alt_entry = calcular_alt_at_star_entry_wp_m(x_sim, h_sim, d_km)
    tiempo_vuelo = calcular_time_from_star_entry_wp_to_iaf_s(x_sim, t_sim, d_km)
    toa = calcular_toa_iaf(tiempo_vuelo)

    filas_calculadas.append({
        "Aircraft model": fila["Aircraft model"],
        "STAR": fila["STAR"],
        "%MLW IAF": pct,
        "Dist. STAR entry WP to IAF [km]": d_km,
        "Alt. at STAR entry WP [m]": alt_entry,
        "Time from STAR entry WP to IAF [s]": tiempo_vuelo,
        "_tiempo_s": tiempo_vuelo,
        "ToA IAF": toa
    })

filas_ordenadas = ordenar_vuelos(filas_calculadas)
filas_finales = calcular_separaciones(filas_ordenadas)

# ==============================================================================
# SORTIDA PER CONSOLA
# ==============================================================================
print("\n" + "=" * 115)
print(f"{'Aircraft model':<16} {'STAR':<10} {'%MLW':<6} {'Dist [km]':<12} {'Alt [m]':<10} {'Time [s]':<10} {'ToA IAF':<10} {'Seq':<5} {'Sep [s]'}")
print("-" * 115)

for r in filas_finales:
    print(
        f"{r['Aircraft model']:<16} "
        f"{r['STAR']:<10} "
        f"{r['%MLW IAF']:<6} "
        f"{r['Dist. STAR entry WP to IAF [km]']:<12.4f} "
        f"{r['Alt. at STAR entry WP [m]']:<10.1f} "
        f"{r['Time from STAR entry WP to IAF [s]']:<10.1f} "
        f"{r['ToA IAF']:<10} "
        f"{r['Seq. Order']:<5} "
        f"{str(r['Sep. to previous at IAF [s]']):<7}"
    )

print("=" * 115 + "\n")