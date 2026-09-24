import math

# ==========================================
# BASE DE DATOS DE AVIONES (Tabla Anexo B)
# Copiada exactamente en sus unidades originales
# ==========================================
aviones={
    'B767-300ER': {
        'MLW_tons': 0.145150E+03,
        'Max_Weight_tons': 0.20410E+03,
        'Max_Payload_tons': 0.46500E+02,
        'S': 0.28350E+03,
        'CD0_app': 0.14000E-01,
        'CD2_app': 0.49000E-01,
        'CD0_clean': 0.17400E-01,
        'CD2_clean': 0.45900E-01,
        'hp_desc': 26418,
        'CT_high': 0.64359E-1,
        'CT_low': 0.55988E-1,
        'CT_app': 0.12475,
        'CT1': .35167E+06,
        'CT2': .44673E+05,
        'CT3': .10129E-09,
        'CF1': .54005E+00,
        'CF2': .55782E+03
    },
    'B777-300': {
        'MLW_tons': 0.237680E+03,
        'Max_Weight_tons': 0.29930E+03,
        'Max_Payload_tons': 0.64900E+02,
        'S': 0.42804E+03,
        'CD0_app': 0.17300E-01,
        'CD2_app': 0.48400E-01,
        'CD0_clean': 0.15700E-01,
        'CD2_clean': 0.42000E-01,
        'hp_desc': 36122,
        'CT_high': 0.44239E-1,
        'CT_low': 0.41065E-1,
        'CT_app': 0.92921E-1,
        'CT1': .42577E+06,
        'CT2': .48987E+05,
        'CT3': .66146E-10,
        'CF1': .87843E+00,
        'CF2': .36897E+04
    },
    'B737': {
        'MLW_tons': 0.51710E+02,
        'Max_Weight_tons': 0.70800E+02,
        'Max_Payload_tons': 0.16920E+02,
        'S': 0.12465E+03,
        'CD0_app': 0.27000E-01,
        'CD2_app': 0.44100E-01,
        'CD0_clean': 0.23500E-01,
        'CD2_clean': 0.44500E-01,
        'hp_desc': 30152,
        'CT_high': 0.36336E-1,
        'CT_low': 0.53395E-1,
        'CT_app': 0.16440,
        'CT1': .14573E+06,
        'CT2': .55638E+05,
        'CT3': .14200E-10,
        'CF1': .94680E+00,
        'CF2': .10000E+15
    },
    'A320-212': {
        'MLW_tons': 0.64500E+02,
        'Max_Weight_tons': 0.77000E+02,
        'Max_Payload_tons': 0.21500E+02,
        'S': 0.12260E+03,
        'CD0_app': 0.24200E-01,
        'CD2_app': 0.46900E-01,
        'CD0_clean': 0.24000E-01,
        'CD2_clean': 0.37500E-01,
        'hp_desc': 12398,
        'CT_high': 0.45711E-1,
        'CT_low': 0.27207E-1,
        'CT_app': 0.13981,
        'CT1': .13605E+06,
        'CT2': .52238E+05,
        'CT3': .26637E-10,
        'CF1': .94000E+00,
        'CF2': .10000E+06
    },
    'A319-131': {
        'MLW_tons': 0.61000E2,
        'Max_Weight_tons': 0.70000E+02,
        'Max_Payload_tons': 0.17000E+02,
        'S': 0.12260E+03,
        'CD0_app': 0.28400E-01,
        'CD2_app': 0.37600E-01,
        'CD0_clean': 0.28000E-01,
        'CD2_clean': 0.31000E-01,
        'hp_desc': 27726,
        'CT_high': 0.83084E-1,
        'CT_low': 0.51765E-1,
        'CT_app': 0.14767,
        'CT1': .13900E+06,
        'CT2': .58900E+05,
        'CT3': .57200E-14,
        'CF1': .68800E+00,
        'CF2': .16700E+04
    }
}


def atmosfera_isa(altitud_ft):
    """
    Calcula la densitat de l'aire (rho) segons el model ISA a una altitud dada en peus (ft).[cite: 1]
    """
    h_m = altitud_ft * 0.3048  # Convertir peus a metres

    # Constants ISA
    T0 = 288.15  # Temperatura a nivell del mar (K)
    L = 0.0065  # Gradient tèrmic (K/m)
    rho0 = 1.225  # Densitat a nivell del mar (kg/m^3)
    g = 9.80665  # Gravetat (m/s^2)
    R = 287.058  # Constant dels gasos (J/(kg·K))

    T = T0 - L * h_m
    rho = rho0 * (T / T0) ** ((g / (L * R)) - 1)

    return rho


def calcular_v_minrod(avion, peso_kg, rho, empuje_T):
    """
    Calcula la velocitat que minimitza el Rate of Descent (ROD)[cite: 3]
    """
    S = avion['S']
    CD0 = avion['CD0_clean']  # Asumim configuració clean
    CD2 = avion['CD2_clean']
    g = 9.80665

    # Aplicem la fórmula demostrada
    numerador = empuje_T + math.sqrt(empuje_T ** 2 + 12 * CD0 * CD2 * (peso_kg * g) ** 2)
    denominador = 3 * CD0 * rho * S

    v_minrod = math.sqrt(numerador / denominador)  # Velocitat en m/s
    return v_minrod


# ==============================================================================
# 3. EMPUJE EN RALENTÍ (IDLE THRUST)
# ==============================================================================
def calcular_empuje_idle(avion, hp_ft, config):
    """
    Calcula el empuje al ralentí T_desc (en Newtons) según las ecuaciones BADA[cite: 1]
    """
    CT1 = avion['CT1']
    CT2 = avion['CT2']
    CT3 = avion['CT3']

    # Empuje máximo a la altitud actual[cite: 1]
    T_max = CT1 * (1.0 - (hp_ft / CT2) + (CT3 * (hp_ft ** 2)))

    # Empuje en ralentí según altitud y configuración[cite: 1]
    if hp_ft > avion['hp_desc']:
        T_desc = avion['CT_high'] * T_max
    else:
        if config == 'clean':
            T_desc = avion['CT_low'] * T_max
        else:  # config == 'approach'
            T_desc = avion['CT_app'] * T_max

    return T_desc


# ==============================================================================
# 4. SIMULADOR DE TRAYECTORIA CDO
# ==============================================================================
def simular_cdo(nombre_avion, porcentaje_peso, altitud_inicial_ft=12000):
    """
    Simula el descenso continuo desde altitud_inicial_ft hasta 5,000 ft[cite: 1]
    """
    avion = aviones[nombre_avion]

    # Conversión: Toneladas de la tabla -> kilogramos dentro del cálculo
    peso_kg = (porcentaje_peso / 100.0) * (avion['MLW_tons'] * 1000.0)
    g = 9.80665

    # Variables iniciales
    h_ft = float(altitud_inicial_ft)
    t = 0.0  # Tiempo acumulado (segundos)
    dt = 1.0  # Paso del bucle (1 segundo)
    distancia_m = 0.0  # Distancia horizontal recorrida (metros)

    # Historial para guardar datos
    historial_t = []
    historial_h = []

    # BUCLE DE SIMULACIÓN HASTA LLEGAR AL IAF (6000[cite: 1]
    while h_ft > 6000:
        # Selección de clave según la altitud (Transición a Approach a 6000 ft)[cite: 1]
        config_key = 'app' if h_ft <= 6000 else 'clean'
        config_idle = 'approach' if h_ft <= 6000 else 'clean'

        # 1. Densidad del aire[cite: 1]
        rho = atmosfera_isa(h_ft)

        # 2. Empuje al ralentí[cite: 1]
        T = calcular_empuje_idle(avion, h_ft, config_idle)

        # 3. Datos aerodinámicos
        S = avion['S']
        CD0 = avion[f'CD0_{config_key}']
        CD2 = avion[f'CD2_{config_key}']

        # 4. Cálculo de V_minROD (Fórmula analítica derivada)[cite: 3]
        num = T + math.sqrt(T ** 2 + 12 * CD0 * CD2 * (peso_kg * g) ** 2)
        den = 3 * CD0 * rho * S
        v_m_s = math.sqrt(num / den)

        # 5. Resistencia (D) y Rate of Descent (ROD)[cite: 1, 3]
        CL = (2 * peso_kg * g) / (rho * (v_m_s ** 2) * S)
        CD = CD0 + CD2 * (CL ** 2)
        D = 0.5 * rho * (v_m_s ** 2) * S * CD

        rod_m_s = ((D - T) * v_m_s) / (peso_kg * g)
        rod_ft_s = rod_m_s / 0.3048  # Convertir m/s a ft/s

        # Guardar distancia horizontal recorrida
        distancia_m += v_m_s * dt

        # Guardar datos en el historial
        historial_t.append(t)
        historial_h.append(h_ft)

        # Actualizar valores para el siguiente segundo
        h_ft -= rod_ft_s * dt
        t += dt

    distancia_NM = distancia_m / 1852.0  # Convertir metros a Millas Náuticas (NM)
    return t, distancia_NM, historial_t, historial_h


# ==============================================================================
# 5. BLOQUE DE PRUEBA (EJECUCIÓN)
# ==============================================================================
if __name__ == '__main__':
    print("--- RESULTADOS DE SIMULACIÓN DESDE 12.000 ft HASTA6000ft6 ft ---")

    # Simulamos todos los aviones al 100% de su MLW
    for nombre in aviones:
        tiempo_s, dist_nm, t_hist, h_hist = simular_cdo(nombre, porcentaje_peso=100, altitud_inicial_ft=12000)
        minutos = tiempo_s / 60.0
        print(f"Avión {nombre:10s} -> Tiempo: {minutos:5.2f} min ({tiempo_s:3.0f} s) | Distancia: {dist_nm:5.2f} NM")