import math
import matplotlib.pyplot as plt

# ==============================================================================
# CLASE: Avion
# ==============================================================================
class Avion:
    def __init__(self, nombre, mlw, s, cd0_limpio, cd2_limpio, cd0_app, cd2_app,
                 hp_desc, ct_alto, ct_bajo, ct_app, ct1, ct2, ct3, cf1, cf2):
        self.nombre = nombre          # Nombre identificador del modelo
        self.mlw = mlw                # Maximum Landing Weight [kg]
        self.s = s                    # Superficie alar [m^2]
        self.cd0_limpio = cd0_limpio  # Coeficiente de resistencia parásita en limpio
        self.cd2_limpio = cd2_limpio  # Coeficiente de resistencia inducida en limpio
        self.cd0_app = cd0_app        # Coeficiente de resistencia parásita en aproximación
        self.cd2_app = cd2_app        # Coeficiente de resistencia inducida en aproximación
        self.hp_desc = hp_desc        # Altitud de transición de empuje de descenso [ft]
        self.ct_alto = ct_alto        # Coeficiente de empuje idle en alta cota (hp > hp_desc)
        self.ct_bajo = ct_bajo        # Coeficiente de empuje idle en baja cota (limpio)
        self.ct_app = ct_app          # Coeficiente de empuje idle en configuración de aproximación
        self.ct1 = ct1                # Parámetro de empuje máximo CT1 [N]
        self.ct2 = ct2                # Parámetro de empuje máximo CT2 [ft]
        self.ct3 = ct3                # Parámetro de empuje máximo CT3 [1/ft^2]
        self.cf1 = cf1                # Coeficiente de consumo de combustible CF1 [kg/(min*kN)]
        self.cf2 = cf2                # Coeficiente de velocidad para consumo CF2 [kt]

# ==============================================================================
# CONSTANTES FÍSICAS, DE ATMÓSFERA ISA Y FACTORES DE CONVERSIÓN
# ==============================================================================
GRAVEDAD = 9.80665        # Aceleración de la gravedad g0 [m/s^2]
R_AIRE = 287.05287        # Constante del gas para el aire [J/(kg*K)]
T0_ISA = 288.15           # Temperatura estándar al nivel del mar [K] (15 °C)
P0_ISA = 101325.0         # Presión atmosférica al nivel del mar [Pa]
GRADIENTE_TEMP = -0.0065  # Variación de temperatura con la altura beta [K/m]
H_TROPOPAUSA = 11000.0    # Altura donde empieza la tropopausa [m]

PIES_A_METROS = 0.3048    # Multiplicar pies por esto para tener metros
METROS_A_PIES = 1.0 / PIES_A_METROS # Multiplicar metros por esto para tener pies
MS_A_NUDOS = 1.94384      # Conversión de m/s a nudos (kt)


# ==============================================================================
# FUNCIONES AUXILIARES
# ==============================================================================
def obtener_densidad_isa(h_m):
    """Calcula la densidad del aire rho en función de la altura usando las constantes ISA."""
    if h_m <= H_TROPOPAUSA:
        T = T0_ISA + GRADIENTE_TEMP * h_m
        exponente = -GRAVEDAD / (GRADIENTE_TEMP * R_AIRE)
        p = P0_ISA * (T / T0_ISA) ** exponente
    else:
        T_trop = T0_ISA + GRADIENTE_TEMP * H_TROPOPAUSA
        exponente = -GRAVEDAD / (GRADIENTE_TEMP * R_AIRE)
        p_trop = P0_ISA * (T_trop / T0_ISA) ** exponente
        T = T_trop
        p = p_trop * math.exp(-GRAVEDAD * (h_m - H_TROPOPAUSA) / (R_AIRE * T))

    return p / (R_AIRE * T)


def obtener_empuje_idle(avion, hp_ft, es_aproximacion):
    """Calcula el empuje idle (T_desc) en Newtons a partir de las formulas BADA."""
    # Formula de T_max corregida (con el signo menos)
    t_max = avion.ct1 * (1.0 - (hp_ft / avion.ct2) + avion.ct3 * (hp_ft ** 2))
    t_max = max(t_max, 0.0)

    # Seleccion del factor reductor CT_desc
    if hp_ft > avion.hp_desc:
        ct = avion.ct_alto
    else:
        ct = avion.ct_app if es_aproximacion else avion.ct_bajo

    return ct * t_max


# ==============================================================================
# FUNCIÓN DE SIMULACIÓN CDO (HACIA ATRÁS DESDE EL IAF)
# ==============================================================================
def simular_cdo(avion, porcentaje_mlw, h_max_m=12200.0, dt=1.0):
    # Condiciones iniciales en IAF (x = 0, h = 6000 ft) corregida con '*'
    masa = (porcentaje_mlw / 100.0) * avion.mlw
    h = 6000.0 * PIES_A_METROS

    x = 0.0
    t = 0.0

    lista_x = [x]
    lista_h = [h]
    lista_t = [t]
    lista_v = []

    while h < h_max_m:
        hp_ft = h * METROS_A_PIES
        es_aproximacion = (hp_ft <= 6000.0)

        # Seleccion de coeficientes segun si lleva flaps (<= 6000 ft) o limpio
        cd0 = avion.cd0_app if es_aproximacion else avion.cd0_limpio
        cd2 = avion.cd2_app if es_aproximacion else avion.cd2_limpio

        rho = obtener_densidad_isa(h)

        # Formula de v_minRoD corregida con '*'
        num = 4.0 * cd2 * ((masa * GRAVEDAD) ** 2)
        den = 3.0 * (rho ** 2) * (avion.s ** 2) * cd0
        v = (num / den) ** 0.25
        lista_v.append(v)

        # Formulas corregidas con '*'
        cl = (2.0 * masa * GRAVEDAD) / (rho * (v ** 2) * avion.s)
        cd = cd0 + cd2 * (cl ** 2)
        resistencia = 0.5 * rho * (v ** 2) * avion.s * cd
        empuje = obtener_empuje_idle(avion, hp_ft, es_aproximacion)

        # B. Formulas de RoD y gamma
        rod = v * (resistencia - empuje) / (masa * GRAVEDAD)
        if rod <= 0:
            rod = 0.1

        sin_gamma = max(min(rod / v, 1.0), -1.0)
        gamma = math.asin(sin_gamma)

        # C. Formula de flujo de combustible FF
        v_nudos = v * MS_A_NUDOS
        eta = avion.cf1 * (1.0 + v_nudos / avion.cf2)
        flujo_combustible = eta * (empuje / 1000.0) / 60.0

        # D. Integracion hacia atras en el tiempo
        h += rod * dt
        x -= v * math.cos(gamma) * dt
        masa += flujo_combustible * dt
        t += dt

        lista_x.append(x)
        lista_h.append(h)
        lista_t.append(t)

    lista_v.append(lista_v[-1])
    return lista_x, lista_h, lista_t, lista_v

# ==============================================================================
# BASE DE DATOS DE AVIONES
# ==============================================================================
lista_aviones = [
    Avion(
        nombre='B767-300ER',
        mlw=145150.0,
        s=283.50,
        cd0_limpio=0.0174,
        cd2_limpio=0.0459,
        cd0_app=0.0140,
        cd2_app=0.0490,
        hp_desc=26418.0,
        ct_alto=0.064359,
        ct_bajo=0.055988,
        ct_app=0.12475,
        ct1=351670.0,
        ct2=44673.0,
        ct3=0.10129e-09,
        cf1=0.54005,
        cf2=557.82
    ),
    Avion(
        nombre='B777-300',
        mlw=237680.0,
        s=428.04,
        cd0_limpio=0.0157,
        cd2_limpio=0.0420,
        cd0_app=0.0173,
        cd2_app=0.0484,
        hp_desc=36122.0,
        ct_alto=0.044239,
        ct_bajo=0.041065,
        ct_app=0.092921,
        ct1=425770.0,
        ct2=48987.0,
        ct3=0.66146e-10,
        cf1=0.87843,
        cf2=3689.7
    ),
    Avion(
        nombre='B737',
        mlw=51710.0,
        s=124.65,
        cd0_limpio=0.0235,
        cd2_limpio=0.0445,
        cd0_app=0.0270,
        cd2_app=0.0441,
        hp_desc=30152.0,
        ct_alto=0.036336,
        ct_bajo=0.053395,
        ct_app=0.16440,
        ct1=145730.0,
        ct2=55638.0,
        ct3=0.14200e-10,
        cf1=0.94680,
        cf2=1.0e14
    ),
    Avion(
        nombre='A320-212',
        mlw=64500.0,
        s=122.60,
        cd0_limpio=0.0240,
        cd2_limpio=0.0375,
        cd0_app=0.0242,
        cd2_app=0.0469,
        hp_desc=12398.0,
        ct_alto=0.045711,
        ct_bajo=0.027207,
        ct_app=0.13981,
        ct1=136050.0,
        ct2=52238.0,
        ct3=0.26637e-10,
        cf1=0.94000,
        cf2=1.0e5
    ),
    Avion(
        nombre='A319-131',
        mlw=61000.0,
        s=122.60,
        cd0_limpio=0.0280,
        cd2_limpio=0.0310,
        cd0_app=0.0284,
        cd2_app=0.0376,
        hp_desc=27726.0,
        ct_alto=0.083084,
        ct_bajo=0.051765,
        ct_app=0.14767,
        ct1=139000.0,
        ct2=58900.0,
        ct3=0.57200e-14,
        cf1=0.68800,
        cf2=1670.0
    )
]

# ==============================================================================
# SCRIPT PRINCIPAL: SÓLO DIBUJA AL EJECUTAR SIMULATOR DIRECTAMENTE
# ==============================================================================
if __name__ == '__main__':
    porcentajes_mlw = [100, 80]

    plt.figure(figsize=(11, 6))

    for avion in lista_aviones:
        for porcentaje in porcentajes_mlw:
            x_pts, h_pts, _, _ = simular_cdo(avion, porcentaje)
            plt.plot(x_pts, h_pts, label=f"{avion.nombre} [{porcentaje}% MLW]", linewidth=1.2)

    plt.title("Perfiles de Descenso Continuo (CDO) al IAF (x = 0 m, h = 6000 ft)")
    plt.xlabel("x [m]")
    plt.ylabel("h [m]")
    plt.xlim(-220000, 5000)
    plt.ylim(1000, 13000)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.show()