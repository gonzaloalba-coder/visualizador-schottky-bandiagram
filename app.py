import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- 1. CONFIGURACIÓN ---
st.set_page_config(page_title="PN Junction Exact Model", layout="wide")
st.title("Unión p-n de Silicio: Modelo de Doble Parábola")
st.markdown("Visualización exacta resolviendo la ecuación de Poisson por tramos (Aproximación de Deplexión).")

# --- 2. PARÁMETROS ---
st.sidebar.header("Propiedades del Material")

# Geometría finita
col1, col2 = st.sidebar.columns(2)
L_p = col1.number_input("Longitud P (nm)", value=150.0, step=10.0)
L_n = col2.number_input("Longitud N (nm)", value=150.0, step=10.0)

# Dopaje
st.sidebar.subheader("Niveles de Dopaje")
Na_exp = st.sidebar.slider("Na (Aceptores) [log cm^-3]", 15.0, 19.0, 17.0, 0.1)
Nd_exp = st.sidebar.slider("Nd (Donores) [log cm^-3]", 15.0, 19.0, 16.5, 0.1)

Na = 10**Na_exp
Nd = 10**Nd_exp

T = 300 # Kelvin

# --- 3. FÍSICA DE SEMICONDUCTORES (MOTOR) ---
def calculate_bands():
    # Constantes Físicas (Unidades CGS para cálculo interno)
    q = 1.602e-19      # Coulombs
    kb = 1.38e-23      # J/K
    kb_eV = 8.617e-5   # eV/K
    eps_0 = 8.854e-14  # F/cm
    eps_r = 11.7       # Silicio
    eps = eps_0 * eps_r
    
    # Parámetros del Silicio (300 K)
    Eg = 1.12          # eV
    Nc = 2.8e19        # cm^-3
    Nv = 1.04e19       # cm^-3
    ni = 1.5e10        # cm^-3

    # 1. Posición de Fermi en el Bulk (Lejos de la unión)
    # Referencia: E_Fermi = 0 eV en todo el dispositivo (equilibrio)
    
    # Lado P (Neutro): Ef - Ev = k*T * ln(Nv/Na)
    # Por tanto: Ev_bulk_p = - (Ef - Ev) = -k*T * ln(Nv/Na)
    delta_p = kb_eV * T * np.log(Nv / Na)
    Ev_p_bulk = -delta_p
    
    # Lado N (Neutro): Ec - Ef = k*T * ln(Nc/Nd)
    # Por tanto: Ev_bulk_n = Ec_bulk_n - Eg = (Ef + delta_n) - Eg
    delta_n = kb_eV * T * np.log(Nc / Nd)
    # Pero necesitamos la diferencia total de potencial (Vbi)
    
    # 2. Potencial Integrado (Built-in Potential)
    # Vbi = (kT/q) * ln(Na*Nd / ni^2)
    Vbi = kb_eV * T * np.log((Na * Nd) / (ni**2))
    
    # Comprobación de consistencia:
    # La diferencia entre las bandas de valencia debe ser Vbi
    Ev_n_bulk = Ev_p_bulk - Vbi

    # 3. Zona de Carga Espacial (Deplexión)
    # Ancho total W
    # Factor 1e7 para convertir cm a nm en el resultado final de W
    term_doping = (1/Na) + (1/Nd)
    W_cm = np.sqrt((2 * eps * Vbi) / q * term_doping)
    W_nm = W_cm * 1e7
    
    # Extensión hacia cada lado (xp y xn)
    xn_nm = W_nm * Na
