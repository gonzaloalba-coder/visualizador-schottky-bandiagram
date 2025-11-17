import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- 1. CONFIGURACIÓN ---
st.set_page_config(page_title="PN Junction Solver", layout="wide")
st.title("Simulador de Unión P-N de Silicio (Método Exacto)")

# --- 2. PARÁMETROS (SIDEBAR) ---
st.sidebar.header("Parámetros del Dispositivo")

# Geometría
col1, col2 = st.sidebar.columns(2)
L_p = col1.number_input("Longitud Lado-P (nm)", value=200.0, step=10.0)
L_n = col2.number_input("Longitud Lado-N (nm)", value=200.0, step=10.0)

# Dopaje (Escala Logarítmica)
st.sidebar.subheader("Niveles de Dopaje")
Na_exp = st.sidebar.slider("Aceptores Na (p-side) [log cm^-3]", 14.0, 19.0, 17.0, 0.1)
Nd_exp = st.sidebar.slider("Donores Nd (n-side) [log cm^-3]", 14.0, 19.0, 16.0, 0.1)

# --- 3. CÁLCULOS FÍSICOS (CORE) ---
def calcular_bandas(na_log, nd_log, lp_nm, ln_nm):
    # 1. Definir Constantes (Silicio a 300K)
    Na = 10**na_log
    Nd = 10**nd_log
    
    q = 1.602e-19       # Carga elemental
    eps_0 = 8.854e-14   # F/cm
    eps_r = 11.7        # Permitividad relativa Si
    eps = eps_0 * eps_r
    k_T_eV = 0.02585    # kT en eV (T=300K)
    ni = 1.5e10         # Concentración intrínseca (cm^-3)
    Eg = 1.12           # Bandgap (eV)

    # 2. Posiciones de los niveles de energía LEJOS de la unión (Bulk)
    # Referencia: E_Fermi = 0 eV en todo el dispositivo
    
    # En lado P: p ~ Na.  p = ni * exp((Ei - Ef)/kT)  => Ei - Ef = kT * ln(Na/ni)
    # Como Ef=0, entonces Ei_p = kT * ln(Na/ni)
    # El Ei en el lado P está "arriba" (positivo)
    Ei_p_bulk = k_T_eV * np.log(Na / ni)
    
    # En lado N: n ~ Nd.  n = ni * exp((Ef - Ei)/kT)  => Ef - Ei = kT * ln(Nd/ni)
    # Como Ef=0, entonces Ei_n = -kT * ln(Nd/ni)
    # El Ei en el lado N está "abajo" (negativo)
    Ei_n_bulk = -k_T_eV * np.log(Nd / ni)
    
    # Potencial Integrado (Built-in) es la diferencia total
    Vbi = Ei_p_bulk - Ei_n_bulk

    # 3. Anchuras de Deplexión (xp y xn)
    # Fórmula estándar de la aproximación de deplexión
    # W = sqrt( 2*eps*Vbi/q * (Na+Nd)/(Na*Nd) )
    # Ojo: Vbi está en eV, para la fórmula de Poisson necesitamos Voltios (que numéricamente es igual)
    # Pero necesitamos q como carga (C) en el denominador si eps está en F/cm.
    
    # Simplificación para unidades: 
    # Factor curvatura = (q * N) / (2 * eps)  [V / cm^2]
    # Trabajaremos todo en cm y luego pasamos a nm
    
    W_cm = np.sqrt( (2 * eps * Vbi) / q * (Na + Nd) / (Na * Nd) )
    W_nm = W_cm * 1e7
    
    xn_nm = W_nm * Na / (Na + Nd)
    xp_nm = W_nm * Nd / (Na + Nd)

    # 4. Generar perfil de Banda Intrínseca Ei(x)
    # Creamos el eje X en nm
    x = np.linspace(-lp_nm, ln_nm, 600)
    Ei = np.zeros_like(x)

    # Factores de curvatura (A/2 * x^2)
    # Poisson: d2V/dx2 = -rho/eps.  Energía E = -qV.  d2E/dx2 = q*rho/eps.
    # Lado P (x < 0): rho = -qNa.  d2E/dx2 = -q^2 Na / eps (Curvatura Negativa/Convexa)
    # Lado N (x > 0): rho = +qNd.  d2E/dx2 = +q^2 Nd / eps (Curvatura Positiva/Cóncava)
    
    # Calculamos constantes en unidades compatibles con nm
    # K = (q^2 * N) / (2 * eps).  Unidades: eV / cm^2.
    # Para pasar a eV / nm^2 hay que dividir por 10^14
    K_p = (q * Na) / (2 * eps) * 1e-14 
    K_n = (q * Nd) / (2 * eps) * 1e-14

    for i, val_x in enumerate(x):
        if val_x < -xp_nm:
            # Zona Neutra P
            Ei[i] = Ei_p_bulk
        elif val_x < 0:
            # Zona Deplexión P (Parábola que baja)
            # E(x) = Ei_bulk - Kp * (x + xp)^2
            Ei[i] = Ei_p_bulk - K_p * (val_x + xp_nm)**2
        elif val_x < xn_nm:
            # Zona Deplexión N (Parábola que sube desde el fondo)
            # E(x) = Ei_bulk_n + Kn * (xn - x)^2
            Ei[i] = Ei_n_bulk + K_n * (xn_nm - val_x)**2
        else:
            # Zona Neutra N
            Ei[i] = Ei_n_bulk

    # 5. Calcular Bandas C y V a partir de Ei
    Ec = Ei + (Eg / 2)
    Ev = Ei - (Eg / 2)
    
    # Nivel de Fermi (Siempre 0)
    Ef = np.zeros_like(x)

    return x, Ec, Ev, Ef, Ei, xp_nm, xn_nm, Vbi, W_nm, Eg

# Ejecutar
x, Ec, Ev, Ef, Ei, xp, xn, Vbi, W, Eg = calcular_bandas(Na_exp, Nd_exp, L_p, L_n)

# --- 4. GRAFICADO (ESTILO SHADING) ---
fig = go.Figure()

# Definir límite inferior para el sombreado
y_min_plot = np.min(Ev) - 0.5
y_max_plot = np.max(Ec) + 0.5

# 1. Pista invisible para el suelo (Base del relleno)
fig.add_trace(go.Scatter(
    x=x, y=np.full_like(x, y_min_plot),
    mode='lines', line=dict(width=0), hoverinfo='skip', showlegend=False
))

# 2. Banda de Valencia (Rellena hasta abajo)
fig.add_trace(go.Scatter(
    x=x, y=Ev,
    name='Banda Valencia (Ev)',
    mode='lines',
    line=dict(color='#1f77b4', width=3),
    fill='tonexty', # Rellena hasta la traza anterior (el suelo)
    fillcolor='rgba(31, 119, 180, 0.3)'
))

# 3. Banda de Conducción
fig.add_trace(go.Scatter(
    x=x, y=Ec,
    name='Banda Conducción (Ec)',
    mode='lines',
    line=dict(color='#d62728', width=3)
))

# 4. Nivel de Fermi
fig.add_trace(go.Scatter(
    x=x, y=Ef,
    name='Nivel de Fermi (Ef)',
    mode='lines',
    line=dict(color='green', width=2, dash='dash')
))

# 5. Nivel Intríseco (Opcional, didáctico)
fig.add_trace(go.Scatter(
    x=x, y=Ei,
    name='Nivel Intríseco (Ei)',
    mode='lines',
    line=dict(color='gray', width=1, dash='dot'),
    visible='legendonly'
))

# Marcas verticales de la unión
fig.add_vline(x=0, line_width=1, line_color="black")
fig.add_vline(x=-xp, line_width=1, line_dash="dot", line_color="gray", annotation_text=f"-xp: {-xp:.1f} nm")
fig.add_vline(x=xn, line_width=1, line_dash="dot", line_color="gray", annotation_text=f"xn: {xn:.1f} nm")

# Layout
fig.update_layout(
    title=f"Diagrama de Bandas en Equilibrio (Vbi = {Vbi:.3f} eV)",
    xaxis_title="Posición (nm) [0 = Unión Metalúrgica]",
    yaxis_title="Energía (eV) [Ef = 0]",
    yaxis=dict(range=[y_min_plot, y_max_plot]),
    xaxis=dict(range=[-L_p, L_n]),
    template="plotly_white",
    height=600
)

st.plotly_chart(fig, use_container_width=True)

# Mostrar datos clave
c1, c2, c3 = st.columns(3)
c1.metric("Potencial Vbi", f"{Vbi:.3f} V")
c2.metric("Ancho Deplexión (W)", f"{W:.1f} nm")
c3.metric("Bandgap", f"{Eg} eV")
