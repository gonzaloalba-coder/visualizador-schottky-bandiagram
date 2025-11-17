import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- 1. CONFIGURACIÓN INICIAL ---
st.set_page_config(page_title="Base p-n Junction", layout="wide")
st.title("Estructura Base: Unión p-n de Silicio (Finito)")

# --- 2. PARÁMETROS DE ENTRADA (GUI) ---
st.sidebar.header("Geometría y Dopaje")

# Definimos longitudes finitas (Eje X)
L_p = st.sidebar.number_input("Longitud lado-p (nm)", value=200.0, step=10.0)
L_n = st.sidebar.number_input("Longitud lado-n (nm)", value=200.0, step=10.0)

# Dopajes
Na_exp = st.sidebar.slider("Na (p-type) [log cm^-3]", 14.0, 19.0, 17.0)
Nd_exp = st.sidebar.slider("Nd (n-type) [log cm^-3]", 14.0, 19.0, 16.0)

Na = 10**Na_exp
Nd = 10**Nd_exp

# --- 3. MOTOR FÍSICO (AQUÍ IRÁN TUS DETALLES) ---
def calculate_pn_structure(Na, Nd, Lp, Ln):
    # --- [DETALLE A COMPLETAR]: Constantes del material (Si) ---
    # Pon aquí tus valores exactos si difieren
    Eg = 1.12        # Bandgap (eV)
    k_T = 0.0259     # kT a 300K (eV)
    ni = 1.5e10      # Intrinsic carrier conc (cm^-3)
    e_si = 11.7 * 8.85e-14 # Permitividad

    # 1. Cálculos de Potencial (Modelo de Deplexión)
    # Vbi: Potencial integrado
    try:
        Vbi = k_T * np.log((Na * Nd) / (ni**2))
    except:
        Vbi = 0.7 # Fallback seguro
        
    # W: Ancho de deplexión total
    # W = sqrt( 2*e*Vbi/q * (1/Na + 1/Nd) ) ... simplificado aquí para la estructura
    # Factor 1e7 para pasar a nm
    q = 1.6e-19
    W = np.sqrt((2 * e_si * Vbi * (Na + Nd)) / (q * Na * Nd)) * 1e7
    
    # Reparto de la zona de deplexión (xp hacia izquierda, xn hacia derecha)
    xn = W * Na / (Na + Nd)
    xp = W * Nd / (Na + Nd)

    # 2. Generación del Eje X (Finito)
    # Creamos un vector desde -Lp hasta +Ln
    # El 0 está en la unión metalúrgica
    x = np.linspace(-Lp, Ln, 500)
    
    # 3. Perfil de Bandas (Ec, Ev, Ef)
    Ec = np.zeros_like(x)
    Ev = np.zeros_like(x)
    
    # Nivel de Fermi (Referencia E_F = 0)
    Ef = np.zeros_like(x) 

    # Posición de las bandas lejos de la unión (Neutral regions)
    # Lado P (izq): Fermi cerca de Ev
    # delta_p = Ef - Ev = kT * ln(Nv/Na) -> Asumimos aproximación simple para estructura
    delta_p = k_T * np.log(1e19/Na) # Placeholder, ajusta con Nv real
    
    # Lado N (der): Fermi cerca de Ec
    # delta_n = Ec - Ef = kT * ln(Nc/Nd)
    delta_n = k_T * np.log(1e19/Nd) # Placeholder, ajusta con Nc real

    # Lógica de bandas (Simplificada para visualización base)
    # En p-neutral: Ev = -delta_p
    # En n-neutral: Ev = -delta_p - Vbi
    # La zona de transición se modela cuadráticamente en la realidad
    
    # --- [LÓGICA DE BUCLE PARA DIBUJAR LA CURVATURA] ---
    for i, val_x in enumerate(x):
        # Zona Neutra P (Izquierda de -xp)
        if val_x < -xp:
            band_bending = 0
        # Zona Deplexión P (Entre -xp y 0)
        elif val_x < 0:
            # Parábola cuadrática
            band_bending = Vbi * (Na / (Na+Nd)) * ((val_x + xp)/xp)**2
        # Zona Deplexión N (Entre 0 y xn)
        elif val_x < xn:
            # Parábola inversa
            band_bending = Vbi - (Vbi * (Nd / (Na+Nd)) * ((xn - val_x)/xn)**2)
        # Zona Neutra N (Derecha de xn)
        else:
            band_bending = Vbi

        # Asignamos energías (Referencia: Lado P plano)
        # Lado P bulk: Ev está a 'delta_p' debajo de Ef(0) si definimos Ef=0 global?
        # Mejor ref: Ev_bulk_p = 0 (arbitrario para visualizar)
        # Ec = Ev + Eg
        
        # Ajuste dinámico: Los niveles bajan qV(x)
        shift = -band_bending
        
        # Definición base lado P (antes de doblar)
        Ev_flat_p = -0.1 # Un poco bajo cero para ver el sombreado
        Ec_flat_p = Ev_flat_p + Eg
        
        Ev[i] = Ev_flat_p + shift
        Ec[i] = Ec_flat_p + shift

    return x, Ec, Ev, Ef, xp, xn, Vbi

# Ejecutamos cálculos
x, Ec, Ev, Ef, xp, xn, Vbi = calculate_pn_structure(Na, Nd, L_p, L_n)

# --- 4. VISUALIZACIÓN (SOMBRADO PERSONALIZADO) ---
fig = go.Figure()

# TRUCO DEL SOMBRADO:
# Para sombrear "hacia abajo" en Plotly sin que pinte todo el infinito,
# creamos una línea invisible en el fondo del gráfico (ej: -5 eV) y rellenamos hasta ella.

Y_BOTTOM_LIMIT = np.min(Ev) - 1.0 # Un poco más abajo del mínimo de Ev

# 1. Línea invisible de suelo (Base para el relleno)
fig.add_trace(go.Scatter(
    x=x, 
    y=np.full_like(x, Y_BOTTOM_LIMIT),
    mode='lines',
    line=dict(width=0),
    showlegend=False,
    hoverinfo='skip'
))

# 2. Banda de Valencia (Ev) - CON SOMBRADO
fig.add_trace(go.Scatter(
    x=x, y=Ev,
    mode='lines',
    name='Valence Band (Ev)',
    line=dict(color='blue', width=3),
    fill='tonexty', # Rellena hasta la traza anterior (la invisible de abajo)
    fillcolor='rgba(0, 0, 255, 0.2)' # Azul semitransparente
))

# 3. Banda de Conducción (Ec)
fig.add_trace(go.Scatter(
    x=x, y=Ec,
    mode='lines',
    name='Conduction Band (Ec)',
    line=dict(color='red', width=3)
))

# 4. Nivel de Fermi (Ef)
fig.add_trace(go.Scatter(
    x=x, y=Ef,
    mode='lines',
    name='Fermi Level (Ef)',
    line=dict(color='green', dash='dash')
))

# Layout
fig.update_layout(
    title=f"Energy Band Diagram (Vbi={Vbi:.2f} eV)",
    xaxis_title="Position (nm)",
    yaxis_title="Energy (eV)",
    xaxis=dict(range=[-L_p, L_n], zeroline=True, zerolinewidth=2, zerolinecolor='black'),
    yaxis=dict(range=[Y_BOTTOM_LIMIT, np.max(Ec)+0.5]),
    template="plotly_white",
    height=500
)

# Líneas verticales para marcar la zona de deplexión
fig.add_vline(x=-xp, line_dash="dot", line_color="gray", annotation_text="-xp")
fig.add_vline(x=xn, line_dash="dot", line_color="gray", annotation_text="xn")
fig.add_vline(x=0, line_width=1, line_color="black") # Intercara física

st.plotly_chart(fig, use_container_width=True)
