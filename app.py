import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- CONFIGURACIÓN DE LA PÁGINA ---
st.set_page_config(page_title="Band Diagram Visualizer", layout="wide")

st.title("Interactive Band Diagram: Metal-Semiconductor Interface")
st.markdown("""
**Author:** Gonzalo Alba | **University of Cádiz** This tool visualizes the Fermi level alignment and band bending in a Metal-Semiconductor junction 
based on the Schottky-Mott theory.
""")

# --- BARRA LATERAL (CONTROLES) ---
st.sidebar.header("Parameters")

# Sliders definidos según tu prompt
T = st.sidebar.slider("Temperature (T) [K]", min_value=100, max_value=600, value=300, step=10)
Na_exp = st.sidebar.slider("Acceptor Conc. (Na) [log cm^-3]", 14.0, 20.0, 16.0, 0.1)
Nd_exp = st.sidebar.slider("Donor Conc. (Nd) [log cm^-3]", 10.0, 18.0, 12.0, 0.1)
Wm = st.sidebar.slider("Metal Workfunction (Wm) [eV]", 3.0, 6.0, 4.5, 0.05)
Xs = st.sidebar.slider("SC Electron Affinity (Xs) [eV]", 3.0, 5.0, 4.05, 0.05)

# Convertir logs a valores lineales
Na = 10**Na_exp
Nd = 10**Nd_exp

# --- MOTOR FÍSICO (Tus Fórmulas de Desmos) ---
def calculate_physics():
    # Constantes
    q = 1.6e-19
    e0 = 8.85e-14
    e1 = 5.0         # Tu valor específico (posiblemente diamante o Si aproximado)
    k_eV = 8.617e-5  # k en eV
    k_J = 1.3806e-23 # k en J
    h = 6.626e-34    # Planck
    m_0 = 9.109e-31  # Masa electrón libre
    m = 0.81 * m_0   # Masa efectiva huecos (aprox Si/Diamante según tu input)
    g = 4.0          # Factor degeneración
    Eg = 1.12        # Bandgap (Si default, cámbialo a 5.47 si es diamante)
    Ea0 = 0.045      # Energía ionización base

    # Ionization energy dependence with T/Na
    # Fórmula: Ea = Ea0 - 5.5e-8 * Na^(1/3)
    Ea = Ea0 - 5.5e-8 * (Na**(1/3))

    # Effective density of states (Nv)
    # Tu fórmula: (2 * ( (2*pi*m*T*k_J) / h^2 )^1.5 ) / 1000000
    # Nota: El /1000000 es para pasar de m^-3 a cm^-3
    Nv = (2 * ((2 * np.pi * m * T * k_J) / (h**2))**1.5) / 1e6

    # Fermi Level (F)
    # Esta es la fórmula compleja que me pasaste
    # F representa la distancia del Nivel de Fermi a la Banda de Conducción (Ec - Ef)
    # o a la de Valencia, según tu sistema de referencia. 
    # En tu Desmos: y_VB = F - Eg. Si y=0 es Fermi, entonces F es la distancia Ec - Fermi.
    
    term1 = k_eV * T * np.log(Nv / (Na - Nd))
    
    # Término dentro de la raíz gigante
    exp_term = (g * np.exp(Ea / (k_eV * T)) * Nd) / Nv
    term_inside_sqrt = 1 + (4 * g * np.exp(Ea / (k_eV * T)) * Na) / Nv + \
                       exp_term * (exp_term - 2)
    
    term2_inner = 0.5 * (np.sqrt(term_inside_sqrt) + 1 + exp_term)
    term2 = k_eV * T * np.log(term2_inner)
    
    F = term1 + term2

    # Workfunction of SC
    Ws = Xs + F 
    
    # Built-in voltage
    Vbi = Ws - Wm
    
    # Schottky barrier height
    Sbh = F + Vbi

    # Depletion width (w) [nm]
    # Tu fórmula: sqrt( (2*|Vbi|*e0*e1) / (q*Na) ) * 10^7
    # El factor 10^7 convierte cm a nm
    w = np.sqrt((2 * abs(Vbi) * e0 * e1) / (q * Na)) * 1e7

    return F, Ws, Vbi, Sbh, w, q, Na, e0, e1, Eg, Xs

# Ejecutar cálculos
F, Ws, Vbi, Sbh, w, q, Na, e0, e1, Eg, Xs = calculate_physics()

# --- GENERACIÓN DE DATOS PARA GRÁFICAS (x, y) ---
# Definir eje X (distancia en nm)
x_metal = np.linspace(-10, 0, 100)
x_sc_dep = np.linspace(0, w, 200) if w > 0 else np.array([0])
x_sc_bulk = np.linspace(w, w + 20, 100)

# Función de Potencial V(z) en la zona de carga espacial
# Tu fórmula: V(z) = (q*Na / (2*e0*e1)) * (z-w)^2 * (1/10^14)
# El 1/10^14 ajusta unidades de nm^2 a cm^2
def V_z(z, w, Vbi):
    potential = (q * Na / (2 * e0 * e1)) * ((z - w)**2) * (1e-14)
    if Vbi < 0:
        return potential # Caso acumulación (simplificado según tu Desmos)
    else:
        return potential # Caso depleción

# Calcular perfiles de energía
# BANDA DE CONDUCCIÓN (Ec)
# Bulk: y = F
# Depletion: -V(z) + F (si Vbi > 0)
Ec_bulk = np.full_like(x_sc_bulk, F)
v_pot = V_z(x_sc_dep, w, Vbi)

if Vbi >= 0:
    Ec_dep = -v_pot + F
else:
    Ec_dep = v_pot + F # Invertir curvatura si Vbi < 0

# BANDA DE VALENCIA (Ev)
# Bulk: y = F - Eg
Ev_bulk = np.full_like(x_sc_bulk, F - Eg)
if Vbi >= 0:
    Ev_dep = -v_pot + F - Eg
else:
    Ev_dep = v_pot + F - Eg

# NIVEL DE VACÍO (Evac)
# Bulk: y = F + Xs
Evac_bulk = np.full_like(x_sc_bulk, F + Xs)
if Vbi >= 0:
    Evac_dep = -v_pot + F + Xs
else:
    Evac_dep = v_pot + F + Xs

# METAL
# Fermi level metal = 0 (Referencia)
# Metal Workfunction level = Wm (en el vacío)
y_metal_fermi = np.zeros_like(x_metal)
y_metal_vac = np.full_like(x_metal, Wm) # Nivel de vacío en el metal está a Wm sobre Fermi

# --- VISUALIZACIÓN CON PLOTLY ---
fig = go.Figure()

# 1. Metal
fig.add_trace(go.Scatter(x=x_metal, y=y_metal_fermi, mode='lines', name='Metal Fermi (Ef)', line=dict(color='gold', width=3, dash='dash')))
fig.add_trace(go.Scatter(x=x_metal, y=y_metal_vac, mode='lines', name='Metal Vacuum', line=dict(color='black', width=2)))
# Relleno visual metal
fig.add_trace(go.Scatter(x=x_metal, y=np.full_like(x_metal, -2), fill='tonexty', mode='none', showlegend=False, fillcolor='rgba(255, 215, 0, 0.1)'))

# 2. Semiconductor (Conduction Band)
x_full_sc = np.concatenate([x_sc_dep, x_sc_bulk])
y_Ec = np.concatenate([Ec_dep, Ec_bulk])
fig.add_trace(go.Scatter(x=x_full_sc, y=y_Ec, mode='lines', name='Conduction Band (Ec)', line=dict(color='red', width=3)))

# 3. Semiconductor (Valence Band)
y_Ev = np.concatenate([Ev_dep, Ev_bulk])
fig.add_trace(go.Scatter(x=x_full_sc, y=y_Ev, mode='lines', name='Valence Band (Ev)', line=dict(color='blue', width=3)))

# 4. Semiconductor (Vacuum Level)
y_Evac = np.concatenate([Evac_dep, Evac_bulk])
fig.add_trace(go.Scatter(x=x_full_sc, y=y_Evac, mode='lines', name='Vacuum Level (Evac)', line=dict(color='gray', dash='dot')))

# 5. Semiconductor (Fermi Level)
fig.add_trace(go.Scatter(x=x_full_sc, y=np.zeros_like(x_full_sc), mode='lines', name='Semi Fermi (Ef)', line=dict(color='green', width=2, dash='dash')))

# Layout del gráfico
fig.update_layout(
    title=f"Energy Band Diagram (Vbi = {Vbi:.2f} eV, Depletion w = {w:.1f} nm)",
    xaxis_title="Position z [nm] (0 = Interface)",
    yaxis_title="Energy [eV] (Ref: Fermi Level = 0)",
    yaxis=dict(range=[-Eg-1, Xs+F+1]), # Auto-zoom razonable
    xaxis=dict(range=[-5, w+15]),
    template="plotly_white",
    height=600
)

# Añadir anotaciones de valores clave
st.plotly_chart(fig, use_container_width=True)

# --- MOSTRAR CÁLCULOS INTERMEDIOS (EDUCATIVO) ---
with st.expander("See Physics Calculations (Math Details)"):
    col1, col2 = st.columns(2)
    with col1:
        st.latex(r"N_v = " + f"{Nv:.2e}" + r" \text{ cm}^{-3}")
        st.latex(r"E_a = " + f"{Ea*1000:.1f}" + r" \text{ meV}")
        st.latex(r"F (E_c - E_F) = " + f"{F:.3f}" + r" \text{ eV}")
    with col2:
        st.latex(r"W_s (\text{Workfunc SC}) = " + f"{Ws:.3f}" + r" \text{ eV}")
        st.latex(r"V_{bi} (\text{Built-in}) = " + f"{Vbi:.3f}" + r" \text{ eV}")
        st.latex(r"w (\text{Depletion}) = " + f"{w:.2f}" + r" \text{ nm}")

# --- EXPORTAR DATOS ---
# Permitir al alumno descargar los datos para hacer su informe
import pandas as pd
df = pd.DataFrame({
    "Position_nm": x_full_sc,
    "Ec_eV": y_Ec,
    "Ev_eV": y_Ev,
    "Evac_eV": y_Evac
})
csv = df.to_csv(index=False).encode('utf-8')
st.download_button("Download Data as CSV", data=csv, file_name="band_diagram_data.csv", mime="text/csv")
