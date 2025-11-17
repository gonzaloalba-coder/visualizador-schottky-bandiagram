import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- CONFIGURACIÓN DE LA PÁGINA ---
st.set_page_config(page_title="Band Diagram Visualizer", layout="wide")

st.title("Interactive Band Diagram: Metal-Semiconductor Interface")
st.markdown("""
**Author:** Gonzalo Alba | **University of Cádiz** This tool visualizes the Fermi level alignment and band bending in a **p-type Semiconductor** junction 
considering incomplete ionization (based on user-provided Desmos model).
""")

# --- BARRA LATERAL (CONTROLES) ---
st.sidebar.header("Parameters")

# Sliders
T = st.sidebar.slider("Temperature (T) [K]", 100, 600, 300, 10)
Na_exp = st.sidebar.slider("Acceptor Conc. (Na) [log cm^-3]", 14.0, 20.0, 17.0, 0.1)
Nd_exp = st.sidebar.slider("Donor Conc. (Nd) [log cm^-3]", 10.0, 18.0, 12.0, 0.1)
Wm = st.sidebar.slider("Metal Workfunction (Wm) [eV]", 3.0, 6.0, 4.5, 0.05)
Xs = st.sidebar.slider("SC Electron Affinity (Xs) [eV]", 3.0, 5.0, 4.05, 0.05)

# Convertir logs a valores lineales
Na = 10**Na_exp
Nd = 10**Nd_exp

# --- ADVERTENCIA DE SEGURIDAD FÍSICA ---
if Nd >= Na:
    st.error(f"⚠️ Math Error: Your formula assumes a p-type semiconductor where Na > Nd. \nCurrent: Na={Na:.1e}, Nd={Nd:.1e}. Please increase Na or decrease Nd.")
    st.stop() # Detiene la ejecución para no romper el código

# --- MOTOR FÍSICO ---
def calculate_physics():
    try:
        # Constantes
        q = 1.6e-19
        e0 = 8.85e-14
        e1 = 5.0         # Tu valor
        k_eV = 8.617e-5  # k en eV
        k_J = 1.3806e-23 # k en J
        h = 6.626e-34    # Planck
        m_0 = 9.109e-31 
        m = 0.81 * m_0   # Masa efectiva huecos
        g = 4.0          # Factor degeneración
        Eg = 1.12        # Bandgap (Si) - Ajustar si es Diamante (5.47)
        Ea0 = 0.045      # Energía ionización base

        # 1. Ionization energy (Ea)
        Ea = Ea0 - 5.5e-8 * (Na**(1/3))

        # 2. Effective density of states (Nv)
        # Factor 1e6 convierte m^-3 a cm^-3
        Nv = (2 * ((2 * np.pi * m * T * k_J) / (h**2))**1.5) / 1e6

        # 3. Fermi Level (F) - The Complex Formula
        # F = kT * ln( Nv / (Na - Nd) ) + Correction_Term
        
        # Pre-check para evitar log(negativo)
        net_doping = Na - Nd
        if net_doping <= 0: return None # Safety catch

        term1 = k_eV * T * np.log(Nv / net_doping)
        
        # Término exponencial y raíz
        # Dividimos por k_eV*T. OJO: Si T es muy baja, esto explota (Overflow).
        # Python aguanta float hasta ~1e308.
        kt_ev = k_eV * T
        exp_factor = np.exp(Ea / kt_ev) 
        
        term_A = (g * exp_factor * Nd) / Nv
        term_B = (4 * g * exp_factor * Na) / Nv
        
        term_inside_sqrt = 1 + term_B + term_A * (term_A - 2)
        
        # Evitar raíz de negativo
        if term_inside_sqrt < 0: term_inside_sqrt = 0 

        term2_inner = 0.5 * (np.sqrt(term_inside_sqrt) + 1 + term_A)
        term2 = kt_ev * np.log(term2_inner)
        
        F = term1 + term2

        # 4. Parámetros derivados
        Ws = Xs + F 
        Vbi = Ws - Wm
        Sbh = F + Vbi
        
        # Depletion width (w)
        # w = sqrt( 2 * |Vbi| * e0 * e1 / q*Na ) * 10^7 (para nm)
        w = np.sqrt((2 * abs(Vbi) * e0 * e1) / (q * Na)) * 1e7

        return F, Ws, Vbi, Sbh, w, q, e0, e1, Eg, Xs
        
    except Exception as e:
        st.error(f"Calculation Error: {e}")
        return None

# Ejecutar cálculos
results = calculate_physics()

if results is None:
    st.stop() # Si hubo error matemático, paramos aquí.

F, Ws, Vbi, Sbh, w, q, e0, e1, Eg, Xs = results

# --- GENERACIÓN DE DATOS (Vectores) ---
# Eje Z (profundidad) en nm
x_metal = np.linspace(-10, 0, 100)

# Depletion Region (0 a w)
if w > 0.1: # Si w es muy pequeño, crea un punto al menos
    x_sc_dep = np.linspace(0, w, 200)
else:
    x_sc_dep = np.array([0, 0.1])

# Bulk Region (w a w+20)
x_sc_bulk = np.linspace(x_sc_dep[-1], x_sc_dep[-1] + 20, 100)

# Potencial V(z)
# V(z) = (q*Na / 2*e0*e1) * (z-w)^2 * conversion
# z en nm. Factor conversión de la fórmula original: 1e-14
def get_potential(z_nm, w_nm, vbi_val):
    # Constante pre-factor
    prefactor = (q * Na) / (2 * e0 * e1) 
    # Ajuste unidades: z^2 es nm^2. Necesitamos cm^2 para cancelar con q*Na (cm^-3)? 
    # Tu fórmula original tenía un factor 1/10^14 manual. Lo mantenemos.
    pot = prefactor * ((z_nm - w_nm)**2) * 1e-14
    return pot

v_pot = get_potential(x_sc_dep, w, Vbi)

# --- PERFILES DE ENERGÍA ---
# Nota sobre signos:
# Si Vbi > 0 (Depletion típica p-type con metal baja función trabajo?): Bandas se doblan abajo?
# Tu formula Desmos: V01(z) = -V(z) + F - Eg  (Para VB)
# Vamos a implementar tu lógica de Desmos tal cual:

# BANDA VALENCIA (Ev)
Ev_bulk_val = F - Eg
Ev_bulk = np.full_like(x_sc_bulk, Ev_bulk_val)

if Vbi >= 0:
    Ev_dep = -v_pot + F - Eg
else:
    Ev_dep = v_pot + F - Eg

# BANDA CONDUCCIÓN (Ec)
Ec_bulk_val = F
Ec_bulk = np.full_like(x_sc_bulk, Ec_bulk_val)

if Vbi >= 0:
    Ec_dep = -v_pot + F
else:
    Ec_dep = v_pot + F

# NIVEL VACÍO (Evac)
Evac_bulk_val = F + Xs
Evac_bulk = np.full_like(x_sc_bulk, Evac_bulk_val)

if Vbi >= 0:
    Evac_dep = -v_pot + F + Xs
else:
    Evac_dep = v_pot + F + Xs

# Concatenar
x_semi = np.concatenate([x_sc_dep, x_sc_bulk])
y_Ev = np.concatenate([Ev_dep, Ev_bulk])
y_Ec = np.concatenate([Ec_dep, Ec_bulk])
y_Evac = np.concatenate([Evac_dep, Evac_bulk])
y_Ef_semi = np.zeros_like(x_semi) # Fermi en 0

# METAL
y_Ef_metal = np.zeros_like(x_metal)
y_Evac_metal = np.full_like(x_metal, Wm) # Nivel vacío metal = Wm (sobre Fermi)

# --- GRÁFICA ---
fig = go.Figure()

# Metal
fig.add_trace(go.Scatter(x=x_metal, y=y_Ef_metal, name="Metal Fermi", line=dict(color="gold", width=3, dash="dash")))
fig.add_trace(go.Scatter(x=x_metal, y=y_Evac_metal, name="Metal Vacuum", line=dict(color="black", width=2)))

# Semiconductor
fig.add_trace(go.Scatter(x=x_semi, y=y_Evac, name="SC Vacuum", line=dict(color="gray", dash="dot")))
fig.add_trace(go.Scatter(x=x_semi, y=y_Ec, name="Conduction Band (Ec)", line=dict(color="red", width=3)))
fig.add_trace(go.Scatter(x=x_semi, y=y_Ev, name="Valence Band (Ev)", line=dict(color="blue", width=3)))
fig.add_trace(go.Scatter(x=x_semi, y=y_Ef_semi, name="SC Fermi (Ef)", line=dict(color="green", width=2, dash="dash")))

# Layout
fig.update_layout(
    title=f"Band Diagram (Vbi={Vbi:.2f} eV, w={w:.1f} nm)",
    xaxis_title="Position z [nm]",
    yaxis_title="Energy [eV]",
    height=600,
    template="plotly_white",
    xaxis=dict(range=[-5, w + 10]),
    yaxis=dict(range=[np.min(y_Ev)-0.5, np.max(y_Evac)+1])
)

st.plotly_chart(fig, use_container_width=True)

# Información extra
st.info(f"**Parameters:** F (Ec-Ef) = {F:.3f} eV | Bandgap = {Eg} eV")
