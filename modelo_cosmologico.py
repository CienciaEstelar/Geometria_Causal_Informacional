# ==========================================
# MODELO COSMOLÓGICO UNIFICADO: ALPs y Energía del Vacío
# ==========================================

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

try:
    from IPython.display import display
    IPYTHON_AVAILABLE = True
except ImportError:
    IPYTHON_AVAILABLE = False

print("=" * 60)
print("MODELO COSMOLÓGICO UNIFICADO - ALPs y Energía del Vacío")
print("=" * 60)

# ==========================================
# CONSTANTES FÍSICAS (valores numéricos actualizados 2025)
# ==========================================
print("\n1. CARGANDO CONSTANTES FÍSICAS...")

h_val = 6.62607015e-34        # Constante de Planck (J s)
hbar_val = h_val / (2 * np.pi) # Constante de Planck reducida (J s)
c_val = 299792458              # Velocidad de la luz (m/s)
G_val = 6.67430e-11           # Constante gravitacional (m³ kg⁻¹ s⁻²)
HO_val_SI = 70.4 * 1000 / 3.086e22  # Constante de Hubble (s⁻¹) [JWST 2025]
rho_vac_observada = 5.36e-10  # Densidad de energía del vacío (J/m³)

print(f"   • h = {h_val:.3e} J·s")
print(f"   • ħ = {hbar_val:.3e} J·s") 
print(f"   • c = {c_val:.3e} m/s")
print(f"   • G = {G_val:.3e} m³/kg·s²")
print(f"   • H₀ = {HO_val_SI:.3e} s⁻¹")
print(f"   • ρ_vac = {rho_vac_observada:.3e} J/m³")

# ==========================================
# PARTE I: DERIVACIÓN SIMBÓLICA DEL FORMALISMO
# ==========================================
print("\n2. DERIVANDO FORMALISMO MATEMÁTICO...")

t, G, xi = sp.symbols('t G xi', positive=True)
a = sp.Function('a')(t)
Phi = sp.Function('Phi')(t)
V = sp.Function('V')(Phi)
R_sym = sp.Symbol('R')

L_EH = R_sym / (16 * sp.pi * G)
L_Phi = -sp.Rational(1, 2) * sp.diff(Phi, t)**2 - V
L_int = -sp.Rational(1, 2) * xi * R_sym * Phi**2
L_total_general = L_EH + L_Phi + L_int

print("   • Lagrangiano general definido")

a_dot = sp.diff(a, t)
a_ddot = sp.diff(a, t, 2)
R_FLRW = 6 * (a_ddot / a + (a_dot / a)**2)

L_total_final = L_total_general.subs(R_sym, R_FLRW)
print("   • Lagrangiano con métrica FLRW obtenido")

eqs_mov = sp.calculus.euler.euler_equations(L_total_final, [Phi, a], t)
print("   • Ecuaciones de movimiento derivadas")

if IPYTHON_AVAILABLE:
    display(L_total_final)
    display(eqs_mov)

# ==========================================
# PARTE II: CONVERGENCIA EN ESCALA DE TRANSICIÓN
# ==========================================
print("\n3. ANALIZANDO CONVERGENCIA EN ω_c...")

omega_p, omega_h, omega_c = sp.symbols('omega_p omega_h omega_c', positive=True)
hbar_sym, c_sym = sp.symbols('hbar c', positive=True)

invariance_eq = sp.Eq(omega_c / omega_h, omega_p / omega_c)
sol_inv = sp.solve(invariance_eq, omega_c)
print("   • Invariancia de escala: ω_c = √(ω_p ω_h)")

rg_eq = sp.Eq(sp.log(omega_p) - sp.log(omega_c), sp.log(omega_c) - sp.log(omega_h))
sol_rg = sp.solve(rg_eq, omega_c)
print("   • Punto fijo RG: ω_c = √(ω_p ω_h)")

rho_vac_paper = hbar_sym * omega_c**4 / c_sym**3
rho_vac_fund = omega_h**2 * c_sym**2 / G
omega_p_sq_def = c_sym**5 / (hbar_sym * G)

consistency_eq = sp.Eq(rho_vac_paper, rho_vac_fund)
omega_c_4_sol = sp.solve(consistency_eq, omega_c**4)[0]
final_eq = sp.Eq(omega_c**4, omega_c_4_sol.subs(G, c_sym**5 / (hbar_sym * omega_p_sq_def)))
print("   • Consistencia de energía: ω_c⁴ = ω_p² ω_h²")

# ==========================================
# PARTE III: CALIBRACIÓN Y PREDICCIÓN FALSABLE
# ==========================================
print("\n4. REALIZANDO CALIBRACIÓN NUMÉRICA...")

omega_p_val = np.sqrt(c_val**5 / (hbar_val * G_val))  # Frecuencia de Planck
omega_h_val = HO_val_SI                              # Frecuencia de Hubble
omega_c_geom_val = np.sqrt(omega_p_val * omega_h_val) # Media geométrica

print(f"   • ω_p = {omega_p_val:.3e} rad/s")
print(f"   • ω_h = {omega_h_val:.3e} rad/s") 
print(f"   • ω_c (geom) = {omega_c_geom_val:.3e} rad/s")

kappa_val = ((2 * np.pi**2 * c_val**3 * rho_vac_observada) / 
             (3 * hbar_val * omega_c_geom_val**4))**(1/4)
nu_c_calib_val = kappa_val * (omega_c_geom_val / (2 * np.pi))

print(f"   • κ = {kappa_val:.4f}")
print(f"   • ν_c = {nu_c_calib_val:.4e} Hz")

masa_alp_kg = (h_val * nu_c_calib_val) / c_val**2
masa_alp_eV = (masa_alp_kg * c_val**2) / 1.60218e-19

print(f"   • m_ALP predicha = {masa_alp_eV:.2e} eV")

# ==========================================
# PARTE IV: SIMULACIÓN NUMÉRICA Y ESTABILIDAD
# ==========================================
print("\n5. EJECUTANDO SIMULACIÓN NUMÉRICA...")

solutions = sp.solve([eqs_mov[0], eqs_mov[1]], (a_ddot, sp.diff(Phi, t, 2)))
a_ddot_expr = solutions[a_ddot]
Phi_ddot_expr = solutions[sp.diff(Phi, t, 2)]

print("   • Aceleraciones despejadas")

m_sym = sp.symbols('m', positive=True)
V_specific = sp.Rational(1, 2) * m_sym**2 * Phi**2
dV_dPhi = sp.diff(V_specific, Phi)

a_ddot_expr_specific = a_ddot_expr.subs({V: V_specific, V.diff(Phi): dV_dPhi})
Phi_ddot_expr_specific = Phi_ddot_expr.subs({V: V_specific, V.diff(Phi): dV_dPhi})

state_vars = [a, a_dot, Phi, sp.diff(Phi, t)]
params = [G, xi, m_sym]

a_ddot_func = sp.lambdify(state_vars + params, a_ddot_expr_specific, modules=['numpy', 'sympy'])
Phi_ddot_func = sp.lambdify(state_vars + params, Phi_ddot_expr_specific, modules=['numpy', 'sympy'])

print("   • Funciones numéricas compiladas")

def derivs(t, y, xi_val, m_val):
    a_val, a_dot_val, Phi_val, Phi_dot_val = y
    
    a_ddot_val = a_ddot_func(a_val, a_dot_val, Phi_val, Phi_dot_val, 
                            G_val, xi_val, m_val)
    Phi_ddot_val = Phi_ddot_func(a_val, a_dot_val, Phi_val, Phi_dot_val, 
                                G_val, xi_val, m_val)
    
    return [a_dot_val, a_ddot_val, Phi_dot_val, Phi_ddot_val]

y0 = [1.0, 0.1, 0.1, 0.0]
t_span = [0, 1000]
args = (0.1, 1e-6)

sol = solve_ivp(derivs, t_span, y0, args=args, method='RK45', 
                dense_output=True, rtol=1e-8)

print("   • Integración numérica completada")

# ==========================================
# PARTE V: VISUALIZACIÓN DE RESULTADOS
# ==========================================
print("\n6. GENERANDO GRÁFICOS...")

t_eval = np.linspace(t_span[0], t_span[1], 5000)
y_sol = sol.sol(t_eval)

plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.plot(t_eval, y_sol[0], 'b-', linewidth=2, label='Factor de escala a(t)')
plt.plot(t_eval, y_sol[2], 'r-', linewidth=2, label='Campo Φ(t)')
plt.yscale('log')
plt.xlabel('Tiempo cósmico (t)')
plt.ylabel('Valor (escala log)')
plt.title('Evolución Cosmológica: a(t) y Φ(t)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(2, 1, 2)
plt.plot(y_sol[2], y_sol[3], 'g-', linewidth=1)
plt.xlabel('Φ(t)')
plt.ylabel('dΦ/dt')
plt.title('Diagrama de Fase del Campo Escalar')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('evolucion_cosmologica.png', dpi=300, bbox_inches='tight')
plt.show()

print("   • Gráficos guardados en 'evolucion_cosmologica.png'")

# ==========================================
# PARTE VI: EXPORTACIÓN A LATEX
# ==========================================
print("\n7. EXPORTANDO ECUACIONES A LATEX...")

with open('L_total.tex', 'w', encoding='utf-8') as f:
    f.write(sp.latex(L_total_final))

with open('eqs_movimiento.tex', 'w', encoding='utf-8') as f:
    f.write("\\text{Ecuación de Klein-Gordon modificada:}\n")
    f.write(sp.latex(eqs_mov[0]) + "\n\n")
    f.write("\\text{Ecuación de Friedmann modificada:}\n")
    f.write(sp.latex(eqs_mov[1]))

print("   • Lagrangiano exportado a 'L_total.tex'")
print("   • Ecuaciones de movimiento exportadas a 'eqs_movimiento.tex'")

# ==========================================
# RESUMEN FINAL
# ==========================================
print("\n" + "=" * 60)
print("RESUMEN DE RESULTADOS:")
print("=" * 60)
print(f"• Frecuencia de corte calibrada: ν_c = {nu_c_calib_val:.4e} Hz")
print(f"• Masa del ALP predicha: m_ALP = {masa_alp_eV:.2e} eV")
print(f"• Factor de escala final: a(t_final) = {y_sol[0, -1]:.4f}")
print(f"• Campo escalar final: Φ(t_final) = {y_sol[2, -1]:.4f}")
print("• Archivos generados:")
print("  - evolucion_cosmologica.png (gráficos)")
print("  - L_total.tex (Lagrangiano)")
print("  - eqs_movimiento.tex (ecuaciones)")
print("=" * 60)
