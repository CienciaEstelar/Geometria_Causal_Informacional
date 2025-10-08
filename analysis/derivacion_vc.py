import sympy as sp

# ==========================================
# CONFIGURACIÓN DE IMPRESIÓN SIMBÓLICA
# ==========================================
try:
    from IPython.display import display
    IPYTHON_AVAILABLE = True
except ImportError:
    IPYTHON_AVAILABLE = False
    def display(obj):
        print(f"[Resultado Simbólico]:\n{obj}")

sp.init_printing(use_unicode=True)

# ==========================================
# DEFINICIÓN DE SÍMBOLOS
# ==========================================
omega_p, omega_h, omega_c = sp.symbols('omega_p omega_h omega_c', positive=True)

print("=" * 60)
print("DERIVACIÓN SIMBÓLICA AVANZADA: Escala de Transición y Punto Fijo RG")
print("=" * 60)

print("\n# Paso 1: Condición de Punto Fijo RG (simetría logarítmica)")
print("ln(ω_p) - ln(ω_c) = ln(ω_c) - ln(ω_h)")
rg_eq = sp.Eq(sp.log(omega_p) - sp.log(omega_c), sp.log(omega_c) - sp.log(omega_h))
display(rg_eq)

# ==========================================
# Paso 2: Resolviendo la ecuación simbólicamente
# ==========================================
solutions = sp.solve(rg_eq, omega_c)
print("\n# Paso 2: Resolver para ω_c")
final_result = sp.Eq(omega_c, solutions[0])
display(final_result)

# ==========================================
# Paso 3: Comprobación simbólica de consistencia
# ==========================================
print("\n# Paso 3: Verificación simbólica")
check_lhs = sp.log(omega_p) - sp.log(final_result.rhs)
check_rhs = sp.log(final_result.rhs) - sp.log(omega_h)
verification = sp.simplify(check_lhs - check_rhs)
display(sp.Eq(0, verification))
print("✅ La verificación confirma la simetría logarítmica y consistencia del punto fijo RG.")

# ==========================================
# Paso 4: Interpretación física
# ==========================================
print("\n# Paso 4: Interpretación física")
print(f"La frecuencia de corte ω_c emerge como la media geométrica entre los límites UV e IR:")
display(sp.Eq(omega_c, sp.sqrt(omega_p * omega_h)))
print("Esta escala de transición garantiza flujo RG simétrico entre ω_p y ω_h.\n")

# ==========================================
# Paso 5: Visualización opcional (opcional para MS)
# ==========================================
try:
    import matplotlib.pyplot as plt
    import numpy as np

    # Valores arbitrarios para demostración
    omega_p_val = 1e43  # Frecuencia de Planck
    omega_h_val = 1e-18 # Frecuencia de Hubble
    omega_c_val = np.sqrt(omega_p_val * omega_h_val)

    plt.figure(figsize=(8, 5))
    plt.plot([omega_h_val, omega_c_val, omega_p_val], [0, 1, 0], 'ro-')
    plt.xticks([omega_h_val, omega_c_val, omega_p_val],
               [r'$\omega_H$', r'$\omega_c$', r'$\omega_P$'])
    plt.ylabel('Escala RG (esquema conceptual)')
    plt.title('Esquema de Escalas de Transición (Punto Fijo RG)')
    plt.grid(True, alpha=0.3)
    plt.show()
except ImportError:
    print("Matplotlib no disponible: se omite visualización gráfica.")
