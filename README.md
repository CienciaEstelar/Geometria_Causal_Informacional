
# 🌌 Geometría Causal-Informacional (GCI)

### Framework de Cosmología Escalar-Tensor (DVT)

Este repositorio contiene el **Dynamic Vacuum Toolkit (DVT)**, enfocado en la implementación de modelos de **Energía Oscura Dinámica**.  
El proyecto principal es la **Geometría Causal-Informacional (GCI)**, que unifica principios de información, gravedad cuántica y cosmología.

---

## 🔎 Descripción del Proyecto

El script principal `modelo_cosmologico.py` es un framework simbólico-numérico de cosmología que realiza:

- **Derivación Simbólica:**  
  Genera las ecuaciones de **Klein-Gordon** y **Friedmann** modificadas a partir del Lagrangiano (usando `sympy`).

- **Calibración:**  
  Fija la frecuencia de corte fundamental `ν_c` y la densidad de energía del vacío `ρ_vac` con alta precisión.

- **Simulación Numérica:**  
  Integra las ecuaciones diferenciales ordinarias (EDOs) para la evolución cosmológica `a(t)` y `Φ(t)`.

- **Predicciones Clave:**  
  Calcula la **masa del axión predicha** (`m_a ≈ 3.61 meV`).

---

## 📂 Estructura del Repositorio

| Carpeta / Archivo | Contenido Principal | Propósito |
|--------------------|--------------------|------------|
| `/` | `modelo_cosmologico.py` | Ejecución central del modelo |
| `analysis/` | Scripts de verificación de consistencia (e.g. `derivacion_vc.py`) | Blindaje teórico — Punto Fijo RG |
| `docs/` | Manuscritos, apéndices y documentación técnica | Referencia teórica |
| `output/` | Ecuaciones (LaTeX), gráficos (`.png`, `.pdf`), y resultados | Resultados de simulación |

---

## ▶️ Instrucciones de Reproducción

### 1. Instalación de Dependencias

Se requiere **Python 3.9+** y las librerías científicas listadas en `requirements.txt`:

```bash
pip install -r requirements.txt
````

---

### 2. Ejecución de la Simulación Principal

Este script realiza la derivación simbólica, calibración y simulación numérica:

```bash
python modelo_cosmologico.py
```

---

### 3. Verificación del Blindaje Teórico

Ejecute este script para verificar simbólicamente la consistencia de la escala `ν_c`
con la simetría del Grupo de Renormalización (**Punto Fijo RG**):

```bash
python analysis/derivacion_vc.py
```

---

## 📈 Resultados

Una vez completada la ejecución, los archivos de salida (`.tex`, `.png`, `.pdf`)
se generarán automáticamente dentro de la carpeta `output/`.

---

## 🧩 Créditos

**Autor:** Juan Galaz
**Proyecto:** *Dynamic Vacuum Toolkit (DVT) — Geometría Causal-Informacional (GCI)*
**Universidad:** Universidad de Santiago de Chile (USACH)

---

## 📜 Licencia

Este proyecto se distribuye bajo licencia **MIT**, promoviendo la libre investigación científica.
