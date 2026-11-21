# SIMTELE-RISK

Simulador educativo de riesgo para un sistema de transporte por cable inspirado en el teleférico de La Paz, Bolivia. Es el proyecto integrador de **Métodos Numéricos**: cuatro casos de estudio, cada uno evaluado en tres condiciones (**NORMAL**, **ESTRÉS** y **MAL CONDICIONADO**), con teoría de errores aplicada en todos los cálculos y un componente de inteligencia artificial contrastado siempre con la solución exacta.

> **Aviso.** Es un proyecto estrictamente académico. Los datos son sintéticos, los parámetros son supuestos y los escenarios son hipotéticos: no describen fallas verificadas de ninguna empresa ni tienen conexión con equipos de control reales. Ningún resultado debe usarse como diagnóstico, norma de seguridad o plan de mantenimiento.

---

## 1. Contenido de la aplicación

La aplicación es **un único archivo** (`index.html`), sin dependencias externas, sin servidor y sin paso de compilación. Todo el cálculo se ejecuta en el navegador.

| Pantalla | Contenido |
|---|---|
| **Inicio** | Cuatro tarjetas de módulos, resumen de teoría de errores y acceso a sistemas lineales. |
| **1 Detención** | Detención de emergencia de cabinas (EDO, integración, raíces, sustituto IA). |
| **2 Deterioro** | Tiempo hasta un umbral de inspección (raíces, calibración con sistemas lineales). |
| **3 Oscilaciones** | Oscilaciones por viento (EDO, integración, mal condicionamiento, IA exploratoria). |
| **4 Sensores** | Reconstrucción de señales con datos faltantes y ruido (interpolación, filtro IA). |
| **Sistemas lineales** | Jacobi, Gauss-Seidel, LU y gradiente conjugado con referencia directa. |

Cada módulo muestra el modelo, los parámetros editables con validación, el selector de escenario, tablas de resultados e iteraciones, gráficos y dos conclusiones **separadas**: una sobre los métodos numéricos y otra sobre el problema aplicado.

---

## 2. Instalación y ejecución

No hay nada que instalar.

**Opción A: abrir el archivo.** Descarga `index.html` y ábrelo con doble clic en un navegador moderno (probado en Chromium). Funciona sin conexión.

**Opción B: servidor local (opcional).**

```bash
python3 -m http.server 8000
# abrir http://localhost:8000
```

**Opción C: GitHub Pages.**

1. Sube `index.html` a la raíz del repositorio (rama `main`).
2. En *Settings → Pages*, elige *Deploy from a branch*, rama `main`, carpeta `/ (root)`.
3. Espera a que se publique y abre la URL que GitHub te indica.

### Estructura del repositorio

```
.
├── index.html                       # aplicación completa (HTML + CSS + JavaScript)
├── README.md                        # este documento
├── supabase_schema.sql              # esquema PostgreSQL para Supabase
├── supabase_permiso_eliminar.sql    # corrección: permiso de eliminar pruebas
├── GUIA_SUPABASE.md                 # crear y cargar la base de datos
└── GUIA_CORRECCIONES_SUPABASE.md    # correcciones y solución de problemas
```

---

## 3. Los cuatro módulos

### Módulo 1: Detención de emergencia de cabinas

**Pregunta:** ¿cómo varían el tiempo y la distancia de detención con la velocidad, el frenado y el paso numérico?

Modelo: `dx/dt = v`, `dv/dt = −a(v)`, con x en m, v en m/s y t en s.

| Escenario | Modelo de frenado | v₀ | a₀ | k | h |
|---|---|---|---|---|---|
| NORMAL | `a = a₀` constante | 5 | 1.5 | 0 | 0.05 |
| ESTRÉS | `a = a₀ + k·v` (fricción viscosa, sobrecarga) | 10 | 1.5 | 0.3 | 0.5 |
| MAL CONDICIONADO | `a = a₀√v` (raíz doble en t_det) | 5 | 0.5 | no aplica | 0.2 |

- **Referencia analítica:** constante, `t_det = v₀/a`, `d = v₀²/(2a)`; viscosa, `v = (v₀+a/k)e^(−kt) − a/k`; raíz cuadrada, `v = (√v₀ − a₀t/2)²`.
- **Métodos:** RK2 (variante de Heun) y RK4; bisección (intervalo [0, 3·t_det]), Newton-Raphson y secante sobre `v(t) = 0`; trapecio y Simpson 1/3 sobre la tabla de RK4; Romberg con refinamiento `n = 2^q`.
- **Salidas:** tabla temporal frente a la referencia, errores por método, raíces, refinamiento de integración, gráfico de v(t), gráfico de convergencia log-log, sensibilidad (+1 % en v₀ y en a₀) y predicción del sustituto IA.
- **Mal condicionado:** `v(t)` toca cero sin cambiar de signo, por lo que la bisección no aplica; Newton converge de forma lineal y un residuo pequeño no garantiza una raíz precisa. Es un subexperimento numérico, no una falla física.

### Módulo 2: Deterioro mecánico y mantenimiento

**Pregunta:** ¿cuándo alcanzaría un indicador hipotético de deterioro un umbral de inspección?

Modelo: `D(t) = D₀ + αt + βt²`, se resuelve `g(t) = D(t) − D_umbral = 0` (D adimensional, t en meses).

| Escenario | D₀ | α | β | D_umbral |
|---|---|---|---|---|
| NORMAL | 0.10 | 0.006 | 0.00008 | 0.65 |
| ESTRÉS | 0.12 | 0.016 | 0.00020 | 0.65 |
| MAL CONDICIONADO | 0.649 | 0.00001 | 0.00000001 | 0.65 |

- **Referencia analítica:** raíz de la cuadrática, con la raíz físicamente admisible. En NORMAL, t* ≈ 53.5 meses.
- **Métodos:** bisección, Newton-Raphson y secante (tolerancia 10⁻¹⁰, máximo 100 iteraciones). Calibración auxiliar de coeficientes con un sistema SPD 4×4 resuelto con Jacobi, Gauss-Seidel, LU y gradiente conjugado.
- **Salidas:** resultados por método, comparación de los tres escenarios, tabla de iteraciones (selector de método), sensibilidad a D₀ (+0.1 %), calibración, dos gráficos y exportación a CSV.

### Módulo 3: Oscilaciones por viento

**Pregunta:** ¿cómo cambian las oscilaciones laterales de una cabina hipotética con la fuerza externa, el amortiguamiento y la rigidez?

Modelo: `m x″ + c x′ + k x = F₀ sen(ωt)`, un grado de libertad (m en kg, c en N·s/m, k en N/m, F₀ en N, ω en rad/s, x en m).

| Escenario | m | c | k | F₀ | ω | x₀ | T | h | Umbral |
|---|---|---|---|---|---|---|---|---|---|
| NORMAL | 1000 | 400 | 4000 | 200 | 1 | 0.10 | 30 | 0.1 | 0.30 |
| ESTRÉS | 1000 | 150 | 4000 | 800 | 1.9 | 0.10 | 30 | 0.1 | 0.30 |
| MAL CONDICIONADO | 1000 | 400 | 4000 | 200 | 1 | 0.10 | 30 | 0.1 | 0.30 + calibración 4×4 |

- **Referencia analítica:** solución exacta de la EDO (subamortiguada, crítica y sobreamortiguada).
- **Métodos:** Heun (RK2) y RK4; integral `J = ∫x(t)² dt` con trapecio, Simpson 1/3 y Romberg (tabla triangular); estudio de refinamiento del paso con orden observado; calibración SPD 4×4 con número de condición κ₂ y amplificación de una perturbación δb.
- **Pestañas:** *Simulación* (métricas, trayectoria, error, tabla temporal), *Validación numérica* (convergencia, cuadratura, sensibilidad ±1 %, mal condicionamiento) e *IA y ciencia de datos* (regresión ridge sobre 100 simulaciones sintéticas, 80 % entrenamiento y 20 % prueba, semilla 2026).
- **Límites:** hasta 5000 pasos por corrida; sin aerodinámica 3D, dinámica del cable, ráfagas estocásticas ni control. El umbral de desplazamiento es didáctico.

### Módulo 4: Reconstrucción de señales de sensores

**Pregunta:** ¿cómo reconstruir registros faltantes y medir la distorsión que produce el ruido?

Datos sintéticos de `y(t) = sin t` en [0, 2π] o `y(t) = e^(−t)` en [0, 5]. Se ocultan nodos interiores como puntos de prueba y se agrega ruido `ε ~ N(0, σ²)` a los nodos conocidos (semilla reproducible).

| Escenario | Señal | N | % ocultos | σ | Nodos |
|---|---|---|---|---|---|
| NORMAL | sin t | 14 | 14 | 0 | equiespaciados |
| ESTRÉS | sin t | 20 | 45 | 0.05 | irregulares |
| MAL CONDICIONADO | sin t | 20 | 15 | 0.01 | equiespaciados (fenómeno de Runge) |

- **Métodos:** Lagrange, Newton por diferencias divididas y spline cúbico natural. El sistema tridiagonal del spline se resuelve con LU y se verifica con Gauss-Seidel (diagonal dominancia estricta).
- **Salidas:** reconstrucción en cada punto ocultado, RMSE y error máximo (en los ocultos y en 401 puntos del intervalo), diferencia Lagrange-Newton, verificación LU contra Gauss-Seidel y dos gráficos (curvas superpuestas y error semilogarítmico).
- **Rangos válidos:** N entre 6 y 40, ocultos entre 0 y 70 %, σ ≥ 0.

### Sistemas lineales (componente transversal)

Se construye `b = A x*` con `x*` conocido para medir el error real. Todas las matrices son simétricas definidas positivas, aptas para gradiente conjugado.

| Escenario | Matriz |
|---|---|
| NORMAL | 4×4 estrictamente diagonal dominante |
| ESTRÉS | 12×12 tridiagonal (2.05 en la diagonal, −1 fuera), dominancia débil |
| MAL CONDICIONADO | Hilbert 7×7, κ∞ ≈ 10⁹ |

Se reportan κ∞(A), iteraciones, convergencia, residuo `‖r‖∞`, error relativo, tiempo y sensibilidad ante una perturbación de 10⁻⁶ en `b`, más un gráfico de residuo por iteración. La referencia directa es LU con pivoteo parcial.

---

## Historial de pruebas (base de datos en Supabase)

Cada módulo (1 a 4) tiene una subpestaña **Historial de pruebas** y una barra de guardado:

- **Guardar como nueva prueba:** guarda manualmente la simulación actual (parámetros, un registro por método, sus iteraciones y el resumen con las conclusiones).
- **Cargar:** devuelve una prueba al panel del módulo, recalcula y muestra una **verificación** (valores guardados frente a un nuevo cálculo) y las iteraciones guardadas.
- **Reemplazar prueba cargada:** guarda la versión modificada como nueva y elimina la anterior.
- **Eliminar:** borra la prueba y, en cascada, sus métodos, iteraciones y resumen.

Configuración: ejecutar `supabase_schema.sql` en el SQL Editor de Supabase (o, si ya estaba creada la base, `supabase_permiso_eliminar.sql`) y seguir `GUIA_SUPABASE.md` y `GUIA_CORRECCIONES_SUPABASE.md`. Las constantes `SBU` y `SBK` (URL y clave pública) están al inicio del bloque «Historial de pruebas» en `index.html`.

## 4. Cobertura de la consigna

| Familia | Métodos | Dónde |
|---|---|---|
| Errores (transversal) | absoluto, relativo %, residual, sensibilidad | todos los módulos |
| Sistemas lineales | Jacobi, Gauss-Seidel, LU, gradiente conjugado, referencia directa | Sistemas lineales, M2, M3, M4 (spline) |
| Interpolación | Lagrange, Newton, spline cúbico | M4 |
| Raíces | bisección, Newton-Raphson, secante | M1, M2 |
| Integración | trapecio, Simpson 1/3, Romberg | M1, M3 |
| EDO | RK2 (Heun), RK4 | M1, M3 |

**Definiciones de error usadas**

- Absoluto: `|e_a| = |v_real − v_aprox|`
- Relativo porcentual: `e_r % = |e_a / v_real| · 100` (se evita dividir por cero)
- Residual: `|R| = |f(x_calc)|`
- Sensibilidad: cambio relativo del resultado al perturbar una entrada.

La aplicación distingue **alta exigencia**, **riesgo físico**, **mal condicionamiento** y **falta de convergencia**, que no son lo mismo.

---

## 5. Inteligencia artificial

Paradigma híbrido: el **motor numérico es la fuente de verdad**; la IA solo anticipa o pre-procesa, y se contrasta siempre con la referencia.

- **Módulo 1, sustituto (*surrogate*):** regresión lineal sobre logaritmos, `ln d = w₀ + w₁ ln v₀ + w₂ ln a₀`, ajustada con 80 simulaciones sintéticas con 2 % de ruido (ecuaciones normales resueltas con LU). Con `a` constante se esperan w₁ ≈ 2 y w₂ ≈ −1; con `a = a₀√v`, w₁ ≈ 1.5. Solo vale dentro del rango entrenado (v₀ de 2 a 12 m/s, a₀ de 0.5 a 3).
- **Módulo 3, exploratoria:** regresión ridge que predice la amplitud máxima a partir de entradas estandarizadas. R² puede ser bajo o negativo; no garantiza nada fuera del dominio sintético.
- **Módulo 4, filtro de atípicos:** mediana móvil de 3 puntos con criterio MAD (`1.4826·MAD`, umbral de 3). Es un detector estadístico simple, **no un autoencoder**; no inventa datos en los nodos ocultos y puede aplanar picos reales.

---

## 6. Guía rápida de uso

1. Entra a un módulo desde las tarjetas de **Inicio** o desde la barra superior.
2. Elige el escenario: los parámetros se rellenan solos y la simulación se ejecuta.
3. Edita los parámetros y pulsa **Ejecutar**. Si un valor no es válido, se muestra un mensaje en rojo.
4. Revisa tablas, gráficos y las conclusiones numérica y del problema.
5. En los módulos 2 y 3 puedes **exportar un CSV** de los resultados.

### Verificación rápida de resultados conocidos

| Caso | Valor esperado |
|---|---|
| M1 NORMAL | t_det = 3.333 s, distancia = 8.333 m (RK2, RK4 y analítica coinciden) |
| M1 ESTRÉS | t_det ≈ 3.662 s, distancia analítica ≈ 15.023 m (RK4 se acerca más que RK2) |
| M1 MAL CONDICIONADO | t_det ≈ 8.944 s, distancia analítica ≈ 14.907 m; la bisección no aplica |
| M2 NORMAL | t* ≈ 53.5 meses |
| Sistemas lineales | κ∞ ≈ 3.14 (NORMAL), ≈ 44.8 (ESTRÉS), ≈ 9.85·10⁸ (MAL CONDICIONADO) |

---

## 7. Limitaciones y alcance de esta versión

- **El historial requiere internet y la base de datos configurada:** los cálculos funcionan sin conexión, pero guardar, listar, cargar y eliminar pruebas usa Supabase (ver la sección «Historial de pruebas»).
- RK2 es la variante de **Heun**. En el módulo 1, Romberg se calcula sobre RK4 con `n = 2^q` intervalos en `[0, t_det]` (no sobre la tabla con paso h).
- El cruce de umbral por raíces sobre la curva interpolada (módulo 4) queda para una fase posterior.
- Los parámetros de los escenarios del módulo 1, del módulo 4 y de sistemas lineales son supuestos elegidos para la demostración.
- Sin pruebas automatizadas: la verificación es manual, comparando con las referencias analíticas de la sección 6.

---

## 8. Créditos y licencia

- **Asignatura:** Métodos Numéricos, proyecto integrador SIMTELE-RISK.
- **Integrantes y mención:** _completar_
- **Docente:** _completar_
- **Licencia:** _completar según las reglas del repositorio_
