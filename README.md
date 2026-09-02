# Modelo de fatiga asimétrico no lineal (Sec. 3.6)

Sustituye `f(t) = 1 + λt/B` (orden-invariante, ver discusión EJOR-D-26-01706)
por un estado de fatiga `G_i` propagado arco a arco:

```
ψ_ij = φ+_ij − ρ·φ-_ij + μ·δ_ij
```

con contribución distinta para subida (`φ+`), bajada (`φ-`, ponderada por
`ρ`) y distancia horizontal (`δ`, ponderada por `μ`). `ρ` se trata como
parámetro de sensibilidad (sweep), no calibrado a un único valor. `μ` se
deriva de la curva de Minetti ya usada para el coste base (no es un
parámetro libre). `λ` tampoco tiene un valor único de literatura
utilizable en las nuevas unidades — se trata como sweep también, con
`1.75e-5` como ancla de orden de magnitud (ver discusión completa en
`paper/orienteering_paper.tex`, Sec. 3.6 y Sec. 7).

El estado real (clip) es `F_j = max(0, F_i + ψ_ij)`; la relajación LP usa
la versión sin clip `G_i` (cota superior válida sobre `F_i`) por las
razones detalladas en el paper.

## Archivos modificados

- `python/core/cost_functions.py` — `derive_mu()`: deriva `μ` de la
  curva de Minetti (coste de 1 metro llano ÷ coste del metro vertical
  más barato alcanzable).
- `python/core/pathfinding.py` — el Dijkstra anisotrópico acumula, sobre
  el mismo árbol de caminos mínimos que usa para `cm`, el desnivel
  positivo (`gain`, φ+), negativo (`loss`, φ-) y la distancia real
  recorrida (`dist`, δ) hasta cada nodo. `compute_cost_matrix` devuelve
  ahora `(cm, gain_m, loss_m, dist_m)`. También corrige un bug de
  unidades preexistente: `gain`/`loss` usaban `step_len` (recuento de
  pasos de grid) en vez de metros reales (`step_len*ds*cell_m`) —
  entendía mal el desnivel real de terreno por un factor `ds*cell_m`
  (16x para Torremocha).
- `python/run_pipeline.py` — exporta `gain`, `loss`, `dist`,
  `rho_default` y `mu_default` en el JSON de cada instancia.
- `benchmark/generate_instances.py` — sintetiza `dist` (distancia
  euclídea del arco) y `mu_default` (copia local de `derive_mu`) para
  las 21 instancias sintéticas de benchmark.
- `cpp/solver_flow_{torremocha,la_muela}.cpp`,
  `cpp/solver_mtz_{torremocha,la_muela}.cpp`,
  `benchmark/benchmark_solver_{flow,mtz,ablation}.cpp` — las 7 variantes
  del solver (ambas formulaciones, ambos terrenos reales, más las
  variantes de benchmark) implementan el modelo completo:
  - `psi_arc(i,j,rho) = gain_ij − rho·loss_ij + mu·dist_ij`.
  - `compute_fatigue_bounds`/`compute_fatigue_lower_bounds`: cotas
    `Ĝ_i`/`Ǧ_i` vía Bellman-Ford en semianillo (max,+)/(min,+), acotado
    a n-1 rondas.
  - `rcost_fatigue_asym`: coste real con clip (`max(0,·)`), usado para
    validar rutas concretas (SA, solución entera del B&C).
  - MTZ: variable de estado `G_col`/`Ghat` + envolvente McCormick
    `u_col` (equivalente a `t_i`/`w_ij` del modelo antiguo).
  - Flow: variables de flujo `g_col` con acoplamiento **de dos lados**
    (`Glow[i]` y `Ghat[i]`, no solo `Ghat[i]` — necesario porque `ψ`
    puede ser negativo, a diferencia del flujo de tiempo original).
  - `main()` acepta `--rho-sweep` (ρ ∈ {0, 0.25, 0.5, 0.75, 1.0}) y
    `--lambda-sweep` (λ ∈ {0, 1e-5, 1.75e-5, 5e-5, 1e-4}) — no ambos a
    la vez. `benchmark_solver_ablation.cpp` no tiene ningún sweep (su
    CLI es para alternar familias de cortes, no para sensibilidad).
  - Cada JSON de salida incluye `fatigue_cost` (modelo corregido) y
    `fatigue_cost_legacy` (modelo original del paper), lado a lado.

## Compilar

No hay toolchain de C++ preinstalado; requiere MSVC (Build Tools) +
HiGHS compilado desde fuente (headers/lib, no solo `highs.dll`):

```bash
# una vez: clonar y compilar HiGHS (CMake + MSVC, Release)
git clone https://github.com/ERGO-Code/HiGHS.git
cmake -B HiGHS/build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DFAST_BUILD=ON HiGHS
cmake --build HiGHS/build --config Release

# compilar cada solver (ejemplo con vcvarsall de MSVC ya cargado)
cl.exe /O2 /MD /std:c++20 /EHsc cpp/solver_flow_torremocha.cpp ^
  /I HiGHS/highs /I HiGHS/build ^
  /link /LIBPATH:HiGHS/build/Release/bin highs.lib /out:solver_flow_torremocha.exe
```

## Cómo probar rápido

```bash
cd python
python run_pipeline.py torremocha --preprocess   # regenera JSONs con gain/loss/dist
cd ../cpp
./solver_flow_torremocha              # una corrida, rho=0.5/mu derivado (defaults del JSON)
./solver_flow_torremocha --rho-sweep     # barrido de sensibilidad en rho
./solver_flow_torremocha --lambda-sweep  # barrido de sensibilidad en lambda
```

Las instancias sintéticas de benchmark (`benchmark/instances/*.json`) ya
incluyen `gain`/`loss`/`dist`/`mu_default` y sirven para probar sin
regenerar datos reales de terreno.

## Pendiente

1. **Datos reales de terreno**: los `op_input_torremocha_*.json` /
   `op_input_la_muela_*.json` no existen en este checkout — hay que
   volver a correr el pipeline de Python (arriba) para regenerarlos con
   `gain`/`loss`/`dist` correctos (incluyendo el fix del bug de unidades).
   Los solvers `cpp/solver_*` nunca se han corrido de punta a punta
   contra terreno real, solo contra las instancias sintéticas.
2. **Recalibración de λ**: sigue sin haber un valor de literatura
   directamente utilizable en las nuevas unidades (metros de desnivel,
   no fracción de presupuesto). Se trata como sweep con `1.75e-5` como
   ancla de orden de magnitud, no como valor calibrado — ver
   `paper/orienteering_paper.tex` Sec. 7 para el detalle de la búsqueda
   de literatura y el razonamiento.
3. **Resultados numéricos del paper** (Sec. 6-7: tablas Torremocha/La
   Muela, comparativa MTZ-vs-Flow, estudio de ablación): siguen
   reflejando el modelo antiguo. Actualizarlos honestamente requiere
   correr la batería experimental completa (horas de cómputo) con datos
   reales regenerados — no hecho todavía.
4. **Validez general de `Ĝ_i`/`Ǧ_i` como cota**: usadas en la práctica y
   verificadas en instancias de prueba, pero no hay una prueba
   rigurosa general de que `Ĝ_i` sea siempre una cota válida sobre el
   estado real para cualquier grafo — no afecta la corrección de las
   rutas finalmente aceptadas (siempre se re-validan exactas antes de
   reportarse), pero sí podría afectar la completitud/optimalidad de la
   búsqueda B&C en casos adversariales no explorados.
