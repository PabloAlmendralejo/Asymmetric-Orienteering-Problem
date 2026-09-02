# Parche: modelo de fatiga asimétrico no lineal (Sec. 3.6)

Sustituye `f(t) = 1 + λt/B` (orden-invariante, ver discusión EJOR-D-26-01706)
por un estado de fatiga `G_i` propagado arco a arco, con contribución
distinta para subida (φ+) y bajada (φ-, ponderada por `ρ`), y con `ρ`
tratado como parámetro de sensibilidad (sweep), no calibrado a un único
valor.

## Archivos modificados

- `python/core/pathfinding.py` — el Dijkstra anisotrópico ahora también
  acumula, sobre el mismo árbol de caminos mínimos que ya usa para `cm`,
  el desnivel positivo (`gain`, φ+) y negativo (`loss`, φ-) hasta cada nodo.
  `compute_cost_matrix` devuelve ahora `(cm, gain_m, loss_m)` en vez de
  solo `cm`.
- `python/run_pipeline.py` — exporta `gain`, `loss` y `rho_default` en el
  JSON de cada instancia, junto a `cm`/`pts`/`bud_raw`/`fatigue_rate` como
  antes.
- `cpp/solver_flow_torremocha.cpp`, `cpp/solver_flow_la_muela.cpp` —
  formulación Flow extendida:
  - `psi_arc(i,j,rho) = phi_plus_ij - rho*phi_minus_ij` (constante por arco).
  - `compute_fatigue_bounds`: cota `Ĝ_i` vía Bellman-Ford en semianillo
    (max,+), acotado a n-1 rondas (válido pese a posibles ciclos de peso
    positivo en el grafo de ψ — ver comentario en el código).
  - `rcost_fatigue_asym`: coste real con clip (`max(0,·)`), usado para
    validar rutas concretas (SA, solución entera del B&C) — sustituye a
    `rcost_fatigue` (modelo antiguo) en `is_feasible_route` y en SA.
  - Nuevas variables de flujo `g_col` + restricciones
    `add_fatigue_flow_coupling` / `add_fatigue_flow_propagation` /
    `add_fatigue_budget_asym`, que sustituyen a `add_fatigue_budget_flow`
    (queda definida pero sin llamar, para comparación A/B si la quieres).
  - `main()` acepta `--rho-sweep`: resuelve cada mapa para
    ρ ∈ {0, 0.25, 0.5, 0.75, 1.0} y escribe `<nombre>_rho<v>.json`.
  - Cada JSON de salida ahora incluye `fatigue_cost` (modelo corregido) y
    `fatigue_cost_legacy` (modelo original del paper), lado a lado.

## Pendiente (no incluido en este parche)

1. **`solver_mtz_*.cpp`**: la formulación MTZ necesita el equivalente
   McCormick de `g_col`/`Ghat` (mismo patrón que ya usas para `t_i`/`w_ij`,
   sustituyendo `t_i → G_i`, `C_ij → psi_ij`). No lo he tocado porque el
   patrón McCormick es distinto al de flujo y merece su propia pasada.
2. **Proposiciones B2/B3 (routing infeasibility, cycle cover)**: falta
   quitar el factor de inflación `(1+λ/2)` de `MRC(Q)`/`MRC(F)` — con el
   bound conservador ya no hace falta, y usarlo ahora sería incorrecto
   (sobreestimaría la cota). No están implementadas como cortes de
   separación en este código (`find_and_add_cover_cuts` es solo B1); si
   los añades, deriva `MRC` sin ese factor.
3. **Recalibración de `fatigue_rate` (λ)**: al quitar la normalización
   `/B` del término de fatiga (ya no tiene sentido dividir un estado en
   metros de desnivel entre un presupuesto en unidades de coste), λ
   necesita recalibrarse — los valores usados en el paper (`λ=0.2`, etc.)
   ya no son directamente comparables.
4. Replicar en `python/config/*.py` si hay algo que asuma la firma antigua
   de `compute_cost_matrix` (no encontré otros call-sites, pero conviene
   revisar antes de correr el pipeline completo).

## Cómo probar rápido

```bash
cd python
python run_pipeline.py torremocha --preprocess   # regenera JSONs con gain/loss
cd ../cpp
# compilar solver_flow_torremocha.cpp como antes (misma toolchain HiGHS)
./solver_flow_torremocha            # una corrida, rho=0.5 (rho_default)
./solver_flow_torremocha --rho-sweep  # barrido de sensibilidad en rho
```
