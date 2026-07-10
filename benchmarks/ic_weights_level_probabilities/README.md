# Snellere en schaalbare level probabilities voor `ic_weights()`

**Onderzoeksvraag:** kan het bepalen van de level probabilities (chi-bar-kwadraat
gewichten) in `ic_weights()` (`R/ext_ic_infer.R`) ook op een andere manier, die
veel sneller is en schaalbaar?

**Antwoord: ja.** De exacte (pmvnorm) methode is *inherent* exponentieel in het
aantal constraints en valt niet fundamenteel te versnellen. Een Monte
Carlo-schatter die elke draw als een standaard **NNLS-probleem**
(Lawson–Hanson, gecompileerd Fortran, package `nnls`) in *constraint-ruimte*
oplost, geeft dezelfde gewichten tot op MC-nauwkeurigheid (~0.001–0.002 bij
R = 1e5), kost ~2 seconden ongeacht g ≤ 20, en schaalt polynomiaal door tot
g = 100+. Zie het prototype `con_weights_nnls()` in `bench_ic_weights.R`.

## 1. Waarom de huidige exacte methode niet schaalbaar is

`ic_weights()` implementeert de methode van Kudô (via ic.infer): voor elke
subsetgrootte k = 1, …, floor((g−2)/2) worden **alle** C(g, k) deelverzamelingen
van constraints doorlopen, met per subset 2–4 orthant-kansen (`pmvnorm`). Het
aantal pmvnorm-aanroepen groeit dus als ~2^(g−1). Gemeten (4 cores, sequentieel):

| g  | exact (pmvnorm) | NNLS-MC (R = 1e5) | max\|verschil\| |
|----|-----------------|--------------------|-----------------|
| 6  | 0.03 s          | 1.7 s              | 0.0011 |
| 8  | 0.22 s          | 1.7 s              | 0.0019 |
| 10 | 1.9 s           | 1.7 s              | 0.0016 |
| 12 | 15 s            | 1.9 s              | 0.0013 |
| 14 | 108 s           | 2.0 s              | 0.0022 |
| 16 | ~15 min (extrapolatie ×7–8 per +2) | ~2 s | — |
| 20 | uren            | 2.5 s              | — |
| 50 | onmogelijk      | 7 s                | — |
| 100| onmogelijk      | 24 s               | — |

De recent toegevoegde future-parallelisatie verzacht dit hooguit met een factor
(#cores); de exponentiële groei blijft. Ook slimmere exacte
orthant-algoritmes (Miwa–Hayter–Kuriki, TVPACK) veranderen niets aan de
C(g, k)-enumeratie zelf — exact rekenen blijft exponentieel.

## 2. Het alternatief: NNLS-simulatie in constraint-ruimte

De level probabilities zijn per definitie

> w_k(W) = P( de projectie van s ~ N(0, W) op de kegel {b ≥ 0} onder de
> metriek W⁻¹ heeft precies k positieve componenten ),  met W = A·V·Aᵀ (g × g).

Herschrijf met de Cholesky-decompositie W = UᵀU en M = U⁻ᵀ, z = M·s ~ N(0, I):

```
min_{b ≥ 0} (s − b)ᵀ W⁻¹ (s − b)   ⇔   min_{b ≥ 0} ‖ z − M·b ‖²
```

Elke MC-draw is dus een **standaard NNLS-probleem met vaste designmatrix M**
en een i.i.d. standaardnormale respons z — geen `rtmvnorm`/`rmvnorm`-sampling,
geen p-dimensionale QP. `nnls::nnls()` (Lawson–Hanson, Fortran) lost dit op in
~16 µs per draw bij g = 8. Tel het aantal positieve componenten, tabuleer:
klaar. Antithetische paren (z, −z) geven gratis variantiereductie.

Voordelen boven de bestaande `con_weights_boot()` (quadprog):

1. **Onafhankelijk van p.** De boot-methode werkt in θ-ruimte (p-dim QP +
   p-dim multivariate sampling); de NNLS-methode alleen in constraint-ruimte
   (g-dim). Gemeten bij g = 8, R = 2e4:

   | p   | boot (µs/draw) | NNLS (µs/draw) |
   |-----|----------------|-----------------|
   | 11  | 36             | 18 |
   | 25  | 53             | 16 |
   | 50  | 127            | 16 |
   | 100 | 497            | 17 |

2. **Nauwkeuriger bij gelijke rekentijd.** Bij g = 8 haalde
   `con_weights_boot` (convergentiecriterium 1e-3, gestopt na 35 000 draws,
   1.24 s) een max. fout van 0.0057 t.o.v. exact; NNLS deed 100 000 draws in
   1.65 s met max. fout 0.0010.

3. **Zelfde doelgrootheid, dus drop-in.** De schatter reproduceert
   `rev(ic_weights(W))` tot op MC-fout (gevalideerd voor g = 4–14). Voor
   meq > 0 geldt dezelfde reductie als nu in `con_weights()`:
   pas de methode toe op `solve(solve(W)[-(1:meq), -(1:meq)])`.

### Nauwkeurigheid voor de GORIC-penalty

De penalty gebruikt PT = 1 + Σ k·w_k. Bij g = 10 (PT_exact = 2.6927, 20
herhalingen):

| R      | sd(PT) | max\|fout\| |
|--------|--------|-------------|
| 1e4    | 0.0050 | 0.011 |
| 5e4    | 0.0035 | 0.008 |
| 1e5    | 0.0020 | 0.004 |

Een fout van ~0.005 op de penalty is verwaarloosbaar op GORIC-schaal
(vergelijkbaar met de MC-fout die `pmvnorm` (GenzBretz) zelf al introduceert).

## 3. Kan `con_weights_boot()` (quadprog) veilig vervangen worden?

Ja — mits de **duale vorm** wordt gebruikt (zie `check_edge_cases.R`). De
primale vorm hierboven vereist `chol(W)` en faalt als W = A·V·Aᵀ singulier is
(meer constraints dan parameters, of redundante constraints). De duale vorm
van hetzelfde projectieprobleem,

```
min_{λ ≥ 0} ‖ z̃ − D·λ ‖²,   D = Lᵀ·Aᵀ (p × g),  V = L·Lᵀ,  z̃ ~ N(0, I_p)
```

is óók een standaard NNLS-probleem, vereist alleen dat V positief-definiet is
(zoals quadprog nu ook al eist), en telt via het aantal positieve λ's exact
hetzelfde als `solve.QP`: dimL = p − #actieve ongelijkheden. Empirisch
gevalideerd:

- **meq = 0, W regulier**: duaal = primaal = exact (max. verschil ~0.001).
- **g > p (W singulier)**: `chol(W)` faalt; duaal en boot geven identieke
  verdelingen (max. verschil 0.0016 = MC-fout van twee onafhankelijke runs).
- **meq > 0**: reductie `solve(solve(W)[−(1:meq), −(1:meq)])` (zoals in
  `con_weights()`) matcht zowel exact als boot; de boot-massa ligt precies op
  levels (p−g):(p−meq) zoals verwacht.

Aandachtspunten bij een echte vervanging:

1. **Output-contract behouden**: boot levert p+1 gewichten (levels 0:p) met
   `attr(method) = "boot"`; `penalty_goric()` en `penalty_complement_goric()`
   indexeren per methode verschillend. De duale telling levert dit formaat
   direct.
2. **meq > 0 in de duale vorm**: equality-multipliers zijn vrij (niet ≥ 0);
   oplosbaar door de equality-kolommen uit D en z̃ weg te projecteren en NNLS
   op de rest te doen, of door de reductie-route te nemen als W regulier is.
3. **`...` naar `rtmvnorm`**: `con_weights_boot()` geeft `...` door aan
   `rtmvnorm` (truncatie via lower/upper). Wordt dat gebruikt, dan is de
   NNLS-schatter met standaardnormale draws niet equivalent — quadprog-pad
   als fallback behouden.
4. **Reproduceerbaarheid**: zelfde seed geeft andere (even geldige) getallen
   dan de oude boot; snapshot-tests met hardgecodeerde waarden moeten worden
   bijgewerkt.
5. `nnls` faalt in de praktijk niet (geen `error.idx`-mechanisme nodig), maar
   de teldrempel (`x > tol`) verdient een vaste, gedocumenteerde tolerantie.

> **Status:** de duale NNLS-motor is geïmplementeerd als standaard engine in
> `con_weights_boot()` (`R/compute_chiBarSquare_weights.R`), met het oude
> rtmvnorm + quadprog-pad als fallback wanneer truncatie-argumenten via `...`
> worden meegegeven. Zie `tests/testthat/test-con_weights_boot.R`.
>
> Naar aanleiding van review is bovendien:
>
> - de **per-draw equivalentie** met `solve.QP` aangetoond met gedeelde
>   trekkingen: 100% overeenstemming in vijf scenario's × 20.000 draws,
>   inclusief meq > 0, bijna-afhankelijke rijen, exact gedupliceerde rijen en
>   veel simultaan actieve restricties (`check_shared_draws.R`; let op de
>   tekenconventie z ↔ −Uᵀz̃ die daar wordt toegelicht);
> - de **stopregel vervangen**: de oude regel vergeleek opeenvolgende
>   cumulatieve gewichten (die convergeren vanzelf, ongeacht de werkelijke
>   MC-fout — empirisch: stop bij crit 1e-3 met werkelijke fout 0.0057).
>   De nieuwe regel stopt zodra de 95%-CI-halfbreedte (1.96 × MC-SE) van elk
>   gewicht onder `convergence_crit` ligt, met de SE geschat uit de
>   between-pair variantie van de antithetische paren (nnls-engine) dan wel
>   de binomiale variantie (quadprog-engine). De default is herijkt naar
>   5e-3 (vergelijkbare runtime, maar nu een echte precisiegarantie —
>   gemeten: claim 0.0050, werkelijke fout 0.0014). Het attribuut `mc_se`
>   rapporteert de SE per gewicht;
> - **antithetic sampling genuanceerd**: de variantiereductie is groot voor
>   de penalty (sd-ratio ~1.6 = variantiefactor ~2.6) maar per afzonderlijk
>   gewicht wisselend; daarom schat de stopregel de SE's empirisch in plaats
>   van onafhankelijkheid aan te nemen;
> - **inputvalidatie** toegevoegd (dimensies, symmetrie, PD via chol, meq/R/
>   chunk_size/convergence_crit), plus een rank-check op de equality-rijen en
>   een guard tegen `n_valid == 0`.

## 4. Aanbeveling

- Houd de exacte route voor **g ≤ 10 à 12** (daar is hij sneller dan MC en
  exact); de analytische shortcuts voor dim 1–3 blijven zoals ze zijn.
- Vervang daarboven — én als snellere motor achter `mix_weights = "boot"` —
  de quadprog-lus door de NNLS-schatter (`con_weights_nnls()` in dit
  benchmarkscript). Het bestaande chunk/convergentie-mechanisme van
  `con_weights_boot()` kan er ongewijzigd omheen. Nieuwe dependency: `nnls`
  (klein, geen verdere dependencies, al decennia op CRAN).
- De NNLS-lus is triviaal te parallelliseren (future), maar dat is pas nodig
  vanaf g ≳ 50 of R ≳ 1e6; single-core is 1e5 draws in ~2 s voor g ≤ 14.
- Verdere opties indien ooit nodig: batched Lawson–Hanson in C/C++ (nog eens
  factor 3–5, scheelt vooral R-loop-overhead) en exact berekende eindpunten
  w_0 en w_g (2 pmvnorm-calls) combineren met MC voor de rest.

## Reproduceren

```r
# vereist: mvtnorm, quadprog, tmvtnorm, nnls, future, future.apply
Rscript benchmarks/ic_weights_level_probabilities/bench_ic_weights.R
```

Gedraaid met R 4.3.3, mvtnorm 1.2-4, nnls, quadprog (Ubuntu 24.04, 4 cores).
