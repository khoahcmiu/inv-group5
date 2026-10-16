# Inventory Lab — Stochastic Inventory Optimization (s, Q)

An educational Python application for solving the **continuous-review (s, Q) inventory model** with **stochastic lead-time demand** and **constant replenishment lead time**.

Inventory Lab implements the EOQ-initialized iterative procedure taught in the *Individual Item with Probabilistic Demand* lecture and illustrated in Examples 1 and 2 of `a.pdf`. It shows the numerical result first, then lets the reader expand the formulas, exact input transformations, and step-by-step calculations.

> **Scope:** This is a teaching implementation of the lecture's continuous approximation with backorders. It is not an exact optimizer for every inventory setting or a simulation of real inventory transactions.

## Contents

- [Features](#features)
- [Run locally](#run-locally)
- [Deploy with Streamlit Community Cloud](#deploy-with-streamlit-community-cloud)
- [Mathematical model](#mathematical-model)
- [Lead-time demand distributions](#lead-time-demand-distributions)
- [Worked examples](#worked-examples)
- [Understanding the results](#understanding-the-results)
- [Accuracy and assumptions](#accuracy-and-assumptions)
- [Tests and project structure](#tests-and-project-structure)
- [Course material](#course-material)

## Features

- **Seven continuous lead-time demand distributions:** Normal, Uniform, Exponential, Gamma, Lognormal, Weibull, and Triangular.
- **Flexible math inputs:** enter integers, decimals, fractions (`14/52`), radicals (`sqrt(18)` or `√18`), and powers (`3^2` or `3²`).
- **Four tabs:**
  - **Results & formulas:** optimized values with expandable derivations and numerical substitutions.
  - **Iterations:** worked iterations **0–10**, plots, and CSV/JSON export. The solver can continue internally until convergence.
  - **Costs & service:** cost components, service indicators, and explanatory notes.
  - **Theory:** named formulas, assumptions, and distribution-specific expressions.
- **Exact input algebra when possible:** SymPy preserves fractions and radicals in the explanation layer.
- **Numerical optimization:** SciPy and Python floating-point arithmetic evaluate inverse distribution functions and expected shortages without rounding solver inputs for display.
- **Readable results:** summary metrics show three decimal places; underlying calculations retain their numerical precision.

## Run locally

### Requirements

- Python **3.11 or later**
- Windows, macOS, or Linux

From the **repository root** (the directory containing `app.py` and `pyproject.toml`):

```bash
python -m pip install -e ".[ui]"
python -m streamlit run app.py
```

Open the local URL printed by Streamlit (usually <http://localhost:8501>).

**First steps:** Choose `Example 1 — a.pdf` or `Example 2 — a.pdf` in the sidebar. The normal-distribution examples load automatically. For modified inputs or a custom dataset, select the relevant distribution, enter your parameters, and click **Calculate**.

Optional developer dependencies and tests:

```bash
python -m pip install -e ".[ui,dev]"
python -m pytest -q
```

You can also use the command-line solver without Streamlit:

```bash
python -m inventory_sq --case ex1 --csv ex1_steps.csv --json ex1_results.json
python -m inventory_sq --case ex2 --csv ex2_steps.csv
```

## Deploy with Streamlit Community Cloud

1. Upload the project to a GitHub repository, keeping `app.py`, `pyproject.toml`, and `src/inventory_sq/` in their current relative locations.
2. Place a file named `requirements.txt` in the repository root with the following contents:

   ```text
   -e .[ui]
   ```

   This installs the local `inventory_sq` package and the Streamlit UI dependencies declared in `pyproject.toml`.
3. Visit [share.streamlit.io](https://share.streamlit.io/), sign in, and create an app from your repository.
4. Choose the repository and branch; use `app.py` as the **main file path**.
5. Deploy and share the generated application URL.

> Do not upload virtual environments, build directories, caches, or private credentials. Streamlit Community Cloud must install the package as well as its third-party dependencies.

## Mathematical model

### Notation

| Symbol | Meaning |
|:--|:--|
| $D$ | Expected demand per accounting period |
| $L$ | Constant lead time, measured in the same time unit |
| $X$ | Random total demand during lead time |
| $\theta=E[X]$ | Expected lead-time demand |
| $\sigma_X$ | Standard deviation of lead-time demand |
| $A$ | Fixed cost per replenishment order |
| $v$ | Purchase cost per unit |
| $i$ | Holding-cost rate per accounting period |
| $h=iv$ | Holding cost per unit per accounting period |
| $p$ | Shortage penalty per unit short |
| $Q$ | Replenishment order quantity |
| $s$ | Reorder point |
| $G(s)$ | Cumulative distribution function of $X$ |
| $n(s)$ | Expected units short per replenishment cycle |

**Consistent time units are essential.** For instance, if $D$ is annual demand, $h$ must be annual holding cost per unit and $L$ must be expressed in years.

### 1. Economic Order Quantity (EOQ) — initial order size

```math
Q_0=\sqrt{\frac{2DA}{h}}
```

### 2. Expected lead-time shortage (loss function)

```math
n(s)=E[(X-s)^+]=\int_s^{\infty}(x-s)g(x)\,dx
```

Here, $g(x)$ is the probability density function of lead-time demand, and $(x-s)^+=\max(x-s,0)$.

### 3. Expected cost per accounting period

```math
Y(Q,s)=\frac{AD}{Q}+h\left(\frac{Q}{2}+s-\theta\right)
+\frac{pD}{Q}n(s)+Dv
```

The four terms represent **ordering**, **approximate holding**, **shortage**, and **purchasing** costs, respectively. Because $Dv$ is constant when unit cost is fixed, it does not affect the optimizing values of $Q$ and $s$.

### 4. First-order condition for the reorder point

```math
\frac{\partial Y}{\partial s}=0
\quad\Longrightarrow\quad
P(X>s)=1-G(s)=\frac{Qh}{pD}
```

### 5. First-order condition for order quantity

```math
\frac{\partial Y}{\partial Q}=0
\quad\Longrightarrow\quad
Q=\sqrt{\frac{2D[A+p\,n(s)]}{h}}
```

### 6. Iterative solution procedure

Start with EOQ. For iteration $i=0,1,2,\ldots$:

```math
\begin{aligned}
Q_0 &= \sqrt{\frac{2DA}{h}},\\[4pt]
\alpha_i &= \frac{Q_i h}{pD},\\[4pt]
s_i &= G^{-1}(1-\alpha_i),\\[4pt]
n_i &= E[(X-s_i)^+],\\[4pt]
Q_{i+1} &= \sqrt{\frac{2D[A+p\,n_i]}{h}}.
\end{aligned}
```

The implementation uses the **inverse survival function** (ISF) for $s_i$, avoiding unnecessary loss of precision from subtracting a small tail probability from 1.

A finite interior solution requires $0<\alpha_i<1$. The solver reports an error if this condition is violated, rather than silently clipping the probability.

## Lead-time demand distributions

The lecture directly illustrates **Normal**, **Uniform**, and **Exponential** cases. The other four distributions are extensions of the same continuous first-order method.

| Distribution | Input parameters | Mean $E[X]$ | Notes |
|:--|:--|:--|:--|
| Normal | Base-period mean $\mu$, SD $\sigma$, conversion $m$, lead time $L$ | $D L$, where $D=m\mu$ | Symmetric, unbounded |
| Uniform | Lower $a$, upper $b$ | $(a+b)/2$ | Bounded and flat |
| Exponential | Scale $\beta$ | $\beta$ | Right-skewed, nonnegative |
| Gamma | Shape $k$, scale $\beta$ | $k\beta$ | Flexible right-skewed model |
| Lognormal | Arithmetic mean $\theta$, SD $\sigma_X$ | $\theta$ | Positive, potentially heavy right tail |
| Weibull | Shape $k$, scale $\beta$ | $\beta\Gamma(1+1/k)$ | Nonnegative, flexible shape |
| Triangular | Lower $a$, mode $c$, upper $b$ | $(a+c+b)/3$ | Bounded with a most-likely value |

All seven models compute the distribution-specific **inverse survival function** and **expected shortage**, rather than substituting a Normal approximation.

For Normal demand, if $k=(s-\theta)/\sigma_X$, the normal loss formula is:

```math
n(s)=\sigma_X\left[\phi(k)-k\bigl(1-\Phi(k)\bigr)\right]
```

For Uniform demand, when $a\le s\le b$:

```math
n(s)=\frac{(b-s)^2}{2(b-a)}
```

For Exponential demand with scale $\beta$, when $s\ge0$:

```math
n(s)=\beta e^{-s/\beta}
```

The application’s **Theory** tab also documents the Gamma, Lognormal, Weibull, and Triangular formulas, including domain-specific cases.

> **Important:** A coherent lead-time demand model ordinarily requires $E[X]=DL$. If the chosen distribution parameters do not satisfy this identity, the app displays a warning; it does not silently adjust user inputs.

## Worked examples

The examples below reproduce the handwritten exercises and retain continuous-valued $Q$ and $s$.

### Example 1 — monthly Normal demand

Monthly demand follows $N(28,8^2)$. With 12 months per year and fixed lead time $L=14/52$ years:

```math
\begin{aligned}
D &= 28\times12=336,\\
\operatorname{Var}(D) &= 12\times8^2=768,\\
\theta &= DL=336\times\frac{14}{52}=\frac{1176}{13},\\
\operatorname{Var}(X) &=768\times\frac{14}{52}=\frac{2688}{13}.
\end{aligned}
```

Given $A=15$, $v=6$, $i=0.30$, and $p=7$:

```math
h=iv=1.8,\qquad Q_0=\sqrt{\frac{2(336)(15)}{1.8}}=20\sqrt{14}
```

| Result | Converged value | Displayed in the app |
|:--|--:|--:|
| Optimal $Q^*$ | 81.32283593 | **81.323** |
| Reorder point $s^*$ | 112.55219958 | **112.552** |
| Safety stock $s^*-\theta$ | 22.09066112 | **22.091** |
| Inventory control cost / year | 186.14429470 | **186.144** |

### Example 2 — Normal demand per cost period

Demand follows $N(100,3^2)$ per period. With $L=2$ periods:

```math
\theta=200,\qquad \operatorname{Var}(X)=18,\qquad\sigma_X=\sqrt{18}
```

Given $A=200$, $v=40$, $i=0.10$, and $p=30$:

```math
h=4,\qquad Q_0=\sqrt{\frac{2(100)(200)}{4}}=100
```

| Result | Converged value | Displayed in the app |
|:--|--:|--:|
| Optimal $Q^*$ | 102.17103581 | **102.171** |
| Reorder point $s^*$ | 204.65597563 | **204.656** |
| Safety stock $s^*-\theta$ | 4.65597563 | **4.656** |

**Source-note correction:** In the handwritten solution, the calculation gives approximately **204.66** for the reorder point; one later copied line says **205.66**. The application follows the calculation, not the transcription error.

## Understanding the results

- **Reorder point:** When inventory position reaches $s$, place an order of size $Q$.
- **Safety stock:** $SS=s-\theta$; this measures inventory above expected demand during the lead time.
- **Ordering cost:** $C_{\mathrm{order}}=AD/Q$.
- **Approximate holding cost:** $C_{\mathrm{hold}}=h(Q/2+s-\theta)$.
- **Shortage cost:** $C_{\mathrm{short}}=(pD/Q)n(s)$.
- **Stockout probability:** $P(X>s)$, the chance that lead-time demand exceeds the reorder point.
- **Approximate fill rate:** $1-n(s)/Q$, the fraction of demand fulfilled immediately on average under this approximation.
- **Stationarity residuals:** Numerical values of $\partial Y/\partial Q$ and $\partial Y/\partial s$ near zero indicate proximity to a stationary point, **not** a guarantee of global optimality.

## Accuracy and assumptions

1. **Continuous approximation:** The method treats order quantity and reorder point as continuous. It does not optimize an integer policy automatically.
2. **Fixed lead time:** The current solver assumes a constant replenishment lead time, continuous review, and backordered shortages.
3. **Time consistency:** Demand rate, lead time, and holding-cost rate must use compatible time units.
4. **Standard deviation vs variance:** The Normal input labeled `σ` accepts **standard deviation**, not variance. The app derives the corresponding lead-time variance.
5. **Distribution validity:** The model does not support discrete demand distributions such as Poisson, Binomial, or Negative Binomial without a separate discrete optimization method.
6. **Numerical accuracy:** Input fractions and radicals can remain symbolic for display; CDF quantiles, loss-function evaluations, and final optima are computed numerically. Three-decimal output **does not** mean intermediate values were rounded.
7. **Convergence:** The current UI computes at least iterations 0–10 (11 updates) and can continue until its stopping tolerance is met. A small stationarity residual is a diagnostic, not a global-optimality certificate.
8. **Interpretation:** Normal demand can assign some probability to negative demand. The lecture's approximate holding expression may also be negative for some input combinations; interpret those cases cautiously.

## Tests and project structure

Run the regression tests:

```bash
python -m pytest -q
```

The tests cover handwritten-example values, distribution loss functions, input handling, and iteration behavior. Browser layout may require a separate manual check after Streamlit upgrades.

```text
inventory-lab/
├── app.py                         # Streamlit user interface
├── pyproject.toml                 # Package metadata and dependencies
├── requirements.txt               # Add for Streamlit Cloud deployment
├── README.md                      # This document
├── src/
│   └── inventory_sq/
│       ├── __init__.py
│       ├── __main__.py
│       ├── cli.py
│       ├── examples.py
│       ├── explanations.py
│       ├── models.py              # Seven continuous distributions
│       └── solver.py              # (s, Q) iterative solution
├── examples/                      # Reproducible outputs and scripts
└── tests/                         # Automated tests
```

## Course material

This project was developed using:

- `a.pdf`, handwritten Examples 1 and 2 (especially pages 1–3).
- `chap2_Singleitem_ProbabilisticDemand.pptx`, continuous-review (s, Q) lectures, especially the loss-function derivation and iterative solution procedure.

Those course files are **not bundled with the application**. They are referenced to explain the teaching context. Four distributions—Gamma, Lognormal, Weibull, and Triangular—are explicitly identified as **extensions** beyond the three distributions illustrated in the lecture.
