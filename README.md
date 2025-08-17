# InfoGravity: Gravidade como Compressão Informacional

> **TL;DR**: Este repositório acopla um termo de **gravidade informacional** $\nabla\Phi_I$ a qualquer solver N-body e traz um **demo 2D** de Poisson modificada $\nabla\cdot(\mu\,\nabla\phi)=4\pi G\rho$. A proposta trata a gravidade como um efeito **entrópico-informacional** (compressão de dados do universo), recuperando Newton no limite clássico e sugerindo correções em baixa aceleração.

---

## Índice

* [Visão Geral](#visão-geral)
* [Fundamentos — a teoria em 10 minutos](#fundamentos--a-teoria-em-10-minutos)
* [Estrutura do Projeto](#estrutura-do-projeto)
* [Instalação Rápida](#instalação-rápida)
* [Uso Rápido (Copy/Paste)](#uso-rápido-copypaste)

  * [C++](#c)
  * [Python](#python)
  * [Delphi (VCL/FM X)](#delphi-vclfm-x)
  * [Demo 2D de Poisson Modificada](#demo-2d-de-poisson-modificada)
* [Parâmetros & Tuning](#parâmetros--tuning)
* [API de Referência](#api-de-referência)

  * [API C++](#api-c)
  * [API Python](#api-python)
  * [API Delphi](#api-delphi)
* [Validação, Testes e Boas Práticas](#validação-testes-e-boas-práticas)
* [Performance & Escalabilidade](#performance--escalabilidade)
* [Roadmap](#roadmap)
* [FAQ](#faq)
* [Referências & Leituras](#referências--leituras)
* [Licença](#licença)

---

## Visão Geral

**InfoGravity** é um kit de integração + demo numérica que implementa o termo de **gravidade emergente por compressão informacional**:

$\Phi_I(\mathbf{x}) = -\alpha\,\ln\!\Big(\frac{\rho_I(\mathbf{x})}{\rho_{I0}}\Big), \quad \; \mathbf{a}_I = -\nabla\Phi_I.$

* **$\rho_I$**: um *proxy* de densidade/ordem informacional (pode vir de barions, temperatura, dispersão de velocidades etc.).
* **$\alpha$**: rigidez "compressiva" que controla a intensidade do acoplamento.
* **$\rho_{I0}$**: escala de referência (normalização) para tornar o log adimensional.

Aceleração total em um N-body fica:
$\mathbf{a}_{\text{total}} = \mathbf{a}_{\text{Newton}} + \mathbf{a}_I.$

Incluímos também um **solver 2D** para a equação de Poisson com coeficiente variável:
$\nabla\cdot\big(\mu(\rho_I)\,\nabla\phi\big) = 4\pi G\,\rho,$
com um $\mu$ simples e suave: $\mu(\rho_I)=1+\eta\,\frac{\rho_I}{\rho_I+\rho_0}$. Isso permite comparar **newtoniano** vs **modificado** em mapas e curvas de rotação sintéticas.

> **Status**: protótipo científico/engenharia — pronto para experimentar, calibrar e falsificar em dados.

---

## Fundamentos — a teoria em 10 minutos

### 1) Axiomas mínimos

* **Pixelização/Holografia**: número de graus de liberdade $\propto$ área (lei de área).
* **Landauer**: apagar/organizar $\Delta I$ nats custa $Q\ge k_B T\,\Delta I$.
* **Segunda Lei da Infodinâmica**: a entropia informacional efetiva $S_I$ tende a **não aumentar** (processo espontâneo busca representações mais compactas).
* **Força Entrópica**: $\mathbf{F}_I = T_{\rm eff}\,\nabla S_I$.
* **Equipartição em telas**: $E=\tfrac{1}{2}N k_B T$.

Dessas peças, recupera-se **Newton** e obtêm-se correções em baixa aceleração (escala $a_0\sim cH_0$).

### 2) Potencial "ZIP" informacional

$\Phi_I = -\alpha\ln(\rho_I/\rho_{I0})$ implica $\nabla\Phi_I = -\alpha\,\nabla\ln\rho_I$. Gradientes de ordem informacional geram aceleração extra. Em regiões de baixa aceleração, o termo **pode simular massa adicional efetiva** (curvas de rotação planas).

### 3) Poisson Modificada

$\mu(\rho_I)>0$ deforma o operador e permite ajustar suavemente a transição entre regime clássico e regime dominado por informação.

> **Cuidado**: isto é **física especulativa com pés na termodinâmica/holografia**. O objetivo é produzir **predições falsificáveis** e comparar com dados (galáxias, aglomerados, lentes fracas etc.).

---

## Estrutura do Projeto

```
infogravity/
├─ cpp/
│  ├─ include/
│  │  └─ info_gravity.hpp        # Header-only: Φ_I + integrador Velocity-Verlet
│  └─ examples/
│     └─ nbody_minimal.cpp       # Exemplo de uso
├─ delphi/
│  ├─ uInfoGravity.pas           # Classe TInfoGravity + tipos TVec3/Provider
│  └─ Example_Integrate.pas      # Esqueleto de acoplamento no seu loop
├─ python/
│  ├─ info_gravity.py            # Classe InfoGravity (plug-in aceleração)
│  └─ demo_poisson2d.ipynb       # Solver 2D: ∇·(μ∇φ)=4πGρ + gráficos e CSV
├─ data/
│  └─ outputs/                   # CSVs e figuras geradas pelo demo 2D
├─ tests/
│  ├─ cpp/
│  ├─ python/
│  └─ delphi/
├─ README.md                     # este arquivo
└─ LICENSE
```

---

## Instalação Rápida

### Requisitos

* **C++**: C++17+, CMake (opcional), qualquer toolchain moderna (GCC/Clang/MSVC).
* **Delphi**: Delphi 12 (Athens) ou superior. Funciona em VCL e FMX.
* **Python**: 3.10+; pacotes: `numpy`, `matplotlib` (para o demo 2D).

### Setup Python rápido

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt  # (numpy, matplotlib)
```

> Se não houver `requirements.txt`, instale manualmente: `pip install numpy matplotlib`.

### Build C++ (opcional com CMake)

```bash
mkdir -p build && cd build
cmake ..
cmake --build . --config Release
```

> O header é *single-include*: também dá para só incluir `cpp/include/info_gravity.hpp` no seu projeto sem CMake.

---

## Uso Rápido (Copy/Paste)

### C++

**Inclua** `info_gravity.hpp` e forneça um *provider* de `rho_I` + `grad(rho_I)`:

```cpp
#include "info_gravity.hpp"

InfoParams P; P.alpha = 0.8; P.rho0 = 1.0; P.epsGrad = 1e-12;
auto rhoProv = [](const Vec3& x, double& rhoI, Vec3& gradRhoI){
  const double sigma2 = 1.0;
  const double r2 = x.x*x.x + x.y*x.y + x.z*x.z;
  rhoI = std::exp(-0.5 * r2 / sigma2);
  gradRhoI = Vec3{ -x.x, -x.y, -x.z } * (rhoI / sigma2);
};
InfoGravity info(P, rhoProv);

// No integrador: a_total = a_newton + a_info
Vec3 a_info = info.acceleration(p[i].x);
```

Se quiser, use o **`NBodyStepper`** já pronto para Velocity-Verlet (veja `cpp/examples/nbody_minimal.cpp`).

### Python

```python
from python.info_gravity import InfoGravity
import math
ig = InfoGravity(alpha=0.8, rho0=1.0)

def rho_gauss(x):
    r2 = x[0]**2 + x[1]**2 + x[2]**2
    rho = math.exp(-0.5*r2)
    grad = (-x[0]*rho, -x[1]*rho, -x[2]*rho)
    return rho, grad

ig.set_rho_provider(rho_gauss)
aI = ig.acceleration((1.0, 0.0, 0.0))
print(aI)
```

### Delphi (VCL/FM X)

```pascal
uses uInfoGravity;

var Info: TInfoGravity; Par: TInfoParams;

procedure RhoGauss(const P: TVec3; out RhoI: Double; out GradRhoI: TVec3);
var R2: Double; begin
  R2 := P.X*P.X + P.Y*P.Y + P.Z*P.Z;
  RhoI := Exp(-0.5*R2);
  GradRhoI := TVec3.Make(-P.X*RhoI, -P.Y*RhoI, -P.Z*RhoI);
end;

begin
  Par.Alpha := 0.8; Par.Rho0 := 1.0; Par.EpsGrad := 1e-12;
  Info := TInfoGravity.Create(Par, RhoGauss);
  // a_I = Info.AccelInfo(Posicao)
end;
```

### Demo 2D de Poisson Modificada

Abra `python/demo_poisson2d.ipynb` e rode as células. Ele:

1. Gera um disco exponencial $\rho_b$ (ou outro perfil) e define $\rho_I$ (default: $\rho_I=\rho_b$).
2. Resolve $\nabla\cdot(\mu\,\nabla\phi)=4\pi G\rho$ em malha regular com `Gauss–Seidel red-black` e *harmonic averaging* de $\mu$ nas faces.
3. Exporta CSV de perfis radiais (\nabla\phi, v\_c(r)) e figuras para `data/outputs/`.

Parâmetros-chave no topo do notebook: `eta`, `rho0`, `max_iters`, `tol`, `dx`, `G`.

---

## Parâmetros & Tuning

* **`alpha`** (Φ\_I): intensidade do termo informacional. Aumenta |a\_I|.
* **`rho0`** (Φ\_I): normalização da densidade informacional para o log. Evita singularidades.
* **`epsGrad`** (Φ\_I): piso numérico para dividir por ρ\_I com segurança.
* **`mu(ρ_I)`** (Poisson 2D): `mu = 1 + eta * (ρ_I / (ρ_I + ρ0))`.

  * `eta`: quanto a modificação pesa (0 ⇒ Newton puro).
  * `rho0`: escala de transição de `mu`.
* **Discretização**: `dx` menor ⇒ mais caro e preciso. Use *multigrid* se escalar malha.

Sugestões práticas:

* Comece com `alpha ∈ [0.3, 1.2]`, `eta ∈ [0.1, 2.0]`.
* Fixe `rho0` próximo à mediana de ρ\_I para boa condição numérica do log e de `mu`.

---

## API de Referência

### API C++

**Tipos**

* `struct Vec3 { double x,y,z; ... }`
* `struct Particle { Vec3 x, v; double m; }`
* `struct InfoParams { double alpha, rho0, epsGrad; }`
* `using RhoIProvider = std::function<void(const Vec3&, double&, Vec3&)>;`

**Principais funções/classes**

* `phiI(alpha, rhoI, rho0) -> double`: retorna Φ\_I.
* `gradPhiI(InfoParams, rhoI, gradRhoI) -> Vec3`: retorna ∇Φ\_I.
* `class InfoGravity { acceleration(x: Vec3)->Vec3; }`.
* `class NBodyStepper { step(vector<Particle>&); }`: Velocity-Verlet minimalista (soma Newton + Info).

**Padrão de uso**

1. Crie `InfoParams`.
2. Crie `RhoIProvider` que compute `rhoI` e `gradRhoI`.
3. `InfoGravity info(P, provider);` e chame `info.acceleration(x)` no seu integrador.

### API Python

* `class InfoGravity(alpha=1.0, rho0=1.0, eps_grad=1e-12)`

  * `set_rho_provider(fn: (Vec3)->(rhoI, grad))`
  * `acceleration(x: Vec3) -> Vec3`
  * `grad_phiI(rhoI, grad_rhoI) -> Vec3`

### API Delphi

**Tipos**

* `TVec3 = record (X,Y,Z: Double)`
* `TRhoIProvider = reference to procedure(const P: TVec3; out RhoI: Double; out GradRhoI: TVec3);`
* `TInfoParams = record (Alpha, Rho0, EpsGrad: Double)`

**Classe**

* `TInfoGravity.Create(Params, RhoProvider)`
* `function GradPhiI(RhoI, GradRhoI): TVec3`
* `function AccelInfo(P: TVec3): TVec3`  // retorna a\_I

**Integração**

* Some `AccelInfo` à aceleração newtoniana no seu loop (Velocity-Verlet, RK4, Leapfrog etc.).

---

## Validação, Testes e Boas Práticas

* **Teste de consistência**: com `alpha=0` e `eta=0`, os resultados devem reproduzir Newton (curvas de rotação e potenciais idênticos ao baseline).
* **Convergência**: refine `dx` e `max_iters` e verifique convergência de `φ` e `|g|`.
* **Energia**: monitore energia total em N-body (tolerância de drift).
* **Unidades**: seja rigoroso. Defina um sistema de unidades consistente (SI ou unidades naturais) e mantenha em todo o pipeline.

---

## Performance & Escalabilidade

* **N-body**: use Barnes–Hut ($\mathcal{O}(N\log N)$) ou FMM se Newton ficar caro. O termo informacional é local (via `rho_provider`) e barato; `rho_I` em grid pode ser amostrado com cache/buckets.
* **Poisson 2D**: o Gauss–Seidel serve para grids pequenos/médios. Para grandes domínios, migre para **multigrid** (V-cycle) ou **precondicionado** (CG/BiCGStab) com operador $\mu$ variável e *harmonic averaging* nas faces.

---

## Roadmap

* [ ] **Notebook** refinado com *multigrid*.
* [ ] **Curvas de rotação reais** (SPARC/THINGS) com ajuste global de `alpha`, `eta`, `rho0`.
* [ ] **Lenteamento fraco simulado** (ray-shooting) sob $\mu(\rho_I)$.
* [ ] **Aglomerados**: comparação com mapas de raio-X/SZ.
* [ ] **Port** GPU (CUDA/OpenCL) para N-body + amostragem de `rho_I` em textura 3D.

---

## FAQ

**Isso substitui a Relatividade Geral?**
Não. A meta é um **modelo efetivo** que recupera Newton e proponha correções testáveis. A extensão totalmente relativística fica para versões futuras.

**Como escolho `ρ_I`?**
Comece com proxies simples ($\rho_I=\rho_b$ suavizada). Em aplicações sérias, inclua estado térmico, dispersões e anisotropias. O importante é consistência e validação empírica.

**Isso é “simulação” de universo-computador?**
Não precisa supor um computador externo. A tese é que **o universo se organiza informacionalmente**, e a gravidade emerge disso.

---

## Referências & Leituras

* Landauer, R. *Irreversibility and Heat Generation in the Computing Process*.
* Bekenstein, J. D.; Hawking, S. W. (entropia de buracos negros, lei de área).
* Unruh, W. (temperatura de Unruh).
* Verlinde, E. (gravidade emergente).

> Esta base inspira; nossa implementação concreta (Φ\_I, μ(ρ\_I)) é um **modelo experimental de trabalho**.

---

## Licença

**MIT** — uso livre com atribuição. Veja `LICENSE`.

---

## Como contribuir

* Issues e PRs são bem-vindos.
* Ao abrir PR, inclua: (i) descrição clara, (ii) resultados antes/depois, (iii) impacto em performance.

---

### Anexos — Dicas práticas

* **Gradient clipping**: para evitar explosões quando ρ\_I→0, use `epsGrad`.
* **Smoothing**: aplique kernel gaussiano em ρ\_I para reduzir ruído numérico.
* **Unidades naturais**: em protótipos, defina `G=1`, `ρ0=1` e escale depois.
* **Exportar CSV**: o notebook salva `radius, g_newton, g_eff, v_newton, v_eff` para plot externos.


