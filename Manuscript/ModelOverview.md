---
title: "Mathematical model overview"
author: "Petr Čapek"
date: \today
link-citations: true
bibliography: "zoteroCapek.bib"
---
# 1. sub-Microbial model
## 1. 1. Overview of basal model structure
The model recognizes two different carbon ($C$) pools of soil microbial biomass - basal microbial biomass ($B$) and growth intermediates ($G$). $G$ is "biomass specific" so it is expressed per unit of $B$ (e.g. $\mu mol~(G)~\mu mol~(B)^{-1}$). Total amount of $G$ per unit of oven-dried soil is thus, $B \times G$. Organic substrate ($S$) is taken up by microbial biomass at the rate defined by an equation:

$U=I_{m} \times B \times f$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S1) 

in which $U$ stands for organic C uptake rate, $I_{m}$ is maximum uptake rate constant, and $f$ is a hyperbolic scaling function $\frac{S}{K_{m} + S}$. $K_{m}$ is then the half-saturation constant. 
Organic C is assimilated into $G$ with fixed efficiency defined by the constant $y_{A}$ (i.e. assimilation efficiency). There are overhead cost of assimilation. The assimilation rate ($A$) is therefore expressed as:

$A=U \times y_{A}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S2).

The overhead costs are covered directly from substrate so the fraction of $U$ (defined as $1-y_{A}$) is respired. 

Once assimilated, $G$ is mobilized and used to produce basal biomass $B$. The mobilization rate ($J$) is defined by equation:

$J=\frac{I_{m} \times y_{A}}{G_{m}} \times G = v \times G$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S3),

in which $G_{m}$ represents the maximum content of $G$. To simplify notation of eq. 3, parameter $v$ derived from combination of $I_{m}$, $y_{A}$, and $G_{m}$ is defined. Energy is required to convert $G$ into $B$. This energy cost is covered from $G$ so analogically to substrate assimilation, only a fraction of $J$ is used for growth. We will denote this fraction of flux $J$ as a specific growth rate ($\mu$). If the energy cost of $B$ formation as well as energy production and yield is constant, $\mu$ can be defined by an equation:

$\mu =\frac{J}{1+g}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S4),

in which $g$ represents cost of biomass production in C-mol units. If respiration is the main process of energy generation, the remaining fraction of $J$ is respired as $CO_{2}$ at the rate $\mu \times g$. Commonly used terms growth yield or carbon use efficiency ($CUE$) are related to $g$. Specifically, $CUE$ equals $\frac{1}{1+g}$.

Pool $B$ decays over time as it ages. This decay rate ($H$) is assumed to be a first-order process defined by equation:

$H =h \times B$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S5),

in which $h$ denotes the "biomass specific" death rate constant. For sake of simplicity, we assume that $B$ requires no maintenance energy costs and that dead microbial biomass forms a stable pool of soil organic C, which cannot be rendered available to active microorganisms in course of three months of incubation (see Material and methods section in the main text for details about the incubation experiment).

The basal sub-Microbial model can be summarized by four differential equations:

$\frac{dS}{dt} =-U = - I_{m} \times B \times f; ~~~ f= \frac{S}{K_{m} + S}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S6),

$\frac{dG}{dt} = A \times\frac{1}{B} - J - \frac{dB}{dt} \times \frac{G}{B} = I_{m} \times f \times y_{A} - v \times G - \frac{dB}{dt} \times \frac{G}{B}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S7),

$\frac{dB}{dt} = \mu \times B - H = B \times (\frac{J}{1 + g} - h)$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S8),

$\frac{dCO_{2}}{dt} = U \times (1 - y_{A}) + B \times \mu \times g$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S9).

The last term of eq. S7 ($\frac{dB}{dt} \times \frac{G}{B}$) defines change of "biomass specific" $G$ due to changing of $B$, which $G$ is normalized to.

## 1. 2. Conversion factors (changing microbial biomass composition)

Total biomass of soil microorganisms in units of organic C per gram of oven-dried soil ($MBC$) is the sum of pools $B$ and $G$, thus:

$MBC = B + B \times G = B \times (1 + G)$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S10).

$MBC$ cannot be measured in soil directly. Instead, different proxy-parameters such as DNA, ATP, chloroform labile organic C (CLC) or phospholipid fatty acids (PLFA) are quantified and converted to $MBC$ using conversion factors. As we have shown in our previous study [@capek_revisiting_2023], the conversion factors are not constant because any proxy-parameter can represent different fractions of $B$ and $G$ (i.e. have different abundance in pool $B$ and $G$). Therefore, any conversion factor ($k_{*}$, symbol * in subscript denoting the proxy-parameter) between any proxy-parameter and $MBC$ can be defined as the weighted mean of an abundance of that proxy-parameter in pool $B$ ($n_{B-*}$) and $G$ ($n_{G-*}$):

$k_{*} = \frac{(n_{B-*} \times B) + (n_{G-*} \times B \times G)}{B + B \times G} = \frac{n_{B-*} + n_{G-*} \times G}{1 + G}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S11).

In this study, we require two conversion factors for two proxy-parameters - CLC ($k_{CLC}$) and PLFA ($k_{PLFA}$). These conversion factors are defined as:

$k_{CLC} = \frac{n_{B-CLC} + n_{G-CLC} \times G}{1 + G}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S12),

$k_{PLFA} = \frac{n_{B-PLFA}}{1 + G}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S13).

Here we use notation $k_{CLC}$ instead of common $k_{ec}$ for consistency. 

Note that eqs. S12 and S13 are not similar because $n_{G-PLFA}$ is assumed to be zero (PLFA is part of $B$ only). In respect to PLFA content of biomass, $k_{PLFA}$ decreases when $G$ increases and vice versa. Thus, the relative abundance of PLFA in $MBC$ is highest in "non-growing" biomass ($ \mu = 0$ when $G=0$). In respect to CLC, the opposite is true because $n_{G-CLC}$ is higher than $n_{B-CLC}$ (approx. two times [@capek_revisiting_2023]). The change of microbial biomass composition with growth rate in respect to any biomass component as defined here indirectly using transient pool $G$ is commonly observed in pure cultures and represents the main mechanism explaining proportion between energy generation by fermentative and respiratory pathway as described below.

## 1. 3. Connection between energy and mass fluxes

To explain why energy is generated by fermentative instead of respiratory pathway at fully aerobic conditions, the dependence of growth rate on energy production needs to be explicitly defined. In next two sections, we first define this dependence for non-fermentable substrate utilization (i.e. substrate that can only be processed by respiratory pathway) and then, we extend the general definition for fermentable substrate. For purpose of our study, we specify fermentable substrate as glucose ($C_{6}H_{12}O_{6}$ hereafter abbreviated as $glu$) and non-fermentable substrate as acetate or acetic acid ($CH_{3}COO^{-}$ or $CH_{3}COOH$ hereafter abbreviated as $acet$). We do not strictly distinguish between dissociated and non-disociated form of acetic acid because both forms exist in soil concurrently and both can be taken up by soil microorganisms. 

### 1. 3. 1. Non-fermentable substrate (acetate)

Energy can be generated from non-fermentable substrate only via respiratory pathway. Acetate is assimilated into Acetyl-CoA so it directly enters the TCA cycle. In contrast to glucose as fermentable substrate, it "skips" the glycolysis so it cannot generate fermentation products (note that we assume that all glucose passes the glycolysis before entering the TCA cycle). The acetate uptake ($U_{acet}$) and assimilation into pool $G$ ($A_{acet}$) follows eqs. S1 and S2, whose parameters $K_{m}$ and $y_A$ are acetate specific ($K_{m, acet}$ and $y_{A, acet}$):

$U_{acet}=I_{m} \times B \times f_{acet}; ~~~~ f_{acet} = \frac{S_{acet}}{K_{m, acet} + S_{acet}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S14),

$A_{acet}=U_{acet} \times y_{A, acet}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S15).

Assimilated acetate is used for energy and biomass production. It has been shown that energy and mass fluxes inside the microbial cells are tightly coupled e.g. [@stouthamer_theoretical_1973]. This coupling is defined by an equation:

$J_{E}=\sigma \times \mu$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S16),

in which $J_{E}$ is the "biomass specific" energy flux (in $\mu mols$ of ATP per unit of $B$ and day) and $\sigma$ denotes energy requirements of biomass (i.e. pool $B$) formation. Flux $J$ (eq. S3) with subscript $acet$ must be divided into $J_{E}$ and $\mu$ according to eq. S16. Therefore, eq. S3 can be rewritten to eq. S17:

$J_{acet}=v_{acet} \times G = \frac{J_{E}}{e_{acet}} + \mu$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S17),

in which $e_{acet}$ defines the molar amount of ATP formed by catalysis of one C-mol mobilized $G$ originating from assimilated acetate. The $v$ with subscript $acet$ reflects acetate specific $y_{A, acet}$. Combining eqs. S16 and S17, specific growth rate can be defined as:

$\mu = \frac{J_{acet}}{1 + \frac{\sigma}{e_{acet}}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S18).

The cost of biomass production defined as $g$ in eq. S4 equals $\frac{\sigma}{e_{acet}}$ in eq. S18. $CUE$ can then be defined as:

$CUE = \frac{1}{1 + \frac{\sigma}{e_{acet}}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S19).

Eq. S19 shows that CUE is 0.5 when $\sigma$ equals $e_{acet}$ (i.e. their ratio equals one), decreases when $\sigma$ (energy cost of biomass formation) decreases or $e_{acet}$ (efficiency of energy gain from assimilated substrate) increases and vice versa. $CO_{2}$ is then produced at rate:

$\frac{dCO_2}{dt} = U_{acet} \times (1 - y_{A, acet}) + B\times \mu \times \frac{\sigma}{e_{acet}} = U_{acet} \times (1 - y_{A, acet}) + B\times \frac{J_{acet}}{1 + \frac{e_{acet}}{\sigma}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S20).

### 1. 3. 2. Fermentable substrate (glucose)

The theory defined below for glucose metabolism is based on study of [@basan_overflow_2015]. Glucose is assimilated into 2 molecules of glyceraldehyde 3-phosphate. Once it passes the glycolysis, it can be either excreted as acetate or enter the TCA cycle. For sake of simplicity, we assume that glucose uptake and assimilation into G follows again eqs. S1 and S2 but with glucose specific $K_{m, gluc}$ and $y_{A, gluc}$:

$U_{gluc}=I_{m} \times B \times f_{gluc}; ~~~~ f_{gluc} = \frac{S_{gluc}}{K_{m, gluc} + S_{gluc}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S21),

$A_{gluc}=U_{gluc} \times y_{A, gluc}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S22).

Because glucose is fermentable substrate, energy can be produced by either fermentative or respiratory pathway. Therefore flux $J_{E}$ has theoretically two components:

$J_{E}=J_{E,r} + J_{E,f}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S23).

$J_{E,r}$ and $J_{E,f}$ denote energy fluxes generated by respiratory and fermentative pathway, respectively. If one of the components is zero, $\mu$ and $CUE$ can be defined analogically to metabolism of non-fermentable substrate (eqs. S18 and S19):

$\mu = \begin{cases} \frac{J_{gluc}}{1 + \frac{\sigma}{e_{r}}};~~~if~~~J_{E,f} = 0
                \\\\\ \frac{J_{gluc}}{1 + \frac{\sigma}{e_{f}}};~~~if~~~J_{E,r} = 0 \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S24).

$CUE = \begin{cases} \frac{1}{1 + \frac{\sigma}{e_{r}}};~~~if~~~J_{E,f} = 0
                \\\\\ \frac{1}{1 + \frac{\sigma}{e_{f}}};~~~if~~~J_{E,r} = 0 \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S25).

In eqs. S24 and S25, $e_{r}$ and $e_{f}$ denote molar amounts of ATP produced per one C-mol of assimilated glucose processed in respiratory and fermentative pathway, respectively. If $J_{E,f}$ is zero, the product of respiratory pathway (i.e. $CO_{2}$) is released at the rate $\frac{J_{E,r}}{e_{r}}$. If $J_{E,r}$ is zero, $CO_{2}$ and acetic acid are produced at rates $\frac{1}{3} \times \frac{J_{E,f}}{e_{f}}$ and $\frac{2}{3} \times \frac{J_{E,f}}{e_{f}}$, respectively (i.e. $C_{6}H_{12}O_{6} + O_{2}\longrightarrow 2CH_{3}COOH+2CO_{2}+2H_{2}O$). Since $e_{f}$ is lower than $e_{r}$, $CUE$ is lower when fermentative pathways is 
used instead of respiratory pathway.

Under fully aerobic conditions (assuming complete absence of anaerobic microsites in soil), both $CO_{2}$ and acetic acid can be produced concurrently. Production of acetic acid increases with increasing growth rate in various microbial species. [@basan_overflow_2015] provided convincing evidence, that the observed relationship between the growth rate and production of acetic acid is caused by changes in proteome allocation along the gradient of $\mu$. Because $\mu$ is function of $G$ in sub-Microbial model, fluxes $J_{E,r}$ and $J_{E,f}$ will be defined as a function of $G$ using the same rational described in [@basan_overflow_2015].

Absolute amount of proteins ($P$ in C-mol units assuming same C:N of all proteins) in soil microbial biomass can be defined using an equation:

$P=B \times (n_{B-proteins} + n_{G-proteins} \times G)$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S26).

Non-growing microbial biomass with $G = 0$ (i.e. $B \times n_{B-proteins}$) contains three different fractions of proteome - proteins of respiratory pathway ($\phi_{r, \theta}$), fermentative pathway ($\phi_{f, \theta}$) and other proteins ($\phi_{0, \theta}$, i.e. mediating growth, defense, production of extracellular enzymes, etc.). The symbol $\theta$ denotes that these fractions sum up to one only at the non-growing conditions.

Once the organic substrate is assimilated into $G$, proteome fractions start to change because none of the proteins in pool $B \times G \times n_{G-proteins}$ belong to respiratory or fermentative pathway. The optimality of proteome allocation lies in the fact that microbial biomass start to allocate more proteins into ribosomes to allow rapid growth when there is enough organic substrate in the environment. For that reason, RNA content of microbial cells typically increases with increasing growth rate at the conditions of sufficient organic C supply [e.g. @herbert_chemical_1961, @Makino2003]. Similarly as PLFA content of $MBC$, proteome fractions of respiratory and fermentative pathways decrease with increasing $G$. Along the gradient of $G$, fractions $\phi_{r}$ and $\phi_{f}$ (without $\theta$ symbol) can be defined as (Fig. S1):

 $\phi_{r}=\phi_{r, \theta} \times \frac{n_{B-proteins}}{n_{B-proteins} + n_{G-proteins} \times G}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S27),

 $\phi_{f}=\phi_{f, \theta} \times \frac{n_{B-proteins}}{n_{B-proteins} + n_{G-proteins} \times G}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S28).

Notice that in contrast to PLFA content, $\phi_{r}$ and $\phi_{f}$ are normalized to unit of total amount of proteins in microbial biomass, not C content of microbial biomass.

![image](FigS1.pdf) 
*Figure S1: Change in proteome fractions responsible for energy generation by respiratory ($\phi_{r}$) and fermentative pathway ($\phi_{f}$) along the gradient of pool of growth intermediates ($G$), which is scaled to its maximum ($G_{m}$). The red dashed line denotes the maximum proteome fractions at the zero growth rate ($\phi_{r, \theta}$ and $\phi_{f, \theta}$).*

Fractions $\phi_{r}$ and $\phi_{f}$ are responsible for fluxes of energy $J_{E,r}$ and $J_{E,f}$, respectively, according to following equations:

$J_{E,r}=\phi_{r} \times \epsilon_{r}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S29),

$J_{E,f}=\phi_{f} \times \epsilon_{f}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S30).

In these equations, $\epsilon_{r}$ and $\epsilon_{f}$ denote energy production rates per respective proteome fractions $\phi_{r}$ and $\phi_{f}$. 

Combining eqs. S27, S28, S29 and S30, and defining a fraction of energy-related proteome $\phi_{E}$ as $\phi_{r,\theta} + \phi_{f, \theta}$ (Fig. S1), fluxes $J_{E,r}$ and $J_{E,f}$ can be defined as a function of $G$ using an equation:

$\phi_{r} + \phi_{f} = \frac{J_{E,r}}{\epsilon_{r}} + \frac{J_{E,f}}{\epsilon_{r}} = 
\phi_{E} \times \frac{n_{B-proteins}}{n_{B-proteins} + n_{G-proteins} \times G}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S31).

Assuming $n_{B-proteins}$ and $n_{G-proteins}$ to be same [@hanegraaf2000], eq. S31 simplifies to:

$\frac{J_{E,r}}{\epsilon_{r}} + \frac{J_{E,f}}{\epsilon_{r}}=\frac{\phi_{E}}{1 + G}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S32).

Eq. S32 can be used to express flux $J_{E,r}$ as function of the flux $J_{E,f}$ and vice versa. Using eqs. S16 and S23, fluxes $J_{E,r}$ and $J_{E,f}$ can be defined in respect to flux of mobilized growth intermediates $J_{gluc}$ using equations:

$J_{E,r} =\frac{-J_{gluc} + {\epsilon_{f} \times \phi_{E} \times (e_{f}^{-1} + \sigma^{-1})} / {(1+G)}}{{\epsilon_{f}}/{\epsilon_{r}} \times (e_{f}^{-1} + \sigma^{-1}) - e_{r}^{-1} - \sigma^{-1}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S33),

$J_{E,f} =\frac{J_{gluc} - {\epsilon_{r} \times \phi_{E} \times (e_{r}^{-1} + \sigma^{-1})} / {(1+G)}}{e_{f}^{-1} + \sigma^{-1} - {\epsilon_{r}}/{\epsilon_{f}} \times (e_{r}^{-1} + \sigma^{-1})}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S34).

Specific growth rate $\mu$ then follows eq. S16. 
Eqs. S33 and S34 cannot be applied across the entire gradient of $G$ because $J_{E,f}$ becomes negative when $G$ is low (slow growth rate) while $J_{E,r}$ becomes negative when $G$ approaches maximum (Fig. S2A). Therefore, there is a range of $G$ content at which $J_{E,r}$ and $J_{E,f}$ changes with $G$ according to eqs. S33 and S34. Outside that range, microbial biomass either respire or ferment only (Fig. S2A). This range can be found between the values of $G$ at which $J_{E,r}$ and $J_{E,f}$ equals zero. Because $J_{gluc}$ equals $v_{gluc} \times G$, positive root of quadratic equation have to found. $J_{E,f}$ is zero or positive when $G$ increases to a value denoted as $G_{f,0}$ defined by equation (Fig. S2A):

$G_{f,0} = \frac{-v + \sqrt(v \times (v + 4 \times \epsilon_{r} \times \phi_{E} \times (e_{r}^{-1} + \sigma^{-1})))}{2 \times v}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S35).

$J_{E,r}$ then becomes zero when $G$ increases further to a value denoted as $G_{r,0}$, which is defined by equation (Fig. S2A):

$G_{r,0} = \frac{v - \sqrt(v \times (v + 4 \times \epsilon_{f} \times \phi_{E} \times (e_{f}^{-1} + \sigma^{-1})))}{-2 \times v}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S36).

If fermentable substrate is used for growth, $\mu$ along the entire gradient of $G$ follows the equation (Fig. S2B):

$\mu = \begin{cases} {J_{E,r} \times \sigma^{-1}}~~~if~~~G \leq G_{f,0}
                \\\\\ {(J_{E,r} + J_{E,f}) \times \sigma^{-1}}~~~if~~~G > G_{f,0}~~and ~~G \leq G_{r,0}
                \\\\\ {J_{E,f} \times \sigma^{-1}}~~~if~~~G > G_{r,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S37).


Rates of $CO_{2}$ and acetic acid production are then (Fig. S2B):

$\frac{dCO_{2}}{dt} = U_{gluc} \times (1 - y_{A, gluc}) +\begin{cases} {B \times J_{E,r} \times e_{r}^{-1}}~~~if~~~G \leq G_{f,0}
                \\\\\ {B \times (J_{E,r} \times e_{r}^{-1} + \frac{1}{3}J_{E,f} \times e_{f}^{-1})}~~~if~~~G > G_{f,0}~~and ~~G \leq G_{r,0}
                \\\\\ {B \times \frac{1}{3}J_{E,f} \times e_{f}^{-1}}~~~if~~~G > G_{r,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S38),

$\frac{dS_{acet}}{dt} = \begin{cases} {0}~~~if~~~G \leq G_{f,0}
                           \\\\\ {B \times \frac{2}{3}J_{E,f} \times e_{f}^{-1}}~~~if~~~G > G_{f,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S39).

The value of $CUE$ depends on the proportion between fluxes $J_{E,r}$ and $J_{E,f}$ according to equation (Fig S2C):

$CUE = \begin{cases} \frac{1}{1 + \frac{\sigma}{e_{r}}}~~~if~~~G \leq G_{f,0}
                \\\\\ {\frac{J_{E,r} + J_{E,f}}{J_{E,r} \times (1 + \frac{\sigma}{e_{r}}) + J_{E,f} \times (1 + \frac{\sigma}{e_{f}})}}~~~if~~~G > G_{f,0}~~and ~~G \leq G_{r,0}
                \\\\\ \frac{1}{1 + \frac{\sigma}{e_{f}}}~~~if~~~G > G_{r,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S40).

Eq. S40 implies that the CUE is entirely determined by respiratory pathway ($J_{E,r}$) when $G$ is low (Fig. S2A). Because $e_{r}$ is lower than $e_{f}$, $CUE$ decreases as $G$ increases and reaches its minimum at $G_{r,0}$ (Fig. S2C). Above $G_{r,0}$, $CUE$ is entirely determined by fermentative pathway. It must be noted that product of fermentative pathway doesn't need to be solely acetic acid but other organic compounds such as lactic acid or ethanol. For the sake of simplicity, we assume that only acetic acid is produced. It also must be noted that $CUE$ can be determined only if $CO_{2}$ production as well as the change of acetic acid concentration over time is measured concurrently. If the acetic acid is not measured and "apparent $CUE$" is calculated from biomass and $CO_{2}$ production, $CUE$ may seem to increase along the gradient of increasing growth rate (Fig. S2C). 

![image](FigS2.pdf) 
*Figure S2: A) Biomass-specific fluxes of energy generated by respiratory (black solid line, $J_{E, r}$) and fermentative (grey solid line, $J_{E, f}$) pathway along the gradient of pool of growth intermediates ($G$), which is scaled to its maximum ($G_{m}$). B) Biomass-specific rates of $CO_2$ (black solid line), acetic acid (grey solid line, $CH_{3}COOH$), and biomass (red solid line, $\mu$) production. C) Change of realized (black solid line) and apparent (black solid line) carbon use efficiency (CUE). Realized CUE is described by eq. S40, while apparent CUE is defined as growth per sum of growth and $CO_{2}$ production. Both CUE definitions assume $y_{A, gluc} = 1$ for simplicity. In all graphs, the two vertical red dashed lines denote two thresholds defined by eqs. S35 and S36. Notice that panel B has two y-axes with different scales.*

### 1. 3. 3. Combination of fermentable and non-fermentable substrate (soil supplemented by glucose)

When soil is supplied by glucose, the product of its metabolism could be acetic acid. When acetic acid is produced, both fermentable and non-fermentable organic substrate becomes present in soil. Both can be consumed and metabolized at the same time. Two issues arising from such situation must be considered. First, microbial biomass may prefer one substrate over the other [@LaCecilia2018], e.g. glucose over acetatic acid/acetate. This can be acknowledged directly in the denominator of functions $f_{gluc}$ and $f_{acet}$ as described in [@LaCecilia2018]. Second, two substrates have different uptake and metabolc pathways so they cannot be mixed in a single pool $G$. For that reason, two different pools of $G$ are defined - $G_{acet}$ and $G_{gluc}$. $G$ is a sum of $G_{acet}$ and $G_{gluc}$.

The soil supplemented by glucose can be therefore described by following equations:

$\frac{dS_{gluc}}{dt} =-U_{gluc} = -I_{m} \times B \times f_{gluc}; ~~~~ f_{gluc} = \frac{S_{gluc}}{K_{m, gluc} \times (1 + \frac{S_{acet}}{K_{m, acet}}) + S_{gluc}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S41),

$\frac{dS_{acet}}{dt} = -U_{acet} + \begin{cases} {0}~~~if~~~G_{gluc} \leq G_{gluc - f,0}
                           \\\\\ {B \times \frac{2}{3}J_{E,f} \times e_{f}^{-1}}~~~if~~~G_{gluc} > G_{gluc - f,0} \end{cases};
                           U_{acet} = I_{m} \times B \times \frac{S_{acet}}{K_{m, acet} \times (1 + \frac{S_{gluc}}{K_{m, gluc}}) + S_{acet}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S42),

$\frac{dG_{gluc}}{dt} = U_{gluc} \times \frac{y_{A, gluc}}{B} - v_{gluc} \times G_{gluc} - \frac{dB}{dt} \times \frac{G_{gluc}}{B}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S43),

$\frac{dG_{acet}}{dt} = U_{acet} \times \frac{y_{A, acet}}{B} - v_{acet} \times G_{acet} - \frac{dB}{dt} \times \frac{G_{acet}}{B}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S44),

$\frac{dB}{dt} = B \times \frac{v_{acet} \times G_{acet}}{1 + \frac{\sigma}{e_{acet}}} - H + 
                                                        \begin{cases} {B \times J_{E,r} \times \sigma^{-1}}~~~if~~~G_{gluc} \leq G_{gluc-f,0}
                                                            \\\\\ {B \times (J_{E,r} + J_{E,f}) \times \sigma^{-1}}~~~if~~~G_{gluc} > G_{gluc-f,0}~~and ~~G_{gluc} \leq G_{gluc-r,0}
                                                            \\\\\ {B \times J_{E,f} \times \sigma^{-1}}~~~if~~~G_{gluc} > G_{gluc-r,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S45),

$\frac{dCO_{2}}{dt} = U_{gluc} \times (1 - y_{A, gluc}) + U_{acet} \times (1 - y_{A, acet}) + B \times \frac{v_{acet} \times G_{acet}}{1 + \frac{e_{acet}}{\sigma}} +
                                                            \begin{cases} {B \times J_{E,r} \times e_{r}^{-1}}~~~if~~~G_{gluc} \leq G_{gluc-f,0}
                                                            \\\\\ {B \times (J_{E,r} \times e_{r}^{-1} + \frac{1}{3}J_{E,f} \times e_{f}^{-1})}~~~if~~~G_{gluc} > G_{gluc-f,0}~~and ~~G_{gluc} \leq G_{gluc-r,0}
                                                            \\\\\ {B \times \frac{1}{3}J_{E,f} \times e_{f}^{-1}}~~~if~~~G_{gluc} > G_{gluc-r,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S46).

Biomass specific fluxes $J_{E, r}$ and $J_{E,f}$ are defined in respect to $G_{gluc}$ as well as $G_{acet}$ because the proteome fractions change along the gradient of $G$, which is defined as $G_{gluc} + G_{acet}$ when fermentable and non-fermentable substrate can be consumed by soil microorganisms at the same time:

$J_{E,r} = \begin{cases} {\frac{v_{gluc} \times G_{gluc}}{\sigma^{-1} + e_{r}^{-1}}}~~~if~~~G_{gluc} \leq G_{gluc-f,0}
             \\\\\ {\frac{-v_{gluc} \times G_{gluc} + {\epsilon_{f} \times \phi_{E} \times (e_{f}^{-1} + \sigma^{-1})} / {(1+G)}}{{\epsilon_{f}}/{\epsilon_{r}} \times (e_{f}^{-1} + \sigma^{-1}) - e_{r}^{-1} - \sigma^{-1}} }~~~if~~~G_{gluc} > G_{gluc-f,0}~~and ~~G_{gluc} \leq G_{gluc-r,0}
             \\\\\ {0}~~~if~~~G_{gluc} > G_{gluc-r,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S47),

$J_{E,f} = \begin{cases} {0}~~~if~~~G_{gluc} \leq G_{gluc-f,0}
             \\\\\ {\frac{v_{gluc} \times G_{gluc} - {\epsilon_{r} \times \phi_{E} \times (e_{r}^{-1} + \sigma^{-1})} / {(1+G)}}{e_{f}^{-1} + \sigma^{-1} - {\epsilon_{r}}/{\epsilon_{f}} \times (e_{r}^{-1} + \sigma^{-1})} }~~~if~~~G_{gluc} > G_{gluc-f,0}~~and ~~G_{gluc} \leq G_{gluc-r,0}
             \\\\\ {\frac{v_{gluc} \times G_{gluc}}{\sigma^{-1} + e_{f}^{-1}}}~~~if~~~G_{gluc} > G_{gluc-r,0} \end{cases}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S48).

Threshold values $G_{gluc-f,0}$ and $G_{gluc-r,0}$ therefore consider both, $G_{gluc}$ and $G_{acet}$:

$G_{gluc-f,0} = \frac{-v_{gluc} \times (1 + G_{acet}) + \sqrt((v_{gluc} \times (1 + G_{acet}))^{2} + 4 \times v_{gluc}\times \epsilon_{r} \times \phi_{E} \times (e_{r}^{-1} + \sigma^{-1}))}{2 \times v_{gluc}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S49),

$G_{gluc-r,0} = \frac{v_{gluc} \times (1 + G_{acet}) - \sqrt((v_{gluc} \times (1 + G_{acet}))^{2} + 4 \times v_{gluc} \times \epsilon_{f} \times \phi_{E} \times (e_{f}^{-1} + \sigma^{-1}))}{-2 \times v_{gluc}}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S50).

Equations S49 and S50 imply that in contrast to separate glucose metabolism, onset of acetic acid production from assimilated glucose appears at lower glucose concentrations because part of pool $G$ can be occupied by non-fermentable substrate.

## 1.4. Model parameters and predictions

Basal sub-Microbial model described by eqs. S6 - S9 has six parameters - $I_{m}$, $K_{m}$, $y_{A}$, $G_{m}$, $g$, and $h$. When the dynamic of some microbial biomass proxy-parameter over time needs to be known, additional parameters $n_{B-*}$ and $n_{G-*}$ have to added. Meaning of the parameters and their possible values that were determined in our previous study are reported in table S1.

***Tab S1***: Parameters of basal and extended sub-Microbial model. Meaning of the parameters, range of their expected values and their units are reported.  

| Parameter | Explanation | Mean value | Minimum | Maximum | Units |
|:----------|:------------|:----------:|:-------:|:-------:|:------:|
| basal sub-Microbial model                                      ||
|$I_{m}$    |Maximum specific substrate uptake rate    |6.49|0.13|16.1|$\mu mol(C)~\mu mol (C~of~pool~B)^{-1} day^{-1}$|
|$K_{m}$    |Half-saturation constant                  |221|0.1|393|$\mu mol(C)~g(DW)^{-1}$                            |
|$y_{A}$    |Assimilation efficiency                   |1|0.24|1|unitless                                              |
|$G_{m}$    |Specific maximum content of $G$ within $B$|3.01|0.21|14.5|$\mu mol(C)~\mu mol (C~of~pool~B)^{-1}$         |
|$g$        |Cost of biomass formation                 |2.07|0.01|4.67|$\mu mol(C)~\mu mol (C~of~pool~B)^{-1}$         |
|$h$        |Specific decay rate                       |-|-|-|$day^{-1}$                                               |
| Conversion factors                                             ||
|${}^{\dagger}n_{B-CLC}$        |Relative amount of chloroform-labile organic C in pool $B$                       |0.24|0.04|0.52|unitless  |
|${}^{\dagger}n_{G-CLC}$        |Relative amount of chloroform-labile organic C in pool $G$                       |0.41|0.34|0.74|unitless  |
|$n_{B-PLFA}$       |Relative amount of PLFA in pool $B$                                              |$8.51 \times 10^{-3}$|-|-|unitless  |
| extended sub-Microbial model (energy-mass dependency)                                     ||
|${}^{\ddagger}K_{m, gluc}$    |Half-saturation constant of glucose uptake   |95.5|1.37|500|$\mu mol(C_{gluc})~g(DW)^{-1}$|
|${}^{\ddagger}K_{m, acet}$    |Half-saturation constant of acetate uptake   |201|0.60|1000|$\mu mol(C_{acet})~g(DW)^{-1}$|
|$\sigma$    |Energy requirement of biomass formation    |1.6|0.88|-|$\mu mol(ATP)~\mu mol (C~of~pool~B)^{-1}$|
|$\psi$  |ATP yield per one reducing equivalent    |2|0|3|unitless|
|$e_{acet}$  |Energy formed from assimilated acetic acid ($0.5 + \frac{3.5}{2} \times \psi$)   |4|0.5|5.75|$\mu mol(ATP)~\mu mol (C~of~pool~B)^{-1}$|
|$e_{r}$  |Energy formed from assimilated glucose by respiratory pathway ($1 + \frac{11}{6} \times \psi$)     |4.67|1|6.5|$\mu mol(ATP)~\mu mol (C~of~pool~B)^{-1}$|
|$e_{f}$  |Energy formed from assimilated glucose by fermentative pathway ($1 + \frac{2}{3} \times \psi$)   |2.33|1|3|$\mu mol(ATP)~\mu mol (C~of~pool~B)^{-1}$|
|$y_{A, gluc}$  |Assimilation efficiency of glucose ($1 - \frac{1}{3} \times e_{r}^{*}{}^{-1}$)    |0.93|0.67|0.95|unitless|
|$y_{A, acet}$  |Assimilation efficiency  of acetic acid/acetate ($1 - \frac{1}{2} \times e_{acet}^{*}{}^{-1}$)    |0.88|0|0.91|unitless|
|${}^{\dagger}\frac{\epsilon_{f}}{\epsilon_{r}}$  |ratio between ATP production rate per unit proteome fraction of fermentative and respiratory pathway     |$\frac{75}{39}$|-|-|unitless|
|${}^{\ddagger}\phi_{E} \times \epsilon_{r}$  |Specific energy production rate per unit fraction of proteome responsible for energy production at zero growth rate     |1.12|0.05|10.0|$\mu mol(ATP)~\mu mol (C~of~pool~B)^{-1} ~ day^{-1}$|
|$\phi_{E} \times \epsilon_{f}$  |$\phi_{E} \times \epsilon_{r} \times \frac{\epsilon_{f}}{\epsilon_{r}}$     |-|-|-|$\mu mol(ATP)~\mu mol (C~of~pool~B)^{-1} ~ day^{-1}$|

*Symbol ${}^{\dagger}$ denotes parameters, which are considered constant*
*Symbol ${}^{\ddagger}$ denotes parameters estimated using data presented in [@santruckova_short-term_2004]*

When dependency of energy and mass fluxes are explicitly specified in order to predict production of acetic acid in soil, additional parameters are required - $K_{m, gluc}$, $K_{m, acet}$, $y_{A, gluc}$, $y_{A, acet}$, $\sigma$, $e_{acet}$, $e_{r}$, $e_{f}$, $\epsilon_{r}$, $\epsilon_{f}$, and $\phi_{E}$. Some of these parameters are reported for *E. Coli* in [@basan_overflow_2015]. The total number of parameters whose value is expected to vary and thus, need to be optimized, can be reduced.

Parameters $\epsilon_{r}$ and $\epsilon_{f}$ can be considered as fixed. Their values has been estimated to be approx. 2.55 and 1.33 $\mu mol(ATP)~\mu mol (C~of~pool~B)^{-1} ~ day^{-1}$ (the production rates are "biomass specific"). There is, however, uncertainty surrounding the conversion between original units normalized to $OD_{600}$ in [@basan_overflow_2015] and C atoms of pool $B$ in our study. For that reason, we use unitless ratio $\epsilon_{f}$ to $\epsilon_{r}$ and its inverse, which are defined in denominators of eqs. S47 and S48 (Table S1). Mean value of parameter $\phi_{E}$ is expected to be 0.19 [@basan_overflow_2015], i.e. 19% of entire proteome at zero growth rate should be responsible for energy production. [@basan_overflow_2015], however, showed that this parameter could vary. Moreover, the $\phi_{E}$ could be manipulated experimentally. To account for variability in $\phi_{E}$, $\phi_{E} \times \epsilon_{r}$ from eq. S48 is defined as one composite model parameter reflecting fraction of proteome invested into the energy yield at zero growth rate. $\phi_{E} \times \epsilon_{f}$  in eq. S47 is then simply $\phi_{E} \times \epsilon_{r} \times \frac{\epsilon_{f}}{\epsilon_{r}}$.

The energy requirement of biomass formation $\sigma$ is also expected to vary. According to [@basan_overflow_2015], $\sigma$ is 1.6  $\mu mol$ ATP per $\mu mol$ of biomass C of *E. Coli*. Theoretical calculations conducted by [@stouthamer_theoretical_1973], however, suggest lower value - 0.88 $\mu mol$ ATP per $\mu mol$ of biomass C. It must be noted, that neither of these estimates consider two pools of microbial biomass so the published estimates may not strictly apply. 

The value of parameters $e_{acet}$, $e_{r}$ and $e_{f}$ can be estimated indirectly from known amount of ATP produced by either substrate-level phosphorylation or ATP synthase as explained in [@basan_overflow_2015]. The amount of ATP produced by substrate-level phosphorylation per one C-mol of glucose and acetate is fixed. In case of glucose, it is one ATP molecule per one C-mol glucose processed by either fermentative and respiratory pathway. For acetate, it is only 0.5 (i.e. 1 mol of ATP in TCA cycle per one mol of acetate). The amount of ATP produced by ATP synthase could be variable because it depends on amount of reducing equivalents, which are produced by respective metabolc pathway, stoichiometry of $e^{-}$ and $H^{+}$ transport [@Unden], membrane resistance (i.e. proton leakage, [@cook_energy-spilling_1994]), and electron acceptors that have been used (under anaerobic conditions in particular). To calculate $e_{acet}$, $e_{r}$ and $e_{f}$, a new parameter, which will be abbreviated as $\psi$ hereafter, is introduced. This parameter defines the ATP yield per mol of reducing equivalent and could range between 0 and 3 (Tab. S1). $\psi$ multiplies the number of catabolic pathway specific amount of produced reducing equivalents per one C-mol of glucose and acetate [@basan_overflow_2015]. $e_{acet}$, $e_{r}$ and $e_{f}$ are defined by following equations: 

$e_{acet}=0.5 + \frac{3.5}{2} \times \psi$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S51),

$e_{r}=1 + \frac{11}{6} \times \psi$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S52),

$e_{f}=1 + \frac{2}{3} \times \psi$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S53).

Eqs. S51 - S52 define the entire range of $e_{acet}$, $e_{r}$ and $e_{f}$ (Tab. S1). Applying mean values for $\sigma$, $e_{acet}$, $e_{r}$ and $e_{f}$ as reported in Table S1, CUE is expected to be 0.71, 0.74 and 0.59 for acetate metabolism, respiratory glucose metabolism and fermentative glucose metabolism, respectively. The theoretical maximum is reached when $\sigma$ is minimal while $e_{r}$ is maximal - i.e. 0.88. Vice versa, the theoretical minimum is 0.24 when acetate is metabolized without any ATP being produced by ATP synthase using proton motive force. 

Assimilation efficiencies $y_{A, gluc}$ and $y_{A, acet}$ could be calculated indirectly using similar considerations as described above acknowledging that 2 molecules of ATP are required to assimilate glucose into 2 molecules of glyceraldehyde 3-phosphate, and 1 molecule of ATP is required to assimilate acetate into acetyl CoA: 

$y_{A, gluc}=1 - \frac{1}{3} \times e_{r}^{*}{}^{-1}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S54),

$y_{A, acet}=1 - \frac{1}{2} \times e_{acet}^{*}{}^{-1}$ $~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ (eq. S55).


For mean $e_{r}$ 4.67 $y_{A, gluc}$ is expected to be 0.93. The theoretical minimum and maximum $y_{A, gluc}$ is given by minimum and maximum values of $e_{r}$, thus 0.67 and 0.95, respectively. For acetate assimilation, mean, minimum and maximum $y_{A, acet}$ is 0.88, 0 and 0.91, respectively. It must be noted, however, that the definition of $y_{A}$ is not entirely accurate because the ATP, which is used for assimilation of glucose or acetate is intracellular, not directly derived from consumed organic substrate. If the original definition is, nevertheless, applied, $y_{A, gluc} = 1$ determined in our previous study, could be theoretically reached only if the phosphorylation of glucose is mediated by polyphosphates, not by ATP directly. 

# 2. Comparison of basal and extended sub-Microbial model

The comparison is performed on data published in study of [@santruckova_short-term_2004]. In this study, ${}^{14}C$ glucose was added to soil and the radioactivity was measured in various organic C pools in soil and respired $CO_{2}$. In addition, glucose concentration in soil was measured directly. This allowed authors to distinguish between residual glucose concentration and the amount of extractractable organic C originating from ${}^{14}C$ glucose metabolism (i.e. not retained in microbial biomass). This "excreted" organic C was defined as fermentation product. 

The basal sub-Microbial model was fitted onto the data ignoring existence of fermentation products. Using 7 optimized parameters ($I_{m}$, $K_{m, gluc}$, $y_{A, gluc}$, $G_{m}$, $g$, $n_{B-CLC}$ and $n_{G-CLC}$), sub-Microbial model was able to simulate measured data with high accuracy. Due to the short-term laboratory incubation, the death rate of soil microbial biomass was considered negligible. Here, we compare the original goodness of fit with goodness of fit of extended sub-Microbial model with 9 parameters ($I_{m}$, $K_{m, gluc}$, $K_{m, acet}$, $G_{m}$, $\sigma$, $\psi$, $\phi_{E} \times \epsilon_{r}$, $n_{B-CLC}$ and $n_{G-CLC}$).     

![image](FigS3.pdf) 
*Figure S3: The dynamic of cumulative $CO_{2}$ production (A), radioactively labelled glucose concentration (B), chloroform-labile organic carbon ($CHCl_3$ flush; C), conversion factor between $CHCl_3$ flush and total soil microbial biomass C ($k_{CLC}$; D), and concentration of fermentation products in soil (E) supplied by glucose. Data are reported in [@santruckova_short-term_2004]. The grey solid line and black dashed line denote simulations of basal and extended sub-Microbial model.*

***Tab S2***: Correspondence between data presented in [@santruckova_short-term_2004], and basal and extended sub-Microbial model (see Fig S3). The goodness of correspondence is expressed as coefficient of determination ($R^{2}$) and it is calculated for five different variables separately.  

| Measured variable | basal sub-Micrbial model | extended sub-Micrbial model |
|:------------------|:-------------------------|:---------------------------:|
|Cumulative respiration|0.996|0.992|
|Glucose concentration|0.979|0.982|
|Chloroform labile organic C|0.996|0.996|
|$k_{CLC}$|0.996|0.991|
|Fermentation products|-|0.804|

***Tab S3***: Parameters of basal and extended sub-Microbial model as optimized using data presented in [@santruckova_short-term_2004]. Mean values and standard deviations are reported for each parameter. Meaning of the parameters and their units are reported in Table S1. 

| Model parameter | basal sub-Micrbial model | extended sub-Micrbial model |
|:------------------|:-------------------------|:---------------------------:|
|$I_{m}$|$6.79 \pm 4.41$|$6.41 \pm 4.01$|
|$K_{m, gluc}$|$123 \pm 187$|$95.5 \pm 116$|
|$K_{m, acet}$|-|$201 \pm 200$|
|$y_{A, gluc}$|$0.99 \pm 0.08$|-|
|$G_{m}$|$1.46 \pm 3.28$|$1.66 \pm 10.0$|
|$g$|$0.34 \pm 0.19$|-|
|$\sigma$|-|$0.92 \pm 1.54$|
|$\psi$|-|$1.30 \pm 0.67$|
|$\phi_{E} \times \epsilon_{r}$|-|$1.12 \pm 3.05$|
|$n_{B-CLC}$|$0.22 \pm 0.06$|$0.22 \pm 0.07$|
|$n_{G-CLC}$|$0.50 \pm 0.11$|$0.45 \pm 0.07$|



Fig. S3 and Table S2 show that extended sub-Microbial model fits all measured variables equally good as original basal sub-Microbial model. In addition, extended model simulates concentration of fermentation products in soil with high accuracy ($R^{2} = 0.8$). When goodness of correspondence of the two models are compared across all measured variables (i.e. including fermentation products) statistically as desribed in [@capek_revisiting_2023], the extended sub-Microbial model performs significantly better ($F_{2, 30} = 4.58, p = 0.022$) regardless of two extra parameters. In fact, the total number of parameters (9) explaining variability in five different measured variables with high accuracy is less than if simple linear regression with two parameters is conducted for each variable. Importantly, estimates of parameters of basal and extended model with the same meaning ($I_{m}$, $G_{m}$, $n_{B-CLC}$ and $n_{G-CLC}$) are indistinguishable. 