# sif4 baseline after reduced

![graph](eval53viz_sif4_large_baseline_after_reduced_simple.png)

- summary nodes: 29
- summary reactions: 146
- drawn nodes: 22
- drawn edges: 71
- colors: gas=blue, surface=orange, bulk/mixed=green

## N1 (orange)

Names: t_HN_NH2(S)

Reactions:
- R350: NH3 + t_HN_SIF(S) => HF + SI(D) + t_HN_NH2(S)
- R351: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S)

## N2 (orange)

Names: 5 merged: t_HN_SIF(S), t_F3SI_NH2(S), t_F2SINH(S), t_H2NFSINH(S), t_HN(FSINH)2(S)

Reactions:
- R350: NH3 + t_HN_SIF(S) => HF + SI(D) + t_HN_NH2(S)
- R351: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S)
- R349: t_F2SINH(S) + t_H2NFSINH(S) => HF + t_HN(FSINH)...

## N3 (orange)

Names: s_HN_NH2(S)

Reactions:
- R356: NH3 + s_HN_SIF(S) => HF + SI(D) + s_HN_NH2(S)
- R357: SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S)

## N4 (orange)

Names: 5 merged: s_HN_SIF(S), s_F3SI_NH2(S), s_F2SINH(S), s_H2NFSINH(S), s_HN(FSINH)2(S)

Reactions:
- R356: NH3 + s_HN_SIF(S) => HF + SI(D) + s_HN_NH2(S)
- R357: SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S)
- R355: s_F2SINH(S) + s_H2NFSINH(S) => HF + s_HN(FSINH)...

## N5 (orange)

Names: k_HN_NH2(S)

Reactions:
- R362: NH3 + k_HN_SIF(S) => HF + SI(D) + k_HN_NH2(S)
- R363: SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)

## N6 (orange)

Names: 5 merged: k_HN_SIF(S), k_F3SI_NH2(S), k_F2SINH(S), k_H2NFSINH(S), k_HN(FSINH)2(S)

Reactions:
- R362: NH3 + k_HN_SIF(S) => HF + SI(D) + k_HN_NH2(S)
- R363: SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)
- R361: k_F2SINH(S) + k_H2NFSINH(S) => HF + k_HN(FSINH)...

## N7 (blue)

Names: N2

Reactions:
- R2, 12, 206: N + NH <=> H + N2 | N + NH2 <=> 2 H + N2 | N + NO <=> N2 + O
- R206, 211, 213, 224: N + NO <=> N2 + O | H + N2O <=> N2 + OH | N2O (+M) <=> N2 + O (+M) | NH + NO <=> N2 + OH
- R7, 8, 211, 224, 230, 234: H + NNH <=> H2 + N2 | NH2 + NNH <=> N2 + NH3 | H + N2O <=> N2 + OH | NH + NO <=> N2 + OH | NNH + O2 <=> HO2 + N2 | CH3 + NNH <=> CH4 + N2
- R230: NNH + O2 <=> HO2 + N2

## N8 (blue)

Names: N

Reactions:
- R3, 220: H + NH <=> H2 + N | NH + OH <=> H2O + N
- R12, 298, 299: N + NH2 <=> 2 H + N2 | CH3 + N <=> H + H2CN | CH3 + N <=> H2 + HCN
- R2, 12, 206: N + NH <=> H + N2 | N + NH2 <=> 2 H + N2 | N + NO <=> N2 + O
- R298, 299: CH3 + N <=> H + H2CN | CH3 + N <=> H2 + HCN
- R207: N + O2 <=> NO + O
- R206: N + NO <=> N2 + O

## N9 (blue)

Names: 10 merged: H2, H, NH, NH2, NNH, N2H2, N2H3, N2H4, HF, NH3

Reactions:
- R351, 357, 363: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S) | SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S) | SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)
- R363: SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)
- R357: SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S)
- R351: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S)
- R65, 70, 76, 91, 92, 211, 222, 224, 230: H + O2 + M <=> HO2 + M | H + O2 <=> O + OH | H + HO2 <=> 2 OH | CH2OH + H <=> CH3 + OH | CH2OH + H <=> CH2(S) + H2O | H + N2O <=> N2 + OH
- R361: k_F2SINH(S) + k_H2NFSINH(S) => HF + k_HN(FSINH)...
- R355: s_F2SINH(S) + s_H2NFSINH(S) => HF + s_HN(FSINH)...
- R349: t_F2SINH(S) + t_H2NFSINH(S) => HF + t_HN(FSINH)...
- R3, 220: H + NH <=> H2 + N | NH + OH <=> H2O + N
- R42, 53, 158, 163, 166, 172, 187, 272, 298, 299, 306, 307: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | CH + CH4 <=> C2H4 + H | CH2 + O2 => CO + H + OH | CH2 + CH3 <=> C2H4 + H | CH2(S) + O2 <=> CO + H + OH
- R12, 298, 299: N + NH2 <=> 2 H + N2 | CH3 + N <=> H + H2CN | CH3 + N <=> H2 + HCN
- R7, 8, 211, 224, 230, 234: H + NNH <=> H2 + N2 | NH2 + NNH <=> N2 + NH3 | H + N2O <=> N2 + OH | NH + NO <=> N2 + OH | NNH + O2 <=> HO2 + N2 | CH3 + NNH <=> CH4 + N2
- R42, 53, 163, 172, 254, 306, 307, 311: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | HCN + O <=> H + NCO | CH3 + O => CO + H + H2
- R194, 195: HCO + H2O <=> CO + H + H2O | HCO + M <=> CO + H + M
- R221, 235: NH + O2 <=> HNO + O | H + NO + M <=> HNO + M
- R91, 92, 109, 111: CH2OH + H <=> CH3 + OH | CH2OH + H <=> CH2(S) + H2O | H + HCCO <=> CH2(S) + CO | CH2CO + H <=> CH3 + CO
- R222: NH + O2 <=> NO + OH
- R109, 111: H + HCCO <=> CH2(S) + CO | CH2CO + H <=> CH3 + CO
- R128, 309: CO + OH <=> CO2 + H | CH3 + OH => CH2O + H2
- R254: HCN + O <=> H + NCO
- R272: CH2 + NO <=> H + HNCO

## N10 (blue)

Names: 2 merged: SIF4, SIF3

Reactions:
- R351, 357, 363: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S) | SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S) | SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)
- R363: SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)
- R357: SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S)
- R351: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S)

## N11 (blue)

Names: 2 merged: O, O2

Reactions:
- R35, 43, 65, 163, 172, 173, 184, 196, 197, 198, 203, 204: H2 + O <=> H + OH | CH4 + O <=> CH3 + OH | H + O2 + M <=> HO2 + M | CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O
- R42, 53, 163, 172, 254, 306, 307, 311: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | HCN + O <=> H + NCO | CH3 + O => CO + H + H2
- R55, 163, 172, 173, 196, 204, 306, 311: C2H2 + O <=> CH2 + CO | CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O | HCO + O2 <=> CO + HO2 | HCCO + O2 <=> 2 CO + OH
- R42, 53, 57, 184, 201, 307: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | C2H4 + O <=> CH3 + HCO | CH3 + O2 <=> CH2O + OH | C2H3 + O2 <=> CH2O + HCO | C2H4 + O <=> CH2CHO + H
- R227: NH2 + O <=> H + HNO
- R206, 213: N + NO <=> N2 + O | N2O (+M) <=> N2 + O (+M)
- R218, 222, 232, 239: NH + O <=> H + NO | NH + O2 <=> NO + OH | NNH + O <=> NH + NO | HNO + O2 <=> HO2 + NO
- R206: N + NO <=> N2 + O
- R75, 116: H + HO2 <=> H2 + O2 | 2 OH <=> H2O + O
- R230: NNH + O2 <=> HO2 + N2
- R254: HCN + O <=> H + NCO

## N12 (blue)

Names: 4 merged: OH, H2O, HO2, H2O2

Reactions:
- R65, 70, 76, 91, 92, 211, 222, 224, 230: H + O2 + M <=> HO2 + M | H + O2 <=> O + OH | H + HO2 <=> 2 OH | CH2OH + H <=> CH3 + OH | CH2OH + H <=> CH2(S) + H2O | H + N2O <=> N2 + OH
- R35, 43, 65, 163, 172, 173, 184, 196, 197, 198, 203, 204: H2 + O <=> H + OH | CH4 + O <=> CH3 + OH | H + O2 + M <=> HO2 + M | CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O
- R211, 224, 278: H + N2O <=> N2 + OH | NH + NO <=> N2 + OH | CH3 + NO <=> H2O + HCN
- R219, 223: NH + OH <=> H + HNO | H2O + NH <=> H2 + HNO
- R163, 172, 173, 184, 278: CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O | CH3 + O2 <=> CH2O + OH | CH3 + NO <=> H2O + HCN
- R91, 92, 196, 204: CH2OH + H <=> CH3 + OH | CH2OH + H <=> CH2(S) + H2O | HCO + O2 <=> CO + HO2 | HCCO + O2 <=> 2 CO + OH
- R75, 116: H + HO2 <=> H2 + O2 | 2 OH <=> H2O + O
- R128, 309: CO + OH <=> CO2 + H | CH3 + OH => CH2O + H2
- R239: HNO + O2 <=> HO2 + NO
- R309: CH3 + OH => CH2O + H2

## N13 (blue)

Names: 2 merged: CO, CO2

Reactions:
- R85, 109, 111, 188, 194, 195, 196, 204: H + HCO <=> CO + H2 | H + HCCO <=> CH2(S) + CO | CH2CO + H <=> CH3 + CO | CH3 + HCO <=> CH4 + CO | HCO + H2O <=> CO + H + H2O | HCO + M <=> CO + H + M
- R55, 163, 172, 173, 196, 204, 306, 311: C2H2 + O <=> CH2 + CO | CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O | HCO + O2 <=> CO + HO2 | HCCO + O2 <=> 2 CO + OH
- R163, 172, 173, 306, 311: CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O | CH3 + O => CO + H + H2 | CH2 + O2 => CO2 + 2 H
- R109, 111: H + HCCO <=> CH2(S) + CO | CH2CO + H <=> CH3 + CO
- R288: H + HNCO <=> CO + NH2
- R168: CH2 + CO (+M) <=> CH2CO (+M)

## N14 (blue)

Names: 13 merged: CH, CH2, CH2(S), CH3, CH4, C2H, C2H2, C2H3, C2H4, C2H5, C2H6, C3H7, C3H8

Reactions:
- R42, 53, 158, 163, 166, 172, 187, 272, 298, 299, 306, 307: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | CH + CH4 <=> C2H4 + H | CH2 + O2 => CO + H + OH | CH2 + CH3 <=> C2H4 + H | CH2(S) + O2 <=> CO + H + OH
- R278, 298, 299: CH3 + NO <=> H2O + HCN | CH3 + N <=> H + H2CN | CH3 + N <=> H2 + HCN
- R42, 53, 147, 153, 168, 183, 184, 201, 307, 309, 312, 315: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | CH3 + HO2 <=> CH3O + OH | CH + O2 <=> HCO + O | CH2 + CO (+M) <=> CH2CO (+M) | CH3 + O2 <=> CH3O + O
- R163, 172, 173, 306, 311: CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O | CH3 + O => CO + H + H2 | CH2 + O2 => CO2 + 2 H
- R163, 172, 173, 184, 278: CH2 + O2 => CO + H + OH | CH2(S) + O2 <=> CO + H + OH | CH2(S) + O2 <=> CO + H2O | CH3 + O2 <=> CH2O + OH | CH3 + NO <=> H2O + HCN
- R91, 92, 109, 111: CH2OH + H <=> CH3 + OH | CH2OH + H <=> CH2(S) + H2O | H + HCCO <=> CH2(S) + CO | CH2CO + H <=> CH3 + CO
- R272: CH2 + NO <=> H + HNCO

## N15 (blue)

Names: 2 merged: NO, N2O

Reactions:
- R206, 211, 213, 224: N + NO <=> N2 + O | H + N2O <=> N2 + OH | N2O (+M) <=> N2 + O (+M) | NH + NO <=> N2 + OH
- R211, 224, 278: H + N2O <=> N2 + OH | NH + NO <=> N2 + OH | CH3 + NO <=> H2O + HCN
- R237, 239: H + HNO <=> H2 + NO | HNO + O2 <=> HO2 + NO
- R207: N + O2 <=> NO + O
- R206, 213: N + NO <=> N2 + O | N2O (+M) <=> N2 + O (+M)
- R218, 222, 232, 239: NH + O <=> H + NO | NH + O2 <=> NO + OH | NNH + O <=> NH + NO | HNO + O2 <=> HO2 + NO
- R222: NH + O2 <=> NO + OH
- R235: H + NO + M <=> HNO + M
- R278: CH3 + NO <=> H2O + HCN
- R272: CH2 + NO <=> H + HNCO

## N16 (blue)

Names: HNO

Reactions:
- R237, 239: H + HNO <=> H2 + NO | HNO + O2 <=> HO2 + NO
- R219, 223: NH + OH <=> H + HNO | H2O + NH <=> H2 + HNO
- R221, 235: NH + O2 <=> HNO + O | H + NO + M <=> HNO + M
- R227: NH2 + O <=> H + HNO
- R235: H + NO + M <=> HNO + M
- R239: HNO + O2 <=> HO2 + NO

## N17 (blue)

Names: 3 merged: HCN, H2CN, HCNN

Reactions:
- R278, 298, 299: CH3 + NO <=> H2O + HCN | CH3 + N <=> H + H2CN | CH3 + N <=> H2 + HCN
- R298, 299: CH3 + N <=> H + H2CN | CH3 + N <=> H2 + HCN
- R278: CH3 + NO <=> H2O + HCN
- R254: HCN + O <=> H + NCO

## N18 (blue)

Names: 3 merged: HCNO, HOCN, HNCO

Reactions:
- R288: H + HNCO <=> CO + NH2
- R289: H + HNCO <=> H2 + NCO
- R272: CH2 + NO <=> H + HNCO

## N19 (blue)

Names: NCO

Reactions:
- R289: H + HNCO <=> H2 + NCO
- R254: HCN + O <=> H + NCO

## N20 (blue)

Names: 10 merged: HCO, CH2O, CH2OH, CH3O, CH3OH, HCCO, CH2CO, HCCOH, CH2CHO, CH3CHO

Reactions:
- R85, 109, 111, 188, 194, 195, 196, 204: H + HCO <=> CO + H2 | H + HCCO <=> CH2(S) + CO | CH2CO + H <=> CH3 + CO | CH3 + HCO <=> CH4 + CO | HCO + H2O <=> CO + H + H2O | HCO + M <=> CO + H + M
- R42, 53, 147, 153, 168, 183, 184, 201, 307, 309, 312, 315: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | CH3 + HO2 <=> CH3O + OH | CH + O2 <=> HCO + O | CH2 + CO (+M) <=> CH2CO (+M) | CH3 + O2 <=> CH3O + O
- R194, 195: HCO + H2O <=> CO + H + H2O | HCO + M <=> CO + H + M
- R42, 53, 57, 184, 201, 307: CH3 + O <=> CH2O + H | C2H2 + O <=> H + HCCO | C2H4 + O <=> CH3 + HCO | CH3 + O2 <=> CH2O + OH | C2H3 + O2 <=> CH2O + HCO | C2H4 + O <=> CH2CHO + H
- R91, 92, 196, 204: CH2OH + H <=> CH3 + OH | CH2OH + H <=> CH2(S) + H2O | HCO + O2 <=> CO + HO2 | HCCO + O2 <=> 2 CO + OH
- R91, 92, 109, 111: CH2OH + H <=> CH3 + OH | CH2OH + H <=> CH2(S) + H2O | H + HCCO <=> CH2(S) + CO | CH2CO + H <=> CH3 + CO
- R168: CH2 + CO (+M) <=> CH2CO (+M)
- R309: CH3 + OH => CH2O + H2

## N21 (green)

Names: SI(D)

Reactions:
- R350: NH3 + t_HN_SIF(S) => HF + SI(D) + t_HN_NH2(S)
- R362: NH3 + k_HN_SIF(S) => HF + SI(D) + k_HN_NH2(S)
- R356: NH3 + s_HN_SIF(S) => HF + SI(D) + s_HN_NH2(S)

## N22 (green)

Names: N(D)

Reactions:
- R351, 357, 363: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S) | SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S) | SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)
- R363: SIF4 + k_HN_NH2(S) => HF + N(D) + k_F3SI_NH2(S)
- R357: SIF4 + s_HN_NH2(S) => HF + N(D) + s_F3SI_NH2(S)
- R351: SIF4 + t_HN_NH2(S) => HF + N(D) + t_F3SI_NH2(S)


SVG: [eval53viz_sif4_large_baseline_after_reduced_simple.svg](eval53viz_sif4_large_baseline_after_reduced_simple.svg)
