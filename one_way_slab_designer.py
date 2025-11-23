# one_way_slab_designer.py
# This is a compact version of the finalized designer including key functions.

from dataclasses import dataclass
import math
from typing import Optional, Dict

UNIT_WEIGHT_BY_FCK = {15:25.0,20:25.0,25:25.0,30:25.0}
XU_MAX_BY_FY = {250:0.53,415:0.48,500:0.46}

def unit_weight_for_fck(fck): 
    keys=sorted(UNIT_WEIGHT_BY_FCK.keys()); return UNIT_WEIGHT_BY_FCK[min(keys,key=lambda k:abs(k-fck))]

def dead_load_self(h_mm, unit_weight=25.0): return (h_mm/1000.0)*unit_weight
def factored_load(dead_kNpm2, live_kNpm2, gamma=1.5): return gamma*(dead_kNpm2+live_kNpm2)
def design_bending_moment(w_kNpm,L_m): return w_kNpm*L_m*L_m/8.0

def mu_limit_coefficient_for_fy(fy):
    nearest=min(sorted(XU_MAX_BY_FY.keys()), key=lambda k:abs(k-fy)); xu_d=XU_MAX_BY_FY[nearest]
    return 0.36*xu_d*(1.0-0.42*xu_d)

def annex_g_solve_Ast(Mu_kNm, fy, fck, b_mm, d_mm):
    Mu_Nmm = Mu_kNm * 1e6
    a = 0.87 * fy * fy / (b_mm * fck)
    b_coeff = -0.87 * fy * d_mm
    c_coeff = Mu_Nmm
    disc = b_coeff*b_coeff - 4.0*a*c_coeff
    if disc < 0: raise ValueError("Annex-G quadratic has negative discriminant")
    r1 = (-b_coeff + math.sqrt(disc))/(2*a); r2 = (-b_coeff - math.sqrt(disc))/(2*a)
    cand=[r for r in (r1,r2) if r>0]
    if not cand: raise ValueError("No positive root for Ast")
    return min(cand)

def shear_capacity_from_As(As_mm2_per_m,b_mm,d_mm,fck):
    if As_mm2_per_m<=0: p_t=0.15
    else: p_t=100.0*As_mm2_per_m/(b_mm*d_mm)
    beta=(0.8*fck)/(6.89*p_t) if p_t>0 else 1e6
    tau_c=0.85*math.sqrt(0.8*fck)*((math.sqrt(1.0+5.0*beta)-1.0)/(6.0*beta))
    tau_c_max=0.625*math.sqrt(fck); tau_c=min(tau_c,tau_c_max)
    return (tau_c*b_mm*d_mm)/1000.0

@dataclass
class OneWaySlabResult:
    short_span_m: float; long_span_m: float; h_mm: float; D_mm: float; d_mm: float
    fck: float; fy: float; unit_weight: float
    dead_self_kNpm2: float; external_dead_kNpm2: float; dead_total_kNpm2: float
    live_kNpm2: float; wu_kNpm: float; Mu_kNm_per_m: float; Mu_lim_kNm_per_m: float
    d_req_from_Mu_lim_mm: Optional[float]; Ast_required_mm2_per_m: float; Ast_min_mm2_per_m: float
    dia_bar_mm: float; spacing_tension_mm: float; Ast_provided_mm2_per_m: float
    dist_spacing_mm: float; Ast_dist_provided_mm2_per_m: float; Vu_kN: float; Vc_kN: float
    shear_ok: bool; L_over_d: float; notes: Dict[str,str]

def run_design(short_span_m,long_span_m,h_mm,fck=20.0,fy=415.0,external_dead_kNpm2=1.5,live_kNpm2=4.0,dia_bar_mm=10.0,dist_bar_mm=8.0,max_agg_mm=20.0,unit_weight_override=None):
    notes={}
    if long_span_m < 2.0*short_span_m: raise ValueError("Not a one-way slab")
    unit_w = unit_weight_override if unit_weight_override is not None else unit_weight_for_fck(fck)
    nominal_cover=20.0
    clear_cover = nominal_cover-5.0 if dia_bar_mm<=12 else nominal_cover
    dead_self = dead_load_self(h_mm, unit_w); dead_total = dead_self + external_dead_kNpm2
    eff_cover = clear_cover + 0.5*dia_bar_mm; d_mm = h_mm - eff_cover
    if d_mm<=0: raise ValueError("Depth too small")
    L_plus_d = short_span_m + d_mm/1000.0; L_c2c = short_span_m + 0.230; L_eff = min(L_plus_d, L_c2c)
    wu_kNpm = factored_load(dead_total, live_kNpm2)
    Mu = design_bending_moment(wu_kNpm, L_eff); Vu = wu_kNpm * L_eff / 2.0
    C_coeff = mu_limit_coefficient_for_fy(fy)
    Mu_Nmm = Mu*1e6; d_req_mm = math.sqrt(Mu_Nmm/(C_coeff*fck*1000.0)); d_req_mm = max(0.0,d_req_mm)
    if d_req_mm > d_mm:
        notes['depth_increase']=f"depth increased to {d_req_mm:.1f} mm"; d_mm=d_req_mm
        L_plus_d = short_span_m + d_mm/1000.0; L_eff = min(L_plus_d, L_c2c); Mu = design_bending_moment(wu_kNpm, L_eff); Vu = wu_kNpm * L_eff / 2.0
    Ast_req = annex_g_solve_Ast(Mu, fy, fck, 1000.0, d_mm)
    D_mm = h_mm; Ast_min = 0.0012*1000.0*D_mm
    if Ast_req < Ast_min: notes['as_min_used']=True; Ast_req = Ast_min
    if dia_bar_mm > D_mm/8.0: raise ValueError("Bar dia exceeds D/8 per IS 26.5.2.2")
    A_bar = math.pi*(dia_bar_mm/2.0)**2; S_calc = (1000.0*A_bar)/Ast_req
    lower_band = 0.85*S_calc; practical = list(range(int(math.ceil(lower_band/5.0)*5), int(math.floor(min(3.0*d_mm,300.0)/5.0)*5)+1,5))
    cand25=[s for s in practical if s%25==0 and lower_band<=s<=S_calc]
    if cand25: spacing_tension_mm=min(cand25)
    else:
        cand5=[s for s in practical if lower_band<=s<=S_calc]
        if cand5: spacing_tension_mm=min(cand5)
        else:
            spacing_tension_mm=min(practical, key=lambda s: abs(s-S_calc)); notes['spacing_warning']='used closest'
    Ast_prov=(1000.0*A_bar)/spacing_tension_mm
    dist_bar_mm=8.0; A_dist_bar=math.pi*(dist_bar_mm/2.0)**2; Ast_dist_req=0.0012*1000.0*h_mm
    S_calc_dist=(1000.0*A_dist_bar)/Ast_dist_req; lower_band_d=0.85*S_calc_dist
    practical_d=list(range(int(math.ceil(lower_band_d/5.0)*5), int(math.floor(min(5.0*d_mm,450.0)/5.0)*5)+1,5))
    cand25_d=[s for s in practical_d if s%25==0 and lower_band_d<=s<=S_calc_dist]
    if cand25_d: spacing_dist_mm=min(cand25_d)
    else:
        cand5_d=[s for s in practical_d if lower_band_d<=s<=S_calc_dist]
        if cand5_d: spacing_dist_mm=min(cand5_d)
        else: spacing_dist_mm=min(practical_d, key=lambda s: abs(s-S_calc_dist)); notes['dist_spacing_warning']='used closest'
    Ast_dist_prov=(1000.0*A_dist_bar)/spacing_dist_mm
    Vc_kN = shear_capacity_from_As(Ast_prov,1000.0,d_mm,fck); shear_ok = Vu <= Vc_kN
    L_over_d = (short_span_m*1000.0)/d_mm
    if L_over_d>20.0: notes['serviceability']=f"L/d = {L_over_d:.2f} > 20"
    res = OneWaySlabResult(short_span_m, long_span_m, h_mm, D_mm, d_mm, fck, fy, unit_w, dead_self, external_dead_kNpm2, dead_total, live_kNpm2, wu_kNpm, Mu, 0.0, d_req_mm if d_req_mm>0 else None, Ast_req, Ast_min, dia_bar_mm, spacing_tension_mm, Ast_prov, spacing_dist_mm, Ast_dist_prov, Vu, Vc_kN, shear_ok, L_over_d, notes)
    return res
