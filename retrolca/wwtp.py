"""Yield losses with an internal wastewater treatment plant (WWTP).

This loss model follows the stoichiometric estimation approach for chemicals
described by Hischier et al. (2005) and Althaus et al. (2007). The reactant
inputs are the stoichiometric amounts divided by the yield; of each reactant
input, the air share is emitted to air as the substance itself, and the rest
of the unreacted amount goes to an internal WWTP:

    input = s / retro_yield
    air   = s / retro_yield * retro_air_share
    WWTP  = s / retro_yield * (1 - retro_yield - retro_air_share)

In the WWTP, depending on the type of substance (fragment of the SMILES code):

* organic: the share retro_wwtp_removal is degraded; its carbon is emitted as
  CO2 (fossil) to air, its halogens, sulfur and phosphorus as halide,
  sulfate and phosphate ions to water, and its nitrogen follows the nitrogen
  rule below. The rest is emitted to water as the substance itself; TOC and
  DOC are its carbon, COD is retro_cod_factor * ThOD, and BOD5 is
  retro_bod_factor * COD.
* inorganic acids (HCl, HBr, HF, HI, H2SO4, HNO3, H3PO4): neutralized with
  sodium hydroxide (an additional input); the anion and sodium go to water.
* ions of salts (e.g. Na+, Cl-, SO4 2-): to water as ions.
* nitrogen compounds without carbon (e.g. ammonia) and the nitrogen of the
  degraded organic part: of the nitrogen, retro_nh4_removal is removed as
  ammonium-N and retro_ntot_removal as total N; the rest is emitted as
  ammonium (1 - retro_nh4_removal) and nitrate (retro_nh4_removal -
  retro_ntot_removal). The removed nitrogen (N2) is not reported.
* boiling point: for organic substances and neutral inorganic substances
  (not acids, ammonia or salts, which are treated as dissolved), each
  generated process gets a local parameter bp_r<n> with the boiling point
  in °C of the reactant (from the flow property Boiling-Point, or a built-in
  value for permanent gases; 9999 if unknown). If it is below retro_bp_gas
  (default -100 °C, i.e. permanent gases like hydrogen), the share
  retro_gas_air_share (default 1) of the unreacted amount is emitted to air
  and only the rest goes to the WWTP. Otherwise, and if the boiling point is
  unknown, the default distribution above is used. The amounts are
  if-formulas of the parameter, so a boiling point entered later in openLCA
  changes them.
* water: nothing.
* other inorganic substances: to water as the substance itself.

All shares are global openLCA parameters and the exchange amounts are
formulas of them.
"""

import logging
from dataclasses import dataclass, field

import olca_schema as o
from rdkit import Chem
from rdkit.Chem import Descriptors

from . import oipc
from .losses import Reactant, YieldModel, _num, formula_of
from .res import unwrap

log = logging.getLogger(__name__)

Y = "retro_yield"
AIR = "retro_air_share"
REM = "retro_wwtp_removal"
NH4_REM = "retro_nh4_removal"
NTOT_REM = "retro_ntot_removal"
COD_F = "retro_cod_factor"
BOD_F = "retro_bod_factor"
BP_GAS = "retro_bp_gas"
GAS_AIR = "retro_gas_air_share"

# boiling point parameter value for an unknown boiling point
BP_UNKNOWN = 9999.0

M_C, M_N, M_S, M_P = 12.011, 14.007, 32.06, 30.974
M_CO2, M_NH4, M_NO3 = 44.009, 18.038, 62.004
M_SO4, M_PO4, M_NA, M_NAOH = 96.06, 94.971, 22.990, 39.997
HALIDES = {"F": 18.998, "Cl": 35.453, "Br": 79.904, "I": 126.904}

# the WWTP inflow of a reactant as formula; `s` is replaced by the amount
_W = f"{{s}} / {Y} * (1 - {Y} - {AIR})"

# names of the elementary flows that are matched by name; the first name
# found in the database is used
FLOW_NAMES: dict[str, tuple[str, list[str]]] = {
    "CO2": (oipc.Compartment.AIR, ["Carbon dioxide, fossil"]),
    "TOC": (oipc.Compartment.WATER, ["TOC, Total Organic Carbon"]),
    "DOC": (oipc.Compartment.WATER, ["DOC, Dissolved Organic Carbon"]),
    "COD": (oipc.Compartment.WATER, ["COD, Chemical Oxygen Demand"]),
    "BOD5": (oipc.Compartment.WATER, ["BOD5, Biological Oxygen Demand"]),
    "NH4": (oipc.Compartment.WATER, ["Ammonium, ion", "Ammonium"]),
    "NO3": (oipc.Compartment.WATER, ["Nitrate"]),
    "Na": (oipc.Compartment.WATER, ["Sodium I", "Sodium, ion", "Sodium"]),
    "K": (oipc.Compartment.WATER, ["Potassium I", "Potassium, ion", "Potassium"]),
    "F": (oipc.Compartment.WATER, ["Fluoride"]),
    "Cl": (oipc.Compartment.WATER, ["Chloride"]),
    "Br": (oipc.Compartment.WATER, ["Bromide"]),
    "I": (oipc.Compartment.WATER, ["Iodide"]),
    "SO4": (oipc.Compartment.WATER, ["Sulfate"]),
    "PO4": (oipc.Compartment.WATER, ["Phosphate"]),
}

# acids: canonical SMILES -> (number of acidic H, anion key, anion mass)
ACIDS = {
    "Cl": (1, "Cl", HALIDES["Cl"]),
    "Br": (1, "Br", HALIDES["Br"]),
    "F": (1, "F", HALIDES["F"]),
    "I": (1, "I", HALIDES["I"]),
    "O=S(=O)(O)O": (2, "SO4", M_SO4),
    "O=[N+]([O-])O": (1, "NO3", M_NO3),
    "O=P(O)(O)O": (3, "PO4", M_PO4),
}

# ions: canonical SMILES of the fragment -> (flow key, ion mass or None if
# the nitrogen rule applies)
IONS = {
    "[Na+]": ("Na", M_NA),
    "[K+]": ("K", 39.098),
    "[F-]": ("F", HALIDES["F"]),
    "[Cl-]": ("Cl", HALIDES["Cl"]),
    "[Br-]": ("Br", HALIDES["Br"]),
    "[I-]": ("I", HALIDES["I"]),
    "O=S(=O)([O-])[O-]": ("SO4", M_SO4),
    "O=[N+]([O-])[O-]": ("NO3", M_NO3),
    "O=P([O-])([O-])[O-]": ("PO4", M_PO4),
}

# boiling points (°C) of common permanent gases, used when no boiling point
# is stored in the flow
GASES = {
    "[H][H]": -252.9,
    "N#N": -195.8,
    "O=O": -183.0,
    "[C-]#[O+]": -191.5,
    "C": -161.5,
    "[He]": -268.9,
    "[Ne]": -246.0,
    "[Ar]": -185.8,
}
NOTHING = {"O", "[OH-]"}  # water, hydroxide


def _canon_keys(d):
    out = {}
    for k, v in d.items():
        m = Chem.MolFromSmiles(k)  # ty: ignore
        out[Chem.MolToSmiles(m) if m else k] = v  # ty: ignore
    return out


ACIDS = _canon_keys(ACIDS)
IONS = _canon_keys(IONS)
GASES = _canon_keys(GASES)
NOTHING = set(_canon_keys({n: None for n in NOTHING}))


def _canon(mol) -> str:
    return Chem.MolToSmiles(mol)  # ty: ignore


def _composition(mol) -> tuple[dict[str, int], float]:
    m = Chem.AddHs(mol)  # ty: ignore
    counts: dict[str, int] = {}
    for a in m.GetAtoms():
        counts[a.GetSymbol()] = counts.get(a.GetSymbol(), 0) + 1
    return counts, Descriptors.MolWt(mol)


def _is_dissolved_nitrogen(canon: str, counts: dict[str, int]) -> bool:
    return canon == "[NH4+]" or (counts.get("C", 0) == 0 and "N" in counts)


def _can_evaporate(mol) -> bool:
    """True for organic substances and neutral inorganic substances that are
    not ions, ammonium or nitrogen compounds without carbon."""
    for frag in Chem.GetMolFrags(mol, asMols=True):  # ty: ignore
        fc = _canon(frag)
        if fc in IONS or fc in NOTHING:
            continue
        counts, _ = _composition(frag)
        if _is_dissolved_nitrogen(fc, counts):
            continue
        if any(a.GetFormalCharge() != 0 for a in frag.GetAtoms()):
            if counts.get("C", 0) == 0:
                continue
        return True
    return False


def thod_of(counts: dict[str, int], mw: float) -> float:
    """Theoretical oxygen demand in g O2 / g (nitrogen ends as NH3)."""
    g = counts.get
    hal = sum(g(x, 0) for x in HALIDES)
    o2 = (
        2 * g("C", 0)
        + 0.5 * (g("H", 0) - hal - 3 * g("N", 0))
        + 3 * g("S", 0)
        + 2.5 * g("P", 0)
        + 0.5 * g("Na", 0)
        - g("O", 0)
    )
    return max(16 * o2 / mw, 0.0)


@dataclass
class _Terms:
    """Formula terms of the summed exchanges of a process."""

    values: dict[str, float] = field(default_factory=dict)
    formulas: dict[str, list[str]] = field(default_factory=dict)

    def add(self, key: str, value: float, formula: str):
        if value <= 0:
            return
        self.values[key] = self.values.get(key, 0.0) + value
        self.formulas.setdefault(key, []).append(formula)


class WwtpModel(YieldModel):
    def __init__(
        self,
        ctx: oipc.IpcContext,
        yield_: float,
        naming,
        flows: list[o.Flow] | None = None,
        naoh_provider: o.TechFlow | None = None,
        air_share: float = 0.002,
        wwtp_removal: float = 0.9,
        nh4_removal: float = 0.7,
        ntot_removal: float = 0.5,
        cod_factor: float = 0.96,
        bod_factor: float = 1.0,
        bp_gas: float = -100.0,
        gas_air_share: float = 1.0,
        pubchem_lookup: bool = True,
    ):
        super().__init__(
            ctx, yield_, 0.0, None, naming, flows, pubchem_lookup
        )
        self.naoh_provider = naoh_provider
        p = self._global_param
        self.air_share = p(
            AIR, air_share,
            "Share of the reactant input emitted to air as the substance "
            "itself in the processes generated by retrolca (purge vents, "
            "fugitive emissions).",
        )
        self.removal = p(
            REM, wwtp_removal,
            "Share of the organic substances degraded in the internal "
            "wastewater treatment; their carbon is emitted as CO2.",
        )
        self.nh4_removal = p(
            NH4_REM, nh4_removal,
            "Removal efficiency for ammonium-N in the internal wastewater "
            "treatment.",
        )
        self.ntot_removal = p(
            NTOT_REM, ntot_removal,
            "Removal efficiency for total N in the internal wastewater "
            "treatment.",
        )
        self.cod_factor = p(
            COD_F, cod_factor,
            "COD = this factor * theoretical oxygen demand (ThOD) of the "
            "organic substances emitted to water.",
        )
        self.bod_factor = p(
            BOD_F, bod_factor,
            "BOD5 = this factor * COD (1 = worst case BOD5 = COD).",
        )
        self.bp_gas = p(
            BP_GAS, bp_gas,
            "Boiling point threshold in °C: unreacted substances with a "
            "lower boiling point (permanent gases) are emitted to air with "
            "the share retro_gas_air_share instead of going to the "
            "internal wastewater treatment.",
        )
        self.gas_air_share = p(
            GAS_AIR, gas_air_share,
            "Share of the unreacted amount of a permanent gas (boiling "
            "point below retro_bp_gas) emitted to air.",
        )
        self.unknown_bp: set[str] = set()
        if self.air_share > 1 - self.yield_value:
            log.warning(
                "The air share (%s) is larger than the unreacted share "
                "(1 - yield = %s)",
                self.air_share,
                1 - self.yield_value,
            )
        self.flow_cache: dict[str, o.Flow | None] = {}
        self.missing_names: set[str] = set()
        self.__warned_naoh = False

    def _uses_emission_share(self) -> bool:
        return False

    def _named_flow(self, key: str) -> o.Flow | None:
        if key in self.flow_cache:
            return self.flow_cache[key]
        comp, names = FLOW_NAMES[key]
        flow = self.index.find_name(comp, names)
        if not flow:
            log.error(
                "No elementary flow found for %s (emission to %s, "
                "unspecified; names tried: %s); this emission is not "
                "modelled",
                key,
                comp,
                ", ".join(names),
            )
            self.missing_names.add(names[0])
        self.flow_cache[key] = flow
        return flow

    def _substance_flow(self, r: Reactant, comp: str) -> o.Flow:
        cas = self._cas_of(r)
        names = self._names_of(r)
        formula = formula_of(r.smiles)
        flow, how = self.index.find_match(
            comp, smiles_code=r.smiles, cas_numbers=cas, names=names,
            formula=formula,
        )
        if flow:
            key = f"{r.smiles}|{comp}"
            if key not in self.matched:
                self.matched[key] = r.flow.name or r.smiles
                log.info(
                    "Emissions of '%s' [%s, CAS %s] to %s mapped to '%s' "
                    "(by %s)",
                    r.flow.name, r.smiles, ", ".join(cas) if cas else "-",
                    comp, flow.name, how,
                )
            return flow
        return self._create_flow(r, cas, formula, comp)

    # --- losses

    def add_losses(self, process: o.Process, reactants: list[Reactant]):
        if not reactants or self.yield_value >= 1:
            return
        terms = _Terms()
        for r in reactants:
            self.__reactant(process, r, terms)
        self.__add_terms(process, terms)

    def __output(self, process, flow, amount, formula, description):
        e = o.new_output(process, flow, amount, self.ctx.kg)
        e.flow_property = self.ctx.mass.to_ref()
        e.amount_formula = formula
        e.description = description
        return e

    def __reactant(self, process: o.Process, r: Reactant, terms: _Terms):
        s_ = r.stoich_amount
        sn = _num(s_)
        name = r.flow.name or r.smiles
        y, a = self.yield_value, self.air_share

        mol = Chem.MolFromSmiles(r.smiles)  # ty: ignore
        if not mol:
            log.warning("Cannot parse %s; no losses modelled", r.smiles)
            return
        canon = _canon(mol)
        if canon in NOTHING:
            return

        # air share and WWTP inflow; for substances that can evaporate, the
        # split depends on the boiling point parameter of the process
        air = s_ / y * a
        air_f = f"{sn} / {Y} * {AIR}"
        w = s_ / y * (1 - y - a)
        w_f = _W.format(s=sn)
        if canon not in ACIDS and _can_evaporate(mol):
            pn, bp = self.__bp_param(process, r, canon)
            g = self.gas_air_share
            if bp < self.bp_gas:
                air, w = s_ / y * (1 - y) * g, s_ / y * (1 - y) * (1 - g)
            cond = f"{pn} < {BP_GAS}"
            air_f = f"{sn} / {Y} * if({cond}; (1 - {Y}) * {GAS_AIR}; {AIR})"
            w_f = (
                f"{sn} / {Y} * if({cond}; (1 - {Y}) * (1 - {GAS_AIR}); "
                f"1 - {Y} - {AIR})"
            )

        # air: the substance itself
        if air > 0:
            flow = self._substance_flow(r, oipc.Compartment.AIR)
            self.__output(
                process, flow, air, air_f,
                f"Unreacted {name}, emission to air (purge vents, "
                f"fugitive emissions, gases)",
            )
        if w <= 0:
            return

        if canon in ACIDS:
            self.__acid(process, r, canon, w, w_f, terms)
            return

        _, mw_total = _composition(mol)
        for frag in Chem.GetMolFrags(mol, asMols=True):  # ty: ignore
            fc = _canon(frag)
            counts, mw = _composition(frag)
            share = mw / mw_total
            if abs(share - 1) < 1e-12:
                fw, fw_f = w, w_f
            else:
                fw, fw_f = w * share, f"{w_f} * {_num(share)}"
            if fc in NOTHING:
                continue
            if fc in IONS:
                key, m_ion = IONS[fc]
                terms.add(key, fw, fw_f)
                continue
            if _is_dissolved_nitrogen(fc, counts):
                n_mass = fw * counts.get("N", 0) * M_N / mw
                n_f = f"{fw_f} * {_num(counts.get('N', 0) * M_N / mw)}"
                self.__nitrogen(n_mass, n_f, terms)
                continue
            if counts.get("C", 0) == 0:
                flow = self._substance_flow(r, oipc.Compartment.WATER)
                self.__output(
                    process, flow, fw, fw_f,
                    f"Unreacted {name} (inorganic), emission to water",
                )
                continue
            self.__organic(process, r, counts, mw, fw, fw_f, terms)

    def __bp_param(
        self, process: o.Process, r: Reactant, canon: str
    ) -> tuple[str, float]:
        """Adds the boiling point parameter of a reactant to the process."""
        from . import boiling

        name = r.flow.name or r.smiles
        bp = boiling.stored_value(r.product_flow)
        if bp is not None:
            source = "flow property Boiling-Point"
        elif canon in GASES:
            bp, source = GASES[canon], "built-in value for permanent gases"
        else:
            bp, source = BP_UNKNOWN, "unknown"
            if r.smiles not in self.unknown_bp:
                self.unknown_bp.add(r.smiles)
                log.info(
                    "No boiling point for '%s'; the default distribution is "
                    "used (process parameter value %s)",
                    name,
                    BP_UNKNOWN,
                )
        process.parameters = process.parameters or []
        pn = f"bp_r{len(process.parameters) + 1}"
        param = o.new_parameter(pn, bp, o.ParameterScope.PROCESS_SCOPE)
        param.is_input_parameter = True
        param.description = (
            f"Boiling point in °C of {name} [{r.smiles}], source: {source}. "
            f"{BP_UNKNOWN:g} = unknown, then the default distribution is "
            f"used. Below {BP_GAS} the unreacted amount is emitted to air."
        )
        process.parameters.append(param)
        return pn, bp

    def __acid(self, process, r, canon, w, w_f, terms: _Terms):
        n_h, key, m_anion = ACIDS[canon]
        mol = Chem.MolFromSmiles(canon)  # ty: ignore
        mw = Descriptors.MolWt(mol)
        terms.add(key, w * m_anion / mw, f"{w_f} * {_num(m_anion / mw)}")
        naoh = w * n_h * M_NAOH / mw
        naoh_f = f"{w_f} * {_num(n_h * M_NAOH / mw)}"
        p = self.naoh_provider
        if not p or not p.flow or not p.provider:
            if not self.__warned_naoh:
                log.warning(
                    "No sodium hydroxide provider found; the neutralization "
                    "of unreacted acids is not modelled"
                )
                self.__warned_naoh = True
            return
        e = o.new_input(process, unwrap(p.flow), naoh, self.ctx.kg)
        e.flow_property = self.ctx.mass.to_ref()
        e.default_provider = p.provider
        e.amount_formula = naoh_f
        e.description = (
            f"Neutralization of unreacted {r.flow.name or r.smiles}"
        )
        terms.add("Na", w * n_h * M_NA / mw, f"{w_f} * {_num(n_h * M_NA / mw)}")

    def __nitrogen(self, n_mass: float, n_f: str, terms: _Terms):
        terms.add(
            "NH4",
            n_mass * (1 - self.nh4_removal) * M_NH4 / M_N,
            f"{n_f} * (1 - {NH4_REM}) * {_num(M_NH4 / M_N)}",
        )
        terms.add(
            "NO3",
            n_mass * (self.nh4_removal - self.ntot_removal) * M_NO3 / M_N,
            f"{n_f} * ({NH4_REM} - {NTOT_REM}) * {_num(M_NO3 / M_N)}",
        )

    def __organic(self, process, r, counts, mw, fw, fw_f, terms: _Terms):
        name = r.flow.name or r.smiles
        rem = self.removal
        deg, deg_f = fw * rem, f"{fw_f} * {REM}"
        rest, rest_f = fw * (1 - rem), f"{fw_f} * (1 - {REM})"
        w_c = counts.get("C", 0) * M_C / mw

        # degraded part
        terms.add("CO2", deg * w_c * M_CO2 / M_C,
                  f"{deg_f} * {_num(w_c * M_CO2 / M_C)}")
        for x, m_x in HALIDES.items():
            if counts.get(x):
                f = counts[x] * m_x / mw
                terms.add(x, deg * f, f"{deg_f} * {_num(f)}")
        if counts.get("S"):
            f = counts["S"] * M_SO4 / mw
            terms.add("SO4", deg * f, f"{deg_f} * {_num(f)}")
        if counts.get("P"):
            f = counts["P"] * M_PO4 / mw
            terms.add("PO4", deg * f, f"{deg_f} * {_num(f)}")
        if counts.get("N"):
            f = counts["N"] * M_N / mw
            self.__nitrogen(deg * f, f"{deg_f} * {_num(f)}", terms)

        # non-degraded part: the substance itself and summary parameters
        if rest <= 0:
            return
        flow = self._substance_flow(r, oipc.Compartment.WATER)
        self.__output(
            process, flow, rest, rest_f,
            f"Unreacted {name}, not degraded in the internal wastewater "
            f"treatment, emission to water",
        )
        terms.add("TOC", rest * w_c, f"{rest_f} * {_num(w_c)}")
        terms.add("DOC", rest * w_c, f"{rest_f} * {_num(w_c)}")
        thod = thod_of(counts, mw)
        cod = rest * thod * self.cod_factor
        cod_f = f"{rest_f} * {_num(thod)} * {COD_F}"
        terms.add("COD", cod, cod_f)
        terms.add("BOD5", cod * self.bod_factor, f"{cod_f} * {BOD_F}")

    def __add_terms(self, process: o.Process, terms: _Terms):
        descriptions = {
            "CO2": "Carbon of the unreacted reactants degraded in the "
            "internal wastewater treatment",
            "NH4": "Nitrogen of the unreacted reactants after the internal "
            "wastewater treatment",
            "NO3": "Nitrogen of the unreacted reactants after the internal "
            "wastewater treatment",
        }
        for key, value in terms.values.items():
            flow = self._named_flow(key)
            if not flow:
                continue
            self.__output(
                process, flow, value, " + ".join(terms.formulas[key]),
                descriptions.get(
                    key,
                    "Unreacted reactants after the internal wastewater "
                    "treatment",
                ),
            )

    def log_summary(self) -> None:
        super().log_summary()
        if self.unknown_bp:
            log.info(
                "%d substance(s) without boiling point; enter it in the "
                "process parameters bp_r<n> (°C) to apply the boiling point "
                "rule",
                len(self.unknown_bp),
            )
        if self.missing_names:
            log.error(
                "Elementary flows not found (emissions not modelled): %s",
                ", ".join(sorted(self.missing_names)),
            )
