import logging as log
import re
from dataclasses import dataclass
from typing import Protocol, cast, override

import olca_ipc as ipc
import olca_schema as o

from . import smiles
from .res import Res, nil, unwrap

# The official openLCA IDs of the quantity types 'mass' and 'chemical amount'.
# see https://github.com/GreenDelta/data/blob/master/refdata/flow_properties.csv
_MASS_ID = "93a60a56-a3c8-11da-a746-0800200b9a66"
_CHEMICAL_AMOUNT_ID = "341fd786-b2ad-4552-a762-5eafcab45dee"


class ProviderSelector(Protocol):
    """Protocol for provider selectors.

    When building process chains and for an input of a chemical a process
    producing that chemical already exists in the database, that process is
    linked as provider for the input. In case there are multiple processes
    producing that chemical, a rule needs to be defined how the linked provider
    is then selected. This is the role of a ProviderSelector then.
    """

    def select(self, providers: list[o.TechFlow]) -> o.TechFlow | None: ...


class DefaultProviderSelector(ProviderSelector):
    def __init__(self, preferred_location: str = "GLO"):
        self.preferred_location = preferred_location

    @override
    def select(self, providers: list[o.TechFlow]) -> o.TechFlow | None:
        if not providers:
            return None
        best: o.TechFlow | None = None
        best_score = 0.0
        for p in providers:
            score = self.__score_of(p)
            if not best or score > best_score:
                best = p
                best_score = score
        return best

    def __score_of(self, p: o.TechFlow) -> float:
        if not p or not p.flow or not p.provider:
            return 0.0
        score = 0.0
        if p.flow.name and p.provider.name:
            market_prefix = "market for " + p.flow.name
            prod_prefix = p.flow.name + " production"
            if p.provider.name.startswith(market_prefix):
                score += 8
            elif p.provider.name.startswith(prod_prefix):
                score += 5

        if p.provider.location:
            if p.provider.location == self.preferred_location:
                score += 10
            elif p.provider.location == "GLO":
                score += 5
            elif p.provider.location == "RoW":
                score += 2
        return score


@dataclass
class IpcContext:
    client: ipc.ProtoClient
    mass: o.FlowProperty
    chem_amount: o.FlowProperty
    kg: o.Unit
    mole: o.Unit
    provider_selector: ProviderSelector

    @staticmethod
    def of(
        client: ipc.ProtoClient,
        provider_selector: ProviderSelector | None = None,
    ) -> Res["IpcContext"]:
        mass = client.get(o.FlowProperty, _MASS_ID)
        if not mass:
            return nil, f"Flow property 'Mass' (id={_MASS_ID}) not found"
        if mass.name != "Mass":
            log.warning(
                "The name of the flow property '%s' (id=%s) should be 'Mass'",
                mass.name,
                _MASS_ID,
            )
        if not mass.unit_group or not mass.unit_group.id:
            return nil, "Flow property for 'Mass' does not link a unit group"

        mass_units = client.get(o.UnitGroup, mass.unit_group.id)
        if not mass_units or not mass_units.units:
            return nil, "No valid unit group for mass units found"
        kg = next((u for u in mass_units.units if u.name == "kg"), None)
        if not kg or not kg.is_ref_unit:
            return nil, "'kg' is not the reference unit of mass"

        chem_amount = client.get(o.FlowProperty, _CHEMICAL_AMOUNT_ID)
        if not chem_amount:
            return nil, (
                f"Flow property 'Chemical amount' (id={_CHEMICAL_AMOUNT_ID}) "
                f"not found"
            )
        if chem_amount.name != "Chemical amount":
            log.warning(
                "The name of the flow property '%s' (id=%s) should be "
                "'Chemical amount'",
                chem_amount.name,
                _CHEMICAL_AMOUNT_ID,
            )
        if not chem_amount.unit_group or not chem_amount.unit_group.id:
            return nil, (
                "Flow property for 'Chemical amount' does not link a unit group"
            )
        chem_amount_units = client.get(o.UnitGroup, chem_amount.unit_group.id)
        if not chem_amount_units or not chem_amount_units.units:
            return nil, "No valid unit group for chemical amount units found"
        mole = next(
            (u for u in chem_amount_units.units if u.name == "mol"), None
        )
        if not mole or not mole.is_ref_unit:
            return nil, "'mol' is not the reference unit of chemical amount"

        return (
            IpcContext(
                client=client,
                mass=mass,
                chem_amount=chem_amount,
                kg=kg,
                mole=mole,
                provider_selector=(
                    provider_selector
                    if provider_selector
                    else DefaultProviderSelector()
                ),
            ),
            nil,
        )

    def mass_prop_of(self, flow: o.Flow) -> o.FlowPropertyFactor | None:
        if not flow or not flow.flow_properties:
            return None
        for f in flow.flow_properties:
            if f.flow_property and f.flow_property.id == self.mass.id:
                return f
        return None

    def chem_prop_of(self, flow: o.Flow) -> o.FlowPropertyFactor | None:
        if not flow or not flow.flow_properties:
            return None
        for f in flow.flow_properties:
            if f.flow_property and f.flow_property.id == self.chem_amount.id:
                return f
        return None

    def molar_mass_of(self, flow: o.Flow | o.Ref) -> float | None:
        """Returns the molar mass in g/mol of the flow for this context.

        Returns ``None`` if the respective flow properties for mass and chemical
        amount are not defined in the flow.
        """
        if not flow:
            return None
        if isinstance(flow, o.Ref):
            f = self.client.get(o.Flow, flow.id)
            if not f:
                return None
            flow = f
        if not flow.flow_properties:
            return None

        mass_prop = self.mass_prop_of(flow)
        chem_prop = self.chem_prop_of(flow)
        if not mass_prop or not chem_prop:
            return None
        if mass_prop.is_ref_flow_property:
            return 1000 / cast(float, chem_prop.conversion_factor)
        if chem_prop.is_ref_flow_property:
            return 1000 * cast(float, mass_prop.conversion_factor)
        return None


@dataclass
class FlowIndex:
    data: dict[str, o.Flow]

    @staticmethod
    def of(ctx: IpcContext, flows: list[o.Flow] | None = None) -> "FlowIndex":
        log.info(
            "Build index of chemical products in openLCA with SMILES codes"
        )
        data = {}
        if flows is None:
            flows = ctx.client.get_all(o.Flow)
        for flow in flows:
            if flow.flow_type == o.FlowType.ELEMENTARY_FLOW:
                continue
            code = smiles.of_flow(flow)
            if not code:
                continue
            mm = ctx.molar_mass_of(flow)
            if not mm:
                continue
            data[code] = flow
        log.info("Created index with %d chemical products", len(data))
        return FlowIndex(data)


_CAS_PATTERN = re.compile(r"^0*(\d{2,7})-(\d{2})-(\d)$")


def normalize_cas(cas: str | None) -> str | None:
    """Normalizes a CAS registry number, e.g. `000811-97-2` to `811-97-2`.

    Leading zeros (as used in some databases) and surrounding white space are
    removed. Returns `None` when the value is not a valid CAS number, also
    when its check digit is wrong.
    """
    if not isinstance(cas, str):
        return None
    m = _CAS_PATTERN.match(cas.strip())
    if not m:
        return None
    a, b, check = m.groups()
    digits = (a + b)[::-1]
    total = sum((i + 1) * int(d) for i, d in enumerate(digits))
    if total % 10 != int(check):
        return None
    return f"{int(a)}-{b}-{check}"


def normalize_name(name: str | None) -> str | None:
    if not isinstance(name, str):
        return None
    n = " ".join(name.lower().split())
    return n if n else None


def normalize_formula(formula: str | None) -> str | None:
    if not isinstance(formula, str):
        return None
    f = "".join(formula.split())
    return f if f else None


class Compartment:
    AIR = "air"
    WATER = "water"


# the categories of the elementary flows to which direct emissions are
# mapped; only the unspecified sub-compartments are used
_COMPARTMENT_PATHS = {
    Compartment.AIR: "elementary flows/emission to air/unspecified",
    Compartment.WATER: "elementary flows/emission to water/unspecified",
}


def compartment_of(flow: o.Flow) -> str | None:
    if not flow or not flow.category:
        return None
    path = "/".join(p.strip() for p in flow.category.lower().split("/"))
    for comp, comp_path in _COMPARTMENT_PATHS.items():
        if path == comp_path:
            return comp
    return None


def compartment_path(comp: str) -> str:
    """Returns the category path of new elementary flows of a compartment."""
    if comp == Compartment.AIR:
        return "Elementary flows/Emission to air/unspecified"
    return "Elementary flows/Emission to water/unspecified"


@dataclass
class ElementaryFlowIndex:
    """An index of elementary flows for mapping chemicals to emissions.

    Elementary flows typically have no SMILES codes, but a CAS number, a
    formula, a name, and synonyms. Flows are matched by CAS number first. The
    name and synonyms are only used when the formula matches too. Elementary
    flows that have a SMILES code in their additional properties (e.g. flows
    generated by retrolca) are matched by that code.
    """

    by_cas: dict[tuple[str, str], o.Flow]
    by_name: dict[tuple[str, str], list[o.Flow]]
    by_smiles: dict[tuple[str, str], o.Flow]

    @staticmethod
    def of(
        ctx: IpcContext, flows: list[o.Flow] | None = None
    ) -> "ElementaryFlowIndex":
        idx = ElementaryFlowIndex({}, {}, {})
        if flows is None:
            flows = ctx.client.get_all(o.Flow)
        for flow in flows:
            if flow.flow_type != o.FlowType.ELEMENTARY_FLOW:
                continue
            comp = compartment_of(flow)
            if not comp:
                continue
            mass_prop = ctx.mass_prop_of(flow)
            if not mass_prop or not mass_prop.is_ref_flow_property:
                continue
            idx.put(flow, comp)
        log.info(
            "Indexed %d elementary flows by CAS number for emissions",
            len(idx.by_cas),
        )
        return idx

    def put(self, flow: o.Flow, comp: str):
        cas = normalize_cas(flow.cas)
        if cas and (cas, comp) not in self.by_cas:
            self.by_cas[(cas, comp)] = flow
        code = smiles.of_flow(flow)
        if code:
            self.by_smiles[(code, comp)] = flow
        names = [flow.name]
        if flow.synonyms:
            names.extend(flow.synonyms.split(";"))
        for name in names:
            n = normalize_name(name)
            if not n:
                continue
            self.by_name.setdefault((n, comp), []).append(flow)

    def find(
        self,
        comp: str,
        smiles_code: str | None = None,
        cas_numbers: list[str] | None = None,
        names: list[str] | None = None,
        formula: str | None = None,
    ) -> o.Flow | None:
        flow, _ = self.find_match(comp, smiles_code, cas_numbers, names, formula)
        return flow

    def find_name(self, comp: str, names: list[str]) -> o.Flow | None:
        """Finds a flow by one of the given names or synonyms (in order)."""
        for name in names:
            n = normalize_name(name)
            if n and (flows := self.by_name.get((n, comp))):
                return flows[0]
        return None

    def find_match(
        self,
        comp: str,
        smiles_code: str | None = None,
        cas_numbers: list[str] | None = None,
        names: list[str] | None = None,
        formula: str | None = None,
    ) -> tuple[o.Flow | None, str | None]:
        """Like `find`, but also returns how the flow was matched."""
        if smiles_code and (f := self.by_smiles.get((smiles_code, comp))):
            return f, "SMILES"
        for cas in cas_numbers or []:
            c = normalize_cas(cas)
            if c and (f := self.by_cas.get((c, comp))):
                return f, f"CAS {c}"
        formula = normalize_formula(formula)
        if not formula:
            return None, None
        for name in names or []:
            n = normalize_name(name)
            if not n:
                continue
            for f in self.by_name.get((n, comp), []):
                if normalize_formula(f.formula) == formula:
                    return f, f"name '{name}' + formula {formula}"
        return None, None


@dataclass
class ProviderIndex:
    data: dict[str, o.TechFlow]

    @classmethod
    def of(
        cls, ctx: IpcContext, flows: list[o.Flow] | None = None
    ) -> "ProviderIndex":
        data = {}
        provider_idx = cls.__index_providers(ctx.client)
        if flows is None:
            flows = ctx.client.get_all(o.Flow)
        for flow in flows:
            smiles_code = smiles.of_flow(flow)
            if not smiles_code:
                continue
            mm = ctx.molar_mass_of(flow)
            if not mm:
                continue
            providers = provider_idx.get(unwrap(flow.id))
            if not providers:
                continue
            provider = ctx.provider_selector.select(providers)
            if provider:
                data[smiles_code] = provider
        return cls(data)

    @staticmethod
    def __index_providers(
        client: ipc.ProtoClient,
    ) -> dict[str, list[o.TechFlow]]:
        provider_idx = {}
        for p in client.get_providers():
            if not p.flow or not p.flow.id:
                continue
            ps = provider_idx.get(p.flow.id)
            if not ps:
                ps = []
                provider_idx[p.flow.id] = ps
            ps.append(p)
        return provider_idx

    def get(self, smiles_code: str) -> o.TechFlow | None:
        return self.data.get(smiles_code)

    def put(self, smiles_code: str, process: o.Process) -> o.TechFlow:
        flow_ref: o.Ref | None = None
        for e in unwrap(process.exchanges):
            if e.is_quantitative_reference:
                flow_ref = e.flow
                break
        provider = o.TechFlow(flow=flow_ref, provider=process.to_ref())
        self.data[smiles_code] = provider
        return provider
