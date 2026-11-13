"""Boiling points of chemicals from PubChem.

The experimental boiling points that PubChem lists for a compound come from
different sources and in different formats, e.g.:

    243 to 244 °F at 760 mmHg (NTP, 1992)
    115.7 °C
    116.50 °C. @ 760.00 mm Hg
    117-118 °C
    116.5 °C @760 [mm Hg]

Only values at normal pressure (760 mmHg, 1 atm, 101.3 kPa) or without a
stated pressure are used; ranges count with their mean; °F and K are
converted to °C. The boiling point is the median of these values.

The result is stored in the additional properties of the flow, so that it can
be checked and corrected in openLCA. A value that is already stored there is
never overwritten, so a value entered by hand always wins.
"""

import datetime
import logging
import re
import statistics
from dataclasses import dataclass

import olca_schema as o
import requests

log = logging.getLogger(__name__)

# keys of the additional properties of a flow
BP_KEY = "Boiling-Point"  # in °C
BP_SOURCE_KEY = "Boiling-Point-Source"
BP_VALUES_KEY = "Boiling-Point-Values"

_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_PUG_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

_NUM = r"[-−]?\d+(?:[.,]\d+)?"
_TEMP = re.compile(
    rf"({_NUM})\s*(?:(?:to|-|–|—)\s*({_NUM})\s*)?"
    r"(?:°|º|deg\.?|degrees?)?\s*(C|F|K)\b",
    re.IGNORECASE,
)
_PRESSURE = re.compile(
    rf"({_NUM})\s*\[?\s*(mm\s*hg|torr|kpa|hpa|mbar|bar|atm|pa)\b",
    re.IGNORECASE,
)

# normal pressure in the respective units
_NORMAL = {
    "mmhg": 760.0,
    "torr": 760.0,
    "kpa": 101.325,
    "hpa": 1013.25,
    "mbar": 1013.25,
    "bar": 1.01325,
    "atm": 1.0,
    "pa": 101325.0,
}


@dataclass
class BoilingPoint:
    value: float  # °C, median of the used values
    values: list[float]
    source: str


def _float(s: str) -> float:
    return float(s.replace("−", "-").replace(",", "."))


def _to_celsius(v: float, unit: str) -> float:
    u = unit.upper()
    if u == "F":
        return (v - 32) * 5 / 9
    if u == "K":
        return v - 273.15
    return v


def _is_normal_pressure(text: str) -> bool:
    """True if the text states no pressure or normal pressure."""
    pressures = _PRESSURE.findall(text)
    if not pressures:
        return True
    for num, unit in pressures:
        normal = _NORMAL.get(re.sub(r"\s+", "", unit.lower()))
        if normal is None:
            return False
        if abs(_float(num) - normal) > 0.02 * normal:
            return False
    return True


def parse_values(text: str) -> list[float]:
    """Extracts the boiling points at normal pressure in °C from a text.

    A text can contain several statements separated by `;`, e.g.
    "56.05 °C at 760 mm Hg; 39.5 °C at 400 mm Hg".
    """
    if not text:
        return []
    values: list[float] = []
    for part in re.split(r";|\n", text):
        low = part.lower()
        if "decompos" in low or "sublim" in low:
            continue
        if not _is_normal_pressure(part):
            continue
        m = _TEMP.search(part)
        if not m:
            continue
        a, b, unit = m.groups()
        v = _float(a)
        if b:
            v = (v + _float(b)) / 2
        values.append(round(_to_celsius(v, unit), 3))
    return values


def _strings_of(node) -> list[str]:
    """Collects the boiling point statements of a PUG View record."""
    out: list[str] = []
    if isinstance(node, dict):
        for info in node.get("Information", []) or []:
            value = info.get("Value", {}) if isinstance(info, dict) else {}
            for s in value.get("StringWithMarkup", []) or []:
                if isinstance(s, dict) and isinstance(s.get("String"), str):
                    out.append(s["String"])
            nums = value.get("Number")
            unit = value.get("Unit")
            if isinstance(nums, list) and isinstance(unit, str):
                for n in nums:
                    out.append(f"{n} {unit}")
        for child in node.get("Section", []) or []:
            out.extend(_strings_of(child))
        if "Record" in node:
            out.extend(_strings_of(node["Record"]))
    return out


def cid_of_smiles(smiles_code: str, timeout=30.0) -> int | None:
    try:
        resp = requests.post(
            f"{_PUG}/compound/smiles/cids/JSON",
            data={"smiles": smiles_code},
            timeout=timeout,
        )
        if resp.status_code != 200:
            return None
        cids = resp.json().get("IdentifierList", {}).get("CID", [])
    except Exception as err:
        log.warning("PubChem CID lookup failed for %s: %s", smiles_code, err)
        return None
    cid = cids[0] if cids else None
    return cid if isinstance(cid, int) and cid > 0 else None


def from_pubchem(smiles_code: str, timeout=30.0) -> BoilingPoint | None:
    """Returns the boiling point of a compound from PubChem, if available."""
    cid = cid_of_smiles(smiles_code, timeout)
    if not cid:
        return None
    try:
        resp = requests.get(
            f"{_PUG_VIEW}/data/compound/{cid}/JSON",
            params={"heading": "Boiling Point"},
            timeout=timeout,
        )
        if resp.status_code != 200:
            return None
        data = resp.json()
    except Exception as err:
        log.warning("PubChem boiling point lookup failed for %s: %s", cid, err)
        return None
    return from_statements(_strings_of(data), f"PubChem CID {cid}")


def from_statements(statements: list[str], origin: str) -> BoilingPoint | None:
    values: list[float] = []
    for s in statements:
        values.extend(parse_values(s))
    if not values:
        return None
    return BoilingPoint(
        value=round(statistics.median(values), 2),
        values=values,
        source=(
            f"{origin}: median of {len(values)} value(s) at 760 mmHg "
            f"or without stated pressure"
        ),
    )


def stored_value(flow: o.Flow | None) -> float | None:
    """Returns the boiling point in °C stored in the flow, if any."""
    if not flow or not flow.other_properties:
        return None
    v = flow.other_properties.get(BP_KEY)
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        try:
            return _float(v)
        except ValueError:
            return None
    return None


def was_checked(flow: o.Flow | None) -> bool:
    return bool(
        flow and flow.other_properties and BP_SOURCE_KEY in flow.other_properties
    )


def store(flow: o.Flow, bp: BoilingPoint | None):
    """Writes the boiling point into the additional properties of the flow.

    If no boiling point was found, only the source is written, so that the
    lookup is not repeated and the value can be entered by hand.
    """
    props = flow.other_properties or {}
    flow.other_properties = props
    if bp:
        props[BP_KEY] = bp.value
        props[BP_VALUES_KEY] = bp.values
        props[BP_SOURCE_KEY] = bp.source
    else:
        today = datetime.date.today().isoformat()
        props[BP_SOURCE_KEY] = (
            f"not found on PubChem ({today}); enter the value in °C under "
            f"'{BP_KEY}' if known"
        )
