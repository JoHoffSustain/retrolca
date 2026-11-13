import hashlib
import json
import sqlite3
from dataclasses import asdict, dataclass
from os import PathLike
from typing import Any, Protocol, override

import cirpy


@dataclass
class NamingInfo:
    name: str
    smiles: str
    formula: str | None
    inchi: str | None
    inchi_key: str | None
    # CAS registry numbers; `None` means not looked up yet, an empty list
    # that the lookup did not return any CAS number
    cas: list[str] | None = None


class NamingService(Protocol):
    id: str

    def get_info(self, smiles: str) -> NamingInfo | None: ...

    def get_name(self, smiles: str) -> str:
        """Returns the name for the given SMILES code.

        It returns the SMILES code if the lookup did not return a result.
        """
        info = self.get_info(smiles)
        return info.name if info else smiles


class CIR(NamingService):
    def __init__(self):
        self.id = "CIR"

    @override
    def get_info(self, smiles: str) -> NamingInfo | None:
        mol = cirpy.Molecule(smiles)
        if not mol or not isinstance(mol.iupac_name, str):
            return None
        return NamingInfo(
            name=mol.iupac_name.lower() if mol.iupac_name else smiles,
            smiles=smiles,
            formula=_str(mol.formula),
            inchi=_str(mol.stdinchi),
            inchi_key=_str(mol.stdinchikey),
            cas=_cas_list(mol),
        )


def _cas_list(mol: cirpy.Molecule) -> list[str]:
    try:
        cas = mol.cas
    except Exception:
        return []
    if isinstance(cas, str):
        return [cas]
    if isinstance(cas, list):
        return [c for c in cas if isinstance(c, str)]
    return []


def _str(v: Any) -> str | None:
    if isinstance(v, str):
        return v
    else:
        return None


class CachingNamingService(NamingService):
    def __init__(self, path: str | PathLike, service: NamingService):
        self.service = service
        self.id = f"{service.id}/cached"
        self.conn = sqlite3.connect(path)
        self.conn.execute(
            """
            CREATE TABLE IF NOT EXISTS tbl_naming_infos (
                cache_key TEXT PRIMARY KEY,
                info_json TEXT NOT NULL
            )
            """
        )
        self.conn.commit()

    @override
    def get_info(self, smiles: str) -> NamingInfo | None:
        key = self.__key_of(smiles)
        cached = self.__lookup(key)
        if cached and cached.smiles != smiles:
            cached = None
        # entries cached before CAS numbers were collected are looked up again
        if cached and cached.cas is not None:
            return cached

        info = self.service.get_info(smiles)
        if not info:
            return cached

        self.__store(key, info)
        return info

    def close(self):
        self.conn.close()
        close = getattr(self.service, "close", None)
        if callable(close):
            close()

    def __key_of(self, smiles: str) -> str:
        raw = f"{self.service.id}||{smiles}".encode("utf-8")
        return hashlib.blake2b(raw, digest_size=32).hexdigest()

    def __lookup(self, key: str) -> NamingInfo | None:
        row = self.conn.execute(
            "SELECT info_json FROM tbl_naming_infos WHERE cache_key = ?",
            (key,),
        ).fetchone()
        if not row:
            return None

        try:
            data = json.loads(row[0])
        except (TypeError, json.JSONDecodeError):
            return None

        if not isinstance(data, dict):
            return None
        try:
            return NamingInfo(**data)
        except TypeError:
            return None

    def __store(self, key: str, info: NamingInfo):
        payload = json.dumps(asdict(info), separators=(",", ":"))
        self.conn.execute(
            """
            INSERT INTO tbl_naming_infos (cache_key, info_json)
            VALUES (?, ?)
            ON CONFLICT(cache_key) DO UPDATE SET info_json = excluded.info_json
            """,
            (key, payload),
        )
        self.conn.commit()
