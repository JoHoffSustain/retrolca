import json
import logging
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import override

import aizynthfinder.aizynthfinder as ai
import requests

from .tool import Reaction, RetroTool
from .res import Res, nil

log = logging.getLogger(__name__)


class ZynthTool(RetroTool):
    def __init__(
        self,
        configfile: str | Path,
        policies: str | list[str] = "uspto",
        return_n: int = 10,
    ):
        """Wraps the AiZynthFinder one-step expansion.

        Args:
            configfile:
                The path to the AiZynthFinder config.yml.
            policies:
                The expansion policy or policies to use, as named in the
                config, e.g. "uspto", "ringbreaker", or
                ["ringbreaker", "uspto"] to combine both. RingBreaker is
                trained on ring-forming reactions (e.g. heterocycle
                syntheses) that the USPTO templates often miss.
            return_n:
                The maximum number of distinct reactant sets returned per
                expansion. The results of the policies are concatenated in
                the given order before this cut, so with several policies
                this should be large enough to keep results of all of them.
        """
        log.info("Loading AiZynthExpander with config: %s", configfile)
        if isinstance(policies, str):
            policies = [policies]
        self.policies = policies
        self.return_n = return_n
        # the policies are part of the ID, as it is used as cache key
        self.id = "zynth" if policies == ["uspto"] else (
            "zynth-" + "+".join(policies)
        )
        self.expander = ai.AiZynthExpander(configfile=str(configfile))
        self.expander.expansion_policy.select(policies)
        self.expander.filter_policy.select("uspto")
        log.info("Selected expansion policies: %s", ", ".join(policies))

    @override
    def expand(self, smiles: str) -> Res[list[Reaction]]:
        log.info("Expand SMILES code: %s", smiles)
        reactions: list[Reaction] = []
        expansion = self.expander.do_expansion(smiles, return_n=self.return_n)
        for group in expansion:
            reaction: Reaction | None = None
            for candidate in group:
                # for the same reactants there could be multiple reactions with
                # different scores
                r = _reaction_of(candidate)
                if not r:
                    continue
                if reaction and reaction.score > r.score:
                    continue
                reaction = r
            if reaction:
                reactions.append(reaction)
        log.info("Found %d reactions", len(reactions))
        for r in reactions:
            log.info(
                "  score=%.4f feasibility=%.4f reactants=%s",
                r.score,
                r.feasibility,
                ".".join(r.smiles),
            )
        return reactions, None


def _reaction_of(r: ai.FixedRetroReaction) -> Reaction | None:
    if not r or not r.metadata or not r.reactants:
        return None
    meta = r.metadata
    score = meta.get("policy_probability", 0.0)
    feasibility = meta.get("feasibility", 0.0)
    smiles = []
    for reactants in r.reactants:
        for mol in reactants:
            smiles.append(mol.smiles)
    if len(smiles) == 0:
        return None
    return Reaction(score=score, feasibility=feasibility, smiles=smiles)


@dataclass
class ZynthConfig:
    endpoint: str
    token: str

    @staticmethod
    def from_file(path: str | PathLike) -> "ZynthConfig":
        with open(path, "r", encoding="utf-8") as f:
            d: dict = json.load(f)
            return ZynthConfig(**d)


class ZynthClient(RetroTool):
    def __init__(self, config: ZynthConfig):
        self.endpoint = config.endpoint.strip().rstrip("/")
        self.session = requests.Session()
        self.session.headers["X-API-Token"] = config.token

    def __enter__(self) -> "ZynthClient":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    @override
    def expand(self, smiles: str) -> Res[list[Reaction]]:
        try:
            url = self.endpoint + "/expand"
            resp = self.session.post(url, json={"smiles": smiles})
            resp.raise_for_status()
            return [Reaction(**d) for d in resp.json()], nil
        except Exception as err:
            return nil, str(err)

    def close(self) -> None:
        self.session.close()
