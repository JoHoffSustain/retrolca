import base64
import io
import logging
from typing import Any

import olca_ipc as ipc
import olca_schema as o
from rdkit import Chem
from rdkit.Chem import Draw

from . import losses, oipc, smiles, tool
from .naming import CIR, NamingService
from .res import Res, nil, unwrap

log = logging.getLogger(__name__)


def _find_provider(
    process_uid: str, providers: list[o.TechFlow]
) -> o.TechFlow | None:
    for p in providers:
        if p.provider and p.provider.id == process_uid:
            return p
    log.warning("Could not find process: %s", process_uid)
    return None


def _is_waste(ctx: oipc.IpcContext, p: o.TechFlow | None) -> bool:
    if not p or not p.provider or not p.flow:
        return False
    flow_type = p.flow.flow_type
    if flow_type is None and p.flow.id:
        # references do not always carry the flow type
        flow = ctx.client.get(o.Flow, p.flow.id)
        flow_type = flow.flow_type if flow else None
        p.flow.flow_type = flow_type
    return flow_type == o.FlowType.WASTE_FLOW


def _provider_name(p: o.TechFlow) -> str:
    return p.provider.name if p.provider and p.provider.name else "?"


class ProcessBuilder:
    def __init__(
        self,
        ctx: oipc.IpcContext,
        tool: tool.RetroTool,
        max_levels=3,
        max_variants=1,
        gen_process: str | None = None,
        bal_process: str | None = None,
        naming: NamingService = CIR(),
        yield_: float | None = 0.95,
        emission_share: float = 0.05,
        waste_process: str | None = None,
    ):
        """Constructs a new process builder.

        Args:
            ctx:
                The IPC context for data exchange with openLCA.
            tool:
                The retrosynthesis tool.
            max_levels:
                The maximum number of levels, or supply-chain depth, to
                generate. At each level, the builder tries to link providers.
                If no provider can be linked earlier, generation continues up
                to this depth. Default is 3.
            max_variants:
                The maximum number of process variants that can be created at
                each level. Default is 1.
            gen_process:
                An optional ID of a generic chemical production process that
                should be linked to the generated processes. This process needs
                to have a single product output measured in mass. It is linked
                as an input to the generated processes. Each generated process
                has an output of 1 kg of the respective product. Thus, the
                generic production process is also linked with 1 kg input.
            bal_process:
                An optional ID of a "balancing" process. This process needs to
                have a reference flow, a waste input or product output, given in
                mass. For a generated process of a reaction, if the mass of the
                reactants is smaller than the mass of the output, the balancing
                process is linked with that mass difference. Typically a waste
                treatment process is used in this case and the mass difference
                is then added as waste output linked to that waste treatment
                process. If a production process is used, the difference is
                linked as avoided product output to the balancing process.
            naming:
                The naming service that should be used to translate SMILES codes
                into names of chemicals. By default CIR is used if not specified
                otherwise.
            yield_:
                The reaction yield, default 0.95. The reactant inputs are the
                stoichiometric amounts divided by the yield. The value is
                stored as global parameter `retro_yield` in openLCA and the
                exchange amounts are formulas of it; if that parameter already
                exists, its value is used. `None` switches the yield losses
                off (100% yield, as in earlier versions).
            emission_share:
                The share of the yield excess (the unreacted reactants) that
                is emitted directly, default 0.05; the rest goes to the waste
                process. Stored as global parameter `retro_em_share`. The
                emissions are mapped to elementary flows (emission to air and
                water, unspecified) by CAS number, or by name and formula. If
                an emission to air exists, the emission is split 50/50 between
                air and water, otherwise it goes to water. If no flow is found,
                a new elementary flow (emission to water) without
                characterisation factors is created.
            waste_process:
                An optional ID of a waste treatment process for the waste part
                of the yield excess. If not given, the balancing process is
                used when it is a waste treatment process.
        """
        self.ctx = ctx
        log.info("Build provider and flow index")
        all_flows = ctx.client.get_all(o.Flow)
        self.providers = oipc.ProviderIndex.of(ctx, all_flows)
        self.flows = oipc.FlowIndex.of(ctx, all_flows)
        self.tool = tool
        self.max_variants = max_variants
        self.max_levels = max_levels

        self.gen_provider: o.TechFlow | None = None
        self.bal_provider: o.TechFlow | None = None
        waste_provider: o.TechFlow | None = None
        if gen_process or bal_process or waste_process:
            providers = ctx.client.get_providers()
            if gen_process:
                self.gen_provider = _find_provider(gen_process, providers)
            if bal_process:
                self.bal_provider = _find_provider(bal_process, providers)
            if waste_process:
                waste_provider = _find_provider(waste_process, providers)
        if waste_provider and not _is_waste(ctx, waste_provider):
            log.error("The waste process must treat a waste flow; ignored")
            waste_provider = None
        if not waste_process and _is_waste(ctx, self.bal_provider):
            waste_provider = self.bal_provider
        self.naming = naming

        self.yield_model: losses.YieldModel | None = None
        if yield_ is not None:
            self.yield_model = losses.YieldModel(
                ctx,
                yield_,
                emission_share,
                waste_provider,
                naming,
                all_flows,
            )

    def expand_process(self, process_id: str, category: str | None):
        process = self.ctx.client.get(o.Process, process_id)
        if not process:
            log.error("could not find process: %s", process_id)
            return
        if not process.exchanges:
            return
        expanded = 0
        log.info("expand process: %s", process.name)
        for e in process.exchanges:
            if (
                not e.is_input
                or e.default_provider
                or e.is_avoided_product
                or not e.flow
            ):
                continue
            flow = self.ctx.client.get(o.Flow, e.flow.id)
            if not flow or flow.flow_type != o.FlowType.PRODUCT_FLOW:
                continue
            smiles_code = smiles.of_flow(flow)
            if not smiles_code:
                continue

            log.info("expand input of: %s", flow.name)
            tech_flow = self.providers.get(smiles_code)
            if not tech_flow:
                self.build(smiles_code, flow.name, category, 0)
                tech_flow = self.providers.get(smiles_code)
            if not tech_flow or not tech_flow.provider:
                continue
            e.default_provider = tech_flow.provider
            expanded += 1

        if expanded == 0:
            log.info("could not expand any inputs in: %s", process.name)
            return

        log.info(
            "update process %s with %d new providers", process.name, expanded
        )
        self.ctx.client.put(process)

    def build(
        self,
        smiles_code: str,
        name: str | None = None,
        category: str | None = None,
        level=0,
        expanded_reactions: frozenset[str] = frozenset(),
    ) -> list[o.Process]:
        """Build process variants for a product SMILES code.

        Args:
            smiles_code:
                The product SMILES code for which retrosynthesis processes
                should be generated.
            name:
                An optional display name for the product flow. If omitted, the
                naming service is used and falls back to the SMILES code.
            category:
                An optional root category path under which generated processes,
                flows, and sources are stored in openLCA. For generated
                processes, subcategories for the respective levels are created
                under that path.
            level:
                The current recursion level. This is managed internally during
                recursive expansion and usually stays at the default value.
            expanded_reactions:
                The reaction SMILES keys that have already been expanded on
                the current path to this product. Reactions that are already
                in this set are not expanded again, which avoids cycles in
                the generated process graph. This is managed internally
                during recursive expansion and usually stays at the default
                value.
        """
        smiles_code = smiles.canonicalize(smiles_code)
        log.info(
            "Create processes for SMILES=%s at level=%d", smiles_code, level
        )
        ref_flow = self.__resolve_product(smiles_code, name, category)
        if not ref_flow:
            return []
        mm_prod = self.ctx.molar_mass_of(ref_flow)
        if not mm_prod:
            return []
        reactions = self.__get_reactions(smiles_code, expanded_reactions)
        if len(reactions) == 0:
            return []

        variant = 0
        processes = []
        for reaction in reactions:
            variant += 1
            reaction_key = smiles.reaction_smiles(reaction.smiles, smiles_code)
            next_expanded = expanded_reactions | {reaction_key}
            process = self.__init_process(
                ref_flow, reaction, variant, smiles_code, category, level
            )

            input_mass = 0.0
            reactants: list[losses.Reactant] = []
            for si in reaction.smiles:
                smiles_i = smiles.canonicalize(si)
                provider = self.__resolve_provider(
                    smiles_i, category, level + 1, next_expanded
                )
                if not provider:
                    log.warning(
                        "Could not resolve an input flow for SMILES=%s",
                        smiles_i,
                    )
                    continue
                flow = unwrap(provider.flow)
                mm_react = self.ctx.molar_mass_of(flow)
                if not mm_react:
                    log.warning(
                        "No molar mass for input flow '%s'; skipping input",
                        flow.name,
                    )
                    continue

                # with n: chemical amount, m: mass, and mm: molar mass
                # as we have a reaction of 1 mol reactant (input) and
                # 1 mol product (output) with mm = m / n
                # n_inp = n_out = 1 mol
                # n_out = m_out / mm_out
                # m_inp = n_inp * mm_inp
                # m_inp = n_out * mm_inp
                # m_inp = (m_out / mm_out) * mm_inp
                # m_inp = m_out * mm_inp / mm_out  | m_out = 1 (kg)
                # m_inp = mm_inp / mm_out
                amount = mm_react / mm_prod
                input_mass += amount
                inp = o.new_input(process, flow, amount, self.ctx.kg)
                inp.flow_property = self.ctx.mass.to_ref()
                inp.default_provider = provider.provider
                if self.yield_model:
                    inp.amount, inp.amount_formula = (
                        self.yield_model.input_amount(amount)
                    )
                    reactants.append(
                        losses.Reactant(
                            smiles=smiles_i,
                            flow=flow,
                            product_flow=self.flows.data.get(smiles_i),
                            stoich_amount=amount,
                        )
                    )

            # the balancing flow covers the by-products of the reaction
            # (e.g. water of a condensation); it is calculated from the
            # stoichiometric amounts only, without the yield excess
            if input_mass > 1.0 and self.bal_provider:
                self.__add_bal_exchange(process, input_mass - 1.0)
            if self.yield_model:
                self.yield_model.add_losses(process, reactants)
                process.other_properties = process.other_properties or {}
                process.other_properties["Retrosynthesis-Yield"] = (
                    self.yield_model.yield_value
                )
                process.other_properties["Retrosynthesis-Emission-Share"] = (
                    self.yield_model.emission_share
                )
            self.__check_mass_balance(process)

            # Store the process only after its inputs have been resolved and
            # register it as the provider for this product afterwards. If the
            # provider was registered earlier, a reaction that links back to
            # an ancestor would be short-circuited by the provider index and
            # the path filter could not cut the cycle.
            self.ctx.client.put(process)
            if variant == 1:
                self.providers.put(smiles_code, process)
            log.info("Stored process %s", process.name)
            processes.append(process)

        if level == 0 and self.yield_model:
            self.yield_model.log_summary()
        return processes

    def __check_mass_balance(self, process: o.Process) -> None:
        """Compares the mass of the inputs and outputs of a process.

        The input of the generic production process is not part of the
        balance. Logs a warning when the balance does not close, which
        typically means that
        the retrosynthesis tool did not return all reactants (e.g. hydrogen
        in reductions) or that no balancing or waste process is set.
        """
        gen_id = (
            self.gen_provider.provider.id
            if self.gen_provider and self.gen_provider.provider
            else None
        )
        m_in, m_out = 0.0, 0.0
        for e in process.exchanges or []:
            if (
                e.is_input
                and gen_id
                and e.default_provider
                and e.default_provider.id == gen_id
            ):
                continue
            amount = e.amount or 0.0
            if e.is_input != bool(e.is_avoided_product):
                m_in += amount
            else:
                m_out += amount
        gap = m_in - m_out
        if abs(gap) <= 1e-6 * max(m_in, m_out, 1.0):
            return
        if gap < 0:
            hint = (
                "the reactants are probably incomplete (e.g. hydrogen or "
                "a reagent not returned by the retrosynthesis tool)"
            )
        else:
            hint = "no balancing or waste process takes up the difference"
        log.warning(
            "Mass balance of '%s' does not close: in %.5f kg, out %.5f kg "
            "(gap %.5f kg); %s",
            process.name,
            m_in,
            m_out,
            abs(gap),
            hint,
        )

    def __get_reactions(
        self, smiles_code: str, expanded_reactions: frozenset[str]
    ) -> list[tool.Reaction]:
        reactions, err = self.tool.expand(smiles_code)
        if err:
            log.info(
                "No retrosynthesis results retrieved for %s: %s",
                smiles_code,
                err,
            )
            return []
        if not isinstance(reactions, list) or len(reactions) == 0:
            log.info("No retrosynthesis results retrieved for: %s", smiles_code)
            return []

        # Reactions that have already been expanded on the current path
        # would close a cycle in the generated process graph and are
        # therefore ignored. The filter is applied before the variants are
        # truncated so that a cycle does not consume a variant slot.
        skipped: list[str] = []
        candidates: list[tool.Reaction] = []
        for reaction in reactions:
            key = smiles.reaction_smiles(reaction.smiles, smiles_code)
            if key in expanded_reactions:
                skipped.append(key)
            else:
                candidates.append(reaction)
        if skipped:
            log.info(
                "Ignoring %d reaction(s) already expanded on the current "
                "path for SMILES=%s: %s",
                len(skipped),
                smiles_code,
                ", ".join(skipped),
            )
        if not candidates:
            log.info(
                "All reactions for SMILES=%s are on the current path; "
                "cutting the cycle",
                smiles_code,
            )
            return []

        candidates.sort(key=lambda r: r.score * r.feasibility, reverse=True)
        if len(candidates) > self.max_variants:
            candidates = candidates[0 : self.max_variants]
        return candidates

    def __init_process(
        self,
        ref_flow: o.Flow,
        reaction: tool.Reaction,
        variant: int,
        product_smiles: str,
        category: str | None,
        level: int,
    ) -> o.Process:
        score = reaction.score * reaction.feasibility
        name = (
            f"production of {ref_flow.name} | variant {variant} "
            f"| c = {score:.4f}"
        )
        process = o.new_process(name)
        process.category = category
        if level > 0 and process.category:
            process.category += f"/level {level}"

        qref = o.new_output(process, ref_flow, 1, self.ctx.kg)
        qref.flow_property = self.ctx.mass.to_ref()
        qref.is_quantitative_reference = True
        process.other_properties = {
            "Retrosynthesis-Score": reaction.score,
            "Retrosynthesis-Feasibility": reaction.feasibility,
        }
        self.__add_reaction_source(process, product_smiles, reaction, category)
        self.__add_gen_input(process)
        return process

    def __resolve_provider(
        self,
        smiles_code: str,
        category: str | None,
        level: int,
        expanded_reactions: frozenset[str] = frozenset(),
    ) -> o.TechFlow | None:
        if p := self.providers.get(smiles_code):
            return p
        if level <= self.max_levels:
            processes = self.build(
                smiles_code,
                category=category,
                level=level,
                expanded_reactions=expanded_reactions,
            )
            if len(processes) > 0:
                return self.providers.put(smiles_code, processes[0])
            log.info(
                "No process generated for SMILES=%s; linking the product "
                "flow as an input without a provider",
                smiles_code,
            )
        # Either the maximum depth is reached or no process could be generated
        # for this chemical (e.g. the retrosynthesis tool returned no
        # reaction). In both cases we still link the product flow as an input
        # without a provider so that the generated process does not silently
        # miss this input.
        flow = self.__resolve_product(smiles_code, category=category)
        return o.TechFlow(flow=flow.to_ref()) if flow else None

    def __resolve_product(
        self,
        smiles_code: str,
        name: str | None = None,
        category: str | None = None,
    ) -> o.Flow | None:
        flow = self.flows.data.get(smiles_code)
        if flow:
            return flow
        log.info("Create flow for SMILES: %s", smiles_code)
        flow, err = self.create_product(smiles_code, name, category)
        if err:
            log.error("Failed to create flow for SMILES: %s", smiles_code)
            return None
        assert flow
        self.flows.data[smiles_code] = flow
        return flow

    def __add_reaction_source(
        self,
        process: o.Process,
        smiles_code: str,
        reaction: tool.Reaction,
        category: str | None,
    ):
        image_data = self.__reaction_image(smiles_code, reaction)
        if not image_data:
            return
        source = o.Source(name=f"Reaction scheme for {process.name}")
        source.category = category
        source.description = (
            "This source was automatically generated for the retrosynthesis "
            f"reaction assigned to the process '{process.name}'. It includes "
            "an image showing the reactant molecules together with the product "
            "molecule."
        )
        source_ref = self.ctx.client.put(source)
        if not source_ref:
            log.error("Failed to create reaction source for %s", process.name)
            return

        uploaded = self.ctx.client.put_source_file(
            source_ref, ipc.FileData("reaction.png", image_data)
        )
        if not uploaded:
            log.error(
                "Failed to upload reaction image for process %s",
                process.name,
            )

        process.process_documentation = o.ProcessDocumentation()
        process.process_documentation.sources = [source_ref]

    def __reaction_image(
        self, product_smiles: str, reaction: tool.Reaction
    ) -> str | None:
        codes = [smiles.canonicalize(code) for code in reaction.smiles]
        codes.append(product_smiles)

        mols = []
        labels = []
        for code in codes:
            mol = Chem.MolFromSmiles(code)
            if not mol:
                log.warning("Could not render molecule for SMILES: %s", code)
                continue
            mols.append(mol)
            labels.append(self.naming.get_name(code))

        if len(mols) == 0:
            return None

        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=len(mols),
            subImgSize=(200, 200),
            legends=labels,
        )
        buffer = io.BytesIO()
        img.save(buffer, format="PNG")
        return base64.b64encode(buffer.getvalue()).decode("ascii")

    def __add_gen_input(self, process: o.Process):
        if not self.gen_provider or not self.gen_provider.flow:
            return
        inp = o.new_input(
            process,
            unwrap(self.gen_provider.flow),
            1,
            self.ctx.kg,
        )
        inp.flow_property = self.ctx.mass.to_ref()
        inp.default_provider = self.gen_provider.provider

    def __add_bal_exchange(self, process: o.Process, amount: float):
        if not self.bal_provider or not self.bal_provider.flow:
            return
        flow = unwrap(self.bal_provider.flow)
        exchange: o.Exchange
        if flow.flow_type == o.FlowType.PRODUCT_FLOW:
            exchange = o.new_input(process, flow, amount, self.ctx.kg)
            exchange.is_avoided_product = True
            exchange.description = (
                "Excess Process Mass (by-products) as avoided product"
            )
        else:
            exchange = o.new_output(process, flow, amount, self.ctx.kg)
            exchange.description = (
                "Excess Process Mass (by-products) to waste treatment"
            )
        exchange.flow_property = self.ctx.mass.to_ref()
        exchange.default_provider = self.bal_provider.provider

    def create_product(
        self,
        smiles_code: str,
        name: str | None = None,
        category: str | None = None,
    ) -> Res[o.Flow]:
        info = self.naming.get_info(smiles_code)
        product: str
        if name:
            product = name
        elif info:
            product = info.name
        else:
            product = smiles_code

        flow = o.new_product(product, self.ctx.mass)
        if not flow or not flow.flow_properties:
            return nil, "Could not create product flow"
        flow.category = category
        flow.description = (
            "This product flow was automatically generated from it's SMILES code. "
            "See also see the additional properties of the flow for more "
            "information."
        )

        # add the chemical amount as flow property
        mw = smiles.mol_weight(smiles_code)
        if not mw:
            return nil, f"Could not calculate the molar mass of: {smiles_code}"
        flow.flow_properties.append(
            o.FlowPropertyFactor(
                conversion_factor=1000 / mw,
                flow_property=self.ctx.chem_amount.to_ref(),
            )
        )

        # additional properties
        props: dict[str, Any] = {}
        flow.other_properties = props
        props["SMILES"] = smiles_code
        props["MolarMass"] = mw
        if info:
            flow.formula = info.formula
            if info.inchi:
                props["InChI-String"] = info.inchi
            if info.inchi_key:
                props["InChI-Key"] = info.inchi_key

        self.ctx.client.put(flow)
        return flow, nil
