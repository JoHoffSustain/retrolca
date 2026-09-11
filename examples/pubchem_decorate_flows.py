"""Decorate openLCA product flows with PubChem data.

This example enriches flows in an openLCA database with PubChem-derived
chemical information and can reuse or persist those decorations as a JSON
dump for later import into another database.
"""

import logging as log
from pathlib import Path

import olca_ipc as ipc

import retrolca as retro
import retrolca.pubchem as pub

DUMP = Path(__file__).parent.parent / "out/pubchem_deco.json"


def main():
    log.basicConfig(level=log.INFO)

    # load the IPC context
    client = ipc.Client()
    ctx, err = retro.IpcContext.of(client)
    if err:
        print(f"Failed to load context: {err}")
        return
    assert ctx

    # if a dump with PubChem data exists, we apply this directly
    if DUMP.exists():
        pub.load_decorations(ctx, DUMP)
        return

    # otherwise, collect chemical information on PubChem, decorate the flows,
    # and store the dump
    pub.IpcFlowDecorator(ctx).try_all(in_path="manufacture of chemicals and chemical products")
    pub.dump_decorations(ctx, DUMP)


if __name__ == "__main__":
    main()
