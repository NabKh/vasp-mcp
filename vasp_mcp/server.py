"""MCP server: lets an MCP client check VASP inputs and read VASP results.

Every tool only reads files. Nothing is written, run or uploaded, and POTCAR
contents are never returned (only the dataset names and cutoffs).
"""

import argparse
import re
from dataclasses import asdict
from pathlib import Path

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations
from vasp_lint import RULES, load, run
from vasp_lint.cli import find_calculations

from . import __version__
from .outputs import FLOAT, summarize

READ_ONLY = ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False)
SKIP_DIRS = {".git", "__pycache__", "node_modules"}
TAIL_BYTES = 256 * 1024

server = MCPServer(
    name="vasp",
    version=__version__,
    instructions=(
        "Tools for VASP calculations on this machine. Use check_vasp_inputs before a job is submitted, "
        "summarize_vasp_run to see whether a finished or crashed job converged and what it gave, and "
        "scan_vasp_project to get an overview of many calculations. Paths are folders on the local disk. "
        "Energies are in eV, forces in eV/Angstrom, pressures in kB. Quote the VASP wiki links returned "
        "with each input finding."
    ),
)
ROOTS = []  # set by --root: folders the tools may read


def _folder(path):
    folder = Path(path).expanduser().resolve()
    if ROOTS and not any(folder == r or r in folder.parents for r in ROOTS):
        raise ToolError(f"{folder} is outside the folders this server may read: {', '.join(map(str, ROOTS))}")
    if folder.is_file():
        folder = folder.parent
    if not folder.is_dir():
        raise ToolError(f"{folder} is not a folder")
    return folder


@server.tool(annotations=READ_ONLY)
def check_vasp_inputs(path: str, ignore: list[str] | None = None) -> dict:
    """Check the VASP inputs (INCAR, POSCAR, POTCAR, KPOINTS) in a folder for common mistakes.

    If the folder has no INCAR, every calculation folder below it is checked.
    Each finding has a rule ID, a severity (error, warning or info), the INCAR
    line, an explanation and the VASP wiki page it is based on.
    """
    folder = _folder(path)
    skip = {s.upper() for s in ignore or []}
    calculations = []
    for calc_dir in find_calculations([folder]):
        calc = load(calc_dir)
        calculations.append(
            {
                "path": str(calc_dir),
                "missing_files": [n for n in ("POSCAR", "POTCAR", "KPOINTS") if getattr(calc, n.lower()) is None],
                "findings": [asdict(f) for f in run(calc, skip)],
            }
        )
    if not calculations:
        raise ToolError(f"no INCAR in or below {folder}")
    return {"calculations": calculations}


@server.tool(annotations=READ_ONLY)
def summarize_vasp_run(path: str) -> dict:
    """Summarize what happened in one VASP calculation folder.

    Reads OUTCAR (or OSZICAR), vasprun.xml and the job's standard output. Reports
    whether the job finished, whether the electronic and ionic loops converged,
    the final energies, the largest force on a free atom, pressure, Fermi level,
    magnetisation, MD temperatures, the band gap per spin channel, VASP warnings,
    known error messages, and a plain list of problems.
    """
    return summarize(_folder(path))


@server.tool(annotations=READ_ONLY)
def read_vasp_inputs(path: str) -> dict:
    """The parsed inputs of one calculation: INCAR tags, structure, POTCAR datasets and k-points."""
    folder = _folder(path)
    if not (folder / "INCAR").is_file():
        raise ToolError(f"no INCAR in {folder}")
    calc = load(folder)
    out = {"path": str(folder), "incar": calc.incar.tags}
    if calc.poscar:
        out["structure"] = {
            "species": calc.poscar.species,
            "counts": calc.poscar.counts,
            "nions": calc.poscar.nions,
            "lattice_A": calc.poscar.lattice,
        }
    if calc.potcar:
        out["potcar"] = [{"symbol": p.symbol, "zval": p.zval, "enmax": p.enmax, "enmin": p.enmin} for p in calc.potcar]
    if calc.kpoints:
        out["kpoints"] = {k: v for k, v in asdict(calc.kpoints).items() if v is not None}
    mesh = calc.kmesh()
    if mesh:
        out["kmesh"] = mesh
    out["has_kpoints_opt"] = calc.has_kpoints_opt
    return out


def _quick_status(folder):
    """Status from the end of OUTCAR only, so that large projects scan fast."""
    entry = {"path": str(folder)}
    outcar = folder / "OUTCAR"
    if not outcar.is_file() or outcar.stat().st_size == 0:
        entry["status"] = "not run" if (folder / "INCAR").is_file() else "no OUTCAR"
        return entry
    with open(outcar, "rb") as f:
        f.seek(max(0, outcar.stat().st_size - TAIL_BYTES))
        tail = f.read().decode(errors="replace")
    entry["status"] = "finished" if "General timing and accounting informations" in tail else "incomplete"
    energies = re.findall(rf"energy\(sigma->0\)\s*=\s*({FLOAT})", tail)
    if energies:
        entry["energy_sigma_to_0_eV"] = float(energies[-1])
    if "reached required accuracy" in tail:
        entry["ionic_converged"] = True
    return entry


@server.tool(annotations=READ_ONLY)
def scan_vasp_project(path: str, max_folders: int = 500) -> dict:
    """List every VASP calculation folder below a path with a quick status.

    For each folder with an INCAR or OUTCAR: finished, incomplete or not run, the
    last energy(sigma->0) and whether a relaxation reached the required accuracy.
    Use summarize_vasp_run on a folder for the full picture.
    """
    root = _folder(path)
    folders = set()
    for name in ("INCAR", "OUTCAR"):
        for p in root.rglob(name):
            parts = p.relative_to(root).parts
            if not any(part.startswith(".") or part in SKIP_DIRS for part in parts):
                folders.add(p.parent)
    folders = sorted(folders)
    entries = [_quick_status(f) for f in folders[:max_folders]]
    counts = {}
    for e in entries:
        counts[e["status"]] = counts.get(e["status"], 0) + 1
    return {"root": str(root), "total": len(folders), "shown": len(entries), "counts": counts, "folders": entries}


@server.tool(annotations=READ_ONLY)
def list_vasp_lint_rules() -> dict:
    """All input checks with their severity and the VASP wiki page each is based on."""
    return {"rules": [{"id": r.id, "severity": r.severity, "title": r.title, "url": r.url} for r in RULES]}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="vasp-mcp", description="MCP server for VASP inputs and outputs.")
    parser.add_argument("--root", action="append", default=[], help="only allow reading below this folder (repeatable)")
    parser.add_argument("--transport", choices=["stdio", "streamable-http"], default="stdio")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args(argv)
    ROOTS[:] = [Path(r).expanduser().resolve() for r in args.root]
    server.run(args.transport)
