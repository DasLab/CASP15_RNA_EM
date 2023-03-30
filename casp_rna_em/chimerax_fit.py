# chimerax_search_fit.py PDB.pdb MRC.mrc out.pdb threhsold

from chimerax.core.commands import run
import sys

run(session, f"open {sys.argv[1]}")
run(session, f"open {sys.argv[2]}")
run(session, f"vol #2 level {sys.argv[4]}")
run(session, f"fit #1 in #2")
run(session, f"save {sys.argv[3]} #1")
run(session, "exit")
