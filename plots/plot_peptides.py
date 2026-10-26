from pathlib import Path

from analytes import YAAF, YAAF_TPSA, YTSF, YTSF_TPSA, yaaf_TPSA, ytsf_TPSA, ytsf_nboc, ytsf_ogalnac, ytsf_ps
from electrolyte import Electrolyte
from sweeps import plot_analyte_sweep_comparison, plot_tpsa_comparison
from solver import solve


FIGS = Path(__file__).parent / "figs"


if __name__ == "__main__":
    for analyte, tpsa in ((YTSF, YTSF_TPSA), (YAAF, YAAF_TPSA), (ytsf_TPSA, YTSF_TPSA), (yaaf_TPSA, YAAF_TPSA)):
        print(f"{analyte.name}: TPSA = {tpsa:.0f} Å², site density = {analyte.site_density:.3e} m^-2") 

    plot_analyte_sweep_comparison(
        [YTSF, YAAF, ytsf_TPSA, yaaf_TPSA],
        ["YTSF baseline", "YAAF baseline", "YTSF TPSA", "YAAF TPSA"],
        Electrolyte(ionic_strength=200.0, permittivity=80.0, c_stern=1.1),
        FIGS / "ytsf_yaaf_tpsa_comparison.png",
        boltzmann_ph=False,
        stern=True,
    )

    plot_analyte_sweep_comparison(
        [ytsf_nboc, ytsf_ogalnac, ytsf_TPSA, ytsf_ps],
        ["YTSF + N-Boc", "YTSF + O-GalNAc", "YTSF (unmodified)", "YTSF + phospho-Ser"],
        Electrolyte(ionic_strength=200.0, permittivity=80.0, c_stern=1.1),
        FIGS / "ytsf_modifications.png",
        boltzmann_ph=False,
        stern=True,
    )

    