from pathlib import Path

from analyte import YAAF, YAAF_TPSA, YTSF, YTSF_TPSA, yaaf_TPSA, ytsf_TPSA
from electrolyte import Electrolyte
from sweeps import plot_analyte_sweep_comparison, plot_tpsa_comparison


if __name__ == "__main__":
    for analyte, tpsa in ((YTSF, YTSF_TPSA), (YAAF, YAAF_TPSA), (ytsf_TPSA, YTSF_TPSA), (yaaf_TPSA, YAAF_TPSA)):
        print(f"{analyte.name}: TPSA = {tpsa:.0f} Å², site density = {analyte.site_density:.3e} m^-2")

    plot_analyte_sweep_comparison(
        [YTSF, YAAF, ytsf_TPSA, yaaf_TPSA],
        ["YTSF baseline", "YAAF baseline", "YTSF TPSA", "YAAF TPSA"],
        Electrolyte(ionic_strength=10.0, permittivity=80.0),
        Path("figs/ytsf_yaaf_tpsa_comparison.png"),
    )

    plot_tpsa_comparison(
        YTSF,
        ytsf_TPSA,
        yaaf_TPSA,
        YTSF_TPSA,
        YAAF_TPSA,
        Electrolyte(ionic_strength=10.0, permittivity=80.0),
        Path("figs/ytsf_yaaf_tpsa_two_panel.png"),
    )