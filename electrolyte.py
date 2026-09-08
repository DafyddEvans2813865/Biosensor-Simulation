from dataclasses import dataclass

from scipy.constants import e as Q          # elementary charge, C
from scipy.constants import k as K_B        # Boltzmann, J/K
from scipy.constants import epsilon_0 as EPS_0
from scipy.constants import N_A             # Avogadro


@dataclass(frozen=True)
class Electrolyte:
	ionic_strength: float  # mol/m^3 (SI); convert mM at the CLI boundary
	temperature: float = 298.15
	kappa_w: float = 80.0
	c_stern: float | None = None

	@property
	def thermal_voltage(self) -> float:
		return K_B * self.temperature / Q

	@property
	def debye_length(self) -> float:
		return (
			self.kappa_w
			* EPS_0
			* K_B
			* self.temperature
			/ (2.0 * Q**2 * self.ionic_strength * N_A)
		) ** 0.5

	@property
	def q0(self) -> float:
		return (
			8.0
			* self.kappa_w
			* EPS_0
			* K_B
			* self.temperature
			* self.ionic_strength
			* N_A
		) ** 0.5
