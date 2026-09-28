from dataclasses import dataclass
from math import log 
from typing import Protocol

from scipy.constants import e as Q

LN10 = log(10.0)

class ChargedSurface(Protocol):
    def charge(self,ph_surface: float) -> float:
        pass #interface

@dataclass(frozen=True)
class AmphotericSurface:
    pka: float 
    pkb: float 
    site_density: float 

    @property
    def ka(self) -> float:
        return 10.0 ** -self.pka

    @property
    def kb(self) -> float:
        return 10.0 ** -self.pkb

    @property
    def pzc(self) -> float:
        return 0.5 * (self.pka + self.pkb)

    def charge(self,ph_surface:float) -> float:
        return amphoteric_charge(self, ph_surface)


def amphoteric_charge(s: AmphotericSurface, ph_surface: float) -> float:
    a = 10.0 ** -ph_surface
    ka, kb = s.ka, s.kb

    #Positive -> SiOH2+ > SiO^-1 -> Below pzc
    #Zero -> at pzc
    #Negative -> SiO^-1 > SiOH2+ -> above pzc
    return Q * s.site_density * (a**2 - ka * kb) / (ka * kb + kb * a + a**2)

def beta_int(s: AmphotericSurface, ph_surface: float) -> float:
    # Intrinsic buffer capacity 
    a = 10.0 ** -ph_surface
    ka, kb = s.ka, s.kb
    num = kb * a**2 + 4 * ka * kb * a + ka * kb**2
    den = (ka * kb + kb * a + a**2) ** 2
    return LN10 * s.site_density * a * num / den #how willing the surface is to change charge - higher means reacts strongly to ph change 

def alpha(beta: float, c_dif: float, thermal_voltage: float) -> float:
    #fraction of the ideal 59 mV/pH response (1-0) 
    return 1.0 / (LN10 * thermal_voltage * c_dif / (Q * beta) + 1.0)

