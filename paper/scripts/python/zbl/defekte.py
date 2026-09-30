"""defekte.py -- Wigner-Seitz-Defektanalyse fuer cascade_serial-Endzustaende.

Jedes Atom wird dem naechsten Gitterplatz des Ausgangsgitters zugeordnet.
Leerplatz  = Gitterplatz ohne Atom.
Zwischengitteratom = zweites (drittes, ...) Atom auf einem Platz, oder ein
Atom, dessen naechster Platz weiter als 0,5 L0 entfernt ist.
Gezaehlt werden nur Plaetze mit mindestens `rand` Zeilen/Spalten Abstand
zum Gitterrand (freie Raender relaxieren und wuerden Scheindefekte liefern).
"""
import numpy as np
from scipy.spatial import cKDTree

DY = np.sqrt(3.0) / 2.0


def gitterplaetze(nx, ny):
    j, i = np.mgrid[0:ny, 0:nx]
    x = i + 0.5 * (j & 1)
    y = j * DY
    return np.column_stack([x.ravel(), y.ravel()]), i.ravel(), j.ravel()


def lies_endzustand(pfad):
    a = np.loadtxt(pfad, delimiter=",", skiprows=1, ndmin=2)
    return a[:, :2]


def ws_analyse(pos, nx, ny, rand=3):
    """Liefert (leerplaetze, zwischengitter) im Innenbereich."""
    plaetze, i, j = gitterplaetze(nx, ny)
    innen = (i >= rand) & (i < nx - rand) & (j >= rand) & (j < ny - rand)
    d, k = cKDTree(plaetze).query(pos)
    belegung = np.bincount(k, minlength=len(plaetze))
    leer = int(np.sum((belegung == 0) & innen))
    zw = int(np.sum(np.clip(belegung - 1, 0, None)[innen]))
    return leer, zw
