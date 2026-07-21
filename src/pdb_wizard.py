#!/usr/bin/env python3
# PDB Wizard v0.5.0 — Amalgamated Single-File Distribution
# Based on pdb_wizard by Adam Hogan (GPLv3) and moltui by Kalman Szenes (MIT)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
from __future__ import annotations

import os as _os
import subprocess as _subprocess
import sys as _sys

_os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
_os.environ.setdefault("MKL_NUM_THREADS", "1")
_os.environ.setdefault("OMP_NUM_THREADS", "1")

def _ensure_installed(package, import_name=None):
    import_name = import_name or package
    try:
        __import__(import_name)
        return True
    except ImportError:
        pass
    print(f"\'{package}\' is not installed.")
    try:
        answer = input(f"Install it now with pip? [Y/n] ").strip().lower()
    except EOFError:
        answer = "n"
    if answer in ("", "y", "yes"):
        _subprocess.check_call([_sys.executable, "-m", "pip", "install", package, "--quiet"])
        return True
    return False

if not _ensure_installed("numpy"):
    _sys.exit(1)

from dataclasses import dataclass
from numpy import array, cos, floor, pi, sin, sqrt
from dataclasses import dataclass, field
from typing import Optional
import copy
import re
import sys
from pathlib import Path
from typing import Optional, TextIO
import struct
from typing import BinaryIO, Optional
from functools import reduce
from math import gcd
from typing import Any, Iterator, TextIO
import math
from typing import Any
from collections import Counter
import json
import urllib.parse
import urllib.request
import asyncio
import os
import argparse
import subprocess

_HAS_TEXTUAL = False
try:
    from textual.app import ComposeResult
    from textual.binding import Binding
    from textual.message import Message
    from textual.widget import Widget
    from textual.widgets import (     Checkbox,     DataTable,     Label,     RadioButton,     RadioSet,     Select,     Static,     TabbedContent,     TabPane, )
    from rich.segment import Segment
    from rich.style import Style
    from textual.app import App, ComposeResult
    from textual.containers import Horizontal, Vertical
    from textual.events import Key
    from textual.screen import ModalScreen
    from textual.strip import Strip
    from textual.widgets import (     Button,     Checkbox,     DataTable,     DirectoryTree,     Footer,     Input,     Label,     ProgressBar,     RadioButton,     RadioSet,     Select,     TabbedContent,     TextArea, )
    _HAS_TEXTUAL = True
except ImportError:
    class _TextualStub(type):
        def __new__(mcls, name, bases, ns, **kw):
            return super().__new__(mcls, name, bases, ns)
        def __init__(cls, name, bases, ns, **kw):
            super().__init__(name, bases, ns)
        def __getattr__(cls, _n):
            return _TextualStub(_n, (), {})
        def __call__(cls, *a, **k):
            return cls
        def __getitem__(cls, _i):
            return cls
    def _textual_stub(_n):
        return _TextualStub(_n, (), {})
    ComposeResult = _textual_stub('ComposeResult')
    Binding = _textual_stub('Binding')
    Message = _textual_stub('Message')
    Widget = _textual_stub('Widget')
    Checkbox = _textual_stub('Checkbox')
    DataTable = _textual_stub('DataTable')
    Label = _textual_stub('Label')
    RadioButton = _textual_stub('RadioButton')
    RadioSet = _textual_stub('RadioSet')
    Select = _textual_stub('Select')
    Static = _textual_stub('Static')
    TabbedContent = _textual_stub('TabbedContent')
    TabPane = _textual_stub('TabPane')
    Segment = _textual_stub('Segment')
    Style = _textual_stub('Style')
    App = _textual_stub('App')
    Horizontal = _textual_stub('Horizontal')
    Vertical = _textual_stub('Vertical')
    Key = _textual_stub('Key')
    ModalScreen = _textual_stub('ModalScreen')
    Strip = _textual_stub('Strip')
    Button = _textual_stub('Button')
    DirectoryTree = _textual_stub('DirectoryTree')
    Footer = _textual_stub('Footer')
    Input = _textual_stub('Input')
    ProgressBar = _textual_stub('ProgressBar')
    TextArea = _textual_stub('TextArea')

# ======================================================================
# Module: constants
# ======================================================================
"""Unified element data merging pdb_wizard and moltui element tables."""



BOHR_TO_ANGSTROM = 0.529177249


def _bond_r_from_atomic_number(z: int) -> float:
    if z == 1:
        return 0.9
    if z <= 5:
        return 1.2
    if z <= 7:
        return 1.7
    if z <= 9:
        return 1.5
    if z <= 13:
        return 1.2
    if z <= 17:
        return 1.8
    return 1.8


@dataclass(frozen=True)
class Element:
    symbol: str
    name: str
    atomic_number: int
    mass: float
    cpk_color: tuple[int, int, int]
    covalent_radius: float
    vdw_radius: float
    bond_r: float


# fmt: off
ELEMENTS: dict[str, Element] = {
    # Period 1
    "H":  Element("H",  "Hydrogen",      1,   1.00797,   (255, 255, 255), 0.31, 1.20, 0.9),
    "He": Element("He", "Helium",        2,   4.0026,    (217, 255, 255), 0.28, 1.40, 1.2),
    # Period 2
    "Li": Element("Li", "Lithium",       3,   6.941,     (204, 128, 255), 1.28, 1.82, 1.2),
    "Be": Element("Be", "Beryllium",     4,   9.01218,   (194, 255,   0), 0.96, 1.98, 1.2),
    "B":  Element("B",  "Boron",         5,  10.81,      (255, 181, 181), 0.84, 1.92, 1.2),
    "C":  Element("C",  "Carbon",        6,  12.011,     (144, 144, 144), 0.76, 1.70, 1.7),
    "N":  Element("N",  "Nitrogen",      7,  14.0067,    ( 48,  80, 248), 0.71, 1.55, 1.7),
    "O":  Element("O",  "Oxygen",        8,  15.9994,    (255,  13,  13), 0.66, 1.52, 1.5),
    "F":  Element("F",  "Fluorine",      9,  18.998403,  (144, 224,  80), 0.57, 1.47, 1.5),
    "Ne": Element("Ne", "Neon",         10,  20.179,     (179, 227, 245), 0.58, 1.54, 1.2),
    # Period 3
    "Na": Element("Na", "Sodium",       11,  22.98977,   (171,  92, 242), 1.66, 2.27, 1.2),
    "Mg": Element("Mg", "Magnesium",    12,  24.305,     (138, 255,   0), 1.41, 1.73, 1.2),
    "Al": Element("Al", "Aluminum",     13,  26.98154,   (191, 166, 166), 1.21, 1.84, 1.2),
    "Si": Element("Si", "Silicon",      14,  28.0855,    (240, 200, 160), 1.11, 2.10, 1.8),
    "P":  Element("P",  "Phosphorus",   15,  30.97376,   (255, 128,   0), 1.07, 1.80, 1.8),
    "S":  Element("S",  "Sulfur",       16,  32.06,      (255, 255,  48), 1.05, 1.80, 1.8),
    "Cl": Element("Cl", "Chlorine",     17,  35.453,     ( 31, 240,  31), 1.02, 1.75, 1.8),
    "Ar": Element("Ar", "Argon",        18,  39.948,     (128, 209, 227), 1.06, 1.88, 1.8),
    # Period 4
    "K":  Element("K",  "Potassium",    19,  39.0983,    (143,  64, 212), 2.03, 2.75, 1.8),
    "Ca": Element("Ca", "Calcium",      20,  40.08,      ( 61, 255,   0), 1.76, 2.31, 1.8),
    "Sc": Element("Sc", "Scandium",     21,  44.9559,    (230, 230, 230), 1.70, 2.15, 1.8),
    "Ti": Element("Ti", "Titanium",     22,  47.9,       (191, 194, 199), 1.60, 2.11, 1.8),
    "V":  Element("V",  "Vanadium",     23,  50.9415,    (166, 166, 171), 1.53, 2.07, 1.8),
    "Cr": Element("Cr", "Chromium",     24,  51.996,     (138, 153, 199), 1.39, 2.06, 1.8),
    "Mn": Element("Mn", "Manganese",    25,  54.938,     (156, 122, 199), 1.50, 2.05, 1.8),
    "Fe": Element("Fe", "Iron",         26,  55.847,     (224, 102,  51), 1.42, 2.04, 1.8),
    "Co": Element("Co", "Cobalt",       27,  58.9332,    (240, 144, 160), 1.38, 2.00, 1.8),
    "Ni": Element("Ni", "Nickel",       28,  58.7,       ( 80, 208,  80), 1.24, 1.97, 1.8),
    "Cu": Element("Cu", "Copper",       29,  63.546,     (200, 128,  51), 1.32, 1.96, 1.8),
    "Zn": Element("Zn", "Zinc",         30,  65.38,      (125, 128, 176), 1.22, 2.01, 1.8),
    "Ga": Element("Ga", "Gallium",      31,  69.72,      (194, 143, 143), 1.22, 1.87, 1.8),
    "Ge": Element("Ge", "Germanium",    32,  72.59,      (102, 143, 143), 1.20, 2.11, 1.8),
    "As": Element("As", "Arsenic",      33,  74.9216,    (189, 128, 227), 1.19, 1.85, 1.8),
    "Se": Element("Se", "Selenium",     34,  78.96,      (255, 161,   0), 1.20, 1.90, 1.8),
    "Br": Element("Br", "Bromine",      35,  79.904,     (166,  41,  41), 1.20, 1.85, 1.8),
    "Kr": Element("Kr", "Krypton",      36,  83.8,       ( 92, 184, 209), 1.16, 2.02, 1.8),
    # Period 5
    "Rb": Element("Rb", "Rubidium",     37,  85.4678,    (112,  46, 176), 2.20, 3.03, 1.8),
    "Sr": Element("Sr", "Strontium",    38,  87.62,      (  0, 255,   0), 1.95, 2.49, 1.8),
    "Y":  Element("Y",  "Yttrium",      39,  88.9059,    (148, 255, 255), 1.90, 2.19, 1.8),
    "Zr": Element("Zr", "Zirconium",    40,  91.22,      (148, 224, 224), 1.75, 2.23, 1.8),
    "Nb": Element("Nb", "Niobium",      41,  92.9064,    (115, 194, 201), 1.64, 2.18, 1.8),
    "Mo": Element("Mo", "Molybdenum",   42,  95.94,      ( 84, 181, 181), 1.54, 2.17, 1.8),
    "Tc": Element("Tc", "Technetium",   43,  98.0,       ( 59, 158, 158), 1.47, 2.16, 1.8),
    "Ru": Element("Ru", "Ruthenium",    44, 101.07,      ( 36, 143, 143), 1.46, 2.13, 1.8),
    "Rh": Element("Rh", "Rhodium",      45, 102.9055,    ( 10, 125, 140), 1.42, 2.10, 1.8),
    "Pd": Element("Pd", "Palladium",    46, 106.4,       (  0, 105, 133), 1.39, 2.10, 1.8),
    "Ag": Element("Ag", "Silver",       47, 107.868,     (192, 192, 192), 1.45, 2.11, 1.8),
    "Cd": Element("Cd", "Cadmium",      48, 112.41,      (255, 217, 143), 1.44, 2.18, 1.8),
    "In": Element("In", "Indium",       49, 114.82,      (166, 117, 115), 1.42, 1.93, 1.8),
    "Sn": Element("Sn", "Tin",          50, 118.69,      (102, 128, 128), 1.39, 2.17, 1.8),
    "Sb": Element("Sb", "Antimony",     51, 121.75,      (158,  99, 181), 1.39, 2.06, 1.8),
    "Te": Element("Te", "Tellurium",    52, 127.6,       (212, 122,   0), 1.38, 2.06, 1.8),
    "I":  Element("I",  "Iodine",       53, 126.9045,    (148,   0, 148), 1.39, 1.98, 1.8),
    "Xe": Element("Xe", "Xenon",        54, 131.3,       ( 66, 158, 176), 1.40, 2.16, 1.8),
    # Period 6
    "Cs": Element("Cs", "Cesium",       55, 132.9054,    ( 87,  23, 143), 2.44, 3.43, 1.8),
    "Ba": Element("Ba", "Barium",       56, 137.33,      (  0, 201,   0), 2.15, 2.68, 1.8),
    "La": Element("La", "Lanthanum",    57, 138.9055,    (112, 212, 255), 2.07, 2.43, 1.8),
    "Ce": Element("Ce", "Cerium",       58, 140.12,      (255, 255, 199), 2.04, 2.42, 1.8),
    "Pr": Element("Pr", "Praseodymium", 59, 140.9077,    (217, 255, 199), 2.03, 2.40, 1.8),
    "Nd": Element("Nd", "Neodymium",    60, 144.24,      (199, 255, 199), 2.01, 2.39, 1.8),
    "Pm": Element("Pm", "Promethium",   61, 145.0,       (163, 255, 199), 1.99, 2.38, 1.8),
    "Sm": Element("Sm", "Samarium",     62, 150.4,       (143, 255, 199), 1.98, 2.36, 1.8),
    "Eu": Element("Eu", "Europium",     63, 151.96,      ( 97, 255, 199), 1.98, 2.35, 1.8),
    "Gd": Element("Gd", "Gadolinium",   64, 157.25,      ( 69, 255, 199), 1.96, 2.34, 1.8),
    "Tb": Element("Tb", "Terbium",      65, 158.9254,    ( 48, 255, 199), 1.94, 2.33, 1.8),
    "Dy": Element("Dy", "Dysprosium",   66, 162.5,       ( 31, 255, 199), 1.92, 2.31, 1.8),
    "Ho": Element("Ho", "Holmium",      67, 164.9304,    (  0, 255, 156), 1.92, 2.30, 1.8),
    "Er": Element("Er", "Erbium",       68, 167.26,      (  0, 230, 117), 1.89, 2.29, 1.8),
    "Tm": Element("Tm", "Thulium",      69, 168.9342,    (  0, 212,  82), 1.90, 2.27, 1.8),
    "Yb": Element("Yb", "Ytterbium",    70, 173.04,      (  0, 191,  56), 1.87, 2.26, 1.8),
    "Lu": Element("Lu", "Lutetium",     71, 174.967,     (  0, 171,  36), 1.87, 2.24, 1.8),
    "Hf": Element("Hf", "Hafnium",      72, 178.49,      ( 77, 194, 255), 1.75, 2.23, 1.8),
    "Ta": Element("Ta", "Tantalum",      73, 180.9479,    ( 77, 166, 255), 1.70, 2.22, 1.8),
    "W":  Element("W",  "Tungsten",      74, 183.85,      ( 33, 148, 214), 1.62, 2.18, 1.8),
    "Re": Element("Re", "Rhenium",       75, 186.207,     ( 38, 125, 171), 1.51, 2.16, 1.8),
    "Os": Element("Os", "Osmium",        76, 190.2,       ( 38, 102, 150), 1.44, 2.16, 1.8),
    "Ir": Element("Ir", "Iridium",       77, 192.22,      ( 23,  84, 135), 1.41, 2.20, 1.8),
    "Pt": Element("Pt", "Platinum",      78, 195.09,      (208, 208, 224), 1.36, 2.13, 1.8),
    "Au": Element("Au", "Gold",          79, 196.9665,    (255, 209,  35), 1.36, 2.14, 1.8),
    "Hg": Element("Hg", "Mercury",       80, 200.59,      (184, 184, 208), 1.32, 2.09, 1.8),
    "Tl": Element("Tl", "Thallium",      81, 204.37,      (166,  84,  77), 1.45, 1.96, 1.8),
    "Pb": Element("Pb", "Lead",          82, 207.2,       ( 87,  89,  97), 1.46, 2.02, 1.8),
    "Bi": Element("Bi", "Bismuth",       83, 208.9804,    (158,  79, 181), 1.48, 2.07, 1.8),
    "Po": Element("Po", "Polonium",      84, 209.0,       (171,  92,   0), 1.40, 1.97, 1.8),
    "At": Element("At", "Astatine",      85, 210.0,       (117,  79,  69), 1.50, 2.02, 1.8),
    "Rn": Element("Rn", "Radon",         86, 222.0,       ( 66, 130, 150), 1.50, 2.20, 1.8),
    # Period 7
    "Fr": Element("Fr", "Francium",      87, 223.0,       ( 66,   0, 102), 2.60, 3.48, 1.8),
    "Ra": Element("Ra", "Radium",        88, 226.0254,    (  0, 125,   0), 2.21, 2.83, 1.8),
    "Ac": Element("Ac", "Actinium",      89, 227.0278,    (112, 171, 250), 2.15, 2.47, 1.8),
    "Th": Element("Th", "Thorium",       90, 232.0381,    (  0, 186, 255), 2.06, 2.45, 1.8),
    "Pa": Element("Pa", "Protactinium",  91, 231.0359,    (  0, 161, 255), 2.00, 2.43, 1.8),
    "U":  Element("U",  "Uranium",       92, 238.029,     (  0, 143, 255), 1.96, 2.41, 1.8),
    "Np": Element("Np", "Neptunium",     93, 237.0482,    (  0, 128, 255), 1.90, 2.39, 1.8),
    "Pu": Element("Pu", "Plutonium",     94, 242.0,       (  0, 107, 255), 1.87, 2.43, 1.8),
    # Beyond Pu (no CPK colors in moltui — use default gray)
    "Am": Element("Am", "Americium",     95, 243.0,       (180, 180, 180), 1.80, 2.44, 1.8),
    "Cm": Element("Cm", "Curium",        96, 247.0,       (180, 180, 180), 1.69, 2.45, 1.8),
    "Bk": Element("Bk", "Berkelium",     97, 247.0,       (180, 180, 180), 1.70, 2.44, 1.8),
    "Cf": Element("Cf", "Californium",   98, 251.0,       (180, 180, 180), 1.70, 2.45, 1.8),
    "Es": Element("Es", "Einsteinium",   99, 252.0,       (180, 180, 180), 1.70, 2.45, 1.8),
    "Fm": Element("Fm", "Fermium",      100, 257.0,       (180, 180, 180), 1.70, 2.45, 1.8),
    "Md": Element("Md", "Mendelevium",  101, 258.0,       (180, 180, 180), 1.70, 2.45, 1.8),
    "No": Element("No", "Nobelium",     102, 250.0,       (180, 180, 180), 1.70, 2.45, 1.8),
    "Lr": Element("Lr", "Lawrencium",   103, 260.0,       (180, 180, 180), 1.70, 2.45, 1.8),
    "Rf": Element("Rf", "Rutherfordium",104, 261.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Db": Element("Db", "Dubnium",      105, 262.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Sg": Element("Sg", "Seaborgium",   106, 263.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Bh": Element("Bh", "Bohrium",      107, 262.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Hs": Element("Hs", "Hassium",      108, 255.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Mt": Element("Mt", "Meitnerium",   109, 256.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Ds": Element("Ds", "Darmstadtium", 110, 281.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Rg": Element("Rg", "Roentgenium",  111, 282.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Cn": Element("Cn", "Copernicium",  112, 285.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Nh": Element("Nh", "Nihonium",     113, 286.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Fl": Element("Fl", "Flerovium",    114, 289.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Mc": Element("Mc", "Moscovium",    115, 290.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Lv": Element("Lv", "Livermorium",  116, 293.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Ts": Element("Ts", "Tennessine",   117, 294.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    "Og": Element("Og", "Oganesson",    118, 294.0,       (180, 180, 180), 1.70, 2.40, 1.8),
    # Placeholders
    "Da": Element("Da", "Dummy",          0,   0.0,       (255,  20, 147), 1.00, 1.70, 1.8),
    "X":  Element("X",  "Unknown",        0,   0.0,       (255,  20, 147), 1.00, 1.70, 1.8),
}
# fmt: on

DEFAULT_ELEMENT = ELEMENTS["X"]

ATOMIC_NUMBER_TO_SYMBOL: dict[int, str] = {
    e.atomic_number: e.symbol for e in ELEMENTS.values() if e.atomic_number > 0
}


def get_element(symbol: str) -> Element:
    key = symbol.strip().capitalize()
    return ELEMENTS.get(key, DEFAULT_ELEMENT)


def get_element_by_number(atomic_number: int) -> Element:
    sym = ATOMIC_NUMBER_TO_SYMBOL.get(atomic_number)
    if sym:
        return ELEMENTS[sym]
    return DEFAULT_ELEMENT

# ======================================================================
# Module: atom
# ======================================================================
"""Unified Atom class compatible with both pdb_wizard and moltui renderers."""


import numpy as np



class Atom:
    __slots__ = (
        "name", "_element", "x", "charge", "alpha", "epsilon",
        "sigma", "c6", "c8", "c10", "id",
    )

    _element_cache: dict[str, "Element"] = {}

    def __init__(self, x: float, y: float, z: float, name: str) -> None:
        self.name = str(name).strip()
        self.x = np.array([float(x), float(y), float(z)])

        cached = Atom._element_cache.get(self.name)
        if cached is not None:
            self._element = cached
        else:
            element_key = "".join(c for c in self.name[:2] if c.isalpha())
            element_key = element_key.lower().capitalize()
            if element_key not in ELEMENTS:
                element_key = element_key[0] if element_key else ""
            self._element = ELEMENTS.get(element_key, DEFAULT_ELEMENT)
            Atom._element_cache[self.name] = self._element

        self.charge = 0.0
        self.alpha = 0.0
        self.epsilon = 0.0
        self.sigma = 0.0
        self.c6 = 0.0
        self.c8 = 0.0
        self.c10 = 0.0
        self.id = 0

    @property
    def element(self) -> Element:
        return self._element

    @property
    def position(self) -> np.ndarray:
        return self.x

    @property
    def bond_r(self) -> float:
        return self._element.bond_r

    @property
    def vdw(self) -> float:
        return self._element.vdw_radius

    @property
    def mass(self) -> float:
        return self._element.mass

    @property
    def atomic_number(self) -> int:
        return self._element.atomic_number

    def __str__(self) -> str:
        return f"Atom instance {self._element.symbol} {self.id}"

# ======================================================================
# Module: pbc
# ======================================================================
"""Periodic Boundary Conditions (PBC) for crystallographic unit cells."""


import numpy as np


class PBC:
    __slots__ = (
        "a", "b", "c", "alpha", "beta", "gamma",
        "volume", "inverse_volume",
        "basis_matrix", "reciprocal_basis_matrix",
    )

    def __init__(
        self, a: float, b: float, c: float,
        alpha: float, beta: float, gamma: float,
    ) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self._compute_matrices()

    def _compute_matrices(self) -> None:
        a, b, c = self.a, self.b, self.c
        alpha, beta, gamma = self.alpha, self.beta, self.gamma

        b00 = a
        b01 = 0.0
        b02 = 0.0
        b10 = b * cos(pi / 180.0 * gamma)
        b11 = b * sin(pi / 180.0 * gamma)
        b12 = 0.0
        b20 = c * cos(pi / 180.0 * beta)
        b21 = ((b * c * cos(pi / 180.0 * alpha)) - (b10 * b20)) / b11
        b22 = sqrt(max(0.0, c * c - b20 * b20 - b21 * b21))

        self.basis_matrix = array([
            [b00, b01, b02],
            [b10, b11, b12],
            [b20, b21, b22],
        ])

        vol = b00 * (b11 * b22 - b12 * b21)
        vol += b01 * (b12 * b20 - b10 * b22)
        vol += b02 * (b10 * b21 - b11 * b20)
        self.volume = vol
        self.inverse_volume = 1.0 / vol if abs(vol) > 1e-30 else 0.0

        iv = self.inverse_volume
        self.reciprocal_basis_matrix = array([
            [iv * (b11 * b22 - b12 * b21), iv * (b02 * b21 - b01 * b22), iv * (b01 * b12 - b02 * b11)],
            [iv * (b12 * b20 - b10 * b22), iv * (b00 * b22 - b02 * b20), iv * (b02 * b10 - b00 * b12)],
            [iv * (b10 * b21 - b11 * b20), iv * (b01 * b20 - b00 * b21), iv * (b00 * b11 - b01 * b10)],
        ])

    def update(
        self, a: float, b: float, c: float,
        alpha: float, beta: float, gamma: float,
    ) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self._compute_matrices()

    def min_image(self, dx: np.ndarray) -> float:
        img = np.matmul(dx, self.reciprocal_basis_matrix)
        img = np.round(img)
        di = np.matmul(img, self.basis_matrix)
        dx_return = dx - di
        return float(np.sqrt(np.dot(dx_return, dx_return)))

    def wrap(self, dx: np.ndarray) -> np.ndarray:
        img = np.matmul(dx, self.reciprocal_basis_matrix)
        img = np.round(img)
        di = np.matmul(img, self.basis_matrix)
        return dx - di

    def wrap_forward(self, dx: np.ndarray) -> np.ndarray:
        img = np.matmul(dx, self.reciprocal_basis_matrix)
        img = floor(img)
        di = np.matmul(img, self.basis_matrix)
        return dx - di

    def get_all_rs_min_image(
        self, system: list,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        i_idx, j_idx = np.triu_indices(len(system), k=1)
        coords = np.array([atom.x for atom in system])
        dx = coords[i_idx] - coords[j_idx]
        frac = np.matmul(dx, self.reciprocal_basis_matrix)
        wrapped = dx - np.matmul(np.round(frac), self.basis_matrix)
        rs = np.linalg.norm(wrapped, axis=1)
        return i_idx, j_idx, rs

    def corners(self) -> np.ndarray:
        b = self.basis_matrix
        verts = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ], dtype=float)
        return verts @ b

    def edges(self) -> list[tuple[int, int]]:
        return [
            (0, 1), (0, 2), (0, 3),
            (1, 4), (1, 5),
            (2, 4), (2, 6),
            (3, 5), (3, 6),
            (4, 7), (5, 7), (6, 7),
        ]

# ======================================================================
# Module: forcefields
# ======================================================================
"""Force field parameter data (OPLS-AA/UFF and PHAHST).

UFF vdW parameters from: Rappe et al., JACS 114, 10024-10035 (1992), Table I.
sigma = x_i / 2^(1/6), epsilon = D_i / k_B  (epsilon in K, sigma in Angstroms).
H, C, N use OPLS-AA values; all others are pure UFF.

Polarizabilities (alpha, Angstrom^3):
  - Non-metals: van Duijnen & Swart, JPCA 102, 2399 (1998).
  - Metal ions: Shannon, J. Appl. Phys. 73, 348-366 (1993), Table I.
    Dielectric polarizabilities alpha_D for common MOF oxidation states.
"""



def _uff(x_i: float, d_i: float, alpha: float = 0.0) -> dict[str, float]:
    """Convert UFF r_min (x_i, Ang) and well depth (d_i, kcal/mol) to sigma/epsilon."""
    sigma = x_i / 1.12246204830937  # 2^(1/6)
    epsilon = d_i / 0.0019872041     # kB in kcal/(mol*K)
    return {"alpha": alpha, "sigma": sigma, "epsilon": epsilon, "c6": 0.0, "c8": 0.0, "c10": 0.0}


def get_forcefield(name: int) -> dict[str, dict[str, float]]:
    ffs: list[dict[str, dict[str, float]]] = []

    # --- OPLS-AA/UFF hybrid ---
    # H, C, N from OPLS-AA; everything else from UFF (Rappe 1992)
    opls_aa_uff: dict[str, dict[str, float]] = {
        # OPLS-AA overrides
        "H":  {"alpha": 0.41380, "sigma": 2.42,    "epsilon": 15.11,    "c6": 0.0, "c8": 0.0, "c10": 0.0},
        "C":  {"alpha": 1.28660, "sigma": 3.55,    "epsilon": 35.25,    "c6": 0.0, "c8": 0.0, "c10": 0.0},
        "N":  {"alpha": 0.97157, "sigma": 3.25,    "epsilon": 85.60,    "c6": 0.0, "c8": 0.0, "c10": 0.0},
        # UFF — Period 1-2
        # alpha sources: non-metals = van Duijnen & Swart 1998; metals = Shannon 1993
        "He": _uff(2.362, 0.056),
        "Li": _uff(2.451, 0.025, alpha=0.029),    # Shannon Li+
        "Be": _uff(2.745, 0.085, alpha=0.008),    # Shannon Be2+
        "B":  _uff(4.083, 0.180, alpha=0.003),    # Shannon B3+
        "O":  _uff(3.500, 0.060, alpha=0.852),    # van Duijnen
        "F":  _uff(3.364, 0.050, alpha=0.444747), # van Duijnen
        "Ne": _uff(3.243, 0.042),
        # Period 3
        "Na": _uff(2.983, 0.030, alpha=0.179),    # Shannon Na+
        "Mg": _uff(3.021, 0.111, alpha=1.32),     # Shannon Mg2+
        "Al": _uff(4.499, 0.505, alpha=0.79),     # Shannon Al3+
        "Si": _uff(4.295, 0.402, alpha=0.87),     # Shannon Si4+
        "P":  _uff(4.147, 0.305, alpha=0.021),    # Shannon P5+
        "S":  _uff(4.035, 0.274, alpha=2.474448), # van Duijnen
        "Cl": _uff(3.947, 0.227, alpha=2.40028),  # van Duijnen
        "Ar": _uff(3.868, 0.185),
        # Period 4
        "K":  _uff(3.812, 0.035, alpha=0.83),     # Shannon K+
        "Ca": _uff(3.399, 0.238, alpha=3.16),     # Shannon Ca2+
        "Sc": _uff(3.295, 0.019, alpha=2.81),     # Shannon Sc3+
        "Ti": _uff(3.175, 0.017, alpha=2.93),     # Shannon Ti4+
        "V":  _uff(3.144, 0.016, alpha=2.92),     # Shannon V5+
        "Cr": _uff(3.023, 0.015, alpha=1.45),     # Shannon Cr3+
        "Mn": _uff(2.961, 0.013, alpha=1.64),     # Shannon Mn2+
        "Fe": _uff(2.912, 0.013, alpha=2.29),     # Shannon Fe3+
        "Co": _uff(2.872, 0.014, alpha=1.65),     # Shannon Co2+
        "Ni": _uff(2.834, 0.015, alpha=1.23),     # Shannon Ni2+
        "Cu": _uff(3.495, 0.005, alpha=2.11),     # Shannon Cu2+
        "Zn": _uff(2.763, 0.124, alpha=2.04),     # Shannon Zn2+
        "Ga": _uff(4.383, 0.415, alpha=1.50),     # Shannon Ga3+
        "Ge": _uff(4.280, 0.379, alpha=1.63),     # Shannon Ge4+
        "As": _uff(4.230, 0.309, alpha=0.47),     # Shannon As5+
        "Se": _uff(4.205, 0.291, alpha=3.29),     # Shannon Se4+
        "Br": _uff(4.189, 0.251, alpha=3.493),    # van Duijnen
        "Kr": _uff(4.141, 0.220),
        # Period 5
        "Rb": _uff(4.114, 0.040, alpha=1.40),     # Shannon Rb+
        "Sr": _uff(3.641, 0.235, alpha=4.24),     # Shannon Sr2+
        "Y":  _uff(3.345, 0.072, alpha=3.81),     # Shannon Y3+
        "Zr": _uff(3.124, 0.069, alpha=3.25),     # Shannon Zr4+
        "Nb": _uff(3.165, 0.059, alpha=3.97),     # Shannon Nb5+
        "Mo": _uff(3.052, 0.056, alpha=3.28),     # Shannon Mo6+
        "Tc": _uff(2.998, 0.048),
        "Ru": _uff(2.963, 0.056, alpha=1.82),     # Shannon Ru4+
        "Rh": _uff(2.929, 0.053, alpha=1.72),     # Shannon Rh3+
        "Pd": _uff(2.899, 0.048, alpha=1.68),     # Shannon Pd2+
        "Ag": _uff(3.148, 0.036, alpha=2.25),     # Shannon Ag+
        "Cd": _uff(2.848, 0.228, alpha=3.40),     # Shannon Cd2+
        "In": _uff(4.463, 0.599, alpha=2.62),     # Shannon In3+
        "Sn": _uff(4.392, 0.567, alpha=2.83),     # Shannon Sn4+
        "Sb": _uff(4.420, 0.449, alpha=1.64),     # Shannon Sb5+
        "Te": _uff(4.470, 0.398, alpha=4.23),     # Shannon Te4+
        "I":  _uff(4.500, 0.339),
        "Xe": _uff(4.404, 0.332),
        # Period 6
        "Cs": _uff(4.517, 0.045, alpha=2.42),     # Shannon Cs+
        "Ba": _uff(3.703, 0.364, alpha=6.40),     # Shannon Ba2+
        "La": _uff(3.522, 0.017, alpha=6.07),     # Shannon La3+
        "Ce": _uff(3.556, 0.013, alpha=6.15),     # Shannon Ce3+
        "Pr": _uff(3.606, 0.010, alpha=5.32),     # Shannon Pr3+
        "Nd": _uff(3.575, 0.010, alpha=5.01),     # Shannon Nd3+
        "Pm": _uff(3.547, 0.009),
        "Sm": _uff(3.520, 0.008, alpha=4.74),     # Shannon Sm3+
        "Eu": _uff(3.493, 0.008, alpha=4.53),     # Shannon Eu3+
        "Gd": _uff(3.368, 0.009, alpha=4.37),     # Shannon Gd3+
        "Tb": _uff(3.451, 0.007, alpha=4.25),     # Shannon Tb3+
        "Dy": _uff(3.428, 0.007, alpha=4.07),     # Shannon Dy3+
        "Ho": _uff(3.409, 0.007, alpha=3.97),     # Shannon Ho3+
        "Er": _uff(3.391, 0.007, alpha=3.81),     # Shannon Er3+
        "Tm": _uff(3.374, 0.006, alpha=3.67),     # Shannon Tm3+
        "Yb": _uff(3.355, 0.228, alpha=3.58),     # Shannon Yb3+
        "Lu": _uff(3.640, 0.041, alpha=3.64),     # Shannon Lu3+
        "Hf": _uff(3.141, 0.072, alpha=3.25),     # Shannon Hf4+
        "Ta": _uff(3.170, 0.081, alpha=4.73),     # Shannon Ta5+
        "W":  _uff(3.069, 0.067, alpha=3.78),     # Shannon W6+
        "Re": _uff(2.954, 0.066, alpha=3.97),     # Shannon Re7+
        "Os": _uff(3.120, 0.037, alpha=2.33),     # Shannon Os4+
        "Ir": _uff(2.840, 0.073, alpha=2.11),     # Shannon Ir4+
        "Pt": _uff(2.754, 0.080, alpha=2.00),     # Shannon Pt4+
        "Au": _uff(3.293, 0.039, alpha=2.47),     # Shannon Au+
        "Hg": _uff(2.705, 0.385),
        "Tl": _uff(4.347, 0.680, alpha=7.28),     # Shannon Tl+
        "Pb": _uff(4.297, 0.663, alpha=6.58),     # Shannon Pb2+
        "Bi": _uff(4.370, 0.518, alpha=6.12),     # Shannon Bi3+
        "Po": _uff(4.709, 0.325),
        "At": _uff(4.750, 0.284),
        "Rn": _uff(4.765, 0.248),
        # Period 7
        "Fr": _uff(4.900, 0.050),
        "Ra": _uff(3.677, 0.404),
        "Ac": _uff(3.478, 0.033),
        "Th": _uff(3.396, 0.026, alpha=4.92),     # Shannon Th4+
        "Pa": _uff(3.424, 0.022),
        "U":  _uff(3.395, 0.022, alpha=4.45),     # Shannon U4+
        "Np": _uff(3.424, 0.019),
        "Pu": _uff(3.424, 0.016),
        "Am": _uff(3.381, 0.014),
        "Cm": _uff(3.326, 0.013),
        "Bk": _uff(3.339, 0.013),
        "Cf": _uff(3.313, 0.013),
        "Es": _uff(3.299, 0.012),
        "Fm": _uff(3.286, 0.012),
        "Md": _uff(3.274, 0.011),
        "No": _uff(3.248, 0.011),
        "Lr": _uff(3.236, 0.011),
    }
    ffs.append(opls_aa_uff)

    # --- PHAHST ---
    phahst: dict[str, dict[str, float]] = {
        "Cu": {"alpha": 0.29252, "sigma": 2.73851, "epsilon": 8.82345, "c6": 6.96956,  "c8": 262.82938, "c10": 13951.49740},
        "C":  {"alpha": 0.71317, "sigma": 3.35929, "epsilon": 4.00147, "c6": 11.88969, "c8": 547.51694, "c10": 27317.97855},
        "O":  {"alpha": 1.68064, "sigma": 3.23867, "epsilon": 3.89544, "c6": 27.70093, "c8": 709.36452, "c10": 19820.89339},
        "H":  {"alpha": 0.02117, "sigma": 1.87446, "epsilon": 3.63874, "c6": 0.16278,  "c8": 5.03239,   "c10": 202.99322},
    }
    ffs.append(phahst)

    return ffs[name]


def apply_ff_to_system(system: list, ff: dict[str, dict[str, float]]) -> list:
    for atom in system:
        el = atom.element.symbol
        if el in ff:
            params = ff[el]
            atom.alpha = params["alpha"]
            atom.sigma = params["sigma"]
            atom.epsilon = params["epsilon"]
            atom.c6 = params["c6"]
            atom.c8 = params["c8"]
            atom.c10 = params["c10"]
    return system

# ======================================================================
# Module: molecule
# ======================================================================
"""Molecule class combining moltui geometry methods with pdb_wizard PBC-aware logic."""



import numpy as np



@dataclass
class Molecule:
    atoms: list[Atom]
    bonds: list[tuple[int, int]] = field(default_factory=list)
    pbc: Optional[PBC] = None

    def center(self) -> np.ndarray:
        if not self.atoms:
            return np.zeros(3)
        positions = np.array([a.x for a in self.atoms])
        return positions.mean(axis=0)

    def radius(self) -> float:
        if not self.atoms:
            return 1.0
        positions = np.array([a.x for a in self.atoms])
        centroid = positions.mean(axis=0)
        distances = np.linalg.norm(positions - centroid, axis=1)
        return float(distances.max()) if len(distances) > 0 else 1.0

    def _adjacency(self) -> dict[int, list[int]]:
        adj: dict[int, list[int]] = {i: [] for i in range(len(self.atoms))}
        for i, j in self.bonds:
            adj[i].append(j)
            adj[j].append(i)
        return adj

    def get_bond_lengths(self) -> list[tuple[int, int, float]]:
        results = []
        for i, j in self.bonds:
            if self.pbc is not None:
                dx = self.atoms[i].x - self.atoms[j].x
                dist = self.pbc.min_image(dx)
            else:
                dist = float(np.linalg.norm(self.atoms[i].x - self.atoms[j].x))
            results.append((i, j, dist))
        return results

    def get_angles(self) -> list[tuple[int, int, int, float]]:
        adj = self._adjacency()
        results = []
        for j, neighbors in adj.items():
            for ni, i in enumerate(neighbors):
                for k in neighbors[ni + 1:]:
                    if self.pbc is not None:
                        v1 = self.pbc.wrap(self.atoms[i].x - self.atoms[j].x)
                        v2 = self.pbc.wrap(self.atoms[k].x - self.atoms[j].x)
                    else:
                        v1 = self.atoms[i].x - self.atoms[j].x
                        v2 = self.atoms[k].x - self.atoms[j].x
                    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
                    angle = float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
                    results.append((i, j, k, angle))
        return results

    def get_dihedrals(self) -> list[tuple[int, int, int, int, float]]:
        adj = self._adjacency()
        results: list[tuple[int, int, int, int, float]] = []
        seen: set[tuple[int, int, int, int]] = set()
        for j, k in self.bonds:
            for i in adj[j]:
                if i == k:
                    continue
                for l_atom in adj[k]:
                    if l_atom == j or l_atom == i:
                        continue
                    key = (i, j, k, l_atom) if i < l_atom else (l_atom, k, j, i)
                    if key in seen:
                        continue
                    seen.add(key)
                    if self.pbc is not None:
                        b1 = self.pbc.wrap(self.atoms[j].x - self.atoms[i].x)
                        b2 = self.pbc.wrap(self.atoms[k].x - self.atoms[j].x)
                        b3 = self.pbc.wrap(self.atoms[l_atom].x - self.atoms[k].x)
                    else:
                        b1 = self.atoms[j].x - self.atoms[i].x
                        b2 = self.atoms[k].x - self.atoms[j].x
                        b3 = self.atoms[l_atom].x - self.atoms[k].x
                    n1 = np.cross(b1, b2)
                    n2 = np.cross(b2, b3)
                    n1_norm = np.linalg.norm(n1) + 1e-10
                    n2_norm = np.linalg.norm(n2) + 1e-10
                    cos_d = np.dot(n1, n2) / (n1_norm * n2_norm)
                    angle = float(np.degrees(np.arccos(np.clip(cos_d, -1.0, 1.0))))
                    results.append((i, j, k, l_atom, angle))
        return results

    def detect_bonds(self, tolerance: float = 1.3, progress_callback=None) -> None:
        self.bonds = []
        n = len(self.atoms)
        if n < 2:
            return

        if self.pbc is not None:
            # Use all-pairs for small systems (vectorized numpy beats Python loop overhead)
            if n < 1500:
                coords = np.array([a.x for a in self.atoms])
                bond_rs = np.array([a.bond_r for a in self.atoms], dtype=float)
                self._detect_bonds_all_pairs(coords, bond_rs, progress_callback)
                return

            coords = np.array([a.x for a in self.atoms])
            bond_rs = np.array([a.bond_r for a in self.atoms], dtype=float)
            cutoff = float(bond_rs.max())

            recip = self.pbc.reciprocal_basis_matrix
            basis = self.pbc.basis_matrix
            frac = coords @ recip
            frac -= np.floor(frac)

            box = np.array([self.pbc.a, self.pbc.b, self.pbc.c])
            n_cells = np.maximum(1, (box / cutoff).astype(int))
            if n_cells.min() < 3:
                self._detect_bonds_all_pairs(coords, bond_rs, progress_callback)
                return

            cell_idx = np.floor(frac * n_cells).astype(int) % n_cells
            grid: dict[tuple, list[int]] = {}
            for ai in range(n):
                key = (int(cell_idx[ai, 0]), int(cell_idx[ai, 1]), int(cell_idx[ai, 2]))
                grid.setdefault(key, []).append(ai)

            # For each atom, only check atoms in its cell + 26 neighbors
            bonds_set: set[tuple[int, int]] = set()
            for ai in range(n):
                ci = (int(cell_idx[ai, 0]), int(cell_idx[ai, 1]), int(cell_idx[ai, 2]))
                candidates = []
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        for dk in (-1, 0, 1):
                            key = (
                                (ci[0] + di) % n_cells[0],
                                (ci[1] + dj) % n_cells[1],
                                (ci[2] + dk) % n_cells[2],
                            )
                            if key in grid:
                                candidates.extend(grid[key])
                if not candidates:
                    continue
                cands = np.array([c for c in candidates if c > ai], dtype=int)
                if len(cands) == 0:
                    continue
                dx = coords[cands] - coords[ai]
                f = dx @ recip
                f -= np.round(f)
                w = f @ basis
                rs = np.sqrt((w * w).sum(axis=1))
                mixed = 0.5 * (bond_rs[cands] + bond_rs[ai])
                mask = rs < mixed
                for j in cands[mask]:
                    bonds_set.add((ai, int(j)))
                if progress_callback is not None and ai % max(1, n // 20) == 0:
                    progress_callback(ai / n)
            self.bonds = sorted(bonds_set)
            if progress_callback is not None:
                progress_callback(1.0)
        else:
            for i in range(n):
                for j in range(i + 1, n):
                    dist = float(np.linalg.norm(self.atoms[i].x - self.atoms[j].x))
                    max_bond = (
                        self.atoms[i].element.covalent_radius
                        + self.atoms[j].element.covalent_radius
                    ) * tolerance
                    if dist < max_bond:
                        self.bonds.append((i, j))

    def _detect_bonds_all_pairs(self, coords, bond_rs, progress_callback=None) -> None:
        """Fallback all-pairs algorithm for small cells where spatial hashing has too few cells."""
        n = len(coords)
        i_all, j_all = np.triu_indices(n, k=1)
        total_pairs = len(i_all)
        chunk_size = max(500000, total_pairs // 20)
        all_bonds_i = []
        all_bonds_j = []
        for start in range(0, total_pairs, chunk_size):
            end = min(start + chunk_size, total_pairs)
            i_chunk = i_all[start:end]
            j_chunk = j_all[start:end]
            dx = coords[i_chunk] - coords[j_chunk]
            frac = np.matmul(dx, self.pbc.reciprocal_basis_matrix)
            wrapped = dx - np.matmul(np.round(frac), self.pbc.basis_matrix)
            rs = np.linalg.norm(wrapped, axis=1)
            mixed = 0.5 * (bond_rs[i_chunk] + bond_rs[j_chunk])
            mask = rs < mixed
            all_bonds_i.append(i_chunk[mask])
            all_bonds_j.append(j_chunk[mask])
            if progress_callback is not None:
                progress_callback(end / total_pairs)
        if all_bonds_i:
            bi = np.concatenate(all_bonds_i)
            bj = np.concatenate(all_bonds_j)
            self.bonds = list(zip(bi.tolist(), bj.tolist()))

    def find_molecules(self) -> list[Molecule]:
        # Use cached bonds when available — detect_bonds uses an O(N) spatial
        # hash, much faster than _find_molecules_pbc's O(N²) pairwise scan.
        if not self.bonds:
            self.detect_bonds()
        return _find_molecules_from_bonds(self.atoms, self.bonds, self.pbc)


def set_atom_ids(system: list[Atom]) -> None:
    for ind, atom in enumerate(system):
        atom.id = ind + 1


def _find_molecules_pbc(atoms: list[Atom], pbc: PBC) -> list[Molecule]:
    i_idxs, j_idxs, rs = pbc.get_all_rs_min_image(atoms)
    bond_rs = np.array([a.bond_r for a in atoms], dtype=float)
    mixed = 0.5 * (bond_rs[i_idxs] + bond_rs[j_idxs])
    bond_mask = mixed > rs
    edge_list = np.stack((i_idxs[bond_mask], j_idxs[bond_mask]), axis=1) if bond_mask.any() else np.empty((0, 2), dtype=int)

    try:
        import graph_tool.all as gt
        g = gt.Graph(g=len(atoms), directed=False)
        g.add_edge_list(edge_list)
        comp, _ = gt.label_components(g)
        mol_labels = np.array(comp.a)
        n_mols = int(np.max(mol_labels)) + 1 if len(mol_labels) > 0 else 0
        mols = []
        for mol_idx in range(n_mols):
            mol_atoms = [atoms[i] for i in range(len(atoms)) if mol_labels[i] == mol_idx]
            mols.append(Molecule(atoms=mol_atoms, pbc=pbc))
        return mols
    except ImportError:
        pass

    set_atom_ids(atoms)
    bonds: dict[int, list[int]] = {}
    for i_idx, j_idx in edge_list:
        bonds.setdefault(int(i_idx), [])
        if int(j_idx) not in bonds[int(i_idx)]:
            bonds[int(i_idx)].append(int(j_idx))
        bonds.setdefault(int(j_idx), [])
        if int(i_idx) not in bonds[int(j_idx)]:
            bonds[int(j_idx)].append(int(i_idx))

    visited: set[int] = set()
    mols_by_idx: list[list[int]] = []

    for idx in range(len(atoms)):
        if idx in visited:
            continue
        if idx not in bonds:
            mols_by_idx.append([idx])
            visited.add(idx)
            continue
        new_mol = [idx]
        visited.add(idx)
        queue = list(bonds.get(idx, []))
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            new_mol.append(current)
            for neighbor in bonds.get(current, []):
                if neighbor not in visited:
                    queue.append(neighbor)
        mols_by_idx.append(new_mol)

    return [Molecule(atoms=[atoms[i] for i in mol_idxs], pbc=pbc) for mol_idxs in mols_by_idx]


def _find_molecules_from_bonds(
    atoms: list[Atom], bonds: list[tuple[int, int]], pbc: PBC | None
) -> list[Molecule]:
    adj: dict[int, list[int]] = {i: [] for i in range(len(atoms))}
    for i, j in bonds:
        adj[i].append(j)
        adj[j].append(i)

    visited: set[int] = set()
    mols: list[Molecule] = []
    for idx in range(len(atoms)):
        if idx in visited:
            continue
        component = []
        queue = [idx]
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            for neighbor in adj[current]:
                if neighbor not in visited:
                    queue.append(neighbor)
        mol_atoms = [atoms[i] for i in component]
        comp_set = set(component)
        old_to_new = {old: new for new, old in enumerate(component)}
        mol_bonds = [
            (old_to_new[i], old_to_new[j])
            for i, j in bonds if i in comp_set and j in comp_set
        ]
        mols.append(Molecule(atoms=mol_atoms, bonds=mol_bonds, pbc=pbc))
    return mols

# ======================================================================
# Module: io
# ======================================================================
"""File I/O for PDB, XYZ, Zmat, CIF formats — readers, writers, trajectory support."""


import io as _io

import numpy as np


# ---------------------------------------------------------------------------
# Single-frame readers
# ---------------------------------------------------------------------------

def _split_packed_floats(s: str) -> list[str]:
    """Split a string of packed floats that may lack spaces (e.g. '100.123-100.456')."""
    import re
    return re.findall(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?", s)


def _parse_pdb_atom(line: str) -> Atom:
    """Parse an ATOM/HETATM line.

    Handles three PDB variants:
    1. Standard PDB: fixed-width columns (name 12:16, coords 30:38/38:46/46:54)
    2. MPMC PDB: wider space-delimited fields (name tokens[2], coords tokens[6:9])
    3. Overflow PDB: serial merges with record type, coords may be packed

    Detection: if tokens[0] is longer than 6 chars (e.g. 'HETATM99999'), the
    serial overflowed and token positions are shifted — use fixed-width + regex.
    Otherwise use token-based parsing which works for both standard and MPMC.
    """
    tokens = line.split()

    if len(tokens[0]) > 6:
        # Overflow: serial merged with record type, use fixed-width name + regex coords
        name = line[12:16].strip()
        floats = _split_packed_floats(line[30:])
        return Atom(float(floats[0]), float(floats[1]), float(floats[2]), name)

    # Token-based: works for both standard PDB and MPMC format
    # Layout: [record, serial, name, resName, chainID, resSeq, x, y, z, ...]
    # MPMC:   [... x, y, z, mass, charge, alpha, epsilon, sigma, ...]
    name = tokens[2]
    x = float(tokens[6])
    y = float(tokens[7])
    z = float(tokens[8])
    atom = Atom(x, y, z, name)
    # Try to read MPMC extended columns
    if len(tokens) > 10:
        try:
            atom.charge = float(tokens[10])
        except (ValueError, IndexError):
            pass
    if len(tokens) > 13:
        try:
            atom.alpha = float(tokens[11])
            atom.epsilon = float(tokens[12])
            atom.sigma = float(tokens[13])
        except (ValueError, IndexError):
            pass
    return atom


def read_pdb(file: TextIO, filename: str = "") -> tuple[list[Atom], Optional[PBC]]:
    lines = file.readlines()
    pbc = None
    system: list[Atom] = []

    for line in lines:
        if len(line) < 4:
            continue
        # Fast record detection without full split
        if line[0] == "A" and line[1:4] == "TOM":
            try:
                atom = _parse_pdb_atom(line)
                if atom.atomic_number > 0:
                    system.append(atom)
            except (ValueError, IndexError):
                pass
            continue
        if line[0] == "H" and line[1:6] == "ETATM":
            try:
                atom = _parse_pdb_atom(line)
                if atom.atomic_number > 0:
                    system.append(atom)
            except (ValueError, IndexError):
                pass
            continue
        if line[0] == "R" and line[:6] == "REMARK":
            tokens = line.split()
            if len(tokens) >= 8 and tokens[1] == "carbasis":
                try:
                    pbc = PBC(
                        float(tokens[2]), float(tokens[3]), float(tokens[4]),
                        float(tokens[5]), float(tokens[6]), float(tokens[7]),
                    )
                except (ValueError, IndexError):
                    pass
            continue
        if line[0] == "C" and line[:6] == "CRYST1":
            try:
                tokens = line.split()
                pbc = PBC(
                    float(tokens[1]), float(tokens[2]), float(tokens[3]),
                    float(tokens[4]), float(tokens[5]), float(tokens[6]),
                )
            except (ValueError, IndexError):
                pass
            continue
        if line[:3] == "END":
            break

    set_atom_ids(system)
    return system, pbc


def read_xyz(file: TextIO, filename: str = "") -> tuple[list[Atom], Optional[PBC]]:
    first_line = file.readline().strip()
    pbc = None
    system: list[Atom] = []

    # First line is atom count
    try:
        n_atoms = int(first_line)
    except ValueError:
        n_atoms = None

    line = file.readline()
    try:
        tokens = line.split()
        if len(tokens) != 6:
            raise ValueError
        pbc = PBC(
            float(tokens[0]), float(tokens[1]), float(tokens[2]),
            float(tokens[3]), float(tokens[4]), float(tokens[5]),
        )
    except ValueError:
        pass

    try:
        if n_atoms is not None:
            # Read exactly n_atoms lines (handles multi-frame files)
            for _ in range(n_atoms):
                line = file.readline()
                if not line:
                    break
                tokens = line.split()
                atom = Atom(tokens[1], tokens[2], tokens[3], tokens[0])
                try:
                    atom.charge = float(tokens[4])
                except (ValueError, IndexError):
                    pass
                system.append(atom)
        else:
            for line in file.readlines():
                if line == "" or line == "\n":
                    continue
                tokens = line.split()
                atom = Atom(tokens[1], tokens[2], tokens[3], tokens[0])
                try:
                    atom.charge = float(tokens[4])
                except (ValueError, IndexError):
                    pass
                system.append(atom)
    except (ValueError, IndexError):
        pass  # stop reading on malformed line

    set_atom_ids(system)
    return system, pbc


# ---------------------------------------------------------------------------
# Z-matrix (from moltui)
# ---------------------------------------------------------------------------

def _zmat_to_cartesian(
    symbols: list[str],
    refs: list[tuple[int, ...]],
    values: list[tuple[float, ...]],
) -> list[np.ndarray]:
    coords: list[np.ndarray] = []
    for i in range(len(symbols)):
        if i == 0:
            coords.append(np.array([0.0, 0.0, 0.0]))
        elif i == 1:
            r = values[i][0]
            coords.append(np.array([r, 0.0, 0.0]))
        elif i == 2:
            r = values[i][0]
            angle = np.radians(values[i][1])
            ref_a, ref_b = refs[i][0], refs[i][1]
            d = coords[ref_a] - coords[ref_b]
            d_norm = d / (np.linalg.norm(d) + 1e-15)
            if abs(d_norm[1]) < 0.9:
                perp = np.cross(d_norm, np.array([0.0, 1.0, 0.0]))
            else:
                perp = np.cross(d_norm, np.array([1.0, 0.0, 0.0]))
            perp /= np.linalg.norm(perp) + 1e-15
            pos = coords[ref_a] + r * (-d_norm * np.cos(angle) + perp * np.sin(angle))
            coords.append(pos)
        else:
            r = values[i][0]
            angle = np.radians(values[i][1])
            dihedral = np.radians(values[i][2])
            ref_a, ref_b, ref_c = refs[i][0], refs[i][1], refs[i][2]

            ab = coords[ref_b] - coords[ref_a]
            ab /= np.linalg.norm(ab) + 1e-15
            bc = coords[ref_c] - coords[ref_b]

            n = ab
            bc_perp = bc - np.dot(bc, n) * n
            bc_perp_norm = np.linalg.norm(bc_perp)
            if bc_perp_norm < 1e-10:
                if abs(n[1]) < 0.9:
                    d2 = np.cross(n, np.array([0.0, 1.0, 0.0]))
                else:
                    d2 = np.cross(n, np.array([1.0, 0.0, 0.0]))
                d2 /= np.linalg.norm(d2)
            else:
                d2 = bc_perp / bc_perp_norm
            d3 = np.cross(n, d2)

            pos = coords[ref_a] + r * (
                -n * np.cos(angle)
                + d2 * np.sin(angle) * np.cos(dihedral)
                + d3 * np.sin(angle) * np.sin(dihedral)
            )
            coords.append(pos)
    return coords


# ---------------------------------------------------------------------------
# CIF reader
# ---------------------------------------------------------------------------

def _strip_cif_uncertainty(s: str) -> float:
    """Strip parenthesized uncertainty and parse as float: '25.832(10)' -> 25.832."""
    return float(re.sub(r"\([^)]*\)", "", s))


def _tokenize_cif_line(line: str) -> list[str]:
    """Split a CIF line into tokens, respecting single- and double-quoted strings."""
    tokens: list[str] = []
    i = 0
    n = len(line)
    while i < n:
        if line[i] in (" ", "\t"):
            i += 1
            continue
        if line[i] in ("'", '"'):
            quote = line[i]
            i += 1
            start = i
            while i < n and line[i] != quote:
                i += 1
            tokens.append(line[start:i])
            i += 1  # skip closing quote
        else:
            start = i
            while i < n and line[i] not in (" ", "\t"):
                i += 1
            tokens.append(line[start:i])
    return tokens


def _parse_symop(expr: str) -> tuple[np.ndarray, float]:
    """Parse one component of a symmetry operation string.

    E.g. '-x+1/2' -> ([-1, 0, 0], 0.5)
    Returns (coefficients_xyz, translation).
    """
    coeffs = np.zeros(3)
    translation = 0.0
    expr = expr.strip().replace(" ", "")

    i = 0
    n = len(expr)
    while i < n:
        # Determine sign
        sign = 1.0
        if expr[i] == "+":
            sign = 1.0
            i += 1
        elif expr[i] == "-":
            sign = -1.0
            i += 1

        if i >= n:
            break

        # Check for x, y, z variable
        if expr[i] in ("x", "y", "z"):
            idx = {"x": 0, "y": 1, "z": 2}[expr[i]]
            coeffs[idx] = sign
            i += 1
        # Check for fraction or integer
        elif expr[i].isdigit() or expr[i] == ".":
            # Read the number (could be fraction like 1/2 or decimal)
            num_start = i
            while i < n and expr[i] not in ("+", "-", "x", "y", "z"):
                i += 1
            num_str = expr[num_start:i]
            if "/" in num_str:
                num, den = num_str.split("/")
                translation += sign * float(num) / float(den)
            else:
                translation += sign * float(num_str)
        else:
            i += 1  # skip unexpected character

    return coeffs, translation


def _parse_symop_full(op_string: str) -> tuple[np.ndarray, np.ndarray]:
    """Parse a full symmetry operation string like 'x, y, z' or '-x+1/2, y, -z+1/2'.

    Returns (rotation_matrix (3x3), translation_vector (3,)).
    """
    parts = op_string.split(",")
    if len(parts) != 3:
        raise ValueError(f"Expected 3 components in symmetry op, got {len(parts)}: {op_string!r}")

    rotation = np.zeros((3, 3))
    translation = np.zeros(3)
    for i, part in enumerate(parts):
        coeffs, trans = _parse_symop(part)
        rotation[i] = coeffs
        translation[i] = trans

    return rotation, translation


def _apply_symops(
    frac_coords: np.ndarray,
    labels: list[str],
    symbols: list[str],
    symops: list[tuple[np.ndarray, np.ndarray]],
    tolerance: float = 0.01,
) -> tuple[np.ndarray, list[str], list[str]]:
    """Apply symmetry operations to generate all atoms in the unit cell.

    Args:
        frac_coords: (N, 3) fractional coordinates of asymmetric unit.
        labels: atom labels for each atom.
        symbols: element symbols for each atom.
        symops: list of (rotation, translation) from _parse_symop_full.
        tolerance: fractional distance below which two atoms are duplicates.

    Returns:
        (all_frac, all_labels, all_symbols) with duplicates removed.
    """
    all_frac: list[np.ndarray] = []
    all_labels: list[str] = []
    all_symbols: list[str] = []

    for rot, trans in symops:
        for j in range(len(frac_coords)):
            new_frac = rot @ frac_coords[j] + trans
            # Wrap into [0, 1)
            new_frac = new_frac % 1.0

            # Check for duplicates
            is_duplicate = False
            for existing in all_frac:
                diff = new_frac - existing
                # Handle periodic boundary: distance across 0/1 boundary
                diff = diff - np.round(diff)
                if np.linalg.norm(diff) < tolerance:
                    is_duplicate = True
                    break
            if not is_duplicate:
                all_frac.append(new_frac)
                all_labels.append(labels[j])
                all_symbols.append(symbols[j])

    return np.array(all_frac), all_labels, all_symbols


def _frac_to_cart(frac_coords: np.ndarray, basis_matrix: np.ndarray) -> np.ndarray:
    """Convert fractional coordinates to Cartesian using the PBC basis matrix.

    The basis_matrix rows are the cell vectors: frac @ basis_matrix = cart.
    """
    return frac_coords @ basis_matrix


def read_cif(file: TextIO, filename: str = "") -> tuple[list[Atom], Optional[PBC]]:
    """Read a CIF file, returning atoms (in Cartesian) and PBC.

    Handles:
    - P1 CIF files with fractional coordinates (most common for MOFs)
    - CIF files with symmetry operations that need expansion
    - Both _atom_site_fract_x and _atom_site_Cartn_x
    """
    lines = file.readlines()

    # --- Phase 1: Parse all key-value pairs and loops ---
    cell_params: dict[str, float] = {}
    symops: list[str] = []
    atom_loop_tags: list[str] = []
    atom_loop_rows: list[list[str]] = []
    _space_group: str = ""

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i].strip()

        # Skip comments and blank lines
        if not line or line.startswith("#"):
            i += 1
            continue

        # Skip data_ block headers
        if line.startswith("data_"):
            i += 1
            continue

        # Handle semicolon text fields (skip them)
        if line.startswith(";"):
            i += 1
            while i < n and not lines[i].strip().startswith(";"):
                i += 1
            i += 1  # skip closing semicolon
            continue

        # Handle loops
        if line.lower() == "loop_":
            i += 1
            loop_tags: list[str] = []
            loop_rows: list[list[str]] = []

            # Read tags
            while i < n:
                tag_line = lines[i].strip()
                if not tag_line or tag_line.startswith("#"):
                    i += 1
                    continue
                if tag_line.startswith("_"):
                    loop_tags.append(tag_line.lower())
                    i += 1
                else:
                    break

            # Read data rows
            while i < n:
                data_line = lines[i].strip()
                if not data_line or data_line.startswith("#"):
                    i += 1
                    continue
                if (data_line.startswith("_") or data_line.lower() == "loop_"
                        or data_line.startswith("data_")):
                    break
                # Handle semicolon text fields in loop data
                if data_line.startswith(";"):
                    # Collect multi-line text as single value
                    text_val = ""
                    i += 1
                    while i < n and not lines[i].strip().startswith(";"):
                        text_val += lines[i].strip() + " "
                        i += 1
                    i += 1  # skip closing semicolon
                    # This text value belongs to the current row context, but for
                    # atom loops this is rare; skip for simplicity
                    continue
                tokens = _tokenize_cif_line(data_line)
                if tokens:
                    loop_rows.append(tokens)
                i += 1

            # Identify what this loop contains
            symop_tag = None
            for tag in loop_tags:
                if tag in (
                    "_symmetry_equiv_pos_as_xyz",
                    "_space_group_symop_operation_xyz",
                ):
                    symop_tag = tag
                    break

            if symop_tag is not None:
                tag_idx = loop_tags.index(symop_tag)
                for row in loop_rows:
                    if tag_idx < len(row):
                        symops.append(row[tag_idx])

            atom_tag = None
            for tag in loop_tags:
                if tag in ("_atom_site_label", "_atom_site_type_symbol"):
                    atom_tag = tag
                    break

            if atom_tag is not None:
                atom_loop_tags = loop_tags
                atom_loop_rows = loop_rows

            continue

        # Handle key-value pairs
        tokens = _tokenize_cif_line(line)
        if len(tokens) >= 2 and tokens[0].startswith("_"):
            tag = tokens[0].lower()
            value = tokens[1]

            cell_tag_map = {
                "_cell_length_a": "a",
                "_cell_length_b": "b",
                "_cell_length_c": "c",
                "_cell_angle_alpha": "alpha",
                "_cell_angle_beta": "beta",
                "_cell_angle_gamma": "gamma",
            }
            if tag in cell_tag_map:
                try:
                    cell_params[cell_tag_map[tag]] = _strip_cif_uncertainty(value)
                except ValueError:
                    pass

            if tag in ("_symmetry_space_group_name_h-m", "_space_group_name_h-m_alt"):
                _space_group = value.strip()  # noqa: F841

        i += 1

    # --- Phase 2: Build PBC from cell parameters ---
    pbc: Optional[PBC] = None
    required_cell = {"a", "b", "c", "alpha", "beta", "gamma"}
    if required_cell.issubset(cell_params):
        pbc = PBC(
            cell_params["a"], cell_params["b"], cell_params["c"],
            cell_params["alpha"], cell_params["beta"], cell_params["gamma"],
        )

    # --- Phase 3: Parse symmetry operations ---
    if not symops:
        # Default: identity only (P1)
        symops = ["x,y,z"]

    parsed_symops = []
    for op_str in symops:
        try:
            parsed_symops.append(_parse_symop_full(op_str))
        except ValueError:
            continue

    if not parsed_symops:
        parsed_symops = [_parse_symop_full("x,y,z")]

    is_p1 = len(parsed_symops) == 1

    # --- Phase 4: Extract atom data ---
    if not atom_loop_tags or not atom_loop_rows:
        return [], pbc

    # Build column index map
    col: dict[str, int] = {}
    for idx, tag in enumerate(atom_loop_tags):
        col[tag] = idx

    has_frac = all(
        k in col for k in ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")
    )
    has_cart = all(
        k in col for k in ("_atom_site_cartn_x", "_atom_site_cartn_y", "_atom_site_cartn_z")
    )

    labels: list[str] = []
    symbols: list[str] = []
    frac_coords: list[np.ndarray] = []
    cart_coords: list[np.ndarray] = []

    for row in atom_loop_rows:
        # Get label and symbol
        label = ""
        symbol = ""
        if "_atom_site_label" in col and col["_atom_site_label"] < len(row):
            label = row[col["_atom_site_label"]]
        if "_atom_site_type_symbol" in col and col["_atom_site_type_symbol"] < len(row):
            symbol = row[col["_atom_site_type_symbol"]]

        # If no type_symbol, derive from label
        if not symbol and label:
            symbol = "".join(c for c in label if c.isalpha())[:2]
            # Normalize: capitalize first, lower rest
            if len(symbol) > 1:
                symbol = symbol[0].upper() + symbol[1].lower()
            else:
                symbol = symbol.upper()
        if not label and symbol:
            label = symbol

        if not label and not symbol:
            continue

        labels.append(label)
        symbols.append(symbol)

        if has_frac:
            try:
                fx = _strip_cif_uncertainty(row[col["_atom_site_fract_x"]])
                fy = _strip_cif_uncertainty(row[col["_atom_site_fract_y"]])
                fz = _strip_cif_uncertainty(row[col["_atom_site_fract_z"]])
                frac_coords.append(np.array([fx, fy, fz]))
            except (ValueError, IndexError):
                frac_coords.append(np.array([0.0, 0.0, 0.0]))
        if has_cart:
            try:
                cx = _strip_cif_uncertainty(row[col["_atom_site_cartn_x"]])
                cy = _strip_cif_uncertainty(row[col["_atom_site_cartn_y"]])
                cz = _strip_cif_uncertainty(row[col["_atom_site_cartn_z"]])
                cart_coords.append(np.array([cx, cy, cz]))
            except (ValueError, IndexError):
                cart_coords.append(np.array([0.0, 0.0, 0.0]))

    if not labels:
        return [], pbc

    # --- Phase 5: Apply symmetry and convert to Cartesian ---
    system: list[Atom] = []

    if has_frac and pbc is not None:
        frac_arr = np.array(frac_coords)

        if not is_p1:
            frac_arr, labels, symbols = _apply_symops(
                frac_arr, labels, symbols, parsed_symops
            )

        # Convert fractional to Cartesian
        cart_arr = _frac_to_cart(frac_arr, pbc.basis_matrix)

        for j in range(len(labels)):
            atom = Atom(cart_arr[j, 0], cart_arr[j, 1], cart_arr[j, 2], symbols[j])
            atom.name = labels[j]
            system.append(atom)
    elif has_cart:
        for j in range(len(labels)):
            atom = Atom(
                cart_coords[j][0], cart_coords[j][1], cart_coords[j][2], symbols[j]
            )
            atom.name = labels[j]
            system.append(atom)
    elif has_frac and pbc is None:
        # Fractional coords but no cell — treat as-is (unusual edge case)
        for j in range(len(labels)):
            atom = Atom(
                frac_coords[j][0], frac_coords[j][1], frac_coords[j][2], symbols[j]
            )
            atom.name = labels[j]
            system.append(atom)

    set_atom_ids(system)
    return system, pbc


def read_zmat(filepath: str | Path) -> tuple[list[Atom], None]:
    filepath = Path(filepath)
    with open(filepath) as f:
        text = f.read()

    sections = text.strip().split("\n\n")
    atom_lines = [l.strip() for l in sections[0].strip().splitlines() if l.strip()]

    variables: dict[str, float] = {}
    for section in sections[1:]:
        for line in section.strip().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                name, val = line.split("=", 1)
                variables[name.strip()] = float(val.strip())

    def _resolve(token: str) -> float:
        try:
            return float(token)
        except ValueError:
            if token.startswith("-") and token[1:] in variables:
                return -variables[token[1:]]
            return variables[token]

    symbols: list[str] = []
    refs: list[tuple[int, ...]] = []
    values: list[tuple[float, ...]] = []

    for i, line in enumerate(atom_lines):
        parts = line.split()
        sym = ""
        for ch in parts[0]:
            if ch.isalpha():
                sym += ch
            else:
                break
        symbols.append(sym)

        if i == 0:
            refs.append(())
            values.append(())
        elif i == 1:
            refs.append((int(parts[1]) - 1,))
            values.append((_resolve(parts[2]),))
        elif i == 2:
            refs.append((int(parts[1]) - 1, int(parts[3]) - 1))
            values.append((_resolve(parts[2]), _resolve(parts[4])))
        else:
            refs.append((int(parts[1]) - 1, int(parts[3]) - 1, int(parts[5]) - 1))
            values.append((_resolve(parts[2]), _resolve(parts[4]), _resolve(parts[6])))

    positions = _zmat_to_cartesian(symbols, refs, values)
    system = [Atom(pos[0], pos[1], pos[2], sym) for sym, pos in zip(symbols, positions)]
    set_atom_ids(system)
    return system, None


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def _resolve_out(out: TextIO | str):
    """Return (file_handle, should_close). Accepts a path or open writable handle."""
    if isinstance(out, str):
        return open(out, "w"), True
    return out, False


def write_xyz(mol: Molecule, out: TextIO | str) -> None:
    fh, close = _resolve_out(out)
    try:
        fh.write(str(len(mol.atoms)))
        if mol.pbc is not None:
            pbc = mol.pbc
            fh.write(f"\n{pbc.a} {pbc.b} {pbc.c} {pbc.alpha} {pbc.beta} {pbc.gamma}\n")
        else:
            fh.write("\n\n")
        for atom in mol.atoms:
            fh.write(f"{atom.element.symbol} {atom.x[0]} {atom.x[1]} {atom.x[2]}\n")
    finally:
        if close:
            fh.close()


def write_standard_pdb(mol: Molecule, out: TextIO | str, skip_mols_step: bool = False) -> None:
    import copy


    if mol.pbc is None:
        return

    # Work on a copy so we don't reorder the caller's molecule
    mol = copy.deepcopy(mol)
    pbc = mol.pbc
    if skip_mols_step:
        mols_list = [mol]
    else:
        mol, mols_list = sort_system(mol, output=False)

    fh, close = _resolve_out(out)
    out = fh
    out.write("MODEL        1\n")
    out.write("COMPND    " + " " * 69 + "\n")
    out.write("AUTHOR    GENERATED BY PDB WIZARD\n")
    out.write(
        f"CRYST1  {round(pbc.a, 3):>7}  {round(pbc.b, 3):7}  {round(pbc.c, 3):7} "
        f"{round(pbc.alpha, 2):>6} {round(pbc.beta, 2):>6} {round(pbc.gamma, 2):>6} P 1           1\n"
    )

    atom_id = 1
    lines: list[str] = []
    for idx, submol in enumerate(mols_list):
        mol_name = "UNK"
        mol_elements = [a.element.symbol for a in submol.atoms]

        if len(mol_elements) == 1:
            mol_name = mol_elements[0].upper()

        if not skip_mols_step:
            submol.atoms.sort(key=lambda a: a.element.symbol)
            mol_elements = [a.element.symbol for a in submol.atoms]
            known = {
                ("H", "H", "O"): "HOH",
                ("H", "H"): "H2",
                ("H", "H", "H", "H", "C"): "MET",
                ("N", "N"): "N2",
                ("H", "H", "C", "C"): "ACE",
                ("H", "H", "H", "H", "C", "C"): "ENE",
                ("H", "H", "H", "H", "H", "H", "C", "C"): "ETH",
                ("Zn",): "ZNA",
                ("Cl", "Cl", "Cl", "Cl", "Zn"): "ZNC",
            }
            mol_name = known.get(tuple(mol_elements), mol_name)

        base_atom = submol.atoms[-1]

        # Pre-compute per-atom name suffix without rebuilding other_elements per atom (O(N²) → O(N)).
        # An atom needs a numeric suffix if its element appears more than once in the submol.
        if not skip_mols_step:
            from collections import Counter
            elem_counts = Counter(a.element.symbol for a in submol.atoms)
            elem_running: dict[str, int] = {}
            for atom in submol.atoms:
                if elem_counts[atom.element.symbol] > 1:
                    elem_running.setdefault(atom.element.symbol, 0)
                    elem_running[atom.element.symbol] += 1
                    atom.name = atom.element.symbol + str(elem_running[atom.element.symbol])
                else:
                    atom.name = atom.element.symbol
        else:
            for atom in submol.atoms:
                atom.name = atom.element.symbol

        for atom in submol.atoms:
            atom.id = atom_id
            dx = pbc.wrap(atom.x - base_atom.x)
            atom.x = base_atom.x + dx

            lines.append(
                f"HETATM {atom.id:>4}  {atom.name:<3} {mol_name:>3} A {idx + 1:>4}    "
                f"{round(atom.x[0], 3):>7} {round(atom.x[1], 3):>7} {round(atom.x[2], 3):>7}"
                f"  1.00  0.00          {atom.element.symbol:>2}\n"
            )
            atom_id += 1
    out.write("".join(lines))

    out.write("END\n")
    set_atom_ids(mol.atoms)
    if close:
        out.close()


def write_mpmc_pdb(
    mol: Molecule,
    filename: TextIO | str,
    write_charges: bool = False,
    write_params: bool = False,
    sorbate_lines: list[str] | None = None,
) -> None:
    import copy


    if mol.pbc is None:
        return

    # Work on a copy so we don't reorder the caller's molecule
    mol = copy.deepcopy(mol)
    pbc = mol.pbc
    mol, _ = sort_system(mol)

    # Box molecule ID: 2 if no sorbate, 3 if sorbate present
    box_mol_id = 3 if sorbate_lines else 2

    fh, close_fh = _resolve_out(filename)
    try:
        parts: list[str] = [
            "MODEL        1\n",
            "COMPND    " + " " * 69 + "\n",
            "AUTHOR    GENERATED BY PDB WIZARD\n",
            (f"CRYST1  {round(pbc.a, 3):>7}  {round(pbc.b, 3):7}  {round(pbc.c, 3):7} "
             f"{round(pbc.alpha, 2):>6} {round(pbc.beta, 2):>6} {round(pbc.gamma, 2):>6} P 1           1\n"),
        ]
        # MOF atoms (molecule 1, frozen)
        for idx, atom in enumerate(mol.atoms):
            parts.append(
                f"ATOM {idx + 1:>6} {atom.name:<4} MOF F    1    "
                f"{round(atom.x[0], 3):>7} {round(atom.x[1], 3):>7} {round(atom.x[2], 3):>7}"
                f" {atom.mass:>9.5f} {atom.charge:>9.5f}"
                f" {atom.alpha:>9.5f} {atom.epsilon:>9.5f} {atom.sigma:>9.5f}"
                f" 0.0 0.0 {atom.c6:>8.4} {atom.c8:>10.4} {atom.c10:>10.2}\n"
            )

        next_atom_id = len(mol.atoms) + 1

        # Sorbate atoms (molecule 2, movable) — before box
        if sorbate_lines:
            parts.extend(sl + "\n" for sl in sorbate_lines)
            next_atom_id += len(sorbate_lines)

        # Box vertices (last molecule, frozen)
        corners = pbc.corners()
        for ind, pos in enumerate(corners):
            parts.append(
                f"ATOM {next_atom_id + ind:>6} X    BOX F {box_mol_id:>4}    "
                f"{round(pos[0], 3):>7} {round(pos[1], 3):>7} {round(pos[2], 3):>7} 0.0 0.0 0.0 0.0 0.0\n"
            )
        parts.extend(f"CONECT {next_atom_id + i - 1:>4} {next_atom_id + j - 1:>4}\n"
                     for i, j in pbc.edges())
        parts.extend(f"REMARK BOX BASIS[{r_idx}]  {row[0]:20.14f} {row[1]:20.14f} {row[2]:20.14f}\n"
                     for r_idx, row in enumerate(pbc.basis_matrix))
        parts.append("END\n")
        fh.write("".join(parts))
    finally:
        if close_fh:
            fh.close()


# ---------------------------------------------------------------------------
# Trajectory support
# ---------------------------------------------------------------------------

def check_xyz_trajectory(filename: str) -> bool:
    try:
        with open(filename) as f:
            n_atoms = int(f.readline())
            f.readline()
            for _ in range(n_atoms):
                if f.readline() == "":
                    return False
            line = f.readline()
            if line == "":
                return False
            int(line)
            f.readline()
            for _ in range(n_atoms):
                if f.readline() == "":
                    return False
        return True
    except (ValueError, OSError):
        return False


def check_pdb_trajectory(filename: str) -> bool:
    try:
        with open(filename) as f:
            n_model = 0
            n_remark_step = 0
            for line in f:
                if line[:6] == "MODEL ":
                    n_model += 1
                elif line.startswith("REMARK step="):
                    n_remark_step += 1
                if n_model > 1 or n_remark_step > 1:
                    return True
        return False
    except OSError:
        return False


def read_xyz_trajectory(
    file: TextIO, progress_callback=None,
) -> tuple[list[Molecule], list[Optional[PBC]]]:
    molecules: list[Molecule] = []
    pbcs: list[Optional[PBC]] = []
    default_pbc = None

    # Get file size for progress reporting
    file_size = 0
    if progress_callback:
        pos = file.tell()
        file.seek(0, 2)
        file_size = file.tell()
        file.seek(pos)

    line = file.readline()

    try:
        while line != "" and line != "\n":
            n_atoms = int(line)
            pbc: Optional[PBC] = None
            system: list[Atom] = []
            line = file.readline()
            try:
                tokens = line.split()
                if len(tokens) != 6:
                    raise ValueError
                pbc = PBC(*[float(t) for t in tokens])
            except ValueError:
                if default_pbc is None:
                    default_pbc = PBC(1000000, 1000000, 1000000, 90, 90, 90)
                pbc = copy.deepcopy(default_pbc)

            for _ in range(n_atoms):
                line = file.readline()
                tokens = line.split()
                atom = Atom(tokens[1], tokens[2], tokens[3], tokens[0])
                try:
                    atom.charge = float(tokens[4])
                except (ValueError, IndexError):
                    pass
                system.append(atom)

            set_atom_ids(system)
            molecules.append(Molecule(atoms=system, pbc=pbc))
            pbcs.append(pbc)

            if progress_callback and file_size > 0:
                progress_callback(file.tell() / file_size)

            line = file.readline()
    except (ValueError, IndexError):
        sys.exit("Error reading XYZ trajectory")

    return molecules, pbcs


def _is_frame_boundary(line: str) -> bool:
    """Check if a line marks the start of a new PDB/PQR frame."""
    return line[:6] == "MODEL " or line.startswith("REMARK step=")


def read_pdb_trajectory(file: TextIO, progress_callback=None) -> tuple[list[Molecule], list[Optional[PBC]]]:
    molecules: list[Molecule] = []
    pbcs: list[Optional[PBC]] = []

    all_lines = file.readlines()
    total_lines = len(all_lines)

    frame_lines: list[str] = []
    for li, line in enumerate(all_lines):
        if _is_frame_boundary(line) and len(frame_lines) > 3:
            buf = _io.StringIO("".join(frame_lines))
            system, pbc = read_pdb(buf)
            if system:
                molecules.append(Molecule(atoms=system, pbc=pbc))
                pbcs.append(pbc)
            frame_lines = []
            if progress_callback and li % 2000 == 0:
                progress_callback(li / total_lines)
        frame_lines.append(line)

    # Last frame
    if len(frame_lines) > 3:
        buf = _io.StringIO("".join(frame_lines))
        system, pbc = read_pdb(buf)
        if system:
            molecules.append(Molecule(atoms=system, pbc=pbc))
            pbcs.append(pbc)

    if progress_callback:
        progress_callback(1.0)

    for idx, pbc in enumerate(pbcs):
        if pbc is None:
            if idx == 0:
                pbcs[idx] = PBC(1000000, 1000000, 1000000, 90, 90, 90)
            else:
                pbcs[idx] = copy.copy(pbcs[0])

    return molecules, pbcs


# ---------------------------------------------------------------------------
# Unified entry points
# ---------------------------------------------------------------------------

def detect_filetype(filepath: str) -> str:
    suffix = Path(filepath).suffix.lower()
    name = Path(filepath).name.upper()
    if suffix in (".pdb", ".ent", ".pqr"):
        return "pdb"
    if suffix == ".xyz":
        return "xyz"
    if suffix in (".zmat", ".zmatrix"):
        return "zmat"
    if suffix == ".cif":
        return "cif"
    if suffix in (".vasp",) or name in ("POSCAR", "CONTCAR"):
        return "poscar"
    if suffix in (".lmp", ".lammps", ".data"):
        return "lammps"
    if suffix == ".log":
        return "gaussian_log"
    if suffix in (".com", ".gjf"):
        return "gaussian_com"
    if suffix in (".sdf", ".mol"):
        return "sdf"
    if suffix == ".mol2":
        return "mol2"
    if suffix == ".dcd":
        return "dcd"
    raise ValueError(f"Unsupported file format: {suffix}")


def read_file(filepath: str) -> Molecule:
    ft = detect_filetype(filepath)
    if ft == "pdb":
        with open(filepath) as f:
            system, pbc = read_pdb(f, filepath)
    elif ft == "xyz":
        with open(filepath) as f:
            system, pbc = read_xyz(f, filepath)
    elif ft == "zmat":
        system, pbc = read_zmat(filepath)
    elif ft == "cif":
        with open(filepath) as f:
            system, pbc = read_cif(f, filepath)
    elif ft == "poscar":
        with open(filepath) as f:
            system, pbc = read_poscar(f)
    elif ft == "lammps":
        with open(filepath) as f:
            system, pbc = read_lammps_data(f)
    elif ft == "gaussian_log":
        with open(filepath) as f:
            system, pbc = read_gaussian_log(f)
    elif ft == "gaussian_com":
        with open(filepath) as f:
            system, pbc = read_gaussian_com(f)
    elif ft == "sdf":
        with open(filepath) as f:
            system, pbc = read_sdf(f)
    elif ft == "mol2":
        with open(filepath) as f:
            system, pbc = read_mol2(f)
    else:
        raise ValueError(f"Unknown filetype: {ft}")

    mol = Molecule(atoms=system, pbc=pbc)
    mol.detect_bonds()
    return mol


def read_file_trajectory(filepath: str) -> list[Molecule] | None:
    """Read a trajectory file if it is one, otherwise return None."""
    ft = detect_filetype(filepath)
    is_traj = False
    if ft == "pdb":
        is_traj = check_pdb_trajectory(filepath)
    elif ft == "xyz":
        is_traj = check_xyz_trajectory(filepath)
    elif ft == "dcd":
        is_traj = check_dcd_trajectory(filepath)
    if not is_traj:
        return None

    if ft == "pdb":
        with open(filepath) as f:
            mols, _ = read_pdb_trajectory(f)
    elif ft == "xyz":
        with open(filepath) as f:
            mols, _ = read_xyz_trajectory(f)
    elif ft == "dcd":
        mols, _ = read_dcd_trajectory(filepath)
    else:
        return None

    for mol in mols:
        mol.detect_bonds()
    return mols


# ---------------------------------------------------------------------------
# Charge files (CP2K RESP and raw columns)
# ---------------------------------------------------------------------------

def parse_resp_charges(lines: list[str]) -> list[float]:
    """Extract partial charges from CP2K RESP output lines.

    CP2K writes a couple of header lines, one charge per atom (the charge is the
    last whitespace-separated token on the line), and a trailing
    ``Total charge of the system: ...`` line. Header lines are skipped because
    their last token is not a float; the ``Total`` line is skipped explicitly
    (its last token *is* a float and would otherwise be read as an atom charge).
    """
    charges: list[float] = []
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        words = line.split()
        if words[0].lower() == "total":
            continue
        try:
            charges.append(float(words[-1]))
        except ValueError:
            continue  # header / non-data line
    return charges


def read_charges_file(
    path: str, skip_first: int = 0, skip_last: int = 0,
) -> list[float]:
    """Read partial charges from a file.

    CP2K ``.resp`` files are auto-detected and parsed via
    :func:`parse_resp_charges` (header/footer/``Total`` lines stripped, last
    column taken). For any other file the last whitespace-separated token of
    each non-empty line is taken as the charge, after trimming ``skip_first``
    leading and ``skip_last`` trailing lines.
    """
    with open(path) as f:
        lines = f.readlines()
    if path.lower().endswith(".resp"):
        return parse_resp_charges(lines)
    end = max(0, len(lines) - skip_last)
    used = lines[skip_first:end]
    charges: list[float] = []
    for raw in used:
        line = raw.strip()
        if not line:
            continue
        charges.append(float(line.split()[-1]))
    return charges

# ======================================================================
# Module: dcd
# ======================================================================
"""CHARMM/NAMD/OpenMM DCD trajectory reader.

DCD is a Fortran-style binary format. Each record is bracketed by 4-byte
int32 length markers. Layout:

    Header record (84 bytes payload):
        4 bytes "CORD" magic
        20 × int32 (icntrl):
            icntrl[0]  = NSET    — number of frames
            icntrl[1]  = ISTART  — start timestep
            icntrl[2]  = NSAVC   — steps between frames
            icntrl[3]  = NSTEP   — total simulation steps
            icntrl[7]  = NDEGF   — degrees of freedom
            icntrl[8]  = NFROZEN — frozen atoms
            icntrl[9]  = DELTA   — time step (float32 stored in int slot)
            icntrl[10] = USE_BOX — 1 if unit-cell info is written
            icntrl[19] = VERSION — CHARMM version (24 = unit cell present)
    Title record (4-byte ntitle + ntitle × 80-byte titles)
    N_atoms record (1 × int32)
    Per frame:
        if USE_BOX: 6 × float64 record (a, gamma, b, beta, alpha, c)
        x[N] as N × float32
        y[N] as N × float32
        z[N] as N × float32

DCD does NOT carry element symbols. Our reader optionally pairs the DCD
with a topology file (PDB/XYZ/PSF) of the same atom count to assign
elements. Without a topology, every atom is loaded as carbon and a flag
is set on the returned molecules so the UI can warn the user.
"""



import numpy as np



def _detect_endian(fh: BinaryIO) -> str:
    """First 4 bytes are the record-length prefix of the 84-byte header.
    Read it in both byte orders; whichever yields 84 is the file's endian."""
    raw = fh.read(4)
    fh.seek(0)
    if len(raw) < 4:
        raise ValueError("DCD file too short")
    if struct.unpack("<i", raw)[0] == 84:
        return "<"
    if struct.unpack(">i", raw)[0] == 84:
        return ">"
    raise ValueError(
        f"DCD header record length is not 84 in either endianness "
        f"(read {raw!r}). File may be corrupt or not a DCD."
    )


def _read_record(fh: BinaryIO, endian: str) -> bytes:
    """Read one Fortran record: <int32 length> <payload> <int32 length>.
    Verifies the suffix length matches the prefix."""
    head = fh.read(4)
    if len(head) < 4:
        raise EOFError("Unexpected EOF reading record header")
    n = struct.unpack(endian + "i", head)[0]
    payload = fh.read(n)
    if len(payload) != n:
        raise EOFError(f"Short read: expected {n} bytes, got {len(payload)}")
    tail = fh.read(4)
    if len(tail) < 4:
        raise EOFError("Unexpected EOF reading record trailer")
    nt = struct.unpack(endian + "i", tail)[0]
    if nt != n:
        raise ValueError(
            f"DCD record bracket mismatch: head={n}, tail={nt} "
            f"(file likely corrupt or written in a different endianness)"
        )
    return payload


def _parse_header(payload: bytes, endian: str) -> dict:
    """Header payload: 4-char magic + 20 int32."""
    if len(payload) != 84:
        raise ValueError(f"Header payload should be 84 bytes, got {len(payload)}")
    magic = payload[:4]
    if magic != b"CORD":
        raise ValueError(f"Expected CORD magic, got {magic!r}")
    icntrl = struct.unpack(endian + "20i", payload[4:84])
    return {
        "n_frames": icntrl[0],
        "i_start": icntrl[1],
        "n_savc": icntrl[2],
        "n_step": icntrl[3],
        "n_degf": icntrl[7],
        "n_frozen": icntrl[8],
        "delta_raw": icntrl[9],  # float32 packed into int32 slot in CHARMM ≥24
        "has_box": bool(icntrl[10]),
        "version": icntrl[19],
    }


def _read_topology_elements(topo_path: str, n_expected: int) -> list[str] | None:
    """Read element symbols from a sibling PDB/XYZ topology file.
    Returns None if it doesn't fit (different atom count) or can't be parsed."""
    suffix = Path(topo_path).suffix.lower()
    try:
        if suffix in (".pdb", ".ent", ".pqr"):
            with open(topo_path) as f:
                atoms, _ = read_pdb(f)
        elif suffix == ".xyz":
            with open(topo_path) as f:
                atoms, _ = read_xyz(f)
        else:
            return None
        if len(atoms) != n_expected:
            return None
        return [a.element.symbol for a in atoms]
    except Exception:
        return None


def _find_sibling_topology(dcd_path: str) -> str | None:
    """Locate a topology file next to the DCD.

    Lookup order:
      1. Same stem with PDB/PQR/ENT/XYZ extension (e.g. traj.pdb next to traj.dcd)
      2. Any .pdb / .pqr / .xyz in the same directory — common when the
         topology was written by openmm (e.g. system.pdb or final.pdb
         next to trajectory.dcd).
    """
    p = Path(dcd_path)
    # 1. Same-stem match
    for ext in (".pdb", ".pqr", ".ent", ".xyz"):
        candidate = p.with_suffix(ext)
        if candidate.exists():
            return str(candidate)
    # 2. Any topology in the same directory
    for ext in (".pdb", ".pqr", ".xyz"):
        siblings = sorted(p.parent.glob(f"*{ext}"))
        if siblings:
            return str(siblings[0])
    return None


def read_dcd_trajectory(
    filepath: str,
    topology_file: str | None = None,
    progress_callback=None,
) -> tuple[list[Molecule], list[Optional[PBC]]]:
    """Read a DCD file. Returns (frames, pbcs) — same shape as
    read_pdb_trajectory and read_xyz_trajectory.

    If `topology_file` is None, searches for a sibling PDB/XYZ next to the
    DCD. Without any topology, every atom is loaded as carbon and the
    returned Molecules carry `mol._dcd_topology_missing = True` so the UI
    can warn the user.
    """
    # Resolve element list from topology
    topology_path = topology_file or _find_sibling_topology(filepath)

    with open(filepath, "rb") as fh:
        endian = _detect_endian(fh)
        header_payload = _read_record(fh, endian)
        hdr = _parse_header(header_payload, endian)

        # Title record (skip — purely informational)
        _read_record(fh, endian)

        # N_atoms record (one int32)
        n_atoms_payload = _read_record(fh, endian)
        if len(n_atoms_payload) != 4:
            raise ValueError(
                f"N_atoms record should be 4 bytes, got {len(n_atoms_payload)}"
            )
        n_atoms = struct.unpack(endian + "i", n_atoms_payload)[0]
        if n_atoms <= 0:
            raise ValueError(f"DCD reports {n_atoms} atoms — invalid")

        # Try to fetch element symbols from topology
        elements: list[str]
        topology_missing = False
        if topology_path:
            els = _read_topology_elements(topology_path, n_atoms)
            if els is not None:
                elements = els
            else:
                elements = ["C"] * n_atoms
                topology_missing = True
        else:
            elements = ["C"] * n_atoms
            topology_missing = True

        # Frames
        frames: list[Molecule] = []
        pbcs: list[Optional[PBC]] = []
        coords_size = n_atoms * 4
        for fi in range(hdr["n_frames"]):
            pbc: Optional[PBC] = None
            if hdr["has_box"]:
                box_payload = _read_record(fh, endian)
                if len(box_payload) != 48:
                    raise ValueError(
                        f"Frame {fi}: box record should be 48 bytes, "
                        f"got {len(box_payload)}"
                    )
                # CHARMM/NAMD order: a, gamma, b, beta, alpha, c
                a, gamma, b, beta, alpha, c = struct.unpack(endian + "6d", box_payload)
                # Some versions store cosines instead of angles (-1..1)
                if -1.0 <= alpha <= 1.0 and -1.0 <= beta <= 1.0 and -1.0 <= gamma <= 1.0:
                    alpha = float(np.degrees(np.arccos(alpha)))
                    beta = float(np.degrees(np.arccos(beta)))
                    gamma = float(np.degrees(np.arccos(gamma)))
                if a > 0 and b > 0 and c > 0:
                    pbc = PBC(a, b, c, alpha, beta, gamma)

            x_payload = _read_record(fh, endian)
            if len(x_payload) != coords_size:
                raise ValueError(
                    f"Frame {fi}: x-record should be {coords_size} bytes, "
                    f"got {len(x_payload)}"
                )
            xs = np.frombuffer(x_payload, dtype=endian + "f4")
            ys = np.frombuffer(_read_record(fh, endian), dtype=endian + "f4")
            zs = np.frombuffer(_read_record(fh, endian), dtype=endian + "f4")

            atoms = [
                Atom(float(xs[i]), float(ys[i]), float(zs[i]), elements[i])
                for i in range(n_atoms)
            ]
            set_atom_ids(atoms)
            mol = Molecule(atoms=atoms, pbc=pbc)
            if topology_missing:
                # Tag the molecule so UI/CLI can surface a warning to users
                mol._dcd_topology_missing = True
            frames.append(mol)
            pbcs.append(pbc)

            if progress_callback and hdr["n_frames"] > 0:
                progress_callback((fi + 1) / hdr["n_frames"])

    return frames, pbcs


def check_dcd_trajectory(filepath: str) -> bool:
    """Return True if the file is a parseable DCD with at least one frame.
    Does not load full data — just probes the header."""
    try:
        with open(filepath, "rb") as fh:
            endian = _detect_endian(fh)
            header_payload = _read_record(fh, endian)
            hdr = _parse_header(header_payload, endian)
            return hdr["n_frames"] >= 1
    except Exception:
        return False

# ======================================================================
# Module: geometry
# ======================================================================
"""Geometry analysis functions (bonds, angles, contacts, overlaps)."""


import numpy as np



def overlap_detector(mol: Molecule) -> Molecule:
    if mol.pbc is None:
        return mol
    pbc = mol.pbc
    system = mol.atoms
    set_atom_ids(system)
    messages: list[str] = []
    overlapping = True
    while overlapping:
        overlapping = False
        for atom in system:
            for atom2 in system:
                if atom.id != atom2.id and not overlapping:
                    r = pbc.min_image(atom.x - atom2.x)
                    if r < 0.05:
                        overlapping = True
                        messages.append(
                            f"Deleting overlapping atoms "
                            f"{atom.element.symbol:>3} {atom.id:>5} --- "
                            f"{atom2.element.symbol:>3} {atom2.id:>5}"
                        )
                        system.remove(atom2)
    set_atom_ids(system)
    mol.atoms = system
    return mol


def get_close_contacts(mol: Molecule) -> list[str]:
    if mol.pbc is None:
        return []
    pbc = mol.pbc
    system = mol.atoms
    set_atom_ids(system)
    i_idx, j_idx, rs = pbc.get_all_rs_min_image(system)
    vdw_radii = np.array([a.vdw for a in system], dtype=float)
    bond_radii = np.array([a.bond_r for a in system], dtype=float)
    mixed_vdw = 0.5 * (vdw_radii[i_idx] + vdw_radii[j_idx])
    mixed_bond = 0.5 * (bond_radii[i_idx] + bond_radii[j_idx])
    mask = (mixed_vdw > rs) & (rs > mixed_bond)
    messages: list[str] = []
    for idx in np.where(mask)[0]:
        i, j, r = int(i_idx[idx]), int(j_idx[idx]), rs[idx]
        el_str = f"{system[i].element.symbol}-{system[j].element.symbol}"
        messages.append(f"{el_str:<5} {i+1:>5} {j+1:>5}   r = {np.round(r, 6)}")
    return messages


def get_bonds_list(mol: Molecule) -> list[str]:
    """Return formatted bond list. Uses cached mol.bonds if present, else
    runs detect_bonds — both vectorized."""
    if mol.pbc is None:
        return []
    if not mol.bonds:
        mol.detect_bonds()
    if not mol.bonds:
        return []
    pbc = mol.pbc
    system = mol.atoms
    set_atom_ids(system)

    bond_arr = np.asarray(mol.bonds, dtype=int)
    coords = np.array([a.x for a in system])
    dx = coords[bond_arr[:, 0]] - coords[bond_arr[:, 1]]
    # Min image vectorized
    frac = dx @ pbc.reciprocal_basis_matrix
    frac -= np.round(frac)
    dx_min = frac @ pbc.basis_matrix
    rs = np.linalg.norm(dx_min, axis=1)

    symbols = [a.element.symbol for a in system]
    ids = [a.id for a in system]
    return [
        f"{symbols[i]+'-'+symbols[j]:<5} {ids[i]:>5} {ids[j]:>5}   r = {round(float(r), 6)}"
        for (i, j), r in zip(bond_arr, rs)
    ]


def get_angles_list(mol: Molecule) -> list[str]:
    """Return formatted bond-angle list. Iterates only over bonded triples
    (atom_i — atom_center — atom_k where both edges are in mol.bonds), so
    the cost is O(degree² * N) instead of O(N³)."""
    if mol.pbc is None:
        return []
    if not mol.bonds:
        mol.detect_bonds()
    if not mol.bonds:
        return []
    pbc = mol.pbc
    system = mol.atoms
    set_atom_ids(system)

    # Build adjacency map: center_idx -> list of bonded neighbor indices
    adj: dict[int, list[int]] = {}
    for a, b in mol.bonds:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    coords = np.array([a.x for a in system])
    symbols = [a.element.symbol for a in system]
    ids = [a.id for a in system]
    recip = pbc.reciprocal_basis_matrix
    basis = pbc.basis_matrix

    def min_image(v):
        f = v @ recip
        f -= np.round(f)
        return f @ basis

    messages: list[str] = []
    # Iterate centers; for each center, all unordered pairs of bonded neighbors form an angle
    for c, neighbors in adj.items():
        if len(neighbors) < 2:
            continue
        nb_arr = np.array(neighbors)
        # vectors center -> each neighbor (min-imaged)
        dx_all = coords[nb_arr] - coords[c]
        dx_all = min_image(dx_all)
        norms = np.linalg.norm(dx_all, axis=1)
        for ai in range(len(neighbors)):
            for bi in range(ai + 1, len(neighbors)):
                v1, v2 = dx_all[ai], dx_all[bi]
                n1, n2 = norms[ai], norms[bi]
                if n1 < 1e-9 or n2 < 1e-9:
                    continue
                cos_a = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
                angle = float(np.degrees(np.arccos(cos_a)))
                i, k = neighbors[ai], neighbors[bi]
                el_str = f"{symbols[i]}-{symbols[c]}-{symbols[k]}"
                messages.append(
                    f"{el_str:<7} {ids[i]:>5} {ids[c]:>5} {ids[k]:>5}   "
                    f"angle = {round(angle, 2):>6} "
                    f"r1, r2 = {round(float(n1), 3):>6}, {round(float(n2), 3):>6}"
                )
    return messages


def get_lone_atoms(mol: Molecule) -> list[Atom]:
    if mol.pbc is None:
        return []
    pbc = mol.pbc
    system = mol.atoms
    set_atom_ids(system)
    n = len(system)
    i_idx, j_idx, rs = pbc.get_all_rs_min_image(system)
    vdw_radii = np.array([a.vdw for a in system], dtype=float)
    mixed_vdw = 0.5 * (vdw_radii[i_idx] + vdw_radii[j_idx])
    has_neighbor = np.zeros(n, dtype=bool)
    contact_mask = rs < mixed_vdw
    has_neighbor[i_idx[contact_mask]] = True
    has_neighbor[j_idx[contact_mask]] = True
    return [system[i] for i in range(n) if not has_neighbor[i]]


def delete_lone_atoms(mol: Molecule) -> Molecule:
    lone = set(id(a) for a in get_lone_atoms(mol))
    if not lone:
        return mol
    keep = [i for i, a in enumerate(mol.atoms) if id(a) not in lone]
    keep_set = set(keep)
    old_to_new = {old: new for new, old in enumerate(keep)}
    mol.atoms = [mol.atoms[i] for i in keep]
    mol.bonds = [
        (old_to_new[a], old_to_new[b])
        for a, b in mol.bonds if a in keep_set and b in keep_set
    ]
    set_atom_ids(mol.atoms)
    return mol


def edit_h_dist(mol: Molecule, second_element: str, distance: float) -> Molecule:
    if mol.pbc is None:
        return mol
    pbc = mol.pbc
    system = mol.atoms
    set_atom_ids(system)
    messages: list[str] = []
    for atom in system:
        for atom2 in system:
            if atom2.id > atom.id:
                is_h_pair = (
                    (atom.element.symbol == "H" and atom2.element.symbol == second_element) or
                    (atom2.element.symbol == "H" and atom.element.symbol == second_element)
                )
                if not is_h_pair:
                    continue
                h_atom = atom if atom.element.symbol == "H" else atom2
                other = atom2 if atom.element.symbol == "H" else atom
                dx = h_atom.x - other.x
                r = pbc.min_image(dx)
                bond_r = 0.5 * (atom.bond_r + atom2.bond_r)
                if r < bond_r:
                    el_str = f"{atom.element.symbol}-{atom2.element.symbol}"
                    messages.append(f"{el_str:<5} {atom.id:>5} {atom2.id:>5}")
                    dx = pbc.wrap(dx)
                    dx *= distance / r
                    h_atom.x = other.x + dx
    return mol

# ======================================================================
# Module: operations
# ======================================================================
"""System operations: wrap, sort, extend axis, formula unit, etc."""



import numpy as np



def progressbar(
    it: list, prefix: str = "", size: int = 60, out: TextIO = sys.stdout,
) -> Iterator:
    count = len(it)
    if count == 0:
        return

    def show(j: int) -> None:
        x = int(size * j / count)
        bar = "#" * x + "." * (size - x)
        print(
            f"{prefix}[{bar}] {j}/{count}",
            end="\r", file=out, flush=True,
        )

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    print("", flush=True, file=out)


def gcd_list(values: list[int]) -> int:
    return reduce(gcd, values)


def formula_unit(mol: Molecule) -> dict[str, int]:
    counts: dict[str, int] = {}
    for atom in mol.atoms:
        sym = atom.element.symbol
        counts[sym] = counts.get(sym, 0) + 1
    return counts


def formula_unit_reduced(mol: Molecule) -> dict[str, int]:
    counts = formula_unit(mol)
    if not counts:
        return counts
    g = gcd_list(list(counts.values()))
    return {k: v // g for k, v in counts.items()}


def wrap_atoms(mol: Molecule, forward: bool = False) -> Molecule:
    if mol.pbc is None:
        return mol
    pbc = mol.pbc
    wrap_fn = pbc.wrap_forward if forward else pbc.wrap
    for atom in mol.atoms:
        atom.x = wrap_fn(atom.x)
    return mol


def sort_system(mol: Molecule, output: bool = True) -> tuple[Molecule, list[Molecule]]:
    if mol.pbc is None:
        return mol, [mol]

    submols = mol.find_molecules()

    for submol in submols:
        submol.atoms.sort(key=lambda a: a.element.symbol)

    submols.sort(key=lambda m: m.atoms[0].element.symbol)
    submols.sort(key=lambda m: len(m.atoms), reverse=True)
    submols.sort(key=lambda m: 0 if max(a.mass for a in m.atoms) > 16 else 1)

    all_atoms: list[Atom] = []
    for submol in submols:
        all_atoms.extend(submol.atoms)
    mol.atoms = all_atoms
    set_atom_ids(mol.atoms)
    mol.detect_bonds()
    return mol, submols


def extend_axis(mol: Molecule, axis: int, times: int) -> Molecule:
    if mol.pbc is None:
        return mol
    pbc = mol.pbc
    new_atoms: list[Atom] = []
    for i in range(1, times + 1):
        for atom in mol.atoms:
            new_atom = copy.deepcopy(atom)
            new_atom.x = atom.x + i * pbc.basis_matrix[axis]
            new_atoms.append(new_atom)
    mol.atoms.extend(new_atoms)

    params = [pbc.a, pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma]
    params[axis] *= (times + 1)
    pbc.update(*params)
    set_atom_ids(mol.atoms)
    return mol


def void_volume(mol: Molecule, n_samples: int = 100000, probe_radius: float = 0.0) -> dict[str, float]:
    """Calculate void volume using Monte Carlo sampling.

    Samples random points in the unit cell and checks if each point is
    further than (vdw_radius + probe_radius) from all atoms, using PBC
    min-image distances.

    Args:
        mol: Molecule with PBC.
        n_samples: Number of random sample points.
        probe_radius: Radius of probe molecule (0 = geometric void,
            1.2 = H2-sized probe, 1.4 = N2-sized probe).

    Returns:
        Dict with void_volume (A^3), total_volume (A^3), void_fraction.
    """
    if mol.pbc is None or len(mol.atoms) == 0:
        return {"void_volume": 0.0, "total_volume": 0.0, "void_fraction": 0.0}

    pbc = mol.pbc
    coords = np.array([a.x for a in mol.atoms])
    vdw_radii = np.array([a.vdw for a in mol.atoms])

    # Generate random fractional coordinates
    rng = np.random.default_rng(42)
    frac_points = rng.random((n_samples, 3))
    cart_points = frac_points @ pbc.basis_matrix

    # Loop over atoms with all sample points vectorized.
    # Uses orthorhombic fast path (no matrix multiply) when possible.
    effective_radii_sq = (vdw_radii + probe_radius) ** 2
    box = np.array([pbc.a, pbc.b, pbc.c])
    is_ortho = (abs(pbc.alpha - 90) < 0.01 and abs(pbc.beta - 90) < 0.01
                and abs(pbc.gamma - 90) < 0.01)

    occupied = np.zeros(n_samples, dtype=bool)

    if is_ortho:
        for ai in range(len(coords)):
            dx = cart_points - coords[ai]
            dx -= np.round(dx / box) * box
            dist_sq = (dx * dx).sum(axis=1)
            occupied |= (dist_sq <= effective_radii_sq[ai])
    else:
        recip = pbc.reciprocal_basis_matrix
        basis = pbc.basis_matrix
        for ai in range(len(coords)):
            dx = cart_points - coords[ai]
            frac_dx = dx @ recip.T
            frac_dx -= np.round(frac_dx)
            wrapped = frac_dx @ basis
            dist_sq = (wrapped * wrapped).sum(axis=1)
            occupied |= (dist_sq <= effective_radii_sq[ai])

    void_frac = (~occupied).sum() / n_samples
    total_vol = pbc.volume
    return {
        "void_volume": void_frac * total_vol,
        "total_volume": total_vol,
        "void_fraction": void_frac,
    }


def surface_area(
    mol: Molecule, n_samples: int = 100000, probe_radius: float = 1.4,
) -> dict[str, float]:
    """Compute solvent-accessible surface area (SASA) via Monte Carlo.

    For each atom, random points are sampled on a sphere of radius
    (vdw_radius + probe_radius).  A point contributes to SASA only if it
    is not inside any *other* atom's (vdw + probe) sphere.  PBC
    min-image convention is used for all distance calculations.

    Args:
        mol: Molecule with PBC.
        n_samples: Total number of sample points distributed across atoms
            (weighted by sphere surface area).
        probe_radius: Probe molecule radius in angstroms (default 1.4,
            typical N2 probe for BET-style calculations).

    Returns:
        Dict with ``surface_area`` (A^2) and ``area_per_volume`` (m^2/g,
        using total_mass * 1.66054e-24 for grams).
    """
    if mol.pbc is None or len(mol.atoms) == 0:
        return {"surface_area": 0.0, "area_per_volume": 0.0}

    pbc = mol.pbc
    coords = np.array([a.x for a in mol.atoms])            # (N, 3)
    vdw_radii = np.array([a.vdw for a in mol.atoms])        # (N,)
    radii = vdw_radii + probe_radius                        # (N,)
    n_atoms = len(mol.atoms)

    # Distribute samples per atom proportional to sphere surface area (4*pi*r^2).
    areas = radii ** 2                                      # proportional to 4*pi*r^2
    weights = areas / areas.sum()
    samples_per_atom = np.round(weights * n_samples).astype(int)
    samples_per_atom = np.maximum(samples_per_atom, 1)

    rng = np.random.default_rng(42)
    total_sa = 0.0

    basis = pbc.basis_matrix                                # (3, 3)
    recip = pbc.reciprocal_basis_matrix                     # (3, 3)

    # Pre-filter via cell list (O(N) instead of O(N^2)). Cell size = max pair
    # cutoff so a query only inspects the home cell and its 26 neighbors.
    r_max = float(radii.max())
    cutoff = 2 * r_max  # generous: any j within r_i+r_j ≤ 2*r_max
    # Map each atom to a wrapped fractional cell index
    frac = (coords @ recip) % 1.0
    n_cells = max(1, int(min(pbc.a, pbc.b, pbc.c) // cutoff))
    cell_idx = np.minimum((frac * n_cells).astype(int), n_cells - 1)
    cell_map: dict[tuple, list[int]] = {}
    for ai, (cx, cy, cz) in enumerate(cell_idx.tolist()):
        cell_map.setdefault((cx, cy, cz), []).append(ai)

    neighbors: list[np.ndarray] = [None] * n_atoms
    for i in range(n_atoms):
        cx, cy, cz = int(cell_idx[i, 0]), int(cell_idx[i, 1]), int(cell_idx[i, 2])
        cand: list[int] = []
        for dx_c in (-1, 0, 1):
            for dy_c in (-1, 0, 1):
                for dz_c in (-1, 0, 1):
                    key = ((cx + dx_c) % n_cells,
                           (cy + dy_c) % n_cells,
                           (cz + dz_c) % n_cells)
                    cand.extend(cell_map.get(key, ()))
        # Filter by actual distance (min-image) and exclude self
        if not cand:
            neighbors[i] = np.empty(0, dtype=int)
            continue
        cand_arr = np.array(cand, dtype=int)
        cand_arr = cand_arr[cand_arr != i]
        if cand_arr.size == 0:
            neighbors[i] = cand_arr
            continue
        dx = coords[cand_arr] - coords[i]
        f = dx @ recip
        f -= np.round(f)
        dx = f @ basis
        dists = np.sqrt((dx ** 2).sum(axis=1))
        cutoff_pair = radii[i] + radii[cand_arr]
        keep_mask = dists < cutoff_pair
        neighbors[i] = cand_arr[keep_mask]

    # Process per atom — only against pre-filtered neighbors
    for i in range(n_atoms):
        n_pts = int(samples_per_atom[i])
        r_i = radii[i]
        nb = neighbors[i]

        if len(nb) == 0:
            total_sa += 4.0 * np.pi * r_i ** 2
            continue

        # Random points on the unit sphere (Marsaglia method)
        pts = rng.standard_normal((n_pts, 3))
        pts /= np.linalg.norm(pts, axis=1, keepdims=True)
        pts = pts * r_i + coords[i]                         # (n_pts, 3)

        nb_coords = coords[nb]                              # (M, 3)
        nb_radii = radii[nb]                                # (M,)

        # Vectorised min-image distance: pts (n_pts,3) vs nb_coords (M,3)
        # Process in sub-batches to keep memory bounded.
        chunk = max(1, 200_000_000 // (len(nb) * 3 * 8))  # ~200 MB limit
        accessible = np.ones(n_pts, dtype=bool)
        for start in range(0, n_pts, chunk):
            end = min(start + chunk, n_pts)
            batch = pts[start:end]                          # (B, 3)
            dx = nb_coords[np.newaxis, :, :] - batch[:, np.newaxis, :]  # (B, M, 3)
            frac_dx = dx @ recip
            dx = dx - np.round(frac_dx) @ basis
            dists = np.sqrt((dx ** 2).sum(axis=2))          # (B, M)
            inside = (dists < nb_radii[np.newaxis, :]).any(axis=1)
            accessible[start:end] &= ~inside

        frac_accessible = accessible.sum() / n_pts
        total_sa += frac_accessible * 4.0 * np.pi * r_i ** 2

    total_mass = sum(a.mass for a in mol.atoms)
    mass_grams = total_mass * 1.66054e-24                   # amu -> grams
    # surface_area is in A^2; convert to m^2: 1 A^2 = 1e-20 m^2
    area_per_volume = (total_sa * 1e-20) / mass_grams if mass_grams > 0 else 0.0

    return {
        "surface_area": total_sa,
        "area_per_volume": area_per_volume,
    }


def pore_size_distribution(
    mol: Molecule,
    n_samples: int = 50000,
    n_bins: int = 50,
    max_pore: float = 20.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute pore size distribution via Monte Carlo sampling.

    Random points are generated uniformly in the unit cell.  For each
    point the distance to the nearest atom *surface* (centre distance
    minus vdw_radius) is recorded as the pore radius at that point.
    Results are histogrammed and normalised so that the histogram sums
    to 1.  PBC min-image convention is used.

    Args:
        mol: Molecule with PBC.
        n_samples: Number of random sample points.
        n_bins: Number of histogram bins.
        max_pore: Upper limit of the histogram (angstroms).

    Returns:
        Tuple ``(bin_centers, histogram)`` as numpy arrays.
    """
    bin_edges = np.linspace(0.0, max_pore, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    if mol.pbc is None or len(mol.atoms) == 0:
        return bin_centers, np.zeros(n_bins)

    pbc = mol.pbc
    coords = np.array([a.x for a in mol.atoms])            # (N, 3)
    vdw_radii = np.array([a.vdw for a in mol.atoms])        # (N,)
    n_atoms = len(mol.atoms)

    basis = pbc.basis_matrix                                # (3, 3)
    recip = pbc.reciprocal_basis_matrix                     # (3, 3)

    rng = np.random.default_rng(42)
    frac_points = rng.random((n_samples, 3))
    cart_points = frac_points @ basis                       # (n_samples, 3)

    # Process sample points in chunks to keep memory bounded
    chunk = max(1, 200_000_000 // (n_atoms * 3 * 8))       # ~200 MB limit
    pore_radii = np.empty(n_samples)

    for start in range(0, n_samples, chunk):
        end = min(start + chunk, n_samples)
        batch = cart_points[start:end]                      # (B, 3)
        dx = coords[np.newaxis, :, :] - batch[:, np.newaxis, :]  # (B, N, 3)
        frac_dx = dx @ recip
        dx = dx - np.round(frac_dx) @ basis                # min-image
        dists = np.sqrt((dx ** 2).sum(axis=2))              # (B, N)
        # Pore radius = distance to nearest atom surface
        surface_dists = dists - vdw_radii[np.newaxis, :]    # (B, N)
        pore_radii[start:end] = surface_dists.min(axis=1)

    # Clamp negatives to zero (point inside an atom)
    pore_radii = np.clip(pore_radii, 0.0, None)

    histogram, _ = np.histogram(pore_radii, bins=bin_edges)
    histogram = histogram.astype(float)
    total = histogram.sum()
    if total > 0:
        histogram /= total

    return bin_centers, histogram


def system_info(mol: Molecule, filename: str) -> dict[str, Any]:
    info: dict[str, Any] = {"filename": filename}
    if mol.pbc is not None:
        pbc = mol.pbc
        info["cell"] = {
            "a": pbc.a, "b": pbc.b, "c": pbc.c,
            "alpha": pbc.alpha, "beta": pbc.beta, "gamma": pbc.gamma,
        }
        info["volume"] = pbc.volume
        total_mass = sum(a.mass for a in mol.atoms)
        info["density"] = total_mass * 1.66054 / pbc.volume
        info["basis_matrix"] = pbc.basis_matrix.tolist()
    info["n_atoms"] = len(mol.atoms)
    info["formula"] = formula_unit(mol)
    info["formula_reduced"] = formula_unit_reduced(mol)
    return info


def print_info(mol: Molecule, filename: str) -> None:
    print("")
    print(r"   ___  ___  ___   __    __ _                  _ ")
    print(r"  / _ \/   \/ __\ / / /\ \ (_)______ _ _ __ __| |")
    print(r" / /_)/ /\ /__\// \ \/  \/ / |_  / _` | '__/ _` |")
    print(r"/ ___/ /_// \/  \  \  /\  /| |/ / (_| | | | (_| |")
    print(r"\/  /___,'\_____/   \/  \/ |_/___\__,_|_|  \__,_|")
    print(f"\nfilename: {filename}")
    if mol.pbc is not None:
        pbc = mol.pbc
        print(
            f"\nCell:\n{round(pbc.a, 3):>7}  {round(pbc.b, 3):7}  {round(pbc.c, 3):7} "
            f"{round(pbc.alpha, 2):>6} {round(pbc.beta, 2):>6} {round(pbc.gamma, 2):>6}\n"
        )
        for row in pbc.basis_matrix:
            print(f"{row[0]:20.14f} {row[1]:20.14f} {row[2]:20.14f}")
        total_mass = sum(a.mass for a in mol.atoms)
        density = total_mass * 1.66054 / pbc.volume
        print(f"Volume: {pbc.volume:10.2f} A^3 Density: {density:10.4f} g/cm^3")
    counts = formula_unit(mol)
    print("\nTotal number of atoms:\n")
    for ele, n in counts.items():
        print(f"{ele} {n}")
    reduced = formula_unit_reduced(mol)
    print("\nFormula unit\n")
    for ele, n in reduced.items():
        print(f"{ele} {n}")


def print_info_movie(n_frames: int, filename: str) -> None:
    print("")
    print(r"   ___  ___  ___   __    __ _                  _ ")
    print(r"  / _ \/   \/ __\ / / /\ \ (_)______ _ _ __ __| |")
    print(r" / /_)/ /\ /__\// \ \/  \/ / |_  / _` | '__/ _` |")
    print(r"/ ___/ /_// \/  \  \  /\  /| |/ / (_| | | | (_| |")
    print(r"\/  /___,'\_____/   \/  \/ |_/___\__,_|_|  \__,_|")
    print("")
    print(f"Trajectory detected\n{n_frames} frames\nExtract a single frame for more detailed options")

# ======================================================================
# Module: renderer
# ======================================================================
"""3D braille rendering engine (adapted from moltui, MO/isosurface code removed)."""



import numpy as np



def rotation_matrix(rx: float, ry: float, rz: float) -> np.ndarray:
    cx, sx = math.cos(rx), math.sin(rx)
    cy, sy = math.cos(ry), math.sin(ry)
    cz, sz = math.cos(rz), math.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


class ImageRenderer:
    def __init__(
        self,
        width: int,
        height: int,
        bg_color: tuple[int, int, int] = (0, 0, 0),
    ):
        self.width = width
        self.height = height
        self.fov = 1.5
        # When True, use parallel-projection (orthographic) instead of
        # perspective. The orthographic scale is set by `ortho_scale` —
        # render_scene picks a value matching the camera distance so that
        # toggling perspective/ortho keeps the molecule roughly the same
        # apparent size.
        self.orthographic = False
        self.ortho_scale = 1.0
        self.light_dir = np.array([0.4, 0.7, -0.6])
        self.light_dir /= np.linalg.norm(self.light_dir)
        self.ambient = 0.35
        self.diffuse_strength = 0.6
        self.specular_strength = 0.4
        self.shininess = 32.0
        self.view_dir = np.array([0.0, 0.0, -1.0])
        self.half_vec = self.light_dir + self.view_dir
        self.half_vec /= np.linalg.norm(self.half_vec)
        self.atom_scale = 0.35
        self.bond_radius = 0.08
        self.bg_color = np.array(bg_color, dtype=np.uint8)
        self.clear()

    def clear(self) -> None:
        self.pixels = np.tile(self.bg_color, (self.height, self.width, 1))
        self.z_buf = np.full((self.height, self.width), float("inf"))

    def _project(self, point: np.ndarray) -> tuple[float, float, float]:
        x, y, z = point
        scale = min(self.width, self.height) / 2
        if self.orthographic:
            # Parallel projection: px,py independent of z. Anything in front
            # of the camera (z > 0) renders; we still need a valid z for the
            # depth buffer, but no perspective divide.
            sx = self.width / 2 + x * self.ortho_scale * scale
            sy = self.height / 2 - y * self.ortho_scale * scale
            return (sx, sy, z if z > 0 else 0.1)
        if z <= 0.1:
            return (float("nan"), float("nan"), z)
        px = x * self.fov / z
        py = y * self.fov / z
        sx = self.width / 2 + px * scale
        sy = self.height / 2 - py * scale
        return (sx, sy, z)

    def render_sphere(
        self,
        center: np.ndarray,
        radius: float,
        color: tuple[int, int, int],
    ) -> None:
        sx, sy, sz = self._project(center)
        if math.isnan(sx):
            return

        scale = min(self.width, self.height) / 2
        if self.orthographic:
            pr = radius * self.ortho_scale * scale
        else:
            pr = radius * self.fov / sz * scale
        if pr < 0.5:
            return

        x_min = max(0, int(sx - pr - 1))
        x_max = min(self.width - 1, int(sx + pr + 1))
        y_min = max(0, int(sy - pr - 1))
        y_max = min(self.height - 1, int(sy + pr + 1))

        if x_min > x_max or y_min > y_max:
            return

        ys = np.arange(y_min, y_max + 1)
        xs = np.arange(x_min, x_max + 1)
        px_grid, py_grid = np.meshgrid(xs, ys)

        dx = (px_grid - sx) / pr
        dy = (py_grid - sy) / pr
        dist_sq = dx * dx + dy * dy
        mask = dist_sq <= 1.0

        dz = np.sqrt(np.maximum(0.0, 1.0 - dist_sq))
        norm_len = np.sqrt(dx * dx + dy * dy + dz * dz) + 1e-10
        nx = dx / norm_len
        ny = -dy / norm_len
        nz = dz / norm_len

        n_dot_l = np.maximum(
            0.0,
            nx * self.light_dir[0] + ny * self.light_dir[1] + nz * self.light_dir[2],
        )
        n_dot_h = np.maximum(
            0.0,
            nx * self.half_vec[0] + ny * self.half_vec[1] + nz * self.half_vec[2],
        )
        specular = np.power(n_dot_h, self.shininess) * self.specular_strength
        diffuse = n_dot_l * self.diffuse_strength
        intensity = np.minimum(1.0, self.ambient + diffuse)

        point_z = sz - radius * dz

        z_slice = self.z_buf[y_min: y_max + 1, x_min: x_max + 1]
        valid = mask & (point_z < z_slice)

        z_slice[valid] = point_z[valid]

        color_arr = np.array(color, dtype=np.float64)
        shaded = np.minimum(
            255, color_arr[None, None, :] * intensity[:, :, None] + 255 * specular[:, :, None]
        ).astype(np.uint8)
        for c in range(3):
            channel = self.pixels[y_min: y_max + 1, x_min: x_max + 1, c]
            channel[valid] = shaded[:, :, c][valid]

    def render_bond(
        self,
        p1: np.ndarray,
        p2: np.ndarray,
        color1: tuple[int, int, int],
        color2: tuple[int, int, int],
    ) -> None:
        sx1, sy1, sz1 = self._project(p1)
        sx2, sy2, sz2 = self._project(p2)
        if math.isnan(sx1) or math.isnan(sx2):
            return

        scale = min(self.width, self.height) / 2
        mid_z = (sz1 + sz2) / 2
        if self.orthographic:
            pr = self.bond_radius * self.ortho_scale * scale
        else:
            pr = self.bond_radius * self.fov / mid_z * scale

        dx = sx2 - sx1
        dy = sy2 - sy1
        length = math.sqrt(dx * dx + dy * dy)
        if length < 1.0:
            return

        # Skip if entirely off-screen
        margin = pr + 2
        if (min(sx1, sx2) > self.width + margin or max(sx1, sx2) < -margin or
                min(sy1, sy2) > self.height + margin or max(sy1, sy2) < -margin):
            return

        nx, ny = -dy / length, dx / length
        half_w = max(1.0, pr)
        steps = int(length * 2) + 1

        hw = int(half_w + 1)
        ts = np.linspace(0, 1, steps + 1)
        offsets = np.arange(-hw, hw + 1, dtype=np.float64)
        d_norm = offsets / half_w
        off_mask = np.abs(d_norm) <= 1.0
        offsets = offsets[off_mask]
        d_norm = d_norm[off_mask]

        cxs = sx1 + dx * ts
        cys = sy1 + dy * ts
        czs = sz1 + (sz2 - sz1) * ts

        all_px = np.round(cxs[:, None] + nx * offsets[None, :]).astype(int)
        all_py = np.round(cys[:, None] + ny * offsets[None, :]).astype(int)

        cyl_nz = np.sqrt(1.0 - d_norm * d_norm)
        cyl_nx_v = nx * d_norm
        cyl_ny_v = -ny * d_norm
        norm_len = np.sqrt(cyl_nx_v**2 + cyl_ny_v**2 + cyl_nz**2) + 1e-10
        cyl_nx_v /= norm_len
        cyl_ny_v /= norm_len
        cyl_nz /= norm_len
        diffuse = np.maximum(
            0.0,
            cyl_nx_v * self.light_dir[0]
            + cyl_ny_v * self.light_dir[1]
            + cyl_nz * self.light_dir[2],
        )
        intensity = np.minimum(1.0, self.ambient + (1.0 - self.ambient) * diffuse)
        pz = czs[:, None] - self.bond_radius * cyl_nz[None, :]
        intensity_2d = np.broadcast_to(intensity[None, :], pz.shape)

        flat_px = all_px.ravel()
        flat_py = all_py.ravel()
        flat_pz = pz.ravel()
        flat_int = intensity_2d.ravel()

        n_steps = len(ts)
        n_off = len(offsets)
        step_idx = np.repeat(np.arange(n_steps), n_off)
        half = n_steps // 2
        is_first_half = step_idx <= half

        valid = (flat_px >= 0) & (flat_px < self.width) & (flat_py >= 0) & (flat_py < self.height)
        flat_px, flat_py, flat_pz, flat_int, is_first_half = (
            flat_px[valid], flat_py[valid], flat_pz[valid], flat_int[valid], is_first_half[valid],
        )

        if len(flat_px) == 0:
            return

        z_current = self.z_buf[flat_py, flat_px]
        z_pass = flat_pz < z_current
        flat_px = flat_px[z_pass]
        flat_py = flat_py[z_pass]
        flat_pz = flat_pz[z_pass]
        flat_int = flat_int[z_pass]
        is_first_half = is_first_half[z_pass]

        if len(flat_px) == 0:
            return

        self.z_buf[flat_py, flat_px] = flat_pz

        c1 = np.array(color1, dtype=np.float64)
        c2 = np.array(color2, dtype=np.float64)
        colors = np.where(is_first_half[:, None], c1[None, :], c2[None, :])
        shaded = np.minimum(255, colors * flat_int[:, None]).astype(np.uint8)
        self.pixels[flat_py, flat_px] = shaded

    @staticmethod
    def _highlight_color() -> tuple[int, int, int]:
        return (255, 255, 50)

    @staticmethod
    def _charge_color(charge: float) -> tuple[int, int, int]:
        """Map charge to color: red (negative) -> white (zero) -> blue (positive)."""
        c = max(-1.0, min(1.0, charge))
        if c >= 0:
            # Positive = blue
            r = int(255 * (1 - c))
            g = int(255 * (1 - c))
            b = 255
        else:
            # Negative = red
            r = 255
            g = int(255 * (1 + c))
            b = int(255 * (1 + c))
        return (r, g, b)

    # Palette for residue / chain / index coloring
    _PALETTE = [
        (50, 130, 255), (255, 100, 80), (80, 220, 100), (255, 200, 50),
        (200, 100, 255), (100, 220, 220), (255, 150, 200), (180, 180, 100),
        (255, 130, 50), (100, 180, 255), (220, 80, 180), (80, 200, 160),
        (200, 200, 80), (150, 100, 200), (80, 160, 255), (255, 180, 130),
    ]

    def _compute_atom_colors(
        self, molecule: Molecule, color_mode: str,
    ) -> list[tuple[int, int, int]]:
        """Compute per-atom colors based on the selected mode."""
        if color_mode == "charge":
            return [self._charge_color(a.charge) for a in molecule.atoms]
        if color_mode == "residue":
            # Group by residue_id if set, otherwise by name runs
            pal = self._PALETTE
            colors = []
            prev_name = None
            group_id = -1
            for a in molecule.atoms:
                rid = getattr(a, "residue_id", 0)
                if rid:
                    colors.append(pal[rid % len(pal)])
                else:
                    if a.name != prev_name:
                        group_id += 1
                        prev_name = a.name
                    colors.append(pal[group_id % len(pal)])
            return colors
        if color_mode == "chain":
            # Color by chain breaks (CA-CA distance > 5A or residue_id resets)
            pal = self._PALETTE
            chain_id = 0
            colors: list[tuple[int, int, int]] = []
            prev_pos = None
            for a in molecule.atoms:
                if prev_pos is not None:
                    d = float(np.linalg.norm(a.position - prev_pos))
                    if d > 8.0:
                        chain_id += 1
                prev_pos = a.position
                colors.append(pal[chain_id % len(pal)])
            return colors
        if color_mode == "index":
            # Rainbow gradient by atom index
            n = max(1, len(molecule.atoms))
            colors = []
            for i in range(len(molecule.atoms)):
                frac = i / n
                # HSV-like rainbow: red -> yellow -> green -> cyan -> blue -> magenta
                h = frac * 6.0
                x = int(255 * (1 - abs(h % 2 - 1)))
                if h < 1:
                    colors.append((255, x, 0))
                elif h < 2:
                    colors.append((x, 255, 0))
                elif h < 3:
                    colors.append((0, 255, x))
                elif h < 4:
                    colors.append((0, x, 255))
                elif h < 5:
                    colors.append((x, 0, 255))
                else:
                    colors.append((255, 0, x))
            return colors
        # Default: element CPK
        return [a.element.cpk_color for a in molecule.atoms]

    def render_molecule(
        self,
        molecule: Molecule,
        rot: np.ndarray,
        camera_distance: float,
        pbc: PBC | None = None,
        pan: tuple[float, float] = (0.0, 0.0),
        highlighted_atoms: set[int] | None = None,
        licorice: bool = False,
        vdw: bool = False,
        color_by_charge: bool = False,
        color_mode: str = "element",
        centroid_override: np.ndarray | None = None,
        ribbon: bool = False,
        show_polyhedra: bool = False,
        density_positions: np.ndarray | None = None,
        density_values: np.ndarray | None = None,
    ) -> None:
        self.clear()
        # color_by_charge is legacy — map to color_mode
        if color_by_charge and color_mode == "element":
            color_mode = "charge"
        centroid = centroid_override if centroid_override is not None else (
            molecule.center() if molecule.atoms else np.zeros(3)
        )
        hl = highlighted_atoms or set()
        has_hl = len(hl) > 0

        # Render unit cell wireframe if PBC available
        if pbc is not None:
            corners = pbc.corners()
            # Center the cell on the atom centroid
            cell_center = np.array([0.5, 0.5, 0.5]) @ pbc.basis_matrix
            cell_offset = centroid - cell_center
            corners = corners + cell_offset
            wire_color = (255, 255, 255)
            transformed_corners = []
            for corner in corners:
                pos = rot @ (corner - centroid)
                pos[0] += pan[0]
                pos[1] += pan[1]
                pos[2] += camera_distance
                transformed_corners.append(pos)
            for i, j in pbc.edges():
                self.render_bond(
                    transformed_corners[i], transformed_corners[j],
                    wire_color, wire_color,
                )

        # Snapshot atoms once — molecule may be mutated concurrently by
        # background threads (extend_axis, sort, etc.)
        snap_atoms = list(molecule.atoms)
        n_snap = len(snap_atoms)
        if n_snap > 0:
            positions = np.array([a.position for a in snap_atoms])
            transformed_arr = (positions - centroid) @ rot.T
            transformed_arr[:, 0] += pan[0]
            transformed_arr[:, 1] += pan[1]
            transformed_arr[:, 2] += camera_distance
            transformed = [transformed_arr[i] for i in range(n_snap)]
        else:
            transformed = []

        # Snapshot bonds too — bonds may reference indices that exist in our snapshot
        snap_bonds = list(molecule.bonds) if molecule.bonds else []

        # Ribbon mode: draw backbone trace through CA atoms
        if ribbon and snap_atoms:
            ca_indices = [i for i, a in enumerate(snap_atoms) if a.name.strip() == "CA"]
            if len(ca_indices) < 2:
                if snap_bonds:
                    neighbors: dict[int, set[int]] = {}
                    for a, b in snap_bonds:
                        if a < n_snap and b < n_snap:
                            neighbors.setdefault(a, set()).add(b)
                            neighbors.setdefault(b, set()).add(a)
                    ca_indices = []
                    for i, a in enumerate(snap_atoms):
                        if a.element.symbol != "C":
                            continue
                        nbrs = neighbors.get(i, set())
                        has_n = any(snap_atoms[j].element.symbol == "N" for j in nbrs if j < n_snap)
                        has_c = any(snap_atoms[j].element.symbol == "C" for j in nbrs if j < n_snap)
                        if has_n and has_c:
                            ca_indices.append(i)
            if len(ca_indices) >= 2:
                _chain_colors = [
                    (50, 130, 255), (255, 100, 80), (80, 220, 100),
                    (255, 200, 50), (200, 100, 255), (100, 220, 220),
                    (255, 150, 200), (180, 180, 100),
                ]
                chain_id = 0
                ca_chains: list[int] = [0]
                for k in range(1, len(ca_indices)):
                    d = np.linalg.norm(
                        snap_atoms[ca_indices[k]].position
                        - snap_atoms[ca_indices[k - 1]].position
                    )
                    if d > 5.0:
                        chain_id += 1
                    ca_chains.append(chain_id)
                old_br = self.bond_radius
                self.bond_radius = 0.25
                for k in range(len(ca_indices) - 1):
                    if ca_chains[k] != ca_chains[k + 1]:
                        continue
                    c = _chain_colors[ca_chains[k] % len(_chain_colors)]
                    p0 = transformed[ca_indices[max(k - 1, 0)]]
                    p1 = transformed[ca_indices[k]]
                    p2 = transformed[ca_indices[k + 1]]
                    p3 = transformed[ca_indices[min(k + 2, len(ca_indices) - 1)]]
                    n_seg = 4
                    prev = p1
                    for s in range(1, n_seg + 1):
                        t = s / n_seg
                        t2 = t * t
                        t3 = t2 * t
                        pt = 0.5 * (
                            (2 * p1)
                            + (-p0 + p2) * t
                            + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
                            + (-p0 + 3 * p1 - 3 * p2 + p3) * t3
                        )
                        self.render_bond(prev, pt, c, c)
                        prev = pt
                self.bond_radius = old_br
        else:
            # Normal atom/bond rendering (non-ribbon)
            # Use snapshot for atom colors to match transformed length
            if color_mode == "charge":
                atom_colors = [self._charge_color(a.charge) for a in snap_atoms]
            elif color_mode in ("residue", "chain", "index"):
                # Build a temporary molecule-like object from snapshot
                _snap_mol = type("M", (), {"atoms": snap_atoms})()
                atom_colors = self._compute_atom_colors(_snap_mol, color_mode)
            else:
                atom_colors = [a.element.cpk_color for a in snap_atoms]

            if not vdw:
                for i, j in snap_bonds:
                    if i >= n_snap or j >= n_snap:
                        continue  # skip bonds that reference atoms outside snapshot
                    c1, c2 = atom_colors[i], atom_colors[j]
                    if has_hl and i in hl and j in hl:
                        c1 = self._highlight_color()
                        c2 = self._highlight_color()
                    self.render_bond(transformed[i], transformed[j], c1, c2)

            atom_order = sorted(
                range(n_snap),
                key=lambda idx: -transformed[idx][2],
            )
            for i in atom_order:
                atom = snap_atoms[i]
                if vdw:
                    radius = atom.element.vdw_radius
                elif licorice:
                    radius = self.bond_radius
                else:
                    radius = atom.element.covalent_radius * self.atom_scale
                color = self._highlight_color() if has_hl and i in hl else atom_colors[i]
                self.render_sphere(transformed[i], radius, color)

        # Draw 3D density blobs
        if density_positions is not None and density_values is not None and len(density_positions) > 0:
            for di in range(len(density_positions)):
                pos = rot @ (density_positions[di] - centroid)
                pos[0] += pan[0]
                pos[1] += pan[1]
                pos[2] += camera_distance
                v = float(density_values[di])
                # Color: low density = blue, high = red
                r_c = int(255 * v)
                b_c = int(255 * (1 - v))
                g_c = int(100 * min(v, 1 - v) * 2)
                radius = 0.2 + v * 0.4
                self.render_sphere(pos, radius, (r_c, g_c, b_c))

        # Draw coordination polyhedra around metal atoms
        if show_polyhedra and snap_atoms and transformed:
            _metals = {
                "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V",
                "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr",
                "Y", "Zr", "Nb", "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
                "Sn", "Cs", "Ba", "La", "Ce", "Hf", "Ta", "W", "Re", "Os",
                "Ir", "Pt", "Au", "Pb", "Bi",
            }
            for ci, atom in enumerate(snap_atoms):
                if atom.element.symbol not in _metals:
                    continue
                cutoff = atom.element.covalent_radius * 2.5 + 0.5
                center_pos = atom.position
                nbr_indices = []
                for ni, natom in enumerate(snap_atoms):
                    if ni == ci:
                        continue
                    diff = natom.position - center_pos
                    if pbc is not None:
                        frac = diff @ pbc.reciprocal_basis_matrix
                        frac -= np.round(frac)
                        diff = frac @ pbc.basis_matrix
                    d = float(np.linalg.norm(diff))
                    if d < cutoff:
                        nbr_indices.append(ni)
                if len(nbr_indices) < 3:
                    continue
                # Draw edges between all neighboring atoms (polyhedron wireframe)
                color = atom.element.cpk_color
                dim = (color[0] // 2, color[1] // 2, color[2] // 2)
                for ai in range(len(nbr_indices)):
                    for bi in range(ai + 1, len(nbr_indices)):
                        self.render_bond(
                            transformed[nbr_indices[ai]],
                            transformed[nbr_indices[bi]],
                            dim, dim,
                        )

        # Draw axis arrows in bottom-left corner (always, including ribbon)
        if pbc is not None:
            self._render_axes(rot, camera_distance, pbc)

    def _render_axes(self, rot: np.ndarray, camera_distance: float,
                     pbc: "PBC") -> None:
        """Draw a/b/c axis arrows in the bottom-left corner."""
        # Normalize basis vectors to equal visual length
        basis = pbc.basis_matrix.copy()
        for i in range(3):
            basis[i] = basis[i] / (np.linalg.norm(basis[i]) + 1e-10)

        # Place the axis origin in screen space: bottom-left
        arrow_len = camera_distance * 0.08
        # Origin in view space (offset to bottom-left)
        ox = -camera_distance * 0.85
        oy = -camera_distance * 0.45
        origin = np.array([ox, oy, camera_distance])

        colors = [(255, 80, 80), (80, 255, 80), (80, 80, 255)]  # a=red, b=green, c=blue
        for i in range(3):
            direction = rot @ basis[i]
            tip = origin + direction * arrow_len
            self.render_bond(origin, tip, colors[i], colors[i])
            # Small sphere at tip
            self.render_sphere(tip, arrow_len * 0.15, colors[i])


def render_scene(
    width: int,
    height: int,
    molecule: Molecule,
    rot: np.ndarray,
    camera_distance: float,
    bg_color: tuple[int, int, int] = (0, 0, 0),
    pbc: PBC | None = None,
    ssaa: int = 2,
    pan: tuple[float, float] = (0.0, 0.0),
    highlighted_atoms: set[int] | None = None,
    licorice: bool = False,
    vdw: bool = False,
    ambient: float | None = None,
    diffuse: float | None = None,
    specular: float | None = None,
    shininess: float | None = None,
    atom_scale: float | None = None,
    bond_radius: float | None = None,
    color_by_charge: bool = False,
    color_mode: str = "element",
    centroid_override: np.ndarray | None = None,
    ribbon: bool = False,
    show_polyhedra: bool = False,
    density_positions: np.ndarray | None = None,
    density_values: np.ndarray | None = None,
    orthographic: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Render with supersampling anti-aliasing.

    Returns (pixels, hit_mask) where pixels is (height, width, 3) uint8
    and hit_mask is (height, width) bool indicating which pixels were drawn.
    """
    r = ImageRenderer(width * ssaa, height * ssaa, bg_color=bg_color)
    if orthographic:
        r.orthographic = True
        # Match the on-screen size of perspective at the chosen camera
        # distance: in perspective, an object at z=camera_distance is
        # scaled by fov/camera_distance; replicate that as ortho_scale.
        r.ortho_scale = r.fov / max(camera_distance, 1e-3)
    if ambient is not None:
        r.ambient = ambient
    if diffuse is not None:
        r.diffuse_strength = diffuse
    if specular is not None:
        r.specular_strength = specular
    if shininess is not None:
        r.shininess = shininess
    if atom_scale is not None:
        r.atom_scale = atom_scale
    if bond_radius is not None:
        r.bond_radius = bond_radius
    r.render_molecule(
        molecule, rot, camera_distance,
        pbc=pbc, pan=pan, highlighted_atoms=highlighted_atoms,
        licorice=licorice, vdw=vdw, color_by_charge=color_by_charge,
        color_mode=color_mode,
        centroid_override=centroid_override,
        ribbon=ribbon, show_polyhedra=show_polyhedra,
        density_positions=density_positions, density_values=density_values,
    )
    hit = np.isfinite(r.z_buf)
    if ssaa == 1:
        return r.pixels, hit
    downsampled = r.pixels.reshape(height, ssaa, width, ssaa, 3).mean(axis=(1, 3)).astype(np.uint8)
    hit_down = hit.reshape(height, ssaa, width, ssaa).any(axis=(1, 3))
    return downsampled, hit_down

# ======================================================================
# Module: charges
# ======================================================================
"""QEq charge equilibration (Rappe & Goddard, J. Phys. Chem. 1991, 95, 3358).

Assigns partial atomic charges by equalizing electronegativity across all atoms
subject to total charge neutrality. Only requires atomic parameters (chi, J)
and interatomic distances.
"""


import numpy as np


# Rappe-Goddard QEq parameters: chi (eV), J (eV)
# From OpenBabel qeq.txt / UFF GMP parameters
# fmt: off
QEQ_PARAMS: dict[str, tuple[float, float]] = {
    "H":  (4.5280, 13.8904),
    "He": (9.66,   29.84),
    "Li": (3.006,   4.772),
    "Be": (4.877,   8.886),
    "B":  (5.11,    9.50),
    "C":  (5.343,  10.126),
    "N":  (6.899,  11.760),
    "O":  (8.741,  13.364),
    "F":  (10.874, 14.948),
    "Ne": (11.04,  21.10),
    "Na": (2.843,   4.592),
    "Mg": (3.951,   7.386),
    "Al": (4.06,    7.18),
    "Si": (4.168,   6.974),
    "P":  (5.463,   8.000),
    "S":  (6.928,   8.972),
    "Cl": (8.564,   9.892),
    "Ar": (9.465,  12.71),
    "K":  (2.421,   3.840),
    "Ca": (3.231,   5.76),
    "Sc": (3.395,   6.16),
    "Ti": (3.47,    6.76),
    "V":  (3.65,    6.82),
    "Cr": (3.415,   7.73),
    "Mn": (3.325,   8.21),
    "Fe": (3.76,    8.28),
    "Co": (4.105,   8.35),
    "Ni": (4.465,   8.41),
    "Cu": (4.20,    8.44),
    "Zn": (5.106,   8.57),
    "Ga": (3.641,   6.32),
    "Ge": (4.051,   6.876),
    "As": (5.188,   7.618),
    "Se": (6.428,   8.262),
    "Br": (7.790,   8.850),
    "Kr": (8.505,  11.43),
    "Rb": (2.331,   3.692),
    "Sr": (3.024,   4.88),
    "Y":  (3.83,    5.62),
    "Zr": (3.40,    7.10),
    "Nb": (3.55,    6.76),
    "Mo": (3.465,   7.51),
    "Tc": (3.29,    7.98),
    "Ru": (3.575,   8.03),
    "Rh": (3.975,   8.01),
    "Pd": (4.32,    8.00),
    "Ag": (4.436,   6.268),
    "Cd": (5.034,   7.914),
    "In": (3.506,   5.792),
    "Sn": (3.987,   6.248),
    "Sb": (4.899,   6.684),
    "Te": (5.816,   7.052),
    "I":  (6.822,   7.524),
    "Xe": (7.595,   9.95),
    "Cs": (2.183,   3.422),
    "Ba": (2.814,   4.792),
    "La": (2.8355,  5.483),
    "Ce": (2.774,   5.384),
    "Pr": (2.858,   5.128),
    "Nd": (2.8685,  5.241),
    "Pm": (2.881,   5.346),
    "Sm": (2.9115,  5.439),
    "Eu": (2.8785,  5.575),
    "Gd": (3.1665,  5.949),
    "Tb": (3.018,   5.668),
    "Dy": (3.0555,  5.743),
    "Ho": (3.127,   5.782),
    "Er": (3.1865,  5.829),
    "Tm": (3.2514,  5.8658),
    "Yb": (3.2889,  5.93),
    "Lu": (2.9629,  4.9258),
    "Hf": (3.70,    6.80),
    "Ta": (5.10,    5.70),
    "W":  (4.63,    6.62),
    "Re": (3.96,    7.84),
    "Os": (5.14,    7.26),
    "Ir": (5.00,    8.00),
    "Pt": (4.79,    8.86),
    "Au": (4.894,   5.172),
    "Hg": (6.27,    8.32),
    "Tl": (3.20,    5.80),
    "Pb": (3.90,    7.06),
    "Bi": (4.69,    7.48),
}
# fmt: on

# Coulomb constant: 14.4 eV*A/e^2
_KE = 14.3996


def compute_qeq_charges(
    mol: Molecule,
    total_charge: float = 0.0,
    progress_callback=None,
) -> np.ndarray:
    """Compute QEq charges for all atoms in the molecule.

    Uses Rappe-Goddard charge equilibration with a cutoff distance for
    the Coulomb interaction to keep memory and compute manageable.

    Args:
        mol: Molecule with atoms (and optionally PBC).
        total_charge: Total system charge (0 for neutral).
        progress_callback: Optional callable(frac: float) for progress updates.

    Returns:
        Array of partial charges (electrons).
    """
    n = len(mol.atoms)
    if n == 0:
        return np.array([])

    # Get parameters
    chi = np.zeros(n)
    J_diag = np.zeros(n)
    for i, atom in enumerate(mol.atoms):
        sym = atom.element.symbol
        if sym in QEQ_PARAMS:
            chi[i], J_diag[i] = QEQ_PARAMS[sym]
        else:
            chi[i], J_diag[i] = 5.0, 10.0

    coords = np.array([a.x for a in mol.atoms])

    if progress_callback:
        progress_callback(0.05)

    if n > 3000:
        raise ValueError(
            f"System too large for QEq ({n} atoms). "
            "Try on a single unit cell before extending."
        )

    b = np.zeros(n + 1)
    for i in range(n):
        b[i] = -chi[i]
    b[n] = total_charge

    if progress_callback:
        progress_callback(0.1)

    # Compute pairwise Coulomb interactions chunked
    i_all, j_all = np.triu_indices(n, k=1)
    total_pairs = len(i_all)
    chunk_size = max(500000, total_pairs // 10)

    # Collect COO data for off-diagonal entries
    row_list = []
    col_list = []
    val_list = []

    for start in range(0, total_pairs, chunk_size):
        end = min(start + chunk_size, total_pairs)
        ic = i_all[start:end]
        jc = j_all[start:end]
        dx = coords[ic] - coords[jc]
        if mol.pbc is not None:
            frac = dx @ mol.pbc.reciprocal_basis_matrix
            dx = dx - np.round(frac) @ mol.pbc.basis_matrix
        rs = np.linalg.norm(dx, axis=1)
        rs = np.maximum(rs, 0.1)
        j_ij = _KE / rs

        # Apply cutoff for memory: ignore very weak interactions
        cutoff = 20.0
        mask = rs < cutoff
        row_list.append(ic[mask])
        row_list.append(jc[mask])
        col_list.append(jc[mask])
        col_list.append(ic[mask])
        val_list.append(j_ij[mask])
        val_list.append(j_ij[mask])

        if progress_callback:
            progress_callback(0.1 + 0.6 * end / total_pairs)

    # Build dense matrix (safe up to ~5000 atoms = 200MB)
    A = np.zeros((n + 1, n + 1))

    # Diagonal
    for i in range(n):
        A[i, i] = J_diag[i]

    # Off-diagonal from COO lists
    if row_list:
        rows = np.concatenate(row_list)
        cols = np.concatenate(col_list)
        vals = np.concatenate(val_list)
        # Direct assignment is safe: each (row,col) appears exactly once
        A[rows, cols] = vals

    # Lagrange multiplier
    A[n, :n] = 1.0
    A[:n, n] = 1.0

    if progress_callback:
        progress_callback(0.8)

    # Solve
    try:
        result = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    charges = result[:n]

    if progress_callback:
        progress_callback(1.0)

    return charges


def apply_qeq_charges(mol: Molecule, progress_callback=None) -> None:
    """Compute and apply QEq charges to all atoms in the molecule."""
    charges = compute_qeq_charges(mol, progress_callback=progress_callback)
    for i, atom in enumerate(mol.atoms):
        atom.charge = float(charges[i])

# ======================================================================
# Module: rdf
# ======================================================================
"""Radial distribution function (RDF) calculation with PBC support."""


import numpy as np



def compute_rdf(
    mol: Molecule,
    type1: str,
    type2: str,
    n_bins: int = 200,
    r_max: float | None = None,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the RDF g(r) between two atom types for a single frame.

    Args:
        mol: Molecule with PBC.
        type1: Element symbol for first atom type.
        type2: Element symbol for second atom type.
        n_bins: Number of histogram bins.
        r_max: Maximum distance (default: half the shortest cell dimension).
        progress_callback: Optional callable(frac: float).

    Returns:
        (r_values, g_r) arrays.
    """
    if mol.pbc is None:
        return np.array([]), np.array([])

    pbc = mol.pbc
    if r_max is None:
        r_max = min(pbc.a, pbc.b, pbc.c) / 2.0

    idx1 = [i for i, a in enumerate(mol.atoms) if a.element.symbol == type1]
    idx2 = [i for i, a in enumerate(mol.atoms) if a.element.symbol == type2]

    if not idx1 or not idx2:
        return np.array([]), np.array([])

    coords = np.array([a.x for a in mol.atoms])
    dr = r_max / n_bins
    hist = np.zeros(n_bins)

    same_type = type1 == type2
    total = len(idx1)

    for ci, i in enumerate(idx1):
        targets = idx2 if not same_type else [j for j in idx2 if j > i]
        if not targets:
            continue
        dx = coords[targets] - coords[i]
        frac = dx @ pbc.reciprocal_basis_matrix
        wrapped = dx - np.round(frac) @ pbc.basis_matrix
        dists = np.linalg.norm(wrapped, axis=1)
        bins = (dists / dr).astype(int)
        valid = (bins >= 0) & (bins < n_bins) & (dists < r_max)
        for b in bins[valid]:
            hist[b] += 1

        if progress_callback and ci % max(1, total // 20) == 0:
            progress_callback(ci / total)

    if progress_callback:
        progress_callback(1.0)

    # Normalize
    r_values = (np.arange(n_bins) + 0.5) * dr
    n1 = len(idx1)
    n2 = len(idx2)
    volume = pbc.volume

    shell_vol = 4.0 * np.pi * r_values**2 * dr
    if same_type:
        # Counted each pair once (j > i). Normalize by N*(N-1)/2 pairs
        # and ideal gas density.
        rho = n1 / volume
        ideal = rho * shell_vol
        g_r = np.zeros(n_bins)
        nonzero = ideal > 0
        # Factor 2: we counted j>i only, but ideal counts full shell
        g_r[nonzero] = 2.0 * hist[nonzero] / (n1 * ideal[nonzero])
    else:
        rho = n2 / volume
        ideal = rho * shell_vol
        g_r = np.zeros(n_bins)
        nonzero = ideal > 0
        g_r[nonzero] = hist[nonzero] / (n1 * ideal[nonzero])

    return r_values, g_r


def compute_rdf_trajectory(
    frames: list[Molecule],
    type1: str,
    type2: str,
    n_bins: int = 200,
    r_max: float | None = None,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the RDF averaged over multiple trajectory frames."""
    if not frames:
        return np.array([]), np.array([])

    g_sum = None
    r_values = None
    n_frames = len(frames)

    for fi, mol in enumerate(frames):
        r, g = compute_rdf(mol, type1, type2, n_bins=n_bins, r_max=r_max)
        if len(r) == 0:
            continue
        if g_sum is None:
            g_sum = g.copy()
            r_values = r
        else:
            g_sum += g
        if progress_callback:
            progress_callback((fi + 1) / n_frames)

    if g_sum is None:
        return np.array([]), np.array([])

    return r_values, g_sum / n_frames

# ======================================================================
# Module: msd
# ======================================================================
"""Mean Squared Displacement (MSD) calculation for diffusion analysis."""


import numpy as np



def compute_msd(
    frames: list[Molecule],
    element: str | None = None,
    name: str | None = None,
    max_lag: int | None = None,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute MSD vs time lag for selected atoms across trajectory frames.

    Uses the direct method: MSD(tau) = <|r(t+tau) - r(t)|^2> averaged over
    all atoms and all time origins t.

    Args:
        frames: List of Molecule objects (trajectory).
        element: Filter to this element (None = all moveable atoms).
        name: Filter to this atom name (None = all matching element).
        max_lag: Maximum lag in frames (default: half the trajectory).
        progress_callback: Optional callable(frac).

    Returns:
        (lag_frames, msd_values) arrays.
    """
    if len(frames) < 2:
        return np.array([]), np.array([])

    # Identify which atoms to track — must be present in all frames
    # Use first frame to get indices
    ref = frames[0]
    if element:
        if name and name != "(all)":
            def atom_filter(a):
                return a.element.symbol == element and a.name.strip() == name
        else:
            def atom_filter(a):
                return a.element.symbol == element
    else:
        # Track all non-framework atoms (moveable = charge != 0 or element is H2-like)
        def atom_filter(a):
            return True

    ref_indices = [i for i, a in enumerate(ref.atoms) if atom_filter(a)]
    if not ref_indices:
        return np.array([]), np.array([])

    # Build position arrays: find matching atoms in each frame by name+position proximity
    # For MPMC trajectories, atom count varies. Track by name matching within element.
    n_frames = len(frames)
    if max_lag is None:
        max_lag = n_frames // 2

    # Collect positions per frame for tracked atoms
    # Use index-based tracking (assumes atom ordering is consistent for matching atoms)
    n_track = len(ref_indices)
    positions = np.zeros((n_frames, n_track, 3))
    valid_count = np.zeros(n_frames, dtype=int)

    for fi, mol in enumerate(frames):
        indices = [i for i, a in enumerate(mol.atoms) if atom_filter(a)]
        n_avail = min(len(indices), n_track)
        for ti in range(n_avail):
            positions[fi, ti] = mol.atoms[indices[ti]].x
        valid_count[fi] = n_avail
        if progress_callback and fi % max(1, n_frames // 20) == 0:
            progress_callback(fi / n_frames * 0.5)

    # PBC unwrap: wrapped trajectories produce garbage MSD past first cell
    # crossing. Reconstruct unwrapped positions by accumulating min-imaged
    # frame-to-frame displacements.
    pbc = ref.pbc
    if pbc is not None:
        recip = pbc.reciprocal_basis_matrix
        basis = pbc.basis_matrix
        unwrapped = np.zeros_like(positions)
        unwrapped[0] = positions[0]
        for fi in range(1, n_frames):
            dx = positions[fi] - positions[fi - 1]
            frac = dx @ recip
            frac -= np.round(frac)
            dx_min = frac @ basis
            unwrapped[fi] = unwrapped[fi - 1] + dx_min
        positions = unwrapped

    # Compute MSD for each lag — only use atoms valid in both frames
    min_valid = int(valid_count.min()) if len(valid_count) > 0 else 0
    lags = np.arange(1, max_lag + 1)
    msd = np.zeros(len(lags))

    for li, lag in enumerate(lags):
        displacements = positions[lag:, :min_valid] - positions[:n_frames - lag, :min_valid]
        sq_disp = np.sum(displacements ** 2, axis=2)
        msd[li] = sq_disp.mean() if sq_disp.size > 0 else 0.0
        if progress_callback and li % max(1, len(lags) // 20) == 0:
            progress_callback(0.5 + li / len(lags) * 0.5)

    if progress_callback:
        progress_callback(1.0)

    return lags.astype(float), msd


def compute_diffusion_constant(
    lags: np.ndarray,
    msd: np.ndarray,
    dt: float = 1.0,
    fit_start: float = 0.25,
    fit_end: float = 0.75,
) -> tuple[float, float]:
    """Compute diffusion constant from MSD via Einstein relation.

    D = slope / (2 * d) where d=3 for 3D diffusion, so D = slope / 6.
    slope is d(MSD)/d(t) in A^2/ps.

    Args:
        lags: Lag values in frames.
        msd: MSD values in A^2.
        dt: Time between frames in picoseconds.
        fit_start: Fraction of lag range to start fitting (skip ballistic regime).
        fit_end: Fraction of lag range to end fitting (skip noisy tail).

    Returns:
        (D_cm2_s, slope_A2_ps): Diffusion constant in cm^2/s and slope in A^2/ps.
    """
    if len(lags) < 4:
        return 0.0, 0.0
    n = len(lags)
    i_start = max(1, int(n * fit_start))
    i_end = max(i_start + 2, int(n * fit_end))
    t = lags[i_start:i_end] * dt  # time in ps
    m = msd[i_start:i_end]
    slope, _ = np.polyfit(t, m, 1)  # A^2/ps
    D_A2_ps = slope / 6.0
    D_cm2_s = D_A2_ps * 1e-16 / 1e-12  # 1 A^2 = 1e-16 cm^2, 1 ps = 1e-12 s
    return D_cm2_s, slope

# ======================================================================
# Module: rmsd
# ======================================================================
"""Root Mean Square Deviation (RMSD) calculation over trajectories."""


import numpy as np



def compute_rmsd(
    frames: list[Molecule],
    element: str | None = None,
    name: str | None = None,
    reference_frame: int = 0,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute RMSD vs frame number relative to a reference frame.

    RMSD(t) = sqrt( (1/N) * sum_i |r_i(t) - r_i(ref)|^2 )

    Args:
        frames: List of Molecule objects (trajectory).
        element: Filter to this element (None = all atoms).
        name: Filter to this atom name (None = all matching element).
        reference_frame: Frame index to use as reference (default 0).
        progress_callback: Optional callable(frac).

    Returns:
        (frame_indices, rmsd_values) arrays.
    """
    if len(frames) < 2:
        return np.array([]), np.array([])

    ref = frames[reference_frame]

    if element:
        if name and name != "(all)":
            def atom_filter(a):
                return a.element.symbol == element and a.name.strip() == name
        else:
            def atom_filter(a):
                return a.element.symbol == element
    else:
        def atom_filter(a):
            return True

    ref_indices = [i for i, a in enumerate(ref.atoms) if atom_filter(a)]
    if not ref_indices:
        return np.array([]), np.array([])

    ref_pos = np.array([ref.atoms[i].x for i in ref_indices])
    n_atoms = len(ref_indices)
    n_frames = len(frames)

    frame_idx = np.arange(n_frames, dtype=float)
    rmsd_vals = np.zeros(n_frames)

    for fi, mol in enumerate(frames):
        indices = [i for i, a in enumerate(mol.atoms) if atom_filter(a)]
        n_avail = min(len(indices), n_atoms)
        if n_avail == 0:
            continue
        pos = np.array([mol.atoms[indices[i]].x for i in range(n_avail)])
        diff = pos[:n_avail] - ref_pos[:n_avail]
        rmsd_vals[fi] = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))

        if progress_callback and fi % max(1, n_frames // 20) == 0:
            progress_callback(fi / n_frames)

    if progress_callback:
        progress_callback(1.0)

    return frame_idx, rmsd_vals

# ======================================================================
# Module: coordination
# ======================================================================
"""Coordination number analysis."""


import numpy as np



def compute_coordination(
    mol: Molecule,
    center_element: str,
    neighbor_element: str,
    r_max: float = 5.0,
    n_bins: int = 200,
    center_name: str | None = None,
    neighbor_name: str | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute coordination number as a function of distance.

    Returns the running coordination number N(r) — the average number of
    neighbor atoms within distance r of each center atom.

    Args:
        mol: Molecule (with or without PBC).
        center_element: Element symbol for the central atoms.
        neighbor_element: Element symbol for the neighbor atoms.
        r_max: Maximum distance in Angstroms.
        n_bins: Number of distance bins.
        center_name: Optional atom name filter for center.
        neighbor_name: Optional atom name filter for neighbor.

    Returns:
        (r_values, coordination_number, avg_at_cutoff):
        r_values and running N(r), plus the average coordination at r_max.
    """
    center_idx = []
    for i, a in enumerate(mol.atoms):
        if a.element.symbol != center_element:
            continue
        if center_name and center_name != "(all)" and a.name.strip() != center_name:
            continue
        center_idx.append(i)

    neighbor_idx = []
    for i, a in enumerate(mol.atoms):
        if a.element.symbol != neighbor_element:
            continue
        if neighbor_name and neighbor_name != "(all)" and a.name.strip() != neighbor_name:
            continue
        neighbor_idx.append(i)

    if not center_idx or not neighbor_idx:
        return np.array([]), np.array([]), 0.0

    center_pos = np.array([mol.atoms[i].x for i in center_idx])
    neighbor_pos = np.array([mol.atoms[i].x for i in neighbor_idx])

    r_values = np.linspace(0, r_max, n_bins)
    dr = r_values[1] - r_values[0]
    histogram = np.zeros(n_bins)

    use_pbc = mol.pbc is not None

    for ci in range(len(center_pos)):
        diffs = neighbor_pos - center_pos[ci]
        if use_pbc:
            # Minimum image convention
            frac = diffs @ mol.pbc.reciprocal_basis_matrix
            frac -= np.round(frac)
            diffs = frac @ mol.pbc.basis_matrix
        dists = np.sqrt(np.sum(diffs ** 2, axis=1))
        # Exclude self (distance ~0) and clamp to r_max
        mask = (dists > 0.1) & (dists < r_max)
        if not mask.any():
            continue
        bins = np.minimum((dists[mask] / dr).astype(int), n_bins - 1)
        histogram += np.bincount(bins, minlength=n_bins)[:n_bins]

    # Average per center atom
    n_center = len(center_idx)
    histogram /= n_center

    # Running coordination number = cumulative sum
    coord_number = np.cumsum(histogram)

    avg_at_cutoff = float(coord_number[-1]) if len(coord_number) > 0 else 0.0

    return r_values, coord_number, avg_at_cutoff


def compute_coordination_trajectory(
    frames: list[Molecule],
    center_element: str,
    neighbor_element: str,
    r_max: float = 5.0,
    n_bins: int = 200,
    center_name: str | None = None,
    neighbor_name: str | None = None,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Average coordination number over trajectory frames."""
    if not frames:
        return np.array([]), np.array([]), 0.0

    r_values = np.linspace(0, r_max, n_bins)
    coord_sum = np.zeros(n_bins)

    for fi, mol in enumerate(frames):
        _, cn, _ = compute_coordination(
            mol, center_element, neighbor_element, r_max, n_bins,
            center_name, neighbor_name,
        )
        if len(cn) == n_bins:
            coord_sum += cn
        if progress_callback and fi % max(1, len(frames) // 20) == 0:
            progress_callback(fi / len(frames))

    coord_avg = coord_sum / len(frames)
    avg_at_cutoff = float(coord_avg[-1]) if len(coord_avg) > 0 else 0.0

    if progress_callback:
        progress_callback(1.0)

    return r_values, coord_avg, avg_at_cutoff

# ======================================================================
# Module: hbonds
# ======================================================================
"""Hydrogen bond detection and analysis."""


import numpy as np



def detect_hbonds(
    mol: Molecule,
    d_cutoff: float = 3.5,
    angle_cutoff: float = 120.0,
) -> list[tuple[int, int, int, float, float]]:
    """Detect hydrogen bonds: D-H...A where D and A are electronegative.

    Args:
        mol: Molecule with bonds detected.
        d_cutoff: Max D-A distance in Angstroms.
        angle_cutoff: Min D-H-A angle in degrees.

    Returns:
        List of (donor_idx, h_idx, acceptor_idx, distance, angle).
    """
    _donors_acceptors = {"N", "O", "F", "S", "Cl"}

    # Build neighbor map
    neighbors: dict[int, list[int]] = {}
    for a, b in mol.bonds:
        neighbors.setdefault(a, []).append(b)
        neighbors.setdefault(b, []).append(a)

    # Find H atoms bonded to donor
    h_donors: list[tuple[int, int]] = []  # (h_idx, donor_idx)
    for i, atom in enumerate(mol.atoms):
        if atom.element.symbol != "H":
            continue
        for nb in neighbors.get(i, []):
            if mol.atoms[nb].element.symbol in _donors_acceptors:
                h_donors.append((i, nb))
                break

    # Find acceptor atoms
    acceptor_indices = [i for i, a in enumerate(mol.atoms) if a.element.symbol in _donors_acceptors]

    if not h_donors or not acceptor_indices:
        return []

    use_pbc = mol.pbc is not None
    hbonds = []

    for h_idx, d_idx in h_donors:
        h_pos = mol.atoms[h_idx].x
        d_pos = mol.atoms[d_idx].x
        for a_idx in acceptor_indices:
            if a_idx == d_idx or a_idx == h_idx:
                continue
            a_pos = mol.atoms[a_idx].x
            diff_da = a_pos - d_pos
            if use_pbc:
                frac = diff_da @ mol.pbc.reciprocal_basis_matrix
                frac -= np.round(frac)
                diff_da = frac @ mol.pbc.basis_matrix
            dist = float(np.linalg.norm(diff_da))
            if dist > d_cutoff:
                continue
            # Angle D-H...A
            dh = d_pos - h_pos
            ha = a_pos - h_pos
            if use_pbc:
                frac = ha @ mol.pbc.reciprocal_basis_matrix
                frac -= np.round(frac)
                ha = frac @ mol.pbc.basis_matrix
            cos_a = np.dot(dh, ha) / (np.linalg.norm(dh) * np.linalg.norm(ha) + 1e-10)
            angle = float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
            if angle >= angle_cutoff:
                hbonds.append((d_idx, h_idx, a_idx, dist, angle))

    return hbonds


def count_hbonds_trajectory(
    frames: list[Molecule],
    d_cutoff: float = 3.5,
    angle_cutoff: float = 120.0,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Count hydrogen bonds per frame over a trajectory.

    Returns:
        (frame_indices, hbond_counts) arrays.
    """
    n = len(frames)
    frame_idx = np.arange(n, dtype=float)
    counts = np.zeros(n)

    for fi, mol in enumerate(frames):
        if not mol.bonds:
            if fi > 0 and frames[0].bonds and len(frames[0].atoms) == len(mol.atoms):
                mol.bonds = frames[0].bonds
            else:
                mol.detect_bonds()
        hb = detect_hbonds(mol, d_cutoff, angle_cutoff)
        counts[fi] = len(hb)
        if progress_callback and fi % max(1, n // 20) == 0:
            progress_callback(fi / n)

    if progress_callback:
        progress_callback(1.0)

    return frame_idx, counts

# ======================================================================
# Module: gyration
# ======================================================================
"""Radius of gyration calculation."""


import numpy as np



def compute_gyration_radius(
    mol: Molecule,
    element: str | None = None,
    name: str | None = None,
) -> float:
    """Compute radius of gyration for selected atoms.

    Rg = sqrt( (1/N) * sum_i |r_i - r_com|^2 )
    """
    if element:
        if name and name != "(all)":
            atoms = [a for a in mol.atoms if a.element.symbol == element and a.name.strip() == name]
        else:
            atoms = [a for a in mol.atoms if a.element.symbol == element]
    else:
        atoms = list(mol.atoms)

    if len(atoms) < 2:
        return 0.0

    positions = np.array([a.x for a in atoms])
    com = positions.mean(axis=0)
    sq_dist = np.sum((positions - com) ** 2, axis=1)
    return float(np.sqrt(sq_dist.mean()))


def compute_gyration_trajectory(
    frames: list[Molecule],
    element: str | None = None,
    name: str | None = None,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute radius of gyration vs frame number.

    Returns:
        (frame_indices, rg_values) arrays.
    """
    if len(frames) < 1:
        return np.array([]), np.array([])

    n = len(frames)
    frame_idx = np.arange(n, dtype=float)
    rg = np.zeros(n)

    for fi, mol in enumerate(frames):
        rg[fi] = compute_gyration_radius(mol, element, name)
        if progress_callback and fi % max(1, n // 20) == 0:
            progress_callback(fi / n)

    if progress_callback:
        progress_callback(1.0)

    return frame_idx, rg

# ======================================================================
# Module: density3d
# ======================================================================
"""3D spatial density histogram for sorbate position analysis."""


import numpy as np



def compute_density_3d(
    frames: list[Molecule],
    element: str,
    n_bins: int = 20,
    name: str | None = None,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray, PBC]:
    """Compute 3D spatial density histogram for an element over trajectory frames.

    Bins atom positions into a 3D grid in fractional coordinates, then
    converts bin centers back to Cartesian for rendering.

    Args:
        frames: Trajectory frames.
        element: Element symbol to histogram.
        n_bins: Number of bins per axis.
        name: Optional atom name filter.
        progress_callback: Optional callable(frac).

    Returns:
        (positions, densities, pbc):
        positions is (N, 3) Cartesian coords of non-zero bins,
        densities is (N,) normalized density values (0-1),
        pbc is the PBC from the first frame.
    """
    if not frames or frames[0].pbc is None:
        return np.zeros((0, 3)), np.zeros(0), None

    pbc = frames[0].pbc
    histogram = np.zeros((n_bins, n_bins, n_bins))

    for fi, mol in enumerate(frames):
        if mol.pbc is None:
            continue
        for atom in mol.atoms:
            if atom.element.symbol != element:
                continue
            if name and name != "(all)" and atom.name.strip() != name:
                continue
            # Convert to fractional coordinates
            frac = atom.x @ mol.pbc.reciprocal_basis_matrix
            frac -= np.floor(frac)  # wrap to [0, 1)
            ix = min(int(frac[0] * n_bins), n_bins - 1)
            iy = min(int(frac[1] * n_bins), n_bins - 1)
            iz = min(int(frac[2] * n_bins), n_bins - 1)
            histogram[ix, iy, iz] += 1

        if progress_callback and fi % max(1, len(frames) // 20) == 0:
            progress_callback(fi / len(frames))

    if progress_callback:
        progress_callback(1.0)

    # Extract non-zero bins
    nonzero = histogram > 0
    if not nonzero.any():
        return np.zeros((0, 3)), np.zeros(0), pbc

    indices = np.argwhere(nonzero)  # (N, 3) array of bin indices
    densities = histogram[nonzero]

    # Normalize to 0-1
    max_dens = densities.max()
    if max_dens > 0:
        densities = densities / max_dens

    # Convert bin centers to Cartesian
    frac_centers = (indices + 0.5) / n_bins  # fractional coordinates
    positions = frac_centers @ pbc.basis_matrix

    return positions, densities, pbc

# ======================================================================
# Module: templates
# ======================================================================
"""Template-based input file generator.

Templates live in src/pdb_wizard/templates/ as plain text files with
{{{varname}}} placeholders. Substitution is intentionally minimal — no
conditionals, no loops, no inheritance. For computed/multi-line content
(like CP2K's per-element &KIND blocks) the generator builds the block in
Python and passes it in as a single string variable.

Users may edit the shipped templates in place to tweak defaults; the
Input Generator UI re-reads them every time the modal opens.
"""



TEMPLATES_DIR = Path(__file__).parent

_PLACEHOLDER = re.compile(r"\{\{\{(\w+)\}\}\}")

# Populated by build.py for the amalgamated single-file dist (which has no
# templates/ folder next to it). Source installs always read from disk so
# users can edit the shipped templates in place.
_BUNDLED_TEMPLATES: dict[str, str] = {}


def load_template(name: str) -> str:
    """Load a template by filename (relative to templates/).

    Prefers an on-disk file under TEMPLATES_DIR so editable installs can
    customise the shipped templates. Falls back to _BUNDLED_TEMPLATES when
    running from the amalgamated single-file dist.
    """
    p = TEMPLATES_DIR / name
    if p.is_file():
        return p.read_text()
    if name in _BUNDLED_TEMPLATES:
        return _BUNDLED_TEMPLATES[name]
    raise FileNotFoundError(f"template not found: {p}")


def find_variables(template_str: str) -> list[str]:
    """Return the unique {{{var}}} names in `template_str`, in first-seen order."""
    seen: dict[str, None] = {}
    for m in _PLACEHOLDER.finditer(template_str):
        seen.setdefault(m.group(1), None)
    return list(seen.keys())


def render(template_str: str, vars: dict, *, strict: bool = True) -> str:
    """Substitute {{{var}}} placeholders with `str(vars[var])`.

    Raises KeyError when strict=True and a placeholder has no value.
    Extra keys in `vars` are silently ignored.
    """
    missing: list[str] = []

    def _sub(m: re.Match) -> str:
        name = m.group(1)
        if name not in vars:
            missing.append(name)
            return m.group(0)
        return str(vars[name])

    out = _PLACEHOLDER.sub(_sub, template_str)
    if missing and strict:
        raise KeyError(
            "unfilled template variables: " + ", ".join(sorted(set(missing)))
        )
    return out


# --- Engine registry (UI metadata) ---


@dataclass
class VarSpec:
    """One user-editable variable in an engine form.

    kind:
        'float'       — SpinBox over reals
        'int'         — SpinBox over ints
        'str'         — single-line Input
        'select'      — Select with `choices`
        'csv_floats'  — Input parsed as comma-sep list of floats (multi_run)
    """
    name: str
    label: str
    kind: str
    default: Any
    choices: tuple[str, ...] | None = None
    min_val: float | None = None
    max_val: float | None = None
    step: float | None = None
    help: str = ""


@dataclass
class EngineSpec:
    """One engine: template + UI form + behavior."""
    key: str
    label: str
    description: str
    template_file: str
    output_filename: str
    user_vars: list[VarSpec] = field(default_factory=list)
    # multi_run=True iff the user picks a 'csv_floats' variable that
    # produces one output directory per value (currently only MPMC isotherm
    # over pressures).
    multi_run: bool = False
    multi_run_var: str = ""   # name of the csv_floats var that drives the loop


ENGINES: dict[str, EngineSpec] = {
    "mpmc_isotherm": EngineSpec(
        key="mpmc_isotherm",
        label="MPMC μVT isotherm",
        description="Grand-canonical Monte Carlo isotherm over multiple pressures.",
        template_file="mpmc_uvt.inp",
        output_filename="mpmc.inp",
        multi_run=True,
        multi_run_var="pressures",
        user_vars=[
            VarSpec("sorbate", "Sorbate model", "select", "h2_bss",
                    choices=("h2_bss", "ch4_trappe", "n2_ttm", "co2_epm2")),
            VarSpec("temperature", "Temperature (K)", "float", 77.0,
                    min_val=1.0, max_val=2000.0),
            VarSpec("pressures", "Pressures (atm, comma-sep)", "csv_floats",
                    "0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0"),
            VarSpec("numsteps", "MC steps per point", "int", 500_000,
                    min_val=1000, max_val=100_000_000, step=100_000),
            VarSpec("corrtime", "Correlation time", "int", 1000,
                    min_val=10, max_val=100_000),
            VarSpec("ensemble", "Ensemble", "select", "uvt",
                    choices=("uvt", "nvt")),
            VarSpec("insert_probability", "Insert probability", "float", 0.667,
                    min_val=0.0, max_val=1.0),
            VarSpec("move_factor", "Move factor", "float", 0.01,
                    min_val=0.0, max_val=1.0),
            VarSpec("rot_factor", "Rotation factor", "float", 0.01,
                    min_val=0.0, max_val=1.0),
        ],
    ),
    "cp2k_cell_opt": EngineSpec(
        key="cp2k_cell_opt",
        label="CP2K cell optimization",
        description="Relax atomic positions and the unit cell with CP2K LBFGS.",
        template_file="cp2k_cell_opt.inp",
        output_filename="input.inp",
        user_vars=[
            VarSpec("project", "Project name", "str", "pdbwizard"),
            VarSpec("cutoff_ry", "Plane-wave cutoff (Ry)", "float", 400.0,
                    min_val=100.0, max_val=2000.0, step=10.0),
            VarSpec("rel_cutoff_ry", "Relative cutoff (Ry)", "float", 60.0,
                    min_val=20.0, max_val=200.0, step=5.0),
            VarSpec("xc_functional", "XC functional", "select", "PBE",
                    choices=("PBE", "BLYP", "PBE0", "BEEFVDW")),
            VarSpec("eps_scf", "SCF tolerance", "float", 1.0e-7,
                    min_val=1.0e-12, max_val=1.0e-3),
            VarSpec("max_scf", "Max SCF iters", "int", 50,
                    min_val=5, max_val=500),
            VarSpec("n_opt_steps", "Cell-opt steps", "int", 200,
                    min_val=1, max_val=5000),
        ],
    ),
    "cp2k_dft_md": EngineSpec(
        key="cp2k_dft_md",
        label="CP2K DFT-MD (NVT)",
        description="Born-Oppenheimer MD with Nose-Hoover thermostat.",
        template_file="cp2k_dft_md.inp",
        output_filename="input.inp",
        user_vars=[
            VarSpec("project", "Project name", "str", "pdbwizard"),
            VarSpec("temperature", "Temperature (K)", "float", 300.0,
                    min_val=1.0, max_val=5000.0),
            VarSpec("timestep_fs", "Timestep (fs)", "float", 0.5,
                    min_val=0.1, max_val=10.0, step=0.1),
            VarSpec("n_md_steps", "MD steps", "int", 10_000,
                    min_val=1, max_val=10_000_000, step=1000),
            VarSpec("md_traj_interval", "Trajectory print every N steps",
                    "int", 10, min_val=1, max_val=10_000),
            VarSpec("cutoff_ry", "Plane-wave cutoff (Ry)", "float", 400.0,
                    min_val=100.0, max_val=2000.0, step=10.0),
            VarSpec("rel_cutoff_ry", "Relative cutoff (Ry)", "float", 60.0,
                    min_val=20.0, max_val=200.0, step=5.0),
            VarSpec("xc_functional", "XC functional", "select", "PBE",
                    choices=("PBE", "BLYP", "PBE0", "BEEFVDW")),
            VarSpec("eps_scf", "SCF tolerance", "float", 1.0e-7,
                    min_val=1.0e-12, max_val=1.0e-3),
            VarSpec("max_scf", "Max SCF iters", "int", 50,
                    min_val=5, max_val=500),
        ],
    ),
    "openmm_npt": EngineSpec(
        key="openmm_npt",
        label="OpenMM NPT MD",
        description="Langevin NPT with Monte-Carlo barostat, amber14/tip3pfb default.",
        template_file="openmm_npt.py",
        output_filename="run.py",
        user_vars=[
            VarSpec("forcefield", "Force field XMLs (comma-sep)", "str",
                    "amber14-all.xml,amber14/tip3pfb.xml"),
            VarSpec("temperature_k", "Temperature (K)", "float", 300.0,
                    min_val=1.0, max_val=2000.0),
            VarSpec("pressure_atm", "Pressure (atm)", "float", 1.0,
                    min_val=0.001, max_val=10000.0),
            VarSpec("timestep_fs", "Timestep (fs)", "float", 2.0,
                    min_val=0.1, max_val=10.0, step=0.5),
            VarSpec("n_steps", "MD steps", "int", 100_000,
                    min_val=10, max_val=100_000_000, step=10_000),
            VarSpec("report_interval", "Report every N steps", "int", 1000,
                    min_val=10, max_val=1_000_000, step=100),
            VarSpec("nonbonded_cutoff_nm", "Nonbonded cutoff (nm)", "float",
                    1.0, min_val=0.5, max_val=5.0, step=0.1),
        ],
    ),
    "openmm_nvt": EngineSpec(
        key="openmm_nvt",
        label="OpenMM NVT MD",
        description="Langevin NVT, amber14/tip3pfb default.",
        template_file="openmm_nvt.py",
        output_filename="run.py",
        user_vars=[
            VarSpec("forcefield", "Force field XMLs (comma-sep)", "str",
                    "amber14-all.xml,amber14/tip3pfb.xml"),
            VarSpec("temperature_k", "Temperature (K)", "float", 300.0,
                    min_val=1.0, max_val=2000.0),
            VarSpec("timestep_fs", "Timestep (fs)", "float", 2.0,
                    min_val=0.1, max_val=10.0, step=0.5),
            VarSpec("n_steps", "MD steps", "int", 100_000,
                    min_val=10, max_val=100_000_000, step=10_000),
            VarSpec("report_interval", "Report every N steps", "int", 1000,
                    min_val=10, max_val=1_000_000, step=100),
            VarSpec("nonbonded_cutoff_nm", "Nonbonded cutoff (nm)", "float",
                    1.0, min_val=0.5, max_val=5.0, step=0.1),
        ],
    ),
}


def engine_defaults(key: str) -> dict[str, Any]:
    """Return {var_name: default} for the engine identified by `key`."""
    spec = ENGINES[key]
    return {v.name: v.default for v in spec.user_vars}

# Bundled template payload (inlined by build.py)
_BUNDLED_TEMPLATES.update({
    'cp2k_cell_opt.inp': '! CP2K cell optimization — generated by pdb_wizard from cp2k_cell_opt.inp\n! Placeholders are triple-brace tokens. The per-element &KIND block\n! (basis + GTH pseudopotential) is built by the generator and slotted in\n! as a single multi-line value.\n\n&GLOBAL\n  PROJECT {{{project}}}\n  RUN_TYPE CELL_OPT\n  PRINT_LEVEL LOW\n&END GLOBAL\n\n&FORCE_EVAL\n  METHOD QUICKSTEP\n  STRESS_TENSOR ANALYTICAL\n  &SUBSYS\n    &CELL\n      A   {{{cell_a}}} 0.000000 0.000000\n      B   0.000000 {{{cell_b}}} 0.000000\n      C   0.000000 0.000000 {{{cell_c}}}\n      ALPHA_BETA_GAMMA {{{cell_alpha}}} {{{cell_beta}}} {{{cell_gamma}}}\n      PERIODIC XYZ\n    &END CELL\n    &TOPOLOGY\n      COORD_FILE_NAME {{{coords_file}}}\n      COORD_FILE_FORMAT XYZ\n      &CENTER_COORDINATES\n      &END\n    &END TOPOLOGY\n{{{kind_blocks}}}  &END SUBSYS\n  &DFT\n    BASIS_SET_FILE_NAME BASIS_MOLOPT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    &MGRID\n      CUTOFF {{{cutoff_ry}}}\n      REL_CUTOFF {{{rel_cutoff_ry}}}\n    &END MGRID\n    &SCF\n      EPS_SCF {{{eps_scf}}}\n      MAX_SCF {{{max_scf}}}\n      &OT\n        MINIMIZER DIIS\n        PRECONDITIONER FULL_SINGLE_INVERSE\n      &END OT\n      &OUTER_SCF\n        EPS_SCF {{{eps_scf}}}\n        MAX_SCF 20\n      &END\n    &END SCF\n    &XC\n      &XC_FUNCTIONAL {{{xc_functional}}}\n      &END\n    &END XC\n  &END DFT\n&END FORCE_EVAL\n\n&MOTION\n  &CELL_OPT\n    OPTIMIZER LBFGS\n    KEEP_ANGLES T\n    MAX_ITER {{{n_opt_steps}}}\n    TYPE DIRECT_CELL_OPT\n  &END CELL_OPT\n  &PRINT\n    &TRAJECTORY\n      &EACH\n        CELL_OPT 1\n      &END\n    &END\n    &CELL\n      &EACH\n        CELL_OPT 1\n      &END\n    &END\n  &END\n&END MOTION\n',
    'cp2k_dft_md.inp': '! CP2K Born-Oppenheimer NVT MD — generated by pdb_wizard from cp2k_dft_md.inp\n! Placeholders are triple-brace tokens. The per-element &KIND block\n! (basis + GTH pseudopotential) is built by the generator and slotted in\n! as a single multi-line value.\n\n&GLOBAL\n  PROJECT {{{project}}}\n  RUN_TYPE MD\n  PRINT_LEVEL LOW\n&END GLOBAL\n\n&FORCE_EVAL\n  METHOD QUICKSTEP\n  STRESS_TENSOR ANALYTICAL\n  &SUBSYS\n    &CELL\n      A   {{{cell_a}}} 0.000000 0.000000\n      B   0.000000 {{{cell_b}}} 0.000000\n      C   0.000000 0.000000 {{{cell_c}}}\n      ALPHA_BETA_GAMMA {{{cell_alpha}}} {{{cell_beta}}} {{{cell_gamma}}}\n      PERIODIC XYZ\n    &END CELL\n    &TOPOLOGY\n      COORD_FILE_NAME {{{coords_file}}}\n      COORD_FILE_FORMAT XYZ\n      &CENTER_COORDINATES\n      &END\n    &END TOPOLOGY\n{{{kind_blocks}}}  &END SUBSYS\n  &DFT\n    BASIS_SET_FILE_NAME BASIS_MOLOPT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    &MGRID\n      CUTOFF {{{cutoff_ry}}}\n      REL_CUTOFF {{{rel_cutoff_ry}}}\n    &END MGRID\n    &SCF\n      EPS_SCF {{{eps_scf}}}\n      MAX_SCF {{{max_scf}}}\n      &OT\n        MINIMIZER DIIS\n        PRECONDITIONER FULL_SINGLE_INVERSE\n      &END OT\n      &OUTER_SCF\n        EPS_SCF {{{eps_scf}}}\n        MAX_SCF 20\n      &END\n    &END SCF\n    &XC\n      &XC_FUNCTIONAL {{{xc_functional}}}\n      &END\n    &END XC\n  &END DFT\n&END FORCE_EVAL\n\n&MOTION\n  &MD\n    ENSEMBLE NVT\n    TIMESTEP {{{timestep_fs}}}\n    STEPS {{{n_md_steps}}}\n    TEMPERATURE {{{temperature}}}\n    &THERMOSTAT\n      TYPE NOSE\n      &NOSE\n        TIMECON 50\n      &END\n    &END\n  &END MD\n  &PRINT\n    &TRAJECTORY\n      FILENAME =trajectory.xyz\n      FORMAT XYZ\n      &EACH\n        MD {{{md_traj_interval}}}\n      &END\n    &END\n    &VELOCITIES OFF\n    &END\n    &RESTART\n      FILENAME =restart\n      &EACH\n        MD 100\n      &END\n    &END\n  &END\n&END MOTION\n',
    'mpmc_uvt.inp': "! MPMC isotherm input — generated by pdb_wizard from template mpmc_uvt.inp\n! Placeholders are written as triple-brace tokens; the Input Generator\n! UI re-renders this file with the values you choose.\n\njob_name              {{{job_name}}}\n\nensemble              {{{ensemble}}}\n\nabcbasis              {{{cell_a}}} {{{cell_b}}} {{{cell_c}}} {{{cell_alpha}}} {{{cell_beta}}} {{{cell_gamma}}}\n\ntemperature           {{{temperature}}}\npressure              {{{pressure}}}\n\nnumsteps              {{{numsteps}}}\ncorrtime              {{{corrtime}}}\ninsert_probability    {{{insert_probability}}}\nmove_factor           {{{move_factor}}}\nrot_factor            {{{rot_factor}}}\n\npqr_input             {{{pqr_input}}}\npqr_output            output.pqr\npqr_restart           restart.pqr\ntraj_output           trajectory.pqr\nenergy_output         energy.dat\n\n! No insert_input: the sorbate is already embedded in pqr_input as a movable\n! molecule, which MPMC clones for insertions. Pointing insert_input at a\n! separate pool triggers a double-free in MPMC's cleanup (crash on exit).\n\nrd_lrc                on\nwrapall               on\n\npop_histogram         on\npop_histogram_output  histogram.dat\n",
    'openmm_npt.py': '"""OpenMM NPT driver — generated by pdb_wizard from openmm_npt.py template.\n\nPlaceholders are triple-brace tokens. The generator fills in force-field\nXMLs, temperature, pressure, step counts, etc. Edit the rendered run.py\nto tune the integrator, platform, or reporter settings.\n\nRun:\n    python {{{script_filename}}}\n\nOutputs:\n    trajectory.dcd  — coords every {{{report_interval}}} steps\n    output.log      — energy / temperature / box (CSV)\n    final.pdb       — last frame\n"""\nfrom openmm import (\n    LangevinMiddleIntegrator,\n    MonteCarloBarostat,\n    Platform,\n)\nfrom openmm.app import (\n    PME,\n    DCDReporter,\n    ForceField,\n    HBonds,\n    Modeller,\n    PDBFile,\n    PDBReporter,\n    Simulation,\n    StateDataReporter,\n)\nfrom openmm.unit import (\n    angstrom,\n    atmosphere,\n    femtoseconds,\n    kelvin,\n    nanometer,\n    picosecond,\n    picoseconds,\n)\nimport sys\n\n\npdb = PDBFile("{{{pdb_filename}}}")\nforcefield = ForceField({{{forcefield_xmls}}})\n\nmodeller = Modeller(pdb.topology, pdb.positions)\n# Add hydrogens / solvent here if needed:\n#   modeller.addHydrogens(forcefield)\n#   modeller.addSolvent(forcefield, padding=1.0 * nanometer)\n\nsystem = forcefield.createSystem(\n    modeller.topology,\n    nonbondedMethod=PME,\n    nonbondedCutoff={{{nonbonded_cutoff_nm}}} * nanometer,\n    constraints=HBonds,\n)\n\n# NPT barostat\nsystem.addForce(MonteCarloBarostat(\n    {{{pressure_atm}}} * atmosphere,\n    {{{temperature_k}}} * kelvin,\n    25,\n))\n\nintegrator = LangevinMiddleIntegrator(\n    {{{temperature_k}}} * kelvin,\n    1.0 / picosecond,\n    {{{timestep_fs}}} * femtoseconds,\n)\n\nplatform = Platform.getPlatformByName("CUDA")\ntry:\n    simulation = Simulation(modeller.topology, system, integrator, platform)\nexcept Exception:\n    platform = Platform.getPlatformByName("CPU")\n    simulation = Simulation(modeller.topology, system, integrator, platform)\n\nsimulation.context.setPositions(modeller.positions)\nsimulation.minimizeEnergy()\n\nsimulation.reporters.append(DCDReporter("trajectory.dcd", {{{report_interval}}}))\nsimulation.reporters.append(PDBReporter("final.pdb", {{{n_steps}}}))\nsimulation.reporters.append(StateDataReporter(\n    sys.stdout, {{{report_interval}}},\n    step=True, time=True, potentialEnergy=True, kineticEnergy=True,\n    totalEnergy=True, temperature=True, volume=True, density=True,\n    speed=True, separator="  ",\n))\nsimulation.reporters.append(StateDataReporter(\n    "output.log", {{{report_interval}}},\n    step=True, time=True, potentialEnergy=True, kineticEnergy=True,\n    totalEnergy=True, temperature=True, volume=True, density=True,\n    separator=",",\n))\n\nsimulation.step({{{n_steps}}})\nprint("done")\n',
    'openmm_nvt.py': '"""OpenMM NVT driver — generated by pdb_wizard from openmm_nvt.py template.\n\nPlaceholders are triple-brace tokens. The {{{nonbonded_setup}}} block is\nbuilt by the generator: PME with a cutoff if the system has PBC, otherwise\nNoCutoff.\n\nRun:\n    python {{{script_filename}}}\n\nOutputs:\n    trajectory.dcd  — coords every {{{report_interval}}} steps\n    output.log      — energy / temperature (CSV)\n    final.pdb       — last frame\n"""\nfrom openmm import (\n    LangevinMiddleIntegrator,\n    Platform,\n)\nfrom openmm.app import (\n    NoCutoff,\n    PME,\n    DCDReporter,\n    ForceField,\n    HBonds,\n    Modeller,\n    PDBFile,\n    PDBReporter,\n    Simulation,\n    StateDataReporter,\n)\nfrom openmm.unit import (\n    angstrom,\n    femtoseconds,\n    kelvin,\n    nanometer,\n    picosecond,\n    picoseconds,\n)\nimport sys\n\n\npdb = PDBFile("{{{pdb_filename}}}")\nforcefield = ForceField({{{forcefield_xmls}}})\n\nmodeller = Modeller(pdb.topology, pdb.positions)\n\nsystem = forcefield.createSystem(\n    modeller.topology,\n{{{nonbonded_setup}}}    constraints=HBonds,\n)\n\nintegrator = LangevinMiddleIntegrator(\n    {{{temperature_k}}} * kelvin,\n    1.0 / picosecond,\n    {{{timestep_fs}}} * femtoseconds,\n)\n\nplatform = Platform.getPlatformByName("CUDA")\ntry:\n    simulation = Simulation(modeller.topology, system, integrator, platform)\nexcept Exception:\n    platform = Platform.getPlatformByName("CPU")\n    simulation = Simulation(modeller.topology, system, integrator, platform)\n\nsimulation.context.setPositions(modeller.positions)\nsimulation.minimizeEnergy()\n\nsimulation.reporters.append(DCDReporter("trajectory.dcd", {{{report_interval}}}))\nsimulation.reporters.append(PDBReporter("final.pdb", {{{n_steps}}}))\nsimulation.reporters.append(StateDataReporter(\n    sys.stdout, {{{report_interval}}},\n    step=True, time=True, potentialEnergy=True, kineticEnergy=True,\n    totalEnergy=True, temperature=True, volume=True, density=True,\n    speed=True, separator="  ",\n))\nsimulation.reporters.append(StateDataReporter(\n    "output.log", {{{report_interval}}},\n    step=True, time=True, potentialEnergy=True, kineticEnergy=True,\n    totalEnergy=True, temperature=True, volume=True, density=True,\n    separator=",",\n))\n\nsimulation.step({{{n_steps}}})\nprint("done")\n',
})

# ======================================================================
# Module: sim_inputs
# ======================================================================
"""Generate input scripts for external simulation packages.

Currently supported:
  * OpenMM (Python script): NPT or NVT MD with TIP3P / amber14 by default
  * CP2K (input file): cell optimization and DFT-MD

All input files are rendered from templates under `src/pdb_wizard/templates/`
using the placeholder syntax `{{{name}}}`. This module is responsible only
for (1) writing the companion coordinate file, (2) computing derived
template variables from the loaded Molecule, and (3) invoking the
template renderer.
"""



try:  # Literal moved into typing in 3.8; typing_extensions backports it for 3.7
    from typing import Literal
except ImportError:  # pragma: no cover - exercised only on Python 3.7
    from typing_extensions import Literal


# ---------------------------------------------------------------------------
# OpenMM
# ---------------------------------------------------------------------------

OpenMMEnsemble = Literal["NPT", "NVT"]


def _ff_python_literal(forcefield: str) -> str:
    """Convert 'amber14-all.xml,amber14/tip3pfb.xml' into a Python argument
    list: `'amber14-all.xml', 'amber14/tip3pfb.xml'` — suitable for splicing
    into `ForceField({...})` in the template."""
    return ", ".join(repr(x.strip()) for x in forcefield.split(",") if x.strip())


def _nvt_nonbonded_setup(has_pbc: bool, cutoff_nm: float) -> str:
    """Build the createSystem() nonbonded kwargs block for NVT.

    With PBC we use PME and pass a cutoff. Without PBC, OpenMM's NoCutoff
    method doesn't accept `nonbondedCutoff`, so we emit just the method.
    """
    if has_pbc:
        return (
            f"    nonbondedMethod=PME,\n"
            f"    nonbondedCutoff={cutoff_nm} * nanometer,\n"
        )
    return "    nonbondedMethod=NoCutoff,\n"


def generate_openmm_script(
    mol: Molecule,
    output_dir: str | Path,
    ensemble: OpenMMEnsemble = "NPT",
    *,
    forcefield: str = "amber14-all.xml,amber14/tip3pfb.xml",
    temperature_K: float = 300.0,
    pressure_atm: float = 1.0,
    timestep_fs: float = 2.0,
    n_steps: int = 100_000,
    report_interval: int = 1000,
    nonbonded_cutoff_nm: float = 1.0,
    pdb_filename: str = "system.pdb",
    script_filename: str = "run.py",
) -> dict:
    """Write an OpenMM Python driver + a PDB topology file.

    Returns a dict with the paths written.
    """
    if ensemble not in ("NPT", "NVT"):
        raise ValueError(f"ensemble must be NPT or NVT, got {ensemble!r}")
    if mol.pbc is None and ensemble == "NPT":
        raise ValueError("NPT simulation requires a periodic box (mol.pbc).")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # PDB topology — skip_mols_step=False so atoms group into named residues
    # (HOH, MET, etc.), which amber-style force fields need to find templates.
    pdb_path = out / pdb_filename
    write_standard_pdb(mol, str(pdb_path), skip_mols_step=False)

    common_vars = {
        "pdb_filename": pdb_filename,
        "script_filename": script_filename,
        "forcefield_xmls": _ff_python_literal(forcefield),
        "temperature_k": temperature_K,
        "timestep_fs": timestep_fs,
        "n_steps": n_steps,
        "report_interval": report_interval,
        "nonbonded_cutoff_nm": nonbonded_cutoff_nm,
    }

    if ensemble == "NPT":
        template = load_template("openmm_npt.py")
        rendered = render(template, {
            **common_vars,
            "pressure_atm": pressure_atm,
        })
    else:  # NVT
        template = load_template("openmm_nvt.py")
        rendered = render(template, {
            **common_vars,
            "nonbonded_setup": _nvt_nonbonded_setup(
                mol.pbc is not None, nonbonded_cutoff_nm,
            ),
        })

    script_path = out / script_filename
    script_path.write_text(rendered)

    return {
        "script": str(script_path),
        "pdb": str(pdb_path),
        "ensemble": ensemble,
        "n_steps": n_steps,
        "temperature_K": temperature_K,
        "pressure_atm": pressure_atm if ensemble == "NPT" else None,
    }


# ---------------------------------------------------------------------------
# CP2K
# ---------------------------------------------------------------------------

CP2KMode = Literal["cell_opt", "dft_md"]


def _formula_summary(mol: Molecule) -> str:
    counts = Counter(a.element.symbol for a in mol.atoms)
    return "".join(f"{el}{n}" if n > 1 else el for el, n in sorted(counts.items()))


def _kind_section(symbol: str, basis: str = "DZVP-MOLOPT-SR-GTH",
                  potential: str = "GTH-PBE") -> str:
    """Build a CP2K &KIND block for one element. Tries to choose a sensible
    Q-value (number of valence electrons) for the GTH pseudopotential; falls
    back to '-q' suffix omitted (CP2K will pick the default)."""
    Q = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7,
        "Ne": 8, "Na": 9, "Mg": 10, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7,
        "Ar": 8, "K": 9, "Ca": 10, "Ti": 12, "V": 13, "Cr": 14, "Mn": 15,
        "Fe": 16, "Co": 17, "Ni": 18, "Cu": 11, "Zn": 12, "Ga": 13, "Br": 7,
        "I": 7, "Pt": 18, "Au": 11,
    }
    pot = potential
    if symbol in Q:
        pot = f"{potential}-q{Q[symbol]}"
    return (
        f"    &KIND {symbol}\n"
        f"      BASIS_SET {basis}\n"
        f"      POTENTIAL {pot}\n"
        f"    &END KIND\n"
    )


def _build_kind_blocks(mol: Molecule) -> str:
    elements_present = sorted({a.element.symbol for a in mol.atoms})
    return "".join(_kind_section(el) for el in elements_present)


def _write_xyz(mol: Molecule, path: Path) -> None:
    lines = [str(len(mol.atoms)), f"{_formula_summary(mol)} from pdb_wizard"]
    for a in mol.atoms:
        lines.append(
            f"{a.element.symbol}  {a.x[0]:.6f}  {a.x[1]:.6f}  {a.x[2]:.6f}"
        )
    path.write_text("\n".join(lines) + "\n")


def generate_cp2k_input(
    mol: Molecule,
    output_dir: str | Path,
    mode: CP2KMode = "cell_opt",
    *,
    project_name: str = "pdbwizard",
    cutoff_Ry: float = 400.0,
    rel_cutoff_Ry: float = 60.0,
    xc_functional: str = "PBE",
    eps_scf: float = 1.0e-7,
    max_scf: int = 50,
    temperature_K: float = 300.0,
    timestep_fs: float = 0.5,
    n_md_steps: int = 10_000,
    n_opt_steps: int = 200,
    md_traj_interval: int = 10,
    coords_filename: str = "system.xyz",
    input_filename: str = "input.inp",
) -> dict:
    """Write a CP2K input file (`.inp`) plus an XYZ coordinate file.

    `mode='cell_opt'` runs CP2K's CELL_OPT to relax both atomic positions
    and the unit cell. `mode='dft_md'` runs Born-Oppenheimer MD in NVT.
    """
    if mode not in ("cell_opt", "dft_md"):
        raise ValueError(f"mode must be cell_opt or dft_md, got {mode!r}")
    if mol.pbc is None:
        raise ValueError("CP2K requires a periodic cell (mol.pbc).")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    xyz_path = out / coords_filename
    _write_xyz(mol, xyz_path)

    pbc = mol.pbc
    common_vars = {
        "project": project_name,
        "cell_a": f"{pbc.a:.6f}",
        "cell_b": f"{pbc.b:.6f}",
        "cell_c": f"{pbc.c:.6f}",
        "cell_alpha": f"{pbc.alpha:.4f}",
        "cell_beta": f"{pbc.beta:.4f}",
        "cell_gamma": f"{pbc.gamma:.4f}",
        "coords_file": coords_filename,
        "kind_blocks": _build_kind_blocks(mol),
        "cutoff_ry": cutoff_Ry,
        "rel_cutoff_ry": rel_cutoff_Ry,
        "xc_functional": xc_functional,
        "eps_scf": eps_scf,
        "max_scf": max_scf,
    }

    if mode == "cell_opt":
        template = load_template("cp2k_cell_opt.inp")
        rendered = render(template, {
            **common_vars,
            "n_opt_steps": n_opt_steps,
        })
    else:  # dft_md
        template = load_template("cp2k_dft_md.inp")
        rendered = render(template, {
            **common_vars,
            "temperature": temperature_K,
            "timestep_fs": timestep_fs,
            "n_md_steps": n_md_steps,
            "md_traj_interval": md_traj_interval,
        })

    inp_path = out / input_filename
    inp_path.write_text(rendered)

    return {
        "input": str(inp_path),
        "xyz": str(xyz_path),
        "mode": mode,
        "project": project_name,
        "elements": sorted({a.element.symbol for a in mol.atoms}),
    }

# ======================================================================
# Module: isotherm
# ======================================================================
"""Isotherm planner — generates MPMC input directories at multiple pressures.

The actual MPMC input file is rendered from `templates/mpmc_uvt.inp` via the
template engine; this module supplies the per-pressure derived variables
(box size, sorbate placement, file paths) and writes the supporting files
(input.pqr, insert.pqr, run.sh).
"""



import numpy as np


# Default MPMC LJ cutoff in Å. The simulation needs min(a,b,c)/2 > this.
_DEFAULT_LJ_CUTOFF = 10.0
# A charge is "set" if |q| exceeds this — anything smaller is rounding noise.
_CHARGE_NONZERO = 1e-9
# An LJ parameter is "set" if epsilon > this — typical organic atom epsilon
# is 30-150 K, hydrogen is ~10 K, so 0.1 K is a safe lower bound for "applied".
_EPSILON_SET = 0.1
# Sorbate–framework overlap distance (Å) below which we warn.
_SORBATE_OVERLAP = 1.5


def validate_for_mpmc(
    mol: Molecule,
    sorbate_pos: tuple[float, float, float] | None = None,
    lj_cutoff: float = _DEFAULT_LJ_CUTOFF,
) -> tuple[list[str], list[str]]:
    """Inspect a molecule for common MPMC-input footguns.

    Returns (errors, warnings). 'errors' block generation outright; 'warnings'
    let the user proceed but flag likely issues.

    Errors:
      - No PBC (MPMC requires a periodic box)
      - All charges effectively zero (zero electrostatics — meaningless)
      - All epsilon/sigma effectively zero (zero LJ — sorbate floats freely)

    Warnings:
      - |net charge| > 1e-3 e (non-neutral box)
      - min(a,b,c)/2 < lj_cutoff (box too small for default LJ cutoff)
      - Sorbate placement overlaps a framework atom (within 1.5 Å)
      - Triclinic cell (α/β/γ deviates from 90° — abcbasis loses info)
    """
    errors: list[str] = []
    warnings: list[str] = []

    if mol.pbc is None:
        errors.append("Molecule has no PBC. MPMC requires a periodic cell — "
                      "use 'Update Unit Cell' to add one.")
        return errors, warnings  # other checks need PBC

    n = len(mol.atoms)
    if n == 0:
        errors.append("Molecule has no atoms.")
        return errors, warnings

    # Charges
    charges = np.array([a.charge for a in mol.atoms])
    n_set = int((np.abs(charges) > _CHARGE_NONZERO).sum())
    if n_set == 0:
        errors.append(
            "All atom charges are zero. Apply a force field or "
            "Generate QEq Charges before planning the isotherm — otherwise "
            "MPMC computes no electrostatics and the uptake is meaningless."
        )
    else:
        net_q = float(charges.sum())
        if abs(net_q) > 1e-3:
            warnings.append(
                f"Net system charge = {net_q:+.4f} e (non-neutral). "
                f"MPMC assumes a neutral simulation cell — add counterions "
                f"or rescale charges before running."
            )

    # LJ parameters
    epsilons = np.array([a.epsilon for a in mol.atoms])
    n_lj_set = int((epsilons > _EPSILON_SET).sum())
    if n_lj_set == 0:
        errors.append(
            "All Lennard-Jones epsilon values are zero. Apply a force field "
            "(OPLS-AA or PHAHST) before planning — without LJ, MPMC sees no "
            "framework potential and the sorbate just floats."
        )
    elif n_lj_set < n // 2:
        warnings.append(
            f"Only {n_lj_set}/{n} atoms have nonzero LJ epsilon. "
            f"Some elements are likely missing from the chosen force field."
        )

    # Box size vs LJ cutoff
    pbc = mol.pbc
    half_min = min(pbc.a, pbc.b, pbc.c) / 2.0
    if half_min < lj_cutoff:
        warnings.append(
            f"Smallest box half-length is {half_min:.2f} Å but MPMC's default "
            f"LJ cutoff is {lj_cutoff:.1f} Å. Either build a supercell "
            f"(Edit > Extend Axis) or reduce the cutoff in mpmc.inp."
        )

    # Triclinic cell
    if (abs(pbc.alpha - 90.0) > 0.5
            or abs(pbc.beta - 90.0) > 0.5
            or abs(pbc.gamma - 90.0) > 0.5):
        warnings.append(
            f"Cell is not orthorhombic (α={pbc.alpha:.2f}, β={pbc.beta:.2f}, "
            f"γ={pbc.gamma:.2f}). 'abcbasis' in mpmc.inp loses skew "
            f"information; use 'pbc_input' with the basis matrix instead."
        )

    # Sorbate placement overlap
    if sorbate_pos is None:
        cx, cy, cz = mol.center()
    else:
        cx, cy, cz = sorbate_pos
    coords = np.array([a.x for a in mol.atoms])
    dx = coords - np.array([float(cx), float(cy), float(cz)])
    # Min-image distance to the sorbate placement
    frac = dx @ pbc.reciprocal_basis_matrix
    frac -= np.round(frac)
    dx = frac @ pbc.basis_matrix
    dists = np.linalg.norm(dx, axis=1)
    n_overlap = int((dists < _SORBATE_OVERLAP).sum())
    if n_overlap > 0:
        warnings.append(
            f"Sorbate placement at ({cx:.2f}, {cy:.2f}, {cz:.2f}) overlaps "
            f"{n_overlap} framework atom{'s' if n_overlap != 1 else ''} "
            f"(within {_SORBATE_OVERLAP} Å). Use 'sorbate_pos=' to pick a "
            f"pore center, or let MPMC's μVT insertion handle placement."
        )

    return errors, warnings


def _pressure_dir_name(pressure: float) -> str:
    p_str = f"{pressure:.4g}".replace(".", "p").replace("-", "m")
    return f"P_{p_str}_atm"


def generate_isotherm(
    mol: Molecule,
    sorbate_model: str,
    temperature: float,
    pressures: list[float],
    output_dir: str,
    numsteps: int = 500000,
    corrtime: int = 1000,
    ensemble: str = "uvt",
    sorbate_pos: tuple[float, float, float] | None = None,
    insert_probability: float = 0.667,
    move_factor: float = 0.01,
    rot_factor: float = 0.01,
) -> list[str]:
    """Generate MPMC input directories for an isotherm.

    Creates one subdirectory per pressure point with:
      - input.pqr   : framework + sorbate, charges + LJ params
      - mpmc.inp    : rendered from templates/mpmc_uvt.inp
      - insert.pqr  : sorbate-only insert pool
      - run.sh      : launch script

    Args:
        mol: MOF structure with PBC.
        sorbate_model: Sorbate model name (from sorbates.py).
        temperature: Temperature in Kelvin.
        pressures: List of pressures in atm.
        output_dir: Base directory for all pressure subdirs.
        numsteps: MC steps per pressure point.
        corrtime: Correlation time for sampling.
        ensemble: "uvt" for grand canonical, "nvt" for canonical.
        sorbate_pos: Optional (x, y, z) for sorbate placement. Auto if None.
        insert_probability, move_factor, rot_factor: MPMC move probabilities.

    Returns:
        List of created directory paths.
    """
    if mol.pbc is None:
        raise ValueError("Molecule must have PBC for MPMC simulation")

    pbc = mol.pbc
    base = Path(output_dir)
    base.mkdir(parents=True, exist_ok=True)

    if sorbate_pos is None:
        center = mol.center()
        sx, sy, sz = float(center[0]), float(center[1]), float(center[2])
    else:
        sx, sy, sz = sorbate_pos

    sorbate_lines = format_sorbate_pqr(
        sorbate_model, sx, sy, sz, mol_id=2, start_atom_id=len(mol.atoms) + 1,
    )

    template = load_template("mpmc_uvt.inp")
    created_dirs: list[str] = []

    for pressure in pressures:
        p_name = _pressure_dir_name(pressure)
        p_dir = base / p_name
        p_dir.mkdir(exist_ok=True)

        pqr_path = str(p_dir / "input.pqr")
        write_mpmc_pdb(
            mol, pqr_path,
            write_charges=True, write_params=True, sorbate_lines=sorbate_lines,
        )

        rendered = render(template, {
            "job_name": f"isotherm_{p_name}",
            "ensemble": ensemble,
            "cell_a": f"{pbc.a:.4f}",
            "cell_b": f"{pbc.b:.4f}",
            "cell_c": f"{pbc.c:.4f}",
            "cell_alpha": f"{pbc.alpha:.2f}",
            "cell_beta": f"{pbc.beta:.2f}",
            "cell_gamma": f"{pbc.gamma:.2f}",
            "temperature": f"{temperature:.1f}",
            "pressure": pressure,
            "numsteps": numsteps,
            "corrtime": corrtime,
            "insert_probability": insert_probability,
            "move_factor": move_factor,
            "rot_factor": rot_factor,
            "pqr_input": "input.pqr",
        })
        (p_dir / "mpmc.inp").write_text(rendered)

        # insert.pqr — just the sorbate atoms + END
        (p_dir / "insert.pqr").write_text("\n".join(sorbate_lines) + "\nEND\n")

        # run.sh
        run_path = p_dir / "run.sh"
        run_path.write_text("#!/bin/bash\nmpmc mpmc.inp\n")
        run_path.chmod(0o755)

        created_dirs.append(str(p_dir))

    # Master run-all script
    run_all = base / "run_all.sh"
    lines = ["#!/bin/bash", f"# Isotherm: {sorbate_model} at {temperature} K", ""]
    for p_dir_str in created_dirs:
        p_name = Path(p_dir_str).name
        lines.append(f"echo 'Running {p_name}...'")
        lines.append(f"cd {p_name} && mpmc mpmc.inp > mpmc.log 2>&1 && cd ..")
    lines.append("")
    lines.append("echo 'All pressure points complete.'")
    run_all.write_text("\n".join(lines) + "\n")
    run_all.chmod(0o755)

    return created_dirs


def parse_isotherm_results(base_dir: str) -> list[tuple[float, float]]:
    """Parse MPMC output to extract isotherm data (pressure, uptake).

    Reads histogram.dat or energy.dat from each pressure subdirectory.

    Returns:
        List of (pressure_atm, avg_molecules) tuples, sorted by pressure.
    """
    base = Path(base_dir)
    results = []

    for p_dir in sorted(base.iterdir()):
        if not p_dir.is_dir() or not p_dir.name.startswith("P_"):
            continue
        p_str = p_dir.name.replace("P_", "").replace("_atm", "").replace("p", ".").replace("m", "-")
        try:
            pressure = float(p_str)
        except ValueError:
            continue

        hist_file = p_dir / "histogram.dat"
        if hist_file.exists():
            try:
                lines = hist_file.read_text().strip().split("\n")
                total_n = 0.0
                total_p = 0.0
                for line in lines:
                    if line.startswith("#"):
                        continue
                    tokens = line.split()
                    if len(tokens) >= 2:
                        n = float(tokens[0])
                        prob = float(tokens[1])
                        total_n += n * prob
                        total_p += prob
                if total_p > 0:
                    avg_n = total_n / total_p
                    results.append((pressure, avg_n))
                    continue
            except Exception:
                pass

        energy_file = p_dir / "energy.dat"
        if energy_file.exists():
            try:
                data = read_energy_dat(str(energy_file))
                if "N" in data and len(data["N"]) > 0:
                    n_vals = data["N"]
                    start = len(n_vals) * 3 // 4
                    avg_n = float(n_vals[start:].mean())
                    results.append((pressure, avg_n))
            except Exception:
                pass

    results.sort(key=lambda x: x[0])
    return results

# ======================================================================
# Module: reduce_cell
# ======================================================================
"""Reduce a supercell back to a single unit cell, with sorbate selection."""


import numpy as np



def identify_framework_and_sorbates(
    mol: Molecule,
    na: int, nb: int, nc: int,
    overlap_tol: float = 0.5,
) -> tuple[list[int], list[list[int]]]:
    """Identify framework atoms (periodic copies) vs sorbate atoms (unique).

    Given a supercell of dimensions na x nb x nc, wrap all atoms into the
    primitive cell. Framework atoms will cluster into groups of na*nb*nc
    overlapping copies. Sorbate atoms are unique — they don't overlap.

    Args:
        mol: Supercell molecule with PBC.
        na, nb, nc: Supercell dimensions along a, b, c.
        overlap_tol: Distance tolerance for considering atoms as copies (A).

    Returns:
        (framework_indices, sorbate_groups):
        framework_indices: indices of one representative atom per framework site
        sorbate_groups: list of groups, each group is a list of atom indices
            belonging to the same sorbate molecule
    """
    if mol.pbc is None or not mol.atoms:
        return [], []

    pbc = mol.pbc
    n_mult = na * nb * nc
    n_atoms = len(mol.atoms)

    # Compute fractional coordinates in the supercell
    coords = np.array([a.x for a in mol.atoms])
    frac = coords @ pbc.reciprocal_basis_matrix

    # Scale fractional coords to primitive cell [0, na) -> [0, 1)
    prim_frac = frac.copy()
    prim_frac[:, 0] *= na
    prim_frac[:, 1] *= nb
    prim_frac[:, 2] *= nc
    # Wrap to [0, 1)
    prim_frac -= np.floor(prim_frac)

    # Primitive cell basis
    prim_basis = pbc.basis_matrix.copy()
    prim_basis[0] /= na
    prim_basis[1] /= nb
    prim_basis[2] /= nc

    # Convert back to Cartesian in primitive cell
    prim_cart = prim_frac @ prim_basis

    # Cluster atoms by position: framework atoms cluster into groups of n_mult
    # Precompute inverse for min-image in primitive cell
    prim_inv = np.linalg.inv(prim_basis)

    np.zeros(n_atoms, dtype=bool)
    groups: list[list[int]] = []

    # Group by element first to reduce O(n^2) comparisons
    elem_indices: dict[str, list[int]] = {}
    for i, a in enumerate(mol.atoms):
        elem_indices.setdefault(a.element.symbol, []).append(i)

    for el, indices in elem_indices.items():
        el_used = set()
        for ii, i in enumerate(indices):
            if i in el_used:
                continue
            group = [i]
            el_used.add(i)
            for j in indices[ii + 1:]:
                if j in el_used:
                    continue
                diff = prim_cart[j] - prim_cart[i]
                frac_diff = diff @ prim_inv
                frac_diff -= np.round(frac_diff)
                cart_diff = frac_diff @ prim_basis
                dist = float(np.linalg.norm(cart_diff))
                if dist < overlap_tol:
                    group.append(j)
                    el_used.add(j)
            groups.append(group)

    # Framework groups have exactly n_mult members (periodic copies)
    # Sorbate groups have fewer
    framework_indices = []
    sorbate_groups = []

    for group in groups:
        if len(group) == n_mult:
            # Framework — keep one representative
            framework_indices.append(group[0])
        else:
            sorbate_groups.append(group)

    # Merge sorbate atoms into molecules using bond connectivity
    # Build neighbor map from existing bonds
    neighbors: dict[int, set[int]] = {}
    for a, b in mol.bonds:
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)

    # BFS to find connected sorbate molecules
    sorbate_atoms = set()
    for sg in sorbate_groups:
        sorbate_atoms.update(sg)

    visited = set()
    sorbate_molecules: list[list[int]] = []
    for seed_group in sorbate_groups:
        for seed in seed_group:
            if seed in visited:
                continue
            # BFS from seed
            mol_atoms = []
            queue = [seed]
            while queue:
                atom = queue.pop(0)
                if atom in visited:
                    continue
                visited.add(atom)
                mol_atoms.append(atom)
                for nb_idx in neighbors.get(atom, set()):
                    if nb_idx in sorbate_atoms and nb_idx not in visited:
                        queue.append(nb_idx)
            if mol_atoms:
                sorbate_molecules.append(sorted(mol_atoms))

    return framework_indices, sorbate_molecules


def reduce_to_primitive(
    mol: Molecule,
    na: int, nb: int, nc: int,
    keep_sorbate_indices: list[int] | None = None,
    overlap_tol: float = 0.5,
) -> Molecule:
    """Reduce a supercell to a single unit cell.

    Args:
        mol: Supercell molecule.
        na, nb, nc: Supercell dimensions.
        keep_sorbate_indices: Which sorbate molecule indices to keep
            (from identify_framework_and_sorbates). None = keep all.
        overlap_tol: Distance tolerance for overlap detection.

    Returns:
        New Molecule with primitive cell PBC and reduced atoms.
    """
    if mol.pbc is None:
        return mol

    framework_idx, sorbate_mols = identify_framework_and_sorbates(
        mol, na, nb, nc, overlap_tol,
    )

    pbc = mol.pbc
    prim_pbc = PBC(
        pbc.a / na, pbc.b / nb, pbc.c / nc,
        pbc.alpha, pbc.beta, pbc.gamma,
    )

    # Collect atoms: framework + selected sorbates
    keep_atoms = list(framework_idx)
    if keep_sorbate_indices is not None:
        for si in keep_sorbate_indices:
            if 0 <= si < len(sorbate_mols):
                keep_atoms.extend(sorbate_mols[si])
    else:
        for sm in sorbate_mols:
            keep_atoms.extend(sm)

    # Build new atom list with positions wrapped into primitive cell
    prim_basis = prim_pbc.basis_matrix
    new_atoms = []
    for idx in sorted(keep_atoms):
        a = mol.atoms[idx]
        # Convert to fractional in primitive cell
        frac = a.x @ prim_pbc.reciprocal_basis_matrix
        frac -= np.floor(frac)
        new_pos = frac @ prim_basis
        new_atom = Atom(float(new_pos[0]), float(new_pos[1]), float(new_pos[2]), a.name)
        new_atom.charge = a.charge
        new_atom.alpha = a.alpha
        new_atom.epsilon = a.epsilon
        new_atom.sigma = a.sigma
        new_atoms.append(new_atom)

    result = Molecule(atoms=new_atoms, pbc=prim_pbc)
    result.detect_bonds()
    return result

# ======================================================================
# Module: energy
# ======================================================================
"""MPMC energy.dat file reader and plotter."""



import numpy as np


def read_energy_dat(filepath: str) -> dict[str, np.ndarray]:
    """Read an MPMC energy.dat file.

    Parses the header line (starting with #) for column names.
    Returns dict mapping column name -> numpy array.
    """
    path = Path(filepath)
    if not path.exists():
        return {}

    lines = path.read_text().splitlines()
    if not lines:
        return {}

    # Parse header for column names
    header_names: list[str] = []
    data_lines: list[list[float]] = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#") or stripped.startswith("!"):
            # Parse header: "#step #energy #coulombic ..." or "# step energy ..."
            tokens = stripped.lstrip("#!").split()
            if tokens:
                header_names = [t.lstrip("#") for t in tokens]
            continue
        tokens = stripped.split()
        try:
            float(tokens[0])
            data_lines.append([float(t) for t in tokens])
        except (ValueError, IndexError):
            continue

    if not data_lines:
        return {}

    arr = np.array(data_lines)
    n_cols = arr.shape[1]

    # If we got header names, use them; otherwise use defaults
    if len(header_names) >= n_cols:
        col_names = header_names[:n_cols]
    elif header_names:
        col_names = header_names + [f"col{i}" for i in range(len(header_names), n_cols)]
    else:
        # Fallback: common MPMC column order
        defaults = ["step", "energy", "coulombic", "rd", "polar", "vdw",
                     "kinetic", "kin_temp", "N", "spin_ratio", "volume", "core_temp"]
        col_names = defaults[:min(n_cols, len(defaults))]
        col_names += [f"col{i}" for i in range(len(col_names), n_cols)]

    result: dict[str, np.ndarray] = {}
    for i, name in enumerate(col_names):
        result[name] = arr[:, i]

    return result

# ======================================================================
# Module: density_profile
# ======================================================================
"""1D density profile along a cell axis."""


import numpy as np



def compute_density_profile(
    mol: Molecule,
    axis: int = 2,
    n_bins: int = 100,
    element: str | None = None,
    name: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute 1D number density profile along a cell axis.

    Args:
        mol: Molecule with PBC.
        axis: 0=a, 1=b, 2=c.
        n_bins: Number of bins along the axis.
        element: Filter to this element (None = all).
        name: Filter to this atom name (None = all).

    Returns:
        (positions_along_axis, density) arrays.
    """
    if mol.pbc is None:
        return np.array([]), np.array([])

    pbc = mol.pbc
    atoms = mol.atoms
    if element:
        if name and name != "(all)":
            atoms = [a for a in atoms if a.element.symbol == element and a.name.strip() == name]
        else:
            atoms = [a for a in atoms if a.element.symbol == element]

    if not atoms:
        return np.array([]), np.array([])

    coords = np.array([a.x for a in atoms])
    # Project onto fractional coordinates along the chosen axis
    frac = coords @ pbc.reciprocal_basis_matrix
    frac_along = frac[:, axis] % 1.0

    # Histogram
    cell_length = [pbc.a, pbc.b, pbc.c][axis]
    bin_edges = np.linspace(0, 1, n_bins + 1)
    hist, _ = np.histogram(frac_along, bins=bin_edges)

    # Convert to number density (atoms per Angstrom)
    bin_width = cell_length / n_bins
    density = hist / bin_width

    positions = (bin_edges[:-1] + bin_edges[1:]) / 2 * cell_length

    return positions, density


def compute_density_profile_trajectory(
    frames: list[Molecule],
    axis: int = 2,
    n_bins: int = 100,
    element: str | None = None,
    name: str | None = None,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Average density profile over trajectory frames."""
    d_sum = None
    positions = None
    for fi, mol in enumerate(frames):
        p, d = compute_density_profile(mol, axis, n_bins, element, name)
        if len(p) == 0:
            continue
        if d_sum is None:
            d_sum = d.copy()
            positions = p
        else:
            d_sum += d
        if progress_callback:
            progress_callback((fi + 1) / len(frames))
    if d_sum is None:
        return np.array([]), np.array([])
    return positions, d_sum / len(frames)

# ======================================================================
# Module: pxrd
# ======================================================================
"""Simulated Powder X-Ray Diffraction (PXRD) pattern calculation.

Computes the diffraction pattern from crystal structure using the Debye
scattering equation with atomic form factors.
"""


import numpy as np


# Cromer-Mann coefficients for X-ray atomic scattering factors
# f(s) = sum_i(a_i * exp(-b_i * s^2)) + c, where s = sin(theta)/lambda
# From International Tables for Crystallography, Vol C
# Format: {element: (a1,b1, a2,b2, a3,b3, a4,b4, c)}
SCATTERING_FACTORS: dict[str, tuple[float, ...]] = {
    "H":  (0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.5519, 0.049879, 2.20159, 0.001305),
    "He": (0.8734, 9.1037, 0.6309, 3.3568, 0.3112, 22.9276, 0.1780, 0.9821, 0.0064),
    "Li": (1.1282, 3.9546, 0.7508, 1.0524, 0.6175, 85.3905, 0.4653, 168.261, 0.0377),
    "Be": (1.5919, 43.6427, 1.1278, 1.8623, 0.5391, 103.483, 0.7029, 0.5420, 0.0385),
    "B":  (2.0545, 23.2185, 1.3326, 1.0210, 1.0979, 60.3498, 0.7068, 0.1403, -0.1932),
    "C":  (2.3100, 20.8439, 1.0200, 10.2075, 1.5886, 0.5687, 0.8650, 51.6512, 0.2156),
    "N":  (12.2126, 0.0057, 3.1322, 9.8933, 2.0125, 28.9975, 1.1663, 0.5826, -11.529),
    "O":  (3.0485, 13.2771, 2.2868, 5.7011, 1.5463, 0.3239, 0.8670, 32.9089, 0.2508),
    "F":  (3.5392, 10.2825, 2.6412, 4.2944, 1.5170, 0.2615, 1.0243, 26.1476, 0.2776),
    "Na": (4.7626, 3.2850, 3.1736, 8.8422, 1.2674, 0.3136, 1.1128, 129.424, 0.676),
    "Mg": (5.4204, 2.8275, 2.1735, 79.2611, 1.2269, 0.3808, 2.3073, 7.1937, 0.8584),
    "Al": (6.4202, 3.0387, 1.9002, 0.7426, 1.5936, 31.5472, 1.9646, 85.0886, 1.1151),
    "Si": (6.2915, 2.4386, 3.0353, 32.3337, 1.9891, 0.6785, 1.5410, 81.6937, 1.1407),
    "P":  (6.4345, 1.9067, 4.1791, 27.1570, 1.7800, 0.5260, 1.4908, 68.1645, 1.1149),
    "S":  (6.9053, 1.4679, 5.2034, 22.2151, 1.4379, 0.2536, 1.5863, 56.1720, 0.8669),
    "Cl": (11.4604, 0.0104, 7.1964, 1.1662, 6.2556, 18.5194, 1.6455, 47.7784, -9.5574),
    "K":  (8.2186, 12.7949, 7.4398, 0.7748, 1.0519, 213.187, 0.8659, 41.6841, 1.4228),
    "Ca": (8.6266, 10.4421, 7.3873, 0.6599, 1.5899, 85.7484, 1.0211, 178.437, 1.3751),
    "Ti": (9.7595, 7.8508, 7.3558, 0.5000, 1.6991, 35.6338, 1.9021, 116.105, 1.2807),
    "V":  (10.2971, 6.8657, 7.3511, 0.4385, 2.0703, 26.8938, 2.0571, 102.478, 1.2199),
    "Cr": (10.6406, 6.1038, 7.3537, 0.3920, 3.3240, 20.2626, 1.4922, 98.7399, 1.1832),
    "Mn": (11.2819, 5.3409, 7.3573, 0.3432, 3.0193, 17.8674, 2.2441, 83.7543, 1.0896),
    "Fe": (11.7695, 4.7611, 7.3573, 0.3072, 3.5222, 15.3535, 2.3045, 76.8805, 1.0369),
    "Co": (12.2841, 4.2791, 7.3409, 0.2784, 4.0034, 13.5359, 2.3488, 71.1692, 1.0118),
    "Ni": (12.8376, 3.8785, 7.2920, 0.2565, 4.4438, 12.1763, 2.3800, 66.3421, 1.0341),
    "Cu": (13.3380, 3.5828, 7.1676, 0.2470, 5.6158, 11.3966, 1.6735, 64.8126, 1.1910),
    "Zn": (14.0743, 3.2655, 7.0318, 0.2333, 5.1652, 10.3163, 2.4100, 58.7097, 1.3041),
    "Br": (17.1789, 2.1723, 5.2358, 16.5796, 5.6377, 0.2609, 3.9851, 41.4328, 2.9557),
    "Ag": (19.2808, 0.6446, 16.6885, 7.4726, 4.8045, 24.6605, 1.0463, 99.8156, 5.179),
    "Cd": (19.2214, 0.5946, 17.6444, 6.9089, 4.4610, 24.7008, 1.6029, 87.4825, 5.0694),
    "I":  (20.1472, 4.3470, 18.9949, 0.3814, 7.5138, 27.7660, 2.2735, 66.8776, 4.0712),
    "Ba": (20.3361, 3.2160, 19.2970, 0.2756, 10.888, 20.2073, 2.6959, 167.202, 2.7731),
    "La": (20.578, 2.94817, 19.599, 0.244475, 11.3727, 18.7726, 3.28719, 133.124, 2.14678),
    "Pt": (27.0059, 1.5129, 17.7639, 8.8117, 15.7131, 0.1243, 5.7837, 38.6103, 11.6883),
    "Au": (16.8819, 0.4611, 18.5913, 8.6216, 25.5582, 1.4826, 5.860, 36.3956, 12.0658),
    "Pb": (31.0617, 0.6902, 13.0637, 2.3576, 18.442, 8.618, 5.9696, 47.2579, 13.4118),
}

# Default for elements not in the table
_DEFAULT_SF = (6.0, 3.0, 3.0, 10.0, 2.0, 0.5, 1.0, 50.0, 0.5)


def _atomic_form_factor(element: str, s: np.ndarray) -> np.ndarray:
    """Compute atomic form factor f(s) for an element.
    s = sin(theta)/lambda in A^-1.
    """
    params = SCATTERING_FACTORS.get(element, _DEFAULT_SF)
    a1, b1, a2, b2, a3, b3, a4, b4, c = params
    s2 = s * s
    return (a1 * np.exp(-b1 * s2) + a2 * np.exp(-b2 * s2) +
            a3 * np.exp(-b3 * s2) + a4 * np.exp(-b4 * s2) + c)


def compute_pxrd(
    mol: Molecule,
    wavelength: float = 1.5406,  # Cu K-alpha in Angstroms
    two_theta_min: float = 5.0,
    two_theta_max: float = 50.0,
    n_points: int = 2000,
    peak_width: float = 0.1,  # Gaussian broadening in degrees
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute simulated PXRD pattern using Bragg's law and structure factors.

    Uses the crystal structure with PBC to compute hkl reflections and
    their intensities from structure factors.

    Args:
        mol: Molecule with PBC (crystal structure).
        wavelength: X-ray wavelength in Angstroms (default Cu K-alpha).
        two_theta_min: Minimum 2-theta in degrees.
        two_theta_max: Maximum 2-theta in degrees.
        n_points: Number of points in the output pattern.
        peak_width: Gaussian peak broadening (FWHM in degrees).
        progress_callback: Optional callable(frac).

    Returns:
        (two_theta, intensity) arrays.
    """
    if mol.pbc is None or len(mol.atoms) == 0:
        return np.array([]), np.array([])

    pbc = mol.pbc
    recip = pbc.reciprocal_basis_matrix  # rows are reciprocal vectors

    # Determine hkl range from max 2-theta
    # d_min = lambda / (2 * sin(theta_max))
    theta_max = np.radians(two_theta_max / 2)
    d_min = wavelength / (2 * np.sin(theta_max))
    1.0 / (2 * d_min)  # sin(theta)/lambda max

    # Generate hkl indices
    # |G_hkl| = |h*a* + k*b* + l*c*|, d = 1/|G|
    h_max = int(np.ceil(pbc.a / d_min)) + 1
    k_max = int(np.ceil(pbc.b / d_min)) + 1
    l_max = int(np.ceil(pbc.c / d_min)) + 1
    # Cap to reasonable range
    h_max = min(h_max, 20)
    k_max = min(k_max, 20)
    l_max = min(l_max, 20)

    coords = np.array([a.x for a in mol.atoms])
    elements = np.array([a.element.symbol for a in mol.atoms])
    unique_elements = list(set(elements.tolist()))

    # Pre-compute fractional coordinates and element-grouped slices
    frac_coords = coords @ recip
    el_groups = {el: frac_coords[elements == el] for el in unique_elements}

    # Generate the full (h,k,l) grid, drop (0,0,0), filter to two_theta range.
    h_arr = np.arange(-h_max, h_max + 1)
    k_arr = np.arange(-k_max, k_max + 1)
    l_arr = np.arange(-l_max, l_max + 1)
    H, K, L = np.meshgrid(h_arr, k_arr, l_arr, indexing="ij")
    hkl = np.stack([H.ravel(), K.ravel(), L.ravel()], axis=1)        # (M, 3)
    hkl = hkl[(hkl != 0).any(axis=1)]
    G_all = hkl.astype(float) @ recip                                # (M, 3)
    G_norm = np.linalg.norm(G_all, axis=1)
    valid = G_norm > 1e-10
    hkl = hkl[valid]
    G_norm = G_norm[valid]
    d_spacing = 1.0 / G_norm
    sin_theta = wavelength / (2 * d_spacing)
    valid = np.abs(sin_theta) <= 1.0
    hkl = hkl[valid]
    sin_theta = sin_theta[valid]
    two_theta_all = 2 * np.degrees(np.arcsin(sin_theta))
    in_range = (two_theta_all >= two_theta_min) & (two_theta_all <= two_theta_max)
    hkl = hkl[in_range]
    two_theta_all = two_theta_all[in_range]
    sin_theta = sin_theta[in_range]
    n_refl = len(hkl)

    # Form factors per element vs s = sin(theta)/lambda — vectorized across reflections
    s_all = sin_theta / wavelength
    f_per_el = {el: _atomic_form_factor(el, s_all) for el in unique_elements}

    # Structure factor: for each element group, sum cos/sin of 2π·(hkl·frac)
    # over its atoms. Then total F = Σ_el f_el(s) * Σ_atoms_in_el e^{2πi(hkl·frac)}
    F_real = np.zeros(n_refl)
    F_imag = np.zeros(n_refl)
    if n_refl > 0:
        for el in unique_elements:
            frac_g = el_groups[el]              # (n_el, 3)
            if len(frac_g) == 0:
                continue
            # phase shape: (n_refl, n_el)
            phase = 2 * np.pi * (hkl.astype(float) @ frac_g.T)
            cos_sum = np.cos(phase).sum(axis=1)
            sin_sum = np.sin(phase).sum(axis=1)
            f_el = f_per_el[el]                # (n_refl,)
            F_real += f_el * cos_sum
            F_imag += f_el * sin_sum

    intensities = F_real ** 2 + F_imag ** 2

    # Lorentz-polarization factor
    theta_rad = np.radians(two_theta_all / 2)
    cos2t = np.cos(np.radians(two_theta_all))
    sin_t = np.sin(theta_rad)
    cos_t = np.cos(theta_rad)
    denom = sin_t * sin_t * cos_t
    lp = np.where(np.abs(denom) > 1e-10, (1 + cos2t ** 2) / denom, 0.0)
    intensities *= np.abs(lp)

    reflections = list(zip(two_theta_all.tolist(), intensities.tolist()))

    if progress_callback:
        progress_callback(1.0)

    if not reflections:
        return np.array([]), np.array([])

    # Build pattern by convolving with Gaussian peaks
    two_theta_arr = np.linspace(two_theta_min, two_theta_max, n_points)
    pattern = np.zeros(n_points)

    sigma = peak_width / (2 * np.sqrt(2 * np.log(2)))  # FWHM to sigma

    for tt, intensity in reflections:
        peak = intensity * np.exp(-0.5 * ((two_theta_arr - tt) / sigma) ** 2)
        pattern += peak

    # Normalize to max = 100
    if pattern.max() > 0:
        pattern = pattern / pattern.max() * 100

    return two_theta_arr, pattern

# ======================================================================
# Module: formats
# ======================================================================
"""Additional file format readers/writers: POSCAR, LAMMPS data, MOL2, CIF write."""



import numpy as np


# ---------------------------------------------------------------------------
# POSCAR / CONTCAR (VASP)
# ---------------------------------------------------------------------------

def read_poscar(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:
    """Read VASP POSCAR/CONTCAR format."""
    comment = file.readline().strip()
    scale = float(file.readline().strip())

    # Lattice vectors
    vecs = []
    for _ in range(3):
        tokens = file.readline().split()
        vecs.append([float(t) * scale for t in tokens[:3]])
    basis = np.array(vecs)

    a = float(np.linalg.norm(basis[0]))
    b = float(np.linalg.norm(basis[1]))
    c = float(np.linalg.norm(basis[2]))

    # Compute angles from dot products
    cos_alpha = np.dot(basis[1], basis[2]) / (b * c)
    cos_beta = np.dot(basis[0], basis[2]) / (a * c)
    cos_gamma = np.dot(basis[0], basis[1]) / (a * b)
    alpha = float(np.degrees(np.arccos(np.clip(cos_alpha, -1, 1))))
    beta = float(np.degrees(np.arccos(np.clip(cos_beta, -1, 1))))
    gamma = float(np.degrees(np.arccos(np.clip(cos_gamma, -1, 1))))

    pbc = PBC(a, b, c, alpha, beta, gamma)

    # Species line (VASP 5+)
    line = file.readline().strip()
    tokens = line.split()
    try:
        int(tokens[0])
        # Old VASP 4 format — no species line, this IS the counts
        species = [comment.split()[i] if i < len(comment.split()) else f"X{i}" for i in range(len(tokens))]
        counts = [int(t) for t in tokens]
    except ValueError:
        species = tokens
        counts = [int(t) for t in file.readline().split()]

    # Selective dynamics or coordinate type
    line = file.readline().strip()
    if line.lower().startswith("s"):
        line = file.readline().strip()

    direct = line.lower().startswith("d")

    system: list[Atom] = []
    for sp, count in zip(species, counts):
        for _ in range(count):
            tokens = file.readline().split()
            coords = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
            if direct:
                coords = coords @ basis
            atom = Atom(coords[0], coords[1], coords[2], sp)
            system.append(atom)

    set_atom_ids(system)
    return system, pbc


def write_poscar(mol: Molecule, out: TextIO) -> None:
    """Write VASP POSCAR format."""
    if mol.pbc is None:
        return

    from collections import OrderedDict
    species_counts: OrderedDict[str, int] = OrderedDict()
    for a in mol.atoms:
        sym = a.element.symbol
        species_counts[sym] = species_counts.get(sym, 0) + 1

    out.write("Generated by PDB Wizard\n")
    out.write("1.0\n")
    for row in mol.pbc.basis_matrix:
        out.write(f"  {row[0]:>16.10f} {row[1]:>16.10f} {row[2]:>16.10f}\n")
    out.write(" ".join(species_counts.keys()) + "\n")
    out.write(" ".join(str(c) for c in species_counts.values()) + "\n")
    out.write("Cartesian\n")
    for a in mol.atoms:
        out.write(f"  {a.x[0]:>16.10f} {a.x[1]:>16.10f} {a.x[2]:>16.10f}\n")


# ---------------------------------------------------------------------------
# LAMMPS data file
# ---------------------------------------------------------------------------

def read_lammps_data(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:
    """Read LAMMPS data file (atomic style)."""
    lines = file.readlines()
    xlo, xhi, ylo, yhi, zlo, zhi = 0.0, 10.0, 0.0, 10.0, 0.0, 10.0
    atom_lines: list[str] = []
    masses: dict[int, float] = {}
    in_atoms = False
    in_masses = False

    _section_keywords = {"Masses", "Atoms", "Bonds", "Angles", "Dihedrals",
                         "Impropers", "Velocities", "Pair Coeffs", "Bond Coeffs"}
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        # Check for section headers
        if any(stripped.startswith(kw) for kw in _section_keywords):
            in_atoms = stripped.startswith("Atoms")
            in_masses = stripped.startswith("Masses")
            continue
        if in_atoms:
            atom_lines.append(stripped)
        elif in_masses:
            tokens = stripped.split()
            if len(tokens) >= 2:
                try:
                    masses[int(tokens[0])] = float(tokens[1])
                except ValueError:
                    pass
        elif "atoms" in stripped:
            try:
                int(stripped.split()[0])
            except ValueError:
                pass
        elif "xlo" in stripped and "xhi" in stripped:
            tokens = stripped.split()
            xlo, xhi = float(tokens[0]), float(tokens[1])
        elif "ylo" in stripped and "yhi" in stripped:
            tokens = stripped.split()
            ylo, yhi = float(tokens[0]), float(tokens[1])
        elif "zlo" in stripped and "zhi" in stripped:
            tokens = stripped.split()
            zlo, zhi = float(tokens[0]), float(tokens[1])

    pbc = PBC(xhi - xlo, yhi - ylo, zhi - zlo, 90.0, 90.0, 90.0)

    # Map masses to elements (approximate)
    mass_to_el: dict[int, str] = {}
    for type_id, mass in masses.items():
        best_el = "X"
        best_diff = 999.0
        for sym, el in ELEMENTS.items():
            if el.mass > 0 and abs(el.mass - mass) < best_diff:
                best_diff = abs(el.mass - mass)
                best_el = sym
        mass_to_el[type_id] = best_el

    # Detect atom style from first data line
    # atomic: id type x y z          (5 tokens)
    # charge: id type charge x y z   (6 tokens)
    # full:   id mol-id type charge x y z  (7 tokens)
    style = "atomic"
    if atom_lines:
        n_tok = len(atom_lines[0].split())
        if n_tok >= 7:
            style = "full"
        elif n_tok >= 6:
            style = "charge"

    system: list[Atom] = []
    for line in atom_lines:
        tokens = line.split()
        try:
            if style == "full":
                # id mol-id type charge x y z
                atom_type = int(tokens[2])
                charge = float(tokens[3])
                x, y, z = float(tokens[4]), float(tokens[5]), float(tokens[6])
            elif style == "charge":
                # id type charge x y z
                atom_type = int(tokens[1])
                charge = float(tokens[2])
                x, y, z = float(tokens[3]), float(tokens[4]), float(tokens[5])
            else:
                # id type x y z
                atom_type = int(tokens[1])
                charge = 0.0
                x, y, z = float(tokens[2]), float(tokens[3]), float(tokens[4])
            el = mass_to_el.get(atom_type, f"T{atom_type}")
            atom = Atom(x, y, z, el)
            atom.charge = charge
            system.append(atom)
        except (ValueError, IndexError):
            pass

    set_atom_ids(system)
    return system, pbc


# ---------------------------------------------------------------------------
# CIF writer
# ---------------------------------------------------------------------------

def write_cif(mol: Molecule, out: TextIO) -> None:
    """Write a CIF file from molecule with PBC."""
    if mol.pbc is None:
        return

    pbc = mol.pbc
    out.write("data_pdb_wizard\n")
    out.write(f"_cell_length_a    {pbc.a:.6f}\n")
    out.write(f"_cell_length_b    {pbc.b:.6f}\n")
    out.write(f"_cell_length_c    {pbc.c:.6f}\n")
    out.write(f"_cell_angle_alpha {pbc.alpha:.4f}\n")
    out.write(f"_cell_angle_beta  {pbc.beta:.4f}\n")
    out.write(f"_cell_angle_gamma {pbc.gamma:.4f}\n")
    out.write("_symmetry_space_group_name_H-M   'P 1'\n")
    out.write("_symmetry_Int_Tables_number       1\n\n")
    out.write("loop_\n")
    out.write("_symmetry_equiv_pos_as_xyz\n")
    out.write("'x, y, z'\n\n")
    out.write("loop_\n")
    out.write("_atom_site_label\n")
    out.write("_atom_site_type_symbol\n")
    out.write("_atom_site_fract_x\n")
    out.write("_atom_site_fract_y\n")
    out.write("_atom_site_fract_z\n")

    # Convert to fractional
    coords = np.array([a.x for a in mol.atoms])
    frac = coords @ pbc.reciprocal_basis_matrix

    label_counts: dict[str, int] = {}
    for i, atom in enumerate(mol.atoms):
        sym = atom.element.symbol
        label_counts[sym] = label_counts.get(sym, 0) + 1
        label = f"{sym}{label_counts[sym]}"
        out.write(f"{label:<6} {sym:<3} {frac[i, 0]:.6f} {frac[i, 1]:.6f} {frac[i, 2]:.6f}\n")


# ---------------------------------------------------------------------------
# Gaussian log/out reader
# ---------------------------------------------------------------------------

def read_gaussian_log(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:
    """Read atoms from the LAST Standard/Input orientation block in a Gaussian log."""

    lines = file.readlines()

    # Find the last orientation block
    last_block_start: int | None = None
    for i, line in enumerate(lines):
        if "Standard orientation" in line or "Input orientation" in line:
            last_block_start = i

    if last_block_start is None:
        return [], None

    # Skip header lines: orientation title, dashes, column headers, dashes
    # Then atom rows until next dashes line
    idx = last_block_start + 1
    # Skip to the line after the second dashed separator
    dash_count = 0
    while idx < len(lines):
        if lines[idx].strip().startswith("---"):
            dash_count += 1
            if dash_count == 2:
                idx += 1
                break
        idx += 1

    system: list[Atom] = []
    while idx < len(lines):
        line = lines[idx].strip()
        if line.startswith("---"):
            break
        tokens = line.split()
        if len(tokens) >= 6:
            atomic_number = int(tokens[1])
            x, y, z = float(tokens[3]), float(tokens[4]), float(tokens[5])
            el = get_element_by_number(atomic_number)
            system.append(Atom(x, y, z, el.symbol))
        idx += 1

    set_atom_ids(system)
    return system, None


# ---------------------------------------------------------------------------
# Gaussian com/gjf reader
# ---------------------------------------------------------------------------

def read_gaussian_com(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:
    """Read atoms from a Gaussian .com/.gjf input file."""
    lines = file.readlines()
    idx = 0

    # Skip Link0 commands (lines starting with %)
    while idx < len(lines) and lines[idx].strip().startswith("%"):
        idx += 1

    # Skip route section (lines starting with #, may span multiple lines)
    if idx < len(lines) and lines[idx].strip().startswith("#"):
        while idx < len(lines) and lines[idx].strip():
            idx += 1

    # Skip blank line after route
    while idx < len(lines) and not lines[idx].strip():
        idx += 1

    # Skip title section (non-blank lines)
    while idx < len(lines) and lines[idx].strip():
        idx += 1

    # Skip blank line after title
    while idx < len(lines) and not lines[idx].strip():
        idx += 1

    # Skip charge and multiplicity line
    if idx < len(lines):
        idx += 1

    # Read atom lines until blank line or EOF
    system: list[Atom] = []
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            break
        tokens = line.split()
        if len(tokens) >= 4:
            element = tokens[0]
            x, y, z = float(tokens[1]), float(tokens[2]), float(tokens[3])
            system.append(Atom(x, y, z, element))
        idx += 1

    set_atom_ids(system)
    return system, None


# ---------------------------------------------------------------------------
# SDF / MOL reader
# ---------------------------------------------------------------------------

def read_sdf(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:
    """Read atoms from an SDF/MOL file."""
    lines = file.readlines()

    # Lines 0-2: header (molecule name, program/timestamp, comment)
    # Line 3: counts line
    if len(lines) < 4:
        return [], None

    counts_line = lines[3]
    n_atoms = int(counts_line[:3].strip())
    # n_bonds = int(counts_line[3:6].strip())  # not used — detect_bonds handles it

    system: list[Atom] = []
    for i in range(4, 4 + n_atoms):
        if i >= len(lines):
            break
        tokens = lines[i].split()
        if len(tokens) >= 4:
            x, y, z = float(tokens[0]), float(tokens[1]), float(tokens[2])
            element = tokens[3]
            system.append(Atom(x, y, z, element))

    set_atom_ids(system)
    return system, None


# ---------------------------------------------------------------------------
# Tripos MOL2 reader
# ---------------------------------------------------------------------------

def read_mol2(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:
    """Read atoms from a Tripos MOL2 file."""
    lines = file.readlines()

    # Find @<TRIPOS>ATOM section
    atom_start: int | None = None
    for i, line in enumerate(lines):
        if line.strip().upper() == "@<TRIPOS>ATOM":
            atom_start = i + 1
            break

    if atom_start is None:
        return [], None

    system: list[Atom] = []
    idx = atom_start
    while idx < len(lines):
        line = lines[idx].strip()
        if line.startswith("@") or not line:
            if line.startswith("@"):
                break
            idx += 1
            continue
        tokens = line.split()
        # Columns: atom_id atom_name x y z atom_type [residue_id residue_name charge]
        if len(tokens) >= 6:
            x, y, z = float(tokens[2]), float(tokens[3]), float(tokens[4])
            atom_type = tokens[5]
            # Element is the part before the dot in atom_type (e.g. "C.3" -> "C")
            element = atom_type.split(".")[0]
            system.append(Atom(x, y, z, element))
        idx += 1

    set_atom_ids(system)
    return system, None


# ---------------------------------------------------------------------------
# LAMMPS data file writer (full style)
# ---------------------------------------------------------------------------

def write_lammps_data(mol: Molecule, out: TextIO) -> None:
    """Write a LAMMPS data file in 'full' atom style."""
    from collections import OrderedDict


    # Build atom type mapping: element symbol -> type id
    type_map: OrderedDict[str, int] = OrderedDict()
    for atom in mol.atoms:
        sym = atom.element.symbol
        if sym not in type_map:
            type_map[sym] = len(type_map) + 1

    n_atoms = len(mol.atoms)
    n_bonds = len(mol.bonds)
    n_atom_types = len(type_map)

    # Assign bond types by element pair
    bond_type_map: dict[tuple[str, str], int] = {}
    for a1, a2 in mol.bonds:
        pair = tuple(sorted([mol.atoms[a1].element.symbol, mol.atoms[a2].element.symbol]))
        if pair not in bond_type_map:
            bond_type_map[pair] = len(bond_type_map) + 1
    n_bond_types = len(bond_type_map) if n_bonds > 0 else 0

    # Assign molecule IDs via connected components
    mol_ids = [0] * n_atoms
    neighbors: dict[int, list[int]] = {}
    for a, b in mol.bonds:
        neighbors.setdefault(a, []).append(b)
        neighbors.setdefault(b, []).append(a)
    current_mol = 0
    visited = [False] * n_atoms
    for seed in range(n_atoms):
        if visited[seed]:
            continue
        current_mol += 1
        stack = [seed]
        while stack:
            node = stack.pop()
            if visited[node]:
                continue
            visited[node] = True
            mol_ids[node] = current_mol
            for nb in neighbors.get(node, []):
                if not visited[nb]:
                    stack.append(nb)

    # Box dimensions
    if mol.pbc is not None:
        xlo, ylo, zlo = 0.0, 0.0, 0.0
        xhi = mol.pbc.a
        yhi = mol.pbc.b
        zhi = mol.pbc.c
    else:
        coords = np.array([a.x for a in mol.atoms])
        margin = 5.0
        xlo = float(coords[:, 0].min()) - margin
        xhi = float(coords[:, 0].max()) + margin
        ylo = float(coords[:, 1].min()) - margin
        yhi = float(coords[:, 1].max()) + margin
        zlo = float(coords[:, 2].min()) - margin
        zhi = float(coords[:, 2].max()) + margin

    # Header
    out.write("LAMMPS data file via pdb_wizard\n\n")
    out.write(f"{n_atoms} atoms\n")
    if n_bonds > 0:
        out.write(f"{n_bonds} bonds\n")
    out.write(f"{n_atom_types} atom types\n")
    if n_bonds > 0:
        out.write(f"{n_bond_types} bond types\n")
    out.write("\n")
    out.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
    out.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
    out.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")

    # Masses section
    out.write("\nMasses\n\n")
    for sym, type_id in type_map.items():
        mass = ELEMENTS[sym].mass if sym in ELEMENTS else 1.0
        out.write(f"{type_id} {mass:.4f} # {sym}\n")

    # Atoms section (full style: id mol-id type charge x y z)
    out.write("\nAtoms\n\n")
    for i, atom in enumerate(mol.atoms):
        atom_id = i + 1
        type_id = type_map[atom.element.symbol]
        charge = atom.charge
        out.write(f"{atom_id} {mol_ids[i]} {type_id} {charge:.6f} "
                  f"{atom.x[0]:.6f} {atom.x[1]:.6f} {atom.x[2]:.6f}\n")

    # Bonds section
    if n_bonds > 0:
        out.write("\nBonds\n\n")
        for i, (a1, a2) in enumerate(mol.bonds):
            bond_id = i + 1
            pair = tuple(sorted([mol.atoms[a1].element.symbol, mol.atoms[a2].element.symbol]))
            bond_type = bond_type_map[pair]
            out.write(f"{bond_id} {bond_type} {a1 + 1} {a2 + 1}\n")

# ======================================================================
# Module: database
# ======================================================================
"""Structure database fetching: COD, RCSB PDB, PubChem, Materials Project."""



# Some servers (notably PubChem) reject Python's default User-Agent.
_USER_AGENT = "pdb-wizard/1.0 (+https://github.com)"
_DEFAULT_TIMEOUT = 30.0


def _http_get(url: str, *, timeout: float = _DEFAULT_TIMEOUT,
              extra_headers: dict | None = None) -> bytes:
    """GET a URL with timeout and User-Agent. Returns the response bytes.
    Raises urllib.error.URLError / HTTPError on failure."""
    headers = {"User-Agent": _USER_AGENT}
    if extra_headers:
        headers.update(extra_headers)
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _http_download(url: str, out_path: str, *,
                   timeout: float = _DEFAULT_TIMEOUT,
                   extra_headers: dict | None = None) -> str:
    """Download a URL to disk with timeout and User-Agent. Returns out_path."""
    payload = _http_get(url, timeout=timeout, extra_headers=extra_headers)
    with open(out_path, "wb") as f:
        f.write(payload)
    return out_path


def fetch_cod(cod_id: str, output_dir: str = ".") -> str:
    """Fetch a CIF file from the Crystallography Open Database by ID.

    Args:
        cod_id: COD entry ID (e.g., '4512072')
        output_dir: Directory to save the file.

    Returns:
        Path to the downloaded CIF file.
    """
    url = f"https://www.crystallography.net/cod/{cod_id}.cif"
    out_path = str(Path(output_dir) / f"COD_{cod_id}.cif")
    return _http_download(url, out_path)


def fetch_rcsb(pdb_id: str, output_dir: str = ".") -> str:
    """Fetch a PDB file from the RCSB Protein Data Bank.

    Args:
        pdb_id: 4-letter PDB code (e.g., '1VII')
        output_dir: Directory to save the file.

    Returns:
        Path to the downloaded PDB file.
    """
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out_path = str(Path(output_dir) / f"{pdb_id.upper()}.pdb")
    return _http_download(url, out_path)


def search_cod(formula: str = "", text: str = "") -> list[dict]:
    """Search COD for structures matching a formula or text query.

    Returns list of dicts with 'id', 'formula', 'title', 'source' keys.
    """
    params = []
    if formula:
        params.append("formula=" + urllib.parse.quote(formula))
    if text:
        params.append("text=" + urllib.parse.quote(text))
    if not params:
        return []

    query = "&".join(params)
    url = "https://www.crystallography.net/cod/result.php?" + query + "&format=json"
    try:
        data = json.loads(_http_get(url, timeout=10))
        results = []
        for entry in data[:50]:
            results.append({
                "source": "cod",
                "id": str(entry.get("file", "")),
                "formula": entry.get("formula", ""),
                "title": entry.get("title", "")[:60],
            })
        return results
    except Exception:
        return []


def search_rcsb(query: str) -> list[dict]:
    """Search RCSB PDB for structures matching a text query.

    Returns list of dicts with 'id', 'formula', 'title', 'source' keys.
    """
    if not query.strip():
        return []
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query},
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 25}},
    }
    try:
        req = urllib.request.Request(
            url,
            data=json.dumps(payload).encode(),
            headers={
                "Content-Type": "application/json",
                "User-Agent": _USER_AGENT,
            },
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
        results = []
        for hit in data.get("result_set", []):
            pdb_id = hit.get("identifier", "")
            results.append({
                "source": "rcsb",
                "id": pdb_id,
                "formula": "",
                "title": pdb_id,
            })
        return results
    except Exception:
        return []


def fetch_pubchem_sdf(name: str, output_dir: str = ".") -> str:
    """Fetch a 3D SDF file from PubChem by compound name.

    Args:
        name: Compound name (e.g., 'aspirin', 'benzene') or numeric CID.
        output_dir: Directory to save the file.

    Returns:
        Path to the downloaded SDF file.
    """
    # PubChem search returns CIDs; allow either name or CID here.
    is_cid = name.isdigit()
    encoded = urllib.parse.quote(name)
    if is_cid:
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            + encoded + "/SDF?record_type=3d"
        )
        out_path = str(Path(output_dir) / f"CID_{name}.sdf")
    else:
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            + encoded + "/SDF?record_type=3d"
        )
        out_path = str(Path(output_dir) / (name.replace(" ", "_") + ".sdf"))
    return _http_download(url, out_path)


def search_pubchem(query: str) -> list[dict]:
    """Search PubChem for compounds by name.

    Returns list of dicts with 'id' (CID), 'formula', 'title', 'source' keys.
    """
    if not query.strip():
        return []
    encoded = urllib.parse.quote(query)
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + encoded + "/property/MolecularFormula,IUPACName/JSON")
    try:
        data = json.loads(_http_get(url, timeout=10))
        results = []
        for prop in data.get("PropertyTable", {}).get("Properties", [])[:25]:
            cid = str(prop.get("CID", ""))
            results.append({
                "source": "pubchem",
                "id": cid,
                "formula": prop.get("MolecularFormula", ""),
                "title": prop.get("IUPACName", "")[:60],
            })
        return results
    except Exception:
        return []


def search_materials_project(formula: str, api_key: str = "") -> list[dict]:
    """Search Materials Project for structures by formula.

    Requires MP_API_KEY environment variable or api_key parameter.
    Returns list of dicts with 'id', 'formula', 'title', 'source' keys.
    """
    import os
    key = api_key or os.environ.get("MP_API_KEY", "")
    if not key or not formula.strip():
        return []
    encoded = urllib.parse.quote(formula.strip())
    url = ("https://api.materialsproject.org/materials/summary/"
           "?formula=" + encoded + "&_limit=25&_fields=material_id,formula_pretty,symmetry")
    try:
        data = json.loads(_http_get(url, timeout=15, extra_headers={"X-API-KEY": key}))
        results = []
        for entry in data.get("data", []):
            mid = entry.get("material_id", "")
            formula_p = entry.get("formula_pretty", "")
            sym = entry.get("symmetry", {})
            sg = sym.get("symbol", "") if isinstance(sym, dict) else ""
            results.append({
                "source": "mp",
                "id": mid,
                "formula": formula_p,
                "title": sg,
            })
        return results
    except Exception:
        return []


def fetch_materials_project_cif(
    material_id: str, api_key: str = "", output_dir: str = ".",
) -> str:
    """Fetch a CIF from Materials Project by material ID (e.g. 'mp-149')."""
    import os
    key = api_key or os.environ.get("MP_API_KEY", "")
    if not key:
        raise ValueError("Set MP_API_KEY env var (free at materialsproject.org)")
    url = ("https://api.materialsproject.org/materials/summary/"
           "?material_ids=" + material_id + "&_fields=structure&_limit=1")
    data = json.loads(_http_get(url, timeout=15, extra_headers={"X-API-KEY": key}))
    entries = data.get("data", [])
    if not entries:
        raise ValueError(f"Material {material_id} not found")
    structure = entries[0].get("structure", {})
    out_path = str(Path(output_dir) / (material_id.replace("/", "_") + ".cif"))
    lattice = structure.get("lattice", {})
    sites = structure.get("sites", [])
    with open(out_path, "w") as f:
        f.write("data_mp\n")
        for p, v in [("a", "a"), ("b", "b"), ("c", "c")]:
            f.write(f"_cell_length_{p}    {lattice.get(v, 1)}\n")
        for p, v in [("alpha", "alpha"), ("beta", "beta"), ("gamma", "gamma")]:
            f.write(f"_cell_angle_{p} {lattice.get(v, 90)}\n")
        f.write("_symmetry_space_group_name_H-M   'P 1'\n\n")
        f.write("loop_\n_atom_site_label\n_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
        for i, site in enumerate(sites):
            species = site.get("species", [{}])
            el = species[0].get("element", "X") if species else "X"
            abc = site.get("abc", [0, 0, 0])
            f.write(f"{el}{i+1} {el} {abc[0]:.6f} {abc[1]:.6f} {abc[2]:.6f}\n")
    return out_path


# ---------------------------------------------------------------------------
# Backend registry — single source of truth for code, label, search, fetch.
# Lowercase keys flow through both the UI Select and the result dict's
# 'source' field so dispatch can never drift between search and fetch.
#
# Each backend declares its supported search modes:
#   - "text":    free-text query (works on every backend)
#   - "formula": chemical formula (e.g. "C9H8O4", "Cu O")
#   - "id":      direct ID — bypasses search entirely, calls fetch with the
#                user's exact string
# UI is responsible for offering only the modes a backend supports.
# ---------------------------------------------------------------------------

def _search_cod_dispatch(query: str, mode: str = "text") -> list[dict]:
    if mode == "formula":
        return search_cod(formula=query)
    return search_cod(text=query)


BACKENDS: dict[str, dict] = {
    "cod": {
        "label": "COD",
        "modes": ("text", "formula", "id"),
        "placeholder": {
            "text": "anything (eg 'cuprite')",
            "formula": "Cu O",
            "id": "1542931",
        },
        "search": _search_cod_dispatch,
        "fetch": fetch_cod,
    },
    "rcsb": {
        "label": "RCSB PDB",
        "modes": ("text", "id"),
        "placeholder": {
            "text": "villin headpiece",
            "id": "1VII",
        },
        "search": lambda q, mode="text": search_rcsb(q),
        "fetch": fetch_rcsb,
    },
    "pubchem": {
        "label": "PubChem",
        "modes": ("text", "id"),
        "placeholder": {
            "text": "aspirin",
            "id": "2244 (CID)",
        },
        "search": lambda q, mode="text": search_pubchem(q),
        "fetch": fetch_pubchem_sdf,
    },
    "mp": {
        "label": "Materials Project",
        "modes": ("formula", "id"),
        "placeholder": {
            "formula": "Si O2",
            "id": "mp-149",
        },
        "search": lambda q, mode="formula": search_materials_project(q),
        # fetch_materials_project_cif takes (material_id, api_key, output_dir)
        # Bind api_key='' so callers pass (id, output_dir) like the others.
        "fetch": lambda mid, output_dir=".": fetch_materials_project_cif(mid, "", output_dir),
    },
}

# ======================================================================
# Module: sorbates
# ======================================================================
"""MPMC sorbate molecule models from the MPMC repository.

Two categories:
  - Nonpolar LJ: Standard Lennard-Jones sigma/epsilon. No polarizability.
  - Polarizable LJ: LJ sigma/epsilon with explicit point polarizabilities.
    Used with many-body polarization turned on in MPMC.

Note: True PHAHST (damped exponential repulsion + C6/C8/C10 dispersion) sorbate
models are not included here — they require parameters specifically fitted for
that potential form. The PHAHST MOF force field parameters are in forcefields.py.
"""




@dataclass
class SorbateSite:
    atomtype: str
    x: float
    y: float
    z: float
    mass: float
    charge: float
    polarizability: float
    epsilon: float  # LJ epsilon (K) or exp-repulsion beta (A^-1)
    sigma: float    # LJ sigma (A) or exp-repulsion rho (A)
    c6: float = 0.0
    c8: float = 0.0
    c10: float = 0.0


@dataclass
class SorbateModel:
    name: str
    label: str
    description: str
    pot_type: str  # "lj", "lj_polar", or "phahst"
    sites: list[SorbateSite] = field(default_factory=list)


def _s(atomtype: str, x: float, y: float, z: float,
       mass: float, charge: float, pol: float, eps: float, sig: float) -> SorbateSite:
    return SorbateSite(atomtype, x, y, z, mass, charge, pol, eps, sig)


def _ph(atomtype: str, x: float, y: float, z: float,
        mass: float, charge: float, pol: float,
        beta: float, rho: float,
        c6: float, c8: float, c10: float) -> SorbateSite:
    return SorbateSite(atomtype, x, y, z, mass, charge, pol, beta, rho, c6, c8, c10)


SORBATE_MODELS: dict[str, SorbateModel] = {
    # =========================================================================
    # Nonpolar LJ models
    # =========================================================================

    # --- Hydrogen ---
    "h2_bss": SorbateModel("h2_bss", "H2", "H2 BSS 5-site", "lj", [
        _s("H2G", 0, 0, 0, 0, -0.7464, 0, 8.8516, 3.2293),
        _s("H2E", 0.371, 0, 0, 1.008, 0.3732, 0, 0, 0),
        _s("H2E", -0.371, 0, 0, 1.008, 0.3732, 0, 0, 0),
        _s("H2N", 0.329, 0, 0, 0, 0, 0, 4.0659, 2.3406),
        _s("H2N", -0.329, 0, 0, 0, 0, 0, 4.0659, 2.3406),
    ]),
    "h2_buch": SorbateModel("h2_buch", "H2", "H2 Buch 1-site", "lj", [
        _s("H2G", 0, 0, 0, 2.016, 0, 0, 34.2, 2.96),
    ]),
    "h2_dl": SorbateModel("h2_dl", "H2", "H2 Darkrim-Levesque 3-site", "lj", [
        _s("H2G", 0, 0, 0, 0, -0.936, 0, 36.7, 2.958),
        _s("H2E", -0.37, 0, 0, 1.008, 0.468, 0, 0, 0),
        _s("H2E", 0.37, 0, 0, 1.008, 0.468, 0, 0, 0),
    ]),
    # --- CO2 ---
    "co2_trappe": SorbateModel("co2_trappe", "CO2", "CO2 TraPPE 3-site", "lj", [
        _s("COG", 0, 0, 0, 12.01, 0.7, 0, 27.0, 2.80),
        _s("COE", 1.16, 0, 0, 16.0, -0.35, 0, 79.0, 3.05),
        _s("COE", -1.16, 0, 0, 16.0, -0.35, 0, 79.0, 3.05),
    ]),
    "co2_epm2": SorbateModel("co2_epm2", "CO2", "CO2 EPM2 Harris-Yung 3-site", "lj", [
        _s("COG", 0, 0, 0, 12.01, 0.6512, 0, 28.129, 2.757),
        _s("COE", 1.149, 0, 0, 16.0, -0.3256, 0, 80.507, 3.033),
        _s("COE", -1.149, 0, 0, 16.0, -0.3256, 0, 80.507, 3.033),
    ]),
    "co2_phast": SorbateModel("co2_phast", "CO2", "CO2 PHAST 5-site", "lj", [
        _s("COG", 0, 0, 0, 12.0107, 0.77106, 0, 8.52238, 3.05549),
        _s("COE", 1.162, 0, 0, 15.9994, -0.38553, 0, 0, 0),
        _s("COE", -1.162, 0, 0, 15.9994, -0.38553, 0, 0, 0),
        _s("CON", 1.091, 0, 0, 0, 0, 0, 76.76607, 2.94473),
        _s("CON", -1.091, 0, 0, 0, 0, 0, 76.76607, 2.94473),
    ]),
    # --- Methane ---
    "ch4_trappe": SorbateModel("ch4_trappe", "CH4", "CH4 TraPPE 1-site", "lj", [
        _s("CHG", 0, 0, 0, 16.0426, 0, 0, 148.0, 3.73),
    ]),
    "ch4_9site": SorbateModel("ch4_9site", "CH4", "CH4 9-site", "lj", [
        _s("CHG", 0, 0, 0, 12.011, -0.5868, 0, 58.53869, 2.22416),
        _s("CHE", 0, 0, 1.099, 1.0079, 0.1467, 0, 0, 0),
        _s("CHE", 1.036, 0, -0.366, 1.0079, 0.1467, 0, 0, 0),
        _s("CHE", -0.518, -0.897, -0.366, 1.0079, 0.1467, 0, 0, 0),
        _s("CHE", -0.518, 0.897, -0.366, 1.0079, 0.1467, 0, 0, 0),
        _s("MOV", 0, 0, 0.816, 0, 0, 0, 16.85422, 2.96286),
        _s("MOV", 0.769, 0, -0.271, 0, 0, 0, 16.85422, 2.96286),
        _s("MOV", -0.385, -0.668, -0.271, 0, 0, 0, 16.85422, 2.96286),
        _s("MOV", -0.385, 0.668, -0.271, 0, 0, 0, 16.85422, 2.96286),
    ]),
    # --- Nitrogen ---
    "n2_trappe": SorbateModel("n2_trappe", "N2", "N2 TraPPE 3-site", "lj", [
        _s("N2G", 0, 0, 0, 0, 0.964, 0, 0, 0),
        _s("N2E", 0.55, 0, 0, 14.0067, -0.482, 0, 36.0, 3.31),
        _s("N2E", -0.55, 0, 0, 14.0067, -0.482, 0, 36.0, 3.31),
    ]),
    "n2_mcquarrie": SorbateModel("n2_mcquarrie", "N2", "N2 McQuarrie 1-site", "lj", [
        _s("N2G", 0, 0, 0, 28.01344, 0, 0, 95.1, 3.7),
    ]),
    # --- Water ---
    "h2o_tip3p": SorbateModel("h2o_tip3p", "H2O", "H2O TIP3P 3-site", "lj", [
        _s("OXY", 0, 0, 0, 16.0, -0.834, 0, 76.42, 3.151),
        _s("HYD", -0.757, -0.586, 0, 1.008, 0.417, 0, 0, 0),
        _s("HYD", 0.757, -0.586, 0, 1.008, 0.417, 0, 0, 0),
    ]),
    "h2o_tip4p": SorbateModel("h2o_tip4p", "H2O", "H2O TIP4P 4-site", "lj", [
        _s("OXY", 0, 0, 0, 16.0, 0, 0, 78.0, 3.154),
        _s("HYD", 0.58588, 0.75695, 0, 1.008, 0.52, 0, 0, 0),
        _s("HYD", 0.58588, -0.75695, 0, 1.008, 0.52, 0, 0, 0),
        _s("M", 0.15, 0, 0, 0, -1.04, 0, 0, 0),
    ]),
    # --- Noble gases ---
    "he": SorbateModel("he", "He", "Helium pore volume probe", "lj", [
        _s("He", 0, 0, 0, 4.002602, 0, 0, 10.22, 2.28),
    ]),
    # --- O2 ---
    "o2_trappe": SorbateModel("o2_trappe", "O2", "O2 TraPPE 3-site", "lj", [
        _s("O", -0.605, 0, 0, 7.9997, -0.113, 0, 49.0, 3.02),
        _s("CoM", 0, 0, 0, 0, 0.226, 0, 0, 0),
        _s("O", 0.605, 0, 0, 7.9997, -0.113, 0, 49.0, 3.02),
    ]),

    # =========================================================================
    # Polarizable LJ models (LJ + explicit point polarizabilities)
    # =========================================================================

    # --- Hydrogen ---
    "h2_bssp": SorbateModel("h2_bssp", "H2", "H2 BSSP 5-site polarizable", "lj_polar", [
        _s("H2G", 0, 0, 0, 0, -0.7464, 0.6938, 12.76532, 3.15528),
        _s("H2E", 0.371, 0, 0, 1.008, 0.3732, 0.00044, 0, 0),
        _s("H2E", -0.371, 0, 0, 1.008, 0.3732, 0.00044, 0, 0),
        _s("H2N", 0.363, 0, 0, 0, 0, 0, 2.16726, 2.37031),
        _s("H2N", -0.363, 0, 0, 0, 0, 0, 2.16726, 2.37031),
    ]),
    # --- CO2 ---
    "co2_phastp": SorbateModel("co2_phastp", "CO2", "CO2 PHASTP 5-site polarizable", "lj_polar", [
        _s("COG", 0, 0, 0, 12.0107, 0.77134, 1.2281, 19.61757, 3.03366),
        _s("COE", 1.162, 0, 0, 15.9994, -0.38567, 0.7395, 0, 0),
        _s("COE", -1.162, 0, 0, 15.9994, -0.38567, 0.7395, 0, 0),
        _s("CON", 1.208, 0, 0, 0, 0, 0, 46.47457, 2.99429),
        _s("CON", -1.208, 0, 0, 0, 0, 0, 46.47457, 2.99429),
    ]),
    "co2_becker": SorbateModel("co2_becker", "CO2", "CO2 Becker 3-site polarizable", "lj_polar", [
        _s("COG", 0, 0, 0, 12.01, 0.7, 0.916, 23.4, 2.8),
        _s("COE", 1.16, 0, 0, 16.0, -0.35, 0.575, 73.08, 3.05),
        _s("COE", -1.16, 0, 0, 16.0, -0.35, 0.575, 73.08, 3.05),
    ]),
    # --- Methane ---
    "ch4_9sitep": SorbateModel("ch4_9sitep", "CH4", "CH4 9-site polarizable", "lj_polar", [
        _s("CHG", 0, 0, 0, 12.011, -0.5868, 1.0987, 45.0973, 2.16247),
        _s("CHE", 0, 0, 1.099, 1.0079, 0.1467, 0.4246, 0, 0),
        _s("CHE", 1.036, 0, -0.366, 1.0079, 0.1467, 0.4246, 0, 0),
        _s("CHE", -0.518, -0.897, -0.366, 1.0079, 0.1467, 0.4246, 0, 0),
        _s("CHE", -0.518, 0.897, -0.366, 1.0079, 0.1467, 0.4246, 0, 0),
        _s("MOV", 0, 0, 0.814, 0, 0, 0, 18.57167, 2.94787),
        _s("MOV", 0.768, 0, -0.270, 0, 0, 0, 18.57167, 2.94787),
        _s("MOV", -0.383, -0.666, -0.270, 0, 0, 0, 18.57167, 2.94787),
        _s("MOV", -0.383, 0.666, -0.270, 0, 0, 0, 18.57167, 2.94787),
    ]),
    # --- Nitrogen ---
    "n2_polar": SorbateModel("n2_polar", "N2", "N2 5-site polarizable", "lj_polar", [
        _s("N2G", 0, 0, 0, 0, 1.04742, 1.4559, 20.6365, 3.42344),
        _s("N2E", 0.549, 0, 0, 14.0067, -0.52371, 0.5138, 0, 0),
        _s("N2E", -0.549, 0, 0, 14.0067, -0.52371, 0.5138, 0, 0),
        _s("N2N", 0.738, 0, 0, 0, 0, 0, 18.12772, 3.15125),
        _s("N2N", -0.738, 0, 0, 0, 0, 0, 18.12772, 3.15125),
    ]),
    # --- Water ---
    "h2o_pol3": SorbateModel("h2o_pol3", "H2O", "H2O POL3 3-site polarizable", "lj_polar", [
        _s("OXY", 0, 0, 0, 15.999, -0.73, 0.528, 78.50225, 3.596),
        _s("HYD", -0.816, -0.577, 0, 1.0079, 0.365, 0.17, 0, 0),
        _s("HYD", 0.816, -0.577, 0, 1.0079, 0.365, 0.17, 0, 0),
    ]),
    # --- Noble gases (polarizable LJ) ---
    "ar_pol": SorbateModel("ar_pol", "Ar", "Argon polarizable (LJ)", "lj_polar", [
        _s("Ar", 0, 0, 0, 39.948, 0, 1.6392212, 128.326802, 3.371914),
    ]),
    "ne_pol": SorbateModel("ne_pol", "Ne", "Neon polarizable (LJ)", "lj_polar", [
        _s("Ne", 0, 0, 0, 20.1797, 0, 0.3913212, 36.824138, 2.785823),
    ]),
    "kr_pol": SorbateModel("kr_pol", "Kr", "Krypton polarizable (LJ)", "lj_polar", [
        _s("Kr", 0, 0, 0, 83.798, 0, 2.5004096, 183.795833, 3.601271),
    ]),
    "xe_pol": SorbateModel("xe_pol", "Xe", "Xenon polarizable (LJ)", "lj_polar", [
        _s("Xe", 0, 0, 0, 131.293, 0, 4.0232578, 237.985247, 3.956802),
    ]),

    # =========================================================================
    # PHAHST models (damped exp repulsion + C6/C8/C10 dispersion)
    # Hogan & Space, JCTC 2020, 16, 7632-7644
    # beta = exp repulsion steepness (A^-1), rho = exp repulsion distance (A)
    # C6/C8/C10 in atomic units (Hartree * Bohr^n)
    # =========================================================================

    # --- Hydrogen ---
    "h2_phahst": SorbateModel("h2_phahst", "H2", "H2 PHAHST 3-site", "phahst", [
        _ph("H2DA", 0, 0, 0, 0, -0.846166, 0, 3.627796, 2.664506, 0, 0, 0),
        _ph("H2H", 0.371, 0, 0, 1.008, 0.423083, 0.34325, 3.100603, 1.859425, 2.884735, 38.97178, 644.95683),
        _ph("H2H", -0.371, 0, 0, 1.008, 0.423083, 0.34325, 3.100603, 1.859425, 2.884735, 38.97178, 644.95683),
    ]),
    # --- Nitrogen ---
    "n2_phahst": SorbateModel("n2_phahst", "N2", "N2 PHAHST 3-site", "phahst", [
        _ph("N2DA", 0, 0, 0, 0, 0.94194, 0, 0, 0, 0, 0, 0),
        _ph("N2N", 0.5507, 0, 0, 14.0067, -0.47103, 0.85092, 3.85368, 3.31513, 17.80503, 416.3235, 11924.913),
        _ph("N2N", -0.5507, 0, 0, 14.0067, -0.47103, 0.85092, 3.85368, 3.31513, 17.80503, 416.3235, 11924.913),
    ]),
    # --- Noble gases ---
    "he_phahst": SorbateModel("he_phahst", "He", "Helium PHAHST", "phahst", [
        _ph("He", 0, 0, 0, 4.0026, 0, 0.2002, 4.68451, 2.38376, 1.40716, 11.13635, 107.964),
    ]),
    "ne_phahst": SorbateModel("ne_phahst", "Ne", "Neon PHAHST", "phahst", [
        _ph("Ne", 0, 0, 0, 20.18, 0, 0.3823, 4.99432, 2.8024, 6.21275, 67.98647, 911.376),
    ]),
    "ar_phahst": SorbateModel("ar_phahst", "Ar", "Argon PHAHST", "phahst", [
        _ph("Ar", 0, 0, 0, 39.948, 0, 1.655, 3.88525, 3.68623, 65.46, 1438.9, 38745.0),
    ]),
    "kr_phahst": SorbateModel("kr_phahst", "Kr", "Krypton PHAHST", "phahst", [
        _ph("Kr", 0, 0, 0, 83.798, 0, 2.497, 3.52894, 4.03477, 130.1, 3981.0, 149225.0),
    ]),
    "xe_phahst": SorbateModel("xe_phahst", "Xe", "Xenon PHAHST", "phahst", [
        _ph("Xe", 0, 0, 0, 131.293, 0, 4.026, 3.24691, 4.47518, 288.4, 11390.0, 551047.0),
    ]),
}


def get_sorbate_names() -> list[tuple[str, str, str]]:
    """Return list of (model_name, description, pot_type) for all sorbates."""
    return [(name, m.description, m.pot_type) for name, m in SORBATE_MODELS.items()]


def format_sorbate_pqr(model_name: str, center_x: float, center_y: float,
                       center_z: float, mol_id: int, start_atom_id: int) -> list[str]:
    """Generate MPMC PQR lines for a sorbate placed at the given center coordinates."""
    model = SORBATE_MODELS[model_name]
    lines = []
    for i, site in enumerate(model.sites):
        atom_id = start_atom_id + i
        x = center_x + site.x
        y = center_y + site.y
        z = center_z + site.z
        lines.append(
            f"ATOM {atom_id:>6} {site.atomtype:<4} {model.label:<3} M {mol_id:>4}    "
            f"{x:>10.6f} {y:>10.6f} {z:>10.6f}"
            f" {site.mass:>9.5f} {site.charge:>9.5f}"
            f" {site.polarizability:>9.5f} {site.epsilon:>9.5f} {site.sigma:>9.5f}"
            f" 0.0 0.0 {site.c6:>10.5f} {site.c8:>12.5f} {site.c10:>14.5f}"
        )
    return lines

# ======================================================================
# Module: panels
# ======================================================================
"""TUI sidebar panels: OperationsPanel, GeometryPanel, and VisualPanel."""




# ---------------------------------------------------------------------------
# OperationsPanel
# ---------------------------------------------------------------------------

# Mirrors the top-bar menus (File, Edit, Analysis) so users have one source
# of truth for every available tool. View-only display toggles live in the
# top-bar View menu and the VisualPanel; they are intentionally NOT in the
# Operations panel because they're not "run-once" actions.
OPERATIONS = [
    # ---- Analysis (mirrors top bar Analysis menu) ----
    ("info",           "System Info",              "analysis"),
    ("void_vol",       "Void Volume",              "analysis"),
    ("surface_area",   "Surface Area",             "analysis"),
    ("pore_size",      "Pore Size Distribution",   "analysis"),
    ("rdf",            "Radial Distribution",      "analysis"),
    ("coordination",   "Coordination Number",      "analysis"),
    ("hbonds",         "Hydrogen Bonds",           "analysis"),
    ("msd",            "Mean Sq. Displacement",    "analysis"),
    ("rmsd",           "RMSD vs Time",             "analysis"),
    ("gyration",       "Radius of Gyration",       "analysis"),
    ("density",        "Density Profile",          "analysis"),
    ("density3d",      "3D Density Map",           "analysis"),
    ("pxrd",           "Powder XRD",               "analysis"),
    ("energy_plot",    "Energy Plot",              "analysis"),
    ("input_generator","Input Generator",          "analysis"),
    ("isotherm_plot",  "Plot Isotherm Results",    "analysis"),
    ("db_search",      "Fetch from Database",      "analysis"),
    # ---- Edit (mirrors top bar Edit menu) ----
    ("undo",           "Undo",                     "modify"),
    ("wrap_center",    "Wrap (centered)",          "modify"),
    ("wrap_forward",   "Wrap (forward)",           "modify"),
    ("sort",           "Sort Atoms",               "modify"),
    ("extend",         "Extend Axis",              "modify"),
    ("del_lone",       "Delete Lone Atoms",        "modify"),
    ("edit_h",         "Edit H Distances",         "modify"),
    ("update_cell",    "Update Unit Cell",         "modify"),
    ("substitute",     "Substitute Element",       "modify"),
    ("reduce_cell",    "Reduce Supercell",         "modify"),
    ("qeq_charges",    "Generate QEq Charges",     "modify"),
    # ---- File (mirrors top bar File menu) ----
    ("open_file",      "Open File...",             "export"),
    ("write_xyz",      "Save as XYZ",              "export"),
    ("write_pdb",      "Save as PDB",              "export"),
    ("write_cif",      "Save as CIF",              "export"),
    ("write_poscar",   "Save as POSCAR",           "export"),
    ("write_lammps",   "Save as LAMMPS Data",      "export"),
    ("write_mpmc",     "Save as MPMC PDB",         "export"),
    ("frame_to_tab",   "Open Frame in New Tab",    "export"),
    ("export_png",     "Export PNG",               "export"),
    ("export_gif",     "Export Rotation GIF",      "export"),
    ("quit",           "Quit",                     "export"),
]


class OperationsPanel(Widget):
    BINDINGS = [
        Binding("enter", "run_selected", "Run", show=True),
    ]
    DEFAULT_CSS = """
    OperationsPanel {
        dock: left;
        width: 30;
        display: none;
        border-right: solid $accent;
    }
    OperationsPanel.visible {
        display: block;
    }
    OperationsPanel DataTable {
        height: 1fr;
    }
    """

    class RunCommand(Message):
        def __init__(self, command: str) -> None:
            super().__init__()
            self.command = command

    def __init__(self) -> None:
        super().__init__()
        self._populated = False

    def compose(self) -> ComposeResult:
        yield DataTable(id="ops-table", cursor_type="row")

    def on_mount(self) -> None:
        if not self._populated:
            self._populate()

    def _populate(self) -> None:
        table = self.query_one("#ops-table", DataTable)
        table.clear(columns=True)
        table.add_columns("", "Command")
        section = ""
        for cmd_id, label, group in OPERATIONS:
            if group != section:
                section = group
                header = {"analysis": "--- Analysis ---", "modify": "--- Edit ---", "export": "--- File ---"}
                table.add_row("", header.get(group, ""), key=f"_header_{group}")
            table.add_row(chr(0x25b6), label, key=cmd_id)
        self._populated = True

    def action_run_selected(self) -> None:
        table = self.query_one("#ops-table", DataTable)
        if table.row_count == 0:
            return
        row_keys = list(table.rows.keys())
        row = max(0, min(table.cursor_row, len(row_keys) - 1))
        key = row_keys[row].value
        if key and not key.startswith("_header_"):
            self.post_message(self.RunCommand(key))

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        if event.row_key.value and not event.row_key.value.startswith("_header_"):
            self.post_message(self.RunCommand(event.row_key.value))


# ---------------------------------------------------------------------------
# GeometryPanel
# ---------------------------------------------------------------------------

class GeometryPanel(Widget):
    BINDINGS = [
        Binding("tab", "next_tab", "Tab", show=True),
        Binding("shift+tab", "prev_tab", "Prev tab", show=False),
        Binding("s", "toggle_sort", "Sort", show=True),
    ]

    _TAB_IDS = [
        "tab-bonds", "tab-angles", "tab-dihedrals",
        "tab-contacts", "tab-molecules", "tab-lone", "tab-coords",
    ]

    DEFAULT_CSS = """
    GeometryPanel {
        dock: right;
        width: 50;
        display: none;
        border-left: solid $accent;
    }
    GeometryPanel.visible {
        display: block;
    }
    GeometryPanel DataTable {
        height: 1fr;
    }
    """

    class HighlightAtoms(Message):
        def __init__(self, atom_indices: tuple[int, ...]) -> None:
            super().__init__()
            self.atom_indices = atom_indices

    class EditCoord(Message):
        def __init__(self, atom_index: int, x: float, y: float, z: float) -> None:
            super().__init__()
            self.atom_index = atom_index
            self.x = x
            self.y = y
            self.z = z

    def __init__(self) -> None:
        super().__init__()
        self._molecule: Molecule | None = None
        self._populating = False
        self._sort_ascending: dict[str, bool] = {t: False for t in self._TAB_IDS}

    def set_molecule(self, molecule: Molecule) -> None:
        self._molecule = molecule
        # Mark all tabs as needing refresh; compute lazily on demand
        self._populated_tabs: set[str] = set()
        if self.is_mounted and self.has_class("visible"):
            self._populate_active_tab()

    def on_mount(self) -> None:
        if not hasattr(self, "_populated_tabs"):
            self._populated_tabs = set()
        if self._molecule is not None and self.has_class("visible"):
            self._populate_active_tab()

    def on_show(self) -> None:
        """Populate the active tab when the panel becomes visible."""
        if self._molecule is not None:
            self._populate_active_tab()

    def _populate_active_tab(self) -> None:
        """Populate only the currently active tab's data."""
        if self._molecule is None:
            return
        try:
            tabs = self.query_one(TabbedContent)
            active = tabs.active
        except Exception:
            return
        self._populate_tab(active)

    def _populate_tab(self, tab_id: str, force: bool = False) -> None:
        """Populate a single tab's data on demand."""
        if self._molecule is None:
            return
        if not force and tab_id in self._populated_tabs:
            return
        self._populated_tabs.add(tab_id)
        self._populating = True
        mol = self._molecule
        if tab_id == "tab-bonds":
            self._populate_bonds(mol)
        elif tab_id == "tab-angles":
            self._populate_angles(mol)
        elif tab_id == "tab-dihedrals":
            self._populate_dihedrals(mol)
        elif tab_id == "tab-contacts":
            self._populate_contacts(mol)
        elif tab_id == "tab-molecules":
            self._populate_submols(mol)
        elif tab_id == "tab-lone":
            self._populate_lone(mol)
        elif tab_id == "tab-coords":
            self._populate_coords(mol)
        self._populating = False

    def action_next_tab(self) -> None:
        tabs = self.query_one(TabbedContent)
        current = tabs.active
        idx = self._TAB_IDS.index(current) if current in self._TAB_IDS else 0
        tabs.active = self._TAB_IDS[(idx + 1) % len(self._TAB_IDS)]

    def action_prev_tab(self) -> None:
        tabs = self.query_one(TabbedContent)
        current = tabs.active
        idx = self._TAB_IDS.index(current) if current in self._TAB_IDS else 0
        tabs.active = self._TAB_IDS[(idx - 1) % len(self._TAB_IDS)]

    def action_toggle_sort(self) -> None:
        if self._molecule is None:
            return
        tabs = self.query_one(TabbedContent)
        active = tabs.active
        selected_row_key = self._current_row_key_value(self._table_for_tab(active))
        asc = self._sort_ascending.get(active, False)
        self._sort_ascending[active] = not asc
        self._populate_tables({active: selected_row_key} if selected_row_key is not None else None)
        table = self._table_for_tab(active)
        table.focus()
        self._emit_current_highlight(table)

    def compose(self) -> ComposeResult:
        with TabbedContent("Bonds", "Angles", "Dihedrals", "Contacts", "Mols", "Lone", "Coords"):
            with TabPane("Bonds", id="tab-bonds"):
                yield DataTable(id="bonds-table", cursor_type="row")
            with TabPane("Angles", id="tab-angles"):
                yield DataTable(id="angles-table", cursor_type="row")
            with TabPane("Dihedrals", id="tab-dihedrals"):
                yield DataTable(id="dihedrals-table", cursor_type="row")
            with TabPane("Contacts", id="tab-contacts"):
                yield DataTable(id="contacts-table", cursor_type="row")
            with TabPane("Mols", id="tab-molecules"):
                yield DataTable(id="molecules-table", cursor_type="row")
            with TabPane("Lone", id="tab-lone"):
                yield DataTable(id="lone-table", cursor_type="row")
            with TabPane("Coords", id="tab-coords"):
                yield DataTable(id="coords-table", cursor_type="row")

    def _atom_label(self, idx: int) -> str:
        if self._molecule is None:
            return str(idx + 1)
        return f"{idx + 1}:{self._molecule.atoms[idx].element.symbol}"

    def _table_for_tab(self, tab_id: str) -> DataTable:
        table_id = {
            "tab-bonds": "#bonds-table",
            "tab-angles": "#angles-table",
            "tab-dihedrals": "#dihedrals-table",
            "tab-contacts": "#contacts-table",
            "tab-molecules": "#molecules-table",
            "tab-lone": "#lone-table",
            "tab-coords": "#coords-table",
        }.get(tab_id, "#bonds-table")
        return self.query_one(table_id, DataTable)

    def _current_row_key_value(self, dt: DataTable) -> str | None:
        if dt.row_count == 0:
            return None
        row_keys = list(dt.rows.keys())
        if not row_keys:
            return None
        row = max(0, min(dt.cursor_row, len(row_keys) - 1))
        return row_keys[row].value

    def _restore_cursor(self, dt: DataTable, row_key: str | None) -> None:
        if row_key is None or dt.row_count == 0:
            return
        row_keys = list(dt.rows.keys())
        for row, key in enumerate(row_keys):
            if key.value == row_key:
                dt.move_cursor(row=row, scroll=False)
                return

    def _populate_tables(self, selected_row_keys: dict[str, str | None] | None = None) -> None:
        """Force-repopulate all tabs (used by sort and explicit refresh)."""
        if self._molecule is None:
            return
        self._populated_tabs = set()
        for tid in self._TAB_IDS:
            self._populate_tab(tid, force=True)

    def _populate_bonds(self, mol) -> None:
        bonds = getattr(mol, '_cached_bond_lengths', None) or mol.get_bond_lengths()
        if self._sort_ascending.get("tab-bonds"):
            bonds.sort(key=lambda x: x[2])
        table = self.query_one("#bonds-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Atom 2", "Length (\u00c5)")
        for i, j, dist in bonds:
            table.add_row(self._atom_label(i), self._atom_label(j), f"{dist:.4f}", key=f"{i}-{j}")

    def _populate_angles(self, mol) -> None:
        angles = getattr(mol, '_cached_angles', None) or mol.get_angles()
        if self._sort_ascending.get("tab-angles"):
            angles.sort(key=lambda x: x[3])
        table = self.query_one("#angles-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Vertex", "Atom 3", "Angle (\u00b0)")
        for i, j, k, angle in angles:
            table.add_row(
                self._atom_label(i), self._atom_label(j), self._atom_label(k),
                f"{angle:.3f}", key=f"{i}-{j}-{k}",
            )

    def _populate_dihedrals(self, mol) -> None:
        dihedrals = getattr(mol, '_cached_dihedrals', None) or mol.get_dihedrals()
        if self._sort_ascending.get("tab-dihedrals"):
            dihedrals.sort(key=lambda x: x[4])
        table = self.query_one("#dihedrals-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Atom 2", "Atom 3", "Atom 4", "Angle (\u00b0)")
        for i, j, k, l_idx, angle in dihedrals:
            table.add_row(
                self._atom_label(i), self._atom_label(j),
                self._atom_label(k), self._atom_label(l_idx),
                f"{angle:.3f}", key=f"{i}-{j}-{k}-{l_idx}",
            )

    def _populate_contacts(self, mol) -> None:
        contacts = getattr(mol, '_cached_contacts', None)
        if contacts is None:
            set_atom_ids(mol.atoms)
            contacts = get_close_contacts(mol)
        table = self.query_one("#contacts-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Atom 2", "Distance (\u00c5)")
        for msg in contacts:
            parts = msg.split()
            el_pair = parts[0]
            id1, id2 = int(parts[1]), int(parts[2])
            r_val = parts[-1]
            table.add_row(
                f"{id1}:{el_pair.split('-')[0]}", f"{id2}:{el_pair.split('-')[1]}",
                r_val, key=f"{id1 - 1}-{id2 - 1}",
            )

    def _populate_submols(self, mol) -> None:
        submols = getattr(mol, '_cached_molecules', None) or mol.find_molecules()
        table = self.query_one("#molecules-table", DataTable)
        table.clear(columns=True)
        table.add_columns("#", "Atoms", "Formula")
        from collections import Counter
        for mi, sm in enumerate(submols):
            sm.atoms.sort(key=lambda a: a.atomic_number, reverse=True)
            counts = Counter(a.element.symbol for a in sm.atoms)
            formula = "".join(f"{el}{n}" if n > 1 else el for el, n in counts.items())
            atom_indices = "-".join(str(mol.atoms.index(a)) for a in sm.atoms)
            table.add_row(str(mi + 1), str(len(sm.atoms)), formula, key=atom_indices)

    def _populate_lone(self, mol) -> None:
        lone = getattr(mol, '_cached_lone_atoms', None)
        if lone is None:
            lone = get_lone_atoms(mol)
        table = self.query_one("#lone-table", DataTable)
        table.clear(columns=True)
        table.add_columns("#", "Element", "x", "y", "z")
        for a in lone:
            idx = mol.atoms.index(a)
            table.add_row(
                str(a.id), a.element.symbol,
                f"{a.x[0]:.4f}", f"{a.x[1]:.4f}", f"{a.x[2]:.4f}",
                key=str(idx),
            )

    def _populate_coords(self, mol) -> None:
        table = self.query_one("#coords-table", DataTable)
        table.clear(columns=True)
        table.add_columns("#", "El", "x", "y", "z")
        for idx, a in enumerate(mol.atoms):
            table.add_row(
                str(a.id), a.element.symbol,
                f"{a.x[0]:.4f}", f"{a.x[1]:.4f}", f"{a.x[2]:.4f}",
                key=str(idx),
            )

    def _emit_current_highlight(self, dt: DataTable) -> None:
        if not self.has_class("visible") or dt.row_count == 0:
            return
        rk = list(dt.rows.keys())[dt.cursor_row]
        if rk.value is not None:
            indices = tuple(int(x) for x in rk.value.split("-"))
            self.post_message(self.HighlightAtoms(indices))

    def on_tabbed_content_tab_activated(self, event: TabbedContent.TabActivated) -> None:
        # Lazily populate the activated tab
        if event.pane.id and self._molecule is not None and self.has_class("visible"):
            self._populate_tab(event.pane.id)
        for dt in event.pane.query(DataTable):
            dt.focus()
            self._emit_current_highlight(dt)
            break

    def on_data_table_row_highlighted(self, event: DataTable.RowHighlighted) -> None:
        if self._populating or not self.has_class("visible"):
            return
        if event.row_key is None or event.row_key.value is None:
            return
        indices = tuple(int(x) for x in event.row_key.value.split("-"))
        self.post_message(self.HighlightAtoms(indices))

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        if event.data_table.id != "coords-table":
            return
        if event.row_key is None or event.row_key.value is None:
            return
        if self._molecule is None:
            return
        idx = int(event.row_key.value)
        a = self._molecule.atoms[idx]
        self.post_message(self.EditCoord(idx, float(a.x[0]), float(a.x[1]), float(a.x[2])))


# ---------------------------------------------------------------------------
# Slider widget
# ---------------------------------------------------------------------------

class _NavRadioSet(RadioSet):
    BINDINGS = [
        Binding("tab", "next_and_toggle", "Next", show=False),
        Binding("shift+tab", "prev_and_toggle", "Prev", show=False),
    ]

    def action_next_and_toggle(self) -> None:
        self.action_next_button()
        self.action_toggle_button()

    def action_prev_and_toggle(self) -> None:
        self.action_previous_button()
        self.action_toggle_button()


class Slider(Static, can_focus=True):
    BINDINGS = [
        Binding("tab", "increase", "Increase", show=False),
        Binding("shift+tab", "decrease", "Decrease", show=False),
    ]
    DEFAULT_CSS = """
    Slider {
        height: 1;
        width: 1fr;
        padding: 0 1;
    }
    Slider:focus {
        background: $accent 30%;
        text-style: bold;
    }
    """

    def __init__(
        self, label: str, value: float = 0.5,
        min_val: float = 0.0, max_val: float = 1.0,
        step: float = 0.05, decimals: int = 2, **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.label = label
        self.value = value
        self.min_val = min_val
        self.max_val = max_val
        self.step = step
        self.decimals = decimals

    class Changed(Message):
        def __init__(self, slider: Slider, value: float) -> None:
            super().__init__()
            self.slider = slider
            self.value = value

    def render(self) -> str:
        bar_width = 10
        frac = (self.value - self.min_val) / max(self.max_val - self.min_val, 1e-9)
        filled = int(frac * bar_width)
        bar = chr(0x2588) * filled + chr(0x2591) * (bar_width - filled)
        prefix = chr(0x25b8) + " " if self.has_focus else "  "
        arrows = " " + chr(0x25c0) + chr(0x25b6) if self.has_focus else ""
        return f"{prefix}{self.label}: {self.value:.{self.decimals}f} [{bar}]{arrows}"

    def _adjust(self, delta: float) -> None:
        new = max(self.min_val, min(self.max_val, self.value + delta))
        if new != self.value:
            self.value = new
            self.refresh()
            self.post_message(self.Changed(self, self.value))

    def action_decrease(self) -> None:
        self._adjust(-self.step)

    def action_increase(self) -> None:
        self._adjust(self.step)


# ---------------------------------------------------------------------------
# VisualPanel (MO isovalue removed)
# ---------------------------------------------------------------------------

class VisualPanel(Widget):
    DEFAULT_CSS = """
    VisualPanel {
        dock: right;
        width: 30;
        display: none;
        border-left: solid $accent;
        padding: 1;
    }
    VisualPanel.visible {
        display: block;
    }
    VisualPanel Label {
        margin-top: 1;
        text-style: bold;
    }
    VisualPanel RadioSet {
        height: auto;
        margin-bottom: 1;
    }
    VisualPanel #visual-help {
        dock: bottom;
        height: auto;
        color: $text-muted;
        margin-top: 1;
    }
    """

    class StyleChanged(Message):
        def __init__(self, licorice: bool, vdw: bool, ribbon: bool = False) -> None:
            super().__init__()
            self.licorice = licorice
            self.vdw = vdw
            self.ribbon = ribbon

    class LightingChanged(Message):
        def __init__(self, ambient: float, diffuse: float, specular: float, shininess: float) -> None:
            super().__init__()
            self.ambient = ambient
            self.diffuse = diffuse
            self.specular = specular
            self.shininess = shininess

    class SizeChanged(Message):
        def __init__(self, atom_scale: float, bond_radius: float) -> None:
            super().__init__()
            self.atom_scale = atom_scale
            self.bond_radius = bond_radius

    def __init__(self) -> None:
        super().__init__()
        self._licorice = False
        self._vdw = False
        self._ribbon = False

    def set_state(
        self, *, licorice: bool, vdw: bool = False, ribbon: bool = False,
        ambient: float, diffuse: float, specular: float, shininess: float,
        atom_scale: float, bond_radius: float,
    ) -> None:
        self._licorice = licorice
        self._vdw = vdw
        self._ribbon = ribbon
        if self.is_mounted:
            self._sync_widgets(
                ambient=ambient, diffuse=diffuse, specular=specular,
                shininess=shininess, atom_scale=atom_scale, bond_radius=bond_radius,
            )

    def _sync_widgets(
        self, *, ambient: float, diffuse: float, specular: float,
        shininess: float, atom_scale: float, bond_radius: float,
    ) -> None:
        radio_set = self.query_one(_NavRadioSet)
        idx = 1 if self._licorice else (2 if self._vdw else (3 if self._ribbon else 0))
        radio_set.query(RadioButton)[idx].value = True
        self.query_one("#slider-atom-scale", Slider).value = atom_scale
        self.query_one("#slider-bond-radius", Slider).value = bond_radius
        self.query_one("#slider-ambient", Slider).value = ambient
        self.query_one("#slider-diffuse", Slider).value = diffuse
        self.query_one("#slider-specular", Slider).value = specular
        self.query_one("#slider-shininess", Slider).value = shininess
        self._update_visibility()
        self.refresh()

    class ThemeChanged(Message):
        def __init__(self, dark: bool) -> None:
            super().__init__()
            self.dark = dark

    class ToggleChanged(Message):
        def __init__(self, toggle: str, value: bool) -> None:
            super().__init__()
            self.toggle = toggle
            self.value = value

    class ColorModeChanged(Message):
        def __init__(self, mode: str) -> None:
            super().__init__()
            self.mode = mode

    def compose(self) -> ComposeResult:
        yield Label("Theme")
        with RadioSet(id="theme-set"):
            yield RadioButton("Dark", value=True, id="radio-dark")
            yield RadioButton("Light", id="radio-light")
        yield Label("Display")
        yield Checkbox("Show bonds", value=True, id="chk-bonds")
        yield Checkbox("Atom numbers", value=False, id="chk-atomnums")
        yield Checkbox("Hide water", value=False, id="chk-hidewater")
        yield Checkbox("Polyhedra", value=False, id="chk-polyhedra")
        yield Label("Color by")
        yield Select(
            [
                ("Element", "element"),
                ("Charge", "charge"),
                ("Residue", "residue"),
                ("Chain", "chain"),
                ("Index", "index"),
            ],
            value="element",
            id="sel-colormode",
        )
        yield Label("Style")
        with _NavRadioSet():
            yield RadioButton("CPK", value=True, id="radio-cpk")
            yield RadioButton("Licorice", id="radio-licorice")
            yield RadioButton("VDW", id="radio-vdw")
            yield RadioButton("Ribbon", id="radio-ribbon")
        yield Label("Sizes", id="label-sizes")
        yield Slider("Atom scale", value=0.35, min_val=0.10, max_val=1.00, id="slider-atom-scale")
        yield Slider("Bond radius", value=0.08, min_val=0.02, max_val=0.30, step=0.02, id="slider-bond-radius")
        yield Label("Lighting")
        yield Slider("Ambient", value=0.35, min_val=0.0, max_val=1.0, id="slider-ambient")
        yield Slider("Diffuse", value=0.60, min_val=0.0, max_val=1.0, id="slider-diffuse")
        yield Slider("Specular", value=0.40, min_val=0.0, max_val=1.0, id="slider-specular")
        yield Slider("Shininess", value=32.0, min_val=1.0, max_val=128.0, step=4.0, id="slider-shininess")
        yield Static("n/p nav; (shift-)tab toggle", id="visual-help")

    def _update_visibility(self) -> None:
        self.query_one("#slider-atom-scale", Slider).display = not self._licorice and not self._vdw
        self.query_one("#slider-bond-radius", Slider).display = not self._vdw
        self.query_one("#label-sizes", Label).display = not self._vdw

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        cid = event.checkbox.id or ""
        if cid == "chk-bonds":
            self.post_message(self.ToggleChanged("bonds", event.value))
        elif cid == "chk-atomnums":
            self.post_message(self.ToggleChanged("atom_numbers", event.value))
        elif cid == "chk-hidewater":
            self.post_message(self.ToggleChanged("hide_water", event.value))
        elif cid == "chk-polyhedra":
            self.post_message(self.ToggleChanged("polyhedra", event.value))

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id == "sel-colormode" and event.value is not None:
            self.post_message(self.ColorModeChanged(str(event.value)))

    def on_radio_set_changed(self, event: RadioSet.Changed) -> None:
        if event.pressed.id in ("radio-dark", "radio-light"):
            self.post_message(self.ThemeChanged(event.pressed.id == "radio-dark"))
        else:
            self._licorice = event.pressed.id == "radio-licorice"
            self._vdw = event.pressed.id == "radio-vdw"
            self._ribbon = event.pressed.id == "radio-ribbon"
            self._update_visibility()
            self.post_message(self.StyleChanged(self._licorice, self._vdw, self._ribbon))

    def on_slider_changed(self, event: Slider.Changed) -> None:
        sid = event.slider.id or ""
        if sid.startswith("slider-atom") or sid.startswith("slider-bond"):
            self.post_message(
                self.SizeChanged(
                    atom_scale=self.query_one("#slider-atom-scale", Slider).value,
                    bond_radius=self.query_one("#slider-bond-radius", Slider).value,
                )
            )
        else:
            self.post_message(
                self.LightingChanged(
                    ambient=self.query_one("#slider-ambient", Slider).value,
                    diffuse=self.query_one("#slider-diffuse", Slider).value,
                    specular=self.query_one("#slider-specular", Slider).value,
                    shininess=self.query_one("#slider-shininess", Slider).value,
                )
            )

# ======================================================================
# Module: viewer
# ======================================================================
"""PDB Wizard TUI viewer (adapted from moltui, MO code removed)."""



if not hasattr(asyncio, "to_thread"):
    # Python 3.8 backport of asyncio.to_thread (added in 3.9). pdb_wizard
    # targets >=3.8 and offloads analysis compute via to_thread throughout.
    import contextvars
    import functools

    async def _to_thread(func, *args, **kwargs):
        loop = asyncio.get_running_loop()
        ctx = contextvars.copy_context()
        call = functools.partial(ctx.run, func, *args, **kwargs)
        return await loop.run_in_executor(None, call)

    asyncio.to_thread = _to_thread

import numpy as np


_BRAILLE_MAP = np.array(
    [
        [0x01, 0x08],
        [0x02, 0x10],
        [0x04, 0x20],
        [0x40, 0x80],
    ],
    dtype=np.uint8,
)


class InputModal(ModalScreen[Optional[str]]):
    """Modal dialog that prompts for a text value with OK/Cancel."""

    # No 'q' binding — the focused Input would capture it as text. Escape only.
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    DEFAULT_CSS = """
    InputModal {
        align: center middle;
    }
    InputModal > Vertical {
        width: 60%;
        min-width: 40;
        max-width: 80;
        height: auto;
        max-height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    InputModal Label {
        margin-bottom: 1;
    }
    InputModal Input {
        margin-bottom: 1;
    }
    InputModal .hint {
        color: $text-muted;
        margin-bottom: 1;
    }
    InputModal Horizontal {
        height: auto;
        align-horizontal: center;
    }
    InputModal Button {
        width: auto;
        min-width: 10;
        margin: 0 1;
    }
    """

    def __init__(self, prompt: str, default: str = "", hint: str = "") -> None:
        super().__init__()
        self._prompt = prompt
        self._default = default
        self._hint = hint

    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label(self._prompt)
            yield Input(value=self._default, id="modal-input")
            if self._hint:
                yield Label(self._hint, classes="hint")
            with Horizontal():
                yield Button("OK", id="input-ok", variant="primary")
                yield Button("Cancel", id="input-cancel")

    def on_mount(self) -> None:
        self.query_one("#modal-input", Input).focus()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self.dismiss(event.value)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "input-ok":
            self.dismiss(self.query_one("#modal-input", Input).value)
        else:
            self.dismiss("")

    def action_cancel(self) -> None:
        self.dismiss("")


class InfoModal(ModalScreen[None]):
    """Modal dialog showing copyable text with an OK button."""

    BINDINGS = [
        Binding("escape", "dismiss_modal", "Close"),
        Binding("enter", "dismiss_modal", "OK"),
    ]
    DEFAULT_CSS = """
    InfoModal {
        align: center middle;
    }
    InfoModal > Vertical {
        width: 70%;
        min-width: 40;
        max-width: 90;
        height: 70%;
        min-height: 12;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    InfoModal Label {
        height: 1;
        margin-bottom: 1;
        text-style: bold;
    }
    InfoModal TextArea {
        height: 1fr;
        margin-bottom: 1;
    }
    InfoModal Button {
        dock: bottom;
        width: auto;
        min-width: 10;
        margin-top: 1;
        align-horizontal: center;
    }
    """

    def __init__(self, title: str, content: str) -> None:
        super().__init__()
        self._title = title
        self._content = content

    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label(self._title)
            yield TextArea(self._content, read_only=True, id="info-text")
            yield Button("OK", id="info-ok", variant="primary")

    def on_mount(self) -> None:
        self.query_one("#info-ok", Button).focus()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        self.dismiss(None)

    def action_dismiss_modal(self) -> None:
        self.dismiss(None)


class ConfirmModal(ModalScreen[bool]):
    """Modal OK/Cancel confirmation dialog."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("enter", "confirm", "OK"),
    ]
    DEFAULT_CSS = """
    ConfirmModal {
        align: center middle;
    }
    ConfirmModal > Vertical {
        width: 60%;
        min-width: 36;
        max-width: 70;
        height: auto;
        max-height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    ConfirmModal Label {
        margin-bottom: 1;
        text-style: bold;
    }
    ConfirmModal .description {
        margin-bottom: 1;
        text-style: none;
    }
    ConfirmModal Horizontal {
        height: auto;
        align-horizontal: center;
    }
    ConfirmModal Button {
        width: auto;
        min-width: 10;
        margin: 0 1;
    }
    """

    def __init__(self, title: str, description: str, yes_no: bool = False) -> None:
        super().__init__()
        self._title = title
        self._description = description
        self._yes_no = yes_no

    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label(self._title)
            yield Label(self._description, classes="description")
            with Horizontal():
                ok_label = "Yes" if self._yes_no else "OK"
                cancel_label = "No" if self._yes_no else "Cancel"
                yield Button(ok_label, id="confirm-ok", variant="primary")
                yield Button(cancel_label, id="confirm-cancel")

    def on_mount(self) -> None:
        self.query_one("#confirm-ok", Button).focus()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        self.dismiss(event.button.id == "confirm-ok")

    def action_confirm(self) -> None:
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)


class FileSaveModal(ModalScreen[Optional[str]]):
    """Modal file browser for choosing a save location."""

    # No 'q' binding — the focused Input/DirectoryTree would capture it. Escape only.
    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]
    DEFAULT_CSS = """
    FileSaveModal {
        align: center middle;
    }
    FileSaveModal > Vertical {
        width: 80%;
        min-width: 45;
        max-width: 120;
        height: 80%;
        min-height: 14;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    FileSaveModal DirectoryTree {
        height: 1fr;
        margin-bottom: 1;
    }
    FileSaveModal Input {
        margin-bottom: 1;
    }
    FileSaveModal Horizontal {
        height: auto;
        align-horizontal: center;
    }
    FileSaveModal Button {
        width: auto;
        min-width: 10;
        margin: 0 1;
    }
    """

    def __init__(self, prompt: str, default: str = "", button_label: str = "Save") -> None:
        super().__init__()
        self._prompt = prompt
        self._default = default
        self._button_label = button_label

    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label(self._prompt)
            yield DirectoryTree(str(Path(self._default).parent or "."), id="file-tree")
            yield Input(value=self._default, id="file-input")
            with Horizontal():
                yield Button(self._button_label, id="file-save", variant="primary")
                yield Button("Cancel", id="file-cancel")

    def on_mount(self) -> None:
        self.query_one("#file-input", Input).focus()

    def on_directory_tree_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        inp = self.query_one("#file-input", Input)
        inp.value = str(event.path)
        inp.focus()

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self.dismiss(event.value)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "file-save":
            self.dismiss(self.query_one("#file-input", Input).value)
        else:
            self.dismiss("")

    def action_cancel(self) -> None:
        self.dismiss("")


class ForceFieldModal(ModalScreen[Optional[str]]):
    """Modal to select a force field with a live parameter preview."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("q", "cancel", "Close"),
    ]
    DEFAULT_CSS = """
    ForceFieldModal {
        align: center middle;
    }
    ForceFieldModal > Vertical {
        width: 75%;
        min-width: 45;
        max-width: 100;
        height: 75%;
        min-height: 16;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    ForceFieldModal Label {
        margin-bottom: 1;
    }
    ForceFieldModal .title {
        text-style: bold;
    }
    ForceFieldModal RadioSet {
        height: auto;
        margin-bottom: 1;
    }
    ForceFieldModal TextArea {
        height: 1fr;
        margin-bottom: 1;
    }
    ForceFieldModal Horizontal {
        height: auto;
        align-horizontal: center;
    }
    ForceFieldModal Button {
        width: auto;
        min-width: 10;
        margin: 0 1;
    }
    """

    def __init__(self, molecule: "Molecule") -> None:
        super().__init__()
        self._molecule = molecule

    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label("Select Force Field", classes="title")
            with RadioSet(id="ff-radio"):
                yield RadioButton("None (skip)", value=True, id="ff-none")
                yield RadioButton("OPLS-AA / UFF", id="ff-oplsaa")
                yield RadioButton("PHAHST", id="ff-phahst")
            yield Label("Parameter preview:")
            yield TextArea("", read_only=True, id="ff-preview")
            with Horizontal():
                yield Button("OK", id="ff-apply", variant="primary")
                yield Button("Cancel", id="ff-cancel")

    def on_mount(self) -> None:
        self._update_preview(-1)

    def on_radio_set_changed(self, event: RadioSet.Changed) -> None:
        idx_map = {"ff-none": -1, "ff-oplsaa": 0, "ff-phahst": 1}
        self._update_preview(idx_map.get(event.pressed.id, -1))

    def _update_preview(self, ff_idx: int) -> None:
        from collections import Counter
        elements = Counter(a.element.symbol for a in self._molecule.atoms)

        if ff_idx < 0:
            lines = ["No force field selected.", "", "Parameters will not be written."]
            self.query_one("#ff-preview", TextArea).text = "\n".join(lines)
            return

        ff = get_forcefield(ff_idx)
        label = ["OPLS-AA/UFF", "PHAHST"][ff_idx]
        is_phahst = ff_idx == 1

        lines = [f"{label} -- {len(ff)} element types available", ""]
        if is_phahst:
            lines.append(" El     #      rho     beta    alpha         C6         C8          C10")
            lines.append("               (A)    (1/A)    (A^3)     (a.u.)     (a.u.)       (a.u.)")
            lines.append("-" * 72)
            for el, count in sorted(elements.items()):
                if el in ff:
                    p = ff[el]
                    sig, eps, alp = p.get("sigma", 0), p.get("epsilon", 0), p.get("alpha", 0)
                    c6, c8, c10 = p.get("c6", 0), p.get("c8", 0), p.get("c10", 0)
                    lines.append(
                        f"{el:>3} {count:>5} {sig:>8.4f} {eps:>8.4f} {alp:>8.5f}"
                        f" {c6:>10.4f} {c8:>10.4f} {c10:>12.4f}"
                    )
                else:
                    lines.append(f"{el:>3} {count:>5}      ---      ---      ---        ---        ---          ---")
        else:
            lines.append(" El     #    sigma      eps    alpha")
            lines.append("               (A)      (K)    (A^3)")
            lines.append("-" * 38)
            for el, count in sorted(elements.items()):
                if el in ff:
                    p = ff[el]
                    sig, eps, alp = p.get("sigma", 0), p.get("epsilon", 0), p.get("alpha", 0)
                    lines.append(
                        f"{el:>3} {count:>5} {sig:>8.4f} {eps:>8.2f} {alp:>8.5f}"
                    )
                else:
                    lines.append(f"{el:>3} {count:>5}      ---      ---      ---")

        self.query_one("#ff-preview", TextArea).text = "\n".join(lines)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "ff-apply":
            if self.query_one("#ff-phahst", RadioButton).value:
                self.dismiss("PHAHST")
            elif self.query_one("#ff-oplsaa", RadioButton).value:
                self.dismiss("OPLSAA")
            else:
                self.dismiss("")  # None selected: skip FF but continue the flow
        else:
            self.dismiss(None)  # Cancel: abort the whole flow

    def action_cancel(self) -> None:
        self.dismiss(None)  # Cancel: abort the whole flow


class TrackSlider(Widget):
    """Clickable/draggable horizontal slider for trajectory scrubbing."""

    DEFAULT_CSS = """
    TrackSlider {
        height: 1;
        width: 1fr;
    }
    """

    class Seeked(Message):
        def __init__(self, fraction: float) -> None:
            super().__init__()
            self.fraction = fraction

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._fraction = 0.0
        self._total = 1
        self._current = 0
        self._dragging = False

    def set_position(self, current: int, total: int) -> None:
        self._current = current
        self._total = max(1, total)
        self._fraction = current / self._total
        self.refresh()

    def render_line(self, y: int) -> Strip:
        w = self.size.width
        if w < 1:
            return Strip.blank(w)
        filled = int(self._fraction * w)
        filled = max(0, min(w, filled))
        empty = w - filled
        segments = []
        if filled > 0:
            segments.append(Segment(" " * filled, Style(bgcolor="rgb(0,180,100)")))
        if empty > 0:
            segments.append(Segment(" " * empty, Style(bgcolor="rgb(40,40,45)")))
        return Strip(segments, w)

    def _seek_to_x(self, x: int) -> None:
        w = self.size.width
        if w <= 0:
            return
        frac = max(0.0, min(1.0, x / w))
        self._fraction = frac
        self.refresh()
        self.post_message(self.Seeked(frac))

    def on_mouse_down(self, event) -> None:
        self._dragging = True
        self.capture_mouse()
        self._seek_to_x(event.x)

    def on_mouse_up(self, event) -> None:
        self._dragging = False
        self.release_mouse()

    def on_mouse_move(self, event) -> None:
        if self._dragging:
            self._seek_to_x(event.x)

    def on_click(self, event) -> None:
        self._seek_to_x(event.x)


class SpinBox(Widget, can_focus=True):
    """Integer input with +/- buttons and mouse wheel support."""

    DEFAULT_CSS = """
    SpinBox {
        height: 3;
        layout: horizontal;
    }
    SpinBox Button {
        width: 4;
        min-width: 4;
        height: 3;
    }
    SpinBox .spin-value {
        width: 1fr;
        height: 3;
        content-align: center middle;
        text-style: bold;
        border: tall $accent;
    }
    SpinBox:focus .spin-value {
        border: tall $accent-lighten-2;
    }
    """

    def __init__(self, value: int = 1, min_val: int = 1, max_val: int = 99, step: int = 1, **kwargs) -> None:
        super().__init__(**kwargs)
        self._value = max(min_val, min(max_val, value))
        self._min = min_val
        self._max = max_val
        self._step = step

    @property
    def value(self) -> int:
        return self._value

    def compose(self) -> ComposeResult:
        yield Button("-", id="spin-dec")
        yield Label(str(self._value), id="spin-display", classes="spin-value")
        yield Button("+", id="spin-inc")

    class Changed(Message):
        def __init__(self, spin_box: "SpinBox", value: int) -> None:
            super().__init__()
            self.spin_box = spin_box
            self.value = value

    def _update(self, new: int) -> None:
        new = max(self._min, min(self._max, new))
        if new != self._value:
            self._value = new
            self.query_one("#spin-display", Label).update(str(self._value))
            self.post_message(self.Changed(self, self._value))

    def on_button_pressed(self, event: Button.Pressed) -> None:
        event.stop()
        if event.button.id == "spin-inc":
            self._update(self._value + self._step)
        elif event.button.id == "spin-dec":
            self._update(self._value - self._step)

    def on_mouse_scroll_up(self, event) -> None:
        self._update(self._value + self._step)

    def on_mouse_scroll_down(self, event) -> None:
        self._update(self._value - self._step)

    def on_key(self, event: Key) -> None:
        if event.key == "up":
            self._update(self._value + self._step)
            event.stop()
        elif event.key == "down":
            self._update(self._value - self._step)
            event.stop()


class ExtendAxisModal(ModalScreen[tuple]):
    """Modal to extend all three axes at once. Returns (na, nb, nc) multipliers."""

    BINDINGS = [Binding("escape", "cancel", "Cancel"), Binding("q", "cancel", "Close")]
    DEFAULT_CSS = """
    ExtendAxisModal {
        align: center middle;
    }
    ExtendAxisModal > Vertical {
        width: 60%;
        min-width: 40;
        max-width: 60;
        height: auto;
        max-height: 80%;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    ExtendAxisModal Label {
        margin-bottom: 0;
    }
    ExtendAxisModal .title {
        text-style: bold;
        margin-bottom: 1;
    }
    ExtendAxisModal SpinBox {
        width: 24;
        height: 3;
        margin: 0 1;
    }
    ExtendAxisModal .axis-label {
        width: 4;
        height: 3;
        content-align: right middle;
        text-style: bold;
    }
    ExtendAxisModal .axis-row {
        height: auto;
        margin-bottom: 0;
        align-horizontal: center;
    }
    ExtendAxisModal Horizontal {
        height: auto;
        align-horizontal: center;
    }
    ExtendAxisModal Button {
        width: auto;
        min-width: 10;
        margin: 0 1;
    }
    """

    def __init__(self, pbc) -> None:
        super().__init__()
        self._pbc = pbc

    def compose(self) -> ComposeResult:
        pbc = self._pbc
        if pbc:
            sub = f"a={pbc.a:.2f}  b={pbc.b:.2f}  c={pbc.c:.2f} A"
        else:
            sub = ""
        with Vertical():
            yield Label("Extend Supercell", classes="title")
            if sub:
                yield Label(sub)
            yield Label("Multipliers (final cell = original x multiplier):")
            with Horizontal(classes="axis-row"):
                yield Label("a:", classes="axis-label")
                yield SpinBox(value=1, min_val=1, max_val=20, id="extend-na")
            with Horizontal(classes="axis-row"):
                yield Label("b:", classes="axis-label")
                yield SpinBox(value=1, min_val=1, max_val=20, id="extend-nb")
            with Horizontal(classes="axis-row"):
                yield Label("c:", classes="axis-label")
                yield SpinBox(value=1, min_val=1, max_val=20, id="extend-nc")
            with Horizontal():
                yield Button("Extend", id="extend-ok", variant="primary")
                yield Button("Cancel", id="extend-cancel")

    def on_mount(self) -> None:
        self.query_one("#extend-na", SpinBox).focus()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "extend-ok":
            self._submit()
        elif event.button.id == "extend-cancel":
            self.dismiss(())

    def _submit(self) -> None:
        na = self.query_one("#extend-na", SpinBox).value
        nb = self.query_one("#extend-nb", SpinBox).value
        nc = self.query_one("#extend-nc", SpinBox).value
        self.dismiss((na, nb, nc))

    def action_cancel(self) -> None:
        self.dismiss(())


class SorbateModal(ModalScreen[Optional[str]]):
    """Modal to select a sorbate molecule to insert into the MPMC PDB."""

    BINDINGS = [Binding("escape", "cancel", "Cancel"), Binding("q", "cancel", "Close")]
    DEFAULT_CSS = """
    SorbateModal {
        align: center middle;
    }
    SorbateModal > Vertical {
        width: 70%;
        min-width: 45;
        max-width: 90;
        height: 70%;
        min-height: 16;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    SorbateModal Label {
        margin-bottom: 1;
    }
    SorbateModal .title {
        text-style: bold;
    }
    SorbateModal DataTable {
        height: 1fr;
        margin-bottom: 1;
    }
    SorbateModal Horizontal {
        height: auto;
        align-horizontal: center;
    }
    SorbateModal Button {
        width: auto;
        min-width: 10;
        margin: 0 1;
    }
    """

    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label("Insert Sorbate Molecule", classes="title")
            yield Label("Select a sorbate to place at the center of the unit cell:")
            yield DataTable(id="sorbate-table", cursor_type="row")
            with Horizontal():
                yield Button("Insert", id="sorb-ok", variant="primary")
                yield Button("Skip", id="sorb-skip")
                yield Button("Cancel", id="sorb-cancel")

    def on_mount(self) -> None:
        table = self.query_one("#sorbate-table", DataTable)
        table.add_columns("Model", "Mol", "Sites", "Description")
        groups = [
            ("lj", "--- LJ (nonpolar) ---"),
            ("lj_polar", "--- LJ (polarizable) ---"),
            ("phahst", "--- PHAHST (exp-6-8-10) ---"),
        ]
        for pot_type, header in groups:
            models = {k: v for k, v in SORBATE_MODELS.items() if v.pot_type == pot_type}
            if models:
                table.add_row("", "", "", header, key=f"_hdr_{pot_type}")
                for name, model in models.items():
                    table.add_row(name, model.label, str(len(model.sites)), model.description, key=name)
        table.focus()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "sorb-ok":
            table = self.query_one("#sorbate-table", DataTable)
            row_keys = list(table.rows.keys())
            if row_keys:
                row = max(0, min(table.cursor_row, len(row_keys) - 1))
                key = row_keys[row].value
                if key and not key.startswith("_hdr_"):
                    self.dismiss(key)
                    return
            self.dismiss("")
        elif event.button.id == "sorb-skip":
            self.dismiss("")  # Skip: no sorbate, but continue the flow
        else:
            self.dismiss(None)  # Cancel: abort the whole flow

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        if event.row_key.value and not event.row_key.value.startswith("_hdr_"):
            self.dismiss(event.row_key.value)

    def action_cancel(self) -> None:
        self.dismiss(None)  # Cancel: abort the whole flow


class ChargesFileModal(ModalScreen[Optional[dict]]):
    """File browser with preview and line-skip controls for loading charges."""

    BINDINGS = [Binding("escape", "cancel", "Cancel"), Binding("q", "cancel", "Close")]
    DEFAULT_CSS = """
    ChargesFileModal {
        align: center middle;
    }
    ChargesFileModal > Vertical {
        width: 90%;
        min-width: 50;
        max-width: 140;
        height: 85%;
        min-height: 16;
        border: thick $accent;
        background: $surface;
        padding: 1 2;
    }
    ChargesFileModal Label {
        margin-bottom: 1;
    }
    ChargesFileModal .title {
        text-style: bold;
    }
    ChargesFileModal .hint {
        color: $text-muted;
    }
    ChargesFileModal #browser-and-preview {
        height: 1fr;
    }
    ChargesFileModal DirectoryTree {
        width: 1fr;
        height: 100%;
        border-right: solid $accent;
    }
    ChargesFileModal #preview-panel {
        width: 1fr;
        height: 100%;
        padding: 0 1;
    }
    ChargesFileModal TextArea {
        height: 1fr;
    }
    ChargesFileModal #skip-controls {
        height: auto;
        margin-top: 1;
    }
    ChargesFileModal #skip-controls SpinBox {
        width: 16;
        margin: 0 1;
    }
    ChargesFileModal #skip-controls Label {
        width: auto;
        margin: 0;
    }
    ChargesFileModal #bottom-buttons {
        height: auto;
        align-horizontal: center;
        margin-top: 1;
    }
    ChargesFileModal Button {
        width: auto;
        min-width: 10;
        margin: 0 1;
    }
    """

    def __init__(self) -> None:
        super().__init__()
        self._file_path = ""
        self._raw_lines: list[str] = []

    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label("Select Charges File", classes="title")
            yield Label(
                "Select a .resp file or a raw column of charges. "
                "Use skip controls to trim header/footer lines.",
                classes="hint",
            )
            with Horizontal(id="browser-and-preview"):
                yield DirectoryTree(".", id="charges-tree")
                with Vertical(id="preview-panel"):
                    yield Label("Preview (lines to be used):", classes="hint")
                    yield TextArea("", read_only=True, id="charges-preview")
            with Horizontal(id="skip-controls"):
                yield Label("Skip first:")
                yield SpinBox(value=0, min_val=0, max_val=9999, id="skip-first")
                yield Label("  Skip last:")
                yield SpinBox(value=0, min_val=0, max_val=9999, id="skip-last")
                yield Label("", id="line-count")
            with Horizontal(id="bottom-buttons"):
                yield Button("OK", id="charges-ok", variant="primary")
                yield Button("Cancel", id="charges-cancel")

    def on_directory_tree_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        self._file_path = str(event.path)
        try:
            self._raw_lines = open(self._file_path).readlines()
        except OSError:
            self._raw_lines = []
        self._update_preview()

    def on_spin_box_changed(self, event: SpinBox.Changed) -> None:
        self._update_preview()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "charges-ok":
            self._submit()
        elif event.button.id == "charges-cancel":
            self.dismiss({})

    def _update_preview(self) -> None:
        # CP2K .resp files are auto-parsed (headers/footer/Total stripped) — the
        # skip controls don't apply, so preview the extracted charges directly.
        if self._file_path.lower().endswith(".resp"):
            charges = parse_resp_charges(self._raw_lines)
            preview = [".resp file — auto-parsed, skip controls ignored", ""]
            preview += [f"{i + 1:>4}: {c:+.6f}" for i, c in enumerate(charges)]
            self.query_one("#charges-preview", TextArea).text = "\n".join(preview)
            self.query_one("#line-count", Label).update(
                f"  {len(charges)} charges (.resp)"
            )
            return

        try:
            skip_first = self.query_one("#skip-first", SpinBox).value
        except Exception:
            skip_first = 0
        try:
            skip_last = self.query_one("#skip-last", SpinBox).value
        except Exception:
            skip_last = 0

        lines = self._raw_lines
        total = len(lines)
        end = max(0, total - skip_last)
        start = min(skip_first, end)
        used = lines[start:end]

        preview_lines = []
        for i, line in enumerate(self._raw_lines):
            prefix = "  " if start <= i < end else "X "
            preview_lines.append(f"{prefix}{i + 1:>4}: {line.rstrip()}")

        self.query_one("#charges-preview", TextArea).text = "\n".join(preview_lines)
        self.query_one("#line-count", Label).update(
            f"  {len(used)}/{total} lines"
        )

    def _submit(self) -> None:
        if not self._file_path:
            return
        try:
            skip_first = self.query_one("#skip-first", SpinBox).value
        except Exception:
            skip_first = 0
        try:
            skip_last = self.query_one("#skip-last", SpinBox).value
        except Exception:
            skip_last = 0
        self.dismiss({
            "path": self._file_path,
            "skip_first": skip_first,
            "skip_last": skip_last,
        })

    def action_cancel(self) -> None:
        self.dismiss({})


class PlotWidget(Widget):
    """Braille-based line plot widget."""

    DEFAULT_CSS = """
    PlotWidget {
        height: 1fr;
        border: solid $accent;
    }
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._x: np.ndarray = np.array([])
        self._y: np.ndarray = np.array([])
        self._x_min = 0.0
        self._x_max = 1.0
        self._y_min = 0.0
        self._y_max = 5.0
        self._title = ""
        self._x_label = ""
        self._overlay_x: np.ndarray | None = None
        self._overlay_y: np.ndarray | None = None
        self._overlay_color = "rgb(255,100,80)"

    def set_data(self, x: np.ndarray, y: np.ndarray,
                 title: str = "", x_label: str = "",
                 y_min: float = 0.0, y_max: float = 5.0) -> None:
        self._x = x
        self._y = y
        self._title = title
        self._x_label = x_label
        self._y_min = y_min
        self._y_max = y_max
        if len(x) > 0:
            self._x_min = float(x.min())
            self._x_max = float(x.max())
        self.refresh()

    def set_overlay(self, x: np.ndarray | None, y: np.ndarray | None,
                    color: str = "rgb(255,100,80)") -> None:
        """Set a second data series to overlay (e.g. experimental PXRD)."""
        self._overlay_x = x
        self._overlay_y = y
        self._overlay_color = color
        self.refresh()

    def clear_overlay(self) -> None:
        self._overlay_x = None
        self._overlay_y = None
        self.refresh()

    def zoom_x(self, factor: float) -> None:
        cx = (self._x_min + self._x_max) / 2
        span = (self._x_max - self._x_min) * factor / 2
        self._x_min = max(0, cx - span)
        self._x_max = cx + span
        self.refresh()

    def render_line(self, y_row: int) -> Strip:
        w, h = self.size.width, self.size.height
        if len(self._x) == 0 or w < 10 or h < 4:
            return Strip.blank(w)

        left_margin = 8
        bottom_rows = 2
        plot_cols = w - left_margin
        plot_rows = h - bottom_rows - 1  # terminal rows for the plot area

        if plot_cols < 2 or plot_rows < 2:
            return Strip.blank(w)

        # Pixel dimensions (braille: 2 wide x 4 tall per cell)
        px_w = plot_cols * 2
        px_h = plot_rows * 4

        if y_row == 0:
            title = self._title[:w]
            pad = max(0, (w - len(title)) // 2)
            return Strip([Segment(" " * pad + title + " " * max(0, w - pad - len(title)))], w)

        if y_row > plot_rows:
            if y_row == plot_rows + 1:
                lo = f"{self._x_min:.1f}"
                hi = f"{self._x_max:.1f}"
                mid = f"{(self._x_min + self._x_max) / 2:.1f}"
                line = " " * left_margin + lo
                gap1 = max(1, plot_cols // 2 - len(lo) - len(mid) // 2)
                gap2 = max(1, plot_cols - len(lo) - gap1 - len(mid) - len(hi))
                line += " " * gap1 + mid + " " * gap2 + hi
                return Strip([Segment(line[:w])], w)
            elif y_row == plot_rows + 2:
                label = self._x_label
                pad = max(0, (w - len(label)) // 2)
                return Strip([Segment(" " * pad + label + " " * max(0, w - pad - len(label)))], w)
            return Strip.blank(w)

        # Plot area: render braille line graph for this row
        cell_row = y_row - 1

        x_range = self._x_max - self._x_min
        y_range = self._y_max - self._y_min
        if x_range <= 0:
            x_range = 1.0
        if y_range <= 0:
            y_range = 1.0

        py_top = cell_row * 4
        py_bot = py_top + 3

        def _rasterize(xdata, ydata, cells):
            prev_py = None
            for px in range(px_w):
                x_val = self._x_min + (px + 0.5) / px_w * x_range
                idx = np.searchsorted(xdata, x_val)
                if idx <= 0:
                    y_val = float(ydata[0])
                elif idx >= len(xdata):
                    y_val = float(ydata[-1])
                else:
                    x0, x1 = float(xdata[idx - 1]), float(xdata[idx])
                    y0, y1 = float(ydata[idx - 1]), float(ydata[idx])
                    t = (x_val - x0) / (x1 - x0) if x1 != x0 else 0
                    y_val = y0 + t * (y1 - y0)
                py = int(px_h - 1 - (y_val - self._y_min) / y_range * (px_h - 1))
                py = max(0, min(px_h - 1, py))
                if prev_py is not None:
                    lo_py = min(prev_py, py)
                    hi_py = max(prev_py, py)
                else:
                    lo_py = hi_py = py
                for fill_py in range(lo_py, hi_py + 1):
                    if py_top <= fill_py <= py_bot:
                        col = px // 2
                        sub_x = px % 2
                        sub_y = fill_py - py_top
                        if sub_y < 3:
                            bit = sub_y + sub_x * 3
                        else:
                            bit = 6 + sub_x
                        cells[col] |= (1 << bit)
                prev_py = py

        braille_cells = np.zeros(plot_cols, dtype=np.uint8)
        _rasterize(self._x, self._y, braille_cells)

        overlay_cells = None
        if self._overlay_x is not None and self._overlay_y is not None and len(self._overlay_x) > 0:
            overlay_cells = np.zeros(plot_cols, dtype=np.uint8)
            _rasterize(self._overlay_x, self._overlay_y, overlay_cells)

        chars = [chr(0x2800 + int(b)) for b in braille_cells]

        # Y-axis label
        if cell_row == 0:
            y_lbl = f"{self._y_max:>7.2f} "
        elif cell_row == plot_rows - 1:
            y_lbl = f"{self._y_min:>7.2f} "
        elif cell_row == plot_rows // 2:
            mid_y = (self._y_max + self._y_min) / 2
            y_lbl = f"{mid_y:>7.2f} "
        else:
            y_lbl = " " * 8

        segments = [Segment(y_lbl, Style(color="rgb(150,150,150)"))]
        if overlay_cells is not None:
            # Interleave: where both have dots, show primary; overlay-only in overlay color
            overlay_chars = [chr(0x2800 + int(b)) for b in overlay_cells]
            for ci in range(plot_cols):
                if braille_cells[ci] and overlay_cells[ci]:
                    # Both: show primary over overlay
                    braille_cells[ci] & overlay_cells[ci]
                    braille_cells[ci] & ~overlay_cells[ci]
                    overlay_cells[ci] & ~braille_cells[ci]
                    # Combine: show primary color (overlay underneath is lost in braille)
                    segments.append(Segment(chars[ci], Style(color="rgb(80,200,255)")))
                elif braille_cells[ci]:
                    segments.append(Segment(chars[ci], Style(color="rgb(80,200,255)")))
                elif overlay_cells[ci]:
                    segments.append(Segment(overlay_chars[ci], Style(color=self._overlay_color)))
                else:
                    segments.append(Segment(chr(0x2800)))
        else:
            segments.append(Segment("".join(chars), Style(color="rgb(80,200,255)")))
        return Strip(segments, w)


# ---------------------------------------------------------------------------
# Base class for analysis screens — eliminates duplicate CSS, button handlers,
# CSV export, frame range slicing, and element/name filter updates.
# ---------------------------------------------------------------------------

_ANALYSIS_CSS = """
    {cls} > Vertical {{ width: 100%; height: 100%; background: $surface; padding: 0 1; }}
    {cls} Select {{ width: 16; margin: 0 1; }}
    {cls} SpinBox {{ width: 18; height: 3; margin: 0 1; }}
    {cls} Button {{ width: auto; min-width: 4; margin: 0 1; }}
    {cls} Label {{ width: auto; margin: 0; }}
    {cls} Checkbox {{ width: auto; margin: 0 1; }}
"""


class _AnalysisScreen(ModalScreen[None]):
    """Base class for analysis screens with common patterns."""

    BINDINGS = [Binding("escape", "close", "Close"), Binding("q", "close", "Close")]

    # Subclasses set these
    _prefix = ""  # e.g. "msd", "rdf" — used for widget IDs
    _csv_header = ""
    _csv_columns: list[str] = []  # attribute names for CSV columns

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__()
        self._molecule = molecule
        self._frames = frames
        self._elements = sorted(set(a.element.symbol for a in molecule.atoms))
        self._compute_task: asyncio.Task | None = None

    def _spawn_task(self, coro) -> asyncio.Task:
        """Track the task so action_close can cancel it on dismiss."""
        if self._compute_task is not None and not self._compute_task.done():
            self._compute_task.cancel()
        self._compute_task = asyncio.create_task(coro)
        return self._compute_task

    def _get_frame_slice(self) -> list:
        """Get frame range from Start/End SpinBoxes if present."""
        frames = self._frames
        if frames and len(frames) > 1:
            try:
                start = self.query_one(f"#{self._prefix}-start", SpinBox).value - 1
                end = self.query_one(f"#{self._prefix}-end", SpinBox).value
                frames = frames[start:end]
            except Exception:
                pass
        return frames

    def _yield_frame_range(self, prefix: str = "") -> ComposeResult:
        """Yield Start/End SpinBox widgets for frame range selection."""
        p = prefix or self._prefix
        if self._frames and len(self._frames) > 1:
            n = len(self._frames)
            yield Label("Start:")
            yield SpinBox(value=1, min_val=1, max_val=n, id=f"{p}-start")
            yield Label("End:")
            yield SpinBox(value=n, min_val=1, max_val=n, id=f"{p}-end")

    def _update_name_select(self, element: str, name_select_id: str) -> None:
        """Update a name Select when element changes."""
        if element == "(all)":
            names = sorted(set(a.name.strip() for a in self._molecule.atoms))
        else:
            names = sorted(set(
                a.name.strip() for a in self._molecule.atoms
                if a.element.symbol == element
            ))
        sel = self.query_one(name_select_id, Select)
        sel.set_options([("(all)", "(all)")] + [(n, n) for n in names])
        sel.value = "(all)"

    def _do_export_csv(self, filepath: str, header: str, *arrays) -> None:
        """Write CSV with header and parallel arrays."""
        if not filepath.strip():
            return
        with open(filepath.strip(), "w") as f:
            f.write(header + "\n")
            for row in zip(*arrays):
                f.write(",".join(f"{v:.6f}" for v in row) + "\n")

    def action_close(self) -> None:
        if self._compute_task is not None and not self._compute_task.done():
            self._compute_task.cancel()
            self._compute_task = None
        self.dismiss(None)


class RdfScreen(_AnalysisScreen):
    """Full-screen RDF with element + molecule filtering, Y-axis control, zoom."""

    _prefix = "rdf"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="RdfScreen") + """
    RdfScreen #rdf-row1 { height: auto; margin-bottom: 1; }
    RdfScreen #rdf-row2 { height: auto; margin-bottom: 1; }
    RdfScreen PlotWidget { height: 1fr; }
    RdfScreen #rdf-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._names = sorted(set(a.name.strip() for a in molecule.atoms))
        self._last_r = np.array([])
        self._last_g = np.array([])
        self._last_title = ""

    def _export_csv(self) -> None:
        if len(self._last_r) == 0:
            self.query_one("#rdf-status", Label).update("No data to export.")
            return
        from pathlib import Path
        stem = Path(self._last_title.replace(" ", "_").replace("/", "-"))
        self.app.push_screen(
            FileSaveModal("Save RDF as CSV:", default=f"rdf_{stem}.csv"),
            callback=lambda fp: self._do_export_csv(fp, "r(A),g(r)", self._last_r, self._last_g),
        )

    def _names_for_element(self, el: str) -> list[str]:
        return sorted(set(a.name.strip() for a in self._molecule.atoms if a.element.symbol == el))

    def compose(self) -> ComposeResult:
        if not self._elements:
            self._elements = ["X"]
        el_opts = [(el, el) for el in self._elements]
        el1 = self._elements[0]
        el2 = self._elements[-1]
        names1 = self._names_for_element(el1)
        names2 = self._names_for_element(el2)
        name_opts1 = [("(all)", "(all)")] + [(n, n) for n in names1]
        name_opts2 = [("(all)", "(all)")] + [(n, n) for n in names2]
        with Vertical():
            with Horizontal(id="rdf-row1"):
                yield Label("Type 1:  El:")
                yield Select(el_opts, value=el1, id="rdf-el1")
                yield Label("  Name:")
                yield Select(name_opts1, value="(all)", id="rdf-name1")
                yield Label("    Type 2:  El:")
                yield Select(el_opts, value=el2, id="rdf-el2")
                yield Label("  Name:")
                yield Select(name_opts2, value="(all)", id="rdf-name2")
            with Horizontal(id="rdf-row2"):
                yield Label("Ymax:")
                yield SpinBox(value=5, min_val=1, max_val=100, id="rdf-ymax")
                yield Label("Xmax:")
                yield SpinBox(value=13, min_val=2, max_val=50, id="rdf-xmax")
                yield Label("Bins:")
                yield SpinBox(value=100, min_val=10, max_val=1000, step=10, id="rdf-bins")
                if self._frames and len(self._frames) > 1:
                    n = len(self._frames)
                    yield Checkbox("Traj", value=True, id="rdf-traj")
                    yield Label("Start:")
                    yield SpinBox(value=1, min_val=1, max_val=n, id="rdf-start")
                    yield Label("End:")
                    yield SpinBox(value=n, min_val=1, max_val=n, id="rdf-end")
                yield Button("Compute", id="rdf-compute", variant="primary")
                yield Button("CSV", id="rdf-csv")
                yield Button("Close", id="rdf-close")
            yield PlotWidget(id="rdf-plot")
            yield Label("Select atom types and click Compute  |  +/- zoom x  |  r reset", id="rdf-status")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        bid = event.button.id
        if bid == "rdf-compute":
            self._spawn_task(self._compute_rdf_async())
        elif bid == "rdf-csv":
            self._export_csv()
        elif bid == "rdf-close":
            self.dismiss(None)

    def on_select_changed(self, event: Select.Changed) -> None:
        sid = event.select.id
        if sid == "rdf-el1":
            self._update_name_select(str(event.value), "#rdf-name1")
        elif sid == "rdf-el2":
            self._update_name_select(str(event.value), "#rdf-name2")

    def on_spin_box_changed(self, event: SpinBox.Changed) -> None:
        sid = event.spin_box.id or ""
        plot = self.query_one("#rdf-plot", PlotWidget)
        if sid == "rdf-ymax":
            # Just update the plot Y range, no recompute
            if len(plot._x) > 0:
                plot._y_max = float(event.value)
                plot.refresh()
        elif sid in ("rdf-xmax", "rdf-bins"):
            # Recompute with new params if we have data
            if len(plot._x) > 0:
                self._spawn_task(self._compute_rdf_async())

    def _get_atom_filter(self, el_id: str, name_id: str):
        """Return a function that tests if an atom matches the filter."""
        el = self.query_one(f"#{el_id}", Select).value
        name = self.query_one(f"#{name_id}", Select).value
        if name == "(all)":
            return lambda a: a.element.symbol == el
        else:
            return lambda a: a.element.symbol == el and a.name.strip() == name

    async def _compute_rdf_async(self) -> None:
        filter1 = self._get_atom_filter("rdf-el1", "rdf-name1")
        filter2 = self._get_atom_filter("rdf-el2", "rdf-name2")

        el1 = self.query_one("#rdf-el1", Select).value
        el2 = self.query_one("#rdf-el2", Select).value
        name1 = self.query_one("#rdf-name1", Select).value
        name2 = self.query_one("#rdf-name2", Select).value
        label1 = el1 if name1 == "(all)" else name1
        label2 = el2 if name2 == "(all)" else name2

        y_max = float(self.query_one("#rdf-ymax", SpinBox).value)
        x_max = float(self.query_one("#rdf-xmax", SpinBox).value)
        n_bins = self.query_one("#rdf-bins", SpinBox).value

        use_traj = False
        start_frame = 0
        end_frame = 0
        try:
            use_traj = self.query_one("#rdf-traj", Checkbox).value and self._frames and len(self._frames) > 1
            if use_traj:
                start_frame = self.query_one("#rdf-start", SpinBox).value - 1  # 0-indexed
                end_frame = self.query_one("#rdf-end", SpinBox).value  # exclusive
        except Exception:
            pass

        status = self.query_one("#rdf-status", Label)
        n_total = (end_frame - start_frame) if use_traj else 1

        # Shared progress state polled by a timer
        self._rdf_progress = [0.0]

        def _update_status() -> None:
            frac = self._rdf_progress[0]
            if use_traj:
                frame_num = int(frac * n_total)
                status.update(f"Computing g(r) {label1}-{label2}... frame {frame_num}/{n_total}")
            else:
                status.update(f"Computing g(r) {label1}-{label2}... {int(frac * 100)}%")

        def _prog(frac: float) -> None:
            self._rdf_progress[0] = frac

        poll_timer = self.set_interval(0.1, _update_status)
        status.update(f"Computing g(r) {label1}-{label2}... 0%")
        await asyncio.sleep(0)

        if use_traj:
            frames_slice = self._frames[start_frame:end_frame]
            r, g = await asyncio.to_thread(
                self._rdf_filtered_traj_frames, frames_slice, filter1, filter2, _prog, x_max, n_bins,
            )
            title = f"g(r) {label1}-{label2} (frames {start_frame+1}-{end_frame})"
        else:
            r, g = await asyncio.to_thread(
                self._rdf_filtered, self._molecule, filter1, filter2, _prog, x_max, n_bins,
            )
            title = f"g(r) {label1}-{label2}"

        poll_timer.stop()

        if len(r) == 0:
            status.update(f"No {label1}-{label2} pairs found")
            return

        # Store for CSV export
        self._last_r = r
        self._last_g = g
        self._last_title = title

        plot = self.query_one("#rdf-plot", PlotWidget)
        plot.set_data(r, g, title=title, x_label="r (A)", y_min=0.0, y_max=y_max)
        peak_r = r[g.argmax()]
        status.update(
            f"max g(r) = {g.max():.2f} at r = {peak_r:.2f} A"
        )

    @staticmethod
    def _rdf_filtered(mol, filter1, filter2, progress_callback=None, r_max=None, n_bins=200):
        """Compute RDF with arbitrary atom filters (fully vectorized)."""
        if mol.pbc is None:
            return np.array([]), np.array([])
        pbc = mol.pbc
        if r_max is None:
            r_max = min(pbc.a, pbc.b, pbc.c) / 2.0
        dr = r_max / n_bins

        idx1 = np.array([i for i, a in enumerate(mol.atoms) if filter1(a)])
        idx2 = np.array([i for i, a in enumerate(mol.atoms) if filter2(a)])
        if len(idx1) == 0 or len(idx2) == 0:
            return np.array([]), np.array([])

        coords = np.array([a.x for a in mol.atoms])
        hist = np.zeros(n_bins)
        same = np.array_equal(np.sort(idx1), np.sort(idx2))

        if same:
            # Self-correlation: use upper triangle
            pairs_i, pairs_j = np.triu_indices(len(idx1), k=1)
            all_i = idx1[pairs_i]
            all_j = idx1[pairs_j]
        else:
            # Cross-correlation: all i x j pairs
            grid_i, grid_j = np.meshgrid(idx1, idx2, indexing="ij")
            all_i = grid_i.ravel()
            all_j = grid_j.ravel()

        # Chunk to avoid memory explosion on huge pair counts
        total_pairs = len(all_i)
        chunk_size = max(200000, total_pairs // 20)

        for start in range(0, total_pairs, chunk_size):
            end = min(start + chunk_size, total_pairs)
            dx = coords[all_i[start:end]] - coords[all_j[start:end]]
            frac = dx @ pbc.reciprocal_basis_matrix
            wrapped = dx - np.round(frac) @ pbc.basis_matrix
            dists = np.linalg.norm(wrapped, axis=1)
            valid = dists < r_max
            bins = (dists[valid] / dr).astype(int)
            bins = bins[(bins >= 0) & (bins < n_bins)]
            np.add.at(hist, bins, 1)
            if progress_callback:
                progress_callback(end / total_pairs)

        r_values = (np.arange(n_bins) + 0.5) * dr
        n1, n2 = len(idx1), len(idx2)
        rho = n2 / pbc.volume
        shell_vol = 4.0 * np.pi * r_values**2 * dr
        ideal = rho * shell_vol
        g_r = np.zeros(n_bins)
        nonzero = ideal > 0
        g_r[nonzero] = hist[nonzero] / (n1 * ideal[nonzero])
        return r_values, g_r

    @staticmethod
    def _rdf_filtered_traj_frames(frames, filter1, filter2, progress_callback=None, r_max=None, n_bins=200):
        g_sum = None
        r_values = None
        n = len(frames)
        for fi, mol in enumerate(frames):
            def _sub_prog(frac, _fi=fi):
                if progress_callback:
                    progress_callback((_fi + frac) / n)

            r, g = RdfScreen._rdf_filtered(mol, filter1, filter2, _sub_prog, r_max, n_bins)
            if len(r) == 0:
                continue
            if g_sum is None:
                g_sum = g.copy()
                r_values = r
            else:
                g_sum += g
        if progress_callback:
            progress_callback(1.0)
        if g_sum is None:
            return np.array([]), np.array([])
        return r_values, g_sum / n


class MsdScreen(_AnalysisScreen):
    """MSD analysis with element selector, frame range, and diffusion constant."""
    _prefix = "msd"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="MsdScreen") + """
    MsdScreen #msd-controls { height: auto; margin-bottom: 1; }
    MsdScreen #msd-row2 { height: auto; margin-bottom: 1; }
    MsdScreen #msd-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._last_lags = np.array([])
        self._last_msd = np.array([])

    def compose(self) -> ComposeResult:
        el_opts = [(el, el) for el in self._elements]
        names = sorted(set(a.name.strip() for a in self._molecule.atoms))
        name_opts = [("(all)", "(all)")] + [(n, n) for n in names]
        with Vertical():
            with Horizontal(id="msd-controls"):
                yield Label("Element:")
                yield Select(el_opts, value=self._elements[-1] if self._elements else "H", id="msd-el")
                yield Label("Name:")
                yield Select(name_opts, value="(all)", id="msd-name")
                yield Button("Compute", id="msd-compute", variant="primary")
                yield Button("CSV", id="msd-csv")
                yield Button("Close", id="msd-close")
            with Horizontal(id="msd-row2"):
                yield Label("dt (ps):")
                yield SpinBox(value=1, min_val=1, max_val=10000, id="msd-dt")
                yield from self._yield_frame_range()
            yield PlotWidget(id="msd-plot")
            yield Label("Select atom type and click Compute", id="msd-status")

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id == "msd-el":
            self._update_name_select(str(event.value), "#msd-name")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "msd-compute":
            self._spawn_task(self._run_msd())
        elif event.button.id == "msd-close":
            self.dismiss(None)
        elif event.button.id == "msd-csv":
            if len(self._last_lags) == 0:
                return
            self.app.push_screen(FileSaveModal("Save MSD as CSV:", default="msd.csv"), callback=self._export_csv)

    async def _run_msd(self) -> None:
        el = str(self.query_one("#msd-el", Select).value)
        name = str(self.query_one("#msd-name", Select).value)
        name_filter = name if name != "(all)" else None
        dt = float(self.query_one("#msd-dt", SpinBox).value)
        frames = self._get_frame_slice()

        status = self.query_one("#msd-status", Label)
        status.update(f"Computing MSD for {el}...")
        await asyncio.sleep(0)
        lags, msd_vals = await asyncio.to_thread(
            compute_msd, frames, element=el, name=name_filter,
        )
        if len(lags) == 0:
            status.update(f"No {el} atoms found")
            return
        self._last_lags = lags
        self._last_msd = msd_vals
        plot = self.query_one("#msd-plot", PlotWidget)
        y_max = float(msd_vals.max()) * 1.1 if msd_vals.max() > 0 else 1.0
        label = el if not name_filter else name_filter
        plot.set_data(lags * dt, msd_vals, title=f"MSD {label}", x_label="time (ps)", y_min=0.0, y_max=y_max)
        if len(lags) > 10:
            D, slope = compute_diffusion_constant(lags, msd_vals, dt=dt)
            status.update(f"D = {D:.4e} cm^2/s  (slope = {slope:.4f} A^2/ps)")
        else:
            status.update("Done")

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "lag_frames,msd_A2", self._last_lags, self._last_msd)


class CoordinationScreen(_AnalysisScreen):
    """Coordination number analysis."""
    _prefix = "coord"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="CoordinationScreen") + """
    CoordinationScreen #coord-controls { height: auto; margin-bottom: 1; }
    CoordinationScreen #coord-row2 { height: auto; margin-bottom: 1; }
    CoordinationScreen #coord-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._last_r = np.array([])
        self._last_cn = np.array([])

    def compose(self) -> ComposeResult:
        el_opts = [(el, el) for el in self._elements]
        el1 = self._elements[0] if self._elements else "X"
        el2 = self._elements[-1] if self._elements else "X"
        with Vertical():
            with Horizontal(id="coord-controls"):
                yield Label("Center:")
                yield Select(el_opts, value=el1, id="coord-el1")
                yield Label("Neighbor:")
                yield Select(el_opts, value=el2, id="coord-el2")
                yield Label("Rmax:")
                yield SpinBox(value=5, min_val=2, max_val=20, id="coord-rmax")
                yield Button("Compute", id="coord-compute", variant="primary")
                yield Button("CSV", id="coord-csv")
                yield Button("Close", id="coord-close")
            with Horizontal(id="coord-row2"):
                if self._frames and len(self._frames) > 1:
                    yield Checkbox("Traj", value=False, id="coord-traj")
                    yield from self._yield_frame_range()
            yield PlotWidget(id="coord-plot")
            yield Label("Select center and neighbor elements, then Compute", id="coord-status")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "coord-compute":
            self._spawn_task(self._run_coord())
        elif event.button.id == "coord-close":
            self.dismiss(None)
        elif event.button.id == "coord-csv":
            if len(self._last_r) == 0:
                return
            self.app.push_screen(FileSaveModal("Save coordination as CSV:", default="coordination.csv"), callback=self._export_csv)

    async def _run_coord(self) -> None:
        el1 = str(self.query_one("#coord-el1", Select).value)
        el2 = str(self.query_one("#coord-el2", Select).value)
        r_max = float(self.query_one("#coord-rmax", SpinBox).value)
        status = self.query_one("#coord-status", Label)

        use_traj = False
        try:
            use_traj = self.query_one("#coord-traj", Checkbox).value and self._frames and len(self._frames) > 1
        except Exception:
            pass

        if use_traj:
            frames = self._get_frame_slice()
            status.update(f"Computing N({el2}) around {el1} over {len(frames)} frames...")
            await asyncio.sleep(0)
            r, cn, avg = await asyncio.to_thread(
                compute_coordination_trajectory, frames, el1, el2, r_max,
            )
        else:
            status.update(f"Computing N({el2}) around {el1}...")
            await asyncio.sleep(0)
            r, cn, avg = await asyncio.to_thread(
                compute_coordination, self._molecule, el1, el2, r_max,
            )

        if len(r) == 0:
            status.update(f"No {el1} or {el2} atoms found")
            return

        self._last_r = r
        self._last_cn = cn
        plot = self.query_one("#coord-plot", PlotWidget)
        y_max = float(cn.max()) * 1.1 if cn.max() > 0 else 1.0
        plot.set_data(r, cn, title=f"N({el2}) around {el1}", x_label="r (A)", y_min=0.0, y_max=y_max)
        status.update(f"Avg coordination at {r_max} A: {avg:.2f} {el2} per {el1}")

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "r_A,coordination_number", self._last_r, self._last_cn)


class RmsdScreen(_AnalysisScreen):
    """RMSD vs time analysis."""
    _prefix = "rmsd"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="RmsdScreen") + """
    RmsdScreen #rmsd-controls { height: auto; margin-bottom: 1; }
    RmsdScreen #rmsd-row2 { height: auto; margin-bottom: 1; }
    RmsdScreen #rmsd-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._last_idx = np.array([])
        self._last_rmsd = np.array([])

    def compose(self) -> ComposeResult:
        el_opts = [("(all)", "(all)")] + [(el, el) for el in self._elements]
        names = sorted(set(a.name.strip() for a in self._molecule.atoms))
        name_opts = [("(all)", "(all)")] + [(n, n) for n in names]
        with Vertical():
            with Horizontal(id="rmsd-controls"):
                yield Label("Element:")
                yield Select(el_opts, value="(all)", id="rmsd-el")
                yield Label("Name:")
                yield Select(name_opts, value="(all)", id="rmsd-name")
                yield Button("Compute", id="rmsd-compute", variant="primary")
                yield Button("CSV", id="rmsd-csv")
                yield Button("Close", id="rmsd-close")
            with Horizontal(id="rmsd-row2"):
                yield Label("dt (ps):")
                yield SpinBox(value=1, min_val=1, max_val=10000, id="rmsd-dt")
                yield from self._yield_frame_range()
            yield PlotWidget(id="rmsd-plot")
            yield Label("Select atom filter and click Compute", id="rmsd-status")

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id == "rmsd-el":
            self._update_name_select(str(event.value), "#rmsd-name")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "rmsd-compute":
            self._spawn_task(self._run_rmsd())
        elif event.button.id == "rmsd-close":
            self.dismiss(None)
        elif event.button.id == "rmsd-csv":
            if len(self._last_idx) == 0:
                return
            self.app.push_screen(FileSaveModal("Save RMSD as CSV:", default="rmsd.csv"), callback=self._export_csv)

    async def _run_rmsd(self) -> None:
        el = str(self.query_one("#rmsd-el", Select).value)
        el_filter = el if el != "(all)" else None
        name = str(self.query_one("#rmsd-name", Select).value)
        name_filter = name if name != "(all)" else None
        dt = float(self.query_one("#rmsd-dt", SpinBox).value)
        frames = self._get_frame_slice()

        status = self.query_one("#rmsd-status", Label)
        label = el_filter or "all"
        status.update(f"Computing RMSD for {label}...")
        await asyncio.sleep(0)
        idx, rmsd_vals = await asyncio.to_thread(
            compute_rmsd, frames, element=el_filter, name=name_filter,
        )
        if len(idx) == 0:
            status.update("No matching atoms found")
            return
        self._last_idx = idx
        self._last_rmsd = rmsd_vals
        plot = self.query_one("#rmsd-plot", PlotWidget)
        time = idx * dt
        y_max = float(rmsd_vals.max()) * 1.1 if rmsd_vals.max() > 0 else 1.0
        plot.set_data(time, rmsd_vals, title=f"RMSD ({label})", x_label="time (ps)", y_min=0.0, y_max=y_max)
        avg_rmsd = float(rmsd_vals[len(rmsd_vals)//2:].mean())
        status.update(f"RMSD avg (2nd half) = {avg_rmsd:.3f} A  |  max = {rmsd_vals.max():.3f} A")

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "frame,rmsd_A", self._last_idx, self._last_rmsd)


class HbondScreen(_AnalysisScreen):
    """Hydrogen bond count over trajectory."""
    _prefix = "hb"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="HbondScreen") + """
    HbondScreen #hb-controls { height: auto; margin-bottom: 1; }
    HbondScreen #hb-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._last_idx = np.array([])
        self._last_counts = np.array([])

    def compose(self) -> ComposeResult:
        with Vertical():
            with Horizontal(id="hb-controls"):
                yield Label("D-A cutoff:")
                yield SpinBox(value=35, min_val=20, max_val=50, id="hb-cutoff")
                yield Label("/10 A")
                yield Label("  Angle min:")
                yield SpinBox(value=120, min_val=90, max_val=170, id="hb-angle")
                yield Label("deg")
                yield from self._yield_frame_range()
                yield Button("Compute", id="hb-compute", variant="primary")
                yield Button("CSV", id="hb-csv")
                yield Button("Close", id="hb-close")
            yield PlotWidget(id="hb-plot")
            yield Label("Click Compute to count H-bonds per frame", id="hb-status")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "hb-compute":
            self._spawn_task(self._run())
        elif event.button.id == "hb-close":
            self.dismiss(None)
        elif event.button.id == "hb-csv":
            if len(self._last_idx) == 0:
                return
            self.app.push_screen(FileSaveModal("Save H-bonds as CSV:", default="hbonds.csv"), callback=self._export_csv)

    async def _run(self) -> None:
        cutoff = self.query_one("#hb-cutoff", SpinBox).value / 10.0
        angle = float(self.query_one("#hb-angle", SpinBox).value)
        status = self.query_one("#hb-status", Label)

        if self._frames and len(self._frames) > 1:
            frames = self._get_frame_slice()
            status.update(f"Counting H-bonds over {len(frames)} frames...")
            await asyncio.sleep(0)
            idx, counts = await asyncio.to_thread(count_hbonds_trajectory, frames, cutoff, angle)
        else:
            status.update("Detecting H-bonds...")
            await asyncio.sleep(0)
            hb = await asyncio.to_thread(detect_hbonds, self._molecule, cutoff, angle)
            idx = np.array([0.0])
            counts = np.array([float(len(hb))])

        if len(idx) == 0:
            status.update("No frames")
            return

        self._last_idx = idx
        self._last_counts = counts
        plot = self.query_one("#hb-plot", PlotWidget)
        y_max = float(counts.max()) * 1.1 if counts.max() > 0 else 1.0
        plot.set_data(idx, counts, title="H-bonds vs frame", x_label="frame", y_min=0.0, y_max=y_max)
        avg = float(counts.mean())
        status.update(f"Avg: {avg:.1f} H-bonds/frame  |  max: {counts.max():.0f}  |  min: {counts.min():.0f}")

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "frame,hbond_count", self._last_idx, self._last_counts)


class GyrationScreen(_AnalysisScreen):
    """Radius of gyration vs time."""
    _prefix = "rg"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="GyrationScreen") + """
    GyrationScreen #rg-controls { height: auto; margin-bottom: 1; }
    GyrationScreen #rg-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._last_idx = np.array([])
        self._last_rg = np.array([])

    def compose(self) -> ComposeResult:
        el_opts = [("(all)", "(all)")] + [(el, el) for el in self._elements]
        with Vertical():
            with Horizontal(id="rg-controls"):
                yield Label("Element:")
                yield Select(el_opts, value="(all)", id="rg-el")
                yield Label("dt (ps):")
                yield SpinBox(value=1, min_val=1, max_val=10000, id="rg-dt")
                yield from self._yield_frame_range()
                yield Button("Compute", id="rg-compute", variant="primary")
                yield Button("CSV", id="rg-csv")
                yield Button("Close", id="rg-close")
            yield PlotWidget(id="rg-plot")
            yield Label("Select filter and click Compute", id="rg-status")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "rg-compute":
            self._spawn_task(self._run())
        elif event.button.id == "rg-close":
            self.dismiss(None)
        elif event.button.id == "rg-csv":
            if len(self._last_idx) == 0:
                return
            self.app.push_screen(FileSaveModal("Save Rg as CSV:", default="rg.csv"), callback=self._export_csv)

    async def _run(self) -> None:
        el = str(self.query_one("#rg-el", Select).value)
        el_filter = el if el != "(all)" else None
        dt = float(self.query_one("#rg-dt", SpinBox).value)
        frames = self._get_frame_slice()

        status = self.query_one("#rg-status", Label)
        label = el_filter or "all"
        status.update(f"Computing Rg for {label}...")
        await asyncio.sleep(0)
        idx, rg = await asyncio.to_thread(compute_gyration_trajectory, frames, el_filter)
        if len(idx) == 0:
            status.update("No data")
            return
        self._last_idx = idx
        self._last_rg = rg
        plot = self.query_one("#rg-plot", PlotWidget)
        time = idx * dt
        y_max = float(rg.max()) * 1.1 if rg.max() > 0 else 1.0
        y_min = float(rg.min()) * 0.9 if rg.min() > 0 else 0.0
        plot.set_data(time, rg, title=f"Rg ({label})", x_label="time (ps)", y_min=y_min, y_max=y_max)
        avg = float(rg.mean())
        status.update(f"Avg Rg = {avg:.3f} A  |  std = {rg.std():.3f} A")

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "frame,rg_A", self._last_idx, self._last_rg)


class ReduceCellScreen(ModalScreen[Optional[Molecule]]):
    """Reduce a supercell to a primitive cell, selecting which sorbates to keep."""
    BINDINGS = [Binding("escape", "close", "Close"), Binding("q", "close", "Close")]
    DEFAULT_CSS = """
    ReduceCellScreen > Vertical { width: 100%; height: 100%; background: $surface; padding: 1 2; }
    ReduceCellScreen #rc-controls { height: auto; margin-bottom: 1; }
    ReduceCellScreen SpinBox { width: 10; height: 3; margin: 0 1; }
    ReduceCellScreen Button { width: auto; min-width: 6; margin: 0 1; }
    ReduceCellScreen Label { width: auto; margin: 0; }
    ReduceCellScreen DataTable { height: 1fr; }
    ReduceCellScreen #rc-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule) -> None:
        super().__init__()
        self._molecule = molecule
        self._framework_idx: list[int] = []
        self._sorbate_mols: list[list[int]] = []
        self._keep: list[bool] = []

    def compose(self) -> ComposeResult:
        with Vertical():
            with Horizontal(id="rc-controls"):
                yield Label("Supercell:")
                yield SpinBox(value=2, min_val=1, max_val=8, id="rc-na")
                yield Label("x")
                yield SpinBox(value=2, min_val=1, max_val=8, id="rc-nb")
                yield Label("x")
                yield SpinBox(value=2, min_val=1, max_val=8, id="rc-nc")
                yield Button("Analyze", id="rc-analyze", variant="primary")
                yield Button("Reduce", id="rc-reduce")
                yield Button("Close", id="rc-close")
            yield DataTable(id="rc-table")
            yield Label("Set supercell dimensions and click Analyze", id="rc-status")

    def on_mount(self) -> None:
        table = self.query_one("#rc-table", DataTable)
        table.add_columns("Keep", "#", "Atoms", "Elements", "Center (A)")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "rc-analyze":
            self._do_analyze()
        elif event.button.id == "rc-reduce":
            self._do_reduce()
        elif event.button.id == "rc-close":
            self.dismiss(None)

    def _do_analyze(self) -> None:
        na = self.query_one("#rc-na", SpinBox).value
        nb = self.query_one("#rc-nb", SpinBox).value
        nc = self.query_one("#rc-nc", SpinBox).value
        status = self.query_one("#rc-status", Label)

        fw, sorbs = identify_framework_and_sorbates(self._molecule, na, nb, nc)
        self._framework_idx = fw
        self._sorbate_mols = sorbs
        self._keep = [True] * len(sorbs)

        table = self.query_one("#rc-table", DataTable)
        table.clear()

        n_mult = na * nb * nc
        status_msg = f"Framework: {len(fw)} unique sites ({len(fw) * n_mult} total atoms)"
        status_msg += f"  |  Sorbate molecules: {len(sorbs)}"

        for i, mol_indices in enumerate(sorbs):
            elements = {}
            positions = []
            for idx in mol_indices:
                a = self._molecule.atoms[idx]
                elements[a.element.symbol] = elements.get(a.element.symbol, 0) + 1
                positions.append(a.x)
            elem_str = " ".join(f"{el}{n}" for el, n in sorted(elements.items()))
            center = np.mean(positions, axis=0) if positions else np.zeros(3)
            center_str = f"{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}"
            table.add_row("[X]", str(i + 1), str(len(mol_indices)), elem_str, center_str)

        status.update(status_msg + "  |  Click rows to toggle keep/discard, then Reduce")

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        """Toggle keep/discard on row click."""
        table = self.query_one("#rc-table", DataTable)
        row_idx = event.cursor_row
        if 0 <= row_idx < len(self._keep):
            self._keep[row_idx] = not self._keep[row_idx]
            mark = "[X]" if self._keep[row_idx] else "[ ]"
            # Update first cell
            row_key = table.get_row_at(row_idx)
            table.update_cell(row_key, table.columns[list(table.columns.keys())[0]].key, mark)

    def _do_reduce(self) -> None:
        if not self._framework_idx and not self._sorbate_mols:
            self.app.notify("Click Analyze first", timeout=2)
            return
        na = self.query_one("#rc-na", SpinBox).value
        nb = self.query_one("#rc-nb", SpinBox).value
        nc = self.query_one("#rc-nc", SpinBox).value
        keep_indices = [i for i, k in enumerate(self._keep) if k]
        result = reduce_to_primitive(self._molecule, na, nb, nc, keep_indices)
        n_kept = len(keep_indices)
        n_total = len(self._sorbate_mols)
        self.app.notify(
            f"Reduced to primitive cell: {len(result.atoms)} atoms, "
            f"kept {n_kept}/{n_total} sorbate molecules",
            timeout=3,
        )
        self.dismiss(result)

    def action_close(self) -> None:
        self.dismiss(None)


class InputGeneratorScreen(ModalScreen[None]):
    """Unified template-driven input generator.

    Renders an input file (and any companion files) from a template under
    `pdb_wizard.templates`. Engines:
      - mpmc_isotherm  : multi-pressure μVT MC tree
      - cp2k_cell_opt  : CP2K cell optimization
      - cp2k_dft_md    : CP2K Born-Oppenheimer NVT MD
      - openmm_npt     : OpenMM NPT MD driver (Python)
      - openmm_nvt     : OpenMM NVT MD driver (Python)

    The screen reads its form definition from `templates.ENGINES`. Switching
    the engine rebuilds the form. Edit the shipped template files in place
    to tweak the rendered output beyond what the form exposes.
    """

    BINDINGS = [Binding("escape", "close", "Close"), Binding("q", "close", "Close")]
    DEFAULT_CSS = """
    InputGeneratorScreen > Vertical { width: 100%; height: 100%; background: $surface; padding: 1 2; }
    InputGeneratorScreen #gen-row-engine { height: auto; margin-bottom: 1; }
    InputGeneratorScreen #gen-row-outdir { height: auto; margin-bottom: 1; }
    InputGeneratorScreen #gen-row-buttons { height: auto; margin-top: 1; margin-bottom: 1; }
    InputGeneratorScreen #gen-form { height: auto; max-height: 16; margin-bottom: 1; }
    InputGeneratorScreen .gen-form-row { height: auto; margin-bottom: 1; }
    InputGeneratorScreen Select { width: 24; margin: 0 1; }
    InputGeneratorScreen #gen-engine { width: 24; }
    InputGeneratorScreen SpinBox { width: 18; height: 3; margin: 0 1; }
    InputGeneratorScreen Input { width: 1fr; margin: 0 1; }
    InputGeneratorScreen Button { width: auto; min-width: 6; margin: 0 1; }
    InputGeneratorScreen Label { width: auto; margin: 0; }
    InputGeneratorScreen .gen-var-label { width: 28; }
    InputGeneratorScreen TextArea { height: 1fr; }
    InputGeneratorScreen #gen-status { height: 1; color: $text-muted; }
    InputGeneratorScreen #gen-description { height: 1; color: $text-muted; }
    """

    def __init__(
        self,
        molecule: Molecule,
        filepath: str = "",
        engine: str = "mpmc_isotherm",
    ) -> None:
        super().__init__()
        if engine not in ENGINES:
            raise ValueError(f"unknown engine key: {engine!r}")
        self._molecule = molecule
        self._filepath = filepath
        self._engine_key = engine
        # widget refs for value extraction, populated by _build_form
        self._form_widgets: dict[str, Widget] = {}
        # cached current values per engine so switching back restores them
        self._values_by_engine: dict[str, dict] = {}

    # ---------------------------------------------------------------
    # Compose
    # ---------------------------------------------------------------
    def compose(self) -> ComposeResult:
        engine_opts = [(spec.label, key) for key, spec in ENGINES.items()]
        with Vertical():
            with Horizontal(id="gen-row-engine"):
                yield Label("Engine:")
                yield Select(engine_opts, value=self._engine_key, id="gen-engine")
                yield Label("", id="gen-description")
                yield Button("Close", id="gen-close")
            with Horizontal(id="gen-row-outdir"):
                yield Label("Output dir:")
                yield Input(value=self._default_outdir(), id="gen-outdir")
            yield Vertical(id="gen-form")
            with Horizontal(id="gen-row-buttons"):
                yield Button("Validate", id="gen-validate")
                yield Button("Show template", id="gen-show-template")
                yield Button("Generate", id="gen-generate", variant="primary")
            yield TextArea(id="gen-preview", read_only=True)
            yield Label("Pick an engine, set values, then Generate.", id="gen-status")

    def _default_outdir(self) -> str:
        """Folder name next to the loaded file (or cwd) for the current engine."""
        base = Path(self._filepath).parent if self._filepath else Path(".")
        folder_for = {
            "mpmc_isotherm": "isotherm",
            "cp2k_cell_opt": "cp2k_cell_opt",
            "cp2k_dft_md": "cp2k_md",
            "openmm_npt": "openmm_npt",
            "openmm_nvt": "openmm_nvt",
        }
        return str(base / folder_for.get(self._engine_key, "sim_input"))

    # ---------------------------------------------------------------
    # Mount + form building
    # ---------------------------------------------------------------
    def on_mount(self) -> None:
        self._build_form()
        self._update_description()
        # If MPMC, run validation immediately so the user sees errors
        if self._engine_key == "mpmc_isotherm":
            self._show_validation()
        else:
            self._set_status("Edit values, then Generate.")
            self._set_preview(
                f"Engine: {self._engine_key}\n"
                f"Output file: {self._engine_spec().output_filename}\n\n"
                f"Click Show template to inspect the source template, "
                f"or Generate to write the rendered files."
            )

    def _engine_spec(self):
        return ENGINES[self._engine_key]

    def _build_form(self) -> None:
        """(Re)build the dynamic form area for the current engine."""
        form = self.query_one("#gen-form", Vertical)
        # Tear down existing children
        for child in list(form.children):
            child.remove()
        self._form_widgets.clear()
        spec = self._engine_spec()
        # Persisted values (if user switched away and back)
        persisted = self._values_by_engine.get(self._engine_key, {})

        for var in spec.user_vars:
            widget_id = f"gen-{var.name}"
            value = persisted.get(var.name, var.default)
            row = Horizontal(classes="gen-form-row")
            form.mount(row)
            row.mount(Label(var.label + ":", classes="gen-var-label"))
            widget: Widget
            if var.kind == "select":
                choices = var.choices or ()
                widget = Select(
                    [(c, c) for c in choices],
                    value=value if value in choices else (choices[0] if choices else ""),
                    id=widget_id,
                )
            elif var.kind == "int":
                widget = SpinBox(
                    value=int(value),
                    min_val=int(var.min_val) if var.min_val is not None else 0,
                    max_val=int(var.max_val) if var.max_val is not None else 99,
                    step=int(var.step) if var.step is not None else 1,
                    id=widget_id,
                )
            else:
                # 'float', 'str', 'csv_floats' all map to Input
                widget = Input(value=str(value), id=widget_id)
            row.mount(widget)
            self._form_widgets[var.name] = widget

    def _update_description(self) -> None:
        try:
            self.query_one("#gen-description", Label).update(
                self._engine_spec().description
            )
        except Exception:
            pass

    # ---------------------------------------------------------------
    # Event handling
    # ---------------------------------------------------------------
    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id != "gen-engine":
            return
        # Snapshot the outgoing engine's values before swapping
        self._snapshot_values()
        new_key = str(event.value)
        if new_key == self._engine_key:
            return
        self._engine_key = new_key
        # Reset outdir to engine default
        self.query_one("#gen-outdir", Input).value = self._default_outdir()
        self._build_form()
        self._update_description()
        if new_key == "mpmc_isotherm":
            self._show_validation()
        else:
            self._set_status("Edit values, then Generate.")
            self._set_preview(
                f"Engine: {new_key}\nOutput: {self._engine_spec().output_filename}"
            )

    def _snapshot_values(self) -> None:
        try:
            self._values_by_engine[self._engine_key] = self._collect_values(strict=False)
        except Exception:
            pass

    def on_button_pressed(self, event: Button.Pressed) -> None:
        bid = event.button.id
        if bid == "gen-close":
            self.dismiss(None)
        elif bid == "gen-validate":
            self._show_validation()
        elif bid == "gen-show-template":
            self._show_template()
        elif bid == "gen-generate":
            self._do_generate()

    def action_close(self) -> None:
        self.dismiss(None)

    # ---------------------------------------------------------------
    # Value collection
    # ---------------------------------------------------------------
    def _collect_values(self, *, strict: bool = True) -> dict:
        """Read each form widget back to its declared type."""
        spec = self._engine_spec()
        out: dict = {}
        for var in spec.user_vars:
            w = self._form_widgets.get(var.name)
            if w is None:
                if strict:
                    raise KeyError(f"missing form widget for {var.name!r}")
                continue
            if var.kind == "int":
                out[var.name] = int(w.value)  # SpinBox.value
            elif var.kind == "float":
                raw = str(w.value).strip()
                if not raw:
                    if strict:
                        raise ValueError(f"{var.label} is empty")
                    continue
                out[var.name] = float(raw)
            elif var.kind == "csv_floats":
                raw = str(w.value)
                parts = [p.strip() for p in raw.split(",") if p.strip()]
                if not parts:
                    if strict:
                        raise ValueError(f"{var.label} is empty")
                    continue
                out[var.name] = [float(p) for p in parts]
            elif var.kind == "select":
                out[var.name] = str(w.value)
            else:  # 'str'
                out[var.name] = str(w.value)
        return out

    # ---------------------------------------------------------------
    # Status / preview helpers
    # ---------------------------------------------------------------
    def _set_status(self, text: str) -> None:
        try:
            self.query_one("#gen-status", Label).update(text)
        except Exception:
            pass

    def _set_preview(self, text: str) -> None:
        try:
            self.query_one("#gen-preview", TextArea).load_text(text)
        except Exception:
            pass

    # ---------------------------------------------------------------
    # MPMC validation
    # ---------------------------------------------------------------
    def _validate_mpmc(self) -> tuple[list[str], list[str]]:
        return validate_for_mpmc(self._molecule)

    def _format_validation(self, errors: list[str], warnings: list[str]) -> str:
        lines: list[str] = []
        if not errors and not warnings:
            lines.append("✓ All checks passed")
            lines.append("")
            lines.append("Charges, LJ params, box size, and sorbate placement all look fine.")
            return "\n".join(lines)
        if errors:
            lines.append(f"✗ {len(errors)} error{'s' if len(errors) != 1 else ''} "
                         f"(must fix before generating):")
            lines.append("")
            for e in errors:
                lines.append(f"  ✗ {e}")
            lines.append("")
        if warnings:
            lines.append(f"⚠ {len(warnings)} warning"
                         f"{'s' if len(warnings) != 1 else ''} "
                         f"(generation will proceed):")
            lines.append("")
            for w in warnings:
                lines.append(f"  ⚠ {w}")
        return "\n".join(lines)

    def _show_validation(self) -> None:
        if self._engine_key == "mpmc_isotherm":
            errors, warnings = self._validate_mpmc()
            self._set_preview(self._format_validation(errors, warnings))
            if errors:
                self._set_status(
                    f"{len(errors)} error{'s' if len(errors) != 1 else ''} — "
                    f"fix before Generate"
                )
            elif warnings:
                self._set_status(
                    f"{len(warnings)} warning{'s' if len(warnings) != 1 else ''} — "
                    f"review then Generate"
                )
            else:
                self._set_status("Ready to generate")
        else:
            # Light-weight check: PBC required for CP2K, NPT
            problems: list[str] = []
            if self._engine_key.startswith("cp2k") and self._molecule.pbc is None:
                problems.append("CP2K requires a periodic box. Add one in Update Unit Cell.")
            if self._engine_key == "openmm_npt" and self._molecule.pbc is None:
                problems.append("OpenMM NPT requires a periodic box.")
            if problems:
                self._set_preview("✗ " + "\n✗ ".join(problems))
                self._set_status("Fix issues before Generate")
            else:
                self._set_preview("✓ No engine-level problems detected.")
                self._set_status("Ready to generate")

    # ---------------------------------------------------------------
    # Show template
    # ---------------------------------------------------------------
    def _show_template(self) -> None:
        try:
            tpl = load_template(self._engine_spec().template_file)
            self._set_preview(
                f"-- {self._engine_spec().template_file} --\n\n{tpl}"
            )
            self._set_status(
                f"Template: {self._engine_spec().template_file} "
                f"(edit the file in src/pdb_wizard/templates/ to change defaults)"
            )
        except Exception as e:
            self._set_status(f"Error loading template: {e}")

    # ---------------------------------------------------------------
    # Generate
    # ---------------------------------------------------------------
    def _do_generate(self) -> None:
        try:
            values = self._collect_values(strict=True)
        except (ValueError, KeyError) as e:
            self._set_status(f"Invalid input: {e}")
            return

        outdir = self.query_one("#gen-outdir", Input).value.strip()
        if not outdir:
            self._set_status("Output dir is required")
            return

        engine = self._engine_key
        try:
            if engine == "mpmc_isotherm":
                self._generate_mpmc(values, outdir)
            elif engine == "cp2k_cell_opt":
                self._generate_cp2k(values, outdir, mode="cell_opt")
            elif engine == "cp2k_dft_md":
                self._generate_cp2k(values, outdir, mode="dft_md")
            elif engine == "openmm_npt":
                self._generate_openmm(values, outdir, ensemble="NPT")
            elif engine == "openmm_nvt":
                self._generate_openmm(values, outdir, ensemble="NVT")
            else:
                self._set_status(f"Unknown engine: {engine}")
        except Exception as e:
            self._set_status(f"Error: {e}")

    def _generate_mpmc(self, values: dict, outdir: str) -> None:
        # Refuse on validation errors (warnings are allowed)
        errors, warnings = self._validate_mpmc()
        if errors:
            self._set_preview(self._format_validation(errors, warnings))
            self._set_status(
                f"Refusing: {len(errors)} error{'s' if len(errors) != 1 else ''} "
                f"would produce a meaningless run"
            )
            self.app.notify(
                "Cannot generate — fix errors shown in the preview",
                severity="error", timeout=5,
            )
            return

        dirs = generate_isotherm(
            self._molecule,
            sorbate_model=values["sorbate"],
            temperature=values["temperature"],
            pressures=values["pressures"],
            output_dir=outdir,
            numsteps=values["numsteps"],
            corrtime=values["corrtime"],
            ensemble=values["ensemble"],
            insert_probability=values["insert_probability"],
            move_factor=values["move_factor"],
            rot_factor=values["rot_factor"],
        )

        lines: list[str] = []
        if warnings:
            lines.append(f"⚠ Proceeded with {len(warnings)} warning"
                         f"{'s' if len(warnings) != 1 else ''}:")
            for w in warnings:
                lines.append(f"  ⚠ {w}")
            lines.append("")
        lines.extend([
            f"Generated {len(dirs)} pressure points:",
            f"  Sorbate: {values['sorbate']}",
            f"  Temperature: {values['temperature']} K",
            f"  Steps: {values['numsteps']}",
            "",
            "Directories:",
        ])
        for d in dirs:
            lines.append(f"  {d}/")
            lines.append("    input.pqr, mpmc.inp, insert.pqr, run.sh")
        lines.extend([
            "",
            f"Master script: {outdir}/run_all.sh",
            "",
            f"To run: cd {outdir} && bash run_all.sh",
        ])
        self._set_preview("\n".join(lines))
        tail = f" ({len(warnings)} warnings)" if warnings else ""
        self._set_status(f"Generated {len(dirs)} directories in {outdir}/{tail}")

    def _generate_cp2k(self, values: dict, outdir: str, mode: str) -> None:
        if self._molecule.pbc is None:
            self._set_status("CP2K requires a periodic box")
            return
        kwargs = dict(
            project_name=values["project"],
            cutoff_Ry=values["cutoff_ry"],
            rel_cutoff_Ry=values["rel_cutoff_ry"],
            xc_functional=values["xc_functional"],
            eps_scf=values["eps_scf"],
            max_scf=values["max_scf"],
        )
        if mode == "dft_md":
            kwargs.update(
                temperature_K=values["temperature"],
                timestep_fs=values["timestep_fs"],
                n_md_steps=values["n_md_steps"],
                md_traj_interval=values["md_traj_interval"],
            )
        else:
            kwargs["n_opt_steps"] = values["n_opt_steps"]
        result = generate_cp2k_input(self._molecule, outdir, mode=mode, **kwargs)
        self._set_preview(
            f"Wrote CP2K {mode} input:\n"
            f"  {result['input']}\n"
            f"  {result['xyz']}\n\n"
            f"Elements: {', '.join(result['elements'])}\n\n"
            f"To run: cd {outdir} && cp2k -i {Path(result['input']).name} -o out.log"
        )
        self._set_status(f"Wrote {result['input']}")

    def _generate_openmm(self, values: dict, outdir: str, ensemble: str) -> None:
        kwargs = dict(
            forcefield=values["forcefield"],
            temperature_K=values["temperature_k"],
            timestep_fs=values["timestep_fs"],
            n_steps=values["n_steps"],
            report_interval=values["report_interval"],
            nonbonded_cutoff_nm=values["nonbonded_cutoff_nm"],
        )
        if ensemble == "NPT":
            kwargs["pressure_atm"] = values["pressure_atm"]
        result = generate_openmm_script(
            self._molecule, outdir, ensemble=ensemble, **kwargs,
        )
        self._set_preview(
            f"Wrote OpenMM {ensemble} driver:\n"
            f"  {result['script']}\n"
            f"  {result['pdb']}\n\n"
            f"To run: cd {outdir} && python {Path(result['script']).name}"
        )
        self._set_status(f"Wrote {result['script']}")


# Backward-compat: existing imports of IsothermPlannerScreen keep working.
# Subclass so the engine arg defaults to MPMC isotherm without callers needing
# to pass it.
class IsothermPlannerScreen(InputGeneratorScreen):
    def __init__(self, molecule: Molecule, filepath: str = "") -> None:
        super().__init__(molecule, filepath, engine="mpmc_isotherm")


class DatabaseSearchScreen(ModalScreen[None]):
    """Unified database search.

    One screen for COD, RCSB PDB, PubChem, and Materials Project.
    Per-source search modes (text / formula / direct ID), source-aware
    placeholder hints, per-source result counts, and a preview line for
    the currently-selected row before fetching."""
    BINDINGS = [Binding("escape", "close", "Close"), Binding("q", "close", "Close")]
    DEFAULT_CSS = """
    DatabaseSearchScreen > Vertical { width: 100%; height: 100%; background: $surface; padding: 1 2; }
    DatabaseSearchScreen #db-row1 { height: auto; margin-bottom: 1; }
    DatabaseSearchScreen #db-row2 { height: auto; margin-bottom: 1; }
    DatabaseSearchScreen Input { width: 1fr; margin: 0 1; }
    DatabaseSearchScreen #db-source { width: 22; margin: 0 1; }
    DatabaseSearchScreen #db-mode { width: 14; margin: 0 1; }
    DatabaseSearchScreen Button { width: auto; min-width: 6; margin: 0 1; }
    DatabaseSearchScreen Label { width: auto; margin: 0; }
    DatabaseSearchScreen DataTable { height: 1fr; }
    DatabaseSearchScreen #db-preview { height: 1; color: $text-muted; margin-top: 1; }
    DatabaseSearchScreen #db-status { height: 1; color: $text-muted; }
    """

    # Modes the "All"-source aggregate accepts. We dispatch the user's query
    # to every backend that supports the chosen mode.
    _ALL_MODES = ("text", "formula")

    def __init__(self) -> None:
        super().__init__()
        self._results: list[dict] = []
        self._compute_task: asyncio.Task | None = None

    def _spawn_task(self, coro) -> asyncio.Task:
        if self._compute_task is not None and not self._compute_task.done():
            self._compute_task.cancel()
        self._compute_task = asyncio.create_task(coro)
        return self._compute_task

    def _modes_for(self, source: str) -> tuple[str, ...]:
        if source == "all":
            return self._ALL_MODES
        return BACKENDS.get(source, {}).get("modes", ("text",))

    def _placeholder_for(self, source: str, mode: str) -> str:
        if source == "all":
            if mode == "formula":
                return "formula across all DBs (eg 'C9H8O4')"
            return "free-text search across all DBs"
        spec = BACKENDS.get(source, {})
        return spec.get("placeholder", {}).get(mode, "query")

    def compose(self) -> ComposeResult:
        db_opts = [("All databases", "all")] + [
            (spec["label"], code) for code, spec in BACKENDS.items()
        ]
        # Default to "all" + "text" — the broadest, most useful start.
        default_source = "all"
        default_mode = "text"
        mode_opts = [(m.capitalize(), m) for m in self._modes_for(default_source)]
        with Vertical():
            with Horizontal(id="db-row1"):
                yield Label("Database:")
                yield Select(db_opts, value=default_source, id="db-source")
                yield Label("Search by:")
                yield Select(mode_opts, value=default_mode, id="db-mode")
                yield Button("Close", id="db-close")
            with Horizontal(id="db-row2"):
                yield Input(
                    placeholder=self._placeholder_for(default_source, default_mode),
                    id="db-query",
                )
                yield Button("Search", id="db-search", variant="primary")
                yield Button("Fetch by ID", id="db-fetch-id")
                yield Button("Fetch selected", id="db-fetch")
            yield DataTable(id="db-results")
            yield Label("", id="db-preview")
            yield Label(
                "Pick a source + mode → type a query → Search. "
                "Or 'Fetch by ID' to skip search.",
                id="db-status",
            )

    def on_mount(self) -> None:
        table = self.query_one("#db-results", DataTable)
        table.add_columns("Source", "ID", "Formula", "Title")
        table.cursor_type = "row"
        # Auto-update the preview when the user moves the cursor
        table.show_cursor = True

    def on_select_changed(self, event: Select.Changed) -> None:
        sid = event.select.id
        if sid == "db-source":
            # Refresh the mode dropdown to match the newly-picked source
            new_source = str(event.value)
            new_modes = self._modes_for(new_source)
            mode_select = self.query_one("#db-mode", Select)
            current_mode = str(mode_select.value)
            mode_select.set_options([(m.capitalize(), m) for m in new_modes])
            # Preserve mode if still valid; otherwise pick first
            mode_select.value = current_mode if current_mode in new_modes else new_modes[0]
            self._update_placeholder()
        elif sid == "db-mode":
            self._update_placeholder()

    def _update_placeholder(self) -> None:
        try:
            source = str(self.query_one("#db-source", Select).value)
            mode = str(self.query_one("#db-mode", Select).value)
            self.query_one("#db-query", Input).placeholder = self._placeholder_for(source, mode)
        except Exception:
            pass

    def on_data_table_row_highlighted(self, event) -> None:
        # Update preview pane
        idx = event.cursor_row
        if idx is None or idx < 0 or idx >= len(self._results):
            return
        r = self._results[idx]
        src = BACKENDS.get(r.get("source", ""), {}).get("label", r.get("source", ""))
        formula = r.get("formula", "")
        title = r.get("title", "")
        bits = [f"[bold]{src}[/bold] {r.get('id', '')}"]
        if formula:
            bits.append(f"formula: {formula}")
        if title:
            bits.append(f"title: {title}")
        try:
            self.query_one("#db-preview", Label).update("  •  ".join(bits))
        except Exception:
            pass

    def on_input_submitted(self, event: Input.Submitted) -> None:
        if event.input.id == "db-query":
            self._spawn_task(self._do_search())

    def on_button_pressed(self, event: Button.Pressed) -> None:
        bid = event.button.id
        if bid == "db-search":
            self._spawn_task(self._do_search())
        elif bid == "db-close":
            self.dismiss(None)
        elif bid == "db-fetch":
            self._fetch_selected()
        elif bid == "db-fetch-id":
            self._fetch_by_id()

    async def _do_search(self) -> None:

        query = self.query_one("#db-query", Input).value.strip()
        source = str(self.query_one("#db-source", Select).value)
        mode = str(self.query_one("#db-mode", Select).value)

        # Direct-ID mode skips search → pop the user straight into a fetch
        if mode == "id":
            if not query:
                self.query_one("#db-status", Label).update(
                    "Enter an ID then click 'Fetch by ID'"
                )
                return
            self._spawn_task(self._do_fetch(source, query))
            return

        if not query:
            return

        status = self.query_one("#db-status", Label)
        if source == "all":
            # Only query backends that support this mode
            codes = [c for c, spec in BACKENDS.items() if mode in spec["modes"]]
        else:
            codes = [source]

        status.update(
            f"Searching {', '.join(BACKENDS[c]['label'] for c in codes)} "
            f"({mode})..."
        )
        await asyncio.sleep(0)

        coros = [asyncio.to_thread(BACKENDS[c]["search"], query, mode) for c in codes]
        per_backend = await asyncio.gather(*coros, return_exceptions=True)

        results: list[dict] = []
        per_source_counts: dict[str, int] = {}
        errors: list[str] = []
        for code, out in zip(codes, per_backend):
            if isinstance(out, Exception):
                errors.append(f"{BACKENDS[code]['label']} ({type(out).__name__})")
                per_source_counts[code] = 0
                continue
            results.extend(out)
            per_source_counts[code] = len(out)

        # Sort results by source label so they group nicely in the table
        results.sort(key=lambda r: (r.get("source", ""), r.get("id", "")))

        self._results = results
        table = self.query_one("#db-results", DataTable)
        table.clear()
        for r in results:
            label = BACKENDS.get(r.get("source", ""), {}).get("label", r.get("source", ""))
            table.add_row(
                label,
                r.get("id", ""),
                r.get("formula", "")[:25],
                r.get("title", "")[:60],
            )

        if source == "all" and len(codes) > 1:
            breakdown = ", ".join(
                f"{BACKENDS[c]['label']}={per_source_counts.get(c, 0)}" for c in codes
            )
            msg = f"{len(results)} results  ({breakdown})"
        else:
            msg = f"{len(results)} results in {BACKENDS[codes[0]]['label']}" if codes else "0 results"
        if errors:
            msg += f"  •  errors: {', '.join(errors)}"
        status.update(msg)
        self.query_one("#db-preview", Label).update("")

    def _fetch_selected(self) -> None:
        table = self.query_one("#db-results", DataTable)
        if table.cursor_row is None or table.cursor_row >= len(self._results):
            self.app.notify("Select a result first", timeout=2)
            return
        result = self._results[table.cursor_row]
        source = result.get("source", "")
        rid = result.get("id", "")
        self._spawn_task(self._do_fetch(source, rid))

    def _fetch_by_id(self) -> None:
        """Direct-ID fetch path (skip search)."""
        source = str(self.query_one("#db-source", Select).value)
        if source == "all":
            self.app.notify("Pick a specific database first (not 'All')", timeout=3)
            return
        rid = self.query_one("#db-query", Input).value.strip()
        if not rid:
            self.app.notify("Enter an ID in the query field", timeout=2)
            return
        self._spawn_task(self._do_fetch(source, rid))

    async def _do_fetch(self, source: str, rid: str) -> None:

        spec = BACKENDS.get(source)
        if spec is None:
            self.query_one("#db-status", Label).update(f"Unknown source: {source!r}")
            return

        label = spec["label"]
        status = self.query_one("#db-status", Label)
        status.update(f"Fetching {label} {rid}...")
        await asyncio.sleep(0)

        try:
            cache_dir = self.app._db_cache_dir()
            path = await asyncio.to_thread(spec["fetch"], rid, cache_dir)
            # Parsing a large CIF is the synchronous slow part — offload it.
            mol = await asyncio.to_thread(read_file, path)
            self.app.open_in_new_tab(mol, path)
            self.app.notify(
                f"Opened {label} {rid} ({len(mol.atoms)} atoms)", timeout=3,
            )
            self.dismiss(None)  # auto-close on success
        except Exception as e:
            status.update(f"Error: {e}")
            self.app.notify(f"Fetch failed: {e}", timeout=4, severity="error")

    def action_close(self) -> None:
        if self._compute_task is not None and not self._compute_task.done():
            self._compute_task.cancel()
            self._compute_task = None
        self.dismiss(None)


class EnergyPlotScreen(_AnalysisScreen):
    """Energy plot screen for MPMC energy.dat files."""
    _prefix = "energy"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="EnergyPlotScreen") + """
    EnergyPlotScreen Select { width: 20; }
    EnergyPlotScreen #energy-controls { height: auto; margin-bottom: 1; }
    EnergyPlotScreen #energy-status { height: 1; color: $text-muted; }
    """

    def __init__(self, data: dict, filepath: str) -> None:
        # Skip _AnalysisScreen.__init__ — no molecule for energy plots
        ModalScreen.__init__(self)
        self._data = data
        self._filepath = filepath
        self._molecule = None
        self._frames = None
        self._elements = []

    def _n_steps(self) -> int:
        """Number of rows in the energy data."""
        for k in ("step", "N", "energy"):
            if k in self._data:
                return len(self._data[k])
        if self._data:
            return len(next(iter(self._data.values())))
        return 0

    def compose(self) -> ComposeResult:
        cols = [k for k in self._data if k != "step"]
        options = [(c, c) for c in cols]
        # Pick a sensible default: N (molecule count) or energy
        default = cols[0] if cols else "N"
        for prefer in ("N", "energy"):
            if prefer in cols:
                default = prefer
                break
        n = self._n_steps()
        with Vertical():
            with Horizontal(id="energy-controls"):
                yield Label("Column:")
                yield Select(options, value=default, id="energy-col")
                if n > 1:
                    yield Label("Start:")
                    yield SpinBox(value=1, min_val=1, max_val=n, id="energy-start")
                    yield Label("End:")
                    yield SpinBox(value=n, min_val=1, max_val=n, id="energy-end")
                yield Button("CSV", id="energy-csv")
                yield Button("Close", id="energy-close")
            yield PlotWidget(id="energy-plot")
            yield Label(f"Loaded {self._filepath}", id="energy-status")

    def on_mount(self) -> None:
        col = str(self.query_one("#energy-col", Select).value)
        self._plot_column(col)

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id == "energy-col" and event.value is not None:
            self._plot_column(str(event.value))

    def on_spin_box_changed(self, event) -> None:
        # Re-plot whenever Start or End changes.
        if event.spin_box.id in ("energy-start", "energy-end"):
            try:
                col = str(self.query_one("#energy-col", Select).value)
            except Exception:
                return
            self._plot_column(col)

    def _step_slice(self) -> slice:
        """Resolve the user's Start/End SpinBoxes into a Python slice.
        Returns the full range if the SpinBoxes don't exist (no data)."""
        n = self._n_steps()
        try:
            start = self.query_one("#energy-start", SpinBox).value - 1  # 1-indexed in UI
            end = self.query_one("#energy-end", SpinBox).value         # inclusive in UI
        except Exception:
            return slice(0, n)
        # Clamp + ensure end > start
        start = max(0, min(start, n - 1))
        end = max(start + 1, min(end, n))
        return slice(start, end)

    def _plot_column(self, col: str) -> None:
        if col not in self._data:
            return
        sl = self._step_slice()
        full_steps = self._data.get("step", np.arange(len(self._data[col])))
        steps = full_steps[sl]
        values = self._data[col][sl]
        if len(values) == 0:
            self.query_one("#energy-status", Label).update("Empty step range")
            return
        plot = self.query_one("#energy-plot", PlotWidget)
        y_min = float(values.min()) * (0.9 if values.min() > 0 else 1.1)
        y_max = float(values.max()) * 1.1 if values.max() > 0 else float(values.max()) * 0.9
        if y_min == y_max:
            y_max = y_min + 1
        title = col
        n_total = self._n_steps()
        if sl.start != 0 or sl.stop != n_total:
            title = f"{col} (steps {sl.start + 1}-{sl.stop})"
        plot.set_data(steps, values, title=title, x_label="step",
                      y_min=y_min, y_max=y_max)
        avg = float(values[-len(values)//4:].mean()) if len(values) > 4 else float(values.mean())
        self.query_one("#energy-status", Label).update(
            f"{col} [{sl.start + 1}-{sl.stop}]: "
            f"last-quarter avg = {avg:.4f}, final = {values[-1]:.4f}"
        )

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "energy-close":
            self.dismiss(None)
        elif event.button.id == "energy-csv":
            self.app.push_screen(
                FileSaveModal("Save energy as CSV:", default="energy.csv"),
                callback=self._export_csv,
            )

    def _export_csv(self, filepath: str) -> None:
        if not filepath.strip():
            return
        sl = self._step_slice()
        with open(filepath.strip(), "w") as f:
            cols = list(self._data.keys())
            f.write(",".join(cols) + "\n")
            for i in range(sl.start, sl.stop):
                row = ",".join(f"{self._data[c][i]:.6f}" for c in cols)
                f.write(row + "\n")
        self.query_one("#energy-status", Label).update(
            f"Exported steps {sl.start + 1}-{sl.stop} to {filepath.strip()}"
        )


class DensityScreen(_AnalysisScreen):
    """Density profile with element/axis selectors."""
    _prefix = "dens"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="DensityScreen") + """
    DensityScreen #dens-controls { height: auto; margin-bottom: 1; }
    DensityScreen #dens-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._last_pos = np.array([])
        self._last_dens = np.array([])

    def compose(self) -> ComposeResult:
        el_opts = [(el, el) for el in self._elements]
        names = sorted(set(a.name.strip() for a in self._molecule.atoms))
        name_opts = [("(all)", "(all)")] + [(n, n) for n in names]
        axis_opts = [("a (x)", "0"), ("b (y)", "1"), ("c (z)", "2")]
        with Vertical():
            with Horizontal(id="dens-controls"):
                yield Label("Element:")
                yield Select(el_opts, value=self._elements[0] if self._elements else "H", id="dens-el")
                yield Label("Name:")
                yield Select(name_opts, value="(all)", id="dens-name")
                yield Label("Axis:")
                yield Select(axis_opts, value="2", id="dens-axis")
                if self._frames and len(self._frames) > 1:
                    yield Checkbox("Traj", value=True, id="dens-traj")
                yield from self._yield_frame_range()
                yield Button("Compute", id="dens-compute", variant="primary")
                yield Button("CSV", id="dens-csv")
                yield Button("Close", id="dens-close")
            yield PlotWidget(id="dens-plot")
            yield Label("Select atom type, axis, and click Compute", id="dens-status")

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id == "dens-el":
            self._update_name_select(str(event.value), "#dens-name")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "dens-compute":
            self._spawn_task(self._run_density())
        elif event.button.id == "dens-close":
            self.dismiss(None)
        elif event.button.id == "dens-csv":
            if len(self._last_pos) == 0:
                return
            self.app.push_screen(FileSaveModal("Save density as CSV:", default="density.csv"), callback=self._export_csv)

    async def _run_density(self) -> None:
        el = str(self.query_one("#dens-el", Select).value)
        name = str(self.query_one("#dens-name", Select).value)
        axis = int(self.query_one("#dens-axis", Select).value)
        name_filter = name if name != "(all)" else None
        use_traj = False
        try:
            use_traj = self.query_one("#dens-traj", Checkbox).value and self._frames and len(self._frames) > 1
        except Exception:
            pass

        # Frame range
        if use_traj:
            frames = self._get_frame_slice()
        else:
            frames = self._frames

        status = self.query_one("#dens-status", Label)
        axis_label = ["a", "b", "c"][axis]
        status.update(f"Computing density profile {el} along {axis_label}...")
        await asyncio.sleep(0)

        if use_traj:
            pos, dens = await asyncio.to_thread(
                compute_density_profile_trajectory, frames, axis, 100, el, name_filter,
            )
            title = f"Density {el} ({axis_label}-axis, {len(frames)} frames)"
        else:
            pos, dens = await asyncio.to_thread(
                compute_density_profile, self._molecule, axis, 100, el, name_filter,
            )
            title = f"Density {el} ({axis_label}-axis)"

        if len(pos) == 0:
            status.update(f"No {el} atoms found")
            return

        self._last_pos = pos
        self._last_dens = dens
        plot = self.query_one("#dens-plot", PlotWidget)
        y_max = float(dens.max()) * 1.1 if dens.max() > 0 else 1.0
        plot.set_data(pos, dens, title=title, x_label=f"position along {axis_label} (A)", y_min=0.0, y_max=y_max)
        status.update(f"max density = {dens.max():.2f} atoms/A at {pos[dens.argmax()]:.2f} A")

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "position_A,density", self._last_pos, self._last_dens)


class PxrdScreen(_AnalysisScreen):
    """Simulated Powder XRD pattern."""
    _prefix = "pxrd"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="PxrdScreen") + """
    PxrdScreen Select { width: 14; }
    PxrdScreen #pxrd-controls { height: auto; margin-bottom: 1; }
    PxrdScreen #pxrd-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule) -> None:
        super().__init__(molecule, None)
        self._last_tt = np.array([])
        self._last_intensity = np.array([])

    def compose(self) -> ComposeResult:
        wl_opts = [("Cu Ka 1.5406", "1.5406"), ("Mo Ka 0.7107", "0.7107"), ("Co Ka 1.7903", "1.7903")]
        with Vertical():
            with Horizontal(id="pxrd-controls"):
                yield Label("Source:")
                yield Select(wl_opts, value="1.5406", id="pxrd-wl")
                yield Label("2T max:")
                yield SpinBox(value=50, min_val=10, max_val=90, id="pxrd-ttmax")
                yield Label("Width:")
                yield SpinBox(value=10, min_val=1, max_val=50, id="pxrd-width")
                yield Button("Compute", id="pxrd-compute", variant="primary")
                yield Button("Load Exp.", id="pxrd-load-exp")
                yield Button("CSV", id="pxrd-csv")
                yield Button("Close", id="pxrd-close")
            yield PlotWidget(id="pxrd-plot")
            yield Label("Select parameters and click Compute  |  Load Exp. to overlay experimental data", id="pxrd-status")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "pxrd-compute":
            self._spawn_task(self._run_pxrd())
        elif event.button.id == "pxrd-close":
            self.dismiss(None)
        elif event.button.id == "pxrd-csv":
            if len(self._last_tt) == 0:
                return
            self.app.push_screen(FileSaveModal("Save PXRD as CSV:", default="pxrd.csv"), callback=self._export_csv)
        elif event.button.id == "pxrd-load-exp":
            self.app.push_screen(
                FileSaveModal("Load experimental PXRD (CSV/XY):", default="", button_label="Open"),
                callback=self._load_experimental,
            )

    async def _run_pxrd(self) -> None:
        wl = float(self.query_one("#pxrd-wl", Select).value)
        tt_max = float(self.query_one("#pxrd-ttmax", SpinBox).value)
        width = self.query_one("#pxrd-width", SpinBox).value / 100.0  # hundredths to degrees

        status = self.query_one("#pxrd-status", Label)
        status.update("Computing PXRD pattern...")

        self._rdf_progress = [0.0]
        def _prog(frac):
            self._rdf_progress[0] = frac
        poll = self.set_interval(0.1, lambda: status.update(
            f"Computing PXRD... {int(self._rdf_progress[0] * 100)}%"
        ))
        await asyncio.sleep(0)

        tt, intensity = await asyncio.to_thread(
            compute_pxrd, self._molecule, wavelength=wl,
            two_theta_max=tt_max, peak_width=width, progress_callback=_prog,
        )

        poll.stop()

        if len(tt) == 0:
            status.update("No reflections found")
            return

        self._last_tt = tt
        self._last_intensity = intensity
        plot = self.query_one("#pxrd-plot", PlotWidget)
        wl_name = "Cu Ka" if abs(wl - 1.5406) < 0.01 else f"{wl:.4f} A"
        plot.set_data(tt, intensity, title=f"PXRD ({wl_name})",
                      x_label="2-theta (deg)", y_min=0.0, y_max=105.0)
        n_peaks = len([r for r in intensity if r > 5])
        status.update(f"{n_peaks} peaks above 5% intensity")

    def _load_experimental(self, filepath: str) -> None:
        if not filepath.strip():
            return
        try:
            data = np.loadtxt(filepath.strip(), delimiter=None, comments="#")
            if data.ndim == 1:
                # Try comma-separated
                data = np.loadtxt(filepath.strip(), delimiter=",", comments="#", skiprows=1)
            if data.ndim != 2 or data.shape[1] < 2:
                self.query_one("#pxrd-status", Label).update("Error: need 2 columns (2theta, intensity)")
                return
            exp_tt = data[:, 0]
            exp_int = data[:, 1]
            # Normalize to 100
            if exp_int.max() > 0:
                exp_int = exp_int / exp_int.max() * 100
            plot = self.query_one("#pxrd-plot", PlotWidget)
            plot.set_overlay(exp_tt, exp_int, color="rgb(255,130,60)")
            self.query_one("#pxrd-status", Label).update(
                f"Experimental data loaded ({len(exp_tt)} points) — cyan=simulated, orange=experimental"
            )
        except Exception as e:
            self.query_one("#pxrd-status", Label).update(f"Error loading: {e}")

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "two_theta_deg,intensity", self._last_tt, self._last_intensity)


class MoleculeView(Widget):
    """Braille-based 3D molecule renderer."""

    can_focus = True

    def __init__(self) -> None:
        super().__init__()
        self.molecule: Molecule | None = None
        self.pbc: PBC | None = None
        self.rot_matrix = rotation_matrix(-0.2, -0.5, 0.0)
        self.camera_distance = 4.0
        self.show_bonds = True
        self.dark_bg = True
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.pan_mode = False
        self.highlighted_atoms: set[int] = set()
        self.show_atom_numbers = False
        self.licorice = False
        self.vdw = False
        self.color_by_charge = False
        self.hide_water = False
        self.ribbon = False
        self.color_mode = "element"
        self.atom_scale = 0.35
        self.bond_radius = 0.08
        self.ambient = 0.50
        self.diffuse = 0.60
        self.specular = 0.40
        self.shininess = 32.0
        self._cached_strips: list[Strip] = []
        self._cached_size: tuple[int, int] = (0, 0)
        self._dragging = False
        self._drag_start_x = 0
        self._drag_start_y = 0
        self._locked_centroid = None
        self._blink_on = True
        self._blink_timer = None
        self._picked_atoms: list[int] = []  # for measurement mode
        self.measure_mode = False
        self._measurements: list[tuple[list[int], str]] = []  # [(atom_indices, label_text)]
        self.show_polyhedra = False
        self.orthographic = False
        self._density_positions: np.ndarray | None = None
        self._density_values: np.ndarray | None = None
        self._density_element: str = ""

    class AtomPicked(Message):
        def __init__(self, atom_index: int, shift: bool) -> None:
            super().__init__()
            self.atom_index = atom_index
            self.shift = shift

    def _pick_atom_at(self, term_x: int, term_y: int) -> int | None:
        """Find the atom index closest to terminal cell (term_x, term_y)."""
        if self.molecule is None or len(self.molecule.atoms) == 0:
            return None
        w, h = self.size.width, self.size.height
        if w == 0 or h == 0:
            return None
        px_w = w * 2
        px_h = h * 4
        # Convert terminal cell to pixel center
        click_px = term_x * 2 + 1
        click_py = term_y * 4 + 2
        fov = 1.5
        scale = min(px_w, px_h) / 2
        centroid = self._locked_centroid if self._locked_centroid is not None else self.molecule.center()
        rot = self.rot_matrix
        best_idx = None
        best_dist_sq = float("inf")
        for i, atom in enumerate(self.molecule.atoms):
            pos = rot @ (atom.position - centroid)
            pos[0] += self.pan_x
            pos[1] += self.pan_y
            pos[2] += self.camera_distance
            if pos[2] <= 0.1:
                continue
            sx = px_w / 2 + pos[0] * fov / pos[2] * scale
            sy = px_h / 2 - pos[1] * fov / pos[2] * scale
            d2 = (sx - click_px) ** 2 + (sy - click_py) ** 2
            if d2 < best_dist_sq:
                best_dist_sq = d2
                best_idx = i
        # Only pick if within reasonable pixel distance
        if best_dist_sq > (30 * 30):
            return None
        return best_idx

    def on_mouse_down(self, event) -> None:
        self._dragging = True
        self._drag_start_x = event.x
        self._drag_start_y = event.y
        self._drag_moved = False
        self.capture_mouse()

    def on_mouse_up(self, event) -> None:
        was_drag = self._drag_moved
        self._dragging = False
        self.release_mouse()
        if not was_drag and self.measure_mode:
            # Only pick atoms in measure mode
            idx = self._pick_atom_at(event.x, event.y)
            if idx is not None:
                self.post_message(self.AtomPicked(idx, event.shift))

    def on_mouse_move(self, event) -> None:
        if not self._dragging:
            return
        dx = event.x - self._drag_start_x
        dy = event.y - self._drag_start_y
        if abs(dx) > 1 or abs(dy) > 1:
            self._drag_moved = True
        self._drag_start_x = event.x
        self._drag_start_y = event.y
        if self.pan_mode:
            self.pan_x += dx * self.camera_distance * 0.01
            self.pan_y -= dy * self.camera_distance * 0.01
            self._clamp_pan()
        else:
            self.rot_matrix = rotation_matrix(-dy * 0.02, -dx * 0.02, 0) @ self.rot_matrix
        self._invalidate_cache()

    def on_mouse_scroll_up(self, event) -> None:
        self.camera_distance = max(1.0, self.camera_distance - 0.5)
        self._invalidate_cache()

    def on_mouse_scroll_down(self, event) -> None:
        self.camera_distance += 0.5
        self._invalidate_cache()

    def set_molecule(self, molecule: Molecule, keep_camera: bool = False) -> None:
        self.molecule = molecule
        self.pbc = molecule.pbc
        if not keep_camera:
            mol_radius = molecule.radius()
            self.camera_distance = max(4.0, mol_radius * 3.0)
            self._locked_centroid = None
        self._invalidate_cache()

    def _clamp_pan(self) -> None:
        if self.molecule is None:
            return
        max_pan = self.molecule.radius() * 0.5
        self.pan_x = max(-max_pan, min(max_pan, self.pan_x))
        self.pan_y = max(-max_pan, min(max_pan, self.pan_y))

    def _invalidate_cache(self) -> None:
        self._cached_size = (0, 0)
        self.refresh()

    def _start_blink(self) -> None:
        if self._blink_timer is None:
            self._blink_on = True
            self._blink_timer = self.set_interval(0.5, self._toggle_blink)

    def _stop_blink(self) -> None:
        if self._blink_timer is not None:
            self._blink_timer.stop()
            self._blink_timer = None
            self._blink_on = True
            self._invalidate_cache()

    def _toggle_blink(self) -> None:
        if not self.highlighted_atoms:
            self._stop_blink()
            return
        self._blink_on = not self._blink_on
        self._invalidate_cache()

    def render_line(self, y: int) -> Strip:
        w, h = self.size.width, self.size.height
        if (w, h) != self._cached_size:
            self._rebuild(w, h)
        if 0 <= y < len(self._cached_strips):
            return self._cached_strips[y]
        return Strip.blank(w)

    def _rebuild(self, cols: int, rows: int) -> None:
        self._cached_size = (cols, rows)

        if self.molecule is None or cols == 0 or rows == 0:
            self._cached_strips = [Strip.blank(cols) for _ in range(rows)]
            return

        px_w = cols * 2
        px_h = rows * 4

        bg = (0, 0, 0) if self.dark_bg else (255, 255, 255)
        rot = self.rot_matrix

        mol = self.molecule
        if self.hide_water:
            # Detect actual water molecules: O bonded to exactly 2 H (and nothing else)
            water_atoms: set[int] = set()
            if mol.bonds:
                neighbors: dict[int, list[int]] = {}
                for a, b in mol.bonds:
                    neighbors.setdefault(a, []).append(b)
                    neighbors.setdefault(b, []).append(a)
                for i, a in enumerate(mol.atoms):
                    if a.element.symbol != "O":
                        continue
                    nbrs = neighbors.get(i, [])
                    h_nbrs = [j for j in nbrs if mol.atoms[j].element.symbol == "H"]
                    if len(h_nbrs) == 2 and len(nbrs) == 2:
                        water_atoms.add(i)
                        water_atoms.update(h_nbrs)
            if water_atoms:
                keep = [i for i in range(len(mol.atoms)) if i not in water_atoms]
                keep_set = set(keep)
                old_to_new = {old: new for new, old in enumerate(keep)}
                new_atoms = [mol.atoms[i] for i in keep]
                new_bonds = [
                    (old_to_new[a], old_to_new[b])
                    for a, b in mol.bonds if a in keep_set and b in keep_set
                ]
                mol = Molecule(atoms=new_atoms, bonds=new_bonds, pbc=mol.pbc)
        if not self.show_bonds:
            mol = Molecule(atoms=mol.atoms, bonds=[], pbc=mol.pbc)

        hl = self.highlighted_atoms if (self.highlighted_atoms and self._blink_on) else None
        pixels, hit = render_scene(
            px_w, px_h, mol, rot, self.camera_distance,
            bg_color=bg, pbc=self.pbc, ssaa=1,
            pan=(self.pan_x, self.pan_y), highlighted_atoms=hl,
            licorice=self.licorice, vdw=self.vdw,
            ambient=self.ambient, diffuse=self.diffuse,
            specular=self.specular, shininess=self.shininess,
            atom_scale=self.atom_scale, bond_radius=self.bond_radius,
            color_by_charge=self.color_by_charge,
            color_mode=self.color_mode,
            centroid_override=self._locked_centroid,
            ribbon=self.ribbon,
            show_polyhedra=self.show_polyhedra,
            density_positions=self._density_positions,
            density_values=self._density_values,
            orthographic=self.orthographic,
        )

        blocks = pixels.reshape(rows, 4, cols, 2, 3)
        is_on = hit.reshape(rows, 4, cols, 2)

        on_count = is_on.sum(axis=(1, 3))
        on_mask = is_on[:, :, :, :, None]
        color_sum = (blocks * on_mask).sum(axis=(1, 3))
        safe_count = np.maximum(on_count, 1)[:, :, None]
        avg_fg = (color_sum / safe_count).astype(np.uint8)

        braille_bits = np.where(is_on, _BRAILLE_MAP[None, :, None, :], 0)
        codepoints = 0x2800 + braille_bits.sum(axis=(1, 3)).astype(np.uint32)

        any_hit = on_count > 0
        bg_style = Style(bgcolor=f"rgb({bg[0]},{bg[1]},{bg[2]})")

        label_cells: dict[tuple[int, int], tuple[str, Style]] = {}
        if self.show_atom_numbers and self.molecule is not None:
            fov = 1.5
            scale = min(px_w, px_h) / 2
            centroid = self.molecule.center()
            label_style = Style(
                color="rgb(255,255,0)" if self.dark_bg else "rgb(0,0,180)",
                bgcolor=f"rgb({bg[0]},{bg[1]},{bg[2]})",
                bold=True,
            )
            for idx, atom in enumerate(self.molecule.atoms):
                pos = rot @ (atom.position - centroid)
                pos[0] += self.pan_x
                pos[1] += self.pan_y
                pos[2] += self.camera_distance
                if pos[2] <= 0.1:
                    continue
                sx = px_w / 2 + pos[0] * fov / pos[2] * scale
                sy = px_h / 2 - pos[1] * fov / pos[2] * scale
                cell_col = int(sx / 2)
                cell_row = int(sy / 4) - 1
                label = str(idx + 1)
                for ci, ch in enumerate(label):
                    c = cell_col + ci
                    if 0 <= cell_row < rows and 0 <= c < cols:
                        label_cells[(cell_row, c)] = (ch, label_style)

        # Axis labels (a, b, c) in bottom-left corner
        if self.pbc is not None:
            fov = 1.5
            scale = min(px_w, px_h) / 2
            basis = self.pbc.basis_matrix.copy()
            for i in range(3):
                basis[i] = basis[i] / (np.linalg.norm(basis[i]) + 1e-10)
            cam = self.camera_distance
            arrow_len = cam * 0.08
            ox = -cam * 0.85
            oy = -cam * 0.45
            origin = np.array([ox, oy, cam])
            axis_labels = ["a", "b", "c"]
            axis_colors = ["rgb(255,80,80)", "rgb(80,255,80)", "rgb(80,80,255)"]
            for i in range(3):
                direction = rot @ basis[i]
                tip = origin + direction * arrow_len * 1.4
                if tip[2] > 0.1:
                    sx = px_w / 2 + tip[0] * fov / tip[2] * scale
                    sy = px_h / 2 - tip[1] * fov / tip[2] * scale
                    cell_col = int(sx / 2)
                    cell_row = int(sy / 4)
                    style = Style(
                        color=axis_colors[i],
                        bgcolor=f"rgb({bg[0]},{bg[1]},{bg[2]})",
                        bold=True,
                    )
                    if 0 <= cell_row < rows and 0 <= cell_col < cols:
                        label_cells[(cell_row, cell_col)] = (axis_labels[i], style)

        # Charge colorbar on right side
        if self.color_by_charge and self.molecule is not None:
            charges = [a.charge for a in self.molecule.atoms]
            if charges:
                max_q = max(abs(min(charges)), abs(max(charges)), 0.01)
                bar_height = min(rows - 4, 16)
                bar_col = cols - 8
                bar_start = (rows - bar_height) // 2
                for bi in range(bar_height):
                    frac = 1.0 - bi / max(bar_height - 1, 1)  # top=+max, bottom=-max
                    charge_val = max_q * (2 * frac - 1)
                    cr, cg, cb = ImageRenderer._charge_color(charge_val)
                    r_row = bar_start + bi
                    if 0 <= r_row < rows:
                        # Color block
                        block_style = Style(bgcolor=f"rgb({cr},{cg},{cb})")
                        label_cells[(r_row, bar_col)] = (" ", block_style)
                        label_cells[(r_row, bar_col + 1)] = (" ", block_style)
                # Labels
                top_label = f"+{max_q:.2f}"
                mid_label = " 0.00"
                bot_label = f"-{max_q:.2f}"
                lbl_style = Style(
                    color="rgb(255,255,255)" if self.dark_bg else "rgb(0,0,0)",
                    bgcolor=f"rgb({bg[0]},{bg[1]},{bg[2]})",
                )
                for ci, ch in enumerate(top_label):
                    c = bar_col + 2 + ci
                    if 0 <= c < cols:
                        label_cells[(bar_start, c)] = (ch, lbl_style)
                mid_row = bar_start + bar_height // 2
                for ci, ch in enumerate(mid_label):
                    c = bar_col + 2 + ci
                    if 0 <= c < cols:
                        label_cells[(mid_row, c)] = (ch, lbl_style)
                for ci, ch in enumerate(bot_label):
                    c = bar_col + 2 + ci
                    if 0 <= c < cols:
                        label_cells[(bar_start + bar_height - 1, c)] = (ch, lbl_style)

        # Pinned measurement labels
        if self._measurements and self.molecule is not None:
            fov = 1.5
            scale = min(px_w, px_h) / 2
            centroid = self._locked_centroid if self._locked_centroid is not None else self.molecule.center()
            meas_style = Style(
                color="rgb(0,255,0)" if self.dark_bg else "rgb(0,120,0)",
                bgcolor=f"rgb({bg[0]},{bg[1]},{bg[2]})",
                bold=True,
            )
            for atom_indices, label_text in self._measurements:
                # Compute midpoint of all measured atoms
                valid = [i for i in atom_indices if i < len(self.molecule.atoms)]
                if not valid:
                    continue
                mid = np.mean([self.molecule.atoms[i].position for i in valid], axis=0)
                pos = rot @ (mid - centroid)
                pos[0] += self.pan_x
                pos[1] += self.pan_y
                pos[2] += self.camera_distance
                if pos[2] <= 0.1:
                    continue
                sx = px_w / 2 + pos[0] * fov / pos[2] * scale
                sy = px_h / 2 - pos[1] * fov / pos[2] * scale
                cell_col = int(sx / 2) + 1
                cell_row = int(sy / 4) - 1
                for ci, ch in enumerate(label_text):
                    c = cell_col + ci
                    if 0 <= cell_row < rows and 0 <= c < cols:
                        label_cells[(cell_row, c)] = (ch, meas_style)

        strips = []
        for row in range(rows):
            segments = []
            prev_style = None
            run_chars: list[str] = []
            for x in range(cols):
                if (row, x) in label_cells:
                    ch, style = label_cells[(row, x)]
                elif self.dark_bg:
                    cp = int(codepoints[row, x])
                    if cp == 0x2800:
                        style = bg_style
                        ch = " "
                    else:
                        fg = avg_fg[row, x]
                        style = Style(
                            color=f"rgb({int(fg[0])},{int(fg[1])},{int(fg[2])})",
                            bgcolor=f"rgb({bg[0]},{bg[1]},{bg[2]})",
                        )
                        ch = chr(cp)
                else:
                    cp = int(codepoints[row, x])
                    if not any_hit[row, x]:
                        style = bg_style
                        ch = " "
                    else:
                        fc = avg_fg[row, x]
                        n_on = int(on_count[row, x])
                        if n_on == 0:
                            style = bg_style
                            ch = " "
                        else:
                            # Blend color with white based on coverage (feathered edges)
                            frac = n_on / 8.0
                            r = int(fc[0] * frac + bg[0] * (1 - frac))
                            g = int(fc[1] * frac + bg[1] * (1 - frac))
                            b = int(fc[2] * frac + bg[2] * (1 - frac))
                            style = Style(bgcolor=f"rgb({r},{g},{b})")
                            ch = " "

                if style == prev_style:
                    run_chars.append(ch)
                else:
                    if run_chars and prev_style is not None:
                        segments.append(Segment("".join(run_chars), prev_style))
                    run_chars = [ch]
                    prev_style = style
            if run_chars and prev_style is not None:
                segments.append(Segment("".join(run_chars), prev_style))
            strips.append(Strip(segments, cols))

        self._cached_strips = strips


def _find_nonoverlap_position(mol: Molecule) -> np.ndarray:
    """Find a position inside the unit cell that maximizes minimum distance to MOF atoms.

    Samples a grid of candidate positions in fractional coordinates, computes
    the minimum distance to any MOF atom for each candidate (using PBC min image),
    and returns the position with the largest minimum distance.
    """
    if mol.pbc is None or len(mol.atoms) == 0:
        return mol.center()

    pbc = mol.pbc
    coords = np.array([a.x for a in mol.atoms])

    # Sample a grid of candidate positions in fractional coordinates
    n = 8  # 8x8x8 = 512 candidates
    fracs = np.linspace(0.1, 0.9, n)
    best_pos = np.array([0.5, 0.5, 0.5]) @ pbc.basis_matrix
    best_min_dist = 0.0

    for fi in fracs:
        for fj in fracs:
            for fk in fracs:
                candidate = np.array([fi, fj, fk]) @ pbc.basis_matrix
                # Compute min-image distances to all atoms
                dx = coords - candidate
                frac_dx = dx @ pbc.reciprocal_basis_matrix
                wrapped = dx - np.round(frac_dx) @ pbc.basis_matrix
                dists = np.linalg.norm(wrapped, axis=1)
                min_dist = dists.min()
                if min_dist > best_min_dist:
                    best_min_dist = min_dist
                    best_pos = candidate

    return best_pos


class MenuBar(Widget):
    """Custom menu bar with File/Edit/View dropdowns and integrated tab bar."""

    DEFAULT_CSS = """
    MenuBar {
        dock: top;
        height: 1;
        background: $panel;
        color: $foreground;
        layout: horizontal;
    }
    MenuBar .menu-btn {
        width: auto;
        height: 1;
        padding: 0 2;
        background: $panel;
        text-style: bold;
    }
    MenuBar .menu-btn:hover {
        background: $accent;
    }
    MenuBar .menu-sep {
        width: 1;
        height: 1;
        background: $accent-darken-2;
    }
    MenuBar .tab-label {
        width: auto;
        height: 1;
        padding: 0 1;
        background: $surface-darken-1;
    }
    MenuBar .tab-label.active {
        background: $accent;
        text-style: bold;
    }
    MenuBar .tab-label:hover {
        background: $accent 50%;
    }
    MenuBar .tab-add {
        width: auto;
        height: 1;
        padding: 0 1;
        background: $surface-darken-1;
    }
    MenuBar .tab-add:hover {
        background: $accent;
    }
    MenuBar .menu-spacer {
        width: 1fr;
        height: 1;
    }
    """

    class MenuAction(Message):
        def __init__(self, action: str) -> None:
            super().__init__()
            self.action = action

    class TabSelected(Message):
        def __init__(self, index: int) -> None:
            super().__init__()
            self.index = index

    class TabClose(Message):
        def __init__(self, index: int, screen_x: int = 0) -> None:
            super().__init__()
            self.index = index
            self.screen_x = screen_x

    class TabAdd(Message):
        pass

    def __init__(self, title: str = "", **kwargs) -> None:
        super().__init__(**kwargs)
        self._title = title
        self._tab_names: list[str] = []
        self._active_tab = 0

    def compose(self) -> ComposeResult:
        yield Label(" File ", id="menu-file", classes="menu-btn")
        yield Label(" Edit ", id="menu-edit", classes="menu-btn")
        yield Label(" View ", id="menu-view", classes="menu-btn")
        yield Label(" Analysis ", id="menu-analysis", classes="menu-btn")
        yield Label(" <Undo ", id="menu-undo", classes="menu-btn")
        yield Label(" Redo> ", id="menu-redo", classes="menu-btn")
        yield Label("", classes="menu-sep")
        # Tab labels will be added dynamically via set_tabs()

    def set_title(self, title: str) -> None:
        self._title = title

    def set_tabs(self, names: list[str], active: int) -> None:
        self._tab_names = names
        self._active_tab = active
        self.call_after_refresh(self._rebuild_tabs)

    async def _rebuild_tabs(self) -> None:
        # Remove old tab labels and add button (must await removal)
        to_remove = [
            child for child in self.children
            if child.has_class("tab-label") or child.has_class("tab-add")
        ]
        for child in to_remove:
            await child.remove()
        for i, name in enumerate(self._tab_names):
            cls = "tab-label active" if i == self._active_tab else "tab-label"
            truncated = (name[:19] + "…") if len(name) > 20 else name
            await self.mount(Label(f" {truncated} ", id=f"tab-{i}", classes=cls))
        await self.mount(Label(" + ", id="tab-add", classes="tab-add"))

    def on_click(self, event) -> None:
        try:
            widget = self.screen.get_widget_at(event.screen_x, event.screen_y)[0]
            if not hasattr(widget, 'id') or not widget.id:
                return
            # Tab clicks
            if widget.id == "tab-add":
                self.post_message(self.TabAdd())
                return
            if widget.id.startswith("tab-") and widget.has_class("tab-label"):
                idx = int(widget.id.replace("tab-", ""))
                if idx == self._active_tab:
                    # Clicking the active tab opens a tiny context menu
                    self.post_message(self.TabClose(idx, screen_x=widget.region.x))
                else:
                    self.post_message(self.TabSelected(idx))
                return
            # Menu clicks
            if widget.id.startswith("menu-") and widget.id != "menu-title-label":
                menu_name = widget.id.replace("menu-", "")
                self.post_message(self.MenuAction(menu_name))
        except Exception:
            pass


class MenuDropdown(ModalScreen[Optional[str]]):
    """Dropdown menu that appears under a menu button."""

    BINDINGS = [
        Binding("escape", "close", ""),
        Binding("q", "close", ""),
    ]
    DEFAULT_CSS = """
    MenuDropdown {
        align: left top;
    }
    MenuDropdown > Vertical {
        width: 35;
        height: auto;
        max-height: 80%;
        background: $surface;
        border: solid $accent;
        margin: 1 0 0 0;
    }
    MenuDropdown .menu-item {
        width: 100%;
        height: 1;
        padding: 0 1;
        background: $surface;
    }
    MenuDropdown .menu-item:hover {
        background: $accent;
    }
    MenuDropdown .menu-sep {
        width: 100%;
        height: 1;
        color: $text-muted;
    }
    """

    def __init__(self, items: list[tuple[str, str]], offset: int = 0) -> None:
        super().__init__()
        self._items = items
        self._offset = offset
        self._already_dismissed = False

    def compose(self) -> ComposeResult:
        with Vertical():
            for action, label in self._items:
                if action == "---":
                    yield Label(chr(0x2500) * 33, classes="menu-sep")
                else:
                    yield Label(label, id=f"mdrop-{action}", classes="menu-item")

    def on_mount(self) -> None:
        v = self.query_one(Vertical)
        v.styles.margin = (1, 0, 0, self._offset)

    def on_click(self, event) -> None:
        # Guard against double-dismiss when click events bubble after screen close
        if not self.is_attached or self._already_dismissed:
            return
        try:
            widget = self.screen.get_widget_at(event.screen_x, event.screen_y)[0]
            if hasattr(widget, 'id') and widget.id and widget.id.startswith("mdrop-"):
                self._already_dismissed = True
                self.dismiss(widget.id.replace("mdrop-", ""))
                return
        except Exception:
            pass
        self._already_dismissed = True
        self.dismiss("")

    def action_close(self) -> None:
        if self._already_dismissed:
            return
        self._already_dismissed = True
        self.dismiss("")


MENU_FILE = [
    ("open_file", "Open File..."),
    ("---", ""),
    ("write_xyz", "Save as XYZ"),
    ("write_pdb", "Save as PDB"),
    ("write_cif", "Save as CIF"),
    ("write_poscar", "Save as POSCAR"),
    ("write_lammps", "Save as LAMMPS Data"),
    ("write_mpmc", "Save as MPMC PDB"),
    ("---", ""),
    ("frame_to_tab", "Open Frame in New Tab"),
    ("export_png", "Export PNG"),
    ("export_gif", "Export Rotation GIF"),
    ("---", ""),
    ("quit", "Quit"),
]

MENU_EDIT = [
    ("undo", "Undo"),
    ("---", ""),
    ("wrap_center", "Wrap (centered)"),
    ("wrap_forward", "Wrap (forward)"),
    ("sort", "Sort Atoms"),
    ("extend", "Extend Axis"),
    ("del_lone", "Delete Lone Atoms"),
    ("edit_h", "Edit H Distances"),
    ("update_cell", "Update Unit Cell"),
    ("substitute", "Substitute Element"),
    ("reduce_cell", "Reduce Supercell"),
    ("---", ""),
    ("qeq_charges", "Generate QEq Charges"),
]

MENU_VIEW = [
    ("toggle_bonds", "Toggle Bonds (b)"),
    ("toggle_style", "Cycle Style (v)"),
    ("toggle_bg", "Toggle Dark/Light (i)"),
    ("toggle_atom_numbers", "Atom Numbers (#)"),
    ("cycle_color_mode", "Cycle Color Mode (C)"),
    ("toggle_hide_water", "Hide Water (w)"),
    ("toggle_measure", "Measure Mode (m)"),
    ("clear_measurements", "Clear Measurements"),
    ("toggle_polyhedra", "Polyhedra (P)"),
    ("toggle_orthographic", "Orthographic Projection (O)"),
    ("---", ""),
    ("toggle_operations", "Operations Panel (o)"),
    ("toggle_geometry", "Geometry Panel (g)"),
    ("toggle_visual", "Visual Settings (V)"),
    ("---", ""),
    ("reset_view", "Reset View (r)"),
]

MENU_ANALYSIS = [
    ("info", "System Info"),
    ("void_vol", "Void Volume"),
    ("surface_area", "Surface Area"),
    ("pore_size", "Pore Size Distribution"),
    ("rdf", "Radial Distribution"),
    ("coordination", "Coordination Number"),
    ("hbonds", "Hydrogen Bonds"),
    ("msd", "Mean Sq. Displacement"),
    ("rmsd", "RMSD vs Time"),
    ("gyration", "Radius of Gyration"),
    ("density", "Density Profile"),
    ("density3d", "3D Density Map"),
    ("pxrd", "Powder XRD"),
    ("energy_plot", "Energy Plot"),
    ("---", ""),
    ("input_generator", "Input Generator"),
    ("isotherm_plot", "Plot Isotherm Results"),
    ("db_search", "Database Search"),
]


class TabState:
    """Per-tab state — isolates each open structure completely."""
    def __init__(self, molecule: Molecule, filepath: str = "", frames: list | None = None):
        self.molecule = molecule
        self.filepath = filepath
        self.frames = frames
        self.frame_idx = 0
        self.undo_stack: list = []
        self.redo_stack: list = []
        # View state
        self.rot_matrix = rotation_matrix(-0.2, -0.5, 0.0)
        self.camera_distance = 4.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.pan_mode = False
        self.dark_bg = True
        self.show_bonds = True
        self.show_atom_numbers = False
        self.color_by_charge = False
        self.hide_water = False
        self.ribbon = False
        self.show_polyhedra = False
        self.orthographic = False
        self.color_mode = "element"
        self.licorice = False
        self.vdw = False
        self.atom_scale = 0.35
        self.bond_radius = 0.08
        self.locked_centroid = None
        self.camera_initialized = False

    @property
    def name(self) -> str:
        from pathlib import Path
        return Path(self.filepath).name if self.filepath else "untitled"



class PdbWizardApp(App):
    ENABLE_COMMAND_PALETTE = False
    CSS = """
    Screen { layout: vertical; }
    #main-content { height: 1fr; }
    MoleculeView { width: 1fr; }
    #status-bar {
        height: 1;
        max-height: 1;
        overflow: hidden;
    }
    #status-bar Label {
        width: auto;
        max-width: 60;
        margin-right: 1;
        color: $text-muted;
    }
    #status-bar ProgressBar {
        width: 1fr;
        display: none;
    }
    #status-bar.busy ProgressBar {
        display: block;
    }
    #status-bar.busy Label {
        color: $text;
    }
    #traj-bar {
        height: 3;
        display: none;
        background: $accent-darken-3;
        dock: bottom;
    }
    #traj-bar.visible {
        display: block;
    }
    #traj-bar Button {
        width: auto;
        min-width: 3;
        height: 3;
        margin: 0;
    }
    #traj-bar TrackSlider {
        width: 1fr;
        height: 3;
        margin: 0 1;
    }
    #traj-bar Label {
        width: auto;
        height: 3;
        content-align: center middle;
        margin: 0 1;
    }
    """

    BINDINGS = [
        Binding("k,up", "rotate_up", "Rot up", show=False),
        Binding("j,down", "rotate_down", "Rot down", show=False),
        Binding("h,left", "rotate_left", "Rot left", show=False),
        Binding("l,right", "rotate_right", "Rot right", show=False),
        Binding("comma", "rotate_cw", ",/. roll", show=False),
        Binding("full_stop", "rotate_ccw", show=False),
        Binding("K", "zoom_in", show=False),
        Binding("J", "zoom_out", show=False),
        Binding("plus,equals_sign", "zoom_in", "Zoom", key_display="+/-"),
        Binding("minus,underscore", "zoom_out", show=False),
        Binding("t", "toggle_mode", "Pan/Rot"),
        Binding("c", "center", show=False),
        Binding("r", "reset_view", show=False),
        Binding("b", "toggle_bonds", show=False),
        Binding("v", "toggle_style", show=False),
        Binding("i", "toggle_bg", show=False),
        Binding("escape", "close_panel", show=False),
        Binding("o", "toggle_operations", show=False),
        Binding("g", "toggle_geometry", show=False),
        Binding("e", "export_png", show=False),
        Binding("number_sign", "toggle_atom_numbers", show=False),
        Binding("C", "toggle_charge_colors", show=False),
        Binding("w", "toggle_hide_water", show=False),
        Binding("m", "toggle_measure", show=False),
        Binding("P", "toggle_polyhedra", show=False),
        Binding("O", "toggle_orthographic", show=False),
        Binding("right_square_bracket", "next_frame", "Frame", key_display="[/]"),
        Binding("left_square_bracket", "prev_frame", show=False),
        Binding("home", "first_frame", show=False),
        Binding("end", "last_frame", show=False),
        Binding("space", "toggle_play", show=False),
        Binding("n", "panel_next", show=False),
        Binding("p", "panel_prev", show=False),
        Binding("V", "toggle_visual", show=False),
        Binding("q", "request_quit", "Quit"),
    ]

    def __init__(self, molecule: Molecule, filepath: str = "", frames: list | None = None,
                 load_trajectory: bool = False):
        super().__init__()
        # Tab system
        self._tabs: list[TabState] = [TabState(molecule, filepath, frames)]
        self._active_tab = 0
        # Convenience accessors (always point to active tab)
        self.molecule = molecule
        self.filepath = filepath
        self._frames = frames
        self._frame_idx = 0
        self._load_trajectory = load_trajectory
        self._update_title()
        self._quit_pending = False
        self._quit_timer = None
        self._play_timer = None
        self._play_speed = 1.0
        self._undo_stack: list = []
        self._redo_stack: list = []
        self._play_interval = 0.15
        # Background tasks tracked so we can cancel pending ones on app exit.
        # Without this, "Task was destroyed but it is pending" warnings flood
        # stderr when the user quits while async work is in flight.
        self._bg_tasks: set[asyncio.Task] = set()
        # Per-session cache for downloaded structures. Created lazily on first
        # fetch; removed on app teardown.
        self._db_cache_dir_path: str | None = None

    def _track_task(self, coro) -> asyncio.Task:
        """Create an asyncio task and track it for cleanup."""
        task = asyncio.create_task(coro)
        self._bg_tasks.add(task)
        task.add_done_callback(self._bg_tasks.discard)
        return task

    def _db_cache_dir(self) -> str:
        """Return a single per-session directory for downloaded structures.
        Created on first call, cleaned up in on_unmount."""
        if self._db_cache_dir_path is None:
            import tempfile
            self._db_cache_dir_path = tempfile.mkdtemp(prefix="pdb_wizard_db_")
        return self._db_cache_dir_path

    async def on_unmount(self) -> None:
        """Cancel any in-flight background tasks and clean caches."""
        for task in list(self._bg_tasks):
            if not task.done():
                task.cancel()
        self._bg_tasks.clear()
        if self._db_cache_dir_path is not None:
            import shutil
            shutil.rmtree(self._db_cache_dir_path, ignore_errors=True)
            self._db_cache_dir_path = None

    def _refresh_tab_bar(self) -> None:
        names = self._disambiguated_tab_names()
        try:
            self.query_one("#app-menu", MenuBar).set_tabs(names, self._active_tab)
        except Exception:
            pass

    def _disambiguated_tab_names(self) -> list[str]:
        """Return tab names with parent directory prefixed for duplicates."""
        from collections import Counter
        base_names = [t.name for t in self._tabs]
        counts = Counter(base_names)
        result = []
        for tab in self._tabs:
            base = tab.name
            if counts[base] > 1 and tab.filepath:
                parent = Path(tab.filepath).parent.name
                if parent:
                    result.append(f"{parent}/{base}")
                else:
                    result.append(base)
            else:
                result.append(base)
        return result

    def _save_current_tab_state(self) -> None:
        """Save current view state back to the active tab."""
        tab = self._tabs[self._active_tab]
        view = self.query_one(MoleculeView)
        tab.molecule = self.molecule
        tab.frames = self._frames
        tab.frame_idx = self._frame_idx
        tab.undo_stack = self._undo_stack
        tab.redo_stack = self._redo_stack
        tab.rot_matrix = view.rot_matrix
        tab.camera_distance = view.camera_distance
        tab.pan_x = view.pan_x
        tab.pan_y = view.pan_y
        tab.pan_mode = view.pan_mode
        tab.dark_bg = view.dark_bg
        tab.show_bonds = view.show_bonds
        tab.show_atom_numbers = view.show_atom_numbers
        tab.color_by_charge = view.color_by_charge
        tab.hide_water = view.hide_water
        tab.ribbon = view.ribbon
        tab.show_polyhedra = view.show_polyhedra
        tab.orthographic = view.orthographic
        tab.color_mode = view.color_mode
        tab.licorice = view.licorice
        tab.vdw = view.vdw
        tab.atom_scale = view.atom_scale
        tab.bond_radius = view.bond_radius
        tab.locked_centroid = view._locked_centroid

    def _restore_tab_state(self, idx: int) -> None:
        """Restore state from tab idx to the view."""
        tab = self._tabs[idx]
        self._active_tab = idx
        self.molecule = tab.molecule
        self.filepath = tab.filepath
        self._frames = tab.frames
        self._frame_idx = tab.frame_idx
        self._undo_stack = tab.undo_stack
        self._redo_stack = tab.redo_stack

        view = self.query_one(MoleculeView)
        view.rot_matrix = tab.rot_matrix
        view.camera_distance = tab.camera_distance
        view.pan_x = tab.pan_x
        view.pan_y = tab.pan_y
        view.pan_mode = tab.pan_mode
        view.dark_bg = tab.dark_bg
        view.show_bonds = tab.show_bonds
        view.show_atom_numbers = tab.show_atom_numbers
        view.color_by_charge = tab.color_by_charge
        view.hide_water = tab.hide_water
        view.ribbon = tab.ribbon
        view.show_polyhedra = tab.show_polyhedra
        view.orthographic = getattr(tab, "orthographic", False)
        view.color_mode = tab.color_mode
        view.licorice = tab.licorice
        view.vdw = tab.vdw
        view.atom_scale = tab.atom_scale
        view.bond_radius = tab.bond_radius
        view._locked_centroid = tab.locked_centroid

        if len(self.molecule.atoms) > 0:
            if tab.camera_initialized:
                view.set_molecule(self.molecule, keep_camera=True)
            else:
                view.set_molecule(self.molecule, keep_camera=False)
                tab.camera_initialized = True
            self.query_one(GeometryPanel).set_molecule(self.molecule)
        self._update_title()
        self._refresh_tab_bar()
        self._setup_traj_bar()
        view._invalidate_cache()

    def open_in_new_tab(self, mol: Molecule, filepath: str = "",
                        frames: list | None = None) -> int:
        """Append a new tab and switch to it. Returns the new tab's index.

        Centralizes the save → append → restore sequence so callers (modal
        screens, file dialogs, DB fetch) don't have to reach into the
        internals of self._tabs and self._active_tab. Stops playback first
        so a running animation on the old tab doesn't bleed into the new one.
        """
        if self._play_timer:
            self._play_timer.stop()
            self._play_timer = None
        self._save_current_tab_state()
        self._tabs.append(TabState(mol, filepath, frames))
        new_idx = len(self._tabs) - 1
        self._restore_tab_state(new_idx)
        return new_idx

    def on_menu_bar_tab_selected(self, event: MenuBar.TabSelected) -> None:
        if event.index == self._active_tab:
            return
        if event.index < 0 or event.index >= len(self._tabs):
            return
        # Stop playback
        if self._play_timer:
            self._play_timer.stop()
            self._play_timer = None
        self._save_current_tab_state()
        self._restore_tab_state(event.index)

    def on_menu_bar_tab_close(self, event: MenuBar.TabClose) -> None:
        if len(self._tabs) <= 1:
            self.notify("Can't close the last tab", timeout=2)
            return
        idx = event.index
        disambig = self._disambiguated_tab_names()
        name = disambig[idx] if idx < len(disambig) else self._tabs[idx].name
        items = [("close_tab", f"Close {name}")]
        self.push_screen(
            MenuDropdown(items, offset=event.screen_x),
            callback=lambda action: self._confirm_close_tab(action, idx),
        )

    def _confirm_close_tab(self, action: str, idx: int) -> None:
        if action != "close_tab":
            return
        if idx < 0 or idx >= len(self._tabs) or len(self._tabs) <= 1:
            return
        # Stop playback
        if self._play_timer:
            self._play_timer.stop()
            self._play_timer = None
        closed_name = self._tabs[idx].name
        self._tabs.pop(idx)
        if self._active_tab >= len(self._tabs):
            self._active_tab = len(self._tabs) - 1
        elif self._active_tab > idx:
            self._active_tab -= 1
        self._restore_tab_state(self._active_tab)
        self.notify(f"Closed {closed_name}", timeout=2)

    def on_menu_bar_tab_add(self, event: MenuBar.TabAdd) -> None:
        self.push_screen(
            FileSaveModal("Open file:", default="", button_label="Open"),
            callback=self._open_new_tab,
        )

    def _open_new_tab(self, filepath: str) -> None:
        if not filepath.strip():
            return
        filepath = filepath.strip()
        try:
            ft = detect_filetype(filepath)
            is_traj = False
            if ft == "pdb":
                is_traj = check_pdb_trajectory(filepath)
            elif ft == "xyz":
                is_traj = check_xyz_trajectory(filepath)
            elif ft == "dcd":
                is_traj = check_dcd_trajectory(filepath)

            if is_traj:
                # Load trajectory in background into a new tab
                # Don't switch yet — stay on current tab until loading finishes
                self._save_current_tab_state()
                empty = Molecule(atoms=[])
                new_tab = TabState(empty, filepath)
                self._tabs.append(new_tab)
                self._refresh_tab_bar()
                self.notify(f"Loading {new_tab.name}...", timeout=2)
                self._track_task(self._load_traj_into_tab(filepath, len(self._tabs) - 1))
            else:
                mol = read_file(filepath)
                self.open_in_new_tab(mol, filepath)
                self.notify(f"Opened {self._tabs[-1].name}", timeout=2)
        except Exception as e:
            self.notify(f"Error: {e}", timeout=3)

    async def _load_traj_into_tab(self, filepath: str, tab_idx: int) -> None:
        """Load a trajectory file into a specific tab in the background."""

        self._show_progress("Loading trajectory...", total=100)
        await asyncio.sleep(0)

        ft = detect_filetype(filepath)
        _progress = [0.0]

        def _prog(frac: float) -> None:
            _progress[0] = frac

        def _poll_progress() -> None:
            pct = int(_progress[0] * 80)
            self._update_progress(pct)
            self.query_one("#status-label", Label).update(
                f"Loading trajectory... {int(_progress[0] * 100)}%"
            )

        timer = self.set_interval(0.1, _poll_progress)

        def _read():
            if ft == "dcd":
                mols, _ = read_dcd_trajectory(filepath, progress_callback=_prog)
                return mols
            with open(filepath) as f:
                if ft == "pdb":
                    mols, _ = read_pdb_trajectory(f, progress_callback=_prog)
                else:
                    mols, _ = read_xyz_trajectory(f, progress_callback=_prog)
            return mols

        frames = await asyncio.to_thread(_read)
        timer.stop()

        # If DCD loaded without a topology file, atoms are all carbon —
        # surface that to the user.
        if frames and getattr(frames[0], "_dcd_topology_missing", False):
            self.notify(
                "DCD has no element info. All atoms loaded as carbon. "
                "Place a sibling .pdb or .xyz with the same atom count "
                "next to the .dcd file.",
                severity="warning", timeout=8,
            )

        self.query_one("#status-label", Label).update("Detecting bonds...")
        self._update_progress(85)
        await asyncio.sleep(0)

        if frames:
            await asyncio.to_thread(frames[0].detect_bonds)

        # Update the tab state
        if tab_idx < len(self._tabs):
            tab = self._tabs[tab_idx]
            tab.frames = frames
            tab.molecule = frames[0] if frames else Molecule(atoms=[])

        self._hide_progress()

        # Switch to the newly loaded tab
        if tab_idx < len(self._tabs):
            self._save_current_tab_state()
            self._restore_tab_state(tab_idx)
        self.notify(f"Loaded {len(frames)} frames", timeout=2)

    def _update_title(self) -> None:
        # Use the disambiguated tab name (includes parent dir if name duplicated)
        try:
            disambig = self._disambiguated_tab_names()
            name = disambig[self._active_tab] if 0 <= self._active_tab < len(disambig) else "PDB Wizard"
        except Exception:
            name = Path(self.filepath).name if self.filepath else "PDB Wizard"
        if self._frames and len(self._frames) > 1:
            title = f"{name}  frame {self._frame_idx + 1}/{len(self._frames)}"
        else:
            title = name
        self.title = title
        try:
            self.query_one("#app-menu", MenuBar).set_title(title)
        except Exception:
            pass

    def action_request_quit(self) -> None:
        if self._quit_pending:
            self.exit()
        else:
            self._quit_pending = True
            self.notify("Press q again to quit", timeout=3)
            self._quit_timer = self.set_timer(3.0, self._cancel_quit)

    def _cancel_quit(self) -> None:
        self._quit_pending = False

    def compose(self) -> ComposeResult:
        yield MenuBar(title=self.title, id="app-menu")
        with Horizontal(id="main-content"):
            yield OperationsPanel()
            yield MoleculeView()
            yield GeometryPanel()
            yield VisualPanel()
        with Horizontal(id="status-bar"):
            yield Label("", id="status-label")
            yield ProgressBar(id="status-progress", total=100, show_eta=False)
        with Horizontal(id="traj-bar"):
            yield Button("|<", id="traj-first")
            yield Button("<", id="traj-prev")
            yield Button("Play", id="traj-play")
            yield Button(">", id="traj-next")
            yield Button(">|", id="traj-last")
            yield Button("Slow", id="traj-slower")
            yield Label("1.0x", id="traj-speed")
            yield Button("Fast", id="traj-faster")
            yield TrackSlider(id="traj-slider")
            yield Label("1/1", id="traj-label")
        yield Footer()

    def on_mount(self) -> None:
        # Open operations panel by default
        ops = self.query_one(OperationsPanel)
        ops.add_class("visible")
        for dt in ops.query(DataTable):
            dt.focus()
            break

        # Initialize tab bar
        self._refresh_tab_bar()

        if self._load_trajectory:
            # Load trajectory in background
            self._track_task(self._load_trajectory_async())
        elif len(self.molecule.atoms) > 0:
            self._show_progress("Loading molecule...")
            view = self.query_one(MoleculeView)
            view.set_molecule(self.molecule)
            self._tabs[self._active_tab].camera_initialized = True
            self._show_progress("Building geometry tables...")
            panel = self.query_one(GeometryPanel)
            panel.set_molecule(self.molecule)
            self._hide_progress()
            self._setup_traj_bar()

    async def _load_trajectory_async(self) -> None:

        self._show_progress("Loading trajectory...", total=100)
        await asyncio.sleep(0)

        ft = detect_filetype(self.filepath)

        def _prog(frac: float) -> None:
            self.call_from_thread(self._update_progress, int(frac * 80))
            self.call_from_thread(
                self.query_one("#status-label", Label).update,
                f"Loading trajectory... {int(frac * 100)}%"
            )

        def _read():
            with open(self.filepath) as f:
                if ft == "pdb":
                    mols, _ = read_pdb_trajectory(f, progress_callback=_prog)
                else:
                    mols, _ = read_xyz_trajectory(f, progress_callback=_prog)
            return mols

        frames = await asyncio.to_thread(_read)

        self.query_one("#status-label", Label).update("Detecting bonds...")
        self._update_progress(85)
        await asyncio.sleep(0)

        # Detect bonds on first frame only (for display), rest on demand
        if frames:
            await asyncio.to_thread(frames[0].detect_bonds)

        self._update_progress(80)

        self._frames = frames
        self._frame_idx = 0
        if frames:
            self.molecule = frames[0]

        view = self.query_one(MoleculeView)
        if self.molecule and len(self.molecule.atoms) > 0:
            view.set_molecule(self.molecule)
            self.query_one(GeometryPanel).set_molecule(self.molecule)

        self._update_progress(90)
        self.query_one("#status-label", Label).update("Ready")
        await asyncio.sleep(0)
        self._hide_progress()
        self._update_title()
        self._setup_traj_bar()

    def _setup_traj_bar(self) -> None:
        traj_bar = self.query_one("#traj-bar")
        if self._frames and len(self._frames) > 1:
            traj_bar.add_class("visible")
            n = len(self._frames)
            idx = self._frame_idx + 1
            self.query_one("#traj-slider", TrackSlider).set_position(idx, n)
            self.query_one("#traj-label", Label).update(f"{idx}/{n}")
        else:
            traj_bar.remove_class("visible")

    def _show_progress(self, label: str, total: int = 100) -> None:
        bar = self.query_one("#status-bar")
        bar.add_class("busy")
        self.query_one("#status-label", Label).update(label)
        pb = self.query_one("#status-progress", ProgressBar)
        pb.update(total=total, progress=0)

    def _update_progress(self, progress: int) -> None:
        self.query_one("#status-progress", ProgressBar).update(progress=progress)

    def _hide_progress(self) -> None:
        bar = self.query_one("#status-bar")
        bar.remove_class("busy")
        # Clear the label so the always-visible bar looks empty when idle.
        self.query_one("#status-label", Label).update("")

    def _run_with_progress(self, label: str, func, *args, **kwargs):
        """Run a function async with a progress bar shown."""
        self._show_progress(label)
        self._track_task(self._run_bg_task(func, *args, **kwargs))

    async def _run_bg_task(self, func, *args, **kwargs):
        try:
            result = await asyncio.to_thread(func, *args, **kwargs)
        finally:
            self._hide_progress()
        return result

    async def _populate_tables_async(self, panel: GeometryPanel, mol: Molecule) -> None:
        """Populate geometry tables one at a time, yielding between each for UI updates."""
        panel._populating = True

        # Use cached data
        bonds = getattr(mol, '_cached_bond_lengths', None) or mol.get_bond_lengths()
        self.query_one("#status-label", Label).update("Populating bonds table...")
        self._update_progress(86)
        await asyncio.sleep(0)
        table = panel.query_one("#bonds-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Atom 2", "Length (\u00c5)")
        for i, j, dist in bonds:
            table.add_row(panel._atom_label(i), panel._atom_label(j), f"{dist:.4f}", key=f"{i}-{j}")

        angles = getattr(mol, '_cached_angles', None) or mol.get_angles()
        self.query_one("#status-label", Label).update("Populating angles table...")
        self._update_progress(88)
        await asyncio.sleep(0)
        table = panel.query_one("#angles-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Vertex", "Atom 3", "Angle (\u00b0)")
        for i, j, k, angle in angles:
            table.add_row(
                panel._atom_label(i), panel._atom_label(j), panel._atom_label(k),
                f"{angle:.3f}", key=f"{i}-{j}-{k}",
            )

        dihedrals = getattr(mol, '_cached_dihedrals', None) or mol.get_dihedrals()
        self.query_one("#status-label", Label).update("Populating dihedrals table...")
        self._update_progress(90)
        await asyncio.sleep(0)
        table = panel.query_one("#dihedrals-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Atom 2", "Atom 3", "Atom 4", "Angle (\u00b0)")
        for i, j, k, l_idx, angle in dihedrals:
            table.add_row(
                panel._atom_label(i), panel._atom_label(j),
                panel._atom_label(k), panel._atom_label(l_idx),
                f"{angle:.3f}", key=f"{i}-{j}-{k}-{l_idx}",
            )

        contacts = getattr(mol, '_cached_contacts', None) or []
        self.query_one("#status-label", Label).update("Populating contacts table...")
        self._update_progress(92)
        await asyncio.sleep(0)
        table = panel.query_one("#contacts-table", DataTable)
        table.clear(columns=True)
        table.add_columns("Atom 1", "Atom 2", "Distance (\u00c5)")
        for msg in contacts:
            parts = msg.split()
            el_pair = parts[0]
            id1, id2 = int(parts[1]), int(parts[2])
            r_val = parts[-1]
            table.add_row(
                f"{id1}:{el_pair.split('-')[0]}", f"{id2}:{el_pair.split('-')[1]}",
                r_val, key=f"{id1 - 1}-{id2 - 1}",
            )

        submols = getattr(mol, '_cached_molecules', None) or mol.find_molecules()
        self.query_one("#status-label", Label).update("Populating molecules table...")
        self._update_progress(94)
        await asyncio.sleep(0)
        table = panel.query_one("#molecules-table", DataTable)
        table.clear(columns=True)
        table.add_columns("#", "Atoms", "Formula")
        for mi, sm in enumerate(submols):
            sm.atoms.sort(key=lambda a: a.atomic_number, reverse=True)
            from collections import Counter
            counts = Counter(a.element.symbol for a in sm.atoms)
            formula = "".join(f"{el}{n}" if n > 1 else el for el, n in counts.items())
            atom_indices = "-".join(str(mol.atoms.index(a)) for a in sm.atoms)
            table.add_row(str(mi + 1), str(len(sm.atoms)), formula, key=atom_indices)

        lone = getattr(mol, '_cached_lone_atoms', None) or []
        self.query_one("#status-label", Label).update("Populating coords table...")
        self._update_progress(96)
        await asyncio.sleep(0)
        table = panel.query_one("#lone-table", DataTable)
        table.clear(columns=True)
        table.add_columns("#", "Element", "x", "y", "z")
        for a in lone:
            idx = mol.atoms.index(a)
            table.add_row(
                str(a.id), a.element.symbol,
                f"{a.x[0]:.4f}", f"{a.x[1]:.4f}", f"{a.x[2]:.4f}",
                key=str(idx),
            )

        self._update_progress(98)
        await asyncio.sleep(0)
        table = panel.query_one("#coords-table", DataTable)
        table.clear(columns=True)
        table.add_columns("#", "El", "x", "y", "z")
        for idx, a in enumerate(mol.atoms):
            table.add_row(
                str(a.id), a.element.symbol,
                f"{a.x[0]:.4f}", f"{a.x[1]:.4f}", f"{a.x[2]:.4f}",
                key=str(idx),
            )

        # Clear caches
        for attr in ('_cached_bond_lengths', '_cached_angles', '_cached_dihedrals',
                     '_cached_contacts', '_cached_molecules', '_cached_lone_atoms'):
            if hasattr(mol, attr):
                delattr(mol, attr)

        panel._populating = False

    @staticmethod
    def _precompute_geometry(mol: Molecule, progress=None) -> None:
        """Pre-compute geometry data in a thread so the panel populate is fast."""
        set_atom_ids(mol.atoms)

        if progress:
            progress(0, "Computing bond lengths...")
        mol._cached_bond_lengths = mol.get_bond_lengths()

        if progress:
            progress(5, "Computing angles...")
        mol._cached_angles = mol.get_angles()

        if progress:
            progress(15, "Computing dihedrals...")
        mol._cached_dihedrals = mol.get_dihedrals()

        if progress:
            progress(30, "Finding close contacts...")
        mol._cached_contacts = get_close_contacts(mol)

        if progress:
            progress(55, "Finding molecules...")
        mol._cached_molecules = mol.find_molecules()

        if progress:
            progress(75, "Finding lone atoms...")
        mol._cached_lone_atoms = get_lone_atoms(mol)

        if progress:
            progress(100, "Done computing geometry")

    def _panel_is_open(self) -> bool:
        return (
            self.query_one(OperationsPanel).has_class("visible")
            or self.query_one(GeometryPanel).has_class("visible")
            or self.query_one(VisualPanel).has_class("visible")
        )

    def check_action(self, action: str, parameters: tuple[object, ...]) -> bool | None:
        if action in ("panel_next", "panel_prev"):
            if not self._panel_is_open():
                return False
        if action in ("next_frame", "prev_frame", "first_frame", "last_frame"):
            if not self._frames or len(self._frames) <= 1:
                return False
        return True

    # --- Rotation / Zoom / Pan ---

    def action_rotate_up(self) -> None:
        view = self.query_one(MoleculeView)
        if view.pan_mode:
            view.pan_y -= view.camera_distance * 0.05
            view._clamp_pan()
        else:
            view.rot_matrix = rotation_matrix(-0.1, 0, 0) @ view.rot_matrix
        view._invalidate_cache()

    def action_rotate_down(self) -> None:
        view = self.query_one(MoleculeView)
        if view.pan_mode:
            view.pan_y += view.camera_distance * 0.05
            view._clamp_pan()
        else:
            view.rot_matrix = rotation_matrix(0.1, 0, 0) @ view.rot_matrix
        view._invalidate_cache()

    def action_rotate_left(self) -> None:
        view = self.query_one(MoleculeView)
        if view.pan_mode:
            view.pan_x += view.camera_distance * 0.05
            view._clamp_pan()
        else:
            view.rot_matrix = rotation_matrix(0, 0.1, 0) @ view.rot_matrix
        view._invalidate_cache()

    def action_rotate_right(self) -> None:
        view = self.query_one(MoleculeView)
        if view.pan_mode:
            view.pan_x -= view.camera_distance * 0.05
            view._clamp_pan()
        else:
            view.rot_matrix = rotation_matrix(0, -0.1, 0) @ view.rot_matrix
        view._invalidate_cache()

    def action_rotate_cw(self) -> None:
        view = self.query_one(MoleculeView)
        view.rot_matrix = rotation_matrix(0, 0, 0.1) @ view.rot_matrix
        view._invalidate_cache()

    def action_rotate_ccw(self) -> None:
        view = self.query_one(MoleculeView)
        view.rot_matrix = rotation_matrix(0, 0, -0.1) @ view.rot_matrix
        view._invalidate_cache()

    def action_zoom_in(self) -> None:
        view = self.query_one(MoleculeView)
        view.camera_distance = max(1.0, view.camera_distance - 0.5)
        view._invalidate_cache()

    def action_zoom_out(self) -> None:
        view = self.query_one(MoleculeView)
        view.camera_distance += 0.5
        view._invalidate_cache()

    def action_toggle_mode(self) -> None:
        view = self.query_one(MoleculeView)
        view.pan_mode = not view.pan_mode
        mode = "PAN" if view.pan_mode else "ROT"
        self.notify(f"Mode: {mode}", timeout=1)

    def action_center(self) -> None:
        view = self.query_one(MoleculeView)
        view.pan_x = 0.0
        view.pan_y = 0.0
        view._invalidate_cache()

    # --- Display toggles ---

    def action_toggle_style(self) -> None:
        view = self.query_one(MoleculeView)
        # Cycle: CPK -> Licorice -> VDW -> Ribbon -> CPK
        if not view.licorice and not view.vdw and not view.ribbon:
            view.licorice = True
            view.vdw = False
            view.ribbon = False
            view.bond_radius = 0.15
        elif view.licorice:
            view.licorice = False
            view.vdw = True
            view.ribbon = False
            view.bond_radius = 0.08
        elif view.vdw:
            view.licorice = False
            view.vdw = False
            view.ribbon = True
            view.bond_radius = 0.08
        else:
            view.licorice = False
            view.vdw = False
            view.ribbon = False
            view.bond_radius = 0.08
        vis = self.query_one(VisualPanel)
        if vis.has_class("visible"):
            vis.set_state(
                licorice=view.licorice, vdw=view.vdw, ribbon=view.ribbon,
                ambient=view.ambient, diffuse=view.diffuse,
                specular=view.specular, shininess=view.shininess,
                atom_scale=view.atom_scale, bond_radius=view.bond_radius,
            )
        view._invalidate_cache()

    def action_toggle_atom_numbers(self) -> None:
        view = self.query_one(MoleculeView)
        view.show_atom_numbers = not view.show_atom_numbers
        view._invalidate_cache()

    def action_toggle_charge_colors(self) -> None:
        self._cycle_color_mode()

    def _cycle_color_mode(self) -> None:
        view = self.query_one(MoleculeView)
        modes = ["element", "charge", "residue", "chain", "index"]
        try:
            idx = modes.index(view.color_mode)
        except ValueError:
            idx = 0
        view.color_mode = modes[(idx + 1) % len(modes)]
        view.color_by_charge = (view.color_mode == "charge")
        self.notify(f"Color: {view.color_mode}", timeout=1)
        view._invalidate_cache()

    def action_toggle_hide_water(self) -> None:
        view = self.query_one(MoleculeView)
        view.hide_water = not view.hide_water
        label = "hidden" if view.hide_water else "visible"
        self.notify(f"Water: {label}", timeout=1)
        view._invalidate_cache()

    def action_toggle_measure(self) -> None:
        view = self.query_one(MoleculeView)
        view.measure_mode = not view.measure_mode
        if view.measure_mode:
            view._picked_atoms = []
            self.notify("Measure mode ON — click atoms to measure distances/angles", timeout=3)
        else:
            view._picked_atoms = []
            view.highlighted_atoms = set()
            view._stop_blink()
            self.notify("Measure mode OFF", timeout=1)
        view._invalidate_cache()

    def action_toggle_polyhedra(self) -> None:
        view = self.query_one(MoleculeView)
        view.show_polyhedra = not view.show_polyhedra
        label = "ON" if view.show_polyhedra else "OFF"
        self.notify(f"Polyhedra: {label}", timeout=1)
        view._invalidate_cache()

    def action_toggle_orthographic(self) -> None:
        view = self.query_one(MoleculeView)
        view.orthographic = not view.orthographic
        label = "Orthographic" if view.orthographic else "Perspective"
        self.notify(f"Projection: {label}", timeout=1)
        view._invalidate_cache()

    def _clear_measurements(self) -> None:
        view = self.query_one(MoleculeView)
        view._measurements.clear()
        self.notify("Measurements cleared", timeout=1)
        view._invalidate_cache()

    def _switch_frame(self, idx: int) -> None:
        if not self._frames or len(self._frames) <= 1:
            return
        idx = max(0, min(idx, len(self._frames) - 1))
        if idx == self._frame_idx:
            return
        self._frame_idx = idx
        self.molecule = self._frames[idx]
        # Reuse bonds from frame 0 if this frame has none (same topology)
        if not self.molecule.bonds and len(self.molecule.atoms) > 0:
            if self._frames[0].bonds and len(self._frames[0].atoms) == len(self.molecule.atoms):
                self.molecule.bonds = self._frames[0].bonds
            else:
                self.molecule.detect_bonds()
        view = self.query_one(MoleculeView)
        # Lock centroid to first frame so MOF doesn't jitter
        if view._locked_centroid is None:
            view._locked_centroid = self._frames[0].center()
        view.set_molecule(self.molecule, keep_camera=True)
        self._update_title()
        n = len(self._frames)
        self.query_one("#traj-slider", TrackSlider).set_position(idx + 1, n)
        self.query_one("#traj-label", Label).update(f"{idx + 1}/{n}")
        # Update geometry panel only when not auto-playing (expensive)
        if self._play_timer is None:
            self.query_one(GeometryPanel).set_molecule(self.molecule)
        view._invalidate_cache()

    def on_track_slider_seeked(self, event: TrackSlider.Seeked) -> None:
        if self._frames:
            idx = int(event.fraction * (len(self._frames) - 1))
            self._switch_frame(idx)

    def _on_traj_button(self, event: Button.Pressed) -> None:
        bid = event.button.id
        if bid == "traj-first":
            self._switch_frame(0)
        elif bid == "traj-prev":
            self._switch_frame(self._frame_idx - 1)
        elif bid == "traj-next":
            self._switch_frame(self._frame_idx + 1)
        elif bid == "traj-last":
            if self._frames:
                self._switch_frame(len(self._frames) - 1)
        elif bid == "traj-play":
            self._toggle_play()
        elif bid == "traj-slower":
            self._play_speed = max(0.125, self._play_speed / 2)
            self._update_play_speed()
        elif bid == "traj-faster":
            self._play_speed = min(16.0, self._play_speed * 2)
            self._update_play_speed()

    def _update_play_speed(self) -> None:
        self.query_one("#traj-speed", Label).update(f"{self._play_speed:.1f}x")
        # Restart timer at new speed if playing
        if self._play_timer is not None:
            self._play_timer.stop()
            interval = self._play_interval / self._play_speed
            self._play_timer = self.set_interval(interval, self._play_step)

    def _toggle_play(self) -> None:
        if not self._frames:
            return
        if self._play_timer is not None:
            self._play_timer.stop()
            self._play_timer = None
            self.query_one("#traj-play", Button).label = "Play"
            # Update geometry panel now that playback stopped
            self.query_one(GeometryPanel).set_molecule(self.molecule)
        else:
            self.query_one("#traj-play", Button).label = "Pause"
            interval = self._play_interval / self._play_speed
            self._play_timer = self.set_interval(interval, self._play_step)

    def _play_step(self) -> None:
        if not self._frames:
            return
        # Skip frames if rendering is slow — advance by play_speed to keep up
        step = max(1, int(self._play_speed))
        next_idx = self._frame_idx + step
        if next_idx >= len(self._frames):
            next_idx = 0  # loop
        self._switch_frame(next_idx)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        bid = event.button.id or ""
        if bid.startswith("traj-"):
            self._on_traj_button(event)

    def action_next_frame(self) -> None:
        self._switch_frame(self._frame_idx + 1)

    def action_prev_frame(self) -> None:
        self._switch_frame(self._frame_idx - 1)

    def action_first_frame(self) -> None:
        self._switch_frame(0)

    def action_last_frame(self) -> None:
        if self._frames:
            self._switch_frame(len(self._frames) - 1)

    def action_toggle_play(self) -> None:
        if self._frames and len(self._frames) > 1:
            self._toggle_play()

    def action_toggle_bonds(self) -> None:
        view = self.query_one(MoleculeView)
        view.show_bonds = not view.show_bonds
        view._invalidate_cache()

    def action_toggle_bg(self) -> None:
        view = self.query_one(MoleculeView)
        view.dark_bg = not view.dark_bg
        self.theme = "textual-dark" if view.dark_bg else "textual-light"
        view._invalidate_cache()

    def action_reset_view(self) -> None:
        view = self.query_one(MoleculeView)
        view.rot_matrix = rotation_matrix(-0.2, -0.5, 0.0)
        view.pan_x = 0.0
        view.pan_y = 0.0
        view.pan_mode = False
        view.highlighted_atoms = set()
        if view.molecule:
            view.camera_distance = max(4.0, view.molecule.radius() * 3.0)
        view._invalidate_cache()

    def _frame_to_tab(self) -> None:
        """Open the current trajectory frame as a standalone structure in a new tab."""
        if not self._frames or not self.molecule or len(self.molecule.atoms) == 0:
            self.notify("No trajectory frame to open", timeout=2)
            return
        import copy
        mol = copy.deepcopy(self.molecule)
        stem = Path(self.filepath).stem if self.filepath else "frame"
        name = f"{stem}_frame{self._frame_idx + 1}"
        self.open_in_new_tab(mol, name)
        self.notify(f"Opened frame {self._frame_idx + 1} as {name}", timeout=2)

    def action_export_png(self) -> None:
        view = self.query_one(MoleculeView)
        if view.molecule is None:
            return
        self.notify("Exporting PNG...", timeout=2)
        self._track_task(self._export_png_async())

    async def _calc_void_volume(self) -> None:
        mol = self.molecule
        self._show_progress("Calculating void volume (geometric)...")
        result_geo = await asyncio.to_thread(void_volume, mol, 100000, 0.0)
        self.query_one("#status-label", Label).update("Calculating void volume (He probe)...")
        self._update_progress(50)
        result_he = await asyncio.to_thread(void_volume, mol, 100000, 1.14)
        self._hide_progress()

        gvv = result_geo.get("void_volume", 0)
        gtv = result_geo.get("total_volume", 0)
        gvf = result_geo.get("void_fraction", 0) * 100
        hvv = result_he.get("void_volume", 0)
        htv = result_he.get("total_volume", 0)
        hvf = result_he.get("void_fraction", 0) * 100
        lines = [
            "Void Volume Analysis (100,000 MC samples)",
            "",
            "Geometric (no probe):",
            f"  Void volume:   {gvv:.2f} A^3",
            f"  Total volume:  {gtv:.2f} A^3",
            f"  Void fraction: {gvf:.1f}%",
            "",
            "He probe (r = 1.14 A):",
            f"  Void volume:   {hvv:.2f} A^3",
            f"  Total volume:  {htv:.2f} A^3",
            f"  Void fraction: {hvf:.1f}%",
        ]
        self.push_screen(InfoModal("Void Volume", "\n".join(lines)))

    async def _calc_surface_area(self) -> None:
        mol = self.molecule
        self._show_progress("Calculating surface area (N2 probe)...")
        result = await asyncio.to_thread(surface_area, mol, 100000, 1.4)
        self._hide_progress()
        sa = result.get("surface_area", 0)
        apv = result.get("area_per_volume", 0)
        lines = [
            "Solvent-Accessible Surface Area",
            "  Probe radius: 1.4 A (N2-sized)",
            f"  Surface area: {sa:.1f} A^2",
            f"  Specific area: {apv:.1f} m^2/g",
        ]
        self.push_screen(InfoModal("Surface Area", "\n".join(lines)))

    async def _calc_pore_size(self) -> None:
        mol = self.molecule
        self._show_progress("Computing pore size distribution...")
        bins, hist = await asyncio.to_thread(pore_size_distribution, mol, 50000, 50, 20.0)
        self._hide_progress()
        if len(bins) == 0:
            self.notify("No PBC for pore size calculation", timeout=3)
            return
        # Find dominant pore size
        peak = bins[hist.argmax()] if hist.max() > 0 else 0
        # Show as a plot
        class PsdScreen(ModalScreen[None]):
            BINDINGS = [Binding("escape", "close", ""), Binding("q", "close", "")]
            DEFAULT_CSS = "PsdScreen > Vertical { width: 100%; height: 100%; background: $surface; padding: 0 1; }"
            def __init__(self, b, h, p):
                super().__init__()
                self._b, self._h, self._p = b, h, p
            def compose(self_inner) -> ComposeResult:
                with Vertical():
                    yield PlotWidget(id="psd-plot")
                    yield Label(f"Dominant pore diameter: {self_inner._p * 2:.2f} A  |  press q to close", id="psd-status")
            def on_mount(self_inner) -> None:
                plot = self_inner.query_one("#psd-plot", PlotWidget)
                y_max = float(self_inner._h.max()) * 1.1 if self_inner._h.max() > 0 else 1.0
                plot.set_data(self_inner._b * 2, self_inner._h, title="Pore Size Distribution",
                              x_label="pore diameter (A)", y_min=0.0, y_max=y_max)
            def action_close(self_inner) -> None:
                self_inner.dismiss(None)
        self.push_screen(PsdScreen(bins, hist, peak))

    async def _export_png_async(self) -> None:
        view = self.query_one(MoleculeView)
        if view.molecule is None:
            return

        mol = view.molecule
        if not view.show_bonds:
            mol = Molecule(atoms=mol.atoms, bonds=[], pbc=mol.pbc)

        export_w, export_h = 1600, 1200
        bg = (0, 0, 0) if view.dark_bg else (255, 255, 255)

        pixels, _ = await asyncio.to_thread(
            render_scene, export_w, export_h, mol, view.rot_matrix,
            view.camera_distance, bg_color=bg, pbc=view.pbc, ssaa=2,
            pan=(view.pan_x, view.pan_y),
            licorice=view.licorice, vdw=view.vdw, ribbon=view.ribbon,
            color_mode=view.color_mode, color_by_charge=view.color_by_charge,
            ambient=0.31, diffuse=0.72, specular=0.42, shininess=96.0,
            atom_scale=view.atom_scale, bond_radius=view.bond_radius,
            show_polyhedra=view.show_polyhedra,
            centroid_override=view._locked_centroid,
            density_positions=view._density_positions,
            density_values=view._density_values,
        )

        try:
            from PIL import Image
        except ImportError:
            self.notify("Pillow required: pip install Pillow", timeout=3)
            return

        stem = Path(self.filepath).stem
        out_path = Path(self.filepath).parent / f"{stem}.png"
        img = Image.fromarray(pixels)
        await asyncio.to_thread(img.save, str(out_path))
        self.notify(f"Saved {out_path}", timeout=3)

    async def _export_gif_async(self) -> None:
        view = self.query_one(MoleculeView)
        if view.molecule is None:
            return
        try:
            from PIL import Image
        except ImportError:
            self.notify("Pillow required: pip install Pillow", timeout=3)
            return

        mol = view.molecule
        if not view.show_bonds:
            mol = Molecule(atoms=mol.atoms, bonds=[], pbc=mol.pbc)

        export_w, export_h = 800, 600
        bg = (0, 0, 0) if view.dark_bg else (255, 255, 255)
        n_frames_gif = 36
        self._show_progress("Rendering GIF...", total=n_frames_gif)

        frames_img = []
        base_rot = view.rot_matrix.copy()
        for i in range(n_frames_gif):
            angle = 2 * np.pi * i / n_frames_gif
            rot = rotation_matrix(0, angle, 0) @ base_rot
            pixels, _ = await asyncio.to_thread(
                render_scene, export_w, export_h, mol, rot,
                view.camera_distance, bg_color=bg, pbc=view.pbc, ssaa=2,
                pan=(view.pan_x, view.pan_y),
                licorice=view.licorice, vdw=view.vdw, ribbon=view.ribbon,
                color_mode=view.color_mode, color_by_charge=view.color_by_charge,
                ambient=0.31, diffuse=0.72, specular=0.42, shininess=96.0,
                atom_scale=view.atom_scale, bond_radius=view.bond_radius,
                show_polyhedra=view.show_polyhedra,
                centroid_override=view._locked_centroid,
                density_positions=view._density_positions,
                density_values=view._density_values,
            )
            frames_img.append(Image.fromarray(pixels))
            self._update_progress(i + 1)

        self._hide_progress()
        stem = Path(self.filepath).stem
        out_path = Path(self.filepath).parent / (stem + ".gif")
        await asyncio.to_thread(
            frames_img[0].save,
            str(out_path),
            save_all=True,
            append_images=frames_img[1:],
            duration=80,
            loop=0,
        )
        self.notify(f"Saved {out_path} ({n_frames_gif} frames)", timeout=3)

    # --- Panel management ---

    def _active_panel_table(self) -> DataTable | None:
        ops = self.query_one(OperationsPanel)
        if ops.has_class("visible"):
            for dt in ops.query(DataTable):
                return dt
        geom = self.query_one(GeometryPanel)
        if geom.has_class("visible"):
            tabs = geom.query_one(TabbedContent)
            pane = tabs.get_pane(tabs.active)
            for dt in pane.query(DataTable):
                return dt
        return None

    def action_panel_next(self) -> None:
        if self.query_one(VisualPanel).has_class("visible"):
            self.screen.focus_next()
            return
        dt = self._active_panel_table()
        if dt is not None:
            dt.action_cursor_down()

    def action_panel_prev(self) -> None:
        if self.query_one(VisualPanel).has_class("visible"):
            self.screen.focus_previous()
            return
        dt = self._active_panel_table()
        if dt is not None:
            dt.action_cursor_up()

    def action_close_panel(self) -> None:
        geom = self.query_one(GeometryPanel)
        vis = self.query_one(VisualPanel)
        if geom.has_class("visible") or vis.has_class("visible"):
            self._close_panels()
            view = self.query_one(MoleculeView)
            view._invalidate_cache()
            view.focus()

    def _close_panels(self) -> None:
        """Close right-side panels. Operations panel is pinned and stays open."""
        view = self.query_one(MoleculeView)
        geom = self.query_one(GeometryPanel)
        vis = self.query_one(VisualPanel)
        if geom.has_class("visible"):
            geom.remove_class("visible")
            view.highlighted_atoms = set()
            view._stop_blink()
        if vis.has_class("visible"):
            vis.remove_class("visible")

    def action_toggle_operations(self) -> None:
        ops = self.query_one(OperationsPanel)
        view = self.query_one(MoleculeView)
        if ops.has_class("visible"):
            ops.remove_class("visible")
            view.focus()
        else:
            ops.add_class("visible")
            for dt in ops.query(DataTable):
                dt.focus()
                break
        view._invalidate_cache()

    def action_toggle_geometry(self) -> None:
        panel = self.query_one(GeometryPanel)
        was_visible = panel.has_class("visible")
        self._close_panels()
        view = self.query_one(MoleculeView)
        if not was_visible:
            panel.add_class("visible")
            # Populate active tab now that the panel is being shown
            panel._populate_active_tab()
            for dt in panel.query(DataTable):
                dt.focus()
                panel._emit_current_highlight(dt)
                break
        else:
            view.focus()
        view._invalidate_cache()

    def action_toggle_visual(self) -> None:
        vis = self.query_one(VisualPanel)
        was_visible = vis.has_class("visible")
        self._close_panels()
        view = self.query_one(MoleculeView)
        if not was_visible:
            vis.set_state(
                licorice=view.licorice, vdw=view.vdw, ribbon=view.ribbon,
                ambient=view.ambient, diffuse=view.diffuse,
                specular=view.specular, shininess=view.shininess,
                atom_scale=view.atom_scale, bond_radius=view.bond_radius,
            )
            vis.add_class("visible")
            for child in vis.query("*"):
                if child.can_focus:
                    child.focus()
                    break
        else:
            view.focus()
        view._invalidate_cache()

    # --- Event handlers ---

    def on_visual_panel_color_mode_changed(self, event: VisualPanel.ColorModeChanged) -> None:
        view = self.query_one(MoleculeView)
        view.color_mode = event.mode
        view.color_by_charge = (event.mode == "charge")
        view._invalidate_cache()

    def on_visual_panel_toggle_changed(self, event: VisualPanel.ToggleChanged) -> None:
        view = self.query_one(MoleculeView)
        if event.toggle == "bonds":
            view.show_bonds = event.value
        elif event.toggle == "atom_numbers":
            view.show_atom_numbers = event.value
        elif event.toggle == "charge_colors":
            view.color_by_charge = event.value
        elif event.toggle == "hide_water":
            view.hide_water = event.value
        elif event.toggle == "polyhedra":
            view.show_polyhedra = event.value
        view._invalidate_cache()

    def on_visual_panel_theme_changed(self, event: VisualPanel.ThemeChanged) -> None:
        view = self.query_one(MoleculeView)
        view.dark_bg = event.dark
        self.theme = "textual-dark" if event.dark else "textual-light"
        view._invalidate_cache()

    def on_visual_panel_style_changed(self, event: VisualPanel.StyleChanged) -> None:
        view = self.query_one(MoleculeView)
        view.licorice = event.licorice
        view.vdw = event.vdw
        view.ribbon = event.ribbon
        view.bond_radius = 0.15 if event.licorice else 0.08
        vis = self.query_one(VisualPanel)
        if vis.has_class("visible"):
            vis.set_state(
                licorice=view.licorice, vdw=view.vdw, ribbon=view.ribbon,
                ambient=view.ambient, diffuse=view.diffuse,
                specular=view.specular, shininess=view.shininess,
                atom_scale=view.atom_scale, bond_radius=view.bond_radius,
            )
        view._invalidate_cache()

    def on_visual_panel_lighting_changed(self, event: VisualPanel.LightingChanged) -> None:
        view = self.query_one(MoleculeView)
        view.ambient = event.ambient
        view.diffuse = event.diffuse
        view.specular = event.specular
        view.shininess = event.shininess
        view._invalidate_cache()

    def on_visual_panel_size_changed(self, event: VisualPanel.SizeChanged) -> None:
        view = self.query_one(MoleculeView)
        view.atom_scale = event.atom_scale
        view.bond_radius = event.bond_radius
        view._invalidate_cache()

    def on_menu_bar_menu_action(self, event: MenuBar.MenuAction) -> None:
        menu_map = {
            "file": (MENU_FILE, 0),
            "edit": (MENU_EDIT, 6),
            "view": (MENU_VIEW, 12),
            "analysis": (MENU_ANALYSIS, 20),
        }
        if event.action == "undo":
            self._do_undo()
        elif event.action == "redo":
            self._do_redo()
        elif event.action in menu_map:
            items, offset = menu_map[event.action]
            self.push_screen(
                MenuDropdown(items, offset=offset),
                callback=self._handle_menu_action,
            )

    def _handle_menu_action(self, action: str) -> None:
        if not action:
            return
        # Map menu actions to existing methods
        action_map = {
            "quit": lambda: self.exit(),
            "frame_to_tab": lambda: self._frame_to_tab(),
            "export_png": lambda: self.action_export_png(),
            "export_gif": lambda: self._track_task(self._export_gif_async()),
            "toggle_bonds": lambda: self.action_toggle_bonds(),
            "toggle_style": lambda: self.action_toggle_style(),
            "toggle_bg": lambda: self.action_toggle_bg(),
            "toggle_atom_numbers": lambda: self.action_toggle_atom_numbers(),
            "toggle_charge_colors": lambda: self.action_toggle_charge_colors(),
            "cycle_color_mode": lambda: self._cycle_color_mode(),
            "toggle_hide_water": lambda: self.action_toggle_hide_water(),
            "toggle_measure": lambda: self.action_toggle_measure(),
            "clear_measurements": lambda: self._clear_measurements(),
            "toggle_polyhedra": lambda: self.action_toggle_polyhedra(),
            "toggle_orthographic": lambda: self.action_toggle_orthographic(),
            "toggle_operations": lambda: self.action_toggle_operations(),
            "toggle_geometry": lambda: self.action_toggle_geometry(),
            "toggle_visual": lambda: self.action_toggle_visual(),
            "reset_view": lambda: self.action_reset_view(),
        }
        if action in action_map:
            action_map[action]()
        else:
            # Dispatch to operations panel handler
            self.on_operations_panel_run_command(OperationsPanel.RunCommand(action))

    def on_molecule_view_atom_picked(self, event: MoleculeView.AtomPicked) -> None:
        view = self.query_one(MoleculeView)
        mol = self.molecule
        idx = event.atom_index
        if not mol or not mol.atoms or idx >= len(mol.atoms):
            return
        atom = mol.atoms[idx]

        if view.measure_mode:
            # Measure mode: clicks accumulate atoms and pin labels
            picked = view._picked_atoms
            if idx not in picked:
                picked.append(idx)
            if len(picked) > 4:
                picked.pop(0)
            view.highlighted_atoms = set(picked)
            view._start_blink()

            if len(picked) == 1:
                a = mol.atoms[picked[0]]
                self.notify(f"Atom {picked[0]+1}: {a.element.symbol} — click more to measure", timeout=3)
            elif len(picked) == 2:
                i, j = picked
                if mol.pbc:
                    d = mol.pbc.min_image(mol.atoms[i].x - mol.atoms[j].x)
                else:
                    d = float(np.linalg.norm(mol.atoms[i].x - mol.atoms[j].x))
                label = f"{d:.3f}A"
                view._measurements.append(([i, j], label))
                self.notify(f"Distance {i+1}-{j+1} = {d:.4f} A (pinned)", timeout=3)
                view._picked_atoms = []  # reset for next measurement
                view.highlighted_atoms = set()
                view._stop_blink()
            elif len(picked) == 3:
                i, j, k = picked
                if mol.pbc:
                    v1 = mol.pbc.wrap(mol.atoms[i].x - mol.atoms[j].x)
                    v2 = mol.pbc.wrap(mol.atoms[k].x - mol.atoms[j].x)
                else:
                    v1 = mol.atoms[i].x - mol.atoms[j].x
                    v2 = mol.atoms[k].x - mol.atoms[j].x
                cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
                angle = float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
                label = f"{angle:.1f}deg"
                view._measurements.append(([i, j, k], label))
                self.notify(f"Angle {i+1}-{j+1}-{k+1} = {angle:.2f} deg (pinned)", timeout=3)
                view._picked_atoms = []
                view.highlighted_atoms = set()
                view._stop_blink()
        elif event.shift:
            # Legacy shift+click measurement (toast only)
            picked = view._picked_atoms
            if idx not in picked:
                picked.append(idx)
            if len(picked) > 3:
                picked.pop(0)
            view.highlighted_atoms = set(picked)
            view._start_blink()

            if len(picked) == 1:
                a = mol.atoms[picked[0]]
                self.notify(f"Atom {picked[0]+1}: {a.element.symbol} ({a.name})", timeout=3)
            elif len(picked) == 2:
                i, j = picked
                if mol.pbc:
                    d = mol.pbc.min_image(mol.atoms[i].x - mol.atoms[j].x)
                else:
                    d = float(np.linalg.norm(mol.atoms[i].x - mol.atoms[j].x))
                self.notify(
                    f"Distance {i+1}:{mol.atoms[i].element.symbol} - "
                    f"{j+1}:{mol.atoms[j].element.symbol} = {d:.4f} A", timeout=5,
                )
            elif len(picked) == 3:
                i, j, k = picked
                if mol.pbc:
                    v1 = mol.pbc.wrap(mol.atoms[i].x - mol.atoms[j].x)
                    v2 = mol.pbc.wrap(mol.atoms[k].x - mol.atoms[j].x)
                else:
                    v1 = mol.atoms[i].x - mol.atoms[j].x
                    v2 = mol.atoms[k].x - mol.atoms[j].x
                cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
                angle = float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
                self.notify(
                    f"Angle {i+1}-{j+1}-{k+1} = {angle:.2f} deg", timeout=5,
                )
        else:
            # Single click: show atom info, highlight it
            view._picked_atoms = [idx]
            view.highlighted_atoms = {idx}
            view._start_blink()
            charge_str = f"  q={atom.charge:.4f}" if abs(atom.charge) > 1e-10 else ""
            self.notify(
                f"Atom {idx+1}: {atom.element.symbol} ({atom.name})"
                f"  [{atom.x[0]:.3f}, {atom.x[1]:.3f}, {atom.x[2]:.3f}]{charge_str}",
                timeout=5,
            )
        view._invalidate_cache()

    def on_geometry_panel_highlight_atoms(self, event: GeometryPanel.HighlightAtoms) -> None:
        view = self.query_one(MoleculeView)
        view.highlighted_atoms = set(event.atom_indices)
        if view.highlighted_atoms:
            view._start_blink()
        else:
            view._stop_blink()
        view._invalidate_cache()

    def on_geometry_panel_edit_coord(self, event: GeometryPanel.EditCoord) -> None:
        idx = event.atom_index
        if not self.molecule or idx >= len(self.molecule.atoms):
            return
        a = self.molecule.atoms[idx]
        self.push_screen(
            InputModal(
                f"Edit coordinates for atom {idx + 1} ({a.element.symbol}):",
                default=f"{event.x:.6f} {event.y:.6f} {event.z:.6f}",
                hint="x y z",
            ),
            callback=lambda val: self._do_edit_coord(idx, val),
        )

    def _do_edit_coord(self, idx: int, value: str) -> None:
        if not value.strip():
            return
        try:
            import numpy as np
            parts = value.strip().split()
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            self.molecule.atoms[idx].x = np.array([x, y, z])
            self._track_task(self._refresh_molecule(f"Updated atom {idx + 1} position"))
        except (ValueError, IndexError):
            self.notify("Invalid input. Need: x y z", timeout=3)

    def on_operations_panel_run_command(self, event: OperationsPanel.RunCommand) -> None:

        mol = self.molecule
        cmd = event.command
        view = self.query_one(MoleculeView)

        if cmd == "info":
            info = system_info(mol, self.filepath)
            fname = info.get("filename", "")
            n_at = info.get("n_atoms", 0)
            lines = [f"Filename:  {fname}", f"Atoms:     {n_at}"]
            if "volume" in info:
                c = info.get("cell", {})
                ca, cb, cc = c.get("a", 0), c.get("b", 0), c.get("c", 0)
                al, be, ga = c.get("alpha", 0), c.get("beta", 0), c.get("gamma", 0)
                vol = info.get("volume", 0)
                dens = info.get("density", 0)
                lines.append(f"Cell:      {ca:.4f}  {cb:.4f}  {cc:.4f}")
                lines.append(f"Angles:    {al:.2f}  {be:.2f}  {ga:.2f}")
                lines.append(f"Volume:    {vol:.2f} A^3")
                lines.append(f"Density:   {dens:.4f} g/cm^3")
                lines.append("")
                for row in info.get("basis_matrix", []):
                    lines.append(f"  {row[0]:>14.8f} {row[1]:>14.8f} {row[2]:>14.8f}")
            fu = formula_unit(mol)
            reduced = formula_unit_reduced(mol)
            lines.append("")
            lines.append("Total:     " + "  ".join(f"{k}{v}" for k, v in fu.items()))
            lines.append("Formula:   " + "  ".join(f"{k}{v}" for k, v in reduced.items()))
            self.push_screen(InfoModal("System Info", "\n".join(lines)))
            return

        elif cmd == "rdf":
            if len(mol.atoms) == 0:
                self.notify("No atoms loaded yet", timeout=2)
                return
            self.push_screen(RdfScreen(mol, frames=self._frames))
            return

        elif cmd == "coordination":
            if len(mol.atoms) == 0:
                self.notify("No atoms loaded", timeout=2)
                return
            self.push_screen(CoordinationScreen(mol, frames=self._frames))
            return

        elif cmd == "hbonds":
            if len(mol.atoms) == 0:
                self.notify("No atoms loaded", timeout=2)
                return
            self.push_screen(HbondScreen(mol, self._frames))
            return

        elif cmd == "msd":
            if not self._frames or len(self._frames) < 3:
                self.notify("Need a trajectory with 3+ frames for MSD", timeout=3)
                return
            self.push_screen(MsdScreen(mol, self._frames))
            return

        elif cmd == "rmsd":
            if not self._frames or len(self._frames) < 3:
                self.notify("Need a trajectory with 3+ frames for RMSD", timeout=3)
                return
            self.push_screen(RmsdScreen(mol, self._frames))
            return

        elif cmd == "gyration":
            if not self._frames or len(self._frames) < 2:
                self.notify("Need a trajectory for Rg", timeout=3)
                return
            self.push_screen(GyrationScreen(mol, self._frames))
            return

        elif cmd == "density":
            if len(mol.atoms) == 0:
                self.notify("No atoms loaded", timeout=2)
                return
            self.push_screen(DensityScreen(mol, self._frames))
            return

        elif cmd == "pxrd":
            if len(mol.atoms) == 0 or mol.pbc is None:
                self.notify("Need a crystal structure with PBC for PXRD", timeout=3)
                return
            self.push_screen(PxrdScreen(mol))
            return

        elif cmd == "density3d":
            if not self._frames or len(self._frames) < 2:
                self.notify("Need a trajectory with PBC for 3D density", timeout=3)
                return
            if mol.pbc is None:
                self.notify("Need PBC for 3D density map", timeout=3)
                return
            self._show_density3d()
            return

        elif cmd in ("input_generator", "isotherm"):
            # "isotherm" kept as legacy command code for any saved keymaps;
            # both open the unified Input Generator (defaults to MPMC isotherm).
            self.push_screen(InputGeneratorScreen(mol, self.filepath))
            return

        elif cmd == "isotherm_plot":
            self.push_screen(
                FileSaveModal("Isotherm results directory:", default="isotherm", button_label="Open"),
                callback=self._plot_isotherm,
            )
            return

        elif cmd == "db_search":
            self.push_screen(DatabaseSearchScreen())
            return

        elif cmd == "energy_plot":
            self._show_energy_plot()
            return

        elif cmd == "undo":
            self._do_undo()
            return

        elif cmd == "void_vol":
            self.notify("Calculating void volume (100k samples)...", timeout=2)
            self._track_task(self._calc_void_volume())
            return

        elif cmd == "surface_area":
            if mol.pbc is None:
                self.notify("Need PBC for surface area", timeout=3)
                return
            self._track_task(self._calc_surface_area())
            return

        elif cmd == "pore_size":
            if mol.pbc is None:
                self.notify("Need PBC for pore size", timeout=3)
                return
            self._track_task(self._calc_pore_size())
            return

        elif cmd == "substitute":
            self.push_screen(
                InputModal("Substitute element:", default="Cu Zn", hint="old new (e.g. Cu Zn)"),
                callback=self._do_substitute,
            )
            return

        elif cmd == "reduce_cell":
            if mol.pbc is None:
                self.notify("Need PBC for supercell reduction", timeout=3)
                return
            self.push_screen(ReduceCellScreen(mol), callback=self._apply_reduced_cell)
            return

        elif cmd == "open_geom":
            self.action_toggle_geometry()
            return

        elif cmd == "wrap_center":
            self.push_screen(
                ConfirmModal(
                    "Wrap Atoms (centered)",
                    "Wrap all atoms into the unit cell centered at the origin.\n"
                    "Fractional coords mapped to (-0.5, 0.5).",
                ),
                callback=lambda ok: self._do_wrap(ok, forward=False),
            )
            return

        elif cmd == "wrap_forward":
            self.push_screen(
                ConfirmModal(
                    "Wrap Atoms (forward)",
                    "Wrap all atoms into the unit cell forward of the origin.\n"
                    "Fractional coords mapped to (0, 1).",
                ),
                callback=lambda ok: self._do_wrap(ok, forward=True),
            )
            return

        elif cmd == "sort":
            self.push_screen(
                ConfirmModal(
                    "Sort Atoms",
                    "Sort atoms by element within each molecule, then sort\n"
                    "molecules by size and composition. Bonds will be\n"
                    "recalculated.",
                ),
                callback=self._do_sort,
            )
            return

        elif cmd == "del_lone":
            self.push_screen(
                ConfirmModal(
                    "Delete Lone Atoms",
                    "Remove all atoms that have no neighbors within\n"
                    "VDW contact distance. This cannot be undone.",
                ),
                callback=self._do_del_lone,
            )
            return

        elif cmd == "extend":
            self.push_screen(
                ExtendAxisModal(mol.pbc),
                callback=self._do_extend,
            )
            return

        elif cmd == "edit_h":
            self.push_screen(
                InputModal(
                    "Edit H bond distances:",
                    default="O 1.0",
                    hint="element and distance, e.g. 'O 1.0' or 'C 1.09'",
                ),
                callback=self._do_edit_h,
            )
            return

        elif cmd == "update_cell":
            if mol.pbc is not None:
                pbc = mol.pbc
                default = f"{pbc.a:.4f} {pbc.b:.4f} {pbc.c:.4f} {pbc.alpha:.2f} {pbc.beta:.2f} {pbc.gamma:.2f}"
            else:
                default = "10.0 10.0 10.0 90.0 90.0 90.0"
            self.push_screen(
                InputModal(
                    "Unit cell parameters:",
                    default=default,
                    hint="a b c alpha beta gamma  •  type 'none' to remove the box",
                ),
                callback=self._do_update_cell,
            )
            return

        elif cmd == "qeq_charges":
            self.push_screen(
                ConfirmModal(
                    "Generate QEq Charges",
                    "Compute Rappe-Goddard charge equilibration\n"
                    "charges for all atoms? This will overwrite\n"
                    "any existing charges.\n\n"
                    f"{len(mol.atoms)} atoms",
                ),
                callback=self._do_qeq_charges,
            )
            return

        elif cmd == "write_xyz":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / f"{stem}_out.xyz")
            self.push_screen(
                FileSaveModal("Save XYZ as:", default=default),
                callback=self._do_write_xyz,
            )
            return

        elif cmd == "write_pdb":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / f"{stem}_out.pdb")
            self.push_screen(
                FileSaveModal("Save standard PDB as:", default=default),
                callback=self._do_write_pdb,
            )
            return

        elif cmd == "write_cif":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / f"{stem}_out.cif")
            self.push_screen(FileSaveModal("Save CIF as:", default=default), callback=self._do_write_cif)
            return

        elif cmd == "write_poscar":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / "POSCAR")
            self.push_screen(FileSaveModal("Save POSCAR as:", default=default), callback=self._do_write_poscar)
            return

        elif cmd == "write_lammps":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / (stem + ".data"))
            self.push_screen(FileSaveModal("Save LAMMPS data as:", default=default), callback=self._do_write_lammps)
            return

        elif cmd == "open_file":
            self.push_screen(
                FileSaveModal(
                    "Open structure file:",
                    default=self.filepath or str(Path.cwd()) + "/",
                    button_label="Open",
                ),
                callback=self._do_open_file,
            )
            return

        elif cmd == "frame_to_tab":
            self._frame_to_tab()
            return

        elif cmd == "export_png":
            self.action_export_png()
            return

        elif cmd == "export_gif":
            self._track_task(self._export_gif_async())
            return

        elif cmd == "write_mpmc":
            self._mpmc_state = {"write_charges": False, "write_ff": False}
            # Check if atoms already have charges
            has_charges = any(abs(a.charge) > 1e-10 for a in mol.atoms)
            if has_charges:
                total_q = sum(a.charge for a in mol.atoms)
                self.push_screen(
                    ConfirmModal(
                        "Charges Detected",
                        f"Atoms already have charges (total: {total_q:.4f} e).\n\n"
                        "Keep existing charges?\n"
                        "Yes = keep, No = choose new charges",
                        yes_no=True,
                    ),
                    callback=self._mpmc_step_keep_charges,
                )
            else:
                self._mpmc_step_charge_source()
            return

        elif cmd == "quit":
            # Direct exit: selecting Quit from the menu is already deliberate,
            # so it doesn't use the bare-'q' "press again" confirmation.
            self.exit()
            return

        view._invalidate_cache()

    async def _refresh_molecule(self, msg: str = "") -> None:
        """Redetect bonds and refresh the view and geometry panel."""
        mol = self.molecule
        self._show_progress("Detecting bonds...", total=100)
        await asyncio.sleep(0)

        def _bond_progress(frac: float) -> None:
            self.call_from_thread(self._update_progress, int(frac * 90))

        await asyncio.to_thread(mol.detect_bonds, progress_callback=_bond_progress)
        self._update_progress(80)
        view = self.query_one(MoleculeView)
        view.set_molecule(mol)
        self.query_one("#status-label", Label).update("Computing geometry tables...")
        await asyncio.sleep(0)
        # Pre-compute the heavy geometry data in a thread
        def _geo_progress(pct: int, label: str) -> None:
            self.call_from_thread(self._update_progress, 80 + pct // 7)
            self.call_from_thread(
                self.query_one("#status-label", Label).update, label
            )

        await asyncio.to_thread(self._precompute_geometry, mol, _geo_progress)
        self._update_progress(95)
        self.query_one("#status-label", Label).update("Populating tables...")
        await asyncio.sleep(0)
        panel = self.query_one(GeometryPanel)
        panel._molecule = mol
        if panel.is_mounted:
            await self._populate_tables_async(panel, mol)
        self._update_progress(100)
        self._hide_progress()
        if msg:
            self.notify(msg, timeout=2)
        view._invalidate_cache()

    def _do_wrap(self, confirmed: bool, forward: bool) -> None:
        if not confirmed:
            return
        self._track_task(self._do_wrap_async(forward))

    async def _do_wrap_async(self, forward: bool) -> None:
        self._save_undo()
        self._show_progress("Wrapping atoms...")
        await asyncio.sleep(0)
        await asyncio.to_thread(wrap_atoms, self.molecule, forward)
        label = "forward" if forward else "centered"
        await self._refresh_molecule(f"Wrapped atoms ({label})")

    def _do_sort(self, confirmed: bool) -> None:
        if not confirmed:
            return
        self._track_task(self._do_sort_async())

    async def _do_sort_async(self) -> None:
        self._save_undo()
        self._show_progress("Sorting atoms...")
        await asyncio.sleep(0)
        await asyncio.to_thread(sort_system, self.molecule, False)
        await self._refresh_molecule("Sorted atoms")

    def _do_del_lone(self, confirmed: bool) -> None:
        if not confirmed:
            return
        self._track_task(self._do_del_lone_async())

    async def _do_del_lone_async(self) -> None:
        self._save_undo()
        before = len(self.molecule.atoms)
        self._show_progress("Finding lone atoms...")
        await asyncio.sleep(0)
        await asyncio.to_thread(delete_lone_atoms, self.molecule)
        after = len(self.molecule.atoms)
        await self._refresh_molecule(f"Deleted {before - after} lone atoms")

    def _do_open_file(self, value: str | None) -> None:
        """Callback for the File > Open... picker."""
        if not value or not value.strip():
            return
        filepath = value.strip()
        # If the buffer is empty (no atoms, fresh launch), reuse the active
        # tab instead of opening a new one — avoids a leftover blank tab.
        if not self.molecule.atoms and not self.filepath:
            try:
                mol = read_file(filepath)
                self.molecule = mol
                self.filepath = filepath
                self._tabs[self._active_tab].molecule = mol
                self._tabs[self._active_tab].filepath = filepath
                view = self.query_one(MoleculeView)
                view.set_molecule(mol, keep_camera=False)
                view._invalidate_cache()
                panel = self.query_one(GeometryPanel)
                panel.set_molecule(mol)
                self._update_title()
                self._refresh_tab_bar()
                self.notify(f"Opened {Path(filepath).name}", timeout=2)
            except Exception as e:
                self.notify(f"Error opening file: {e}", timeout=4, severity="error")
            return
        # Otherwise open in a new tab.
        self._open_new_tab(filepath)

    def _do_extend(self, value: tuple) -> None:
        if not value or len(value) != 3:
            return
        na, nb, nc = value
        if na == 1 and nb == 1 and nc == 1:
            return  # no-op
        self._track_task(self._do_extend_async(na, nb, nc))

    async def _do_extend_async(self, na: int, nb: int, nc: int) -> None:
        self._save_undo()
        self._show_progress(f"Extending {na}x{nb}x{nc}...")
        await asyncio.sleep(0)

        def _extend_all(mol):
            # All three axes inside one thread — main loop only sees the
            # final state, never the intermediate (renderer doesn't render
            # the partial 2x1x1 / 2x2x1 frames).
            if na > 1:
                extend_axis(mol, 0, na - 1)
            if nb > 1:
                extend_axis(mol, 1, nb - 1)
            if nc > 1:
                extend_axis(mol, 2, nc - 1)
            return mol

        await asyncio.to_thread(_extend_all, self.molecule)
        await self._refresh_molecule(f"Extended to {na}x{nb}x{nc}")

    def _confirm_overwrite(self, filepath: str, callback) -> None:
        """If filepath exists, prompt to overwrite. Otherwise call callback directly."""
        if Path(filepath).exists():
            self.push_screen(
                ConfirmModal(
                    "Overwrite?",
                    f"{filepath}\nalready exists. Overwrite?",
                    yes_no=True,
                ),
                callback=lambda ok: callback(filepath) if ok else None,
            )
        else:
            callback(filepath)

    def _do_write_xyz(self, value: str) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), self._write_xyz_file)

    def _write_xyz_file(self, filepath: str) -> None:
        try:
            with open(filepath, "w") as f:
                write_xyz(self.molecule, f)
            self.notify(f"Wrote {filepath}", timeout=3)
        except OSError as e:
            self.notify(f"Error: {e}", timeout=3)

    def _do_write_pdb(self, value: str) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), self._write_pdb_file)

    def _write_pdb_file(self, filepath: str) -> None:
        try:
            with open(filepath, "w") as f:
                write_standard_pdb(self.molecule, f)
            self.notify(f"Wrote {filepath}", timeout=3)
        except OSError as e:
            self.notify(f"Error: {e}", timeout=3)

    def _do_write_cif(self, value: str) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), self._write_cif_file)

    def _write_cif_file(self, filepath: str) -> None:
        try:
            with open(filepath, "w") as f:
                write_cif(self.molecule, f)
            self.notify(f"Wrote {filepath}", timeout=3)
        except OSError as e:
            self.notify(f"Error: {e}", timeout=3)

    def _do_write_poscar(self, value: str) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), self._write_poscar_file)

    def _write_poscar_file(self, filepath: str) -> None:
        try:
            with open(filepath, "w") as f:
                write_poscar(self.molecule, f)
            self.notify(f"Wrote {filepath}", timeout=3)
        except OSError as e:
            self.notify(f"Error: {e}", timeout=3)

    def _do_write_lammps(self, value: str) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), self._write_lammps_file)

    def _write_lammps_file(self, filepath: str) -> None:
        try:
            with open(filepath, "w") as f:
                write_lammps_data(self.molecule, f)
            self.notify(f"Wrote {filepath}", timeout=3)
        except OSError as e:
            self.notify(f"Error: {e}", timeout=3)

    def _apply_reduced_cell(self, result: Molecule | None) -> None:
        if result is None:
            return
        self._save_undo()
        self.molecule = result
        self._tabs[self._active_tab].molecule = result
        view = self.query_one(MoleculeView)
        view.set_molecule(result, keep_camera=False)
        view._invalidate_cache()
        panel = self.query_one(GeometryPanel)
        panel.set_molecule(result)
        self._update_title()

    def _do_substitute(self, value: str) -> None:
        if not value.strip():
            return
        parts = value.strip().split()
        if len(parts) != 2:
            self.notify("Enter: old_element new_element (e.g. Cu Zn)", timeout=3)
            return
        old_el, new_el = parts[0], parts[1]
        new_elem = get_element(new_el)
        if new_elem.symbol == "X":
            self.notify(f"Unknown element: {new_el}", timeout=3)
            return
        self._save_undo()
        count = 0
        for atom in self.molecule.atoms:
            if atom.element.symbol == old_el:
                atom._element = new_elem
                atom.name = new_elem.symbol
                count += 1
        if count == 0:
            self.notify(f"No {old_el} atoms found", timeout=2)
            return
        view = self.query_one(MoleculeView)
        view.set_molecule(self.molecule, keep_camera=True)
        view._invalidate_cache()
        panel = self.query_one(GeometryPanel)
        panel.set_molecule(self.molecule)
        self.notify(f"Substituted {count} {old_el} -> {new_el}", timeout=3)

    def _do_edit_h(self, value: str) -> None:
        if not value.strip():
            return
        try:
            parts = value.strip().split()
            element = parts[0].strip().capitalize()
            distance = float(parts[1])
            self._track_task(self._do_edit_h_async(element, distance))
        except (ValueError, IndexError):
            self.notify("Invalid input. Use: O 1.0, C 1.09, etc.", timeout=3)

    async def _do_edit_h_async(self, element: str, distance: float) -> None:
        self._save_undo()
        self._show_progress(f"Editing {element}-H distances...")
        await asyncio.sleep(0)
        await asyncio.to_thread(edit_h_dist, self.molecule, element, distance)
        await self._refresh_molecule(f"Set {element}-H distances to {distance} A")

    def _do_qeq_charges(self, confirmed: bool) -> None:
        if not confirmed:
            return
        self._track_task(self._do_qeq_async())

    async def _do_qeq_async(self) -> None:
        mol = self.molecule
        self._show_progress("Computing QEq charges...", total=100)
        await asyncio.sleep(0)

        def _progress(frac: float) -> None:
            self.call_from_thread(self._update_progress, int(frac * 100))

        try:
            await asyncio.to_thread(apply_qeq_charges, mol, progress_callback=_progress)
        except Exception as e:
            self._hide_progress()
            self.notify(f"QEq failed: {e}", timeout=5)
            return
        self._hide_progress()

        # Show summary
        charges = np.array([a.charge for a in mol.atoms])
        lines = [
            f"QEq charges computed for {len(mol.atoms)} atoms",
            "",
            f"Total charge:  {charges.sum():.6f} e",
            f"Min charge:    {charges.min():.6f} e",
            f"Max charge:    {charges.max():.6f} e",
            f"Mean |charge|: {np.abs(charges).mean():.6f} e",
            "",
            "Per-element averages:",
        ]
        from collections import defaultdict
        by_el: dict[str, list[float]] = defaultdict(list)
        for a in mol.atoms:
            by_el[a.element.symbol].append(a.charge)
        for el in sorted(by_el):
            arr = by_el[el]
            lines.append(f"  {el:<3} ({len(arr):>4}):  avg = {np.mean(arr):>8.5f}  range = [{min(arr):.5f}, {max(arr):.5f}]")
        self.push_screen(InfoModal("QEq Charges", "\n".join(lines)))

    # MSD and Density are now launched as full screens (MsdScreen, DensityScreen)

    def _show_density3d(self) -> None:
        """Prompt for element and compute 3D density map overlay."""
        elements = sorted(set(a.element.symbol for a in self.molecule.atoms))
        items = [(el, f"3D density: {el}") for el in elements]
        self.push_screen(
            MenuDropdown(items, offset=20),
            callback=self._do_density3d,
        )

    def _do_density3d(self, element: str) -> None:
        if not element:
            return
        self._track_task(self._compute_density3d(element))

    async def _compute_density3d(self, element: str) -> None:
        self._show_progress(f"Computing 3D density for {element}...", total=100)
        _prog = [0.0]
        def _cb(f):
            _prog[0] = f
        timer = self.set_interval(0.1, lambda: self._update_progress(int(_prog[0] * 100)))
        await asyncio.sleep(0)

        positions, densities, pbc = await asyncio.to_thread(
            compute_density_3d, self._frames, element, n_bins=15, progress_callback=_cb,
        )
        timer.stop()
        self._hide_progress()

        if len(positions) == 0:
            self.notify(f"No {element} atoms found in trajectory", timeout=3)
            return

        view = self.query_one(MoleculeView)
        view._density_positions = positions
        view._density_values = densities
        view._density_element = element
        view._invalidate_cache()
        self.notify(f"3D density: {len(positions)} bins for {element} (toggle with View menu)", timeout=3)

    def _plot_isotherm(self, dirpath: str) -> None:
        if not dirpath.strip():
            return
        import os
        dirpath = dirpath.strip()
        if not os.path.isdir(dirpath):
            self.notify(f"Not a directory: {dirpath}", timeout=3)
            return
        results = parse_isotherm_results(dirpath)
        if not results:
            self.notify("No completed pressure points found", timeout=3)
            return
        pressures = np.array([r[0] for r in results])
        uptakes = np.array([r[1] for r in results])

        class IsothermPlotScreen(_AnalysisScreen):
            _prefix = "iso"
            DEFAULT_CSS = _ANALYSIS_CSS.format(cls="IsothermPlotScreen") + """
            IsothermPlotScreen #iso-status { height: 1; color: $text-muted; }
            """
            def __init__(self, p, u, path):
                ModalScreen.__init__(self)
                self._molecule = None
                self._frames = None
                self._elements = []
                self._compute_task = None
                self._p, self._u, self._path = p, u, path
            def compose(self_inner) -> ComposeResult:
                with Vertical():
                    with Horizontal():
                        yield Button("CSV", id="iso-csv")
                        yield Button("Close", id="iso-close")
                    yield PlotWidget(id="iso-plot")
                    yield Label("", id="iso-status")
            def on_mount(self_inner) -> None:
                plot = self_inner.query_one("#iso-plot", PlotWidget)
                y_max = float(self_inner._u.max()) * 1.1 if self_inner._u.max() > 0 else 1.0
                plot.set_data(self_inner._p, self_inner._u,
                              title="Adsorption Isotherm",
                              x_label="pressure (atm)", y_min=0.0, y_max=y_max)
                n = len(self_inner._p)
                self_inner.query_one("#iso-status", Label).update(
                    f"{n} pressure points from {self_inner._path}"
                )
            def on_button_pressed(self_inner, event: Button.Pressed) -> None:
                if event.button.id == "iso-close":
                    self_inner.dismiss(None)
                elif event.button.id == "iso-csv":
                    self_inner.app.push_screen(
                        FileSaveModal("Save isotherm CSV:", default="isotherm.csv"),
                        callback=lambda fp: self_inner._do_export_csv(
                            fp, "pressure_atm,avg_molecules", self_inner._p, self_inner._u),
                    )

        self.push_screen(IsothermPlotScreen(pressures, uptakes, dirpath))

    def _show_energy_plot(self) -> None:
        from pathlib import Path
        # Look for energy.dat in same directory as the loaded file
        base = Path(self.filepath).parent
        candidates = [base / "energy.dat", base / "energy.csv"]
        energy_file = None
        for c in candidates:
            if c.exists():
                energy_file = str(c)
                break
        if energy_file is None:
            # Ask user for file
            self.push_screen(
                FileSaveModal("Open energy file:", default=str(base / "energy.dat"), button_label="Open"),
                callback=self._load_energy_plot,
            )
        else:
            self._load_energy_plot(energy_file)

    def _load_energy_plot(self, filepath: str) -> None:
        if not filepath.strip():
            return
        data = read_energy_dat(filepath.strip())
        if not data:
            self.notify(f"Could not read {filepath}", timeout=3)
            return
        self.push_screen(EnergyPlotScreen(data, filepath.strip()))

    def _save_undo(self) -> None:
        import copy
        self._undo_stack.append(copy.deepcopy(self.molecule))
        if len(self._undo_stack) > 10:
            self._undo_stack.pop(0)
        self._redo_stack.clear()

    def _do_undo(self) -> None:
        if not self._undo_stack:
            self.notify("Nothing to undo", timeout=2)
            return
        import copy
        self._redo_stack.append(copy.deepcopy(self.molecule))
        self.molecule = self._undo_stack.pop()
        view = self.query_one(MoleculeView)
        view.set_molecule(self.molecule)
        self.query_one(GeometryPanel).set_molecule(self.molecule)
        view._invalidate_cache()
        self.notify("Undone", timeout=1)

    def _do_redo(self) -> None:
        if not self._redo_stack:
            self.notify("Nothing to redo", timeout=2)
            return
        import copy
        self._undo_stack.append(copy.deepcopy(self.molecule))
        self.molecule = self._redo_stack.pop()
        view = self.query_one(MoleculeView)
        view.set_molecule(self.molecule)
        self.query_one(GeometryPanel).set_molecule(self.molecule)
        view._invalidate_cache()
        self.notify("Redone", timeout=1)

    # Density profile is now launched as DensityScreen

    def _do_update_cell(self, value: str) -> None:
        if value is None or not value.strip():
            return
        text = value.strip().lower()
        # Sentinel words to remove the box entirely (turn off PBC)
        if text in ("none", "off", "remove", "delete", "no", "0"):
            if self.molecule.pbc is None:
                self.notify("This system already has no PBC", timeout=2)
                return
            self.molecule.pbc = None
            self._track_task(self._refresh_molecule("PBC removed"))
            return
        try:
            parts = value.strip().split()
            a, b, c = float(parts[0]), float(parts[1]), float(parts[2])
            alpha, beta, gamma = float(parts[3]), float(parts[4]), float(parts[5])
            if self.molecule.pbc is not None:
                self.molecule.pbc.update(a, b, c, alpha, beta, gamma)
            else:
                self.molecule.pbc = PBC(a, b, c, alpha, beta, gamma)
            self._track_task(
                self._refresh_molecule(f"Cell: {a:.3f} {b:.3f} {c:.3f} {alpha:.1f} {beta:.1f} {gamma:.1f}")
            )
        except (ValueError, IndexError):
            self.notify(
                "Invalid input. Need: a b c alpha beta gamma  (or 'none' to remove)",
                timeout=4,
            )

    def _mpmc_step_keep_charges(self, keep: bool) -> None:
        if keep:
            self._mpmc_state["write_charges"] = True
            self._mpmc_step_ask_ff()
        else:
            # Zero out charges first
            for a in self.molecule.atoms:
                a.charge = 0.0
            self._mpmc_step_charge_source()

    def _mpmc_step_charge_source(self) -> None:
        self.push_screen(
            ConfirmModal(
                "Charge Source",
                "How would you like to assign charges?\n\n"
                "Yes = Load from file (.resp or raw)\n"
                "No = Generate QEq charges automatically",
                yes_no=True,
            ),
            callback=self._mpmc_step_charge_choice,
        )

    def _mpmc_step_charge_choice(self, from_file: bool) -> None:
        if from_file:
            self._mpmc_step_charges_file_prompt()
        else:
            # Run QEq
            self._mpmc_state["write_charges"] = True
            self._track_task(self._mpmc_run_qeq())

    async def _mpmc_run_qeq(self) -> None:
        self._show_progress("Computing QEq charges...", total=100)
        await asyncio.sleep(0)

        def _progress(frac: float) -> None:
            self.call_from_thread(self._update_progress, int(frac * 100))

        try:
            await asyncio.to_thread(apply_qeq_charges, self.molecule, progress_callback=_progress)
            self._hide_progress()
            total_q = sum(a.charge for a in self.molecule.atoms)
            self.notify(f"QEq charges applied (total: {total_q:.4f} e)", timeout=2)
        except Exception as e:
            self._hide_progress()
            self.notify(f"QEq failed: {e}", timeout=5)
        self._mpmc_step_ask_ff()

    def _mpmc_step_charges_file_prompt(self) -> None:
        self._mpmc_state["write_charges"] = True
        self.push_screen(
            ChargesFileModal(),
            callback=self._mpmc_step_charges_file,
        )

    def _mpmc_step_charges_file(self, result: dict) -> None:
        if not result:
            return  # cancelled
        mol = self.molecule
        try:
            charges = read_charges_file(
                result["path"],
                skip_first=result.get("skip_first", 0),
                skip_last=result.get("skip_last", 0),
            )

            n = len(mol.atoms)
            if len(charges) == n:
                for ind, atom in enumerate(mol.atoms):
                    atom.charge = charges[ind]
            elif len(charges) > 0 and n % len(charges) == 0:
                for ind, atom in enumerate(mol.atoms):
                    atom.charge = charges[ind % len(charges)]
                self.notify(
                    f"Applied {len(charges)} charges recursively "
                    f"({n // len(charges)}x)", timeout=3,
                )
            else:
                self.notify(
                    f"Charge count ({len(charges)}) doesn't match "
                    f"atom count ({n})", timeout=5,
                )
                return
            self.notify(f"Loaded {len(charges)} charges", timeout=2)
        except (OSError, ValueError) as e:
            self.notify(f"Error reading charges: {e}", timeout=5)
            return
        self._mpmc_step_ask_ff()

    def _mpmc_step_ask_ff(self) -> None:
        self.push_screen(
            ForceFieldModal(self.molecule),
            callback=self._mpmc_step_ff_apply,
        )

    def _mpmc_step_ff_apply(self, ff_name: str | None) -> None:
        if ff_name is None:
            return  # Cancel: abort the MPMC export
        if not ff_name:
            self._mpmc_state["write_ff"] = False
        else:
            ff_map = {"OPLSAA": 0, "PHAHST": 1}
            ff_idx = ff_map.get(ff_name, 0)
            apply_ff_to_system(self.molecule.atoms, get_forcefield(ff_idx))
            label = ["OPLS-AA/UFF", "PHAHST"][ff_idx]
            self.notify(f"Applied {label} force field", timeout=2)
        self._mpmc_step_sorbate()

    def _mpmc_step_sorbate(self) -> None:
        self.push_screen(
            SorbateModal(),
            callback=self._mpmc_step_sorbate_result,
        )

    def _mpmc_step_sorbate_result(self, model_name: str | None) -> None:
        if model_name is None:
            return  # Cancel: abort the MPMC export
        self._mpmc_state["sorbate"] = model_name if model_name else None
        self._mpmc_step_filename()

    def _mpmc_step_filename(self) -> None:
        stem = Path(self.filepath).stem
        default = str(Path(self.filepath).parent / f"{stem}_out_mpmc.pdb")
        self.push_screen(
            FileSaveModal("Save MPMC PDB as:", default=default),
            callback=self._mpmc_step_write,
        )

    def _mpmc_step_write(self, value: str) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), self._write_mpmc_file)

    def _write_mpmc_file(self, filepath: str) -> None:


        sorbate_lines = None
        sorbate_name = self._mpmc_state.get("sorbate")
        if sorbate_name:
            mol = self.molecule
            center = _find_nonoverlap_position(mol)
            n_atoms = len(mol.atoms)
            sorbate_lines = format_sorbate_pqr(
                sorbate_name, center[0], center[1], center[2],
                mol_id=2, start_atom_id=n_atoms + 1,
            )

        try:
            write_mpmc_pdb(
                self.molecule, filepath,
                write_charges=self._mpmc_state["write_charges"],
                write_params=self._mpmc_state["write_ff"],
                sorbate_lines=sorbate_lines,
            )
            self.notify(f"Wrote {filepath}", timeout=3)
        except OSError as e:
            self.notify(f"Error: {e}", timeout=3)

# ======================================================================
# Module: menu
# ======================================================================
"""Classic text menu interface (from original pdb_wizard)."""





def _list_coords(mol: Molecule) -> None:
    for atom in mol.atoms:
        print(f"{atom.element.symbol} {atom.x}")


def _vmd_preview(mol: Molecule) -> None:
    write_mpmc_pdb(mol, "pdb_wizard.tmp.pdb")
    os.system("vmd pdb_wizard.tmp.pdb")
    os.system("rm pdb_wizard.tmp.pdb")


def _prompt_cell() -> tuple[float, float, float, float, float, float]:
    while True:
        try:
            a = float(input("Enter cell information\na>     "))
            b = float(input("b>     "))
            c = float(input("c>     "))
            alpha = float(input("alpha> "))
            beta = float(input("beta>  "))
            gamma = float(input("gamma> "))
            return a, b, c, alpha, beta, gamma
        except ValueError:
            print("!!! Error converting input to float !!!\n")


def _menu_update_pbc(pbc: PBC) -> PBC:
    pbc.update(*_prompt_cell())
    return pbc


def _prompt_axis_times(pbc: PBC) -> tuple[int, int] | None:
    """Prompt for an axis (0/1/2 or x/y/z) and a replication count, showing the
    current cell. Returns (axis, times), or None if the user cancels."""
    print(
        f"\nCurrent cell:\n{round(pbc.a, 3):>7}  {round(pbc.b, 3):7}  {round(pbc.c, 3):7} "
        f"{round(pbc.alpha, 2):>6} {round(pbc.beta, 2):>6} {round(pbc.gamma, 2):>6}\n"
    )
    for row in pbc.basis_matrix:
        print(f"{row[0]:20.14f} {row[1]:20.14f} {row[2]:20.14f}")
    while True:
        try:
            axis_in = input(
                "\nWhat axis would you like to extend? 0, 1, 2 or x, y, z or q(uit)\n\n> "
            )
            if axis_in.lower() in ("q", "quit"):
                return None
            axis_map = {"x": 0, "y": 1, "z": 2}
            # Don't use dict.get(..., int(axis_in)): the default is evaluated
            # eagerly, so a letter axis would raise before the lookup.
            if axis_in.lower() in axis_map:
                axis = axis_map[axis_in.lower()]
            else:
                axis = int(axis_in)
            if axis < 0 or axis > 2:
                raise ValueError
            times = int(input("\nHow many times would you like to extend it?\n\n> "))
            if times < 1:
                raise ValueError
            return axis, times
        except ValueError:
            print("!!! Error converting input to int or x, y, z !!!")


def _prompt_frame(n_frames: int) -> int | None:
    """Prompt for a 1-based frame number, returning the 0-based index (or None
    if the user cancels)."""
    while True:
        try:
            raw = input(f"\nWhich frame? (1-{n_frames}, q to cancel)\n\n> ")
            if raw.lower() in ("q", "quit"):
                return None
            idx = int(raw) - 1
            if idx < 0 or idx >= n_frames:
                raise ValueError
            return idx
        except ValueError:
            print(f"!!! Error: enter a frame number from 1 to {n_frames} !!!")


def _menu_extend_axis(mol: Molecule) -> Molecule:
    if mol.pbc is None:
        print("No PBC data available")
        return mol
    res = _prompt_axis_times(mol.pbc)
    if res is None:
        return mol
    axis, times = res
    return extend_axis(mol, axis, times)


def _write_mpmc_options(mol: Molecule) -> None:
    if mol.pbc is None:
        return
    while True:
        try:
            wc = input(
                "\nWould you like to read in charges?\n"
                "('yes', 'y', 1 or 'no', 'n', 0)\n\n> "
            )
            write_charges = wc.lower() in ("yes", "y", "1")
            if wc.lower() not in ("yes", "y", "1", "no", "n", "0"):
                raise ValueError
            break
        except ValueError:
            print("!!! Error reading input !!!")

    if write_charges:
        while True:
            charges_filename = input(
                "\nEnter a resp file or a valid column of raw charges\n"
                "charges file name > "
            )
            try:
                charges = read_charges_file(charges_filename)
                n = len(mol.atoms)
                if len(charges) == n:
                    for ind, atom in enumerate(mol.atoms):
                        atom.charge = charges[ind]
                elif n % len(charges) == 0:
                    for ind, atom in enumerate(mol.atoms):
                        atom.charge = charges[ind % len(charges)]
                else:
                    raise ValueError
                break
            except (TypeError, ValueError, FileNotFoundError) as e:
                print(f"!!! Error reading charges: {e} !!!")

    while True:
        try:
            wf = input(
                "\nWould you like to automatically apply a forcefield?\n"
                "('yes', 'y', 1 or 'no', 'n', 0)\n\n> "
            )
            write_ff = wf.lower() in ("yes", "y", "1")
            if wf.lower() not in ("yes", "y", "1", "no", "n", "0"):
                raise ValueError
            break
        except ValueError:
            print("!!! Error reading input !!!")

    if write_ff:
        while True:
            try:
                ff_in = input(
                    "\nWhich force field?\n"
                    "valid answers are 'OPLSAA' (0) or 'PHAHST' (1)\n\n> "
                )
                ff_map = {"OPLSAA": 0, "PHAHST": 1}
                ff_idx = ff_map.get(ff_in, int(ff_in))
                if ff_idx not in (0, 1):
                    raise ValueError
                apply_ff_to_system(mol.atoms, get_forcefield(ff_idx))
                break
            except ValueError:
                print("!!! Error reading input !!!")

    filename = input("\noutput filename > ")
    write_mpmc_pdb(mol, filename, write_charges=write_charges, write_params=write_ff)


def _menu_geom_analysis(mol: Molecule) -> Molecule:
    while True:
        option = 0
        try:
            option = int(input(
                "\nWhat would you like to do?\n\n"
                "1 = list bonds\n"
                "2 = list close vdw contacts\n"
                "3 = list angles\n"
                "4 = list lone atoms\n"
                "5 = delete lone atoms\n"
                "6 = list coordinates\n"
                "7 = edit hydrogen bond distances\n"
                "8 = preview with VMD\n"
                "9 = back to main menu\n\n> "
            ))
        except ValueError:
            print("!!! Error converting input to int !!!")
        if option == 1:
            print("\nBonded atoms:\n")
            for msg in get_bonds_list(mol):
                print(msg)
        elif option == 2:
            print("\nClose contacts:\n")
            for msg in get_close_contacts(mol):
                print(msg)
        elif option == 3:
            print("\nAngles:\n")
            for msg in get_angles_list(mol):
                print(msg)
        elif option == 4:
            lone = get_lone_atoms(mol)
            if not lone:
                print("\nNo lone atoms found\n")
            else:
                print("\nLone atoms:\n")
                for atom in lone:
                    print(f"{atom.element.symbol:>3} {atom.id:>5} {atom.x}")
        elif option == 5:
            mol = delete_lone_atoms(mol)
        elif option == 6:
            _list_coords(mol)
        elif option == 7:
            el = input("\nLook for hydrogens bonded with which element?\n\n> ")
            dist = float(input(f"\nWhat distance shall {el}-H bonds be set to?\n\n> "))
            mol = edit_h_dist(mol, el, dist)
        elif option == 8:
            _vmd_preview(mol)
        elif option == 9:
            return mol
    return mol


def _main_loop_single(mol: Molecule, filename: str) -> None:
    mol = overlap_detector(mol)

    while True:
        print_info(mol, filename)
        try:
            option = int(input(
                "\nWhat would you like to do?\n\n"
                "1 = geometry analysis\n"
                "2 = extend along axis\n"
                "3 = wrap atoms from (0, 0, 0) to (1, 1, 1)\n"
                "4 = wrap atoms from (-1/2, -1/2, -1/2) to (1/2, 1/2, 1/2)\n"
                "5 = update cell dimensions\n"
                "6 = write .xyz\n"
                "7 = write MPMC .pdb\n"
                "8 = write standardized .pdb\n"
                "9 = quit\n\n> "
            ))
        except ValueError:
            print("!!! Error converting input to int !!!")
            continue
        if option == 1:
            mol = _menu_geom_analysis(mol)
        elif option == 2:
            mol = _menu_extend_axis(mol)
        elif option == 3:
            wrap_atoms(mol, forward=True)
            print("\nWrapped atoms forward of origin")
        elif option == 4:
            wrap_atoms(mol)
            print("\nWrapped atoms around origin")
        elif option == 5:
            if mol.pbc is not None:
                mol.pbc = _menu_update_pbc(mol.pbc)
        elif option == 6:
            out_filename = input("\noutput filename > ")
            with open(out_filename, "w") as out:
                write_xyz(mol, out)
            print(f"wrote {out_filename}")
        elif option == 7:
            _write_mpmc_options(mol)
        elif option == 8:
            out_filename = input("\noutput filename > ")
            with open(out_filename, "w") as out:
                write_standard_pdb(mol, out)
            print(f"wrote {out_filename}")
        elif option == 9:
            break
        else:
            print("\nInvalid option!")


def _main_loop_movie(molecules: list[Molecule], filename: str) -> None:
    while True:
        print_info_movie(len(molecules), filename)
        n = len(molecules)
        try:
            option = int(input(
                "\nWhat would you like to do?\n\n"
                "1 = geometry analysis (single frame)\n"
                "2 = extend along axis (all frames)\n"
                "3 = wrap atoms from (0, 0, 0) to (1, 1, 1) (all frames)\n"
                "4 = wrap atoms from (-1/2, -1/2, -1/2) to (1/2, 1/2, 1/2) (all frames)\n"
                "5 = update cell dimensions (all frames)\n"
                "6 = write .xyz (all frames)\n"
                "7 = write MPMC .pdb (single frame)\n"
                "8 = write standardized .pdb (all frames)\n"
                "9 = quit\n\n> "
            ))
        except ValueError:
            print("!!! Error converting input to int !!!")
            continue
        if option == 1:
            idx = _prompt_frame(n)
            if idx is not None:
                molecules[idx] = _menu_geom_analysis(molecules[idx])
        elif option == 2:
            if molecules[0].pbc is None:
                print("No PBC data available")
            else:
                res = _prompt_axis_times(molecules[0].pbc)
                if res is not None:
                    axis, times = res
                    molecules = [extend_axis(m, axis, times) for m in molecules]
                    print(f"\nExtended all {len(molecules)} frames")
        elif option == 3:
            for m in molecules:
                wrap_atoms(m, forward=True)
            print(f"\nWrapped atoms forward of origin in all {n} frames")
        elif option == 4:
            for m in molecules:
                wrap_atoms(m)
            print(f"\nWrapped atoms around origin in all {n} frames")
        elif option == 5:
            if molecules[0].pbc is None:
                print("No PBC data available")
            else:
                vals = _prompt_cell()
                for m in molecules:
                    if m.pbc is not None:
                        m.pbc.update(*vals)
                print(f"\nUpdated cell dimensions in all {n} frames")
        elif option == 6:
            out_filename = input("\noutput filename > ")
            with open(out_filename, "w") as out:
                for m in molecules:
                    write_xyz(m, out)
            print(f"wrote {n} frames to {out_filename}")
        elif option == 7:
            idx = _prompt_frame(n)
            if idx is not None:
                _write_mpmc_options(molecules[idx])
        elif option == 8:
            out_filename = input("\noutput filename > ")
            with open(out_filename, "w") as out:
                for m in molecules:
                    write_standard_pdb(m, out)
            print(f"wrote {n} frames to {out_filename}")
        elif option == 9:
            return
        else:
            print("\nInvalid option!")


def run_classic(filename: str) -> None:
    ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""

    is_trajectory = False
    if ext == "xyz":
        is_trajectory = check_xyz_trajectory(filename)
    elif ext in ("pdb", "ent", "pqr"):
        is_trajectory = check_pdb_trajectory(filename)
    elif ext == "cif":
        sys.exit("Use Mercury to save a .cif as a .xyz or .pdb first")
    else:
        sys.exit("Unable to determine if input file is xyz or pdb (please rename)")

    if is_trajectory:
        with open(filename) as f:
            if ext == "xyz":
                molecules, pbcs = read_xyz_trajectory(f)
            else:
                molecules, pbcs = read_pdb_trajectory(f)
        _main_loop_movie(molecules, filename)
    else:
        with open(filename) as f:
            if ext == "xyz":
                system, pbc = read_xyz(f, filename)
            else:
                system, pbc = read_pdb(f, filename)

        if pbc is None:
            while True:
                try:
                    print(f"Cell information not found in {filename}")
                    a = float(input("a>     "))
                    b = float(input("b>     "))
                    c = float(input("c>     "))
                    alpha = float(input("alpha> "))
                    beta = float(input("beta>  "))
                    gamma = float(input("gamma> "))
                    break
                except ValueError:
                    print("!!! Error converting input to float !!!\n")
            pbc = PBC(a, b, c, alpha, beta, gamma)

        mol = Molecule(atoms=system, pbc=pbc)
        mol.detect_bonds()
        _main_loop_single(mol, filename)

# ======================================================================
# Module: __main__
# ======================================================================
"""PDB Wizard entry point — defaults to interactive TUI viewer, --classic for text menu.

Batch CLI flags (--info, --rdf, --msd, --convert, --surface, --psd) run
analysis from the command line and exit without launching the TUI.
"""



# Prevent OpenBLAS/MKL threading conflicts with Python asyncio threads
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")


def _ensure_installed(package: str, import_name: str | None = None) -> bool:
    """Try to import a package; if missing, offer to pip install it."""
    import_name = import_name or package
    try:
        __import__(import_name)
        return True
    except ImportError:
        pass
    print(f"'{package}' is not installed.")
    try:
        answer = input("Install it now with pip? [Y/n] ").strip().lower()
    except EOFError:
        answer = "n"
    if answer in ("", "y", "yes"):
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "--quiet"])
        return True
    return False


# ------------------------------------------------------------------
# Batch-mode handlers
# ------------------------------------------------------------------

def _check_file(filepath: str) -> None:
    """Exit with error if file doesn't exist."""
    if not os.path.exists(filepath):
        print(f"Error: file not found: {filepath}", file=sys.stderr)
        sys.exit(1)


def _batch_info(filepath: str) -> None:
    """Print system info and exit."""

    _check_file(filepath)
    mol = read_file(filepath)
    print_info(mol, filepath)


def _batch_rdf(filepath: str, element1: str, element2: str, csv_path: str | None) -> None:
    """Compute RDF and optionally write CSV."""
    _check_file(filepath)

    frames = read_file_trajectory(filepath)
    if frames is not None:
        r, g = compute_rdf_trajectory(frames, element1, element2)
    else:
        mol = read_file(filepath)
        r, g = compute_rdf(mol, element1, element2)

    if len(r) == 0:
        print(f"No RDF data for {element1}-{element2} (check elements or PBC).",
              file=sys.stderr)
        sys.exit(1)

    if csv_path:
        import numpy as np
        np.savetxt(csv_path, np.column_stack([r, g]),
                   header="r(A),g(r)", delimiter=",", comments="")
        print(f"RDF written to {csv_path}")
    else:
        for ri, gi in zip(r, g):
            print(f"{ri:.4f}  {gi:.6f}")


def _batch_msd(filepath: str, element: str | None, csv_path: str | None) -> None:
    """Compute MSD and optionally write CSV."""
    _check_file(filepath)

    frames = read_file_trajectory(filepath)
    if frames is None:
        print("MSD requires a trajectory file with multiple frames.",
              file=sys.stderr)
        sys.exit(1)

    lags, msd = compute_msd(frames, element=element)

    if len(lags) == 0:
        print("No MSD data (check element or trajectory).", file=sys.stderr)
        sys.exit(1)

    if csv_path:
        import numpy as np
        np.savetxt(csv_path, np.column_stack([lags, msd]),
                   header="lag(frames),MSD(A^2)", delimiter=",", comments="")
        print(f"MSD written to {csv_path}")
    else:
        for li, mi in zip(lags, msd):
            print(f"{li:.0f}  {mi:.6f}")


def _batch_convert(filepath: str, output: str) -> None:
    """Convert between molecular file formats."""
    _check_file(filepath)
    from pathlib import Path


    mol = read_file(filepath)
    suffix = Path(output).suffix.lower()

    with open(output, "w") as out:
        if suffix == ".xyz":
            write_xyz(mol, out)
        elif suffix in (".pdb", ".ent", ".pqr"):
            write_standard_pdb(mol, out, skip_mols_step=True)
        elif suffix == ".cif":
            write_cif(mol, out)
        elif suffix in (".vasp", ".poscar"):
            write_poscar(mol, out)
        elif suffix in (".lmp", ".lammps", ".data"):
            write_lammps_data(mol, out)
        else:
            print(f"Unsupported output format: {suffix}", file=sys.stderr)
            sys.exit(1)

    print(f"Converted {filepath} -> {output}")


def _batch_surface(filepath: str) -> None:
    """Compute and print surface area."""
    _check_file(filepath)

    mol = read_file(filepath)
    result = surface_area(mol)
    print(f"Surface area: {result['surface_area']:.2f} A^2")
    print(f"Area per volume: {result['area_per_volume']:.2f} m^2/g")


def _batch_psd(filepath: str) -> None:
    """Compute and print pore size distribution."""
    _check_file(filepath)

    mol = read_file(filepath)
    centers, hist = pore_size_distribution(mol)

    if hist.sum() == 0:
        print("No pore size data (check PBC or atom count).", file=sys.stderr)
        sys.exit(1)

    print("pore_radius(A)  probability")
    for c, h in zip(centers, hist):
        print(f"{c:.4f}  {h:.6f}")


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

def main() -> None:
    if not _ensure_installed("numpy"):
        sys.exit(1)

    parser = argparse.ArgumentParser(
        prog="pdb_wizard",
        description="Molecular structure processor and interactive 3D terminal viewer",
    )
    parser.add_argument(
        "file", nargs="?",
        help="molecular structure file (.pdb, .xyz, .ent, .pqr, .zmat). "
             "Optional: launch with no file to start with an empty buffer "
             "and pick one later via File > Open.",
    )
    parser.add_argument(
        "--classic", action="store_true",
        help="use the classic text menu instead of the interactive TUI viewer",
    )

    # Batch-mode analysis flags
    batch = parser.add_argument_group("batch analysis (run and exit, no TUI)")
    batch.add_argument(
        "--info", action="store_true",
        help="print system info and exit",
    )
    batch.add_argument(
        "--rdf", action="store_true",
        help="compute radial distribution function (requires -e1 and -e2)",
    )
    batch.add_argument(
        "--msd", action="store_true",
        help="compute mean squared displacement (requires -e for element, trajectory file)",
    )
    batch.add_argument(
        "--convert", action="store_true",
        help="convert file to another format (requires -o for output path)",
    )
    batch.add_argument(
        "--surface", action="store_true",
        help="compute solvent-accessible surface area",
    )
    batch.add_argument(
        "--psd", action="store_true",
        help="compute pore size distribution",
    )

    # Element selectors
    batch.add_argument("-e", dest="element", default=None, help="element symbol (for --msd)")
    batch.add_argument("-e1", dest="element1", default=None, help="first element (for --rdf)")
    batch.add_argument("-e2", dest="element2", default=None, help="second element (for --rdf)")

    # Output options
    batch.add_argument("-o", "--output", default=None, help="output file path (for --convert)")
    batch.add_argument("--csv", default=None, help="CSV output file (for --rdf, --msd)")

    args = parser.parse_args()

    # ---- batch mode dispatch ----
    # Batch flags require a file argument.
    batch_flags = (args.info, args.rdf, args.msd, args.convert,
                   args.surface, args.psd)
    if any(batch_flags) and not args.file:
        parser.error("Batch flags require a file argument.")

    if args.info:
        _batch_info(args.file)
        sys.exit(0)

    if args.rdf:
        if not args.element1 or not args.element2:
            parser.error("--rdf requires -e1 and -e2")
        _batch_rdf(args.file, args.element1, args.element2, args.csv)
        sys.exit(0)

    if args.msd:
        _batch_msd(args.file, args.element, args.csv)
        sys.exit(0)

    if args.convert:
        if not args.output:
            parser.error("--convert requires -o / --output")
        _batch_convert(args.file, args.output)
        sys.exit(0)

    if args.surface:
        _batch_surface(args.file)
        sys.exit(0)

    if args.psd:
        _batch_psd(args.file)
        sys.exit(0)

    # ---- interactive mode (existing) ----
    if args.classic:
        if not args.file:
            parser.error("Classic menu requires a file argument.")
        run_classic(args.file)
    else:
        if not _ensure_installed("textual"):
            # textual unavailable (missing and not installed) — fall back to the
            # classic text menu rather than exiting, so the tool still works.
            if not args.file:
                parser.error(
                    "The interactive TUI needs 'textual'. Install it, or pass a "
                    "file to use the classic text menu."
                )
            print(
                "'textual' unavailable — falling back to the classic text menu "
                "(use --classic to select it directly).",
                file=sys.stderr,
            )
            run_classic(args.file)
            return

        if not args.file:
            # No-file launch: empty buffer, user picks via File > Open
            empty = Molecule(atoms=[])
            app = PdbWizardApp(molecule=empty, filepath="")
        else:
            ft = detect_filetype(args.file)
            is_traj = False
            if ft == "pdb":
                is_traj = check_pdb_trajectory(args.file)
            elif ft == "xyz":
                is_traj = check_xyz_trajectory(args.file)

            if is_traj:
                # Launch app empty, load trajectory in background
                empty = Molecule(atoms=[])
                app = PdbWizardApp(molecule=empty, filepath=args.file, load_trajectory=True)
            else:
                molecule = read_file(args.file)
                app = PdbWizardApp(molecule=molecule, filepath=args.file)
        app.run()


if __name__ == '__main__':
    main()
