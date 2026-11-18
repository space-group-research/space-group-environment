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
    try:
        if not _sys.stdin.isatty():
            return False  # never ask a script: the question would swallow its input
    except (AttributeError, ValueError):
        return False
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
from typing import Optional
from numpy import array, cos, floor, pi, sin, sqrt
from dataclasses import dataclass, field
import copy
import logging
import re
from pathlib import Path
from typing import Callable, Optional, TextIO
import struct
from typing import BinaryIO, Optional
import json
from typing import List, Optional
from typing import Callable, List, Optional
from typing import List, Optional, TextIO
import sys
from functools import reduce
from math import gcd
from typing import Any, Iterator, TextIO
import math
from typing import Callable
from numpy import cos, pi, sin, sqrt
from typing import Any
from collections import Counter
from typing import Callable, Dict, List, Optional, Sequence
from typing import Any, Dict, List, Optional, Tuple
from typing import Optional, TextIO
import urllib.parse
import urllib.request
import asyncio
import os
from typing import Callable, Dict, List, Optional
import concurrent.futures
import threading
import shutil
import tempfile
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
        "sigma", "c6", "c8", "c10", "id", "omega", "gwp_alpha", "c9",
        "mpmc_flag", "mpmc_mol", "mpmc_res", "file_mass",
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
        # MPMC columns after sigma that the tool does not use but must not lose.
        self.omega = 0.0
        self.gwp_alpha = 0.0
        self.c9 = 0.0
        self.id = 0
        # Set only when read from an MPMC file: F/M flag, molecule ID, residue
        # name and the mass column (dummy sites like M or LP are massless).
        self.mpmc_flag: Optional[str] = None
        self.mpmc_mol: Optional[int] = None
        self.mpmc_res: Optional[str] = None
        self.file_mass: Optional[float] = None

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

    def perpendicular_widths(self) -> np.ndarray:
        """Distance between each pair of opposite cell faces (Å), for a, b, c.

        Half the smallest of these is the largest cutoff for which the
        rounding minimum image (used here and in MPMC) finds every pair.
        For a skewed cell it is shorter than min(a, b, c) / 2.
        """
        m = self.basis_matrix
        areas = np.array([
            np.linalg.norm(np.cross(m[1], m[2])),
            np.linalg.norm(np.cross(m[2], m[0])),
            np.linalg.norm(np.cross(m[0], m[1])),
        ])
        return abs(self.volume) / areas

    def shortest_lattice_vector(self, max_coef: int = 5) -> float:
        """Length of the shortest nonzero lattice vector (Å).

        Mirrors MPMC's pbc_cutoff() (src/energy/pbc.c, MAX_VECT_COEF = 5):
        when mpmc.inp has no pbc_cutoff, MPMC uses half of this value.
        """
        r = np.arange(-max_coef, max_coef + 1)
        coefs = np.array(np.meshgrid(r, r, r, indexing="ij")).reshape(3, -1).T
        coefs = coefs[np.any(coefs != 0, axis=1)]
        vecs = coefs @ self.basis_matrix
        return float(np.sqrt((vecs * vecs).sum(axis=1)).min())

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


def cell_problem(a: float, b: float, c: float, alpha: float, beta: float, gamma: float):
    """Why these six numbers do not describe a unit cell, or None if they do."""
    values = (a, b, c, alpha, beta, gamma)
    if not all(np.isfinite(v) for v in values):
        return "Cell parameters must be finite numbers"
    if min(a, b, c) <= 0:
        return "Cell lengths must be positive"
    if not all(0 < x < 180 for x in (alpha, beta, gamma)):
        return "Cell angles must be between 0 and 180 degrees"
    ca, cb, cg = (np.cos(np.radians(x)) for x in (alpha, beta, gamma))
    # (V / abc)^2; zero or negative when the three angles cannot close a cell
    if 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg <= 1e-10:
        return "These angles do not make a cell (zero or negative volume)"
    return None

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


def apply_ff_to_system_report(
    system: list, ff: dict[str, dict[str, float]],
) -> tuple[list, dict[str, int]]:
    """Apply force-field parameters; return (system, {missing element: atom count}).

    Atoms whose element is absent from the force field get all parameters set
    to zero, as the original pdb_wizard did, and are counted in the returned
    dict so the caller can warn the user instead of silently writing zeros.
    Movable (M) molecules read from an MPMC file are skipped: they carry their
    own model parameters (e.g. a water model), which a framework force field
    must not overwrite.
    """
    missing: dict[str, int] = {}
    for atom in system:
        if atom.mpmc_flag == "M":
            continue
        el = atom.element.symbol
        if el in ff:
            params = ff[el]
            atom.alpha = params["alpha"]
            atom.sigma = params["sigma"]
            atom.epsilon = params["epsilon"]
            atom.c6 = params["c6"]
            atom.c8 = params["c8"]
            atom.c10 = params["c10"]
        else:
            atom.alpha = atom.sigma = atom.epsilon = 0.0
            atom.c6 = atom.c8 = atom.c10 = 0.0
            missing[el] = missing.get(el, 0) + 1
    return system, missing


def apply_ff_to_system(system: list, ff: dict[str, dict[str, float]]) -> list:
    system, _ = apply_ff_to_system_report(system, ff)
    return system


def missing_ff_messages(missing: dict[str, int]) -> list[str]:
    """Human-readable warnings for apply_ff_to_system_report's missing dict."""
    return [
        f"!!! {n} atom{'s' if n != 1 else ''} of element {el} not found in forcefield, "
        f"parameters set to all zeros !!!"
        for el, n in sorted(missing.items())
    ]

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

            # grid cells must be at least `cutoff` thick; in a skewed cell their thickness
            # comes from the perpendicular widths, which are shorter than the edges
            box = self.pbc.perpendicular_widths()
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


_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Single-frame readers
# ---------------------------------------------------------------------------

def _split_packed_floats(s: str) -> list[str]:
    """Split a string of packed floats that may lack spaces (e.g. '100.123-100.456')."""
    import re
    return re.findall(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?", s)


# An MPMC atom record carries x y z mass charge alpha epsilon sigma (and
# optionally omega gwp_alpha c6 c8 c10) after the six leading fields, so it has
# at least 14 whitespace tokens. A standard PDB line has at most 12 (record,
# serial, name, resName, chainID, resSeq, x, y, z, occupancy, B-factor,
# element, formal charge) and a pdb2pqr .pqr line 10-11.
MPMC_MIN_TOKENS = 14


def _is_decimal(tok: str) -> bool:
    if "." not in tok:
        return False
    try:
        float(tok)
        return True
    except ValueError:
        return False


# PDB coordinates are %8.3f: when two overflow into each other the only
# unambiguous way to split "-1234.5002345.250" is on exactly three decimals.
_PDB_COORD_RE = re.compile(r"[+-]?\d+\.\d{3}")


def _expand_numeric_tokens(tokens: list[str]) -> list[str]:
    """Split tokens that pack several numbers together ('-1000.000-1000.000').

    Packed tokens are split on the %8.3f pattern first (the only way packing
    arises), then on a general float pattern. Non-numeric tokens are kept
    as-is so a run of coordinates can still be located relative to them.
    """
    out: list[str] = []
    for tok in tokens:
        if _is_decimal(tok) or not any(ch.isdigit() for ch in tok):
            out.append(tok)
            continue
        parts = _PDB_COORD_RE.findall(tok)
        if len(parts) > 1 and "".join(parts) == tok:
            out.extend(parts)
            continue
        parts = _split_packed_floats(tok)
        if len(parts) > 1 and "".join(parts) == tok:
            out.extend(parts)
        else:
            out.append(tok)
    return out


def _coords_from_tokens(tokens: list[str]) -> tuple[float, float, float, list[str]]:
    """(x, y, z, tokens after z) from the first run of three decimal numbers
    after the residue fields.

    Used when fixed columns are not reliable: free-form spacing, or values too
    wide for their 8-character field that ran into the neighbouring one.
    Residue sequence numbers are integers, so requiring a decimal point skips
    them.
    """
    toks = _expand_numeric_tokens(tokens[4:])
    for i in range(len(toks) - 2):
        if _is_decimal(toks[i]) and _is_decimal(toks[i + 1]) and _is_decimal(toks[i + 2]):
            return float(toks[i]), float(toks[i + 1]), float(toks[i + 2]), toks[i + 3:]
    raise ValueError("no x y z coordinates found")


def _is_mpmc_record(tokens: list[str]) -> bool:
    """MPMC atom record: 14+ tokens with x y z mass charge alpha epsilon sigma
    all numeric. A standard PDB line with a segment ID, element and formal
    charge can also reach 14 tokens, but its token 11 is the segment ID.
    """
    if len(tokens) < MPMC_MIN_TOKENS:
        return False
    try:
        for tok in tokens[6:14]:
            float(tok)
    except ValueError:
        return False
    return True


# Marker lines/fields that identify a whole file as MPMC format even when
# some atom records are short (MPMC tolerates missing trailing columns).
_MPMC_FILE_MARKERS = ("REMARK BOX BASIS",)
_MPMC_FLAGS = ("F", "M", "A", "S", "T")  # frozen / movable / adiabatic / spectre / target


def _is_mpmc_dummy_site(atom: Atom) -> bool:
    """A massless model site from an MPMC file: water M/LP charge sites, H2
    centers of mass (H2G), and so on. Shown as element X, but expected."""
    return atom.file_mass is not None and atom.file_mass < 0.5


_SITE_MASS_TOLERANCE = 0.02  # amu


def site_element(atom: Atom):
    """The chemical element an atom really is, or None for a massless model site.

    Sites of MPMC sorbate models are named after the model, not the element ('COG' and
    'COE' are the carbon and oxygens of CO2, 'CHE' a methane hydrogen), so ``atom.element``,
    which is read from the name, can be wrong for them. Their mass column is reliable: it
    is matched to the periodic table here. Framework atoms, atoms without an MPMC mass,
    and united-atom sites whose mass is no single element's keep ``atom.element``."""
    if atom.file_mass is None or atom.mpmc_flag != "M":
        return atom.element  # framework atoms are named after their element
    if atom.file_mass < 0.5:
        return None
    best = min((el for el in ELEMENTS.values() if el.atomic_number > 0),
               key=lambda el: abs(el.mass - atom.file_mass))
    if abs(best.mass - atom.file_mass) <= _SITE_MASS_TOLERANCE:
        return best
    return atom.element


def _mpmc_mass(atom: Atom) -> float:
    """Mass to write for a movable site: the file's own value when it came from
    MPMC (dummy sites are massless). Framework atoms keep element masses."""
    return atom.file_mass if atom.file_mass is not None else atom.mass


def _mpmc_omega_fields(atom: Atom) -> str:
    """The omega and gwp_alpha columns: '0.0 0.0' as always when both are zero, else the
    values read from the file."""
    if not atom.omega and not atom.gwp_alpha:
        return "0.0 0.0"
    return f"{atom.omega:.5f} {atom.gwp_alpha:.5f}"


def _mpmc_c9_field(atom: Atom) -> str:
    """The optional c9 column (after c10), written only when the file had a value."""
    return f" {atom.c9:.5f}" if atom.c9 else ""


def export_sites(mol: Molecule) -> list:
    """(atom, element symbol) for every atom a non-MPMC file can hold.

    MPMC sorbate sites are named after their model ('COG' is a carbon), so the element
    comes from ``site_element``; massless model sites (water M, H2 centre) are left out.
    """
    out = []
    for atom in mol.atoms:
        element = site_element(atom)
        if element is not None:
            out.append((atom, element.symbol))
    return out


def lattice_to_standard(basis) -> tuple:
    """(PBC, M) for three lattice vectors (rows) in any orientation.

    The PBC always has a along x and b in the xy plane. M is None when the vectors are
    already in that orientation; otherwise ``x @ M`` moves Cartesian positions given in the
    vectors' frame into the PBC's frame (same fractional coordinates).
    """
    basis = np.array(basis, dtype=float)
    lengths = [float(np.linalg.norm(v)) for v in basis]
    if min(lengths) <= 0.0 or abs(np.linalg.det(basis)) < 1e-12:
        raise ValueError("Lattice vectors are zero or coplanar")

    def angle(u, v) -> float:
        cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))

    pbc = PBC(lengths[0], lengths[1], lengths[2], angle(basis[1], basis[2]),
              angle(basis[0], basis[2]), angle(basis[0], basis[1]))
    if np.abs(basis - pbc.basis_matrix).max() < 1e-6:
        return pbc, None
    return pbc, np.linalg.inv(basis) @ pbc.basis_matrix


def assign_framework_charges(mol: Molecule, charges: list[float]) -> str | None:
    """Put a list of charges on the framework atoms, in order.

    Movable (M) molecules read from an MPMC file keep their model charges, so
    the charges file only has to cover the framework. A list that divides the
    framework count is repeated (charges for one cell of a supercell). Returns
    None on success or a message saying why nothing was assigned.
    """
    framework = [a for a in mol.atoms if a.mpmc_flag != "M"]
    n = len(framework)
    if not charges or n % len(charges) != 0:
        return f"Charge count ({len(charges)}) doesn't match framework atom count ({n})"
    for ind, atom in enumerate(framework):
        atom.charge = charges[ind % len(charges)]
    return None


def mpmc_movable_count(mol: Molecule) -> int:
    """Number of movable (M) molecules read from an MPMC file and still in ``mol``."""
    return len({a.mpmc_mol for a in mol.atoms if a.mpmc_flag == "M"})


def next_mpmc_mol_id(mol: Molecule) -> int:
    """Molecule id for the next sorbate added to ``mol``: one past the highest movable id
    (the framework is molecule 1). Counting molecules is not enough: supercells and
    deletions leave gaps, and a reused id merges two sorbates into one MPMC molecule."""
    ids = [a.mpmc_mol for a in mol.atoms if a.mpmc_flag == "M" and isinstance(a.mpmc_mol, int)]
    return max(ids) + 1 if ids else 2


def _file_is_mpmc(lines: list[str]) -> bool:
    for line in lines:
        if line.startswith(_MPMC_FILE_MARKERS):
            return True
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            tokens = line.split()
            if len(tokens) > 5 and tokens[3] == "BOX" and tokens[4] in _MPMC_FLAGS:
                return True
            if _is_mpmc_record(tokens):
                return True
    return False


def _parse_mpmc_atom(tokens: list[str]) -> Atom:
    """MPMC dialect: whitespace tokens, free field widths.

    Layout: ATOM serial name resName flag molID x y z mass charge alpha epsilon
    sigma [omega gwp_alpha c6 c8 c10]. MPMC itself reads this with
    ``sscanf("%s %s ...")``, so columns carry no meaning here.
    """
    name = tokens[2]
    atom = Atom(float(tokens[6]), float(tokens[7]), float(tokens[8]), name)
    if len(tokens) > 5 and tokens[4] in _MPMC_FLAGS:
        atom.mpmc_flag = tokens[4]
        atom.mpmc_res = tokens[3]
        if tokens[5].isdigit():
            atom.mpmc_mol = int(tokens[5])
    try:
        atom.file_mass = float(tokens[9])
    except (ValueError, IndexError):
        pass
    try:
        atom.charge = float(tokens[10])
    except (ValueError, IndexError):
        pass
    try:
        atom.alpha = float(tokens[11])
        atom.epsilon = float(tokens[12])
        atom.sigma = float(tokens[13])
    except (ValueError, IndexError):
        pass
    try:
        atom.omega = float(tokens[14])
        atom.gwp_alpha = float(tokens[15])
    except (ValueError, IndexError):
        pass
    if len(tokens) > 18:
        try:
            atom.c6 = float(tokens[16])
            atom.c8 = float(tokens[17])
            atom.c10 = float(tokens[18])
        except (ValueError, IndexError):
            pass
    if len(tokens) > 19:
        try:
            atom.c9 = float(tokens[19])
        except ValueError:
            pass
    return atom


def _parse_formal_charge(field: str) -> Optional[float]:
    """PDB columns 79-80: '2+', '1-' (or the lenient '+2' / '-1')."""
    f = field.strip()
    if len(f) != 2:
        return None
    if f[0].isdigit() and f[1] in "+-":
        return float(f[0]) * (1.0 if f[1] == "+" else -1.0)
    if f[1].isdigit() and f[0] in "+-":
        return float(f[1]) * (1.0 if f[0] == "+" else -1.0)
    return None


def _coordinate_columns_aligned(line: str) -> bool:
    """True if the x/y/z fields (cols 31-54) are bounded by blanks.

    A value too wide for its 8-character field runs into the neighbouring
    column; the fixed-column slice would then still parse, but to a wrong
    (truncated or sign-flipped) number. Column 30 is blank by specification
    and column 55 is the first character of the occupancy field.
    """
    if len(line) > 29 and line[29] != " ":
        return False
    if len(line) > 54 and line[54] not in " \r\n":
        return False
    return True


def _element_column(line: str) -> str:
    """Element symbol from PDB columns 77-78, tolerating an off-by-one shift.

    The symbol is right-justified in columns 77-78 (0-based 76:78), but many
    writers put it one column early or late. Take the alphabetic run touching
    those columns; give up (return '') if it is longer than two letters, which
    means it merged with a segment identifier.
    """
    n = len(line)
    for start in (76, 77, 75):
        if start < n and line[start].isalpha():
            left = start
            while left > 72 and line[left - 1].isalpha():
                left -= 1
            right = start
            while right + 1 < n and line[right + 1].isalpha():
                right += 1
            tok = line[left:right + 1]
            return tok if len(tok) <= 2 else ""
    return ""


def _parse_standard_atom(line: str, tokens: list[str], pqr: bool = False) -> Atom:
    """Standard PDB (fixed columns) or pdb2pqr .pqr dialect.

    Coordinates come from columns 31-54 as the format specifies; a line whose
    columns are misaligned falls back to the first run of three decimals.
    The element column (77-78) overrides name-based inference, so 'CA' is an
    alpha carbon and 'HO' a hydrogen, not calcium and holmium. The charge
    comes from the formal-charge columns (79-80) or, for pdb2pqr files, from
    the token after z — never from the B-factor.
    """
    name = line[12:16].strip() or tokens[2]
    after_z: Optional[list[str]] = None
    try:
        if not _coordinate_columns_aligned(line):
            raise ValueError("coordinate fields overflow their columns")
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
    except ValueError:
        x, y, z, after_z = _coords_from_tokens(tokens)
        name = tokens[2]
    atom = Atom(x, y, z, name)

    elem_col = _element_column(line)
    if elem_col and not pqr:
        if elem_col.upper() == "D":
            elem_col = "H"  # deuterium: no separate entry in the element table
        el = get_element(elem_col)
        if el.symbol != "X":
            atom._element = el

    if pqr:
        # pdb2pqr / APBS: ... x y z charge radius
        if after_z is None:
            try:
                _, _, _, after_z = _coords_from_tokens(tokens)
            except ValueError:
                after_z = None
        if after_z:
            try:
                atom.charge = float(after_z[0])
            except ValueError:
                pass
    else:
        fc = _parse_formal_charge(line[78:80])
        if fc is not None:
            atom.charge = fc
    return atom


def _parse_pdb_atom(line: str, pqr: bool = False, mpmc: bool = False) -> Atom:
    """Parse an ATOM/HETATM line, detecting the dialect per line.

    1. Serial overflow (``HETATM99999``): fixed-width name + regex coordinates.
    2. MPMC (14+ numeric-tailed tokens, or ``mpmc=True`` because the file
       carries MPMC markers): whitespace tokens, extended force-field columns.
    3. Standard PDB / pdb2pqr: fixed columns with a token fallback.

    ``pqr`` says the file has a .pqr suffix, which only matters for the
    pdb2pqr charge/radius convention; MPMC files also use .pqr but are
    recognised by their content.
    """
    tokens = line.split()

    if 6 < len(tokens[0]) <= 11 and tokens[0][6:].isdigit():
        # 'HETATM10000': the serial fills its five columns, nothing has moved.
        tokens = [tokens[0][:6], tokens[0][6:]] + tokens[1:]

    if len(tokens[0]) > 6:
        name = line[12:16].strip()
        nums = [t for t in _expand_numeric_tokens(line[30:].split()) if _is_decimal(t)]
        if len(nums) < 3:
            raise ValueError("no x y z coordinates found")
        return Atom(float(nums[0]), float(nums[1]), float(nums[2]), name)

    if mpmc or _is_mpmc_record(tokens):
        return _parse_mpmc_atom(tokens)

    return _parse_standard_atom(line, tokens, pqr=pqr)


def _is_placeholder_cryst1(values: list) -> bool:
    """'CRYST1 1.000 1.000 1.000 90.00 90.00 90.00': what the PDB format prescribes for a
    structure that is not a crystal. It is not a 1 A cell."""
    return (all(abs(v - 1.0) < 1e-6 for v in values[:3])
            and all(abs(v - 90.0) < 1e-6 for v in values[3:6]))


def read_pdb(file: TextIO, filename: str = "",
             notes: Optional[list] = None) -> tuple[list[Atom], Optional[PBC]]:
    """Read atoms and cell from a PDB / MPMC PDB / pdb2pqr PQR file.

    MPMC box-corner pseudo-atoms (element X in residue BOX) are dropped. Other
    atoms with an unrecognised element are kept as element 'X' and reported
    through the module logger, instead of vanishing silently.
    """
    system: list[Atom] = []
    pbc: Optional[PBC] = None
    pqr = filename.lower().endswith(".pqr")
    unknown_names: dict[str, int] = {}
    skipped: list[int] = []
    box_basis: dict[int, list] = {}
    altlocs: set = set()

    lines = file.readlines()
    # A file with MPMC markers (box corners, basis remarks, or any full MPMC
    # record) is parsed as MPMC throughout, so short records that only carry
    # mass and charge are not mistaken for standard PDB occupancy/B-factor.
    mpmc_file = _file_is_mpmc(lines)
    for lineno, line in enumerate(lines, 1):
        if len(line) < 4:
            continue
        # Fast record detection without full split
        is_atom = line[0] == "A" and line[1:4] == "TOM"
        is_hetatm = line[0] == "H" and line[1:6] == "ETATM"
        if is_atom or is_hetatm:
            try:
                atom = _parse_pdb_atom(line, pqr=pqr, mpmc=mpmc_file)
            except (ValueError, IndexError):
                skipped.append(lineno)
                continue
            if atom.atomic_number == 0:
                tokens = line.split()
                if len(tokens) > 3 and tokens[3] == "BOX":
                    continue  # MPMC box corner, not an atom
                if not _is_mpmc_dummy_site(atom):
                    unknown_names[atom.name] = unknown_names.get(atom.name, 0) + 1
            if (not mpmc_file and len(line) > 54 and line[16] != " " and line[20] == " "
                    and _coordinate_columns_aligned(line)):
                altlocs.add(line[16])
            system.append(atom)
            continue
        if line[0] == "R" and line[:6] == "REMARK":
            tokens = line.split()
            if line.startswith("REMARK BOX BASIS[") and len(tokens) >= 6:
                try:
                    box_basis[int(tokens[2][6])] = [float(t) for t in tokens[-3:]]
                except (ValueError, IndexError):
                    pass
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
                values = [float(t) for t in tokens[1:7]]
                if len(values) == 6 and not _is_placeholder_cryst1(values):
                    pbc = PBC(*values)
            except (ValueError, IndexError):
                pass
            continue
        if line[:3] == "END":
            break

    if sorted(box_basis) == [0, 1, 2]:
        # MPMC's exact cell. Its vectors may point anywhere (mpmc.inp basis1/2/3); when
        # they are not a-along-x, b-in-xy the atoms are turned with the cell so both agree.
        try:
            exact, turn = lattice_to_standard([box_basis[0], box_basis[1], box_basis[2]])
        except ValueError:
            turn = None
        if turn is not None:
            pbc = exact
            for atom in system:
                atom.x = atom.x @ turn
    if notes is not None and len(altlocs) > 1:
        notes.append(
            f"The PDB has alternate locations ({', '.join(sorted(altlocs))}). Every one was "
            f"kept as an atom, so disordered groups appear more than once. Remove the "
            f"extra copies before using it.")

    set_atom_ids(system)
    where = f" in {filename}" if filename else ""
    if skipped:
        shown = ", ".join(str(n) for n in skipped[:5])
        more = f" (+{len(skipped) - 5} more)" if len(skipped) > 5 else ""
        _log.warning(f"Skipped {len(skipped)} unparsable atom line(s){where}: {shown}{more}")
    if unknown_names:
        names = ", ".join(f"{n} x{c}" for n, c in sorted(unknown_names.items()))
        _log.warning(f"Unrecognised element names{where} kept as 'X': {names}")
    return system, pbc


def _xyz_atom(tokens: list[str]) -> Atom:
    """Atom from an XYZ record; the symbol column may hold an atomic number."""
    symbol = tokens[0]
    if symbol.isdigit():
        el = get_element_by_number(int(symbol))
        atom = Atom(tokens[1], tokens[2], tokens[3], el.symbol)
        atom._element = el
        return atom
    return Atom(tokens[1], tokens[2], tokens[3], symbol)


def _looks_like_xyz_atom(line: str) -> bool:
    tokens = line.split()
    if len(tokens) < 4:
        return False
    try:
        float(tokens[1]), float(tokens[2]), float(tokens[3])
    except ValueError:
        return False
    return True


def read_xyz(file: TextIO, filename: str = "",
             notes: Optional[list] = None) -> tuple[list[Atom], Optional[PBC]]:
    """Read the first frame of an XYZ file. ``notes`` (a list) receives a sentence when the
    atom count in the header and the atom lines disagree; the atoms read are returned
    either way."""
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

    problem = ""
    try:
        if n_atoms is not None:
            # Read exactly n_atoms lines (handles multi-frame files)
            for _ in range(n_atoms):
                line = file.readline()
                if not line:
                    problem = (f"The XYZ header says {n_atoms} atoms but the file ends after "
                               f"{len(system)}.")
                    break
                tokens = line.split()
                if not _looks_like_xyz_atom(line):
                    problem = (f"XYZ atom line {len(system) + 1} could not be read "
                               f"({line.strip()!r}); only the {len(system)} atoms before it "
                               f"of {n_atoms} were loaded.")
                    break
                atom = _xyz_atom(tokens)
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
                atom = _xyz_atom(tokens)
                try:
                    atom.charge = float(tokens[4])
                except (ValueError, IndexError):
                    pass
                system.append(atom)
    except (ValueError, IndexError):
        pass  # stop reading on malformed line

    if n_atoms is not None and not problem and len(system) == n_atoms:
        rest = file.readline()
        while rest and not rest.strip():
            rest = file.readline()
        if rest and _looks_like_xyz_atom(rest):
            problem = (f"The XYZ header says {n_atoms} atoms but more atom lines follow; "
                       f"only the first {n_atoms} were loaded.")
    if problem:
        _log.warning(problem)
        if notes is not None:
            notes.append(problem)
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
    expr = expr.strip().replace(" ", "").lower()  # 'X,Y,Z' is as valid as 'x,y,z'

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
    basis: Optional[np.ndarray] = None,
    max_distance: float = 0.5,
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
                if np.linalg.norm(diff) < tolerance and (
                        basis is None or np.linalg.norm(diff @ basis) < max_distance):
                    # In a large cell 0.01 fractional is more than a bond length, so two
                    # images are also required to be closer than max_distance (A).
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


def read_cif(file: TextIO, filename: str = "",
             notes: Optional[list] = None) -> tuple[list[Atom], Optional[PBC]]:
    """Read a CIF file, returning atoms (in Cartesian) and PBC.

    Handles:
    - P1 CIF files with fractional coordinates (most common for MOFs)
    - CIF files with symmetry operations that need expansion
    - Both _atom_site_fract_x and _atom_site_Cartn_x

    ``notes`` (a list) receives a sentence for each thing the reader could not honour: a
    space group named without its symmetry operators (the atoms are then only the
    asymmetric unit) and partially occupied or disordered sites (all kept as full atoms).
    """
    lines = file.readlines()
    _sg_number = ""

    # --- Phase 1: Parse all key-value pairs and loops ---
    cell_params: dict[str, float] = {}
    symops: list[str] = []
    atom_loop_tags: list[str] = []
    atom_loop_rows: list[list[str]] = []
    _space_group: str = ""
    extra_blocks = False

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i].strip()

        # Skip comments and blank lines
        if not line or line.startswith("#"):
            i += 1
            continue

        # data_ block headers: the first block that has atoms is the structure
        if line.startswith("data_"):
            if atom_loop_rows:
                extra_blocks = True
                break
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
            loop_tokens: list[str] = []

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
                    loop_tokens.append(text_val.strip())
                    continue
                tokens = _tokenize_cif_line(data_line)
                if tokens:
                    loop_rows.append(tokens)
                    loop_tokens.extend(tokens)
                i += 1

            # A row may be wrapped over several lines: when the values divide evenly by
            # the number of tags, rows are cut by count instead of by line.
            n_tags = len(loop_tags)
            if (n_tags and loop_tokens and len(loop_tokens) % n_tags == 0
                    and any(len(row) != n_tags for row in loop_rows)):
                loop_rows = [loop_tokens[k:k + n_tags]
                             for k in range(0, len(loop_tokens), n_tags)]

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

            if any(tag.startswith("_atom_site.") for tag in loop_tags):
                raise ValueError(
                    "This is an mmCIF (PDBx) file, which is not supported. Download the "
                    "structure as .pdb instead.")

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
        if len(tokens) == 1 and tokens[0].startswith("_"):
            # the value may be on the next line
            j = i + 1
            while j < n and (not lines[j].strip() or lines[j].strip().startswith("#")):
                j += 1
            nxt = lines[j].strip() if j < n else ""
            if nxt and not nxt.startswith(("_", ";", "data_")) and nxt.lower() != "loop_":
                tokens += _tokenize_cif_line(nxt)[:1]
                i = j
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
                _space_group = value.strip().strip("'\"")
            if tag in ("_symmetry_int_tables_number", "_space_group_it_number"):
                _sg_number = value.strip()

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
    if notes is not None and is_p1:
        name = _space_group.replace(" ", "").replace("_", "").upper()
        named_higher = bool(name) and name not in ("P1", "?", ".")
        numbered_higher = _sg_number.isdigit() and int(_sg_number) > 1
        if named_higher or numbered_higher:
            which = _space_group if named_higher else f"number {_sg_number}"
            notes.append(
                f"The CIF names space group {which} but lists no symmetry operators, so it "
                f"was read as P1: the atoms are only the asymmetric unit, not the full cell. "
                f"Re-export it with symmetry operators, or as P1 with all atoms (Mercury, "
                f"VESTA), before using it.")

    # --- Phase 4: Extract atom data ---
    if notes is not None and extra_blocks:
        notes.append("The CIF holds more than one data block. Only the first one with atoms "
                     "was read.")

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
    no_coordinates = 0

    for row in atom_loop_rows:
        # Get label and symbol
        label = ""
        symbol = ""
        if "_atom_site_label" in col and col["_atom_site_label"] < len(row):
            label = row[col["_atom_site_label"]]
        if "_atom_site_type_symbol" in col and col["_atom_site_type_symbol"] < len(row):
            symbol = row[col["_atom_site_type_symbol"]]

        # If no type_symbol, derive from label: its leading letters are the element
        # ('C1A' is carbon 1 part A, not calcium)
        if not symbol and label:
            lead = re.match(r"[A-Za-z]*", label).group(0)
            symbol = (lead or "".join(c for c in label if c.isalpha()))[:2]
            # Normalize: capitalize first, lower rest
            if len(symbol) > 1:
                symbol = symbol[0].upper() + symbol[1].lower()
            else:
                symbol = symbol.upper()
        if not label and symbol:
            label = symbol

        if not label and not symbol:
            continue

        # A site without usable coordinates ('?', '.', a short row) is left out: it used
        # to become an atom at the origin.
        try:
            if has_frac:
                frac_xyz = np.array([_strip_cif_uncertainty(row[col[k]]) for k in (
                    "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z")])
            if has_cart:
                cart_xyz = np.array([_strip_cif_uncertainty(row[col[k]]) for k in (
                    "_atom_site_cartn_x", "_atom_site_cartn_y", "_atom_site_cartn_z")])
        except (ValueError, IndexError):
            no_coordinates += 1
            continue

        labels.append(label)
        symbols.append(symbol)
        if has_frac:
            frac_coords.append(frac_xyz)
        if has_cart:
            cart_coords.append(cart_xyz)

    if notes is not None and no_coordinates:
        notes.append(f"{no_coordinates} site{'s' if no_coordinates != 1 else ''} in the CIF "
                     f"had no coordinates and {'were' if no_coordinates != 1 else 'was'} "
                     f"left out.")

    if not labels:
        return [], pbc

    if notes is not None:
        partial = 0
        groups = set()
        for row in atom_loop_rows:
            idx = col.get("_atom_site_occupancy")
            if idx is not None and idx < len(row):
                try:
                    if _strip_cif_uncertainty(row[idx]) < 0.999:
                        partial += 1
                except ValueError:
                    pass
            idx = col.get("_atom_site_disorder_group")
            if idx is not None and idx < len(row) and row[idx] not in (".", "?", "0", ""):
                groups.add(row[idx])
        if partial or groups:
            what = []
            if partial:
                what.append(f"{partial} site{'s' if partial != 1 else ''} with occupancy "
                            f"below 1")
            if groups:
                what.append(f"{len(groups)} disorder group{'s' if len(groups) != 1 else ''}")
            notes.append(
                f"The CIF has {' and '.join(what)}. Every site was kept as a full atom, so "
                f"alternative positions overlap and the composition is wrong. Remove the "
                f"extra sites before using it.")

    # --- Phase 5: Apply symmetry and convert to Cartesian ---
    system: list[Atom] = []

    if has_frac and pbc is not None:
        frac_arr = np.array(frac_coords)

        if not is_p1:
            frac_arr, labels, symbols = _apply_symops(
                frac_arr, labels, symbols, parsed_symops, basis=pbc.basis_matrix
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
        sites = export_sites(mol)  # sorbate sites by their real element, no massless ones
        fh.write(str(len(sites)))
        if mol.pbc is not None:
            pbc = mol.pbc
            fh.write(f"\n{pbc.a} {pbc.b} {pbc.c} {pbc.alpha} {pbc.beta} {pbc.gamma}\n")
        else:
            fh.write("\n\n")
        for atom, symbol in sites:
            fh.write(f"{symbol} {atom.x[0]} {atom.x[1]} {atom.x[2]}\n")
    finally:
        if close:
            fh.close()


class CellRequiredError(ValueError):
    """Raised by writers for formats that cannot be written without a unit cell."""


def require_cell(mol: Molecule, what: str) -> None:
    if mol.pbc is None:
        raise CellRequiredError(
            f"{what} requires a unit cell (PBC); this structure has none. "
            "Set one with Update Cell first."
        )


def write_standard_pdb(mol: Molecule, out: TextIO | str, skip_mols_step: bool = False) -> None:
    """Write a standardized PDB (CRYST1, HETATM records, END).

    With ``skip_mols_step`` the atoms are written in their current order and
    positions; otherwise molecules are found, sorted, named and re-imaged whole.
    A structure without a cell is written without a CRYST1 record (and without
    molecule re-imaging) — a valid non-periodic PDB, as OpenMM NoCutoff setups
    need. Earlier versions silently wrote nothing in that case.
    """
    import copy


    # Work on a copy so we don't reorder the caller's molecule
    mol = copy.deepcopy(mol)
    pbc = mol.pbc
    sites = export_sites(mol)
    if len(sites) != len(mol.atoms) or any(a.element.symbol != sym for a, sym in sites):
        # MPMC sorbate sites: write the element their mass says (the name 'COG' reads as
        # cobalt) and leave massless model sites out. Only this copy is changed.
        for atom, symbol in sites:
            atom._element = get_element(symbol)
        mol.atoms = [atom for atom, _ in sites]
        mol.bonds = []
    if skip_mols_step:
        mols_list = [mol]
    else:
        mol, mols_list = sort_system(mol, output=False)
    mols_list = [m for m in mols_list if m.atoms]  # nothing to write for an empty one

    fh, close = _resolve_out(out)
    out = fh
    out.write("MODEL        1\n")
    out.write("COMPND    " + " " * 69 + "\n")
    out.write("AUTHOR    GENERATED BY PDB WIZARD\n")
    if pbc is not None:
        out.write(
            f"CRYST1{round(pbc.a, 3):>9}{round(pbc.b, 3):>9}{round(pbc.c, 3):>9} "
            f"{round(pbc.alpha, 2):>6} {round(pbc.beta, 2):>6} {round(pbc.gamma, 2):>6} "
            "P 1           1\n"
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
            if not skip_mols_step and pbc is not None:
                # Keep each molecule whole: image every atom next to its base atom.
                dx = pbc.wrap(atom.x - base_atom.x)
                atom.x = base_atom.x + dx

            # Exact PDB columns: serial 7-11, name 13-16 (short names start at
            # column 14), resName 18-20, chain 22, resSeq 23-26, x/y/z 31-54,
            # occupancy 55-60, B-factor 61-66, element 77-78.
            name = atom.name[:4]
            name_field = f" {name:<3}" if len(name) < 4 else name
            lines.append(
                # serial and residue number wrap so that every column stays in place
                f"HETATM{atom.id % 100000:>5} {name_field} {mol_name:>3} "
                f"A{(idx + 1) % 10000:>4}    "
                f"{atom.x[0]:8.3f}{atom.x[1]:8.3f}{atom.x[2]:8.3f}"
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
    """Write an MPMC-format PDB (what MPMC calls a .pqr).

    Every atom record is written complete — mass, charge, alpha, epsilon,
    sigma, c6, c8, c10 — using the values currently stored on the atoms, so the
    file is always loadable by MPMC. ``write_charges`` / ``write_params`` are
    kept for API compatibility and record whether the caller assigned charges /
    force-field parameters; they do not change the record layout (the original
    tool wrote ``xxxNAMExxx`` placeholders instead, which MPMC cannot read).
    Raises CellRequiredError when ``mol.pbc`` is None.
    """
    import copy


    require_cell(mol, "Writing an MPMC PDB")

    # Work on a copy so we don't reorder the caller's molecule
    mol = copy.deepcopy(mol)
    pbc = mol.pbc
    # Movable (M) molecules read from an MPMC file stay movable; everything
    # else is framework and gets sorted as before.
    movable = [a for a in mol.atoms if a.mpmc_flag == "M"]
    mol = Molecule(atoms=[a for a in mol.atoms if a.mpmc_flag != "M"], pbc=pbc)
    if mol.atoms:
        mol, _ = sort_system(mol)
    groups: dict[object, list[Atom]] = {}
    for a in movable:
        groups.setdefault(a.mpmc_mol, []).append(a)

    # Box molecule ID: one past the last movable/sorbate molecule (2 if none).
    # MPMC skips BOX records, but distinct IDs keep the file unambiguous.
    box_mol_id = 2 + len(groups)
    for sl in sorbate_lines or ():
        tok = sl.split()
        if len(tok) > 5 and tok[5].isdigit():
            box_mol_id = max(box_mol_id, int(tok[5]) + 1)

    fh, close_fh = _resolve_out(filename)
    try:
        parts: list[str] = [
            "MODEL        1\n",
            "COMPND    " + " " * 69 + "\n",
            "AUTHOR    GENERATED BY PDB WIZARD\n",
            (f"CRYST1{round(pbc.a, 3):>9}{round(pbc.b, 3):>9}{round(pbc.c, 3):>9} "
             f"{round(pbc.alpha, 2):>6} {round(pbc.beta, 2):>6} {round(pbc.gamma, 2):>6} P 1           1\n"),
        ]
        # MOF atoms (molecule 1, frozen)
        for idx, atom in enumerate(mol.atoms):
            parts.append(
                f"ATOM {idx + 1:>6} {atom.name:<4} MOF F    1    "
                f"{atom.x[0]:7.3f} {atom.x[1]:7.3f} {atom.x[2]:7.3f}"
                f" {atom.mass:>9.5f} {atom.charge:>9.5f}"
                f" {atom.alpha:>9.5f} {atom.epsilon:>9.5f} {atom.sigma:>9.5f}"
                # omega and gwp_alpha (unused) then the PHAHST dispersion
                # coefficients at full precision. The original wrote c10 with
                # only two significant digits ({:>10.2}), e.g. 13951.5 -> 1.4e+04,
                # which MPMC happily read as a 0.35% error.
                f" {_mpmc_omega_fields(atom)} {atom.c6:>10.5f} {atom.c8:>11.4f} {atom.c10:>12.3f}"
                f"{_mpmc_c9_field(atom)}\n"
            )

        next_atom_id = len(mol.atoms) + 1

        # Movable molecules from the input file (molecules 2, 3, ...), each made
        # whole again in case a wrap split it across a cell face.
        for k, atoms in enumerate(groups.values()):
            ref = atoms[0].x
            for atom in atoms:
                x = ref + pbc.wrap(atom.x - ref)
                parts.append(
                    f"ATOM {next_atom_id:>6} {atom.name:<4} {(atom.mpmc_res or 'MOL'):<3} M"
                    f" {k + 2:>4}    "
                    f"{x[0]:>10.6f} {x[1]:>10.6f} {x[2]:>10.6f}"
                    f" {_mpmc_mass(atom):>9.5f} {atom.charge:>9.5f}"
                    f" {atom.alpha:>9.5f} {atom.epsilon:>9.5f} {atom.sigma:>9.5f}"
                    f" {_mpmc_omega_fields(atom)} {atom.c6:>10.5f} {atom.c8:>12.5f}"
                    f" {atom.c10:>14.5f}{_mpmc_c9_field(atom)}\n"
                )
                next_atom_id += 1

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
        parts.extend(f"CONECT {next_atom_id + i:>4} {next_atom_id + j:>4}\n"
                     for i, j in pbc.edges())
        # MPMC's read_pqr_box needs the "=": it checks that the 4th token is "=".
        parts.extend(f"REMARK BOX BASIS[{r_idx}] = "
                     f"{row[0]:20.14f} {row[1]:20.14f} {row[2]:20.14f}\n"
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


def _is_atom_record(line: str) -> bool:
    return line[:4] == "ATOM" or line[:6] == "HETATM"


def _is_frame_end(line: str) -> bool:
    """END / ENDMDL record (VMD and many converters separate frames with END)."""
    rec = line[:6].rstrip()
    return rec in ("END", "ENDMDL")


def check_pdb_trajectory(filename: str) -> bool:
    """True if the file holds more than one frame.

    Frames may be introduced by ``MODEL`` or ``REMARK step=`` records, or
    simply separated by ``END``/``ENDMDL`` with further ATOM records after it.
    """
    try:
        with open(filename) as f:
            n_model = 0
            n_remark_step = 0
            atoms_in_block = False
            block_closed = False
            for line in f:
                if line[:6] == "MODEL ":
                    n_model += 1
                elif line.startswith("REMARK step="):
                    n_remark_step += 1
                elif _is_atom_record(line):
                    if block_closed:
                        return True
                    atoms_in_block = True
                elif _is_frame_end(line) and atoms_in_block:
                    block_closed = True
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

    while True:
        while line != "" and not line.strip():
            line = file.readline()  # blank lines between frames or at the end
        if line == "":
            break
        try:
            n_atoms = int(line)
        except ValueError:
            problem = (f"XYZ trajectory: expected an atom count at the start of frame "
                       f"{len(molecules) + 1}, found {line.strip()!r}")
            if not molecules:
                raise ValueError(problem) from None
            _log.warning(problem + "; reading stopped there")
            break
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
                default_pbc = fallback_cell()
            pbc = copy.deepcopy(default_pbc)

        for _ in range(n_atoms):
            line = file.readline()
            tokens = line.split()
            if not _looks_like_xyz_atom(line):
                break
            atom = _xyz_atom(tokens)
            try:
                atom.charge = float(tokens[4])
            except (ValueError, IndexError):
                pass
            system.append(atom)
        if len(system) < n_atoms:
            # a job still running or killed: keep the complete frames
            _log.warning(f"XYZ trajectory: frame {len(molecules) + 1} is incomplete "
                         f"({len(system)} of {n_atoms} atoms); it and anything after it "
                         f"were not loaded")
            break

        set_atom_ids(system)
        molecules.append(Molecule(atoms=system, pbc=pbc))
        pbcs.append(pbc)

        if progress_callback and file_size > 0:
            progress_callback(file.tell() / file_size)

        line = file.readline()

    return molecules, pbcs


def _is_frame_boundary(line: str) -> bool:
    """Check if a line marks the start of a new PDB/PQR frame."""
    return line[:6] == "MODEL " or line.startswith("REMARK step=")


# Edge length (Å) of the box assumed for trajectories that carry no cell. It is
# large enough that minimum-image never wraps, i.e. the system is treated as
# non-periodic. Callers should warn the user when they see it.
NO_CELL_FALLBACK_LENGTH = 1000000.0


def fallback_cell() -> PBC:
    """A huge cubic box standing in for 'no cell' in trajectory frames."""
    return PBC(NO_CELL_FALLBACK_LENGTH, NO_CELL_FALLBACK_LENGTH, NO_CELL_FALLBACK_LENGTH,
               90.0, 90.0, 90.0)


def is_fallback_cell(pbc: Optional[PBC]) -> bool:
    """True if ``pbc`` is the stand-in box from fallback_cell()."""
    if pbc is None:
        return False
    return (abs(pbc.a - NO_CELL_FALLBACK_LENGTH) < 1.0
            and abs(pbc.b - NO_CELL_FALLBACK_LENGTH) < 1.0
            and abs(pbc.c - NO_CELL_FALLBACK_LENGTH) < 1.0)


def read_pdb_trajectory(file: TextIO, progress_callback=None) -> tuple[list[Molecule], list[Optional[PBC]]]:
    """Read every frame of a multi-frame PDB/PQR.

    Frames start at ``MODEL``/``REMARK step=`` records or after an
    ``END``/``ENDMDL`` that closes a block of atoms. Frames without a cell get
    the shared fallback_cell() so minimum-image code works (non-periodic).
    """
    molecules: list[Molecule] = []
    pbcs: list[Optional[PBC]] = []

    all_lines = file.readlines()
    total_lines = len(all_lines)

    frame_lines: list[str] = []
    atoms_in_frame = 0

    def _flush() -> None:
        nonlocal frame_lines, atoms_in_frame
        if atoms_in_frame > 0:
            buf = _io.StringIO("".join(frame_lines))
            system, pbc = read_pdb(buf)
            if system:
                molecules.append(Molecule(atoms=system, pbc=pbc))
                pbcs.append(pbc)
            frame_lines = []
        else:
            # A cell written before MODEL (CRYST1 per frame, or once at the top as VMD
            # and MDAnalysis do) belongs to the frame that follows; don't drop it.
            frame_lines = [ln for ln in frame_lines if ln.startswith(("CRYST1", "REMARK BOX"))]
        atoms_in_frame = 0

    for li, line in enumerate(all_lines):
        if _is_frame_boundary(line):
            _flush()
            if progress_callback and li % 2000 == 0:
                progress_callback(li / total_lines)
        frame_lines.append(line)
        if _is_atom_record(line):
            atoms_in_frame += 1
        elif _is_frame_end(line):
            _flush()

    _flush()  # last frame

    if progress_callback:
        progress_callback(1.0)

    for idx, pbc in enumerate(pbcs):
        if pbc is None:
            pbcs[idx] = fallback_cell() if idx == 0 else copy.copy(pbcs[0])
            molecules[idx].pbc = pbcs[idx]

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
    if suffix == ".traj":
        return "ase_traj"
    if suffix == ".out":
        return "crystal_out"
    raise ValueError(f"Unsupported file format: {suffix}")


def read_file(
    filepath: str,
    remove_overlaps: bool = True,
    on_message: Optional[Callable[[str], None]] = None,
) -> Molecule:
    """Read a single-structure file into a Molecule.

    Like the original pdb_wizard, duplicate atoms sitting on top of each other
    (closer than ``geometry.OVERLAP_CUTOFF`` Å under minimum image) are deleted
    on load when the file carries a cell. Each deletion is reported through
    ``on_message``; when no callback is given the messages go to the
    ``pdb_wizard.io`` logger at WARNING level (stderr by default) so the
    deletion is never silent. Pass ``remove_overlaps=False`` to keep every atom.
    """
    ft = detect_filetype(filepath)
    load_notes: list[str] = []  # what the reader had to assume; kept on the Molecule
    if ft == "pdb":
        with open(filepath) as f:
            system, pbc = read_pdb(f, filepath, load_notes)
    elif ft == "xyz":
        with open(filepath) as f:
            system, pbc = read_xyz(f, filepath, load_notes)
    elif ft == "zmat":
        system, pbc = read_zmat(filepath)
    elif ft == "cif":
        with open(filepath) as f:
            system, pbc = read_cif(f, filepath, load_notes)
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
    elif ft == "crystal_out":
        with open(filepath, errors="replace") as f:
            system, pbc = read_crystal_output(f)
    elif ft == "dcd":
        raise ValueError(
            "DCD files are trajectories; open them with read_file_trajectory "
            "(or the TUI / --msd / --rdf)."
        )
    elif ft == "ase_traj":
        raise ValueError(
            "ASE .traj files are trajectories; open them with read_file_trajectory "
            "(or the viewers / --msd / --rdf)."
        )
    else:
        raise ValueError(f"Unknown filetype: {ft}")

    mol = Molecule(atoms=system, pbc=pbc)
    if load_notes:
        mol.load_notes = load_notes  # read by checks.load_notes and shown by the front ends
    if remove_overlaps:
        mol, messages = remove_overlapping_atoms(mol)
        report = on_message if on_message is not None else _log.warning
        for message in messages:
            report(message)
    if load_notes:
        mol.load_notes = load_notes
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
    elif ft == "ase_traj":
        is_traj = check_ase_trajectory(filepath)
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
    elif ft == "ase_traj":
        mols, _ = read_ase_trajectory(filepath)
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


def _looks_like_resp(lines: list[str]) -> bool:
    """A CP2K RESP charge file whatever it is called: a 'RESP charges' title, or atom lines
    that begin with the word RESP."""
    head = [ln.strip().lower() for ln in lines[:6] if ln.strip()]
    return any(ln.startswith("resp") for ln in head)


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
    if path.lower().endswith(".resp") or _looks_like_resp(lines):
        return parse_resp_charges(lines)
    while skip_last and lines and not lines[-1].strip():
        lines.pop()  # a trailing blank line is not the footer
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


def _sibling_topologies(dcd_path: str) -> list:
    """Every candidate topology next to the DCD, in the order of _find_sibling_topology:
    same-stem files first, then the other .pdb / .pqr / .xyz files of the folder."""
    p = Path(dcd_path)
    found = [str(p.with_suffix(ext)) for ext in (".pdb", ".pqr", ".ent", ".xyz")
             if p.with_suffix(ext).exists()]
    for ext in (".pdb", ".pqr", ".xyz"):
        for sibling in sorted(p.parent.glob(f"*{ext}")):
            if str(sibling) not in found:
                found.append(str(sibling))
    return found


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
    # (the first sibling with the right atom count, not just the first sibling)
    topology_paths = [topology_file] if topology_file else _sibling_topologies(filepath)

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
        elements: list[str] = ["C"] * n_atoms
        topology_missing = True
        for topology_path in topology_paths:
            els = _read_topology_elements(topology_path, n_atoms)
            if els is not None:
                elements = els
                topology_missing = False
                break

        # Frames
        frames: list[Molecule] = []
        pbcs: list[Optional[PBC]] = []
        coords_size = n_atoms * 4
        # A header frame count of 0 (never updated) means: read what is there.
        n_expected = hdr["n_frames"] if hdr["n_frames"] > 0 else None
        fi = -1
        while n_expected is None or fi + 1 < n_expected:
            fi += 1
            frame_start = fh.tell()
            try:
                # A job still running or killed leaves a short last frame, or fewer
                # frames than the header says: keep the complete ones.
                for _ in range(4 if hdr["has_box"] else 3):
                    _read_record(fh, endian)
            except EOFError:
                break
            fh.seek(frame_start)
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

    if not frames and hdr["n_frames"] > 0:
        raise ValueError("The DCD file holds no complete frame")
    return frames, pbcs


def check_dcd_trajectory(filepath: str) -> bool:
    """Return True if the file is a parseable DCD with at least one frame.
    Does not load full data — just probes the header."""
    try:
        with open(filepath, "rb") as fh:
            endian = _detect_endian(fh)
            header_payload = _read_record(fh, endian)
            hdr = _parse_header(header_payload, endian)
            if hdr["n_frames"] >= 1:
                return True
            # frame count never written to the header: is there data after it?
            _read_record(fh, endian)  # titles
            _read_record(fh, endian)  # atom count
            return len(fh.read(4)) == 4
    except Exception:
        return False

# ======================================================================
# Module: ase_traj
# ======================================================================
"""Reader for ASE trajectory files (``.traj``), without needing ASE.

ASE writes trajectories in its ULM container: a 48-byte header (the magic ``- of Ulm``, a
16-byte tag, then the format version, the number of items and the position of the offset
table as little-endian 64-bit integers), a table of item offsets, and for each frame a
length-prefixed JSON document. Large arrays are not in the JSON: a key ending in a dot
(``positions.``) holds ``{"ndarray": [shape, dtype, offset]}`` pointing at raw
little-endian data elsewhere in the file. The first frame carries what does not change
(atomic numbers, periodic flags); later frames repeat only what does (positions, cell).

Read here: elements, positions, the cell (when the system is periodic) and, when the
trajectory was written with a calculator that reports them, per-atom charges.
"""



import numpy as np


ASE_TRAJ_MAGIC = b"- of Ulm"
_ASE_HEADER = struct.Struct("<8s16sqqq")


def _ase_header(f) -> tuple:
    """(tag, version, number of items, offset of the offset table), or raise ValueError."""
    raw = f.read(_ASE_HEADER.size)
    if len(raw) < _ASE_HEADER.size or raw[:8] != ASE_TRAJ_MAGIC:
        raise ValueError("not an ASE trajectory (the file does not start with '- of Ulm')")
    _, tag, version, n_items, table = _ASE_HEADER.unpack(raw)
    if n_items < 0 or table < 0:
        raise ValueError("ASE trajectory header is damaged")
    return tag.decode("ascii", "replace").strip(), version, n_items, table


def check_ase_trajectory(filepath: str) -> bool:
    """True when the file is an ASE trajectory with at least one frame."""
    try:
        with open(filepath, "rb") as f:
            return _ase_header(f)[2] > 0
    except (OSError, ValueError):
        return False


def _ase_array(f, ref) -> np.ndarray:
    shape, dtype, offset = ref["ndarray"]
    count = int(np.prod(shape)) if shape else 1
    dt = np.dtype(dtype).newbyteorder("<")
    f.seek(offset)
    data = f.read(count * dt.itemsize)
    if len(data) < count * dt.itemsize:
        raise ValueError("ASE trajectory ends in the middle of an array (still being written?)")
    return np.frombuffer(data, dtype=dt, count=count).reshape(shape)


def _ase_value(f, item: dict, key: str):
    """The value stored under ``key`` (inline) or ``key.`` (a reference to raw data)."""
    if key + "." in item:
        return _ase_array(f, item[key + "."])
    if key in item:
        return np.asarray(item[key])
    return None


def _ase_item(f, offset: int) -> dict:
    f.seek(offset)
    raw = f.read(8)
    if len(raw) < 8:
        raise ValueError("ASE trajectory ends before a frame (still being written?)")
    (length,) = struct.unpack("<q", raw)
    data = f.read(length)
    if length < 0 or len(data) < length:
        raise ValueError("ASE trajectory ends in the middle of a frame (still being written?)")
    return json.loads(data.decode("utf-8"))


def _ase_cell(cell, periodic) -> tuple:
    """(PBC or None, 3x3 lattice or None) for an ASE cell. A cell counts only when the
    system is periodic in some direction and the three vectors span a volume."""
    if cell is None or periodic is None or not np.any(periodic):
        return None, None
    lattice = np.asarray(cell, dtype=float).reshape(3, 3)
    if abs(np.linalg.det(lattice)) < 1e-9:
        return None, None
    a, b, c = (float(np.linalg.norm(v)) for v in lattice)

    def angle(u, v) -> float:
        cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))

    pbc = PBC(a, b, c, angle(lattice[1], lattice[2]), angle(lattice[0], lattice[2]),
              angle(lattice[0], lattice[1]))
    return pbc, lattice


def read_ase_trajectory(filepath: str, progress_callback=None) -> tuple:
    """Read an ASE ``.traj`` file. Returns (frames, cells), like the other trajectory
    readers: frames without a periodic cell get the stand-in box and ``None`` in ``cells``.

    Frames that were only partly written (a run still in progress) are left out; the
    complete ones before them are returned. Raises ValueError when no frame can be read.
    """
    frames: List[Molecule] = []
    cells: List[Optional[PBC]] = []
    with open(filepath, "rb") as f:
        _, _, n_items, table = _ase_header(f)
        f.seek(table)
        raw = f.read(8 * n_items)
        if len(raw) < 8 * n_items:
            raise ValueError("ASE trajectory ends before its table of frames")
        offsets = struct.unpack("<%dq" % n_items, raw)

        numbers = None
        periodic = None
        symbols: List[str] = []
        for index, offset in enumerate(offsets):
            try:
                item = _ase_item(f, offset)
                frame_numbers = _ase_value(f, item, "numbers")
                if frame_numbers is not None:
                    numbers = np.asarray(frame_numbers, dtype=int)
                    symbols = [get_element_by_number(int(z)).symbol for z in numbers]
                if "pbc" in item or "pbc." in item:
                    periodic = np.asarray(_ase_value(f, item, "pbc"), dtype=bool)
                positions = _ase_value(f, item, "positions")
                if numbers is None or positions is None:
                    raise ValueError("ASE trajectory frame has no atoms")
                positions = np.asarray(positions, dtype=float).reshape(-1, 3)
                if len(positions) != len(numbers):
                    raise ValueError("ASE trajectory frame does not match its atom list")
                pbc, lattice = _ase_cell(_ase_value(f, item, "cell"), periodic)
                charges = _ase_value(f, item.get("calculator.", item.get("calculator", {})),
                                     "charges")
            except ValueError:
                if frames:  # a run that is still writing: keep what is complete
                    break
                raise
            if pbc is not None:
                # ASE allows any orientation of the cell; this program keeps a along x and
                # b in the xy plane, so the atoms are carried over by fractional coordinates
                positions = (positions @ np.linalg.inv(lattice)) @ pbc.basis_matrix
            atoms = [Atom(float(x), float(y), float(z), symbol)
                     for (x, y, z), symbol in zip(positions, symbols)]
            if charges is not None and len(charges) == len(atoms):
                for atom, q in zip(atoms, np.asarray(charges, dtype=float)):
                    atom.charge = float(q)
            set_atom_ids(atoms)
            frames.append(Molecule(atoms=atoms, pbc=pbc if pbc is not None else fallback_cell()))
            cells.append(pbc)
            if progress_callback is not None:
                progress_callback((index + 1) / max(n_items, 1))
    if not frames:
        raise ValueError("ASE trajectory has no frames")
    return frames, cells

# ======================================================================
# Module: trajectory
# ======================================================================
"""Multi-frame files (PDB, XYZ, DCD) for the viewers."""





def is_trajectory_file(filepath: str) -> bool:
    """True for multi-frame PDB/XYZ files and readable DCD and ASE .traj files."""
    ft = detect_filetype(filepath)
    if ft == "pdb":
        return check_pdb_trajectory(filepath)
    if ft == "xyz":
        return check_xyz_trajectory(filepath)
    if ft == "dcd":
        return check_dcd_trajectory(filepath)
    if ft == "ase_traj":
        return check_ase_trajectory(filepath)
    return False


def read_trajectory(filepath: str,
                    progress: Optional[Callable[[float], None]] = None) -> List[Molecule]:
    """All frames of a PDB, XYZ, DCD or ASE .traj trajectory."""
    ft = detect_filetype(filepath)
    if ft == "dcd":
        mols, _ = read_dcd_trajectory(filepath, progress_callback=progress)
        return mols
    if ft == "ase_traj":
        mols, _ = read_ase_trajectory(filepath, progress_callback=progress)
        return mols
    with open(filepath) as f:
        if ft == "pdb":
            mols, _ = read_pdb_trajectory(f, progress_callback=progress)
        else:
            mols, _ = read_xyz_trajectory(f, progress_callback=progress)
    return mols


def trajectory_caveats(frames: Optional[list]) -> List[str]:
    """Warnings about anything the reader had to assume or fake."""
    if not frames:
        return []
    first = frames[0]
    notes = []
    # DCD without a topology file: every atom is a carbon placeholder.
    if getattr(first, "_dcd_topology_missing", False):
        notes.append("DCD has no element info. All atoms loaded as carbon. "
                     "Place a sibling .pdb or .xyz with the same atom count "
                     "next to the .dcd file.")
    # No cell in the file: the reader substituted a huge non-periodic box.
    if is_fallback_cell(first.pbc):
        notes.append(f"Trajectory has no cell information; assuming a "
                     f"{NO_CELL_FALLBACK_LENGTH:.0f} Å box (non-periodic). "
                     "Use Update Cell if the system is periodic.")
    return notes


def frame_bonds(frames: List[Molecule], idx: int) -> None:
    """Give frame ``idx`` bonds: frame 0's when the topology matches, else detect them."""
    mol = frames[idx]
    if mol.bonds or not mol.atoms:
        return
    if frames[0].bonds and len(frames[0].atoms) == len(mol.atoms):
        mol.bonds = frames[0].bonds
    else:
        mol.detect_bonds()

# ======================================================================
# Module: geometry
# ======================================================================
"""Geometry analysis functions (bonds, angles, contacts, overlaps)."""


import numpy as np


# Two atoms closer than this (Å, minimum image) are treated as duplicates.
OVERLAP_CUTOFF = 0.05
# Upper bound on (rows x atoms) distance entries evaluated per chunk.
_OVERLAP_CHUNK_BUDGET = 4_000_000


def find_overlapping_pairs(
    mol: Molecule, cutoff: float = OVERLAP_CUTOFF,
) -> list[tuple[int, int]]:
    """Return index pairs (i, j), i < j, of atoms within ``cutoff`` of each other.

    Uses the minimum-image convention, so requires ``mol.pbc``. Pairs are
    sorted by (i, j). Distances are evaluated in row chunks so the memory
    footprint stays bounded for large frameworks.
    """
    pbc = mol.pbc
    n = len(mol.atoms)
    if pbc is None or n < 2:
        return []
    coords = np.array([a.x for a in mol.atoms], dtype=float)
    recip = pbc.reciprocal_basis_matrix
    basis = pbc.basis_matrix
    cut2 = cutoff * cutoff
    # Keep each chunk's (chunk x n x 3) intermediates to a few million floats.
    chunk = max(1, _OVERLAP_CHUNK_BUDGET // n)
    pairs: list[tuple[int, int]] = []
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        dx = coords[start:stop, None, :] - coords[None, :, :]
        frac = dx @ recip
        wrapped = dx - np.round(frac) @ basis
        r2 = np.einsum("ijk,ijk->ij", wrapped, wrapped)
        ii, jj = np.nonzero(r2 < cut2)
        ii = ii + start
        keep = ii < jj
        pairs.extend(zip(ii[keep].tolist(), jj[keep].tolist()))
    return pairs


def _same_movable_molecule(a: Atom, b: Atom) -> bool:
    return (a.mpmc_flag == "M" and b.mpmc_flag == "M" and a.mpmc_mol is not None
            and a.mpmc_mol == b.mpmc_mol)


def _pairs_within(mol: Molecule, cutoff: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(i, j, r) for minimum-image pairs i < j closer than ``cutoff``, in (i, j) order.

    Evaluated in row chunks, so memory stays bounded where an all-pairs matrix would take
    N^2 (gigabytes for a few thousand atoms).
    """
    pbc = mol.pbc
    n = len(mol.atoms)
    empty = (np.zeros(0, dtype=int), np.zeros(0, dtype=int), np.zeros(0))
    if pbc is None or n < 2:
        return empty
    coords = np.array([a.x for a in mol.atoms], dtype=float)
    recip, basis = pbc.reciprocal_basis_matrix, pbc.basis_matrix
    chunk = max(1, _OVERLAP_CHUNK_BUDGET // n)
    out_i, out_j, out_r = [], [], []
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        dx = coords[start:stop, None, :] - coords[None, :, :]
        wrapped = dx - np.round(dx @ recip) @ basis
        r = np.linalg.norm(wrapped, axis=2)
        ii, jj = np.nonzero(r < cutoff)
        keep = ii + start < jj
        out_i.append(ii[keep] + start)
        out_j.append(jj[keep])
        out_r.append(r[ii[keep], jj[keep]])
    if not out_i:
        return empty
    return np.concatenate(out_i), np.concatenate(out_j), np.concatenate(out_r)


def remove_overlapping_atoms(
    mol: Molecule, cutoff: float = OVERLAP_CUTOFF,
) -> tuple[Molecule, list[str]]:
    """Delete duplicate atoms that sit on top of each other.

    For each cluster of overlapping atoms the first one in file order is kept
    and the later ones are removed, matching the behaviour of the original
    pdb_wizard. Returns the molecule and one human-readable message per
    deleted atom. Nothing is printed here; callers decide how to report.
    """
    if mol.pbc is None:
        return mol, []
    system = mol.atoms
    set_atom_ids(system)
    pairs = find_overlapping_pairs(mol, cutoff)
    if not pairs:
        return mol, []

    removed: set[int] = set()
    messages: list[str] = []
    for i, j in pairs:
        if i in removed or j in removed:
            continue
        if _same_movable_molecule(system[i], system[j]):
            continue  # sites of one sorbate model (BSS H2 sits 0.04 A apart) are not copies
        removed.add(j)
        a, b = system[i], system[j]
        messages.append(
            f"Deleting overlapping atoms\n"
            f"{a.element.symbol:>3} {a.id:>5} --- {b.element.symbol:>3} {b.id:>5}\n"
            f" {a.x}\n {b.x}"
        )

    mol.atoms = [atom for k, atom in enumerate(system) if k not in removed]
    set_atom_ids(mol.atoms)
    return mol, messages


def overlap_detector(mol: Molecule) -> Molecule:
    """Classic-menu entry point: remove overlapping atoms and print what was deleted."""
    mol, messages = remove_overlapping_atoms(mol)
    for message in messages:
        print(message)
    return mol


def get_close_contacts(mol: Molecule) -> list[str]:
    if mol.pbc is None:
        return []
    system = mol.atoms
    set_atom_ids(system)
    vdw_radii = np.array([a.vdw for a in system], dtype=float)
    bond_radii = np.array([a.bond_r for a in system], dtype=float)
    # a contact is closer than the mean vdW radius, so nothing beyond the largest matters
    i_idx, j_idx, rs = _pairs_within(mol, float(vdw_radii.max()) if len(system) else 0.0)
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
    system = mol.atoms
    set_atom_ids(system)
    n = len(system)
    vdw_radii = np.array([a.vdw for a in system], dtype=float)
    i_idx, j_idx, rs = _pairs_within(mol, float(vdw_radii.max()) if n else 0.0)
    mixed_vdw = 0.5 * (vdw_radii[i_idx] + vdw_radii[j_idx])
    has_neighbor = np.zeros(n, dtype=bool)
    contact_mask = rs < mixed_vdw
    has_neighbor[i_idx[contact_mask]] = True
    has_neighbor[j_idx[contact_mask]] = True
    return [system[i] for i in range(n) if not has_neighbor[i]]


def remove_lone_atoms(mol: Molecule) -> tuple[Molecule, list[Atom]]:
    """Delete atoms with no vdW-contact neighbour; return (mol, deleted atoms).

    Requires a cell: without ``mol.pbc`` nothing is removed and the returned
    list is empty, so callers can distinguish "no lone atoms" from "no cell".
    """
    lone_atoms = get_lone_atoms(mol)
    if not lone_atoms:
        return mol, []
    lone = set(id(a) for a in lone_atoms)
    keep = [i for i, a in enumerate(mol.atoms) if id(a) not in lone]
    keep_set = set(keep)
    old_to_new = {old: new for new, old in enumerate(keep)}
    mol.atoms = [mol.atoms[i] for i in keep]
    mol.bonds = [
        (old_to_new[a], old_to_new[b])
        for a, b in mol.bonds if a in keep_set and b in keep_set
    ]
    set_atom_ids(mol.atoms)
    return mol, lone_atoms


def delete_lone_atoms(mol: Molecule) -> Molecule:
    """Delete lone atoms, returning only the molecule (see remove_lone_atoms)."""
    mol, _ = remove_lone_atoms(mol)
    return mol


def adjust_h_distances(
    mol: Molecule, second_element: str, distance: float,
) -> tuple[Molecule, list[str]]:
    """Set every bonded ``second_element``-H distance to ``distance`` Å.

    Returns the molecule and one message per adjusted pair, in the format the
    original pdb_wizard printed (``"O-H       12    13"``). Requires a cell:
    with ``mol.pbc`` None nothing is changed and no messages are returned.
    """
    if mol.pbc is None:
        return mol, []
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
    return mol, messages


def edit_h_dist(mol: Molecule, second_element: str, distance: float) -> Molecule:
    """Adjust X-H distances, returning only the molecule (see adjust_h_distances)."""
    mol, _ = adjust_h_distances(mol, second_element, distance)
    return mol


# Metals that get coordination polyhedra in the viewers.
POLYHEDRON_METALS = frozenset({
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Cs", "Ba", "La", "Ce", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Pb", "Bi",
})


def water_atom_indices(mol: Molecule) -> set[int]:
    """Indices of water molecules' atoms: an O bonded to exactly two H and nothing else.

    Uses ``mol.bonds`` as they are (call ``detect_bonds`` first).
    """
    neighbors: dict[int, list[int]] = {}
    for a, b in mol.bonds or []:
        neighbors.setdefault(a, []).append(b)
        neighbors.setdefault(b, []).append(a)
    water: set[int] = set()
    for i, atom in enumerate(mol.atoms):
        if atom.element.symbol != "O":
            continue
        nbrs = neighbors.get(i, [])
        h_nbrs = [j for j in nbrs if mol.atoms[j].element.symbol == "H"]
        if len(h_nbrs) == 2 and len(nbrs) == 2:
            water.add(i)
            water.update(h_nbrs)
    return water


def coordination_polyhedra(atoms: list, pbc=None) -> list[tuple[int, np.ndarray]]:
    """(metal index, neighbour positions) for each metal with three or more neighbours.

    Neighbours are within 2.5 covalent radii + 0.5 Å of the metal and placed at their
    minimum image around it, so a metal near a cell face gets a whole polyhedron.
    """
    if not atoms:
        return []
    positions = np.array([a.position for a in atoms])
    out = []
    for ci, atom in enumerate(atoms):
        if atom.element.symbol not in POLYHEDRON_METALS:
            continue
        cutoff = atom.element.covalent_radius * 2.5 + 0.5
        diff = positions - positions[ci]
        if pbc is not None:
            frac = diff @ pbc.reciprocal_basis_matrix
            frac -= np.round(frac)
            diff = frac @ pbc.basis_matrix
        d = np.linalg.norm(diff, axis=1)
        near = np.where(d < cutoff)[0]
        near = near[near != ci]
        if len(near) >= 3:
            out.append((ci, positions[ci] + diff[near]))
    return out

# ======================================================================
# Module: checks
# ======================================================================
"""Guard rails: things worth knowing about a structure or trajectory before simulating it.

``check_structure`` looks for the mistakes that give a wrong simulation without any error
message: a CIF read as its asymmetric unit, disordered sites kept as full atoms, missing
hydrogens, atoms on top of each other, no charges, a cell too small for the cutoff, and so
on. It never changes the structure. The terminal viewer, the desktop GUI and the classic
menu all show the same report (``structure_check_text``).
"""



import numpy as np


CHECK_ERROR, CHECK_WARNING, CHECK_NOTE = "error", "warning", "note"
_CHECK_MARKS = {CHECK_ERROR: "✗", CHECK_WARNING: "⚠", CHECK_NOTE: "•"}

CHECK_CUTOFF = 10.0  # A, the cutoff the cell is judged against (half the smallest width)
_CLOSE_CONTACT = 0.7  # A, closer than any bond: overlapping or disordered atoms
# Elements usually open-shell in MOF nodes; CP2K needs UKS and a multiplicity for them.
OPEN_SHELL_METALS = frozenset({"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Mo", "Ru",
                               "Rh", "W", "Re", "Os", "Ir", "Ce", "Gd", "Eu", "Tb", "Dy"})
_ORGANIC_NEEDS_H = 4  # carbons before "no hydrogens at all" is worth a warning


@dataclass
class CheckFinding:
    level: str  # CHECK_ERROR, CHECK_WARNING or CHECK_NOTE
    title: str
    detail: str = ""


def load_notes(mol: Optional[Molecule]) -> List[str]:
    """What the file reader had to assume (set by ``read_file``; empty for most files)."""
    return list(getattr(mol, "load_notes", None) or [])


def open_shell_metals(mol: Molecule) -> List[str]:
    """Open-shell metals really present: sorbate sites are judged by their mass, so the
    'COG'/'COE' sites of a CO2 model are not mistaken for cobalt."""
    symbols = {el.symbol for el in (site_element(a) for a in mol.atoms) if el is not None}
    return sorted(symbols & OPEN_SHELL_METALS)


def lone_atom_summary(mol: Molecule) -> str:
    """What Delete Lone Atoms would remove, e.g. '48 Cu, 2 O' ('' when nothing)."""
    if mol.pbc is None:
        return ""
    counts: dict = {}
    for atom in get_lone_atoms(mol):
        counts[atom.element.symbol] = counts.get(atom.element.symbol, 0) + 1
    return ", ".join(f"{n} {el}" for el, n in sorted(counts.items(), key=lambda kv: -kv[1]))


def lone_atom_warning(mol: Molecule) -> str:
    """A caution when Delete Lone Atoms would remove metal centres ('' otherwise)."""
    if mol.pbc is None:
        return ""
    metals = sorted({a.element.symbol for a in get_lone_atoms(mol) if a.element.symbol
                     not in ("H", "C", "N", "O", "F", "Cl", "Br", "I", "S", "P", "He", "Ne",
                             "Ar", "Kr", "Xe", "B", "Si")})
    if not metals:
        return ""
    return (f"This includes metal atoms ({', '.join(metals)}). Metal-ligand bonds are longer "
            f"than the contact distance this test uses, so framework metals can look "
            f"'lone'. Deleting them would break the framework.")


_BIO_ELEMENTS = frozenset({"H", "C", "N", "O", "S", "P", "Se", "Na", "K", "Cl", "Mg", "Ca"})


def openmm_template_warning(mol: Molecule, forcefield: str) -> str:
    """Why a stock biomolecular force field will not load this structure ('' if it may)."""
    if not any(key in forcefield.lower() for key in ("amber", "charmm", "tip", "spce")):
        return ""  # the user's own force-field file
    other = sorted({a.element.symbol for a in mol.atoms} - _BIO_ELEMENTS)
    names = {a.name.strip().upper() for a in mol.atoms}
    protein = {"CA", "N", "C", "O"} <= names
    water_only = {a.element.symbol for a in mol.atoms} <= {"H", "O"}
    if not other and (protein or water_only):
        return ""
    what = f"contains {', '.join(other)}" if other else "is not a protein or water"
    return (f"The force field ({forcefield}) only has templates for proteins, nucleic acids, "
            f"water and common ions, and this structure {what}. OpenMM will stop with "
            f"'No template found for residue'. Give your own force-field XML for a framework.")


def supercell_for_cutoff(mol: Molecule, cutoff: float = CHECK_CUTOFF) -> tuple:
    """Cells needed along a, b, c for half the smallest width to reach ``cutoff``."""
    widths = mol.pbc.perpendicular_widths()
    return tuple(int(np.ceil(2.0 * cutoff / w - 1e-9)) for w in widths)


def mc_trajectory_reason(frames: Optional[list]) -> str:
    """Why a trajectory looks like Monte Carlo output ('' if it does not).

    Monte Carlo frames have no time between them, and in grand-canonical runs molecules
    appear and disappear, so atom i is not the same molecule from frame to frame."""
    if not frames or len(frames) < 2:
        return ""
    counts = {len(m.atoms) for m in frames}
    if len(counts) > 1:
        return (f"the number of atoms changes between frames ({min(counts)} to "
                f"{max(counts)}), as in a grand-canonical (μVT) run")
    first = frames[0]
    flags = {getattr(a, "mpmc_flag", "") for a in first.atoms}
    if "M" in flags and "F" in flags:
        return "the frames are MPMC output (frozen framework 'F' and movable 'M' molecules)"
    return ""


def mc_trajectory_warning(frames: Optional[list], what: str) -> str:
    reason = mc_trajectory_reason(frames)
    if not reason:
        return ""
    return (f"⚠ This looks like a Monte Carlo trajectory: {reason}. Monte Carlo has no time, "
            f"so {what} here is not a physical result.")


def check_structure(mol: Optional[Molecule], frames: Optional[list] = None,
                    cutoff: float = CHECK_CUTOFF) -> List[CheckFinding]:
    """Everything worth a look before simulating ``mol``. Most serious first."""
    out: List[CheckFinding] = []
    if mol is None or not mol.atoms:
        return [CheckFinding(CHECK_ERROR, "No atoms loaded")]
    atoms = mol.atoms
    symbols = [a.element.symbol for a in atoms]
    n = len(atoms)

    for note in load_notes(mol):
        out.append(CheckFinding(CHECK_WARNING, "From the file reader", note))

    # massless sites of sorbate models (TIP4P's M, TIP5P's LP) are not unknown atoms
    unknown = sum(s == "X" and site_element(a) is not None for s, a in zip(symbols, atoms))
    if unknown:
        out.append(CheckFinding(CHECK_WARNING, f"{unknown} atoms have an unknown element",
                                "Their masses and force-field parameters will be wrong."))

    # ---- cell
    has_cell = mol.pbc is not None and not is_fallback_cell(mol.pbc)
    if not has_cell:
        out.append(CheckFinding(
            CHECK_WARNING, "No unit cell",
            "Periodic simulations (MPMC, CP2K) need one. Set it with Update Unit Cell."))
    else:
        half = float(mol.pbc.perpendicular_widths().min()) / 2.0
        if half < cutoff:
            na, nb, nc = supercell_for_cutoff(mol, cutoff)
            out.append(CheckFinding(
                CHECK_WARNING, f"Cell is small for a {cutoff:g} Å cutoff",
                f"Half the smallest cell width is {half:.2f} Å. A {na} × {nb} × {nc} supercell "
                f"({n * na * nb * nc:,} atoms) reaches {cutoff:g} Å."))
        else:
            out.append(CheckFinding(CHECK_NOTE, f"Cell allows a cutoff up to {half:.2f} Å"))

    # ---- atoms on top of each other
    if has_cell and n > 1:
        i, j, d = _pairs_within(mol, _CLOSE_CONTACT)
        # sites of one rigid sorbate model (and its massless sites) sit close by design
        keep = np.array([atoms[p].mass > 0 and atoms[q].mass > 0
                         and not _same_movable_molecule(atoms[p], atoms[q])
                         for p, q in zip(i, j)], dtype=bool)
        i, j, d = i[keep], j[keep], d[keep]
        if len(i):
            k = int(np.argmin(d))
            out.append(CheckFinding(
                CHECK_ERROR, f"{len(i)} pairs of atoms are closer than {_CLOSE_CONTACT} Å",
                f"Closest: {symbols[i[k]]}{i[k] + 1}-{symbols[j[k]]}{j[k] + 1} at {d[k]:.2f} Å. "
                f"This usually means disordered sites, a duplicated molecule, or the wrong "
                f"cell."))

    # ---- hydrogens
    n_c = symbols.count("C")
    if n_c >= _ORGANIC_NEEDS_H and "H" not in symbols:
        out.append(CheckFinding(
            CHECK_WARNING, "No hydrogen atoms",
            f"There are {n_c} carbons and no hydrogens. Crystal structures often leave "
            f"them out; add them before assigning charges or a force field."))

    # ---- unbonded atoms
    if has_cell:
        summary = lone_atom_summary(mol)
        if summary:
            out.append(CheckFinding(
                CHECK_NOTE, f"Atoms with no close neighbour: {summary}",
                lone_atom_warning(mol) or "These may be leftover solvent or counter-ions."))

    # ---- charges and force-field parameters (framework atoms only)
    frame_atoms = [a for a in atoms if getattr(a, "mpmc_flag", "") != "M"] or atoms
    charges = np.array([a.charge for a in frame_atoms])
    if not np.any(np.abs(charges) > 1e-9):
        out.append(CheckFinding(
            CHECK_WARNING, "No charges",
            "Every charge is zero. Load fitted charges or generate QEq charges."))
    else:
        net = float(np.array([a.charge for a in atoms]).sum())
        if abs(net) > 1e-3:
            out.append(CheckFinding(
                CHECK_WARNING, f"Net charge {net:+.4f} e",
                "A periodic simulation cell should be neutral. Check for missing "
                "counter-ions or rescale the charges."))
        else:
            out.append(CheckFinding(CHECK_NOTE, f"Charges present, net {net:+.4f} e"))
    eps = np.array([a.epsilon for a in frame_atoms])
    n_eps = int((eps > 0).sum())
    if n_eps == 0:
        out.append(CheckFinding(CHECK_WARNING, "No force-field parameters",
                                "Every ε is zero. Apply a force field before writing MPMC "
                                "input."))
    elif n_eps < len(frame_atoms):
        missing = sorted({a.element.symbol for a in frame_atoms if a.epsilon <= 0
                          and a.mass > 0})
        if missing:
            out.append(CheckFinding(
                CHECK_WARNING, f"No force-field parameters for {', '.join(missing)}",
                "The chosen force field has no entry for these elements."))

    # ---- things to remember later
    metals = open_shell_metals(mol)
    if metals:
        out.append(CheckFinding(
            CHECK_NOTE, f"Open-shell metals: {', '.join(metals)}",
            "For CP2K add UKS and a MULTIPLICITY; the generated input is closed-shell."))
    reason = mc_trajectory_reason(frames)
    if reason:
        out.append(CheckFinding(
            CHECK_NOTE, "Monte Carlo trajectory",
            f"{reason[0].upper()}{reason[1:]}. MSD and RMSD have no physical meaning here; "
            f"RDF, density and coordination are fine."))

    order = {CHECK_ERROR: 0, CHECK_WARNING: 1, CHECK_NOTE: 2}
    out.sort(key=lambda f: order[f.level])
    return out


def format_findings(findings: List[CheckFinding]) -> str:
    n_err = sum(f.level == CHECK_ERROR for f in findings)
    n_warn = sum(f.level == CHECK_WARNING for f in findings)
    if n_err or n_warn:
        parts = []
        if n_err:
            parts.append(f"{n_err} problem{'s' if n_err != 1 else ''}")
        if n_warn:
            parts.append(f"{n_warn} warning{'s' if n_warn != 1 else ''}")
        lines = ["Structure check: " + ", ".join(parts), ""]
    else:
        lines = ["Structure check: nothing suspicious found", ""]
    for f in findings:
        lines.append(f"{_CHECK_MARKS[f.level]} {f.title}")
        if f.detail:
            lines.append(f"    {f.detail}")
    lines += ["", "This is a set of heuristics, not a proof that the structure is right."]
    return "\n".join(lines)


def structure_check_text(mol: Optional[Molecule], frames: Optional[list] = None) -> str:
    return format_findings(check_structure(mol, frames))

# ======================================================================
# Module: crystal
# ======================================================================
"""CRYSTAL23 input and output, and the symmetry search that makes the input small.

Writing: ``write_crystal_d12`` produces a complete ``.d12`` (the group's cell-optimisation
template from common-scripts) whose geometry block follows the CRYSTAL manual's "Geometry
input for crystalline compounds":

    CRYSTAL
    0 0 0                 IFLAG IFHR IFSO: group by number, hexagonal axes, origin
                          setting (1 = first setting, for groups that have two)
    225                   space group number
    26.343                the minimal set of cell parameters for that crystal system
    7                     atoms in the asymmetric unit
    29 0.2850 0.2850 0.0  atomic number and fractional coordinates

Without a symmetry reduction the structure is written in P1 with every atom (always valid,
but the most expensive way to run CRYSTAL). With one, only the asymmetric unit of the
conventional cell is written. The reduction uses spglib, whose standard setting is the one
CRYSTAL assumes for a group given by number (unique axis b, hexagonal axes for
rhombohedral groups, centred conventional cell) except for the origin: spglib uses the
first origin setting where there are two and CRYSTAL defaults to the second, so IFSO=1 is
written for those groups.

Symmetry: ``find_symmetry`` runs spglib at a ladder of tolerances and ranks the space
groups it finds by how far the atoms have to move to fit each one exactly.

Reading: ``read_crystal_output`` takes the last geometry printed in a CRYSTAL ``.out``
(the optimised one after OPTGEOM) as a P1 structure in the primitive cell.
"""



import numpy as np


try:  # optional: only the symmetry search needs it
    import spglib as _spglib
except ImportError:  # pragma: no cover - exercised by the no-spglib tests via monkeypatch
    _spglib = None

SPGLIB_HINT = ("The symmetry search needs spglib (pip install spglib). Without it the "
               "structure can still be written in P1.")
SYM_TOLERANCES = (0.001, 0.003, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5)  # A
SYM_SUGGEST_DRIFT = 0.1  # A: the suggested group is the most symmetric one within this
_SYM_MATCH_CHUNK = 2_000_000

# The group's cell-optimisation settings (common-scripts/CRYSTAL/geo-opt.d12).
CRYSTAL_D12_TAIL = """OPTGEOM
FULLOPTG
ENDOPT
BASISSET
POB-TZVP-REV2
DFT
HSE06
XLGRID
ENDDFT
SHRINK
4 8
TOLINTEG
7 7 7 7 14
TOLDEE
7
MAXCYCLE
200
FMIXING
70
DIIS
END
"""

SPACE_GROUP_231 = ("Space group 231 is not in the International Tables: there are only 230. "
                   "Number 231 is the Space group at NC State, point group Brian, with no "
                   "known symmetry operations except showing up to group meeting.")


class SymmetryUnavailable(RuntimeError):
    """spglib is not installed."""


@dataclass
class SymmetryCandidate:
    """One space group found by the search."""

    number: int
    symbol: str
    tolerance: float  # the smallest tolerance (A) at which it was found
    max_drift: float  # the furthest any atom moves to fit the group exactly (A)
    rms_drift: float
    n_operations: int  # symmetry operations of the cell that was searched
    n_asymmetric: int  # atoms CRYSTAL needs in the input
    suggested: bool = False


@dataclass
class SymmetryReduction:
    """The asymmetric unit of the conventional cell, ready for a CRYSTAL input."""

    number: int
    symbol: str
    tolerance: float
    cell: tuple  # a, b, c, alpha, beta, gamma of the conventional cell
    numbers: List[int] = field(default_factory=list)  # atomic numbers
    positions: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))  # fractional
    n_conventional: int = 0  # atoms in the conventional cell
    max_drift: float = 0.0
    origin_flag: int = 0  # CRYSTAL's IFSO: 1 when the cell is in the first origin setting


# ----------------------------------------------------------------------
# Symmetry search
# ----------------------------------------------------------------------

def has_spglib() -> bool:
    return _spglib is not None


def space_group_name(number: int) -> str:
    """'F m -3 m'-style symbol for a space group number; 231 has its own entry."""
    if number == 231:
        return SPACE_GROUP_231
    if not 1 <= number <= 230:
        raise ValueError("Space groups are numbered 1 to 230")
    if _spglib is None:
        raise SymmetryUnavailable(SPGLIB_HINT)
    for hall in range(1, 531):
        info = _spglib.get_spacegroup_type(hall)
        if _field(info, "number") == number:
            return str(_field(info, "international_short"))
    raise ValueError(f"No space group {number}")  # pragma: no cover


def _field(obj, name: str):
    """spglib returned dicts before 2.5 and objects after."""
    return obj[name] if isinstance(obj, dict) else getattr(obj, name)


def _real_atoms(mol: Molecule) -> list:
    """(atom, atomic number) for every atom with a nucleus. Massless model sites (water M,
    H2 centres) are left out, and sorbate sites get their element from their mass, not
    their MPMC name (see ``io.site_element``)."""
    out = []
    for atom in mol.atoms:
        element = site_element(atom)
        if element is not None and element.symbol != "X" and element.atomic_number > 0:
            out.append((atom, element.atomic_number))
    return out


def _spglib_cell(mol: Molecule) -> tuple:
    if mol.pbc is None or is_fallback_cell(mol.pbc):
        raise ValueError("A symmetry search needs a unit cell")
    atoms = _real_atoms(mol)
    if not atoms:
        raise ValueError("No atoms to search")
    cart = np.array([a.x for a, _ in atoms], dtype=float)
    frac = cart @ mol.pbc.reciprocal_basis_matrix
    frac -= np.floor(frac)
    numbers = [z for _, z in atoms]
    return np.array(mol.pbc.basis_matrix, dtype=float), frac, numbers


def _dataset(cell: tuple, tolerance: float):
    if _spglib is None:
        raise SymmetryUnavailable(SPGLIB_HINT)
    return _spglib.get_symmetry_dataset(cell, symprec=tolerance)


def _drifts(cell: tuple, ds) -> np.ndarray:
    """How far (A) each input atom is from its place in spglib's idealised cell.

    spglib gives x_std = P x + p for the input atoms; the idealised structure has an atom
    of the same element at (or a lattice vector from) that point, moved by the drift."""
    _, frac, numbers = cell
    P = np.array(_field(ds, "transformation_matrix"), dtype=float)
    p = np.array(_field(ds, "origin_shift"), dtype=float)
    moved = frac @ P.T + p
    std_lat = np.array(_field(ds, "std_lattice"), dtype=float)
    std_pos = np.array(_field(ds, "std_positions"), dtype=float)
    std_types = np.array(_field(ds, "std_types"))
    numbers = np.array(numbers)
    out = np.zeros(len(moved))
    for z in np.unique(numbers):
        rows = np.where(numbers == z)[0]
        targets = std_pos[std_types == z]
        if len(targets) == 0:  # cannot happen for a consistent dataset
            out[rows] = np.inf
            continue
        step = max(1, _SYM_MATCH_CHUNK // max(len(targets), 1))
        for start in range(0, len(rows), step):
            part = rows[start:start + step]
            d = moved[part][:, None, :] - targets[None, :, :]
            d -= np.round(d)
            dist = np.linalg.norm(d @ std_lat, axis=2)
            out[part] = dist.min(axis=1)
    return out


def _asymmetric_unit(ds, tolerance: float) -> tuple:
    """(atomic numbers, fractional positions) of one atom per orbit of the idealised
    conventional cell, and the size of that cell."""
    std_lat = np.array(_field(ds, "std_lattice"), dtype=float)
    std_pos = np.array(_field(ds, "std_positions"), dtype=float)
    std_types = [int(z) for z in _field(ds, "std_types")]
    inner = _dataset((std_lat, std_pos, std_types), min(tolerance, 1e-3))
    if inner is None or int(_field(inner, "number")) != int(_field(ds, "number")):
        inner = _dataset((std_lat, std_pos, std_types), tolerance)
    if inner is None or int(_field(inner, "number")) != int(_field(ds, "number")):
        raise ValueError("The idealised cell did not keep its space group; try another "
                         "tolerance")
    reps = sorted(set(int(i) for i in _field(inner, "equivalent_atoms")))
    pos = std_pos[reps] - np.floor(std_pos[reps] + 1e-9)
    pos[np.abs(pos) < 5e-10] = 0.0
    return [std_types[i] for i in reps], pos, len(std_types)


def _cell_parameters(lattice: np.ndarray) -> tuple:
    a, b, c = (float(np.linalg.norm(v)) for v in lattice)

    def angle(u, v) -> float:
        cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))

    return (a, b, c, angle(lattice[1], lattice[2]), angle(lattice[0], lattice[2]),
            angle(lattice[0], lattice[1]))


def find_symmetry(mol: Molecule, tolerances=SYM_TOLERANCES) -> List[SymmetryCandidate]:
    """The space groups spglib finds as the tolerance is loosened, most symmetric last.

    Each group is listed once, at the smallest tolerance that finds it, with the atomic
    drift needed to make the structure fit it exactly."""
    cell = _spglib_cell(mol)
    found: List[SymmetryCandidate] = []
    seen = set()
    for tol in tolerances:
        ds = _dataset(cell, tol)
        if ds is None:
            continue
        number = int(_field(ds, "number"))
        if number in seen:
            continue
        seen.add(number)
        drift = _drifts(cell, ds)
        try:
            n_asym = len(_asymmetric_unit(ds, tol)[0])
        except ValueError:
            continue
        found.append(SymmetryCandidate(
            number, str(_field(ds, "international")), tol, float(drift.max()),
            float(np.sqrt((drift ** 2).mean())), len(_field(ds, "rotations")), n_asym))
    ok = [c for c in found if c.max_drift <= SYM_SUGGEST_DRIFT]
    if ok:
        max(ok, key=lambda c: (c.n_operations, -c.max_drift)).suggested = True
    return found


def reduce_symmetry(mol: Molecule, tolerance: float) -> SymmetryReduction:
    """The asymmetric unit of the conventional cell at this tolerance."""
    cell = _spglib_cell(mol)
    ds = _dataset(cell, tolerance)
    if ds is None:
        raise ValueError(f"No symmetry found at {tolerance:g} Å")
    numbers, positions, n_conv = _asymmetric_unit(ds, tolerance)
    lattice = np.array(_field(ds, "std_lattice"), dtype=float)
    return SymmetryReduction(int(_field(ds, "number")), str(_field(ds, "international")),
                             tolerance, _cell_parameters(lattice), numbers, positions, n_conv,
                             float(_drifts(cell, ds).max()),
                             _origin_flag(int(_field(ds, "number")), str(_field(ds, "choice"))))


def _origin_flag(number: int, choice: str) -> int:
    """CRYSTAL's IFSO for spglib's standard setting, after checking the setting is one
    CRYSTAL assumes for a group given by number.

    spglib's ``choice`` is '1' or '2' for the groups with two origins (it standardises to
    the first; CRYSTAL's default, IFSO=0, is the second, so IFSO=1 is written), 'H' for
    rhombohedral groups on hexagonal axes, and starts with 'b' for monoclinic unique
    axis b."""
    choice = choice.strip()
    if 3 <= number <= 15 and not choice.startswith("b"):
        raise ValueError(f"Monoclinic setting '{choice}' is not the unique-axis-b setting "
                         f"CRYSTAL expects")
    if choice == "R":
        raise ValueError("Rhombohedral axes are not supported; expected hexagonal axes")
    return 1 if choice == "1" else 0


def reduction_for_group(mol: Molecule, number: int,
                        tolerances=SYM_TOLERANCES) -> SymmetryReduction:
    """The reduction for a space group picked from ``find_symmetry``'s list."""
    if number == 231:
        raise ValueError(SPACE_GROUP_231)
    for cand in find_symmetry(mol, tolerances):
        if cand.number == number:
            return reduce_symmetry(mol, cand.tolerance)
    raise ValueError(f"Space group {number} was not found for this structure at any "
                     f"tolerance up to {max(tolerances):g} Å")


def symmetry_report(mol: Molecule) -> str:
    """The ranking shown by Find Symmetry in every front end."""
    if _spglib is None:
        return SPGLIB_HINT
    try:
        found = find_symmetry(mol)
    except ValueError as e:
        return str(e)
    n = len(_real_atoms(mol))
    lines = [f"Space groups found for {n} atoms, loosening the tolerance step by step:", "",
             f"{'group':>5}  {'symbol':<12} {'tolerance':>9}  {'max drift':>9}  "
             f"{'rms drift':>9}  {'atoms in input':>14}"]
    for c in found:
        lines.append(f"{c.number:>5}  {c.symbol:<12} {c.tolerance:>7.3f} Å  "
                     f"{c.max_drift:>7.4f} Å  {c.rms_drift:>7.4f} Å  {c.n_asymmetric:>14}"
                     + ("   <- suggested" if c.suggested else ""))
    lines += ["",
              "tolerance: the smallest search tolerance that finds the group.",
              "drift: how far atoms move when the structure is made to fit it exactly.",
              f"suggested: the most symmetric group within {SYM_SUGGEST_DRIFT:g} Å of drift. "
              f"A large drift means the",
              "symmetry is only approximate (or a guest breaks it); a higher group that "
              "needs it is not real.",
              "", "Use the group number when saving a CRYSTAL23 .d12 to write only its "
              "asymmetric unit."]
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------

def crystal_lattice_line(number: int, cell: tuple) -> str:
    """The minimal set of cell parameters CRYSTAL reads for a space group number."""
    a, b, c, alpha, beta, gamma = cell

    def fmt(*values) -> str:
        return " ".join(f"{v:.6f}" for v in values)

    if number >= 195:  # cubic
        return fmt(a)
    if number >= 75:  # tetragonal, trigonal and hexagonal (hexagonal axes)
        return fmt(a, c)
    if number >= 16:  # orthorhombic
        return fmt(a, b, c)
    if number >= 3:  # monoclinic, unique axis b
        return fmt(a, b, c, beta)
    return fmt(a, b, c, alpha, beta, gamma)


def crystal_geometry(mol: Molecule, reduction: Optional[SymmetryReduction] = None) -> str:
    """The geometry block, from ``CRYSTAL`` to the last atom (no END)."""
    unknown = []
    for atom in mol.atoms:
        element = site_element(atom)  # None: a massless model site, left out by design
        if element is not None and (element.symbol == "X" or element.atomic_number <= 0):
            unknown.append(atom.name)
    if unknown:
        raise ValueError(
            "A CRYSTAL input needs an atomic number for every atom; these atom names have "
            f"no known element: {', '.join(sorted(set(unknown)))}")
    if reduction is None:
        if mol.pbc is None or is_fallback_cell(mol.pbc):
            raise ValueError("A CRYSTAL input needs a unit cell")
        atoms = _real_atoms(mol)
        if not atoms:
            raise ValueError("No atoms to write")
        pbc = mol.pbc
        frac = np.array([a.x for a, _ in atoms]) @ pbc.reciprocal_basis_matrix
        frac -= np.floor(frac + 1e-9)
        number = 1
        cell = (pbc.a, pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma)
        numbers = [z for _, z in atoms]
    else:
        number, cell = reduction.number, reduction.cell
        numbers, frac = reduction.numbers, reduction.positions
    origin = 0 if reduction is None else reduction.origin_flag
    lines = ["CRYSTAL", f"0 0 {origin}", str(number), crystal_lattice_line(number, cell),
             str(len(numbers))]
    for z, (x, y, zc) in zip(numbers, frac):
        lines.append(f"{z:<3d} {x + 0.0:15.10f} {y + 0.0:15.10f} {zc + 0.0:15.10f}")
    return "\n".join(lines)


def write_crystal_d12(mol: Molecule, out: TextIO,
                      reduction: Optional[SymmetryReduction] = None,
                      title: str = "") -> None:
    """A complete CRYSTAL23 cell-optimisation input (HSE06/POB-TZVP-REV2, the group's
    template). ``reduction`` None writes every atom in P1."""
    if not title:
        what = "P1" if reduction is None else f"space group {reduction.number}"
        title = f"PDB Wizard - cell optimization, HSE06/POB-TZVP-REV2 ({what})"
    out.write(title.splitlines()[0] + "\n")
    out.write(crystal_geometry(mol, reduction) + "\n")
    out.write(CRYSTAL_D12_TAIL)


def crystal_reduction_choice(mol: Molecule, answer: str) -> Optional[SymmetryReduction]:
    """Turn what the user typed at "space group?" into a reduction.

    Blank, 1 or P1 means no reduction (None). A number means that group from the Find
    Symmetry list. Raises ValueError with the message to show otherwise."""
    text = answer.strip().lower().replace(" ", "")
    if text in ("", "1", "p1"):
        return None
    if not text.isdigit():
        raise ValueError("Enter a space group number from Find Symmetry, or 1 for P1")
    number = int(text)
    if number == 231:
        raise ValueError(SPACE_GROUP_231)
    if not 1 <= number <= 230:
        raise ValueError("Space groups are numbered 1 to 230")
    if _spglib is None:
        raise ValueError(SPGLIB_HINT)
    return reduction_for_group(mol, number)


# ----------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------

def read_crystal_output(file: TextIO) -> tuple:
    """(atoms, PBC) of the last geometry in a CRYSTAL output: the optimised structure
    after OPTGEOM, or the input one for a single point. It is the primitive cell with all
    of its atoms, so the result is a complete P1 structure."""
    lines = file.read().splitlines()
    cell = None
    atoms_at = None
    for i, line in enumerate(lines):
        if "ATOMS IN THE ASYMMETRIC UNIT" in line and "ATOMS IN THE UNIT CELL" in line:
            # the cell printed just above belongs to this block
            for j in range(i - 1, max(i - 12, -1), -1):
                if "ALPHA" in lines[j] and "BETA" in lines[j] and "GAMMA" in lines[j]:
                    try:
                        values = [float(v) for v in lines[j + 1].split()[:6]]
                    except (ValueError, IndexError):
                        break
                    if len(values) == 6:
                        cell, atoms_at = values, i
                    break
    if cell is None or atoms_at is None:
        raise ValueError("No CRYSTAL geometry found (expected 'ATOMS IN THE ASYMMETRIC "
                         "UNIT' after a lattice-parameter table)")
    pbc = PBC(*cell)
    atoms: List[Atom] = []
    for line in lines[atoms_at + 1:]:
        tokens = line.split()
        if len(tokens) == 7 and tokens[0].isdigit() and tokens[1] in ("T", "F"):
            try:
                frac = np.array([float(t) for t in tokens[4:7]])
                z = int(tokens[2]) % 100  # "conventional" atomic numbers carry a basis tag
            except ValueError:
                break
            symbol = tokens[3].capitalize()
            if get_element(symbol).symbol == "X" or get_element(symbol).atomic_number != z:
                symbol = _symbol_for(z)
            xyz = frac @ pbc.basis_matrix
            atom = Atom(xyz[0], xyz[1], xyz[2], symbol)
            atoms.append(atom)
        elif atoms and not line.strip().startswith("*"):
            break
    if not atoms:
        raise ValueError("The CRYSTAL geometry block has no atoms")
    set_atom_ids(atoms)
    return atoms, pbc


def _symbol_for(atomic_number: int) -> str:
    return get_element_by_number(atomic_number).symbol

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


def site_mass(atom: Atom) -> float:
    """Mass of an atom: the model's own value for a movable MPMC site (massless sites are
    0, a united-atom site keeps its group mass), the element mass otherwise."""
    if atom.mpmc_flag == "M" and atom.file_mass is not None:
        return float(atom.file_mass)
    return atom.mass


def real_atoms(mol: Molecule) -> list:
    """(atom, element) for every atom that is a real atom. Sites of MPMC sorbate models
    get their true element ('COG' is carbon, not cobalt) and massless sites are left out."""
    out = []
    for atom in mol.atoms:
        element = site_element(atom)
        if element is not None:
            out.append((atom, element))
    return out


def formula_unit(mol: Molecule) -> dict[str, int]:
    counts: dict[str, int] = {}
    for _atom, element in real_atoms(mol):
        sym = element.symbol
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
    mol_ids = [a.mpmc_mol for a in mol.atoms if a.mpmc_flag == "M" and a.mpmc_mol is not None]
    id_step = max(mol_ids, default=0)
    for i in range(1, times + 1):
        for atom in mol.atoms:
            new_atom = copy.deepcopy(atom)
            new_atom.x = atom.x + i * pbc.basis_matrix[axis]
            if new_atom.mpmc_flag == "M" and new_atom.mpmc_mol is not None:
                new_atom.mpmc_mol += i * id_step  # each copy is its own molecule
            new_atoms.append(new_atom)
    mol.atoms.extend(new_atoms)

    params = [pbc.a, pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma]
    params[axis] *= (times + 1)
    pbc.update(*params)
    set_atom_ids(mol.atoms)
    # The cached bond list only covers the original atoms; drop it so the
    # geometry listings recompute instead of showing a stale, short list.
    mol.bonds = []
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
    if mol.pbc is None:
        return {"void_volume": 0.0, "total_volume": 0.0, "void_fraction": 0.0}

    pbc = mol.pbc
    sites = real_atoms(mol)
    if not sites:  # nothing in the cell: all of it is void
        return {"void_volume": pbc.volume, "total_volume": pbc.volume, "void_fraction": 1.0}
    coords = np.array([a.x for a, _ in sites])
    vdw_radii = np.array([el.vdw_radius for _, el in sites])

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
            frac_dx = dx @ recip
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
    sites = real_atoms(mol)
    if mol.pbc is None or len(sites) == 0:
        return {"surface_area": 0.0, "area_per_volume": 0.0}

    pbc = mol.pbc
    coords = np.array([a.x for a, _ in sites])              # (N, 3)
    vdw_radii = np.array([el.vdw_radius for _, el in sites])  # (N,)
    radii = vdw_radii + probe_radius                        # (N,)
    n_atoms = len(sites)

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
    # cells must be at least `cutoff` thick: in a skewed cell that is set by the
    # perpendicular widths, not by the edge lengths
    n_cells = max(1, int(float(pbc.perpendicular_widths().min()) // cutoff))
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

    total_mass = sum(site_mass(a) for a in mol.atoms)
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

    sites = real_atoms(mol)
    if mol.pbc is None or len(sites) == 0:
        return bin_centers, np.zeros(n_bins)

    pbc = mol.pbc
    coords = np.array([a.x for a, _ in sites])              # (N, 3)
    vdw_radii = np.array([el.vdw_radius for _, el in sites])  # (N,)
    n_atoms = len(sites)

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
        total_mass = sum(site_mass(a) for a in mol.atoms)
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
        total_mass = sum(site_mass(a) for a in mol.atoms)
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


def classify_periodic_bonds(positions, bonds, pbc, tol: float = 1e-6):
    """Split bonds into normal vs. periodic (wrapped across the cell).

    Bonds are detected with the minimum-image convention, so a bond may connect
    two atoms sitting on opposite faces of the cell. Drawing a straight line
    between their real coordinates streaks across the whole structure; instead
    such bonds should be drawn as stubs toward each atom's nearest periodic
    image.

    Returns ``(is_periodic, d_min)`` where ``is_periodic[k]`` is True when bond
    ``k`` wraps, and ``d_min[k]`` is its minimum-image displacement
    ``pos_j - pos_i`` (real space). With no PBC, all bonds are non-periodic.
    """
    n = len(bonds)
    if n == 0:
        return np.zeros(0, dtype=bool), np.zeros((0, 3))
    n_atoms = len(positions)
    bond_i = np.clip(np.array([b[0] for b in bonds]), 0, max(n_atoms - 1, 0))
    bond_j = np.clip(np.array([b[1] for b in bonds]), 0, max(n_atoms - 1, 0))
    d_raw = positions[bond_j] - positions[bond_i]
    if pbc is None:
        return np.zeros(n, dtype=bool), d_raw
    d_min = pbc.wrap(d_raw)
    is_periodic = np.linalg.norm(d_raw - d_min, axis=1) > tol
    return is_periodic, d_min


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
                # Periodic bonds connect atoms on opposite faces of the cell.
                # Draw them as two short stubs (one from each atom toward its
                # partner's nearest image) instead of a line straight across.
                stub_view = None
                is_periodic = None
                if snap_bonds:
                    is_periodic, d_min = classify_periodic_bonds(
                        positions, snap_bonds, pbc
                    )
                    stub_view = (0.5 * d_min) @ rot.T  # half-bond, view space
                for k, (i, j) in enumerate(snap_bonds):
                    if i >= n_snap or j >= n_snap:
                        continue  # skip bonds that reference atoms outside snapshot
                    c1, c2 = atom_colors[i], atom_colors[j]
                    if has_hl and i in hl and j in hl:
                        c1 = self._highlight_color()
                        c2 = self._highlight_color()
                    if is_periodic is not None and is_periodic[k]:
                        self.render_bond(
                            transformed[i], transformed[i] + stub_view[k], c1, c1
                        )
                        self.render_bond(
                            transformed[j], transformed[j] - stub_view[k], c2, c2
                        )
                    else:
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
        if (density_positions is not None and density_values is not None
                and len(density_positions) > 0):
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
            for ci, nbr_world in coordination_polyhedra(snap_atoms, pbc):
                nbr_view = []
                for p in nbr_world:
                    pt = rot @ (p - centroid)
                    pt[0] += pan[0]
                    pt[1] += pan[1]
                    pt[2] += camera_distance
                    nbr_view.append(pt)
                # Draw edges between all neighboring atoms (polyhedron wireframe)
                color = snap_atoms[ci].element.cpk_color
                dim = (color[0] // 2, color[1] // 2, color[2] // 2)
                for ai in range(len(nbr_view)):
                    for bj in range(ai + 1, len(nbr_view)):
                        self.render_bond(nbr_view[ai], nbr_view[bj], dim, dim)

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

# GMP orbital radii R (Å) from Rappé & Goddard; with chi and J they are the UFF "Xi",
# "Hard" (= J/2) and "Radius" columns (RDKit/Open Babel UFF.prm).
# fmt: off
QEQ_RADII: dict[str, float] = {
    "H": 0.371, "He": 1.3, "Li": 1.557, "Be": 1.24, "B": 0.822, "C": 0.759, "N": 0.715,
    "O": 0.669, "F": 0.706, "Ne": 1.768, "Na": 2.085, "Mg": 1.5, "Al": 1.201, "Si": 1.176,
    "P": 1.102, "S": 1.047, "Cl": 0.994, "Ar": 2.108, "K": 2.586, "Ca": 2.0, "Sc": 1.75,
    "Ti": 1.607, "V": 1.47, "Cr": 1.402, "Mn": 1.533, "Fe": 1.393, "Co": 1.406, "Ni": 1.398,
    "Cu": 1.434, "Zn": 1.4, "Ga": 1.211, "Ge": 1.189, "As": 1.204, "Se": 1.224, "Br": 1.141,
    "Kr": 2.27, "Rb": 2.77, "Sr": 2.415, "Y": 1.998, "Zr": 1.758, "Nb": 1.603, "Mo": 1.53,
    "Tc": 1.5, "Ru": 1.5, "Rh": 1.509, "Pd": 1.544, "Ag": 1.622, "Cd": 1.6, "In": 1.404,
    "Sn": 1.354, "Sb": 1.404, "Te": 1.38, "I": 1.333, "Xe": 2.459, "Cs": 2.984, "Ba": 2.442,
    "La": 2.071, "Ce": 1.925, "Pr": 2.007, "Nd": 2.007, "Pm": 2.0, "Sm": 1.978, "Eu": 2.227,
    "Gd": 1.968, "Tb": 1.954, "Dy": 1.934, "Ho": 1.925, "Er": 1.915, "Tm": 2.0, "Yb": 2.158,
    "Lu": 1.896, "Hf": 1.759, "Ta": 1.605, "W": 1.538, "Re": 1.6, "Os": 1.7, "Ir": 1.866,
    "Pt": 1.557, "Au": 1.618, "Hg": 1.6, "Tl": 1.53, "Pb": 1.444, "Bi": 1.514,
}
# fmt: on

# Coulomb constant: 14.4 eV*A/e^2
_KE = 14.3996
# Slater exponent scaling: zeta = LAMBDA (2n + 1) / (2 R)  (Rappé & Goddard 1991)
_LAMBDA = 0.4913
_DEFAULT_PARAMS = (5.0, 10.0, 1.2)  # chi, J, R for elements without parameters

# Short-range shielding is tabulated out to this distance (Å); beyond it J_AB = K/R to
# better than 1e-7 eV for every element pair.
_SHIELD_CUTOFF = 10.0
_SHIELD_STEP = 0.002


def _principal_n(atomic_number: int) -> int:
    """Principal quantum number of the valence shell (the Slater ns orbital used)."""
    for n, last in ((1, 2), (2, 10), (3, 18), (4, 36), (5, 54)):
        if atomic_number <= last:
            return n
    return 6


def _slater_density_ft(k: np.ndarray, zeta: float, n: int) -> np.ndarray:
    """Fourier transform of a normalized Slater ns density |phi|^2 ~ r^(2n-2) e^(-2 zeta r)."""
    a = 2.0 * zeta
    k = np.asarray(k, dtype=float)
    out = np.ones_like(k)
    nz = k > 1e-12
    kk = k[nz]
    out[nz] = a ** (2 * n + 1) * np.imag((a - 1j * kk) ** (-2 * n)) / (2 * n * kk)
    return out


def slater_coulomb(zeta_a: float, n_a: int, zeta_b: float, n_b: int,
                   r: np.ndarray) -> np.ndarray:
    """Coulomb integral (1/Å) between two Slater ns densities centred ``r`` Å apart.

    J(R) = (2/pi) Int_0^inf rhoA(k) rhoB(k) sin(kR)/(kR) dk. Tends to 1/R at long range
    and stays finite at R = 0 (5 zeta / 8 for two identical 1s densities).
    """
    k_max = 80.0 * max(zeta_a, zeta_b)  # rhoA * rhoB < 1e-11 beyond this
    k = np.linspace(0.0, k_max, 60001)
    w = _slater_density_ft(k, zeta_a, n_a) * _slater_density_ft(k, zeta_b, n_b)
    r = np.atleast_1d(np.asarray(r, dtype=float))
    out = np.empty_like(r)
    dk = k[1] - k[0]
    for i, R in enumerate(r):
        kr = k * R
        sinc = np.where(kr > 1e-12, np.sin(kr) / np.where(kr > 1e-12, kr, 1.0), 1.0)
        f = w * sinc
        out[i] = (2.0 / np.pi) * dk * (f.sum() - 0.5 * (f[0] + f[-1]))
    return out


def _erfc(x: np.ndarray) -> np.ndarray:
    """Complementary error function for numpy arrays (relative error < 1.2e-7).

    Numerical Recipes' erfcc; the package must run with numpy alone.
    """
    z = np.abs(x)
    t = 1.0 / (1.0 + 0.5 * z)
    poly = (-1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (
        -0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (
            -0.82215223 + t * 0.17087277)))))))))
    ans = t * np.exp(-z * z + poly)
    return np.where(x >= 0, ans, 2.0 - ans)


class _Shielding:
    """Tabulated short-range part s_AB(R) = 1/R - J_AB(R) (1/Å) for each pair of orbitals.

    J_AB(R) is the Coulomb integral between the two Slater densities; s vanishes
    (below 1e-7) well inside _SHIELD_CUTOFF.
    """

    def __init__(self) -> None:
        self.grid = np.arange(_SHIELD_STEP, _SHIELD_CUTOFF + _SHIELD_STEP / 2, _SHIELD_STEP)
        self._tables: dict = {}

    def _entry(self, a: tuple, b: tuple) -> tuple:
        key = (a, b) if a <= b else (b, a)
        if key not in self._tables:
            (za, na), (zb, nb) = key
            # J is smooth: integrate on a coarse grid, then interpolate onto the fine one
            coarse = np.concatenate([np.linspace(0.0, 2.0, 81),
                                     np.linspace(2.05, _SHIELD_CUTOFF, 160)])
            j = slater_coulomb(za, na, zb, nb, coarse)
            s_coarse = 1.0 / coarse[1:] - j[1:]
            s_fine = np.interp(self.grid, coarse[1:], s_coarse)
            self._tables[key] = (s_fine, float(j[0]))
        return self._tables[key]

    def s(self, a: tuple, b: tuple, r: np.ndarray) -> np.ndarray:
        """s_AB at distances r > 0 (zero beyond the table)."""
        s_fine, _ = self._entry(a, b)
        return np.interp(r, self.grid, s_fine, right=0.0)

    def j0(self, a: tuple, b: tuple) -> float:
        """J_AB(0) (1/Å): two orbitals on top of each other."""
        return self._entry(a, b)[1]


def _cell_widths(basis: np.ndarray) -> np.ndarray:
    vol = abs(np.linalg.det(basis))
    return np.array([vol / np.linalg.norm(np.cross(basis[(i + 1) % 3], basis[(i + 2) % 3]))
                     for i in range(3)])


def _coulomb_matrix(coords: np.ndarray, orbitals: list, pbc=None,
                    eta: Optional[float] = None, progress=None) -> np.ndarray:
    """Shielded Coulomb matrix (1/Å) between the atoms' Slater orbitals.

    Off-diagonal: J_ij summed over all periodic images (Ewald for the 1/R part, plus the
    short-range shielding). Diagonal: only the interaction with the atom's own images
    (zero without a cell); the QEq hardness J_i is added by the caller.
    """
    n = len(coords)
    shield = _Shielding()
    kinds = sorted(set(orbitals))
    kind_of = np.array([kinds.index(o) for o in orbitals])
    A = np.zeros((n, n))
    periodic = pbc is not None
    if periodic:
        eta = 0.35 if eta is None else eta
        r_real = max(4.0 / eta, _SHIELD_CUTOFF)  # erfc(4) < 2e-8 and s ~ 0 beyond
        basis = pbc.basis_matrix
        frac = coords @ np.linalg.inv(basis)
        widths = _cell_widths(basis)
        reach = np.ceil(r_real / widths + 0.5).astype(int)
        shifts = np.mgrid[-reach[0]:reach[0] + 1, -reach[1]:reach[1] + 1,
                          -reach[2]:reach[2] + 1].reshape(3, -1).T.astype(float)
    else:
        eta = 0.0
        r_real = np.inf
        shifts = np.zeros((1, 3))

    # ---- real space: erfc(eta r)/r - s(r) over images (1/r - s(r) without a cell) ----
    block = max(1, 2_000_000 // max(n, 1))
    two_eta_rpi = 2.0 * eta / np.sqrt(np.pi)
    for b0 in range(0, n, block):
        rows = np.arange(b0, min(b0 + block, n))
        if periodic:
            df = frac[None, :, :] - frac[rows][:, None, :]
            df -= np.round(df)
        else:
            dx = coords[None, :, :] - coords[rows][:, None, :]
        for shift in shifts:
            if periodic:
                d = (df + shift) @ basis
            else:
                d = dx
            r = np.sqrt(np.einsum("ijk,ijk->ij", d, d))
            self_term = (rows[:, None] == np.arange(n)[None, :]) & np.all(shift == 0)
            near = (r < r_real) & ~self_term
            if not near.any():
                continue
            ii, jj = np.nonzero(near)
            rr = r[ii, jj]
            val = np.empty_like(rr)
            zero = rr < 1e-8  # two sites on top of each other (e.g. MPMC dummy sites)
            pos = ~zero
            coul = 1.0 / rr[pos] if not periodic else _erfc(eta * rr[pos]) / rr[pos]
            val[pos] = coul
            ki, kj = kind_of[rows[ii]], kind_of[jj]
            for ka in range(len(kinds)):
                for kb in range(len(kinds)):
                    m = (ki == ka) & (kj == kb)
                    if not m.any():
                        continue
                    mp = m & pos
                    if mp.any():
                        val[mp] -= shield.s(kinds[ka], kinds[kb], rr[mp])
                    mz = m & zero
                    if mz.any():
                        val[mz] = shield.j0(kinds[ka], kinds[kb]) - two_eta_rpi
            np.add.at(A, (rows[ii], jj), val)
        if progress:
            progress(min(b0 + block, n) / n)

    if not periodic:
        return A

    # ---- reciprocal space and the self term ------------------------------------------
    vol = abs(np.linalg.det(basis))
    recip = 2.0 * np.pi * np.linalg.inv(basis).T  # rows: reciprocal vectors (with 2 pi)
    k_cut = 9.0 * eta  # exp(-k^2 / 4 eta^2) < 2e-9
    m_max = np.ceil(k_cut * np.linalg.norm(basis, axis=1) / (2.0 * np.pi)).astype(int)
    m = np.mgrid[-m_max[0]:m_max[0] + 1, -m_max[1]:m_max[1] + 1,
                 -m_max[2]:m_max[2] + 1].reshape(3, -1).T
    # half space: k and -k contribute equally
    half = (m[:, 0] > 0) | ((m[:, 0] == 0) & (m[:, 1] > 0)) | (
        (m[:, 0] == 0) & (m[:, 1] == 0) & (m[:, 2] > 0))
    kvec = m[half] @ recip
    k2 = np.einsum("ij,ij->i", kvec, kvec)
    keep = (k2 > 0) & (k2 < k_cut * k_cut)
    kvec, k2 = kvec[keep], k2[keep]
    weight = 2.0 * (4.0 * np.pi / vol) * np.exp(-k2 / (4.0 * eta * eta)) / k2
    chunk = max(1, 4_000_000 // max(n, 1))
    for c0 in range(0, len(kvec), chunk):
        phase = coords @ kvec[c0:c0 + chunk].T
        w = np.sqrt(weight[c0:c0 + chunk])
        C, S = np.cos(phase) * w, np.sin(phase) * w
        A += C @ C.T + S @ S.T
    A[np.diag_indices(n)] -= two_eta_rpi
    return A


def compute_qeq_charges(
    mol: Molecule,
    total_charge: float = 0.0,
    progress_callback=None,
    eta: Optional[float] = None,
) -> np.ndarray:
    """QEq charges (Rappé & Goddard, J. Phys. Chem. 1991, 95, 3358).

    Minimizes E(q) = sum_i (chi_i q_i + J_i q_i^2 / 2) + sum_{i<j} J_ij q_i q_j subject to
    sum q = total_charge. J_ij is the Coulomb integral between Slater ns orbitals with
    exponents zeta = 0.4913 (2n + 1) / (2 R) (Rappé-Goddard GMP radii R); it is shielded
    at bonding distances and tends to 14.4/R at long range. Periodic cells use an Ewald
    sum over all images. Rappé-Goddard's charge-dependent hydrogen hardness is not
    included.

    Earlier versions used a bare 14.4/R between bonded atoms, which exceeds the atoms'
    own hardness at bond lengths and inverted the signs (water came out O +0.85).

    Args:
        mol: Molecule with atoms (a cell, if any, is used periodically).
        total_charge: Total system charge (0 for neutral).
        progress_callback: Optional callable(frac: float) for progress updates.
        eta: Ewald splitting parameter (1/Å); the charges do not depend on it.

    Returns:
        Array of partial charges (electrons).
    """

    n = len(mol.atoms)
    if n == 0:
        return np.array([])
    if n > 3000:
        raise ValueError(
            f"System too large for QEq ({n} atoms). "
            "Try on a single unit cell before extending."
        )

    chi = np.zeros(n)
    hard = np.zeros(n)
    orbitals = []
    for i, atom in enumerate(mol.atoms):
        sym = atom.element.symbol
        c, j = QEQ_PARAMS.get(sym, _DEFAULT_PARAMS[:2])
        radius = QEQ_RADII.get(sym, _DEFAULT_PARAMS[2])
        n_shell = _principal_n(atom.element.atomic_number or 6)
        chi[i], hard[i] = c, j
        orbitals.append((round(_LAMBDA * (2 * n_shell + 1) / (2.0 * radius), 10), n_shell))

    pbc = mol.pbc if mol.pbc is not None and not is_fallback_cell(mol.pbc) else None
    coords = np.array([a.x for a in mol.atoms], dtype=float)

    if progress_callback:
        progress_callback(0.05)

    def _prog(frac: float) -> None:
        if progress_callback:
            progress_callback(0.05 + 0.85 * frac)

    J = _KE * _coulomb_matrix(coords, orbitals, pbc, eta, _prog)
    J[np.diag_indices(n)] += hard

    A = np.zeros((n + 1, n + 1))
    A[:n, :n] = J
    A[n, :n] = A[:n, n] = 1.0
    b = np.concatenate([-chi, [total_charge]])
    try:
        result = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    if progress_callback:
        progress_callback(1.0)
    return result[:n]


def apply_qeq_charges(mol: Molecule, progress_callback=None) -> None:
    """Compute and apply QEq charges to the framework atoms.

    Movable (M) molecules read from an MPMC file keep their model charges and
    are left out of the equalization.
    """
    framework = [a for a in mol.atoms if a.mpmc_flag != "M"]
    if len(framework) != len(mol.atoms):
        mol = Molecule(atoms=framework, pbc=mol.pbc)
    charges = compute_qeq_charges(mol, progress_callback=progress_callback)
    for i, atom in enumerate(mol.atoms):
        atom.charge = float(charges[i])

# ======================================================================
# Module: rdf
# ======================================================================
"""Radial distribution function (RDF) calculation with PBC support."""


import numpy as np



def _select(mol: Molecule, element: str, name: str | None) -> np.ndarray:
    return np.array([i for i, a in enumerate(mol.atoms)
                     if a.element.symbol == element
                     and (name is None or a.name.strip() == name)], dtype=int)


def half_cell_width(pbc) -> float:
    """Half the smallest distance between opposite cell faces: the largest distance for
    which the minimum image finds every pair. min(a, b, c) / 2 for an orthogonal cell."""
    half = float(pbc.perpendicular_widths().min()) / 2.0
    edge = min(pbc.a, pbc.b, pbc.c) / 2.0
    return edge if not np.isfinite(half) or abs(half - edge) <= 1e-9 * edge else half


# Most pair distances an image sum may evaluate. A placeholder cell (CRYST1 1 1 1 around
# thousands of atoms) would need far more; it keeps the plain minimum image, as before.
_IMAGE_SUM_MAX_WORK = 5e9


def image_shifts(pbc, r_max: float, n_pairs: int = 1) -> np.ndarray:
    """Lattice translations (in cell units, integer rows) needed besides the minimum image
    to reach every pair within ``r_max``; the zero shift is not included. Empty when
    ``r_max`` is within half the smallest cell width (or the cell is degenerate, or the
    sum over ``n_pairs`` pairs would be unreasonably long)."""
    none = np.zeros((0, 3))
    widths = pbc.perpendicular_widths()
    if not np.all(np.isfinite(widths)) or widths.min() <= 0:
        return none
    if r_max <= half_cell_width(pbc) * (1.0 + 1e-12):
        return none
    reach = np.ceil(r_max / widths + 0.5).astype(int)
    if float(np.prod(2.0 * reach + 1.0)) * n_pairs > _IMAGE_SUM_MAX_WORK:
        return none
    grid = np.mgrid[-reach[0]:reach[0] + 1, -reach[1]:reach[1] + 1,
                    -reach[2]:reach[2] + 1].reshape(3, -1).T.astype(float)
    return grid[np.any(grid != 0, axis=1)]


def compute_rdf(
    mol: Molecule,
    type1: str,
    type2: str,
    n_bins: int = 200,
    r_max: float | None = None,
    progress_callback=None,
    name1: str | None = None,
    name2: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the RDF g(r) between two atom types for a single frame.

    Args:
        mol: Molecule with PBC.
        type1: Element symbol for first atom type.
        type2: Element symbol for second atom type.
        n_bins: Number of histogram bins.
        r_max: Maximum distance (default: half the smallest cell width). Beyond that
            the pairs are summed over periodic images.
        progress_callback: Optional callable(frac: float).
        name1, name2: Optional atom names (e.g. "OW") narrowing each type.

    Returns:
        (r_values, g_r) arrays. g(r) tends to 1 at long range for an ideal gas, both for
        a type with itself (each pair counted once, hence the factor 2) and for two types.
    """
    if mol.pbc is None:
        return np.array([]), np.array([])

    pbc = mol.pbc
    if r_max is None:
        r_max = half_cell_width(pbc)

    idx1 = _select(mol, type1, name1)
    idx2 = _select(mol, type2, name2)

    if not len(idx1) or not len(idx2):
        return np.array([]), np.array([])

    coords = np.array([a.x for a in mol.atoms])
    dr = r_max / n_bins
    hist = np.zeros(n_bins)

    same_type = len(idx1) == len(idx2) and bool(np.array_equal(idx1, idx2))
    total = len(idx1)
    targets = coords[idx2]
    # rows of the i x j distance matrix per block, keeping a block near 2M pairs
    block = max(1, 2_000_000 // len(idx2))
    # beyond half the cell width one image per pair is not enough
    shifts = image_shifts(pbc, r_max, len(idx1) * len(idx2)) @ pbc.basis_matrix

    for b0 in range(0, total, block):
        rows = np.arange(b0, min(b0 + block, total))
        dx = targets[None, :, :] - coords[idx1[rows]][:, None, :]
        frac = dx @ pbc.reciprocal_basis_matrix
        wrapped = dx - np.round(frac) @ pbc.basis_matrix
        dists = np.linalg.norm(wrapped, axis=2)
        if same_type:  # each pair once: j after i
            keep = np.arange(len(idx2))[None, :] > rows[:, None]
        else:  # overlapping selections (e.g. all O vs OW): never an atom with itself
            keep = idx2[None, :] != idx1[rows][:, None]
        d = dists[keep]
        bins = (d / dr).astype(int)
        valid = (bins >= 0) & (bins < n_bins) & (d < r_max)
        hist += np.bincount(bins[valid], minlength=n_bins)[:n_bins]
        for shift in shifts:
            # further images of every atom, its own included. For one type with itself
            # (i, j, +L) and (j, i, -L) are the same pair, so each counts a half.
            d = np.linalg.norm(wrapped + shift, axis=2).ravel()
            d = d[d < r_max]
            hist += np.bincount((d / dr).astype(int), minlength=n_bins)[:n_bins] * (
                0.5 if same_type else 1.0)

        if progress_callback:
            progress_callback(rows[-1] / total)

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
    name1: str | None = None,
    name2: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the RDF averaged over multiple trajectory frames."""
    if not frames:
        return np.array([]), np.array([])

    g_sum = None
    r_values = None
    n_frames = len(frames)

    for fi, mol in enumerate(frames):
        def _frame_progress(frac, _fi=fi):
            if progress_callback:
                progress_callback((_fi + frac) / n_frames)

        r, g = compute_rdf(mol, type1, type2, n_bins=n_bins, r_max=r_max,
                           progress_callback=_frame_progress, name1=name1, name2=name2)
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

    # Atoms that leave through a cell face come back on the other side. Follow them:
    # count the cell crossings frame to frame (minimum image) and undo them, so a
    # wrapped trajectory gives the displacement actually travelled.
    unwrapped = []
    last = None
    n_seen = 0
    crossings = np.zeros((n_atoms, 3))
    for mol in frames:
        indices = [i for i, a in enumerate(mol.atoms) if atom_filter(a)]
        n_avail = min(len(indices), n_atoms)
        pos = np.array([mol.atoms[indices[i]].x for i in range(n_avail)]).reshape(-1, 3)
        pbc = mol.pbc if mol.pbc is not None else ref.pbc
        if pbc is not None and n_avail:
            if last is None:
                last = np.zeros((n_atoms, 3))
            n_known = min(n_avail, n_seen)  # atoms that were there in an earlier frame
            crossings[:n_known] += np.round(
                (pos[:n_known] - last[:n_known]) @ pbc.reciprocal_basis_matrix)
            last[:n_avail] = pos
            n_seen = max(n_seen, n_avail)
            if crossings[:n_avail].any():
                pos = pos - crossings[:n_avail] @ pbc.basis_matrix
        unwrapped.append(pos)
    if reference_frame != 0:
        ref_pos = unwrapped[reference_frame]

    for fi in range(n_frames):
        pos = unwrapped[fi]
        n_avail = min(len(pos), len(ref_pos))
        if n_avail == 0:
            continue
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
    # beyond half the cell width one image per pair is not enough
    shifts = (image_shifts(mol.pbc, r_max, len(center_idx) * len(neighbor_idx))
              @ mol.pbc.basis_matrix) if use_pbc else ()

    for ci in range(len(center_pos)):
        diffs = neighbor_pos - center_pos[ci]
        if use_pbc:
            # Minimum image convention
            frac = diffs @ mol.pbc.reciprocal_basis_matrix
            frac -= np.round(frac)
            diffs = frac @ mol.pbc.basis_matrix
        dists = np.sqrt(np.sum(diffs ** 2, axis=1))
        if len(shifts):
            dists = np.concatenate([dists] + [np.sqrt(np.sum((diffs + shift) ** 2, axis=1))
                                              for shift in shifts])
        # Exclude self (distance ~0) and clamp to r_max
        mask = (dists > 0.1) & (dists < r_max)
        if not mask.any():
            continue
        # entry i of the running sum is the number closer than r_values[i]
        bins = np.minimum((dists[mask] / dr).astype(int) + 1, n_bins - 1)
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
                # the donor nearest to its H (a wrapped molecule can be split by a cell
                # face), and the same image of the acceptor that the distance used
                frac = dh @ mol.pbc.reciprocal_basis_matrix
                frac -= np.round(frac)
                dh = frac @ mol.pbc.basis_matrix
                ha = dh + diff_da
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



def _gyration_bonds(mol: Molecule, reuse=None) -> list:
    """The molecule's bonds: its own list, else ``reuse`` (an earlier frame's) when every
    one of those is still a bond here, else detected now."""
    if mol.bonds:
        return mol.bonds
    if reuse:
        pairs = np.asarray(reuse, dtype=int)
        if pairs.max() < len(mol.atoms):
            xyz = np.array([a.x for a in mol.atoms])
            radii = np.array([a.bond_r for a in mol.atoms], dtype=float)
            d = xyz[pairs[:, 1]] - xyz[pairs[:, 0]]
            d -= np.round(d @ mol.pbc.reciprocal_basis_matrix) @ mol.pbc.basis_matrix
            if np.all(np.linalg.norm(d, axis=1) < 0.5 * (radii[pairs[:, 0]] + radii[pairs[:, 1]])):
                return reuse
    probe = Molecule(atoms=mol.atoms, pbc=mol.pbc)
    probe.detect_bonds()
    return probe.bonds


def _gyration_whole_positions(mol: Molecule, bonds: list):
    """Positions with molecules that a cell face has split put back together, or None
    when no bond crosses a face (nothing to do).

    Each bonded group is rebuilt from its first atom along the bonds (minimum image).
    A group bonded to its own periodic image, such as a framework, cannot be made whole
    and keeps its coordinates. The molecule itself is not changed.
    """
    if not bonds:
        return None
    pbc = mol.pbc
    xyz = np.array([a.x for a in mol.atoms])
    pairs = np.asarray(bonds, dtype=int)
    if pairs.max() >= len(xyz):
        return None
    raw = xyz[pairs[:, 1]] - xyz[pairs[:, 0]]
    image = np.round(raw @ pbc.reciprocal_basis_matrix)
    if not image.any():
        return None
    step = raw - image @ pbc.basis_matrix  # bond vectors, nearest image
    adj: dict = {}
    for k, (i, j) in enumerate(pairs.tolist()):
        adj.setdefault(i, []).append((j, k, 1.0))
        adj.setdefault(j, []).append((i, k, -1.0))
    out = xyz.copy()
    seen: set = set()
    for k0 in np.nonzero(image.any(axis=1))[0]:
        start = int(pairs[k0, 0])
        if start in seen:
            continue
        group = [start]
        seen.add(start)
        placed = {start: xyz[start]}
        queue = [start]
        while queue:
            i = queue.pop()
            for j, k, sign in adj[i]:
                if j not in placed:
                    placed[j] = placed[i] + sign * step[k]
                    seen.add(j)
                    group.append(j)
                    queue.append(j)
        members = set(group)
        closed = all(np.allclose(placed[int(pairs[k, 1])] - placed[int(pairs[k, 0])], step[k],
                                 atol=1e-6)
                     for i in group for _j, k, sign in adj[i] if sign > 0 and _j in members)
        if closed:  # a finite molecule; otherwise it is bonded to its own image
            first = min(group)  # keep the lowest-numbered atom where it is
            offset = xyz[first] - placed[first]
            for i in group:
                out[i] = placed[i] + offset
    return out


def compute_gyration_radius(
    mol: Molecule,
    element: str | None = None,
    name: str | None = None,
    bonds: list | None = None,
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
    if mol.pbc is not None:
        whole = _gyration_whole_positions(mol, _gyration_bonds(mol, bonds))
        if whole is not None:
            index = {id(a): i for i, a in enumerate(mol.atoms)}
            positions = whole[[index[id(a)] for a in atoms]]
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

    bonds = None
    for fi, mol in enumerate(frames):
        if mol.pbc is not None:
            bonds = _gyration_bonds(mol, bonds)
        rg[fi] = compute_gyration_radius(mol, element, name, bonds)
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



import numpy as np



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
    # Berendsen et al., J. Phys. Chem. 91, 6269 (1987): r_OH 1.0 A, HOH 109.47 deg
    "h2o_spce": SorbateModel("h2o_spce", "H2O", "H2O SPC/E 3-site", "lj", [
        _s("OXY", 0, 0, 0, 15.9994, -0.8476, 0, 78.2, 3.166),
        _s("HYD", 0.81650, -0.57735, 0, 1.008, 0.4238, 0, 0, 0),
        _s("HYD", -0.81650, -0.57735, 0, 1.008, 0.4238, 0, 0, 0),
    ]),
    # Mahoney & Jorgensen, J. Chem. Phys. 112, 8910 (2000): r_OL 0.70 A, LOL 109.47 deg
    "h2o_tip5p": SorbateModel("h2o_tip5p", "H2O", "H2O TIP5P 5-site", "lj", [
        _s("OXY", 0, 0, 0, 15.9994, 0, 0, 80.51, 3.12),
        _s("HYD", 0.75695, -0.58588, 0, 1.008, 0.241, 0, 0, 0),
        _s("HYD", -0.75695, -0.58588, 0, 1.008, 0.241, 0, 0, 0),
        _s("LP", 0, 0.40415, 0.57155, 0, -0.241, 0, 0, 0),
        _s("LP", 0, 0.40415, -0.57155, 0, -0.241, 0, 0, 0),
    ]),
    # --- Noble gases ---
    "he": SorbateModel("he", "He", "Helium pore volume probe", "lj", [
        _s("He", 0, 0, 0, 4.002602, 0, 0, 10.22, 2.28),
    ]),
    # --- O2 ---
    # O mass 15.9994: the MPMC repo's file has 7.9997 per O, half the molecule's mass
    "o2_trappe": SorbateModel("o2_trappe", "O2", "O2 TraPPE 3-site", "lj", [
        _s("O", -0.605, 0, 0, 15.9994, -0.113, 0, 49.0, 3.02),
        _s("CoM", 0, 0, 0, 0, 0.226, 0, 0, 0),
        _s("O", 0.605, 0, 0, 15.9994, -0.113, 0, 49.0, 3.02),
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
        # LJ sites at 1.114 Å (Mullen et al., JCTC 2013, Table 2; Pham et al. 2017); the
        # MPMC repo's file has 1.208
        _s("CON", 1.114, 0, 0, 0, 0, 0, 46.47457, 2.99429),
        _s("CON", -1.114, 0, 0, 0, 0, 0, 46.47457, 2.99429),
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
        # The MPMC repo's n2_polar LJ sites; 18.128 K / 3.151 Å at 0.738 Å, used here
        # before, are the non-polarizable model's
        _s("N2N", 0.7445, 0, 0, 0, 0, 0, 16.142, 3.16141),
        _s("N2N", -0.7445, 0, 0, 0, 0, 0, 16.142, 3.16141),
    ]),
    # --- Water ---
    "h2o_pol3": SorbateModel("h2o_pol3", "H2O", "H2O POL3 3-site polarizable", "lj_polar", [
        _s("OXY", 0, 0, 0, 15.999, -0.73, 0.528, 78.50225, 3.596),
        _s("HYD", -0.816, -0.577, 0, 1.0079, 0.365, 0.17, 0, 0),
        _s("HYD", 0.816, -0.577, 0, 1.0079, 0.365, 0.17, 0, 0),
    ]),
    # Dang & Chang, J. Chem. Phys. 106, 8149 (1997): polarizable M site 0.215 A from O
    "h2o_dc97": SorbateModel("h2o_dc97", "H2O", "H2O DC97 4-site polarizable", "lj_polar", [
        _s("OXY", 0, 0, 0, 15.9994, 0, 0, 91.84, 3.234),
        _s("HYD", 0.75695, -0.58588, 0, 1.008, 0.519, 0, 0, 0),
        _s("HYD", -0.75695, -0.58588, 0, 1.008, 0.519, 0, 0, 0),
        _s("M", 0, -0.215, 0, 0, -1.038, 1.444, 0, 0),
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
    return [
        _sorbate_atom_line(model, site, start_atom_id + i, mol_id,
                           (center_x + site.x, center_y + site.y, center_z + site.z))
        for i, site in enumerate(model.sites)
    ]


def _sorbate_atom_line(model: SorbateModel, site: SorbateSite, atom_id: int,
                       mol_id: int, xyz) -> str:
    """One movable (M) MPMC atom record for a sorbate site at cartesian ``xyz``."""
    x, y, z = xyz
    return (
        f"ATOM {atom_id:>6} {site.atomtype:<4} {model.label:<3} M {mol_id:>4}    "
        f"{x:>10.6f} {y:>10.6f} {z:>10.6f}"
        f" {site.mass:>9.5f} {site.charge:>9.5f}"
        f" {site.polarizability:>9.5f} {site.epsilon:>9.5f} {site.sigma:>9.5f}"
        f" 0.0 0.0 {site.c6:>10.5f} {site.c8:>12.5f} {site.c10:>14.5f}"
    )


# ======================================================================
# Module: insert
# ======================================================================
"""Pack rigid sorbate molecules into a periodic cell (a small packmol for MPMC).

Molecules get a random position and a uniformly random orientation. A trial is
kept only if no site lands too close to a framework atom or to an already
placed molecule, measured with the minimum image in the actual (triclinic)
cell, so skewed and monoclinic cells work. "Too close" for a site pair is

    d < max(min_dist, scale * (r_i + r_j))

where r is sigma/2 for sites with LJ parameters, the element vdW radius for
framework atoms that have no sigma yet, and 0 for charge-only sites (H, M, LP),
which are then only kept min_dist away.
"""


def read_sorbate_template(path: str | Path) -> SorbateModel:
    """First movable (M) molecule of an MPMC .pqr, as a SorbateModel.

    Lets users insert any model MPMC can read, not just the built-in ones.
    Raises ValueError when the file has no M molecule.
    """
    sites: list[SorbateSite] = []
    label, mol_id = "MOL", None
    for line in Path(path).read_text().splitlines():
        tok = _expand_numeric_tokens(line.split())
        if len(tok) < 14 or tok[0] not in ("ATOM", "HETATM") or tok[4] != "M":
            continue
        if mol_id is None:
            mol_id, label = tok[5], tok[3]
        elif tok[5] != mol_id:
            break
        c = [float(v) for v in tok[16:19]] if len(tok) > 18 else [0.0, 0.0, 0.0]
        sites.append(SorbateSite(tok[2], float(tok[6]), float(tok[7]), float(tok[8]),
                                 float(tok[9]), float(tok[10]), float(tok[11]),
                                 float(tok[12]), float(tok[13]), *c))
    if not sites:
        raise ValueError(f"{path}: no movable (M) molecule found")
    pot = "phahst" if any(s.c6 for s in sites) else (
        "lj_polar" if any(s.polarizability for s in sites) else "lj")
    return SorbateModel(Path(path).stem, label, f"from {Path(path).name}", pot, sites)


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Uniformly distributed rotation matrix (Shoemake 1992)."""
    u1, u2, u3 = rng.random(3)
    a, b = sqrt(1.0 - u1), sqrt(u1)
    x, y = a * sin(2 * pi * u2), a * cos(2 * pi * u2)
    z, w = b * sin(2 * pi * u3), b * cos(2 * pi * u3)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def insert_sorbates(
    mol: Molecule,
    model: SorbateModel | str,
    n: int | None = None,
    *,
    scale: float = 0.85,
    min_dist: float = 1.0,
    max_failures: int = 5000,
    seed: int | None = None,
    progress_callback: Callable[[int], None] | None = None,
) -> list[np.ndarray]:
    """Place up to ``n`` copies of ``model`` in ``mol``'s cell without overlaps.

    ``n=None`` fills the cell: insertion stops after ``max_failures`` rejected
    trials in a row. With a number it also stops there, so fewer than ``n``
    may fit; check ``len()`` of the result. Returns one (n_sites, 3) array of
    cartesian site positions per placed molecule. ``mol`` is not modified.
    Raises CellRequiredError when ``mol`` has no cell.
    """
    require_cell(mol, "Inserting sorbates")
    if isinstance(model, str):
        model = SORBATE_MODELS[model]
    pbc = mol.pbc
    basis, recip = pbc.basis_matrix, pbc.reciprocal_basis_matrix
    rng = np.random.default_rng(seed)

    sites = np.array([[s.x, s.y, s.z] for s in model.sites], dtype=float)
    masses = np.array([s.mass for s in model.sites], dtype=float)
    center = (masses @ sites / masses.sum()) if masses.sum() > 0 else sites.mean(axis=0)
    sites -= center
    site_r = np.array([s.sigma / 2 if s.sigma > 0 and model.pot_type != "phahst" else 0.0
                       for s in model.sites])

    fw_xyz = np.array([a.x for a in mol.atoms], dtype=float).reshape(-1, 3)
    fw_r = np.array([a.sigma / 2 if a.sigma > 0 else a.vdw for a in mol.atoms], dtype=float)
    # Sample centers in the same periodic image as the framework (e.g. a cell
    # stored from -1/2 to 1/2), so molecules show up inside it in the viewer.
    # Molecules are kept whole: never wrap their atoms one by one.
    window = (np.mean(fw_xyz @ recip, axis=0) - 0.5) if len(fw_xyz) else np.zeros(3)

    placed: list[np.ndarray] = []
    others_xyz, others_r = fw_xyz, fw_r
    failures = 0
    while (n is None or len(placed) < n) and failures < max_failures:
        trial = sites @ _random_rotation(rng).T + (window + rng.random(3)) @ basis
        if len(others_xyz):
            d = trial[:, None, :] - others_xyz[None, :, :]
            d -= np.round(d @ recip) @ basis
            dist = np.sqrt((d * d).sum(axis=2))
            limit = np.maximum(min_dist, scale * (site_r[:, None] + others_r[None, :]))
            if (dist < limit).any():
                failures += 1
                continue
        placed.append(trial)
        others_xyz = np.vstack([others_xyz, trial])
        others_r = np.concatenate([others_r, site_r])
        failures = 0
        if progress_callback is not None:
            progress_callback(len(placed))
    return placed


def format_inserted_sorbates(model: SorbateModel | str, placements: list[np.ndarray],
                             start_atom_id: int, start_mol_id: int = 2) -> list[str]:
    """MPMC atom records for ``insert_sorbates`` output, one molecule ID per copy.

    Pass the result to ``write_mpmc_pdb(..., sorbate_lines=...)``.
    """
    if isinstance(model, str):
        model = SORBATE_MODELS[model]
    lines, atom_id = [], start_atom_id
    for k, xyz in enumerate(placements):
        if len(xyz) != len(model.sites):
            raise ValueError(f"placement {k} has {len(xyz)} sites, "
                             f"{model.name} has {len(model.sites)}")
        for site, pos in zip(model.sites, xyz):  # no strict=: must run on Python 3.7
            lines.append(_sorbate_atom_line(model, site, atom_id, start_mol_id + k, pos))
            atom_id += 1
    return lines


def sorbate_atoms(model: SorbateModel, placements, first_mol_id: int) -> list:
    """Atoms for packed sorbate molecules, flagged movable (M) with consecutive molecule
    ids so an MPMC export keeps them movable. ``placements`` is what insert_sorbates
    returns: one (n_sites, 3) array per molecule."""
    atoms = []
    for k, xyz in enumerate(placements):
        for site, pos in zip(model.sites, xyz):
            atom = Atom(pos[0], pos[1], pos[2], site.atomtype)
            atom.charge, atom.alpha = site.charge, site.polarizability
            atom.epsilon, atom.sigma = site.epsilon, site.sigma
            atom.c6, atom.c8, atom.c10 = site.c6, site.c8, site.c10
            atom.mpmc_flag, atom.mpmc_mol = "M", first_mol_id + k
            atom.mpmc_res, atom.file_mass = model.label, site.mass
            atoms.append(atom)
    return atoms


def parse_sorbate_count(text: str):
    """'fill' -> None (pack until full), a positive whole number -> that number.

    Anything else raises ValueError (unicode digits such as '²' included).
    """
    raw = text.strip().lower()
    if raw == "fill":
        return None
    if raw.isascii() and raw.isdigit() and int(raw) > 0:
        return int(raw)
    raise ValueError("Enter a positive number or 'fill'")

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

    Source installs read the on-disk file under TEMPLATES_DIR so editable installs can
    customise the shipped templates. The amalgamated single-file dist uses only its
    _BUNDLED_TEMPLATES: there TEMPLATES_DIR is whatever folder the file was copied to, and
    an unrelated file of the same name must not be taken for a template.
    """
    p = TEMPLATES_DIR / name
    if _BUNDLED_TEMPLATES:
        if name in _BUNDLED_TEMPLATES:
            return _BUNDLED_TEMPLATES[name]
        raise FileNotFoundError(f"template not found: {name}")
    if p.is_file():
        return p.read_text()
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
                    choices=tuple(SORBATE_MODELS)),  # every built-in model
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
    'cp2k_cell_opt.inp': '! CP2K cell optimization — generated by pdb_wizard from cp2k_cell_opt.inp\n! Placeholders are triple-brace tokens. The per-element &KIND block\n! (basis + GTH pseudopotential) is built by the generator and slotted in\n! as a single multi-line value.\n\n&GLOBAL\n  PROJECT {{{project}}}\n  RUN_TYPE CELL_OPT\n  PRINT_LEVEL LOW\n&END GLOBAL\n\n&FORCE_EVAL\n  METHOD QUICKSTEP\n  STRESS_TENSOR ANALYTICAL\n  &SUBSYS\n    &CELL\n      A   {{{cell_vec_a}}}\n      B   {{{cell_vec_b}}}\n      C   {{{cell_vec_c}}}\n      PERIODIC XYZ\n    &END CELL\n    &TOPOLOGY\n      COORD_FILE_NAME {{{coords_file}}}\n      COORD_FILE_FORMAT XYZ\n      &CENTER_COORDINATES\n      &END\n    &END TOPOLOGY\n{{{kind_blocks}}}  &END SUBSYS\n  &DFT\n    BASIS_SET_FILE_NAME BASIS_MOLOPT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    &MGRID\n      CUTOFF {{{cutoff_ry}}}\n      REL_CUTOFF {{{rel_cutoff_ry}}}\n    &END MGRID\n    &SCF\n      EPS_SCF {{{eps_scf}}}\n      MAX_SCF {{{max_scf}}}\n      &OT\n        MINIMIZER DIIS\n        PRECONDITIONER FULL_SINGLE_INVERSE\n      &END OT\n      &OUTER_SCF\n        EPS_SCF {{{eps_scf}}}\n        MAX_SCF 20\n      &END\n    &END SCF\n    &XC\n      &XC_FUNCTIONAL {{{xc_functional}}}\n      &END\n    &END XC\n  &END DFT\n&END FORCE_EVAL\n\n&MOTION\n  &CELL_OPT\n    OPTIMIZER LBFGS\n    KEEP_ANGLES T\n    MAX_ITER {{{n_opt_steps}}}\n    TYPE DIRECT_CELL_OPT\n  &END CELL_OPT\n  &PRINT\n    &TRAJECTORY\n      &EACH\n        CELL_OPT 1\n      &END\n    &END\n    &CELL\n      &EACH\n        CELL_OPT 1\n      &END\n    &END\n  &END\n&END MOTION\n',
    'cp2k_dft_md.inp': '! CP2K Born-Oppenheimer NVT MD — generated by pdb_wizard from cp2k_dft_md.inp\n! Placeholders are triple-brace tokens. The per-element &KIND block\n! (basis + GTH pseudopotential) is built by the generator and slotted in\n! as a single multi-line value.\n\n&GLOBAL\n  PROJECT {{{project}}}\n  RUN_TYPE MD\n  PRINT_LEVEL LOW\n&END GLOBAL\n\n&FORCE_EVAL\n  METHOD QUICKSTEP\n  STRESS_TENSOR ANALYTICAL\n  &SUBSYS\n    &CELL\n      A   {{{cell_vec_a}}}\n      B   {{{cell_vec_b}}}\n      C   {{{cell_vec_c}}}\n      PERIODIC XYZ\n    &END CELL\n    &TOPOLOGY\n      COORD_FILE_NAME {{{coords_file}}}\n      COORD_FILE_FORMAT XYZ\n      &CENTER_COORDINATES\n      &END\n    &END TOPOLOGY\n{{{kind_blocks}}}  &END SUBSYS\n  &DFT\n    BASIS_SET_FILE_NAME BASIS_MOLOPT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n    &MGRID\n      CUTOFF {{{cutoff_ry}}}\n      REL_CUTOFF {{{rel_cutoff_ry}}}\n    &END MGRID\n    &SCF\n      EPS_SCF {{{eps_scf}}}\n      MAX_SCF {{{max_scf}}}\n      &OT\n        MINIMIZER DIIS\n        PRECONDITIONER FULL_SINGLE_INVERSE\n      &END OT\n      &OUTER_SCF\n        EPS_SCF {{{eps_scf}}}\n        MAX_SCF 20\n      &END\n    &END SCF\n    &XC\n      &XC_FUNCTIONAL {{{xc_functional}}}\n      &END\n    &END XC\n  &END DFT\n&END FORCE_EVAL\n\n&MOTION\n  &MD\n    ENSEMBLE NVT\n    TIMESTEP {{{timestep_fs}}}\n    STEPS {{{n_md_steps}}}\n    TEMPERATURE {{{temperature}}}\n    &THERMOSTAT\n      TYPE NOSE\n      &NOSE\n        TIMECON 50\n      &END\n    &END\n  &END MD\n  &PRINT\n    &TRAJECTORY\n      FILENAME =trajectory.xyz\n      FORMAT XYZ\n      &EACH\n        MD {{{md_traj_interval}}}\n      &END\n    &END\n    &VELOCITIES OFF\n    &END\n    &RESTART\n      FILENAME =restart\n      &EACH\n        MD 100\n      &END\n    &END\n  &END\n&END MOTION\n',
    'mpmc_uvt.inp': "! MPMC isotherm input — generated by pdb_wizard from template mpmc_uvt.inp\n! Placeholders are written as triple-brace tokens; the Input Generator\n! UI re-renders this file with the values you choose.\n\njob_name              {{{job_name}}}\n\nensemble              {{{ensemble}}}\n\nabcbasis              {{{cell_a}}} {{{cell_b}}} {{{cell_c}}} {{{cell_alpha}}} {{{cell_beta}}} {{{cell_gamma}}}\n\ntemperature           {{{temperature}}}\npressure              {{{pressure}}}\n\nnumsteps              {{{numsteps}}}\ncorrtime              {{{corrtime}}}\ninsert_probability    {{{insert_probability}}}\nmove_factor           {{{move_factor}}}\nrot_factor            {{{rot_factor}}}\n\npqr_input             {{{pqr_input}}}\npqr_output            output.pqr\npqr_restart           restart.pqr\ntraj_output           trajectory.pqr\nenergy_output         energy.dat\n\n! No insert_input: the sorbate is already embedded in pqr_input as a movable\n! molecule, which MPMC clones for insertions. Pointing insert_input at a\n! separate pool triggers a double-free in MPMC's cleanup (crash on exit).\n\nrd_lrc                on\nwrapall               on\n\n{{{model_keywords}}}\n\npop_histogram         on\npop_histogram_output  histogram.dat\n",
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


def _cp2k_atoms(mol: Molecule) -> list:
    """(symbol, atom) for the real atoms. Sites of MPMC sorbate models get their true
    element ('COG' is carbon, not cobalt) and massless sites are left out: CP2K needs
    atoms, not model sites."""
    out = []
    for atom in mol.atoms:
        element = site_element(atom)
        if element is not None:
            out.append((element.symbol, atom))
    return out


def _formula_summary(mol: Molecule) -> str:
    counts = Counter(symbol for symbol, _ in _cp2k_atoms(mol))
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
    elements_present = sorted({symbol for symbol, _ in _cp2k_atoms(mol)})
    return "".join(_kind_section(el) for el in elements_present)


def _write_xyz(mol: Molecule, path: Path) -> None:
    atoms = _cp2k_atoms(mol)
    lines = [str(len(atoms)), f"{_formula_summary(mol)} from pdb_wizard"]
    for symbol, a in atoms:
        lines.append(
            f"{symbol}  {a.x[0]:.6f}  {a.x[1]:.6f}  {a.x[2]:.6f}"
        )
    path.write_text("\n".join(lines) + "\n")


def _cell_vector_vars(pbc) -> dict:
    """{'cell_vec_a': 'ax ay az', ...}: the rows of the cell's basis matrix, in Å."""
    def fmt(row):
        return " ".join(f"{round(float(v), 6) + 0.0:.6f}" for v in row)

    a, b, c = pbc.basis_matrix
    return {"cell_vec_a": fmt(a), "cell_vec_b": fmt(b), "cell_vec_c": fmt(c)}


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

    # The cell vectors themselves, in the frame of the XYZ coordinates. CP2K's
    # ALPHA_BETA_GAMMA only applies together with ABC, so writing orthogonal A/B/C plus
    # angles (as before) gave a rectangular box for hexagonal or monoclinic cells.
    common_vars = {
        "project": project_name,
        **_cell_vector_vars(mol.pbc),
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
        "elements": sorted({symbol for symbol, _ in _cp2k_atoms(mol)}),
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


# LJ cutoff (Å) the box is checked against. Minimum image needs half the
# smallest perpendicular cell width to be at least this.
_DEFAULT_LJ_CUTOFF = 10.0
# A charge is "set" if |q| exceeds this — anything smaller is rounding noise.
_CHARGE_NONZERO = 1e-9
# An LJ parameter is "set" if epsilon > this — typical organic atom epsilon
# is 30-150 K, hydrogen is ~10 K, so 0.1 K is a safe lower bound for "applied".
_EPSILON_SET = 0.1
# Sorbate–framework overlap distance (Å) below which we warn.
_SORBATE_OVERLAP = 1.5


# Polarization and PHAHST settings from the MPMC tutorials, which call them "generally
# acceptable for most situations" (tutorials/02_bssp_h2_mof-5, 03_phahst_h2_hkust-1).
_POLARIZATION_KEYWORDS = """polarization          on
polar_damp_type       exponential
polar_damp            2.1304
polar_wolf            on
polar_wolf_alpha      0.13
polar_iterative       on
polar_palmo           on
polar_gs_ranked       on
polar_gamma           1.03
polar_max_iter        4"""
_PHAHST_KEYWORDS = """disp_expansion        on
extrapolate_disp_coeffs on
damp_dispersion       on
cavity_autoreject_repulsion 10000"""


def _model_flags(mol: Molecule, sorbate_model=None) -> dict:
    """What the framework and sorbate parameters need MPMC to switch on."""
    model = SORBATE_MODELS.get(sorbate_model) if isinstance(sorbate_model, str) else sorbate_model
    sites = model.sites if model is not None else []
    framework = [a for a in mol.atoms if a.mpmc_flag != "M"]
    return {
        "polar": any(a.alpha > 0 for a in mol.atoms) or any(s.polarizability > 0 for s in sites),
        "phahst_framework": any(max(a.c6, a.c8, a.c10) > 0 for a in framework),
        "phahst_sorbate": any(max(s.c6, s.c8, s.c10) > 0 for s in sites),
        "lj_sorbate": bool(sites) and not any(max(s.c6, s.c8, s.c10) > 0 for s in sites),
        "sorbate_charged": any(abs(s.charge) > _CHARGE_NONZERO for s in sites),
        "sorbate_known": model is not None,
    }


def mpmc_model_keywords(mol: Molecule, sorbate_model=None) -> str:
    """The mpmc.inp lines the PQR parameters need.

    Without them MPMC ignores the polarizabilities and reads PHAHST beta/rho (columns
    13-14) as Lennard-Jones epsilon/sigma. The cutoff is half the smallest perpendicular
    cell width: the same as MPMC's default (half the shortest lattice vector) for an
    orthogonal cell, and inside the minimum-image limit for a skewed one.
    """
    flags = _model_flags(mol, sorbate_model)
    cutoff = float(mol.pbc.perpendicular_widths().min()) / 2.0
    lines = [f"pbc_cutoff            {np.floor(cutoff * 1e4) / 1e4:.4f}"]
    if flags["phahst_framework"] or flags["phahst_sorbate"]:
        lines += ["", _PHAHST_KEYWORDS]
    if flags["polar"]:
        lines += ["", _POLARIZATION_KEYWORDS]
    return "\n".join(lines)


def validate_for_mpmc(
    mol: Molecule,
    sorbate_pos: tuple[float, float, float] | None = None,
    lj_cutoff: float = _DEFAULT_LJ_CUTOFF,
    sorbate_model=None,
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
      - Half the smallest perpendicular width < lj_cutoff (cell too thin)
      - MPMC's automatic cutoff exceeds the minimum-image limit (skewed cell)
      - Sorbate placement overlaps a framework atom (within 1.5 Å)

    Triclinic cells need no warning: 'abcbasis a b c α β γ' builds the same
    basis (a along x, b in the xy plane) as PBC._compute_matrices.
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

    # Charges (an uncharged framework is fine with an uncharged sorbate, e.g. TraPPE CH4)
    flags = _model_flags(mol, sorbate_model)
    charges = np.array([a.charge for a in mol.atoms])
    n_set = int((np.abs(charges) > _CHARGE_NONZERO).sum())
    if n_set == 0 and flags["sorbate_known"] and not flags["sorbate_charged"]:
        pass
    elif n_set == 0:
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

    if flags["phahst_framework"] != flags["phahst_sorbate"] and flags["sorbate_known"]:
        which = "framework" if flags["phahst_framework"] else "sorbate"
        warnings.append(
            f"Only the {which} has PHAHST (C6/C8/C10) parameters. MPMC's disp_expansion "
            f"applies to every pair, so the other side's LJ epsilon/sigma would be read as "
            f"exponential-repulsion beta/rho. Use PHAHST (or LJ) for both."
        )

    # Box size vs LJ cutoff. Minimum image is only exact up to half the
    # smallest perpendicular width, which for a skewed cell is less than
    # min(a, b, c) / 2.
    pbc = mol.pbc
    max_cutoff = float(pbc.perpendicular_widths().min()) / 2.0
    if max_cutoff < lj_cutoff:
        warnings.append(
            f"Cell is too thin for a {lj_cutoff:.1f} Å cutoff: half the "
            f"smallest perpendicular width is {max_cutoff:.2f} Å. Either build "
            f"a supercell (Edit > Extend Axis) or set pbc_cutoff ≤ "
            f"{max_cutoff:.2f} in mpmc.inp."
        )

    # Without pbc_cutoff, MPMC uses half the shortest lattice vector, which
    # exceeds the minimum-image limit in any non-orthogonal cell. Only pairs
    # in the shell between the two radii are affected, so stay quiet for
    # mild skew (<5%); a hexagonal cell (γ = 120°) is ~15% over.
    mpmc_auto = pbc.shortest_lattice_vector() / 2.0
    if mpmc_auto > 1.05 * max_cutoff:
        warnings.append(
            f"Cell is skewed: MPMC's automatic cutoff ({mpmc_auto:.2f} Å, half "
            f"the shortest lattice vector) is larger than the minimum-image "
            f"limit ({max_cutoff:.2f} Å), so some pairs would be missed. Set "
            f"pbc_cutoff ≤ {max_cutoff:.2f} in mpmc.inp."
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
        sorbate_model, sx, sy, sz, mol_id=next_mpmc_mol_id(mol),
        start_atom_id=len(mol.atoms) + 1,
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
            # one decimal as always, unless that would not be the temperature asked for
            "temperature": (f"{temperature:.1f}" if float(f"{temperature:.1f}") == temperature
                            else repr(float(temperature))),
            "pressure": pressure,
            "numsteps": numsteps,
            "corrtime": corrtime,
            "insert_probability": insert_probability,
            "move_factor": move_factor,
            "rot_factor": rot_factor,
            "pqr_input": "input.pqr",
            "model_keywords": mpmc_model_keywords(mol, sorbate_model),
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
        if hist_file.exists() and not _is_opendx(hist_file):
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


def _is_opendx(path: Path) -> bool:
    """MPMC's histogram.dat is an OpenDX density grid, not an N histogram."""
    try:
        with open(path, errors="replace") as f:
            head = f.read(400)
    except OSError:
        return False
    return "OpenDX" in head or "gridpositions" in head


_LOG_NAMES = ("mpmc.log", "runlog", "run.log", "output.log")
_LOG_VALUE = re.compile(r"OUTPUT:\s*(N|wt %|wt % \(ME\)|qst)\s*=\s*([-+0-9.eE]+)"
                        r"(?:\s*\+-\s*([-+0-9.eE]+))?")
_LOG_STEP = re.compile(r"Completed step\s+(\d+)/(\d+)")


def read_mpmc_log(p_dir: Path) -> dict:
    """MPMC's own final averages from the run log in a pressure folder.

    Looks for mpmc.log (what run_all.sh writes), a few common names, then any *.log or
    slurm-*.out. Returns {} when there is no log, else a dict with ``file``, ``finished``,
    ``step``/``steps`` and, when printed, ``N``, ``wt``, ``wt_me`` (value, error) and
    ``qst``.
    """
    candidates = [p_dir / name for name in _LOG_NAMES]
    candidates += sorted(p_dir.glob("*.log")) + sorted(p_dir.glob("slurm-*.out"))
    for path in candidates:
        if not path.is_file():
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        if "OUTPUT:" not in text:
            continue
        out: dict = {"file": path.name, "finished": "simulation exiting successfully" in text}
        keys = {"N": "N", "wt %": "wt", "wt % (ME)": "wt_me", "qst": "qst"}
        for name, value, err in _LOG_VALUE.findall(text):  # the last block wins
            try:
                out[keys[name]] = (float(value), float(err) if err else None)
            except ValueError:
                continue
        steps = _LOG_STEP.findall(text)
        if steps:
            out["step"], out["steps"] = int(steps[-1][0]), int(steps[-1][1])
        return out
    return {}


def isotherm_log_table(base_dir: str) -> str:
    """MPMC's own uptake, error bars, wt% and Qst for each pressure folder, as text, with a
    line for every run that has not finished. '' when no folder has a log."""
    base = Path(base_dir)
    if not base.is_dir():
        return ""
    rows = []
    for p_dir in sorted(base.iterdir()):
        if not p_dir.is_dir() or not p_dir.name.startswith("P_"):
            continue
        p_str = p_dir.name.replace("P_", "").replace("_atm", "").replace("p", ".").replace("m", "-")
        try:
            pressure = float(p_str)
        except ValueError:
            continue
        log = read_mpmc_log(p_dir)
        if log:
            rows.append((pressure, log))
    if not rows:
        return ""
    rows.sort(key=lambda r: r[0])

    def pm(log: dict, key: str, fmt: str) -> str:
        if key not in log:
            return "-"
        value, err = log[key]
        return f"{value:{fmt}}" + (f" ± {err:{fmt}}" if err is not None else "")

    lines = ["MPMC's own averages (whole run, from the log files):", "",
             f"{'P (atm)':>10}  {'N (molecules)':>20}  {'wt %':>18}  {'Qst (kJ/mol)':>12}"]
    unfinished = []
    for pressure, log in rows:
        lines.append(f"{pressure:>10g}  {pm(log, 'N', '.3f'):>20}  {pm(log, 'wt', '.4f'):>18}  "
                     f"{pm(log, 'qst', '.3f'):>12}")
        if not log["finished"]:
            done = f" (step {log['step']}/{log['steps']})" if "step" in log else ""
            unfinished.append(f"⚠ {pressure:g} atm has not finished{done}; its numbers will "
                              f"still change.")
    lines += [""] + unfinished
    lines.append("The plot shows the average of N over the last quarter of energy.dat; MPMC's "
                 "average covers the whole run, including equilibration, so the two differ "
                 "until the run is long.")
    return "\n".join(lines)

# ======================================================================
# Module: reduce_cell
# ======================================================================
"""Reduce a supercell back to a single unit cell, with sorbate selection."""



import numpy as np


# Two movable molecules whose sites all coincide within this (A) once folded into the
# small cell are periodic copies of one molecule.
_REDUCE_COPY_TOL = 0.05


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

    # Group by element first to reduce O(n^2) comparisons. Movable (M) molecules are
    # sorbates by definition and are handled below: their sites can sit within the
    # tolerance of each other (PHAHST H2, the CON/COE sites of CO2) and would otherwise
    # be mistaken for periodic copies of one framework atom.
    elem_indices: dict[str, list[int]] = {}
    for i, a in enumerate(mol.atoms):
        if a.mpmc_flag != "M":
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

    n_frozen = sum(len(indices) for indices in elem_indices.values())
    n_periodic = n_mult * len(framework_indices)
    if n_mult > 1 and n_frozen and 2 * n_periodic < n_frozen:
        raise ValueError(
            f"this is not a {na} x {nb} x {nc} supercell: only {n_periodic} of its "
            f"{n_frozen} framework atoms repeat {n_mult} times under the sub-cell "
            f"translations. Check the numbers of cells along a, b and c.")

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

    # Movable molecules, one entry per molecule id. Periodic copies of the same
    # molecule (a supercell made with Extend Axis after the sorbates were added) land
    # on top of each other in the small cell: keep the first of each.
    movable: dict = {}
    for i, a in enumerate(mol.atoms):
        if a.mpmc_flag == "M":
            movable.setdefault(a.mpmc_mol if a.mpmc_mol is not None else ("atom", i), []).append(i)
    kept: list[list[int]] = []
    for indices in movable.values():
        duplicate = False
        for other in kept:
            if len(other) != len(indices):
                continue
            diff = (prim_cart[indices] - prim_cart[other]) @ prim_inv
            diff -= np.round(diff)
            if (np.linalg.norm(diff @ prim_basis, axis=1) < _REDUCE_COPY_TOL).all() and all(
                    mol.atoms[i].name == mol.atoms[j].name for i, j in zip(indices, other)):
                duplicate = True
                break
        if not duplicate:
            kept.append(indices)
    sorbate_molecules.extend(kept)

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
    prim_recip = prim_pbc.reciprocal_basis_matrix
    # A movable molecule moves as a whole (by the lattice vector that brings its first
    # site into the cell), so it stays the rigid model MPMC expects.
    first_site: dict = {}
    for sm in sorbate_mols:
        if all(mol.atoms[i].mpmc_flag == "M" for i in sm):
            for i in sm:
                first_site[i] = sm[0]
    new_atoms = []
    for idx in sorted(keep_atoms):
        a = mol.atoms[idx]
        if idx in first_site:
            origin = mol.atoms[first_site[idx]].x
            rel = (a.x - origin) @ prim_recip
            rel -= np.round(rel)
            frac = origin @ prim_recip
            frac = frac - np.floor(frac) + rel
        else:
            # Convert to fractional in primitive cell
            frac = a.x @ prim_recip
            frac -= np.floor(frac)
        new_pos = frac @ prim_basis
        # a full copy: charge and LJ parameters, and also C6/C8/C10, the movable flag,
        # molecule id, residue name and model mass
        new_atom = copy.deepcopy(a)
        new_atom.x = np.array([float(new_pos[0]), float(new_pos[1]), float(new_pos[2])])
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
            if tokens and not header_names:  # later comment lines are not the header
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

    # MPMC flushes the file in blocks, so a running job's last line can be cut short.
    data_lines = [row for row in data_lines if len(row) == len(data_lines[0])]

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
    f = (a1 * np.exp(-b1 * s2) + a2 * np.exp(-b2 * s2) +
         a3 * np.exp(-b3 * s2) + a4 * np.exp(-b4 * s2) + c)
    if element not in SCATTERING_FACTORS:
        # No tabulated coefficients: the generic curve, scaled so that f(0) is the
        # element's number of electrons (it was 12.5 for every such element).
        z = ELEMENTS[element].atomic_number if element in ELEMENTS else 0
        if z > 0:
            f = f * (z / (a1 + a2 + a3 + a4 + c))
    return f


def pxrd_untabulated_elements(mol: Molecule) -> list:
    """Elements in ``mol`` that have no tabulated form factor (a scaled generic curve is
    used for them, so their peak intensities are approximate)."""
    symbols = {el.symbol for el in (site_element(a) for a in mol.atoms) if el is not None}
    return sorted(sym for sym in symbols
                  if sym not in SCATTERING_FACTORS and sym in ELEMENTS
                  and ELEMENTS[sym].atomic_number > 0)


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
    two_theta_arr, pattern = pxrd_raw_pattern(mol, wavelength, two_theta_min, two_theta_max,
                                              n_points, peak_width, progress_callback)

    # Normalize to max = 100
    if len(pattern) and pattern.max() > 0:
        pattern = pattern / pattern.max() * 100

    return two_theta_arr, pattern


def pxrd_raw_pattern(
    mol: Molecule,
    wavelength: float = 1.5406,
    two_theta_min: float = 5.0,
    two_theta_max: float = 50.0,
    n_points: int = 2000,
    peak_width: float = 0.1,
    progress_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    """compute_pxrd's pattern before it is scaled to a maximum of 100: the summed
    Gaussian-broadened intensities on ``linspace(two_theta_min, two_theta_max, n_points)``.
    Empty arrays without a cell, atoms or reflections in range. Patterns of trajectory
    frames are averaged in this form (intensities, not amplitudes) and scaled afterwards.
    """
    if mol.pbc is None or len(mol.atoms) == 0:
        return np.array([]), np.array([])

    two_theta_all, intensities = pxrd_reflections(mol, wavelength, two_theta_min,
                                                  two_theta_max)
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

    return two_theta_arr, pattern


def pxrd_reflections(mol: Molecule, wavelength: float = 1.5406, two_theta_min: float = 5.0,
                     two_theta_max: float = 50.0, with_hkl: bool = False):
    """(2-theta, intensity) of every reflection in range; with ``with_hkl`` also the indices.

    Intensities are |F|^2 times the Lorentz-polarization factor.
    """
    pbc = mol.pbc
    # pbc.reciprocal_basis_matrix is inv(basis): its COLUMNS are the reciprocal vectors
    # a*, b*, c* (a_i . a*_j = delta_ij), so G = h a* + k b* + l c* = hkl @ inv(basis).T.
    # Using its rows put peaks in the wrong place for any non-orthogonal cell.
    inv_basis = pbc.reciprocal_basis_matrix
    recip_vectors = inv_basis.T

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

    # real atoms only, with their true elements: sites of MPMC sorbate models are named
    # after the model ('COG' is carbon) and some are massless
    sites = [(a, el) for a, el in ((a, site_element(a)) for a in mol.atoms) if el is not None]
    if not sites:
        empty = np.zeros(0)
        return (empty, empty, np.zeros((0, 3), dtype=int)) if with_hkl else (empty, empty)
    coords = np.array([a.x for a, _ in sites])
    elements = np.array([el.symbol for _, el in sites])
    unique_elements = list(set(elements.tolist()))

    # Pre-compute fractional coordinates and element-grouped slices
    frac_coords = coords @ inv_basis
    el_groups = {el: frac_coords[elements == el] for el in unique_elements}

    # Generate the full (h,k,l) grid, drop (0,0,0), filter to two_theta range.
    h_arr = np.arange(-h_max, h_max + 1)
    k_arr = np.arange(-k_max, k_max + 1)
    l_arr = np.arange(-l_max, l_max + 1)
    H, K, L = np.meshgrid(h_arr, k_arr, l_arr, indexing="ij")
    hkl = np.stack([H.ravel(), K.ravel(), L.ravel()], axis=1)        # (M, 3)
    hkl = hkl[(hkl != 0).any(axis=1)]
    G_all = hkl.astype(float) @ recip_vectors                        # (M, 3)
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
            # phase shape: (n_refl, n_el), built in slices to bound the memory
            cos_sum = np.empty(n_refl)
            sin_sum = np.empty(n_refl)
            step = max(1, 4_000_000 // len(frac_g))
            for r0 in range(0, n_refl, step):
                phase = 2 * np.pi * (hkl[r0:r0 + step].astype(float) @ frac_g.T)
                cos_sum[r0:r0 + step] = np.cos(phase).sum(axis=1)
                sin_sum[r0:r0 + step] = np.sin(phase).sum(axis=1)
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

    if with_hkl:
        return two_theta_all, intensities, hkl
    return two_theta_all, intensities

# ======================================================================
# Module: analyses
# ======================================================================
"""Analysis runs shared by the terminal viewer and the desktop GUI.

Each ``run_*`` function takes plain parameters (elements, atom names, frame range, ...),
calls the numerical routine and returns an :class:`AnalysisResult` with the curve, its
labels, the one-line summary the viewers show and the CSV header. Keeping this out of the
UI code means both front ends compute and report exactly the same thing.
"""



import numpy as np


ALL_ATOMS = "(all)"
PXRD_SOURCES = (("Cu Ka 1.5406", 1.5406), ("Mo Ka 0.7107", 0.7107), ("Co Ka 1.7903", 1.7903))


@dataclass
class AnalysisResult:
    """One analysis outcome. ``x``/``y`` are empty when there is nothing to plot."""

    x: np.ndarray = field(default_factory=lambda: np.array([]))
    y: np.ndarray = field(default_factory=lambda: np.array([]))
    title: str = ""
    x_label: str = ""
    summary: str = ""
    csv_header: str = ""
    text: str = ""
    y_min: Optional[float] = 0.0
    y_max: Optional[float] = None
    extra: Dict = field(default_factory=dict)
    csv_data: Optional[tuple] = None  # columns for CSV when they differ from (x, y)

    @property
    def empty(self) -> bool:
        return len(self.x) == 0

    def csv_arrays(self) -> tuple:
        return self.csv_data if self.csv_data is not None else (self.x, self.y)


def _opt(value: Optional[str]) -> Optional[str]:
    """'(all)' and '' mean no filter."""
    return None if value in (None, "", ALL_ATOMS) else value


def _grow(values: np.ndarray) -> float:
    return float(values.max()) * 1.1 if len(values) and values.max() > 0 else 1.0


def frame_slice(frames: Optional[list], start: int = 1,
                end: Optional[int] = None) -> Optional[list]:
    """Frames ``start``..``end`` (1-based, inclusive), as the viewers' Start/End boxes mean."""
    if not frames or len(frames) < 2:
        return frames
    return frames[max(start, 1) - 1:len(frames) if end is None else end]


def elements_of(mol: Molecule) -> List[str]:
    return sorted({a.element.symbol for a in mol.atoms})


def names_of(mol: Molecule, element: Optional[str] = None) -> List[str]:
    """Atom names, optionally only those of one element (for the Name filters)."""
    el = _opt(element)
    return sorted({a.name.strip() for a in mol.atoms if el is None or a.element.symbol == el})


def unavailable(key: str, mol: Optional[Molecule], frames: Optional[list]) -> Optional[str]:
    """Why an analysis cannot run on this structure (None if it can).

    The terminal viewer's Analysis-menu rules and messages, plus one more: void volume
    needs a cell (the terminal viewer reports zeros without one).
    """
    n_frames = len(frames) if frames else 0
    # A trajectory without cell data carries a huge stand-in box; it is not a cell.
    no_cell = mol is None or mol.pbc is None or is_fallback_cell(mol.pbc)
    if key in ("energy_plot", "isotherm_plot", "db_search"):
        return None
    if key in ("msd", "rmsd") and n_frames < 3:
        return f"Need a trajectory with 3+ frames for {key.upper()}"
    if key == "gyration" and n_frames < 2:
        return "Need a trajectory for Rg"
    if key == "density3d":
        if n_frames < 2:
            return "Need a trajectory with PBC for 3D density"
        if no_cell:
            return "Need PBC for 3D density map"
    if key == "pxrd" and (mol is None or not mol.atoms or no_cell):
        return "Need a crystal structure with PBC for PXRD"
    if mol is None:
        return "No atoms loaded"
    if key in ("info", "input_generator", "check", "symmetry"):
        return None
    if not mol.atoms:
        return "No atoms loaded yet" if key == "rdf" else "No atoms loaded"
    if no_cell:
        need = {"surface_area": "Need PBC for surface area", "pore_size": "Need PBC for pore size",
                "void_vol": "Need PBC for void volume"}
        if key in need:
            return need[key]
    return None


# ----------------------------------------------------------------------
# Single-structure reports
# ----------------------------------------------------------------------

def system_info_text(mol: Molecule, filepath: str = "") -> str:
    info = system_info(mol, filepath)
    lines = [f"Filename:  {info.get('filename', '')}", f"Atoms:     {info.get('n_atoms', 0)}"]
    if "volume" in info:
        c = info.get("cell", {})
        lines.append(f"Cell:      {c.get('a', 0):.4f}  {c.get('b', 0):.4f}  {c.get('c', 0):.4f}")
        lines.append(
            f"Angles:    {c.get('alpha', 0):.2f}  {c.get('beta', 0):.2f}  {c.get('gamma', 0):.2f}")
        lines.append(f"Volume:    {info.get('volume', 0):.2f} A^3")
        lines.append(f"Density:   {info.get('density', 0):.4f} g/cm^3")
        lines.append("")
        for row in info.get("basis_matrix", []):
            lines.append(f"  {row[0]:>14.8f} {row[1]:>14.8f} {row[2]:>14.8f}")
    lines.append("")
    lines.append("Total:     " + "  ".join(f"{k}{v}" for k, v in formula_unit(mol).items()))
    lines.append("Formula:   " + "  ".join(f"{k}{v}" for k, v in formula_unit_reduced(mol).items()))
    return "\n".join(lines)


def void_volume_report(geo: Dict, he: Dict, n_samples: int = 100000,
                       probe: float = 1.14) -> str:
    """The Void Volume report from the geometric and He-probe results."""
    lines = [f"Void Volume Analysis ({n_samples:,} MC samples)", ""]
    for label, r in (("Geometric (no probe):", geo), (f"He probe (r = {probe:g} A):", he)):
        lines += [label,
                  f"  Void volume:   {r.get('void_volume', 0):.2f} A^3",
                  f"  Total volume:  {r.get('total_volume', 0):.2f} A^3",
                  f"  Void fraction: {r.get('void_fraction', 0) * 100:.1f}%", ""]
    return "\n".join(lines[:-1])


def void_volume_text(mol: Molecule, n_samples: int = 100000, probe: float = 1.14) -> str:
    """Geometric and He-probe void volume, computed and formatted."""
    return void_volume_report(void_volume(mol, n_samples, 0.0), void_volume(mol, n_samples, probe),
                              n_samples, probe)


def surface_area_text(mol: Molecule, n_samples: int = 100000, probe: float = 1.4) -> str:
    r = surface_area(mol, n_samples, probe)
    return "\n".join([
        "Solvent-Accessible Surface Area",
        f"  Probe radius: {probe:g} A" + (" (N2-sized)" if probe == 1.4 else ""),
        f"  Surface area: {r.get('surface_area', 0):.1f} A^2",
        f"  Specific area: {r.get('area_per_volume', 0):.1f} m^2/g",
    ])


def check_text(mol: Optional[Molecule], frames: Optional[list] = None) -> str:
    """The Check Structure report (see ``checks``)."""
    return structure_check_text(mol, frames)


def symmetry_text(mol: Optional[Molecule]) -> str:
    """The Find Symmetry report: space groups ranked by the atomic drift they need."""
    if mol is None or not mol.atoms:
        return "No atoms loaded"
    return symmetry_report(mol)


def charge_report(mol: Molecule) -> str:
    """Summary after QEq: totals, extremes and per-element averages."""
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
    by_el: Dict[str, List[float]] = {}
    for a in mol.atoms:
        by_el.setdefault(a.element.symbol, []).append(a.charge)
    for el in sorted(by_el):
        arr = by_el[el]
        lines.append(f"  {el:<3} ({len(arr):>4}):  avg = {np.mean(arr):>8.5f}  "
                     f"range = [{min(arr):.5f}, {max(arr):.5f}]")
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Curves
# ----------------------------------------------------------------------

def run_pore_size(mol: Molecule, n_samples: int = 50000, n_bins: int = 50,
                  max_pore: float = 20.0) -> AnalysisResult:
    if mol.pbc is None:
        return AnalysisResult(summary="No PBC for pore size calculation")
    bins, hist = pore_size_distribution(mol, n_samples, n_bins, max_pore)
    peak = bins[hist.argmax()] if hist.max() > 0 else 0
    return AnalysisResult(bins * 2, hist, "Pore Size Distribution", "pore diameter (A)",
                          f"Dominant pore diameter: {peak * 2:.2f} A",
                          "pore_diameter_A,fraction", y_max=_grow(hist))


def _empty_window(frames: Optional[list], start: int,
                  end: Optional[int]) -> Optional[AnalysisResult]:
    """A result explaining an empty Start/End window (e.g. Start after End), else None."""
    if frames and len(frames) > 1 and not frame_slice(frames, start, end):
        last = len(frames) if end is None else end
        return AnalysisResult(summary=f"No frames between {start} and {last} "
                                      f"(the trajectory has {len(frames)})")
    return None


def run_rdf(mol: Molecule, el1: str, el2: str, name1: str = ALL_ATOMS, name2: str = ALL_ATOMS,
            r_max: float = 13.0, n_bins: int = 100, frames: Optional[list] = None,
            start: int = 1, end: Optional[int] = None, use_traj: bool = True,
            progress: Optional[Callable[[float], None]] = None) -> AnalysisResult:
    n1, n2 = _opt(name1), _opt(name2)
    label1, label2 = n1 or el1, n2 or el2
    traj = use_traj and frames and len(frames) > 1
    if traj:
        empty = _empty_window(frames, start, end)
        if empty:
            return empty
        sub = frame_slice(frames, start, end)
        r, g = compute_rdf_trajectory(sub, el1, el2, n_bins, r_max, progress, name1=n1, name2=n2)
        first = max(start, 1)
        title = f"g(r) {label1}-{label2} (frames {first}-{first + len(sub) - 1})"
    else:
        r, g = compute_rdf(mol, el1, el2, n_bins, r_max, progress, name1=n1, name2=n2)
        title = f"g(r) {label1}-{label2}"
    if len(r) == 0:
        return AnalysisResult(summary=f"No {label1}-{label2} pairs found")
    return AnalysisResult(r, g, title, "r (A)",
                          f"max g(r) = {g.max():.2f} at r = {r[g.argmax()]:.2f} A", "r(A),g(r)")


def run_coordination(mol: Molecule, center: str, neighbor: str, r_max: float = 5.0,
                     frames: Optional[list] = None, start: int = 1, end: Optional[int] = None,
                     use_traj: bool = False) -> AnalysisResult:
    if use_traj and frames and len(frames) > 1:
        empty = _empty_window(frames, start, end)
        if empty:
            return empty
        r, cn, avg = compute_coordination_trajectory(frame_slice(frames, start, end), center,
                                                     neighbor, r_max)
    else:
        r, cn, avg = compute_coordination(mol, center, neighbor, r_max)
    if len(r) == 0:
        return AnalysisResult(summary=f"No {center} or {neighbor} atoms found")
    return AnalysisResult(r, cn, f"N({neighbor}) around {center}", "r (A)",
                          f"Avg coordination at {r_max} A: {avg:.2f} {neighbor} per {center}",
                          "r_A,coordination_number", y_max=_grow(cn))


def run_hbonds(mol: Molecule, cutoff: float = 3.5, angle: float = 120.0,
               frames: Optional[list] = None, start: int = 1,
               end: Optional[int] = None) -> AnalysisResult:
    if frames and len(frames) > 1:
        empty = _empty_window(frames, start, end)
        if empty:
            return empty
        idx, counts = count_hbonds_trajectory(frame_slice(frames, start, end), cutoff, angle)
        text = ""
    else:
        hb = detect_hbonds(mol, cutoff, angle)
        idx, counts = np.array([0.0]), np.array([float(len(hb))])
        text = "\n".join(["donor      H  acceptor   D...A (A)   angle (deg)"] + [
            f"{d + 1:>5} {h + 1:>6} {a + 1:>9}   {dist:9.3f}   {ang:11.1f}"
            for d, h, a, dist, ang in hb])
    if len(idx) == 0:
        return AnalysisResult(summary="No frames")
    return AnalysisResult(idx, counts, "H-bonds vs frame", "frame",
                          f"Avg: {counts.mean():.1f} H-bonds/frame  |  max: {counts.max():.0f}"
                          f"  |  min: {counts.min():.0f}",
                          "frame,hbond_count", text=text, y_max=_grow(counts))


def run_msd(frames: list, element: str, name: str = ALL_ATOMS, dt: float = 1.0, start: int = 1,
            end: Optional[int] = None) -> AnalysisResult:
    nm = _opt(name)
    empty = _empty_window(frames, start, end)
    if empty:
        return empty
    lags, msd = compute_msd(frame_slice(frames, start, end), element=element, name=nm)
    if len(lags) == 0:
        return AnalysisResult(summary=f"No {element} atoms found")
    summary = "Done"
    D = None
    if len(lags) > 10:
        D, slope = compute_diffusion_constant(lags, msd, dt=dt)
        summary = f"D = {D:.4e} cm^2/s  (slope = {slope:.4f} A^2/ps)"
    return AnalysisResult(lags * dt, msd, f"MSD {nm or element}", "time (ps)", summary,
                          "lag_frames,msd_A2", y_max=_grow(msd), extra={"D": D},
                          csv_data=(lags, msd),
                          text=mc_trajectory_warning(frames, "a diffusion constant"))


def run_rmsd(frames: list, element: str = ALL_ATOMS, name: str = ALL_ATOMS, dt: float = 1.0,
             start: int = 1, end: Optional[int] = None) -> AnalysisResult:
    el, nm = _opt(element), _opt(name)
    empty = _empty_window(frames, start, end)
    if empty:
        return empty
    idx, vals = compute_rmsd(frame_slice(frames, start, end), element=el, name=nm)
    if len(idx) == 0:
        return AnalysisResult(summary="No matching atoms found")
    avg = float(vals[len(vals) // 2:].mean())
    return AnalysisResult(idx * dt, vals, f"RMSD ({el or 'all'})", "time (ps)",
                          f"RMSD avg (2nd half) = {avg:.3f} A  |  max = {vals.max():.3f} A",
                          "frame,rmsd_A", y_max=_grow(vals), csv_data=(idx, vals),
                          text=mc_trajectory_warning(frames, "an RMSD against time"))


def run_gyration(frames: list, element: str = ALL_ATOMS, dt: float = 1.0, start: int = 1,
                 end: Optional[int] = None) -> AnalysisResult:
    el = _opt(element)
    empty = _empty_window(frames, start, end)
    if empty:
        return empty
    idx, rg = compute_gyration_trajectory(frame_slice(frames, start, end), el)
    if len(idx) == 0:
        return AnalysisResult(summary="No data")
    return AnalysisResult(idx * dt, rg, f"Rg ({el or 'all'})", "time (ps)",
                          f"Avg Rg = {rg.mean():.3f} A  |  std = {rg.std():.3f} A", "frame,rg_A",
                          y_min=float(rg.min()) * 0.9 if rg.min() > 0 else 0.0, y_max=_grow(rg),
                          csv_data=(idx, rg))


def run_density_profile(mol: Molecule, element: str, name: str = ALL_ATOMS, axis: int = 2,
                        frames: Optional[list] = None, start: int = 1,
                        end: Optional[int] = None, use_traj: bool = True,
                        n_bins: int = 100) -> AnalysisResult:
    nm = _opt(name)
    axis_label = "abc"[axis]
    traj = use_traj and frames and len(frames) > 1
    if traj:
        empty = _empty_window(frames, start, end)
        if empty:
            return empty
        sub = frame_slice(frames, start, end)
        pos, dens = compute_density_profile_trajectory(sub, axis, n_bins, element, nm)
        title = f"Density {element} ({axis_label}-axis, {len(sub)} frames)"
    else:
        pos, dens = compute_density_profile(mol, axis, n_bins, element, nm)
        title = f"Density {element} ({axis_label}-axis)"
    if len(pos) == 0:
        return AnalysisResult(summary=f"No {element} atoms found")
    return AnalysisResult(pos, dens, title, f"position along {axis_label} (A)",
                          f"max density = {dens.max():.2f} atoms/A at {pos[dens.argmax()]:.2f} A",
                          "position_A,density", y_max=_grow(dens))


def _pxrd_frame_average(frames: list, wavelength: float, two_theta_max: float,
                        peak_width: float,
                        progress: Optional[Callable[[float], None]]) -> tuple:
    """(2-theta, pattern, frames used): the mean of the frames' un-normalised patterns,
    each with its own cell, scaled so the highest point is 100. Frames without a real
    cell are left out."""
    used: list = []
    x, total = np.array([]), None
    for i, frame in enumerate(frames):
        if frame.pbc is not None and not is_fallback_cell(frame.pbc):
            used.append(frame)
            fx, fy = pxrd_raw_pattern(frame, wavelength=wavelength,
                                      two_theta_max=two_theta_max, peak_width=peak_width)
            if len(fy):
                x, total = fx, (fy if total is None else total + fy)
        if progress:
            progress((i + 1) / len(frames))
    if total is None:
        return np.array([]), np.array([]), used
    mean = total / len(used)  # frames with a cell but no reflections count as zero
    if mean.max() > 0:
        mean = mean / mean.max() * 100
    return x, mean, used


def run_pxrd(mol: Molecule, wavelength: float = 1.5406, two_theta_max: float = 50.0,
             peak_width: float = 0.10,
             progress: Optional[Callable[[float], None]] = None,
             frames: Optional[list] = None, start: int = 1, end: Optional[int] = None,
             use_traj: bool = False) -> AnalysisResult:
    wl_name = "Cu Ka" if abs(wavelength - 1.5406) < 0.01 else f"{wavelength:.4f} A"
    title, note, generic = f"PXRD ({wl_name})", "", None
    if use_traj and frames and len(frames) > 1:
        empty = _empty_window(frames, start, end)
        if empty:
            return empty
        sub = frame_slice(frames, start, end)
        first = max(start, 1)
        span = f"frames {first}-{first + len(sub) - 1}"
        tt, intensity, used = _pxrd_frame_average(sub, wavelength, two_theta_max, peak_width,
                                                  progress)
        if not used:
            return AnalysisResult(summary=f"No unit cell in {span}: "
                                          f"PXRD needs a crystal structure")
        title = f"PXRD ({wl_name}, {span})"
        if len(used) < len(sub):
            note = (f"  |  {len(sub) - len(used)} of {len(sub)} frames have no unit cell "
                    f"and were left out")
        generic = sorted(set().union(*(pxrd_untabulated_elements(f) for f in used)))
    else:
        tt, intensity = compute_pxrd(mol, wavelength=wavelength, two_theta_max=two_theta_max,
                                     peak_width=peak_width, progress_callback=progress)
    if len(tt) == 0:
        return AnalysisResult(summary="No reflections found")
    y = np.asarray(intensity)  # peaks: local maxima of the curve, not points above 5%
    n_peaks = int(((y[1:-1] > 5) & (y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:])).sum())
    summary = f"{n_peaks} peaks above 5% intensity"
    if generic is None:
        generic = pxrd_untabulated_elements(mol)
    if generic:
        summary += (f"  |  no tabulated form factor for {', '.join(generic)}: "
                    f"intensities approximate")
    return AnalysisResult(tt, intensity, title, "2-theta (deg)",
                          summary + note, "two_theta_deg,intensity",
                          y_max=105.0)


def read_experimental_pxrd(path: str) -> tuple:
    """Two-column (2-theta, intensity) pattern, scaled to a maximum of 100.

    Columns may be separated by commas, semicolons, tabs or spaces. '#' comments, blank
    lines and text header lines are skipped; extra columns are ignored.
    """
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            fields = [f for f in line.replace(",", " ").replace(";", " ").split() if f]
            try:
                values = [float(f) for f in fields[:2]]
            except ValueError:
                continue  # a header line
            if len(values) == 2:
                rows.append(values)
    if not rows:
        raise ValueError("need 2 columns (2theta, intensity)")
    data = np.array(rows)
    tt, inten = data[:, 0], data[:, 1]
    if inten.max() > 0:
        inten = inten / inten.max() * 100
    return tt, inten


def run_density_3d(frames: list, element: str, n_bins: int = 15,
                   progress: Optional[Callable[[float], None]] = None) -> AnalysisResult:
    positions, values, _ = compute_density_3d(frames, element, n_bins=n_bins,
                                              progress_callback=progress)
    if len(positions) == 0:
        return AnalysisResult(summary=f"No {element} atoms found in trajectory")
    return AnalysisResult(summary=f"3D density: {len(positions)} bins for {element}",
                          extra={"positions": positions, "values": values, "element": element})


# ----------------------------------------------------------------------
# Files: MPMC energy.dat and isotherm folders
# ----------------------------------------------------------------------

def find_energy_file(structure_path: str) -> Optional[str]:
    """energy.dat or energy.csv next to the loaded structure, if present."""
    base = Path(structure_path).parent if structure_path else Path(".")
    for cand in (base / "energy.dat", base / "energy.csv"):
        if cand.exists():
            return str(cand)
    return None


def energy_columns(data: Dict) -> tuple:
    """(plottable columns, default column): everything but 'step'; N, then energy, preferred."""
    cols = [k for k in data if k != "step"]
    default = cols[0] if cols else "N"
    for prefer in ("N", "energy"):
        if prefer in cols:
            default = prefer
            break
    return cols, default


def energy_n_steps(data: Dict) -> int:
    for k in ("step", "N", "energy"):
        if k in data:
            return len(data[k])
    return len(next(iter(data.values()))) if data else 0


def energy_slice(data: Dict, start: int = 1, end: Optional[int] = None) -> slice:
    """1-based inclusive Start/End clamped to the data, as a slice."""
    n = energy_n_steps(data)
    if n == 0:
        return slice(0, 0)
    s = max(0, min(start - 1, n - 1))
    e = n if end is None else end
    return slice(s, max(s + 1, min(e, n)))


def run_energy(data: Dict, column: str, start: int = 1,
               end: Optional[int] = None) -> AnalysisResult:
    if column not in data:
        return AnalysisResult(summary=f"No column {column!r}")
    sl = energy_slice(data, start, end)
    steps = data.get("step", np.arange(len(data[column])))[sl]
    values = data[column][sl]
    if len(values) == 0:
        return AnalysisResult(summary="Empty step range")
    finite = values[np.isfinite(values)]
    if not len(finite):
        return AnalysisResult(summary=f"{column}: no finite values in steps "
                                      f"{sl.start + 1}-{sl.stop}")
    lo, hi = float(finite.min()), float(finite.max())
    y_min = lo * (0.9 if lo > 0 else 1.1)
    y_max = hi * 1.1 if hi > 0 else hi * 0.9
    if y_min == y_max:
        y_max = y_min + 1
    title = column
    if sl.start != 0 or sl.stop != energy_n_steps(data):
        title = f"{column} (steps {sl.start + 1}-{sl.stop})"
    tail = values[-len(values) // 4:] if len(values) > 4 else values
    avg = float(np.nanmean(tail)) if np.isfinite(tail).any() else float("nan")
    return AnalysisResult(steps, values, title, "step",
                          f"{column} [{sl.start + 1}-{sl.stop}]: "
                          f"last-quarter avg = {avg:.4f}, final = {values[-1]:.4f}",
                          ",".join(data.keys()), y_min=y_min, y_max=y_max,
                          extra={"slice": sl})


def write_energy_csv(data: Dict, path: str, sl: slice) -> None:
    """All columns over the chosen steps (the energy plot's CSV export)."""
    cols = list(data.keys())
    with open(path, "w") as f:
        f.write(",".join(cols) + "\n")
        for i in range(sl.start, sl.stop):
            f.write(",".join(f"{data[c][i]:.6f}" for c in cols) + "\n")


def run_isotherm(directory: str) -> AnalysisResult:
    results = parse_isotherm_results(directory)
    if not results:
        return AnalysisResult(summary="No completed pressure points found")
    p = np.array([r[0] for r in results])
    u = np.array([r[1] for r in results])
    return AnalysisResult(p, u, "Adsorption Isotherm", "pressure (atm)",
                          f"{len(p)} pressure points from {directory}",
                          "pressure_atm,avg_molecules", y_max=_grow(u),
                          text=isotherm_log_table(directory))


def write_csv(path: str, header: str, *arrays: Sequence[float]) -> None:
    """Parallel arrays as CSV with six decimals (the viewers' CSV export)."""
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in zip(*arrays):
            f.write(",".join(f"{v:.6f}" for v in row) + "\n")

# ======================================================================
# Module: inputgen
# ======================================================================
"""Simulation-input generation shared by the terminal viewer and the desktop GUI.

The forms come from ``templates.ENGINES``; this module turns form values into files and
produces the status line and preview text both front ends show.
"""




_FOLDERS = {
    "mpmc_isotherm": "isotherm",
    "cp2k_cell_opt": "cp2k_cell_opt",
    "cp2k_dft_md": "cp2k_md",
    "openmm_npt": "openmm_npt",
    "openmm_nvt": "openmm_nvt",
}


@dataclass
class GenOutcome:
    """What to show after an action. ``preview`` None leaves the preview as it was."""

    ok: bool
    status: str
    preview: Optional[str] = None


def _plural_s(n: int, word: str) -> str:
    return f"{n} {word}{'s' if n != 1 else ''}"


def default_engine_outdir(filepath: str, engine: str) -> str:
    """Folder next to the loaded file (or the current one), named after the engine."""
    base = Path(filepath).parent if filepath else Path(".")
    return str(base / _FOLDERS.get(engine, "sim_input"))


def parse_form_value(var: VarSpec, raw: Any) -> Any:
    """A form entry converted to the variable's type; empty numbers raise ValueError."""
    if var.kind == "int":
        return int(raw)
    if var.kind == "float":
        text = str(raw).strip()
        if not text:
            raise ValueError(f"{var.label} is empty")
        return float(text)
    if var.kind == "csv_floats":
        parts = [p.strip() for p in str(raw).split(",") if p.strip()]
        if not parts:
            raise ValueError(f"{var.label} is empty")
        return [float(p) for p in parts]
    return str(raw)


def format_validation(errors: List[str], warnings: List[str]) -> str:
    if not errors and not warnings:
        return ("✓ All checks passed\n\n"
                "Charges, LJ params, box size, and sorbate placement all look fine.")
    lines: List[str] = []
    if errors:
        lines += [f"✗ {_plural_s(len(errors), 'error')} (must fix before generating):", ""]
        lines += [f"  ✗ {e}" for e in errors] + [""]
    if warnings:
        lines += [f"⚠ {_plural_s(len(warnings), 'warning')} (generation will proceed):", ""]
        lines += [f"  ⚠ {w}" for w in warnings]
    return "\n".join(lines)


def validate_mpmc_inputs(mol: Molecule,
                         sorbate: Optional[str] = None) -> Tuple[List[str], List[str]]:
    return validate_for_mpmc(mol, sorbate_model=sorbate)


def check_engine(mol: Molecule, engine: str,
                 values: Optional[Dict[str, Any]] = None) -> GenOutcome:
    """The Validate button: MPMC's full checks, or the periodic-box check for the others.

    ``values`` are the form's current values (the MPMC checks depend on the sorbate)."""
    if engine not in ENGINES:
        return GenOutcome(False, f"Unknown engine: {engine}")
    if engine == "mpmc_isotherm":
        errors, warnings = validate_mpmc_inputs(mol, (values or {}).get("sorbate"))
        preview = format_validation(errors, warnings)
        if errors:
            return GenOutcome(False, f"{_plural_s(len(errors), 'error')} — fix before Generate",
                              preview)
        if warnings:
            return GenOutcome(True, f"{_plural_s(len(warnings), 'warning')} — review then Generate",
                           preview)
        return GenOutcome(True, "Ready to generate", preview)
    problems = []
    if engine.startswith("cp2k") and mol.pbc is None:
        problems.append("CP2K requires a periodic box. Add one in Update Unit Cell.")
    if engine == "openmm_npt" and mol.pbc is None:
        problems.append("OpenMM NPT requires a periodic box.")
    if problems:
        return GenOutcome(False, "Fix issues before Generate", "✗ " + "\n✗ ".join(problems))
    cautions = engine_cautions(mol, engine, values)
    if cautions:
        n = _plural_s(len(cautions), "warning")
        return GenOutcome(True, f"{n} — review then Generate", "⚠ " + "\n\n⚠ ".join(cautions))
    return GenOutcome(True, "Ready to generate", "✓ No engine-level problems detected.")


def engine_cautions(mol: Molecule, engine: str,
                    values: Optional[Dict[str, Any]] = None) -> List[str]:
    """Things the CP2K and OpenMM templates do not handle for this structure. They do not
    stop generation; they say what to change in the generated files."""
    out: List[str] = []
    if engine.startswith("cp2k"):
        metals = open_shell_metals(mol)
        if metals:
            out.append(
                f"{', '.join(metals)} usually has unpaired electrons, and the generated input "
                f"is closed-shell. Add UKS to &DFT and set MULTIPLICITY (and check the "
                f"charge) before running, or the SCF will converge to the wrong state or "
                f"not at all.")
    if engine.startswith("openmm"):
        ff = str((values or {}).get("forcefield") or engine_defaults(engine)["forcefield"])
        warning = openmm_template_warning(mol, ff)
        if warning:
            out.append(warning)
    return out


def engine_template_text(engine: str) -> GenOutcome:
    if engine not in ENGINES:
        return GenOutcome(False, f"Unknown engine: {engine}")
    spec = ENGINES[engine]
    try:
        tpl = load_template(spec.template_file)
    except Exception as e:  # missing or unreadable template file
        return GenOutcome(False, f"Error loading template: {e}")
    return GenOutcome(True, f"Template: {spec.template_file} "
                         f"(edit the file in src/pdb_wizard/templates/ to change defaults)",
                   f"-- {spec.template_file} --\n\n{tpl}")


def generate_engine_inputs(mol: Molecule, engine: str, values: Dict[str, Any],
                           outdir: str) -> GenOutcome:
    """Write the input files for ``engine`` into ``outdir``."""
    if not outdir:
        return GenOutcome(False, "Output dir is required")
    try:
        if engine == "mpmc_isotherm":
            return _generate_mpmc(mol, values, outdir)
        if engine in ("cp2k_cell_opt", "cp2k_dft_md"):
            return _generate_cp2k(mol, values, outdir, engine.split("_", 1)[1])
        if engine in ("openmm_npt", "openmm_nvt"):
            return _generate_openmm(mol, values, outdir, engine[-3:].upper())
    except Exception as e:  # reported in the status line
        return GenOutcome(False, f"Error: {e}")
    return GenOutcome(False, f"Unknown engine: {engine}")


def _generate_mpmc(mol: Molecule, values: Dict[str, Any], outdir: str) -> GenOutcome:
    errors, warnings = validate_mpmc_inputs(mol, values.get("sorbate"))
    if errors:  # refuse; warnings are allowed
        return GenOutcome(False, f"Refusing: {_plural_s(len(errors), 'error')} "
                              f"would produce a meaningless run",
                       format_validation(errors, warnings))
    dirs = generate_isotherm(
        mol, sorbate_model=values["sorbate"], temperature=values["temperature"],
        pressures=values["pressures"], output_dir=outdir, numsteps=values["numsteps"],
        corrtime=values["corrtime"], ensemble=values["ensemble"],
        insert_probability=values["insert_probability"], move_factor=values["move_factor"],
        rot_factor=values["rot_factor"],
    )
    lines: List[str] = []
    if warnings:
        lines.append(f"⚠ Proceeded with {_plural_s(len(warnings), 'warning')}:")
        lines += [f"  ⚠ {w}" for w in warnings] + [""]
    lines += [f"Generated {len(dirs)} pressure points:", f"  Sorbate: {values['sorbate']}",
              f"  Temperature: {values['temperature']} K", f"  Steps: {values['numsteps']}",
              "", "Directories:"]
    for d in dirs:
        lines += [f"  {d}/", "    input.pqr, mpmc.inp, insert.pqr, run.sh"]
    lines += ["", f"Master script: {outdir}/run_all.sh", "",
              f"To run: cd {outdir} && bash run_all.sh"]
    tail = f" ({len(warnings)} warnings)" if warnings else ""
    return GenOutcome(True, f"Generated {len(dirs)} directories in {outdir}/{tail}",
                      "\n".join(lines))


def _caution_block(mol: Molecule, engine: str, values: Dict[str, Any]) -> str:
    cautions = engine_cautions(mol, engine, values)
    return ("\n\n⚠ " + "\n\n⚠ ".join(cautions)) if cautions else ""


def _generate_cp2k(mol: Molecule, values: Dict[str, Any], outdir: str, mode: str) -> GenOutcome:
    if mol.pbc is None:
        return GenOutcome(False, "CP2K requires a periodic box")
    kwargs = dict(project_name=values["project"], cutoff_Ry=values["cutoff_ry"],
                  rel_cutoff_Ry=values["rel_cutoff_ry"], xc_functional=values["xc_functional"],
                  eps_scf=values["eps_scf"], max_scf=values["max_scf"])
    if mode == "dft_md":
        kwargs.update(temperature_K=values["temperature"], timestep_fs=values["timestep_fs"],
                      n_md_steps=values["n_md_steps"], md_traj_interval=values["md_traj_interval"])
    else:
        kwargs["n_opt_steps"] = values["n_opt_steps"]
    result = generate_cp2k_input(mol, outdir, mode=mode, **kwargs)
    return GenOutcome(True, f"Wrote {result['input']}",
                   f"Wrote CP2K {mode} input:\n  {result['input']}\n  {result['xyz']}\n\n"
                   f"Elements: {', '.join(result['elements'])}\n\n"
                   f"To run: cd {outdir} && cp2k -i {Path(result['input']).name} -o out.log"
                   + _caution_block(mol, f"cp2k_{mode}", values))


def _generate_openmm(mol: Molecule, values: Dict[str, Any], outdir: str,
                     ensemble: str) -> GenOutcome:
    kwargs = dict(forcefield=values["forcefield"], temperature_K=values["temperature_k"],
                  timestep_fs=values["timestep_fs"], n_steps=values["n_steps"],
                  report_interval=values["report_interval"],
                  nonbonded_cutoff_nm=values["nonbonded_cutoff_nm"])
    if ensemble == "NPT":
        kwargs["pressure_atm"] = values["pressure_atm"]
    result = generate_openmm_script(mol, outdir, ensemble=ensemble, **kwargs)
    return GenOutcome(True, f"Wrote {result['script']}",
                   f"Wrote OpenMM {ensemble} driver:\n  {result['script']}\n  {result['pdb']}\n\n"
                   f"To run: cd {outdir} && python {Path(result['script']).name}"
                   + _caution_block(mol, f"openmm_{ensemble.lower()}", values))

# ======================================================================
# Module: formats
# ======================================================================
"""Additional file format readers/writers: POSCAR, LAMMPS data, MOL2, CIF write."""



import numpy as np


# ---------------------------------------------------------------------------
# POSCAR / CONTCAR (VASP)
# ---------------------------------------------------------------------------

def read_poscar(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:
    """Read VASP POSCAR/CONTCAR format.

    The cell is returned in the tool's orientation (a along x, b in the xy plane); when
    the file's vectors point elsewhere the atoms are turned with them. The scaling factor
    applies to the lattice and to Cartesian positions; a negative one is the cell volume
    and three values scale x, y and z.
    """
    try:
        return _read_poscar_body(file)
    except IndexError:
        raise ValueError("The POSCAR ends before all of its lattice vectors, species "
                         "counts and atom positions were read") from None


def _read_poscar_body(file: TextIO) -> tuple[list[Atom], Optional[PBC]]:

    comment = file.readline().strip()
    scale_tokens = file.readline().split()
    if len(scale_tokens) not in (1, 3):
        raise ValueError("POSCAR line 2 must hold one scaling factor or three")
    scales = [float(t) for t in scale_tokens]

    # Lattice vectors
    vecs = []
    for _ in range(3):
        tokens = file.readline().split()
        vecs.append([float(t) for t in tokens[:3]])
    raw = np.array(vecs)
    if len(scales) == 3:
        scale = np.array(scales)  # one factor per Cartesian axis
    elif scales[0] < 0:
        scale = (abs(scales[0]) / abs(np.linalg.det(raw))) ** (1.0 / 3.0)  # target volume
    else:
        scale = scales[0]
    basis = raw * scale

    pbc, turn = lattice_to_standard(basis)

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
            else:
                coords = coords * scale
            if turn is not None:
                coords = coords @ turn
            atom = Atom(coords[0], coords[1], coords[2], sp)
            system.append(atom)

    set_atom_ids(system)
    return system, pbc


def write_poscar(mol: Molecule, out: TextIO) -> None:
    """Write VASP POSCAR format. Raises CellRequiredError without a cell."""
    require_cell(mol, "Writing a POSCAR")

    from collections import OrderedDict

    sites = export_sites(mol)
    species_counts: OrderedDict[str, int] = OrderedDict()
    for _, sym in sites:
        species_counts[sym] = species_counts.get(sym, 0) + 1

    out.write("Generated by PDB Wizard\n")
    out.write("1.0\n")
    for row in mol.pbc.basis_matrix:
        out.write(f"  {row[0]:>16.10f} {row[1]:>16.10f} {row[2]:>16.10f}\n")
    out.write(" ".join(species_counts.keys()) + "\n")
    out.write(" ".join(str(c) for c in species_counts.values()) + "\n")
    out.write("Cartesian\n")
    # VASP gives the first N1 positions to the first species, the next N2 to the second:
    # the positions have to be grouped the way the counts are.
    for species in species_counts:
        for a, sym in sites:
            if sym == species:
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
    mass_symbols: dict[int, str] = {}
    in_atoms = False
    in_masses = False
    style_hint = ""
    tilt = None

    for line in lines:
        # a '#' starts a comment anywhere on a line (the hint after 'Atoms' is kept)
        stripped, _, comment = line.partition("#")
        stripped = stripped.strip()
        comment = comment.strip()
        if not stripped:
            continue
        # Section headers ('Masses', 'Atoms', 'Angle Coeffs', ...) start with a letter;
        # header and data lines start with a number.
        if stripped[0].isalpha():
            in_atoms = stripped.startswith("Atoms")
            in_masses = stripped.startswith("Masses")
            if in_atoms and comment:
                style_hint = comment.split()[0].lower()
            continue
        if in_atoms:
            atom_lines.append(stripped)
        elif in_masses:
            tokens = stripped.split()
            if len(tokens) >= 2:
                try:
                    masses[int(tokens[0])] = float(tokens[1])
                except ValueError:
                    continue
                if comment and comment.split()[0] in ELEMENTS and comment.split()[0] != "X":
                    mass_symbols[int(tokens[0])] = comment.split()[0]
        elif "xy" in stripped and "xz" in stripped and "yz" in stripped:
            tokens = stripped.split()
            tilt = (float(tokens[0]), float(tokens[1]), float(tokens[2]))
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

    if tilt is not None and any(abs(t) > 1e-12 for t in tilt):
        # triclinic box: a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz)
        pbc, _ = lattice_to_standard([[xhi - xlo, 0.0, 0.0], [tilt[0], yhi - ylo, 0.0],
                                      [tilt[1], tilt[2], zhi - zlo]])
    else:
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
        mass_to_el[type_id] = mass_symbols.get(type_id, best_el)

    # Atom style: the 'Atoms # style' hint when the file has one, else the token count of
    # the first data line. Three image flags may follow the coordinates.
    # atomic: id type x y z          (5 tokens, 8 with image flags)
    # charge: id type charge x y z   (6 tokens, 9 with image flags)
    # full:   id mol-id type charge x y z  (7 tokens, 10 with image flags)
    style = "atomic"
    if atom_lines:
        n_tok = len(atom_lines[0].split())
        if n_tok == 8:
            style = "atomic"
        elif n_tok == 9:
            style = "charge"
        elif n_tok >= 7:
            style = "full"
        elif n_tok >= 6:
            style = "charge"
    if style_hint in ("atomic", "charge", "full"):
        style = style_hint
    elif style_hint in ("molecular", "bond", "angle"):
        style = "molecular"

    system: list[Atom] = []
    for line in atom_lines:
        tokens = line.split()
        try:
            if style == "full":
                # id mol-id type charge x y z
                atom_type = int(tokens[2])
                charge = float(tokens[3])
                x, y, z = float(tokens[4]), float(tokens[5]), float(tokens[6])
            elif style == "molecular":
                # id mol-id type x y z
                atom_type = int(tokens[2])
                charge = 0.0
                x, y, z = float(tokens[3]), float(tokens[4]), float(tokens[5])
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
    """Write a CIF file. Raises CellRequiredError without a cell."""
    require_cell(mol, "Writing a CIF")

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
    sites = export_sites(mol)
    coords = np.array([a.x for a, _ in sites]).reshape(-1, 3)
    frac = coords @ pbc.reciprocal_basis_matrix

    label_counts: dict[str, int] = {}
    for i, (atom, sym) in enumerate(sites):
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
    vectors: list = []  # 'Tv' lines: translation vectors of a periodic input
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            break
        tokens = line.replace(",", " ").split()
        if len(tokens) >= 4:
            element = tokens[0]
            xyz = tokens[1:4]
            if len(tokens) >= 5 and _is_int_literal(tokens[1]) and _is_float_literal(tokens[4]):
                xyz = tokens[2:5]  # 'C -1 x y z': the integer is the freeze code
            x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
            if element.lower() == "tv":
                vectors.append([x, y, z])
            elif element.isdigit():
                system.append(Atom(x, y, z, get_element_by_number(int(element)).symbol))
            else:
                system.append(Atom(x, y, z, element))
        idx += 1

    pbc = None
    if len(vectors) == 3:
        pbc, turn = lattice_to_standard(vectors)
        if turn is not None:
            for atom in system:
                atom.x = atom.x @ turn
    set_atom_ids(system)
    return system, pbc


def _is_int_literal(token: str) -> bool:
    return token.lstrip("+-").isdigit()


def _is_float_literal(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


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
    if "V3000" in counts_line:
        raise ValueError("This is a V3000 SDF/MOL file, which is not supported. Save it "
                         "as V2000 (the default of most programs) and open that.")
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
            # Element is the part before the dot in atom_type (e.g. "C.3" -> "C").
            # antechamber writes lower-case GAFF types instead ('ca' is an aromatic
            # carbon, 'os' an ether oxygen): their first letter is the element, except
            # for cl and br.
            element = atom_type.split(".")[0]
            if atom_type.islower():
                element = (atom_type[:2] if atom_type[:2] in ("cl", "br")
                           else atom_type[0]).capitalize()
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


    # Atoms a LAMMPS file can hold: sorbate sites by their real element, massless model
    # sites left out (see io.export_sites). Bonds are renumbered to match.
    sites = export_sites(mol)
    atoms = [atom for atom, _ in sites]
    symbols = [sym for _, sym in sites]
    if len(atoms) == len(mol.atoms):
        bonds = list(mol.bonds)
    else:
        kept = {id(atom): k for k, atom in enumerate(atoms)}
        new_index = [kept.get(id(atom)) for atom in mol.atoms]
        bonds = [(new_index[a], new_index[b]) for a, b in mol.bonds
                 if new_index[a] is not None and new_index[b] is not None]
    unknown = sorted({atom.name for atom, sym in sites
                      if sym not in ELEMENTS or ELEMENTS[sym].mass <= 0})
    if unknown:
        raise ValueError(
            "A LAMMPS data file needs an element with a mass for every atom; these atom "
            f"names have none: {', '.join(unknown)}")

    # Build atom type mapping: element symbol -> type id
    type_map: OrderedDict[str, int] = OrderedDict()
    for sym in symbols:
        if sym not in type_map:
            type_map[sym] = len(type_map) + 1

    n_atoms = len(atoms)
    n_bonds = len(bonds)
    n_atom_types = len(type_map)

    # Assign bond types by element pair
    bond_type_map: dict[tuple[str, str], int] = {}
    for a1, a2 in bonds:
        pair = tuple(sorted([symbols[a1], symbols[a2]]))
        if pair not in bond_type_map:
            bond_type_map[pair] = len(bond_type_map) + 1
    n_bond_types = len(bond_type_map) if n_bonds > 0 else 0

    # Assign molecule IDs via connected components
    mol_ids = [0] * n_atoms
    neighbors: dict[int, list[int]] = {}
    for a, b in bonds:
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

    # Box dimensions. LAMMPS: a = (xhi-xlo, 0, 0), b = (xy, yhi-ylo, 0), c = (xz, yz,
    # zhi-zlo), which is the orientation of pbc.basis_matrix.
    tilt = None
    if mol.pbc is not None:
        xlo, ylo, zlo = 0.0, 0.0, 0.0
        basis = mol.pbc.basis_matrix
        xhi = float(basis[0][0])
        yhi = float(basis[1][1])
        zhi = float(basis[2][2])
        if max(abs(basis[1][0]), abs(basis[2][0]), abs(basis[2][1])) > 1e-9:
            tilt = (float(basis[1][0]), float(basis[2][0]), float(basis[2][1]))
    else:
        coords = np.array([a.x for a in atoms]).reshape(-1, 3)
        if not len(coords):
            coords = np.zeros((1, 3))
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
    if tilt is not None:
        out.write(f"{tilt[0]:.6f} {tilt[1]:.6f} {tilt[2]:.6f} xy xz yz\n")

    # Masses section
    out.write("\nMasses\n\n")
    for sym, type_id in type_map.items():
        mass = ELEMENTS[sym].mass if sym in ELEMENTS else 1.0
        out.write(f"{type_id} {mass:.4f} # {sym}\n")

    # Atoms section (full style: id mol-id type charge x y z)
    out.write("\nAtoms\n\n")
    for i, atom in enumerate(atoms):
        atom_id = i + 1
        type_id = type_map[symbols[i]]
        charge = atom.charge
        out.write(f"{atom_id} {mol_ids[i]} {type_id} {charge:.6f} "
                  f"{atom.x[0]:.6f} {atom.x[1]:.6f} {atom.x[2]:.6f}\n")

    # Bonds section
    if n_bonds > 0:
        out.write("\nBonds\n\n")
        for i, (a1, a2) in enumerate(bonds):
            bond_id = i + 1
            pair = tuple(sorted([symbols[a1], symbols[a2]]))
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


def _db_file_stem(name: str) -> str:
    """A compound name as a file name: spaces and path separators become '_' ('N/A' used
    to fail and '../x' to leave the output folder)."""
    stem = name.replace(" ", "_").replace("/", "_").replace("\\", "_").lstrip(".")
    return stem or "structure"


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


def search_cod(formula: str = "", text: str = "", strict: bool = False) -> list[dict]:
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
                "formula": entry.get("formula") or "",
                "title": (entry.get("title") or "")[:60],  # COD sends null for no title
            })
        return results
    except Exception as e:
        if strict and isinstance(e, OSError):  # network down: search_backends reports it
            raise
        return []


def search_rcsb(query: str, strict: bool = False) -> list[dict]:
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
    except Exception as e:
        if strict and isinstance(e, OSError):
            raise
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
        out_path = str(Path(output_dir) / (_db_file_stem(name) + ".sdf"))
    return _http_download(url, out_path)


def search_pubchem(query: str, strict: bool = False) -> list[dict]:
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
                "formula": prop.get("MolecularFormula") or "",
                "title": (prop.get("IUPACName") or "")[:60],
            })
        return results
    except Exception as e:
        if strict and isinstance(e, OSError):
            raise
        return []


def search_materials_project(formula: str, api_key: str = "",
                             strict: bool = False) -> list[dict]:
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
    except Exception as e:
        if strict and isinstance(e, OSError):
            raise
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
        # The occupancy is written so that the CIF reader can say when a site is shared
        # by several species (only the first one is kept here).
        f.write("_atom_site_occupancy\n")
        for i, site in enumerate(sites):
            species = site.get("species", [{}])
            el = species[0].get("element", "X") if species else "X"
            occupancy = species[0].get("occu", 1) if species else 1
            abc = site.get("abc", [0, 0, 0])
            f.write(f"{el}{i+1} {el} {abc[0]:.6f} {abc[1]:.6f} {abc[2]:.6f} {occupancy}\n")
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
        return search_cod(formula=query, strict=True)
    return search_cod(text=query, strict=True)


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
        "search": lambda q, mode="text": search_rcsb(q, strict=True),
        "fetch": fetch_rcsb,
    },
    "pubchem": {
        "label": "PubChem",
        "modes": ("text", "id"),
        "placeholder": {
            "text": "aspirin",
            "id": "2244 (CID)",
        },
        "search": lambda q, mode="text": search_pubchem(q, strict=True),
        "fetch": fetch_pubchem_sdf,
    },
    "mp": {
        "label": "Materials Project",
        "modes": ("formula", "id"),
        "placeholder": {
            "formula": "Si O2",
            "id": "mp-149",
        },
        "search": lambda q, mode="formula": search_materials_project(q, strict=True),
        # fetch_materials_project_cif takes (material_id, api_key, output_dir)
        # Bind api_key='' so callers pass (id, output_dir) like the others.
        "fetch": lambda mid, output_dir=".": fetch_materials_project_cif(mid, "", output_dir),
    },
}


# Modes the "All databases" search accepts; the query goes to every backend with that mode.
ALL_MODES = ("text", "formula")


def modes_for(source: str) -> tuple:
    if source == "all":
        return ALL_MODES
    return BACKENDS.get(source, {}).get("modes", ("text",))


def placeholder_for(source: str, mode: str) -> str:
    if source == "all":
        if mode == "formula":
            return "formula across all DBs (eg 'C9H8O4')"
        return "free-text search across all DBs"
    return BACKENDS.get(source, {}).get("placeholder", {}).get(mode, "query")


def backends_for(source: str, mode: str) -> list:
    if source == "all":
        return [c for c, spec in BACKENDS.items() if mode in spec["modes"]]
    return [source]


def search_backends(source: str, mode: str, query: str) -> tuple:
    """Search one backend or all of them (in parallel). Returns (results, status message).

    A backend that fails is reported in the message; the others' results are kept.
    """
    from concurrent.futures import ThreadPoolExecutor
    codes = backends_for(source, mode)
    with ThreadPoolExecutor(max_workers=max(1, len(codes))) as pool:
        futures = {c: pool.submit(BACKENDS[c]["search"], query, mode) for c in codes}
    results: list = []
    counts: dict = {}
    errors: list = []
    for code in codes:
        try:
            out = futures[code].result()
        except Exception as e:  # network or parse failure of one backend
            errors.append(f"{BACKENDS[code]['label']} ({type(e).__name__})")
            counts[code] = 0
            continue
        results.extend(out)
        counts[code] = len(out)
    results.sort(key=lambda r: (r.get("source", ""), r.get("id", "")))
    if source == "all" and len(codes) > 1:
        breakdown = ", ".join(f"{BACKENDS[c]['label']}={counts.get(c, 0)}" for c in codes)
        msg = f"{len(results)} results  ({breakdown})"
    else:
        msg = f"{len(results)} results in {BACKENDS[codes[0]]['label']}" if codes else "0 results"
    if errors:
        msg += f"  •  errors: {', '.join(errors)}"
    return results, msg

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
    ("check",          "Check Structure",          "analysis"),
    ("symmetry",       "Find Symmetry",            "analysis"),
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
    ("pack_sorbates",  "Pack Sorbates",            "modify"),
    ("qeq_charges",    "Generate QEq Charges",     "modify"),
    # ---- File (mirrors top bar File menu) ----
    ("open_file",      "Open File...",             "export"),
    ("write_xyz",      "Save as XYZ",              "export"),
    ("write_pdb",      "Save as PDB",              "export"),
    ("write_cif",      "Save as CIF",              "export"),
    ("write_poscar",   "Save as POSCAR",           "export"),
    ("write_lammps",   "Save as LAMMPS Data",      "export"),
    ("write_crystal",  "Save as CRYSTAL23 .d12",   "export"),
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
            allow_blank=False,
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
    """Modal to extend all three axes at once. Returns (na, nb, nc), the number
    of extra copies per axis, matching the classic menu's "How many times"."""

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
            yield Label("Times to extend each axis (1 doubles it, as in the classic menu):")
            with Horizontal(classes="axis-row"):
                yield Label("a:", classes="axis-label")
                yield SpinBox(value=0, min_val=0, max_val=19, id="extend-na")
            with Horizontal(classes="axis-row"):
                yield Label("b:", classes="axis-label")
                yield SpinBox(value=0, min_val=0, max_val=19, id="extend-nb")
            with Horizontal(classes="axis-row"):
                yield Label("c:", classes="axis-label")
                yield SpinBox(value=0, min_val=0, max_val=19, id="extend-nc")
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
    """Modal to select a sorbate molecule to insert into the MPMC PDB.

    ``standalone`` (Edit > Pack Sorbates) swaps Skip for "From .pqr file", which
    dismisses with ``FROM_FILE``.
    """

    FROM_FILE = "__file__"

    def __init__(self, standalone: bool = False) -> None:
        super().__init__()
        self._standalone = standalone

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
            yield Label("Select a sorbate to insert (you choose how many next):")
            yield DataTable(id="sorbate-table", cursor_type="row")
            with Horizontal():
                yield Button("Insert", id="sorb-ok", variant="primary")
                if self._standalone:
                    yield Button("From .pqr file", id="sorb-file")
                else:
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
        elif event.button.id == "sorb-file":
            self.dismiss(self.FROM_FILE)
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
        except (OSError, UnicodeDecodeError):
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
        if x is not None and y is not None and len(x) == len(y) and len(x) > 1:
            x, y = np.asarray(x), np.asarray(y)
            if np.any(np.diff(x) < 0):
                # The plot looks samples up by x, so they must ascend (a pattern may be
                # stored from high to low angle).
                order = np.argsort(x, kind="stable")
                x, y = x[order], y[order]
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
                if not np.isfinite(y_val):
                    prev_py = None  # a gap: nothing drawn, and no line across it
                    continue
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
        """Track the task so action_close can cancel it on dismiss.

        A failure inside the task is reported in the screen's status label and
        as a toast; without this an exception would only reach stderr at
        garbage-collection time and the status would stay at "Computing...".
        """
        if self._compute_task is not None and not self._compute_task.done():
            self._compute_task.cancel()
        label = _task_label(coro)
        self._compute_task = asyncio.create_task(coro)
        self._compute_task.add_done_callback(
            lambda t, label=label: self._on_compute_done(t, label)
        )
        return self._compute_task

    def _on_compute_done(self, task: asyncio.Task, label: str) -> None:
        exc = _task_failure(task)
        if exc is None:
            return
        msg = f"{label} failed: {exc}"
        try:
            self.query_one(f"#{self._prefix}-status", Label).update(f"Error: {exc}")
        except Exception:
            pass
        try:
            self.app.notify(msg, timeout=6, severity="error")
        except Exception:
            pass

    def _get_frame_slice(self) -> list:
        """Get frame range from Start/End SpinBoxes if present."""
        return frame_slice(self._frames, *self._frame_range())

    def _frame_range(self) -> tuple:
        """(start, end), 1-based inclusive, from the Start/End SpinBoxes; (1, None) without."""
        if self._frames and len(self._frames) > 1:
            try:
                return (self.query_one(f"#{self._prefix}-start", SpinBox).value,
                        self.query_one(f"#{self._prefix}-end", SpinBox).value)
            except Exception:
                pass
        return 1, None

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
        path = filepath.strip()
        try:
            write_csv(path, header, *arrays)
        except OSError as e:
            # Raised inside a modal dismiss callback: an uncaught error here
            # would take down the whole app.
            self._safe_notify(f"Could not write {path}: {e}", severity="error")
            return
        self._safe_notify(f"Wrote {path}")

    def _safe_notify(self, message: str, severity: str = "information") -> None:
        """notify() via the app, tolerating a screen that is not mounted."""
        try:
            self.app.notify(message, timeout=6 if severity == "error" else 3,
                            severity=severity)
        except Exception:
            pass

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
                yield Select(el_opts, value=el1, id="rdf-el1", allow_blank=False)
                yield Label("  Name:")
                yield Select(name_opts1, value="(all)", id="rdf-name1", allow_blank=False)
                yield Label("    Type 2:  El:")
                yield Select(el_opts, value=el2, id="rdf-el2", allow_blank=False)
                yield Label("  Name:")
                yield Select(name_opts2, value="(all)", id="rdf-name2", allow_blank=False)
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
            self.action_close()

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

        el1 = str(self.query_one("#rdf-el1", Select).value)
        el2 = str(self.query_one("#rdf-el2", Select).value)
        name1 = str(self.query_one("#rdf-name1", Select).value)
        name2 = str(self.query_one("#rdf-name2", Select).value)
        label1 = el1 if name1 == "(all)" else name1
        label2 = el2 if name2 == "(all)" else name2

        y_max = float(self.query_one("#rdf-ymax", SpinBox).value)
        x_max = float(self.query_one("#rdf-xmax", SpinBox).value)
        n_bins = self.query_one("#rdf-bins", SpinBox).value

        use_traj = False
        try:
            use_traj = bool(self.query_one("#rdf-traj", Checkbox).value
                            and self._frames and len(self._frames) > 1)
        except Exception:
            pass
        start, end = self._frame_range() if use_traj else (1, None)

        status = self.query_one("#rdf-status", Label)
        n_total = (end - start + 1) if use_traj else 1

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

        try:
            res = await asyncio.to_thread(
                run_rdf, self._molecule, el1, el2, name1, name2, x_max, n_bins,
                self._frames, start, end, use_traj, _prog,
            )
        finally:
            poll_timer.stop()

        if res.empty:
            status.update(res.summary)
            return

        # Store for CSV export
        self._last_r = res.x
        self._last_g = res.y
        self._last_title = res.title

        plot = self.query_one("#rdf-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=0.0, y_max=y_max)
        status.update(res.summary)

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
                yield Select(el_opts, value=self._elements[-1] if self._elements else "H",
                             id="msd-el", allow_blank=False)
                yield Label("Name:")
                yield Select(name_opts, value="(all)", id="msd-name", allow_blank=False)
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
            self.action_close()
        elif event.button.id == "msd-csv":
            if len(self._last_lags) == 0:
                return
            self.app.push_screen(FileSaveModal("Save MSD as CSV:", default="msd.csv"), callback=self._export_csv)

    async def _run_msd(self) -> None:
        el = str(self.query_one("#msd-el", Select).value)
        name = str(self.query_one("#msd-name", Select).value)
        dt = float(self.query_one("#msd-dt", SpinBox).value)
        status = self.query_one("#msd-status", Label)
        status.update(f"Computing MSD for {el}...")
        await asyncio.sleep(0)
        res = await asyncio.to_thread(run_msd, self._frames, el, name, dt, *self._frame_range())
        if res.empty:
            status.update(res.summary)
            return
        self._last_lags, self._last_msd = res.csv_arrays()
        plot = self.query_one("#msd-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=0.0,
                      y_max=res.y_max)
        status.update(res.summary)
        if res.text:  # Monte Carlo trajectory: the curve is not diffusion
            self.app.notify(res.text, severity="warning", timeout=15)

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
                yield Select(el_opts, value=el1, id="coord-el1", allow_blank=False)
                yield Label("Neighbor:")
                yield Select(el_opts, value=el2, id="coord-el2", allow_blank=False)
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
            self.action_close()
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
            use_traj = bool(self.query_one("#coord-traj", Checkbox).value
                            and self._frames and len(self._frames) > 1)
        except Exception:
            pass
        if use_traj:
            start, end = self._frame_range()
            n = len(self._get_frame_slice())
            status.update(f"Computing N({el2}) around {el1} over {n} frames...")
        else:
            start, end = 1, None
            status.update(f"Computing N({el2}) around {el1}...")
        await asyncio.sleep(0)
        res = await asyncio.to_thread(
            run_coordination, self._molecule, el1, el2, r_max, self._frames, start, end, use_traj,
        )
        if res.empty:
            status.update(res.summary)
            return
        self._last_r = res.x
        self._last_cn = res.y
        plot = self.query_one("#coord-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=0.0,
                      y_max=res.y_max)
        status.update(res.summary)

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
                yield Select(el_opts, value="(all)", id="rmsd-el", allow_blank=False)
                yield Label("Name:")
                yield Select(name_opts, value="(all)", id="rmsd-name", allow_blank=False)
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
            self.action_close()
        elif event.button.id == "rmsd-csv":
            if len(self._last_idx) == 0:
                return
            self.app.push_screen(FileSaveModal("Save RMSD as CSV:", default="rmsd.csv"), callback=self._export_csv)

    async def _run_rmsd(self) -> None:
        el = str(self.query_one("#rmsd-el", Select).value)
        name = str(self.query_one("#rmsd-name", Select).value)
        dt = float(self.query_one("#rmsd-dt", SpinBox).value)
        status = self.query_one("#rmsd-status", Label)
        status.update(f"Computing RMSD for {el if el != '(all)' else 'all'}...")
        await asyncio.sleep(0)
        res = await asyncio.to_thread(run_rmsd, self._frames, el, name, dt, *self._frame_range())
        if res.empty:
            status.update(res.summary)
            return
        self._last_idx, self._last_rmsd = res.csv_arrays()
        plot = self.query_one("#rmsd-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=0.0,
                      y_max=res.y_max)
        status.update(res.summary)
        if res.text:  # Monte Carlo trajectory: the curve has no time axis
            self.app.notify(res.text, severity="warning", timeout=15)

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
            self.action_close()
        elif event.button.id == "hb-csv":
            if len(self._last_idx) == 0:
                return
            self.app.push_screen(FileSaveModal("Save H-bonds as CSV:", default="hbonds.csv"), callback=self._export_csv)

    async def _run(self) -> None:
        cutoff = self.query_one("#hb-cutoff", SpinBox).value / 10.0
        angle = float(self.query_one("#hb-angle", SpinBox).value)
        status = self.query_one("#hb-status", Label)
        if self._frames and len(self._frames) > 1:
            status.update(f"Counting H-bonds over {len(self._get_frame_slice())} frames...")
        else:
            status.update("Detecting H-bonds...")
        await asyncio.sleep(0)
        res = await asyncio.to_thread(run_hbonds, self._molecule, cutoff, angle, self._frames,
                                      *self._frame_range())
        if res.empty:
            status.update(res.summary)
            return
        self._last_idx = res.x
        self._last_counts = res.y
        plot = self.query_one("#hb-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=0.0,
                      y_max=res.y_max)
        status.update(res.summary)

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
                yield Select(el_opts, value="(all)", id="rg-el", allow_blank=False)
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
            self.action_close()
        elif event.button.id == "rg-csv":
            if len(self._last_idx) == 0:
                return
            self.app.push_screen(FileSaveModal("Save Rg as CSV:", default="rg.csv"), callback=self._export_csv)

    async def _run(self) -> None:
        el = str(self.query_one("#rg-el", Select).value)
        dt = float(self.query_one("#rg-dt", SpinBox).value)
        status = self.query_one("#rg-status", Label)
        status.update(f"Computing Rg for {el if el != '(all)' else 'all'}...")
        await asyncio.sleep(0)
        res = await asyncio.to_thread(run_gyration, self._frames, el, dt, *self._frame_range())
        if res.empty:
            status.update(res.summary)
            return
        self._last_idx, self._last_rg = res.csv_arrays()
        plot = self.query_one("#rg-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=res.y_min,
                      y_max=res.y_max)
        status.update(res.summary)

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
        table.cursor_type = "row"  # Enter / click toggles a whole row
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

        try:
            fw, sorbs = identify_framework_and_sorbates(self._molecule, na, nb, nc)
        except ValueError as e:  # not an na x nb x nc supercell
            self._framework_idx, self._sorbate_mols, self._keep = [], [], []
            self.query_one("#rc-table", DataTable).clear()
            status.update(str(e))
            return
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
            table.add_row("[X]", str(i + 1), str(len(mol_indices)), elem_str, center_str,
                          key=str(i))

        status.update(status_msg + "  |  Click rows to toggle keep/discard, then Reduce")

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        """Toggle keep/discard on row click."""
        table = self.query_one("#rc-table", DataTable)
        row_idx = event.cursor_row
        if 0 <= row_idx < len(self._keep):
            self._keep[row_idx] = not self._keep[row_idx]
            mark = "[X]" if self._keep[row_idx] else "[ ]"
            keep_column = next(iter(table.columns))
            table.update_cell(event.row_key, keep_column, mark)

    def _do_reduce(self) -> None:
        if not self._framework_idx and not self._sorbate_mols:
            self.app.notify("Click Analyze first", timeout=2)
            return
        na = self.query_one("#rc-na", SpinBox).value
        nb = self.query_one("#rc-nb", SpinBox).value
        nc = self.query_one("#rc-nc", SpinBox).value
        keep_indices = [i for i, k in enumerate(self._keep) if k]
        try:
            result = reduce_to_primitive(self._molecule, na, nb, nc, keep_indices)
        except ValueError as e:  # the sizes were changed after Analyze
            self.query_one("#rc-status", Label).update(str(e))
            return
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
                yield Select(engine_opts, value=self._engine_key, id="gen-engine",
                             allow_blank=False)
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
        return default_engine_outdir(self._filepath, self._engine_key)

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
                    allow_blank=not choices,
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
        out: dict = {}
        for var in self._engine_spec().user_vars:
            w = self._form_widgets.get(var.name)
            if w is None:
                if strict:
                    raise KeyError(f"missing form widget for {var.name!r}")
                continue
            try:
                out[var.name] = parse_form_value(var, w.value)
            except ValueError:
                if strict:
                    raise
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
    def _show_validation(self) -> None:
        self._show_outcome(check_engine(self._molecule, self._engine_key,
                                        self._collect_values(strict=False)))

    def _show_outcome(self, outcome) -> None:
        if outcome.preview is not None:
            self._set_preview(outcome.preview)
        self._set_status(outcome.status)

    # ---------------------------------------------------------------
    # Show template
    # ---------------------------------------------------------------
    def _show_template(self) -> None:
        self._show_outcome(engine_template_text(self._engine_key))

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
        outcome = generate_engine_inputs(self._molecule, self._engine_key, values, outdir)
        self._show_outcome(outcome)
        if outcome.status.startswith("Refusing"):
            self.app.notify(
                "Cannot generate — fix errors shown in the preview",
                severity="error", timeout=5,
            )

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
        label = _task_label(coro)
        self._compute_task = asyncio.create_task(coro)
        self._compute_task.add_done_callback(
            lambda t, label=label: self._on_compute_done(t, label)
        )
        return self._compute_task

    def _on_compute_done(self, task: asyncio.Task, label: str) -> None:
        exc = _task_failure(task)
        if exc is None:
            return
        try:
            self.query_one("#db-status", Label).update(f"Error: {exc}")
        except Exception:
            pass
        try:
            self.app.notify(f"{label} failed: {exc}", timeout=6, severity="error")
        except Exception:
            pass

    def _modes_for(self, source: str) -> tuple[str, ...]:
        return modes_for(source)

    def _placeholder_for(self, source: str, mode: str) -> str:
        return placeholder_for(source, mode)

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
                yield Select(db_opts, value=default_source, id="db-source", allow_blank=False)
                yield Label("Search by:")
                yield Select(mode_opts, value=default_mode, id="db-mode", allow_blank=False)
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
            self.action_close()
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
        codes = backends_for(source, mode)
        status.update(
            f"Searching {', '.join(BACKENDS[c]['label'] for c in codes)} "
            f"({mode})..."
        )
        await asyncio.sleep(0)
        results, msg = await asyncio.to_thread(search_backends, source, mode, query)
        self._results = results
        table = self.query_one("#db-results", DataTable)
        table.clear()
        for r in results:
            label = BACKENDS.get(r.get("source", ""), {}).get("label", r.get("source", ""))
            table.add_row(
                label,
                str(r.get("id") or ""),
                str(r.get("formula") or "")[:25],
                str(r.get("title") or "")[:60],
            )
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
            load_msgs: list[str] = []
            mol = await asyncio.to_thread(read_file, path, True, load_msgs.append)
            # The download sits in a per-session cache that is deleted on exit. The tab
            # keeps only the file name, so Save and Export default to the working directory.
            self.app.open_in_new_tab(mol, Path(path).name)
            self.app.notify(
                f"Opened {label} {rid} ({len(mol.atoms)} atoms)", timeout=3,
            )
            self.app.notify_overlaps_removed(load_msgs)
            self.app.notify_load_notes(mol)
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
        self._compute_task = None

    def _n_steps(self) -> int:
        """Number of rows in the energy data."""
        return energy_n_steps(self._data)

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
                yield Select(options, value=default, id="energy-col",
                             allow_blank=not options)
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
        try:
            start = self.query_one("#energy-start", SpinBox).value
            end = self.query_one("#energy-end", SpinBox).value
        except Exception:
            return energy_slice(self._data)
        return energy_slice(self._data, start, end)

    def _plot_column(self, col: str) -> None:
        if col not in self._data:
            return
        sl = self._step_slice()
        res = run_energy(self._data, col, sl.start + 1, sl.stop)
        if res.empty:
            self.query_one("#energy-status", Label).update(res.summary)
            return
        plot = self.query_one("#energy-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label,
                      y_min=res.y_min, y_max=res.y_max)
        self.query_one("#energy-status", Label).update(res.summary)

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
        path = filepath.strip()
        try:
            write_energy_csv(self._data, path, sl)
        except OSError as e:
            self.query_one("#energy-status", Label).update(f"Error: {e}")
            self.app.notify(f"Could not write {path}: {e}", timeout=6, severity="error")
            return
        self.query_one("#energy-status", Label).update(
            f"Exported steps {sl.start + 1}-{sl.stop} to {path}"
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
                yield Select(el_opts, value=self._elements[0] if self._elements else "H",
                             id="dens-el", allow_blank=False)
                yield Label("Name:")
                yield Select(name_opts, value="(all)", id="dens-name", allow_blank=False)
                yield Label("Axis:")
                yield Select(axis_opts, value="2", id="dens-axis", allow_blank=False)
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
            self.action_close()
        elif event.button.id == "dens-csv":
            if len(self._last_pos) == 0:
                return
            self.app.push_screen(FileSaveModal("Save density as CSV:", default="density.csv"), callback=self._export_csv)

    async def _run_density(self) -> None:
        el = str(self.query_one("#dens-el", Select).value)
        name = str(self.query_one("#dens-name", Select).value)
        axis = int(self.query_one("#dens-axis", Select).value)
        use_traj = False
        try:
            use_traj = bool(self.query_one("#dens-traj", Checkbox).value
                            and self._frames and len(self._frames) > 1)
        except Exception:
            pass
        start, end = self._frame_range() if use_traj else (1, None)
        status = self.query_one("#dens-status", Label)
        status.update(f"Computing density profile {el} along {'abc'[axis]}...")
        await asyncio.sleep(0)
        res = await asyncio.to_thread(run_density_profile, self._molecule, el, name, axis,
                                      self._frames, start, end, use_traj)
        if res.empty:
            status.update(res.summary)
            return
        self._last_pos = res.x
        self._last_dens = res.y
        plot = self.query_one("#dens-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=0.0,
                      y_max=res.y_max)
        status.update(res.summary)

    def _export_csv(self, filepath: str) -> None:
        self._do_export_csv(filepath, "position_A,density", self._last_pos, self._last_dens)


class PxrdScreen(_AnalysisScreen):
    """Simulated Powder XRD pattern."""
    _prefix = "pxrd"
    DEFAULT_CSS = _ANALYSIS_CSS.format(cls="PxrdScreen") + """
    PxrdScreen Select { width: 14; }
    PxrdScreen #pxrd-controls { height: auto; margin-bottom: 1; }
    PxrdScreen #pxrd-row2 { height: auto; margin-bottom: 1; }
    PxrdScreen #pxrd-status { height: 1; color: $text-muted; }
    """

    def __init__(self, molecule: Molecule, frames: list | None = None) -> None:
        super().__init__(molecule, frames)
        self._last_tt = np.array([])
        self._last_intensity = np.array([])

    def compose(self) -> ComposeResult:
        wl_opts = [("Cu Ka 1.5406", "1.5406"), ("Mo Ka 0.7107", "0.7107"), ("Co Ka 1.7903", "1.7903")]
        with Vertical():
            with Horizontal(id="pxrd-controls"):
                yield Label("Source:")
                yield Select(wl_opts, value="1.5406", id="pxrd-wl", allow_blank=False)
                yield Label("2T max:")
                yield SpinBox(value=50, min_val=10, max_val=90, id="pxrd-ttmax")
                yield Label("Width:")
                yield SpinBox(value=10, min_val=1, max_val=50, id="pxrd-width")
                yield Button("Compute", id="pxrd-compute", variant="primary")
                yield Button("Load Exp.", id="pxrd-load-exp")
                yield Button("CSV", id="pxrd-csv")
                yield Button("Close", id="pxrd-close")
            if self._frames and len(self._frames) > 1:
                with Horizontal(id="pxrd-row2"):  # own row: row 1 still fits 120 columns
                    yield Checkbox("Traj", value=True, id="pxrd-traj")
                    yield from self._yield_frame_range()
            yield PlotWidget(id="pxrd-plot")
            yield Label("Select parameters and click Compute  |  Load Exp. to overlay experimental data", id="pxrd-status")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "pxrd-compute":
            self._spawn_task(self._run_pxrd())
        elif event.button.id == "pxrd-close":
            self.action_close()
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
        use_traj = False
        try:
            use_traj = bool(self.query_one("#pxrd-traj", Checkbox).value
                            and self._frames and len(self._frames) > 1)
        except Exception:
            pass
        start, end = self._frame_range() if use_traj else (1, None)
        n_total = len(self._get_frame_slice() or []) if use_traj else 1
        status = self.query_one("#pxrd-status", Label)
        status.update("Computing PXRD pattern...")
        self._rdf_progress = [0.0]

        def _prog(frac):
            self._rdf_progress[0] = frac

        def _update_status() -> None:
            frac = self._rdf_progress[0]
            if use_traj:
                status.update(f"Computing PXRD... frame {int(frac * n_total)}/{n_total}")
            else:
                status.update(f"Computing PXRD... {int(frac * 100)}%")

        poll = self.set_interval(0.1, _update_status)
        await asyncio.sleep(0)
        try:
            res = await asyncio.to_thread(run_pxrd, self._molecule, wl, tt_max, width, _prog,
                                          self._frames, start, end, use_traj)
        finally:
            poll.stop()
        if res.empty:
            status.update(res.summary)
            return
        self._last_tt = res.x
        self._last_intensity = res.y
        plot = self.query_one("#pxrd-plot", PlotWidget)
        plot.set_data(res.x, res.y, title=res.title, x_label=res.x_label, y_min=0.0,
                      y_max=res.y_max)
        status.update(res.summary)

    def _load_experimental(self, filepath: str) -> None:
        if not filepath.strip():
            return
        try:
            exp_tt, exp_int = read_experimental_pxrd(filepath.strip())
        except Exception as e:
            self.query_one("#pxrd-status", Label).update(f"Error loading: {e}")
            return
        plot = self.query_one("#pxrd-plot", PlotWidget)
        plot.set_overlay(exp_tt, exp_int, color="rgb(255,130,60)")
        self.query_one("#pxrd-status", Label).update(
            f"Experimental data loaded ({len(exp_tt)} points) — "
            "cyan=simulated, orange=experimental"
        )

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
        self._drag_moved = False
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
        self._atom_count = 0  # of the structure last given to set_molecule

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
        # A release can arrive without a press on this widget: the button went down on
        # another widget or before the app was ready, or the terminal was left reporting
        # the mouse by an earlier program. That is not a click here, so it is ignored.
        if not self._dragging:
            return
        was_drag = self._drag_moved
        self._dragging = False
        self._drag_moved = False
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
        n = len(molecule.atoms)
        if n != self._atom_count:
            # Another frame or an edit changed the atom count: indices picked before no
            # longer name the same atoms (or any atom).
            self._picked_atoms = []
            self.highlighted_atoms = {i for i in self.highlighted_atoms if i < n}
            self._atom_count = n
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
            water_atoms = water_atom_indices(mol)
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
            centroid = (self._locked_centroid if self._locked_centroid is not None
                        else self.molecule.center())
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


def _mpmc_wizard_molecule(app):
    """The molecule the MPMC export wizard was started for. The wizard runs across several
    dialogs and a background step, and the active tab can change in between. Steps called
    without a wizard (or before this was recorded) use the active molecule."""
    mol = getattr(app, "_mpmc_state", {}).get("mol")
    return app.molecule if mol is None else mol


def _task_label(coro) -> str:
    """Human-readable name for a coroutine, e.g. '_do_wrap_async' -> 'wrap'."""
    name = getattr(coro, "__name__", None) or "task"
    for prefix in ("_do_", "_calc_", "_compute_", "_run_", "_load_", "_export_", "_"):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    for suffix in ("_async",):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name.replace("_", " ") or "task"


def _task_failure(task: asyncio.Task) -> BaseException | None:
    """Return the exception a finished task raised, or None (also for cancels)."""
    if task.cancelled():
        return None
    try:
        return task.exception()
    except asyncio.CancelledError:
        return None


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
            await self.mount(Label(f" {truncated} ", id=f"tab-{i}", classes=cls, markup=False))
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
    ("write_crystal", "Save as CRYSTAL23 .d12"),
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
    ("pack_sorbates", "Pack Sorbates"),
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
    ("check", "Check Structure"),
    ("symmetry", "Find Symmetry"),
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
                 load_trajectory: bool = False, load_messages: list[str] | None = None):
        super().__init__()
        # Messages produced while reading the startup file (e.g. overlapping
        # atoms deleted by read_file). Shown as a toast once the UI is mounted.
        self._startup_load_messages: list[str] = list(load_messages or [])
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
        # The edit being applied in the background (worker thread plus refresh), if any.
        self._edit_task: asyncio.Task | None = None

    def _edit_running(self) -> bool:
        """True, with a toast, while an earlier edit is still being applied. A second edit
        or an undo started now would race it on the same molecule."""
        task = self._edit_task
        if task is None or task.done():
            return False
        self.notify("The previous edit is still running. Try again when it has finished.",
                    timeout=4, severity="warning")
        return True

    def _start_edit(self, coro) -> asyncio.Task:
        self._edit_task = self._track_task(coro)
        return self._edit_task

    @staticmethod
    def _has_real_cell(mol) -> bool:
        """A unit cell from the file or the user, not the huge stand-in box given to
        trajectories that carry no cell."""
        return mol.pbc is not None and not is_fallback_cell(mol.pbc)

    def notify_load_notes(self, mol) -> None:
        """Warn about anything the file reader could not honour (e.g. a CIF read as its
        asymmetric unit); Check Structure repeats it."""
        for note in load_notes(mol):
            self.notify(note, timeout=20, severity="warning")

    def notify_overlaps_removed(self, messages: list[str]) -> None:
        """Toast a summary of atoms deleted by read_file's overlap check.

        The classic menu prints each deletion; in the TUI stdout is not
        visible, so the user must be told here or the atoms vanish silently.
        """
        n = len(messages)
        if n == 0:
            return
        plural = "s" if n != 1 else ""
        self.notify(
            f"Deleted {n} overlapping atom{plural} (closer than {OVERLAP_CUTOFF} Å)",
            timeout=6, severity="warning",
        )

    def _track_task(self, coro) -> asyncio.Task:
        """Create an asyncio task and track it for cleanup.

        A failure inside the task is surfaced as an error toast and the busy
        progress bar is cleared. Otherwise the exception would only be printed
        to stderr at garbage-collection time ("Task exception was never
        retrieved") and the status bar would stay stuck on "...".
        """
        label = _task_label(coro)
        task = asyncio.create_task(coro)
        self._bg_tasks.add(task)
        task.add_done_callback(self._bg_tasks.discard)
        task.add_done_callback(lambda t, label=label: self._on_tracked_task_done(t, label))
        return task

    def _on_tracked_task_done(self, task: asyncio.Task, label: str) -> None:
        exc = _task_failure(task)
        if exc is None:
            return
        try:
            self._hide_progress()
        except Exception:
            pass
        try:
            self.notify(f"{label} failed: {exc}", timeout=6, severity="error")
        except Exception:
            pass

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
        else:
            # An empty tab (a trajectory still loading): do not leave the previous tab's
            # structure on screen.
            view.molecule = None
            view.pbc = None
            view.highlighted_atoms = set()
            view._picked_atoms = []
            view._atom_count = 0
            panel = self.query_one(GeometryPanel)
            panel._molecule = None
            panel._populated_tabs = set()
            for dt in panel.query(DataTable):
                dt.clear()
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

    @staticmethod
    def _is_trajectory_file(filepath: str) -> bool:
        """True for multi-frame PDB/XYZ files and readable DCD files."""
        return is_trajectory_file(filepath)

    def _open_new_tab(self, filepath: str) -> None:
        if not filepath.strip():
            return
        filepath = filepath.strip()
        try:
            is_traj = self._is_trajectory_file(filepath)

            if is_traj:
                # Load trajectory in background into a new tab
                # Don't switch yet — stay on current tab until loading finishes
                self._save_current_tab_state()
                empty = Molecule(atoms=[])
                new_tab = TabState(empty, filepath)
                self._tabs.append(new_tab)
                self._refresh_tab_bar()
                self.notify(f"Loading {new_tab.name}...", timeout=2)
                self._track_task(self._load_traj_into_tab(filepath, new_tab))
            else:
                load_msgs: list[str] = []
                mol = read_file(filepath, on_message=load_msgs.append)
                self.open_in_new_tab(mol, filepath)
                self.notify(f"Opened {self._tabs[-1].name}", timeout=2)
                self.notify_overlaps_removed(load_msgs)
                self.notify_load_notes(mol)
        except Exception as e:
            self.notify(f"Error: {e}", timeout=3)

    async def _load_traj_into_tab(self, filepath: str, tab) -> None:
        """Load a trajectory file into a specific tab in the background.

        ``tab`` is the TabState itself: tabs may be closed or opened while the file loads,
        so a position in the tab bar would go stale. (An index is still accepted.)
        """
        if isinstance(tab, int):
            tab = self._tabs[tab] if 0 <= tab < len(self._tabs) else None
        self._show_progress("Loading trajectory...", total=100)
        await asyncio.sleep(0)

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
            return read_trajectory(filepath, _prog)

        try:
            frames = await asyncio.to_thread(_read)
        finally:
            timer.stop()

        self._warn_trajectory_caveats(frames)

        self.query_one("#status-label", Label).update("Detecting bonds...")
        self._update_progress(85)
        await asyncio.sleep(0)

        if frames:
            await asyncio.to_thread(frames[0].detect_bonds)

        self._hide_progress()

        tab_idx = next((i for i, t in enumerate(self._tabs) if t is tab), None)
        if tab_idx is None:
            return  # the tab was closed while its file was loading

        # Update the tab state
        tab.frames = frames
        tab.molecule = frames[0] if frames else Molecule(atoms=[])

        # Switch to the newly loaded tab. When loading into the already-active
        # (empty) tab, saving the current state first would clobber the frames
        # we just stored with the view's stale None.
        if tab_idx == self._active_tab:
            self._restore_tab_state(tab_idx)
        elif len(self.screen_stack) > 1:
            # A dialog is open for the structure on screen: its answer must not land on
            # another tab, so the new one waits for a click.
            self.notify(f"Loaded {len(frames)} frames into the {tab.name} tab", timeout=4)
            return
        else:
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
            yield Label("1x", id="traj-speed")
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

        if self._startup_load_messages:
            self.notify_overlaps_removed(self._startup_load_messages)
            self._startup_load_messages = []
        self.notify_load_notes(self.molecule)

    async def _load_trajectory_async(self) -> None:
        tab = self._tabs[self._active_tab]  # the tab this file belongs to
        self._show_progress("Loading trajectory...", total=100)
        await asyncio.sleep(0)


        def _prog(frac: float) -> None:
            self.call_from_thread(self._update_progress, int(frac * 80))
            self.call_from_thread(
                self.query_one("#status-label", Label).update,
                f"Loading trajectory... {int(frac * 100)}%"
            )

        def _read():
            return read_trajectory(self.filepath, _prog)

        frames = await asyncio.to_thread(_read)
        self._warn_trajectory_caveats(frames)

        self.query_one("#status-label", Label).update("Detecting bonds...")
        self._update_progress(85)
        await asyncio.sleep(0)

        # Detect bonds on first frame only (for display), rest on demand
        if frames:
            await asyncio.to_thread(frames[0].detect_bonds)

        self._update_progress(80)

        if tab is not self._tabs[self._active_tab]:
            # Another file was opened in a new tab while this one loaded: the frames belong
            # to the startup tab, not to whatever is on screen now.
            self._hide_progress()
            if any(t is tab for t in self._tabs):
                tab.frames = frames
                tab.frame_idx = 0
                if frames:
                    tab.molecule = frames[0]
                self.notify(f"Loaded {len(frames)} frames into the {tab.name} tab", timeout=4)
            return

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

    def _warn_trajectory_caveats(self, frames: list | None) -> None:
        """Toast anything the trajectory reader had to assume or fake."""
        for note in trajectory_caveats(frames):
            self.notify(note, severity="warning", timeout=8)

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
        self.query_one("#traj-speed", Label).update(f"{self._play_speed:g}x")
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
        self.push_screen(InfoModal("Void Volume", void_volume_report(result_geo, result_he)))

    async def _calc_surface_area(self) -> None:
        mol = self.molecule
        self._show_progress("Calculating surface area (N2 probe)...")
        text = await asyncio.to_thread(surface_area_text, mol, 100000, 1.4)
        self._hide_progress()
        self.push_screen(InfoModal("Surface Area", text))

    async def _calc_pore_size(self) -> None:
        mol = self.molecule
        self._show_progress("Computing pore size distribution...")
        res = await asyncio.to_thread(run_pore_size, mol, 50000, 50, 20.0)
        self._hide_progress()
        if res.empty:
            self.notify(res.summary, timeout=3)
            return

        class PsdScreen(ModalScreen[None]):
            BINDINGS = [Binding("escape", "close", ""), Binding("q", "close", "")]
            DEFAULT_CSS = ("PsdScreen > Vertical { width: 100%; height: 100%; "
                           "background: $surface; padding: 0 1; }")

            def __init__(self, result):
                super().__init__()
                self._res = result

            def compose(self_inner) -> ComposeResult:
                with Vertical():
                    yield PlotWidget(id="psd-plot")
                    yield Label(f"{self_inner._res.summary}  |  press q to close", id="psd-status")

            def on_mount(self_inner) -> None:
                r = self_inner._res
                self_inner.query_one("#psd-plot", PlotWidget).set_data(
                    r.x, r.y, title=r.title, x_label=r.x_label, y_min=0.0, y_max=r.y_max)

            def action_close(self_inner) -> None:
                self_inner.dismiss(None)

        self.push_screen(PsdScreen(res))

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

        stem = Path(self.filepath).stem or "untitled"
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
        stem = Path(self.filepath).stem or "untitled"
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
        if self._edit_running():
            return
        try:
            import numpy as np
            parts = value.strip().split()
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            self._save_undo()
            self.molecule.atoms[idx].x = np.array([x, y, z])
            self._start_edit(self._refresh_molecule(f"Updated atom {idx + 1} position"))
        except (ValueError, IndexError):
            self.notify("Invalid input. Need: x y z", timeout=3)

    def on_operations_panel_run_command(self, event: OperationsPanel.RunCommand) -> None:
        mol = self.molecule
        cmd = event.command
        view = self.query_one(MoleculeView)

        if cmd == "info":
            self.push_screen(InfoModal("System Info", system_info_text(mol, self.filepath)))
            return

        elif cmd == "check":
            self.push_screen(InfoModal("Check Structure", check_text(mol, self._frames)))
            return

        elif cmd == "symmetry":
            self.push_screen(InfoModal("Find Symmetry", symmetry_text(mol)))
            return

        elif cmd == "write_crystal":
            if not self._require_cell("Writing a CRYSTAL23 input"):
                return
            self.push_screen(
                InputModal("Space group for the CRYSTAL23 input:", default="1",
                           hint="1 writes every atom in P1. A number from Find Symmetry "
                                "writes only that group's asymmetric unit."),
                callback=lambda v: self._crystal_step_group(v, mol),
            )
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
            if len(mol.atoms) == 0 or not self._has_real_cell(mol):
                self.notify("Need a crystal structure with PBC for PXRD", timeout=3)
                return
            self.push_screen(PxrdScreen(mol, self._frames))
            return

        elif cmd == "density3d":
            if not self._frames or len(self._frames) < 2:
                self.notify("Need a trajectory with PBC for 3D density", timeout=3)
                return
            if not self._has_real_cell(mol):
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
            if len(mol.atoms) == 0 or not self._has_real_cell(mol):
                self.notify("Need PBC for void volume", timeout=3)
                return
            self.notify("Calculating void volume (100k samples)...", timeout=2)
            self._track_task(self._calc_void_volume())
            return

        elif cmd == "surface_area":
            if not self._has_real_cell(mol):
                self.notify("Need PBC for surface area", timeout=3)
                return
            self._track_task(self._calc_surface_area())
            return

        elif cmd == "pore_size":
            if not self._has_real_cell(mol):
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

        elif cmd == "pack_sorbates":
            if not self._require_cell("Packing sorbates"):
                return
            self.push_screen(SorbateModal(standalone=True), callback=self._pack_step_model)
            return

        elif cmd == "reduce_cell":
            if not self._has_real_cell(mol):
                self.notify("Need PBC for supercell reduction", timeout=3)
                return
            self.push_screen(ReduceCellScreen(mol), callback=self._apply_reduced_cell)
            return

        elif cmd == "open_geom":
            self.action_toggle_geometry()
            return

        elif cmd == "wrap_center":
            if not self._require_cell("Wrapping atoms"):
                return
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
            if not self._require_cell("Wrapping atoms"):
                return
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
            if not self._require_cell("Sorting atoms into molecules"):
                return
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
            if not self._require_cell("Finding lone atoms"):
                return
            found = lone_atom_summary(mol)
            caution = lone_atom_warning(mol)
            self.push_screen(
                ConfirmModal(
                    "Delete Lone Atoms",
                    "Remove all atoms that have no neighbors within\n"
                    "VDW contact distance. Undo restores them.\n\n"
                    + (f"Would delete: {found}" if found else "No lone atoms found.")
                    + (f"\n\n⚠ {caution}" if caution else ""),
                ),
                callback=self._do_del_lone,
            )
            return

        elif cmd == "extend":
            if not self._require_cell("Extending the cell"):
                return
            self.push_screen(
                ExtendAxisModal(mol.pbc),
                callback=self._do_extend,
            )
            return

        elif cmd == "edit_h":
            if not self._require_cell("Editing H distances"):
                return
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

        elif cmd in ("write_pdb", "write_lammps") and len(mol.atoms) == 0:
            self.notify("Nothing to save: no atoms loaded", timeout=3, severity="warning")
            return

        elif cmd == "write_xyz":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / f"{stem}_out.xyz")
            self.push_screen(
                FileSaveModal("Save XYZ as:", default=default),
                callback=lambda v: self._do_write_xyz(v, mol),
            )
            return

        elif cmd == "write_pdb":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / f"{stem}_out.pdb")
            self.push_screen(
                FileSaveModal("Save standard PDB as:", default=default),
                callback=lambda v: self._do_write_pdb(v, mol),
            )
            return

        elif cmd == "write_cif":
            if not self._require_cell("Writing a CIF"):
                return
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / f"{stem}_out.cif")
            self.push_screen(FileSaveModal("Save CIF as:", default=default),
                             callback=lambda v: self._do_write_cif(v, mol))
            return

        elif cmd == "write_poscar":
            if not self._require_cell("Writing a POSCAR"):
                return
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / "POSCAR")
            self.push_screen(FileSaveModal("Save POSCAR as:", default=default),
                             callback=lambda v: self._do_write_poscar(v, mol))
            return

        elif cmd == "write_lammps":
            stem = Path(self.filepath).stem
            default = str(Path(self.filepath).parent / (stem + ".data"))
            self.push_screen(FileSaveModal("Save LAMMPS data as:", default=default),
                             callback=lambda v: self._do_write_lammps(v, mol))
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
            if not self._require_cell("Writing an MPMC PDB"):
                return
            # The wizard mutates charges and force-field parameters on the
            # live molecule, so make the whole export undoable in one step.
            self._save_undo()
            self._mpmc_state = {"write_charges": False, "write_ff": False,
                                "mol": mol, "filepath": self.filepath}
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

    def _require_cell(self, what: str) -> bool:
        """Toast and return False when the active molecule has no unit cell.

        Core operations such as wrap/extend/sort and the PDB/CIF/POSCAR/MPMC
        writers are no-ops without a cell; guarding here keeps the TUI from
        reporting success for something that did nothing.
        """
        if self._has_real_cell(self.molecule):
            return True
        self.notify(
            f"{what} needs a unit cell (PBC). This structure has none — "
            "set one with Edit > Update Cell first.",
            timeout=5, severity="warning",
        )
        return False

    def _discard_last_undo(self, stack=None) -> None:
        """Drop the most recent undo entry (an operation turned out to be a no-op).
        ``stack`` is the undo stack the entry was pushed on (default: the active tab's)."""
        stack = self._undo_stack if stack is None else stack
        if stack:
            stack.pop()

    def _write_text_file(self, filepath: str, writer, what: str) -> bool:
        """Render with ``writer(handle)`` into memory, then write the file.

        Buffering means a writer that raises (e.g. CellRequiredError) never
        leaves a truncated or empty file behind. Returns True on success.
        """
        import io as _io_mod
        buf = _io_mod.StringIO()
        try:
            writer(buf)
        except ValueError as e:
            self.notify(f"{what} failed: {e}", timeout=6, severity="error")
            return False
        try:
            with open(filepath, "w") as f:
                f.write(buf.getvalue())
        except OSError as e:
            self.notify(f"Error: {e}", timeout=4, severity="error")
            return False
        self.notify(f"Wrote {filepath}", timeout=3)
        return True

    async def _refresh_molecule(self, msg: str = "", mol=None) -> None:
        """Redetect bonds and refresh the view and geometry panel.

        ``mol`` is the molecule that was edited (default: the active one). If another tab
        is on screen by the time a step finishes, the view and panel are left alone.
        """
        mol = self.molecule if mol is None else mol
        self._show_progress("Detecting bonds...", total=100)
        await asyncio.sleep(0)

        def _bond_progress(frac: float) -> None:
            self.call_from_thread(self._update_progress, int(frac * 90))

        await asyncio.to_thread(mol.detect_bonds, progress_callback=_bond_progress)
        self._update_progress(80)
        view = self.query_one(MoleculeView)
        if mol is self.molecule:
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
        if mol is self.molecule:
            panel._molecule = mol
            if panel.is_mounted:
                await self._populate_tables_async(panel, mol)
        else:
            for attr in ('_cached_bond_lengths', '_cached_angles', '_cached_dihedrals',
                         '_cached_contacts', '_cached_molecules', '_cached_lone_atoms'):
                if hasattr(mol, attr):
                    delattr(mol, attr)
        self._update_progress(100)
        self._hide_progress()
        if msg:
            self.notify(msg, timeout=2)
        view._invalidate_cache()

    def _do_wrap(self, confirmed: bool, forward: bool) -> None:
        if not confirmed or self._edit_running():
            return
        self._start_edit(self._do_wrap_async(forward))

    async def _do_wrap_async(self, forward: bool) -> None:
        mol = self.molecule
        self._save_undo()
        self._show_progress("Wrapping atoms...")
        await asyncio.sleep(0)
        await asyncio.to_thread(wrap_atoms, mol, forward)
        label = "forward" if forward else "centered"
        await self._refresh_molecule(f"Wrapped atoms ({label})", mol)

    def _do_sort(self, confirmed: bool) -> None:
        if not confirmed or self._edit_running():
            return
        self._start_edit(self._do_sort_async())

    async def _do_sort_async(self) -> None:
        mol = self.molecule
        self._save_undo()
        self._show_progress("Sorting atoms...")
        await asyncio.sleep(0)
        await asyncio.to_thread(sort_system, mol, False)
        await self._refresh_molecule("Sorted atoms", mol)

    def _do_del_lone(self, confirmed: bool) -> None:
        if not confirmed or self._edit_running():
            return
        self._start_edit(self._do_del_lone_async())

    async def _do_del_lone_async(self) -> None:
        mol = self.molecule
        self._save_undo()
        undo = self._undo_stack
        self._show_progress("Finding lone atoms...")
        await asyncio.sleep(0)
        _, deleted = await asyncio.to_thread(remove_lone_atoms, mol)
        if not deleted:
            self._discard_last_undo(undo)
            self._hide_progress()
            self.notify("No lone atoms found", timeout=3)
            return
        n = len(deleted)
        await self._refresh_molecule(f"Deleted {n} lone atom{'s' if n != 1 else ''}", mol)

    def _do_open_file(self, value: str | None) -> None:
        """Callback for the File > Open... picker."""
        if not value or not value.strip():
            return
        filepath = value.strip()
        # If the buffer is empty (no atoms, fresh launch), reuse the active
        # tab instead of opening a new one — avoids a leftover blank tab.
        if not self.molecule.atoms and not self.filepath:
            try:
                if self._is_trajectory_file(filepath):
                    # Same background path as "open in new tab", but into the
                    # empty active tab so no blank tab is left behind.
                    tab = self._tabs[self._active_tab]
                    tab.filepath = filepath
                    self.filepath = filepath
                    self._refresh_tab_bar()
                    self.notify(f"Loading {Path(filepath).name}...", timeout=2)
                    self._track_task(self._load_traj_into_tab(filepath, tab))
                    return
                load_msgs: list[str] = []
                mol = read_file(filepath, on_message=load_msgs.append)
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
                self.notify_overlaps_removed(load_msgs)
                self.notify_load_notes(mol)
            except Exception as e:
                self.notify(f"Error opening file: {e}", timeout=4, severity="error")
            return
        # Otherwise open in a new tab.
        self._open_new_tab(filepath)

    def _do_extend(self, value: tuple) -> None:
        if not value or len(value) != 3:
            return
        na, nb, nc = value
        if na == 0 and nb == 0 and nc == 0:
            return  # no-op
        if self._edit_running():
            return
        self._start_edit(self._do_extend_async(na, nb, nc))

    async def _do_extend_async(self, na: int, nb: int, nc: int) -> None:
        """Add na, nb, nc extra copies along a, b, c (classic-menu counting)."""
        size = f"{na + 1}x{nb + 1}x{nc + 1}"
        target = self.molecule
        self._save_undo()
        self._show_progress(f"Extending to {size}...")
        await asyncio.sleep(0)

        def _extend_all(mol):
            # All three axes inside one thread — main loop only sees the
            # final state, never the intermediate (renderer doesn't render
            # the partial 2x1x1 / 2x2x1 frames).
            if na > 0:
                extend_axis(mol, 0, na)
            if nb > 0:
                extend_axis(mol, 1, nb)
            if nc > 0:
                extend_axis(mol, 2, nc)
            return mol

        await asyncio.to_thread(_extend_all, target)
        await self._refresh_molecule(f"Extended to {size}", target)

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

    # The writers take the molecule the Save dialog was opened for (``mol``); the active
    # tab can change while the dialog is up. Without it they write the active molecule.

    def _do_write_xyz(self, value: str, mol=None) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), lambda p: self._write_xyz_file(p, mol))

    def _write_xyz_file(self, filepath: str, mol=None) -> None:
        mol = self.molecule if mol is None else mol
        try:
            with open(filepath, "w") as f:
                write_xyz(mol, f)
            self.notify(f"Wrote {filepath}", timeout=3)
        except OSError as e:
            self.notify(f"Error: {e}", timeout=3)

    def _do_write_pdb(self, value: str, mol=None) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), lambda p: self._write_pdb_file(p, mol))

    def _write_pdb_file(self, filepath: str, mol=None) -> None:
        mol = self.molecule if mol is None else mol
        self._write_text_file(
            filepath, lambda f: write_standard_pdb(mol, f), "Writing PDB",
        )

    def _do_write_cif(self, value: str, mol=None) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), lambda p: self._write_cif_file(p, mol))

    def _write_cif_file(self, filepath: str, mol=None) -> None:
        mol = self.molecule if mol is None else mol
        self._write_text_file(filepath, lambda f: write_cif(mol, f), "Writing CIF")

    def _do_write_poscar(self, value: str, mol=None) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), lambda p: self._write_poscar_file(p, mol))

    def _write_poscar_file(self, filepath: str, mol=None) -> None:
        mol = self.molecule if mol is None else mol
        self._write_text_file(
            filepath, lambda f: write_poscar(mol, f), "Writing POSCAR",
        )

    def _crystal_step_group(self, value: str | None, mol=None) -> None:
        if value is None:
            return
        mol = self.molecule if mol is None else mol
        try:
            reduction = crystal_reduction_choice(mol, value)
        except ValueError as e:
            self.notify(str(e), timeout=12, severity="warning")
            return
        stem = Path(self.filepath).stem
        default = str(Path(self.filepath).parent / f"{stem}.d12")
        self.push_screen(
            FileSaveModal("Save CRYSTAL23 input as:", default=default),
            callback=lambda path: self._do_write_crystal(path, reduction, mol),
        )

    def _do_write_crystal(self, value: str, reduction, mol=None) -> None:
        if not value or not value.strip():
            return
        mol = self.molecule if mol is None else mol
        self._confirm_overwrite(
            value.strip(),
            lambda path: self._write_text_file(
                path, lambda f: write_crystal_d12(mol, f, reduction),
                "Writing CRYSTAL23 input"),
        )

    def _do_write_lammps(self, value: str, mol=None) -> None:
        if not value.strip():
            return
        self._confirm_overwrite(value.strip(), lambda p: self._write_lammps_file(p, mol))

    def _write_lammps_file(self, filepath: str, mol=None) -> None:
        mol = self.molecule if mol is None else mol
        self._write_text_file(
            filepath, lambda f: write_lammps_data(mol, f), "Writing LAMMPS data",
        )

    def _apply_reduced_cell(self, result: Molecule | None) -> None:
        if result is None or self._edit_running():
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
        if self._edit_running():
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
            self._discard_last_undo()
            self.notify(f"No {old_el} atoms found", timeout=2)
            return
        # Bond radii changed with the element: re-detect bonds and rebuild
        # the geometry tables like every other mutating operation does.
        self._start_edit(
            self._refresh_molecule(f"Substituted {count} {old_el} -> {new_el}")
        )

    def _do_edit_h(self, value: str) -> None:
        if not value.strip():
            return
        try:
            parts = value.strip().split()
            element = parts[0].strip().capitalize()
            distance = float(parts[1])
            if not np.isfinite(distance) or distance <= 0:
                self.notify("The distance must be a positive number, e.g. O 1.0", timeout=3)
                return
            if self._edit_running():
                return
            self._start_edit(self._do_edit_h_async(element, distance))
        except (ValueError, IndexError):
            self.notify("Invalid input. Use: O 1.0, C 1.09, etc.", timeout=3)

    async def _do_edit_h_async(self, element: str, distance: float) -> None:
        mol = self.molecule
        self._save_undo()
        undo = self._undo_stack
        self._show_progress(f"Editing {element}-H distances...")
        await asyncio.sleep(0)
        _, messages = await asyncio.to_thread(
            adjust_h_distances, mol, element, distance,
        )
        if not messages:
            self._discard_last_undo(undo)
            self._hide_progress()
            self.notify(f"No bonded {element}-H pairs found; nothing changed", timeout=4)
            return
        n = len(messages)
        await self._refresh_molecule(
            f"Set {n} {element}-H distance{'s' if n != 1 else ''} to {distance} A", mol
        )

    def _do_qeq_charges(self, confirmed: bool) -> None:
        if not confirmed or self._edit_running():
            return
        self._start_edit(self._do_qeq_async())

    async def _do_qeq_async(self) -> None:
        mol = self.molecule
        self._save_undo()
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
        self.push_screen(InfoModal("QEq Charges", charge_report(mol)))

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

        try:
            positions, densities, pbc = await asyncio.to_thread(
                compute_density_3d, self._frames, element, n_bins=15, progress_callback=_cb,
            )
        finally:
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
        res = run_isotherm(dirpath)
        if res.empty:
            self.notify(res.summary, timeout=3)
            return
        pressures, uptakes = res.x, res.y

        class IsothermPlotScreen(_AnalysisScreen):
            _prefix = "iso"
            DEFAULT_CSS = _ANALYSIS_CSS.format(cls="IsothermPlotScreen") + """
            IsothermPlotScreen #iso-status { height: 1; color: $text-muted; }
            """
            def __init__(self, p, u, path, log_text=""):
                ModalScreen.__init__(self)
                self._molecule = None
                self._frames = None
                self._elements = []
                self._compute_task = None
                self._p, self._u, self._path = p, u, path
                self._log_text = log_text  # MPMC's own averages and unfinished runs
            def compose(self_inner) -> ComposeResult:
                with Vertical():
                    with Horizontal():
                        yield Button("CSV", id="iso-csv")
                        if self_inner._log_text:
                            yield Button("MPMC averages", id="iso-log")
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
                elif event.button.id == "iso-log":
                    self_inner.app.push_screen(
                        InfoModal("MPMC averages", self_inner._log_text))
                elif event.button.id == "iso-csv":
                    self_inner.app.push_screen(
                        FileSaveModal("Save isotherm CSV:", default="isotherm.csv"),
                        callback=lambda fp: self_inner._do_export_csv(
                            fp, "pressure_atm,avg_molecules", self_inner._p, self_inner._u),
                    )

        self.push_screen(IsothermPlotScreen(pressures, uptakes, dirpath, res.text))
        if "has not finished" in res.text:
            self.notify("Some pressure points have not finished; see MPMC averages.",
                        severity="warning", timeout=8)

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
        try:
            data = read_energy_dat(filepath.strip())
        except (OSError, ValueError, IndexError) as e:  # a folder, not text, ragged rows
            self.notify(f"Could not read {filepath}: {e}", timeout=5, severity="error")
            return
        if not data:
            self.notify(f"Could not read {filepath}", timeout=3)
            return
        if not [k for k in data if k != "step"]:
            self.notify(f"Nothing to plot in {filepath}: it has only a step column",
                        timeout=5, severity="warning")
            return
        self.push_screen(EnergyPlotScreen(data, filepath.strip()))

    def _save_undo(self, frame_cells: bool = False, mol=None, undo=None, redo=None) -> None:
        """Snapshot the molecule for Undo. ``frame_cells`` also keeps every trajectory
        frame's cell (for an edit that changes them all). ``mol``, ``undo`` and ``redo``
        name another tab's molecule and stacks (default: the active tab's)."""
        import copy
        snap = copy.deepcopy(self.molecule if mol is None else mol)
        if self._frames and (mol is None or mol is self.molecule):
            snap._undo_frame_idx = self._frame_idx
        if frame_cells and self._frames:
            snap._undo_frame_cells = [copy.deepcopy(f.pbc) for f in self._frames]
        undo = self._undo_stack if undo is None else undo
        undo.append(snap)
        if len(undo) > 10:
            undo.pop(0)
        (self._redo_stack if redo is None else redo).clear()

    def _swap_in_snapshot(self, source: list, target: list) -> None:
        """Make the top of ``source`` the molecule; push the current one on ``target``."""
        import copy
        snap = source.pop()
        idx = snap.__dict__.pop("_undo_frame_idx", None)
        if (idx is not None and self._frames and idx != self._frame_idx
                and 0 <= idx < len(self._frames)):
            # The snapshot is of another frame than the one on screen: go back to it, or
            # the frame on screen would be overwritten with that frame's atoms.
            self._frame_idx = idx
            self.molecule = self._frames[idx]
        current = copy.deepcopy(self.molecule)
        if self._frames:
            current._undo_frame_idx = self._frame_idx
        cells = snap.__dict__.pop("_undo_frame_cells", None)
        if cells is not None and self._frames and len(cells) == len(self._frames):
            current._undo_frame_cells = [copy.deepcopy(f.pbc) for f in self._frames]
            for frame, cell in zip(self._frames, cells):
                frame.pbc = cell
        target.append(current)
        self.molecule = snap
        if self._frames and 0 <= self._frame_idx < len(self._frames):
            # The shown frame is this molecule: without this the edit comes back as soon
            # as the user steps to another frame and returns.
            self._frames[self._frame_idx] = snap

    def _do_undo(self) -> None:
        if not self._undo_stack:
            self.notify("Nothing to undo", timeout=2)
            return
        if self._edit_running():
            return
        self._swap_in_snapshot(self._undo_stack, self._redo_stack)
        view = self.query_one(MoleculeView)
        view.set_molecule(self.molecule)
        self.query_one(GeometryPanel).set_molecule(self.molecule)
        if self._frames:
            self._update_title()
            self._setup_traj_bar()
        view._invalidate_cache()
        self.notify("Undone", timeout=1)

    def _do_redo(self) -> None:
        if not self._redo_stack:
            self.notify("Nothing to redo", timeout=2)
            return
        if self._edit_running():
            return
        self._swap_in_snapshot(self._redo_stack, self._undo_stack)
        view = self.query_one(MoleculeView)
        view.set_molecule(self.molecule)
        self.query_one(GeometryPanel).set_molecule(self.molecule)
        if self._frames:
            self._update_title()
            self._setup_traj_bar()
        view._invalidate_cache()
        self.notify("Redone", timeout=1)

    # Density profile is now launched as DensityScreen

    def _do_update_cell(self, value: str) -> None:
        if value is None or not value.strip():
            return
        text = value.strip().lower()
        if self._edit_running():
            return
        # A trajectory has one cell per frame; the new cell applies to all of them (as in
        # the classic menu), or the analyses over frames would use the old ones.
        others = [f for f in (self._frames or []) if f is not self.molecule]
        # Sentinel words to remove the box entirely (turn off PBC)
        if text in ("none", "off", "remove", "delete", "no", "0"):
            if self.molecule.pbc is None:
                self.notify("This system already has no PBC", timeout=2)
                return
            self._save_undo(frame_cells=True)
            self.molecule.pbc = None
            for frame in others:
                frame.pbc = None
            self._start_edit(self._refresh_molecule("PBC removed"))
            return
        try:
            parts = value.strip().split()
            a, b, c = float(parts[0]), float(parts[1]), float(parts[2])
            alpha, beta, gamma = float(parts[3]), float(parts[4]), float(parts[5])
            problem = cell_problem(a, b, c, alpha, beta, gamma)
            if problem:
                self.notify(f"{problem}. The cell was not changed.", timeout=5,
                            severity="warning")
                return
            self._save_undo(frame_cells=True)
            if self.molecule.pbc is not None:
                self.molecule.pbc.update(a, b, c, alpha, beta, gamma)
            else:
                self.molecule.pbc = PBC(a, b, c, alpha, beta, gamma)
            for frame in others:
                frame.pbc = PBC(a, b, c, alpha, beta, gamma)
            self._start_edit(
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
            # Existing charges stay until a new source actually succeeds, so
            # cancelling the next step doesn't leave the system with all-zero
            # charges. Both the file loader and QEq overwrite every framework atom.
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
            if self._edit_running():
                return
            self._mpmc_state["write_charges"] = True
            self._start_edit(self._mpmc_run_qeq())

    async def _mpmc_run_qeq(self) -> None:
        mol = _mpmc_wizard_molecule(self)
        self._show_progress("Computing QEq charges...", total=100)
        await asyncio.sleep(0)

        def _progress(frac: float) -> None:
            self.call_from_thread(self._update_progress, int(frac * 100))

        try:
            await asyncio.to_thread(apply_qeq_charges, mol, progress_callback=_progress)
            self._hide_progress()
            total_q = sum(a.charge for a in mol.atoms)
            self.notify(f"QEq charges applied (total: {total_q:.4f} e)", timeout=2)
        except Exception as e:
            # No charges were made, so the export stops here rather than writing a file
            # whose charge column was never filled.
            self._hide_progress()
            self.notify(f"QEq failed: {e}. Nothing was written.", timeout=8, severity="error")
            return
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
        mol = _mpmc_wizard_molecule(self)
        try:
            charges = read_charges_file(
                result["path"],
                skip_first=result.get("skip_first", 0),
                skip_last=result.get("skip_last", 0),
            )

            n = sum(1 for a in mol.atoms if a.mpmc_flag != "M")
            problem = assign_framework_charges(mol, charges)
            if problem:
                self.notify(problem, timeout=5)
                return
            if len(charges) != n:
                self.notify(
                    f"Applied {len(charges)} charges recursively "
                    f"({n // len(charges)}x)", timeout=3,
                )
            self.notify(f"Loaded {len(charges)} charges", timeout=2)
        except (OSError, ValueError) as e:
            self.notify(f"Error reading charges: {e}", timeout=5)
            return
        self._mpmc_step_ask_ff()

    def _mpmc_step_ask_ff(self) -> None:
        self.push_screen(
            ForceFieldModal(_mpmc_wizard_molecule(self)),
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
            mol = _mpmc_wizard_molecule(self)
            _, missing = apply_ff_to_system_report(mol.atoms, get_forcefield(ff_idx))
            label = ["OPLS-AA/UFF", "PHAHST"][ff_idx]
            self.notify(f"Applied {label} force field", timeout=2)
            for msg in missing_ff_messages(missing):
                self.notify(msg.strip("! "), timeout=8, severity="warning")
        self._mpmc_step_sorbate()

    # --- Edit > Pack Sorbates: add movable sorbate molecules to the structure ---

    def _pack_step_model(self, choice: str | None) -> None:
        if not choice:
            return  # Cancel
        if choice == SorbateModal.FROM_FILE:
            self.push_screen(
                FileSaveModal("Sorbate template (MPMC .pqr):", default="", button_label="Open"),
                callback=self._pack_step_template,
            )
            return
        self._pack_ask_count(SORBATE_MODELS[choice])

    def _pack_step_template(self, path: str | None) -> None:
        if not path or not path.strip():
            return
        try:
            model = read_sorbate_template(path.strip())
        except (OSError, ValueError) as e:
            self.notify(f"Could not read the template: {e}", severity="error", timeout=6)
            return
        self._pack_ask_count(model)

    def _pack_ask_count(self, model) -> None:
        self.push_screen(
            InputModal(f"How many {model.label} molecules?", default="1",
                       hint="A number, or 'fill' to pack until nothing else fits. "
                            "Positions and orientations are random, with no overlaps "
                            "(minimum image in the real cell)."),
            callback=lambda value: self._pack_step_count(model, value),
        )

    def _pack_step_count(self, model, value: str | None) -> None:
        if not value or not value.strip():
            return  # Cancel / Escape
        try:
            n = parse_sorbate_count(value)
        except ValueError as e:
            self.notify(str(e), severity="warning")
            self._pack_ask_count(model)
            return
        if self._edit_running():
            return
        self._start_edit(self._pack_sorbates_async(model, n))

    async def _pack_sorbates_async(self, model, n: int | None) -> None:
        mol = self.molecule
        undo, redo = self._undo_stack, self._redo_stack
        self._show_progress(f"Packing {model.label}...")
        await asyncio.sleep(0)
        placed = await asyncio.to_thread(insert_sorbates, mol, model, n)
        if not placed:
            self._hide_progress()
            self.notify("No room for even one molecule with these settings.",
                        severity="error", timeout=5)
            return
        self._save_undo(mol=mol, undo=undo, redo=redo)
        mol.atoms.extend(sorbate_atoms(model, placed, next_mpmc_mol_id(mol)))
        full = "" if n is None or len(placed) == n else f" (asked for {n}; the cell is full)"
        await self._refresh_molecule(
            f"Packed {len(placed)} {model.label}{full}. Save as MPMC PDB for MPMC.", mol)

    def _mpmc_step_sorbate(self) -> None:
        self.push_screen(
            SorbateModal(),
            callback=self._mpmc_step_sorbate_result,
        )

    def _mpmc_step_sorbate_result(self, model_name: str | None) -> None:
        if model_name is None:
            return  # Cancel: abort the MPMC export
        self._mpmc_state["sorbate"] = model_name if model_name else None
        self._mpmc_state["sorbate_n"] = 1
        if not model_name:
            self._mpmc_step_filename()
            return
        self.push_screen(
            InputModal("How many molecules?", default="1",
                       hint="A number, or 'fill' to pack until nothing else fits. "
                            "Positions and orientations are random, with no overlaps "
                            "(minimum image in the real cell)."),
            callback=self._mpmc_step_sorbate_count,
        )

    def _mpmc_step_sorbate_count(self, value: str | None) -> None:
        if not value or not value.strip():
            return  # Cancel / Escape (InputModal returns ""): abort the MPMC export
        try:
            self._mpmc_state["sorbate_n"] = parse_sorbate_count(value)
        except ValueError as e:
            self.notify(str(e), severity="warning")
            self._mpmc_step_sorbate_result(self._mpmc_state["sorbate"])
            return
        self._mpmc_step_filename()

    def _mpmc_step_filename(self) -> None:
        filepath = getattr(self, "_mpmc_state", {}).get("filepath", self.filepath)
        stem = Path(filepath).stem
        default = str(Path(filepath).parent / f"{stem}_out_mpmc.pqr")
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
        mol = _mpmc_wizard_molecule(self)
        if sorbate_name:
            n_want = self._mpmc_state.get("sorbate_n", 1)
            placed = insert_sorbates(mol, sorbate_name, n_want)
            if not placed:
                self.notify("No room for even one molecule; nothing written",
                            severity="error")
                return
            if n_want is None or len(placed) < n_want:
                self.notify(f"Placed {len(placed)} molecules (cell full)", timeout=5)
            sorbate_lines = format_inserted_sorbates(
                sorbate_name, placed, start_atom_id=len(mol.atoms) + 1,
                start_mol_id=next_mpmc_mol_id(mol),
            )

        self._write_text_file(
            filepath,
            lambda f: write_mpmc_pdb(
                mol, f,
                write_charges=self._mpmc_state["write_charges"],
                write_params=self._mpmc_state["write_ff"],
                sorbate_lines=sorbate_lines,
            ),
            "Writing MPMC PDB",
        )

# ======================================================================
# Module: menu_tools
# ======================================================================
"""Classic text menu: the analyses, edits, exports and generators the viewers have.

Everything here is reached through menu numbers that the original pdb_wizard never used
(11 and up in the main menus, 10 and up in the geometry menu), so input scripts written
for the original numbering keep working. The numerical work is done by the same shared
modules the terminal viewer and the desktop GUI call (``analyses``, ``inputgen``,
``database``, ...), so all three front ends report the same numbers.
"""



import numpy as np


# Where each viewer feature (the terminal viewer's Operations-panel keys) lives in the
# classic menu. tests/test_classic_parity.py fails if a viewer feature is missing here.
CLASSIC_FEATURES: Dict[str, str] = {
    "info": "main 11 > 1",
    "void_vol": "main 11 > 2",
    "surface_area": "main 11 > 3",
    "pore_size": "main 11 > 4",
    "rdf": "main 11 > 5",
    "coordination": "main 11 > 6",
    "hbonds": "main 11 > 7",
    "msd": "main 11 > 8",
    "rmsd": "main 11 > 9",
    "gyration": "main 11 > 10",
    "density": "main 11 > 11",
    "density3d": "main 11 > 12",
    "pxrd": "main 11 > 13",
    "energy_plot": "main 11 > 14",
    "check": "main 11 > 15",
    "symmetry": "main 13 > 7",
    "write_crystal": "main 13 > 7",
    "input_generator": "main 14",
    "isotherm_plot": "main 15",
    "db_search": "main 16",
    "undo": "main 17",
    "wrap_center": "main 4",
    "wrap_forward": "main 3",
    "sort": "main 12 > 1",
    "extend": "main 2",
    "del_lone": "main 1 > 5",
    "edit_h": "main 1 > 7",
    "update_cell": "main 5",
    "substitute": "main 12 > 2",
    "reduce_cell": "main 12 > 3",
    "pack_sorbates": "main 12 > 4",
    "qeq_charges": "main 12 > 5",
    "open_file": "main 18",
    "write_xyz": "main 6",
    "write_pdb": "main 8",
    "write_cif": "main 13 > 1",
    "write_poscar": "main 13 > 2",
    "write_lammps": "main 13 > 3",
    "write_mpmc": "main 7",
    "export_png": "main 13 > 4",
    "export_gif": "main 13 > 5",
    "quit": "main 9",
}
# Viewer features with no meaning in a text menu.
CLASSIC_NOT_APPLICABLE = frozenset({"frame_to_tab"})

CM_ANALYSES = [
    ("info", "system info"),
    ("void_vol", "void volume"),
    ("surface_area", "surface area"),
    ("pore_size", "pore size distribution"),
    ("rdf", "radial distribution function"),
    ("coordination", "coordination number"),
    ("hbonds", "hydrogen bonds"),
    ("msd", "mean squared displacement"),
    ("rmsd", "RMSD vs time"),
    ("gyration", "radius of gyration"),
    ("density", "density profile"),
    ("density3d", "3D density map"),
    ("pxrd", "powder XRD"),
    ("energy_plot", "energy plot"),
    ("check", "check structure"),
]

_CM_QUIT = ("q", "quit")
CM_UNDO_MAX_ATOMS = 300_000  # atoms (over all frames) above which undo is not recorded


class _CmCancel(Exception):
    """The user typed q at a prompt."""


# ----------------------------------------------------------------------
# Prompts
# ----------------------------------------------------------------------

def _cm_ask(label: str, default: object = None) -> str:
    """One line of input; blank returns the default, q cancels the whole action."""
    hint = f" [{default}]" if default not in (None, "") else ""
    raw = input(f"{label}{hint} > ").strip()
    if raw.lower() in _CM_QUIT:
        raise _CmCancel
    return raw if raw else ("" if default is None else str(default))


def _cm_ask_number(label: str, default: float, cast: Callable = float, low: float = None):
    while True:
        raw = _cm_ask(label, default)
        try:
            value = cast(raw)
            if isinstance(value, float) and not math.isfinite(value):
                raise ValueError  # nan and inf are not numbers anyone means
            if low is not None and value < low:
                raise ValueError
            return value
        except ValueError:
            print(f"!!! Error: enter a number{'' if low is None else f' of at least {low:g}'} !!!")


def _cm_ask_choice(label: str, choices: List[str], default: str) -> str:
    """One of ``choices`` (case-insensitive)."""
    lookup = {c.lower(): c for c in choices}
    while True:
        raw = _cm_ask(f"{label} ({', '.join(choices)})", default)
        if raw.lower() in lookup:
            return lookup[raw.lower()]
        print(f"!!! Error: choose one of {', '.join(choices)} !!!")


def _cm_ask_yes(label: str, default: bool = False) -> bool:
    while True:
        raw = _cm_ask(f"{label} (y/n)", "y" if default else "n").lower()
        if raw in ("yes", "y", "1"):
            return True
        if raw in ("no", "n", "0"):
            return False
        print("!!! Error reading input !!!")


def _cm_ask_option(prompt: str) -> int:
    """A menu number; 0 for anything that is not a number."""
    try:
        return int(input(prompt))
    except ValueError:
        print("!!! Error converting input to int !!!")
        return 0


def _cm_ask_element(label: str, mol: Molecule, default: str, allow_all: bool = False) -> str:
    elements = elements_of(mol)
    while True:
        raw = _cm_ask(f"{label} ({', '.join(elements)}"
                      f"{', or all' if allow_all else ''})", default)
        if allow_all and raw.lower() in ("all", ALL_ATOMS):
            return ALL_ATOMS
        symbol = get_element(raw).symbol
        if symbol in elements:
            return symbol
        print(f"!!! No {raw} atoms in this structure !!!")


def _cm_ask_name(label: str = "Atom name filter (blank = all)") -> str:
    return _cm_ask(label) or ALL_ATOMS


def _cm_ask_range(frames: Optional[list]) -> tuple:
    """(start, end) frames, 1-based inclusive; (1, None) without a trajectory."""
    if not frames or len(frames) < 2:
        return 1, None
    n = len(frames)
    start = _cm_ask_number("First frame", 1, int, 1)
    end = _cm_ask_number("Last frame", n, int, 1)
    return start, end


# ----------------------------------------------------------------------
# Showing results
# ----------------------------------------------------------------------

def cm_ascii_plot(x, y, width: int = 64, height: int = 12) -> str:
    """A small text plot of y against x (finite points only)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) == 0:
        return ""
    x_lo, x_hi = float(x.min()), float(x.max())
    y_lo, y_hi = float(y.min()), float(y.max())
    x_span = (x_hi - x_lo) or 1.0
    y_span = (y_hi - y_lo) or 1.0
    grid = [[" "] * width for _ in range(height)]
    cols = np.minimum(((x - x_lo) / x_span * (width - 1)).round().astype(int), width - 1)
    rows = np.minimum(((y - y_lo) / y_span * (height - 1)).round().astype(int), height - 1)
    for c, r in zip(cols, rows):
        grid[height - 1 - r][c] = "*"
    lines = []
    for i, row in enumerate(grid):
        label = f"{y_hi:>10.4g}" if i == 0 else f"{y_lo:>10.4g}" if i == height - 1 else " " * 10
        lines.append(f"{label} |{''.join(row).rstrip()}")
    lines.append(" " * 11 + "+" + "-" * width)
    lines.append(" " * 12 + f"{x_lo:<.4g}".ljust(width - 10) + f"{x_hi:>10.4g}")
    return "\n".join(lines)


def _cm_show(res: AnalysisResult) -> None:
    """Print a result the way the viewers show it: title, curve, summary, then offer CSV."""
    if res.empty:
        print(f"\n{res.summary}\n")
        return
    print(f"\n{res.title}\n")
    if len(res.x) > 1:
        print(cm_ascii_plot(res.x, res.y))
        print(" " * 12 + res.x_label)
    print(f"\n{res.summary}")
    if res.text:
        print(f"\n{res.text}")
    path = _cm_ask("\nCSV filename (blank to skip)")
    if path:
        try:
            write_csv(path, res.csv_header, *res.csv_arrays())
            print(f"wrote {path}")
        except OSError as e:
            print(f"!!! Error: {e} !!!")


# ----------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------

def _cm_run_analysis(key: str, mol: Molecule, frames: Optional[list], filename: str) -> None:
    problem = unavailable(key, mol, frames)
    if problem:
        print(f"\n{problem}\n")
        return
    traj = bool(frames and len(frames) > 1)
    elements = elements_of(mol) if mol is not None and mol.atoms else []
    first = elements[0] if elements else ""
    if key == "info":
        print("\n" + system_info_text(mol, filename) + "\n")
    elif key == "check":
        print("\n" + check_text(mol, frames) + "\n")
    elif key == "void_vol":
        print("\n" + void_volume_text(mol) + "\n")
    elif key == "surface_area":
        print("\n" + surface_area_text(mol) + "\n")
    elif key == "pore_size":
        _cm_show(run_pore_size(mol))
    elif key == "rdf":
        el1 = _cm_ask_element("First element", mol, first)
        el2 = _cm_ask_element("Second element", mol, elements[1] if len(elements) > 1 else first)
        name1 = _cm_ask_name("First atom name filter (blank = all)")
        name2 = _cm_ask_name("Second atom name filter (blank = all)")
        r_max = _cm_ask_number("Maximum r (A)", 13.0, float, 0.1)
        n_bins = _cm_ask_number("Number of bins", 100, int, 1)
        start, end = _cm_ask_range(frames)
        _cm_show(run_rdf(mol, el1, el2, name1, name2, r_max, n_bins, frames, start, end))
    elif key == "coordination":
        center = _cm_ask_element("Centre element", mol, first)
        neighbor = _cm_ask_element("Neighbour element", mol,
                                   elements[1] if len(elements) > 1 else first)
        r_max = _cm_ask_number("Maximum r (A)", 5.0, float, 0.1)
        use_traj = traj and _cm_ask_yes("Average over the trajectory?")
        start, end = _cm_ask_range(frames) if use_traj else (1, None)
        _cm_show(run_coordination(mol, center, neighbor, r_max, frames, start, end, use_traj))
    elif key == "hbonds":
        cutoff = _cm_ask_number("Donor-acceptor cutoff (A)", 3.5, float, 0.1)
        angle = _cm_ask_number("Minimum D-H...A angle (deg)", 120.0, float, 0.0)
        start, end = _cm_ask_range(frames)
        _cm_show(run_hbonds(mol, cutoff, angle, frames, start, end))
    elif key == "msd":
        element = _cm_ask_element("Element", mol, first)
        name = _cm_ask_name()
        dt = _cm_ask_number("Time between frames (ps)", 1.0, float, 0.0)
        start, end = _cm_ask_range(frames)
        _cm_show(run_msd(frames, element, name, dt, start, end))
    elif key == "rmsd":
        element = _cm_ask_element("Element", mol, "all", allow_all=True)
        name = _cm_ask_name()
        dt = _cm_ask_number("Time between frames (ps)", 1.0, float, 0.0)
        start, end = _cm_ask_range(frames)
        _cm_show(run_rmsd(frames, element, name, dt, start, end))
    elif key == "gyration":
        element = _cm_ask_element("Element", mol, "all", allow_all=True)
        dt = _cm_ask_number("Time between frames (ps)", 1.0, float, 0.0)
        start, end = _cm_ask_range(frames)
        _cm_show(run_gyration(frames, element, dt, start, end))
    elif key == "density":
        element = _cm_ask_element("Element", mol, first)
        name = _cm_ask_name()
        axis = "abc".index(_cm_ask_choice("Axis", ["a", "b", "c"], "c"))
        start, end = _cm_ask_range(frames)
        _cm_show(run_density_profile(mol, element, name, axis, frames, start, end))
    elif key == "density3d":
        element = _cm_ask_element("Element", mol, first)
        n_bins = _cm_ask_number("Bins per axis", 15, int, 1)
        res = run_density_3d(frames, element, n_bins)
        print(f"\n{res.summary}")
        if res.extra:
            path = _cm_ask("\nCSV filename for the grid (blank to skip)")
            if path:
                pos = res.extra["positions"]
                write_csv(path, "x_A,y_A,z_A,density", pos[:, 0], pos[:, 1], pos[:, 2],
                          res.extra["values"])
                print(f"wrote {path}")
    elif key == "pxrd":
        wavelength = _cm_ask_number("Wavelength (A; Cu Ka 1.5406, Mo Ka 0.7107, Co Ka 1.7903)",
                                    1.5406, float, 0.01)
        tt_max = _cm_ask_number("Maximum 2-theta (deg)", 50.0, float, 5.1)
        width = _cm_ask_number("Peak width (deg)", 0.10, float, 0.001)
        use_traj = traj and _cm_ask_yes("Average over the trajectory?", True)
        start, end = _cm_ask_range(frames) if use_traj else (1, None)
        _cm_show(run_pxrd(mol, wavelength, tt_max, width, frames=frames, start=start, end=end,
                          use_traj=use_traj))
    elif key == "energy_plot":
        _cm_energy(filename)


def _cm_energy(filename: str) -> None:
    path = _cm_ask("MPMC energy file", find_energy_file(filename))
    data = read_energy_dat(path) if path else {}
    if not data:
        print(f"\nNo energy data in {path or 'that file'}\n")
        return
    columns, default = energy_columns(data)
    column = _cm_ask_choice("Column", columns, default)
    n = energy_n_steps(data)
    start = _cm_ask_number("First step", 1, int, 1)
    end = _cm_ask_number("Last step", n, int, 1)
    res = run_energy(data, column, start, end)
    if res.empty:
        print(f"\n{res.summary}\n")
        return
    print(f"\n{res.title}\n")
    print(cm_ascii_plot(res.x, res.y))
    print(" " * 12 + res.x_label)
    print(f"\n{res.summary}")
    out = _cm_ask("\nCSV filename (blank to skip)")
    if out:
        write_energy_csv(data, out, res.extra["slice"])
        print(f"wrote {out}")


def cm_menu_analysis(mol: Molecule, frames: Optional[list], filename: str) -> None:
    """The Analysis menu. ``frames`` is the whole trajectory, or None for one structure."""
    listing = "".join(f"{i} = {label}\n" for i, (_, label) in enumerate(CM_ANALYSES, 1))
    back = len(CM_ANALYSES) + 1
    while True:
        option = _cm_ask_option(f"\nWhich analysis?\n\n{listing}{back} = back to main menu\n\n> ")
        if option == back:
            return
        if 1 <= option <= len(CM_ANALYSES):
            try:
                _cm_run_analysis(CM_ANALYSES[option - 1][0], mol, frames, filename)
            except _CmCancel:
                print("\nCancelled")
            except (OSError, ValueError) as e:
                print(f"!!! Error: {e} !!!")
        elif option:
            print("\nInvalid option!")


def cm_menu_isotherm() -> None:
    try:
        directory = _cm_ask("Isotherm folder (the one holding the pressure folders)", ".")
        _cm_show(run_isotherm(directory))
    except _CmCancel:
        print("\nCancelled")
    except (OSError, ValueError) as e:
        print(f"!!! Error: {e} !!!")


# ----------------------------------------------------------------------
# Edit
# ----------------------------------------------------------------------

def cm_prompt_sorbate_model():
    """Pick a built-in sorbate model or an MPMC .pqr template; None if cancelled."""
    names = list(SORBATE_MODELS)
    print("\nSorbate models:")
    for i, name in enumerate(names):
        print(f"{i:>3} = {name:<14} {SORBATE_MODELS[name].description}")
    while True:
        choice = input("\nModel number or name, or a path to an MPMC .pqr template "
                       "(q to cancel)\n\n> ").strip()
        if choice.lower() in _CM_QUIT:
            return None
        try:
            if choice.isdigit() and int(choice) < len(names):
                return SORBATE_MODELS[names[int(choice)]]
            if choice in SORBATE_MODELS:
                return SORBATE_MODELS[choice]
            return read_sorbate_template(choice)
        except (OSError, ValueError) as e:
            print(f"!!! {e} !!!")


def cm_prompt_sorbate_count() -> Optional[int]:
    """How many molecules to place; None means fill the cell."""
    while True:
        raw = input("\nHow many molecules? (a number, or 'fill' to pack until "
                    "nothing else fits)\n\n> ").strip().lower()
        try:
            return parse_sorbate_count(raw)
        except ValueError:
            print("!!! Enter a positive number or 'fill' !!!")


def cm_place_sorbates(mol: Molecule, model, n: Optional[int]) -> list:
    placed = insert_sorbates(mol, model, n,
                             progress_callback=lambda k: print(f"\rplaced {k}", end=""))
    print(f"\nPlaced {len(placed)} {model.label}"
          + (f" (asked for {n}; the cell is full)" if n is not None and len(placed) < n else ""))
    return placed


def _cm_sort(mol: Molecule) -> Molecule:
    if mol.pbc is None:
        print("\nSorting atoms into molecules needs a unit cell")
        return mol
    sort_system(mol, False)
    print("\nSorted atoms")
    return mol


def _cm_substitute(mol: Molecule) -> Molecule:
    parts = _cm_ask("Substitute which element with which? (old new, e.g. Cu Zn)").split()
    if len(parts) != 2:
        print("!!! Enter: old_element new_element (e.g. Cu Zn) !!!")
        return mol
    old, new = parts
    new_el = get_element(new)
    if new_el.symbol == "X":
        print(f"!!! Unknown element: {new} !!!")
        return mol
    count = 0
    for atom in mol.atoms:
        if atom.element.symbol == old:
            atom._element = new_el
            atom.name = new_el.symbol
            count += 1
    if not count:
        print(f"\nNo {old} atoms found")
        return mol
    mol.bonds = []
    mol.detect_bonds()
    print(f"\nReplaced {count} {old} atom{'s' if count != 1 else ''} with {new_el.symbol}")
    return mol


def _cm_reduce(mol: Molecule) -> Molecule:
    if mol.pbc is None:
        print("\nNeed PBC for supercell reduction")
        return mol
    print("\nHow many unit cells does the supercell hold along each axis?")
    na = _cm_ask_number("a", 2, int, 1)
    nb = _cm_ask_number("b", 2, int, 1)
    nc = _cm_ask_number("c", 2, int, 1)
    try:
        framework, sorbates = identify_framework_and_sorbates(mol, na, nb, nc)
    except ValueError as e:  # not an na x nb x nc supercell
        print(f"\n!!! {e} !!!")
        return mol
    print(f"\nFramework: {len(framework)} unique sites ({len(framework) * na * nb * nc} total "
          f"atoms)  |  Sorbate molecules: {len(sorbates)}")
    keep = list(range(len(sorbates)))
    if sorbates:
        print("\n  #  atoms  formula      centre")
        for i, indices in enumerate(sorbates, 1):
            counts: Dict[str, int] = {}
            for idx in indices:
                symbol = mol.atoms[idx].element.symbol
                counts[symbol] = counts.get(symbol, 0) + 1
            centre = np.mean([mol.atoms[idx].x for idx in indices], axis=0)
            formula = " ".join(f"{el}{n}" for el, n in sorted(counts.items()))
            print(f"{i:>3}  {len(indices):>5}  {formula:<12} "
                  f"{centre[0]:.1f}, {centre[1]:.1f}, {centre[2]:.1f}")
        while True:
            raw = _cm_ask("\nSorbate molecules to keep (all, none, or numbers like 1,3,4)",
                          "all").lower()
            try:
                if raw == "all":
                    break
                keep = [] if raw == "none" else sorted(
                    {int(p) - 1 for p in raw.replace(",", " ").split()})
                if any(k < 0 or k >= len(sorbates) for k in keep):
                    raise ValueError
                break
            except ValueError:
                print(f"!!! Enter all, none, or numbers from 1 to {len(sorbates)} !!!")
    try:
        result = reduce_to_primitive(mol, na, nb, nc, keep)
    except ValueError as e:
        print(f"\n!!! {e} !!!")
        return mol
    notes = load_notes(mol)
    if notes:
        result.load_notes = notes  # what the reader had to assume still applies
    print(f"\nReduced to primitive cell: {len(result.atoms)} atoms, "
          f"kept {len(keep)}/{len(sorbates)} sorbate molecules")
    return result


def _cm_pack(mol: Molecule) -> Molecule:
    if mol.pbc is None:
        print("\n!!! Packing sorbates needs a unit cell (update cell dimensions first) !!!")
        return mol
    model = cm_prompt_sorbate_model()
    if model is None:
        return mol
    placed = cm_place_sorbates(mol, model, cm_prompt_sorbate_count())
    if not placed:
        print("No room for even one molecule.")
        return mol
    mol.atoms.extend(sorbate_atoms(model, placed, next_mpmc_mol_id(mol)))
    set_atom_ids(mol.atoms)
    mol.bonds = []
    print("Added them to the structure. Write an MPMC .pdb (main menu 7) to use them in MPMC.")
    return mol


def _cm_qeq(mol: Molecule) -> Molecule:
    if not mol.atoms:
        print("\nNo atoms loaded")
        return mol
    print("\nComputing QEq charges...")
    try:
        apply_qeq_charges(mol)
    except Exception as e:  # too many atoms, singular matrix, ...
        print(f"!!! QEq failed: {e} !!!")
        return mol
    print("\n" + charge_report(mol) + "\n")
    return mol


def cm_menu_edit(mol: Molecule) -> Molecule:
    """The Edit menu. Returns the (possibly new) structure."""
    actions = {1: _cm_sort, 2: _cm_substitute, 3: _cm_reduce, 4: _cm_pack, 5: _cm_qeq}
    while True:
        option = _cm_ask_option(
            "\nWhat would you like to do?\n\n"
            "1 = sort atoms into molecules\n"
            "2 = substitute an element\n"
            "3 = reduce a supercell to one unit cell\n"
            "4 = pack sorbates into the structure\n"
            "5 = generate QEq charges\n"
            "6 = back to main menu\n\n> "
        )
        if option == 6:
            return mol
        if option in actions:
            try:
                mol = actions[option](mol)
            except _CmCancel:
                print("\nCancelled")
        elif option:
            print("\nInvalid option!")


# ----------------------------------------------------------------------
# Undo
# ----------------------------------------------------------------------

def cm_signature(mol: Molecule) -> tuple:
    """Everything an edit can change, for telling whether an action did anything."""
    pbc = mol.pbc
    cell = None if pbc is None else (pbc.a, pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma)
    return (cell, tuple((a.element.symbol, a.name, float(a.x[0]), float(a.x[1]), float(a.x[2]),
                         a.charge, a.alpha, a.epsilon, a.sigma, a.c6, a.c8, a.c10)
                        for a in mol.atoms))


def _cm_snapshot(state):
    """An independent copy of a structure (or a list of them) for undo.

    copy.deepcopy takes about 12 microseconds per atom, which an old script on a large
    trajectory pays at every wrap or cell update; copying the atoms' fields directly is
    some fifteen times faster and copies the same things."""
    if isinstance(state, list):
        return [_cm_snapshot(m) for m in state]
    new = copy.copy(state)
    for key, value in vars(state).items():
        if key != "atoms":
            setattr(new, key, copy.deepcopy(value))  # cell, bonds, load notes, ...
    atoms = []
    for atom in state.atoms:
        if type(atom) is not Atom:
            atoms.append(copy.deepcopy(atom))
            continue
        twin = Atom.__new__(Atom)
        for field in Atom.__slots__:
            setattr(twin, field, getattr(atom, field))
        twin.x = atom.x.copy()
        atoms.append(twin)
    new.atoms = atoms
    return new


class CmHistory:
    """Undo for the classic menu: ``before`` an action, ``after`` it, ``restore`` to go back.

    ``before(state, index)`` records only ``state[index]`` (one frame of a trajectory)."""

    def __init__(self, depth: int = 20):
        self._depth = depth
        self._stack: list = []
        self._pending = None
        self._unrecorded = False  # the last change was too large to record

    def before(self, state, index: Optional[int] = None) -> None:
        if self._pending is not None:  # already recorded for this action
            return
        target = state if index is None else state[index]
        # Copying a very large structure would stall the menu, so those are not recorded.
        # Older steps are dropped too: undoing one of them would also revert this change.
        frames = target if isinstance(target, list) else [target]
        if sum(len(m.atoms) for m in frames) > CM_UNDO_MAX_ATOMS:
            self._stack = []
            self._unrecorded = True
            return
        self._pending = (index, _cm_snapshot(target), self._sig(target))

    def after(self, state) -> None:
        if self._pending is not None:
            index, snapshot, signature = self._pending
            target = state if index is None else state[index]
            if signature != self._sig(target):
                self._stack.append((index, snapshot))
                del self._stack[:-self._depth]
                self._unrecorded = False
        self._pending = None

    def undo(self):
        """The previous whole state, or None when there is nothing to undo."""
        return self._stack.pop()[1] if self._stack else None

    def restore(self, state):
        """``state`` with the last change undone, or None when there is nothing to undo."""
        if not self._stack:
            return None
        index, snapshot = self._stack.pop()
        if index is None:
            return snapshot
        state[index] = snapshot
        return state

    def nothing_to_undo(self) -> str:
        """Why ``restore`` returned None."""
        if self._unrecorded:
            return (f"Undo is not available: changes to more than {CM_UNDO_MAX_ATOMS:,} atoms "
                    f"are not recorded")
        return "Nothing to undo"

    @staticmethod
    def _sig(state) -> tuple:
        if isinstance(state, list):
            return tuple(cm_signature(m) for m in state)
        return cm_signature(state)


# ----------------------------------------------------------------------
# Measuring
# ----------------------------------------------------------------------

def cm_measure(mol: Molecule, ids: List[int]) -> str:
    """Distance (2 atoms), angle (3) or dihedral (4) between 1-based atom numbers,
    using the minimum image when there is a cell."""
    if len(ids) not in (2, 3, 4) or len(set(ids)) != len(ids):
        raise ValueError("Enter 2, 3 or 4 different atom numbers")
    if any(i < 1 or i > len(mol.atoms) for i in ids):
        raise ValueError(f"Atom numbers run from 1 to {len(mol.atoms)}")
    atoms = [mol.atoms[i - 1] for i in ids]
    label = "-".join(f"{a.element.symbol}{i}" for a, i in zip(atoms, ids))

    def vec(a, b) -> np.ndarray:
        d = b.x - a.x
        return mol.pbc.wrap(d) if mol.pbc is not None else d

    if len(ids) == 2:
        return f"{label}   r = {np.linalg.norm(vec(atoms[0], atoms[1])):.4f} A"
    if len(ids) == 3:
        v1, v2 = vec(atoms[1], atoms[0]), vec(atoms[1], atoms[2])
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        return f"{label}   angle = {np.degrees(np.arccos(np.clip(cos, -1, 1))):.2f} deg"
    # unsigned (0-180), as Molecule.get_dihedrals: the dihedrals list and the GUI agree
    b1, b2, b3 = (vec(atoms[0], atoms[1]), vec(atoms[1], atoms[2]), vec(atoms[2], atoms[3]))
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    cos_d = np.dot(n1, n2) / ((np.linalg.norm(n1) + 1e-10) * (np.linalg.norm(n2) + 1e-10))
    angle = np.degrees(np.arccos(np.clip(cos_d, -1.0, 1.0)))
    return f"{label}   dihedral = {angle:.2f} deg"


def cm_menu_measure(mol: Molecule) -> None:
    while True:
        raw = input("\nAtom numbers: 2 for a distance, 3 for an angle, 4 for a dihedral "
                    "(q to go back)\n\n> ").strip()
        if raw.lower() in _CM_QUIT or not raw:
            return
        try:
            print("\n" + cm_measure(mol, [int(p) for p in raw.replace(",", " ").split()]))
        except ValueError as e:
            message = str(e)
            print(f"!!! {message if 'atom' in message.lower() else 'Enter atom numbers'} !!!")


def cm_dihedrals_list(mol: Molecule) -> List[str]:
    if not mol.bonds:
        mol.detect_bonds()
    atoms = mol.atoms
    return [
        f"{'-'.join(atoms[n].element.symbol for n in (i, j, k, m)):<11} "
        f"{i + 1:>5} {j + 1:>5} {k + 1:>5} {m + 1:>5}   dihedral = {round(angle, 4)}"
        for i, j, k, m, angle in mol.get_dihedrals()
    ]


# ----------------------------------------------------------------------
# Other file formats and pictures
# ----------------------------------------------------------------------

def _cm_render(mol: Molecule, width: int, height: int, turn: float = 0.0) -> np.ndarray:
    """The viewers' default view, turned about the vertical axis by ``turn`` radians."""
    if not mol.bonds:
        mol.detect_bonds()
    rot = rotation_matrix(0, turn, 0) @ rotation_matrix(-0.2, -0.5, 0.0)
    pixels, _ = render_scene(width, height, mol, rot, 4.0, bg_color=(0, 0, 0), pbc=mol.pbc,
                             ssaa=2, ambient=0.31, diffuse=0.72, specular=0.42, shininess=96.0)
    return pixels


def _cm_write_picture(mol: Molecule, path: str, animated: bool) -> None:
    try:
        from PIL import Image
    except ImportError:
        print("!!! Pillow required: pip install Pillow !!!")
        return
    if not animated:
        Image.fromarray(_cm_render(mol, 1600, 1200)).save(path)
        print(f"wrote {path}")
        return
    n = 36
    images = []
    for i in range(n):
        images.append(Image.fromarray(_cm_render(mol, 800, 600, 2 * np.pi * i / n)))
        print(f"\rrendered {i + 1}/{n}", end="")
    images[0].save(path, save_all=True, append_images=images[1:], duration=80, loop=0)
    print(f"\nwrote {path} ({n} frames)")


def _cm_write_crystal(mol: Molecule, stem: str) -> None:
    """Print the symmetry ranking, ask which group to write, write the .d12."""
    if mol.pbc is None:
        print("\n!!! A CRYSTAL23 input needs a unit cell !!!")
        return
    print("\n" + symmetry_text(mol) + "\n")
    try:
        while True:
            answer = _cm_ask("Space group to write (1 = P1 with every atom)", 1)
            try:
                reduction = crystal_reduction_choice(mol, answer)
                break
            except ValueError as e:
                print(f"\n{e}\n")
        path = _cm_ask("\noutput filename", stem + ".d12")
        with open(path, "w") as out:
            write_crystal_d12(mol, out, reduction)
        what = "P1" if reduction is None else (f"space group {reduction.number} "
                                               f"({reduction.symbol}), "
                                               f"{len(reduction.numbers)} atoms")
        print(f"wrote {path} in {what}")
    except _CmCancel:
        print("\nCancelled")
    except (OSError, ValueError) as e:
        print(f"!!! Error: {e} !!!")


def cm_menu_write_other(mol: Molecule, filename: str) -> None:
    stem = os.path.splitext(filename)[0]
    writers = {1: (write_cif, ".cif"), 2: (write_poscar, ".vasp"), 3: (write_lammps_data, ".data")}
    while True:
        option = _cm_ask_option(
            "\nWhat would you like to write?\n\n"
            "1 = .cif\n"
            "2 = VASP POSCAR\n"
            "3 = LAMMPS data file\n"
            "4 = PNG picture\n"
            "5 = rotating GIF\n"
            "6 = back to main menu\n"
            "7 = CRYSTAL23 .d12 (lists the space groups it fits first)\n\n> "
        )
        if option == 6:
            return
        if option == 7:
            _cm_write_crystal(mol, stem)
            continue
        if option not in (1, 2, 3, 4, 5):
            if option:
                print("\nInvalid option!")
            continue
        if not mol.atoms:
            print("\nNo atoms loaded")
            continue
        try:
            if option in writers:
                writer, _ = writers[option]
                path = input("\noutput filename > ").strip()
                with open(path, "w") as out:
                    writer(mol, out)
                print(f"wrote {path}")
            else:
                path = _cm_ask("\noutput filename", stem + (".png" if option == 4 else ".gif"))
                _cm_write_picture(mol, path, animated=option == 5)
        except _CmCancel:
            print("\nCancelled")
        except (OSError, ValueError) as e:
            print(f"!!! Error: {e} !!!")


# ----------------------------------------------------------------------
# Simulation inputs
# ----------------------------------------------------------------------

def _cm_format_default(value: object) -> str:
    if isinstance(value, (list, tuple)):
        return ", ".join(f"{v:g}" for v in value)
    return str(value)


def _cm_check_range(var, value) -> None:
    """Whole-number fields have the limits the viewers' spin boxes enforce."""
    if var.kind != "int":
        return
    lo = int(var.min_val) if var.min_val is not None else 0
    hi = int(var.max_val) if var.max_val is not None else 99
    if not lo <= value <= hi:
        raise ValueError(f"{var.label} must be between {lo} and {hi}")


def cm_menu_inputs(mol: Molecule, filename: str) -> None:
    """Generate MPMC, CP2K or OpenMM inputs from the same forms the viewers show."""
    keys = list(ENGINES)
    listing = "".join(f"{i} = {ENGINES[k].label}\n" for i, k in enumerate(keys, 1))
    back = len(keys) + 1
    option = _cm_ask_option(f"\nWhich inputs?\n\n{listing}{back} = back to main menu\n\n> ")
    if not 1 <= option <= len(keys):
        if option not in (0, back):
            print("\nInvalid option!")
        return
    engine = keys[option - 1]
    spec = ENGINES[engine]
    print(f"\n{spec.description}\nPress Enter to keep a default; q cancels.\n")
    values: Dict[str, object] = {}
    try:
        for var in spec.user_vars:
            while True:
                if var.kind == "select" and var.choices:
                    values[var.name] = _cm_ask_choice(var.label, list(var.choices), var.default)
                    break
                raw = _cm_ask(var.label, _cm_format_default(var.default))
                try:
                    value = parse_form_value(var, raw)
                    _cm_check_range(var, value)
                    values[var.name] = value
                    break
                except ValueError as e:
                    print(f"!!! Error: {e} !!!")
        check = check_engine(mol, engine, values)
        print(f"\n{check.preview}\n\n{check.status}")
        if not check.ok:
            return
        outdir = _cm_ask("\nOutput folder", default_engine_outdir(filename, engine))
    except _CmCancel:
        print("\nCancelled")
        return
    outcome = generate_engine_inputs(mol, engine, values, outdir)
    if outcome.preview:
        print(f"\n{outcome.preview}")
    print(f"\n{outcome.status}\n")


# ----------------------------------------------------------------------
# Databases
# ----------------------------------------------------------------------

def cm_menu_database() -> Optional[str]:
    """Search COD, RCSB, PubChem or the Materials Project and download a structure.

    Returns the downloaded file's path when the user wants to open it, else None."""
    sources = ["all"] + list(BACKENDS)
    try:
        source = _cm_ask_choice("Database", sources, "all")
        modes = list(modes_for(source))
        mode = _cm_ask_choice("Search by", modes, modes[0]) if len(modes) > 1 else modes[0]
        query = _cm_ask(f"Search for (e.g. {placeholder_for(source, mode)})")
        if not query:
            return None
        if mode == "id":
            picked = {"source": source, "id": query}
        else:
            print("\nSearching...")
            results, message = search_backends(source, mode, query)
            print(message)
            if not results:
                return None
            shown = results[:50]
            print()
            for i, r in enumerate(shown, 1):
                print(f"{i:>3}  {BACKENDS[r['source']]['label']:<17} {str(r.get('id') or ''):<12} "
                      f"{str(r.get('formula') or ''):<16} {str(r.get('title') or '')}")
            if len(results) > len(shown):
                print(f"... and {len(results) - len(shown)} more; narrow the search to see them")
            number = _cm_ask_number("\nDownload which one? (number)", 1, int, 1)
            if number > len(shown):
                print("\nInvalid option!")
                return None
            picked = shown[number - 1]
        outdir = _cm_ask("Download into folder", ".")
        print("\nDownloading...")
        path = BACKENDS[picked["source"]]["fetch"](picked["id"], outdir)
        print(f"wrote {path}")
        return path if _cm_ask_yes("Open it now? (replaces the current structure)") else None
    except _CmCancel:
        print("\nCancelled")
    except Exception as e:  # network errors, bad ids, unwritable folder
        print(f"!!! Error: {e} !!!")
    return None

# ======================================================================
# Module: menu
# ======================================================================
"""Classic text menu interface (from original pdb_wizard)."""





class _CmEditScope:
    """Set by the main menus while geometry analysis or the MPMC writer runs, so those can
    record an undo step just before they change the structure (and only then)."""

    about_to_change = None


def _cm_about_to_change() -> None:
    if _CmEditScope.about_to_change is not None:
        _CmEditScope.about_to_change()


def _cm_write_failed(error: OSError) -> None:
    """An output file could not be written: say so and go back to the menu."""
    print(f"!!! Error: could not write the file: {error} !!!")


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
                _cm_about_to_change()
                problem = assign_framework_charges(mol, charges)
                if problem:
                    raise ValueError(problem)
                break
            except (TypeError, ValueError, FileNotFoundError) as e:
                print(f"!!! Error reading charges: {e} !!!")
            except OSError as e:  # a folder, no permission: back to the main menu
                print(f"!!! Error reading charges: {e} !!!")
                return

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
                key = ff_in.strip().upper()
                if key in ff_map:
                    ff_idx = ff_map[key]
                elif key.isdigit():
                    ff_idx = int(key)
                else:
                    raise ValueError
                if ff_idx not in (0, 1):
                    raise ValueError
                _cm_about_to_change()
                _, missing = apply_ff_to_system_report(mol.atoms, get_forcefield(ff_idx))
                for msg in missing_ff_messages(missing):
                    print(msg)
                break
            except ValueError:
                print("!!! Error reading input !!!")

    filename = input("\noutput filename > ")
    try:
        write_mpmc_pdb(mol, filename, write_charges=write_charges, write_params=write_ff)
    except OSError as e:
        _cm_write_failed(e)
        return
    print(f"\nwrote {filename}\n")


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
                "9 = back to main menu\n"
                "10 = list dihedrals\n"
                "11 = measure a distance, angle or dihedral\n\n> "
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
            caution = lone_atom_warning(mol)
            _cm_about_to_change()
            mol, deleted = remove_lone_atoms(mol)
            if deleted and caution:
                print(f"\n!!! Warning: {caution} (main menu 17 undoes this) !!!")
            if not deleted:
                print("\nNo lone atoms found\n")
            else:
                print("\nDeleting lone atoms:\n")
                for atom in deleted:
                    print(f"{atom.element.symbol:>3} {atom.id:>5} {atom.x}")
        elif option == 6:
            _list_coords(mol)
        elif option == 7:
            while True:
                try:
                    el = input(
                        "\nLook for hydrogens bonded with which element? "
                        "(e.g. C, O, N, etc)\n\n> "
                    ).strip()
                    if not el or get_element(el).symbol == "X":
                        raise ValueError
                    el = get_element(el).symbol
                    dist = float(input(
                        f"\nWhat distance (in angstroms) shall {el}-H bonds be set to?\n\n> "
                    ))
                    break
                except ValueError:
                    print("!!! Error finding element or reading distance !!!")
            _cm_about_to_change()
            mol, messages = adjust_h_distances(mol, el, dist)
            if not messages:
                print(f"\nNo {el}-H bonds found\n")
            for msg in messages:
                print(msg)
        elif option == 8:
            _vmd_preview(mol)
        elif option == 9:
            return mol
        elif option == 10:
            print("\nDihedrals:\n")
            for msg in cm_dihedrals_list(mol):
                print(msg)
        elif option == 11:
            cm_menu_measure(mol)
    return mol


_NEW_OPTIONS = (
    "11 = analysis{scope}\n"
    "12 = edit structure{one}\n"
    "13 = write other formats and pictures{one}\n"
    "14 = generate simulation inputs{one}\n"
    "15 = isotherm results\n"
    "16 = fetch a structure from a database\n"
    "17 = undo the last change\n"
    "18 = open another file\n\n> "
)


def _prompt_open_file() -> str | None:
    """A readable structure file to switch to, or None."""
    path = input("\nfile to open (q to cancel) > ").strip()
    if not path or path.lower() in ("q", "quit"):
        return None
    if not os.path.isfile(path):
        print(f"!!! Error: {path} not found !!!")
        return None
    try:
        detect_filetype(path)
    except ValueError as e:
        print(f"!!! Error: {e} !!!")
        return None
    return _cm_preload(path)


# Structures read by "open another file" / the database menu before the current one is
# given up, keyed by path; _run_classic_file takes them from here.
_CM_PRELOADED: dict = {}


def _cm_preload(path: str) -> str | None:
    """Read ``path`` now. Returns it when it loaded (the menu then switches to it), or
    prints why not and returns None so the menu stays on the current structure."""
    try:
        _CM_PRELOADED[path] = _load_classic_file(path)
    except SystemExit as e:  # the readers' own "could not read" exits
        print(f"!!! Error: {e.code if isinstance(e.code, str) else 'could not read ' + path} !!!")
        return None
    except (OSError, ValueError) as e:
        print(f"!!! Error: could not read {path}: {e} !!!")
        return None
    return path


def _main_loop_single(mol: Molecule, filename: str) -> str | None:
    """The single-structure menu. Returns a file to open next, or None to quit."""
    mol = overlap_detector(mol)
    history = CmHistory()

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
                "9 = quit\n"
                "10 = insert sorbates and write MPMC .pdb\n"
                + _NEW_OPTIONS.format(scope="", one="")
            ))
        except ValueError:
            print("!!! Error converting input to int !!!")
            continue
        if option in (2, 3, 4, 5, 12):
            history.before(mol)
        elif option in (1, 7):  # recorded only if (and just before) they change something
            _CmEditScope.about_to_change = lambda: history.before(mol)
        if option == 1:
            try:
                mol = _menu_geom_analysis(mol)
            finally:
                _CmEditScope.about_to_change = None
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
                mol.bonds = []  # bonds depend on the cell; recompute lazily
        elif option == 6:
            out_filename = input("\noutput filename > ")
            try:
                with open(out_filename, "w") as out:
                    write_xyz(mol, out)
                print(f"wrote {out_filename}")
            except OSError as e:
                _cm_write_failed(e)
        elif option == 7:
            try:
                _write_mpmc_options(mol)
            finally:
                _CmEditScope.about_to_change = None
        elif option == 8:
            out_filename = input("\noutput filename > ")
            try:
                with open(out_filename, "w") as out:
                    write_standard_pdb(mol, out)
                print(f"wrote {out_filename}")
            except OSError as e:
                _cm_write_failed(e)
        elif option == 9:
            return None
        elif option == 10:
            _menu_insert_sorbates(mol)
        elif option == 11:
            cm_menu_analysis(mol, None, filename)
        elif option == 12:
            mol = cm_menu_edit(mol)
        elif option == 13:
            cm_menu_write_other(mol, filename)
        elif option == 14:
            cm_menu_inputs(mol, filename)
        elif option == 15:
            cm_menu_isotherm()
        elif option == 16:
            path = cm_menu_database()
            if path and _cm_preload(path):
                return path
        elif option == 17:
            previous = history.restore(mol)
            if previous is None:
                print("\n" + history.nothing_to_undo())
            else:
                mol = previous
                print("\nUndid the last change")
        elif option == 18:
            path = _prompt_open_file()
            if path:
                return path
        else:
            print("\nInvalid option!")
        if option in (1, 2, 3, 4, 5, 7, 12):
            history.after(mol)


def _menu_insert_sorbates(mol: Molecule) -> None:
    """Classic-menu flow: pick a model, pack N copies (or fill), write MPMC .pdb."""
    if mol.pbc is None:
        print("\n!!! Inserting sorbates needs a unit cell (update cell dimensions first) !!!")
        return
    model = cm_prompt_sorbate_model()
    if model is None:
        return
    n = cm_prompt_sorbate_count()
    placed = cm_place_sorbates(mol, model, n)
    if not placed:
        return
    lines = format_inserted_sorbates(model, placed, start_atom_id=len(mol.atoms) + 1,
                                     start_mol_id=next_mpmc_mol_id(mol))
    filename = input("\noutput filename > ")
    try:
        write_mpmc_pdb(mol, filename, sorbate_lines=lines)
    except OSError as e:
        _cm_write_failed(e)
        return
    print(f"\nwrote {filename}\n")


def _main_loop_movie(molecules: list[Molecule], filename: str) -> str | None:
    """The trajectory menu. Returns a file to open next, or None to quit."""
    history = CmHistory(depth=3)
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
                "9 = quit\n"
                "10 = insert sorbates and write MPMC .pdb (single frame)\n"
                + _NEW_OPTIONS.format(scope=" (whole trajectory)", one=" (single frame)")
            ))
        except ValueError:
            print("!!! Error converting input to int !!!")
            continue
        if option in (2, 3, 4, 5):
            history.before(molecules)
        if option == 1:
            idx = _prompt_frame(n)
            if idx is not None:
                # only this frame can change, and only if the user edits it
                _CmEditScope.about_to_change = lambda: history.before(molecules, idx)
                try:
                    molecules[idx] = _menu_geom_analysis(molecules[idx])
                finally:
                    _CmEditScope.about_to_change = None
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
                        m.bonds = []
                print(f"\nUpdated cell dimensions in all {n} frames")
        elif option == 6:
            out_filename = input("\noutput filename > ")
            try:
                with open(out_filename, "w") as out:
                    for m in molecules:
                        write_xyz(m, out)
                print(f"wrote {n} frames to {out_filename}")
            except OSError as e:
                _cm_write_failed(e)
        elif option == 7:
            idx = _prompt_frame(n)
            if idx is not None:
                _CmEditScope.about_to_change = lambda: history.before(molecules, idx)
                try:
                    _write_mpmc_options(molecules[idx])
                finally:
                    _CmEditScope.about_to_change = None
        elif option == 8:
            out_filename = input("\noutput filename > ")
            try:
                with open(out_filename, "w") as out:
                    for m in molecules:
                        write_standard_pdb(m, out)
                print(f"wrote {n} frames to {out_filename}")
            except OSError as e:
                _cm_write_failed(e)
        elif option == 9:
            return None
        elif option == 10:
            idx = _prompt_frame(n)
            if idx is not None:
                _menu_insert_sorbates(molecules[idx])
        elif option == 11:
            cm_menu_analysis(molecules[0], molecules, filename)
        elif option == 12:
            idx = _prompt_frame(n)
            if idx is not None:
                history.before(molecules, idx)
                molecules[idx] = cm_menu_edit(molecules[idx])
        elif option == 13:
            idx = _prompt_frame(n)
            if idx is not None:
                cm_menu_write_other(molecules[idx], filename)
        elif option == 14:
            idx = _prompt_frame(n)
            if idx is not None:
                cm_menu_inputs(molecules[idx], filename)
        elif option == 15:
            cm_menu_isotherm()
        elif option == 16:
            path = cm_menu_database()
            if path and _cm_preload(path):
                return path
        elif option == 17:
            previous = history.restore(molecules)
            if previous is None:
                print("\n" + history.nothing_to_undo())
            else:
                molecules = previous
                print("\nUndid the last change")
        elif option == 18:
            path = _prompt_open_file()
            if path:
                return path
        else:
            print("\nInvalid option!")
        if option in (1, 2, 3, 4, 5, 7, 12):
            history.after(molecules)


def run_classic(filename: str) -> None:
    """Run the menu on ``filename``, then on any file the user opens from it."""
    next_file: str | None = filename
    try:
        while next_file:
            next_file = _run_classic_file(next_file)
    except EOFError:
        # stdin closed mid-prompt (a script without its final "9"): no traceback, still code 1
        print("\nInput ended before the menu was quit; exiting.", file=sys.stderr)
        sys.exit(1)


def _run_classic_file(filename: str) -> str | None:
    if filename in _CM_PRELOADED:
        kind, data = _CM_PRELOADED.pop(filename)
    else:
        try:
            kind, data = _load_classic_file(filename)
        except (OSError, UnicodeDecodeError) as e:
            sys.exit(f"Could not read {filename}: {e}")
    if kind == "movie":
        return _main_loop_movie(data, filename)
    return _main_loop_single(data, filename)


def _load_classic_file(filename: str) -> tuple:
    """Read ``filename`` (asking for a cell if it has none). Returns ("movie", frames) or
    ("single", molecule)."""
    ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""

    is_trajectory = False
    other_format = False
    file_notes: list[str] = []
    if ext == "xyz":
        is_trajectory = check_xyz_trajectory(filename)
    elif ext in ("pdb", "ent", "pqr"):
        is_trajectory = check_pdb_trajectory(filename)
    else:
        # Formats the original menu could not read (.cif, POSCAR, LAMMPS data, .dcd, ...)
        try:
            other_format = detect_filetype(filename) not in ("dcd", "ase_traj")
        except ValueError:
            sys.exit("Unable to determine if input file is xyz or pdb (please rename)")
        is_trajectory = not other_format

    if is_trajectory:
        if ext == "xyz":
            with open(filename) as f:
                molecules, pbcs = read_xyz_trajectory(f)
        elif ext in ("pdb", "ent", "pqr"):
            with open(filename) as f:
                molecules, pbcs = read_pdb_trajectory(f)
        else:
            molecules = read_trajectory(filename)
        return "movie", molecules

    if other_format:
        try:
            loaded = read_file(filename, remove_overlaps=False)
        except (OSError, ValueError) as e:
            sys.exit(f"Could not read {filename}: {e}")
        file_notes = load_notes(loaded)
        for note in file_notes:
            print(f"\n!!! Warning: {note} !!!\n")
        system, pbc = loaded.atoms, loaded.pbc
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
    if file_notes:
        mol.load_notes = file_notes  # so Check Structure repeats them
    mol.detect_bonds()
    return "single", mol

# ======================================================================
# Module: gui_tools
# ======================================================================
"""Desktop GUI tools: plots, background tasks, analysis windows, edit dialogs, the input
generator, database search, the geometry tables and the trajectory bar.

Every computation goes through the same shared modules as the terminal viewer
(``analyses``, ``inputgen``, ``database``, ``trajectory``), so the two front ends offer the
same features and produce the same numbers. Like gui.py, this imports without PySide6.
"""



import numpy as np


try:
    from PySide6 import QtCore, QtGui, QtWidgets
    _HAS_QT_TOOLS = True
except ImportError:  # optional extra, see gui.py
    _HAS_QT_TOOLS = False

if _HAS_QT_TOOLS:
    _QWidgetT = QtWidgets.QWidget
    _QDialogT = QtWidgets.QDialog
else:  # keeps the module importable without Qt
    _QWidgetT = _QDialogT = object

_PLOT_COLORS = [(61, 174, 233), (255, 130, 60)]


def _gui_real_cell(mol: Optional[Molecule]) -> bool:
    """A cell that is really there: not None and not the 10^6 Å stand-in box that readers
    give trajectories without cell data."""
    return mol is not None and mol.pbc is not None and not is_fallback_cell(mol.pbc)


def _mono_font():
    return QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont)


# ======================================================================
# Background work
# ======================================================================

class GuiTask:
    """Run ``fn`` on a worker thread; ``on_done`` / ``on_error`` run on the GUI thread.

    A QTimer polls the future, so no Qt object crosses threads. ``wait()`` blocks until the
    task has finished and its callback has run (for scripts and tests).
    """

    _live: set = set()  # running tasks stay referenced, so a second task never drops the first

    def __init__(self, fn: Callable, on_done: Callable, on_error: Callable) -> None:
        # one daemon thread per task: a slow task (a dead network, a big trajectory) neither
        # holds up the next one nor keeps the process alive after the window is closed
        self._future: concurrent.futures.Future = concurrent.futures.Future()
        threading.Thread(target=self._work, args=(fn,), daemon=True).start()
        self._on_done, self._on_error = on_done, on_error
        self.finished = False
        GuiTask._live.add(self)
        self._timer = QtCore.QTimer()
        self._timer.setInterval(40)
        self._timer.timeout.connect(self._poll)
        self._timer.start()

    def _work(self, fn: Callable) -> None:
        try:
            result = fn()
        except BaseException as e:  # SystemExit too: _poll hands it to on_error
            self._future.set_exception(e)
        else:
            self._future.set_result(result)

    def detach(self) -> None:
        """Drop the callbacks (their window is closing); the work finishes unseen."""
        self._on_done = self._on_error = lambda *_: None

    def _poll(self) -> None:
        if self.finished or not self._future.done():
            return
        self._timer.stop()
        self.finished = True
        GuiTask._live.discard(self)
        try:
            result = self._future.result()
        except BaseException as e:  # reported to the user by on_error; a worker that
            callback, arg = self._on_error, e  # calls sys.exit() must not end the GUI
        else:
            callback, arg = self._on_done, result
        try:
            callback(arg)
        except RuntimeError as e:
            if "already deleted" not in str(e):  # its window was destroyed meanwhile: fine
                raise

    def wait(self, timeout: float = 600.0) -> None:
        concurrent.futures.wait([self._future], timeout=timeout)
        self._poll()


# ======================================================================
# Plot
# ======================================================================

def _nice_ticks(lo: float, hi: float, n: int = 6) -> List[float]:
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        return [lo]
    raw = (hi - lo) / n
    mag = 10 ** np.floor(np.log10(raw))
    step = next((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= raw), raw)
    start = np.ceil(lo / step) * step
    return [float(v) for v in np.arange(start, hi + step * 1e-3, step)]


class QtPlot(_QWidgetT):
    """Line plot drawn with QPainter, like the terminal viewer's plot.

    The wheel (or +/-) zooms the x axis around the cursor, a double click (or r) resets
    it, and hovering shows the nearest point. One overlay curve can be added (e.g. an
    experimental PXRD pattern).
    """

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setMinimumSize(420, 260)
        self.setMouseTracking(True)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.x = np.array([])
        self.y = np.array([])
        self.overlay: Optional[tuple] = None
        self.title = self.x_label = ""
        self.y_min: Optional[float] = 0.0
        self.y_max: Optional[float] = None
        self.x_range: Optional[tuple] = None  # zoom window, None = everything
        self._hover: Optional[int] = None
        self._frame = QtCore.QRectF()

    def set_data(self, x, y, title: str = "", x_label: str = "",
                 y_min: Optional[float] = 0.0, y_max: Optional[float] = None) -> None:
        self.x, self.y = np.asarray(x, float), np.asarray(y, float)
        self.title, self.x_label = title, x_label
        self.y_min, self.y_max = y_min, y_max
        self.x_range = None
        self.overlay = None
        self._hover = None
        self.update()

    def set_result(self, res: AnalysisResult) -> None:
        self.set_data(res.x, res.y, res.title, res.x_label, res.y_min, res.y_max)

    def set_overlay(self, x, y, label: str = "") -> None:
        self.overlay = (np.asarray(x, float), np.asarray(y, float), label)
        self.update()

    def clear(self) -> None:
        self.set_data([], [])

    # ---- geometry -----------------------------------------------------------------------

    def limits(self) -> tuple:
        xs = [a for a in (self.x, self.overlay[0] if self.overlay else None)
              if a is not None and len(a)]
        if not xs:
            return 0.0, 1.0, 0.0, 1.0
        x0 = min(float(np.nanmin(a)) for a in xs)
        x1 = max(float(np.nanmax(a)) for a in xs)
        if self.x_range is not None:
            x0, x1 = self.x_range
        ys = [a for a in (self.y, self.overlay[1] if self.overlay else None)
              if a is not None and len(a)]
        y0 = self.y_min if self.y_min is not None else min(float(np.nanmin(a)) for a in ys)
        y1 = self.y_max if self.y_max is not None else max(float(np.nanmax(a)) for a in ys)
        if x1 <= x0:
            x0, x1 = x0 - 0.5, x1 + 0.5
        if y1 <= y0:
            y1 = y0 + 1.0
        return x0, x1, y0, y1

    def zoom_x(self, factor: float, center: Optional[float] = None) -> None:
        if not len(self.x):
            return
        x0, x1, _, _ = self.limits()
        c = (x0 + x1) / 2 if center is None else center
        lo, hi = c - (c - x0) / factor, c + (x1 - c) / factor
        xs = [self.x] + ([self.overlay[0]] if self.overlay is not None else [])
        full0 = min(float(np.nanmin(a)) for a in xs if len(a))
        full1 = max(float(np.nanmax(a)) for a in xs if len(a))
        self.x_range = None if hi - lo >= full1 - full0 else (max(lo, full0), min(hi, full1))
        self.update()

    def reset_zoom(self) -> None:
        self.x_range = None
        self.update()

    # ---- painting -----------------------------------------------------------------------

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt API)
        p = QtGui.QPainter(self)
        try:  # an exception with the painter still active crashes Qt on the next paint
            self._paint(p)
        finally:
            p.end()

    def _paint(self, p) -> None:
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        pal = self.palette()
        p.fillRect(self.rect(), pal.color(QtGui.QPalette.Base))
        text = pal.color(QtGui.QPalette.Text)
        fm = p.fontMetrics()
        left, right, top, bottom = 72, 18, 30, 46
        frame = QtCore.QRectF(left, top, max(self.width() - left - right, 10),
                              max(self.height() - top - bottom, 10))
        self._frame = frame
        x0, x1, y0, y1 = self.limits()

        def px(x, y):
            return (frame.left() + (x - x0) / (x1 - x0) * frame.width(),
                    frame.bottom() - (y - y0) / (y1 - y0) * frame.height())

        grid = QtGui.QColor(text)
        grid.setAlpha(38)
        xt, yt = _nice_ticks(x0, x1), _nice_ticks(y0, y1)
        p.setPen(QtGui.QPen(grid, 1))
        for t in xt:
            xx = px(t, y0)[0]
            p.drawLine(QtCore.QPointF(xx, frame.top()), QtCore.QPointF(xx, frame.bottom()))
        for t in yt:
            yy = px(x0, t)[1]
            p.drawLine(QtCore.QPointF(frame.left(), yy), QtCore.QPointF(frame.right(), yy))
        p.setPen(QtGui.QPen(text, 1))
        p.drawRect(frame)
        for t in xt:
            label = f"{t:g}"
            p.drawText(QtCore.QPointF(px(t, y0)[0] - fm.horizontalAdvance(label) / 2,
                                      frame.bottom() + 16), label)
        for t in yt:
            label = f"{t:.4g}"
            p.drawText(QtCore.QPointF(frame.left() - 8 - fm.horizontalAdvance(label),
                                      px(x0, t)[1] + 4), label)
        p.drawText(QtCore.QRectF(0, 6, self.width(), 20), QtCore.Qt.AlignCenter, self.title)
        p.drawText(QtCore.QRectF(left, self.height() - 22, frame.width(), 18),
                   QtCore.Qt.AlignCenter, self.x_label)
        if not len(self.x):
            p.drawText(frame, QtCore.Qt.AlignCenter, "Set the options and press Compute")
            return

        p.setClipRect(frame)
        curves = [(self.x, self.y, QtGui.QColor(*_PLOT_COLORS[0]))]
        if self.overlay is not None:
            curves.append((self.overlay[0], self.overlay[1], QtGui.QColor(*_PLOT_COLORS[1])))
        for xs, ys, color in curves:
            path = QtGui.QPainterPath()
            for i, (xi, yi) in enumerate(zip(xs, ys)):
                a, b = px(xi, yi)
                if i == 0:
                    path.moveTo(a, b)
                else:
                    path.lineTo(a, b)
            p.setPen(QtGui.QPen(color, 2))
            p.setBrush(QtCore.Qt.NoBrush)
            p.drawPath(path)
            if len(xs) == 1:
                p.setBrush(color)
                p.drawEllipse(QtCore.QPointF(*px(xs[0], ys[0])), 4, 4)
        p.setClipping(False)
        if self.overlay is not None and self.overlay[2]:
            p.setPen(QtGui.QColor(*_PLOT_COLORS[1]))
            p.drawText(QtCore.QPointF(frame.right() - 8 - fm.horizontalAdvance(self.overlay[2]),
                                      frame.top() + 16), self.overlay[2])
        if self._hover is not None and self._hover < len(self.x):
            i = self._hover
            hx, _ = px(self.x[i], self.y[i])
            p.setPen(QtGui.QPen(pal.color(QtGui.QPalette.Highlight), 1, QtCore.Qt.DashLine))
            p.drawLine(QtCore.QPointF(hx, frame.top()), QtCore.QPointF(hx, frame.bottom()))
            p.setPen(text)
            p.drawText(QtCore.QPointF(frame.left() + 6, frame.top() + 16),
                       f"x = {self.x[i]:.4g}   y = {self.y[i]:.4g}")

    # ---- input --------------------------------------------------------------------------

    def _x_at(self, pos) -> float:
        x0, x1, _, _ = self.limits()
        return x0 + (pos.x() - self._frame.left()) / max(self._frame.width(), 1) * (x1 - x0)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        pos = event.position()
        if not len(self.x) or not self._frame.contains(pos):
            self._hover = None
        else:
            dist = np.abs(self.x - self._x_at(pos))
            self._hover = int(np.nanargmin(dist)) if np.isfinite(dist).any() else None
        self.update()

    def wheelEvent(self, event) -> None:  # noqa: N802
        factor = 1.25 if event.angleDelta().y() > 0 else 1 / 1.25
        self.zoom_x(factor, self._x_at(event.position()))

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        self.reset_zoom()

    def keyPressEvent(self, event) -> None:  # noqa: N802
        key = event.key()
        if key in (QtCore.Qt.Key_Plus, QtCore.Qt.Key_Equal):
            self.zoom_x(1.25)
        elif key == QtCore.Qt.Key_Minus:
            self.zoom_x(1 / 1.25)
        elif key == QtCore.Qt.Key_R:
            self.reset_zoom()
        else:
            super().keyPressEvent(event)


# ======================================================================
# Analysis windows
# ======================================================================

@dataclass
class Field:
    """One option of an analysis window.

    kind: 'element' (default 'first' or 'last'), 'element_all', 'name:<element field>',
    'int', 'float', 'choice' (``choices`` = ((label, value), ...)) or 'dir'.
    """

    key: str
    label: str
    kind: str
    default: object = None
    lo: float = 0
    hi: float = 0
    step: float = 1
    choices: tuple = ()


@dataclass
class AnalysisSpec:
    title: str
    fields: list
    run: Callable  # (mol, frames, params, progress) -> AnalysisResult
    frame_range: bool = False  # Start/End boxes for trajectories
    traj_default: Optional[bool] = None  # a "Use trajectory" box with this default


def _rng(p: Dict) -> tuple:
    return p.get("start", 1), p.get("end")


ANALYSES: Dict[str, AnalysisSpec] = {
    "rdf": AnalysisSpec(
        "Radial Distribution",
        [Field("el1", "Type 1 element", "element", "first"),
         Field("name1", "Type 1 name", "name:el1", ALL_ATOMS),
         Field("el2", "Type 2 element", "element", "last"),
         Field("name2", "Type 2 name", "name:el2", ALL_ATOMS),
         Field("xmax", "r max (Å)", "int", 13, 2, 50),
         Field("bins", "Bins", "int", 100, 10, 1000, 10),
         Field("ymax", "Y max (display)", "int", 5, 1, 100)],
        lambda mol, fr, p, prog: run_rdf(mol, p["el1"], p["el2"], p["name1"], p["name2"],
                                         float(p["xmax"]), p["bins"], fr, *_rng(p),
                                         p.get("traj", False), prog),
        frame_range=True, traj_default=True),
    "coordination": AnalysisSpec(
        "Coordination Number",
        [Field("center", "Center", "element", "first"),
         Field("neighbor", "Neighbor", "element", "last"),
         Field("rmax", "r max (Å)", "int", 5, 2, 20)],
        lambda mol, fr, p, prog: run_coordination(mol, p["center"], p["neighbor"],
                                                  float(p["rmax"]), fr, *_rng(p),
                                                  p.get("traj", False)),
        frame_range=True, traj_default=False),
    "hbonds": AnalysisSpec(
        "Hydrogen Bonds",
        [Field("cutoff", "D–A cutoff (Å)", "float", 3.5, 2.0, 5.0, 0.1),
         Field("angle", "Min D–H···A angle (°)", "int", 120, 90, 170)],
        lambda mol, fr, p, prog: run_hbonds(mol, p["cutoff"], float(p["angle"]), fr, *_rng(p)),
        frame_range=True),
    "msd": AnalysisSpec(
        "Mean Squared Displacement",
        [Field("element", "Element", "element", "last"),
         Field("name", "Name", "name:element", ALL_ATOMS),
         Field("dt", "dt (ps per frame)", "int", 1, 1, 10000)],
        lambda mol, fr, p, prog: run_msd(fr, p["element"], p["name"], float(p["dt"]), *_rng(p)),
        frame_range=True),
    "rmsd": AnalysisSpec(
        "RMSD vs Time",
        [Field("element", "Element", "element_all", ALL_ATOMS),
         Field("name", "Name", "name:element", ALL_ATOMS),
         Field("dt", "dt (ps per frame)", "int", 1, 1, 10000)],
        lambda mol, fr, p, prog: run_rmsd(fr, p["element"], p["name"], float(p["dt"]),
                                          *_rng(p)),
        frame_range=True),
    "gyration": AnalysisSpec(
        "Radius of Gyration",
        [Field("element", "Element", "element_all", ALL_ATOMS),
         Field("dt", "dt (ps per frame)", "int", 1, 1, 10000)],
        lambda mol, fr, p, prog: run_gyration(fr, p["element"], float(p["dt"]), *_rng(p)),
        frame_range=True),
    "density": AnalysisSpec(
        "Density Profile",
        [Field("element", "Element", "element", "first"),
         Field("name", "Name", "name:element", ALL_ATOMS),
         Field("axis", "Axis", "choice", 2, choices=(("a (x)", 0), ("b (y)", 1), ("c (z)", 2)))],
        lambda mol, fr, p, prog: run_density_profile(mol, p["element"], p["name"], p["axis"],
                                                     fr, *_rng(p), p.get("traj", False)),
        frame_range=True, traj_default=True),
    "pxrd": AnalysisSpec(
        "Powder XRD",
        [Field("wavelength", "Source", "choice", 1.5406, choices=PXRD_SOURCES),
         Field("tt_max", "2θ max (°)", "int", 50, 10, 90),
         Field("width", "Peak width (0.01°)", "int", 10, 1, 50)],
        lambda mol, fr, p, prog: run_pxrd(mol, p["wavelength"], float(p["tt_max"]),
                                          p["width"] / 100.0, prog, fr, *_rng(p),
                                          p.get("traj", False)),
        frame_range=True, traj_default=True),
    "pore_size": AnalysisSpec(
        "Pore Size Distribution", [],
        lambda mol, fr, p, prog: run_pore_size(mol, 50000, 50, 20.0)),
    "isotherm_plot": AnalysisSpec(
        "Isotherm Results", [Field("directory", "Isotherm folder", "dir", "isotherm")],
        lambda mol, fr, p, prog: run_isotherm(p["directory"])),
}


def _couple_range(start, end) -> None:
    """Keep a Start spin box at or below its End box, so a frame range is never empty."""
    start.valueChanged.connect(end.setMinimum)
    end.valueChanged.connect(start.setMaximum)
    end.setMinimum(start.value())
    start.setMaximum(end.value())


class AnalysisDialog(_QDialogT):
    """Options on top, Compute / CSV / Save-plot buttons, the plot and a status line.

    The same windows as the terminal viewer's analysis screens, driven by ``ANALYSES``.
    """

    def __init__(self, key: str, mol: Optional[Molecule], frames: Optional[list] = None,
                 parent=None) -> None:
        super().__init__(parent)
        spec = ANALYSES[key]
        self.key, self.spec, self.mol, self.frames = key, spec, mol, frames
        self.doc = None  # the window's document; MainWindow keeps mol/frames in step with it
        self._overlay: Optional[tuple] = None
        self.result: Optional[AnalysisResult] = None
        self.task: Optional[GuiTask] = None
        self._progress = [0.0]
        self.setWindowTitle(spec.title)
        self.resize(940, 600)
        lay = QtWidgets.QVBoxLayout(self)
        grid = QtWidgets.QGridLayout()
        self.inputs: Dict[str, QtWidgets.QWidget] = {}
        cells = []
        for f in spec.fields:
            w = self._make_input(f)
            self.inputs[f.key] = w
            cells.append((f.label, w))
        n_frames = len(frames) if frames else 0
        self.traj_box = None
        if n_frames > 1 and spec.traj_default is not None:
            self.traj_box = QtWidgets.QCheckBox("Use trajectory")
            self.traj_box.setChecked(spec.traj_default)
            cells.append(("", self.traj_box))
        self.start_spin = self.end_spin = None
        if n_frames > 1 and spec.frame_range:
            self.start_spin, self.end_spin = QtWidgets.QSpinBox(), QtWidgets.QSpinBox()
            for spin, value in ((self.start_spin, 1), (self.end_spin, n_frames)):
                spin.setRange(1, n_frames)
                spin.setValue(value)
            _couple_range(self.start_spin, self.end_spin)
            cells += [("Start frame", self.start_spin), ("End frame", self.end_spin)]
        # Element/name pairs read best two to a row; the frame range starts its own row.
        ncols = 2 if any(f.kind.startswith("name:") for f in spec.fields) else 3
        row = col = 0
        for label, w in cells:
            if w is self.start_spin and col:
                row, col = row + 1, 0
            if label:
                grid.addWidget(QtWidgets.QLabel(label), row, 2 * col)
            grid.addWidget(w, row, 2 * col + 1)
            col += 1
            if col == ncols:
                row, col = row + 1, 0
        for c in range(ncols):
            grid.setColumnStretch(2 * c + 1, 1)
        lay.addLayout(grid)
        buttons = QtWidgets.QHBoxLayout()
        self.compute_button = QtWidgets.QPushButton("Compute")
        self.compute_button.setDefault(True)
        self.compute_button.clicked.connect(self.run)
        self.csv_button = QtWidgets.QPushButton("CSV…")
        self.csv_button.clicked.connect(lambda: self.export_csv())
        self.png_button = QtWidgets.QPushButton("Save plot…")
        self.png_button.clicked.connect(lambda: self.export_png())
        buttons.addWidget(self.compute_button)
        if key == "pxrd":
            self.exp_button = QtWidgets.QPushButton("Load experimental…")
            self.exp_button.clicked.connect(lambda: self.load_experimental())
            buttons.addWidget(self.exp_button)
        buttons.addStretch(1)
        buttons.addWidget(self.csv_button)
        buttons.addWidget(self.png_button)
        lay.addLayout(buttons)
        self.plot = QtPlot()
        lay.addWidget(self.plot, 1)
        self.text = QtWidgets.QPlainTextEdit()
        self.text.setReadOnly(True)
        self.text.setFont(_mono_font())
        self.text.setMaximumHeight(140)
        self.text.hide()
        lay.addWidget(self.text)
        self.status = QtWidgets.QLabel("Set the options and press Compute")
        self.status.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        lay.addWidget(self.status)
        self._progress_timer = QtCore.QTimer(self)
        self._progress_timer.setInterval(100)
        self._progress_timer.timeout.connect(
            lambda: self.status.setText(f"Computing… {int(self._progress[0] * 100)}%"))
        if "ymax" in self.inputs:
            self.inputs["ymax"].valueChanged.connect(lambda _: self._apply_ymax())
        self._set_export_enabled(False)

    # ---- inputs -------------------------------------------------------------------------

    def _make_input(self, f: Field):
        mol = self.mol
        if f.kind in ("element", "element_all") or f.kind.startswith("name:"):
            w = QtWidgets.QComboBox()
            if f.kind.startswith("name:"):
                el_box = self.inputs[f.kind.split(":", 1)[1]]
                w.addItems([ALL_ATOMS] + (names_of(mol, el_box.currentText()) if mol else []))
                el_box.currentTextChanged.connect(lambda el, box=w: self._refill_names(box, el))
            else:
                items = elements_of(mol) if mol is not None else []
                if f.kind == "element_all":
                    items = [ALL_ATOMS] + items
                w.addItems(items or ["X"])
                if f.default == "last":
                    w.setCurrentIndex(w.count() - 1)
                elif f.default not in ("first", None) and f.default in items:
                    w.setCurrentText(f.default)
            return w
        if f.kind == "int":
            w = QtWidgets.QSpinBox()
            w.setRange(int(f.lo), int(f.hi))
            w.setSingleStep(int(f.step))
            w.setValue(int(f.default))
            return w
        if f.kind == "float":
            w = QtWidgets.QDoubleSpinBox()
            w.setRange(f.lo, f.hi)
            w.setSingleStep(f.step)
            w.setDecimals(2)
            w.setValue(float(f.default))
            return w
        if f.kind == "choice":
            w = QtWidgets.QComboBox()
            for label, value in f.choices:
                w.addItem(label, value)
            w.setCurrentIndex(max(w.findData(f.default), 0))
            return w
        w = QtWidgets.QLineEdit(str(f.default or ""))  # 'dir'
        act = w.addAction(self.style().standardIcon(QtWidgets.QStyle.SP_DirOpenIcon),
                          QtWidgets.QLineEdit.TrailingPosition)
        act.triggered.connect(lambda _=False, edit=w: self._browse_dir(edit))
        return w

    def _refill_names(self, box, element: str) -> None:
        keep = box.currentText()
        box.clear()
        box.addItems([ALL_ATOMS] + (names_of(self.mol, element) if self.mol is not None else []))
        if box.findText(keep) >= 0:
            box.setCurrentText(keep)

    def set_structure(self, mol: Optional[Molecule], frames: Optional[list]) -> None:
        """Follow the document after an edit, undo or frame change: new atoms, same options."""
        self.mol, self.frames = mol, frames
        elements = elements_of(mol) if mol is not None else []
        for f in self.spec.fields:
            if f.kind not in ("element", "element_all"):
                continue
            box = self.inputs[f.key]
            keep = box.currentText()
            items = ([ALL_ATOMS] if f.kind == "element_all" else []) + elements
            box.blockSignals(True)
            box.clear()
            box.addItems(items or ["X"])
            if box.findText(keep) >= 0:
                box.setCurrentText(keep)
            elif f.default == "last":
                box.setCurrentIndex(box.count() - 1)
            box.blockSignals(False)
            box.currentTextChanged.emit(box.currentText())  # refresh linked name boxes
        n = len(frames) if frames else 0
        if self.end_spin is not None and n > 1:
            for spin in (self.start_spin, self.end_spin):
                spin.setMaximum(n)

    def _snapshot(self) -> tuple:
        """Copies the worker can read while edits change the live structure. Only the
        frame being shown can be edited, so the other frames are shared, not copied."""
        mol = copy.deepcopy(self.mol) if self.mol is not None else None
        frames = self.frames
        if frames:
            frames = [mol if f is self.mol else f for f in frames]
        return mol, frames

    def _browse_dir(self, edit) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a folder", edit.text())
        if path:
            edit.setText(path)

    def params(self) -> Dict:
        out: Dict = {}
        for f in self.spec.fields:
            w = self.inputs[f.key]
            if isinstance(w, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)):
                out[f.key] = w.value()
            elif isinstance(w, QtWidgets.QComboBox):
                out[f.key] = w.currentData() if f.kind == "choice" else w.currentText()
            else:
                out[f.key] = w.text().strip()
        out["traj"] = bool(self.traj_box is not None and self.traj_box.isChecked())
        if self.start_spin is not None:
            out["start"], out["end"] = self.start_spin.value(), self.end_spin.value()
        return out

    # ---- running ------------------------------------------------------------------------

    def run(self) -> GuiTask:
        params = self.params()
        mol, frames = self._snapshot()
        fn = self.spec.run
        self._progress[0] = 0.0
        self.compute_button.setEnabled(False)
        self.status.setText("Computing…")
        self._progress_timer.start()

        def progress(frac: float) -> None:
            self._progress[0] = frac

        self.task = GuiTask(lambda: fn(mol, frames, params, progress), self._done, self._failed)
        return self.task

    def _done(self, res: AnalysisResult) -> None:
        self._progress_timer.stop()
        self.compute_button.setEnabled(True)
        self.result = res
        self.status.setText(res.summary)
        self.text.setPlainText(res.text)
        self.text.setVisible(bool(res.text))
        if res.empty:
            self.plot.clear()
            self._set_export_enabled(False)
            return
        self.plot.set_result(res)
        if self._overlay is not None:  # keep a loaded experimental pattern across recomputes
            self.plot.set_overlay(*self._overlay)
        self._apply_ymax()
        self._set_export_enabled(True)

    def _failed(self, error: Exception) -> None:
        self._progress_timer.stop()
        self.compute_button.setEnabled(True)
        self.status.setText(f"Error: {error}")

    def closeEvent(self, event) -> None:  # noqa: N802
        if self.task is not None:
            self.task.detach()
        self._progress_timer.stop()
        super().closeEvent(event)

    def _apply_ymax(self) -> None:
        if "ymax" in self.inputs and len(self.plot.x):
            self.plot.y_max = float(self.inputs["ymax"].value())
            self.plot.update()

    def _set_export_enabled(self, on: bool) -> None:
        self.csv_button.setEnabled(on)
        self.png_button.setEnabled(on)

    # ---- export -------------------------------------------------------------------------

    def export_csv(self, path: Optional[str] = None) -> Optional[str]:
        if self.result is None or self.result.empty:
            return None
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save CSV", f"{self.key}.csv", "CSV (*.csv)")
        if not path:
            return None
        try:
            write_csv(path, self.result.csv_header, *self.result.csv_arrays())
        except OSError as e:
            self.status.setText(f"Could not write {path}: {e}")
            return None
        self.status.setText(f"Wrote {path}")
        return path

    def export_png(self, path: Optional[str] = None) -> Optional[str]:
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save plot", f"{self.key}.png", "PNG (*.png)")
        if not path:
            return None
        if not self.plot.grab().save(path):
            self.status.setText(f"Could not write {path}")
            return None
        self.status.setText(f"Wrote {path}")
        return path

    def load_experimental(self, path: Optional[str] = None) -> bool:
        if not path:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Experimental PXRD (2θ, intensity)", "",
                "Patterns (*.xy *.csv *.dat *.txt);;All files (*)")
        if not path:
            return False
        try:
            tt, inten = read_experimental_pxrd(path)
        except Exception as e:  # malformed file
            self.status.setText(f"Error loading: {e}")
            return False
        self._overlay = (tt, inten, "experimental")
        self.plot.set_overlay(*self._overlay)
        self.status.setText(f"Experimental data loaded ({len(tt)} points) — "
                            "blue = simulated, orange = experimental")
        return True


class EnergyDialog(_QDialogT):
    """MPMC energy.dat plot: pick a column and a step range; CSV writes every column."""

    def __init__(self, data: Dict, path: str, parent=None) -> None:
        super().__init__(parent)
        self.data, self.path = data, path
        self.result: Optional[AnalysisResult] = None
        self.setWindowTitle("Energy Plot")
        self.resize(900, 560)
        lay = QtWidgets.QVBoxLayout(self)
        row = QtWidgets.QHBoxLayout()
        cols, default = energy_columns(data)
        self.column = QtWidgets.QComboBox()
        self.column.addItems(cols)
        self.column.setCurrentText(default)
        n = energy_n_steps(data)
        self.start, self.end = QtWidgets.QSpinBox(), QtWidgets.QSpinBox()
        for spin, value in ((self.start, 1), (self.end, n)):
            spin.setRange(1, max(n, 1))
            spin.setValue(value)
        _couple_range(self.start, self.end)
        for label, w in (("Column", self.column), ("Start", self.start), ("End", self.end)):
            row.addWidget(QtWidgets.QLabel(label))
            row.addWidget(w)
        row.addStretch(1)
        self.csv_button = QtWidgets.QPushButton("CSV…")
        self.csv_button.clicked.connect(lambda: self.export_csv())
        row.addWidget(self.csv_button)
        lay.addLayout(row)
        self.plot = QtPlot()
        lay.addWidget(self.plot, 1)
        self.status = QtWidgets.QLabel(f"Loaded {path}")
        lay.addWidget(self.status)
        self.column.currentTextChanged.connect(lambda _: self.replot())
        self.start.valueChanged.connect(lambda _: self.replot())
        self.end.valueChanged.connect(lambda _: self.replot())
        self.replot()

    def replot(self) -> None:
        res = run_energy(self.data, self.column.currentText(), self.start.value(),
                         self.end.value())
        self.result = res
        self.status.setText(res.summary)
        if res.empty:
            self.plot.clear()
        else:
            self.plot.set_result(res)

    def export_csv(self, path: Optional[str] = None) -> Optional[str]:
        if self.result is None or self.result.empty:
            return None
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save CSV", "energy.csv",
                                                            "CSV (*.csv)")
        if not path:
            return None
        sl = self.result.extra["slice"]
        try:
            write_energy_csv(self.data, path, sl)
        except OSError as e:
            self.status.setText(f"Error: {e}")
            return None
        self.status.setText(f"Exported steps {sl.start + 1}-{sl.stop} to {path}")
        return path


class TextReport(_QDialogT):
    """A read-only monospace report (System Info, Void Volume, QEq summary, ...)."""

    def __init__(self, title: str, text: str, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(580, 400)
        lay = QtWidgets.QVBoxLayout(self)
        self.text = QtWidgets.QPlainTextEdit(text)
        self.text.setReadOnly(True)
        self.text.setFont(_mono_font())
        lay.addWidget(self.text)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Close)
        copy = buttons.addButton("Copy", QtWidgets.QDialogButtonBox.ActionRole)
        copy.clicked.connect(lambda: QtWidgets.QApplication.clipboard().setText(text))
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)


# ======================================================================
# Edit dialogs
# ======================================================================

class CellDialog(_QDialogT):
    """Update the cell (a b c α β γ) or remove it."""

    def __init__(self, pbc, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Update Unit Cell")
        form = QtWidgets.QFormLayout(self)
        values = ((pbc.a, pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma) if pbc is not None
                  else (10.0, 10.0, 10.0, 90.0, 90.0, 90.0))
        labels = ("a (Å)", "b (Å)", "c (Å)", "α (°)", "β (°)", "γ (°)")
        self.spins = []
        self._original = tuple(float(v) for v in values)
        self._edited = [False] * 6
        for k, (label, v) in enumerate(zip(labels, values)):
            spin = QtWidgets.QDoubleSpinBox()
            if k < 3:
                spin.setRange(0.0001, 1.0e6)
            else:
                spin.setRange(0.001, 179.999)
            spin.setDecimals(6)
            spin.setValue(v)
            spin.valueChanged.connect(lambda _v, i=k: self._edited.__setitem__(i, True))
            self.spins.append(spin)
            form.addRow(label, spin)
        self.remove = QtWidgets.QCheckBox("Remove the cell (no periodic boundaries)")
        self.remove.setEnabled(pbc is not None)
        self.remove.toggled.connect(lambda on: [s.setEnabled(not on) for s in self.spins])
        form.addRow(self.remove)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def values(self) -> Optional[tuple]:
        """(a, b, c, alpha, beta, gamma), or None to remove the cell."""
        if self.remove.isChecked():
            return None
        # untouched fields return the cell's own values, not the spin boxes' rounding
        return tuple(s.value() if edited else orig
                     for s, edited, orig in zip(self.spins, self._edited, self._original))


class ReduceCellDialog(_QDialogT):
    """Reduce a supercell: set its size, Analyze, untick sorbates to drop, Reduce."""

    def __init__(self, mol: Molecule, parent=None) -> None:
        super().__init__(parent)
        self.mol = mol
        self.result_mol: Optional[Molecule] = None
        self.framework: list = []
        self.sorbates: list = []
        self.setWindowTitle("Reduce Supercell")
        self.resize(660, 460)
        lay = QtWidgets.QVBoxLayout(self)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Supercell:"))
        self.spins = []
        for k in range(3):
            spin = QtWidgets.QSpinBox()
            spin.setRange(1, 8)
            spin.setValue(2)
            spin.valueChanged.connect(lambda _v: self._sizes_changed())
            self.spins.append(spin)
            row.addWidget(spin)
            if k < 2:
                row.addWidget(QtWidgets.QLabel("×"))
        self.analyze_button = QtWidgets.QPushButton("Analyze")
        self.analyze_button.clicked.connect(self.analyze)
        self.reduce_button = QtWidgets.QPushButton("Reduce")
        self.reduce_button.clicked.connect(self.reduce)
        self.reduce_button.setEnabled(False)
        row.addStretch(1)
        row.addWidget(self.analyze_button)
        row.addWidget(self.reduce_button)
        lay.addLayout(row)
        self.table = QtWidgets.QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Keep", "#", "Atoms", "Elements", "Center (Å)"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        lay.addWidget(self.table, 1)
        self.status = QtWidgets.QLabel("Set the supercell dimensions and click Analyze")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)

    def sizes(self) -> tuple:
        return tuple(s.value() for s in self.spins)

    def _sizes_changed(self) -> None:
        if self.reduce_button.isEnabled():  # the analysis was for another size
            self.reduce_button.setEnabled(False)
            self.status.setText("Supercell size changed: click Analyze again")

    def analyze(self) -> None:
        na, nb, nc = self.sizes()
        if self.mol.pbc is None:
            self.status.setText("Reducing a supercell needs a periodic cell")
            self.reduce_button.setEnabled(False)
            return
        try:
            self.framework, self.sorbates = identify_framework_and_sorbates(self.mol, na, nb, nc)
        except ValueError as e:  # not an na x nb x nc supercell
            self.framework, self.sorbates = [], []
            self.table.setRowCount(0)
            self.status.setText(str(e))
            self.reduce_button.setEnabled(False)
            return
        self.table.setRowCount(len(self.sorbates))
        for i, idx in enumerate(self.sorbates):
            counts: dict = {}
            for j in idx:
                el = self.mol.atoms[j].element.symbol
                counts[el] = counts.get(el, 0) + 1
            center = np.mean([self.mol.atoms[j].x for j in idx], axis=0) if idx else np.zeros(3)
            keep = QtWidgets.QTableWidgetItem()
            keep.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            keep.setCheckState(QtCore.Qt.Checked)
            self.table.setItem(i, 0, keep)
            cells = (str(i + 1), str(len(idx)),
                     " ".join(f"{e}{n}" for e, n in sorted(counts.items())),
                     f"{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}")
            for col, text in enumerate(cells, 1):
                self.table.setItem(i, col, QtWidgets.QTableWidgetItem(text))
        self.status.setText(
            f"Framework: {len(self.framework)} unique sites "
            f"({len(self.framework) * na * nb * nc} total atoms)  |  "
            f"Sorbate molecules: {len(self.sorbates)}  |  Untick sorbates to drop, then Reduce")
        self.reduce_button.setEnabled(True)

    def keep(self) -> List[int]:
        return [i for i in range(self.table.rowCount())
                if self.table.item(i, 0).checkState() == QtCore.Qt.Checked]

    def reduce(self) -> Optional[Molecule]:
        if not self.reduce_button.isEnabled() or (not self.framework and not self.sorbates):
            self.status.setText("Click Analyze first")
            return None
        try:
            self.result_mol = reduce_to_primitive(self.mol, *self.sizes(), self.keep())
        except ValueError as e:
            self.status.setText(str(e))
            self.reduce_button.setEnabled(False)
            return None
        self.accept()
        return self.result_mol


# ======================================================================
# Input generator
# ======================================================================

class InputGeneratorDialog(_QDialogT):
    """MPMC isotherm, CP2K and OpenMM inputs from the shared engine templates."""

    def __init__(self, mol: Molecule, filepath: str = "", parent=None,
                 engine: str = "mpmc_isotherm") -> None:
        super().__init__(parent)
        self.mol, self.filepath = mol, filepath
        self.engine = engine
        self._values_by_engine: Dict[str, Dict] = {}
        self.widgets: Dict[str, QtWidgets.QWidget] = {}
        self.setWindowTitle("Input Generator")
        self.resize(780, 660)
        lay = QtWidgets.QVBoxLayout(self)
        top = QtWidgets.QHBoxLayout()
        self.engine_box = QtWidgets.QComboBox()
        for key, spec in ENGINES.items():
            self.engine_box.addItem(spec.label, key)
        self.engine_box.setCurrentIndex(self.engine_box.findData(engine))
        self.description = QtWidgets.QLabel()
        self.description.setWordWrap(True)
        top.addWidget(QtWidgets.QLabel("Engine"))
        top.addWidget(self.engine_box)
        top.addWidget(self.description, 1)
        lay.addLayout(top)
        out = QtWidgets.QHBoxLayout()
        self.outdir = QtWidgets.QLineEdit()
        browse = QtWidgets.QPushButton("…")
        browse.setFixedWidth(32)
        browse.clicked.connect(self._browse)
        out.addWidget(QtWidgets.QLabel("Output folder"))
        out.addWidget(self.outdir, 1)
        out.addWidget(browse)
        lay.addLayout(out)
        self.form_host = QtWidgets.QWidget()
        self.form = QtWidgets.QFormLayout(self.form_host)
        lay.addWidget(self.form_host)
        buttons = QtWidgets.QHBoxLayout()
        self.validate_button = QtWidgets.QPushButton("Validate")
        self.template_button = QtWidgets.QPushButton("Show template")
        self.generate_button = QtWidgets.QPushButton("Generate")
        self.generate_button.setDefault(True)
        self.validate_button.clicked.connect(self.validate)
        self.template_button.clicked.connect(self.show_template)
        self.generate_button.clicked.connect(self.generate)
        buttons.addWidget(self.validate_button)
        buttons.addWidget(self.template_button)
        buttons.addStretch(1)
        buttons.addWidget(self.generate_button)
        lay.addLayout(buttons)
        self.preview = QtWidgets.QPlainTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setFont(_mono_font())
        lay.addWidget(self.preview, 1)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        lay.addWidget(self.status)
        self.engine_box.currentIndexChanged.connect(
            lambda _: self.set_engine(self.engine_box.currentData()))
        self._load_engine()

    def _browse(self) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Output folder",
                                                          self.outdir.text())
        if path:
            self.outdir.setText(path)

    def set_engine(self, key: str) -> None:
        if key == self.engine:
            return
        self._values_by_engine[self.engine] = self.values(strict=False)
        self.engine = key
        self.engine_box.blockSignals(True)
        self.engine_box.setCurrentIndex(self.engine_box.findData(key))
        self.engine_box.blockSignals(False)
        self._load_engine()

    def _load_engine(self) -> None:
        spec = ENGINES[self.engine]
        self.description.setText(spec.description)
        self.outdir.setText(default_engine_outdir(self.filepath, self.engine))
        while self.form.rowCount():
            self.form.removeRow(0)
        self.widgets.clear()
        saved = self._values_by_engine.get(self.engine, {})
        for var in spec.user_vars:
            value = saved.get(var.name, var.default)
            if var.kind == "select":
                w = QtWidgets.QComboBox()
                for choice in var.choices or ():
                    model = SORBATE_MODELS.get(choice) if var.name == "sorbate" else None
                    w.addItem(f"{choice}  –  {model.description}" if model else choice, choice)
                w.setCurrentIndex(max(w.findData(value), 0))
            elif var.kind == "int":
                w = QtWidgets.QSpinBox()
                w.setRange(int(var.min_val) if var.min_val is not None else 0,
                           int(var.max_val) if var.max_val is not None else 99)
                w.setSingleStep(int(var.step) if var.step is not None else 1)
                w.setValue(int(value))
            else:  # float, str and csv_floats are typed in, as in the terminal viewer
                if isinstance(value, list):
                    value = ", ".join(repr(float(v)) for v in value)  # exact, not rounded
                w = QtWidgets.QLineEdit(str(value))
            if var.help:
                w.setToolTip(var.help)
            self.widgets[var.name] = w
            self.form.addRow(var.label, w)
        if self.engine == "mpmc_isotherm":
            self.validate()
        else:
            self.status.setText("Edit values, then Generate.")
            self.preview.setPlainText(
                f"Engine: {self.engine}\nOutput file: {spec.output_filename}\n\n"
                "Click Show template to inspect the source template, "
                "or Generate to write the rendered files.")

    def values(self, strict: bool = True) -> Dict:
        out: Dict = {}
        for var in ENGINES[self.engine].user_vars:
            w = self.widgets[var.name]
            if isinstance(w, QtWidgets.QSpinBox):
                raw = w.value()
            elif isinstance(w, QtWidgets.QComboBox):
                raw = w.currentData()
            else:
                raw = w.text()
            try:
                out[var.name] = parse_form_value(var, raw)
            except ValueError:
                if strict:
                    raise
        return out

    def set_structure(self, mol: Molecule, frames: Optional[list] = None) -> None:
        """Follow the document after edits (validation and generation use the new atoms)."""
        self.mol = mol

    def _show(self, outcome) -> None:
        if outcome.preview is not None:
            self.preview.setPlainText(outcome.preview)
        self.status.setText(outcome.status)

    def _engine_mol(self) -> Molecule:
        """The structure the engines see: the stand-in cell of a cell-less trajectory is
        handed over as no cell, so they refuse it as they refuse a structure without one."""
        mol = self.mol
        if mol is not None and mol.pbc is not None and not _gui_real_cell(mol):
            mol = copy.copy(mol)
            mol.pbc = None
        return mol

    def validate(self) -> None:
        self._show(check_engine(self._engine_mol(), self.engine, self.values(strict=False)))

    def show_template(self) -> None:
        self._show(engine_template_text(self.engine))

    def generate(self):
        try:
            values = self.values(strict=True)
        except (ValueError, KeyError) as e:
            self.status.setText(f"Invalid input: {e}")
            return None
        outcome = generate_engine_inputs(self._engine_mol(), self.engine, values,
                                         self.outdir.text().strip())
        self._show(outcome)
        return outcome


# ======================================================================
# Database search
# ======================================================================

class DatabaseDialog(_QDialogT):
    """COD, RCSB PDB, PubChem and Materials Project: search, then fetch into a new tab.

    ``on_fetched(path, label, rid)`` is called with the downloaded file.
    """

    def __init__(self, cache_dir: str, on_fetched: Callable, parent=None) -> None:
        super().__init__(parent)
        self.cache_dir, self.on_fetched = cache_dir, on_fetched
        self.results: list = []
        self.task: Optional[GuiTask] = None
        self.setWindowTitle("Fetch from Database")
        self.resize(820, 520)
        lay = QtWidgets.QVBoxLayout(self)
        row1 = QtWidgets.QHBoxLayout()
        self.source = QtWidgets.QComboBox()
        self.source.addItem("All databases", "all")
        for code, spec in BACKENDS.items():
            self.source.addItem(spec["label"], code)
        self.mode = QtWidgets.QComboBox()
        row1.addWidget(QtWidgets.QLabel("Database"))
        row1.addWidget(self.source)
        row1.addWidget(QtWidgets.QLabel("Search by"))
        row1.addWidget(self.mode)
        row1.addStretch(1)
        lay.addLayout(row1)
        row2 = QtWidgets.QHBoxLayout()
        self.query = QtWidgets.QLineEdit()
        self.search_button = QtWidgets.QPushButton("Search")
        self.fetch_id_button = QtWidgets.QPushButton("Fetch by ID")
        self.fetch_button = QtWidgets.QPushButton("Fetch selected")
        row2.addWidget(self.query, 1)
        for b in (self.search_button, self.fetch_id_button, self.fetch_button):
            row2.addWidget(b)
        lay.addLayout(row2)
        self.table = QtWidgets.QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Source", "ID", "Formula", "Title"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        lay.addWidget(self.table, 1)
        self.preview = QtWidgets.QLabel("")
        self.status = QtWidgets.QLabel(
            "Pick a source + mode → type a query → Search. Or 'Fetch by ID' to skip search.")
        lay.addWidget(self.preview)
        lay.addWidget(self.status)
        self.source.currentIndexChanged.connect(lambda _: self._source_changed())
        self.mode.currentIndexChanged.connect(lambda _: self._update_placeholder())
        self.search_button.clicked.connect(self.search)
        self.query.returnPressed.connect(self.search)
        self.fetch_id_button.clicked.connect(self.fetch_by_id)
        self.fetch_button.clicked.connect(self.fetch_selected)
        self.table.doubleClicked.connect(lambda _: self.fetch_selected())
        self.table.currentCellChanged.connect(lambda r, *_: self._preview(r))
        self._source_changed()

    def closeEvent(self, event) -> None:  # noqa: N802
        if self.task is not None:
            self.task.detach()
        super().closeEvent(event)

    def reject(self) -> None:  # Esc does not send a close event
        if self.task is not None:
            self.task.detach()
        super().reject()

    def _source_changed(self) -> None:
        current = self.mode.currentData()
        self.mode.blockSignals(True)
        self.mode.clear()
        for m in modes_for(self.source.currentData()):
            self.mode.addItem(m.capitalize(), m)
        self.mode.setCurrentIndex(max(self.mode.findData(current), 0))
        self.mode.blockSignals(False)
        self._update_placeholder()

    def _update_placeholder(self) -> None:
        self.query.setPlaceholderText(
            placeholder_for(self.source.currentData(), self.mode.currentData()))

    def _preview(self, row: int) -> None:
        if not 0 <= row < len(self.results):
            return
        r = self.results[row]
        label = BACKENDS.get(r.get("source", ""), {}).get("label", r.get("source", ""))
        bits = [f"<b>{label}</b> {r.get('id') or ''}"]
        if r.get("formula"):
            bits.append(f"formula: {r['formula']}")
        if r.get("title"):
            bits.append(f"title: {r['title']}")
        self.preview.setText("  •  ".join(bits))

    def search(self) -> Optional[GuiTask]:
        source, mode = self.source.currentData(), self.mode.currentData()
        query = self.query.text().strip()
        if mode == "id":  # direct-ID mode skips the search
            if not query:
                self.status.setText("Enter an ID then click 'Fetch by ID'")
                return None
            return self._fetch(source, query)
        if not query:
            return None
        labels = ", ".join(BACKENDS[c]["label"] for c in backends_for(source, mode))
        self.status.setText(f"Searching {labels} ({mode})...")
        self.task = GuiTask(lambda: search_backends(source, mode, query), self._show_results,
                            lambda e: self.status.setText(f"Error: {e}"))
        return self.task

    def _show_results(self, out: tuple) -> None:
        self.results, msg = out
        self.table.setRowCount(len(self.results))
        for i, r in enumerate(self.results):
            label = BACKENDS.get(r.get("source", ""), {}).get("label", r.get("source", ""))
            cells = (label, str(r.get("id") or ""), str(r.get("formula") or "")[:25],
                     str(r.get("title") or "")[:60])
            for j, text in enumerate(cells):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(text))
        self.status.setText(msg)
        self.preview.setText("")

    def fetch_selected(self) -> Optional[GuiTask]:
        row = self.table.currentRow()
        if not 0 <= row < len(self.results):
            self.status.setText("Select a result first")
            return None
        r = self.results[row]
        if not r.get("id"):
            self.status.setText("This result has no ID to fetch")
            return None
        return self._fetch(r.get("source", ""), str(r["id"]))

    def fetch_by_id(self) -> Optional[GuiTask]:
        source = self.source.currentData()
        if source == "all":
            self.status.setText("Pick a specific database first (not 'All')")
            return None
        rid = self.query.text().strip()
        if not rid:
            self.status.setText("Enter an ID in the query field")
            return None
        return self._fetch(source, rid)

    def _fetch(self, source: str, rid: str) -> Optional[GuiTask]:
        spec = BACKENDS.get(source)
        if spec is None:
            self.status.setText(f"Unknown source: {source!r}")
            return None
        if source == "all":
            self.status.setText("Pick a specific database first (not 'All')")
            return None
        self.status.setText(f"Fetching {spec['label']} {rid}...")
        cache = self.cache_dir

        def done(path: str) -> None:
            try:
                opened = self.on_fetched(path, spec["label"], rid)
            except Exception as e:  # e.g. a download that is not a readable structure
                self.status.setText(f"Error opening {rid}: {e}")
                return
            if opened is False:  # the callback reported the problem itself
                self.status.setText(f"Could not open {spec['label']} {rid}")
                return
            self.accept()

        self.task = GuiTask(lambda: spec["fetch"](rid, cache), done,
                            lambda e: self.status.setText(f"Error: {e}"))
        return self.task


# ======================================================================
# Geometry tables
# ======================================================================

GEOMETRY_TABS = ["Bonds", "Angles", "Dihedrals", "Contacts", "Mols", "Lone", "Coords"]
_SORTABLE = {"Bonds", "Angles", "Dihedrals"}


def _label(mol: Molecule, i: int) -> str:
    return f"{i + 1}:{mol.atoms[i].element.symbol}"


def geometry_rows(mol: Optional[Molecule], tab: str, ascending: bool = False) -> tuple:
    """(column headers, [(cells, atom indices)]) for one geometry tab.

    The same tables as the terminal viewer's geometry panel.
    """
    if mol is None or not mol.atoms:
        return [], []
    if tab in _SORTABLE and not mol.bonds:
        mol.detect_bonds()
    if tab == "Bonds":
        data = mol.get_bond_lengths()
        if ascending:
            data.sort(key=lambda x: x[2])
        return (["Atom 1", "Atom 2", "Length (Å)"],
                [((_label(mol, i), _label(mol, j), f"{d:.4f}"), (i, j)) for i, j, d in data])
    if tab == "Angles":
        data = mol.get_angles()
        if ascending:
            data.sort(key=lambda x: x[3])
        return (["Atom 1", "Vertex", "Atom 3", "Angle (°)"],
                [((_label(mol, i), _label(mol, j), _label(mol, k), f"{a:.3f}"), (i, j, k))
                 for i, j, k, a in data])
    if tab == "Dihedrals":
        data = mol.get_dihedrals()
        if ascending:
            data.sort(key=lambda x: x[4])
        return (["Atom 1", "Atom 2", "Atom 3", "Atom 4", "Angle (°)"],
                [((_label(mol, i), _label(mol, j), _label(mol, k), _label(mol, m), f"{a:.3f}"),
                  (i, j, k, m)) for i, j, k, m, a in data])
    if tab == "Contacts":
        set_atom_ids(mol.atoms)
        rows = []
        for msg in get_close_contacts(mol):
            parts = msg.split()
            el1, el2 = parts[0].split("-")
            id1, id2 = int(parts[1]), int(parts[2])
            rows.append(((f"{id1}:{el1}", f"{id2}:{el2}", parts[-1]), (id1 - 1, id2 - 1)))
        return ["Atom 1", "Atom 2", "Distance (Å)"], rows
    if tab == "Mols":
        index = {id(a): k for k, a in enumerate(mol.atoms)}
        rows = []
        for mi, sm in enumerate(mol.find_molecules()):
            atoms = sorted(sm.atoms, key=lambda a: a.atomic_number, reverse=True)
            counts: dict = {}
            for a in atoms:
                counts[a.element.symbol] = counts.get(a.element.symbol, 0) + 1
            formula = "".join(f"{el}{n}" if n > 1 else el for el, n in counts.items())
            rows.append(((str(mi + 1), str(len(atoms)), formula),
                         tuple(index[id(a)] for a in atoms if id(a) in index)))
        return ["#", "Atoms", "Formula"], rows
    if tab == "Lone":
        set_atom_ids(mol.atoms)
        index = {id(a): k for k, a in enumerate(mol.atoms)}
        return (["#", "Element", "x", "y", "z"],
                [((str(a.id), a.element.symbol, f"{a.x[0]:.4f}", f"{a.x[1]:.4f}",
                   f"{a.x[2]:.4f}"), (index[id(a)],)) for a in get_lone_atoms(mol)])
    if tab == "Coords":
        return (["#", "El", "x", "y", "z"],
                [((str(k + 1), a.element.symbol, f"{a.x[0]:.4f}", f"{a.x[1]:.4f}",
                   f"{a.x[2]:.4f}"), (k,)) for k, a in enumerate(mol.atoms)])
    raise ValueError(f"unknown geometry tab {tab!r}")


class GuiGeometryPanel(_QWidgetT):
    """Bonds, angles, dihedrals, close contacts, molecules, lone atoms and coordinates.

    Selecting a row highlights its atoms in the view; editing a Coords cell moves the atom.
    A table is filled only when its tab is shown.
    """

    def __init__(self, on_highlight: Callable[[tuple], None],
                 on_edit_coord: Callable[[int, tuple], None], parent=None) -> None:
        super().__init__(parent)
        self.on_highlight, self.on_edit_coord = on_highlight, on_edit_coord
        self.mol: Optional[Molecule] = None
        self._stale = set(GEOMETRY_TABS)
        self._ascending: Dict[str, bool] = {}
        self._indices: Dict[str, list] = {}
        self._filling = False
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.tabs = QtWidgets.QTabWidget()
        self.tables: Dict[str, QtWidgets.QTableWidget] = {}
        for name in GEOMETRY_TABS:
            t = QtWidgets.QTableWidget(0, 0)
            t.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            t.verticalHeader().setVisible(False)
            t.horizontalHeader().setStretchLastSection(True)
            t.currentCellChanged.connect(lambda r, *_, n=name: self._row_selected(n, r))
            t.itemChanged.connect(lambda item, n=name: self._item_changed(n, item))
            self.tables[name] = t
            self.tabs.addTab(t, name)
        self.tabs.currentChanged.connect(lambda _: self.refresh_visible())
        lay.addWidget(self.tabs)
        row = QtWidgets.QHBoxLayout()
        self.sort_button = QtWidgets.QPushButton("Sort by value")
        self.sort_button.setCheckable(True)
        self.sort_button.toggled.connect(self._toggle_sort)
        self.hint = QtWidgets.QLabel("")
        row.addWidget(self.sort_button)
        row.addWidget(self.hint, 1)
        lay.addLayout(row)

    def current_tab(self) -> str:
        return self.tabs.tabText(self.tabs.currentIndex())

    def set_molecule(self, mol: Optional[Molecule]) -> None:
        self.mol = mol
        self._stale = set(GEOMETRY_TABS)
        self.refresh_visible()

    def refresh_visible(self, force: bool = False) -> None:
        name = self.current_tab()
        self.sort_button.setEnabled(name in _SORTABLE)
        self.sort_button.blockSignals(True)
        self.sort_button.setChecked(self._ascending.get(name, False))
        self.sort_button.blockSignals(False)
        self.hint.setText("Double-click x, y or z to move an atom" if name == "Coords" else "")
        if not self.isVisible() or (name not in self._stale and not force):
            return
        headers, rows = geometry_rows(self.mol, name, self._ascending.get(name, False))
        t = self.tables[name]
        self._filling = True
        t.clear()
        t.setColumnCount(len(headers))
        t.setHorizontalHeaderLabels(headers)
        t.setRowCount(len(rows))
        for r, (cells, _) in enumerate(rows):
            for c, text in enumerate(cells):
                item = QtWidgets.QTableWidgetItem(text)
                if name != "Coords" or c < 2:
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
                t.setItem(r, c, item)
        self._indices[name] = [idx for _, idx in rows]
        self._filling = False
        self._stale.discard(name)

    def _toggle_sort(self, on: bool) -> None:
        self._ascending[self.current_tab()] = on
        self.refresh_visible(force=True)

    def _row_selected(self, name: str, row: int) -> None:
        idx = self._indices.get(name, [])
        if not self._filling and 0 <= row < len(idx):
            self.on_highlight(idx[row])

    def _item_changed(self, name: str, item) -> None:
        if self._filling or name != "Coords" or self.mol is None:
            return
        row, t = item.row(), self.tables[name]
        try:
            xyz = tuple(float(t.item(row, c).text()) for c in (2, 3, 4))
        except (ValueError, AttributeError):
            xyz = ()
        if len(xyz) != 3 or not np.all(np.isfinite(xyz)):
            self.refresh_visible(force=True)  # not a usable number: put the old value back
            return
        self.on_edit_coord(self._indices[name][row][0], xyz)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self.refresh_visible()


# ======================================================================
# Trajectory bar
# ======================================================================

class GuiFrameBar(_QWidgetT):
    """First / previous / play / next / last, a slider, the frame count and the speed."""

    BASE_INTERVAL_MS = 100

    def __init__(self, on_frame: Callable[[int], None], parent=None) -> None:
        super().__init__(parent)
        self.on_frame = on_frame
        self.n = 1
        self.speed = 1.0
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(6, 2, 6, 2)
        sp = QtWidgets.QStyle
        self.buttons: Dict[str, QtWidgets.QToolButton] = {}
        for name, icon, tip in (("first", sp.SP_MediaSkipBackward, "First frame (Home)"),
                                ("prev", sp.SP_MediaSeekBackward, "Previous frame ([)"),
                                ("play", sp.SP_MediaPlay, "Play / pause (Space)"),
                                ("next", sp.SP_MediaSeekForward, "Next frame (])"),
                                ("last", sp.SP_MediaSkipForward, "Last frame (End)")):
            b = QtWidgets.QToolButton()
            b.setIcon(self.style().standardIcon(icon))
            b.setToolTip(tip)
            b.setAutoRaise(True)
            self.buttons[name] = b
            lay.addWidget(b)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        lay.addWidget(self.slider, 1)
        self.label = QtWidgets.QLabel()
        self.label.setMinimumWidth(80)
        lay.addWidget(self.label)
        self.buttons["slower"] = QtWidgets.QToolButton()
        self.buttons["slower"].setText("−")
        self.buttons["slower"].setToolTip("Slower")
        self.speed_label = QtWidgets.QLabel("1x")
        self.buttons["faster"] = QtWidgets.QToolButton()
        self.buttons["faster"].setText("+")
        self.buttons["faster"].setToolTip("Faster")
        for w in (self.buttons["slower"], self.speed_label, self.buttons["faster"]):
            lay.addWidget(w)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(self.BASE_INTERVAL_MS)
        self.timer.timeout.connect(self.step)
        self.buttons["first"].clicked.connect(lambda: self.go(0))
        self.buttons["last"].clicked.connect(lambda: self.go(self.n - 1))
        self.buttons["prev"].clicked.connect(lambda: self.go(self.slider.value() - 1))
        self.buttons["next"].clicked.connect(lambda: self.go(self.slider.value() + 1))
        self.buttons["play"].clicked.connect(self.toggle_play)
        self.buttons["slower"].clicked.connect(lambda: self.set_speed(self.speed / 2))
        self.buttons["faster"].clicked.connect(lambda: self.set_speed(self.speed * 2))
        self.slider.valueChanged.connect(self._changed)
        self.setVisible(False)

    @property
    def index(self) -> int:
        return self.slider.value()

    @property
    def playing(self) -> bool:
        return self.timer.isActive()

    def set_count(self, n: int, current: int = 0) -> None:
        self.stop(notify=False)
        self.n = max(n, 1)
        self.slider.blockSignals(True)
        self.slider.setRange(0, self.n - 1)
        self.slider.setValue(min(current, self.n - 1))
        self.slider.blockSignals(False)
        self._label()
        self.setVisible(self.n > 1)

    def go(self, idx: int) -> None:
        self.slider.setValue(max(0, min(idx, self.n - 1)))

    def _label(self) -> None:
        self.label.setText(f"{self.slider.value() + 1}/{self.n}")

    def _changed(self, value: int) -> None:
        self._label()
        self.on_frame(value)

    def set_speed(self, speed: float) -> None:
        self.speed = float(min(16.0, max(0.125, speed)))
        self.speed_label.setText(f"{self.speed:g}x")
        self.timer.setInterval(max(1, int(self.BASE_INTERVAL_MS / self.speed)))

    def toggle_play(self) -> None:
        if self.timer.isActive():
            self.stop()
        elif self.n > 1:
            self.timer.start()
            self.buttons["play"].setIcon(
                self.style().standardIcon(QtWidgets.QStyle.SP_MediaPause))

    def stop(self, notify: bool = True) -> None:
        was = self.timer.isActive()
        self.timer.stop()
        self.buttons["play"].setIcon(self.style().standardIcon(QtWidgets.QStyle.SP_MediaPlay))
        if was and notify:  # listeners refresh what they skip while playing
            self.on_frame(self.slider.value())

    def step(self) -> None:
        self.slider.setValue((self.slider.value() + 1) % self.n)

# ======================================================================
# Module: gui
# ======================================================================
"""Desktop GUI for PDB Wizard (Qt / PySide6).

``pdb_wizard --gui [file]`` opens a window with a 3D view, the structure's details, an atom
table, geometry tables and tool panels. It has the terminal viewer's features: every edit
and analysis, the input generator, database search, trajectories, tabs, view styles and
PNG/GIF export, plus sorbate packing, force fields and charges. Computations go through the
same shared modules as the terminal viewer (see gui_tools.py).

The view is drawn with QPainter (depth-sorted, shaded spheres), not OpenGL, so it works on
any display, over SSH X-forwarding, and headless (QT_QPA_PLATFORM=offscreen).

PySide6 is optional (``pip install "pdb-wizard[gui]"``). Without it this module still
imports, so the single-file build keeps working; run_gui() then explains what to install.
"""


import io as _io

import numpy as np


try:
    from PySide6 import QtCore, QtGui, QtWidgets
    _HAS_QT = True
except ImportError:  # optional extra; see run_gui()
    _HAS_QT = False

if _HAS_QT:
    _QWidget = QtWidgets.QWidget
    _QMainWindow = QtWidgets.QMainWindow
    _QTableModel = QtCore.QAbstractTableModel
else:  # keeps the module (and the single-file build) importable without Qt
    _QWidget = _QMainWindow = _QTableModel = object

GUI_INSTALL_HINT = 'The desktop GUI needs PySide6:  pip install "pdb-wizard[gui]"'
GUI_BROKEN_QT_HINT = (
    "This PySide6 release is broken on this Python version: every Qt call that returns "
    "nothing loses a reference to None, and Python stops with 'Fatal Python error: "
    "none_dealloc' while the window is drawing (seen with PySide6 6.12.0 on Python 3.10).\n"
    'Install a working release:  pip install "PySide6!=6.12.0"   (or use Python 3.12+)'
)


def qt_loses_none_references(calls: int = 64) -> bool:
    """True when this Qt binding drops a reference to None on each call of a method that
    returns nothing. Python 3.12+ cannot be hurt by it (None is immortal there)."""
    if not _HAS_QT or sys.version_info >= (3, 12):
        return False
    probe = QtCore.QObject()
    before = sys.getrefcount(None)
    for _ in range(calls):
        probe.setObjectName("probe")
    return before - sys.getrefcount(None) >= calls // 2

# writers by output suffix; MPMC's dialect is .pqr
_WRITERS = {
    ".xyz": write_xyz,
    ".pdb": write_standard_pdb,
    ".pqr": write_mpmc_pdb,
    ".cif": write_cif,
    ".vasp": write_poscar,
    ".poscar": write_poscar,
    ".lmp": write_lammps_data,
    ".data": write_lammps_data,
    ".d12": write_crystal_d12,  # P1; MainWindow.save_crystal writes a symmetry-reduced one
}
_CELL_FORMATS = (".pqr", ".cif", ".vasp", ".poscar", ".d12")  # cannot be written without a cell
_SAVE_FILTER = ("MPMC (*.pqr);;PDB (*.pdb);;XYZ (*.xyz);;CIF (*.cif);;"
                "POSCAR (*.vasp *.poscar);;LAMMPS data (*.lmp *.data);;CRYSTAL23 input (*.d12)")
_OPEN_FILTER = ("Structures (*.pdb *.pqr *.ent *.cif *.xyz *.vasp *.poscar POSCAR CONTCAR "
                "*.lmp *.data *.sdf *.mol2 *.log *.com *.gjf *.zmat *.dcd *.traj *.out);;All files (*)")
_FORCE_FIELDS = [("None", None), ("OPLS-AA / UFF", 0), ("PHAHST", 1)]
_DUMMY_RADIUS = 0.15  # A, drawn size of massless model sites (water M/LP, H2 centres)

# The terminal viewer's style and colour cycles
STYLES = [("cpk", "Ball and stick"), ("licorice", "Licorice"), ("vdw", "Space filling"),
          ("ribbon", "Ribbon (protein backbone)")]
COLOR_MODES = [("element", "Element"), ("charge", "Charge"), ("residue", "Residue"),
               ("chain", "Chain"), ("index", "Index")]
BACKGROUNDS = [("theme", "Window theme"), ("dark", "Dark"), ("light", "Light")]
_DEFAULT_BOND_RADIUS = 0.11
CONTROLS_HELP = """\
Mouse (in the 3D view)
  Left drag            rotate
  Right or middle drag pan (left drag also pans in pan mode)
  Wheel                zoom
  Click an atom        select it (shown in the Structure panel and the Atoms table)
  Shift+click atoms    2 atoms give a distance, 3 an angle, 4 a dihedral (in the status bar)

Keys (click the 3D view first)
  Arrows or h j k l    rotate (pan in pan mode)
  + / -  or  K / J     zoom in / out
  , / .                roll
  r                    reset the view
  c                    centre the structure
  t                    pan mode on/off
  v                    next style
  C                    next colour mode
  i                    dark/light background
  b                    bonds
  #                    atom numbers
  w                    hide water
  P                    polyhedra
  O                    orthographic / perspective
  m                    measure mode: click two atoms to pin their distance

Trajectories (click the 3D view or the frame bar first)
  [ / ]                previous / next frame
  Home / End           first / last frame
  Space                play / pause

Window
  Ctrl+O open    Ctrl+S save as    Ctrl+W close tab
  Ctrl+Z undo    Ctrl+Shift+Z or Ctrl+Y redo
  Ctrl+D database    Ctrl+I input generator    Ctrl+E export image
  Ctrl+G geometry panel    F1 this help
"""
_LICORICE_BOND_RADIUS = 0.15


def save_structure(mol: Molecule, path: str) -> None:
    """Write ``mol`` in the format given by the file suffix (``.pqr`` = MPMC)."""
    writer = _WRITERS.get(Path(path).suffix.lower())
    if writer is None:
        raise ValueError(f"unknown output format: {Path(path).suffix or path}")
    with open(path, "w") as fh:
        writer(mol, fh)


# ======================================================================
# 3D view
# ======================================================================

def _rotation(axis: int, angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    if axis == 0:
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == 1:
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def _mix(a: "QtGui.QColor", b: "QtGui.QColor", t: float) -> "QtGui.QColor":
    t = float(np.clip(t, 0.0, 1.0)) if np.isfinite(t) else 0.0
    return QtGui.QColor(int(a.red() * (1 - t) + b.red() * t),
                        int(a.green() * (1 - t) + b.green() * t),
                        int(a.blue() * (1 - t) + b.blue() * t))


def _scaled(c: "QtGui.QColor", f: float) -> "QtGui.QColor":
    return QtGui.QColor(*(int(np.clip(v * f, 0, 255)) for v in (c.red(), c.green(), c.blue())))


def _ca_chains(mol: Molecule) -> List[List[int]]:
    """Backbone traces for the ribbon style: CA atoms (or C bonded to N and C), split
    into chains where consecutive CAs are more than 5 Å apart (as the terminal viewer)."""
    atoms = mol.atoms
    ca = [i for i, a in enumerate(atoms) if a.name.strip() == "CA"]
    if len(ca) < 2 and mol.bonds:
        nbrs: dict = {}
        for a, b in mol.bonds:
            nbrs.setdefault(a, set()).add(b)
            nbrs.setdefault(b, set()).add(a)
        ca = [i for i, a in enumerate(atoms) if a.element.symbol == "C"
              and any(atoms[j].element.symbol == "N" for j in nbrs.get(i, ()))
              and any(atoms[j].element.symbol == "C" for j in nbrs.get(i, ()))]
    if len(ca) < 2:
        return []
    chains, cur = [], [ca[0]]
    for prev, i in zip(ca, ca[1:]):
        if np.linalg.norm(atoms[i].position - atoms[prev].position) > 5.0:
            chains.append(cur)
            cur = []
        cur.append(i)
    chains.append(cur)
    return [c for c in chains if len(c) >= 2]


class StructureView(_QWidget):
    """3D structure view drawn with QPainter.

    Left drag rotates (pans in pan mode), right or middle drag pans, the wheel zooms and a
    click selects the atom under the cursor (a double click is just two clicks; r resets
    the view). In measure mode two clicks pin a distance; shift-click up to four atoms for a
    distance, angle or dihedral.
    The terminal viewer's keys work when the view has focus (b, v, i, #, C, w, m, P, O, t,
    c, r, arrows, +/-, ',' and '.').
    """

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setMinimumSize(320, 240)
        self.setMouseTracking(True)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.mol: Optional[Molecule] = None
        self.show_bonds = True
        self.show_cell = True
        self.depth_cue = True
        self.style_name = "cpk"
        self.color_mode = "element"
        self.background = "theme"
        self.show_atom_numbers = False
        self.hide_water = False
        self.show_polyhedra = False
        self.orthographic = True
        self.pan_mode = False
        self.measure_mode = False
        self.atom_scale = 1.0
        self.bond_radius = _DEFAULT_BOND_RADIUS
        self.ambient, self.diffuse, self.specular, self.shininess = 0.50, 0.60, 0.40, 32.0
        self.selected: Optional[int] = None
        self.highlighted: set = set()
        self.picked: List[int] = []
        self.measurements: List[tuple] = []  # (atom indices, label)
        self.density: Optional[tuple] = None  # (positions, values)
        self.show_density = True
        self.on_hover: Optional[Callable[[Optional[int]], None]] = None
        self.on_select: Optional[Callable[[Optional[int]], None]] = None
        self.on_message: Optional[Callable[[str], None]] = None
        self.on_option: Optional[Callable[[], None]] = None  # a key changed a view option
        self._rot = np.eye(3)
        self._zoom = 20.0
        self._pan = QtCore.QPointF(0, 0)
        self._center = np.zeros(3)
        self._radius = 10.0
        self._fit_pending = False
        self._drag_pos: Optional[QtCore.QPointF] = None
        self._drag_button = None
        self._screen = np.zeros((0, 2))
        self._radii_px = np.zeros(0)
        self._depth = np.zeros(0)
        self._visible = np.zeros(0, dtype=bool)
        self._bonds: list = []
        self._periodic = np.zeros(0, dtype=bool)
        self._d_min = np.zeros((0, 3))
        self._cache: dict = {}

    # ---- data ---------------------------------------------------------------------------

    def set_molecule(self, mol: Optional[Molecule], reset_view: bool = True,
                     redetect: bool = True) -> None:
        self.mol = mol
        self.selected = None
        self.picked = []
        self.highlighted = set()
        n = len(mol.atoms) if mol is not None else 0
        if reset_view:
            self.measurements = []
        else:  # frames of one trajectory: keep pins whose atoms still exist
            self.measurements = [m for m in self.measurements if all(i < n for i in m[0])]
        self.refresh_bonds(redetect)
        if reset_view:
            self.reset_view()
        self.update()

    def refresh_bonds(self, redetect: bool = True) -> None:
        """Re-detect bonds (after an edit) and rebuild the drawing caches."""
        self._cache = {}
        mol = self.mol
        if mol is None or not mol.atoms:
            self._bonds, self._periodic, self._d_min = [], np.zeros(0, dtype=bool), np.zeros((0, 3))
            return
        if redetect or not mol.bonds:
            mol.detect_bonds()
        dummy = {i for i, a in enumerate(mol.atoms) if _is_mpmc_dummy_site(a)}
        self._bonds = [(i, j) for i, j in mol.bonds if i not in dummy and j not in dummy]
        pos = np.array([a.x for a in mol.atoms])
        self._periodic, self._d_min = classify_periodic_bonds(pos, self._bonds, mol.pbc)

    def invalidate(self) -> None:
        self._cache = {}
        self.update()

    def has_cell(self) -> bool:
        """A real cell to draw. Trajectories without one get a 10^6 Å stand-in box from
        the reader; drawing or fitting the view to that would shrink the atoms to a dot."""
        return _gui_real_cell(self.mol)

    def cell_corners(self) -> np.ndarray:
        """Cell corners around the atoms: the box is placed in whichever periodic image
        the structure is stored (0..1, -1/2..1/2, ...), snapped to half cells."""
        pbc = self.mol.pbc
        frac = np.array([a.x for a in self.mol.atoms]) @ pbc.reciprocal_basis_matrix
        origin = np.round((frac.mean(axis=0) - 0.5) * 2) / 2
        return pbc.corners() + origin @ pbc.basis_matrix

    def reset_view(self) -> None:
        """Look down c with a slight tilt, structure centred and filling the view."""
        self._rot = _rotation(0, -0.35) @ _rotation(1, 0.45)
        self._pan = QtCore.QPointF(0, 0)
        mol = self.mol
        if mol is None or not mol.atoms:
            return
        pts = np.array([a.x for a in mol.atoms])
        if self.has_cell():
            corners = self.cell_corners()
            self._center = corners.mean(axis=0)
            pts = np.vstack([corners, pts])
        else:
            self._center = pts.mean(axis=0)
        self._radius = max(float(np.linalg.norm(pts - self._center, axis=1).max()), 2.0)
        if self.width() < 50 or self.height() < 50:
            self._fit_pending = True  # not laid out yet: fit on the first paint
            return
        self._fit_pending = False
        self._zoom = 0.47 * min(self.width(), self.height()) / self._radius
        self.update()

    def center_view(self) -> None:
        """Undo panning (the terminal viewer's 'c')."""
        self._pan = QtCore.QPointF(0, 0)
        self.update()

    def rotate(self, d_yaw: float, d_pitch: float) -> None:
        self._rot = _rotation(0, d_pitch) @ _rotation(1, d_yaw) @ self._rot
        self.update()

    def roll(self, angle: float) -> None:
        self._rot = _rotation(2, angle) @ self._rot
        self.update()

    def zoom(self, factor: float) -> None:
        self._zoom = float(np.clip(self._zoom * factor, 0.2, 500.0))
        self.update()

    def view_state(self) -> tuple:
        return self._rot.copy(), self._zoom, QtCore.QPointF(self._pan), self._center.copy()

    def set_view_state(self, state: tuple) -> None:
        self._rot, self._zoom, self._pan, self._center = state
        self._fit_pending = False
        self.update()

    # ---- options ------------------------------------------------------------------------

    def set_style(self, name: str) -> None:
        self.style_name = name
        self.bond_radius = _LICORICE_BOND_RADIUS if name == "licorice" else _DEFAULT_BOND_RADIUS
        self.invalidate()

    def cycle_style(self) -> str:
        names = [k for k, _ in STYLES]
        self.set_style(names[(names.index(self.style_name) + 1) % len(names)])
        return self.style_name

    def cycle_color_mode(self) -> str:
        modes = [k for k, _ in COLOR_MODES]
        self.color_mode = modes[(modes.index(self.color_mode) + 1) % len(modes)]
        self.invalidate()
        return self.color_mode

    def set_measure_mode(self, on: bool) -> None:
        self.measure_mode = on
        self.picked = []
        self.highlighted = set()
        self.update()

    def clear_measurements(self) -> None:
        self.measurements = []
        self.picked = []
        self.update()

    def set_density(self, positions, values) -> None:
        self.density = None if positions is None else (np.asarray(positions),
                                                       np.asarray(values, float))
        self.update()

    # ---- projection ---------------------------------------------------------------------

    def _camera(self) -> float:
        return 3.0 * self._radius

    def _project(self, xyz: np.ndarray) -> tuple:
        """Screen points, depth (larger is nearer) and per-point perspective scale."""
        p = (np.asarray(xyz, dtype=float).reshape(-1, 3) - self._center) @ self._rot.T
        if self.orthographic:
            f = np.ones(len(p))
        else:
            cam = self._camera()
            f = cam / np.maximum(cam - p[:, 2], 0.1 * cam)
        cx = self.width() / 2 + self._pan.x()
        cy = self.height() / 2 + self._pan.y()
        return (np.column_stack([cx + p[:, 0] * self._zoom * f, cy - p[:, 1] * self._zoom * f]),
                p[:, 2], f)

    def atom_radius(self, atom: Atom) -> float:
        """Drawn radius in A for the current style."""
        if _is_mpmc_dummy_site(atom):
            return _DUMMY_RADIUS
        if self.style_name == "licorice":
            return self.bond_radius
        if self.style_name == "vdw":
            return atom.element.vdw_radius * self.atom_scale
        return 0.28 * atom.element.vdw_radius * self.atom_scale

    def atom_at(self, x: float, y: float) -> Optional[int]:
        """Index of the front-most atom drawn under the point, or None."""
        if not len(self._screen):
            return None
        if self.mol is None or len(self._screen) != len(self.mol.atoms):
            return None  # the atoms changed since the last paint
        d = np.hypot(self._screen[:, 0] - x, self._screen[:, 1] - y)
        hit = np.where((d <= np.maximum(self._radii_px, 3.0)) & self._visible)[0]
        if not len(hit):
            return None
        return int(hit[np.argmax(self._depth[hit])])

    # ---- cached per-structure data ------------------------------------------------------

    def _cached(self, key: str, fn: Callable):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    def atom_colors(self) -> List["QtGui.QColor"]:
        mol = self.mol

        def build():
            if self.color_mode == "element":
                rgb = [a.element.cpk_color for a in mol.atoms]
            else:  # the terminal renderer's palettes, so both viewers colour alike
                rgb = ImageRenderer(1, 1)._compute_atom_colors(mol, self.color_mode)
            return [QtGui.QColor(*c) for c in rgb]

        return self._cached(f"colors:{self.color_mode}", build)

    def atom_color(self, atom: Atom) -> "QtGui.QColor":
        return self.atom_colors()[self.mol.atoms.index(atom)]

    def _hidden(self) -> set:
        if not self.hide_water:
            return set()
        return self._cached("water", lambda: water_atom_indices(self.mol))

    def _background(self) -> "QtGui.QColor":
        if self.background == "dark":
            return QtGui.QColor(14, 16, 20)
        if self.background == "light":
            return QtGui.QColor(250, 250, 250)
        return self.palette().color(QtGui.QPalette.Base)

    def _text_color(self) -> "QtGui.QColor":
        bg = self._background()
        return QtGui.QColor(30, 30, 30) if bg.lightness() > 140 else QtGui.QColor(225, 228, 232)

    def _fog(self, color: "QtGui.QColor", t: float) -> "QtGui.QColor":
        """Blend far objects (t near 0) toward the background."""
        if not self.depth_cue:
            return color
        return _mix(color, self._background(), 0.55 * (1.0 - t))

    def _shaded(self, center, r: float, base: "QtGui.QColor"):
        """Sphere shading from the ambient / diffuse / specular / shininess settings."""
        lit = _scaled(base, self.ambient + self.diffuse)
        g = QtGui.QRadialGradient(center, r, center + QtCore.QPointF(-0.35 * r, -0.4 * r))
        g.setColorAt(0.0, _mix(lit, QtGui.QColor(255, 255, 255), self.specular * 1.3))
        g.setColorAt(float(np.clip(0.9 / np.sqrt(self.shininess / 8.0), 0.08, 0.9)), lit)
        g.setColorAt(1.0, _scaled(base, self.ambient * 1.08))
        return g

    # ---- painting -----------------------------------------------------------------------

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt API)
        painter = QtGui.QPainter(self)
        try:  # an exception with the painter still active crashes Qt on the next paint
            self._paint(painter)
        finally:
            painter.end()

    def _paint(self, painter) -> None:
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        bg = self._background()
        grad = QtGui.QRadialGradient(self.width() / 2, self.height() / 2,
                                     0.75 * max(self.width(), self.height()))
        grad.setColorAt(0.0, bg.lighter(115) if bg.lightness() < 128 else bg)
        grad.setColorAt(1.0, bg.darker(115))
        painter.fillRect(self.rect(), grad)
        mol = self.mol
        if self._fit_pending and mol is not None and mol.atoms:
            self.reset_view()
        if mol is None or not mol.atoms:
            painter.setPen(self.palette().color(QtGui.QPalette.PlaceholderText))
            painter.drawText(self.rect(), QtCore.Qt.AlignCenter,
                             "Open a structure (File → Open or drop a file here)")
            self._screen = np.zeros((0, 2))
            return
        pos = np.array([a.x for a in mol.atoms])
        screen, depth, persp = self._project(pos)
        radii = np.array([self.atom_radius(a) for a in mol.atoms])
        hidden = self._hidden()
        visible = np.ones(len(mol.atoms), dtype=bool)
        if hidden:
            visible[list(hidden)] = False
        self._screen, self._depth, self._visible = screen, depth, visible
        self._radii_px = radii * self._zoom * persp
        zmin, zmax = float(depth.min()), float(depth.max())
        span = max(zmax - zmin, 1e-6)
        colors = self.atom_colors()

        if self.show_cell and self.has_cell():
            self._paint_cell(painter, mol)
        if self.density is not None and self.show_density:
            self._paint_density(painter)
        chains = _ca_chains(mol) if self.style_name == "ribbon" else []
        if chains:
            self._paint_ribbon(painter, chains, screen)
        else:
            if self.show_bonds and self._bonds and self.style_name != "vdw":
                self._paint_bonds(painter, pos, screen, depth, persp, zmin, span, colors,
                                  visible)
            self._paint_atoms(painter, mol, screen, depth, zmin, span, colors, visible)
        if self.show_polyhedra:
            self._paint_polyhedra(painter, mol)
        self._paint_marks(painter, mol, screen)
        if self.show_atom_numbers:
            self._paint_numbers(painter, screen, depth, visible)
        self._paint_measurements(painter, mol, screen)
        self._paint_axes(painter, mol)

    def _paint_atoms(self, painter, mol, screen, depth, zmin, span, colors, visible) -> None:
        outline = QtGui.QPen(QtGui.QColor(0, 0, 0, 90), 0.8)
        for i in np.argsort(depth):  # far to near
            if not visible[i]:
                continue
            atom = mol.atoms[i]
            r = float(self._radii_px[i])
            c = QtCore.QPointF(screen[i, 0], screen[i, 1])
            base = self._fog(colors[i], (depth[i] - zmin) / span)
            if _is_mpmc_dummy_site(atom):
                painter.setPen(QtCore.Qt.NoPen)
                painter.setBrush(QtGui.QColor(base.red(), base.green(), base.blue(), 150))
                painter.drawEllipse(c, max(r, 1.5), max(r, 1.5))
                continue
            painter.setBrush(self._shaded(c, r, base))
            painter.setPen(outline if r > 4 else QtCore.Qt.NoPen)
            painter.drawEllipse(c, r, r)

    def _paint_cell(self, painter, mol) -> None:
        corners = self._project(self.cell_corners())[0]
        color = self._text_color()
        color.setAlpha(110)
        painter.setPen(QtGui.QPen(color, 1.2, QtCore.Qt.DashLine))
        for i, j in mol.pbc.edges():
            painter.drawLine(QtCore.QPointF(*corners[i]), QtCore.QPointF(*corners[j]))

    def _paint_bonds(self, painter, pos, screen, depth, persp, zmin, span, colors,
                     visible) -> None:
        segments = []  # (depth, p, q, atom index, perspective)
        for k, (i, j) in enumerate(self._bonds):
            if not (visible[i] and visible[j]):
                continue
            if self._periodic[k]:  # stubs toward each atom's nearest image
                half = self._d_min[k] / 2
                ends = self._project(np.array([pos[i] + half, pos[j] - half]))[0]
                segments.append((depth[i], screen[i], ends[0], i))
                segments.append((depth[j], screen[j], ends[1], j))
            else:
                mid = (screen[i] + screen[j]) / 2
                segments.append((depth[i], screen[i], mid, i))
                segments.append((depth[j], screen[j], mid, j))
        segments.sort(key=lambda s: s[0])
        for z, p, q, idx in segments:
            width = max(1.5, 2 * self.bond_radius * self._zoom * float(persp[idx]))
            color = self._fog(colors[idx], (z - zmin) / span).darker(115)
            painter.setPen(QtGui.QPen(color, width, QtCore.Qt.SolidLine, QtCore.Qt.RoundCap))
            painter.drawLine(QtCore.QPointF(*p), QtCore.QPointF(*q))

    def _paint_ribbon(self, painter, chains, screen) -> None:
        palette = [(50, 130, 255), (255, 100, 80), (80, 220, 100), (255, 200, 50),
                   (200, 100, 255), (100, 220, 220), (255, 150, 200), (180, 180, 100)]
        width = max(2.0, 0.5 * self._zoom)
        for k, chain in enumerate(chains):
            pts = screen[chain]
            path = QtGui.QPainterPath(QtCore.QPointF(*pts[0]))
            for a in range(len(pts) - 1):  # Catmull-Rom through the CA positions
                p0, p1 = pts[max(a - 1, 0)], pts[a]
                p2, p3 = pts[a + 1], pts[min(a + 2, len(pts) - 1)]
                c1, c2 = p1 + (p2 - p0) / 6, p2 - (p3 - p1) / 6
                path.cubicTo(QtCore.QPointF(*c1), QtCore.QPointF(*c2), QtCore.QPointF(*p2))
            painter.setPen(QtGui.QPen(QtGui.QColor(*palette[k % len(palette)]), width,
                                      QtCore.Qt.SolidLine, QtCore.Qt.RoundCap,
                                      QtCore.Qt.RoundJoin))
            painter.setBrush(QtCore.Qt.NoBrush)
            painter.drawPath(path)

    def _paint_polyhedra(self, painter, mol) -> None:
        polys = self._cached("polyhedra", lambda: coordination_polyhedra(mol.atoms, mol.pbc))
        colors = self.atom_colors()
        for ci, nbrs in polys:
            pts = self._project(nbrs)[0]
            c = colors[ci]
            fill = QtGui.QColor(c.red(), c.green(), c.blue(), 40)
            edge = QtGui.QColor(c.red() // 2, c.green() // 2, c.blue() // 2, 200)
            painter.setPen(QtCore.Qt.NoPen)
            painter.setBrush(fill)
            painter.drawPolygon(QtGui.QPolygonF(_convex_hull(pts)))
            painter.setPen(QtGui.QPen(edge, 1.0))
            for a in range(len(pts)):
                for b in range(a + 1, len(pts)):
                    painter.drawLine(QtCore.QPointF(*pts[a]), QtCore.QPointF(*pts[b]))

    def _paint_density(self, painter) -> None:
        positions, values = self.density
        if not len(values):
            return
        pts, depth, persp = self._project(positions)
        vmax = float(values.max()) or 1.0
        for i in np.argsort(depth):
            t = float(values[i] / vmax)
            color = _mix(QtGui.QColor(60, 120, 255), QtGui.QColor(255, 70, 60), t)
            color.setAlpha(int(40 + 140 * t))
            r = (2.0 + 6.0 * t) * float(persp[i])
            painter.setPen(QtCore.Qt.NoPen)
            painter.setBrush(color)
            painter.drawEllipse(QtCore.QPointF(*pts[i]), r, r)

    def _paint_marks(self, painter, mol, screen) -> None:
        rings = [(i, self.palette().color(QtGui.QPalette.Highlight))
                 for i in ([self.selected] if self.selected is not None else [])]
        marked = sorted(self.highlighted | set(self.picked))
        rings += [(i, QtGui.QColor(255, 210, 60)) for i in marked]
        painter.setBrush(QtCore.Qt.NoBrush)
        for i, color in rings:
            if i < len(mol.atoms):
                painter.setPen(QtGui.QPen(color, 2.5))
                r = float(self._radii_px[i]) + 4
                painter.drawEllipse(QtCore.QPointF(screen[i, 0], screen[i, 1]), r, r)

    def _paint_numbers(self, painter, screen, depth, visible) -> None:
        painter.setPen(self._text_color())
        order = [i for i in np.argsort(depth) if visible[i]][-3000:]  # nearest 3000
        for i in order:
            painter.drawText(QtCore.QPointF(screen[i, 0] + 0.6 * self._radii_px[i],
                                            screen[i, 1] - 0.6 * self._radii_px[i]), str(i + 1))

    def _paint_measurements(self, painter, mol, screen) -> None:
        if not self.measurements:
            return
        accent = QtGui.QColor(255, 210, 60)
        for idx, _ in self.measurements:
            if any(i >= len(mol.atoms) for i in idx):
                continue
            label = f"{self.distance(*idx):.3f} Å" if len(idx) == 2 else ""
            pts = [QtCore.QPointF(*screen[i]) for i in idx]
            painter.setPen(QtGui.QPen(accent, 1.5, QtCore.Qt.DashLine))
            for a, b in zip(pts, pts[1:]):
                painter.drawLine(a, b)
            mid = sum((p for p in pts), QtCore.QPointF(0, 0)) / len(pts)
            painter.setPen(accent)
            painter.drawText(mid + QtCore.QPointF(6, -6), label)

    def _paint_axes(self, painter, mol) -> None:
        if not self.has_cell():
            return
        origin = QtCore.QPointF(42, self.height() - 42)
        colors = [QtGui.QColor(230, 80, 80), QtGui.QColor(90, 190, 90), QtGui.QColor(90, 140, 240)]
        for k, (name, vec) in enumerate(zip("abc", mol.pbc.basis_matrix)):
            d = (vec / np.linalg.norm(vec)) @ self._rot.T
            tip = origin + QtCore.QPointF(26 * d[0], -26 * d[1])
            painter.setPen(QtGui.QPen(colors[k], 2.2, QtCore.Qt.SolidLine, QtCore.Qt.RoundCap))
            painter.drawLine(origin, tip)
            painter.drawText(tip + QtCore.QPointF(3 * np.sign(d[0]) - 3, 4), name)

    # ---- measuring ----------------------------------------------------------------------

    def _vec(self, i: int, j: int) -> np.ndarray:
        d = self.mol.atoms[i].x - self.mol.atoms[j].x
        return self.mol.pbc.wrap(d) if self.mol.pbc is not None else d

    def distance(self, i: int, j: int) -> float:
        d = self.mol.atoms[i].x - self.mol.atoms[j].x
        if self.mol.pbc is not None:
            return self.mol.pbc.min_image(d)
        return float(np.linalg.norm(d))

    def angle(self, i: int, j: int, k: int) -> float:
        v1, v2 = self._vec(i, j), self._vec(k, j)
        cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        return float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))

    def dihedral(self, i: int, j: int, k: int, m: int) -> float:
        """Dihedral i-j-k-m in degrees, unsigned (0-180), minimum image when there is a
        cell: the formula of Molecule.get_dihedrals, so it matches the Geometry panel."""
        b1, b2, b3 = self._vec(j, i), self._vec(k, j), self._vec(m, k)
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        cos_d = np.dot(n1, n2) / ((np.linalg.norm(n1) + 1e-10) * (np.linalg.norm(n2) + 1e-10))
        return float(np.degrees(np.arccos(np.clip(cos_d, -1.0, 1.0))))

    def clear_marks(self) -> None:
        """Forget selection, picks, highlights and pinned measurements (after an edit the
        indices may point at other atoms, or past the end)."""
        self.selected = None
        self.picked = []
        self.highlighted = set()
        self.measurements = []

    def pick(self, idx: int, shift: bool = False) -> str:
        """Measure-mode / shift-click picking, as in the terminal viewer. Returns the
        message to show."""
        n = len(self.mol.atoms) if self.mol is not None else 0
        self.picked = [i for i in self.picked if i < n]
        if not 0 <= idx < n:
            return ""
        a = self.mol.atoms[idx]
        if self.measure_mode:
            if idx not in self.picked:
                self.picked.append(idx)
            if len(self.picked) == 1:
                msg = f"Atom {idx + 1}: {a.element.symbol} — click another to measure"
            else:
                i, j = self.picked[:2]
                d = self.distance(i, j)
                self.measurements.append(((i, j), f"{d:.3f} Å"))
                msg = f"Distance {i + 1}-{j + 1} = {d:.4f} Å (pinned)"
                self.picked = []
        else:
            if idx not in self.picked:
                self.picked.append(idx)
            del self.picked[:-4]  # a fifth atom drops the first
            p = self.picked
            numbers = "-".join(str(i + 1) for i in p)
            if len(p) == 1:
                msg = f"Atom {p[0] + 1}: {a.element.symbol} ({a.name})"
            elif len(p) == 2:
                msg = f"Distance {numbers} = {self.distance(*p):.4f} Å"
            elif len(p) == 3:
                msg = f"Angle {numbers} = {self.angle(*p):.2f}°"
            else:
                msg = f"Dihedral {numbers} = {self.dihedral(*p):.2f}°"
        self.update()
        return msg

    # ---- input --------------------------------------------------------------------------

    def mousePressEvent(self, event) -> None:  # noqa: N802
        self._drag_pos = event.position()
        self._drag_button = event.button()
        self._pressed_at = event.position()
        self.setFocus()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        p = event.position()
        if self._drag_pos is not None:
            dx, dy = p.x() - self._drag_pos.x(), p.y() - self._drag_pos.y()
            if self._drag_button == QtCore.Qt.LeftButton and not self.pan_mode:
                self.rotate(dx * 0.01, dy * 0.01)
            else:
                self._pan += QtCore.QPointF(dx, dy)
                self.update()
            self._drag_pos = p
        elif self.on_hover is not None:
            self.on_hover(self.atom_at(p.x(), p.y()))

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        p = event.position()
        moved = (p - getattr(self, "_pressed_at", p)).manhattanLength()
        if moved < 4 and event.button() == QtCore.Qt.LeftButton:
            idx = self.atom_at(p.x(), p.y())
            shift = bool(event.modifiers() & QtCore.Qt.ShiftModifier)
            if idx is not None and (self.measure_mode or shift):
                msg = self.pick(idx, shift)
                if self.on_message is not None:
                    self.on_message(msg)
            else:
                if not shift:
                    self.picked = []
                self.selected = idx
                if self.on_select is not None:
                    self.on_select(self.selected)
            self.update()
        self._drag_pos = None

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        # Qt sends press, release, double click, release: treat the double click as the
        # second press, so it selects or picks like any click and never moves the camera
        self.mousePressEvent(event)

    def wheelEvent(self, event) -> None:  # noqa: N802
        self.zoom(1.12 if event.angleDelta().y() > 0 else 1 / 1.12)

    def keyPressEvent(self, event) -> None:  # noqa: N802
        key, text = event.key(), event.text()
        step = 0.1
        moves = {QtCore.Qt.Key_Left: (-step, 0), QtCore.Qt.Key_Right: (step, 0),
                 QtCore.Qt.Key_Up: (0, -step), QtCore.Qt.Key_Down: (0, step),
                 QtCore.Qt.Key_H: (-step, 0), QtCore.Qt.Key_L: (step, 0),
                 QtCore.Qt.Key_K: (0, -step), QtCore.Qt.Key_J: (0, step)}
        toggles = {"b": "show_bonds", "#": "show_atom_numbers", "w": "hide_water",
                   "P": "show_polyhedra", "O": "orthographic", "t": "pan_mode"}
        if key in moves and text not in ("K", "J"):
            dx, dy = moves[key]
            if self.pan_mode:
                self._pan += QtCore.QPointF(-dx * 300, -dy * 300)
                self.update()
            else:
                self.rotate(dx, dy)
        elif key in (QtCore.Qt.Key_Plus, QtCore.Qt.Key_Equal) or text == "K":
            self.zoom(1.12)
        elif key == QtCore.Qt.Key_Minus or text == "J":
            self.zoom(1 / 1.12)
        elif text == ",":
            self.roll(0.1)
        elif text == ".":
            self.roll(-0.1)
        elif text == "r":
            self.reset_view()
        elif text == "c":
            self.center_view()
        elif text == "v":
            self.cycle_style()
        elif text == "C":
            self.cycle_color_mode()
        elif text == "i":
            self.background = "light" if self._background().lightness() < 128 else "dark"
            self.update()
        elif text == "m":
            self.set_measure_mode(not self.measure_mode)
        elif text in toggles:
            setattr(self, toggles[text], not getattr(self, toggles[text]))
            self.invalidate()
        else:
            super().keyPressEvent(event)
            return
        if self.on_option is not None:
            self.on_option()

    def resizeEvent(self, event) -> None:  # noqa: N802
        old, new = event.oldSize(), event.size()
        if old.width() > 50 and old.height() > 50 and not self._fit_pending:
            self._zoom *= min(new.width(), new.height()) / max(min(old.width(), old.height()), 1)
        super().resizeEvent(event)


def _qt_alive(obj) -> bool:
    """False once Qt has deleted the C++ side of a widget (closed with WA_DeleteOnClose)."""
    try:
        obj.isVisible()
    except RuntimeError:
        return False
    return True


def _cell_problem(a: float, b: float, c: float, alpha: float, beta: float,
                  gamma: float) -> Optional[str]:
    """Why these cell parameters are unusable, or None."""
    values = (a, b, c, alpha, beta, gamma)
    if not all(np.isfinite(v) for v in values):
        return "Cell parameters must be finite numbers"
    if min(a, b, c) <= 0:
        return "Cell lengths must be positive"
    if not all(0 < x < 180 for x in (alpha, beta, gamma)):
        return "Cell angles must be between 0 and 180 degrees"
    ca, cb, cg = (np.cos(np.radians(x)) for x in (alpha, beta, gamma))
    if 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg <= 1e-10:
        return "These angles do not make a cell (zero or negative volume)"
    return None


def _convex_hull(pts: np.ndarray) -> list:
    """Screen-space convex hull (monotone chain) as QPointFs, for polyhedron shading."""
    p = sorted(map(tuple, np.asarray(pts, float)))
    if len(p) < 3:
        return [QtCore.QPointF(*q) for q in p]

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for q in p:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], q) <= 0:
            lower.pop()
        lower.append(q)
    for q in reversed(p):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], q) <= 0:
            upper.pop()
        upper.append(q)
    return [QtCore.QPointF(*q) for q in lower[:-1] + upper[:-1]]


# ======================================================================
# Atom table
# ======================================================================

class AtomTableModel(_QTableModel):
    HEADERS = ["#", "Element", "Name", "Molecule", "x", "y", "z", "q (e)", "α (Å³)",
               "ε (K)", "σ (Å)"]

    def __init__(self) -> None:
        super().__init__()
        self.mol: Optional[Molecule] = None

    def set_molecule(self, mol: Optional[Molecule]) -> None:
        self.beginResetModel()
        self.mol = mol
        self.endResetModel()

    def rowCount(self, parent=None) -> int:  # noqa: N802
        return 0 if self.mol is None else len(self.mol.atoms)

    def columnCount(self, parent=None) -> int:  # noqa: N802
        return len(self.HEADERS)

    def headerData(self, section, orientation, role=None):  # noqa: N802
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self.HEADERS[section]
        return None

    def data(self, index, role=None):
        if role not in (QtCore.Qt.DisplayRole, QtCore.Qt.TextAlignmentRole) or self.mol is None:
            return None
        col = index.column()
        if role == QtCore.Qt.TextAlignmentRole:
            return int(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter) if col >= 4 else None
        a = self.mol.atoms[index.row()]
        if col == 0:
            return index.row() + 1
        if col == 1:
            return a.element.symbol
        if col == 2:
            return a.name
        if col == 3:
            return f"M {a.mpmc_mol} ({a.mpmc_res})" if a.mpmc_flag == "M" else "framework"
        values = [a.x[0], a.x[1], a.x[2], a.charge, a.alpha, a.epsilon, a.sigma]
        return f"{values[col - 4]:.4f}" if col in (4, 5, 6, 7) else f"{values[col - 4]:.3f}"


# ======================================================================
# Documents (tabs)
# ======================================================================

class GuiDocument:
    """One open structure or trajectory: what a tab holds."""

    def __init__(self, mol: Optional[Molecule] = None, path: Optional[str] = None,
                 frames: Optional[list] = None, name: str = "") -> None:
        self.mol = mol
        self.path = path
        self.frames = frames
        self.frame_idx = 0
        self.undo: list = []  # (frame index, Molecule), newest last
        self.redo: list = []  # states undo replaced, same layout; a new edit clears it
        self.name = name or (Path(path).name if path else "untitled")
        self.view_state: Optional[tuple] = None
        self.density: Optional[tuple] = None
        self.measurements: list = []
        self.loading = False  # a trajectory is still being read into this tab

    @property
    def empty(self) -> bool:
        return self.mol is None or not self.mol.atoms


def _doc_attr(name: str):
    return property(lambda self: getattr(self.doc, name),
                    lambda self, value: setattr(self.doc, name, value))


# ======================================================================
# Main window
# ======================================================================

class MainWindow(_QMainWindow):
    """PDB Wizard's desktop window."""

    mol = _doc_attr("mol")
    path = _doc_attr("path")
    frames = _doc_attr("frames")
    _undo = _doc_attr("undo")
    _redo = _doc_attr("redo")

    def __init__(self, path: Optional[str] = None) -> None:
        super().__init__()
        self.docs: List[GuiDocument] = [GuiDocument()]
        self.doc_index = 0
        self._redo_stash: Optional[tuple] = None  # (document, redo list) during an edit
        self.dialogs: list = []  # keeps open tool windows alive
        self.load_task: Optional[GuiTask] = None
        self.task: Optional[GuiTask] = None
        self._db_cache: Optional[str] = None
        self.setAcceptDrops(True)
        self.resize(1360, 860)
        self.view = StructureView()
        self.view.on_hover = self._hovered
        self.view.on_select = self._selected
        self.view.on_message = lambda msg: self.statusBar().showMessage(msg, 8000)
        self.view.on_option = self._sync_view_controls
        self.tab_bar = QtWidgets.QTabBar()
        self.tab_bar.setTabsClosable(True)
        self.tab_bar.setDocumentMode(True)
        self.tab_bar.setExpanding(False)
        self.tab_bar.addTab("untitled")
        self.tab_bar.currentChanged.connect(self._switch_doc)
        self.tab_bar.tabCloseRequested.connect(self.close_tab)
        self.frame_bar = GuiFrameBar(self.set_frame)
        central = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(central)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.tab_bar)
        lay.addWidget(self.view, 1)
        lay.addWidget(self.frame_bar)
        self.setCentralWidget(central)
        self._build_info_dock()
        self._build_tools_dock()
        self._build_table_dock()
        self._build_geometry_dock()
        self._build_actions()
        self.statusBar().showMessage("Ready")
        self._update_title()
        self._update_supercell_preview()
        if path:
            self.open_path(path)

    @property
    def doc(self) -> GuiDocument:
        return self.docs[self.doc_index]

    # ---- layout -------------------------------------------------------------------------

    def _build_info_dock(self) -> None:
        dock = QtWidgets.QDockWidget("Structure", self)
        dock.setObjectName("structure")
        body = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(body)
        self.info_form = QtWidgets.QFormLayout()
        self.info_labels = {}
        for key in ("File", "Formula", "Atoms", "Sorbates", "a, b, c (Å)",
                    "α, β, γ (°)", "Volume (Å³)", "Density (g/cm³)",
                    "Net charge (e)", "Frame"):
            label = QtWidgets.QLabel("—")
            label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            label.setWordWrap(True)
            self.info_labels[key] = label
            self.info_form.addRow(f"<b>{key}</b>", label)
        lay.addLayout(self.info_form)
        lay.addWidget(QtWidgets.QLabel("<b>Elements</b>"))
        self.legend = QtWidgets.QListWidget()
        self.legend.setFrameShape(QtWidgets.QFrame.NoFrame)
        lay.addWidget(self.legend, 1)
        self.selection_label = QtWidgets.QLabel("Click an atom to inspect it.")
        self.selection_label.setWordWrap(True)
        lay.addWidget(self.selection_label)
        dock.setWidget(body)
        dock.setMinimumWidth(250)
        self.info_dock = dock
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)

    def _build_tools_dock(self) -> None:
        dock = QtWidgets.QDockWidget("Tools", self)
        dock.setObjectName("tools")
        tools = QtWidgets.QTabWidget()
        tools.setDocumentMode(True)
        tools.addTab(self._view_page(), "View")

        # supercell
        page = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(page)
        self.extend_spins = []
        for axis in "abc":
            spin = QtWidgets.QSpinBox()
            spin.setRange(1, 20)
            spin.setToolTip("Unit cells along this axis in the result (1 leaves it as it is)")
            spin.valueChanged.connect(lambda _: self._update_supercell_preview())
            self.extend_spins.append(spin)
            form.addRow(f"Cells along {axis}", spin)
        self.supercell_preview = QtWidgets.QLabel()
        form.addRow(self.supercell_preview)
        btn = QtWidgets.QPushButton("Build supercell")
        btn.clicked.connect(self.build_supercell)
        form.addRow(btn)
        wrap_row = QtWidgets.QHBoxLayout()
        for label, forward in (("Wrap (centred)", False), ("Wrap (0–1)", True)):
            b = QtWidgets.QPushButton(label)
            b.clicked.connect(lambda _=False, f=forward: self.wrap(f))
            wrap_row.addWidget(b)
        form.addRow(wrap_row)
        reduce_btn = QtWidgets.QPushButton("Reduce supercell…")
        reduce_btn.clicked.connect(lambda: self.reduce_supercell())
        form.addRow(reduce_btn)
        tools.addTab(page, "Supercell")

        # sorbates
        page = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(page)
        self.sorbate_combo = QtWidgets.QComboBox()
        for name, model in SORBATE_MODELS.items():
            self.sorbate_combo.addItem(f"{model.description}", name)
        self.sorbate_combo.addItem("From an MPMC .pqr file…", "__file__")
        idx = self.sorbate_combo.findData("h2_bss")
        self.sorbate_combo.setCurrentIndex(max(idx, 0))
        form.addRow("Model", self.sorbate_combo)
        self.sorbate_count = QtWidgets.QSpinBox()
        self.sorbate_count.setRange(1, 100000)
        self.sorbate_count.setValue(10)
        self.sorbate_fill = QtWidgets.QCheckBox("Fill the cell")
        self.sorbate_fill.toggled.connect(lambda on: self.sorbate_count.setEnabled(not on))
        count_row = QtWidgets.QHBoxLayout()
        count_row.addWidget(self.sorbate_count)
        count_row.addWidget(self.sorbate_fill)
        form.addRow("Molecules", count_row)
        self.sorbate_scale = QtWidgets.QDoubleSpinBox()
        self.sorbate_scale.setRange(0.5, 1.2)
        self.sorbate_scale.setSingleStep(0.05)
        self.sorbate_scale.setValue(0.85)
        self.sorbate_scale.setToolTip("Closest allowed approach, as a fraction of the LJ σ")
        form.addRow("Overlap scale", self.sorbate_scale)
        btn = QtWidgets.QPushButton("Pack sorbates")
        btn.clicked.connect(lambda: self.pack_sorbates())
        form.addRow(btn)
        hint = QtWidgets.QLabel("Random positions and orientations without overlaps. Relax "
                                "them in MPMC (NVT) before starting μVT.")
        hint.setWordWrap(True)
        form.addRow(hint)
        tools.addTab(page, "Sorbates")

        # force field and charges
        page = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(page)
        self.ff_combo = QtWidgets.QComboBox()
        for label, idx in _FORCE_FIELDS:
            self.ff_combo.addItem(label, idx)
        self.ff_combo.setCurrentIndex(1)
        ff_btn = QtWidgets.QPushButton("Apply force field")
        ff_btn.clicked.connect(self.apply_force_field)
        form.addRow("Force field", self.ff_combo)
        form.addRow(ff_btn)
        q_btn = QtWidgets.QPushButton("Load charges file…")
        q_btn.clicked.connect(lambda: self.load_charges())
        qeq_btn = QtWidgets.QPushButton("QEq charges")
        qeq_btn.clicked.connect(self.qeq_charges)
        form.addRow(q_btn)
        form.addRow(qeq_btn)
        note = QtWidgets.QLabel("Only framework atoms are changed; packed sorbates keep their "
                                "model parameters.")
        note.setWordWrap(True)
        form.addRow(note)
        tools.addTab(page, "Charges && FF")

        self.tools = tools
        dock.setWidget(tools)
        dock.setMinimumWidth(300)
        self.tools_dock = dock
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, dock)

    def _view_page(self) -> "QtWidgets.QWidget":
        page = QtWidgets.QScrollArea()
        page.setWidgetResizable(True)
        page.setFrameShape(QtWidgets.QFrame.NoFrame)
        body = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(body)
        v = self.view
        self.style_combo = QtWidgets.QComboBox()
        for key, label in STYLES:
            self.style_combo.addItem(label, key)
        self.style_combo.currentIndexChanged.connect(
            lambda _: self.set_style(self.style_combo.currentData()))
        form.addRow("Style", self.style_combo)
        self.color_combo = QtWidgets.QComboBox()
        for key, label in COLOR_MODES:
            self.color_combo.addItem(label, key)
        self.color_combo.currentIndexChanged.connect(
            lambda _: self._set_view("color_mode", self.color_combo.currentData(), True))
        form.addRow("Colour by", self.color_combo)
        self.bg_combo = QtWidgets.QComboBox()
        for key, label in BACKGROUNDS:
            self.bg_combo.addItem(label, key)
        self.bg_combo.currentIndexChanged.connect(
            lambda _: self._set_view("background", self.bg_combo.currentData()))
        form.addRow("Background", self.bg_combo)
        self.view_checks = {}
        for attr, label in (("show_bonds", "Bonds"), ("show_cell", "Unit cell"),
                            ("show_atom_numbers", "Atom numbers"), ("hide_water", "Hide water"),
                            ("show_polyhedra", "Polyhedra"), ("orthographic", "Orthographic"),
                            ("depth_cue", "Depth cueing"), ("measure_mode", "Measure mode")):
            chk = QtWidgets.QCheckBox(label)
            chk.setChecked(getattr(v, attr))
            chk.toggled.connect(lambda on, a=attr: self._toggle_view(a, on))
            self.view_checks[attr] = chk
            form.addRow(chk)
        # kept for scripts and tests written against the first version
        self.chk_bonds = self.view_checks["show_bonds"]
        self.chk_cell = self.view_checks["show_cell"]
        self.chk_fog = self.view_checks["depth_cue"]
        self.sliders = {}
        for attr, label, lo, hi, scale in (("atom_scale", "Atom size", 0.3, 2.0, 100),
                                           ("bond_radius", "Bond radius (Å)", 0.02, 0.3, 100),
                                           ("ambient", "Ambient", 0.0, 1.0, 100),
                                           ("diffuse", "Diffuse", 0.0, 1.0, 100),
                                           ("specular", "Specular", 0.0, 1.0, 100),
                                           ("shininess", "Shininess", 1.0, 128.0, 1)):
            s = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            s.setRange(int(lo * scale), int(hi * scale))
            s.setValue(int(round(getattr(v, attr) * scale)))
            s.valueChanged.connect(lambda val, a=attr, k=scale: self._set_view(a, val / k))
            self.sliders[attr] = s
            form.addRow(label, s)
        self.size_slider = self.sliders["atom_scale"]
        row = QtWidgets.QHBoxLayout()
        for label, slot in (("Reset view", self.view.reset_view),
                            ("Clear measurements", self.view.clear_measurements)):
            b = QtWidgets.QPushButton(label)
            b.clicked.connect(slot)
            row.addWidget(b)
        form.addRow(row)
        page.setWidget(body)
        return page

    def _build_table_dock(self) -> None:
        dock = QtWidgets.QDockWidget("Atoms", self)
        dock.setObjectName("atoms")
        self.table_model = AtomTableModel()
        self.table = QtWidgets.QTableView()
        self.table.setModel(self.table_model)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setDefaultSectionSize(self.table.fontMetrics().height() + 8)
        self.table.clicked.connect(lambda idx: self._selected(idx.row(), from_table=True))
        dock.setWidget(self.table)
        self.table_dock = dock
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        self.resizeDocks([dock], [190], QtCore.Qt.Vertical)

    def _build_geometry_dock(self) -> None:
        dock = QtWidgets.QDockWidget("Geometry", self)
        dock.setObjectName("geometry")
        self.geometry = GuiGeometryPanel(self._highlight, self.edit_coord)
        dock.setWidget(self.geometry)
        self.geometry_dock = dock
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, dock)
        self.tabifyDockWidget(self.table_dock, dock)
        self.table_dock.raise_()

    def _build_actions(self) -> None:
        style = self.style()
        bar = self.addToolBar("Main")
        bar.setObjectName("main")
        bar.setMovable(False)
        menus = {name: self.menuBar().addMenu(name) for name in
                 ("&File", "&Edit", "&View", "&Tools", "&Analysis", "&Help")}
        self.actions_by_name = {}

        def action(menu, text, slot, shortcut=None, icon=None, toolbar=False, checkable=False):
            act = QtGui.QAction(text, self)
            if icon is not None:
                act.setIcon(style.standardIcon(icon))
            if shortcut:
                act.setShortcut(QtGui.QKeySequence(shortcut))
            act.setCheckable(checkable)
            act.triggered.connect(slot)
            menus[menu].addAction(act)
            if toolbar:
                bar.addAction(act)
            self.actions_by_name[text.split("\t")[0].replace("&", "").rstrip("…")] = act
            return act

        sp = QtWidgets.QStyle
        f, e, v, t, a = "&File", "&Edit", "&View", "&Tools", "&Analysis"
        action(f, "&Open…", self.open_dialog, "Ctrl+O", sp.SP_DialogOpenButton, True)
        action(f, "Save &as…", self.save_dialog, "Ctrl+S", sp.SP_DialogSaveButton, True)
        menus[f].addSeparator()
        action(f, "Open frame in new &tab", self.frame_to_tab)
        action(f, "Fetch from &database…", lambda: self.open_tool("db_search"), "Ctrl+D")
        menus[f].addSeparator()
        action(f, "Export &image…", lambda: self.export_image(), "Ctrl+E")
        action(f, "Export rotation &GIF…", lambda: self.export_gif())
        menus[f].addSeparator()
        action(f, "&Close tab", lambda: self.close_tab(self.doc_index), "Ctrl+W")
        action(f, "&Quit", self.close, "Ctrl+Q")

        self.undo_action = action(e, "&Undo", self.undo, "Ctrl+Z", sp.SP_ArrowBack, True)
        # "Re&do": R is taken in this menu by Reduce supercell
        self.redo_action = action(e, "Re&do", self.redo, None, sp.SP_ArrowForward, True)
        self.redo_action.setShortcuts([QtGui.QKeySequence("Ctrl+Shift+Z"),
                                       QtGui.QKeySequence("Ctrl+Y")])
        self.undo_action.setEnabled(False)
        self.redo_action.setEnabled(False)
        menus[e].addSeparator()
        action(e, "&Build supercell…", self._show_supercell_tab)
        action(e, "&Reduce supercell…", lambda: self.reduce_supercell())
        action(e, "Wrap (&centred)", lambda: self.wrap(False))
        action(e, "Wrap (&forward)", lambda: self.wrap(True))
        action(e, "Update unit c&ell…", lambda: self.update_cell())
        menus[e].addSeparator()
        action(e, "&Sort atoms", self.sort_atoms)
        action(e, "Substitute ele&ment…", lambda: self.substitute_element())
        action(e, "Delete &lone atoms", lambda: self.delete_lone_atoms(ask=True))
        action(e, "Edit &H distances…", lambda: self.edit_h_distances())
        bar.addSeparator()

        action(v, "&Reset view", self.view.reset_view, "R", sp.SP_BrowserReload, True)
        action(v, "C&entre\tC", self.view.center_view)
        action(v, "Zoom &in", lambda: self.view.zoom(1.2), "Ctrl++")
        action(v, "Zoom &out", lambda: self.view.zoom(1 / 1.2), "Ctrl+-")
        menus[v].addSeparator()
        action(v, "Cycle &style\tV", lambda: self._key_option(self.view.cycle_style))
        action(v, "Cycle &colour mode\tShift+C",
               lambda: self._key_option(self.view.cycle_color_mode))
        menus[v].addSeparator()
        # The same options as the View tab; the keys work while the 3D view has the focus.
        self.view_actions = {}
        for attr, label in (("show_bonds", "&Bonds\tB"), ("show_cell", "&Unit cell"),
                            ("show_atom_numbers", "Atom &numbers\t#"),
                            ("hide_water", "Hide &water\tW"),
                            ("show_polyhedra", "&Polyhedra\tShift+P"),
                            ("orthographic", "Ortho&graphic\tShift+O"),
                            ("depth_cue", "Dept&h cueing"),
                            ("measure_mode", "&Measure mode\tM")):
            act = action(v, label, lambda on, a=attr: self._menu_toggle_view(a, on),
                         checkable=True)
            act.setChecked(bool(getattr(self.view, attr)))
            self.view_actions[attr] = act
        action(v, "Clear measurements", self.view.clear_measurements)
        menus[v].addSeparator()
        self.density_action = action(v, "Show 3D &density", lambda: None, checkable=True)
        self.density_action.setChecked(True)
        self.density_action.toggled.connect(self._toggle_density)
        self.density_action.setEnabled(False)
        action(v, "Clear 3D density", lambda: self.set_density(None, None))
        menus[v].addSeparator()
        for dock in (self.info_dock, self.tools_dock, self.table_dock, self.geometry_dock):
            menus[v].addAction(dock.toggleViewAction())
        self.geometry_dock.toggleViewAction().setShortcut(QtGui.QKeySequence("Ctrl+G"))

        action(t, "&Pack sorbates…", lambda: self._show_tools_tab(2))
        action(t, "&Force field and charges…", lambda: self._show_tools_tab(3))
        action(t, "Load &charges…", lambda: self.load_charges())
        action(t, "Generate &QEq charges", self.qeq_charges)
        menus[t].addSeparator()
        action(t, "&Input generator…", lambda: self.open_tool("input_generator"), "Ctrl+I")

        for key, label in (("info", "System &info"), ("check", "Chec&k structure"),
                           ("symmetry", "Find s&ymmetry"),
                           ("void_vol", "&Void volume"),
                           ("surface_area", "&Surface area"),
                           ("pore_size", "&Pore size distribution"), (None, None),
                           ("rdf", "&Radial distribution"),
                           ("coordination", "&Coordination number"),
                           ("hbonds", "&Hydrogen bonds"), ("msd", "&Mean squared displacement"),
                           ("rmsd", "RMSD vs &time"), ("gyration", "Radius of &gyration"),
                           ("density", "&Density profile"), ("density3d", "&3D density map"),
                           ("pxrd", "Powder &XRD"), (None, None),
                           ("energy_plot", "&Energy plot (energy.dat)…"),
                           ("isotherm_plot", "Plot is&otherm results…")):
            if key is None:
                menus[a].addSeparator()
            else:
                action(a, label, lambda _=False, k=key: self.open_tool(k))
        action("&Help", "&Keyboard and mouse", self.show_controls, "F1")
        action("&Help", "&About PDB Wizard", self.about)

        # the terminal viewer's trajectory keys
        for keys, slot in (("]", lambda: self.frame_bar.go(self.frame_bar.index + 1)),
                           ("[", lambda: self.frame_bar.go(self.frame_bar.index - 1)),
                           ("Home", lambda: self.frame_bar.go(0)),
                           ("End", lambda: self.frame_bar.go(self.frame_bar.n - 1)),
                           ("Space", self.frame_bar.toggle_play)):
            # for the 3D view and the frame bar only: window-wide they would take Space,
            # Home and End away from the buttons, check boxes and tables in the docks
            for owner in (self.view, self.frame_bar):
                sc = QtGui.QShortcut(QtGui.QKeySequence(keys), owner)
                sc.setContext(QtCore.Qt.WidgetWithChildrenShortcut)
                sc.activated.connect(slot)

    # ---- documents ----------------------------------------------------------------------

    def _store_view(self) -> None:
        d = self.doc
        d.view_state = self.view.view_state()
        d.density = self.view.density
        d.measurements = list(self.view.measurements)

    def _switch_doc(self, index: int) -> None:
        if index < 0 or index >= len(self.docs) or index == self.doc_index:
            return
        self._store_view()
        self.frame_bar.stop(notify=False)
        self.doc_index = index
        self._restore_doc()

    def _restore_doc(self) -> None:
        """Show the current document with its saved camera, measurements and density."""
        d = self.doc
        self.view.set_molecule(d.mol, reset_view=d.view_state is None, redetect=False)
        if d.view_state is not None:
            self.view.set_view_state(d.view_state)
        self.view.measurements = list(d.measurements)
        self.view.set_density(*(d.density if d.density else (None, None)))
        self.density_action.setEnabled(d.density is not None)
        self.frame_bar.set_count(len(d.frames) if d.frames else 1, d.frame_idx)
        self._refresh_panels()
        self._update_history_actions()

    def new_tab(self, mol: Optional[Molecule], path: Optional[str] = None,
                name: str = "", frames: Optional[list] = None) -> GuiDocument:
        """Open ``mol`` in a tab: the current one if it is empty, else a new one."""
        d = GuiDocument(mol, path, frames, name)
        if self.doc.empty and not self.doc.loading:
            self.docs[self.doc_index] = d
            self.tab_bar.setTabText(self.doc_index, d.name)
            self._show_doc()
        else:
            self._store_view()
            self.docs.append(d)
            self.tab_bar.addTab(d.name)
            self.tab_bar.setCurrentIndex(len(self.docs) - 1)  # -> _switch_doc
        return d

    def _show_doc(self) -> None:
        d = self.doc
        self.view.set_molecule(d.mol, reset_view=True, redetect=False)
        self.view.set_density(None, None)
        self.density_action.setEnabled(False)
        self.frame_bar.set_count(len(d.frames) if d.frames else 1, d.frame_idx)
        self._refresh_panels()
        self._update_history_actions()

    def close_tab(self, index: int) -> None:
        if not 0 <= index < len(self.docs):
            return
        self.frame_bar.stop(notify=index != self.doc_index)  # the shown tab stays: refresh it
        self._close_windows_of(self.docs[index])
        if len(self.docs) == 1:
            self.docs[0] = GuiDocument()
            self.tab_bar.setTabText(0, "untitled")
            self._show_doc()
            return
        current = self.doc_index
        if index != current:
            self._store_view()
        del self.docs[index]
        target = current - 1 if index < current else min(current, len(self.docs) - 1)
        self.tab_bar.blockSignals(True)
        self.tab_bar.removeTab(index)
        self.tab_bar.setCurrentIndex(target)
        self.tab_bar.blockSignals(False)
        self.doc_index = target
        if index == current:
            self._restore_doc()

    # ---- state --------------------------------------------------------------------------

    def set_molecule(self, mol: Optional[Molecule], reset_view: bool = True) -> None:
        self.mol = mol
        if self.frames and mol is not None:
            self.frames[self.doc.frame_idx] = mol
        self.view.set_molecule(mol, reset_view=reset_view)
        self._refresh_panels()

    def _refresh_panels(self) -> None:
        self.table_model.set_molecule(self.mol)
        self.geometry.set_molecule(self.mol)
        self._refresh_info()
        self._update_title()
        self._update_supercell_preview()
        for d in self._live_windows():  # tool windows work on their tab's current atoms
            if getattr(d, "doc", None) is self.doc and hasattr(d, "set_structure"):
                d.set_structure(self.doc.mol, self.doc.frames)

    def _changed(self, message: str) -> None:
        """After an edit: redraw everything and say what happened."""
        self.view.clear_marks()
        self._selected(None)
        self.view.refresh_bonds()
        self.view.update()
        self._refresh_panels()
        self.statusBar().showMessage(message, 8000)

    def _update_history_actions(self) -> None:
        """Undo and Redo follow the shown tab's history."""
        self.undo_action.setEnabled(bool(self._undo))
        self.redo_action.setEnabled(bool(self._redo))

    def _checkpoint(self) -> None:
        # remember which frame was edited: undo must restore that frame, not the shown one
        self._undo.append((self.doc.frame_idx, copy.deepcopy(self.mol)))
        dropped = self._undo[:-20]
        del self._undo[:-20]
        # a new edit ends the redo history; it is kept aside (with the undo steps the cap
        # just dropped) until the edit is known to change something, since an edit that
        # finds nothing to do calls _discard_checkpoint
        self._redo_stash = (self.doc, self._redo, dropped)
        self._redo = []
        self._update_history_actions()

    def _checkpoint_cells(self) -> None:
        """_checkpoint for an edit of the cell, which a trajectory applies to every frame:
        the undo step also holds each frame's cell, so one undo restores them all."""
        self._checkpoint()
        if self.frames:
            self._undo[-1] += ([copy.deepcopy(f.pbc) for f in self.frames],)

    def _discard_checkpoint(self) -> None:
        if self._undo:
            self._undo.pop()
        stash, self._redo_stash = self._redo_stash, None
        if stash is not None and stash[0] is self.doc:
            self._redo = stash[1]
            self._undo[:0] = stash[2]
        self._update_history_actions()

    def _step_history(self, source: list, target: list, message: str) -> None:
        """Undo (from the undo list onto the redo list) or redo (the other way round):
        show the saved state on the frame it was taken from, and keep the state it replaces
        in ``target`` so the step can be taken back."""
        entry = source.pop()
        frame, mol = entry[0], entry[1]
        cells = entry[2] if len(entry) > 2 else None  # a cell edit: every frame's cell
        if self.frames and frame < len(self.frames):
            self.frame_bar.stop(notify=False)
            self.doc.frame_idx = frame
            self.frame_bar.set_count(len(self.frames), frame)
        frame = self.doc.frame_idx
        back = (frame, self.frames[frame] if self.frames else self.mol)
        if cells is not None and self.frames and len(cells) == len(self.frames):
            back += ([copy.deepcopy(f.pbc) for f in self.frames],)
            for k, pbc in enumerate(cells):
                if k != frame:  # that frame is replaced below by ``mol``, saved with its cell
                    self.frames[k].pbc = copy.deepcopy(pbc)
        target.append(back)
        del target[:-20]
        self._redo_stash = None
        self.set_molecule(mol, reset_view=False)
        self.view.clear_marks()
        self._update_history_actions()
        self.statusBar().showMessage(message, 4000)

    def undo(self) -> None:
        if self._undo:
            self._step_history(self._undo, self._redo, "Undone")

    def redo(self) -> None:
        if self._redo:
            self._step_history(self._redo, self._undo, "Redone")

    def _need_structure(self, what: str, cell: bool = False) -> bool:
        if self.mol is None or not self.mol.atoms:
            self._warn(f"{what} needs a structure. Open one first.")
            return False
        if cell and not _gui_real_cell(self.mol):
            none = ("This structure has none" if self.mol.pbc is None else
                    "This trajectory has none (its huge box is only a stand-in)")
            self._warn(f"{what} needs a unit cell (PBC). {none} — "
                       "set one with Edit > Update unit cell first.")
            return False
        return True

    def _set_view(self, attr: str, value, invalidate: bool = False) -> None:
        setattr(self.view, attr, value)
        if invalidate:
            self.view.invalidate()
        else:
            self.view.update()

    def _toggle_view(self, attr: str, on: bool) -> None:
        if attr == "measure_mode":
            self.view.set_measure_mode(on)
            if on:
                self.statusBar().showMessage(
                    "Measure mode: click two atoms to pin their distance", 6000)
        else:
            self._set_view(attr, on, invalidate=True)

    def _menu_toggle_view(self, attr: str, on: bool) -> None:
        self._toggle_view(attr, on)
        self._sync_view_controls()

    def set_style(self, name: str) -> None:
        self.view.set_style(name)
        self._sync_view_controls()
        if name == "ribbon" and self.mol is not None and self.mol.atoms:
            if not _ca_chains(self.mol):
                self.statusBar().showMessage(
                    "No protein backbone (CA atoms) here; showing ball and stick", 6000)

    def _key_option(self, fn: Callable) -> None:
        fn()
        self._sync_view_controls()

    def _sync_view_controls(self) -> None:
        """Make the View tab match the view after a key or menu changed it."""
        v = self.view
        for combo, value in ((self.style_combo, v.style_name), (self.color_combo, v.color_mode),
                             (self.bg_combo, v.background)):
            combo.blockSignals(True)
            combo.setCurrentIndex(combo.findData(value))
            combo.blockSignals(False)
        for attr, chk in self.view_checks.items():
            for box in (chk, getattr(self, "view_actions", {}).get(attr)):
                if box is not None:
                    box.blockSignals(True)
                    box.setChecked(bool(getattr(v, attr)))
                    box.blockSignals(False)
        s = self.sliders["bond_radius"]
        s.blockSignals(True)
        s.setValue(int(round(v.bond_radius * 100)))
        s.blockSignals(False)

    def _update_title(self) -> None:
        name = self.doc.name if not self.doc.empty or self.doc.path else "no structure"
        frames = self.frames
        extra = (f" (frame {self.doc.frame_idx + 1}/{len(frames)})"
                 if frames and len(frames) > 1 else "")
        self.setWindowTitle(f"PDB Wizard — {name}{extra}")
        self.tab_bar.setTabText(self.doc_index, self.doc.name)

    def _refresh_info(self) -> None:
        labels = self.info_labels
        self.legend.clear()
        if self.mol is None:
            for label in labels.values():
                label.setText("—")
            return
        info = system_info(self.mol, self.path or "")
        labels["File"].setText(Path(self.path).name if self.path else f"({self.doc.name})")
        counts: dict = {}
        for a in self.mol.atoms:
            if not _is_mpmc_dummy_site(a):
                counts[a.element.symbol] = counts.get(a.element.symbol, 0) + 1
        labels["Formula"].setText(" ".join(f"{el}<sub>{n}</sub>" for el, n in counts.items()))
        labels["Atoms"].setText(str(info["n_atoms"]))
        labels["Sorbates"].setText(str(mpmc_movable_count(self.mol)))
        labels["Net charge (e)"].setText(f"{sum(a.charge for a in self.mol.atoms):+.4f}")
        frames = self.frames
        labels["Frame"].setText(f"{self.doc.frame_idx + 1} / {len(frames)}"
                                if frames and len(frames) > 1 else "—")
        if "cell" in info:
            c = info["cell"]
            labels["a, b, c (Å)"].setText(f"{c['a']:.3f}, {c['b']:.3f}, {c['c']:.3f}")
            labels["α, β, γ (°)"].setText(
                f"{c['alpha']:.2f}, {c['beta']:.2f}, {c['gamma']:.2f}")
            labels["Volume (Å³)"].setText(f"{info['volume']:.2f}")
            labels["Density (g/cm³)"].setText(f"{info['density']:.4f}")
        else:
            for key in ("a, b, c (Å)", "α, β, γ (°)", "Volume (Å³)", "Density (g/cm³)"):
                labels[key].setText("no cell")
        counts = {}
        for a in self.mol.atoms:
            key = "model sites" if _is_mpmc_dummy_site(a) else a.element.symbol
            counts[key] = counts.get(key, 0) + 1
        for key, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            pix = QtGui.QPixmap(14, 14)
            pix.fill(QtCore.Qt.transparent)
            p = QtGui.QPainter(pix)
            p.setRenderHint(QtGui.QPainter.Antialiasing)
            color = (QtGui.QColor(255, 20, 147) if key == "model sites"
                     else QtGui.QColor(*next(a for a in self.mol.atoms
                                             if a.element.symbol == key).element.cpk_color))
            p.setBrush(color)
            p.setPen(QtGui.QPen(color.darker(160)))
            p.drawEllipse(1, 1, 12, 12)
            p.end()
            self.legend.addItem(QtWidgets.QListWidgetItem(QtGui.QIcon(pix), f"{key}   {n}"))

    def _hovered(self, idx: Optional[int]) -> None:
        if idx is None or self.mol is None:
            self.statusBar().clearMessage()
            return
        self.statusBar().showMessage(self._describe(idx))

    def _selected(self, idx: Optional[int], from_table: bool = False) -> None:
        self.view.selected = idx
        self.view.update()
        if idx is None or self.mol is None:
            self.selection_label.setText("Click an atom to inspect it.")
            return
        self.selection_label.setText(self._describe(idx))
        if not from_table:
            self.table.selectRow(idx)
            self.table.scrollTo(self.table_model.index(idx, 0))

    def _highlight(self, indices: tuple) -> None:
        self.view.highlighted = set(indices)
        self.view.update()

    def _describe(self, idx: int) -> str:
        a = self.mol.atoms[idx]
        where = f"sorbate {a.mpmc_mol}" if a.mpmc_flag == "M" else "framework"
        return (f"#{idx + 1} {a.name} ({a.element.symbol}, {where})  "
                f"x={a.x[0]:.3f} y={a.x[1]:.3f} z={a.x[2]:.3f}  q={a.charge:+.4f} e  "
                f"ε={a.epsilon:.2f} K  σ={a.sigma:.3f} Å")

    def _warn(self, text: str) -> None:
        QtWidgets.QMessageBox.warning(self, "PDB Wizard", text)

    def _live_windows(self) -> list:
        self.dialogs = [d for d in self.dialogs if _qt_alive(d)]
        return list(self.dialogs)

    def _show_window(self, dialog, doc: Optional[GuiDocument] = None) -> object:
        """Show a tool window without blocking. It is freed when closed; ``doc`` ties it to
        a tab (it follows that tab's edits and closes with it)."""
        self._live_windows()
        dialog.doc = doc
        self.dialogs.append(dialog)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
        dialog.show()
        return dialog

    def _close_windows_of(self, doc: GuiDocument) -> None:
        for d in self._live_windows():
            if getattr(d, "doc", None) is doc:
                d.close()

    def report(self, title: str, text: str) -> object:
        return self._show_window(TextReport(title, text, self))

    def _busy(self, message: str, fn: Callable, done: Callable) -> GuiTask:
        """Run ``fn`` in the background with a wait cursor and a status message."""
        self.statusBar().showMessage(message)
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)

        def finish(result):
            QtWidgets.QApplication.restoreOverrideCursor()
            self.statusBar().clearMessage()
            done(result)

        def failed(e):
            QtWidgets.QApplication.restoreOverrideCursor()
            self.statusBar().clearMessage()
            self._warn(f"{message.rstrip('.…')} failed:\n{e}")

        self.task = GuiTask(fn, finish, failed)
        return self.task

    def _snapshot(self) -> Optional[Molecule]:
        """A copy of the shown structure for background work (edits may change the live one)."""
        return copy.deepcopy(self.mol) if self.mol is not None else None

    # ---- files --------------------------------------------------------------------------

    def open_dialog(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open structure", "", _OPEN_FILTER)
        if path:
            self.open_path(path)

    def open_path(self, path: str) -> bool:
        """Open a structure (or a trajectory, loaded in the background) in a tab."""
        if not os.path.isfile(path):
            self._warn(f"No such file: {path}")
            return False
        try:
            trajectory = is_trajectory_file(path)
        except Exception:  # an unreadable file: let read_file report it below
            trajectory = False
        if trajectory:
            return self._open_trajectory(path)
        messages: list = []
        try:
            mol = read_file(path, on_message=messages.append)
        except Exception as e:  # any reader error goes to the user, not the console
            self._warn(f"Could not read {Path(path).name}:\n{e}")
            return False
        if not mol.atoms:
            self._warn(f"No atoms found in {Path(path).name}")
            return False
        self.new_tab(mol, path)
        n = len(messages)
        note = (f" (deleted {n} overlapping atom{'s' if n != 1 else ''}, closer than "
                f"{OVERLAP_CUTOFF} Å)" if n else "")
        self.statusBar().showMessage(f"Opened {Path(path).name}: {len(mol.atoms)} atoms{note}",
                                     8000)
        self._warn_load_notes(mol, Path(path).name)
        return True

    def _warn_load_notes(self, mol: Molecule, name: str) -> None:
        """Say what the file reader could not honour (e.g. a CIF read as its asymmetric
        unit). Analysis > Check structure repeats it."""
        notes = load_notes(mol)
        if notes:
            self._warn(f"{name}:\n\n" + "\n\n".join(notes))

    def _open_trajectory(self, path: str) -> bool:
        d = self.new_tab(Molecule(atoms=[]), path, name=f"{Path(path).name} (loading)")
        d.loading = True
        progress = [0.0]

        def read():
            frames = read_trajectory(path, lambda f: progress.__setitem__(0, f))
            if frames:
                frames[0].detect_bonds()
            return frames

        def finished() -> bool:
            timer.stop()
            d.loading = False
            self.statusBar().clearMessage()
            return d in self.docs  # False if its tab was closed meanwhile

        def done(frames):
            if not finished():
                return
            d.frames, d.frame_idx = frames, 0
            d.mol = frames[0] if frames else Molecule(atoms=[])
            d.name = Path(path).name
            d.view_state = None  # fit the camera to the frames, not the empty placeholder
            if d is self.doc:
                self._show_doc()
            self.tab_bar.setTabText(self.docs.index(d), d.name)
            self.statusBar().showMessage(f"Loaded {len(frames)} frames", 6000)
            for note in trajectory_caveats(frames):
                self._warn(note)

        def failed(e):
            if finished():
                self.close_tab(self.docs.index(d))  # no empty "(loading)" tab left behind
            self._warn(f"Could not read {Path(path).name}:\n{e}")

        timer = QtCore.QTimer(self)
        timer.setInterval(100)
        timer.timeout.connect(lambda: self.statusBar().showMessage(
            f"Loading trajectory... {int(progress[0] * 100)}%"))
        timer.start()
        self.load_task = GuiTask(read, done, failed)
        return True

    def set_frame(self, idx: int) -> None:
        frames = self.frames
        if not frames or len(frames) <= 1:
            return
        idx = max(0, min(idx, len(frames) - 1))
        playing = self.frame_bar.playing
        if idx == self.doc.frame_idx and not playing and self.mol is frames[idx]:
            self._refresh_panels()
            return
        self.doc.frame_idx = idx
        frame_bonds(frames, idx)
        self.mol = frames[idx]
        self.view.set_molecule(self.mol, reset_view=False, redetect=False)
        if not playing:  # tables are expensive; refresh them when playback stops
            self._refresh_panels()
        else:
            self._update_title()

    def frame_to_tab(self) -> Optional[GuiDocument]:
        """Open the current trajectory frame as a standalone structure in a new tab."""
        if not self.frames or self.mol is None or not self.mol.atoms:
            self.statusBar().showMessage("No trajectory frame to open", 4000)
            return None
        stem = Path(self.path).stem if self.path else "frame"
        name = f"{stem}_frame{self.doc.frame_idx + 1}"
        d = self.new_tab(copy.deepcopy(self.mol), None, name)
        self.statusBar().showMessage(f"Opened frame as {name}", 4000)
        return d

    def save_dialog(self) -> None:
        if not self._need_structure("Saving"):
            return
        start = str(Path(self.path).with_suffix(".pqr")) if self.path else "structure.pqr"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save structure", start, _SAVE_FILTER)
        if path and path.lower().endswith(".d12"):
            group = self._ask_space_group()
            if group is not None:
                self.save_crystal(path, group)
        elif path:
            self.save_path(path)

    def _ask_space_group(self) -> Optional[str]:
        """Which space group to write a CRYSTAL23 input in (None if cancelled)."""
        text, ok = QtWidgets.QInputDialog.getText(
            self, "CRYSTAL23 input",
            "Space group number (1 = P1, every atom).\n"
            "Analysis > Find symmetry lists the groups this structure fits.", text="1")
        return text if ok else None

    def save_crystal(self, path: str, group: str = "1") -> bool:
        """Write a CRYSTAL23 .d12: P1 for group 1, else that group's asymmetric unit."""
        if not self._need_structure("Writing a CRYSTAL23 input", cell=True):
            return False
        try:
            reduction = crystal_reduction_choice(self.mol, group)
            with open(path, "w") as fh:
                write_crystal_d12(self.mol, fh, reduction)
        except (OSError, ValueError) as e:
            self._warn(f"Could not write {Path(path).name}:\n{e}")
            return False
        what = "P1" if reduction is None else f"space group {reduction.number} " \
            f"({reduction.symbol}), {len(reduction.numbers)} atoms"
        self.statusBar().showMessage(f"Saved {path} in {what}", 8000)
        return True

    def save_path(self, path: str) -> bool:
        if (Path(path).suffix.lower() in _CELL_FORMATS and self.mol is not None
                and self.mol.pbc is not None and not _gui_real_cell(self.mol)):
            # the writers refuse a structure without a cell; the stand-in box is not one
            return self._need_structure(f"Saving {Path(path).name}", cell=True)
        try:
            save_structure(self.mol, path)
        except Exception as e:
            self._warn(f"Could not save {Path(path).name}:\n{e}")
            return False
        self.statusBar().showMessage(f"Saved {path}", 8000)
        return True

    def export_image(self, path: Optional[str] = None) -> Optional[str]:
        if not path:
            start = str(Path(self.path).with_suffix(".png")) if self.path else "structure.png"
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export image", start,
                                                            "PNG (*.png)")
        if not path:
            return None
        if not self.view.grab().save(path):
            self._warn(f"Could not write {path} (check the folder and the extension)")
            return None
        self.statusBar().showMessage(f"Saved image {path}", 6000)
        return path

    def export_gif(self, path: Optional[str] = None, n_frames: int = 36) -> Optional[str]:
        """A full turn about the vertical axis, 80 ms per frame (as the terminal viewer)."""
        if not self._need_structure("A GIF"):
            return None
        try:
            from PIL import Image
        except ImportError:
            self._warn("Pillow required: pip install Pillow")
            return None
        if not path:
            start = str(Path(self.path).with_suffix(".gif")) if self.path else "structure.gif"
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export rotation GIF", start,
                                                            "GIF (*.gif)")
        if not path:
            return None
        n_frames = max(1, int(n_frames))
        base = self.view._rot.copy()
        images = []
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            for k in range(n_frames):
                self.view._rot = _rotation(1, 2 * np.pi * k / n_frames) @ base
                buf = QtCore.QBuffer()
                buf.open(QtCore.QIODevice.WriteOnly)
                self.view.grab().save(buf, "PNG")
                images.append(Image.open(_io.BytesIO(bytes(buf.data()))).convert("RGB"))
        finally:
            self.view._rot = base
            self.view.update()
            QtWidgets.QApplication.restoreOverrideCursor()
        try:
            images[0].save(path, save_all=True, append_images=images[1:], duration=80, loop=0)
        except (OSError, ValueError) as e:
            self._warn(f"Could not write {path}:\n{e}")
            return None
        self.statusBar().showMessage(f"Saved {path} ({n_frames} frames)", 6000)
        return path

    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:  # noqa: N802
        for url in event.mimeData().urls():
            self.open_path(url.toLocalFile())

    def closeEvent(self, event) -> None:  # noqa: N802
        self.frame_bar.stop(notify=False)
        if self.task is not None and not self.task.finished:
            self.task.detach()  # its report must not pop up after the window has gone
            self.task = None
            QtWidgets.QApplication.restoreOverrideCursor()
        for d in self._live_windows():
            d.close()
        if self._db_cache and os.path.isdir(self._db_cache):
            shutil.rmtree(self._db_cache, ignore_errors=True)
        super().closeEvent(event)

    # ---- edits --------------------------------------------------------------------------

    def _show_tools_tab(self, index: int) -> None:
        self.tools_dock.show()
        self.tools_dock.raise_()
        self.tools.setCurrentIndex(index)

    def _show_supercell_tab(self) -> None:
        self._show_tools_tab(1)

    def build_supercell(self) -> None:
        if not self._need_structure("A supercell", cell=True):
            return
        counts = [spin.value() for spin in self.extend_spins]
        if all(n == 1 for n in counts):
            self.statusBar().showMessage("Set at least one axis above 1", 5000)
            return
        self._checkpoint()
        for axis, n in enumerate(counts):
            if n > 1:
                extend_axis(self.mol, axis, n - 1)  # extend_axis counts the extra copies
        self.view.reset_view()
        size = "×".join(str(n) for n in counts)
        for spin in self.extend_spins:
            spin.setValue(1)
        self._changed(f"Built a {size} supercell: {len(self.mol.atoms)} atoms")

    def _update_supercell_preview(self) -> None:
        counts = [spin.value() for spin in self.extend_spins]
        mol = self.mol
        if mol is None or not mol.atoms or not self.view.has_cell():
            self.supercell_preview.setText("Needs a structure with a unit cell")
            return
        n = len(mol.atoms) * counts[0] * counts[1] * counts[2]
        size = " × ".join(str(c) for c in counts)
        note = " (no change)" if n == len(mol.atoms) else ""
        self.supercell_preview.setText(f"Result: {size} cells, {n:,} atoms{note}")

    def wrap(self, forward: bool) -> None:
        if not self._need_structure("Wrapping atoms", cell=True):
            return
        self._checkpoint()
        wrap_atoms(self.mol, forward=forward)
        self._changed("Wrapped atoms into " + ("0–1" if forward else "−½–½"))

    def sort_atoms(self) -> None:
        if not self._need_structure("Sorting atoms into molecules", cell=True):
            return
        self._checkpoint()
        sort_system(self.mol, False)
        self._changed("Sorted atoms")

    def _confirm(self, title: str, text: str) -> bool:
        box = QtWidgets.QMessageBox
        return box.question(self, title, text, box.Yes | box.No, box.No) == box.Yes

    def delete_lone_atoms(self, ask: bool = False) -> int:
        """Delete atoms with no close neighbour. ``ask`` (the menu) shows what would go
        first, with a caution when that includes framework metals."""
        if not self._need_structure("Finding lone atoms", cell=True):
            return 0
        if ask:
            found = lone_atom_summary(self.mol)
            if not found:
                self.statusBar().showMessage("No lone atoms found", 5000)
                return 0
            caution = lone_atom_warning(self.mol)
            if not self._confirm("Delete lone atoms", f"This would delete: {found}."
                                 + (f"\n\n{caution}" if caution else "") + "\n\nDelete them?"):
                return 0
        self._checkpoint()
        _, deleted = remove_lone_atoms(self.mol)
        if not deleted:
            self._discard_checkpoint()
            self.statusBar().showMessage("No lone atoms found", 5000)
            return 0
        n = len(deleted)
        self._changed(f"Deleted {n} lone atom{'s' if n != 1 else ''}")
        return n

    def edit_h_distances(self, element: Optional[str] = None,
                         distance: Optional[float] = None) -> int:
        if not self._need_structure("Editing H distances", cell=True):
            return 0
        if element is None or distance is None:
            others = [e for e in elements_of(self.mol) if e != "H"] or ["O"]
            element, ok = QtWidgets.QInputDialog.getItem(
                self, "Edit H distances", "Set every bonded X–H distance where X is:", others,
                others.index("O") if "O" in others else 0, False)
            if not ok:
                return 0
            distance, ok = QtWidgets.QInputDialog.getDouble(
                self, "Edit H distances", f"{element}–H distance (Å):", 1.0, 0.5, 3.0, 3)
            if not ok:
                return 0
        element = element.strip().capitalize()
        self._checkpoint()
        _, messages = adjust_h_distances(self.mol, element, distance)
        if not messages:
            self._discard_checkpoint()
            self.statusBar().showMessage(f"No bonded {element}-H pairs found; nothing changed",
                                         6000)
            return 0
        n = len(messages)
        self._changed(f"Set {n} {element}-H distance{'s' if n != 1 else ''} to {distance} A")
        return n

    def update_cell(self, values="ask") -> bool:
        """Set the cell to (a, b, c, alpha, beta, gamma), or remove it with None."""
        if not self._need_structure("Updating the cell"):
            return False
        if values == "ask":
            dlg = CellDialog(self.mol.pbc, self)
            if dlg.exec() != QtWidgets.QDialog.Accepted:
                return False
            values = dlg.values()
        if values is None:
            if self.mol.pbc is None:
                self.statusBar().showMessage("This system already has no PBC", 4000)
                return False
            self._checkpoint_cells()
            for f in self.frames or [self.mol]:
                f.pbc = None
            self._changed("PBC removed")
            return True
        a, b, c, alpha, beta, gamma = (float(v) for v in values)
        problem = _cell_problem(a, b, c, alpha, beta, gamma)
        if problem:
            self._warn(problem)
            return False
        self._checkpoint_cells()
        for f in self.frames or [self.mol]:  # a trajectory: every frame, not just this one
            if f.pbc is not None:
                f.pbc.update(a, b, c, alpha, beta, gamma)
            else:
                f.pbc = PBC(a, b, c, alpha, beta, gamma)
        self._changed(f"Cell: {a:.3f} {b:.3f} {c:.3f} {alpha:.1f} {beta:.1f} {gamma:.1f}")
        return True

    def substitute_element(self, old: Optional[str] = None, new: Optional[str] = None) -> int:
        if not self._need_structure("Substituting"):
            return 0
        if old is None or new is None:
            present = elements_of(self.mol)
            old, ok = QtWidgets.QInputDialog.getItem(self, "Substitute element",
                                                     "Replace every:", present, 0, False)
            if not ok:
                return 0
            new, ok = QtWidgets.QInputDialog.getText(self, "Substitute element",
                                                     f"Replace {old} with (element symbol):")
            if not ok or not new.strip():
                return 0
        new_elem = get_element(new)
        if new_elem.symbol == "X":
            self._warn(f"Unknown element: {new}")
            return 0
        self._checkpoint()
        count = 0
        for atom in self.mol.atoms:
            if atom.element.symbol == old:
                atom._element = new_elem
                atom.name = new_elem.symbol
                count += 1
        if count == 0:
            self._discard_checkpoint()
            self.statusBar().showMessage(f"No {old} atoms found", 4000)
            return 0
        self._changed(f"Substituted {count} {old} -> {new_elem.symbol}")
        return count

    def reduce_supercell(self, dialog: Optional[ReduceCellDialog] = None) -> bool:
        """Reduce Supercell. Pass an already-filled dialog to skip showing it (tests)."""
        if not self._need_structure("Supercell reduction", cell=True):
            return False
        if dialog is None:
            dialog = ReduceCellDialog(self.mol, self)
            if dialog.exec() != QtWidgets.QDialog.Accepted:
                return False
        result = dialog.result_mol
        if result is None:
            return False
        kept, total = len(dialog.keep()), len(dialog.sorbates)
        self._checkpoint()
        self.set_molecule(result, reset_view=True)
        self.statusBar().showMessage(
            f"Reduced to primitive cell: {len(result.atoms)} atoms, "
            f"kept {kept}/{total} sorbate molecules", 8000)
        return True

    def edit_coord(self, idx: int, xyz: tuple) -> None:
        if self.mol is None or not 0 <= idx < len(self.mol.atoms):
            return
        if len(xyz) != 3 or not np.all(np.isfinite(xyz)):
            self._warn("Coordinates must be three finite numbers")
            return
        self._checkpoint()
        self.mol.atoms[idx].x = np.array(xyz, dtype=float)
        self._changed(f"Updated atom {idx + 1} position")

    def pack_sorbates(self, model_path: Optional[str] = None) -> int:
        if not self._need_structure("Packing sorbates", cell=True):
            return 0
        key = self.sorbate_combo.currentData()
        if key == "__file__":
            if not model_path:
                model_path, _ = QtWidgets.QFileDialog.getOpenFileName(
                    self, "Sorbate template (MPMC .pqr)", "", "MPMC (*.pqr *.pdb)")
            if not model_path:
                return 0
            try:
                model = read_sorbate_template(model_path)
            except (OSError, ValueError) as e:
                self._warn(str(e))
                return 0
        else:
            model = SORBATE_MODELS[key]
        n = None if self.sorbate_fill.isChecked() else self.sorbate_count.value()
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            placed = insert_sorbates(self.mol, model, n, scale=self.sorbate_scale.value())
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        if not placed:
            self._warn("No room for even one molecule with these settings.")
            return 0
        self._checkpoint()
        first_id = next_mpmc_mol_id(self.mol)
        self.mol.atoms.extend(sorbate_atoms(model, placed, first_id))
        full = "" if n is None or len(placed) == n else f" (asked for {n}; the cell is full)"
        self._changed(f"Packed {len(placed)} {model.label}{full}. Save as .pqr for MPMC.")
        return len(placed)

    def apply_force_field(self) -> None:
        if not self._need_structure("A force field"):
            return
        idx = self.ff_combo.currentData()
        if idx is None:
            return
        self._checkpoint()
        _, missing = apply_ff_to_system_report(self.mol.atoms, get_forcefield(idx))
        self._changed(f"Applied {self.ff_combo.currentText()} to the framework")
        if missing:
            self._warn("\n".join(m.strip("! ") for m in missing_ff_messages(missing)))

    def load_charges(self, path: Optional[str] = None) -> bool:
        if not self._need_structure("Charges"):
            return False
        if not path:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Charges (.resp or one per line)")
        if not path:
            return False
        try:
            charges = read_charges_file(path)
        except (OSError, ValueError) as e:
            self._warn(f"Could not read charges:\n{e}")
            return False
        self._checkpoint()
        problem = assign_framework_charges(self.mol, charges)
        if problem:
            self._discard_checkpoint()
            self._warn(problem)
            return False
        self._changed(f"Loaded {len(charges)} charges onto the framework")
        return True

    def qeq_charges(self) -> None:
        if not self._need_structure("QEq charges"):
            return
        self._checkpoint()
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            apply_qeq_charges(self.mol)
        except Exception as e:  # e.g. an element without QEq parameters
            self._discard_checkpoint()
            self._warn(f"QEq failed: {e}")
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._changed("QEq charges assigned to the framework")
        self.report("QEq Charges", charge_report(self.mol))

    # ---- analysis and tools -------------------------------------------------------------

    def open_tool(self, key: str):
        """Open an Analysis-menu entry. Returns its window (or the background task)."""
        why = unavailable(key, self.mol, self.frames)
        if why:
            self.statusBar().showMessage(why, 6000)
            self._warn(why)
            return None
        mol, frames, doc = self.mol, self.frames, self.doc
        if key == "info":
            return self.report("System Info", system_info_text(mol, self.path or ""))
        if key == "check":
            return self.report("Check Structure", check_text(mol, frames))
        if key == "symmetry":
            return self.report("Find Symmetry", symmetry_text(mol))
        if key in ("void_vol", "surface_area"):
            mol = self._snapshot()
        if key == "void_vol":
            return self._busy("Calculating void volume (100k samples)...",
                              lambda: void_volume_text(mol),
                              lambda text: self.report("Void Volume", text))
        if key == "surface_area":
            return self._busy("Calculating surface area (N2 probe)...",
                              lambda: surface_area_text(mol),
                              lambda text: self.report("Surface Area", text))
        if key == "density3d":
            return self.density_3d()
        if key == "energy_plot":
            return self.energy_plot()
        if key == "input_generator":
            return self._show_window(InputGeneratorDialog(mol, self.path or "", self), doc)
        if key == "db_search":
            if self._db_cache is None:
                self._db_cache = tempfile.mkdtemp(prefix="pdb_wizard_db_")
            return self._show_window(DatabaseDialog(self._db_cache, self._open_fetched, self))
        dialog = AnalysisDialog(key, mol, frames, self)
        self._show_window(dialog, doc if key != "isotherm_plot" else None)
        if not dialog.spec.fields:  # nothing to choose: compute straight away
            dialog.run()
        return dialog

    def _open_fetched(self, path: str, label: str, rid: str) -> bool:
        """Open a downloaded file in a tab. False (after a warning) if it is not a structure."""
        messages: list = []
        try:
            mol = read_file(path, True, messages.append)
        except Exception as e:  # a download that is not a readable structure
            self._warn(f"Could not read {label} {rid}:\n{e}")
            return False
        if not mol.atoms:
            self._warn(f"No atoms found in {label} {rid}")
            return False
        # no path: the download sits in a temp folder that is deleted when the window
        # closes, so nothing (Save as, the input generator) may default to that folder
        self.new_tab(mol, None, Path(path).name)
        self.statusBar().showMessage(f"Opened {label} {rid} ({len(mol.atoms)} atoms)", 6000)
        self._warn_load_notes(mol, f"{label} {rid}")
        return True

    def energy_plot(self, path: Optional[str] = None):
        """MPMC energy.dat next to the structure, or one you pick."""
        path = path or find_energy_file(self.path or "")
        if not path:
            start = str(Path(self.path).parent / "energy.dat") if self.path else "energy.dat"
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open energy file", start,
                                                            "Energy (*.dat *.csv);;All files (*)")
        if not path:
            return None
        try:
            data = read_energy_dat(path)
        except (OSError, ValueError) as e:  # a folder, not text, rows of unequal length
            self._warn(f"Could not read {path}:\n{e}")
            return None
        if not data:
            self._warn(f"Could not read {path}")
            return None
        return self._show_window(EnergyDialog(data, path, self))

    def density_3d(self, element: Optional[str] = None) -> Optional[GuiTask]:
        if element is None:
            items = elements_of(self.mol)
            element, ok = QtWidgets.QInputDialog.getItem(self, "3D density map", "Element:",
                                                         items, 0, False)
            if not ok:
                return None
        doc = self.doc
        frames = list(self.frames)
        current = next((k for k, f in enumerate(frames) if f is self.mol), -1)
        if current >= 0:  # the shown frame may be edited meanwhile; the rest cannot
            frames[current] = self._snapshot()

        def done(res):
            if "positions" not in res.extra:  # no atoms of that element
                self.statusBar().showMessage(res.summary, 6000)
                return
            if doc not in self.docs:
                return
            doc.density = (np.asarray(res.extra["positions"]),
                           np.asarray(res.extra["values"], float))
            if doc is self.doc:  # the map belongs to the tab it was computed for
                self.view.set_density(*doc.density)
                self.density_action.setEnabled(True)
            self.statusBar().showMessage(f"{res.summary} (View > Show 3D density)", 8000)

        return self._busy(f"Computing 3D density for {element}...",
                          lambda: run_density_3d(frames, element, 15), done)

    def set_density(self, positions, values) -> None:
        self.view.set_density(positions, values)
        self.doc.density = None if positions is None else self.view.density
        self.density_action.setEnabled(positions is not None)

    def _toggle_density(self, on: bool) -> None:
        self.view.show_density = on
        self.view.update()

    def show_controls(self):
        """Help > Keyboard and mouse: everything the 3D view responds to."""
        return self.report("Keyboard and mouse", CONTROLS_HELP)

    def about(self) -> None:
        QtWidgets.QMessageBox.about(
            self, "About PDB Wizard",
            "<h3>PDB Wizard</h3><p>Structure handling for porous-materials simulations.</p>"
            "<p>Space group, NC State University. GPL-3.0-or-later.</p>")


# ======================================================================
# Start-up
# ======================================================================

def _apply_theme(app) -> None:
    """Use the desktop's own style (Breeze on KDE); otherwise a dark Fusion theme."""
    desktop = os.environ.get("XDG_CURRENT_DESKTOP", "").lower()
    if "kde" in desktop or os.environ.get("KDE_FULL_SESSION"):
        return
    app.setStyle("Fusion")
    pal, c, P = QtGui.QPalette(), QtGui.QColor, QtGui.QPalette
    text, panel = c(225, 228, 232), c(48, 52, 60)
    colors = {
        P.Window: c(37, 40, 46), P.WindowText: text, P.Base: c(28, 31, 36),
        P.AlternateBase: c(42, 46, 53), P.Text: text, P.Button: panel, P.ButtonText: text,
        P.Highlight: c(61, 174, 233), P.HighlightedText: c(255, 255, 255),
        P.ToolTipBase: panel, P.ToolTipText: text, P.PlaceholderText: c(140, 146, 155),
    }
    for role, color in colors.items():
        pal.setColor(role, color)
    for role in (P.Text, P.ButtonText, P.WindowText):
        pal.setColor(P.Disabled, role, c(120, 124, 130))
    app.setPalette(pal)


def run_gui(path: Optional[str] = None) -> int:
    """Open the desktop window (blocking). Returns the exit code."""
    if not _HAS_QT:
        print(GUI_INSTALL_HINT)
        return 1
    if qt_loses_none_references():
        print(GUI_BROKEN_QT_HINT)
        return 1
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(["pdb_wizard"])
    app.setApplicationName("PDB Wizard")
    _apply_theme(app)
    window = MainWindow(path)
    window.show()
    return app.exec()

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


def _stdin_is_terminal() -> bool:
    """False when input is piped or redirected (a script is driving the program)."""
    try:
        return bool(sys.stdin.isatty())
    except (AttributeError, ValueError):  # stdin closed or replaced
        return False


def _ensure_installed(package: str, import_name: str | None = None) -> bool:
    """Try to import a package; if missing, offer to pip install it."""
    import_name = import_name or package
    try:
        __import__(import_name)
        return True
    except ImportError:
        pass
    if not _stdin_is_terminal():
        # Never ask a script: the question would swallow its first line of input.
        return False
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


def _read_or_exit(filepath: str):
    """read_file with a clean one-line error instead of a traceback."""
    try:
        mol = read_file(filepath)
    except (ValueError, OSError) as e:
        print(f"Error: could not read {filepath}: {e}", file=sys.stderr)
        sys.exit(1)
    for note in load_notes(mol):  # e.g. a CIF read as its asymmetric unit only
        print(f"Warning: {note}", file=sys.stderr)
    return mol


def _read_trajectory_or_exit(filepath: str):
    """read_file_trajectory with a clean one-line error instead of a traceback."""
    try:
        return read_file_trajectory(filepath)
    except (ValueError, OSError) as e:
        print(f"Error: could not read {filepath}: {e}", file=sys.stderr)
        sys.exit(1)


def _batch_element(symbol: str | None) -> str | None:
    """An element symbol as the structure stores it: 'o' and 'CU' mean O and Cu."""
    if symbol and symbol.strip().capitalize() in ELEMENTS:
        return symbol.strip().capitalize()
    return symbol


def _save_csv_or_exit(csv_path: str, columns, header: str) -> None:
    import numpy as np
    try:
        np.savetxt(csv_path, np.column_stack(columns), header=header, delimiter=",",
                   comments="")
    except OSError as e:
        print(f"Error: could not write {csv_path}: {e}", file=sys.stderr)
        sys.exit(1)


def _batch_info(filepath: str) -> None:
    """Print system info and exit."""

    _check_file(filepath)
    mol = _read_or_exit(filepath)
    print_info(mol, filepath)


def _warn_fallback_cell(frames) -> None:
    """Tell the user when a trajectory carried no cell (a huge box was assumed)."""
    if frames and is_fallback_cell(frames[0].pbc):
        print(
            f"Warning: trajectory has no cell information; assuming a "
            f"{NO_CELL_FALLBACK_LENGTH:.0f} A cubic box (non-periodic). "
            "Periodic quantities such as g(r) normalisation will be meaningless.",
            file=sys.stderr,
        )


def _batch_rdf(filepath: str, element1: str, element2: str, csv_path: str | None) -> None:
    """Compute RDF and optionally write CSV."""
    _check_file(filepath)

    element1, element2 = _batch_element(element1), _batch_element(element2)
    frames = _read_trajectory_or_exit(filepath)
    _warn_fallback_cell(frames)
    if frames is not None:
        r, g = compute_rdf_trajectory(frames, element1, element2)
    else:
        mol = _read_or_exit(filepath)
        r, g = compute_rdf(mol, element1, element2)

    if len(r) == 0:
        print(f"No RDF data for {element1}-{element2} (check elements or PBC).",
              file=sys.stderr)
        sys.exit(1)

    if csv_path:
        _save_csv_or_exit(csv_path, [r, g], "r(A),g(r)")
        print(f"RDF written to {csv_path}")
    else:
        for ri, gi in zip(r, g):
            print(f"{ri:.4f}  {gi:.6f}")


def _batch_msd(filepath: str, element: str | None, csv_path: str | None) -> None:
    """Compute MSD and optionally write CSV."""
    _check_file(filepath)

    element = _batch_element(element)
    frames = _read_trajectory_or_exit(filepath)
    if frames is None:
        print("MSD requires a trajectory file with multiple frames.",
              file=sys.stderr)
        sys.exit(1)
    _warn_fallback_cell(frames)

    lags, msd = compute_msd(frames, element=element)
    caution = mc_trajectory_warning(frames, "a diffusion constant")
    if caution:
        print(caution, file=sys.stderr)

    if len(lags) == 0:
        print("No MSD data (check element or trajectory).", file=sys.stderr)
        sys.exit(1)

    if csv_path:
        _save_csv_or_exit(csv_path, [lags, msd], "lag(frames),MSD(A^2)")
        print(f"MSD written to {csv_path}")
    else:
        for li, mi in zip(lags, msd):
            print(f"{li:.0f}  {mi:.6f}")


# Formats whose writers cannot work without a unit cell. A standard PDB can be
# written without CRYST1, so .pdb/.ent are not listed.
_CELL_REQUIRED_SUFFIXES = (".pqr", ".cif", ".vasp", ".poscar")


def _batch_convert(filepath: str, output: str) -> None:
    """Convert between molecular file formats.

    ``.pdb``/``.ent`` write a standardized PDB; ``.pqr`` writes an MPMC-format
    PDB (charges and force-field columns preserved); ``.xyz``, ``.cif``,
    ``.vasp``/``.poscar`` and LAMMPS data are also supported. Formats that
    need a unit cell refuse to write when the structure has none, instead of
    leaving an empty file behind.
    """
    _check_file(filepath)
    from pathlib import Path


    suffix = Path(output).suffix.lower()
    supported = (".xyz", ".pdb", ".ent", ".lmp", ".lammps", ".data") + _CELL_REQUIRED_SUFFIXES
    if suffix not in supported:
        print(f"Unsupported output format: {suffix}", file=sys.stderr)
        sys.exit(1)

    n_frames = _trajectory_frame_count(filepath)
    if n_frames > 1 and os.path.exists(output) and os.path.samefile(filepath, output):
        print(
            f"Error: {filepath} is a trajectory of {n_frames} frames and --convert writes only "
            "the first; writing it onto itself would lose the rest. Choose another output "
            "file.",
            file=sys.stderr,
        )
        sys.exit(1)

    mol = _read_or_exit(filepath)

    if suffix in _CELL_REQUIRED_SUFFIXES and mol.pbc is None:
        print(
            f"Error: {suffix} output requires a unit cell, but {filepath} has none. "
            "Convert to .xyz or add a cell to the input first.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not mol.atoms and suffix != ".xyz":  # .xyz has always written a 0-atom file
        print(f"Error: could not write {output}: the structure has no atoms", file=sys.stderr)
        sys.exit(1)

    existed = os.path.exists(output)
    if existed:
        print(f"Warning: overwriting existing file {output}", file=sys.stderr)

    try:
        if suffix == ".pqr":
            write_mpmc_pdb(mol, output)
        else:
            with open(output, "w") as out:
                if suffix == ".xyz":
                    write_xyz(mol, out)
                elif suffix in (".pdb", ".ent"):
                    write_standard_pdb(mol, out, skip_mols_step=True)
                elif suffix == ".cif":
                    write_cif(mol, out)
                elif suffix in (".vasp", ".poscar"):
                    write_poscar(mol, out)
                else:
                    write_lammps_data(mol, out)
    except (OSError, ValueError, IndexError) as e:
        # e.g. a missing folder, or a structure with no atoms: no half-written file is left
        if not existed and os.path.isfile(output):
            os.remove(output)
        reason = str(e) if isinstance(e, OSError) or mol.atoms else "the structure has no atoms"
        print(f"Error: could not write {output}: {reason}", file=sys.stderr)
        sys.exit(1)

    if n_frames > 1:
        print(
            f"Warning: {filepath} has {n_frames} frames; only the first was converted.",
            file=sys.stderr,
        )
    print(f"Converted {filepath} -> {output}")


def _trajectory_frame_count(filepath: str) -> int:
    """Frames in a multi-frame PDB or XYZ file (1 for anything else)."""
    try:
        kind = detect_filetype(filepath)
        if kind == "xyz" and check_xyz_trajectory(filepath):
            with open(filepath) as f:
                return max(len(read_xyz_trajectory(f)[0]), 1)
        if kind == "pdb" and check_pdb_trajectory(filepath):
            with open(filepath) as f:
                return max(len(read_pdb_trajectory(f)[0]), 1)
    except (OSError, ValueError, IndexError, SystemExit):
        pass
    return 1


def _batch_surface(filepath: str) -> None:
    """Compute and print surface area."""
    _check_file(filepath)

    mol = _read_or_exit(filepath)
    result = surface_area(mol)
    print(f"Surface area: {result['surface_area']:.2f} A^2")
    print(f"Area per volume: {result['area_per_volume']:.2f} m^2/g")


def _batch_psd(filepath: str) -> None:
    """Compute and print pore size distribution."""
    _check_file(filepath)

    mol = _read_or_exit(filepath)
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
    parser.add_argument(
        "--gui", action="store_true",
        help="open the desktop window (needs PySide6: pip install \"pdb-wizard[gui]\")",
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
        help="convert file to another format (requires -o for output path; "
             ".pqr output is MPMC format)",
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
    if args.gui:
        sys.exit(run_gui(args.file))
    elif args.classic:
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
            if not os.path.exists(args.file):
                parser.error(f"file not found: {args.file}")
            try:
                ft = detect_filetype(args.file)
            except ValueError as e:
                parser.error(str(e))
            if not _stdin_is_terminal():
                # A script is piping answers in (as into the original pdb_wizard.py): the
                # viewer cannot run on a pipe, so the classic menu reads them.
                print(
                    "Input is not a terminal — falling back to the classic text menu "
                    "(use --classic to select it directly).",
                    file=sys.stderr,
                )
                run_classic(args.file)
                return
            is_traj = False
            if ft == "pdb":
                is_traj = check_pdb_trajectory(args.file)
            elif ft == "xyz":
                is_traj = check_xyz_trajectory(args.file)
            elif ft == "dcd":
                is_traj = check_dcd_trajectory(args.file)
                if not is_traj:
                    parser.error(f"{args.file} is not a readable DCD trajectory")
            elif ft == "ase_traj":
                is_traj = check_ase_trajectory(args.file)
                if not is_traj:
                    parser.error(f"{args.file} is not a readable ASE trajectory")

            if is_traj:
                # Launch app empty, load trajectory in background
                empty = Molecule(atoms=[])
                app = PdbWizardApp(molecule=empty, filepath=args.file, load_trajectory=True)
            else:
                load_messages: list = []
                try:
                    molecule = read_file(args.file, on_message=load_messages.append)
                except (ValueError, OSError) as e:
                    parser.error(f"could not read {args.file}: {e}")
                app = PdbWizardApp(molecule=molecule, filepath=args.file,
                                   load_messages=load_messages)
        app.run()


if __name__ == '__main__':
    main()
