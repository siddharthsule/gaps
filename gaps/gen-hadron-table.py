#!/usr/bin/env python3

"""
Generate database.{h,cpp} (cpu-shower) and database.{cuh,cu} (gpu-shower): the
mass [GeV], proper decay length [mm] and decay modes of every state an event
can carry, keyed by signed pid.

The states come from the PDG numbering, with the name, mass and lifetime
`particle` and `pdg` give them; the modes from the EvtGen files `decaylanguage`
ships. A state with no mass, a mode no decay can follow, and a state left with
no way out are all dropped, as is a mode naming a state hadron_weight gives no
weight.

Run it with the py-shower-had venv, which has the three packages. The tables
are written beside this script, so it can be run from anywhere:

    ~/Code/py-shower-had/.venv/bin/python gaps/gen-hadron-table.py
"""

import collections
import functools
import os
import sys

try:
    import decaylanguage
    import pdg
    from decaylanguage import DecFileParser
    from particle import InvalidParticle, Particle, ParticleNotFound
    from particle.exceptions import MatchingIDNotFound
    from pdg.errors import PdgAmbiguousValueError, PdgNoDataError
except ImportError:
    sys.exit("needs 'decaylanguage', 'particle' and 'pdg': "
             "pip install decaylanguage particle pdg")

# The script sits in the source root, beside cpu-shower/ and gpu-shower/
SRC = os.path.dirname(os.path.abspath(__file__))

# The two trees written. The CUDA tables live on the device alone.
VARIANTS = (
    ("cpu-shower", "h", "cpp", "extern const", "inline"),
    ("gpu-shower", "cuh", "cu", "extern __device__ const",
     "__device__ inline"),
)

# The decay files in the order they are tried, and the proper decay length at
# or above which a state is a final state particle [mm]
DEC_NAMES = ("DECAY_LHCB.DEC", "DECAY_BELLE2.DEC")
CTAU_MAX = 1000.0

# hbar c [GeV mm] and the speed of light [mm / s]
HBARC = 1.973269804e-13
C = 2.99792458e11

# ctau sentinels, and the cap on children, above which one_to_n_decay mostly
# exhausts its trials
CTAU_STABLE = "1e30"
CTAU_NO_WIDTH = "-1."
MAX_CHILDREN = 10


def sixFigures(value):
    """Six significant figures, all the precision the PDG quotes and all the
    tables carry."""
    return float("%.6g" % value)

# ------------------------------------------------------------------------------
# The states the event can never hold, and so can never decay into


def isForbidden(pid):
    """True where hadron_weight gives the state no weight."""
    a = abs(pid)
    n1 = (a // 1000) % 10
    n2 = (a // 100) % 10
    n3 = (a // 10) % 10
    nJ = a % 10

    # Mesons, leptons and the photon: the pseudoscalar quarkonia alone
    if n1 == 0:
        return a in (441, 551, 100441, 100551)

    # Baryons
    return (a > 9999 or nJ > 4 or (nJ == 4 and n2 < n3)
            or n1 < n2 or n1 < n3)

# ------------------------------------------------------------------------------
# The ids the table is walked out from

# Mesons, by spin: pseudoscalar, vector, tensor, spin-3 and spin-4
spin0Mesons = [111,
               211, 221,
               311, 321, 331,
               411, 421, 431, 441,
               511, 521, 531, 541, 551]

spin1Mesons = [113,
               213, 223,
               313, 323, 333,
               413, 423, 433, 443,
               513, 523, 533, 543, 553]

spin2Mesons = [115,
               215, 225,
               315, 325, 335,
               415, 425, 435, 445,
               515, 525, 535, 545, 555]

spin3Mesons = [117,
               217, 227,
               317, 327, 337,
               417, 427, 437, 447,
               517, 527, 537, 547, 557]

spin4Mesons = [119,
               219, 229,
               319, 329, 339,
               419, 429, 439, 449,
               519, 529, 539, 549, 559]

# Baryons, spin 1/2 and spin 3/2. 1112, 2222, 3332, 4442 and 5552 are left out
# by Pauli, and a swapped id is only added when n1 > n2 > n3.
spinHalfBaryons = [2112,
                   2212,
                   3112,
                   3212, 3122, 3222,
                   3312, 3322,
                   4112,
                   4212, 4122, 4222,
                   4312, 4132, 4322, 4232, 4332,
                   4412, 4422, 4432,
                   5112,
                   5212, 5122, 5222,
                   5312, 5132, 5322, 5232, 5332,
                   5412, 5142, 5422, 5242, 5432, 5342, 5442,
                   5512, 5522, 5532, 5542]

spin3HalfBaryons = [1114,
                    2114,
                    2214, 2224,
                    3114,
                    3214, 3224,
                    3314, 3324, 3334,
                    4114,
                    4214, 4224,
                    4314, 4324, 4334,
                    4414, 4424, 4434, 4444,
                    5114,
                    5214, 5224,
                    5314, 5324, 5334,
                    5414, 5424, 5434, 5444,
                    5514, 5524, 5534, 5544, 5554]


def startIds():
    """Every meson under the orbital and radial excitation prefixes
    10000 nL + 100000 nR and the prefixes 9000000 + 10000 n of the states
    outside that scheme, and every baryon in its ground state."""
    start = []

    for fullList in (spin0Mesons, spin1Mesons, spin2Mesons, spin3Mesons,
                     spin4Mesons):
        for pid in fullList:
            for nL in (0, 1, 2, 3):
                for nR in (0, 1, 2, 3):
                    start.append(100000 * nR + 10000 * nL + pid)
            for n in range(10):
                start.append(9000000 + 10000 * n + pid)

    for fullList in (spinHalfBaryons, spin3HalfBaryons):
        start.extend(fullList)

    return start

# ------------------------------------------------------------------------------
# Names, masses and proper decay lengths. `particle` stores masses in MeV and
# decay lengths in mm, `pdg` masses and widths in GeV and lifetimes in s.

api = pdg.connect()

# Observed states `pdg` holds with no MC id, by the name it gives them
pdgNames = {
    417: "D_3^*(2750)+",
    427: "D_3^*(2750)0",
    437: "D_s3^*(2860)+",
    447: "psi_3(3842)",
    4422: "Xi_cc()++",
    5312: "Xi_b^'(5935)-",
    5314: "Xi_b(5955)-",
    5324: "Xi_b(5945)0",
}


@functools.lru_cache(maxsize=None)
def lookup(pid):
    """(name, mass [GeV], ctau [mm]) of a positive id, from `particle` first
    and `pdg` for what it leaves out, or None where neither holds the id. A
    mass or ctau neither gives comes back None."""
    found, name, mass, ctau = False, None, None, None

    try:
        p = Particle.from_pdgid(pid)
    except (ParticleNotFound, InvalidParticle):
        pass
    else:
        found, name, ctau = True, p.name, p.ctau
        if p.mass is not None:
            mass = p.mass / 1000.0
        elif pid in (12, 14, 16):
            mass = 0.0

    try:
        q = api.get_particle_by_mcid(pid)
    except ValueError:
        q = (api.get_particle_by_name(pdgNames[pid])
             if pid in pdgNames else None)

    if q is not None:
        found, name = True, name or q.name

        if mass is None:
            try:
                mass = q.mass
            except (PdgNoDataError, PdgAmbiguousValueError):
                pass

        if ctau is None:
            try:
                lifetime = q.lifetime if q.has_lifetime_entry else None
            except (PdgNoDataError, PdgAmbiguousValueError):
                lifetime = None
            try:
                width = q.width if q.has_width_entry else None
            except (PdgNoDataError, PdgAmbiguousValueError):
                width = None

            if lifetime:
                ctau = lifetime * C
            elif width:
                ctau = HBARC / width

    return (name, mass, ctau) if found else None


@functools.lru_cache(maxsize=None)
def selfConjugate(pid):
    """True for a state that is its own antiparticle: one `particle` holds
    with no separate antiparticle id, or one it does not hold that is a meson
    of a quark and its own antiquark."""
    try:
        Particle.from_pdgid(pid)
    except (ParticleNotFound, InvalidParticle):
        a = abs(pid)
        return (a // 1000) % 10 == 0 and (a // 100) % 10 == (a // 10) % 10
    try:
        Particle.from_pdgid(-pid)
    except (ParticleNotFound, InvalidParticle):
        return True
    return False


@functools.lru_cache(maxsize=None)
def nameOf(pid):
    """The name of a signed id, with a trailing ~ on an antiparticle
    `particle` does not hold."""
    try:
        return Particle.from_pdgid(pid).name
    except (ParticleNotFound, InvalidParticle):
        row = lookup(abs(pid))
        name = row[0] if row is not None else str(abs(pid))
        return name if pid > 0 else name + "~"

# ------------------------------------------------------------------------------
# The decay files, each parsed once with the names of the states it gives modes
# for

sources = []
for decName in DEC_NAMES:
    parser = DecFileParser(os.path.join(decaylanguage.__path__[0], "data",
                                        decName))
    parser.parse()
    sources.append((decName, parser,
                    frozenset(parser.list_decay_mother_names())))

# States `particle` gives no EvtGen name, by the names the decay files give
# the particle and its antiparticle
evtgenNames = {
    4422: ("Xi_cc++", "anti-Xi_cc--"),
    5212: ("Sigma_b0", "anti-Sigma_b0"),
    5214: ("Sigma_b*0", "anti-Sigma_b*0"),
    5312: ("Xi'_b-", "anti-Xi'_b+"),
    5314: ("Xi_b*-", "anti-Xi_b*+"),
    5324: ("Xi_b*0", "anti-Xi_b*0"),
    100551: ("eta_b(2S)", None),
}

# The id behind every EvtGen name a mode can hold
pidByName = {}
for entry in Particle.all():
    try:
        pidByName.setdefault(entry.evtgen_name, int(entry.pdgid))
    except MatchingIDNotFound:
        pass
for pid, (evtName, antiName) in evtgenNames.items():
    pidByName.setdefault(evtName, pid)
    if antiName is not None:
        pidByName.setdefault(antiName, -pid)

# The name behind every alias each decay file defines
aliases = {decName: parser.dict_aliases() for decName, parser, _ in sources}


def evtgenName(pid):
    """The EvtGen name of an id, or None where nothing gives one."""
    if pid in evtgenNames:
        return evtgenNames[pid][0]
    try:
        return Particle.from_pdgid(pid).evtgen_name
    except (ParticleNotFound, InvalidParticle, MatchingIDNotFound):
        return None


def broken(mode):
    """True for a mode no decay can follow: switched off, with no children,
    or with a child that has no id or is a quark, gluon or diquark."""
    bf, children = mode
    return bf <= 0 or not children or any(
        d is None or abs(d) <= 8 or abs(d) == 21
        or (1000 <= abs(d) <= 9999 and (abs(d) // 10) % 10 == 0)
        for d in children)


def decayModes(pid):
    """Every mode one decay file gives the state, as (bf, child ids), a child
    with no id held as None. The file is the first in DEC_NAMES giving the
    state an unbroken mode, or where none does, the first listing it at
    all."""
    evtName = evtgenName(pid)
    fallback = ()
    for decName, parser, parents in sources:
        if evtName not in parents:
            continue
        modes = []
        for mode in parser.build_decay_chains(
                evtName, stable_particles=parents)[evtName]:
            names = [aliases[decName].get(d, d) for d in mode["fs"]]
            modes.append((mode["bf"], tuple(pidByName.get(d) for d in names)))
        if not all(broken(mode) for mode in modes):
            return tuple(modes)
        fallback = fallback or tuple(modes)
    return fallback

# ------------------------------------------------------------------------------
# The table, walked out from the start ids


def isFinal(ctau):
    """True for a state living at or past CTAU_MAX."""
    return ctau is not None and ctau >= CTAU_MAX


def buildTable(cuts):
    """{pid: (name, mass, ctau, modes)} over positive ids. A state with no
    mass stays out, a final state carries no mode, a mode survives when its
    children have masses that leave it open, and every child a surviving mode
    reaches joins in turn."""
    table = {}
    queue = sorted(set(startIds()), reverse=True)
    seen = set(queue)

    while queue:
        pid = queue.pop()
        row = lookup(pid)
        if row is None:
            cuts["start ids neither package holds"] += 1
            continue
        name, mass, ctau = row
        if mass is None:
            cuts["states with no mass"] += 1
            continue

        kept = []
        for mode in (() if isFinal(ctau) else decayModes(pid)):
            if broken(mode):
                cuts["modes broken"] += 1
                continue

            bf, children = mode
            childRows = [lookup(abs(d)) for d in children]
            if any(r is None or r[1] is None for r in childRows):
                cuts["modes with a child that has no mass"] += 1
                continue
            if sum(r[1] for r in childRows) > mass:
                cuts["modes closed at the PDG masses"] += 1
                continue

            kept.append(mode)

        table[pid] = (name, mass, ctau, tuple(kept))
        for child in (abs(d) for _, children in kept for d in children):
            if child not in seen:
                seen.add(child)
                queue.append(child)

    return table


def pruneDeadEnds(table, cuts):
    """Drop every mode reaching a short-lived state with no mode, which would
    sit in the event forever, and the states that leaves empty. Dropping a
    mode can empty another state, so this runs to a fixed point."""
    changed = True
    while changed:
        dead = set(p for p, row in table.items()
                   if not row[3] and not isFinal(row[2]))
        changed = False
        for p, (name, mass, ctau, modes) in list(table.items()):
            kept = tuple(mode for mode in modes
                         if not any(abs(d) in dead for d in mode[1]))
            if len(kept) == len(modes):
                continue
            cuts["modes reaching a state that cannot decay"] += (len(modes)
                                                                 - len(kept))
            table[p] = (name, mass, ctau, kept)
            changed = True

    for p in [p for p, row in table.items()
              if not row[3] and not isFinal(row[2])]:
        cuts["states left with no mode"] += 1
        del table[p]


def buildEntries(table):
    """(pid, name, mass, ctau, modes) for every signed id, sorted as the
    binary search needs. A state that is not its own antiparticle is written
    with its antiparticle, which carries the charge conjugate modes."""
    entries = []
    for pid, (name, mass, ctau, modes) in table.items():
        entries.append((pid, name, mass, ctau, modes))
        if not selfConjugate(pid):
            entries.append((-pid, nameOf(-pid), mass, ctau, tuple(
                (bf, tuple(d if selfConjugate(d) else -d for d in children))
                for bf, children in modes)))
    entries.sort(key=lambda entry: entry[0])
    return entries

# ------------------------------------------------------------------------------
# Write it out


def ctauLiteral(ctau):
    """The lifetime, with a sentinel for a final state and for no width."""
    if ctau is None:
        return CTAU_NO_WIDTH
    if ctau == float("inf"):
        return CTAU_STABLE
    return "%.6e" % sixFigures(ctau)


HEADER = '''#ifndef %(guard)s
#define %(guard)s

#include "base.%(h)s"

// %(rule)s
// Hadron masses [GeV], lifetimes [mm] and decay modes, keyed by
// signed pid
//
// GENERATED FILE - do not edit by hand.

// ctau of a final state, and of a state with no PDG width (K0,
// K0bar); neither is a lifetime to compare against ctau_max.
inline constexpr double ctau_stable = %(ctau_stable)s;
inline constexpr double ctau_no_width = %(ctau_no_width)s;

// The most children a mode may name
inline constexpr int max_decay_products = %(max_children)d;

// Struct to hold hadron information
struct hadron_row {
  int pid;         // signed, as the PDG numbers it
  double mass;     // GeV
  double ctau;     // mm, or one of the sentinels above
  int mode_begin;  // first row of its modes in decay_database
  int n_modes;     // how many, zero for a state that decays nowhere
};

// One child: K0/K0bar mixing into K_S/K_L, with all the momentum
struct decay_row {
  double br;                         // branching fraction, unnormalised
  int n_children;                    // 1 to max_decay_products
  int children[max_decay_products];  // their pids, the rest unused
};

inline constexpr int n_hadron_rows = %(n_hadron_rows)d;
inline constexpr int n_decay_rows = %(n_decay_rows)d;

// %(rule)s
// The tables, defined in database.%(c)s

%(qualifier)s hadron_row hadron_database[n_hadron_rows];
%(qualifier)s decay_row decay_database[n_decay_rows];

// %(rule)s
// Row of a state in hadron_database, or < 0 where the table does
// not hold it.

%(function)s int hadron_index(int pid) {
  int lo = 0;
  int hi = n_hadron_rows - 1;

  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    if (hadron_database[mid].pid == pid) return mid;
    if (hadron_database[mid].pid < pid) {
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }

  return -1;
}

// %(rule)s
// Mass of a state [GeV], or < 0 where the table does not hold
// it.

%(function)s double hadron_mass(int pid) {
  int i = hadron_index(pid);
  return (i < 0) ? -1.0 : hadron_database[i].mass;
}

#endif  // %(guard)s
'''

SOURCE = '''#include "database.%(h)s"

// %(rule)s
// GENERATED FILE - do not edit by hand.

%(qualifier)s hadron_row hadron_database[n_hadron_rows] = {
%(hadron_rows)s
};

%(qualifier)s decay_row decay_database[n_decay_rows] = {
%(decay_rows)s
};
'''


def hadronRows(masses, names):
    """One line per state, in pid order."""
    return "\n".join(
        "    {%d, %.6f, %s, %d, %d},  // %s"
        % (pid, mass, ctauLiteral(ctau), begin, count, names[pid])
        for pid, mass, ctau, begin, count in masses)


def decayRows(decays, names):
    """One line per mode, under a comment naming its parent. Slots past
    n_children are left unwritten, as clang-format has them."""
    lines, parent = [], None
    for pid, bf, children in decays:
        if pid != parent:
            parent = pid
            lines.append("    // %s (%d)" % (names[pid], pid))
        lines.append("    {%.6e, %d, {%s}},"
                     % (sixFigures(bf), len(children),
                        ", ".join(str(d) for d in children)))
    return "\n".join(lines)


def write(path, text):
    with open(path, "w") as handle:
        handle.write(text)
    return os.path.relpath(path, SRC)


def main():
    cuts = collections.Counter()
    table = buildTable(cuts)
    pruneDeadEnds(table, cuts)
    entries = buildEntries(table)
    names = {pid: name for pid, name, _, _, _ in entries}

    # --------------------------------------------------------------------------
    # Flatten the states and their modes, dropping the ones the shower cannot
    # carry

    masses = []   # (pid, mass, ctau, mode_begin, n_modes)
    decays = []   # (pid, bf, children)

    for pid, _, mass, ctau, modes in entries:
        kept = 0
        for bf, children in modes:
            if any(isForbidden(d) for d in children):
                cuts["modes naming a forbidden state"] += 1
                continue
            if len(children) > MAX_CHILDREN:
                cuts["modes with too many children"] += 1
                continue
            decays.append((pid, bf, children))
            kept += 1

        masses.append((pid, mass, ctau, len(decays) - kept, kept))

    # --------------------------------------------------------------------------
    # What the lookups rely on, checked here rather than by a static_assert so
    # that the one check covers the CUDA copy too

    # hadron_index binary searches the table, so it must be sorted by pid
    for i in range(1, len(masses)):
        if masses[i - 1][0] >= masses[i][0]:
            sys.exit("hadron_database is not sorted by pid at row %d" % i)

    # hadron_mass looks a signed pid up as given, so every child must be a row
    # and every antiparticle must have its particle
    for pid, _, children in decays:
        for child in children:
            if child not in names:
                sys.exit("%d decays to %d, which the table does not hold"
                         % (pid, child))

    for pid in names:
        if pid < 0 and -pid not in names:
            sys.exit("%d is listed without %d" % (pid, -pid))

    # --------------------------------------------------------------------------
    # Write both trees from the one set of rows

    # The same rows in both trees, which differ only in their qualifiers
    common = {
        "rule": "-" * 75,
        "ctau_stable": CTAU_STABLE,
        "ctau_no_width": CTAU_NO_WIDTH,
        "max_children": MAX_CHILDREN,
        "n_hadron_rows": len(masses),
        "n_decay_rows": len(decays),
        "hadron_rows": hadronRows(masses, names),
        "decay_rows": decayRows(decays, names),
    }

    for tree, h, c, qualifier, function in VARIANTS:
        fields = dict(common, guard="database_%s_" % h, h=h, c=c,
                      qualifier=qualifier, function=function)
        print(write(os.path.join(SRC, tree, "base", "include",
                                 "database.%s" % h), HEADER % fields))
        print(write(os.path.join(SRC, tree, "base", "src",
                                 "database.%s" % c), SOURCE % fields))

    final = sum(1 for row in masses if row[2] == float("inf"))
    print("%d states (%d final, %d with no mode), %d decay modes"
          % (len(masses), final, sum(1 for row in masses if row[4] == 0),
             len(decays)))
    for reason, n in cuts.most_common():
        print("    cut %d, %s" % (n, reason))


if __name__ == "__main__":
    main()
