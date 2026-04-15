#!/usr/bin/env python3
"""Patch the vendored SuperReads source for modern Linux toolchains."""

from __future__ import annotations

import argparse
from pathlib import Path


MARKER = "Micropeptidome SuperReads compatibility patch"


def replace_once(path: Path, old: str, new: str) -> bool:
    text = path.read_text()
    if new in text:
        return False
    if old not in text:
        raise SystemExit(f"Expected text not found in {path}: {old!r}")
    path.write_text(text.replace(old, new, 1))
    return True


def patch_configure_ac(global_dir: Path) -> bool:
    path = global_dir / "configure.ac"
    text = path.read_text()
    if f"{MARKER}: disable unused SWIG/Perl bindings" in text:
        return False

    old = """AC_SUBST([MAYBE_SWIG], [swig])
# Conf Perl
PERL_EXT_LIB='$(libdir)/perl'
AX_PERL_EXT
AM_CONDITIONAL([PERL_BINDING], [true])"""
    new = f"""AC_SUBST([MAYBE_SWIG], [])
# Conf Perl
# {MARKER}: disable unused SWIG/Perl bindings that fail with modern Perl headers.
PERL_EXT_LIB=''
:
AM_CONDITIONAL([PERL_BINDING], [false])"""

    if old not in text:
        raise SystemExit(f"Expected SWIG/Perl configure block not found in {path}")
    path.write_text(text.replace(old, new, 1))
    return True


def find_class_end(lines: list[str], class_index: int) -> int:
    depth = 0
    seen_open = False
    for index in range(class_index, len(lines)):
        depth += lines[index].count("{")
        if "{" in lines[index]:
            seen_open = True
        depth -= lines[index].count("}")
        if seen_open and depth <= 0 and lines[index].lstrip().startswith("};"):
            return index
    raise SystemExit("Could not find end of numeric_limits<__int128> specialization")


def patch_int128(global_dir: Path) -> bool:
    path = global_dir / "jellyfish" / "include" / "jellyfish" / "int128.hpp"
    text = path.read_text()
    if f"{MARKER}: libstdc++ provides __int128 numeric_limits" in text:
        return False

    lines = text.splitlines(keepends=True)
    first_class = next(
        i for i, line in enumerate(lines)
        if "class numeric_limits<__int128>" in line
    )
    guard_start = first_class
    while guard_start > 0 and lines[guard_start - 1].strip() in {"", "template<>"}:
        guard_start -= 1

    second_class = next(
        i for i, line in enumerate(lines)
        if "class numeric_limits<unsigned __int128>" in line
    )
    guard_end = find_class_end(lines, second_class) + 1

    lines.insert(
        guard_start,
        (
            f"#if !defined(__GLIBCXX_TYPE_INT_N_0)\n"
            f"// {MARKER}: libstdc++ provides __int128 numeric_limits on newer GCC.\n"
        ),
    )
    lines.insert(guard_end + 1, "#endif\n")
    path.write_text("".join(lines))
    return True


def patch_misc_hpp(global_dir: Path) -> bool:
    path = global_dir / "SuperReadsR" / "include" / "misc.hpp"
    text = path.read_text()
    if "#include <cstdint>" in text:
        return False
    return replace_once(path, "#include <signal.h>\n", "#include <signal.h>\n#include <cstdint>\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("superreads_dir", type=Path)
    args = parser.parse_args()

    global_dir = args.superreads_dir / "global-1"
    if not global_dir.is_dir():
        raise SystemExit(f"Missing SuperReads global-1 directory: {global_dir}")

    changed = [
        patch_configure_ac(global_dir),
        patch_int128(global_dir),
        patch_misc_hpp(global_dir),
    ]
    print(f"Applied SuperReads compatibility patches: {sum(changed)} file(s) changed")


if __name__ == "__main__":
    main()
