#!/usr/bin/env python3
"""Guard against regex syntax that breaks Nextflow/Groovy parsing."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


def test_nf_files_do_not_use_trailing_regex_i_flag():
    def has_trailing_i_flag(text: str) -> bool:
        idx = text.find("/i")
        while idx != -1:
            next_pos = idx + 2
            if next_pos >= len(text):
                return True
            next_char = text[next_pos]
            if not (next_char.isalnum() or next_char == "_"):
                return True
            idx = text.find("/i", idx + 2)
        return False

    for nf_path in REPO_ROOT.rglob("*.nf"):
        text = nf_path.read_text(encoding="utf-8")
        assert not has_trailing_i_flag(text), (
            f"{nf_path} contains /.../i regex syntax; use inline (?i) instead"
        )
