#!/usr/bin/env python3
"""Check manuscript references, the completed LaTeX log, and PDF text geometry.

This is an automated document-integrity check, not a proof check or visual review.
The literal-command scan suits this repository's self-contained LaTeX source;
it does not expand arbitrary TeX macros. The compilation log supplies a second
check on references resolved by TeX itself.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import re
import sys
import unicodedata

import fitz


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STEM = ROOT / "papers" / "bounded_batch_discovery"


def normalized(text: str) -> str:
    return " ".join(unicodedata.normalize("NFKC", text).split())


def displayed_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def fingerprint(path: Path) -> dict:
    content = path.read_bytes()
    return {"path": displayed_path(path), "size_bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest()}


def uncomment(source: str) -> str:
    """Mask TeX comments while preserving positions, including escaped percent."""
    lines = []
    for line in source.splitlines(keepends=True):
        for index, char in enumerate(line):
            if char != "%":
                continue
            backslashes = 0
            previous = index - 1
            while previous >= 0 and line[previous] == "\\":
                backslashes += 1
                previous -= 1
            if backslashes % 2 == 0:
                line = line[:index] + "".join(
                    "\n" if c == "\n" else " " for c in line[index:])
                break
        lines.append(line)
    return "".join(lines)


def command_keys(source: str, commands: str, optional: bool = False) -> list[dict]:
    options = r"(?:\s*\[[^\]]*\])*" if optional else ""
    pattern = re.compile(
        rf"\\(?:{commands})\*?{options}\s*\{{([^{{}}]*)\}}", re.MULTILINE)
    result = []
    for match in pattern.finditer(source):
        for key in match.group(1).split(","):
            if key.strip():
                result.append({"key": key.strip(),
                               "line": source.count("\n", 0, match.start()) + 1})
    return result


def inspect_source(path: Path, errors: list[str]) -> dict:
    source = uncomment(path.read_text(encoding="utf-8"))
    labels = command_keys(source, "label")
    references = command_keys(source, "ref|eqref|pageref|autoref|cref|Cref")
    bibliography = command_keys(source, "bibitem", optional=True)
    citations = command_keys(source, "cite[a-zA-Z]*|nocite", optional=True)
    label_counts = Counter(item["key"] for item in labels)
    bibliography_counts = Counter(item["key"] for item in bibliography)
    reference_keys = {item["key"] for item in references}
    citation_keys = {item["key"] for item in citations}
    duplicates = sorted(k for k, n in label_counts.items() if n > 1)
    duplicate_bib = sorted(k for k, n in bibliography_counts.items() if n > 1)
    missing_refs = sorted(reference_keys - label_counts.keys())
    missing_cites = sorted(citation_keys - bibliography_counts.keys() - {"*"})
    if duplicates:
        errors.append("Duplicate LaTeX labels: " + ", ".join(duplicates))
    if duplicate_bib:
        errors.append("Duplicate bibliography keys: " + ", ".join(duplicate_bib))
    if missing_refs:
        errors.append("References to undefined labels: " + ", ".join(missing_refs))
    if missing_cites:
        errors.append("Citations without bibliography entries: " + ", ".join(missing_cites))
    return {
        "path": displayed_path(path),
        "scope": "Literal label/reference/citation/bibitem commands after masking comments.",
        "label_count": len(labels), "reference_count": len(references),
        "bibliography_entry_count": len(bibliography), "citation_count": len(citations),
        "duplicate_labels": duplicates, "duplicate_bibliography_keys": duplicate_bib,
        "missing_reference_targets": missing_refs, "missing_citation_targets": missing_cites,
        "unused_labels": sorted(label_counts.keys() - reference_keys),
        "uncited_bibliography_keys": ([] if "*" in citation_keys else
                                     sorted(bibliography_counts.keys() - citation_keys)),
        "labels": labels, "references": references,
        "bibliography_entries": bibliography, "citations": citations,
    }


def inspect_log(path: Path, errors: list[str]) -> dict:
    log = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "undefined_or_duplicate_references":
            r"(?:Reference|Citation) .+ undefined|There were undefined references|"
            r"multiply[- ]defined|multiply defined",
        "rerun_required": r"Rerun to get|Please rerun LaTeX|Label\(s\) may have changed",
        "tex_errors": r"^!|LaTeX Error:|Emergency stop|Fatal error occurred|"
                      r"Undefined control sequence",
        "missing_characters": r"Missing character:",
        "overfull_boxes": r"Overfull \\[hv]box",
        "underfull_boxes": r"Underfull \\[hv]box",
        "other_warnings": r"(?:LaTeX|Package|Class) .*[Ww]arning:|LaTeX Warning:",
    }
    matches = {}
    for category, pattern in patterns.items():
        expression = re.compile(pattern, re.IGNORECASE)
        matches[category] = [
            {"line": number, "text": line.strip()}
            for number, line in enumerate(log.splitlines(), start=1)
            if expression.search(line)
        ]
    for category in ("undefined_or_duplicate_references", "rerun_required",
                     "tex_errors", "missing_characters", "overfull_boxes"):
        if matches[category]:
            errors.append(f"LaTeX log contains {len(matches[category])} {category.replace('_', ' ')}.")
    output = re.search(r"Output written on\s+(.+?)\s+\((\d+)\s+pages?\b", log, re.DOTALL)
    if output is None:
        errors.append("LaTeX log has no completed PDF-output record.")
    return {
        "path": displayed_path(path),
        "completed_pdf_output": output is not None,
        "reported_page_count": int(output.group(2)) if output else None,
        "warning_policy": "Underfull boxes and other nonfatal warnings are recorded for review; "
                          "undefined references, required reruns, TeX errors, missing glyphs, "
                          "and overfull boxes fail validation.",
        "counts": {key: len(value) for key, value in matches.items()},
        **matches,
    }


def inspect_pdf(path: Path, header: str, side: float, vertical: float,
                errors: list[str]) -> dict:
    pages = []
    outliers = []
    expected_header = normalized(header)
    with fitz.open(path) as document:
        if document.is_encrypted:
            raise ValueError("The PDF is encrypted.")
        if len(document) == 0:
            errors.append("PDF contains no pages.")
        for number, page in enumerate(document, start=1):
            text = page.get_text("text")
            rect = page.rect
            header_text = normalized(page.get_text(
                "text", clip=fitz.Rect(0, 0, rect.width, 80)))
            header_present = expected_header in header_text
            if number > 1 and not header_present:
                errors.append(f"PDF page {number} lacks the expected running header.")
            if "??" in text:
                errors.append(f"PDF page {number} contains an unresolved '??' marker.")
            if "\ufffd" in text:
                errors.append(f"PDF page {number} contains Unicode replacement characters.")
            if not text.strip():
                errors.append(f"PDF page {number} contains no extractable text.")
            page_outliers = 0
            for block in page.get_text("dict")["blocks"]:
                for line in block.get("lines", []):
                    for span in line.get("spans", []):
                        if not span["text"].strip():
                            continue
                        x0, y0, x1, y1 = span["bbox"]
                        outside = []
                        if x0 < side:
                            outside.append("left")
                        if x1 > rect.width - side:
                            outside.append("right")
                        if y0 < vertical:
                            outside.append("top")
                        if y1 > rect.height - vertical:
                            outside.append("bottom")
                        if outside:
                            page_outliers += 1
                            outliers.append({"page": number,
                                             "bbox_pt": [round(v, 3) for v in span["bbox"]],
                                             "outside": outside, "text": span["text"]})
            pages.append({"page": number, "width_pt": round(rect.width, 3),
                          "height_pt": round(rect.height, 3),
                          "text_character_count": len(text),
                          "running_header_expected": number > 1,
                          "running_header_present": header_present,
                          "margin_outlier_count": page_outliers})
        metadata = document.metadata
    if outliers:
        errors.append(f"PDF has {len(outliers)} text spans outside the configured safe margins.")
    return {
        "path": displayed_path(path), "page_count": len(pages), "metadata": metadata,
        "expected_running_header": header,
        "safe_margins_pt": {"left": side, "right": side,
                            "top": vertical, "bottom": vertical},
        "pages": pages, "margin_outliers": outliers,
        "visual_review": "Not performed by this script; rendered-page review is separate.",
        "limitations": "Text extraction and bounding boxes cannot detect every rendering defect "
                       "or establish mathematical correctness.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_STEM.with_suffix(".tex"))
    parser.add_argument("--pdf", type=Path, default=DEFAULT_STEM.with_suffix(".pdf"))
    parser.add_argument("--log", type=Path, default=DEFAULT_STEM.with_suffix(".log"))
    parser.add_argument("--output", type=Path,
                        default=ROOT / "reports" / "manuscript_validation.json")
    parser.add_argument("--header-text", default="Bounded-batch species discovery")
    parser.add_argument("--side-margin", type=float, default=64.0,
                        help="Minimum left/right text margin in PDF points (default: 64).")
    parser.add_argument("--vertical-margin", type=float, default=25.0,
                        help="Minimum top/bottom text margin in PDF points (default: 25).")
    args = parser.parse_args()
    if args.side_margin < 0 or args.vertical_margin < 0:
        parser.error("Margins must be nonnegative.")
    errors: list[str] = []
    report = {
        "schema_version": 1,
        "description": "Automated manuscript integrity and text-geometry checks; "
                       "not mathematical proof or visual review.",
        "environment": {"python_version": platform.python_version(),
                        "pymupdf_version": version("PyMuPDF")},
    }
    operations = (
        ("source", args.source, lambda: inspect_source(args.source, errors)),
        ("latex_log", args.log, lambda: inspect_log(args.log, errors)),
        ("pdf", args.pdf, lambda: inspect_pdf(args.pdf, args.header_text,
                                            args.side_margin, args.vertical_margin, errors)),
    )
    for section, path, operation in operations:
        if not path.is_file():
            errors.append(f"Missing required file: {displayed_path(path)}")
            report[section] = {"path": displayed_path(path), "available": False}
            continue
        try:
            report[section] = operation()
        except (OSError, UnicodeError, ValueError, RuntimeError) as exc:
            errors.append(f"Cannot inspect {displayed_path(path)}: {exc}")
            report[section] = {"path": displayed_path(path), "inspection_failed": True}
    log_pages = report["latex_log"].get("reported_page_count")
    pdf_pages = report["pdf"].get("page_count")
    if log_pages is not None and pdf_pages is not None and log_pages != pdf_pages:
        errors.append(f"LaTeX log reports {log_pages} pages but PDF has {pdf_pages}.")
    files = [args.source, args.pdf, args.log, *sorted((ROOT / "scripts").glob("*.py")),
             ROOT / "Makefile", ROOT / "requirements.txt"]
    report["files"] = [fingerprint(path) for path in files if path.is_file()]
    report["errors"] = errors
    report["status"] = "failed" if errors else "passed"
    destination = args.output.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8")
    print(f"Manuscript validation {report['status']}: {len(errors)} error(s).")
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    print(f"Report: {destination}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
