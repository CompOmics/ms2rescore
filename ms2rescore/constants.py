"""Shared constants for ms2rescore."""

import re

# Regex pattern to strip charge state suffix (e.g., "/2") from peptidoform strings
CHARGE_PATTERN = re.compile(r"(/\d+$)")
