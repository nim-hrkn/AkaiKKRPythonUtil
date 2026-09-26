# coding: utf-8
# Copyright (c) 2021-2026 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.
"""``begin_option`` block of AkaiKKR: key registry, validation, writing and reading.

The registry follows ``akaikkr_common/source/m_optn.f`` of AkaiKKR 2022.0721
(see docs/akaikkr_option_keys.md). specx reads the block token by token:
a scalar is ``key=`` followed by the value as a *separate* token (``mse= 5``;
``mse=5`` stops specx with ``unknown token``), an array is ``begin_<name>`` ...
``end_<name>``. Values are kept as strings (80 characters) inside specx and
converted where they are used; when a given key is used specx echoes one line
``optnwrt:<name> <value>`` to the output.
"""
import io
import os
import re
import warnings
from dataclasses import dataclass
from typing import Iterable, Optional, Tuple

from .Error import KKRUnknownOptionError, KKROptionValueError

OPTION_KEYS_VERSION = "2022.0721"
"""AkaiKKR version whose m_optn.f the registry follows."""

_LEN_C = 80  # character(LEN_C) in m_optn.f
_CODES = ("akaikkr", "akaikkr_cnd", "cpa2021v01")
_TRUE = ("t", "true", ".true.")
_FALSE = ("f", "false", ".false.")


@dataclass(frozen=True)
class OptionKey:
    """one key of the begin_option block.

    Attributes:
        name: canonical name (the Fortran variable name).
        aliases: other spellings accepted by specx.
        kind: "int" | "float" | "str" | "bool" | "list".
        default: value used by specx when the key is absent. None if it is derived
            (mse from ewidth) or not applicable (klabel). A dict maps code -> default
            when the default differs between builds.
        codes: builds in which the key has an effect.
        echo: True if specx writes ``optnwrt:<echo_name> <value>`` when the key is used.
        description: one line.
        echo_name: name in the echo line (the Fortran variable name of the caller;
            differs from name only for cemesr_ref, echoed as ``ref``).
    """
    name: str
    aliases: Tuple[str, ...]
    kind: str
    default: object
    codes: Tuple[str, ...]
    echo: bool
    description: str
    echo_name: str = ""

    def __post_init__(self):
        if not self.echo_name:
            object.__setattr__(self, "echo_name", self.name)

    def default_for(self, code=None):
        """default value for a code ("akaikkr", ...); code=None returns the akaikkr one."""
        if isinstance(self.default, dict):
            return self.default.get(code or "akaikkr")
        return self.default


_ALL = _CODES
_CND = ("akaikkr_cnd",)

OPTION_KEYS = {k.name: k for k in [
    OptionKey("mse", ("number_emesh",), "int", None, _ALL, True,
              "number of points of the complex energy path; default 2*int(ewidth*35/2)+1 (<=201), 201 for dos/spc/mcd"),
    OptionKey("tol", ("thresh_scf", "thresh_go"), "float", 1e-6, _ALL, True,
              "SCF convergence threshold on the potential residual cnvq"),
    OptionKey("ng", ("ndegree_cheb",), "int", 21, _ALL, True,
              "number of terms of the Chebyshev expansion of the energy integration (<= ngmx=21)"),
    OptionKey("dex", (), "float", 0.05, _ALL, True,
              "symmetry-breaking field of the first iterations (EF shifted by +-dex; also used by kick)"),
    OptionKey("rkick", (), "float", 1, _ALL, True,
              "number of kick iterations (same role as the digit of magtyp=kick3)"),
    OptionKey("critic", (), "float", -1.0, _ALL, True,
              "log10(cnvq) below which Tchebyshev mixing is replaced by simple mixing with 10*pmix (<=0.1)"),
    OptionKey("cemesr_ref", (), "float", {"akaikkr": 0.75, "akaikkr_cnd": 0.5, "cpa2021v01": 0.75}, _ALL, True,
              "fraction of the real-axis energy window (dos/spc) below EF: kef=(mse-1)*ref+1", echo_name="ref"),
    OptionKey("spmain_bnd2", (), "str", "a", _ALL, True,
              "band-energy end-point correction: a Im(e*z), b ef*Im(z), c Re(e)*Im(z)"),
    OptionKey("spckkr_dmpc0", (), "float", 1.0, _ALL, True,
              "initial damping of the CPA iteration in the spectral-function calculation"),
    OptionKey("spckkr_itrmx", (), "int", 100, _ALL, True,
              "maximum number of CPA iterations in the spectral-function calculation"),
    OptionKey("ie", ("ndirection_cnd",), "int", 3, _CND, True,
              "current direction index (1,2,3=x,y,z) of the conductivity (akaikkr_cnd; read but unused in akaikkr)"),
    OptionKey("cpaitr_show", (), "bool", False, _CND, True,
              "show the progress of the CPA iteration (akaikkr_cnd)"),
    OptionKey("cpaitr_tol", (), "float", 1e-8, _CND, True,
              "convergence threshold of the CPA iteration (akaikkr_cnd)"),
    OptionKey("ddos", (), "bool", False, _CND, True,
              "write the energy derivative of the DOS to ddos.txt (Seebeck; akaikkr_cnd)"),
    OptionKey("tempmu", (), "float", 1e-6, _CND, True,
              "temperature of the Seebeck calculation go=sbk (akaikkr_cnd)"),
    OptionKey("klabel", (), "list", None, _ALL, False,
              "array of k-point labels (begin_klabel ... end_klabel); kept and shown only"),
]}
"""canonical name -> OptionKey (AkaiKKR 2022.0721)."""

_ALIAS = {}
_ECHO_NAME = {}
for _k in OPTION_KEYS.values():
    _ALIAS[_k.name] = _k.name
    for _a in _k.aliases:
        _ALIAS[_a] = _k.name
    if _k.echo:
        _ECHO_NAME[_k.echo_name] = _k.name


def canonical_name(key):
    """canonical name of an option key; a trailing '=' is ignored.

    Raises:
        KKRUnknownOptionError: unknown key.
    """
    k = str(key).strip()
    if k.endswith("="):
        k = k[:-1]
    if k not in _ALIAS:
        raise KKRUnknownOptionError("unknown option key: {}".format(key))
    return _ALIAS[k]


def is_known_option(key):
    k = str(key).strip()
    if k.endswith("="):
        k = k[:-1]
    return k in _ALIAS


def option_key(key):
    """OptionKey of a (possibly aliased) key."""
    return OPTION_KEYS[canonical_name(key)]


def list_option_keys(code=None):
    """OptionKey list, optionally only those effective in a code ("akaikkr", "akaikkr_cnd", "cpa2021v01")."""
    if code is None:
        return list(OPTION_KEYS.values())
    if code not in _CODES:
        raise ValueError("unknown code {}, must be one of {}".format(code, _CODES))
    return [k for k in OPTION_KEYS.values() if code in k.codes]


def _to_bool(value):
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    if s in _TRUE:
        return True
    if s in _FALSE:
        return False
    raise KKROptionValueError("not a logical value: {!r}".format(value))


def _convert(name, value, kind):
    """typed value of a string/number (used by the parsers)."""
    try:
        if kind == "int":
            return int(str(value).strip())
        if kind == "float":
            return float(str(value).strip().replace("d", "e").replace("D", "E"))
        if kind == "bool":
            return _to_bool(value)
        if kind == "list":
            return [str(v) for v in value]
        return str(value).strip()
    except (TypeError, ValueError) as e:
        raise KKROptionValueError("option {}: {}".format(name, e)) from e


def _to_card_string(name, value, kind):
    """string written to the inputcard for a value."""
    if kind == "list":
        if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
            raise KKROptionValueError("option {} must be a list, got {!r}".format(name, value))
        return [str(v) for v in value]
    if kind == "bool":
        return "T" if _to_bool(value) else "F"
    if kind == "int":
        if isinstance(value, bool):
            raise KKROptionValueError("option {} must be an integer, got {!r}".format(name, value))
        try:
            iv = int(value)
        except (TypeError, ValueError) as e:
            raise KKROptionValueError("option {} must be an integer, got {!r}".format(name, value)) from e
        if isinstance(value, float) and iv != value:
            raise KKROptionValueError("option {} must be an integer, got {!r}".format(name, value))
        return str(iv)
    if kind == "float":
        if isinstance(value, bool):
            raise KKROptionValueError("option {} must be a number, got {!r}".format(name, value))
        try:
            float(value)
        except (TypeError, ValueError) as e:
            raise KKROptionValueError("option {} must be a number, got {!r}".format(name, value)) from e
        return str(value).strip()
    s = str(value).strip()
    if not s or any(c.isspace() for c in s):
        raise KKROptionValueError("option {} must be a single token, got {!r}".format(name, value))
    return s


def normalize_option(option, code=None, strict=True):
    """validate an option dict and return a new dict ready for make_inputcard.

    Keys are replaced by their canonical names, scalar values by the strings written
    to the inputcard (bool -> "T"/"F"), list values by lists of strings. The input
    dict is not modified.

    Args:
        option (dict): {key: value}; key may be an alias, with or without trailing '='.
        code (str, optional): "akaikkr", "akaikkr_cnd" or "cpa2021v01". A warning is
            issued for keys that have no effect in that code. Defaults to None.
        strict (bool, optional): raise on unknown keys. If False, unknown keys are kept
            as given with a warning. Defaults to True.

    Raises:
        KKRUnknownOptionError: unknown key (strict).
        KKROptionValueError: value not convertible to the key's kind, or longer than 80 characters.

    Returns:
        dict: canonical name -> string (or list of strings).
    """
    if option is None:
        return {}
    if code is not None and code not in _CODES:
        raise ValueError("unknown code {}, must be one of {}".format(code, _CODES))
    result = {}
    for key, value in option.items():
        kraw = str(key).strip()
        if kraw.endswith("="):
            kraw = kraw[:-1]
        if not is_known_option(kraw):
            if strict:
                raise KKRUnknownOptionError("unknown option key: {}".format(key))
            warnings.warn("unknown option key {!r} is passed to specx as is".format(key))
            if isinstance(value, (list, tuple)):
                result[kraw] = [str(v) for v in value]
            else:
                result[kraw] = str(value)
            continue
        ok = OPTION_KEYS[_ALIAS[kraw]]
        if ok.name in result:
            raise KKROptionValueError("option {} is given twice ({})".format(ok.name, key))
        if code is not None and code not in ok.codes:
            warnings.warn("option {} has no effect in {}".format(ok.name, code))
        s = _to_card_string(ok.name, value, ok.kind)
        for token in (s if isinstance(s, list) else [s]):
            if len(token) > _LEN_C:
                raise KKROptionValueError(
                    "option {}: value longer than {} characters".format(ok.name, _LEN_C))
        result[ok.name] = s
    return result


def format_option_card(option):
    """lines of the begin_option ... end_option block (empty list if option is empty).

    The format is the one AkaikkrJob.make_inputcard has always written: an empty
    line, begin_option, one space + "key= value" per scalar, " begin_key" / " v1 v2 ..." /
    " end_key" per list, end_option, an empty line.
    """
    result = []
    if not option:
        return result
    result.append("")
    result.append("begin_option")
    for key, value in option.items():
        if isinstance(value, list):
            value = list(map(str, value))
            result.append(" begin_{}".format(key))
            result.append(" " + " ".join(value))
            result.append(" end_{}".format(key))
        else:
            value = str(value)
            result.append(" " + " ".join([key + "=", value]))
    result.append("end_option")
    result.append("")
    return result


def _tokens(lines):
    """whitespace tokens of lines, '#' comment lines dropped (like xtoken of specx)."""
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        for t in s.split():
            yield t


def _parse_tokens(tokens, typed, strict):
    """tokens after begin_option -> dict; consumes up to end_option."""
    result = {}
    it = iter(tokens)
    while True:
        try:
            token = next(it)
        except StopIteration:
            raise KKRUnknownOptionError('failed to read "end_option", but found EOF')
        if token == "end_option":
            return result
        if token.startswith("begin_"):
            name = token[len("begin_"):]
            endmark = "end_" + name
            values = []
            while True:
                try:
                    t = next(it)
                except StopIteration:
                    raise KKRUnknownOptionError(
                        'failed to find the token "{}", but found EOF'.format(endmark))
                if t == endmark:
                    break
                values.append(t)
            if name in OPTION_KEYS and OPTION_KEYS[name].kind == "list":
                result[name] = values
            elif strict:
                raise KKRUnknownOptionError("unknown array label {}".format(name))
            else:
                warnings.warn("unknown array label {} (ignored by specx)".format(name))
                result[name] = values
            continue
        if token.endswith("="):
            kraw = token[:-1]
            try:
                value = next(it)
            except StopIteration:
                raise KKRUnknownOptionError(
                    'failed to read the value of {}, but found EOF'.format(token))
            if is_known_option(kraw):
                ok = OPTION_KEYS[_ALIAS[kraw]]
                result[ok.name] = _convert(ok.name, value, ok.kind) if typed else value
            elif strict:
                raise KKRUnknownOptionError("unknown token: {}".format(token))
            else:
                warnings.warn("unknown option key {} (specx stops on it)".format(token))
                result[kraw] = value
            continue
        # e.g. "mse=5": specx stops with 'unknown token'
        raise KKRUnknownOptionError("unknown token: {}".format(token))


def parse_option_block(lines, typed=True, strict=True):
    """parse a begin_option ... end_option block with the token rules of specx.

    Args:
        lines (Iterable[str] or str): lines (or one string) starting at or before begin_option.
        typed (bool, optional): convert values to int/float/bool by the key's kind.
            Defaults to True.
        strict (bool, optional): raise on unknown scalar keys and unknown array labels.
            If False they are kept (as strings) with a warning. Defaults to True.

    Raises:
        KKRUnknownOptionError: no begin_option, no end_option, "key=value" without a
            space, or an unknown key (strict).

    Returns:
        dict: canonical name -> value.
    """
    if isinstance(lines, str):
        lines = lines.splitlines()
    tokens = _tokens(lines)
    for token in tokens:
        if token == "begin_option":
            break
    else:
        raise KKRUnknownOptionError("begin_option not found")
    return _parse_tokens(tokens, typed, strict)


def _lines_of(inputcard):
    if isinstance(inputcard, str):
        if os.path.isfile(inputcard):
            with open(inputcard) as f:
                return f.read().splitlines()
        return inputcard.splitlines()
    if isinstance(inputcard, io.TextIOBase):
        inputcard.seek(0)
        return inputcard.read().splitlines()
    return list(inputcard)


def read_inputcard_option_blocks(inputcard, typed=True, strict=True):
    """all begin_option blocks of an inputcard, in order (spc inputs may have two).

    Args:
        inputcard (str, io.TextIOBase, Iterable[str]): file path, text, file handle or lines.

    Returns:
        list[dict]: one dict per block; [] if there is none.
    """
    lines = _lines_of(inputcard)
    tokens = _tokens(lines)
    blocks = []
    for token in tokens:
        if token == "begin_option":
            blocks.append(_parse_tokens(tokens, typed, strict))
    return blocks


def read_inputcard_option(inputcard, typed=True, strict=True):
    """begin_option of an inputcard as one dict; a later block overrides an earlier one
    (as specx does, since optnrd_rd writes into the same optparam). {} if there is no block.
    """
    result = {}
    for block in read_inputcard_option_blocks(inputcard, typed=typed, strict=strict):
        result.update(block)
    return result


_ECHO = re.compile(r"^\s*optnwrt:(\w+)\s+(\S+)")
_ERROR_PATTERNS = ("unknown token:", 'failed to read "end_option"', "failed to find the token")


def parse_option_echo(lines, typed=True):
    """option values echoed by specx (``optnwrt:<name> <value>``, written when the value is used).

    Only keys that were given *and* used by the run appear (e.g. cemesr_ref only in
    dos/spc runs; klabel is never echoed). The echo name of cemesr_ref is ``ref``;
    the result uses canonical names. The last echo of a key wins.

    Returns:
        dict: canonical name -> value; {} if nothing was echoed.
    """
    if isinstance(lines, str):
        lines = lines.splitlines()
    result = {}
    for line in lines:
        m = _ECHO.match(line)
        if not m:
            continue
        name, value = m.group(1), m.group(2)
        name = _ECHO_NAME.get(name, name)
        if name in OPTION_KEYS and typed:
            result[name] = _convert(name, value, OPTION_KEYS[name].kind)
        else:
            result[name] = value
    return result


def find_option_error(lines):
    """the line reporting an option error of specx (unknown token / missing end_option), or None."""
    if isinstance(lines, str):
        lines = lines.splitlines()
    for line in lines:
        s = line.strip()
        if any(p in s for p in _ERROR_PATTERNS):
            return s
    return None


__all__ = ["OPTION_KEYS", "OPTION_KEYS_VERSION", "OptionKey", "canonical_name", "is_known_option",
           "option_key", "list_option_keys", "normalize_option", "format_option_card",
           "parse_option_block", "read_inputcard_option", "read_inputcard_option_blocks",
           "parse_option_echo", "find_option_error"]
