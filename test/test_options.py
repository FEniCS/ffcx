# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging
import pprint
import re
from unittest import mock

import pytest

import ffcx.options


def test_options_not_formatted_when_info_disabled():
    """Option values must not be pretty-printed unless INFO is emitted.

    ``pformat`` accounts for nearly all of ``get_options``, which runs on
    every JIT cache hit, so formatting a discarded record is pure waste.
    """
    with mock.patch.object(pprint, "pformat", wraps=pprint.pformat) as pformat:
        ffcx.options.get_options({"verbosity": logging.WARNING})
    pformat.assert_not_called()


def test_options_formatted_when_info_enabled(caplog):
    """The option values are still logged when INFO is enabled."""
    with caplog.at_level(logging.INFO, logger="ffcx"):
        options = ffcx.options.get_options({"verbosity": logging.INFO})
    assert "Final option values" in caplog.text
    assert f"'verbosity': {logging.INFO}" in caplog.text
    assert options["verbosity"] == logging.INFO


def test_unknown_option_key_raises():
    """An unrecognised option key must be rejected, not silently accepted.

    An unknown key has no effect on the generated code, but it does
    change the JIT cache signature, so the form is recompiled and the
    option appears to have been applied.
    """
    with pytest.raises(ValueError, match="totally_made_up"):
        ffcx.options.get_options({"totally_made_up": "xyz"})

    # 'quadrature_degree' is a legacy FFC form compiler parameter; it is
    # set through integral metadata, not through the options
    with pytest.raises(ValueError, match="quadrature_degree"):
        ffcx.options.get_options({"quadrature_degree": 4})


def test_known_option_keys_accepted():
    """Every default option key must be accepted."""
    defaults = {k: v[1] for k, v in ffcx.options.FFCX_DEFAULT_OPTIONS.items()}
    assert ffcx.options.get_options(defaults) == defaults


def test_unknown_option_key_in_config_file_raises():
    """Keys read from a configuration file are checked too."""
    with mock.patch.object(ffcx.options, "_load_options", return_value=({"bad_key": 1}, {})):
        with pytest.raises(ValueError, match=re.escape("ffcx_options.json")):
            ffcx.options.get_options()
