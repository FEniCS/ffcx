# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging
import pprint
from unittest import mock

import ffcx.options


def test_options_not_formatted_when_info_disabled():
    """Option values must not be pretty-printed unless INFO is emitted.

    ``pprint.pformat`` accounts for essentially all of ``get_options``, and
    ``get_options`` runs on every JIT cache hit, so formatting it for a log
    record that is discarded is pure overhead.
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
