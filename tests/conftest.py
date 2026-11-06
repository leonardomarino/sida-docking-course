"""Test configuration.

In CI (SIDA_CI=1) a skipped test is a failure. Many tests here skip when the structure cache or a program is missing;
on a fresh runner that would turn the whole suite into skips and a green tick that proves nothing about the install.
The one skip that is legitimate in CI is listed below: the comparison with the receptors prepared by the Linux-only
`sida-docking` environment, which are not in the repository.
"""
import os

import pytest

ALLOWED_SKIPS_IN_CI = ("validated prepared receptors not present",)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if os.environ.get("SIDA_CI") and report.skipped:
        reason = report.longrepr[2] if isinstance(report.longrepr, tuple) else str(report.longrepr)
        if not any(allowed in reason for allowed in ALLOWED_SKIPS_IN_CI):
            report.outcome = "failed"
            report.longrepr = "skipped in CI, which would hide a broken install: %s" % reason
