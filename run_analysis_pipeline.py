#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Entry point script that launches the MfsFlow analysis pipeline via the CLI.

This script provides a convenient way to run the MfsFlow pipeline directly
from the command line without installing the package. It simply delegates
to the main CLI entry point.
"""

from mfsflow.cli import main


if __name__ == "__main__":
    main()
