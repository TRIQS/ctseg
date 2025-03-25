# Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

import sys
from code import InteractiveInterpreter


def main():
    """
    Print lines of input along with output.
    """
    source_lines = (line.rstrip() for line in sys.stdin)
    console = InteractiveInterpreter()
    source = ''
    try:
        while True:
            source = next(source_lines)
            # Allow the user to ignore specific lines of output.
            if not source.endswith('# ignore'):
                print('>>>', source)
            more = console.runsource(source)
            while more:
                next_line = next(source_lines)
                print('...', next_line)
                source += '\n' + next_line
                more = console.runsource(source)
    except StopIteration:
        if more:
            print('... ')
            more = console.runsource(source + '\n')


if __name__ == '__main__':
    main()
