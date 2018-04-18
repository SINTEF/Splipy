#!/usr/bin/env python3

from subprocess import run, CalledProcessError
from tempfile import TemporaryDirectory
from os.path import abspath, dirname, isdir, join
from os import listdir, remove
from shutil import copytree, rmtree
import sys


if __name__ == '__main__':
    if len(sys.argv) > 1:
        remote = sys.argv[1]
    else:
        remote = input('Name of remote to push documentation to: ')

    src = abspath(dirname(__file__))

    # Ensure that local branch gh-pages exists
    try:
        run(['git', 'rev-parse', 'gh-pages'], cwd=src, check=True)
    except CalledProcessError:
        run(['git', 'branch', 'gh-pages', '{}/gh-pages'.format(remote)], check=True)

    # Ensure that local branch gh-pages is updated
    run(['git', 'fetch', remote, 'gh-pages:gh-pages'], check=True)

    with TemporaryDirectory() as tgt:

        # Clone to a temporary directory and check out the gh-pages branch
        run(['git', 'clone', src, tgt], check=True)
        run(['git', 'checkout', 'gh-pages'], cwd=tgt, check=True)

        # Copy files from doc/_build/html over
        build = join(src, 'doc', '_build', 'html')
        for c in listdir(build):
            run(['cp', '-R', join(build, c), tgt])

        # Add them all, show git status and commit
        run(['git', 'add', '-A'], cwd=tgt, check=True)
        run(['git', 'status'], cwd=tgt, check=True)
        run(['git', 'commit', '--allow-empty', '-m', 'Update documentation'], cwd=tgt, check=True)

        # Push to local repository from the temp repository
        run(['git', 'push', 'origin', 'gh-pages'], cwd=tgt, check=True)

        # Push to the remote
        run(['git', 'push', remote, 'gh-pages'], cwd=src, check=True)
