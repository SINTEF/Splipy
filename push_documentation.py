#!/usr/bin/env python3

from subprocess import run
from tempfile import TemporaryDirectory
from os.path import abspath, dirname, isdir, join
from os import listdir, remove
from shutil import copytree, rmtree


if __name__ == '__main__':
    src = abspath(join(dirname(__file__), '..'))
    with TemporaryDirectory() as tgt:
        run(['git', 'clone', src, tgt])
        run(['git', 'checkout', 'gh-pages'], cwd=tgt)
        for c in listdir(tgt):
            if c != '.git':
                c = join(tgt, c)
                if isdir(c):
                    rmtree(c)
                else:
                    remove(c)
        build = join(src, 'docs', 'build', 'html')
        for c in listdir(build):
            run(['cp', '-R', join(build, c), join(tgt, c)])
        run(['git', 'add', '-A'], cwd=tgt)
        run(['git', 'status'], cwd=tgt)
        run(['git', 'commit', '--allow-empty', '-m', 'Update documentation'], cwd=tgt)
        run(['git', 'push', 'origin', 'gh-pages'], cwd=tgt)
        run(['git', 'push', 'origin', 'gh-pages'], cwd=src)
